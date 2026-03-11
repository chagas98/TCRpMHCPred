import pandas as pd
import numpy as np
import json
import random
from itertools import product
import Levenshtein
from datetime import datetime


class TCRNegativeSampler:
    def __init__(self, input_path, cdr_ranges_path, proportion, readable_df_path=None):
        self.input_path = input_path
        self.proportion = proportion

        with open(cdr_ranges_path, 'r') as f:
            self.cdr_ranges = json.load(f)

        self.df = pd.read_csv(self.input_path, dtype=str)
        self.TCRab_cols = [col for col in self.df.columns if col.endswith("a") or col.endswith("b")]
        self.pMHC_cols = [col for col in self.df.columns if col.endswith("p") or col.endswith("mhc")]
        self.pep_cols = [col for col in self.df.columns if col.endswith("pep")]
        self.TCR_ID_col = [col for col in self.df.columns if col.startswith("TCR_ID")]
        self.readable_cols = []

        if readable_df_path is not None:
            print(f"Adding readable DataFrame from {readable_df_path}")
            self._add_readable_df(readable_df_path)
        
        self._prepare_data()

    def _add_readable_df(self, readable_df_path):
        readable_df = pd.read_csv(readable_df_path, dtype=str)
        if 'TCR_ID' in readable_df.columns:
            self.df = self.df.merge(readable_df, on='TCR_ID', how='left')
            self.pmhc_readable_cols = ['peptide', 'MHCseq']
            self.tcr_readable_cols = ['TRA', 'TRB']
            self.extra_readable_cols = ['assigned_allele', 'MHCa']
        else:
            raise ValueError("Readable DataFrame must contain 'TCR_ID' column for merging.")
    def _prepare_data(self):
        start_a, end_a = self.cdr_ranges['alpha']['CDR3a']
        start_b, end_b = self.cdr_ranges['beta']['CDR3b']
        # Get the column names for CDR3a and CDR3b based on their positions
        cdr3a_cols = self.df.columns[self.df.columns.get_loc(f'{start_a}a'):self.df.columns.get_loc(f'{end_a}a') + 1].tolist()
        cdr3b_cols = self.df.columns[self.df.columns.get_loc(f'{start_b}b'):self.df.columns.get_loc(f'{end_b}b') + 1].tolist()

        print(f"CDR3a columns: {cdr3a_cols}")
        print(f"CDR3b columns: {cdr3b_cols}")

        self.df['TCR'] = self.df[self.TCRab_cols].agg(''.join, axis=1)
        self.df['PeptideMHC'] = self.df[self.pMHC_cols].agg(''.join, axis=1)
        self.df['TCR_partial'] = self.df[cdr3a_cols + cdr3b_cols].agg(''.join, axis=1)
        self.df['pep_partial'] = self.df[self.pep_cols].agg(''.join, axis=1)

        pep_order = self.df['pep_partial'].value_counts().index
        self.df['pep_partial'] = pd.Categorical(self.df['pep_partial'], categories=pep_order, ordered=True)
        self.df = self.df.sort_values('pep_partial').reset_index(drop=True)

    def _precompute_peptide_distances(self):
        unique_peptides = self.df['pep_partial'].unique()
        pep_pairing = pd.DataFrame(product(unique_peptides, repeat=2), columns=['col1', 'col2'])
        pep_pairing['ld'] = pep_pairing.apply(lambda row: Levenshtein.distance(row['col1'], row['col2']), axis=1)
        return pep_pairing

    def generate_negatives(self):
        pep_pairing = self._precompute_peptide_distances()
        available_partial_tcrs = list(np.repeat(self.df['TCR_partial'].values, self.proportion))
        list_new_partial_tcrs = []

        random.seed(42)  # For reproducibility
        np.random.seed(42)  # For reproducibility
        for i, row in self.df.iterrows():
            current_pep = row['pep_partial']
            current_pep_complete = row['peptide']
            current_tcr = row['TCR_partial']
            valid_peps = pep_pairing.query("col1 == @current_pep and ld >= 3")['col2'].unique()
            possible_rows = self.df[self.df['pep_partial'].isin(valid_peps)]
            possible_tcrs = possible_rows[possible_rows['TCR_partial'].isin(available_partial_tcrs)]['TCR_partial'].tolist()

            new_partial_tcr_list = []
            pool = possible_tcrs if len(possible_tcrs) >= self.proportion else possible_rows['TCR_partial'].tolist()
            
            while len(new_partial_tcr_list) < self.proportion:
                sampled_tcr = np.random.choice(pool)
                if Levenshtein.distance(current_tcr, sampled_tcr) >= 3 and sampled_tcr not in new_partial_tcr_list:
                    new_partial_tcr_list.append(sampled_tcr)
                    if sampled_tcr in available_partial_tcrs:
                        available_partial_tcrs.remove(sampled_tcr)

            list_new_partial_tcrs.append(new_partial_tcr_list)

        return list_new_partial_tcrs

    def build_dataset(self, list_new_partial_tcrs):
        flat_tcrs = [tcr for sublist in list_new_partial_tcrs for tcr in sublist]
        new_tcr_rows = []

        if len(self.pmhc_readable_cols) > 0:
            final_cols = self.TCR_ID_col + self.TCRab_cols + self.tcr_readable_cols
            pmhc_data_readable = self.df[self.pmhc_readable_cols + self.extra_readable_cols].loc[self.df.index.repeat(self.proportion)].reset_index(drop=True)
        else:
            pmhc_data_readable = None

        for tcr in flat_tcrs:
            match_row = self.df[self.df['TCR_partial'] == tcr].iloc[0]
            new_tcr_rows.append(match_row[final_cols].values)

        new_dt_tcr = pd.DataFrame(new_tcr_rows, columns=final_cols)
        pmhc_data_x = self.df[self.pMHC_cols].loc[self.df.index.repeat(self.proportion)].reset_index(drop=True)
        new_dt_full_neg = pd.concat([new_dt_tcr.reset_index(drop=True), pmhc_data_x, pmhc_data_readable], axis=1)
        # Number the _neg for each repeat of TCR_ID
        new_dt_full_neg[self.TCR_ID_col[0]] = (
            new_dt_full_neg.groupby(self.TCR_ID_col[0]).cumcount()
            .add(1)
            .astype(str)
            .radd(new_dt_full_neg[self.TCR_ID_col[0]] + '_neg')
        )
        new_dt_full_neg['class'] = 0

        positive_dt = self.df[final_cols + self.pMHC_cols + self.pmhc_readable_cols + self.extra_readable_cols].copy()
        positive_dt['class'] = 1

        full_dataset = pd.concat([positive_dt, new_dt_full_neg], ignore_index=True)

        if len(self.pmhc_readable_cols) > 0:
            save_cols =  self.tcr_readable_cols + self.pmhc_readable_cols + self.extra_readable_cols
            full_dataset_readable = full_dataset[['TCR_ID'] + save_cols + ['class']].copy()
            full_dataset.drop(columns=save_cols, inplace=True)

        cdate = datetime.now()
        full_dataset.to_csv(f'featured_dataset_train_NEG_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv', index=False)
        full_dataset_readable.to_csv(f'readable_dataset_train_NEG_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv', index=False)

