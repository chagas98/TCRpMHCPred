from Stitchr import stitchrfunctions as fxn
from multiprocessing import Pool

import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
import re
from os import cpu_count
import os
import timeit

class RunThimble:
    def __init__(self, input_path, output_path='combined_dataset_fullTCR.csv', tmp_path='tmp', thimble_path='thimble', retries=False):

        self.input_path = input_path
        self.input_df = pd.read_csv(input_path)
        self.output_path = output_path
        self.tmp_out_files = tmp_path
        self.thimble_path = thimble_path
        self.retries = retries
        
        os.makedirs(self.tmp_out_files, exist_ok=True)
        
    def format_for_thimble(self, input_path, output_path='thimble.tsv', missing_strategy='remove', chain='A'):
        
        df = pd.read_csv(input_path)
        
        if 'TCR_ID' not in df.columns:
            raise ValueError("Input DataFrame must contain 'TCR_ID' column.")
        
        rename_cols = {
            'TCR_ID': 'TCR_name',
            'CDR3A': 'TRA_CDR3',
            'CDR3B': 'TRB_CDR3',
            'epitope': 'peptide',
            'VA': 'TRAV',
            'VB': 'TRBV',
            'JA': 'TRAJ',
            'JB': 'TRBJ',
            'sp': 'sp',
            'MHCa': 'MHCa',
            'Score': 'Score',
            'Reference': 'Reference'
        }

        df.rename(columns=rename_cols, inplace=True)

        thimble_columns = ['TCR_name', 'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3', 'TRAC', 'TRBC', 'TRA_leader',
                        'TRB_leader',
                        'Linker', 'Link_order', 'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq',
                        'TRB_3_prime_seq']
        
        for column in thimble_columns:
            if column not in df.columns:
                df[column] = np.nan
        for column in ['TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3']:
            if column[2] != chain:
                df[column] = np.nan
            elif missing_strategy != 'remove':
                df[column] = df[column].map(lambda x: x.replace('unknown', '%'))
            else:
                df = df[df[column] != 'unknown']
        
        df = df[thimble_columns]
        print(f'colums for thimble: {df.columns}')
        print(f'rows for thimble: {df.shape[0]}')

        df.to_csv(output_path, sep='\t', index=False)

    def run_thimble_schain(self, thimble_input_path: str, output_path: str, chain: str, 
                           retry_F: bool = True, retry_CF: bool = True, retry_W: bool = True):
        """ 
        Run Thimble for a specific chain (A or B).
        Args:
            thimble_input_path (str): Path to the formatted Thimble input file.
            output_path (str): Path to save the Thimble output file.
            chain (str): Chain to run Thimble for ('A' or 'B').
            retry (bool): Whether to retry Thimble with terminal Ws for failed CDR3s.
        Returns:
            dict: A dictionary mapping TCR names to their amino acid sequences for the specified chain.
            dict: A dictionary mapping TCR names to their new CDR3 sequences if retry is True.
        """
        
        subprocess.call(f'{self.thimble_path} -in {thimble_input_path} -o {output_path} -s human', shell=True)

        outputs = pd.read_csv(output_path, sep='\t')
        fixed = outputs[outputs[f'TR{chain}_aa'].notna()]
        output = dict(fixed[['TCR_name', f'TR{chain}_aa']].values)
        
        if retry_F:
            print(f'Successfully ran Thimble for chain {chain} for {fixed.shape[0]} sequences. Running With terminal Fs for failures.')
            retry_path = thimble_input_path.replace(".tsv", "_redo_F.tsv")
            input_file = pd.read_csv(thimble_input_path, sep='\t')
            to_redo = set(input_file['TCR_name'].unique()).difference(set(fixed['TCR_name'].unique()))
            to_redo_df = input_file[input_file['TCR_name'].isin(to_redo)].copy()
            to_redo_df[f'TR{chain}_CDR3'] = to_redo_df[f'TR{chain}_CDR3'].map(lambda x: x + 'F')
            to_redo_df.to_csv(retry_path, sep='\t', index=False)

            output_path = output_path.replace(".tsv", "_F.tsv")
            subprocess.call(f'{self.thimble_path} -in {retry_path} -o {output_path} -s human', shell=True)
            outputs_F = pd.read_csv(output_path, sep='\t')
            output.update(dict(outputs_F[outputs_F[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_aa']].values))
            print(f'Found {outputs_F[f"TR{chain}_aa"].notna().sum()} sequences with terminal F.')

        if retry_CF:
            outputs = pd.read_csv(output_path, sep='\t')
            fixed = outputs[outputs[f'TR{chain}_aa'].notna()]
            print(f'Successfully ran Thimble for chain {chain} for {fixed.shape[0]} sequences. Running with start C and end Fs for failures.')
            
            retry_path = thimble_input_path.replace(".tsv", "_redo_CF.tsv")
            input_file = pd.read_csv(thimble_input_path, sep='\t')
            to_redo = set(input_file['TCR_name'].unique()).difference(set(fixed['TCR_name'].unique()))
            to_redo_df = input_file[input_file['TCR_name'].isin(to_redo)].copy()
            to_redo_df[f'TR{chain}_CDR3'] = to_redo_df[f'TR{chain}_CDR3'].map(lambda x: 'C' + x + 'F')
            to_redo_df.to_csv(retry_path, sep='\t', index=False)

            output_path = output_path.replace(".tsv", "_CF.tsv")
            subprocess.call(f'{self.thimble_path} -in {retry_path} -o {output_path} -s human', shell=True)
            outputs = pd.read_csv(output_path, sep='\t')
            output.update(dict(outputs[outputs[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_aa']].values))
            new_cdr3bs = dict(outputs[outputs[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_CDR3']].values)
            print(f'Found {outputs[f"TR{chain}_aa"].notna().sum()} sequences with start C and terminal F.')

        if retry_W:

            outputs = pd.read_csv(output_path, sep='\t')
            fixed = outputs[outputs[f'TR{chain}_aa'].notna()]
            print(f'Successfully ran Thimble for chain {chain} for {fixed.shape[0]} sequences. Running With terminal Ws for failures.')
            
            retry_path = thimble_input_path.replace(".tsv", "_redo_W.tsv")
            input_file = pd.read_csv(thimble_input_path, sep='\t')
            to_redo = set(input_file['TCR_name'].unique()).difference(set(fixed['TCR_name'].unique()))
            to_redo_df = input_file[input_file['TCR_name'].isin(to_redo)].copy()
            to_redo_df[f'TR{chain}_CDR3'] = to_redo_df[f'TR{chain}_CDR3'].map(lambda x: x + 'W')
            to_redo_df.to_csv(retry_path, sep='\t', index=False)

            output_path = output_path.replace(".tsv", "_W.tsv")
            subprocess.call(f'{self.thimble_path} -in {retry_path} -o {output_path} -s human', shell=True)
            outputs = pd.read_csv(output_path, sep='\t')
            output.update(dict(outputs[outputs[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_aa']].values))
            new_cdr3bs = dict(outputs[outputs[f'TR{chain}_aa'].notna()][['TCR_name', f'TR{chain}_CDR3']].values)
            print(f'Found {outputs[f"TR{chain}_aa"].notna().sum()} sequences with terminal W.')

        else:
            new_cdr3bs = {}
        
        return output, new_cdr3bs
    
    
    def run(self):
        """ 
        Run Thimble for both chains A and B.
        """
        thimble_outputs = {}
        for chain in ['A', 'B']:
            
            print(f'Formatting chain {chain} data.')
            
            thimble_format_path = os.path.join(self.tmp_out_files, f'thimble_chain{chain}.tsv')

            self.format_for_thimble(input_path       = self.input_path,
                                    output_path      = thimble_format_path, 
                                    missing_strategy = 'remove', 
                                    chain            = chain)
            
            print(f'Running Thimble for chain {chain}.')
            
            thimble_results_path = os.path.join(self.tmp_out_files, f'thimble_chain{chain}_results.tsv')

            if self.retries:
                retry_F = True
                retry_CF = True
                retry_W = True
            else:
                retry_F = False
                retry_CF = False
                retry_W = False

            thimble_outputs[chain] = self.run_thimble_schain(thimble_input_path = thimble_format_path,
                                                             output_path        = thimble_results_path,
                                                             chain              = chain,
                                                             retry_F            = retry_F,
                                                             retry_CF           = retry_CF,
                                                             retry_W            = retry_W)
            
            # create a new column in the input DataFrame for the TCR sequences
            self.input_df[f'TR{chain}'] = self.input_df['TCR_ID'].map(thimble_outputs[chain][0])

            # exclude NA sequences
            self.input_df = self.input_df[(self.input_df[f'TR{chain}'].notna()) & (self.input_df[f'TR{chain}'] != '')]

            print(f'LENTGH VALUE: Chain {chain} processed. Found {self.input_df[f"TR{chain}"].notna().sum()} sequences.')
        return  self.input_df
    
    @property
    def result(self):
        """
        Return the DataFrame with TCR sequences for both chains A and B.
        """
        if 'TRB' not in self.input_df.columns or 'TRA' not in self.input_df.columns:
            raise ValueError("Run Thimble first to get TCR sequences.")
        return self.input_df


class MHCSeqAssigner:
    """ Class to assign MHC sequences to TCRs based on their alleles.
    It reads a FASTA file containing MHC sequences and matches them to TCR alleles.
    The FASTA file should contain sequences in the format 'HLA-A*01:01:01:01'.
    The class provides methods to parse the FASTA file, expand alleles, find closest alleles, and assign sequences to a DataFrame of TCR samples.
    """

    def __init__(self, fasta_path, tmp_path='tmp/'):
        self.fasta_df = self._parse_fasta(fasta_path)
        self.tmp_path = tmp_path
        
    def _parse_fasta(self, fasta_file):
        records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            match = re.search(r'\w+\*\d{1,2}(?::\d{1,3}){1,3}', record.description)
            allele = match.group(0) if match else None

            if allele:
                records.append({
                    'Full_Allele': allele,
                    'mhc_seq': str(record.seq)
                })

        df = pd.DataFrame(records).dropna()
        return df

    def _expand_allele(self, allele, max_fields=4):
        parts = allele.split(":")
        expansions = []
        for i in range(max_fields - len(parts) + 1):
            new_fields = parts + ["01"] * i
            expansions.append(":".join(new_fields))
        return expansions

    def _find_closest_allele(self, allele, fasta_dict):
        parts = allele.split(":")
        expansions = [":".join(parts + ["01"] * i) for i in range(4 - len(parts) + 1)]

        for exp_a in expansions:
            pattern = re.compile("^" + re.escape(exp_a) + r"(?:[:].*)?$")
            for full_allele, seq in fasta_dict.items():
                if pattern.match(full_allele):
                    return {"Core_Allele": allele, "Full_Allele": full_allele, "mhc_seq": seq}
        print(f"No match found for allele: {allele}")
        return {"Core_Allele": allele, "Full_Allele": None, "mhc_seq": None}


    def multiprocess_finder(self, alleles):
        """
        Find closest alleles for a list of alleles using multiprocessing.
        """
        num_workers = max(1, int(0.8 * cpu_count()))
        if num_workers > 25:
            num_workers = 25
        print(f"Using {num_workers} CPU cores for multiprocessing.")

        fasta_dict = dict(zip(self.fasta_df["Full_Allele"], self.fasta_df["mhc_seq"]))

        with Pool(num_workers) as pool:
            # Passa tuplas de (allele, fasta_dict) para cada item individual
            results_nested = pool.starmap(
                self._find_closest_allele,
                [(allele, fasta_dict) for allele in alleles])
            
        return results_nested
    
    def assign_sequences(self, input_df: pd.DataFrame, description="all", keep_cols=None):
        """ Assign MHC sequences to TCR samples based on their alleles.
        Args:
            samples_df (pd.DataFrame): DataFrame containing TCR samples with a column 'MHCa' for MHC alleles.
            description (str): Description of the dataset, used to determine the output columns.
        Returns:
            pd.DataFrame: DataFrame with TCR samples and their assigned MHC sequences.
        """
        
        df = input_df.copy()
        df["Core_Allele"] = df["MHCa"].str.replace("HLA-", "", regex=False)
        
        # Find closest alleles using multiprocessing
        matches = self.multiprocess_finder(df['Core_Allele'].unique())
        match_df = pd.DataFrame(matches)

        # Create a mapping from Core_Allele to Full_Allele and mhc_seq
        df['Full_Allele'] = df['Core_Allele'].map(match_df.set_index('Core_Allele')['Full_Allele'])
        df['mhc_seq'] = df['Core_Allele'].map(match_df.set_index('Core_Allele')['mhc_seq']) 

        # Drop the helper column after merging
        df.drop(columns=["Core_Allele"], inplace=True)

        # Drop none or empty sequences
        df = df[df['mhc_seq'].notna() & (df['mhc_seq'] != '')]


        # Base columns
        base_cols = ["TCR_ID", "TRA", "TRB", "epitope", "MHCa", "Full_Allele", "mhc_seq", "class", "Score", "Reference"]
        if keep_cols:
            base_cols += keep_cols

        # Keep only columns present in df
        present_cols = [col for col in base_cols if col in df.columns]
        df = df[present_cols]

        print(f"LENGTH VALUE: Assigned MHC sequences to {df.shape[0]} samples.")

        # Define renaming based on description
        rename_all = {
            "TCR_ID": "TCR_ID", "TRA": "TRA", "TRB": "TRB", "epitope": "peptide",
            "MHCa": "MHCa", "Full_Allele": "assigned_allele", "mhc_seq": "MHCseq",
            "Score": "Score", "Reference": "Reference"
        }

        rename_test = rename_all | {"class": "class"}

        rename_map = rename_all if description == "train" else rename_test

        # Only rename columns that exist in df
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        self.result_df = df.copy()
        
        return self.result_df

    @property
    def result(self):
        """
        Return the DataFrame with MHC sequences.
        """
        if 'MHCseq' not in self.result_df.columns:
            raise ValueError("Run MHCSeqAssigner first to get MHC sequences.")
        return self.result_df