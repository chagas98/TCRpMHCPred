
"""
Pandas DataFrame datasets for different types of input representations
"""
import pandas as pd
import os
import numpy as np
import requests
from datetime import datetime
import json
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from Bio import SeqIO
from anarci import number
import subprocess
from io import StringIO


_AHO_CDR_RANGES = {
    "CDR1": list(range(25, 41)),
    "CDR2": list(range(58, 78)),
    "CDR3": list(range(109, 138)),
    "CDR2.5": list(range(83, 89))
}

class VDJdbdataset():
    def __init__(self, path, save_dir='.', transform=None, score_threshold=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir
        self.score_threshold = score_threshold

        # Load the dataset
        #self.data = pd.read_csv(os.path.join(self.path, f'{self.split}.csv'))
        self.data = pd.read_csv(self.path, sep='\t', dtype=str)
        self.data = self.preprocess_vdjdb()
        if self.score_threshold is not None:
            self.data = self.data[self.data['Score'] >= self.score_threshold]
    
    def preprocess_vdjdb(self, species_filter: list = ['HomoSapiens']):
        """
        Preprocess VDJdb data.

        Parameters:
            score_threshold (float): Minimum score to include TCR pairs.
            species_filter (list): List of species to retain (e.g., ['HomoSapiens']).

        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()
        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')

        # Keep relevant columns
        df.columns = [x.replace(' ', '.') for x in df.columns]
        df = df[['complex.id', 'Gene', 'CDR3', 'V', 'J', 'Species', 'MHC.A', 'Epitope', 'Score', 'Reference']]

        # Drop rows without complex.id or Score
        df = df.dropna(subset=['complex.id', 'Score'])

        # Filter for complex.id appearing at least twice (paired TCRs)
        valid_ids = df['complex.id'].value_counts()
        valid_ids = valid_ids[valid_ids >= 2].index
        df = df[df['complex.id'].isin(valid_ids)]

        # Add row number within complex.id group
        df['row_num'] = df.groupby('complex.id').cumcount() + 1

        # Pivot to wide format
        pivot_cols = ['CDR3', 'V', 'J', 'Species', 'MHC.A', 'Epitope', 'Score', 'Reference']
        df_wide = df.pivot(index='complex.id', columns='row_num', values=pivot_cols)

        # Flatten MultiIndex columns
        df_wide.columns = [f"{col}_{i}" for col, i in df_wide.columns]
        df_wide.reset_index(inplace=True)

        # Filter by species
        df_wide = df_wide[
            df_wide['Species_1'].isin(species_filter) &
            df_wide['Species_2'].isin(species_filter)
        ]

        # Filter equal MHC types
        df_wide = df_wide[df_wide['MHC.A_1'] == df_wide['MHC.A_2']]
        df_wide = df_wide[df_wide['MHC.A_1'].notna() & (df_wide['MHC.A_1'] != '') & (df_wide['MHC.A_1'] != 'nan')]
        df_wide = df_wide[df_wide['MHC.A_2'].notna() & (df_wide['MHC.A_2'] != '') & (df_wide['MHC.A_2'] != 'nan')]
        df_wide['Score'] = df_wide[['Score_1', 'Score_2']].min(axis=1)
        df_wide = df_wide[df_wide['Score'].notna() & (df_wide['Score'] != '')]
        df_wide = df_wide[df_wide['Epitope_1'] == df_wide['Epitope_2']]
        df_wide = df_wide[(df_wide['Epitope_1'].str.len() >= 8) & (df_wide['Epitope_1'].str.len() <= 14)]
        
        # Reorganize
        df_out = pd.DataFrame({
            'CDR3A': df_wide['CDR3_1'],
            'CDR3B': df_wide['CDR3_2'],
            'epitope': df_wide['Epitope_1'],
            'VA': df_wide['V_1'],
            'VB': df_wide['V_2'],
            'JA': df_wide['J_1'],
            'JB': df_wide['J_2'],
            'sp': df_wide['Species_1'],
            'MHCa': df_wide['MHC.A_1'],
            'Score': df_wide['Score'],
            'Reference': df_wide['Reference_1']
        })

        # Keep only sequences that start with C and end with F
        mask = (
            df_out['CDR3A'].str.startswith('C') & df_out['CDR3A'].str.endswith('F') &
            df_out['CDR3B'].str.startswith('C') & df_out['CDR3B'].str.endswith('F')
        )
        df_out = df_out[mask]

        # Select MHC Class I
        df_out = df_out[df_out['MHCa'].str.contains(r'^HLA-[ABCDEFG]\*', na=False)]

        # Clean reference column
        df_out['Reference'] = df_out['Reference'].str.replace('PMID:', '')
        
        # Filter out empty V genes
        df_out = df_out[df_out['VA'].notna() & (df_out['VA'] != '')]

        # Remove duplicates
        df_out = df_out.drop_duplicates(subset=df_out.columns[:-2])
        df_out = df_out.dropna(subset=['CDR3A', 'CDR3B', 'VA', 'VB', 'JA', 'JB', 'MHCa', 'epitope'])

        # Normalize species names
        df_out['sp'] = df_out['sp'].replace({'HomoSapiens': 'human', 'MusMusculus': 'mouse'})

        # Normalize HLA types
        hla_map = {
            'HLA-A*02'  : 'HLA-A*02:01',
            'HLA-A*11'  : 'HLA-A*11:01',
            'HLA-B*07'  : 'HLA-B*07:01',
            'HLA-B*08'  : 'HLA-B*08:01',
            'HLA-B*27 ' : 'HLA-B*27'
        }
        df_out['MHCa'] = df_out['MHCa'].replace(hla_map)

        # Save to CSV
        df_out['TCR_ID'] = [f"VDJdb{id}" for id in range(len(df_out))]
        df_out = df_out[['TCR_ID', 'CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'sp', 'MHCa', 'Score', 'Reference']]
        cdate = datetime.now()
        
        vdjdb_file = os.path.join(self.save_dir, '01_raw', f'vdjdb_{species_filter[0]}_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv')
        df_out.to_csv(vdjdb_file, index=False)

        return df_out

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
    def __repr__(self):
        return f"VDJdbdataset(num_samples={len(self.data)})"
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data
    

class IEDBdataset():
    def __init__(self, path, save_dir='.', transform=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir

        # Load the dataset
        self.data = pd.read_csv(self.path, dtype=str)
        self.data = self.preprocess_iedb()

    def preprocess_iedb(self, species_filter: list = ['HomoSapiens']):
        """
        Preprocess IEDB data similarly to the provided R script.

        Parameters:
            species_filter (list): List of species to retain (e.g., ['HomoSapiens'])

        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()

        # Filter for alphabeta receptor type
        # Find the column matching the pattern: it starts with "Receptor" and ends with "Type"
        col_reptype = df.columns[df.columns.str.contains(r'^Receptor.*Type$')][0]
        df = df[df[col_reptype] == 'alphabeta']

        # Filter rows with same species in both chains
        col_org1 = df.columns[df.columns.str.contains(r'^Chain 1.*Organism IRI$')][0]
        col_org2 = df.columns[df.columns.str.contains(r'^Chain 2.*Organism IRI$')][0]
        df = df[df[col_org1] == df[col_org2]]

        # Replace taxon IDs with species names
        species_map = {
            "http://purl.obolibrary.org/obo/NCBITaxon_9606": "HomoSapiens",
            "http://purl.obolibrary.org/obo/NCBITaxon_10090": "MusMusculus"
        }

        df[col_org1] = df[col_org1].map(species_map).fillna('other')
        df[col_org2] = df[col_org2].map(species_map).fillna('other')

        # Filter by species
        df = df[
            (df[col_org1].isin(species_filter)) &
            (df[col_org2].isin(species_filter))
        ]

        # Fill missing calculated genes with curated genes
        for chain in ['Chain 1', 'Chain 2']:
            for gene_type in ['V Gene', 'J Gene']:
                calc_col = df.columns[df.columns.str.contains(f'^{chain}.*Calculated {gene_type}$')][0]
                curated_col = df.columns[df.columns.str.contains(f'^{chain}.*Curated {gene_type}$')][0]
                df[calc_col] = df[calc_col].fillna(df[curated_col])
                df.loc[df[calc_col] == '', calc_col] = df[curated_col]

        # Build the output DataFrame
        col_cdr3a = df.columns[df.columns.str.contains(r'^Chain 1.*CDR3 Curated$')][0]
        col_cdr3b = df.columns[df.columns.str.contains(r'^Chain 2.*CDR3 Curated$')][0]
        col_epitope = df.columns[df.columns.str.contains(r'^Epitope.*Name$')][0]
        col_va = df.columns[df.columns.str.contains(r'^Chain 1.*Calculated V Gene$')][0]
        col_vb = df.columns[df.columns.str.contains(r'^Chain 2.*Calculated V Gene$')][0]
        col_ja = df.columns[df.columns.str.contains(r'^Chain 1.*Calculated J Gene$')][0]
        col_jb = df.columns[df.columns.str.contains(r'^Chain 2.*Calculated J Gene$')][0]
        col_mhc = df.columns[df.columns.str.contains(r'^Assay.*MHC Allele Names$')][0]
        col_ref = df.columns[df.columns.str.contains(r'^Reference.*IEDB IRI$')][0]
        df = df[(df[col_epitope].str.len() >= 8) & (df[col_epitope].str.len() <= 14)]

        df_out = pd.DataFrame({
            'CDR3A': df[col_cdr3a],
            'CDR3B': df[col_cdr3b],
            'epitope': df[col_epitope],
            'VA': df[col_va],
            'VB': df[col_vb],
            'JA': df[col_ja],
            'JB': df[col_jb],
            'sp': df[col_org1],
            'MHCa': df[col_mhc],
            'Score': 0,  # Placeholder for score, not provided in IEDB
            'Reference': df[col_ref]
        })

        # Remove rows with any missing/empty string values
        df_out = df_out[(df_out != '').all(axis=1)]

        # Keep only CDR3 sequences starting with C and ending with F
        mask = (
            df_out['CDR3A'].str.startswith('C') & df_out['CDR3A'].str.endswith('F') &
            df_out['CDR3B'].str.startswith('C') & df_out['CDR3B'].str.endswith('F')
        )
        df_out = df_out[mask]

        # Convert species names
        df_out['sp'] = df_out['sp'].replace({'HomoSapiens': 'human', 'MusMusculus': 'mouse'})

        # Remove invalid or ambiguous MHC alleles
        invalid_hla = ["HLA-A*02:01, HLA-A*24:02", "HLA-B*35:01, HLA-B*35:08", "HLA-Cw3"]
        df_out = df_out[~df_out['MHCa'].isin(invalid_hla)]

        hla_corrections = {
            'HLA-B35': 'HLA-B*35:01',
            'HLA-A2': 'HLA-A*02',
            'HLA-A*02': 'HLA-A*02:01',
            'HLA-B7': 'HLA-B*07:01',
            'HLA-B8': 'HLA-B*08:01'
        }
        df_out['MHCa'] = df_out['MHCa'].replace(hla_corrections)

        # Manual VA corrections
        #TRAV names correction 
        #IMGT convention: https://www.imgt.org/textes/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=human&group=TRAV
        df_out = df_out[~df_out['VA'].isin(['TRAV36*01', 'TRAV38-2*01', 'TRDAV1*01'])]
        manual_va_map = {
            'TRAV1': 'TRAV1-1*01',
            'TRAV1*01': 'TRAV1-1*01',
            'TRAV14/DV4': 'TRAV14/DV4*01',
            'TRAV23/DV6': 'TRAV23/DV6*01',
            'TRAV29/DV5': 'TRAV29/DV5*01',
            'TRAV36/DV7': 'TRAV36/DV7*01',
            'TRAV38-2/DV8': 'TRAV38-2/DV8*01'
        }
        df_out['VA'] = df_out['VA'].replace(manual_va_map)

        # Pad missing *01 in VA
        for i in range(1, 46):
            for suffix in ['', '-1', '-2', '-3', '-4', '-5', '-6', '-7', '-8', '-9']:
                old = f'TRAV{i}{suffix}'
                new = f'{old}*01'
                df_out['VA'] = df_out['VA'].replace(old, new)

        # Pad missing *01 in VB
        for i in range(1, 31):
            for suffix in ['', '-1', '-2', '-3', '-4', '-5', '-6', '-7', '-8', '-9']:
                old = f'TRBV{i}{suffix}'
                new = f'{old}*01'
                df_out['VB'] = df_out['VB'].replace(old, new)
        for i in [3, 4, 5, 6, 7, 8, 10, 11, 12, 20, 21, 22, 23, 24, 25, 29]:
            df_out['VB'] = df_out['VB'].replace(f'TRBV{i}', f'TRBV{i}-1*01')

        # Pad missing *01 in JA
        for i in range(1, 61):
            df_out['JA'] = df_out['JA'].replace(f'TRAJ{i}', f'TRAJ{i}*01')

        # Pad missing *01 in JB
        df_out['JB'] = df_out['JB'].replace({
            'TRBJ1': 'TRBJ1-1*01',
            'TRBJ2': 'TRBJ2-1*01'
        })
        for i in range(1, 3):
            for suffix in ['-1', '-2', '-3', '-4', '-5', '-6', '-7']:
                old = f'TRBJ{i}{suffix}'
                new = f'{old}*01'
                df_out['JB'] = df_out['JB'].replace(old, new)

        # Drop duplicates
        df_out = df_out.drop_duplicates(subset=df_out.columns[:-2])
        df_out = df_out.dropna(subset=['CDR3A', 'CDR3B', 'VA', 'VB', 'JA', 'JB', 'MHCa', 'epitope'])
        
        # Keep only CDR3 sequences starting with C and ending with F
        mask = (
            df_out['CDR3A'].str.startswith('C') & df_out['CDR3A'].str.endswith('F') &
            df_out['CDR3B'].str.startswith('C') & df_out['CDR3B'].str.endswith('F')
        )
        df_out = df_out[mask]

        # Select MHC Class I
        df_out = df_out[df_out['MHCa'].str.contains(r'^HLA-[ABCDEFG]\*', na=False)] 
        
        # Save
        df_out['TCR_ID'] = [f"IEDB{id}" for id in range(len(df_out))]
        df_out = df_out[['TCR_ID', 'CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'sp', 'MHCa', 'Score', 'Reference']]
        output_path = os.path.join(self.save_dir, f'iedb_{species_filter[0]}.csv')
        df_out.to_csv(output_path, index=False)
        return df_out

    
    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
    
    def __repr__(self):
        return f"IEDBdataset(num_samples={len(self.data)})"
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data
    
class McPASTCRdataset():
    def __init__(self, path, save_dir='.', transform=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir

        # Load the dataset
        #self.data = pd.read_csv(os.path.join(self.path, f'{self.split}.csv'))
        self.data = pd.read_csv(self.path, dtype=str)
        self.data = self.preprocess_mcpas()
    
    @staticmethod
    def compute_min_score(row):
        """
        Compute the minimum score based on NGS and Single.cell columns.
        """
        
        ngs_true = str(row['NGS']).lower() == 'yes'
        singlecell_true = str(row['Single.cell']).lower() == 'yes'
        if singlecell_true:
            return 3
        elif ngs_true and not singlecell_true:
            return 1
        else:
            return 0
        
    def preprocess_mcpas(self, species_filter: list = ['Human'], singlecell: bool = True):
        """
        Preprocess McPAS-TCR data.
        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()
        
        # Keep relevant columns
        df = df[['CDR3.alpha.aa', 'CDR3.beta.aa', 'Species', 'NGS', \
                 'Single.cell', 'Epitope.peptide', \
                 'MHC', 'TRAV', 'TRAJ', 'TRBV', 'TRBJ', 'PubMed.ID']]


        df['Score'] = df.apply(self.compute_min_score, axis=1)
        df = df[df['Species'].isin(species_filter)]
        df = df[df['MHC'].notna() & (df['MHC'] != '')]
        df = df[df['CDR3.alpha.aa'].notna() & (df['CDR3.alpha.aa'] != '')]
        df = df[df['CDR3.beta.aa'].notna() & (df['CDR3.beta.aa'] != '')]
        df = df[df['TRAV'].notna() & (df['TRAV'] != '')]
        df = df[df['TRBV'].notna() & (df['TRBV'] != '')]
        df = df[df['TRAJ'].notna() & (df['TRAJ'] != '')]
        df = df[df['TRBJ'].notna() & (df['TRBJ'] != '')]
        df = df[df['MHC'].str.contains(r'^HLA-[ABCDEFG]\*', na=False)] 
        df = df[(df['Epitope.peptide'].str.len() >= 8) & (df['Epitope.peptide'].str.len() <= 14)]

        df_out = pd.DataFrame({
            'CDR3A': df['CDR3.alpha.aa'],
            'CDR3B': df['CDR3.beta.aa'],
            'epitope': df['Epitope.peptide'],
            'VA': df['TRAV'],
            'VB': df['TRBV'],
            'JA': df['TRAJ'],
            'JB': df['TRBJ'],
            'sp': df['Species'],
            'MHCa': df['MHC'],
            'Score': df['Score'],
            'Reference': df['PubMed.ID']
        })

        # Remove duplicates
        df_out.drop_duplicates(subset=df_out.columns[:-2], inplace=True)
        df_out.reset_index(drop=True, inplace=True)

        # Keep only CDR3 sequences starting with C and ending with F
        mask = (
            df_out['CDR3A'].str.startswith('C') & df_out['CDR3A'].str.endswith('F') &
            df_out['CDR3B'].str.startswith('C') & df_out['CDR3B'].str.endswith('F')
        )
        df_out = df_out[mask]
        
        # Save to CSV
        df_out['TCR_ID'] = [f"MCPAS{id}" for id in range(len(df_out))]
        df_out = df_out[['TCR_ID', 'CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'sp', 'MHCa', 'Score', 'Reference']]
        output_path = os.path.join(self.save_dir, 'mcpas_HomoSapiens.csv')
        df_out.to_csv(output_path, index=False)

        return df_out
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data
    def __len__(self):
        return len(self.data)
    def __getitem__(self, idx):     
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
    def __repr__(self):
        return f"McPASTCRdataset(num_samples={len(self.data)})"
    

class NeoTCRdataset():
    def __init__(self, path, save_dir='.', transform=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir

        # Load the dataset
        self.data = pd.read_excel(self.path, dtype=str)
        self.data = self.preprocess_neotcr()

    def preprocess_neotcr(self):
        """
        Preprocess NeoTCR data.
        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()

        # Keep relevant columns

        df = df[['Neoepitope', 'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3', 'HLA Allele', 'PubMed ID']]
        
        df = df[df['HLA Allele'].notna() & (df['HLA Allele'] != '')]
        df = df[df['TRAV'].notna() & (df['TRAV'] != '') & (df['TRAV'] != 'n.a.')]
        df = df[df['TRBV'].notna() & (df['TRBV'] != '') & (df['TRBV'] != 'n.a.')]
        df = df[df['TRA_CDR3'].notna() & (df['TRA_CDR3'] != '') & (df['TRA_CDR3'] != 'n.a.')]
        df = df[df['TRB_CDR3'].notna() & (df['TRB_CDR3'] != '') & (df['TRB_CDR3'] != 'n.a.')]
        df = df[df['Neoepitope'].notna() & (df['Neoepitope'] != '') & (df['Neoepitope'] != 'n.a.')]
        df = df[df['HLA Allele'].notna() & (df['HLA Allele'] != '') & (df['HLA Allele'] != 'n.a.') & (~df['HLA Allele'].str.contains(','))]
        df['PubMed ID'] = df['PubMed ID'].str.replace('PMID:', '')
        df = df[df['HLA Allele'].str.contains(r'^HLA-[ABCDEFG]\*', na=False)] 
        df = df[(df['Neoepitope'].str.len() >= 8) & (df['Neoepitope'].str.len() <= 14)]

        df_out = pd.DataFrame({
            'CDR3A': df['TRA_CDR3'],
            'CDR3B': df['TRB_CDR3'],
            'epitope': df['Neoepitope'],
            'VA': df['TRAV'],
            'VB': df['TRBV'],
            'JA': df['TRAJ'],
            'JB': df['TRBJ'],
            'sp': 'human',  # Assuming all NeoTCR data is from human
            'MHCa': df['HLA Allele'],
            'Score': 0,  # Placeholder for score, not provided in NeoTCR
            'Reference': df['PubMed ID']
        })

        # Remove duplicates
        df_out.drop_duplicates(subset=df_out.columns[:-2], inplace=True)
        df_out.reset_index(drop=True, inplace=True)
        
        # Keep only CDR3 sequences starting with C and ending with F
        mask = (
            df_out['CDR3A'].str.startswith('C') & df_out['CDR3A'].str.endswith('F') &
            df_out['CDR3B'].str.startswith('C') & df_out['CDR3B'].str.endswith('F')
        )
        df_out = df_out[mask]

        # Save to CSV
        df_out['TCR_ID'] = [f"NeoTCR{id}" for id in range(len(df_out))]
        df_out = df_out[['TCR_ID', 'CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'sp', 'MHCa', 'Score', 'Reference']]
        output_path = os.path.join(self.save_dir, 'neotcr_human.csv')
        df_out.to_csv(output_path, index=False)
        return df_out
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data

class PosTRAITdataset(): #TODO: verify VA/VB/JA/JB alleles, verify F/W endings
    def __init__(self, path, save_dir='.', transform=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir

        # Load the dataset
        self.data = pd.read_excel(self.path, dtype=str, sheet_name=0)
        self.data = self.preprocess_postrait()

    @staticmethod
    def vdjdb_score(row):
        """
        Compute the VDJdb score based on sequencing and identification methods.
        Parameters:
            row (pd.Series): A row from the DataFrame containing sequencing and identification methods.
        Returns:
            int: The computed score based on the methods.
        """
        method_sequencing = row.get('Sequencing_methods', '')
        method_identification = row.get('Identification_methods', '')
        has_structure_id = row.get('Structure', '') is not None and row.get('Structure', '') != 'NA' and row.get('Structure', '') != ''
        
        score = {
            'sequencing': 0,
            'identification': 0,
            'verification': 0
        }

        # Normalize input to lowercase and handle None
        sequencing = (str(method_sequencing) or '').lower()
        identification = (str(method_identification) or '').lower()

        # === 1. TCR Sequence Confidence ===
        if 'singlecell' in sequencing:
            score['sequencing'] = 3
        elif 'sanger' in sequencing:
            score['sequencing'] = 1  # frequency unknown, assume 1
        elif 'amplicon-seq' in sequencing:
            score['sequencing'] = 0  # no freq, assume low confidence
        else:
            score['sequencing'] = 1  # conservative default

        # === 2. Identification Confidence ===
        if any(kw in identification for kw in ['sort', 'stain', 'tetramer', 'pentamer', 'dextramer', 'streptamer']):
            score['identification'] = 1  # sort-based assumed
        elif 'culture' in identification or 'ctl' in identification:
            score['identification'] = 1  # culture-based assumed
        elif 'limiting' in identification:
            score['identification'] = 1  # assume some confidence
        else:
            score['identification'] = 0  # unknown or low-confidence methods

        # === 3. Verification Confidence ===
        if has_structure_id or 'pdb' in identification:
            score['verification'] = 3
        elif 'stim' in identification:
            score['verification'] = 2
        elif any(kw in identification for kw in ['stain', 'sort', 'multimer', 'dextramer', 'tetramer']):
            score['verification'] = 1
        else:
            score['verification'] = 0

        # === Override sequencing if verification was done ===
        if score['verification'] > 0:
            score['sequencing'] = 3

        score_selected = max(score.values())
        return score_selected
    
    def preprocess_postrait(self, species_filter: list = ['HomoSapiens']):
        """
        Preprocess PosTRAIT data.
        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()

        # Keep relevant columns
        df = df[['CDR3α', 'CDR3β', 'Epitope', 'TRAV', 'TRBV', 'TRAJ', 'TRBJ', \
                 'Species', 'MHC_A', 'Sequencing_methods','Identification_methods', 'PMID']]

        df = df[df['Species'].isin(species_filter)]
        df = df[df['MHC_A'].notna() & (df['MHC_A'] != '')]
        df = df[df['CDR3α'].notna() & (df['CDR3α'] != '')]
        df = df[df['CDR3β'].notna() & (df['CDR3β'] != '')]
        df = df[df['TRAV'].notna() & (df['TRAV'] != '')]
        df = df[df['TRBV'].notna() & (df['TRBV'] != '')]
        df = df[df['TRAJ'].notna() & (df['TRAJ'] != '')]
        df = df[df['TRBJ'].notna() & (df['TRBJ'] != '')]
        df = df[df['MHC_A'].str.contains(r'^HLA-[ABCDEFG]\*', na=False)]
        df['PMID'] = df['PMID'].str.replace('PMID:', '')
        df['Score'] = df.apply(lambda x: self.vdjdb_score(x), axis=1)
        df['Species'] = df['Species'].replace({'HomoSapiens': 'human', 'MusMusculus': 'mouse'})
        df = df[(df['Epitope'].str.len() >= 8) & (df['Epitope'].str.len() <= 14)]

        # Normalize HLA types
        hla_map = {
            'HLA-A*02': 'HLA-A*02:01',
            'HLA-A*11': 'HLA-A*11:01',
            'HLA-B*07': 'HLA-B*07:01',
            'HLA-B*08': 'HLA-B*08:01',
            'HLA-B*27 ': 'HLA-B*27'
        }
        df['MHC_A'] = df['MHC_A'].replace(hla_map)

        
        df_out = pd.DataFrame({
            'CDR3A': df['CDR3α'],
            'CDR3B': df['CDR3β'],
            'epitope': df['Epitope'],
            'VA': df['TRAV'],
            'VB': df['TRBV'],
            'JA': df['TRAJ'],
            'JB': df['TRBJ'],
            'sp': df['Species'],
            'MHCa': df['MHC_A'],
            'Score': df['Score'],
            'Reference': df['PMID']
        })

        # Rename columns
        df_out.columns = ['CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'sp', 'MHCa', 'Score', 'Reference']

        # Keep only CDR3 sequences starting with C and ending with F
        mask = (
            df_out['CDR3A'].str.startswith('C') & df_out['CDR3A'].str.endswith('F') &
            df_out['CDR3B'].str.startswith('C') & df_out['CDR3B'].str.endswith('F')
        )
        df_out = df_out[mask]

        # Filter out rows with missing values and drop duplicates
        df_out.dropna(subset=df_out.columns[:-2], inplace=True)
        df_out.drop_duplicates(subset=df_out.columns[:-2], inplace=True)
        df_out.reset_index(drop=True, inplace=True)
        
        # Save to CSV
        df_out['TCR_ID'] = [f"POSTRAIT{id}" for id in range(len(df_out))]
        df_out = df_out[['TCR_ID', 'CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'sp', 'MHCa', 'Score', 'Reference']]
        output_path = os.path.join(self.save_dir, 'postrait_human.csv')
        df_out.to_csv(output_path, index=False)
        
        #print(df_out['VA'].unique())
        #print(df_out['VB'].unique())
        #print(df_out['JA'].unique())
        #print(df_out['JB'].unique())
        #print(df_out['MHCa'].unique())
        return df_out
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data
    def __len__(self):
        return len(self.data)
    def __getitem__(self, idx):
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
    def __repr__(self):
        return f"PosTRAITdataset(num_samples={len(self.data)})"  

class IMMREP25dataset():
    def __init__(self, path, outfile, map_cols: dict, save_dir='.', transform=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir
        self.outfile = outfile
        self.map_cols = map_cols

        # Load the dataset
        self.data = pd.read_csv(self.path, dtype=str)
        self.data = self.preprocess_test()

    def preprocess_test(self):
        """
        Preprocess test data.
        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()
        df.rename(columns=self.map_cols, inplace=True)
        
        # Ensure all required columns are present
        required_cols = ['CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'MHCa', 'class']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")
        
        # Add HLA- prefix to all values in MHCa if not already present
        df['MHCa'] = df['MHCa'].apply(lambda x: x if str(x).startswith('HLA-') else f'HLA-{x}')
        df['TCR_ID'] = [f"TEST{id}" for id in range(len(df))]
        df = df[['TCR_ID'] + required_cols]
        df.to_csv(f'{self.save_dir}/{self.outfile}', index=False)
        return df
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data
    
    def __len__(self):
        return len(self.data)
    def __getitem__(self, idx):     
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
    def __repr__(self):
        return f"TESTdataset(num_samples={len(self.data)})"

class TESTdataset():
    def __init__(self, path, outfile, map_cols: dict, save_dir='.', transform=None):
        self.path = path
        self.transform = transform
        self.save_dir = save_dir
        self.outfile = outfile
        self.map_cols = map_cols

        # Load the dataset
        self.data = pd.read_csv(self.path, dtype=str)
        self.data = self.preprocess_test()

    def preprocess_test(self):
        """
        Preprocess test data.
        Returns:
            pd.DataFrame: Processed DataFrame.
        """
        # Load the raw data
        df = self.data.copy()
        df.rename(columns=self.map_cols, inplace=True)
        
        # Ensure all required columns are present
        required_cols = ['CDR3A', 'CDR3B', 'epitope', 'VA', 'VB', 'JA', 'JB', 'MHCa', 'class']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")
        
        # Add HLA- prefix to all values in MHCa if not already present
        df['MHCa'] = df['MHCa'].apply(lambda x: x if str(x).startswith('HLA-') else f'HLA-{x}')
        df['TCR_ID'] = [f"TEST{id}" for id in range(len(df))]
        df = df[['TCR_ID'] + required_cols]
        df.to_csv(f'{self.save_dir}/{self.outfile}', index=False)
        return df
    
    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        return self.data
    
    def __len__(self):
        return len(self.data)
    def __getitem__(self, idx):     
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
    def __repr__(self):
        return f"TESTdataset(num_samples={len(self.data)})"

class TCR3dataset():
    def __init__(self, pdb_dir, suffix, blast=True, pep_len_max=25, transform=None, scheme='aho', from_fasta=False, ref_mhc=None, selection_ids = None):
        self.pdb_dir = pdb_dir
        self.suffix = suffix
        self.transform = transform
        self.pep_len_max = pep_len_max
        self.scheme = scheme
        self.from_fasta = from_fasta
        self.ref_mhc = ref_mhc
        self.blast = blast
        self.selection_ids = [id.upper() for id in selection_ids]

        # Load the dataset
        if self.blast:
            self._get_blastdb_path()
        
        self.preprocess_tc3d()

    def _get_noncanonical_aa_(self, entity_count=None):
        # Download PDB file from RCSB PDB
        if entity_count is not None and entity_count > 0:

            for entity_id in range(1, entity_count + 1):
            
                pdb_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{self.pdb_id.upper()}/{entity_id}"
                polymer_file = requests.get(pdb_url)

                if polymer_file.status_code != 200:
                    print(f"Failed to retrieve polymer entity {entity_id} for {self.pdb_id}.pdb")
                    return False

                polymer_data = polymer_file.json()
                strand_id = polymer_data.get('entity_poly', {}).get('pdbx_strand_id', None)

                if strand_id is None:
                    print(f"Failed to retrieve Strand ID for {self.pdb_id}.pdb")
                    return False

                sequence = polymer_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code', None)

                if sequence is None:
                    print(f"Failed to retrieve sequence for {self.pdb_id}.pdb")
                    return False
                
                non_can = ['(', ')', '[', ']']
                if any(char in sequence for char in non_can):
                    return True
                
        else:
            print(f"No protein entities found in {self.pdb_id}.pdb")
            return False

    def _get_pdb_metadata(self):
        '''
        This function returns the release date (Y, M, D) and the resolution of the PDB
        '''
        try:

            # Make a request to the RCSB PDB REST API
            api_url = f'https://data.rcsb.org/rest/v1/core/entry/{self.pdb_id.upper()}'
            response = requests.get(api_url)

            # Check if the request was successful
            response.raise_for_status()
            data = response.json()

        except requests.exceptions.RequestException as e:
            data = {}
            print(f"Request failed: {e} for {self.pdb_id}")

        # Extract release date and resolution
        release_date_str = data.get('rcsb_accession_info', {}).get('initial_release_date', '')
        resolution = data.get('rcsb_entry_info', {}).get('resolution_combined', 0)
        pubmed = data.get('citation', {})[0].get('pdbx_database_id_pub_med', '')
        doi = data.get('citation', {})[0].get('pdbx_database_id_doi', '')
        reference = pubmed if pubmed else doi
        pol_ent_count = data.get('rcsb_entry_info', {}).get('polymer_entity_count_protein', 0)
        non_canonical = self._get_noncanonical_aa_(pol_ent_count)

        if release_date_str:
            release_date = datetime.strptime(release_date_str.split("T")[0], '%Y-%m-%d').strftime('%Y-%m-%d')
        else:
            release_date = 'Unknown'
        
        if isinstance(resolution, list):
            resolution = resolution[0]

        return {'release_date': release_date, 'resolution': resolution, 'non_canonical': non_canonical, 'Reference': reference}

        
    def _get_chain_sequence(self, structure, chain_id, start=None, end=None):
        try:
            chain = structure[0][chain_id]
            residues = [
                res for res in chain
                if res.id[0] == ' ' and
                   (start is None or start <= res.id[1] <= end)
            ]
            aa_seq = "".join([seq1(res.get_resname()) for res in residues])

            return aa_seq
        except KeyError:
            return None

    def _parse_structure(self):
    
        structure = self.parser.get_structure(self.pdb_id, self.pdb_file)
        
        alpha = self._get_chain_sequence(structure, "D")
        beta = self._get_chain_sequence(structure, "E")
        epitope = self._get_chain_sequence(structure, "C")
        mhc = self._get_chain_sequence(structure, "A")

        return {"TRA": alpha, "TRB": beta, "epitope": epitope, "MHCseq": mhc}
    
    def _run_anarci(self, seq, organism):
        try:

            results, chain = number(seq, scheme=self.scheme, allowed_species=organism)
            
            if not results:
                return None, None, None 
            
            seq_numbered = [(pos[0], aa) for pos, aa in results if aa != '-']
            amino_acids = "".join([seq[1] for seq in seq_numbered])

            numbers = [seq[0] for seq in seq_numbered]
            return amino_acids, numbers, chain

        except Exception as e:
            print(f"Error annotating sequence: {e}")
            return None, None, None
 
    def _get_blastdb_path(self, venv_name='blastenv'):
        """
        Get the path to the local BLAST database for the specified organism.
        """
        command = ['conda', 'run', '-n', venv_name,
                   'makeblastdb', '-in', self.ref_mhc, '-dbtype', 'prot']
        try:
            subprocess.run(command, check=True)
            print(f"BLAST database created successfully at {self.ref_mhc}")
            return self.ref_mhc
        except subprocess.CalledProcessError as e:
            print(f"Error creating BLAST database: {e}")
            return None
    
    def _get_mismatch_(self, row):
        qseq = row['qseq']
        sseq = row['sseq']
        qstart = int(row['qstart'])
        mismatches = []
        for i, (q_base, s_base) in enumerate(zip(qseq, sseq)):
            if q_base != s_base and q_base != '-' and s_base != '-':
                # mismatch position in query coordinates (1-based)
                pos = qstart + i
                mismatches.append(f'{s_base}{pos}{q_base}')

        return ','.join(mismatches)

    def _check_blast_alleles(self, df, window_size = 5):
        """ Check and parse BLAST results to identify MHC alleles."""
        
        print(df.head(20))
        df['slen']   = df['send'] - df['sstart'] + 1
        df['allele'] = df['salltitles'].apply(lambda x: x.split(' ')[1] if ' ' in x else None)
        df['first']  = df['allele'].apply(lambda x: x.split(':')[0] if x in x else None)
        df['second'] = df['allele'].apply(lambda x: x.split(':')[1] if len(x.split(':')) > 1 and ':' in x else None)
        df           = df[(df['slen'] > 100) & (df['gaps'] == 0)]

        allele     = None
        mismatches = None
        sequence   = None

        print(f"Checking gaps")
        if any(df.gaps == 0):
            print(f"Checking gaps == 0")
            df_wo_mismatch = df[(df['mismatch'] == 0)]
            if not df_wo_mismatch.empty:
                print(f"Checking mismatches")

                first_list = df_wo_mismatch['first'].dropna().unique()
                second_list = df_wo_mismatch['second'].dropna().unique()
                sequence = df_wo_mismatch['qseq'].iloc[0]
                print(f"First list: {first_list}, Second list: {second_list}")
                mismatches = None
                
                if len(second_list) == 1 and len(first_list) == 1:
                    print(f"Checking single allele")
                    allele = f'{first_list[0]}:{second_list[0]}'

                elif len(first_list) == 1 and len(second_list) < 4:
                    print(f"Checking multiple alleles")
                    alleles = [f'{first_list[0]}:{second_list[0]}' if second_list.size > 0 else f'{first_list[0]}']
                    allele = ','.join(alleles)

                elif len(first_list) == 1:
                    print(f"Checking single allele with multiple mismatches")
                    allele = first_list[0]

                elif len(first_list) > 1:
                    print(f"Checking multiple alleles with multiple mismatches")
                    allele = ','.join(first_list)

                sequence = df_wo_mismatch['qseq'].iloc[0]

            elif df_wo_mismatch.empty and any(df.pident >= 90):
                print(f"Checking pident")
                high_pident = max(df.pident)
                df_high_pident = df[(df.pident == high_pident)]

                first_list = df_high_pident['first'].dropna().unique()
                second_list = df_high_pident['second'].dropna().unique()
                sequence = df_high_pident['qseq'].iloc[0]
                print(f"High pident: {high_pident}, First list: {first_list}, Second list: {second_list}")
                
                df_high_pident.loc[:, 'mismatch'] = df_high_pident.apply(self._get_mismatch_, axis=1)

                if  len(first_list) == 1 and len(second_list) == 1:
                    print(f"Checking single allele")
                    allele = f'{first_list[0]}:{second_list[0]}' if second_list[0] is not None else f'{first_list[0]}'
                    print(allele)
                    if len(df_high_pident) > 1:
                        mismatches = df_high_pident[df_high_pident['allele'].str.contains(allele, regex=False)]['mismatch'].iloc[0]
                        sequence = df_high_pident['qseq'].iloc[0]
                    else:
                        mismatches = df_high_pident['mismatch'].iloc[0]
                        sequence = df_high_pident['qseq'].iloc[0]

                    mismatches = ','.join(list(set(mismatches.split(','))))

                elif len(first_list) == 1 and (len(second_list) < 4):
                    print(f"Checking multiple alleles (gaps != 0)")
                    alleles = [f'{first_list[0]}:{second_list[0]}' if second_list.size > 0 else f'{first_list[0]}']
                    print(second_list)
                    mismatch_list = []
                    for allele in alleles:
                        if len(df_high_pident) > 1:
                            mismatches = df_high_pident[df_high_pident['allele'].str.contains(allele, regex=False)]['mismatch'].iloc[0]
                            #mismatches = df_high_pident[df_high_pident['allele'].str.contains(allele)]['mismatch'].iloc[0]
                        else:
                            mismatches = df_high_pident['mismatch'].iloc[0]

                        mismatch_list.append(mismatches)

                    allele = ','.join(alleles)
                    mismatches = ','.join(list(set(mismatch_list)))
                
                elif len(first_list) == 1 and len(second_list) >= 4:
                    print(f"Checking single allele with multiple mismatches")
                    allele = first_list[0]
                    mismatches = df_high_pident[df_high_pident['allele'].str.contains(allele, regex=False)]['mismatch']
                    mismatches = [mismatch.split(',') for mismatch in mismatches]
                    mismatches = ','.join(list(set([item for sublist in mismatches for item in sublist])))

                elif len(first_list) > 1:
                    print(f"Checking multiple alleles with multiple mismatches")
                    allele = ','.join(first_list)
                    mismatches = 'Inf'

        print(f"Allele: {allele}, Mismatches: {mismatches}, Sequence: {sequence}")
        return {'allele': allele, 'mismatches': mismatches, 'sequence': sequence}


    def _run_localblast_(self, sequence, seq_name, venv_name='blastenv'):
        # Create a temporary FASTA file for the query sequence
        query_file = f"query_tmp.fasta"
        with open(query_file, 'w') as f:
            f.write(f">{seq_name}\n{sequence}\n")

        #Run BLAST search
        command = ['conda', 'run', '-n', venv_name,
                   'blastp', '-query', query_file, '-db', self.ref_mhc,
                   '-outfmt', '6 qseqid salltitles qstart qend sstart send slen pident nident mismatch gaps qseq sseq', 
                   '-max_target_seqs', '100']
        
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            blast_output = result.stdout.strip()
            blast_df = pd.read_csv(StringIO(blast_output), sep='\t', header=None,
                                    names=['qseqid', 'salltitles', 'qstart', 'qend', 'sstart', 'send',
                                           'slen', 'pident', 'nident', 'mismatch', 'gaps', 'qseq', 'sseq'])

            print(blast_df[['qseqid', 'salltitles', 'qstart', 'qend', 'slen', 'pident', 'nident', 'mismatch', 'gaps']].head(50))
        except subprocess.CalledProcessError as e:
            print(f"Error running BLAST search: {e}")

        finally:
            # Clean up the temporary query file
            if os.path.exists(query_file):
                os.remove(query_file)
        
        return self._check_blast_alleles(blast_df)

    def _assign_mhc_(self, sequence):
        """
        Assign MHC sequence based on the reference MHC allele.
        """
        # read the multiple fasta file with MHC alleles

        for record in SeqIO.parse(self.ref_mhc, "fasta"):
            if sequence in str(record.seq):
                return str(record.seq), record.description.split(' ')[1]
        return sequence, None  # Return the original sequence if not found

    def _select_cdr_seq(self, variable_domain_seq, numbers, cdr_name):
        cdr_ranges = _AHO_CDR_RANGES.get(cdr_name)
        if not cdr_ranges:
            raise ValueError(f"Unknown CDR name: {cdr_name}")

        start, end = cdr_ranges[0], cdr_ranges[-1]
        selected_seq = ''.join(
            aa for aa, num in zip(variable_domain_seq, numbers)
            if start <= num <= end
        )
        return selected_seq

    def _parse_fasta(self):
        response = {}
        mhc_candidates = []
        # Filter and process only relevant fasta files
        for file in (f for f in self.fasta_files if self.pdb_id.lower() in f.lower()):
            fasta_path = os.path.join(self.pdb_dir, 'sequences', file)
            record = next(SeqIO.parse(fasta_path, "fasta"))
            sequence = str(record.seq)

            # If sequence length fits epitope criteria, use it and skip further processing.
            if 6 <= len(sequence) < 25:
                print(f"Sequence length for {file} is in the range [6, 25]: {len(sequence)}")
                response['epitope'] = sequence
                continue

            # Run anarci for sequences that not have X
            if 'X' not in sequence:
                record_id_lower = record.id.lower()
                organism = "mouse" if "mus musculus" in record_id_lower else "human"
                
                variable_domain_seq, numbers, chain = self._run_anarci(sequence, organism)

                if variable_domain_seq and chain in ['A', 'B', 'D']:
                    print(f"Variable domain sequence found: {file}")
                    
                    chain = 'A' if chain in ['A', 'D'] else chain  # Treat 'D' as 'A' for TRA
                    response[f'TR{chain.upper()}'] = variable_domain_seq
                    response[f'TR{chain.upper()}_num'] = numbers
                    response[f'CDR1{chain.upper()}'] = self._select_cdr_seq(variable_domain_seq, numbers, "CDR1")
                    response[f'CDR2{chain.upper()}'] = self._select_cdr_seq(variable_domain_seq, numbers, "CDR2")
                    response[f'CDR3{chain.upper()}'] = self._select_cdr_seq(variable_domain_seq, numbers, "CDR3")

                    continue
            
            if len(sequence) > 120:
                mhc_candidates.append(sequence)

        if 'epitope' in response.keys():
            sequence = max(mhc_candidates, key=len)
            mhcseq, allele = self._assign_mhc_(sequence)
            response['MHCseq_ref'] = mhcseq
            response['MHCseq'] = sequence
            response['allele'] = allele
            print(f"Assigned MHC sequence: {mhcseq} with allele {allele}")
        else:
            struct_info = self._parse_structure()
            response['MHCseq_ref'] = struct_info['MHCseq']
            response['MHCseq']     = struct_info['MHCseq']
            response['epitope']    = struct_info['epitope']
            response['allele']     = None
            print(f"Using MHC sequence from structure.")


        # Second try to assign MHC alleles - > RUN BLAST
        if self.blast and not response.get('allele'):
            blast_results = self._run_localblast_(response['MHCseq'], f"{self.pdb_id}")

            response['allele_blast'] = blast_results.get('allele')
            response['mismatches'] = blast_results.get('mismatches')
            response['MHCseq_ref'] = blast_results.get('sequence') or response['MHCseq_ref']

        # Ensure all expected keys are present in the response.
        for key in ['TRA', 'TRB', 'epitope', 'MHCseq', 'MHCseq_ref', 'allele', 'allele_blast', 'mismatches']:
            response.setdefault(key, None)

        return response

    def preprocess_tc3d(self):
        self.parser = PDBParser(QUIET=True)
        self.rows = []
        self.fasta_files = os.listdir(os.path.join(self.pdb_dir, 'sequences'))

        for pdb_file in os.listdir(self.pdb_dir):
            if pdb_file.endswith(self.suffix):
                self.pdb_id = pdb_file[:-len(self.suffix)]

                if self.selection_ids is not None and self.pdb_id.upper() not in self.selection_ids:
                    print(f"Skipping {self.pdb_id} as it is not in the selection list.")
                    continue
                
                self.pdb_file = os.path.join(self.pdb_dir, pdb_file)

                # If from_fasta is True, use the fasta data, otherwise use the PDB data
                if self.from_fasta:
                    structure_data = self._parse_fasta()
                else:
                    structure_data = self._parse_structure()

                # Get metadata from PDB: resolution and release date
                metadata = self._get_pdb_metadata()

                self.rows.append({
                    "TCR_ID"         : f"PDB{self.pdb_id}",
                    "TRA"            : structure_data["TRA"],
                    "TRB"            : structure_data["TRB"],
                    "CDR1A"          : structure_data.get("CDR1A", None),
                    "CDR2A"          : structure_data.get("CDR2A", None),
                    "CDR3A"          : structure_data.get("CDR3A", None),
                    "CDR1B"          : structure_data.get("CDR1B", None),
                    "CDR2B"          : structure_data.get("CDR2B", None),
                    "CDR3B"          : structure_data.get("CDR3B", None),
                    "TRA_num"        : structure_data.get("TRA_num", None),
                    "TRB_num"        : structure_data.get("TRB_num", None),
                    "peptide"        : structure_data["epitope"],
                    "MHCseq"         : structure_data["MHCseq"],
                    "MHCseq_ref"     : structure_data["MHCseq_ref"],
                    "allele"         : structure_data["allele"],
                    "allele_blast"   : structure_data["allele_blast"],
                    "mismatches"     : structure_data["mismatches"],
                    "Score"          : 3,
                    "Release_date"   : metadata['release_date'],
                    "Resolution"     : round(metadata['resolution'], 2),
                    "NonCanonical"   : metadata['non_canonical'],
                    "PDB_ID"         : self.pdb_id,
                    "Reference"      : metadata['Reference']
                })


        print(f"Total PDB files processed: {len(self.rows)}")
        print(self.rows)

    def to_df(self):
        """
        Convert the dataset to a pandas DataFrame.
        """
        self.data = pd.DataFrame(self.rows)
        return self.data

class ConcatDataset():
    def __init__(self, datasets: list, transform=None, deduplicate_epitope=None, labels: list = None):
        self.datasets = datasets
        self.transform = transform
        self.deduplicate = deduplicate_epitope
        self.labels = labels

        if self.deduplicate is not None:
            self.deduplicate_epitope()

    def deduplicate_epitope(self):
        """
        Remove duplicate entries based on peptide column.
        """

        index_guide = next((i for i, dataset in enumerate(self.datasets) if dataset.to_df().equals(self.deduplicate.to_df())), None)
        self.datasets_filtered = [None] * len(self.datasets)
        ref = self.deduplicate.to_df()

        for i, dataset in enumerate(self.datasets):

            # If the dataset is the guide, keep it as is
            if i == index_guide:

                self.datasets_filtered[i] = dataset.to_df().copy()
                continue
            
            # Ensure that df and ref have different epitopes
            df = dataset.to_df()
            mask = df['epitope'].isin(ref['epitope'])
            df = df[~mask]
            self.datasets_filtered[i] = df.copy()
    
    def add_dataset(self, dataset, label=None):
        """
        Add a new dataset to the concatenated dataset.
        """
        dataset_df = dataset.to_df()
        self.datasets_filtered.append(dataset_df)
        self.labels.append(label)

    def concat(self):
        """
        Concatenate all datasets into a single DataFrame.
        """
        if self.deduplicate is None:
            self.datasets_filtered = [dataset.to_df() for dataset in self.datasets]
            
        for dataset, label in zip(self.datasets_filtered, self.labels):
            dataset['set_name'] = label

        self.data = pd.concat(self.datasets_filtered, ignore_index=True)
        
        # Reset index
        self.data.reset_index(drop=True, inplace=True)

        return self.data
    
    def to_df(self):
        """
        Convert the concatenated dataset to a pandas DataFrame.
        """
        if not hasattr(self, 'data'):
            self.data = self.concat()
        
        return self.data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        item = self.data.iloc[idx]
        if self.transform:
            item = self.transform(item)
        return item

    def __repr__(self):
        return f"ConcatDataset(num_samples={len(self.data)})"