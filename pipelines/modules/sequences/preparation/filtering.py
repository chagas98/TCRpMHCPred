#!/usr/bin/env python3

import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import argparse
import pickle
import json
import timeit
import ast

class Filtering:
    def __init__(self, 
                 input_df: pd.DataFrame,
                 selected_factors: str | None = None, 
                 exclude_contacts: str = None, 
                 max_pep_length: int = 14, 
                 just10x: bool = None,
                 filter_rows: bool = True,
                 tmp_path: str = 'tmp/', 
                 description: str = 'train'):
        """
        Initialize the Filtering class with parameters for filtering data.
        :param input_df: DataFrame containing the input data.
        :param exclude_contacts: Dictionary containing positions to exclude for peptide, TCR alpha, TCR beta, MHC, and alleles.
        :param max_pep_length: Maximum length of peptide sequences to include.
        :param just10x: Boolean indicating whether to filter for just 10x data.
        :param description: Description of the dataset.
        :param outfile: Output file path for saving the filtered data.
        """

        self.input_df = input_df
        
        with open(exclude_contacts) as f:
            self.exclude_contacts = json.load(f)
        
        self.max_pep_length   = max_pep_length
        self.just10x          = just10x
        self.filter_rows      = filter_rows
        self.tmp_path         = tmp_path
        self.description      = description
        
        # read a list from selected_factors.txt
        if selected_factors is not None:
            with open(selected_factors) as f:
                self.selected_factors = f.read().splitlines()
        else:
            self.selected_factors = None

        print('PARAMETERS DATA')
        print('------------------------------------------------------------------------------')
        print('Positions to exclude:', self.exclude_contacts)     
        print('Max peptide length:',   self.max_pep_length)
        print('Selected factors:',     self.selected_factors)
        print('Filter rows:',          self.filter_rows)
        print('Just 10x:',             self.just10x)
        print('Description:',          self.description)
        print('------------------------------------------------------------------------------\n')

    def filtering_just10x(self, data, col='just10x'):
        # Filter rows where the specified column contains True and drop the column
        filtered_data = data.loc[data[col]]
        filtered_data = filtered_data.drop(columns=[col], errors='ignore')
        return filtered_data

    def filtering_PepLen(self, data, col, len_max):

        pep_columns = [pep for pep in data.columns if pep.endswith('pep')]
        excl_columns = []

        if isinstance(len_max, int):
            data = data[data[col].str.len() <= len_max]
            col_name_include = [f'{n}pep' for n in range(1,len_max+1)]
            excl_columns = [col for col in pep_columns if col not in col_name_include]
            print(excl_columns)
        
        if len(excl_columns) > 0:
            print(excl_columns)
            data.drop(columns=excl_columns, inplace=True)

        return data

    def filtering_positions(self, data: pd.DataFrame, 
                            selected_columns: list | None,
                            suffix: str | None) -> pd.DataFrame:
        
        if len(selected_columns) == 0 or selected_columns is None:
            return data

        if not isinstance(selected_columns, list):
            raise ValueError("selected_columns must be a list or None")

        selected_columns = [f'{col}{suffix}' for col in selected_columns]
        # Ensure we only use valid column names
        selected_columns = [col for col in selected_columns if col in data.columns]

        if not selected_columns:
            return data  # No valid columns to filter on, return the original DataFrame

        # Filter rows where all selected columns contain 'X'
        #filtered_data = data[data[selected_columns].eq('X').all(axis=1)].drop(columns=selected_columns)
        filtered_data = data.drop(columns=selected_columns, errors='ignore')

        return filtered_data

    def filtering_mhc_misclass(self, data, col, select_alleles):
        if len(select_alleles) == 0 or select_alleles is None:
            return data

        if isinstance(select_alleles, list):
            pattern = '|'.join(select_alleles)
            data = data[~data[col].str.contains(pattern, regex=True, na=False)]
        return data
    
    def filtering_anarci_errors(self, data, selected_columns):
        return data[~data[selected_columns].apply(lambda x: x.str.contains('*', regex=False)).any(axis=1)]

    def select_columns_final(self, data):

        if self.description != 'train':
            # Filter columns based on the description
            pattern = r'^(\d+[a-b]+|\d+pep|\d+mhc|class)' # add 'class' to the pattern
            selected_columns = data.filter(regex=pattern)

        elif self.description == 'train' and self.filter_rows:
            pattern = r'^(\d+[a-b]+|\d+pep|\d+mhc)'
            selected_columns = data.filter(regex=pattern)
        
        elif self.description == 'train' and not self.filter_rows:
            pattern = r'^(\d+[a-b]+|\d+pep|\d+mhc|class)' # add 'class' to the pattern
            selected_columns = data.filter(regex=pattern)

        # Save the selected columns to a text file
        columns = selected_columns.columns.tolist()
        with open(f'selected_positions.txt', 'w') as f:
            for col in columns:
                f.write(f"{col}\n")
            
        # Preserve the original order of columns
        ordered_columns = [col for col in data.columns if col in selected_columns.columns]
        ordered_columns = ['TCR_ID'] + ordered_columns  # Ensure 'TCR_ID' is the first column
        return data[ordered_columns]

    def run(self):

        print('---------------------START------------------------------------')
        print('Before start filter...')
        print("Samples:", len(self.input_df))
        print('Shape:', self.input_df.shape)
        print('Columns:', self.input_df.columns.tolist())
        
        if self.filter_rows:
            if self.just10x and self.description == 'train':
                input_df = self.filtering_just10x(self.input_df)
                print(f'\nAfter filter by just10x: {len(self.input_df)}')
                print('Shape:', self.input_df.shape)
            else:
                print('\nNo filter by just10x\n')
                input_df = self.input_df.copy()

            filtered_data = self.filtering_PepLen(input_df, 'peptide', self.max_pep_length)
            
            print(f'After filter by peptide length and columns: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)

            if self.description == "train":
                alleles_to_exclude = self.exclude_contacts.get('alleles', [])
                filtered_data = self.filtering_mhc_misclass(filtered_data, 'MHCa', alleles_to_exclude)

            print(f'After filter by mhc misclassified: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)

            filtered_data = self.filtering_anarci_errors(filtered_data, ['TRA', 'TRB'])

            print(f'After filter by anarci errors: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)

        else:
            print('\nNo filter by rows!\n')
            input_df = self.input_df.copy()
            pep_pos_to_exclude = self.exclude_contacts.get('peptide', [])
            TRA_pos_to_exclude = self.exclude_contacts.get('alpha', [])
            TRB_pos_to_exclude = self.exclude_contacts.get('beta', [])
            mhc_pos_to_exclude = self.exclude_contacts.get('mhc', [])

            if self.selected_factors is not None:
                columns_to_keep = ['TCR_ID'] + self.selected_factors

                #ADD COLUMNS NOT IN SELECTED FACTORS
                diff_columns =  set(columns_to_keep) - set(self.input_df.columns)
                if diff_columns:
                    print(f'Adding columns not in selected factors: {diff_columns}')
                    
                    #fill with X not in selected factors
                    for col in diff_columns:
                        if col.endswith(('a','b','pep','mhc')):
                            input_df[col] = 'X'
                
                input_df = input_df[columns_to_keep]

                print(f'After select columns: {len(input_df)}')
                print('Shape:', input_df.shape)
                
            else:
                print('No selected factors provided, using all columns.')
            

            filtered_data = self.filtering_positions(input_df, TRA_pos_to_exclude, 'a')

            print(f'After filter by large alpha cdr: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)

            filtered_data = self.filtering_positions(filtered_data, TRB_pos_to_exclude, 'b')

            print(f'After filter by large beta cdr: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)

            filtered_data = self.filtering_positions(filtered_data, mhc_pos_to_exclude, 'mhc')

            print(f'After filter by mhc positions: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)


            filtered_data = self.filtering_positions(filtered_data, pep_pos_to_exclude, 'pep')
            print(f'After filter by peptide positions: {len(filtered_data)}')
            print('Shape:', filtered_data.shape)

        self.output_df = self.select_columns_final(filtered_data)

        #filtered_data_final.to_csv(outfile, index=False)
        print(f'---------------------DONE------------------------------------')

    def result(self, outfile:str = None):
        
        if 'TRA' in self.output_df.columns:
            raise ValueError("Some columns are still present in the DataFrame. Please check the filtering steps.")
        
        if outfile is not None:
            self.output_df.to_csv(outfile, index=False)

        return self.output_df