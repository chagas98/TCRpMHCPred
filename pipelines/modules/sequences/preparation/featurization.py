import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import os
import argparse
import pickle
import json
import timeit

def join_aho_numbers(aho_numbers):
    with open(aho_numbers[0]) as f:
        data1 = json.load(f)
    
    with open(aho_numbers[1]) as f: 
        data2 = json.load(f)

    data1.update(data2)

    print('AHO numbers updated.')
    return data1

class GetPositionsTCR():
    def __init__(self, dataset: pd.DataFrame, contact_to_include: dict[list] = None, ref_anarci: dict[str] = None):
        
        self.dataset = dataset
        self.ref_anarci = ref_anarci
        
        if isinstance(contact_to_include, dict):

            self.targets = contact_to_include.keys()
            self.contact_num = contact_to_include

            if isinstance(contact_to_include.keys(), str):
                self.targets = [contact_to_include.keys()]
            elif isinstance(contact_to_include.keys(), list):
                self.targets = contact_to_include.keys()
        elif all(contact_to_include.values() == "None"):
            self.targets = contact_to_include.keys()
            self.contact_num = None

        self._filterNA_anarci_seq_()
        self._create_empty_indexed_vector_()
                 
    def _get_len_max_(self, cols: list[str] | str):
        df = self.dataset
        lengths = {}

        for col  in cols:
            lengths[col] = df[col].apply(lambda x: len(x)).max()
            
            indices_len_max = df[df[col].str.len() == lengths[col]].index
            if len(indices_len_max) > 1:
                print(f'Multiple indices with max length for {col}')
                TCR_cols = self.ref_anarci.keys()
                data_max = self.dataset.iloc[indices_len_max].drop_duplicates(subset=TCR_cols)
                
                if len(data_max) == 1:
                    print(f'Same TCR with max length for {col}')
                    print('************************************************')
                
                else:
                    raise Exception(f'Different TCRs with max length for {col}')

        return lengths

    def _check_ref_indices_(self, tcr_chain): 
        # essa função vai organizar os indices de contato do anarci, para que sejam ordenados e possamos criar um vetor vazio indexado
        for region in self.contact_num[tcr_chain].keys():    
            self.contact_num[tcr_chain][region] = np.sort(self.contact_num[tcr_chain][region]).tolist()
    
    def _filterNA_anarci_seq_(self):
        
        for tcr, positions_num in self.ref_anarci.items():
            
            old_len = len(self.dataset)
            
            self.dataset = self.dataset[self.dataset['TCR_ID'].isin(positions_num.keys())]

            new_len = len(self.dataset[~pd.isna(self.dataset['TCR_ID'])])

            print(f'\n\nFrom {old_len} of samples. {old_len - new_len} not have anarci assignment in {tcr}.')
            print(f'New N samples: {new_len}\n\n')  

        return self.dataset
                
    def _create_empty_indexed_vector_(self): # cria um vetor vazio indexado com os indices do aho pra se incluir
        self.indexed_vectors = {}
        print(self.contact_num)
        for tcr_chain in self.contact_num.keys():
            
            # sorting the indices to include
            self._check_ref_indices_(tcr_chain)

            self.indexed_vectors[tcr_chain] = {}
            for region, indices_to_include in self.contact_num[tcr_chain].items():
                
                length = len(indices_to_include)
                dtype = np.dtype([('ind', np.int32), ('aa', 'U10')])
                empty_vector = np.zeros(length, dtype=dtype)
                empty_vector['ind'] = indices_to_include
                self.indexed_vectors[tcr_chain][region] = empty_vector
                
    def _create_indexed_seq_vector_(indices: list[int], seq: list[str] = None): # cria um vetor da sequencia do dataset indexado

        if seq is None or len(indices) != len(seq):
            print('len(indices) and len(seq) are not equal')
            return None
        
        dtype = np.dtype([('ind', np.int32), ('aa', 'U10')])
        vector = np.zeros(len(seq), dtype=dtype)
        
        vector['ind'], vector['aa'] = np.array(indices), np.array(seq)
        
        return vector

    def _assign_index_seq_(self, empty_vector, tcr_chain, tcr_id):
        seq = self.ref_anarci[tcr_chain][tcr_id]
        seq = np.array(seq, dtype=object)
        
        # Cria um dicionário para acesso rápido
        seq_dict = {int(ind): aa for ind, aa in seq}

        # Preenche os aminoácidos com base nos índices
        aa_filled = [seq_dict.get(ind, 'X') for ind in empty_vector['ind']]
        empty_vector['aa'] = aa_filled
        
        return aa_filled 


    def get_indexed_seq_vector(self, tcr_chain, keep_regions: bool = False, keep_gaps: bool = False):
    
        for region in self.contact_num[tcr_chain].keys():
            
            print('Assigning index to', region, 'of', tcr_chain)
            
            emptyvec = self.indexed_vectors[tcr_chain][region]
            
            chain = tcr_chain[-1].lower()
            df = self.dataset.apply( 
                lambda x: self._assign_index_seq_(emptyvec, tcr_chain, x['TCR_ID']), axis=1
                )
            # join the amino acids into a single string
            
            if keep_regions:
                new_cols = [region]
                if keep_gaps:
                    df = df.apply(lambda x: ''.join(x))
                else:
                    df = df.apply(lambda x: ''.join(x).replace('X', ''))

            else:
                new_cols = [f"{num}{chain}" for num in emptyvec['ind']]

            df = df.apply(pd.Series)
            df.columns = new_cols

            self.dataset = pd.concat([self.dataset, df], axis=1)
        return self.dataset
    
    def save_dataset(self):
        return self.dataset

class GetPositionsPep():
    def __init__(self, dataset: pd.DataFrame, peptide_col: str, max_len: int = None):
        self.dataset = dataset
        self.peptide_col = peptide_col
        self.max_len = max_len

    def extract_positions(self, peptide, len_max):
        sequence = [char for char in peptide]
        #del sequence[self.to_exclude[0]], sequence[self.to_exclude[1]]
        
        required_length = len_max - len(sequence)
        if required_length > 0:
            sequence.extend(['X'] * required_length)
        return sequence
    
    def generate_embed(self):
        
        df_pep_len_max = self.dataset[self.peptide_col].str.len().max()
        len_max = self.max_len if self.max_len is not None else df_pep_len_max
        df_pep_len_max = max(df_pep_len_max, len_max)
            
        df = self.dataset[self.peptide_col].apply(lambda x: self.extract_positions(x, len_max))
        df = df.apply(pd.Series)
        df.columns = [f'{i}pep' for i in range(1, df_pep_len_max+1)]
        
        self.dataset = pd.concat([self.dataset, df], axis=1)

        return self.dataset
    
    def save_dataset(self):
        return self.dataset

class GetPositionsMHC():
    def __init__(self, dataset: pd.DataFrame, to_include: list, mhc_col: str, 
                 second_mhc_col: str = None, alignment_ref: str = None):
        self.dataset = dataset
        self.to_include = to_include
        self.alignment_ref = alignment_ref
        self.mhc_col = mhc_col
        self.sec_mhc_col = second_mhc_col
        
        if alignment_ref != "None":
            self.alignment = AlignIO.read(str(alignment_ref), "fasta")

    def _extract_mhc_seq_(self, allele):
        selections = []

        for record in self.alignment:
            if record.id == allele[self.mhc_col]:
                selections.append(record)

            elif self.sec_mhc_col != None and record.id == allele[self.sec_mhc_col]:
                selections.append(record)
    
        try:
            if len(selections) == 1:
                return str(selections[0].seq)
            else:
                return None
        except:
            return None
    
    def extract_positions(self, mhc_seq):
        
        if len(self.to_include) != 0:
            if mhc_seq is not None:
                sequence = [char for i, char in enumerate(mhc_seq) if i + 1 in self.to_include]
                return sequence        
            else:
                return [None] * len(self.to_include)

        elif len(self.to_include) == 0:
            return []
        
        elif self.to_include == "None":
            return [char for char in mhc_seq]
        

    def generate_embed(self):

        self.dataset['MHC_trimmed'] = self.dataset.apply(lambda x: self._extract_mhc_seq_(x), axis=1)
        df = self.dataset['MHC_trimmed'].apply(lambda x: self.extract_positions(x))
        print(self.dataset['MHC_trimmed'])
        print(self.dataset['MHC_trimmed'].isna().sum())
        len_max = self.dataset['MHC_trimmed'].apply(lambda x: len(x) if isinstance(x, str) else 0).max()
        df = df.apply(pd.Series)

        if len(self.to_include) != 0:
            df.columns = [f'{i}mhc' for i in self.to_include]
        
        elif len(self.to_include) == 0:
            df.columns = [f'{i}mhc' for i in range(1, len_max+1)]
        
        elif self.to_include == "None":
            df.columns = []
        
        self.dataset = pd.concat([self.dataset, df], axis=1)

        return self.dataset
    
    def save_dataset(self):
        return self.dataset


def pos_mhc_interface(path, threshold):
    data = pd.read_csv(path, header=None, sep=' ')
    data = data[data.iloc[:, 2] > threshold]
    print(data.iloc[:, 0].tolist())
    return data.iloc[:, 0].tolist()

def just10x_assignment(dataset: pd.DataFrame, ref_dataset: str):
    # Load reference values
    ref_values = pd.read_csv(ref_dataset)[['Reference_IEDB', 'Ref', 'PMID', 'VDJdb_alias', 'Exclude']]
    
    # Ensure is treated as boolean and filter
    ref_values['Exclude'] = ref_values['Exclude'].astype(bool)
    excluded_refs = ref_values[ref_values['Exclude']]

    #classify as true/false reference from dataset
    exclude_refs = excluded_refs[['Reference_IEDB', 'Ref', 'PMID', 'VDJdb_alias']].values.flatten()
    dataset['just10x'] = dataset['Reference'].isin(exclude_refs)
    return dataset


class GetPositionsSeq:
    def __init__(self, dataset: pd.DataFrame, 
                 anarci_positions_selected: list = None,
                 anarci_positions: str = None,
                 mhc_interface_contacts_path: str = None,
                 mhc_alignment_path: str = None,
                 mhc_contacts_threshold: int = None,
                 just10x_table_path: str = None,
                 max_len: int = None,
                 description: str = "train",
                 tmp_path: str = 'tmp/', 
                 outfile: str = 'combined_datasets_encoded.csv'):
        
        self.dataset = dataset
        self.mhc_interface_contacts_path = mhc_interface_contacts_path
        self.mhc_alignment_path = mhc_alignment_path
        self.mhc_contacts_threshold = mhc_contacts_threshold
        self.just10x_table = just10x_table_path
        self.max_len = max_len
        self.tmp_path = tmp_path
        self.description = description
        self.outfile = outfile

        # ANARCI data
        # positions that occur in the dataset
        if anarci_positions_selected is None:
            raise ValueError("anarci_contacts_path not provided")
        else:
            self.anarci_positions_selected = anarci_positions_selected
        
        # positions that are indexed in the pickle files
        if anarci_positions is None:
            raise ValueError("anarci_positions_path must be provided")
        else:
            self.anarci_positions = anarci_positions
        
        
        print('PARAMETERS DATA')
        print('------------------------------------------------------------------------------')
        print('TRA anarci data:', len(anarci_positions['TRA'].keys()))
        print('TRB anarci data:', len(anarci_positions['TRB'].keys()))
        print('Contacts to include:',anarci_positions_selected)
        print('MHC interface unique contacts:', mhc_interface_contacts_path)
        print('MHC alignment fasta file:', mhc_alignment_path)
        print('MHC interface contacts threshold:', mhc_contacts_threshold)
        print('Just10x classification:', just10x_table_path)
        print('Max length of peptides:', max_len)
        print('Description:', description)
        print('------------------------------------------------------------------------------')
    
    def run(self):

        tic = timeit.default_timer()

        if self.description == "train":
            print('Just10x classification will be applied to the dataset.')
            dataset = just10x_assignment(self.dataset, self.just10x_table)
        else:
            print('Just10x classification will not be applied to the dataset.')
            dataset = self.dataset.copy()

        print('------------------------------------------------------------------------------')
        print('INITIAL DATA LENTGH:', len(dataset))
        print('Shape:', dataset.shape)
        print('------------------------------------------------------------------------------')
        
        run_positionsTCR = GetPositionsTCR(dataset            = dataset,
                                           contact_to_include = self.anarci_positions_selected,
                                           ref_anarci         = self.anarci_positions)
        run_positionsTCR.get_indexed_seq_vector('TRA')
        run_positionsTCR.get_indexed_seq_vector('TRB')
        self.features_df = run_positionsTCR.save_dataset()

        print('After GetPositionsTCR:', len(self.features_df))
        print('Shape:', self.features_df.shape)

        run_positionsPep = GetPositionsPep(dataset     = self.features_df,
                                           max_len     = self.max_len,
                                           peptide_col = 'peptide')
        run_positionsPep.generate_embed()
        self.features_df = run_positionsPep.save_dataset()
        
        print('After GetPositionsPep:', len(self.features_df))
        print('Shape:', self.features_df.shape)

        mhc_index_to_include = pos_mhc_interface(self.mhc_interface_contacts_path,
                                                 self.mhc_contacts_threshold)

        run_positionsMHC = GetPositionsMHC(dataset        = self.features_df,
                                           to_include     = mhc_index_to_include,
                                           mhc_col        = 'assigned_allele',
                                           second_mhc_col = 'MHCa',
                                           alignment_ref  = self.mhc_alignment_path)
        
        run_positionsMHC.generate_embed()
        self.features_df = run_positionsMHC.save_dataset()

        print('After GetPositionsMHC:', len(self.features_df))
        print('Shape:', self.features_df.shape)

        # End time (toc)
        toc = timeit.default_timer()
        print(f"Execution time: {toc - tic:.6f} seconds")
    
    @property
    def result(self):
        return self.features_df