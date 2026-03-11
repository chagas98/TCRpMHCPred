import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from anarci import number
from os import cpu_count
import multiprocessing
import argparse
import timeit
import pickle
import json
import os

def load_pickle_dict(file_path, tmp_path):
    """
    Load a dictionary from a pickle file.
    
    Args:
        file_path (str): Path to the pickle file.
        
    Returns:
        dict: The loaded dictionary.
    """
    with open(file_path, "rb") as file:
        data_dict = pickle.load(file)
    
    name = file_path.split("/")[-1].replace(".pkl", "")
    
    try:

        with open(f"{tmp_path}check_result_{name}.txt", "w") as f:
            for key, value in data_dict.items():
                print(f"{key}: {''.join(value['aa'])}", file=f)

    except Exception as e:
        print(f"Error loading pickle file: {e}")

class ANARCIndexer:
    def __init__(self, input_df:pd.DataFrame, 
                 selections_include_path:str, 
                 selections_exclude_path:str, 
                 tmp_path:str, 
                 description:str,
                 scheme:str = 'aho',
                 keep_variable_domain:bool = False):
        """
        Initialize the ANARCIndexer with input DataFrame, selections path, and description.
        
        :param input_df: DataFrame containing TCR sequences.
        :param selections_path: Path to the JSON file containing CDR ranges and residues selections.
        :param description: Description of the dataset class.
        """
        self.input_df = input_df
        self.tmp_path = tmp_path
        self.description = description
        self.keep_variable_domain = keep_variable_domain
        self.scheme = scheme

        # Load the JSON file
        if selections_include_path is not None and selections_exclude_path is not None:
            with open(selections_include_path, "r") as f:
                to_include = json.load(f)

            # Convert lists back to tuples
            self.alpha_to_include = {key: tuple(value) for key, value in to_include["alpha"].items()}
            self.beta_to_include = {key: tuple(value) for key, value in to_include["beta"].items()}

            # Load the JSON file
            with open(selections_exclude_path, "r") as f:
                to_exclude = json.load(f)

            # Convert lists back to tuples
            self.alpha_to_exclude = {key: tuple(value) for key, value in to_exclude["alpha"].items()}
            self.beta_to_exclude = {key: tuple(value) for key, value in to_exclude["beta"].items()}

            # Print to verify
            print("alpha_to_include =", self.alpha_to_include)
            print("beta_to_include =", self.beta_to_include)

    def get_unique_seqs(self, df, cols):
        print(f'Original Length: {len(df)}')
        df_unique = df.drop_duplicates(subset=cols)
        print(f'Deduplicated: {len(df_unique)}')
        return df_unique

    def get_anarci_number(self, tcr_seq):
        try:

            results, _ = number(tcr_seq, scheme=self.scheme, allowed_species='human')
            seq_numbered = [(pos[0],aa) for pos, aa in results if aa != '-']
            #indices, residues = zip(*seq_numbered)
            
            dtype = np.dtype([('ind', np.int32), ('aa', 'U10')])
            vector = np.array(seq_numbered, dtype=dtype)
            
            return vector
        
        except Exception as e:
            print(f"Error annotating sequence: {e}")
            return None


    def annotate_cdr_positions(self, tcr_pair):

        tcrA, tcrB, tcr_id, index = tcr_pair

        if index%1000 == 0:
            print(f'Get numbers from {index}')
        try:
            gapA = tcrA.find('*')
            gapB = tcrB.find('*')
            if gapA == -1 and gapB == -1:
                tcrA_numbered = self.get_anarci_number(tcrA)
                tcrB_numbered = self.get_anarci_number(tcrB)

                return tcr_id, tcrA_numbered , tcrB_numbered   # if no results found
            
            else:
                print('Find *')
                return tcr_id, None, None

        except Exception as e:
            print(f"Error annotating sequence: {e}")
            return None, None
        
    def numbers_counter(self, tcr_numbered):
        all_num = np.concatenate(tcr_numbered)['ind']
        unique, counts = np.unique(all_num, return_counts=True)
        return unique, counts

    def gather_all_counts(self, numbered_seqs):
        """ Gather counts and CDR positions from numbered sequences.
        :param numbered_seqs: Dictionary containing numbered TCR sequences.
        :return: List of tuples containing DataFrame with counts, CDR positions, and name.
        """

        list_of_counts =[]
        print(type(numbered_seqs))
        for name, values in numbered_seqs.items():

            try:                
                numbered_tcr, to_include, to_exclude = values
                all_num = np.concatenate(np.array(numbered_tcr))['ind']
                unique, counts = np.unique(all_num, return_counts=True)
                df = pd.DataFrame({'res': unique, 'counts': counts})

                print(f'Including residues: {to_include}')
                if to_include is not None:
                    cdr_positions={}
                    for cdr_name, pos in to_include.items():
                        
                        if cdr_name not in ['TRA', 'TRB']:
                            cdr = df.loc[df['res'].isin(pos)]['res'].tolist()
                            
                            #add outer from pos 
                            #missing = list(set(pos) - set(cdr))
                            #cdr += missing  ---------------------------> # TODO this is not needed, outer residues are not included in the CDRs
                            cdr_positions[cdr_name] = sorted(cdr)
                        
                        else:
                            cdr = df['res'].tolist()
                            
                            #add outer from pos 
                            missing = list(set(pos) - set(cdr))
                            cdr = sorted(cdr + missing)                       
                            cdr_positions[cdr_name] = cdr
                
                # exclude residues are superior in comparison to the included residues
                print("Excluding residues:", to_exclude)
                if to_exclude is not None:
                    if all(key in to_include for key in to_exclude):
                        for cdr_name, pos in to_exclude.items():
                            cdr_prev = cdr_positions[cdr_name]                  
                            cdr_updated = [res for res in cdr_prev if res not in pos]
                            cdr_positions[cdr_name] =  cdr_updated
      
                    else:
                        raise ValueError(f"{name} does not have the same keys in to_include and to_exclude.")
                        
                list_of_counts.append((df, cdr_positions, name))

            except Exception as e:
                print(f"Error gather all counts: {e}")

        return list_of_counts

    def multiprocessing_annotation(self, data_unique):
        """
        Annotate TCR sequences in parallel using multiprocessing.
        
        :param data_unique: DataFrame containing unique TCR sequences.
        :return: Tuple of Series containing annotated TCR sequences.
        """
        
        num_workers = max(1, int(0.9 * cpu_count()))
        if num_workers > 25:
            num_workers = 25
        print(f"Using {num_workers} CPU cores for multiprocessing.")
        index_df = range(len(data_unique))
        with multiprocessing.Pool() as pool:
            results = pool.map(self.annotate_cdr_positions, 
                               zip(data_unique['TRA'], data_unique['TRB'], data_unique['TCR_ID'], index_df))
        
        id, tcr_num_A, tcr_num_B = zip(*results)
        return id, tcr_num_A, tcr_num_B
    
    def run(self):

        tic=timeit.default_timer()

        # unique sequences is almost the same as the input_df, lets use self.input_df directly
        #data_unique = self.get_unique_seqs(self.input_df, ['TRA', 'TRB']).reset_index()

        # Extract TCR numbered sequences
        self.input_df.reset_index(inplace=True, drop=True)
        tcr_id, tcr_numbered_A, tcr_numbered_B =  self.multiprocessing_annotation(self.input_df)
        
        if self.keep_variable_domain:
            print('Keep variable domains')
            self._variable_domains = self.input_df.copy()
            
            print('Lengths:', len(tcr_numbered_A), len(tcr_numbered_B), len(tcr_id))
            for id, tcr_A, tcr_B in zip(tcr_id, tcr_numbered_A, tcr_numbered_B):
                self._variable_domains.loc[self.input_df['TCR_ID'] == id, 'TRA'] = ''.join(tcr_A['aa']) if tcr_A is not None else None
                self._variable_domains.loc[self.input_df['TCR_ID'] == id, 'TRB'] = ''.join(tcr_B['aa']) if tcr_B is not None else None
            print('Variable domains saved')
            return self._variable_domains
               
        _tcr_numbered_A_ = pd.Series(tcr_numbered_A).dropna()
        _tcr_numbered_B_ = pd.Series(tcr_numbered_B).dropna()

        numbered_seqs = {f'alpha_{self.description}': (_tcr_numbered_A_, self.alpha_to_include, self.alpha_to_exclude),
                         f'beta_{self.description}': (_tcr_numbered_B_, self.beta_to_include, self.beta_to_exclude)}

        # Gather counts and CDR positions
        print('Gathering counts and CDR positions...')
        results_counts = self.gather_all_counts(numbered_seqs)

        tcr_num_idA = {
            tcr_id[i]: tcr_numbered_A[i] for i in range(len(tcr_numbered_A)) if tcr_numbered_A[i] is not None
        }
        tcr_num_idB = {
            tcr_id[i]: tcr_numbered_B[i] for i in range(len(tcr_numbered_B)) if tcr_numbered_B[i] is not None
        }

        self._positions_selected = {}
        for final_data in results_counts:

            final_data[0].to_csv(f'{self.tmp_path}ANARCI_freq_{final_data[2]}.csv')

            with open(f'{self.tmp_path}ANARCI_cdr_{final_data[2]}.json', 'w', encoding='utf-8') as f:
                json.dump(final_data[1], f, ensure_ascii=False, indent=4)
                f.close()
            
            if 'alpha' in final_data[2]:
                self._positions_selected['TRA'] = final_data[1]
            elif 'beta' in final_data[2]:
                self._positions_selected['TRB'] = final_data[1]
            else:
                raise ValueError(f"Unknown TCR type in {final_data[2]}")
            
        with open(f'{self.tmp_path}ANARCI_alpha_{self.description}.pkl', 'wb') as f:
            pickle.dump(tcr_num_idA, f)
            f.close()
        
        load_pickle_dict(f'{self.tmp_path}ANARCI_alpha_{self.description}.pkl', tmp_path=self.tmp_path) # for checking results

        with open(f'{self.tmp_path}ANARCI_beta_{self.description}.pkl', 'wb') as f:
            pickle.dump(tcr_num_idB, f)
            f.close()

        load_pickle_dict(f'{self.tmp_path}ANARCI_beta_{self.description}.pkl', tmp_path=self.tmp_path) # for checking results
        
        print('-----------------------------------------------------------')
        toc=timeit.default_timer()
        print(f'Time: {toc-tic}')

        self._anarci_positions = {'TRA': tcr_num_idA, 'TRB': tcr_num_idB}
        
    @property
    def anarci_positions(self):
        return self._anarci_positions

    @property
    def positions_selected(self):
        return self._positions_selected
    
    @property
    def variable_domains(self):
        return self._variable_domains