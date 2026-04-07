import sys
import os
from urllib import response
import pandas as pd
import json
import argparse
import datetime
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from anarci import number
from modules.sequences.preparation.Seqs import RunThimble, MHCSeqAssigner
from modules.sequences.preparation.indexer import ANARCIndexer
from modules.sequences.preparation.featurization import GetPositionsSeq
from modules.sequences.preparation.filtering import Filtering
from modules.sequences.datasets import IMMREP25dataset, ConcatDataset
from datetime import datetime

_AHO_CDR_RANGES = {
    "CDR1": list(range(25, 41)),
    "CDR2": list(range(58, 78)),
    "CDR3": list(range(106, 140)),
    "CDR2.5": list(range(83, 89))
}


def get_datasets(location: str, map_datasets: dict, description: str):
    """ Load and concatenate datasets for TCR analysis.
    
    This function loads multiple TCR datasets, processes them, and returns a concatenated DataFrame.
    It also handles the mapping of columns for the test dataset.
    
    Returns:
        pd.DataFrame: Concatenated DataFrame of all datasets.
    """

    # Map Test Dataset Columns
    map_cols = {
        'tcra_cdr1': 'CDR1A',
        'tcra_cdr2': 'CDR2A',
        'tcra_cdr3': 'CDR3A',
        'tcra_v': 'VA',
        'tcra_j': 'JA',
        'tcra_trimmed': 'TRA',
        'tcrb_cdr1': 'CDR1B',
        'tcrb_cdr2': 'CDR2B',
        'tcrb_cdr3': 'CDR3B',
        'tcrb_v': 'VB',
        'tcrb_j': 'JB',
        'tcrb_trimmed': 'TRB',
        'peptide': 'epitope',
        'hla': 'assigned_allele',
        'hla_sequence': 'MHCa_sequence',
        'hla_trimmed': 'MHCseq',
        'label': 'class'
    }
    
    location_db = f'{location}/immrep25'
    immrep25_file = f'{location_db}/00_immrep25/{map_datasets["immrep25"]}'
    dataImmRep25 = pd.read_csv(immrep25_file, sep='\t')
    dataImmRep25.rename(columns=map_cols, inplace=True)
    print(dataImmRep25.head())
    dataImmRep25 = dataImmRep25[map_cols.values()]
    dataImmRep25.sort_values(by='epitope', inplace=True)
    dataImmRep25.reset_index(drop=True, inplace=True)

    dataImmRep25['Reference'] = 'immrep25'
    dataImmRep25['Score'] = 3
        
    print(f"Concatenated dataset contains {dataImmRep25.shape[0]} rows and {dataImmRep25.shape[1]} columns.")
    print(50*"-")

    return dataImmRep25, location_db

def _select_cdr_seq(variable_domain_seq, numbers, cdr_name):
    cdr_ranges = _AHO_CDR_RANGES.get(cdr_name)
    if not cdr_ranges:
        raise ValueError(f"Unknown CDR name: {cdr_name}")

    start, end = cdr_ranges[0], cdr_ranges[-1]
    selected_seq = ''.join(
        aa for aa, num in zip(variable_domain_seq, numbers)
        if start <= num <= end
    )
    return selected_seq

def _run_anarci(seq):
    try:

        results, chain = number(seq, scheme='aho', allowed_species='human')
        
        if not results:
            return None, None, None 
        
        seq_numbered = [(pos[0], aa) for pos, aa in results if aa != '-']
        amino_acids = "".join([seq[1] for seq in seq_numbered])

        numbers = [seq[0] for seq in seq_numbered]
        return amino_acids, numbers, chain

    except Exception as e:
        print(f"Error annotating sequence: {e}")
        return None, None, None

def run_indexer(df):

    rows = []
    for i, row in df.iterrows():

        if (i+1) % 100 == 0:
                print(f"Processed {i+1} sequences...")

        if pd.isna(row['TRA']) and pd.isna(row['TRB']):
            print(f"Skipping row {i} due to missing TRA and TRB sequences.")
            continue

        for chain in ['A', 'B']:

            variable_domain_seq, nums, chain = _run_anarci(seq=row[f'TR{chain}'])
            
            if variable_domain_seq and chain in ['A', 'B', 'D']:

                chain = 'A' if chain in ['A', 'D'] else chain  # Treat 'D' as 'A' for TRA

                row[f'TR{chain.upper()}'] = variable_domain_seq
                row[f'TR{chain.upper()}_num'] = nums
                row[f'CDR1{chain.upper()}'] = _select_cdr_seq(variable_domain_seq, nums, "CDR1")
                row[f'CDR2{chain.upper()}'] = _select_cdr_seq(variable_domain_seq, nums, "CDR2")
                row[f'CDR3{chain.upper()}'] = _select_cdr_seq(variable_domain_seq, nums, "CDR3")
        
            else:
                print(f"Skipping row {i} due to missing TR{chain} sequence.")
                continue

        # gather rows
        rows.append(row)
    
    df_indexed = pd.DataFrame(rows)

    print("Completed ANARCI indexing for all sequences.")
    
    return df_indexed

def extract_mhc_variableregion(allele, aligned_mhc_records):
        selections = []

        for record in aligned_mhc_records:
            if record.id == allele['assigned_allele']:
                selections.append(record)

            elif record.id == allele['MHCa']:
                selections.append(record)
    
        try:
            if len(selections) == 1:
                return str(selections[0].seq)
            else:
                return None
        except:
            return None
   
def split_basedon_score(df):

    unique_scores = df['Score'].dropna().unique()
    
    all_ids = []
    for score in unique_scores:
        print(f"Score found in dataset: {score}")
        print(df[df['Score'] == score].head())
        
        df_score = df[df['Score'] == score]
        
        cdate = datetime.now()
        filename_path = f'databases/fullSeqs_dataset_score{score}_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv'
        df_score.to_csv(filename_path, index=False)

        ids = df_score['TCR_ID'].tolist()
        all_ids.extend(ids)

        #check for duplicates
        if len(ids) != len(set(ids)):
            print(f"Warning: Duplicate TCR_IDs found for score {score}!")

    # check for repeat ids across scores
    if len(all_ids) != len(set(all_ids)):
        print(f"Warning: Duplicate TCR_IDs found across scores!")

    print(f"Total unique TCR_IDs across all scores: {len(set(all_ids))}")

    return 

def preparation(df, location_db=None):
    """ Main function to run the preparation pipeline.
    This function orchestrates the entire process of preparing the dataset
    by running Thimble, MHC assignment, ANARCI indexing
    """
    
    # Temporary directory for intermediate files
    tmp =  'tmp/'
    os.makedirs(tmp, exist_ok=True)
    
    print("Starting ANARCI indexing...")
    fullSeqs_indexed_df = run_indexer(df)
    print("ANARCI indexing completed.")

    #sort columns
    fullSeqs_indexed_df['TCR_ID'] = fullSeqs_indexed_df.index.to_series().apply(lambda x: f'immrep{x}')
    column_order = ['TCR_ID', 
                    'TRA','CDR1A', 'CDR2A', 'CDR3A', 'TRA_num',
                    'TRB','CDR1B', 'CDR2B', 'CDR3B', 'TRB_num',
                    'epitope', 'assigned_allele', 'MHCa_sequence', 'MHCseq',
                    'Reference', 'Score', 'class']
    fullSeqs_indexed_df = fullSeqs_indexed_df[column_order]
    print("Final DataFrame has shape:", fullSeqs_indexed_df.shape)
    print(fullSeqs_indexed_df.columns)

    cdate = datetime.now()
    fullSeqs_indexed_df.to_csv(f"{location_db}/02_processed/fullSeqs_dataset_immrep25_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv", index=False)

    # Split based on score
    #split_basedon_score(fullSeqs_indexed_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run get_datasets script to load and concatenate TCR datasets.')
    parser.add_argument("-i", "--datasets_list", nargs='+', help="choose datasets to include", type=str)
    parser.add_argument("-f", "--datasets_dir", help="location of public datasets", type=str, default='../databases')
    parser.add_argument("-m", "--map_datasets", help="map dataset file names", type=json.loads, 
                        default='{"immrep25": "immrep2025_for_release.tsv"}')
    parser.add_argument("-d", "--description", help="description of the dataset", type=str, default='train')
    parser.add_argument("-l", "--max_pep_length", help="Maximum peptide length for filtering", type=int, default=30)
    parser.add_argument("-r", "--ref_mhc", help="Reference MHC sequences in FASTA format", type=str, default='./assets/hla_prot_includeMusMusculus.fasta')
    parser.add_argument("-a", "--aligned_mhc", help="Aligned MHC sequences in FASTA format", type=str, default='./assets/A02010101_aligned_trimmed24_180_addGinit.fasta')
    parser.add_argument("--thimble_path", help="Path to Thimble executable", type=str, default='thimble')
    parser.add_argument("--outdir", help="Output directory", type=str, default='output')
    args = parser.parse_args()

    print("-" * 50 +
        f"\nRunning get_datasets with parameters:\n"
        f"Datasets List: {args.datasets_list}\n"
        f"Location: {args.datasets_dir}\n"
        f"Map Datasets: {args.map_datasets}\n"
        f"Description: {args.description}\n"
        f"Max Peptide Length: {args.max_pep_length}\n"
        f"Thimble Path: {args.thimble_path}\n"
        f"Aligned MHC: {args.aligned_mhc}\n"
        f"Reference MHC: {args.ref_mhc}\n"
        f"Output Directory: {args.outdir}\n" +
        "-" * 50)

    df, location_db = get_datasets(location  = args.datasets_dir,
                                   map_datasets  = args.map_datasets,
                                   description   = args.description)

    preparation(location_db=location_db,
                df=df)
    
        

#python3 process_sequences.py -i vdjdb --outdir /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb