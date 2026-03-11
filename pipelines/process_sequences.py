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
from modules.sequences.datasets import VDJdbdataset, McPASTCRdataset, IEDBdataset, NeoTCRdataset, PosTRAITdataset, TESTdataset, ConcatDataset
from datetime import datetime

_AHO_CDR_RANGES = {
    "CDR1": list(range(25, 41)),
    "CDR2": list(range(58, 78)),
    "CDR3": list(range(109, 138)),
    "CDR2.5": list(range(83, 89))
}


def get_datasets(datasets_list: list, location: str, map_datasets: dict, description: str):
    """ Load and concatenate datasets for TCR analysis.
    
    This function loads multiple TCR datasets, processes them, and returns a concatenated DataFrame.
    It also handles the mapping of columns for the test dataset.
    
    Returns:
        pd.DataFrame: Concatenated DataFrame of all datasets.
    """

    # Map Test Dataset Columns
    map_cols = {
        'CDR3a': 'CDR3A',
        'CDR3b': 'CDR3B',
        'peptide': 'epitope',
        'Va': 'VA',
        'Vb': 'VB',
        'Ja': 'JA',
        'Jb': 'JB',
        'MHCa': 'MHCa',
        'Score': 'Score',
        'class': 'class'
    }
    
    datasets = {}
    if 'vdjdb' in datasets_list:
        location_db = f'{location}/VDJdb'
        vdjdb_file = f'{location_db}/00_vdjdb/{map_datasets["vdjdb"]}'
        dataVDJ = VDJdbdataset(path=vdjdb_file, score_threshold=0, save_dir=location_db)
        datasets['vdjdb'] = dataVDJ
    if 'mcpas' in datasets_list:
        location_db = f'{location}/McPAS'
        mcpas_file = f'{location_db}/{map_datasets["mcpas"]}'
        dataMcPAS = McPASTCRdataset(path=mcpas_file)
        datasets['mcpas'] = dataMcPAS
    if 'iedb' in datasets_list:
        location_db = f'{location}/IEDB'
        iedb_file = f'{location_db}/{map_datasets["iedb"]}'
        dataIEDB = IEDBdataset(path=iedb_file)
        datasets['iedb'] = dataIEDB
    if 'neo' in datasets_list:
        location_db = f'{location}/Neo'
        neo_file = f'{location_db}/{map_datasets["neo"]}'
        dataNeo = NeoTCRdataset(path=neo_file)
        datasets['neo'] = dataNeo
    if 'trait' in datasets_list:
        location_db = f'{location}/PosTrait'
        trait_file = f'{location_db}/{map_datasets["trait"]}'
        dataTrait = PosTRAITdataset(path=trait_file)
        datasets['trait'] = dataTrait
    if 'test' in datasets_list:
        location_db = f'{location}/Test'
        test_file = f'{location_db}/{map_datasets["test"]}'
        dataTest = TESTdataset(path=test_file, map_cols=map_cols, outfile=f'test_{description}.csv')
        datasets['test'] = dataTest

    concat_data = ConcatDataset(datasets = list(datasets.values()),
                                labels   = list(datasets.keys()))
    df = concat_data.to_df()
    df.drop_duplicates(subset=['CDR3A','CDR3B','epitope','VA','VB','JA','JB'], inplace=True)
    
    print(f"Concatenated dataset contains {df.shape[0]} rows and {df.shape[1]} columns.")
    print(50*"-")

    cdate = datetime.now()
    filename_path = f'{location_db}/02_processed/combined_dataset_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv'
    df.to_csv(filename_path, index=False)

    return df, filename_path, location_db

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

def preparation(dataset_file, thimble_path, ref_mhc, aligned_mhc, description='train', retries_thimble=False, location_db=None):
    """ Main function to run the preparation pipeline.
    This function orchestrates the entire process of preparing the dataset
    by running Thimble, MHC assignment, ANARCI indexing
    """
    # Temporary directory for intermediate files
    tmp =  'tmp/'
    os.makedirs(tmp, exist_ok=True)

    thimble = RunThimble(input_path  = dataset_file,
                        thimble_path = thimble_path,
                        tmp_path     = tmp,
                        retries      = retries_thimble)
    thimble.run()
    
    # MHC assigner
    print("Starting MHC assignment...")
    mhcassigner = MHCSeqAssigner(fasta_path  = ref_mhc,
                                 tmp_path    = tmp)
    mhcassigner.assign_sequences(input_df    = thimble.result, 
                                 description = description)
    fullSeqs_df = mhcassigner.result
    print("MHC assignment completed.")

    print("Extracting MHC variable regions...")
    mhc_alignment = AlignIO.read(str(aligned_mhc), "fasta")
    fullSeqs_df['MHCseq'] = fullSeqs_df.apply(lambda x: extract_mhc_variableregion(x, mhc_alignment), axis=1)
    print("MHC variable region extraction completed.")


    print("Full sequences DataFrame created with shape:", fullSeqs_df.shape)
    print(fullSeqs_df.head())
    
    print("Starting ANARCI indexing...")
    fullSeqs_indexed_df = run_indexer(fullSeqs_df)
    print("ANARCI indexing completed.")

    print("DataFrame after ANARCI indexing has shape:", fullSeqs_indexed_df.shape)
    fullSeqs_indexed_df.drop_duplicates(subset=['TRA', 'TRB', 'peptide', 'MHCseq'], inplace=True)

    print("Final DataFrame after removing duplicates has shape:", fullSeqs_indexed_df.shape)
    print(fullSeqs_indexed_df.columns)

    cdate = datetime.now()
    fullSeqs_indexed_df.to_csv(f"{location_db}/02_processed/fullSeqs_dataset_{description}_{cdate.year}{cdate.month:02d}{cdate.day:02d}.csv", index=False)

    # Split based on score
    #split_basedon_score(fullSeqs_indexed_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run get_datasets script to load and concatenate TCR datasets.')
    parser.add_argument("-i", "--datasets_list", nargs='+', help="choose datasets to include", type=str)
    parser.add_argument("-f", "--datasets_dir", help="location of public datasets", type=str, default='../databases')
    parser.add_argument("-m", "--map_datasets", help="map dataset file names", type=json.loads, 
                        default='{"vdjdb": "vdjdb_human_monkey_mouse_TRA_TRB_paired_MHCI_score0_110326.tsv", \
                                "mcpas": "McPAS_TCR_updatedSep2022_10062025.csv", \
                                "iedb": "iedb_raw_linear_onlypositive_MHCI_TCRpaired.csv", \
                                "neo": "NeoTCR data-20221220.xlsx", \
                                "trait": "20250312-TRAIT_search_download.xlsx"}')
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

    df, dataset_file, location_db = get_datasets(datasets_list = args.datasets_list,
                                    location  = args.datasets_dir,
                                    map_datasets  = args.map_datasets,
                                    description   = args.description)

    preparation(dataset_file = dataset_file,
                thimble_path = args.thimble_path,
                ref_mhc      = args.ref_mhc,
                aligned_mhc  = args.aligned_mhc,
                description  = args.description,
                location_db   = location_db)

#python3 process_sequences.py -i vdjdb --outdir /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb