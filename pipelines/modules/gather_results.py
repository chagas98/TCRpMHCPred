#!/usr/bin/env python3
import sys
import os
import pandas as pd
import json
import argparse
import glob
from datetime import datetime
from pathlib import Path

cdate = datetime.now()
DATE_STRING = f"{cdate.year}{cdate.month:02d}{cdate.day:02d}"

def find_tcr_models(input_df, path_to_tcr_models):
    tcr_models = []
    for index, row in input_df.iterrows():
        tcr_id = row['TCR_ID']
        unique_tcr_id = row['unique_tcr_id']
        model_file_path = os.path.join(path_to_tcr_models, unique_tcr_id, 'ranked_0.pdb')

        if model_file_path.endswith('.pdb') and os.path.isfile(model_file_path):
            row['filepath_a'] = os.path.abspath(model_file_path)
        else:
            row['filepath_a'] = None
        
        row['label'] = 1
        row['source'] = 'colabfold_pmhc/tcrmodel2_tcr'
        tcr_models.append(row)
    
    tcr_models_df = pd.DataFrame(tcr_models)
    return tcr_models_df

def gather_pmhc_models(path_to_pmhc_models, pattern_mhc = "*_relaxed_rank_001_alphafold2_multimer_v3_model_*"):

    source_dir = Path(path_to_pmhc_models)
    pattern = pattern_mhc

    matching_files = [f for f in source_dir.rglob(pattern) if f.is_file()]

    return matching_files

def find_pmhc_models(input_df,dir_path):
    pmhc_models = []
    for index, row in input_df.iterrows():
        unique_pmhc_id = row['unique_pmhc_id']
        model_name_pattern = os.path.join(dir_path, str(unique_pmhc_id) + '*')
        
        if not glob.glob(model_name_pattern):
            row['filepath_b'] = None
        else:
            model_file_path = glob.glob(model_name_pattern)[0]
            row['filepath_b'] = os.path.abspath(model_file_path)

        row['label'] = 1
        row['source'] = 'colabfold_pmhc/tcrmodel2_tcr'
        pmhc_models.append(row)


    pmhc_models_df = pd.DataFrame(pmhc_models)
    return pmhc_models_df



if __name__ == "__main__":
    #set the args
    parser = argparse.ArgumentParser(description='Run get_datasets script to load and concatenate TCR datasets.')
    parser.add_argument("-i", "--input", help="dataset of all samples", type=str)
    parser.add_argument("-t", "--match_ids_tcr", help="match ids tcr file", type=str)
    parser.add_argument("-p", "--match_ids_pmhc", help="match ids pmhc file", type=str)
    parser.add_argument("-l", "--path_to_tcr", help="location of tcr models", type=str)
    parser.add_argument("-r", "--path_to_pmhc", help="location of pmhc models", type=str)
    args = parser.parse_args()

    input_df = pd.read_csv(args.input)
    match_ids_tcr = pd.read_csv(args.match_ids_tcr)
    match_ids_pmhc = pd.read_csv(args.match_ids_pmhc)

    # merge
    combined_dataset = input_df.merge(match_ids_tcr, on='TCR_ID', how='left').merge(match_ids_pmhc, on='TCR_ID', how='left')

    dataset_df = find_tcr_models(combined_dataset, args.path_to_tcr)
    dataset_df = find_pmhc_models(dataset_df, args.path_to_pmhc)

    dataset_df.to_csv(f'AF_vdjdb_score3_{DATE_STRING}.csv', index=False)


#python3 AF_prepare_data.py -i /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb/datasets/fullSeqs_dataset_score3_20251021.csv -t /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb/score_3/assets/match_ids_20251022.csv -p /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb/score_3/assets/match_ids_mhc_20251022.csv -l /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb/score_3/TCR -r /home/samuel.assis/MatchImm/3_StructPred/predictions/VDJdb/score_3/pMHC_renum
