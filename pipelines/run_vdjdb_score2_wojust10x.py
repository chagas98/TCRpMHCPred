#!/usr/bin/env python3

from modules.TCRpreparation import TCRpreparation, SlurmTCRModel2 
from modules.pMHCpreparation import pMHCpreparation, SlurmColabFold
from modules.renumber_pMHC import renumbering_pMHC
from modules.gather_results import find_tcr_models, find_pmhc_models, gather_pmhc_models
import os
import pandas as pd
import argparse as arg
from datetime import datetime


# Get current date for file naming
cdate = datetime.now()
DATE_STRING = f"{cdate.year}{cdate.month:02d}{cdate.day:02d}"
BASEPATH = "/home/samuel.assis/MatchImm/TCRpMHCPred"

class RunResources:
    def __init__(self, name:str, pipeline_name:str, input_path:str, processed_path:str, stage:str, direct_path:str = ""):
        self.name = name
        self.pipeline_name = pipeline_name
        self.input_path = input_path
        self.processed_path = processed_path
        self.direct_path = direct_path

        # Create Run environment directories
        self.runs_dir     = f"{BASEPATH}/runs/{self.pipeline_name}"
        self.id_path      = f"{self.runs_dir}/01_identification"
        self.pred_path    = f"{self.runs_dir}/02_predictions"
        self.results_path = f"{self.runs_dir}/03_results"
        self.assets_path  = f"{self.runs_dir}/04_assets"
        self.logs_path    = f"{self.runs_dir}/05_logs"

        if stage == "prep":
            self.set_resources()
        
        elif stage == "gather":
            self.data_path = self.direct_path

    def set_resources(self):
        # Define paths and create directories
        self.data_path = self.processed_path

        os.makedirs(self.runs_dir, exist_ok=True)
        os.makedirs(self.id_path, exist_ok=True)
        os.makedirs(self.pred_path, exist_ok=True)
        os.makedirs(self.results_path, exist_ok=True)
        os.makedirs(self.assets_path, exist_ok=True)
        os.makedirs(self.logs_path, exist_ok=True)

        # Remove just10x references from VDJdb scored dataset
        just10x_ref_path = "assets/immrep_just10x_exclusion.csv"

        raw_df = pd.read_csv(self.input_path)
        reference = pd.read_csv(just10x_ref_path)
        just10x_df = reference[reference['Exclude'] == 1]

        just10x_df['PMID'] = just10x_df['PMID'].replace({'PMID: ': ''})

        filtered_df = raw_df[~raw_df['Reference'].isin(just10x_df['PMID'])]
        filtered_df = filtered_df[~filtered_df['Reference'].isin(just10x_df['VDJdb_alias'])]
        filtered_df.reset_index(drop=True, inplace=True)
        filtered_df.to_csv(self.data_path, index=False)
        

        # statistics
        print(f"Original dataset size: {raw_df.shape[0]}")
        print(f"Filtered dataset size: {filtered_df.shape[0]}")

class GatherResults:
    def __init__(self, name:str, data_path:str, match_ids_tcr:str, match_ids_pmhc:str, run_dir:str, pmhc_pred_path:str, tcr_pred_path:str):
        self.name = name
        self.data_path = data_path
        self.match_ids_tcr = match_ids_tcr
        self.match_ids_pmhc = match_ids_pmhc
        self.run_dir = run_dir
        self.pmhc_pred_path = pmhc_pred_path
        self.tcr_pred_path = tcr_pred_path
        self.gather()

    def gather(self):
        print("Gathering results...")
        
        data_df = pd.read_csv(self.data_path)
        match_ids_tcr = pd.read_csv(self.match_ids_tcr)
        match_ids_pmhc = pd.read_csv(self.match_ids_pmhc)

        # merge
        combined_dataset = data_df.merge(
            match_ids_tcr, on='TCR_ID', how='left').merge(
                match_ids_pmhc, on='TCR_ID', how='left')

        # Find TCR models
        tcr_models_dir = os.path.join(self.tcr_pred_path, 'TCR_raw')
        data_tcr_df = find_tcr_models(combined_dataset, tcr_models_dir)

        # Find pMHC models
        pmhc_models_dir = os.path.join(self.pmhc_pred_path, 'pMHC_raw')
        pmhc_models_paths = gather_pmhc_models(pmhc_models_dir)
        
        # Renumber pMHC models
        pMHC_renum_path = os.path.join(self.pmhc_pred_path, 'pMHC_renum')
        renumbering_pMHC(pmhc_paths=pmhc_models_paths,
                         output_dir=pMHC_renum_path,
                         ref_pdb='assets/reference_1ao7.trunc.fit.pdb')

        data_all_df = find_pmhc_models(data_tcr_df, pMHC_renum_path)
        
        
        data_all_df.to_csv(f'AF_{name}_{DATE_STRING}.csv', index=False)


if __name__ == "__main__":

    parser = arg.ArgumentParser(description="Run TCRpMHCpred")
    parser.add_argument('--stage', type=str, default='gather', help='Stage of the pipeline to run')

    args = parser.parse_args()
    stage = args.stage
    pipeline_name="vdjdb_score2_wojust10x"
    name="score2_wojust10x"
    input_path=f"{BASEPATH}/databases/VDJdb/02_processed/fullSeqs_dataset_score2_20251021.csv"
    processed_path=f"{BASEPATH}/databases/VDJdb/02_processed/fullSeqs_dataset_{name}_{DATE_STRING}.csv"
    direct_path = f"{BASEPATH}/databases/VDJdb/02_processed/fullSeqs_dataset_score2_wojust10x_20251205.csv"
    date_to_gather = "20251205"
    
    print("Preparing data...")
    resources = RunResources(name=name, 
                                pipeline_name=pipeline_name, 
                                input_path=input_path,
                                processed_path=processed_path,
                                stage=stage,
                                direct_path=direct_path)
    if stage == "prep":

        # Run TCRmodel2
        TCRpreparation(
            input_path=resources.data_path,
            run_dir=resources.runs_dir)

        SlurmTCRModel2(
            pipeline_name=pipeline_name,
            run_dir=resources.runs_dir,
            input_model_path=os.path.join(resources.id_path, f'unique_tcr_samples_{DATE_STRING}.csv'))

        # Run ColabFold for pMHCs
        pMHCpreparation(
            input_path=resources.data_path,
            run_dir=resources.runs_dir)

        SlurmColabFold(
            pipeline_name=pipeline_name,
            run_dir=resources.runs_dir,
            input_model_path=os.path.join(resources.pred_path, f'unique_pmhc_samples_{DATE_STRING}.csv'),
            colabfold_path='/shared/groups/proteindesign/singularity_containers/colabfold_155.sif')

    if stage == "gather":
        print("Gathering results...")

        GatherResults(
            name=name,
            data_path=resources.data_path,
            match_ids_tcr=os.path.join(resources.id_path, f'match_ids_tcr_{date_to_gather}.csv'),
            match_ids_pmhc=os.path.join(resources.id_path, f'match_ids_pmhc_{date_to_gather}.csv'),
            run_dir=resources.runs_dir,
            pmhc_pred_path=os.path.join(resources.pred_path, 'pMHC'),
            tcr_pred_path=os.path.join(resources.pred_path, 'TCR'))
