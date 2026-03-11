#!/usr/bin/env python3
import pandas as pd
import os
import subprocess
import shutil
from pyMatchImm.benchmark.preparer import DatasetPreparer
import time

#NetTCR allows just 12-mer peptides
# and uses a specific scheme for TCR sequences, so we adapt the DatasetPreparer accordingly.

class NetTCRrun(DatasetPreparer):
    def __init__(self, input_path: str, selections_include_path: str, selections_exclude_path: str,
                 keep_regions: bool = True, keep_gaps: bool = False, model_dir: str = 'NetTCR-2.2',
                 tmp_path: str = 'tmp/', description: str = 'benchmark', scheme: str = 'imgt', 
                 env_name: str = 'nettcr_2_2_env'):
        
        self.input_path = input_path
        self.selections_include_path = selections_include_path
        self.selections_exclude_path = selections_exclude_path
        self.tmp_path = tmp_path
        self.description = description
        self.scheme = scheme
        self.env_name = env_name
        self.keep_regions = keep_regions
        self.keep_gaps = keep_gaps
        self.model_dir = model_dir

        # Initialize the base class (DatasetPreparer) with the input path
        os.makedirs(tmp_path, exist_ok=True)
        super().__init__(input_path)

    def _get_dataset_(self):
        print("Running NetTCRrun...")
        self.dataset = self.TCRseq(selections_exclude_path = self.selections_exclude_path,
                                   selections_include_path = self.selections_include_path,
                                   keep_regions            = self.keep_regions,
                                   keep_gaps               = self.keep_gaps,
                                   tmp_path                = self.tmp_path,
                                   description             = self.description,
                                   scheme                  = self.scheme)
        print("NetTCRrun completed successfully.")

        #select columns
        self.dataset = self.dataset[['TCR_ID', 'A1','A2','A3','B1','B2','B3','peptide','MHCa', 'class']]
        self.dataset.rename(columns={'MHCa': 'allele',
                                     'class': 'Target'}, inplace=True)

        self.dataset.to_csv(os.path.join(self.tmp_path, 'nettcr_dataset.csv'), index=False)
        print("Selected columns:", self.dataset.columns.tolist())
        print("Dataset shape:", self.dataset.shape)

    def netTCRrun(self):
        
        # Ensure the dataset is prepared
        self._get_dataset_()
        
        # define paths for the model and dataset
        model_name = "t.4.v.3"
        model_type = "pan"
        predict_path = os.path.join(self.model_dir, 'src', 'predict.py')
        dataset_path = os.path.join(self.tmp_path, 'nettcr_dataset.csv')
        outdir_path = os.path.join(self.model_dir, 'models', 'nettcr_2_2_pan')        
        outfile_path = os.path.join(outdir_path, f"{model_name}_prediction.csv") # fixed output file of nettcr
        
        # Ensure output directory exists
        self.dataset.to_csv(dataset_path, index=False)
        
        # check paths
        assert os.path.exists(predict_path), "Predict script not found at {}".format(predict_path)
        assert os.path.exists(dataset_path), "Dataset file not found at {}".format(dataset_path)
        assert os.path.exists(outdir_path), "Output directory not found at {}".format(outdir_path)

        # Run the NetTCR prediction script
        command = [
            "conda", "run", "-n", self.env_name,
            'python3', predict_path,
            '--test_data', dataset_path,
            '--model_name', model_name,
            '--model_type', model_type,
            '--outdir', outdir_path
        ]
        
        # Execute the command and capture its output
        print("Running command:", ' '.join(command))
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        print("Command Output:\n", stdout)
        if stderr:
            print("Command Errors:\n", stderr)
        
        if process.returncode != 0:
            raise RuntimeError(f"NetTCR-2.2 prediction failed with return code {process.returncode}")
        else:
            print("NetTCR-2.2 prediction completed successfully.")

            while not os.path.exists(outfile_path):
                time.sleep(3)

        # Save the results to a new directory
        self.results = pd.read_csv(outfile_path)
        self.results.rename(columns={'class': 'Target'}, inplace=True)
        
        os.makedirs('NetTCRmodel', exist_ok=True)
        new_path = os.path.join('NetTCRmodel', "predictions_NetTCR.csv")
        self.results.to_csv(new_path, index=False)
        shutil.rmtree(self.tmp_path)

        print("NetTCR prediction completed successfully.")

    @property
    def processed_length(self):
        """Return the length of the processed dataset."""
        return len(self.dataset)
    
    def __call__(self):
        """Execute the NetTCR prediction and return the results."""
        self.netTCRrun()
        return self.results

        

class TULIPrun(DatasetPreparer):
    def __init__(self, input_path: str, selections_include_path: str, selections_exclude_path: str,
                 keep_regions: bool = True, keep_gaps: bool = False, model_dir: str = 'TULIP',
                 tmp_path: str = 'tmp/', description: str = 'benchmark', scheme: str = 'aho', 
                 env_name: str = 'tulip_env', outdir_path: str = '.'):
        
        print("Initializing TULIPrun...")
        self.input_path              = input_path
        self.selections_include_path = selections_include_path
        self.selections_exclude_path = selections_exclude_path
        self.tmp_path                = tmp_path
        self.description             = description
        self.scheme                  = scheme
        self.env_name                = env_name
        self.model_dir               = model_dir
        self.keep_regions            = keep_regions
        self.keep_gaps               = keep_gaps
        self.outdir_path             = outdir_path

        # Initialize the base class (DatasetPreparer) with the input path
        os.makedirs(tmp_path, exist_ok=True)
        super().__init__(input_path)
            
    def _get_dataset_(self):
        print("Running TULIPrun...")
        self.dataset = self.TCRseq(selections_exclude_path = self.selections_exclude_path,
                                   selections_include_path = self.selections_include_path,
                                   keep_regions            = self.keep_regions,
                                   keep_gaps               = self.keep_gaps,
                                   tmp_path                = self.tmp_path,
                                   description             = self.description,
                                   scheme                  = self.scheme)
        print("TULIPrun completed successfully.")

        #select columns
        self.dataset = self.dataset[['TCR_ID', 'CDR3a', 'CDR3b','MHCa', 'peptide', 'class']]
        self.dataset.rename(columns={'MHCa': 'MHC',
                                     'class': 'binder'}, inplace=True)
        self.dataset.to_csv(os.path.join(self.tmp_path, 'tulip_dataset.csv'), index=False)
        print("Selected columns:", self.dataset.columns.tolist())

    def TULIPrun(self):
        
        # Ensure the dataset is prepared
        self._get_dataset_()
        os.makedirs('TULIPmodel', exist_ok=True)
        
        predict_path = os.path.join(self.model_dir, 'predict.py')
        dataset_path = os.path.join('TULIPmodel', 'tulip_dataset.csv')
        
        # Ensure output directory exists
        self.dataset.to_csv(dataset_path, index=False)
        
        # check paths
        assert os.path.exists(predict_path), "Predict script not found at {}".format(predict_path)
        assert os.path.exists(dataset_path), "Dataset file not found at {}".format(dataset_path)
        assert os.path.exists('TULIPmodel'), "Output directory not found at {}".format('TULIPmodel')

        # Run the NetTCR prediction script
        command = [
            "conda", "run", "-n", self.env_name,
            'python3', predict_path,
            '--test_dir', dataset_path,
            '--output', 'TULIPmodel/'
        ]

        # Execute the command and capture its output
        print("Running command:", ' '.join(command))
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        print("Command Output:\n", stdout)
        if stderr:
            print("Command Errors:\n", stderr)
        
        if process.returncode != 0:
            raise RuntimeError(f"TULIP prediction failed with return code {process.returncode}")
        else:
            print("TULIP prediction completed successfully.")
            # Wait until auc_scores.csv exists
            auc_file = os.path.join('TULIPmodel', 'auc_scores.csv')
            while not os.path.exists(auc_file):
                time.sleep(3)
            
        # Save the results to a new directory
        self.results = pd.read_csv(os.path.join('TULIPmodel', 'auc_scores.csv'))
        
        os.rename(auc_file, os.path.join('TULIPmodel', 'predictions_TULIP.csv'))
        shutil.rmtree(self.tmp_path)

        print("TULIP prediction completed successfully.")
    
    @property
    def processed_length(self):
        """Return the length of the processed dataset."""
        return len(self.dataset)
    
    def __call__(self):
        """Run the TULIP prediction and return the results."""
        self.TULIPrun()
        return self.results
        