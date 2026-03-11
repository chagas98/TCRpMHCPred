#!/usr/bin/env python3
import pandas as pd
import os
import subprocess
import shutil
import pyMatchImm.models.evaluation as eval
from pyMatchImm.benchmark.models  import NetTCRrun, TULIPrun


class BenchmarkModels:
    """Class to run different models for benchmarking."""

    def __init__(self, test_readable: pd.DataFrame = None):

        self.test_readable_path = test_readable
        self.models = {
            'NetTCR': NetTCRrun,
            'TULIP': TULIPrun
        }
        self.outputs_list = {}
        self.processed_length = {}

    def run_model(self, model_name,
                      model_dir='trained',
                      keep_regions=True, keep_gaps=False, tmp_path='tmp/',
                      description='benchmark', scheme='aho'):
        
        if model_name not in self.models:
            raise ValueError(f"Model {model_name} is not supported.")
        
        # get this file path
        base_path = os.path.dirname(os.path.abspath(__file__))
        match model_name:
            case 'NetTCR':
                selections_include_path = os.path.join(base_path, 'assets', 'NetTCR22_include.json')
                selections_exclude_path = os.path.join(base_path, 'assets', 'NetTCR22_exclude.json')
            case 'TULIP':
                selections_include_path = os.path.join(base_path, 'assets', 'TULIP_include.json')
                selections_exclude_path = os.path.join(base_path, 'assets', 'TULIP_exclude.json')
            case _:
                raise ValueError(f"Model {model_name} is not supported.")
        
        modelRun = self.models[model_name]
        runner   = modelRun(
                input_path              = self.test_readable_path,
                selections_include_path = selections_include_path,
                selections_exclude_path = selections_exclude_path,
                keep_regions            = keep_regions,
                keep_gaps               = keep_gaps,
                model_dir               = model_dir,
                tmp_path                = tmp_path,
                description             = description,
                scheme                  = scheme
            )
        # Run the model
        self.outputs_list[model_name] = runner()
        self.processed_length[model_name] = runner.processed_length

    
    def evaluation(self, experiment_name: str = '', test_set_name: str = '', threshold: int = 0.5, description: str = 'benchmark'):
        """Evaluate the models' predictions."""
        if not self.outputs_list:
            raise ValueError("No model outputs to evaluate.")
        
        results = {}
        for model_name, predictions in self.outputs_list.items():
            
            if model_name not in ['TULIP']:
                # Perform evaluation using the evaluation module
                eval_results = eval.results_analysis(df_predictions  = predictions,
                                                    experiment_name = experiment_name, 
                                                    test_set_name   = test_set_name,
                                                    threshold       = threshold,
                                                    description     = description,
                                                    model_type      = model_name,
                                                    save_dir        = f"{model_name}model")
            elif model_name == 'TULIP':
                # For TULIP, we assume the predictions are already in the correct format
                eval_results = eval.results_analysis_tulip(df_predictions  = predictions,
                                                           experiment_name = experiment_name, 
                                                           test_set_name   = test_set_name,
                                                           description     = description,
                                                           model_type      = model_name,
                                                           df_length       = self.processed_length[model_name],
                                                           save_dir        = f"{model_name}model")

            else:
                raise ValueError(f"Model {model_name} is not supported for evaluation.")    
                
            results[model_name] = eval_results

            # calculate auc_per_peptide
            

        results_df = pd.concat(results.values(), keys=results.keys())
        results_df.reset_index(level=0, drop=True, inplace=True)

        results_df.to_csv(f'results_benchmark_{experiment_name}.csv', index=False)

        return results_df
    
    

