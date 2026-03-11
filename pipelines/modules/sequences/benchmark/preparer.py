#!/usr/bin/env python3
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import json
import os
import shutil
from pyMatchImm.preparation.indexer import ANARCIndexer
from pyMatchImm.preparation.featurization import GetPositionsTCR


class DatasetPreparer:
    def __init__(self, input_path:str):
        
        self.input_df = pd.read_csv(input_path)

    def TCRseq(self, selections_include_path: str, selections_exclude_path: str,
               keep_regions: bool = True, keep_gaps: bool = False, 
               tmp_path: str = 'tmp/', description: str = 'benchmark', scheme: str = 'aho'):
    
        indexer = ANARCIndexer(input_df                = self.input_df,
                               selections_include_path = selections_include_path,
                               selections_exclude_path = selections_exclude_path,
                               tmp_path                = tmp_path,
                               description             = description,
                               scheme                  = scheme)
        indexer.run()

        run_positionsTCR = GetPositionsTCR(dataset            = self.input_df,
                                           contact_to_include = indexer.positions_selected,
                                           ref_anarci         = indexer.anarci_positions)

        run_positionsTCR.get_indexed_seq_vector('TRA', keep_regions=True, keep_gaps=False)
        run_positionsTCR.get_indexed_seq_vector('TRB', keep_regions=True, keep_gaps=False)
        self.features_df = run_positionsTCR.save_dataset()
        
        print("Dataset prepared with features:", self.features_df.columns.tolist())

        return self.features_df
    
    