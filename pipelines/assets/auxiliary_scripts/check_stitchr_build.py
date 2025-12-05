#!/usr/bin/env python3
import os
import pandas as pd


def stitchr_commands(input_data):
    df = pd.read_csv(input_data)
    df = df[df['Score'] == 2]
    for index, row in df.iterrows():
        tcr_id = row['TCR_ID']
        score = row['Score']
        tra_cdr3 = row['CDR3A']
        trb_cdr3 = row['CDR3B']

        tra_v_gene = row['VA']
        tra_j_gene = row['JA']
        trb_v_gene = row['VB']
        trb_j_gene = row['JB']
        print(f"# TCR ID: {tcr_id} | Score: {score}")
        print(f"stitchr  --v {tra_v_gene} --j {tra_j_gene} --cdr3 {tra_cdr3}")
        print(f"stitchr  --v {trb_v_gene} --j {trb_j_gene} --cdr3 {trb_cdr3}")

if __name__ == "__main__":
    input_data = "/home/samuel.assis/MatchImm/TCRpMHCPred/databases/VDJdb/02_processed/combined_dataset_20251021.csv"
    stitchr_commands(input_data)