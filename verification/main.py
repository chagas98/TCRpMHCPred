from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Dict, Any
import pandas as pd
from Bio.PDB import PDBParser, Polypeptide

# DS_train = final dataset with all data


@dataclass
class SummaryChecker:
    """Container for comparison outputs."""
    id: str
    duplicated_ids: bool
    missing_ids: bool
    CDR3_existence: bool
    peptide_matching: bool
    MHCa_matching: bool
    MHCseq_matching: bool
    pMHC_structure_matching: bool
    TCR_structure_matching: bool


def extract_pdb_sequences(pdb_file: str) -> dict[str, str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    sequences = {}

    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                # Keep only standard amino acids
                if Polypeptide.is_aa(residue, standard=True):
                    try:
                        seq.append(Polypeptide.protein_letters_3to1[residue.resname])
                    except KeyError:
                        pass
            if seq:
                sequences[chain.id] = "".join(seq)
        break  # use only first model

    return sequences



class DataVerifier:
    def __init__(
    self,
    DS_train: pd.DataFrame,
    DS_raw: pd.DataFrame,
    pMHC_struct_path: str,
    TCR_struct_path: str,
    MHC_refseq_path: str,
    id_col: str = 'TCR_ID') -> None:
    
        self.DS_train = DS_train.copy()
        self.DS_raw = DS_raw.copy()
        self.pMHC_struct_path = pMHC_struct_path
        self.TCR_struct_path = TCR_struct_path
        self.MHC_refseq_path = MHC_refseq_path
        self.id_col = id_col


    def check_duplicates(self) -> bool:
        """Check for duplicate IDs in the training dataset."""
        duplicated_train = self.DS_train[self.id_col].duplicated().any()
        duplicated_raw = self.DS_raw[self.id_col].duplicated().any()
        if duplicated_train or duplicated_raw:
            print(f"Duplicate IDs found in {'training' if duplicated_train else 'raw'} dataset.")
            return False
        return True
    
    def get_id_sets(self) -> tuple[List[str], List[str]]:
        """Return unique ID sets from both datasets."""
        train_ids = self.DS_train[self.id_col].dropna().unique()
        raw_ids = self.DS_raw[self.id_col].dropna().unique()
        return train_ids, raw_ids

    def check_cdr3_existence(self, train_ids, raw_ids) -> bool:
        """Check if CDR3 sequences exist for all IDs in both datasets."""
        train_cdr3 = self.DS_train[self.DS_train[self.id_col].isin(train_ids)][['TCR_ID', 'TRA', 'TRB']]
        raw_cdr3 = self.DS_raw[self.DS_raw[self.id_col].isin(raw_ids)][['TCR_ID', 'CDR3A', 'CDR3B']]
        
        merged = pd.merge(train_cdr3, raw_cdr3, on='TCR_ID', suffixes=('_train', '_raw'))
        maskA = merged.apply(lambda row: str(row['CDR3A']) in str(row['TRA']), axis=1)
        maskB = merged.apply(lambda row: str(row['CDR3B']) in str(row['TRB']), axis=1)
        missing_cdr3 = merged[~(maskA & maskB)]
        if not missing_cdr3.empty:
            print(f"Missing CDR3 sequences found for IDs: {missing_cdr3['TCR_ID'].tolist()}")
            return False
        return True

    def check_peptide_matching(self, train_ids, raw_ids) -> bool:
        """Check if peptides in the training dataset match those in the raw dataset."""
        train_peptides = self.DS_train[self.DS_train[self.id_col].isin(train_ids)][['TCR_ID', 'peptide']]
        raw_peptides = self.DS_raw[self.DS_raw[self.id_col].isin(raw_ids)][['TCR_ID', 'epitope']]
        merged = pd.merge(train_peptides, raw_peptides, on='TCR_ID', suffixes=('_train', '_raw'))
        mismatches = merged[merged['peptide'] != merged['epitope']]
        if not mismatches.empty:
            print(f"Peptide mismatches found for IDs: {mismatches['TCR_ID'].tolist()}")
            return False
        return True
    
    def check_MHCa_matching(self, train_ids, raw_ids) -> bool:
        """Check if MHC alleles in the training dataset match those in the raw dataset."""
        train_MHCa = self.DS_train[self.DS_train[self.id_col].isin(train_ids)][['TCR_ID', 'MHCa', 'assigned_allele']]
        raw_MHCa = self.DS_raw[self.DS_raw[self.id_col].isin(raw_ids)][['TCR_ID', 'MHCa']]
        train_MHCa['MHCa'] = train_MHCa['MHCa'].str.replace('HLA-', '', regex=False)
        raw_MHCa['MHCa'] = raw_MHCa['MHCa'].str.replace('HLA-', '', regex=False)
        
        merged = pd.merge(train_MHCa, raw_MHCa, on='TCR_ID', suffixes=('_train', '_raw'))
        
        mask = merged.apply(
            lambda row: str(row["MHCa_raw"]) in str(row["assigned_allele"]),
            axis=1
        )
        mismatches_assigned = merged[~mask]

        mismatches_allele = merged[merged['MHCa_train'] != merged['MHCa_raw']]
        if not mismatches_allele.empty:
            print(f"MHC allele mismatches found for IDs: {mismatches_allele['TCR_ID'].tolist()}")
            return False
        
        if not mismatches_assigned.empty:
            print(f"Assigned allele mismatches found for IDs: {mismatches_assigned['TCR_ID'].tolist()}")
            return False
        
        return True
    
    def get_mhcrefseq(self) -> Dict[str, str]:
        from Bio import SeqIO
        mhc_dict = {}   
        for seq_record in SeqIO.parse(self.MHC_refseq_path, "fasta"):
            mhc_dict[seq_record.id] = str(seq_record.seq)

        return mhc_dict

    def check_MHCseq_matching(self, train_ids, raw_ids) -> bool:
        """Check if MHC sequences in the training dataset match those in the raw dataset."""
        train_MHCseq = self.DS_train[self.DS_train[self.id_col].isin(train_ids)][['TCR_ID', 'assigned_allele', 'MHCseq']]
        mismatches_Ids = []
        for idx, row in train_MHCseq.iterrows():
            allele = row['assigned_allele']
            mhc_seq = row['MHCseq']
            if pd.isna(allele) or pd.isna(mhc_seq):
                continue
            ref_seq = self.mhc_refseq_dict.get(allele)
            if ref_seq is None:
                print(f"Reference sequence not found for allele: {allele}")
                continue
            if mhc_seq != ref_seq:
                print(f"MHC sequence mismatch for ID {row['TCR_ID']}: expected {ref_seq}, got {mhc_seq}")
                mismatches_Ids.append(row['TCR_ID'])
        if mismatches_Ids:
            print(f"MHC sequence mismatches found for IDs: {mismatches_Ids}")
            return False
        return True
    
    def check_pMHC_structure_matching(self, train_ids, raw_ids) -> bool:
        """Check if pMHC structures in the training dataset match those in the raw dataset."""
        train_pMHC = self.DS_train[self.DS_train[self.id_col].isin(train_ids)][['TCR_ID', 'peptide', 'MHCseq', 'filepath_b']]
        raw_pMHC = self.DS_raw[self.DS_raw[self.id_col].isin(raw_ids)][['TCR_ID', 'epitope']]

        mismatches_ids = []
        for idx, row in train_pMHC.iterrows():
            train_id = row['TCR_ID']
            train_peptide = row['peptide']
            train_MHCseq = row['MHCseq']
            train_filepath = row['filepath_b']

            if pd.isna(train_filepath):
                continue

            struct_sequences = extract_pdb_sequences(train_filepath)
            if 'A' not in struct_sequences or 'C' not in struct_sequences:
                print(f"Missing chains in structure for ID {train_id}")
                continue
            struct_peptide = struct_sequences['C']
            struct_MHCseq = struct_sequences['A']

            if train_peptide != struct_peptide:
                print(f"Peptide mismatch for ID {train_id}: expected {train_peptide}, got {struct_peptide}")
                mismatches_ids.append(train_id)
            if train_MHCseq != struct_MHCseq:
                print(f"MHC sequence mismatch for ID {train_id}: expected {train_MHCseq}, got {struct_MHCseq}")
                mismatches_ids.append(train_id)

        if mismatches_ids:
            print(f"pMHC structure mismatches found for IDs: {mismatches_ids}")
            return False
        return True

    def check_TCR_structure_matching(self, train_ids, raw_ids) -> bool:
        """Check if TCR structures in the training dataset match those in the raw dataset."""
        train_TCR = self.DS_train[self.DS_train[self.id_col].isin(train_ids)][['TCR_ID', 'TRA', 'TRB', 'filepath_a']]

        for idx, row in train_TCR.iterrows():
            train_id = row['TCR_ID']
            train_TRA = row['TRA']
            train_TRB = row['TRB']
            train_filepath = row['filepath_a']

            if pd.isna(train_filepath):
                continue

            struct_sequences = extract_pdb_sequences(train_filepath)
            if 'D' not in struct_sequences or 'E' not in struct_sequences:
                print(f"Missing chains in structure for ID {train_id}")
                continue
            struct_TRA = struct_sequences['D']
            struct_TRB = struct_sequences['E']

            if train_TRA != struct_TRA:
                print(f"TRA mismatch for ID {train_id}: expected {train_TRA}, got {struct_TRA}")
                return False
            if train_TRB != struct_TRB:
                print(f"TRB mismatch for ID {train_id}: expected {train_TRB}, got {struct_TRB}")
                return False
            return True
    

    def run(self) -> None:
        train_ids, raw_ids = self.get_id_sets()
        self.mhc_refseq_dict = self.get_mhcrefseq()
        if not self.check_duplicates():
            print("Duplicate IDs found. Please resolve duplicates before proceeding.")
            return
        if not self.check_cdr3_existence(train_ids, raw_ids):
            print("Missing CDR3 sequences found. Please resolve missing CDR3s before proceeding.")
            return
        if not self.check_peptide_matching(train_ids, raw_ids):
            print("Peptide mismatches found. Please resolve mismatches before proceeding.")
            return
        if not self.check_MHCa_matching(train_ids, raw_ids):
            print("MHC allele mismatches found. Please resolve mismatches before proceeding.")
            return
        if not self.check_MHCseq_matching(train_ids, raw_ids):
            print("MHC sequence mismatches found. Please resolve mismatches before proceeding.")
            return
        if not self.check_pMHC_structure_matching(train_ids, raw_ids):
            print("pMHC structure mismatches found. Please resolve mismatches before proceeding.")
            return
        if not self.check_TCR_structure_matching(train_ids, raw_ids):
            print("TCR structure mismatches found. Please resolve mismatches before proceeding.")
            return
        print("All checks passed successfully.")
        
if __name__ == "__main__":
    
    DS_raw = pd.read_csv('/home/samuel.assis/MatchImm/TCRpMHCPred/databases/VDJdb/01_raw/vdjdb_HomoSapiens.csv')
    
    ref_mhc_path = '/home/samuel.assis/MatchImm/Public_Datasets/mhc/A02010101_aligned_trimmed24_180_addGinit.fasta'

    DS_train_score2 = pd.read_csv('/home/samuel.assis/MatchImm/MatchImmNet/data/01-raw/AF_vdjdb_score2_wojust10x_20260307.csv')
    TCR_struct_path_score2 = '/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score2_wojust10x/02_predictions/TCR/TCR_raw'
    pMHC_struct_path_score2 = '/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score2_wojust10x/02_predictions/pMHC/pMHC_renum'
    verifier = DataVerifier(DS_train_score2, DS_raw, pMHC_struct_path_score2, TCR_struct_path_score2, ref_mhc_path)
    verifier.run()


    DS_train_score3 = pd.read_csv('/home/samuel.assis/MatchImm/MatchImmNet/data/01-raw/AF_vdjdb_score3_20251212.csv')
    TCR_struct_path_score3 = '/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score_3/02_predictions/TCR/TCR_raw'
    pMHC_struct_path_score3 = '/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score_3/02_predictions/pMHC/pMHC_renum'
    verifier = DataVerifier(DS_train_score3, DS_raw, pMHC_struct_path_score3, TCR_struct_path_score3, ref_mhc_path)
    verifier.run()