# pMHCs diferentes tiveram numeração duplicada e o slurm array coletou apenas um dos arquivos fasta,
# então esse script exclui pMHC que não foram preditos

# - verifica quais pMHC foram preditos (score2_wojust10x) e exclui os que não foram preditos do score2_wojust10x.csv
# - gera fastas para com os pMHC não preditos e com nova numeração (continuando a numeração dos pMHCs preditos) para serem preditos no score2_wojust10x

import pandas as pd
import os
from collections import Counter

pMHC_fastas_dir_path = "/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score2_wojust10x/04_assets/pMHC_fastas"
predictions_path = "/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score2_wojust10x/02_predictions/pMHC/pMHC_raw"
all_pMHC_path = "/home/samuel.assis/MatchImm/TCRpMHCPred/runs/vdjdb_score2_wojust10x/01_identification/unique_pmhc_samples_20251205.csv"

all_pMHC = pd.read_csv(all_pMHC_path)['unique_pmhc_id'].tolist()
pMHC_fasta_list = [f.split('/')[-1].split('.')[0] for f in os.listdir(pMHC_fastas_dir_path) if f.endswith('.fasta')]
predictions_list = [f.split('/')[-1] for f in os.listdir(predictions_path)]

# Verificando pMHC
print(len(all_pMHC))
print(len(pMHC_fasta_list))
print(len(predictions_list))

#Verificando pMHC preditos duplicados
predictions_alleles = ['_'.join(f.split('_')[0:-1]) for f in predictions_list]
predictions_alleles_count = Counter(predictions_alleles) # count
duplicated_predictions_alleles = {allele: count for allele, count in predictions_alleles_count.items() if count > 1} # filter only duplicated

print(f"pMHC preditos duplicados (sem considerar numeração): {len(duplicated_predictions_alleles)}")

#Verificando pMHC fastas duplicados
pMHC_fasta_alleles = ['_'.join(f.split('_')[0:-1]) for f in pMHC_fasta_list]
pMHC_fasta_alleles_count = Counter(pMHC_fasta_alleles) # count
duplicated_pMHC_fasta_alleles = {allele: count for allele, count in pMHC_fasta_alleles_count.items() if count > 1} # filter only duplicated

print(f"pMHC fastas duplicados (sem considerar numeração): {len(duplicated_pMHC_fasta_alleles)}")

#Verificando fastas que nao tem pMHC predito
fastas_without_predictions_indexed = set(pMHC_fasta_list) - set(predictions_list)
print(f"Fastas sem predição: {len(fastas_without_predictions_indexed)}")

fastas_without_predictions = set(pMHC_fasta_alleles) - set(predictions_alleles)
print(f"Fastas sem predição (sem considerar numeração): {len(fastas_without_predictions)}")


# deduplicate 
fastas_without_predictions_alleles = set(['_'.join(f.split('_')[0:-1]) for f in fastas_without_predictions_indexed])
print("Alelos duplicados nos fastas e que não tem predição (sem considerar numeração):")
alleles_to_predict = set(fastas_without_predictions_alleles) - set(predictions_alleles) # check if any of the fastas without predictions (without considering numbering) are not in the predictions (without considering numbering)

alleles_to_predict = {'B08010101_MAALPRLIAF', 'B08010101_QAKWRLQTL', 'E01030101_VMTTVLATL', 'B2701_LRVMMLAPF', 'A02010101_LLMEGVPKSL', 'A02010101_VLFGLGFAI', 'B08010101_NLIFKWIL', 'A29020101_EVLPFFLFF', 'A02010101_VVMSWAPPV', 'A01010101_HSNLNDATY', 'B08010101_ILKGKFQTA', 'A02010101_CINGVCWTV', 'A02010101_CLGGLLTMV', 'B35080101_LPEPLPQGQGTAY', 'A02010101_LVGPTPVNI', 'A02010101_AIIRILQQL', 'A02010101_GMDYHNGHL', 'A03010101_AAFKRSCLK', 'A11010101_AVFDRKSDAK', 'A02010101_VLDLFQGQL', 'A02010101_KIVALGINAV', 'A02010101_ALTAVAEEV', 'A02010101_YLYDRLLRV', 'A02010101_GVYALIAGA', 'B08010101_TLKKMREI', 'A02010101_HLVEALYLV', 'B08010101_RAKFKQLL', 'A02010101_ALWMRLLPLL', 'C03030101_FVVPYMIYLL', 'A01010101_YSEHPTFTSQY', 'B44030101_TEDEHFEFY', 'A02010101_CVNGSCFTV', 'B08010101_LLNVKLAL', 'B08010101_YAYAKWKL', 'B08010101_ELRRKMMYM', 'A02010101_VMNILLQYV', 'A03010101_WLIRETQPITK', 'B44020101_GEEDGAGGHSL', 'B08010101_HPRYFNQL', 'A02010101_YLLPRRGPRL', 'A02010101_RLLPLLALL', 'A02010101_YLLEMLWRL', 'B35080101_LPEPLPQGQLTGY', 'A02010101_ILSAHVATA', 'A02010101_LALWGPDPAA', 'A01010101_ALASCMGLIY', 'B08010101_WMRLLPLL', 'B07020101_GARGVGKSAL'}

#create new fastas to predict
ref_pmhc_list = pd.read_csv(all_pMHC_path)
for index, row in ref_pmhc_list.iterrows():

    if row['unique_pmhc_id'] in alleles_to_predict:
        
        fasta_path = os.path.join(pMHC_fastas_dir_path, f"{row['unique_pmhc_id']}_new_{index:03d}.fasta")
        with open(fasta_path, 'w') as fasta_file:
            fasta_file.write(f">{row['unique_pmhc_id']}\n{row['MHCseq']}:{row['peptide']}\n")
        fasta_file.close()


