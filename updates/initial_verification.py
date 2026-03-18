import pandas as pd

df_old = pd.read_csv("/home/samuel.assis/MatchImm/TCRpMHCPred/databases/VDJdb/02_processed/combined_dataset_20251021.csv")
df_new = pd.read_csv("/home/samuel.assis/MatchImm/TCRpMHCPred/databases/VDJdb/02_processed/combined_dataset_20260311.csv")


score_mismatches = []
for i, row_old in df_old.iterrows():
    tcr_id_old = row_old['TCR_ID']
    cdr3a_old = row_old['CDR3A']
    cdr3b_old = row_old['CDR3B']
    va_old = row_old['VA']
    vb_old = row_old['VB']
    peptide_old = row_old['epitope']
    reference_old = row_old['Reference']

    matches = df_new[(df_new['CDR3A'] == cdr3a_old) &
                     (df_new['CDR3B'] == cdr3b_old) &
                     (df_new['VA'] == va_old) &
                     (df_new['VB'] == vb_old) &
                     (df_new['epitope'] == peptide_old) &
                     (df_new['Reference'] == reference_old)]
    if len(matches) == 0:
        print(f"No match found for TCR_ID {tcr_id_old} in the new dataset.")

    elif len(matches) > 1:
        print(f"Multiple matches found for TCR_ID {tcr_id_old} in the new dataset.")
    
    else:
        row_new = matches.iloc[0]
        tcr_id_new = row_new['TCR_ID']

        if row_new['Score'] != row_old['Score']:
            score_mismatches.append((tcr_id_old, tcr_id_new, row_old['Score'], row_new['Score'], row_old['Reference']))
            print(f"Score mismatch for TCR_ID {tcr_id_old}: old score {row_old['Score']} vs new score {row_new['Score']}")

#save
score_mismatches_df = pd.DataFrame(score_mismatches, columns=['TCR_ID_old', 'TCR_ID_new', 'Score_old', 'Score_new', 'Reference_old'])
score_mismatches_df.to_csv("mismatches_score.csv", index=False)