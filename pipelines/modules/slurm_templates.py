#!/usr/bin/env python3

import os

class SlurmTCRModel2:
    def __init__(self, 
                 output_dir:str, 
                 input_model_path:str, 
                 tcrmodel_path:str,
                 ignore_pdbs_string:str, 
                 relax_structures:str, 
                 max_template_date:str, 
                 databases_path:str):

        self.output_dir = output_dir
        self.input_model_path = input_model_path
        self.tcrmodel_path = tcrmodel_path
        self.ignore_pdbs_string = ignore_pdbs_string
        self.relax_structures = relax_structures
        self.max_template_date = max_template_date
        self.databases_path = databases_path

        os.makedirs(os.path.join(self.output_dir, '05_logs', 'TCR'), exist_ok=True)

        self.generate_slurm_template()

        assets_dir = os.path.join(self.output_dir, 'assets')
        with open(f'{assets_dir}/tcrmodel_sbatch.sh', 'w') as file:
            file.write(self.template)

    def generate_slurm_template(self, dir_slurm:str=None):

        self.template = f'''
#!/bin/bash
#SBATCH --array=1-800%1
#SBATCH --job-name=tcrmodel2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=short-gpu-big
#SBATCH --gres=gpu:a100:1
#SBATCH --output={self.output_dir}/05_logs/TCR/slurm.%A_%a.out
#SBATCH --error={self.output_dir}/05_logs/TCR/slurm.%A_%a.error

#
######## END SLURM HEADER
#

mkdir -p {self.output_dir}/TCR

echo $SLURM_ARRAY_TASK_ID

value=$((SLURM_ARRAY_TASK_ID + 0))

readarray -t table < <(tail -n +2 {self.input_model_path} | cut -d, -f1-4)

# assuming $value starts at 1
line="${{table[$((value-1))]}}"

sample_id=$(echo "$line" | cut -d, -f1)
TCRA=$(echo "$line" | cut -d, -f2)
TCRB=$(echo "$line" | cut -d, -f3)

echo "Processing sample ID: $sample_id"
echo "TCRA Sequence: $TCRA"
echo "TCRB Sequence: $TCRB"


singularity run --nv -B /public/alphafold_db_20231114:/database {self.tcrmodel_path} --pwd /opt/tcrmodel2  python /opt/tcrmodel2/run_tcrmodel2_ub_tcr.py \
--job_id=$sample_id \
--output_dir={self.output_dir}/TCR \
--tcra_seq=$TCRA \
--tcrb_seq=$TCRB \
--ori_db=/database/ \
--tp_db={self.databases_path} \
{self.keep_pkl} \
--relax_structures={self.relax_structures} \
--max_template_date={self.max_template_date}'''

        return self.template