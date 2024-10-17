#!/bin/bash

#SBATCH -J combinevcf
#SBATCH -p general
#SBATCH -o %j.txt
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nsabry@iu.edu
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH -A r01003

module load python/3.10

# Set variables for job
JOB_ID=${SLURM_JOB_ID}
START_INDEX=$1
END_INDEX=$2
OUTPUT_PATH=$3

# Call the Python script with the additional parameters
python PGScore.py --job_id ${JOB_ID} --start_index ${START_INDEX} --end_index ${END_INDEX} --output_file ${OUTPUT_PATH}
