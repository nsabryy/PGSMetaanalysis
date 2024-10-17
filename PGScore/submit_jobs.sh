#!/bin/bash

# Set the total number of PGS IDs and the batch size
START=1261
TOTAL_IDS=1270
BATCH_SIZE=1
OUTPUT_DIR="/N/project/compgen/PGSCalc/scoring_results/score_analysis"

# Loop over the ranges and submit a job for each range
for ((i=START; i<=TOTAL_IDS; i+=BATCH_SIZE)); do
  start_index=$i
  end_index=$((i + BATCH_SIZE - 1))

  # Ensure the last batch doesn't exceed TOTAL_IDS
  if [ $end_index -gt $TOTAL_IDS ]; then
    end_index=$TOTAL_IDS
  fi

  # Set the output file for this range
  output_file="${OUTPUT_DIR}/score_metadata_${start_index}_to_${end_index}.tsv"

  # Submit the SLURM job with the calculated indices and output path
  sbatch run_metadata.sh $start_index $end_index $output_file
done