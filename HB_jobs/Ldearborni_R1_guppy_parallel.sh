#!/bin/bash
#SBATCH --partition=96x24gpu4
#SBATCH --job-name=Ld_R1_guppy_parallel
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_guppy_%A_%a.out
#SBATCH --error=Ld_R1_guppy_%A_%a.err
#SBATCH --ntasks=1                   # Only one task for the parallel execution
#SBATCH --cpus-per-task=4            # Adjust CPU cores based on available resources
#SBATCH --gres=gpu:4                 # Request 4 GPUs
#SBATCH --mem=68GB                   # Shared memory for all tasks

# Move to the working directory
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/20240606_1038_MN32724_FAU44241_e56fb14d

# Load guppy and parallel
module load guppy
module load parallel

# Get the list of input files (glob pattern)
INPUT_DIR="pod5s"
FILES=($(ls ${INPUT_DIR}/*.pod5))

# Number of files
NUM_FILES=${#FILES[@]}

# Calculate the number of files per task
FILES_PER_TASK=4

# Distribute files across tasks
for i in $(seq 0 $((NUM_FILES / FILES_PER_TASK))); do
  START_INDEX=$((i * FILES_PER_TASK))
  END_INDEX=$(((i + 1) * FILES_PER_TASK - 1))

  # Ensure we don't go out of bounds
  if [ $END_INDEX -ge $NUM_FILES ]; then
    END_INDEX=$((NUM_FILES - 1))
  fi

  # Slice the files to process
  FILES_TO_PROCESS="${FILES[@]:$START_INDEX:$((END_INDEX - START_INDEX + 1))}"

  # Run `guppy_basecaller` in parallel across 4 tasks
  parallel -j 4 --delay 0.1 --no-notice "guppy_basecaller \
    --input_path ${INPUT_DIR} \
    --save_path guppy_out_parallel/gpu_{#} \
    --flowcell FLO-MIN112 \
    --kit SQK-LSK112 \
    --compress_fastq \
    --device cuda:{#} " ::: ${FILES_TO_PROCESS[@]}
done
