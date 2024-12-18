#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=Ld_R1_pod5_colibri
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_pod5_colibri.out
#SBATCH --error=Ld_R1_pod5_colibri.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=68GB

# Move directories
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/20240606_1038_MN32724_FAU44241_e56fb14d

# Load pod5
module load miniconda3

conda activate pod5

# Find all Fast5 files and split into 5 batches
files=(./fast5/*.fast5)
total_files=${#files[@]}
batch_size=$((total_files / 5))

echo "Total Fast5 files: $total_files"
echo "Batch size: $batch_size"

# Loop through 5 batches
for i in {0..4}; do
    start=$((i * batch_size))
    end=$((start + batch_size - 1))

    # For the last batch, make sure it includes any remaining files
    if [ $i -eq 4 ]; then
        end=$((total_files - 1))
    fi

    batch_files="${files[@]:$start:$((end - start + 1))}"

    # Run pod5 conversion for the current batch
    echo "Processing batch $((i + 1))..."
    pod5 convert fast5 $batch_files --output Ldearborni_R1_converted_colibri_$((i + 1)).pod5

    # Check if the conversion was successful for this batch
    if [ -f "Ldearborni_R1_converted_colibri_$((i + 1)).pod5" ]; then
        echo "Batch $((i + 1)) conversion successful"
    else
        echo "Batch $((i + 1)) conversion failed" >&2
        exit 1
    fi
done

echo "All batches processed successfully."





