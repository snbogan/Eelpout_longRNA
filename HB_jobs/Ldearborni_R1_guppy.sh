#!/bin/bash
#SBATCH --partition=96x24gpu4
#SBATCH --job-name=Ld_R1_guppy
#SBATCH --time=0-72:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_guppy.out
#SBATCH --error=Ld_R1_guppy.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=68GB

# Move directories
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/20240606_1038_MN32724_FAU44241_e56fb14d

# Load dorado
module load guppy

# hac basecall on pod5 directory
guppy_basecaller --input_path pod5s/ --save_path guppy_out --flowcell FLO-MIN112 --kit SQK-LSK112 --compress_fastq
