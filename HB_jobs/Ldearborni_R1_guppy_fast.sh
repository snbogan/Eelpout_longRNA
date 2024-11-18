#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=Ld_R1_guppy_fast
#SBATCH --time=0-72:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_guppy_fast.out
#SBATCH --error=Ld_R1_guppy_fast.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=68GB

# Move directories
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/20240606_1038_MN32724_FAU44241_e56fb14d

# Load dorado
module load guppy

# hac basecall on pod5 directory
guppy_basecaller --input_path pod5s/ --save_path guppy_fast_out --flowcell FLO-MIN112 --kit SQK-LSK112 --compress_fastq --config dna_r9.4.1_450bps_hac.cfg --device cpu
