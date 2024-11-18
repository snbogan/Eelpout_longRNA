#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=Ld_R1_fastqc
#SBATCH --time=0-6:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_fastqc.out
#SBATCH --error=Ld_R1_fastqc.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB

# Move directories
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/20240606_1038_MN32724_FAU44241_e56fb14d/guppy_out_parallel

# Load fastqc
module load fastqc

# Run on passed reads
fastqc -o fastqc/pass Ld_liver_rna_ont_r1_gpu1_pass.fastq.gz

# Run on passed reads
fastqc -o fastqc/fail Ld_liver_rna_ont_r1_gpu1_fail.fastq.gz









