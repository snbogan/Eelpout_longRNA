#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=Ld_R1_trimmed_featureCounts
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_trimmed_featureCounts.out
#SBATCH --error=Ld_R1_trimmed_featureCounts.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB

# Move to read count directory
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/read_counts/

# Load subread
module load miniconda3
conda activate subread

featureCounts -T 8 -a Ldear_hifionly_primary_AFPIII.gff -o Ld_R1_liver_cDNA_trimmed_unique_AFP_counts.txt -t gene -g gene_id /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/alignments/Ld_R1_liver_pass_cutadapt_minimap2_cDNA_sorted_unique.bam

featureCounts -T 8 -a Ldear_hifionly_primary_AFPIII.gff -o Ld_R1_liver_cDNA_trimmed_AFP_counts.txt -t gene -g gene_id /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/alignments/Ld_R1_liver_pass_cutadapt_minimap2_cDNA.sorted.bam
