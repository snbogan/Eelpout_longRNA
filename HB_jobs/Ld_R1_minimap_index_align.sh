#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=Ld_R1_minimap_index_align
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=Ld_R1_minimap_index_align.out
#SBATCH --error=Ld_R1_minimap_index_align.err
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=68GB

# Move directories
cd /hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/alignments

# Load minimap2 and samtools
module load minimap2
module load samtools

# Create inputs, outputs
REFERENCE="/hb/home/snbogan/PolarFish/Long_AFP/new_genomes/Ldear_hifionly.asm.bp.p_ctg.fa"

READS="/hb/home/snbogan/PolarFish/RNAseq/Ldear_run1/20240606_1038_MN32724_FAU44241_e56fb14d/guppy_out_parallel/Ld_liver_rna_ont_r1_gpu1_pass.fastq.gz"

OUTPUT="Ld_R1_liver_pass_minimap2.sam"

# Index L dearborni primary reference genome
echo "Indexing the reference genome..."
minimap2 -d Ld_primary_reference.mmi "$REFERENCE"

# Align with minimap2
echo "Aligning reads to the reference genome..."
minimap2 -ax splice -uf -k14 -t 8 Ld_primary_reference.mmi "$READS" > "$OUTPUT"

# Bam convert and sort sam file
echo "Converting SAM to BAM, sorting, and indexing..."
samtools view -@ 8 -bS "$OUTPUT" | samtools sort -@ 8 -o Ld_R1_liver_pass_minimap2.sorted.bam
samtools index Ld_R1_liver_pass_minimap2.sorted.bam









