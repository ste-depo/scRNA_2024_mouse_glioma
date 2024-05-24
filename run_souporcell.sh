#!/bin/bash 
#SBATCH --job-name=souporcell_pipeline
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

source ~/.bashrc

# create reference vcf using M2_veh sample

singularity exec --bind ${GENOMES_FOLDER}/refdata-gex-mm10-2020-A/fasta:/fasta,InputData/M2_veh/outs:/mnt \
	  ${PATH_TO_SOUPORCELL_SIF}/souporcell_latest.sif souporcell_pipeline.py \
	  -i /mnt/possorted_genome_bam.bam -b /mnt/filtered_feature_bc_matrix/barcodes.tsv.gz \
	  -f /fasta/genome.fa -t 16 -o /mnt/soup_or_cell_ref -k 2

# run souporcell on individual samples using the reference vcf created using M2_veh

for FOLDER in InputData/*/outs;
do
	singularity exec --bind ${GENOMES_FOLDER}/refdata-gex-mm10-2020-A/fasta:/fasta,$FOLDER:/mnt,InputData/M2_veh/outs/soup_or_cell_ref:/ref \
	  ${PATH_TO_SOUPORCELL_SIF}/souporcell_latest.sif souporcell_pipeline.py \
	  --common_variants /ref/souporcell_merged_sorted_vcf.vcf.gz \
	  -i /mnt/possorted_genome_bam.bam -b /mnt/filtered_feature_bc_matrix/barcodes.tsv.gz \
	  -f /fasta/genome.fa -t 16 -o /mnt/soup_or_cell -k 2
done
