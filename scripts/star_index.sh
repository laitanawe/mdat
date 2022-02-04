#!/usr/bin/env bash
#SBATCH -J star_index #Job Name
#SBATCH	-N 1 #Number of nodes
#SBATCH	-q batch #QOS
#SBATCH	-t 72:00:00 #Time
#SBATCH	-c 24 #CPUs
#SBATCH	--mem 64g #memory
#SBATCH	-o './%x.%j.out' #output file; %x is Job Name, %j is Job number
#SBATCH -e './%x.%j.err' #error file; %x is Job Name, %j is Job number

module load singularity

cpus=64
echo "cpu: $cpu"

ram=64
echo "ram: $ram"

dir_out=stardb
echo "dir_out: $dir_out"

genome_fasta=$1
echo "genome_fasta: $genome_fasta"

container="/container_path/rnaseq.sif"
echo "container: $container"

singularity exec $container STAR \
   --runThreadN $cpus \
   --runMode genomeGenerate \
   --genomeDir $dir_out \
   --genomeFastaFiles $genome_fasta


echo "STAR indexing process completed for the reference!"
