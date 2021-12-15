#!/usr/bin/env bash
#SBATCH -J fastqc #Job Name
#SBATCH	-N 1 #Number of nodes
#SBATCH	-q batch #QOS
#SBATCH	-t 72:00:00 #Time
#SBATCH	-c 24 #CPUs
#SBATCH	--mem 64g #memory
#SBATCH	-o './%x.%j.out' #output file; %x is Job Name, %j is Job number
#SBATCH -e './%x.%j.err' #error file; %x is Job Name, %j is Job number

module load singularity

cpu=24
echo "cpu:$cpu"

ram=64
echo "ram:$ram"

genome_index_prefix='/path/to/databases/GRCm39/bowtie2/GRCm39'
echo "genome_index: $genome_index_prefix"

file_in=$1
echo "file_in: $file_in"

output_sam="$($file_in/\.fastq/.sam)"
echo "output_sam: $output_sam"

container="/container_path/fastqc.sif"
echo "container: $container"

singularity exec $container bowtie2 -threads $cpu --local -x $genome_index_prefix \
-U $file_in -S $output_sam
echo "bowtie2 alignment process completed for $file_in!"


## From the command line
## for f in fastq_dir/*.fastq; do echo "processing: "$f; sbatch ./bowtie2_align.sh $f; done
