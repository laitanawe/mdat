#!/usr/bin/env bash
#SBATCH -J bwt2_index #Job Name
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

output_index_dir_name='/path/to/databases/GRCm39/bowtie2/GRCm39'
echo "genome_index: $genome_index_prefix"
# /path/to/databases/GRCm39/ contains the reference genome: GRCm39_genomic.fna, .gff, .gtf
# /path/to/databases/GRCm39/bowtie2/GRCm39*.* are the index files

file_in=$1
echo "file_in: $file_in"

container="/container_path/rnaseq.sif"
echo "container: $container"

singularity exec $container bowtie2-build --threads $cpu $file_in $output_index_dir_name
echo "bowtie2 index build process completed for $file_in!"


## From the command line
## for f in fastq_dir/*.fastq; do echo "processing: "$f; sbatch ./bowtie2_align.sh $f; done
