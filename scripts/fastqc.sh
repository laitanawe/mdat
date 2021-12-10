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

file_in=$1
echo "file_in: $file_in"

container="/container_path/fastqc.sif"
echo "container: $container"

singularity exec $container fastqc -t $cpu $file_in
echo "fastqc process completed for $file_in!"


## From the command line
## for f in fastq_dir/*.fastq; do echo "processing: "$f; ./fastqc.sh $f; done
