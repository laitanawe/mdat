#!/usr/bin/env bash
#SBATCH -J indexbam #Job Name
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

input_sortedbam=$1
echo "input_bam: $input_sortedbam"

container="/container_path/rnaseq.sif"
echo "container: $container"

singularity exec $container samtools index -@$cpu $input_sortedbam
echo "samtools index bam process completed for $input_sortedbam!"

## From the command line
## for f in bam_dir/*.sort.bam; do echo "processing: "$f; sbatch ./samtools_indexbam.sh $f; done
