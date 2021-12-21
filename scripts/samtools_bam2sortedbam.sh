#!/usr/bin/env bash
#SBATCH -J bam2sortedbam #Job Name
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

input_bam=$1
echo "input_bam: $input_bam"

output_name=${input_bam/.bam/}
echo "output_name: $output_name"

container="/container_path/rnaseq.sif"
echo "container: $container"

singularity exec $container samtools sort -@$cpu -O BAM -o "$output_name".sort.bam $input_bam
echo "samtools bam to sorted bam process completed for $input_bam!"

## From the command line
## for f in bam_dir/*.bam; do echo "processing: "$f; sbatch ./samtools_bam2sortedbam.sh $f; done
