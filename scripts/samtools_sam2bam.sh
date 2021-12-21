#!/usr/bin/env bash
#SBATCH -J sam2bam #Job Name
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

input_sam=$1
echo "input_sam: $input_sam"

output_name=${input_sam/.sam/}
echo "output_name: $output_name"

container="/container_path/rnaseq.sif"
echo "container: $container"

singularity exec $container samtools view \
  -@$cpu -q 10 -b -o "$output_name".bam $input_sam
echo "samtools sam to bam process completed for $input_sam!"




## From the command line
## for f in sam_dir/*.sam; do echo "processing: "$f; sbatch ./samtools_sam2bam.sh $f; done
