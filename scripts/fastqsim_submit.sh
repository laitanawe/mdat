#!/usr/bin/env bash
#SBATCH -J fastqsim #Job Name
#SBATCH -N 1 #Number of nodes
#SBATCH -q batch #QOS
#SBATCH -t 72:00:00 #Time
#SBATCH -c 64 #CPUs
#SBATCH --mem 64g #memory
#SBATCH -o './%x.%j.out' #output file; %x is Job Name, %j is Job number
#SBATCH -e './%x.%j.err' #error file; %x is Job Name, %j is Job number

module load singularity

cpu=24
echo "cpu:$cpu"

ram=64
echo "ram:$ram"

file_in=$1
echo "file_in: $file_in"

file_ref=$2
echo "file_ref: $file_ref"

file_bam=$3
echo "file_bam: $file_bam"

container="/containers/fastqsim.sif"
echo "container: $container"

singularity exec $container ./FASTQcharacterize.sh -input $file_in -reference $file_ref -bam $file_bam
echo "FASTQSim process completed for $file_in!"

# Location:
# Command to Run from terminal
# for f in ../*.merged.fastq; do sbatch fastqsim_laura_submit.sh $f ../genome_ref/mm10.fa ${f/.fastq/.bam}; done
# You can change ../*.merged.fastq to match the path to your fastqs. This example assumes that your .bam and .fastq files are in the same path.
# This also assumes that your reference genome is in ../genome_ref/mm10.fa (You can always correct the path if different)
