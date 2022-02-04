#!/usr/bin/env bash
#SBATCH -J featurecounts #Job Name
#SBATCH	-N 1 #Number of nodes
#SBATCH	-q batch #QOS
#SBATCH	-t 72:00:00 #Time
#SBATCH	-c 24 #CPUs
#SBATCH	--mem 64g #memory
#SBATCH	-o './%x.%j.out' #output file; %x is Job Name, %j is Job number
#SBATCH -e './%x.%j.err' #error file; %x is Job Name, %j is Job number

cpus=64
echo "cpu:$cpus"

ram=64
echo "ram:$ram"

out_file=$1
echo "out_file: $out_file"

## The next line file_array works for the mapped reads
file_array=$(echo $(ls ../bams/*.out.bam) ) ## read lines from standard input into an array variable and this approach removes trailing new lines, for paired end reads, DON'T pre-sort bams b4 gene counting
#file_array=$(ls ../bams/*.bam) ## read lines from standard input into an array variable and this approach contains trailing new lines
echo "file_array: $file_array"

# These lines don't work:
#files[*]=$2
#files=$2
#files[*]=$(ls ../bams/*.bam)
#echo "files: $files"
#echo "files: $files"

mapq_min=10
echo "mapq_min: $mapq_min"

min_overlap=5
echo "min_overlap: $min_overlap"

genome_gtf=/path/to/data/reference/gencode_hs39/gencode.v39.primary_assembly.annotation.gtf

container="/container_path/rnaseq.sif"
echo "container: $container"

module load singularity
singularity exec $container featureCounts \
   -T $cpus -Q $mapq_min --minOverlap $min_overlap --primary --ignoreDup \
   -a $genome_gtf -o $out_file \
   $file_array
#   ${files[*]}

echo "featureCounts process completed!"

## Usage:
## ls ../bams/*.bam | xargs /scripts/featurecounts.sh counts.txt 
