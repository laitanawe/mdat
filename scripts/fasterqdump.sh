#!/usr/bin/env bash
#SBATCH -J fasterqdump #Job Name
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

accession_list=$1
echo "accession_list: $accession_list"

outdir=$2
echo "outdir: $outdir"

container="/projects/researchit/crf/containers/data/defs/sratoolkit.sif"
echo "container: $container"

## E.g: If you have a list of accessions, create fastq dir, you can test if SRA Toolkit works multithread 64cores:
## cat accession_list.txt | xargs singularity exec mycontainer.sif fasterq-dump	--outdir . -e 64

cat $accession_list | xargs singularity exec $container fasterq-dump --outdir $outdir -e $ram
