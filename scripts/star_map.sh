#!/usr/bin/env bash
#SBATCH -J star_map #Job Name
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

prefix_out=$1
echo "prefix_out: $prefix_out"

reads_f=$2
echo "reads_forward: $reads_f"

reads_r=$3
echo "reads_reverse: $reads_r"

genome_fasta=/path/to/data/reference/gencode_hs39/GRCh38.primary_assembly.genome.fa
genome_gtf=/path/to/data/reference/gencode_hs39/gencode.v39.primary_assembly.annotation.gtf
genome_idx=/path/to/data/rnaseq_samples/stardb

container="/container_path/rnaseq.sif"
echo "container: $container"

module load singularity
singularity exec $container STAR \
   --runThreadN $cpus --genomeDir $genome_idx --sjdbGTFfile $genome_gtf \
   --quantMode GeneCounts --outSAMtype BAM Unsorted \
   --outFileNamePrefix $prefix_out --readFilesIn $reads_f $reads_r


echo "STAR mapping process completed for the reference!"


## Usage:
## sn=0; for f_forward in ../fastqs/SRR*_1.fastq; do f_reverse="${f_forward/\_1/_2}"; ((sn++)); sbatch /scripts/star_map.sh mapped_${sn}'_' $f_forward $f_reverse; done
## for f_forward in SRR*_1.fastq; do f_reverse="${f_forward/\_1/_2}"; sbatch /scripts/star_map.sh ${f_forward/\_1.fastq/}'_' $f_forward $f_reverse; done
## No Job Scheduling:
## for f_forward in SRR*_1.fastq; do f_reverse="${f_forward/\_1/_2}"; /scripts/star_map.sh ${f_forward/\_1.fastq/}'_' $f_forward $f_reverse; done
