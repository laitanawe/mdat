#!/usr/bin/env bash
#SBATCH -J bowtie2
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 16
#SBATCH --mem 16g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

set -uo pipefail
module load singularity

cpus=16
container=/containers/chipseq.sif

[ $# -eq 3 ] || {
  echo "USAGE: sbatch bowtie2.sh <genome_index> <reads.fastq[.gz]> <output.sam>" >&2
  exit 11
}
genome_index="$1"
reads_fastq="$2"
output_sam="$3"

echo "genome_index:$genome_index"
echo "reads_fastq:$reads_fastq"
echo "output_sam:$output_sam"
echo "container:$container"
echo "cpus:$cpus"
echo "pwd:$(pwd)"

## print out software version:

result=$(singularity exec "$container" bowtie2 --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:bowtie2_version: $result" >&2
  exit 21
}
echo "bowtie2_version:$result"

## actual work:

echo "begin:$(date +'%Y%m%d%H%M%S')"

result=$(singularity exec "$container" bowtie2 \
  -p $cpus \
  --local \
  -x "$genome_index" \
  -U "$reads_fastq" \
  -S "$output_sam" 2>&1)

[ $? -eq 0 ] || {
  echo "ERROR:bowtie2: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "bowtie2_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
