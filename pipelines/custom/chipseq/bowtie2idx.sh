#!/usr/bin/env bash
#SBATCH -J bowtie2idx
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 24
#SBATCH --mem 16g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

set -uo pipefail
module load singularity

cpus=24
container=/containers/chipseq.sif
seed=101

[ $# -eq 2 ] || {
  echo "USAGE: sbatch bowtie2idx.sh <genome.fasta> <genome_index>" >&2
  echo "  Creates a bowtie2 genome index from a genome fasta file." >&2
  exit 11
}
genome_fasta="$1"
genome_index="$2"

echo "genome_fasta:$genome_fasta"
echo "genome_index:$genome_index"
echo "container:$container"
echo "cpus:$cpus"
echo "seed:$seed"
echo "pwd:$(pwd)"

## print out software version:

result=$(singularity exec "$container" bowtie2-build --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:bowtie2-build_version: $result" >&2
  exit 21
}
echo "bowtie2-build_version:$result"

## actual work:

echo "begin:$(date +'%Y%m%d%H%M%S')"

result=$(singularity exec "$container" bowtie2-build \
  --threads $cpus \
  --seed $seed \
  "$genome_fasta" "$genome_index" 2>&1)

[ $? -eq 0 ] || {
  echo "ERROR:bowtie2-build: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "bowtie2-build_result:$result"

echo "finish: $(date +'%Y%m%d%H%M%S')"
