#!/usr/bin/env bash
#SBATCH -J stardb
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 24
#SBATCH --mem 64g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

cpus=24
container='/containers/rnaseq.sif'

[ $# -eq 2 ] || {
  echo "Usage: sbatch stardb.sh <genome.fasta> <dir_out>" >&2
  exit 11
}
genome_fasta="$1"
dir_out="$2"

echo "genome_fasta:$genome_fasta"
echo "dir_out:$dir_out"
echo "cpus:$cpus"
echo "container:$container"

module load singularity

version=$(singularity exec "$container" STAR --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:star_version: $version" >&2
  exit 21
}
echo "star_version:$version"

echo "begin:$(date +'%Y%m%d%H%M%S')"

result=$(singularity exec "$container" STAR \
  --runThreadN $cpus --runMode genomeGenerate \
  --genomeDir "$dir_out" --genomeFastaFiles "$genome_fasta" 2>&1)

[ $? -eq 0 ] || {
  echo "ERROR:star: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "star_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
