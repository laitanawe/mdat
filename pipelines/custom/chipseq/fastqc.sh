#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH -q batch
#SBATCH -c 2
#SBATCH -t 72:00:00
#SBATCH --mem 8g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

set -uo pipefail
module load singularity

container=/containers/chipseq.sif
cpus=2

[ $# -eq 1 ] || {
  echo "USAGE: sbatch fastqc.sh <reads.fastq>" >&2
  exit 11
}
reads_fastq="$1"

echo "reads_fastq:$reads_fastq"
echo "container:$container"
echo "cpus:$cpus"
echo "pwd:$(pwd)"

## print out software version:

result=$(singularity exec "$container" fastqc --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastqc_version: $result" >&2
  exit 21
}
echo "fastqc_version:$result"

## actual work:

echo "begin:$(date +'%Y%m%d%H%M%S')"

result=$(singularity exec "$container" fastqc -t $cpus "$reads_fastq" 2>&1)

[ $? -eq 0 ] || {
  echo "ERROR:fastqc: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "fastqc_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
