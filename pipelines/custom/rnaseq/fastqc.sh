#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH -q batch
#SBATCH -c 2
#SBATCH -t 72:00:00
#SBATCH --mem 8g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

module load singularity
container=/containers/chipseq.sif
cpus=2

[ $# -eq 1 ] || {
  echo "USAGE: sbatch fastqc.sh <reads.fastq>" >&2
  exit 11
}
reads="$1"

echo "reads:$reads"
echo "container:$container"
echo "cpus:$cpus"

## print out software version:

rslt=$(singularity exec "$container" fastqc --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastqc_version: $rslt" >&2
  exit 21
}
echo "fastqc_version:$rslt"

## actual work:

echo "begin: $(date +'%Y%m%d%H%M%S')"

rslt=$(singularity exec "$container" fastqc -t $cpus "$reads" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastqc: $rslt" >&2
  exit 31
}
[ -n "$rslt" ] && echo "fastqc_result:$rslt"

echo "finish: $(date +'%Y%m%d%H%M%S')"
