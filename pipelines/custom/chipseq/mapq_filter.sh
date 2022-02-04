#!/usr/bin/env bash
#SBATCH -J mapq_filter
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 8
#SBATCH --mem 8g          ## >=1g/cpu
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

set -uo pipefail
module load singularity

cpus=8
container=/containers/chipseq.sif
mapq_cutoff=10

[ $# -eq 2 ] || {
  echo "Usage: sbatch mapq_filter.sh <input.sam> <output.bam>" >&2
  echo "  filters <input.sam> for alignments w/ mapq >= $mapq_cutoff and converts output to .bam"
  exit 11
}
input_sam="$1"
output_bam="$2"

echo "input_sam:$input_sam"
echo "output_bam:$output_bam"
echo "mapq_cutoff:$mapq_cutoff"
echo "container:$container"
echo "cpus:$cpus"
echo "pwd:$(pwd)"

result=$(samtools --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:samtools_version: $result" >&2
  exit 21
}
echo "samtools_version:$result"

uid="$(hostname).$(date +'%s').$$"
echo "uid:$uid"

echo "begin_filter:$(date +'%Y%m%d%H%M%S')"
result=$(singularity exec "$container" samtools view \
  -@$cpus -q $mapq_cutoff -b -o "$output_bam.$uid" "$input_sam" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:filter: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "filter_result:$result"

echo "begin_sort:$(date +'%Y%m%d%H%M%S')"
result=$(singularity exec "$container" samtools sort \
  -@$cpus -o "$output_bam" "$output_bam.$uid" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:sort: $result" >&2
  exit 33
}
[ -n "$result" ] && echo "sort_result:$result"

echo "removing:'$output_bam.$uid'"
rm "$output_bam.$uid"

echo "begin_index:$(date +'%Y%m%d%H%M%S')"
result=$(singularity exec "$container" samtools index \
  -@$cpus "$output_bam" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:index: $result" >&2
  exit 35
}
[ -n "$result" ] && echo "index_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
