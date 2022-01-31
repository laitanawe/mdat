#!/usr/bin/env bash
#SBATCH -J feature_counts
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 24
#SBATCH --mem 16g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

module load singularity

cpus=24
mapq_min=10
min_overlap=5
gtf='/db/GRCh38/release-103/gtf/Homo_sapiens.GRCh38.103.gtf'
out_file='feature_counts.txt'
container='/containers/rnaseq.sif'

echo "begin:$(date +'%Y%m%d%H%M%S')"
echo "cpus:$cpus"
echo "mapq_min:$mapq_min"
echo "min_overlap:$min_overlap"
echo "gtf:$gtf"
echo "out_file:$out_file"
echo "container:$container"

## input alignment files:
files=($(ls *.bam))
echo "nfiles:${#files[@]}"
echo "files:${files[*]}"

## the work:
result=$(singularity exec /containers/rnaseq.sif \
  featureCounts \
    -T $cpus \
    -Q $mapq_min \
    --minOverlap $min_overlap \
    --primary \
    --ignoreDup \
    -a "$gtf" \
    -o "$out_file" \
    ${files[*]} 2>&1)

[ $? -eq 0 ] || {
  echo "ERROR:counts: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "counts_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
