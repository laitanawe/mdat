#!/usr/bin/env bash
#SBATCH -J htseq_counts
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 24
#SBATCH --mem 32g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

module load singularity

ncpus=24
mapq_min=10
mode='union'
nonunique='none'
secondary='ignore'
supplemental='ignore'

gtf='/db/GRCh38/release-103/gtf/Homo_sapiens.GRCh38.103.gtf'
out_file='htseq_counts.txt'
container='/containers/rnaseq.sif'

echo "begin:$(date +'%Y%m%d%H%M%S')"
echo "ncpus:$ncpus"
echo "mapq_min:$mapq_min"
echo "mode:$mode"
echo "nonunique:$nonunique"
echo "secondary:$secondary"
echo "supplemental:$supplemental"
echo "gtf:$gtf"
echo "out_file:$out_file"
echo "container:$container"

## input alignment files:
files=($(ls *.bam))
echo "nfiles:${#files[@]}"
echo "files:${files[*]}"

## header row:
echo $(IFS=':' ; echo 'gene:'; echo "${files[*]}") | sed 's/Aligned.out.bam//g' | tr ':' $'\t' > "$out_file"

## the work:
result=$(singularity exec "$container" \
  htseq-count \
    -n $ncpus \
    -a $mapq_min \
    -m $mode \
    --nonunique=$nonunique \
    --secondary-alignments=$secondary \
    --supplementary-alignments=$supplemental \
    ${files[*]} "$gtf" >> "$out_file" 2>&1)

[ $? -eq 0 ] || {
  echo "ERROR:htseq: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "htseq_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
