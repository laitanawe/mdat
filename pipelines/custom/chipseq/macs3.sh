#!/usr/bin/env bash
#SBATCH -J macs3
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 2
#SBATCH --mem 8g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

## <chip.list> and <sham.list> are text files containing list of full paths to bam files (one per line):
[ $# -eq 3 ] || {
  echo "Usage: sbatch macs3.sh <chip_list> <sham_list> <prefix_out>" >&2
  exit 11
}
chip_list="$1"       ## real immunoprecipitation data
sham_list="$2"       ## sham immunoprecip data
prefix_out="$3"      ## prefix for output files

genome_size='hs'     ## presets: hs:human[2.7e9]; mm:mouse[1.9e9]; ce:worm[9e7]; dm:fly[1.2e8]. Enter haploid size otherwise.
fdr_cutoff=0.05      ## false discovery rate cutoff
save_bedgraph=0      ## 0: save; 1: don't save
input_format='BAM'   ## type of files listed in chip_list and sham_list (bam works better than sam)
container='/containers/chipseq.sif'

module load singularity

echo "chip_list:$chip_list"
echo "sham_list:$sham_list"
echo "prefix_out:$prefix_out"
echo "genome_size:$genome_size"
echo "fdr_cutoff:$fdr_cutoff"
echo "save_bedgraph:$save_bedgraph"
echo "input_format:$input_format"
echo "container:$container"
echo "pwd:$(pwd)"

rslt=$(singularity exec "$container" macs3 --version 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:macs3_version:$rslt" >&2
  exit 13
}
echo "macs3_version:$rslt"

## make array of chip file paths:
chip_files=()
while read f; do
  f=$(echo "$f" | sed 's/^ *//' | sed 's/ *$//')
  [ -z "$f" ] && continue
  chip_files+=("$f")
done < "$chip_list"
echo "chip_files:${chip_files[@]}"

## make array of sham file paths:
sham_files=()
while read f; do
  f=$(echo "$f" | sed 's/^ *//' | sed 's/ *$//')
  [ -z "$f" ] && continue
  sham_files+=("$f")
done < "$sham_list"
echo "sham_files:${sham_files[@]}"

## set -B flag (or not):
flag_B=''
[ "$save_bedgraph" == 0 ] && flag_B='-B'
echo "flag_B:$flag_B"

## get to work:

echo "begin:$(date +'%Y%m%d%H%M%S')"

result=$(singularity exec "$container" macs3 callpeak $flag_B -t ${chip_files[@]} -c ${sham_files[@]}  \
       -n "$prefix_out" -f $input_format -g $genome_size -q $fdr_cutoff 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR: macs3: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "macs3_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
