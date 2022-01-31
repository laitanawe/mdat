#!/usr/bin/env bash
#SBATCH -J star
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 12
#SBATCH --mem 64g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

ncpus=12                   ## number of threads (should match SBATCH -c setting!!!)
genome_idx='/db/GRCh38/release-103/dna/star'   ## star index directory
genome_gtf='/db/GRCh38/release-103/gtf/Homo_sapiens.GRCh38.103.gtf'
container='/containers/rnaseq.sif'
multimap_nmax=1            ## max number of genomic locations read maps to (else read unmapped)
mismatch_pmax=0.04         ## max proportion of read residues not matching template (else read unmapped)

[ $# -eq 2 ] || [ $# -eq 3 ] || {
  echo "Usage: sbatch star.sh <prefix_out> <reads_forward.fastq> [<reads_reverse.fastq>]" >&2
  exit 11
}
prefix_out="$1"            ## output prefix (can include path)
reads_f="$2"               ## forward strand reads
reads_r="$3"               ## optional reverse strand read for each forward read

echo "prefix_out:$prefix_out"
echo "reads_f:$reads_f"
echo "reads_r:$reads_r"
echo "genome_idx:$genome_idx"
echo "genome_gtf:$genome_gtf"
echo "multimap_nmax:$multimap_nmax"
echo "mismatch_pmax:$mismatch_pmax"
echo "ncpus:$ncpus"
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
  --runThreadN $ncpus --genomeDir "$genome_idx" --sjdbGTFfile "$genome_gtf" \
  --quantMode GeneCounts --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$prefix_out" --readFilesIn "$reads_f" "$reads_r" >&2)

[ $? -eq 0 ] || {
  echo "ERROR:star: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "star_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
