#!/usr/bin/env bash
#SBATCH -J sra
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem 8g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

## takes NCBI SRA id (starts w/ 'SRA') and outputs corresponding .fastq.gz file

export PATH="$PATH:/opt/sratoolkit/bin"

[ $# -eq 1 ] || {
  echo "USAGE: sbatch sra.sh <sra_id>" >&2
  echo "  creates .fastq.gz file corresponding to <sra_id>" >&2
  exit 11
}
sra_id="$1"
echo "sra_id:$sra_id"

echo "begin_prefetch:$(date +'%Y%m%d%H%M%S')"
result=$(prefetch -v "$sra_id" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:prefetch: $result" >&2
  exit 21
}
[ -n "$result" ] && echo "prefetch_result:$result"

echo "begin_fastq:$(date +'%Y%m%d%H%M%S')"
cd "$sra_id"
result=$(fastq-dump "${sra_id}.sra" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastq: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "fastq_result:$result"

rm "${sra_id}.sra"
gzip "${sra_id}.fastq"

echo "finish:$(date +'%Y%m%d%H%M%S')"
