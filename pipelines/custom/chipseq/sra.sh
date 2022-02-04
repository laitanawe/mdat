#!/usr/bin/env bash
#SBATCH -J sra
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH -c 2
#SBATCH --mem 2g
#SBATCH -o '%x.%j.out'
#SBATCH -e '%x.%j.err'

set -uo pipefail
module load singularity

## takes NCBI SRA id (starts w/ 'SRA') and outputs corresponding .fastq.gz file

container='/containers/chipseq.sif'

[ $# -eq 1 ] || {
  echo "USAGE: sbatch sra.sh <sra_id>" >&2
  echo "  creates .fastq.gz file corresponding to <sra_id>" >&2
  exit 11
}
sra_id="$1"
echo "sra_id:$sra_id"
echo "container:$container"
echo "pwd:$(pwd)"

result=$(singularity exec chipseq.sif prefetch --version | grep prefetch 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:prefetch_version: $result" >&2
  exit 21
}
echo "prefetch_version:$result"

result=$(singularity exec chipseq.sif fastq-dump --version | grep fastq 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastq-dump_version: $result" >&2
  exit 23
}
echo "fastq-dump_version:$result"

echo "begin_prefetch:$(date +'%Y%m%d%H%M%S')"
result=$(prefetch -v "$sra_id" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:prefetch: $result" >&2
  exit 31
}
[ -n "$result" ] && echo "prefetch_result:$result"

echo "begin_fastq:$(date +'%Y%m%d%H%M%S')"
cd "$sra_id"
result=$(fastq-dump "${sra_id}.sra" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:fastq: $result" >&2
  exit 33
}
[ -n "$result" ] && echo "fastq_result:$result"

echo "removing .sra"
rm "${sra_id}.sra"

echo "begin_gzip:$(date +'%Y%m%d%H%M%S')"
result=$(gzip "${sra_id}.fastq" 2>&1)
[ $? -eq 0 ] || {
  echo "ERROR:gzip: $result" >&2
  exit 35
}
[ -n "$result" ] && echo "gzip_result:$result"

echo "finish:$(date +'%Y%m%d%H%M%S')"
