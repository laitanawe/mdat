#!/usr/bin/env bash
#Read count quantification
#Script to count the number of reads aligned to the reference genome using HTSeq.
#module load htseq
mycounts="counts"
echo "Creating directory for Raw Gene Expression Quantification ongoing ... $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"
mkdir -pv $mycounts
#cd $mycounts
#cd /srv/data/my_shared_data_folder/ace2covid/results/feature-counts
# process paired-end data
#BAM_DIR="/srv/data/my_shared_data_folder/ace2covid/results/bam-files"
BAM_DIR="../../../stubdata/sortedbams_coord_bwa2"
bam_file_name=$(basename "$BAM_DIR" .sorted.bam)
#GTF_FILE="/srv/data/my_shared_data_folder/ace2covid/data/ref-index/GCF_000001405.39_GRCh38.p13_genomic.gtf"
GTF_FILE="../../../databases/human_gencode/GCA_000001405.29_GRCh38.p14_genomic.gtf.gz"
n=0
for bam_file in ${BAM_DIR}/*.coordsort.bam; do b=$(basename $bam_file); acc=$(echo $b | sed 's/.coordsort.bam//');
      ((n++))
      echo; \
      echo "processing BAM File $n: "$b;\
      htseq-count \
            -f bam \
            -r pos \
            -s no \
            -t exon \
            -i gene \
            $bam_file \
            $GTF_FILE \
            > $mycounts/${acc}.counts.txt
      done
      echo
echo "Raw Gene Expression Quantification process completed for all .BAM files now, $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"
