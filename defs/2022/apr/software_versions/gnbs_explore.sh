#!/usr/bin/env bash

#Data Download Tools: These packages come with the VM
  wget --version | head -1
  zip --version | head -2 | tail -1 | sed 's/^This is //'
  perl --version | head -2 | tail -1 | sed 's/^This is p/P/'
  python2 --version | head -2 | tail -1 | sed 's/^ *//'
  python3 --version | head -2 | tail -1 | sed 's/^ *//'

  #Bionformatics Tools:
  alias trimmomatic='java -jar /usr/share/java/trimmomatic-0.36.jar'
  echo -n "trimmomatic: "; java -jar /usr/share/java/trimmomatic-0.36.jar -version | head -n 1
  fastqc --version
  echo -n "STAR: "; STAR --version
  R --version | head -n 1
  echo -n "Bioconductor:"; Rscript -e 'BiocManager::version()' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the Bioconductor vers$
  echo -n "ggplot2:"; Rscript -e 'packageVersion("ggplot2")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the ggplot2 version
#  echo -n "GAPIT:"; Rscript -e 'packageVersion("GAPIT3")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
  echo -n "BLINK:"; Rscript -e 'packageVersion("BLINK")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
  echo -n "Cutadapt: "; cutadapt --version
  samtools --version | head -n 1
  hisat2 --version | head -1
  subread-align -version 2>&1 | sed '/^$/d'
  featureCounts -v 2>&1 | sed '/^$/d'
  echo -n "edgeR:"; Rscript -e 'packageVersion("edgeR")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the edgeR version
  echo -n "limma:"; Rscript -e 'packageVersion("limma")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the limma version
  echo -n "Nextflow: "; sudo nextflow -version | head -3 | tail -1 | sed 's/^ *//'
#  htseq-count -h | tail -1 | sed 's/^Public License v3. Part of the //'
# GATK
  gatk --version 2>&1 | head -3 # | sed 's/[^Genome]*//'
  echo -n "snakemake: "; snakemake --version
  echo -n "BWA: "; bwa 2>&1 | head -3 | tail -1
  echo -n "trim_galore: "; trim_galore --version | tail -4 | head -1 | sed 's/^ *//'
