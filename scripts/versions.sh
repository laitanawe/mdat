##!/usr/bin/env bash

#Data Download Tools: These packages come with the VM
wget --version | head -1
zip --version | head -2 | tail -1 | sed 's/^This is //'
perl --version | head -2 | tail -1 | sed 's/^This is p/P/'
python2 --version | head -2 | tail -1 | sed 's/^ *//'
python3 --version | head -2 | tail -1 | sed 's/^ *//'

#Bionformatics Tools:
R --version | head -n 1
echo -n "Bioconductor:"; Rscript -e 'BiocManager::version()' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
fasterq-dump --version | head -2 | tail -1
echo -n "STAR: "; STAR --version
fastqc --version
echo -n "trimmomatic: "; java -jar /usr/share/java/trimmomatic-0.36.jar -version | head -n 1
echo -n "trim_galore: "; trim_galore --version | tail -4 | head -1 | sed 's/^ *//'
python3.8 --version | head -2 | tail -1 | sed 's/^ *//'
echo -n "Numpy: "; python3 -c 'import numpy as np; print(np.__version__)'
echo -n "Matplotlib: "; python3 -c 'import matplotlib as mat; print(mat.__version__)'
echo -n "Scipy: "; python3 -c 'import scipy as sc; print(sc.__version__)' 2>&1 | tail -1
echo -n "Pandas: "; python3 -c 'import pandas as pd; print(pd.__version__)' 2>&1 | tail -1
echo -n "BioPython: "; python3 -c 'import Bio as bio; print(bio.__version__)'
echo -n "Seaborn: "; python3 -c 'import seaborn as sea; print(sea.__version__)'
echo -n "scikit-learn: "; python3 -c 'import sklearn as sk; print(sk.__version__)'

blastx -h | tail -3 | head -1 | sed 's/[^B]*//'
multiqc --version
echo -n "BLINK:"; Rscript -e 'packageVersion("BLINK")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "Cutadapt: "; cutadapt --version
echo -n "edgeR:"; Rscript -e 'packageVersion("edgeR")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "limma:"; Rscript -e 'packageVersion("limma")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "DESeq2:"; Rscript -e 'packageVersion("DESeq2")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "ggplot2:"; Rscript -e 'packageVersion("ggplot2")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "GAPIT:"; Rscript -e 'packageVersion("GAPIT3")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
samtools --version | head -n 1
hisat2 --version 2>&1 | head -1
subread-align -version 2>&1 | sed '/^$/d'
featureCounts -v 2>&1 | sed '/^$/d'
htseq-count -h | tail -1 | sed 's/^Public License v3. Part of the //'
echo -n "BWA "; bwa 2>&1 | head -3 | tail -1
echo -n "snakemake: "; snakemake --version
echo -n "Nextflow: "; nextflow -version | head -3 | tail -1 | sed 's/^ *//'
# GATK
gatk --version 2>&1 | tail -3 # | sed 's/[^Genome]*//'
docker --version
echo -n "Tidyverse:"; Rscript -e 'packageVersion("tidyverse")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
bowtie2 --version 2>&1 | head -1
tophat2 --version
cufflinks --version 2>&1 | head -2 | tail -1
echo -n "RforProteomics:"; Rscript -e 'packageVersion("RforProteomics")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
java -jar /opt/snpEff/snpEff.jar -version

hmmalign -h | head -2 | tail -1 | sed 's/^# //'
echo -n "DECIPHER:"; Rscript -e 'packageVersion("DECIPHER")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
spades --version
echo -n "castor:"; Rscript -e 'packageVersion("castor")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "soapdenovo2: "; soapdenovo2-127mer --version 2>&1 | tail -10 | head -1
echo -n "rRDP:"; Rscript -e 'packageVersion("rRDP")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "phyloseq:"; Rscript -e 'packageVersion("phyloseq")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "dada2:"; Rscript -e 'packageVersion("dada2")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'

bcftools --version | head -2
echo -n "DEXSeq:"; Rscript -e 'packageVersion("DEXSeq")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "DiffBind:"; Rscript -e 'packageVersion("DiffBind")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "ChIPQC:"; Rscript -e 'packageVersion("ChIPQC")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
bedtools --version
picard --version | sed 's/[^ ]*//'
echo -n "cython: "; python3 -c 'import cython as cyt; print(cyt.__version__)'
echo -n "cykhash: "; python3 -c 'import cykhash as cyk; print(cyk.__version__)'
macs2 --version
echo -n "py2bit: "; python3 -c 'import py2bit as py2b; print(py2b.__version__)'
echo -n "pyBigWig: "; python3 -c 'import pyBigWig as pyB; print(pyB.__version__)'
echo -n "pysam: "; python3 -c 'import pysam as pys; print(pys.__version__)'
idr --version | head -1

echo -n "KODAMA:"; Rscript -e 'packageVersion("KODAMA")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "RNASeqR:"; Rscript -e 'packageVersion("RNASeqR")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "mixOmics:"; Rscript -e 'packageVersion("mixOmics")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "dupRadar:"; Rscript -e 'packageVersion("dupRadar")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
echo -n "gsalib:"; Rscript -e 'packageVersion("gsalib")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
