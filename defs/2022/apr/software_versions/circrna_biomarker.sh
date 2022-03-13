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
  subread-align -version 2>&1 | sed '/^$/d'
  R --version | head -1
  echo -n "Bioconductor:"; Rscript -e 'BiocManager::version()' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the Bioconductor version (formatted)
  echo -n "edgeR:"; Rscript -e 'packageVersion("edgeR")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the edgeR version
  echo -n "limma:"; Rscript -e 'packageVersion("limma")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the limma version
