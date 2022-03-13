#!/usr/bin/env bash

#Data Download Tools: These packages come with the VM
  wget --version | head -1
  zip --version | head -2 | tail -1 | sed 's/^This is //'
  perl --version | head -2 | tail -1 | sed 's/^This is p/P/'
  python --version | head -2 | tail -1 | sed 's/^ *//'
  python3 --version | head -2 | tail -1 | sed 's/^ *//'
  fasterq-dump --version | head -2 | tail -1

  #Taxonomic Profiling Tools:
  R --version | head -n 1
  echo -n "Bioconductor:"; Rscript -e 'BiocManager::version()' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the Bioconductor version (formatted)
  echo -n "dada2:"; Rscript -e 'packageVersion("dada2")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
  echo -n "ggplot2:"; Rscript -e 'packageVersion("ggplot2")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
  echo -n "phyloseq:"; Rscript -e 'packageVersion("phyloseq")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
  echo -n "rRDP:"; Rscript -e 'packageVersion("rRDP")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
  echo -n "DECIPHER:"; Rscript -e 'packageVersion("DECIPHER")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'

  #Functional Analysis Tools:
  hmmalign -h | head -2 | tail -1 | sed 's/^# //'
  echo -n "castor:"; Rscript -e 'packageVersion("castor")' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g'
 # echo -n "gappa: "; gappa --version
 # epa-ng --version
  picrust2_pipeline.py --version
  head -n3 /opt/MinPath/MinPath.py | tail -1 | sed 's/MinPath/\/opt\/MinPath\/MinPath.py/'| sed 's/current //' | sed 's/^# //'
