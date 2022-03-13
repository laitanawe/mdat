#!/usr/bin/env bash

#Data Download Tools: These packages come with the VM
  wget --version | head -1
  zip --version | head -2 | tail -1 | sed 's/^This is //'
  perl --version | head -2 | tail -1 | sed 's/^This is p/P/'
  python2 --version | head -2 | tail -1 | sed 's/^ *//'
  python3 --version | head -2 | tail -1 | sed 's/^ *//'
  esearch --help | head -1

#Bionformatics Tools:
  R --version | head -n 1
  echo -n "Bioconductor:"; Rscript -e 'BiocManager::version()' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the Bioconductor version (formatted)
