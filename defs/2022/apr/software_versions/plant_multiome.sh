#!/usr/bin/env bash

#Data Download Tools: These packages come with the VM
  wget --version | head -1
  zip --version | head -2 | tail -1 | sed 's/^This is //'
  perl --version | head -2 | tail -1 | sed 's/^This is p/P/'
  python2 --version | head -2 | tail -1 | sed 's/^ *//'
  python3 --version | head -2 | tail -1 | sed 's/^ *//'
  esearch --help | head -1

  echo -n "Tensorflow: "; python3 -c 'import tensorflow as tf; print(tf.__version__)' 2>&1 | tail -1
  echo -n "scikit-learn: "; python -c 'import sklearn as sk; print(sk.__version__)'
  echo -n "Scipy: "; python -c 'import scipy as sc; print(sc.__version__)' 2>&1 | tail -1
  echo -n "Pandas: "; python -c 'import pandas as pd; print(pd.__version__)' 2>&1 | tail -1
  echo -n "BioPython: "; python -c 'import Bio as bio; print(bio.__version__)'
  echo -n "Seaborn: "; python -c 'import seaborn as sea; print(sea.__version__)'
  echo -n "Matplotlib: "; python -c 'import matplotlib as mat; print(mat.__version__)'
  blastx -h | tail -3 | head -1 | sed 's/^ *Translated Query-Protein Subject B/B/'

#Bionformatics Tools:
  R --version | head -n 1
  echo -n "Bioconductor:"; Rscript -e 'BiocManager::version()' | sed 's/\[1]//' | sed 's/‘//g' | sed 's/’//g' ## Display the Bioconductor version (formatted)
