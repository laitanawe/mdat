#!/usr/bin/env bash

#Data Download Tools: These packages come with the VM
  wget --version | head -1
  zip --version | head -2 | tail -1 | sed 's/^This is //'
  perl --version | head -2 | tail -1 | sed 's/^This is p/P/'
  python2 --version | head -2 | tail -1 | sed 's/^ *//'
  python3 --version | head -2 | tail -1 | sed 's/^ *//'
  echo -n "Tensorflow: "; python3 -c 'import tensorflow as tf; print(tf.__version__)' 2>&1 | tail -1
  echo -n "scikit-learn: "; python3 -c 'import sklearn as sk; print(sk.__version__)'
  echo -n "Scipy: "; python3 -c 'import scipy as sc; print(sc.__version__)' 2>&1 | tail -1
  echo -n "Pandas: "; python3 -c 'import pandas as pd; print(pd.__version__)' 2>&1 | tail -1
