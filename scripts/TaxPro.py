import argparse
from posixpath import join
import subprocess
import threading
from pathlib import Path
import sys
import os
from Bio import Entrez
from collections import Counter
import pandas as pd
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scripts import assembly, taxa_ann, taxonomy_frequency, ploting

parser = argparse.ArgumentParser(description='Give row data')
parser.add_argument('-f', '--forward', help='forward sequence')
parser.add_argument('-r', '--reverse', help=' reverse sequnece')
# parser.add_argument('-s','--sra', help='SRA Run Accession Number')
parser.add_argument('-o', '--output', help='Output folder')
arguments = parser.parse_args()
arguments = arguments.__dict__
file_path = os.path.realpath(__file__)
file_path = os.path.abspath(file_path)
file_path = os.path.dirname(file_path)
forward = arguments['forward']
reverse = arguments['reverse']
# sra = arguments['sra']
output = arguments['output']

if output == None:
    raise KeyError ("Output folder should be mentioned")
elif os.path.isdir(output) == False:
    raise TypeError("The output folder is not supported")
else:
    output = output.rstrip('/')
    output = os.path.abspath(output)
    if forward == None or reverse == None:
        pass
    else:
        try:
            assembly = threading.Thread(target=assembly.assemble_pe, args=(file_path,forward,reverse, output))
            assembly.start()
            print('assembly started')
            print('######################################################################################')
            assembly.join()
            print("Assembly finished")
            print('#######################################################################################')
            annotation = threading.Thread (target = taxa_ann.prodigal, args = (file_path,f"{output}/assembly", output))
            annotation.start()
            print("Taxonomy assignment started")
            print("#########################################################################################")
            annotation.join()
            print("Assignment finished")
            print ('#########################################################################################')
            Path(f"{output}/taxa_stat").mkdir(parents=True, exist_ok=True)
            taxa = threading.Thread (target = taxonomy_frequency.tax_to_freq, args = (f"{output}/prodigal/contigs.gff", f"{output}/taxa_stat/taxa_freq.csv"))
            taxa.start()
            taxa.join()
            data = pd.read_csv(f"{output}/taxa_stat/taxa_freq.csv",sep='\t')
            plots = threading.Thread (target = ploting.ploting, args = (data,f"{output}/taxa_stat"))
            plots.start()
            plots.join()
        except:
            raise ChildProcessError("Error in a step")
