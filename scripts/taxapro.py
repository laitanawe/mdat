import sys
from Bio import Entrez
from collections import Counter
import pandas as pd
import glob
import os

###########################################
def get_tax_id(species):
    """to get data from NCBI taxomomy, we need to have the taxid.  We can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]
###############################################
def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)
###############################################

## Usage: python taxonomy_code.py *_1.fastq*_outputs/contigs.gff *_1.fastq*_outputs/contigs.tsv laitanawe@yahoo.com
## Usage: python taxapro.py (assuming that metagenome was assembled using microbetaxasars2 and taxapro is in same directory as fastqs)
Entrez.email = "laitanawe@yahoo.com"
in_file_list = sorted(glob.glob('*_1.fastq*_outputs/contigs.gff'))
print("InFileList: " + str(in_file_list))
out_file_list = [i[:-3] + 'tsv' for i in in_file_list]
print("\nOutFileList: " + str(out_file_list))

print("\nNumber of Whole Genomes to characterize: " + str(len(in_file_list)))

gffcounter = 0
for in_file in in_file_list:
    out_file = out_file_list[gffcounter]
    print('\nCurrently Reading Input File ..... ' + in_file)
    print('\nOutput File Name to Write: ' + out_file)
    with open (in_file,'r', encoding='utf-8') as inpt:
        sps = []
        content = inpt.readlines()
        print('Extracting nodes ... ')
        for i in content:
            if '# Model Data:' in i:
                sp = i.split ('|')[1].split('|')[0].replace('_',' ')
                sps.append(sp)
    print('Counting nodes .... ')
    counter = dict(Counter(sps))
    total = 0
    for i in counter:
        try:
            taxid = get_tax_id(i)
            taxdata = get_tax_data(taxid)
            total = total + counter[i]
        except:
            pass
    print ('Retrieving taxonomy ... ')
    tax = []
    for i in counter:
        try:
            taxid = get_tax_id(i)
            taxdata = get_tax_data(taxid)
            lineage = {d['Rank']:d['ScientificName'] for d in
            taxdata[0]['LineageEx'] if d['Rank'] in ['superkingdom','phylum','class','order', 'family','genus','species']}
            lineage['Strain']=i
            lineage['#Hits']=counter[i]
            lineage['Frequency (%)']=(counter[i]/total)*100
            tax.append(lineage)
        except:
            pass
    print('Total found sp:',len(tax))
    print('Total nodes: ',total)
    print('Writing output file >',out_file)
    df = pd.DataFrame(tax)
    df.to_csv(out_file,sep='\t',index=False)
    gffcounter += 1

print("Taxonomic and genus-level characterization completed for the entire WGS dataset.")
