import glob
import os
import gc

def assemble (forward,reverse,output):
    return os.system('metaspades.py -1 '+forward+' -2 '+reverse+' -t 16 -m 64 -o '+output)
def prodigal (output):
    return os.system('prodigal -a '+output+'/contigs.aa.fast -d '+output+'/contigs.nuc.fasta -i '+output+'/contigs.fasta -f gff -p meta > '+output+'/contigs.gff')
def meta (output):
    return os.system('metaphlan '+output+'/contigs.nuc.fasta --input_type fasta -o '+output+'/profiled_metagenome.txt')

# Python microbetaxasars2.py (assuming all fastqs are in same dir as python microbetaxasars2)
#Read in all forward and reverse sequences -- CHANGE FILEPATH IF NEEDED#
forward_list=sorted(glob.glob('*_1.fastq*'))
reverse_list=sorted(glob.glob('*_2.fastq*'))
output_list=[]
for files in forward_list:
	sample=files.split('_1.fastq.gz')[0]
	output_list.append(sample+'_outputs')

#start WGS analysis
for i in range(len(forward_list)):
	#assemble the forward and reverse sequences
	forward=forward_list[i]
	reverse=reverse_list[i]
	output=output_list[i]
	print('STARTING METAGENOME ASSEMBLY WITH METASPADES...')
	assemble(forward,reverse,output)
	print('METAGENOME ASSEMBLY COMPLETE...')
	gc.collect()

	#annotate using prodigal
	print('STARTING METAGENOME ANNOTATION WITH PRODIGAL...')
	prodigal(output)
	print('METAGENOME ANNOTATION COMPLETE...')
	gc.collect()

	#taxonomy characterization needs tot be done
	gc.collect()
	print ('TAXONOMY CHARACTERIZATION NEEDS TO BE DONE!')
