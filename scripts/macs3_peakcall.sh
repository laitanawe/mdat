#!/usr/bin/env bash
#SBATCH -J macs3_peakcall #Job Name
#SBATCH	-N 1 #Number of nodes
#SBATCH	-q batch #QOS
#SBATCH	-t 72:00:00 #Time
#SBATCH	-c 24 #CPUs
#SBATCH	--mem 64g #memory
#SBATCH	-o './%x.%j.out' #output file; %x is Job Name, %j is Job number
#SBATCH -e './%x.%j.err' #error file; %x is Job Name, %j is Job number

module load singularity

cpu=24
echo "cpu:$cpu"

ram=64
echo "ram:$ram"

input_chipfile=$1
echo "input_chipfile: $input_chipfile"

container="/container_path/chipseq.sif"
echo "container: $container"

shamfile="/path/testdata/shamfile.cfg"

singularity exec $container macs3 callpeak -B -t $(cat $input_chipfile) -c $(cat $shamfile) \
 -n "macs3" -f BAM -g "mm" -q 0.05 --outdir macs3_outfiles
echo "macs3 peak calling process completed for sorted bams in $input_chipfile!"

## From the command line
## sbatch ./macs3_peakcall.sh path_to_chipfile.cfg
## singularity exec $container macs3 callpeak -B -t chipfile1.sort.bam chipfile2.sort.bam -c shamfile.sort.bam \
## -m "macs3" -f BAM -g 'mm' -q 0.05 --outdir macs3_outfiles
## Use -g 'hs' for human genome, 'mm' for mouse, 'ds' for drosophila, 'ce' for C. elegans
