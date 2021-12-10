#!/usr/bin/env bash
#SBATCH -J gunzip #Job Name
#SBATCH	-N 1 #Number of nodes
#SBATCH	-q batch #QOS
#SBATCH	-t 72:00:00 #Time
#SBATCH	-c 24 #CPUs
#SBATCH	--mem 64g #memory
#SBATCH	-o './%x.%j.out' #output file; %x is Job Name, %j is Job number
#SBATCH -e './%x.%j.err' #error file; %x is Job Name, %j is Job number

filein=$1
echo "filein:$filein"
gunzip $filein

## From the command line, it unzips each file in the URL list
## for f in `ls $path_to_zipped_files`; do echo "unzipping $f" && ./gunzip.sh $f; done

## You can also use sbatch for job submission
## for f in `ls $path_to_zipped_files`; do echo "unzipping $f" && sbatch ./gunzip.sh $f; done
