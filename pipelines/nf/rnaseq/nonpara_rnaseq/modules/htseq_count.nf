#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//This script performs Raw Gene Expression Quantification (RGEQ) on sorted BAM files

//This module should be included into the main script in order to do RGEQ

process htseq_count {

//input: tuple val(mycounts), path(my_pattern)
//input: val mcounts //, tuple val(my_pattern)

input: val(my_pattern)
output: tuple val(my_pattern), path('*.counts.txt') //stdout //, tuple val(mycounts), tuple val(my_pattern), path('*.txt')

   shell:

    '''
    echo "Creating directory for Raw Gene Expression Quantification ongoing ... $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"
    echo "We will now create the directory !{params.counts_dir} for RGEQ ..."
    mkdir -pv test_!{params.counts_dir}

    #redirect output to standard-out
    echo !{params.counts_dir} "is a val and was initially my first argument in the anonymous workflow, I'm in '$PWD'"
    echo !{my_pattern} "is a val and is part of my second argument in the anonymous workflow, I'm in '$PWD'"

    n=0
    for sortedbam_file in !{my_pattern}; do echo; 
      ##++n works
      ((++n)) 
      ##n++ does not work. Why?
      #((n++)); 
      b=$(basename $sortedbam_file); acc=$(echo $b | sed 's/.coordsort.bam//');
      echo "I'm in "$PWD;
      echo "Now processing sample "$n": "$acc;
      htseq-count \
            -f bam \
            -r pos \
            -s no \
            -t exon \
            -i gene \
            $sortedbam_file \
            !{params.GTF_FILE} \
            > ${acc}.!{params.counts_dir}.txt 
      done
    matrix_list=$(ls -Ah *.!{params.counts_dir}.txt > !{params.counts_dir.toLowerCase()}.gather.txt)
    echo
    echo "List of counts: "$matrix_list
    echo !{params.BAM_DIR} "is our BAM directory!"
    echo !{params.GTF_FILE} "is our GTF file!"
    echo
    echo "Raw Gene Expression Quantification process completed for all .BAM files now, $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"

    #redirect output to a file your_groovy_code_goes_within_braces
    echo "!{params.counts_dir} my friend, I'm in '$PWD'" > !{params.counts_dir.toLowerCase()}.echo.txt
    '''
 }
