#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//Include / Import this module to the main script

//This script performs Raw Gene Expression Quantification on sorted BAM files

process htseq_count {

//tuple val(word), path('*.out')

//input: tuple val(mycounts), path(my_pattern)
input: val mycounts //, tuple val(my_pattern)
input: val(my_pattern)
output: tuple val(my_pattern), path('*.txt') //stdout //, tuple val(mycounts), tuple val(my_pattern), path('*.txt')
//mycounts = 'counts'

   shell:

    '''
    echo "Creating directory for Raw Gene Expression Quantification ongoing ... $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"
    echo "We will now create the directory !{mycounts} for RGEQ ..."
    mkdir -pv !{mycounts}

    #redirect output to standard-out
    echo !{mycounts} "is a val and is my first argument in the anonymous workflow, I'm in '$PWD'"
    echo !{my_pattern} "is a val and is my second argument in the anonymous workflow, I'm in '$PWD'"

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
            > !{mycounts}/${acc}.counts.txt 
      done
    matrix_list=$(ls -Ah !{mycounts}/*.counts.txt > !{mycounts.toLowerCase()}.gather.txt)
    echo
    echo "List of counts: "$matrix_list
    echo !{params.BAM_DIR} "is our BAM directory!"
    echo !{params.GTF_FILE} "is our GTF file!"
    echo
    echo "Raw Gene Expression Quantification process completed for all .BAM files now, $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"

    #redirect output to a file your_groovy_code_goes_within_braces
    echo "!{mycounts} my friend, I'm in '$PWD'" > !{mycounts.toLowerCase()}.echo.txt
    '''
 }
