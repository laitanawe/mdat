#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//Include / Import the ffg modules to this main script
// include { htseq_count } from './modules/htseq_count.nf'
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

    #redirect output to a file your_groovy_code_goes_within_braces
    #echo "!{mycounts} my friend, I'm in '$PWD'" > !{mycounts.toLowerCase()}.txt
    echo "!{mycounts} my friend, I'm in '$PWD'" > !{mycounts.toLowerCase()}.echo.txt
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

    echo
    echo !{params.BAM_DIR} "is our BAM directory!"
    echo !{params.GTF_FILE} "is our GTF file!"
    echo
    echo "Raw Gene Expression Quantification process completed for all .BAM files now, $(date +%a) $(date +'%Y-%m-%d %H:%M:%S')"

    '''
 
 }

workflow {

  ch0 = channel.of('count2', params.BAM_DIR)
  ch0.subscribe({ println("ch0: $it") })

  ch1 = channel.of(params.counts_dir, params.bam_files)
  ch1.subscribe({ println("ch1: $it") })

  ch2 = htseq_count(params.counts_dir, params.bam_files)
  ch2.subscribe({ println("ch2: $it") })
}
