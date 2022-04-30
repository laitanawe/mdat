#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//This pipeline performs Raw Gene Expression Quantification (RGEQ)

//Include the ffg modules into the main script:

include { htseq_count } from './modules/htseq_count.nf'

workflow {

  ch0 = channel.of('count2', params.BAM_DIR)
  ch0.subscribe({ println("ch0: $it") })

  ch1 = channel.of(params.counts_dir, params.bam_files)
  ch1.subscribe({ println("ch1: $it") })

  ch_bam = Channel.fromPath("${params.dir_bam}*.coordsort.bam")
  ch_bam.subscribe({ println("ch_bam: $it") })

  ch2 = htseq_count(ch_bam)
  ch2.subscribe({ println("ch2: $it") })
}
