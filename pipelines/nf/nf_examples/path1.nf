#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*********************************************************************************
This script demonstrates a problematic situation that can occur when a workflow
  branches and the output of the branches needs to be combined so that all
  outputs for a given input channel item can be brought back together.

In this case, we use the path output channel where files are generated that contain the word "come" in three different Nigerian languages (Yoruba, Hausa and Igbo) and English.
  Then several types of checksum files are generated for each language file. This
  script shows how a naive approach can result in ambiguities when trying to
  combine the language file channel with the checksum file channel so that
  the language file is correctly grouped with the corresponding checksum files.

*********************************************************************************/

process say_it {

  input: val word
  output: path 'come.out'

  shell:
    '''
    echo "!{word} my friend, I'm in '$PWD'" > come.out
    '''
}

process check_sums {

  input: path file_in
  output: path '*.{md5,sha256,sha512}'

  shell:
    '''
    # can nest 1x or 2x quoted strings w/i triple single quoted multi-line string.
    echo "file_in: '!{file_in}'"
    md5sum '!{file_in}' > '!{file_in}.md5'
    sha256sum '!{file_in}' > '!{file_in}.sha256'
    sha512sum '!{file_in}' > '!{file_in}.sha512'
    '''
}

workflow {

  ch_in = channel.of('come', 'wa', 'zo', 'bia')
  ch_in.subscribe({ println("ch_in: $it\n") })

  ch_come = say_it(ch_in)
  ch_come.subscribe({ println("ch_come: $it\n") })

  ch_sums = check_sums(ch_come)
  ch_sums.subscribe({ println("ch_sums: $it\n") })
}
