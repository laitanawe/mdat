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
  //example redirect output to a file.
  output: path('*.out')
  //example output to standard-out
  //output: stdout

  shell:

//If you're redirecting a channel out item to a file/path for each input channel item, during process execution, a path is created for come.out under separate work subdirectories when the command is executed in bash, and the work subdirs will also have the ffg. hidden files: .command.begin, .command.err, .command.log, .command.out, .command.run, .command.sh, .exitcode
//Note: Anything shell script that should go to standard out actually goes to .command.out
    '''
    #redirect output to a file your_groovy_code_goes_within_braces
    echo "!{word} my friend, I'm in '$PWD'" > !{word.toLowerCase()}.out
    #redirect output to standard-out
    echo "!{word} my friend, I'm in '$PWD'"
    '''
}

process check_sums {

  input: path(file_in)
  output: path('*.{md5,sha256,sha512}')

  shell:
    '''
    # can nest 1x or 2x quoted strings w/i triple single quoted multi-line string.
    echo "file_in: '!{file_in}'"
    md5sum '!{file_in}' > '!{file_in}.md5'
    sha256sum '!{file_in}' > '!{file_in}.sha256'
    sha512sum '!{file_in}' > '!{file_in}.sha512'
    '''
}

process check_sums_intermediate {

//a tuple is an immutable list where no addition or removal of items are allowed.
  input: val(word)
  output: tuple val(word), path('*.out')

  shell:
    '''
    #redirect output to a file your_groovy_code_goes_within_braces
    echo "!{word} my friend, I'm in '$PWD'" > !{word.toLowerCase()}.out
    #redirect output to standard-out
    echo "!{word} my friend, I'm in '$PWD'"
    '''
}

process check_sums_advanced {

  input: tuple val(word), path(file_in)
  output: tuple val(word), path('*.{txt,md5,sha256,sha512}')

  shell:
    '''
    echo "file_in: '!{file_in}'"
    md5sum '!{file_in}' > '!{file_in}.md5'
    sha256sum '!{file_in}' > '!{file_in}.sha256'
    sha512sum '!{file_in}' > '!{file_in}.sha512'
    echo "!{word} my friend, I'm in '$PWD'" > !{word.toLowerCase()}.txt
    '''
}

workflow {

//within workflow, use $it to refer to an item because there's not input variable name.
//within process, use the input variable name declared e.g. !{word} for say_it or !{file_in} for check_sums
  ch_in = channel.of('Come', 'Wa', 'Zo', 'Bia')
  ch_in.subscribe({ println("ch_in: $it") })

  ch_come = say_it(ch_in)
  ch_come.subscribe({ println("ch_come: $it") })

  ch_sums = check_sums(ch_come)
  ch_sums.subscribe({ println("ch_sums: $it\n") })

  ch_sums_intermediate = check_sums_intermediate(ch_in)
  ch_sums_intermediate.subscribe({ println("ch_sums_intermediate: $it\n") })

  ch_sums_advanced = check_sums_advanced(ch_sums_intermediate)
  ch_sums_advanced.subscribe({ println("ch_sums_advanced: $it\n") })

  // join channels on (by default) first element (grouping key) in each tuple:

  ch_join = ch_sums_advanced.join(ch_sums_intermediate)
  ch_join.subscribe({ println("ch_join: $it\n") })

  // ch_join tuples are tuple(val, list(path), path); here, convert to
  //   tuple(val, list(path)), by taking third element (the last path) and
  //   moving it into the second element (the list of paths). This will
  //   simplify definition, use and manipulation of downstream channels:

  // grab the key get(0), clone the list, join one item get(2) to it.
  ch_reformat = ch_join.map({
    key = it.get(0) // get the first item i.e. the common input/output variable between the two joined channels.
    val = it.get(1).clone() // make shallow copy of second item (a list). Cloning 2nd item (the list).
    val.add(it.get(2)) // We're appending this path/value as an item of the cloned list.
    return tuple(key, val) // Return a tuple of [key,[value]]
  })
  ch_reformat.subscribe({ println("ch_reformat: $it\n") })

}
