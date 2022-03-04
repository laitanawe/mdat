#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


list_of_num = [3,7,4,8]

counter = 0
for ( num : list_of_num ) {
counter += 1
println ("Item $counter is $num!")
}

shell:

  // '$' variable interpolation:
  println("Interpolate: $word")

baseDir = "$PWD"

params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"

/* Let us create a channel of paths that will be used for mapping .fq files to a reference genome

*/
workflow {
ch_input = []

}
