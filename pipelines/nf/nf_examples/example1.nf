#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


list_of_num = [3,7,4,8]

counter = 0
for ( num : list_of_num ) {
counter += 1
println ("Item $counter is $num!")
}


process my_fun_process {

//Within the process, we can specify the input type and the name of the variable. The same applies to output.


  shell:

    // '$' variable interpolation:
    println("Interpolate: ")

  //baseDir = "$PWD"

  //params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"

}

/* Let us create a channel of paths that will be used for mapping .fq files to a reference genome

*/
workflow {

// create a channel of values which are strings or channel from paths which are tupples.
ch_input1 = channel.of('Drib09.48-.55','Drib10.00-.25','Drib11.16','DribDown12.49','DribUp1.35')
// We can inspect each item of our channel as follows:
// Each value in the channel gets assigned to the variable $it.
ch_input1.subscribe({ println("ch_input1: $it") })


//ch_input2 = channel.from( ${params.reads} )

// Process each item in parallel. The order will depend on runtime.
ch_output1 = my_fun_process(ch_input1)
ch_output1.subscribe({ println("ch_output1: $it") })

}
