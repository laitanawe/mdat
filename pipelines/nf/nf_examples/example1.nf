#!/usr/bin/env nextflow

list_of_num = [3,7,4,8]

counter = 0
for ( num : list_of_num ) {
counter += 1
println ("Item $counter is $num!")
}
