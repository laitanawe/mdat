#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

list_of_strings = ['Dribble09.48-.55','Dribble10.00-.25','Dribble11.16','Down12.49','Up1.35','Dribble2.02-.38-.55','Dribble3.08-.50']
list_of_num = [3,7,4,8]

counter_s = 0
for ( string : list_of_strings ) {
counter_s += 1
println ("String $counter_s is $string!")
}

println()

counter = 0
for ( num : list_of_num ) {
counter += 1
println ("Item $counter is $num!")
}
