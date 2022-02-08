#!/usr/bin/env bash
#Short runtime fastq -> fasta conversion: (from rtliu (https://www.biostars.org/u/6324/))
sed -n '1~4s/^@/>/p;2~4p' my.fastq > my.fasta
