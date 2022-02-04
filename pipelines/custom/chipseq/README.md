# ChIP-seq basics

---

### How to run a basic ChIP-seq analysis

---

### Table of contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Check sequence quality](#check-sequence-quality)
- [Generate mapping index](#generate-mapping-index)
- [Map reads](#map-reads)
- [Filter mappings](#filter-mappings)
- [Find peaks](#find-peaks)

---

### Introduction

This primer describes the steps involved in a basic chromatin immunoprecipitation sequencing (ChIP-seq) analysis.

What's ChIP-seq?

[Table of contents](#table-of-contents)

---

### Overview

The steps we will go over are:

  1) Check quality of read sequences in fastq files using FASTQC.
  2) Generate a Bowtie2 mapping index from a genomic sequence fasta file.
  3) Map reads to the genome using the genome index.
  4) Filter poor quality and ambiguous mappings.
  5) Find clusters of read mappings (peaks).

[Table of contents](#table-of-contents)

---

### Requirements

This primer describes the steps involved in a basic bulk ChIP-seq analysis.

[Table of contents](#table-of-contents)

---

### Inputs

This primer describes the steps involved in a basic bulk ChIP-seq analysis.

Header lines:

```
gunzip -c SRR1635435.fastq.gz > tmp.fastq

## get header line: NR is total number of records seen thus far:
## InstrumentID:RunNumber:FlowcellID:LaneNumber:TileNumber:XCoord:YCoord
## HWI-ST1378  :51       :C14U9ACXX :6         :1101      :2432  :2046

$ awk 'NR % 4 == 1' tmp.fastq | head
@SRR1635435.1 HWI-ST1378:51:C14U9ACXX:6:1101:2627:2077 length=50
@SRR1635435.2 HWI-ST1378:51:C14U9ACXX:6:1101:2934:2087 length=50
@SRR1635435.3 HWI-ST1378:51:C14U9ACXX:6:1101:3235:2052 length=50
@SRR1635435.4 HWI-ST1378:51:C14U9ACXX:6:1101:3190:2157 length=50
@SRR1635435.5 HWI-ST1378:51:C14U9ACXX:6:1101:3237:2229 length=50
@SRR1635435.6 HWI-ST1378:51:C14U9ACXX:6:1101:3471:2044 length=50
@SRR1635435.7 HWI-ST1378:51:C14U9ACXX:6:1101:3388:2044 length=50
@SRR1635435.8 HWI-ST1378:51:C14U9ACXX:6:1101:3348:2061 length=50
@SRR1635435.9 HWI-ST1378:51:C14U9ACXX:6:1101:3479:2092 length=50
@SRR1635435.10 HWI-ST1378:51:C14U9ACXX:6:1101:3444:2146 length=50

## An Illumina flow cell has 8 lanes; each lane has 2 columns;
##   each column has up to 50 tiles; each tile imaged 4 times per
##   cycle, or one image per base. In newer sequencers (e.g. Novaseq
##   etc), continuous movie taken of whole flow cell, so no tiles.
```

Quality scores:

```
## ascii = 33 + Q; Q = -10 * log10(P.wrong)
## ascii = 33 - 10 * log10(P.wrong)
## P.wrong = 10^((33 - ascii) / 10)

x=$(printf '%d' "'B")

$ Rscript -e "10^((33 - $x) / 10)"
[1] 0.0005011872

x=$(printf '%d' "'#")

$ Rscript -e "10^((33 - $x) / 10)"
[1] 0.6309573

x=$(printf '%d' "'1")

$ Rscript -e "10^((33 - $x) / 10)"
[1] 0.02511886

x=$(printf '%d' "'H")

$ Rscript -e "10^((33 - $x) / 10)"
[1] 0.0001258925
```

[Table of contents](#table-of-contents)

---

### Check sequence quality

We will use the FastQC program to check for commonly encountered sequence quality issues. The concept of k-mers and their use.

```
singularity exec "$container" fastqc -t $cpus "$reads_fastq"
```

Execute batch script:

```
find . -name '*.fastq.gz' | while read f; do
  echo -n "$f: "
  sbatch ~/dev/crf_chipseq_basics/fastqc.sh "$f"
done
```

[Table of contents](#table-of-contents)

---

### Generate mapping index

For our example data, we can map our reads to the human genome using Bowtie2. Mapping programs typically need the target sequence (typically a genome or transcriptome) converted into an index structure prior to mapping the reads. The index encodes the locations of each distinct k-mer present within the genome. This is used during the mapping stage, which processes each read by splitting it into k-mers, and looking for clustered locations of matching k-mers in the genome by doing a k-mer lookup using the index. Subsequently the mapper tries to stitch together any clusters of matching k-mers in the genome by extending the read's alignment into regions between the k-mer matches, rewarding matches and penalizing mismatches/gaps, and summarizing the quality of each potential match within the genome.

```
singularity exec "$container" bowtie2-build \
  --threads $cpus \
  --seed $seed \
  "$genome_fasta" "$genome_index"
```

Execute batch script:

```
## contains genome fasta file Homo_sapiens.GRCh38.dna.primary_assembly.fa
d=/db/GRCh38/release-103/dna

## make a directory for bowtie2 index files:
mkdir $d/bowtie2

## build index:

sbatch ~/bowtie2idx.sh \
  $d/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  $d/bowtie2/grch38r103

## view results:


```

[Table of contents](#table-of-contents)

---

### Map reads

Map some reads here.

```
singularity exec "$container" bowtie2 \
  -p $cpus \
  --local \
  -x "$genome_index" \
  -U "$reads_fastq" \
  -S "$output_sam"
```

Mapping statistics. Expect mapping rates of 90% or higher, though lower may work out fine. This contrasts to RNA-seq experiments, where mapping rates are typically about 10-20% lower.

```
for f in *.sam; do
  singularity exec /projects/researchit/crf/containers/chipseq.sif samtools flagstat $f
done
```

[Table of contents](#table-of-contents)

---

### Filter mappings

Filter some mappings using samtools. some stuff on mapq. need output to be bam. samtools can do both in one step:

```
samtools view -@$cpus -q $mapq_cutoff -o "$output_bam" "$input_sam"
```

[Table of contents](#table-of-contents)

---

### Find peaks

Find some peaks with MACS3.

Effective genome size:
  `https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html`
  for human should be between 2.65Gb and 2.9Gb, depending on read/kmer size and genome
  version. For GRCh38 and 50mers w/o ambiguous mappings: 2701495761 bases.

If accepting non-unique mappings, set to number.bases - number.Ns; so just
  the total number of non-N bases.

If not accepting non-unique mappings, set to number of uniquely mappable positions,
  which is usually estimated as the number.unique.kmers * kmer.size; where
    kmer.size is set to read.length;

`$flag_B` sets the `-B` option, which outputs a bedgraph`

```
singularity exec "$container" macs3 callpeak $flag_B \
  -t ${chip_files[@]} \
  -c ${sham_files[@]}  \
  -n "$prefix_out" \
  -f $input_format \
  -g $genome_size \
  -q $fdr_cutoff
```

[Table of contents](#table-of-contents)

---
