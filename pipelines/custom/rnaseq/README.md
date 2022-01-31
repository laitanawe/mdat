# RNA-seq basics

# WORK IN PROGRESS

---

### How to run a basic bulk RNA-seq analysis on Sumner

---

### Table of contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Check sequence quality](#check-sequence-quality)
- [Generate mapping index](#generate-mapping-index)
- [Map reads](#map-reads)
- [Generate expression matrix](#generate-expression-matrix)
- [Screen for outlier samples](#screen-for-outlier-samples)
- [edgeR](#edger)
- [limma trend](#limma-trend)
- [limma voom](#limma-voom)
- [DESeq2](#deseq2)

---

### Introduction

This primer describes the steps involved in a basic bulk RNA-seq analysis.

[Table of contents](#table-of-contents)

---

### Overview

The steps we will go over are:
  1) Check quality of read sequences in fastq files using FASTQC.
  2) Generate a STAR mapping index from a genomic sequence fasta file.
  3) Map read sequences in fastq files to STAR genome index, producing .bam output files.
  4) Generate a raw expression matrix from the STAR mapping output.
  5) Screen for outlier samples.
  6) Do a differential expression analysis using edgeR.
  7) Do a differential expression analysis using limma trend.
  8) Do a differential expression analysis using limma voom.
  9) Do a differential expression analysis using DESeq2.

[Table of contents](#table-of-contents)

---
### Requirements

The scripts in this repository for the first four steps (up to and including expression matrix generation) assume you have at least 24 threads and 64 GB of RAM available on a compute node. Those scripts are further written to take advantage of the SLURM workload manager as it is set up. With slight modification, the scripts can be run in other Linux environments. Steps 5 through 9 can be carried out or more typically (except for very large datasets) on your laptop. If you decide to do those steps on your laptop, you will need to install R along with the Bioconductor packages `edgeR`, `limma`, and `DESeq2`. The dependencies for all the steps (including R) are present in the Singularity container:

```
/containers/rnaseq.sif
```

A list of software versions in the container can be obtained by running:

```
module load singularity
singularity run /containers/rnaseq.sif
```

If you have R on your laptop, you can install the required R packages with the following commands:

```
Install commands here
```

[Table of contents](#table-of-contents)

---
### Inputs

The normal input files are fastq formatted files containing reads and fasta file containing genomic sequence. Also need a metadata file connecting the fastq files to experimental variables. Fastq file may be gzipped. Metadata file typically tab-delimited. Example input data and metadata file in:

```
/data/rnaseq/
```

[Table of contents](#table-of-contents)

---
### Check sequence quality

We will use the FastQC program to check for commonly encountered sequence quality issues. The concept of k-mers and their use.

```
singularity exec "$container" fastqc -t $cpus "$reads"
```

[Table of contents](#table-of-contents)

---
### Generate mapping index

For our example data, we can map our reads to the human genome using STAR. Mapping programs typically need the target sequence (typically a genome or transcriptome) converted into an index structure prior to mapping the reads. The index encodes the locations of each distinct k-mer present within the genome. This is used during the mapping stage, which processes each read by splitting it into k-mers, and looking for clustered locations of matching k-mers in the genome by doing a k-mer lookup using the index. Subsequently the mapper tries to stitch together any clusters of matching k-mers in the genome by extending the read's alignment into regions between the k-mer matches, rewarding matches and penalizing mismatches/gaps, and summarizing the quality of each potential match within the genome.

```
singularity exec "$container" STAR \
  --runThreadN $cpus \
  --runMode genomeGenerate \
  --genomeDir "$dir_out" \
  --genomeFastaFiles "$genome_fasta"
```

[Table of contents](#table-of-contents)

---
### Map reads

We will map our reads to

```
singularity exec "$container" STAR \
  --runThreadN $ncpus \
  --genomeDir "$genome_idx" \
  --sjdbGTFfile "$genome_gtf" \
  --quantMode GeneCounts \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$prefix_out" \
  --readFilesIn "$reads_f" "$reads_r"
```

[Table of contents](#table-of-contents)

---
### Generate expression matrix

We will generate an expression matrix using featureCounts program from the Subreads R package. Available as a stand-alone program, so do not need to invoke R to run it.

```
singularity exec /containers/rnaseq.sif \
  featureCounts \
    -T $cpus \
    -Q $mapq_min \
    --minOverlap $min_overlap \
    --primary \
    --ignoreDup \
    -a "$gtf" \
    -o "$out_file" \
    ${files[*]}
```

[Table of contents](#table-of-contents)

---

### Screen for outlier samples

In-progress

```
x <- 1:1000
```

[Table of contents](#table-of-contents)

---


### edgeR

We will first do differential expression analysis using the R package edgeR.


Set up data:

```
setwd("C:/Users/user/data/rnaseq")

## read in counts:
x <- read.table("feature_counts.txt", skip=1, header=T, sep="\t", as.is=T)

## gene metadata in first 6 columns:
genes <- x[, 1:6]

## keep track of Geneids this way after delete Geneid column:
rownames(x) <- x$Geneid

## get rid of columns (including Geneids) that are not counts:
x <- x[, 7:ncol(x)]

## clean up sample names and allow indexing columns by meta$sra_id:
names(x) <- sub("Aligned.out.bam", "", names(x))

## read in sample metadata:
meta <- read.table("metadata.txt", header=T, sep="\t", as.is=T)
rownames(meta) <- meta$sra_id
meta <- meta[names(x), ]

## isolate data of interest:
i <- meta$disease %in% c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma')
meta <- meta[i, ]
rownames(meta) <- NULL

## set 'control' as first level of factor:
meta$disease <- factor(meta$disease,
  levels=c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma'))

## order columns of x so match rows of meta:
x <- x[, meta$sra_id]
```

2-way analysis:

```
y <- meta$disease == 'endometrioidadenocarcinoma'

> table(y, useNA='ifany')
FALSE  TRUE
   11    21

dat <- edgeR::DGEList(counts=x, group=y)
min.samples <- 3
min.counts <- 10

i.keep <- apply(x, 1, function(v) sum(v >= min.counts, na.rm=T) >= min.samples)

> table(i.keep, useNA='ifany')
FALSE  TRUE
38518 22148

dat <- dat[i.keep, , keep.lib.sizes=F]
dat <- edgeR::calcNormFactors(dat)
dat <- edgeR::estimateDisp(dat)
tst <- edgeR::exactTest(dat)            ## only works for 2 groups!!!

edgeR::topTags(tst)
```

3-way analysis:

```
dat <- edgeR::DGEList(counts=x, group=meta$disease)
min.samples <- 3
min.counts <- 10
i.keep <- apply(x, 1, function(v) sum(v >= min.counts, na.rm=T) >= min.samples)

> table(i.keep, useNA='ifany')
FALSE  TRUE
38518 22148

y <- factor(meta$disease,
  levels=c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma'))

des <- model.matrix(~ y)

dat <- edgeR::calcNormFactors(dat)
dat <- edgeR::estimateDisp(dat, design=des)
fit <- edgeR::glmQLFit(dat, design=des)

## are there any group differences?:
tst <- edgeR::glmQLFTest(fit, coef=2:3)
edgeR::topTags(tst)

## compare first level to second:
tst <- edgeR::glmQLFTest(fit, coef=2)
edgeR::topTags(tst)

## compare first level to third:
tst <- edgeR::glmQLFTest(fit, coef=3)
edgeR::topTags(tst)

## compare second level to third:
tst <- edgeR::glmQLFTest(fit, contrast=c(0,-1,1))
edgeR::topTags(tst)
```

[Table of contents](#table-of-contents)

---

### Limma trend

We will next do differential expression analysis using the R package limma and the 'trend' approach.

```
setwd("C:/Users/user/data/rnaseq")

x <- read.table("feature_counts.txt", skip=1, header=T, sep="\t", as.is=T)
genes <- x[, 1:6]
rownames(x) <- x$Geneid
x <- x[, 7:ncol(x)]
names(x) <- sub("Aligned.out.bam", "", names(x))

meta <- read.table("metadata.txt", header=T, sep="\t", as.is=T)
rownames(meta) <- meta$sra_id
meta <- meta[names(x), ]

i <- meta$disease %in% c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma')
meta <- meta[i, ]
rownames(meta) <- NULL

x <- x[, meta$sra_id]

## ensure that control condition is first factor level:

meta$disease <- factor(meta$disease,
  levels=c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma'))

dat <- edgeR::DGEList(counts=x)
min.samples <- 3
min.counts <- 10
i.keep <- apply(x, 1, function(v) sum(v >= min.counts, na.rm=T) >= min.samples)

dat <- dat[i.keep, , keep.lib.sizes=F]
dat <- edgeR::calcNormFactors(dat)

log.cpm <- edgeR::cpm(dat, log=T, prior.count=3)
des <- model.matrix(~meta$disease)

fit <- limma::lmFit(log.cpm, design=des)
fit <- limma::eBayes(fit, trend=TRUE)

## for 1 coef: moderated t-test; for >1 coef: moderated F-test:

> limma::topTable(fit, coef=1, number=5)
                   logFC  AveExpr        t      P.Value    adj.P.Val         B
ENSG00000166710 14.82690 14.82366 151.2523 7.881485e-63 1.745591e-58 103.52432
ENSG00000210082 17.27087 17.34214 144.4361 6.319143e-62 6.997818e-58 103.00644
ENSG00000211459 14.68722 14.77767 120.4374 2.297503e-58 1.696170e-54 100.59242
ENSG00000198888 15.38447 15.46347 113.0165 4.041692e-57 2.237885e-53  99.59527
ENSG00000198727 13.36827 13.46238 110.6800 1.036383e-56 4.084926e-53  99.24942

> limma::topTable(fit, coef=2, number=5)
                    logFC   AveExpr         t     P.Value adj.P.Val         B
ENSG00000176200  2.196273 0.3606955  3.400193 0.001416266 0.9998574 -4.506011
ENSG00000005882 -2.268117 2.2841232 -3.290980 0.001940582 0.9998574 -4.511111
ENSG00000145982 -1.483557 3.5725335 -3.245564 0.002208909 0.9998574 -4.513221
ENSG00000256618  2.452429 1.3115048  3.216633 0.002397766 0.9998574 -4.514561
ENSG00000137878 -1.363261 6.0449029 -3.168569 0.002745678 0.9998574 -4.516780

> limma::topTable(fit, coef=3, number=5)
                    logFC    AveExpr         t     P.Value adj.P.Val         B
ENSG00000204710 -2.524127 3.22903165 -3.394811 0.001438586 0.9005333 -4.518016
ENSG00000187714  1.483744 0.03530998  3.011060 0.004249252 0.9005333 -4.533387
ENSG00000256618  2.450653 1.31150478  2.990140 0.004499121 0.9005333 -4.534209
ENSG00000279267 -2.150084 2.07019454 -2.700728 0.009701401 0.9005333 -4.545356
ENSG00000166311 -1.391009 3.79976215 -2.687312 0.010042778 0.9005333 -4.545861

> limma::topTable(fit, coef=2:3, number=5)
                meta.diseasemucinousadenocarcinoma meta.diseaseserousadenocarcinoma   AveExpr        F     P.Value adj.P.Val
ENSG00000256618                         2.45242875                        2.4506527 1.3115048 7.997283 0.001061944 0.9710095
ENSG00000176200                         2.19627325                       -0.2654143 0.3606955 6.395720 0.003580032 0.9710095
ENSG00000204710                         0.04046983                       -2.5241275 3.2290316 6.064302 0.004642337 0.9710095
ENSG00000005882                        -2.26811693                        0.2100621 2.2841232 5.900531 0.005284262 0.9710095
ENSG00000095637                         1.88740074                       -1.4410516 2.3588357 5.878838 0.005375995 0.9710095
```

[Table of contents](#table-of-contents)

---

### Limma voom

We will next do differential expression analysis using the R package limma and the 'voom' approach.

```
rm(list=ls())

setwd("C:/Users/user/data/rnaseq")

x <- read.table("feature_counts.txt", skip=1, header=T, sep="\t", as.is=T)
genes <- x[, 1:6]
rownames(x) <- x$Geneid
x <- x[, 7:ncol(x)]
names(x) <- sub("Aligned.out.bam", "", names(x))

meta <- read.table("metadata.txt", header=T, sep="\t", as.is=T)
rownames(meta) <- meta$sra_id
meta <- meta[names(x), ]

i <- meta$disease %in% c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma')
meta <- meta[i, ]
rownames(meta) <- NULL

x <- x[, meta$sra_id]

## ensure that control condition is first factor level:

meta$disease <- factor(meta$disease,
  levels=c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma'))

dat <- edgeR::DGEList(counts=x)
min.samples <- 3
min.counts <- 10
i.keep <- apply(x, 1, function(v) sum(v >= min.counts, na.rm=T) >= min.samples)

> table(i.keep, useNA='ifany')
FALSE  TRUE
38518 22148

dat <- dat[i.keep, , keep.lib.sizes=F]
dat <- edgeR::calcNormFactors(dat)
des <- model.matrix(~meta$disease)

## mean vs. variance plot:
v <- limma::voom(dat, design=des, plot=T)

fit <- limma::lmFit(v, design=des)
fit <- limma::eBayes(fit)

## for 1 coef: moderated t-test; for >1 coef: moderated F-test:

> limma::topTable(fit, coef=1, number=5)
                   logFC  AveExpr        t      P.Value    adj.P.Val        B
ENSG00000075624 14.24467 14.08730 61.30493 9.980136e-47 1.105200e-42 94.74364
ENSG00000204592 12.69477 12.57254 61.37031 9.495292e-47 1.105200e-42 94.69762
ENSG00000210082 17.27324 17.34214 60.32781 2.113970e-46 1.560673e-42 94.15766
ENSG00000166710 14.82770 14.82363 59.34913 4.537087e-46 2.512185e-42 93.42557
ENSG00000087086 11.69968 11.55681 58.98954 6.025411e-46 2.669016e-42 92.99029

> limma::topTable(fit, coef=2, number=5)
                    logFC   AveExpr         t     P.Value adj.P.Val         B
ENSG00000071553 -1.459662  5.480389 -3.241886 0.002178705 0.9997236 -4.461950
ENSG00000075624 -1.133079 14.087303 -2.281022 0.027093744 0.9997236 -4.472515
ENSG00000182551 -1.291163  6.773107 -2.585222 0.012873379 0.9997236 -4.489661
ENSG00000140396 -1.844712  5.155285 -3.005510 0.004233042 0.9997236 -4.495302
ENSG00000167996 -1.146207 12.223549 -2.115996 0.039636790 0.9997236 -4.496916

> limma::topTable(fit, coef=3, number=5)
                   logFC  AveExpr        t      P.Value adj.P.Val         B
ENSG00000127951 2.654252 3.922963 3.628703 0.0006970853 0.6453993 -4.426928
ENSG00000139318 2.039336 4.104504 3.402679 0.0013668992 0.6453993 -4.446585
ENSG00000179344 1.794213 3.563195 3.628161 0.0006982295 0.6453993 -4.446663
ENSG00000163131 1.905036 6.512962 2.898385 0.0056695075 0.6453993 -4.453638
ENSG00000182718 1.507547 5.057839 3.117518 0.0031001287 0.6453993 -4.455375

> limma::topTable(fit, coef=2:3, number=5)
                meta.diseasemucinousadenocarcinoma meta.diseaseserousadenocarcinoma     AveExpr         F      P.Value adj.P.Val
ENSG00000053918                         -2.1854512                         2.493383  0.61068439 12.921674 3.332356e-05 0.7380502
ENSG00000120738                         -3.1621053                         1.910497  0.86931824  9.279837 3.982218e-04 0.7707887
ENSG00000214063                         -3.1963715                         0.603909  0.71337640  9.106378 4.512020e-04 0.7707887
ENSG00000197168                          2.9678828                        -1.766876 -0.03622534  9.051885 4.693229e-04 0.7707887
ENSG00000128268                         -0.8406509                        -3.132850 -0.87368405  8.904265 5.223256e-04 0.7707887

```

[Table of contents](#table-of-contents)

---

### DESeq2

Finally, we will do differential expression analysis using the R package DESeq2.

```
rm(list=ls())

setwd("C:/Users/user/data/rnaseq")

x <- read.table("feature_counts.txt", skip=1, header=T, sep="\t", as.is=T)
genes <- x[, 1:6]
rownames(x) <- x$Geneid
x <- x[, 7:ncol(x)]
names(x) <- sub("Aligned.out.bam", "", names(x))

meta <- read.table("metadata.txt", header=T, sep="\t", as.is=T)
rownames(meta) <- meta$sra_id
meta <- meta[names(x), ]

i <- meta$disease %in% c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma')
meta <- meta[i, ]
rownames(meta) <- NULL

meta$disease <- factor(meta$disease,
  levels=c('endometrioidadenocarcinoma', 'mucinousadenocarcinoma', 'serousadenocarcinoma'))

x <- x[, meta$sra_id]

dat <- DESeq2::DESeqDataSetFromMatrix(countData=x, colData=meta, design=~disease)

## ensure that control condition is first factor level; the actual hypothesis testing
##   is first level vs last!!! Does not seem to have an equivalent of omnibus test:

> levels(dat@colData@listData$disease)
[1] "endometrioidadenocarcinoma" "mucinousadenocarcinoma"     "serousadenocarcinoma"      

## equivalent to estimateSizeFactors() -> estimateDispersions() -> nbinomWaldTest():
dat <- DESeq2::DESeq(dat)

## 'independent filtering' conducted by default:

res <- results(dat, alpha=0.1)
res <- res[order(res$pvalue), ]

## png("ma_20211216a.png", width=720, height=480)
plotMA(res, ylim=c(-8, 8))
## dev.off()

out <- data.frame(res@listData)
out <- cbind(gene_id=res@rownames, out)

> dim(out)
[1] 60666     7

> head(out)
          gene_id  baseMean log2FoldChange    lfcSE      stat       pvalue      padj
1 ENSG00000128268  4.059533      -5.503071 1.331413 -4.133255 3.576610e-05 0.7065867
2 ENSG00000266718 15.756166      -7.059842 1.725149 -4.092309 4.270998e-05 0.7065867
3 ENSG00000279207  8.545335      -4.359412 1.091790 -3.992903 6.526917e-05 0.7065867
4 ENSG00000287336 13.341784      -4.403653 1.172537 -3.755662 1.728838e-04 0.7065867
5 ENSG00000004939 40.742947      -3.997542 1.066019 -3.749973 1.768536e-04 0.7065867
6 ENSG00000261272 15.043408      -5.836269 1.598729 -3.650568 2.616613e-04 0.7065867

> summary(out$padj)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
  0.707   0.767   0.962   0.878   0.962   1.000    7663

res <- results(dat, contrast=c('disease', 'endometrioidadenocarcinoma', 'mucinousadenocarcinoma'), alpha=0.1)
summary(res$padj)

res <- results(dat, contrast=c('disease', 'endometrioidadenocarcinoma', 'serousadenocarcinoma'), alpha=0.1)
summary(res$padj)

res <- results(dat, contrast=c('disease', 'mucinousadenocarcinoma', 'serousadenocarcinoma'), alpha=0.1)
summary(res$padj)
```

[Table of contents](#table-of-contents)

---
