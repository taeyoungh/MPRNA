# MPRNA (Massively Parallel RNa Assay)

## 0. Prerequsite

### Softwares
[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)\
[R](https://www.r-project.org/)\
[R package: optparse](https://cran.r-project.org/web/packages/optparse/index.html)\
Additionally, download all the programs in the folder "Programs" in this github.

### Oligo Pool Design file
Contact Rinn Lab.

## 1. Fastq

Check quality by running fastqc.

## 2. Alignment (mapping reads by barcode)

This step takes a sequencing read file and maps sequencing reads to the barcodes in an oligo pool.

`zcat SAMPLE.fastq.gz | barcodeAligner.py -i stdin -o ALIGNMENT_OUTPUT.txt -g OLIGO_POOL.fa`\
-o output file name\
-g oligo pool design file (fasta format)

**Example**\
`zcat maxi1.fastq.gz | barcodeAligner.py -i stdin -o maxi1_alignment.txt -g oligoPool.fa`

## 3. Counting reads per barcodes

This step generates COUNT_OUTPUT_PREFIX.Rdata for downstream analysis.

`barcodeCounter.R -n SAMPLE -f ALIGNMENT_OUTPUT.txt -g OLIGO_POOL.fa -t COUNTING_THRESHOLD -o COUNT_OUTPUT_PREFIX`\
-n sample name to be used in output\
-f output file of barcodeAligner.py\
-g oligo pool design file (fasta format)\
-t the threshold value of mismatched nucleotides between a read and a design for a read to be counted, for example -t 2 counts reads whose mismatches less than or equal to 2.\
-o output_prefix

**Example**\
`barcodeCounter.R -n maxi1 -f maxi1_alignment.txt -g oligoPool.fasta -t 2 -o maxi1_count`

**Note** barcodeCounter.R can take multiple samples at a time for the same oligo pool design. In the case, all the arguments except for -o and -g will be comma-separated. For example,\
`barcodeCounter.R -n SAMPLE1,SAMPLE2 -f ALIGNMENT1_OUTPUT.txt,ALIGNMENT2_OUTPUT.txt -g OligoPool.fa -t COUNTING_THRESHOLD1,COUNTING_THRESHOLD2 -o COUNT_OUTPUT_PREFIX`

## 4. Representation of oligos

This step is done with R. See "representation.md".

