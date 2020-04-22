# MPRNA (Massively Parallel RNa Assay)

## 0. Prerequsite

### Softwares
[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Additionally, download all the manual programs in the folder "Programs" in this github.

### Oligo Pool Design file
Contact Rinn Lab.

## 1. Fastq

Check quality by running fastqc.

## 2. Alignment / barcode-mapping

`zcat SAMPLE.fastq.gz | barcodeAligner2.py -i stdin -o ALIGNMENT_OUTPUT.txt -g OligoPool.fa`
-g oligo pool design file (fasta)

## 3. Read-counting per barcodes

barcodeCount.R generates COUNT_OUTPUT_PREFIX.Rdata for downstream analysis.

`barcodeCounter.R -n SAMPLE -f ALIGNMENT_OUTPUT.txt -b BARCODE_NUMBER -t COUNTING_THRESHOLD -o COUNT_OUTPUT_PREFIX`\
-n sample name\
-f output file of barcodeAligner2.py\
-b the number of barcodes per a design\
-t the threshold value of matched nucleotides between a read and a design for a read to be counted

**Note** barodeCounter.R can take multiple samples at a time. In the case, all the arguments except for -o will be comma-separated. For example,\
`barcodeCounter.R -n SAMPLE1,SAMPLE2 -f ALIGNMENT1_OUTPUT.txt,ALIGNMENT2_OUTPUT.txt -b BARCODE_NUMBER1, BARCODE_NUMBER2 -t COUNTING_THRESHOLD1,COUNTING_THRESHOLD2 -o COUNT_OUTPUT_PREFIX`


