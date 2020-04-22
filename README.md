# MPRNA (Massively Parallel RNa Assay)

## 1. Fastq

## 2. Alignment / barcode-mapping

`zcat SAMPLE.fastq.gz | barcodeAligner2.py -i stdin -o ALIGNMENT_OUTPUT.txt -g OligoPool.fa`

## 3. Read-counting per barcodes

`barcodeCounter.R -n SAMPLE -f ALIGNMENT_OUTPUT.txt -b BARCODE_NUMBER -t COUNTING_THRESHOLD -o COUNT_OUTPUT_PREFIX`

barcodeCount.R generates COUNT_OUTPUT_PREFIX.Rdata for downstream analysis.

barodeCounter.R can take multiple samples at a time. In the case, all the arguments except for -o will be comma-separated. For example,\
`barcodeCounter.R -n SAMPLE1,SAMPLE2 -f ALIGNMENT1_OUTPUT.txt,ALIGNMENT2_OUTPUT.txt -b BARCODE_NUMBER1, BARCODE_NUMBER2 -t COUNTING_THRESHOLD1,COUNTING_THRESHOLD2 -o COUNT_OUTPUT_PREFIX`


