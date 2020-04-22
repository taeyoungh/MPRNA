#! /usr/bin/Rscript --vanilla

suppressMessages(library("optparse"))

############################
# Usage
############################

#srun -n1 --mem=10G --pty bash

#barcodeCounter.R -n SAMPLE_NAME1,SAMPLE_NAME2 \
#                 -f SAMPLE_NAME1.alignment.txt,SAMPLE_NAME2.alignment.txt \
#                 -b BARCODE_NUMBER_PER_DESIGN,BARCODE_NUMBER_PER_DESIGN \
#                 -t NUMBER_OF_NUCLEOTIDES_THRESHOLD_COUNTING,NUMBER_OF_NUCLEOTIDES_THRESHOLD_COUNTING \
#                 -o OUTPUT_FILENAME_PREFIX

############################
# 1. Command line arguments
############################

# argument: -n sampleName1,sampleName2
# annotation file is expected to have a header.
# each line should have 6 columns of gene_name, chrom, txStart, txEnd, exonStarts, exonEnds, strand

# argument: output
# default is standard output.

option_list <- list(make_option(c("-n", "--sampleNames"), type="character", default=NULL, help="Comma-separated list of sample names"),
					          make_option(c("-f", "--alignmentFiles"), type="character", default=NULL, help="Comma-separated list of alignment file names"),
					          make_option(c("-b", "--barcodeNum"), type="character", default=NULL, help="The number of barcodes per a design"),
                    make_option(c("-t", "--threshold"), type="character", default=NULL, help="Match threshold for counting"),
					          make_option(c("-o", "--output"), type="character", default="count", help="Output file name"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$sampleNames) | is.null(opt$alignmentFiles) | is.null(opt$barcodeNum) | is.null(opt$threshold)) {
	stop("One of required arguments is missing\n", call.=FALSE)
}

############################
# 2. Preprocessing
############################

sampleList <- unlist(strsplit(opt$sampleNames, split=","))
alignmentFiles <- unlist(strsplit(opt$alignmentFiles, split=","))
barcodeNum <- setNames(as.numeric(unlist(strsplit(opt$barcodeNum, split=","))), sampleList)
threshold <- setNames(as.numeric(unlist(strsplit(opt$threshold, split=","))), sampleList)

############################
# 3. Read alignments
############################
align <- list()

for (i in 1:length(sampleList)) {
  cat(paste0("Reading ", alignmentFiles[i], "\n"))
	temp <- read.table(alignmentFiles[i], header=F, sep="\t", colClasses = c("NULL", "character", "character", "numeric")) # don't read qname.
	colnames(temp) <- c("Barcode", "DesignID", "MatchNum")
	align[[sampleList[i]]] <- temp
}

############################
# 4. QC
############################
align.qc <- list()

### Total number of reads
align.qc[["Total"]] <- sapply(align, nrow)
cat("==Total number of reads==\n")
print(align.qc[["Total"]])

### Number of unmapped reads
align.qc[["Unmapped"]] <- sapply(align, function(x) nrow(subset(x, DesignID=="Unmapped")))
cat("==Proportion of unmapped reads==\n")
print(align.qc[["Unmapped"]] / align.qc[["Total"]])

### Number of mapped
temp <- lapply(align, function(x) as.data.frame(table(x$MatchNum)))
alignNum.mapped <- do.call(rbind, temp)
colnames(alignNum.mapped) <- c("MatchNum", "Number")
alignNum.mapped$MatchNum <- as.numeric(as.character(alignNum.mapped$MatchNum))
alignNum.mapped$Sample <- sapply(strsplit(rownames(alignNum.mapped), split="[.]"), "[[", 1)

alignNum.mapped$Prop <- NA
for (s in sampleList) {
	idx <- which(alignNum.mapped$Sample==s)
	alignNum.mapped$Prop[idx] <- alignNum.mapped$Number[idx]/align.qc[["Total"]][s]
}

align.qc[["Mapped"]] <- alignNum.mapped

cat("==Proportion of mapped reads==\n")
for (s in sampleList) {
  cat(paste0(s,"\n"))
  print(subset(align.qc[["Mapped"]], Sample==s & MatchNum>=threshold[s]))
}

############################
# 5. Count matrix
############################
cat("\nGenerating count matrix\n")
countMat <- function(s) {
	temp <- subset(align[[s]], MatchNum>=threshold[s])
	temp <- split(temp, temp$DesignID)
	count <- t(sapply(temp, function(x) {y <- table(x$Barcode); c(y, rep(0,barcodeNum[s]-length(y)))}))
	colnames(count) <- paste("Barcode", 1:barcodeNum[s], sep=":")
	return(count)
}

count <- lapply(sampleList, countMat)
names(count) <- sampleList

############################
# 6. Output
############################

save(list=c("align.qc", "count"), file=paste0(opt$output, ".Rdata"))
