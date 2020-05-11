#! /usr/bin/Rscript --vanilla
suppressMessages(library("optparse"))

# ====================================
# Version
# ====================================

# [Author]
# Taeyoung Hwang, Ph.D.,
# taeyoung.hwang@colorado.edu

# [Last update]
# 2020-05-10
# This is a cleaned version of barcodeCounter_v2.R (RinnLab in-house tool).

############################
# 1. Command line arguments
############################

option_list <- list(make_option(c("-n", "--sampleNames"), type="character", default=NULL, help="Comma-separated list of sample names"),
					          make_option(c("-f", "--alignmentFiles"), type="character", default=NULL, help="Comma-separated list of alignment file names"),
					          make_option(c("-g", "--oligoPool"), type="character", default=NULL, help="Fasta file of oligo pool"),
                    make_option(c("-t", "--threshold"), type="character", default=NULL, help="Mismatch threshold for counting"),
					          make_option(c("-o", "--output"), type="character", default="count", help="Output file name"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$sampleNames) | is.null(opt$alignmentFiles) | is.null(opt$oligoPool) | is.null(opt$threshold)) {
	stop("One of required arguments is missing\n", call.=FALSE)
}

############################
# 2. Preprocessing
############################

sampleList <- unlist(strsplit(opt$sampleNames, split=","))
alignmentFiles <- unlist(strsplit(opt$alignmentFiles, split=","))
threshold <- setNames(as.numeric(unlist(strsplit(opt$threshold, split=","))), sampleList)

############################
# 3. Read alignments
############################

# oligoPool
temp <- read.table(opt$oligoPool, sep="\n", header=F, stringsAsFactors=F)[,1]
temp <- strsplit(substr(temp[seq(1, length(temp), 2)], 2, 1000), split="_")
oligoPool <- data.frame(DesignID = sapply(temp, "[[", 1), BarcodeID = sapply(strsplit(sapply(temp, "[[", 2), split="[.]"), "[[", 1), stringsAsFactors=F)

# Alignment
align <- list()

for (i in 1:length(sampleList)) {
  cat(paste0("Reading ", alignmentFiles[i], "\n"))
	temp <- read.table(alignmentFiles[i], header=F, sep="\t", colClasses = c("NULL", "NULL", "character", "character", "numeric")) # don't read qname and barcodeSeq.
	colnames(temp) <- c("BarcodeID", "DesignID", "MismatchNum")
	align[[sampleList[i]]] <- temp
}

############################
# 4. QC
############################
align.qc <- list()

### Total number of reads
align.qc[["Total"]] <- sapply(align, nrow)
cat("\n==Total number of reads==\n")
print(align.qc[["Total"]])

### Number of unmapped reads
align.qc[["Unmapped"]] <- sapply(align, function(x) nrow(subset(x, DesignID=="Unmapped")))
cat("==Proportion of unmapped reads==\n")
print(align.qc[["Unmapped"]] / align.qc[["Total"]])

### Number of mapped reads
temp <- lapply(align, function(x) as.data.frame(table(x$MismatchNum)))
alignNum.mapped <- do.call(rbind, temp)
colnames(alignNum.mapped) <- c("MismatchNum", "Number")
alignNum.mapped$MismatchNum <- as.numeric(as.character(alignNum.mapped$MismatchNum))
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
  print(subset(align.qc[["Mapped"]], Sample==s & MismatchNum<=threshold[s]+5 & MismatchNum>=0))
	cat("\n")
}

############################
# 5. Count matrix
############################
cat("\nGenerating count matrix\n")

BARCODE_NUM <- length(unique(oligoPool$BarcodeID))

countMat <- function(s) {
	filtered <- subset(align[[s]], MismatchNum<=threshold[s] & MismatchNum>=0)

	# Count
	filtered <- split(filtered, filtered$DesignID)
	counted <- t(sapply(filtered, function(x) {temp <- setNames(rep(0,BARCODE_NUM), 1:BARCODE_NUM); y <- table(x$BarcodeID); temp[names(y)]<-y; return(temp) }))
	colnames(counted) <- paste("Barcode", 1:BARCODE_NUM, sep=":")

	# Augment missing oligo into count matrix
	temp <- unique(oligoPool$DesignID[! (oligoPool$DesignID %in% rownames(counted))])
	temp <- matrix(0, nrow=length(temp), ncol=BARCODE_NUM, dimnames = list(temp, colnames(counted)))
	out <- rbind(counted, temp)

	return(out)
}

count <- lapply(sampleList, countMat)
names(count) <- sampleList

############################
# 6. Output
############################

save(list=c("align.qc", "count"), file=paste0(opt$output, ".Rdata"))
