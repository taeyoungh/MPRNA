#! /usr/bin/python2

# ====================================
# Version
# ====================================

# [Author]
# Taeyoung Hwang, Ph.D.,
# taeyoung.hwang@colorado.edu

# [Last update]
# 2020-05-10
# This is a cleaned version of barcodeAligner_v3.py (RinnLab in-house tool).

# ====================================
# 0. Parameters and Functions
# ====================================

BARCODE_SIZE = 10 # nts

import string
import sys
import getopt

def revComp(seq):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = seq.translate(complements)[::-1]
	return rcseq

# ====================================
# 1. Input and ouput
# ====================================

opts,args = getopt.getopt(sys.argv[1:],"hg:i:o:",["help", "genome=", "input=", "output="])

for opt,val in opts:
	if opt=="-h":
		print("barcodeAligner.py -i INPUT.fastq/stdin -o OUTPUT.txt/stdout -g OligoPool.fa")
	if opt in ("-i","--input"):
		if val=="stdin":
			fin=sys.stdin
		else:
			fin=open(val)
	if opt in ("-o","--output"):
		if val=="stdout":
			fout=sys.stdout
		else:
			fout=open(val, "w")
	if opt in ("-g","--genome"):
		genomeFile=open(val)

ferr=sys.stderr

# =======================================
# 2. Read Oligo pool design file.
# =======================================

oligoPool = dict()

lineRaw= genomeFile.readline()
while lineRaw:

	line=lineRaw.strip()

	if line[0] == ">":

		# header
		temp = line[1:].split("_") # format is ">DesignID_0.AAACAAAGAA": DesignID_BarcodeID.Barcode
		designID = "".join(temp[:-1]) # to consider the case that designID has "_"
		barcodeID = temp[-1].split(".")[0]

		# sequence
		line = genomeFile.readline()
		temp = line.strip()
		seq = temp[16:(-17-BARCODE_SIZE)]
		barcode = temp[(-17-BARCODE_SIZE):-17]

		oligoPool[barcode] = [designID, barcodeID, seq]

	else:
		ferr.write("Error in oligoPool: "+lineRaw+"\n")
		sys.exit()

	lineRaw= genomeFile.readline()

genomeFile.close()

# ==============================================
# 3. Read and align reads by barcode
# ==============================================

lineNum=0

for lineRaw in fin:

	lineNum = lineNum + 1
	line = lineRaw.strip()

	temp = lineNum % 4
	if temp == 1 : # 1st line of fastq
		qname = line.split(" ")[0]
		continue
	elif temp == 2 : # 2nd line of fastq
		seq = line
		continue
	elif temp == 3 : # 3rd line of fastq
		continue
	else : # 4th line of fastq
		qual = line

	seq_rc = revComp(seq)
	barcodeRead = seq_rc[-BARCODE_SIZE:]
	designRead = seq_rc[:-BARCODE_SIZE]

	if barcodeRead in oligoPool :
		[designID, barcodeID, designSeq] = oligoPool[barcodeRead]
		mismatchNum = sum([1 for i in range(0, min(len(designSeq), len(designRead))) if designSeq[-i-1] != designRead[-i-1]])
		out = [qname, barcodeRead, barcodeID, designID, str(mismatchNum)]
	else :
		out = [qname, barcodeRead, "Unmapped", "Unmapped", str(-1)]

	fout.write("\t".join(out)+"\n")

ferr.write("The number of reads: " + str(lineNum/4) + "\n")

fin.close()
fout.close()
ferr.close()
