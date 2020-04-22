#! /usr/bin/python2

# ====================================
# Usage and Update
# ====================================

'''
Design file: -g DESIGN_FILE.fa
A header format is ">DesignID_BarcodeID.Barcode", for example, "">DesignID_0.AAACAAAGAA".
'''

# ====================================
# PARAMETERS
# ====================================

import string
import sys
import getopt

BARCODE_SIZE = 10

# ====================================
# Functions
# ====================================

def revComp(seq):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = seq.translate(complements)[::-1]
	return rcseq


# ====================================
# I/O
# ====================================

opts,args = getopt.getopt(sys.argv[1:],"hg:i:o:",["help", "genome=", "input=", "output="])

for opt,val in opts:
	if opt=="-h":
		print("barcodeAligner2.py -i INPUT.fastq/stdin -o OUTPUT.txt/stdout -g OligoPool.fa")
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
# 1. Read Oligo pool design file.
# =======================================

oligoPool = dict()

lineRaw= genomeFile.readline()
while lineRaw:

	line=lineRaw.strip()

	if line[0] == ">":

		# header
		temp = line[1:].split("_") # format is ">DesignID_0.AAACAAAGAA": DesignID_BarcodeID.Barcode
		designID = "".join(temp[:-1]) # to consider the case that designID has "_"
		#barcode = temp[-1].split(".")[1] # let's get barcode from sequence, not from header.

		# sequence
		line = genomeFile.readline()
		temp = line.strip()
		# temp[:16] # 5' adapter: ACT GGC CGC TTC ACTG
		# temp[-17:] # 3' adapter: AGA TCG GAA GAG CGT CG
		seq = temp[16:(-17-BARCODE_SIZE)]
		barcode = temp[(-17-BARCODE_SIZE):-17]

		oligoPool[barcode] = [designID, seq]

	else:
		ferr.write("Error in oligoPool: "+lineRaw+"\n")
		sys.exit()

	lineRaw= genomeFile.readline()

genomeFile.close()

# ==============================================
# 2. Read fastq and align sequence by barcode
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
	tileRead = seq_rc[:-BARCODE_SIZE]

	if barcodeRead in oligoPool :
		[designID, tile] = oligoPool[barcodeRead]
		matchNum = sum([1 for i in range(0, min(len(tile), len(tileRead))) if tile[-i-1] == tileRead[-i-1]])
		out = [qname, barcodeRead, designID, str(matchNum)]
	else :
		out = [qname, barcodeRead, "Unmapped", "-1"]

	fout.write("\t".join(out)+"\n")

ferr.write("The number of reads: " + str(lineNum/4) + "\n")

fin.close()
fout.close()
ferr.close()
