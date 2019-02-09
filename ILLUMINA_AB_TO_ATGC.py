import datetime
dttime=datetime.datetime.now().strftime ("%Y%m%d")

'''
A python 2 script to parse illumina matrix genotypes convert AB calls into ATGC 
Author : Aditya Ambati ambati@stanford.edu

'''

def NonAmbiguous(SNP, IllStrand):
	''' This fucntion will process nonambigous snps in AB format to ATGC format again refer to https://www.illumina.com/documents/products/technotes/technote_topbot.pdf'''
	if SNP in ['[A/G]', '[G/A]']:
		assert IllStrand == 'TOP'
		A, B = 'A', 'G'
	elif SNP in ['[A/C]', '[C/A]']:
		assert IllStrand == 'TOP'
		A, B = 'A', 'C'
	elif SNP in ['[T/G]', '[G/T]']:
		assert IllStrand == 'BOT'
		A, B = 'T', 'G'
	elif SNP in ['[T/C]', '[T/A]']:
		assert IllStrand == 'BOT'
		A, B = 'T', 'C'
	return A, B


def ConvertAB(SNP, IllStrand):
	''' This fucntion will process ambigous snps in AB format to ATGC format again refer to https://www.illumina.com/documents/products/technotes/technote_topbot.pdf'''
	ATs = ['[A/T]', '[T/A]']
	GCs= ['[G/C]', '[C/G]']
	if SNP in ATs and IllStrand == 'TOP': ## this will get us if it is A/T snp and next is to get TOP /BOT
		A, B = 'A','T'
	elif SNP in ATs and IllStrand == 'BOT':
		B, A = 'A','T'
	elif SNP in GCs and IllStrand == 'TOP':
		A, B = 'C', 'G'
	elif SNP in GCs and IllStrand == 'BOT':
		A, B = 'G', 'C'
	else:
		A, B = NonAmbiguous(SNP, IllStrand)
	return [A, B]


def ParseReference(RefFile):
	'''this function takes in a illumina reference file and parses the TOP BOT, SNP Id and Refstrna dinformation and then makes a dict with Chr+':'+pos+':'+rsid as key
	and the genotype alleles and strand info as the values, further this converts AB calls to ATGC calls as mentioned in https://www.illumina.com/documents/products/technotes/technote_topbot.pdf
	'''
	from collections import defaultdict
	RsDic = defaultdict(lambda:defaultdict(str))
	SNPCount = 0
	for n, line in enumerate(open(RefFile, 'r')):
		if n > 0:
			ParseLine= line.strip().split(',')
			Chr, Pos, Rsid, IllStrand, SNP, RefStrand= ParseLine[9], ParseLine[10], ParseLine[1], ParseLine[2], ParseLine[3], ParseLine[20]#, ParseLine[17]
			RsKey = Chr+':'+Pos+':'+Rsid
			if SNP not in ['[D/I]', '[I/D]']: ## Drop deletions we will let the imputation handle these 
				GetAB=ConvertAB(SNP = SNP, IllStrand=IllStrand)
				SNPCount += 1
				RsDic[RsKey]['IllStrand'] = IllStrand
				RsDic[RsKey]['RefStrand'] = RefStrand
				RsDic[RsKey]['SNP'] = SNP
				RsDic[RsKey]['A'] = GetAB[0]
				RsDic[RsKey]['B'] = GetAB[1]
	return RsDic



def ProcessGeno(FileIn, FlipFile, TfamFile, TpedFile, RefFile):
	'''This is the main engine of the code, processing the matrix of genotypes exported out of bead studio
	this does the following
		1. Reads the reference illumina file and makes a dictionary item containing the snp, strand information so on
		2. operates on the first line parsing the sample IDs and writing out a Tfam plink file refer to plink documentation
		3. operates on rest of the lines picking out every 3 item in single line as a genotype call, further the rsid in the line is queried against the refdic and AB genotypes are convertd into ATGC
		4. Further if the function encounters a rsid that has been assigned in  negative strand ('-'), then this will write out the rsids so these can be flipped using plink
	'''
	SNPCount =0
	Discard = 0
	print 'PROCESSING THE REFERENCE FILE {}'.format(RefFile)
	RefRsDic = ParseReference(RefFile=RefFile)
	with open(FileIn) as Report:
		for n, line in enumerate(Report):
			ParseLine = line.strip().split(',')
			if n == 0:
				for sample in range(3, len(ParseLine), 6):
					TfamFile.write(' '.join([ParseLine[sample].split('.')[0], ParseLine[sample].split('.')[0], str(0), str(0), str(0), str(-9)])+'\n')
				TfamFile.close()
			else:
				Chr = ParseLine[1]
				Rsid = ParseLine[0]#.replace('exm-', '')
				PosCoord = '0'
				Pos = ParseLine[2]
				RsKey = Chr+':'+Pos+':'+Rsid
				if RsKey in RefRsDic:
					SNPCount += 1
					TpedFile.write(' '.join([Chr, Rsid, PosCoord, Pos])+' ')
					GetCurDic = RefRsDic.get(RsKey)
					if GetCurDic.get('RefStrand') == '-':
						FlipFile.write(Rsid+'\n')

					for i in range(3, len(ParseLine), 6):
						if ParseLine[i] == 'NC':
							GenoCall = '00'
						else:
							GenoCall=GetCurDic.get(ParseLine[i][0])+GetCurDic.get(ParseLine[i][1])
						TpedFile.write(GenoCall[0]+' '+GenoCall[1]+' ')
					TpedFile.write('\n')
				else:
					Discard += 1
		print 'PROCESSED {} GENOTYPES FROM  {} AND DISCARDED {} GENOTYPES'.format(SNPCount, FileIn, Discard)
		TpedFile.close()
		FlipFile.close()


def main():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-RefFile', help='illumina file with reference snp information', required=True)
	parser.add_argument('-FileIn', help='illumina bead studio exported genotype matrix file containing AB genotypes', required=True)
	parser.add_argument('-OutFile', help='A prefix for the outfile', required=True)
	args = parser.parse_args()
	RefFile = args.RefFile
	FileIn = args.FileIn
	TfamFile, TpedFile, FlipFile = [open(args.OutFile+i, 'w') for i in ['.tfam', '.tped', '.txt']]
	#print args
	ProcessGeno(FileIn = FileIn, FlipFile=FlipFile, TpedFile=TpedFile, TfamFile=TfamFile, RefFile=RefFile)

if __name__ == '__main__':main()


# # ##### Argument Parser #####
# FileIn = 'GSE83709_Matrix_processed_Q35.csv'
# FlipFile = open('FLIPFILE.txt', 'w')
# TfamFile = open('test.tfam', 'w')
# TpedFile = open('test.tped', 'w')
# RefFile='HumanOmniExpressExome-8-v1-1-C.csv'







# AnnotTest = open('HumanOmniExpressExome-8-v1-1-C.csv', 'r')
# RsDic = {}
# SNPCount = 0
# for n, line in enumerate(AnnotTest):
# 	if n > 0:
# 		ParseLine= line.strip().split(',')
# 		#Rsid, IllStrand, SNP, RefStrand, Genomic= ParseLine[1], ParseLine[2], ParseLine[3], ParseLine[20], ParseLine[17]
# 		Rsid, IllStrand, SNP, RefStrand= ParseLine[1], ParseLine[2], ParseLine[3], ParseLine[20]#, ParseLine[17]
# 		if SNP not in ['[D/I]', '[I/D]']:
# 			print ConvertAB(SNP = SNP, IllStrand=IllStrand), Rsid, IllStrand, SNP, RefStrand
# 			SNPCount += 1
# 		else:
# 			pass
# 		if SNP in []
# 		SNPA, SNPB = SNP[1], SNP[3]
# 		if IllStrand == 'TOP':
# 			if SNPA == 'A':
# 				A = SNPA

