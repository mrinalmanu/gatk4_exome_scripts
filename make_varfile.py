import pandas as pd
import sys
import os

filename=sys.argv[1]

def read_vcf(filename):
	
	########### READ VCF  ################################################################################
	if filename.split(".")[-1]=="gz":
		openfile=gzip.open(filename)
	else:
		openfile=filename
			
	df_vcf = pd.read_csv(openfile,
					sep='	',
					comment='#',
					header=None,
					names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR'])
	#######################################################################################################
	
	return df_vcf

if os.path.exists(filename) == True:
	df_vcf = read_vcf(filename)
	varfile=df_vcf[['CHROM','POS','REF','REF']]
	varfile.to_csv(filename.replace(".vcf",".varfile"), index= False, sep='\t', header=None)
else:
	print >> sys.stderr, "MAKE VARFILE: vcf file not found."
