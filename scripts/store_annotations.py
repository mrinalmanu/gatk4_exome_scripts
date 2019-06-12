import sys
import os
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
# from Bio.Alphabet import generic_dna, generic_protein
# from Bio.Alphabet import IUPAC
# from Bio.SeqUtils import seq3

import pandas as pd
import gzip
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import OrderedDict
import math
from natsort import natsorted, ns
import argparse
import string
import random
from config_file import *
import subprocess


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
	return ''.join(random.choice(chars) for _ in range(size))


def check_path(path):                		## Check if file exists
			if os.path.exists(path) == True:
				return True
			else:
				return False
	
def read_annotated_vcf(filename):
	########### READ VCF  ################################################################################
	df_vcf = pd.read_csv(filename,
					sep='	',
					comment='#',
					header=None,
					names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR'])
	#######################################################################################################
	if len(df_vcf)==0:
		print "no variance"
		return
#	if len(df_vcf)> 300:
#		print "long vcf"
#		quit()	
	############################ DETERMINE TUMOR-NORMAL TYPE #################################################
	
	df_vcf["CHROM"]=df_vcf["CHROM"].astype(str)
	df_vcf["POSVAR"]=df_vcf["CHROM"].map(str)+"_"+df_vcf["POS"].map(str)+"_"+df_vcf["REF"]+"_"+df_vcf["ALT"]

	###########   SAMPLE INFO   ##########
	df_vcf['SOURCE_FILE']=filename
	
	####### VEP ANNOTATIONS #################### (always run VEP with --pick option)
	df_vcf['CONSEQ']=df_vcf['INFO'].str.split('|').str[1]
	df_vcf['IMPACT']=df_vcf['INFO'].str.split('|').str[2]
	df_vcf['SYMBOL']=df_vcf['INFO'].str.split('|').str[3]
	df_vcf['SIFT']=df_vcf['INFO'].str.split('|').str[27]
	df_vcf['PolyPhen']=df_vcf['INFO'].str.split('|').str[28]
	
	df_vcf['gnomAD_AFR|AMR|ASJ|EAS|FIN|']= df_vcf['INFO'].str.split('|').str[38]+"|"+ \
	df_vcf['INFO'].str.split('|').str[39]+"|"+ \
	df_vcf['INFO'].str.split('|').str[40]+"|"+ \
	df_vcf['INFO'].str.split('|').str[41]+"|"+ \
	df_vcf['INFO'].str.split('|').str[42]
	
	df_vcf['pubmed']=df_vcf['INFO'].str.split('|').str[49]
	df_vcf['rs_no']=df_vcf['INFO'].str.split('|').str[17]
	df_vcf['transcript_ID']=df_vcf['INFO'].str.split('|').str[6]
	df_vcf['codons']=df_vcf['INFO'].str.split('|').str[16]
	df_vcf = df_vcf[df_vcf['CHROM'].isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM','X','Y','M','MT','mt'])]
	df_vcf=df_vcf[['CHROM','POS','rs_no','REF','ALT','POSVAR','SYMBOL','codons','CONSEQ','IMPACT','SIFT','PolyPhen','gnomAD_AFR|AMR|ASJ|EAS|FIN|','pubmed','transcript_ID']]
	return df_vcf


def read_unannotated_file(filename):
	df_csv = pd.read_csv(filename, sep=None, engine='python')
	if "POSVAR" not in df_csv:
		print "'POSVAR' column not in file. Check format of file '"+filename+"'."
		quit()
	return df_csv 

def amino_change(row):
	chrom= row['CHROM']
	codons= row['codons']
	
	if codons in ["","-","Nan","NaN","NA","N/A"]:
		return ""
	if chrom in ['chrM','chrMT','M','MT']:
		translate_table = "Vertebrate Mitochondrial"
	else:
		translate_table = "Standard"
			
	ref_codon= codons.split("/")[0]
	alt_codon= codons.split("/")[1]

	translated_codon_ref = seq3(Seq(ref_codon).translate(table= translate_table, to_stop=True))	
	translated_codon_alt = seq3(Seq(alt_codon).translate(table= translate_table, to_stop=True))	
	return translated_codon_ref+"->"+translated_codon_alt

def create_dummy_vcf(df_not_in_db, dummy_ID):
	df_posvar = pd.DataFrame(df_not_in_db["POSVAR"])
	df_posvar["#CHROM"] = df_posvar["POSVAR"].str.split("_").str[0]
	df_posvar["POS"] = df_posvar["POSVAR"].str.split("_").str[1]
	df_posvar["ID"] = "dummy_ID"
	df_posvar["REF"] = df_posvar["POSVAR"].str.split("_").str[2]
	df_posvar["ALT"] = df_posvar["POSVAR"].str.split("_").str[3]
	df_posvar["QUAL"] = "99"
	df_posvar["FILTER"] = "PASS"
	df_posvar["INFO"] = "NoInfo"
	df_posvar["FORMAT"] = "NoFormat"
	df_posvar["SAMPLE_1"] = "0/1"

	df_posvar = df_posvar[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE_1']]
	df_posvar.to_csv(output_folder+"dummy_vcf_"+dummy_ID+".vcf", index=False, sep="\t")
	
def checkContainer(container_name):
		containers = subprocess.check_output(['docker','ps','-a']).decode(encoding="437")
		if container_name in containers.split():
			subprocess.call(["docker","rm",container_name])
		else:
			return


def run_in_docker(cmd,t_n, image, threads_needed, stdout=None, stderr=None):
		
		container_name = "db_vep"
		checkContainer(container_name)
		dcmd = ["docker", "run","--name",container_name,
						"-v", "{}:{}".format(output_folder, output_folder),
						"-v", "{}:{}".format(reference_folder, reference_folder),
						"-u","0",
						image]
		dcmd += cmd
	
		########### RUN COMMAND IN DOCKER #################
		os.system(" ".join(dcmd))
		checkContainer(container_name)
		
		###################### END OF RUN COMMANDS IN DOCKER CONTAINERS ###################


###### annotation in Variant Effect Predictor #####
def annotateVep(vcf_file):
	annotate_cmd = ["vep",
							"--af", 
							"--af_1kg", 
							"--af_esp",
							"--af_gnomad", 
							"--appris", 
							"--biotype", 
							"--check_existing", 
							"--pick", 
							"--distance 5000", 
							"--merged", 
							"--polyphen b", 
							"--pubmed", 
							"--regulatory", 
							"--sift b", 
							"--species homo_sapiens", 
							"--symbol", 
							"--tsl", 
							"--cache", 
							"--offline", 
							"--dir_cache", vep_reference_folder_hg19, 
							"--fork", max_nr_threads, 
							"--input_file",vcf_file, 
							"--force_overwrite", 
							"--fasta", reference_fasta_hg19, 
							"-o", vcf_file.replace(".vcf",".annotated.vcf"), 
							"--vcf"]
	run_in_docker(annotate_cmd,"asa", vep_docker_image, max_nr_threads)


def create_empty_db():
	empty_df = pd.DataFrame(columns=['CHROM','POS','rs_no','REF','ALT','POSVAR','SYMBOL','codons','CONSEQ','IMPACT','SIFT','PolyPhen','gnomAD_AFR|AMR|ASJ|EAS|FIN|','pubmed','transcript_ID'])
	return empty_df
	
def import_annotations_pck(pck_file, df_annotations, rewrite):
	df_pck = read_annotated_pck(pck_file)
	if rewrite==True:
		df_not_in_pck = df_annotations[~df_annotations["POSVAR"].isin(df_pck["POSVAR"])]
		df_annotations = pd.concat([df_not_in_pck,df_pck])
		print len(df_pck)," variants read from ",pck_file,
		print len(df_annotations)-len(df_not_in_pck)," variants in database were updated from imported data (--force_update was invoked)."

	else:	
		df_vcf_novel = df_vcf[~df_vcf["POSVAR"].isin(df_annotations["POSVAR"])]
		df_annotations = pd.concat([df_annotations,df_vcf_novel])
		print len(df_vcf)," variants read from ",vcf_file,
		print len(df_vcf_novel), "out of them are novel."
	return df_annotations

def import_annotations_vcf(vcf_file, df_annotations, rewrite):
	df_vcf = read_annotated_vcf(vcf_file)
	if df_vcf == None:
		return create_empty_db()
	if rewrite==True:
		df_not_in_vcf = df_annotations[~df_annotations["POSVAR"].isin(df_vcf["POSVAR"])]
		df_annotations = pd.concat([df_not_in_vcf,df_vcf])
		print len(df_vcf)," variants read from ",vcf_file,
		print len(df_annotations)-len(df_not_in_vcf)," variants in database were updated from imported data (--force_update was invoked)."
	else:	
		df_vcf_novel = df_vcf[~df_vcf["POSVAR"].isin(df_annotations["POSVAR"])]
		df_annotations = pd.concat([df_annotations,df_vcf_novel])
		print len(df_vcf)," variants read from ",vcf_file,
		print len(df_vcf_novel), "out of them are novel."
	return df_annotations
	
def main_func():
	########### PARSE COMMAND LINE ARGUMENTS #########################
	# Instantiate the parser
	parser = argparse.ArgumentParser(description='***Update annotations database and annotate csv/vcf files.***')
	
	# Optional argument
	parser.add_argument('--annot_db', type=str, help='choose [DEFAULT] to use DB from config file OR enter full path to pickle file with stored annotations.')
	parser.add_argument('--batch_import', type=str, help='List of annotated vcf files to import annotations from.')
	parser.add_argument('--to_annotate', type=str, help='List of csv or vcf files to be annotated.')
	parser.add_argument('--reannotate_database',type=str , help='If set to TRUE, runs annotation on all variants in database.')
	parser.add_argument('--import_annotations', type=str, help='Path to file. Import annotations from pickle or vcf.')
	parser.add_argument('--batch_name', type=str, help='Any string. ')
	parser.add_argument('--run_annotator', type=str, help='Run annotator on variants not present in database.')
	parser.add_argument('--force_update', type=str, help='If set to TRUE, imported variants overwrite variants already present in database.')                 
	parser.add_argument('--save_updated', type=str, help='If TRUE, updated annotations will be saved to loaded database. If [PATH], updated annotations will be saved to [PATH].')
                                              
	args = parser.parse_args()
	############ END OF PARSER #######################################
	
	if args.force_update in ["True","true","TRUE","yes","y","YES"]:
		rewrite=True
	else:
		rewrite=False
	
	###### read pickled database #########
	if args.annot_db !=None:
		if args.annot_db in ["DEFAULT","Default","default"]:
			print "Loading default database: ",default_annotations_db
			df_annotations = pd.read_pickle(default_annotations_db)
			print str(len(df_annotations))+" variants found in "+default_annotations_db
		elif check_path(args.annot_db)==True:
			print "Loading database: ",args.annot_db
			df_annotations = pd.read_pickle(args.annot_db)
			print str(len(df_annotations))+" variants found in "+args.annot_db
		else:
			print "Annotations file not found, exiting."
			quit()
	else:
		print "Database file not specified, creating empty DB."
		df_annotations = create_empty_db()
			
	##### import annotations from pickle or vcf #############
	if args.import_annotations != None:
		if args.import_annotations.split(".")[-1]=="vcf":
			print "Importing annotations from vcf: ",args.import_annotations 
			df_annotations = import_annotations_vcf(args.import_annotations, df_annotations, rewrite)
		if args.import_annotations.split(".")[-1]=="pck":
			df_annotations = import_annotations_pickle(args.import_annotations, df_annotations, rewrite)
					
	##### batch import ##################################################
	if args.batch_import != None:
		if check_path(args.batch_import)==False:
			print "ERROR: batch import list not found."
			quit()
		
		filelist = open(args.batch_import).read().splitlines()
		for onefile in filelist:
			if check_path(onefile)==False:
				print "File "+onefile+" not found."
				continue
			if onefile.split(".")[-1]=="vcf":
				print "Importing annotations from vcf: ",onefile
				df_annotations = import_annotations_vcf(onefile, df_annotations, rewrite)
			if onefile.split(".")[-1]=="pck":
				print "Importing annotations from pickle: ",onefile
				df_annotations = import_annotations_pickle(onefile, df_annotations, rewrite)
						
	###### re-annotate database #########
	if args.reannotate_database in ["True","true","TRUE","yes","y","YES"]:
		reannotate = raw_input('You choose to annotate whole database. This may take long time. Continue? [yes/no]')
		if reannotate in ['yes','y','YES','Y']:
			dummy_ID = id_generator()
			create_dummy_vcf(df_annotations, dummy_ID)
			annotateVep(output_folder+"dummy_vcf_"+dummy_ID+".vcf")
			import_annotations_vcf(output_folder+"dummy_vcf_"+dummy_ID+".annotated.vcf", df_annotations, rewrite) 
		
	##### read unannotated file ###########		
	if args.to_annotate !=None:
		unannotated_file = args.to_annotate
		if check_path(unannotated_file)==False:
			print "File "+args.to_annotate+" not found."
			quit()
		df_csv = read_unannotated_file(unannotated_file)	
		df_not_in_db = df_csv[~df_csv["POSVAR"].isin(df_annotations["POSVAR"])]
		## run annotator ##
		if args.run_annotator in ["True","true","TRUE","yes","y","YES"]:
			dummy_ID = id_generator()
			create_dummy_vcf(df_not_in_db, dummy_ID)
			annotateVep(output_folder+"dummy_vcf_"+dummy_ID+".vcf")
			import_annotations_vcf(output_folder+"dummy_vcf_"+dummy_ID+".annotated.vcf", df_annotations, rewrite)
		else:
			if len(df_not_in_db) !=0:
				print "WARNING: "+str(len(df_not_in_db))+" variants in "+unannotated_file+" not found in database and cannot be annotated."
				print "Run --run_annotator if you wish to annotate all variants."
		
		merge = pd.merge(df_csv, df_annotations,on='POSVAR', how='left')	
	
		#### MANAGE SAVING OF ANNOTATIONS TO PICKLE ##############################################################
		if (args.annot_db !=None) and (args.save_updated ==None):
			rewrite_db = raw_input('Do you want to save updated annotations into '+args.annot_db+'? [yes/no]')
			if rewrite_db in ['yes','y','YES','Y']:
				if args.annot_db in ["DEFAULT","Default","default"]:
					print "Saving to default database: ",default_annotations_db
					df_annotations.to_pickle(default_annotations_db)
				else:
					print "Saving to database: ",args.annot_db
					df_annotations.to_pickle(args.annot_db)
					
		elif args.save_updated != None:
			if args.save_updated in ["True","true","TRUE","yes","y","YES"]:			
				if args.annot_db in ["DEFAULT","Default","default"]:
					print "Saving to default database: ",default_annotations_db
					df_annotations.to_pickle(default_annotations_db)
				else:
					print "Saving to database: ",args.annot_db
					df_annotations.to_pickle(args.annot_db)
			else:
				print "Saving to database: ",args.save_updated
				df_annotations.to_pickle(args.save_updated)
					
		else:
			save_db = raw_input('New annotations were generated. Do you want to save them? [yes/no]')
			if save_db in ['yes','y','YES','Y']:
				db_path = raw_input('Enter file name with full path (extension must be .pck) OR write DEFAULT to save to default database.')
				if db_path in ["DEFAULT","Default","default"]:
					print "Saving to default database: ",default_annotations_db
					df_annotations.to_pickle(default_annotations_db)
				else:
					print "Saving to database: ",db_path
					df_annotations.to_pickle(db_path)
		#### END OF: MANAGE SAVING OF ANNOTATIONS TO PICKLE ##############################################################
	
	
	
		
#	df_annotations['AMINO_CHANGE']=df_master.apply(amino_change, axis=1)
#	df_annotations = df_master.sort_values(by=['POSVAR'])
	
	#df_annotations.to_csv("test.csv", index=False)
	df_annotations.to_pickle("ann_db.pck")

if __name__== "__main__":
	main_func()
	





