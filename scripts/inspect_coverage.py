import sys
import os

import pandas as pd
import gzip
import numpy as np
from collections import OrderedDict
import math
from natsort import natsorted, ns
from multiprocessing import Process
from natsort import natsorted
from config_file import *
import argparse
import time
import datetime

from store_annotations import create_empty_db
from store_annotations import id_generator
from store_annotations import create_dummy_vcf
from store_annotations import annotateVep
from store_annotations import import_annotations_vcf

def get_time():
	time_stamp=time.time()
	value = datetime.datetime.fromtimestamp(time_stamp)
	return value.strftime('%Y-%m-%d %H:%M:%S ')

def check_path(path):                		## Check if file exists
		if os.path.exists(path) == True:
			return True
		else:
			return False

def create_samples_paths(samples):
	
	samples_paths={}
	for sample in samples:
		sample_name = sample[0]
		pipeline = sample[1]
		
		final_vcf_naming = vcf_naming_dict[pipeline]
		recalibrated_bam_naming = bam_naming_dict[pipeline]
		
		if check_path(output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+final_vcf_naming)==True:
			vcf_path= output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+final_vcf_naming
		elif check_path(output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+final_vcf_naming+".gz")==True:
			vcf_path= output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+final_vcf_naming+".gz"
		else:
			vcf_path="N/A"
			print get_time()+ sample_name+" Error: vcf file nof found in output folder"
			print get_time()+ output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+final_vcf_naming
			continue
		
		if check_path(output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+recalibrated_bam_naming+"_tumor.bam")==True:
			tumor_bam_path= output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+recalibrated_bam_naming+"_tumor.bam"
		else:
			tumor_bam_path="N/A"
			print get_time()+ sample_name+" Error: recalibrated tumor bam file not found in output folder"
			continue
			
			
		if check_path(output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+recalibrated_bam_naming+"_normal.bam")==True:
			normal_bam_path= output_folder+"/"+sample_name+"/output/"+pipeline+"/"+sample_name+recalibrated_bam_naming+"_normal.bam"
		else:
			normal_bam_path="N/A"
			print get_time()+ sample_name+" Error: recalibrated normal bam file not found in output folder"
			continue
			
		samples_paths.update({sample_name:[vcf_path, tumor_bam_path, normal_bam_path]})
	return samples_paths



def read_vcf(vcf_file, germline_somatic, sample_name, tissue):
	
	
	########### READ VCF  ################################################################################
	df_vcf = pd.read_csv(vcf_file,
					sep='	',
					comment='#',
					header=None,
					names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR'],engine='python')
	#######################################################################################################
	if len(df_vcf)==0:
		print get_time()+ "no variance"
		return
	############################ DETERMINE TUMOR-NORMAL TYPE #################################################
	
	#print get_time()+ germline_somatic
	if germline_somatic=="Normal":
		df_vcf["GENOTYPE"]=df_vcf["NORMAL"]
	if germline_somatic=="Tumor":
		df_vcf["GENOTYPE"]=df_vcf["TUMOR"]	
	
	df_vcf["CHROM"]=df_vcf["CHROM"].astype(str)
	df_vcf["POSVAR"]=df_vcf["CHROM"].map(str)+"_"+df_vcf["POS"].map(str)+"_"+df_vcf["REF"]+"_"+df_vcf["ALT"]
	df_vcf["CHROMPOS"]=df_vcf["CHROM"].map(str)+"_"+df_vcf["POS"].map(str)
	
	### check vcf format ###
	mutect_format="GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB"
	varscan_format="GT:GQ:DP:RD:AD:FREQ:DP4"
	
	if set(["F1R2","F2R1","MBQ","MFRL","AD"]).issubset(df_vcf.iloc[0]['FORMAT'].split(":")):
		caller_format = "mutect"
	elif set(["FREQ","DP4","AD","RD"]).issubset(df_vcf.iloc[0]['FORMAT'].split(":")):
		caller_format = "varscan"
	else:
		print get_time()+"ERROR: vcf from unsupported variant caller (supported: varscan2 T-N, mutect2 T-N)"
		return 	
	########################
	
	### read allelic depths #####################################
	if caller_format =="mutect":
		print get_time()+"Detected format: mutect2"
		AO_index=df_vcf.iloc[0]['FORMAT'].split(':').index("AD")
		df_vcf['ALT_SUPPORT']=df_vcf['GENOTYPE'].str.split(":").str[AO_index].str.split(',').str[1].map(int)
		df_vcf['REF_SUPPORT']=df_vcf['GENOTYPE'].str.split(":").str[AO_index].str.split(',').str[0].map(int)
		df_vcf['ALT_SUPPORT']= df_vcf['ALT_SUPPORT'].astype(int)
		df_vcf['REF_SUPPORT']= df_vcf['REF_SUPPORT'].astype(int)
	#	df_vcf['READ_DEPTH']= df_vcf['INFO'].str.split(";").str[0].str.split("=").str[1].apply(pd.to_numeric)
		df_vcf['DEPTH']=df_vcf['ALT_SUPPORT']+df_vcf['REF_SUPPORT']
	if caller_format =="varscan":
		print get_time()+"Detected format: mutect2"
		AD_index=df_vcf.iloc[0]['FORMAT'].split(':').index("AD")
		RD_index=df_vcf.iloc[0]['FORMAT'].split(':').index("RD")
		df_vcf['ALT_SUPPORT']=df_vcf['GENOTYPE'].str.split(":").str[AD_index].map(int)
		df_vcf['REF_SUPPORT']=df_vcf['GENOTYPE'].str.split(":").str[RD_index].map(int)
		df_vcf['ALT_SUPPORT']= df_vcf['ALT_SUPPORT'].astype(int)
		df_vcf['REF_SUPPORT']= df_vcf['REF_SUPPORT'].astype(int)
	#	df_vcf['READ_DEPTH']= df_vcf['INFO'].str.split(";").str[0].str.split("=").str[1].apply(pd.to_numeric)
		df_vcf['DEPTH']=df_vcf['ALT_SUPPORT']+df_vcf['REF_SUPPORT']
	##############################################################	
		
	df_vcf['PATIENT']=sample_name

	tumor_normal= germline_somatic
	patient= sample_name#filename.split('_')[1]
	df_vcf['TUMOR_NORMAL']= tumor_normal
	## filter out alternative contigs ##
	df_vcf = df_vcf[df_vcf['CHROM'].isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM','X','Y','M','MT','mt'])]
	df_vcf=df_vcf[['PATIENT','TUMOR_NORMAL','FILTER','CHROM','CHROMPOS','POS','REF','ALT','POSVAR','DEPTH','ALT_SUPPORT','REF_SUPPORT']]

	return df_vcf

def bam_coverage(patient, bed_file, tumor_bam, normal_bam):
	
	df_bed = pd.read_csv(bed_file,
						sep='	',
						comment='#',
						header=None,
						names=['CHROM','POS'],engine='python')
	df_bed['POS'] = df_bed['POS'].map(str)
	df_bed['CHROMPOS']=df_bed['CHROM'].map(str)+"_"+df_bed['POS']
	df_bed = df_bed[['CHROMPOS']]
	
	df_tumor_cov = pd.read_csv(tumor_bam+".cov",
						sep='	',
						comment='#',
						header=None,
						names=['CHROM','POS','TUMOR_COV'],engine='python')
	df_tumor_cov['TUMOR_COV'] = df_tumor_cov['TUMOR_COV'].astype(int)				
	df_tumor_cov['POS'] = df_tumor_cov['POS'].map(str)
	df_tumor_cov['CHROMPOS']=df_tumor_cov['CHROM']+"_"+df_tumor_cov['POS']
	df_tumor_cov = pd.merge(df_bed, df_tumor_cov, how="left",  on= ["CHROMPOS"])
	df_tumor_cov = df_tumor_cov.set_index('CHROMPOS')
	df_tumor_cov = df_tumor_cov.fillna(0)

	
	df_normal_cov = pd.read_csv(normal_bam+".cov",
						sep='	',
						comment='#',
						header=None,
						names=['CHROM','POS','NORMAL_COV'],engine='python')
	df_normal_cov['NORMAL_COV'] = df_normal_cov['NORMAL_COV'].astype(int)				
	df_normal_cov['POS'] = df_normal_cov['POS'].map(str)
	df_normal_cov['CHROMPOS']=df_normal_cov['CHROM']+"_"+df_normal_cov['POS']
	df_normal_cov = pd.merge(df_bed, df_normal_cov, how="left",  on= ["CHROMPOS"])
	df_normal_cov = df_normal_cov.set_index('CHROMPOS')
	df_normal_cov = df_normal_cov.fillna(0)
	
	df_tumor_normal = pd.concat([df_tumor_cov, df_normal_cov], axis=1,join_axes=[df_tumor_cov.index])
	df_tumor_normal['CHROMPOS']=df_tumor_normal.index
	df_tumor_normal = df_tumor_normal[['CHROMPOS','TUMOR_COV','NORMAL_COV']]
	df_tumor_normal = df_tumor_normal.fillna(0)
	df_tumor_normal['PATIENT'] = patient
	
	return df_tumor_normal

def pass_fail(filter_field):
	if filter_field =='PASS':
		return 'CALL:PASS'
	else:
		return 'CALL:FAIL'	

def submit_coverage(bamfile, bed_file, base_quality_threshold, mapping_quality_threshold):
	os.system("samtools depth "+bamfile+" -q "+base_quality_threshold+" -Q "+mapping_quality_threshold+" -b "+bed_file+" > "+bamfile+".cov")
		
def main_func(bed_file, base_quality_threshold, mapping_quality_threshold):
	########### PARSE COMMAND LINE ARGUMENTS #########################
	# Instantiate the parser
	parser = argparse.ArgumentParser(description='***Create lookup table from vcf files.***')
	
	# Optional argument
	parser.add_argument('--sample_list', type=str, help='Txt file w/list of samples and pipelines.')
	parser.add_argument('--output', type=str, help='Name of output table.')
	parser.add_argument('--add_annotations', type=str, help='If TRUE, calls will be annotated')
	parser.add_argument('--annot_db', type=str, help='choose [DEFAULT] to use DB from config file OR enter full path to pickle file with stored annotations.')
	parser.add_argument('--run_annotator', type=str, help='Run annotator on variants not present in database.')
	parser.add_argument('--save_updated', type=str, help='If TRUE, updated annotations will be saved to loaded database. If [PATH], updated annotations will be saved to [PATH].')
	parser.add_argument('--force_update', type=str, help='If set to TRUE, imported variants overwrite variants already present in database.')                 

	
	args = parser.parse_args()
	
	if args.force_update in ["True","true","TRUE","yes","y","YES"]:
		rewrite=True
	else:
		rewrite=False
	
	bed_file = output_folder+args.output+'.bed'
	sample_list=open(args.sample_list).read().splitlines()
	samples=[]
	for line in sample_list:
		samples.append(line.split())	
	samples_paths = create_samples_paths(samples)
	
	## CREATE DATAFRAME FROM VCF FILES ####################
	first=True
	print get_time()+ "Parsing vcf files..."
	for sample, paths in samples_paths.iteritems():
		sample_name=sample
		vcffile=paths[0]
		tumor_bam = paths[1]
		normal_bam = paths[2]
		
		
		tissue="tiss"
		print get_time()+ vcffile, "Normal vcf..."
		df_vcf = read_vcf(vcffile,'Normal',sample_name, tissue)
		if first == True:
			df_master=df_vcf
			first=False

		else:
			df_master=pd.concat([df_master,df_vcf])

		print get_time()+ vcffile, "Tumor vcf..."
		df_vcf = read_vcf(vcffile,'Tumor', sample_name, tissue)
		df_master=pd.concat([df_master,df_vcf])
	df_master=df_master[df_master["TUMOR_NORMAL"]== 'Tumor']
	df_master = df_master.sort_values(by=['POSVAR'])
	df_master=df_master[['PATIENT','TUMOR_NORMAL','FILTER','CHROM','CHROMPOS','POS','REF','ALT','POSVAR','DEPTH','ALT_SUPPORT','REF_SUPPORT']]
	df_master['PASS_FAIL']=df_master['FILTER'].apply(pass_fail)
	
	# ### CREATE BED FILE ##############################
	# print get_time()+ "Creating bed file with union of calls..."
	# all_chrompos = df_master["CHROMPOS"].tolist()
	# all_chrompos = natsorted(list(set(all_chrompos)))
	# with open(bed_file, 'w') as the_file:
		# for chrompos in all_chrompos:
			# chrompos = chrompos.split('_')
			# the_file.write(chrompos[0]+"\t"+chrompos[1]+'\n') 
	
	# #### CREATE COVERAGE FILES FOR TUMOR AND NORMAL
	# print get_time()+ "Generating coverage files..."
	# proc=[]
	# for sample, paths in samples_paths.iteritems():
		
		# tumor_bam = paths[1]
		# normal_bam = paths[2]
		
		# p = Process(target = submit_coverage, args=[tumor_bam, bed_file, base_quality_threshold, mapping_quality_threshold])
		# p.start()
		# proc.append(p)
		
		# p = Process(target = submit_coverage, args=[normal_bam, bed_file, base_quality_threshold, mapping_quality_threshold])
		# p.start()
		# proc.append(p)

	# for p in proc:
		# p.join()
	########################################################
	
	########## READ BAM COVERAGE DATA #######################
	print get_time()+ "Parsing coverage data..."
	first=True
	for sample, paths in samples_paths.iteritems():
		bed_file = output_folder+args.output+".bed"
		tumor_bam = paths[1]
		normal_bam = paths[2]
		
		df_sample = bam_coverage(sample, bed_file, tumor_bam, normal_bam)
		if first == True:
			df_tumor_normal = df_sample
			first=False
		else:
			df_tumor_normal = pd.concat([df_sample, df_tumor_normal])
	df_tumor_normal = df_tumor_normal[['CHROMPOS','PATIENT','TUMOR_COV','NORMAL_COV']]
	#df_tumor_normal.to_csv("tn.csv")
	################################################################
	
	### MERGE VCF CALLS AND COVERAGE DATA ##########################
	print get_time()+ "Merging calls with coverage information..."
	all_posvar= natsorted(df_master['POSVAR'].unique())
	df_posvar= pd.DataFrame({'POSVAR':all_posvar})
	df_posvar['CHROM'] = df_posvar['POSVAR'].str.split("_").str[0].astype(str)
	df_posvar['POS'] = df_posvar['POSVAR'].str.split("_").str[1].astype(str)
	df_posvar["CHROMPOS"]= df_posvar["CHROM"]+"_"+df_posvar["POS"]
	df_posvar=df_posvar[["POSVAR","CHROMPOS"]]

	df_posvar = pd.merge(df_posvar, df_tumor_normal,  on= ["CHROMPOS"])
	df_master = pd.merge(df_posvar, df_master, on=["POSVAR",'PATIENT'], how='left')
	df_master["PASS_FAIL"] = df_master["PASS_FAIL"].fillna("NO_CALL")
	df_master["REF_SUPPORT"] = df_master["REF_SUPPORT"].fillna(0)
	df_master["ALT_SUPPORT"] = df_master["ALT_SUPPORT"].fillna(0)
	df_master["DEPTH"] = df_master["DEPTH"].fillna(0)
	df_master['CALL_STATUS']=df_master['PASS_FAIL'].map(str)+";ALT:"+df_master['ALT_SUPPORT'].map(str)+",DEP:"+df_master['DEPTH'].map(str)#"#+"COV|T:"+df_master['TUMOR_COV'].map(str)+",N:"+df_master['NORMAL_COV'].map(str)
	df_master['TN_COV']="T:"+df_master['TUMOR_COV'].map(str)+",N:"+df_master['NORMAL_COV'].map(str)
	
	df_pivot_detail = df_master[['POSVAR','PATIENT','TUMOR_COV','NORMAL_COV','REF_SUPPORT','ALT_SUPPORT','DEPTH','PASS_FAIL']]
	df_pivot_detail= df_pivot_detail[df_pivot_detail[['POSVAR','PATIENT']].duplicated() == False]
	df_pivot_detail = df_pivot_detail.pivot(index='POSVAR',columns='PATIENT')
	df_pivot_detail.to_csv(args.output+"_detailed.csv", index=True) 
	
	
	df_master=df_master[['POSVAR','PATIENT','TN_COV','CALL_STATUS']]
	df_master['CALL_STATUS'] = df_master['CALL_STATUS'].fillna("NO_CALL")
	df_master['SAMPLE'] = df_master["CALL_STATUS"]+"|COV(bqt:"+base_quality_threshold+",mqt:"+mapping_quality_threshold+"):"+df_master["TN_COV"]
#	df_master['SAMPLE'] = df_master['SAMPLE'].str.replace(".0","")
	df_master['SAMPLE'] = df_master['SAMPLE'].str.replace(";ALT:0.0,DEP:0.0","")
	
	df_pivot= df_master[['POSVAR','PATIENT','SAMPLE']]
	df_pivot= df_pivot[df_pivot[['POSVAR','PATIENT']].duplicated() == False]
	df_pivot= df_pivot.pivot(index='POSVAR',columns='PATIENT')
	df_pivot.to_csv(args.output+"_summary.csv", index=True)
	
	################################################################################
	## ADD ANNOTATIONS #############################################################
	################################################################################
	
	
	if args.add_annotations in ["True","true","TRUE","yes","y","YES"]:
		print get_time()+ "Adding annotation information..."
		###### read pickled database #########
		if args.annot_db !=None:
			if args.annot_db in ["DEFAULT","Default","default"]:
				print get_time()+ "Loading default database: ",default_annotations_db
				df_annotations = pd.read_pickle(default_annotations_db)
				print get_time()+ str(len(df_annotations))+" variants found in "+default_annotations_db
			elif check_path(args.annot_db)==True:
				print get_time()+ "Loading database: ",args.annot_db
				df_annotations = pd.read_pickle(args.annot_db)
				print get_time()+ str(len(df_annotations))+" variants found in "+args.annot_db
			else:
				print get_time()+ "Annotations file not found, exiting."
				quit()
		else:
			print get_time()+ "Database file not specified, creating empty DB."
			df_annotations = create_empty_db()
		
		### CHECK FOR VARIANTS NOT IN DB
		df_pivot['POSVAR']= df_pivot.index
		df_not_in_db = df_pivot[~df_pivot["POSVAR"].isin(df_annotations["POSVAR"])]
		## run annotator ##
		if args.run_annotator in ["True","true","TRUE","yes","y","YES"]:
			print get_time()+ "Running annotation engine..."
			dummy_ID = id_generator()
			create_dummy_vcf(df_not_in_db, dummy_ID)
			annotateVep(output_folder+"dummy_vcf_"+dummy_ID+".vcf")
			import_annotations_vcf(output_folder+"dummy_vcf_"+dummy_ID+".annotated.vcf", df_annotations, rewrite)
		else:
			if len(df_not_in_db) !=0:
				print get_time()+ "WARNING: "+str(len(df_not_in_db))+" variants were not found in database and cannot be annotated."
				print get_time()+ "Run [--run_annotator TRUE] if you wish to annotate all variants."
		## MERGE LOOKUP TABLE W/ANNOTATIONS 
		df_pivot_annotated = pd.merge(df_annotations, df_pivot, on=["POSVAR"], how="right")
		df_pivot_annotated.to_pickle(args.output+"_annotated.pck")
		df_pivot_annotated.to_csv(args.output+"_annotated.csv", index=False)
		
		######## MANAGE SAVING OF ANNOTATIONS TO PICKLE ##############################################################
		if (args.annot_db !=None) and (args.save_updated ==None):
			rewrite_db = raw_input('Do you want to save updated annotations into '+args.annot_db+'? [yes/no]')
			if rewrite_db in ['yes','y','YES','Y']:
				if args.annot_db in ["DEFAULT","Default","default"]:
					print get_time()+ "Saving to default database: ",default_annotations_db
					df_annotations.to_pickle(default_annotations_db)
				else:
					print get_time()+ "Saving to database: ",args.annot_db
					df_annotations.to_pickle(args.annot_db)
					
		elif args.save_updated != None:
			if args.save_updated in ["True","true","TRUE","yes","y","YES"]:			
				if args.annot_db in ["DEFAULT","Default","default"]:
					print get_time()+ "Saving to default database: ",default_annotations_db
					df_annotations.to_pickle(default_annotations_db)
				else:
					print get_time()+ "Saving to database: ",args.annot_db
					df_annotations.to_pickle(args.annot_db)
			else:
				print get_time()+ "Saving to database: ",args.save_updated
				df_annotations.to_pickle(args.save_updated)
		else:
			save_db = raw_input('New annotations were generated. Do you want to save them? [yes/no]')
			if save_db in ['yes','y','YES','Y']:
				db_path = raw_input('Enter file name with full path (extension must be .pck) OR write DEFAULT to save to default database.')
				if db_path in ["DEFAULT","Default","default"]:
					print get_time()+ "Saving to default database: ",default_annotations_db
					df_annotations.to_pickle(default_annotations_db)
				else:
					print get_time()+ "Saving to database: ",db_path
					df_annotations.to_pickle(db_path)
	#################################################################################################################################
	#################################################################################################################################
						



path=""

germline_frq_dict={}
germline_depth_dict={}

medium_high_terms=["protein_altering_variant","missense_variant","inframe_deletion","inframe_insertion","transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification"]
vep_lof_terms=["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification","missense_variant"]
################# END OF DEFINITIONS #################################


################# START MAIN FUNCTION ##########################
if __name__== "__main__":
	main_func(bed_file, base_quality_threshold, mapping_quality_threshold)

	





