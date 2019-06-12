### GATK4 somatic pipeline for tumor-normal pairs
### Last update: May 1, 2019

import sys
import subprocess
import logging
import os
from gatk4_docker import gatk_docker
from config_file import *

def check_path(path):                		## Check if file exists
		if os.path.exists(path) == True:
			return True
		else:
			return False

def checkContainer(container_name):
		containers = subprocess.check_output(['docker','ps','-a']).decode(encoding="437")
		if container_name in containers.split():
			subprocess.call(["docker","rm",container_name])
		else:
			return
			
class GATK_varscan_TN_v1():
	def __init__(self, sample_name, genome_build, f1, f2,tf1,tf2, cleanup,lib_ID, pl_ID, pu_ID, docker_images_dict):

		self.sample_name= sample_name
		self.fastq1_tumor = tf1
		self.fastq2_tumor = tf2
		self.fastq1 = f1
		self.fastq2 = f2
		if (f1 =="") or (f2 =="") or (tf1 =="") or (tf2 ==""):
			print >>sys.stderr, t_n_sample_name+": FASTQ file is missing, aborting."
			quit()
			
		self.threads = max_nr_threads
		self.ram = "50000"
		self.open_files = []
		self.genome_build = genome_build
		self.cleanup = cleanup
		self.docker_images_dict = docker_images_dict
		self.lib_ID = lib_ID
		self.pl_ID = pl_ID
		self.pu_ID = pu_ID
		self.output_folder = output_folder+"RUN/"+self.sample_name+"/output/"
		
		if check_path(self.output_folder)==False:
			print "MKDIR", self.output_folder			
			os.system("mkdir -p "+self.output_folder)		
		os.system("chmod 777 "+self.output_folder)			

		# logger setup
		LOG_FORMAT = "%(levelname)s\t%(asctime)s\t%(module)s\t%(funcName)s\t%(message)s"
		logging.basicConfig(filename = None,
												level = logging.DEBUG,
												format = LOG_FORMAT)
		self.logger = logging.getLogger()
		
	###################### RUN COMMANDS IN DOCKER CONTAINERS ###################
	def run_in_docker(self, cmd,t_n, image, stdout=None, stderr=None):
		""" Run a command inside docker container"""
	
		container_name = self.sample_name+"_"+t_n+"_"+cmd[0]
		checkContainer(container_name)
		dcmd = ["docker","run","--name",container_name,
						"-v", "{}:{}".format(self.input_folder, self.input_folder),
						"-v", "{}:{}".format(output_folder, output_folder),
						"-v", "{}:{}".format(reference_folder, reference_folder),
						image]
		dcmd += cmd
	
		if stdout is not None:
			stdout = open(stdout,"w")
			self.open_files.append(stdout)
		if stderr is not None:
			stderr = open(stderr,"w")
			self.open_files.append(stderr)
		
		############ LOCK DOWN NEEDED THREADS #############
		threads_needed = int(threads_needed)
		wait_go = "wait"
		while wait_go !='go':
			time.sleep(10)
			wait_go = check_threads(batch_ID, container_name,"start",(-1)*threads_needed, max_threads)	
		
		########### RUN COMMAND IN DOCKER #################
		self.logger.info("Running: "+" ".join(dcmd))
		os.system(" ".join(dcmd))
		checkContainer(container_name)
		
		########### RELEASE LOCK DOWN CORES ###############
		wait_go = "wait"
		while wait_go !='go':
			time.sleep(10)
			wait_go = check_threads(batch_ID, container_name,"finish", threads_needed, max_threads)	
	###################### END OF RUN COMMANDS IN DOCKER CONTAINERS ###################
	
	def run_pipeline(self):
		#def run pipeline(self)
		
		### Genome build selection ###
		if self.genome_build == "GRCh38":
			bwa_index = bwa_index_GRCh38
			reference_fasta = reference_fasta_GRCh38
			dbsnp_vcf = dbsnp_vcf_GRCh38
		if self.genome_build == "hg19":
			bwa_index = bwa_index_hg19
			reference_fasta = reference_fasta_hg19
			dbsnp_vcf = dbsnp_vcf_hg19
		
		def align_BWA(t_n):
			if t_n =="tumor":
				fastq1 = self.fastq1_tumor
				fastq2 = self.fastq2_tumor
			if t_n =="normal":
				fastq1 = self.fastq1
				fastq2 = self.fastq2	
			# Align to reference genome (bwa) and sort (samtools)
		
			self.run_in_docker(["python align_bwa_sort_samtools.py "
				+reference_fasta+" "
				+fastq1+" "
				+fastq2+" "
				+self.threads+" "
				+self.output_folder+t_n_sample_name+" "
				+unsorted_bam],
				t_n,
				"bwa_samtools131")
			
		def markDuplicates(t_n):
			# Mark duplicates (GATK)
			parameters_dict = {
								"input":unsorted_bam,
								"output":markDuplicates_bam,
								"mark_dupl_metrics":markDuplicates_metrics,
								"optical_duplicate_pixel_dist":"2500",
								"assume_sort_order":"queryname",
								"clear_DT":"false",
								"add_pg_tag_to_reads":"false"
									}
			markDuplicates_log = self.output_folder+t_n_sample_name+".MarkDuplicates.log"
			gatk_docker("gatk_mark_duplicates",parameters_dict,
						markDuplicates_log, self.ram,self.docker_images_dict["gatk"])
		
		def addReadGroups(t_n):
			# Add read groups (GATK)
			parameters_dict = {
								"input":markDuplicates_bam,
								"output":ReadGroups_bam,
								"RGLB":self.lib_ID,
								"RGPL":self.pl_ID,
								"RGPU":self.pu_ID,
								"RGSM":self.sample_name+"_r_"+t_n
								}
			ReadGroups_log = self.output_folder+t_n_sample_name+".ReadGroups."+t_n+".log"
			gatk_docker("gatk_add_read_groups", parameters_dict,
						ReadGroups_log, self.ram,self.docker_images_dict["gatk"])
	
		def buildRecalibrator(t_n):
			# Base recalibrator (GATK)
			parameters_dict = {
								"input":sorted_bam,
								"output":BaseRecalibrator_metrics,
								"reference":reference_fasta,
								"known-sites":dbsnp_vcf,
								"use-original-qualities":"true"
								}
			BaseRecalibrator_log = self.output_folder+t_n_sample_name+".BaseRecalibrator."+t_n+".log"
			gatk_docker("gatk_build_recalibrator", parameters_dict,
						BaseRecalibrator_log, self.ram,self.docker_images_dict["gatk"])

			
		def applyRecalibrator(t_n):
			# Base recalibrator - applying model (GATK)
			parameters_dict = {
								"input":sorted_bam,
								"output":BaseRecalibrator_bam,
								"bqsr":BaseRecalibrator_metrics,
								"use-original-qualities":"true",
								"static-quantized-quals":"10",
								"static-quantized-quals":"20",
								"static-quantized-quals":"30"
								}
			ApplyBQSR_log = self.output_folder+t_n_sample_name+".ApplyBQSR."+t_n+".log"
			gatk_docker("gatk_apply_recalibrator", parameters_dict,
						ApplyBQSR_log, self.ram,self.docker_images_dict["gatk"])
		
		def variantCalling_mutect2():
			# Haplotype caller
			vcf = self.output_folder+self.sample_name+".mutect2.vcf.gz"
			parameters_dict = {
								"input_tumor":self.output_folder+self.sample_name+"_r_tumor.bam",
								"input_normal":self.output_folder+self.sample_name+"_r_normal.bam",
								"normal-sample":self.sample_name+"_r_normal",
								"tumor-sample":self.sample_name+"_r_tumor",
								"output":vcf,
								"reference":reference_fasta,
							#	"--disable-read-filter":""
							#	"dbsnp":
							#	"emit-reference-confidence":"true",
							#	"max-alternate-alleles":"3"
								}
			mutect2_log = self.output_folder+self.sample_name+".Mutect2.log"
			gatk_docker("gatk_mutect2", parameters_dict,
						mutect2_log, self.ram,self.docker_images_dict["gatk"])
		def variantCalling_varscan2():
			varscan_cmd = ["java", "-jar" ,"VarScan.jar","somatic", self.output_folder+self.sample_name+"_r_normal.bam.pileup",self.output_folder+self.sample_name+"_r_tumor.bam.pileup", self.output_folder+self.sample_name,"--output-vcf"]
			index_log= self.output_folder+self.sample_name+".varscan2.log"
			self.run_in_docker(varscan_cmd,t_n,self.docker_images_dict["varscan"])
			
		def validateSam(t_n):
			# Validate alignment file integrity
			parameters_dict = {
								"input":self.output_folder+t_n_sample_name+".bwa."+t_n+".bam",
								"MODE":"SUMMARY"
								}
			validate_log = self.output_folder+t_n_sample_name+".validateSamFile."+t_n+".log"
			gatk_docker(gatk_validate_sam, parameters_dict,
						HaplotypeCaller_log, self.ram,self.docker_images_dict["gatk"])
		
		def sortBam(unsorted_bam, sorted_bam, t_n):
			if check_path(unsorted_bam) == True:
				sort_cmd = ['samtools','sort', unsorted_bam]
				sort_cmd += ['-o', sorted_bam]
				sort_log= self.output_folder+self.sample_name+".samtoolssort."+t_n+".log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(sort_cmd,t_n, self.docker_images_dict["samtools"])
			else:
				raise RuntimeError("unsorted bam file not found")
				quit()
		
		
		def indexBam(path, t_n):
			if check_path(path) == True:
				index_cmd = ['samtools','index',path]
				index_log= self.output_folder+self.sample_name+".samtoolsindex."+t_n+".log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(index_cmd,t_n,self.docker_images_dict["samtools"])
			else:
				raise RuntimeError("BaseRecalibrator bam file not found")
				quit()
		
		def processVcf(input_vcf, filename):
			cmd= ['/bin/bash -c "zcat '+input_vcf+' | vt decompose -s - | vt normalize -q - -n -r '+reference_fasta+' 2> /dev/null " | bgzip -c > '+filename+'normalized.vcf.gz']
			self.run_in_docker(cmd, "asa","vt")
		
		def createPileups(BaseRecalibrator_bam, t_n):
			if check_path(BaseRecalibrator_bam) == True:
				pileup_cmd = ['samtools','mpileup','-b',self.output_folder+"list.txt",'-f',reference_fasta,'-o',BaseRecalibrator_bam+".pileup"]
				pileup_log= self.output_folder+self.sample_name+".samtoolsmpileup."+t_n+".log"
				os.system("echo "+ BaseRecalibrator_bam+" > "+self.output_folder+"list.txt")
				self.run_in_docker(pileup_cmd,t_n,self.docker_images_dict["samtools"])
			else:
				raise RuntimeError("BaseRecalibrator bam file not found")
				quit()
		def create_varfile(somatic_vcf):
			cmd = "python make_varfile.py "+somatic_vcf
			os.system(cmd)
			
		def modify_pileup(tumor_pileup):
			os.system("awk 'BEGIN{OFS ="+r'"\t"'+";ORS="+r'"\n"'+"}{if ($4>0) print $1,$2,$3,$4,$5,$6}' "+tumor_pileup+" > "+tumor_pileup+".mod")
		
		def generate_readcounts(modified_pileup, varfile):
			t_n ="varscan" 
			varscan_cmd = ["java", "-jar" ,"VarScan.jar","readcounts"]
			varscan_cmd += [modified_pileup,"--min-coverage 0", "--min-base-qual 15", "--output-file", self.sample_name+"readcounts", "--variants-file", varfile]
			index_log= self.output_folder+self.sample_name+".varscan2.log"
			self.run_in_docker(varscan_cmd,t_n,self.docker_images_dict["varscan"])	
		
		
		def write_vcf(read_counts, sample_name):
			cmd = "python write_vcf.py "+modified_pileup+" "+sample_name	
			os.system(cmd)
		
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
							"--dir_cache", vep_reference_folder, 
							"--fork", self.threads, 
							"--input_file",vcf_file, 
							"--force_overwrite", 
							"--fasta", reference_fasta, 
							"-o", vcf_file.replace(".vcf.gz",".annotated.vcf.gz"), 
							"--vcf"]
								
			self.run_in_docker(annotate_cmd,"asa", "vep")		
		#***********************************
		# Pipeline workflow ###
		#***********************************
		
		#### pre-processing for tumor and normal separately ###
		for t_n in ["tumor","normal"]:
			t_n_sample_name= self.sample_name+"_"+t_n
			### Pipeline auxiliary files ###
			unsorted_bam = self.output_folder+self.sample_name+".bwa."+t_n+".bam"
			sorted_bam = self.output_folder+self.sample_name+".bwa.sorted."+t_n+".bam"
			sort_err = self.output_folder+self.sample_name+".sort."+t_n+".err"
			markDuplicates_bam = self.output_folder+self.sample_name+".MarkDuplicates."+t_n+".bam"
			markDuplicates_metrics = self.output_folder+self.sample_name+".MarkDuplicates-metrics."+t_n+".txt"
			ReadGroups_bam = self.output_folder+self.sample_name+".ReadGroups."+t_n+".bam"
			BaseRecalibrator_metrics = self.output_folder+self.sample_name+".BaseRecalibrator-metrics."+t_n+".txt"
			BaseRecalibrator_bam = self.output_folder+self.sample_name+"_r_"+t_n+".bam"
			tumor_pileup = self.output_folder+self.sample_name+"_r_tumor.bam.pileup"
			somatic_vcf = self.output_folder+self.sample_name+".snp.vcf"
			varfile = somatic_vcf.replace(".vcf",".varfile")
			modified_pileup = tumor_pileup+".mod"
			########################################
			
			# ALIGNMENT BWA (Check if alignment file exists)
			if check_path(self.output_folder+self.sample_name+".bwa."+t_n+".bam")==False:
				self.logger.info("align")
				align_BWA(t_n)
			
			else:
				self.logger.info(t_n+" sample:alignment done, skipping...")
			print ""
			
			# MARK DUPLICATES (Check if duplicated reads were marked)
			if check_path(self.output_folder+self.sample_name+".MarkDuplicates."+t_n+".bam")==False:
				self.logger.info("mark_dupl")
				markDuplicates(t_n)
			else:	
				self.logger.info(t_n+" sample:duplicates marked, skipping...")
			print ""
			
			# ADD READGROUPS (Check if readgroups were added)
			if check_path(self.output_folder+self.sample_name+".ReadGroups."+t_n+".bam")==False:
				self.logger.info("add_readgroups")
				addReadGroups(t_n)
			else:	
				self.logger.info(t_n+" sample: readgroups ok, skipping...")
			print ""
			
			#Sort Bam
	#		self.logger.info("sort")
	#		sortBam(ReadGroups_bam, sorted_bam, t_n)
	#		print ""
			
			# BASE RECALIBRATOR (Check if recalibrator model was built)
			if check_path(self.output_folder+self.sample_name+".BaseRecalibrator-metrics."+t_n+".txt")==False:
				self.logger.info("build_recalib")
				buildRecalibrator(t_n)
			else:	
				self.logger.info(t_n+" sample: recalibrator model already build, skipping...")	
			print ""
			
			# APPLY BQSR (Check if recalibrator model was applied)
			if check_path(BaseRecalibrator_bam)==False:
				self.logger.info("apply_recalib")		
				applyRecalibrator(t_n)
			else:	
				self.logger.info(t_n+" sample: recalibrator model already applied, skipping...")		
			print ""
			 
			# Index bam before variant calling
			if check_path(BaseRecalibrator_bam+".bai")==False:
				self.logger.info(t_n+" sample: indexing")		
				indexBam(BaseRecalibrator_bam, t_n)
			else:
				self.logger.info(t_n+" sample: indexed, skipping...")		
			print ""
			
			# CREATE PILEUPS
			if check_path(BaseRecalibrator_bam+".pileup")==False:
				self.logger.info(t_n+" sample: creating pileups")		
				createPileups(BaseRecalibrator_bam, t_n)
			else:
				self.logger.info(t_n+" sample: pileups created, skipping...")	
			
			#### END OF pre-processing for tumor and normal separately ###
	
		### TUMOR- NORMAL VARIANT CALLING ###########################################################
		# VARIANT CALLING VARSCAN (Check if variant calling was performed)
		
		if check_path(self.output_folder+self.sample_name+".snp.vcf")==False:
			self.logger.info("VarScan2 variant calling")
			variantCalling_varscan2()
		else:
			self.logger.info(t_n+" sample: variant calling done, skipping...")	
		
		# PROCESS VCF
		vcf_file = self.output_folder+self.sample_name+".mutect2.vcf.gz"
		if check_path(vcf_file)==True:
			processVcf(vcf_file, vcf_file.replace("vcf.gz",""))
			
		# FILTER CALLS
		
		
		
		# CREATE VARFILE
		create_varfile(somatic_vcf)
		
		# MODIFY PILEUP
		if check_path(modified_pileup)==False:
			modify_pileup(tumor_pileup)
		else:
			self.logger.info("pileup modified, skipping...")		
		
		#GENERATE READCOUNTS
		generate_readcounts(modified_pileup, varfile)
		
		#WRITE VCF		
		write_vcf(tumor_pileup+".mod", self.sample_name)
		
		## ANNOTATE CALLS WITH VEP (Check if already annotated)
		if check_path(self.output_folder+self.sample_name+".mutect2.normalized.filtered.vep.vcf.gz")==False:		
			self.logger.info(self.sample_name+" annotating with VEP")
			annotateVep(self.output_folder+self.sample_name+".w_alt_alleles.vcf")
		else:
			self.logger.info(self.sample_name+": annotated, skipping...")
		
		
		
			
		# Remove intermediate GATK-produced bams
		if self.cleanup=="YES":
			try:
				subprocess.Popen("rm "+self.output_folder+t_n_sample_name+".MarkDuplicates."+t_n+".bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+t_n_sample_name+".ReadGroup."+t_n+".bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+t_n_sample_name+".BaseRecalibrator."+t_n+".bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
				
		self.logger.info("sample "+self.sample_name+"_"+t_n+" pipeline finished.")	
	
