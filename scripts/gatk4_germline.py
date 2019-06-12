### GATK4 germline pipeline
### Last update: May 1, 2019

import sys
import subprocess
import logging
import os
from gatk4_docker import gatk_docker
from multiprocessing import Process
from config_file import *
from lock_module import *
import time
import random

def check_path(path):                		## Check if file exists
			if os.path.exists(path) == True:
				return True
			else:
				return False
				
def checkContainer(container_name):
		containers = subprocess.check_output(['docker','ps','-a']).decode(encoding="437")
		if container_name in containers.split():
			subprocess.call(["docker","rm","-f",container_name])
		else:
			return				

class GATK4_germline():
	def __init__(self, sample_name, genome_build, f1, f2, tf1, tf2, cleanup,lib_ID, pl_ID, pu_ID, docker_images_dict):
		self.sample_name = sample_name
########################################################################################################################

# We need to feed fastq samples

		self.fastq1_tumor = tf1
		self.fastq2_tumor = tf2
		self.fastq1 = f1
		self.fastq2 = f2

#########################################################################################################################
		if (f1 =="") or (f2 ==""):
			print (self.sample_name+": ERROR: FASTQ file is missing, aborting.")
			quit()
		if (tf1 !="") or (tf2 !=""):
			print (self.sample_name+": ERROR: tumor FASTQ file specified in germline pipeline, aborting.")
			quit()
		self.input_folder=input_folder
		self.threads = max_threads_process
		self.paralellized_threads = max_paralellized_threads
		self.ram = "50000"
		self.open_files = []
		self.genome_build = genome_build
		self.cleanup = cleanup
		self.docker_images_dict = docker_images_dict
		self.lib_ID = lib_ID
		self.pl_ID = pl_ID
		self.pu_ID = pu_ID
		
		self.output_folder = output_folder+"/"+self.sample_name+"/GATK4_germline/output/"
		if check_path(self.output_folder)==False:
			print ("MKDIR", self.output_folder)			
			os.system("mkdir -p "+self.output_folder)		
		os.system("chmod 777 "+self.output_folder)	
		
		# logger setup
		LOG_FORMAT = "%(levelname)s\t%(asctime)s\t%(module)s\t%(funcName)s\t%(message)s"
		logging.basicConfig(filename = None,
												level = logging.DEBUG,
												format = LOG_FORMAT)
		self.logger = logging.getLogger()
	
	
	def run_in_docker(self, cmd, image, threads_needed, stdout=None, stderr=None):
		""" Run a command inside docker container"""
		if len(cmd)<2:
			container_name = self.sample_name+"_vt"+"_"+''.join(random.choice('abcdefghijklmnoprstuvxyz') for _ in range(3))
		else:
			container_name = self.sample_name+"_"+cmd[0]+"_"+''.join(random.choice('abcdefghijklmnoprstuvxyz') for _ in range(3))
		checkContainer(container_name)
	
		
		dcmd = ["docker", "run","--name",container_name,
						"-v", "{}:{}".format(self.input_folder, self.input_folder),
						"-v", "{}:{}".format(output_folder, output_folder),
						"-v", "{}:{}".format(reference_folder, reference_folder),
						self.docker_images_dict[image]
						]
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
			wait_go = check_threads(batch_ID, container_name,"start",(-1)*threads_needed, max_nr_threads)	
		
		########### RUN COMMAND IN DOCKER #################
		self.logger.info("Running: "+" ".join(dcmd))
		os.system(" ".join(dcmd))
		checkContainer(container_name)
		
		########### RELEASE LOCK DOWN CORES ###############
		wait_go = "wait"
		while wait_go !='go':
			time.sleep(10)
			wait_go = check_threads(batch_ID, container_name,"finish", threads_needed, max_nr_threads)	
	###################### END OF RUN COMMANDS IN DOCKER CONTAINERS ###################
			
	def run_pipeline(self):
		### Genome build selection ###
		if self.genome_build == "HG1K37":
			bwa_index = bwa_index_HG1K37
			reference_fasta = reference_fasta_HG1K37
			dbsnp_vcf = dbsnp_vcf_HG1K37
		#if self.genome_build == "hg19":
		#	bwa_index = bwa_index_hg19
		#	reference_fasta = reference_fasta_hg19
		#	dbsnp_vcf = dbsnp_vcf_hg19
		
		### Pipeline auxiliary files ###
		unsorted_bam = self.output_folder+self.sample_name+".bwa.bam"
		sorted_bam = self.output_folder+self.sample_name+".bwa.sorted.bam"
		sort_err = self.output_folder+self.sample_name+".sort.err"
		markDuplicates_bam = self.output_folder+self.sample_name+".MarkDuplicates.bam"
		markDuplicates_metrics = self.output_folder+self.sample_name+".MarkDuplicates-metrics.txt"
		ReadGroups_bam = self.output_folder+self.sample_name+".ReadGroups.bam"
		BaseRecalibrator_metrics = self.output_folder+self.sample_name+".BaseRecalibrator-metrics.txt"
		BaseRecalibrator_bam = self.output_folder+self.sample_name+".BaseRecalibrator.bam"
		vcf = self.output_folder+self.sample_name+".vcf.gz"
		########################################
		
		
		
		def align_BWA():
			# Align to reference genome (bwa) and sort (samtools)
			print ""
			print "STARTING BWA ALIGNMENT"
			print ""
			self.run_in_docker(["python align_bwa_sort_samtools.py "
				+reference_fasta+" "
				+self.fastq1+" "
				+self.fastq2+" "
				+self.threads+" "
				+self.output_folder+self.sample_name+" "
				+unsorted_bam],
				"bwa_samtools131",
				self.threads)
			
		def markDuplicates(input_bam):
			# Mark duplicates (GATK)
			parameters_dict = {
								"input":input_bam,
								"output":input_bam.replace("bam","MarkDuplicates.bam"),
								"mark_dupl_metrics":markDuplicates_metrics,
								"optical_duplicate_pixel_dist":"2500",
								"assume_sort_order":"queryname",
								"clear_DT":"false",
								"add_pg_tag_to_reads":"false"
									}
			markDuplicates_log = self.output_folder+self.sample_name+".MarkDuplicates.log"
			gatk_docker("gatk_mark_duplicates",parameters_dict,
						markDuplicates_log, self.ram,self.docker_images_dict["gatk"])
		
		def addReadGroups():
			# Add read groups (GATK)
			print ""
			print "STARTING ADDING READGROUPS"
			print ""
			parameters_dict = {
								"input":unsorted_bam,
								"output":ReadGroups_bam,
								"RGLB":self.lib_ID,
								"RGPL":self.pl_ID,
								"RGPU":self.pu_ID,
								"RGSM":self.sample_name
								}
			ReadGroups_log = self.output_folder+self.sample_name+".ReadGroups.log"
			gatk_docker("gatk_add_read_groups", parameters_dict,
						ReadGroups_log, self.ram,self.docker_images_dict["gatk"])
	
		def buildRecalibrator():
			# Base recalibrator (GATK)
			print ""
			print "STARTING BUILDING RECALIBRATOR"
			print ""
			parameters_dict = {
								"input":sorted_bam,
								"output":BaseRecalibrator_metrics,
								"reference":reference_fasta,
								"known-sites":dbsnp_vcf,
								"use-original-qualities":"true"
								}
			BaseRecalibrator_log = output_folder+self.sample_name+".BaseRecalibrator.log"
			gatk_docker("gatk_build_recalibrator", parameters_dict,
						BaseRecalibrator_log, self.ram,self.docker_images_dict["gatk"])

			
		def applyRecalibrator():
			# Base recalibrator - applying model (GATK)
			print ""
			print "STARTING APPLY RECALIBRATOR"
			print ""
			parameters_dict = {
								"input":sorted_bam,
								"output":BaseRecalibrator_bam,
								"bqsr":BaseRecalibrator_metrics,
								"use-original-qualities":"true",
								"static-quantized-quals":"10",
								"static-quantized-quals":"20",
								"static-quantized-quals":"30"
								}
			ApplyBQSR_log = output_folder+self.sample_name+".ApplyBQSR.log"
			gatk_docker("gatk_apply_recalibrator", parameters_dict,
						ApplyBQSR_log, self.ram,self.docker_images_dict["gatk"])
		
		def variantCalling_HaplotypeCaller(input_bam):
			# Haplotype caller
			parameters_dict = {
								"input":input_bam,
								"output":input_bam.replace(".bam",".vcf.gz"),
								"reference":reference_fasta
							#	"emit-reference-confidence":"true",
							#	"max-alternate-alleles":"3"
								}
			HaplotypeCaller_log = input_bam+".HaplotypeCaller.log"
			gatk_docker("gatk_haplotype_caller", parameters_dict,
						HaplotypeCaller_log, self.ram,self.docker_images_dict["gatk"])

			
		def validateSam():
			# Validate alignment file integrity
			parameters_dict = {
								"input":output_folder+self.sample_name+".bwa.bam",
								"MODE":"SUMMARY"
								}
			validate_log = output_folder+self.sample_name+".validateSamFile.log"
			gatk_docker(gatk_validate_sam, parameters_dict,
						HaplotypeCaller_log, self.ram,self.docker_images_dict["gatk"])

		def indexBam(path):
			if check_path(path) == True:
				index_cmd = ['samtools','index',path]
				index_log= output_folder+self.sample_name+".samtoolsindex.log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(index_cmd,"samtools",1)
			else:
				raise RuntimeError("BaseRecalibrator bam file not found")
				quit()
		def sortBam(bam_to_sort, sorted_bam, t_n):
			if check_path(bam_to_sort) == True:
				sort_cmd = ['samtools','sort', bam_to_sort]
				sort_cmd += ['-o', sorted_bam]
				sort_cmd += ['-@',self.threads]
				sort_log= self.output_folder+self.sample_name+".samtoolssort.log"
				#self.run_in_docker(index_cmd, stderr=index_log)
				self.run_in_docker(sort_cmd,"samtools",self.threads)
			else:
				raise RuntimeError("unsorted bam file not found")
				quit()		
		
		def processVcf(input_vcf, filename):
			cmd= ['/bin/bash -c "zcat '+input_vcf+' | vt decompose -s - | vt normalize -q - -n -r '+reference_fasta+' 2> /dev/null " | bgzip -c > '+filename+'normalized.vcf.gz']
			self.run_in_docker(cmd,"vt",1)
		
		def sortVcf(vcf_file):
			# index vcf file
			index_cmd=["bcftools","sort","-O z","-o",vcf_file.replace(".vcf.gz",".sort.vcf.gz"),vcf_file]
			self.run_in_docker(index_cmd,"bcftools",1)
		
		def indexVcf(vcf_file):
			# index vcf file
			index_cmd=["bcftools","index","-t","-f",vcf_file]
			self.run_in_docker(index_cmd,"bcftools",1)
		
		def extractChromosome(bam_file, chromosome):
			index_cmd = ['samtools','view', bam_file, "chr"+chromosome, "-b","-o", bam_file.replace("bam","chr"+chromosome+".bam")]
			index_cmd += ['-@', self.paralellized_threads]
			index_log= self.output_folder+self.sample_name+".samtools_split.log"
			#self.run_in_docker(index_cmd, stderr=index_log)
			self.run_in_docker(index_cmd,"samtools",self.paralellized_threads)
		
		def joinChromosomeVcfs():
			joint_vcf = self.output_folder+self.sample_name+".HaplotypeCaller.vcf.gz"
			partial_vcf_list = []
			for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
				partial_vcf_list.append(BaseRecalibrator_bam.replace("bam","chr"+chromosome+".MarkDuplicates.vcf.gz"))
			self.run_in_docker(["bcftools","concat","-o",joint_vcf,"-O","z"," ".join(partial_vcf_list)], "bcftools", self.threads)	
					
		#***********************************
		# Pipeline workflow ###
		#***********************************
		
		# Check if alignment file exists
		if check_path(self.output_folder+self.sample_name+".bwa.bam")==False:
			self.logger.info(self.sample_name+" sample: align")
			align_BWA()
		
#		else:
#			self.logger.info("alignment file present, checking integrity")
#			validateSam()
		
		# Check if readgroups were added
		if check_path(self.output_folder+self.sample_name+".ReadGroups.bam")==False:
			self.logger.info(self.sample_name+" sample: Add readgroups")
			addReadGroups()
		else:	
			self.logger.info("readgroups ok, skipping...")
		#Sort Bam
		self.logger.info("sort")
		sortBam(ReadGroups_bam, sorted_bam, "foo")
		print ("")
		
		# Check if recalibrator model was built
		if check_path(self.output_folder+self.sample_name+".BaseRecalibrator-metrics.txt")==False:
			buildRecalibrator()
		else:	
			self.logger.info("recalibrator model already build, skipping...")	

		# Check if recalibrator model was applied
		if check_path(self.output_folder+self.sample_name+".BaseRecalibrator.bam")==False:
			self.logger.info(self.sample_name+" sample: Base recalibrator apply")
			applyRecalibrator()
		else:	
			self.logger.info("recalibrator model already applied, skipping...")		
		 
		# Index bam before variant calling
		indexBam(BaseRecalibrator_bam)
		
		# Split BAM to chromosomes for parallel variant calling
		
		proc = []
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			p = Process(target = extractChromosome, args=[BaseRecalibrator_bam, chromosome])
			p.start()
			proc.append(p)
		for p in proc:
			p.join()	
		
		proc = []	
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			p = Process(target = indexBam, args=[BaseRecalibrator_bam.replace("bam","chr"+chromosome+".bam")])
			p.start()
			proc.append(p)
		for p in proc:
			p.join()
			
		proc = []
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			p = Process(target = markDuplicates, args=[BaseRecalibrator_bam.replace("bam","chr"+chromosome+".bam")])
			p.start()
			proc.append(p)
		for p in proc:
			p.join()
		
		proc = []	
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			p = Process(target = indexBam, args=[BaseRecalibrator_bam.replace("bam","chr"+chromosome+".MarkDuplicates.bam")])
			p.start()
			proc.append(p)
		for p in proc:
			p.join()
		
		proc = []
		for chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']:
			p = Process(target = variantCalling_HaplotypeCaller, args=[BaseRecalibrator_bam.replace("bam","chr"+chromosome+".MarkDuplicates.bam")])
			p.start()
			proc.append(p)
		
		for p in proc:
			p.join()
		
		# # Join chromosome vcf together
		
		joinChromosomeVcfs()
			
		# PROCESS VCF in VT (check if already processed)
		if check_path(self.output_folder+self.sample_name+".HaplotypeCaller.normalized.vcf.gz")==False:
			vcf_file = self.output_folder+self.sample_name+".HaplotypeCaller.vcf.gz"
			if check_path(vcf_file)==True:
				processVcf(vcf_file, vcf_file.replace("vcf.gz",""))
		else:
			self.logger.info(self.sample_name+": vcf processing already done, skipping...")
		
		sortVcf(self.output_folder+self.sample_name+".HaplotypeCaller.normalized.vcf.gz")
		
		# Index normalized vcf (bcftools index vcf)
		if check_path(self.output_folder+self.sample_name+".HaplotypeCaller.normalized.sort.vcf.gz.tbi")==False:
			self.logger.info(self.sample_name+" sample: indexing")		
			indexVcf(self.output_folder+self.sample_name+".HaplotypeCaller.normalized.sort.vcf.gz")
		else:
			self.logger.info(self.sample_name+" sample: normalized vcf indexed, skipping...")		
			print ("")
		
		# Remove intermediate GATK-produced bams
		if self.cleanup=="YES":
			try:
				subprocess.Popen("rm "+self.output_folder+self.sample_name+".MarkDuplicates.bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+self.sample_name+".ReadGroup.bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
			try:
				subprocess.Popen("rm "+self.output_folder+self.sample_name+".BaseRecalibrator.bam",shell=True,stderr=subprocess.PIPE)
			except:
				pass
				
		self.logger.info("sample "+self.sample_name+" pipeline finished.")	
