This is the tree structure:

.
├── 2_dry_run.sh
├── EXAMPLE_SAMPLE_LIST_TEXT_FILE.txt
├── gatk4_docker.py
├── make_varfile.py
├── preparing_environment
│   ├── Dockerfile
│   ├── pull_images.sh
│   └── run_ngs_container.sh
├── raw_data
├── README_mrinal.txt
├── run_pipeline.py
├── scripts
│   ├── config.py
│   ├── gatk4_germline.py
│   ├── gatk_varscan.py
│   ├── __init__.py
│   ├── inspect_coverage.py
│   ├── lock_module.py
│   └── store_annotations.py
└── write_vcf.py


Since we can use docker it's better to set up docker environment

Before launching the original script, we should do a dry run to see that everything goes smoothly.

Finally, we can launch the 3_real_run.sh script.

###############################################################################################################################
## Step 1:

_______________________________________________

docker image build -t [container name] --no-cache 

bash pull_images.sh

_______________________________________________

This will set up docker environment with the given set of dependencies.

############################################################################
## Step 2:

_______________________________________________


run_ngs_container.sh

_______________________________________________


This is the base container script, it should mount the following:
    
    • docker socker (-v /var/run/docker.sock:/var/run/docker.sock)
    • location where input data will be stored (container will mount this location as /home/input)
    • location where output data will be stored (container will mount this location as /home/output)
    • location with reference files (container will mount this location as /home/Reference)


mandatory flags:
-r: Reference dir 
-i: Input docker mount point folder
-o: Output docker mount point folder
-I: Base container image
-n: Container name (username)*

############################################################################
## Step 3:

Create sample list
Create a comma-separated text file with mandatory header. Next rows contain samples, one sample per row.
NOTE: input and output location is mounted at /home/input (output). Keep this in mind when specifying location of FASTQ files and specify only subfolders! Example: If your FASTQ is in /home/analysis/FASTQ/sample1/sample1.fastq, and your mounted input location is /home/analysis, then set your input location in sample file as /home/input/FASTQ/sample1/sample1.fastq
Header fields are:
SAMPLE_ID,		- any text string (no white spaces)
PIPELINE,		- pipeline_ID [gatk4_germline]
GENOME_BUILD,		- [HG1K37]
FASTQ1,			- location of normal fastq1, if tumor-only run: leave empty
FASTQ2,			- location of normal fastq2, if tumor-only run: leave empty
FASTQ1_TUMOR,		- location of tumor fastq1, if germline run: leave empty
FASTQ2_TUMOR,		- location of tumor fastq2, if germline run: leave empty
CLEANUP,		- set YES is wish to remove intermediate BAM
LIB_ID,			- library_ID any text string (not important)
PL_ID,			- PL string (not important)
PU_ID			- PU string (not important)


An example entry:

SAMPLE_ID,PIPELINE,GENOME_BUILD,FASTQ1,FASTQ2,FASTQ1_TUMOR,FASTQ2_TUMOR,CLEANUP LIB_ID,PL_ID,PU_ID
Sample1,gatk4_germline,HG1K37,NORMAL_R1.fastq,NORMAL_R2.fastq,TUMOUR_R1.fastq,TUMOUR_R2.fastq,NO,lib1,Illumina,XXX

############################################################################
## Step 4:

It contains names of indexes and reference genome fasta files in Reference folder. 

We have to make sure these point to existing files in your mounted Reference directory. Maximum threads available can be set by changing of variable "max_nr_threads".

############################################################################
## Step 5:

_______________________________________________

snakemake -ps $PWD/run_pipeline.py sample_list.txt

_______________________________________________

Script has three parts: first parses docker_images_versions.txt, second parses sample list and creates an instance of class of respective pipeline for every sample in sample list. Third part checks available cores and submits samples for processing by invoking run_pipeline() function.




