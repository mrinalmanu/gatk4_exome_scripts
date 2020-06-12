# gatk4_exome_scripts
gatk4 variant analysis in snakemake

Please go thorugh README_extended.txt for detailed description of the pipeline.

*YAB_corrections* contains the corrected scripts. Courtesy of Yury Barbitoff from JetBrains. This folder contains the entire OLDER implementation of the pipeline.

*Pure Snakemake* implementation contains a simpler workflow, free of pure python. This is the NEWER implementation. I also added some scripts that can be used to fetch samples from the directory. They can be used to create a sample list, which can be added to *config.json* for a corrosponding sample. In this way, a sankefile and config for the given batch are constructed. Format: SNAKEFILE_batchname, and corrosponding CONFIG_batchname.json.

**Objectives**

1. Annotate file with VEP. We need the most comprehensive type of
annotation, including:
a. Variant effect per se
b. gnomAD allele frequency data
c. Clinical significance (ClinVar)

2. Load the annotated file into Hail.

3. Remove samples with low call-rate and low average GQ (after inspecting
the plots)

4. Remove variant sites with low call rate (after inspecting the plots)

Details of steps 1 through 4.

Tasks: Load the final annotated file into Hail. Remove all samples with mean
GQ &lt; 30
1. After removal of low GQ samples redo the call rate analysis.

2. Kinship analysis - it’s better these statistics using variant sites with the
highest call rates. (use variants with call rates &gt; 0.9, and AN field (one of
the INFO fields) - the number of samples * 2. - and subset the dataset (the
matrix table) using the AN field value.

3. After calculating the relatedness matrix, remove related samples.

4. Do a PCA analysis using genotype data (also use variants with call rates &gt;
0.9). To check whether the population is homogeneous or there are some
subgroups.

5. Analyze the correlation between allele frequencies in our dataset and
gnomAD_AF (it is best to draw several scatterplots - REX vs gnomAD;
REX vs gnomAD_NFE, REX vs gnomad_AFR, and others)

5. Investigate the relatedness metrics (kin statistics, overrepresentation
analysis)

6. Remove related samples

7. Make PCA analysis of the data (project individuals into the 2d space to
see whether they will group together into subpopulations). This should be
done using only sites with the maximum call rates (0.95 or more -
depending on the plots)
8. Investigate the frequency of all ClinVar pathogenic variants that are
linked to autosomal-recessive disorders according to OMIM (I’ll put
some data to GD for you to obtain lists of AR genes). For each such SNP,
compute the binomial overrepresentation p-value similar to the one we
calculated in original paper.

Description: We did a binomial test where - gnomad_NFE_AF is our success
probability. AN is our sample size, and AC is the number of successes;
we need to calculate the binomial probability of this event for all variants.
After calculation we can select only ClinVar pathogenic variants that are
associated with autosomal recessive disease.

9. Deploy ExAC-type variant server (see repository) - we need the very
basic things (only genes, dbSNP, and our variants - no coverage data)


