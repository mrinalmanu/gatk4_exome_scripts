# gatk4_exome_scripts
gatk4 variant analysis in snakemake

Please go thorugh README_extended.txt for detailed description of the pipeline.

*YAB_corrections* contains the corrected scripts. Courtesy of Yury Barbitoff from JetBrains.

*Pure Snakemake* implementation contains a simpler workflow, free of pure python. I also added some scripts that can be used to fetch samples from the directory. They can be used to create a sample list, which can be added to *config.json* for a corrosponding sample. In this way, a sankefile and config for the given batch are constructed. Format: SNAKEFILE_batchname, and corrosponding CONFIG_batchname.json.

Now I am working on some *hail* script for analysis of variant data.
