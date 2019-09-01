reference="A.fa"
enrichment="1.f"
PWD = "."
paddingcohort = "???"
SNPVQSLOD = "???"
INDELVQSLOD="???"
BUNDLE = " "
DBSNP = " "
CLINVAR = " "
LACALAF = " "
LOCALAFDB = " "
GWASCAT = " "
SNPSIFT = " "
KGENOMES = " "
DBNSFP = " "
SNPEFF = " "
run="1"
sample="2"

R_list = ["R1","R2"]
run_list = ["1","2"]
samples_list = ["1","2"]
lanes_list = ["1","2"]

rule all:
    input:
        expand("fastqs/run_{run}.sample_{sample}.lane_{lane}.{R}.fastq.gz", R = R_list, run = run_list, sample = samples_list, lane = lanes_list) 
        
rule run_bwamem:
    input: 
        expand("fastqs/run_{run}.sample_{sample}.lane_{lane}.{R}.fastq.gz", R = R_list, run = run_list, sample = samples_list, lane = lanes_list) 
    output:
        expand("bams/run_{run}.sample_{sample}.lane_{lane}.bam", run = run_list, sample = samples_list,lane = lanes_list) 
    params:
        rg = "@RG\tID:run_{run}.sample_{sample}.lane_{lane}\tSM:S{sample}\tLB:1\tPL:illumina"
    shell:
        #"bwa mem -M -t 12 -R '{params.rg}' {input.fastq} {input.ref} | samtools view -Sb -> {output}"
        "mv {input} {output}"
rule sort_bams:
    input: rules.run_bwamem.output
    output: "bams/run_{run}.sample_{sample}.lane_{lane}.sorted.bam"
    shell:
        "samtools sort -T run_{run}.sample_{sample}.lane_{lane}.sorted -o {output} {input}"

rule rewrite_bams:
    input: rules.sort_bams.output
    output: "bams/run_{run}.sample_{sample}.lane_{lane}.bam"
    shell:
        "mv {input} {output}"

rule merge_bams:
    input: rules.rewrite_bams.output
    output: "bams/run_{run}.sample_{sample}.bam"
    shell:
        "samtools merge {output} {input}"
        "rm *.lane_*"

rule mark_dups:
    input: "bams/run_{run}.sample_{sample}.bam"
    output: "bams/run_{run}.sample_{sample}.dedup.bam"
    shell:
        "Java -Xmx2g -jar Picard MarkDuplicates I={input} O={output} ASSUME_SORTED=true TMP_DIR=./tmp"

rule collect_hsm:
    input: rules.mark_dups.output
    output: "/hs_metrics/run_{run}.sample_{sample}.HS.metrics"
    shell:
        "Java -Xmx2g -jar Picard CalculateHsMetrics BAIT_INTERVALS={enrichment} TARGET_INTERVALS={enrichment} INPUT={input} OUTPUT=/hs_metrics/{output}.HSMetrics NEAR_DISTANCE=0"

rule index_bams:
    output: "bams/run_{run}.sample_{sample}.bam"
    shell:
        "Samtools index {input}"
        "touch .status"

rule realign_indels:
    input: "bams/run_{run}.sample_{sample}.dedup.bam"
    output: "bams/run_{run}.sample_{sample}.target.intervals"
    shell:
        "Java -Xmx8g -jar gatk -T RealignerTargetCreator -R reference -I {input} -L {enrichment} -o {output}"
        
rule realign_indels2:
    input: 'bams/run_{run}.sample_{sample}.dedup.bam'
    output: 'bams/run_{run}.sample_{sample}.realigned.bam'
    shell:
        "JAVA -Xmx8g -jar gatk -T IndelRealigner -R {reference} -I {input} -L {enrichment} -targetintervals {input}.target.intervals -known {BUNDLE}/Mills_and_1000G_gold_standard.indels.b37.vcf -o {output}"

rule remove:
    input: rules.realign_indels.output
    output: ".status"
    shell:
        "mv *.intervals ../logs"
        "mkdir dedupped"
        "mv *.dedup.* dedupped"
        "touch .status"

rule recal_tables:
    input: "bams/run_{run}.sample_{sample}.realigned.bam"
    output: "bams/run_{run}.sample_{sample}.recal.table"
    shell:
        "Java -Xmx6g -jar gatk -T BaseRecalibrator -R {reference} -I {input} -knownSites {BUNDLE}/{DBSNP}_138.b37.vcf -knownSites {BUNDLE}/Mills_and_1000G_gold_standard.indels.b37.vcf -L {enrichment} -o {output}"

rule baserecall:
    input: #"bams/run_{run}.sample_{sample}.realigned.bam"
    output: "bams/run_{run}.sample_{sample}.recal.bam"
    shell:
        "Java -Xmx6g -jar gatk -T PrintReads -R {reference} -I {input} -BQSR bams/run_{run}.sample_{sample}.recal.table -o {output}"

rule mv:
    input: rules.recal_tables.output
    output: ".status"
    shell:
        "mv {input} realigned/"

rule run_gatkhs:
    input: rules.baserecall.input
    output: "vcfs/run_{run}.sample_{sample}.g.vcf"
    shell:
        "Java -Xmx6g -jar gatk -T HaplotypeCaller -R {reference} -L {enrichment} -I {input} -o {output} -ERC GVCF"

rule mk_gvcfs:
    input: "vcfs/run_{run}.sample_{sample}.g.vcf"
    output: "gvcfs.list"
    shell:
        "mkdir gvcfs "
        "mv {input} gvcfs"
        #"ls {PWD}/gvcfs/{input} | cat {paddingcohort} - > {output}

rule genotype_gvcfs:
    input: rules.mk_gvcfs.output
    output: "run_{run}.sample_{sample}.raw.vcf"
    shell:
        "Java -Xmx64g -jar gatk -T GenotypeGVCFs -R {reference} -V {input} -o {output}"

rule modelling_snps:
    input: rules.genotype_gvcfs.output
    output: "run_{run}.sample_{sample}.snp.recal"
    shell:
        "Java -Xmx64g -jar gatk -T VariantRecalibrator "
            "-R {reference} "
            "-input {input} "
            "-resourse:hapmap,known=false,training=true,truth=true,prior=15.0 {BUNDLE}/hapmap_3.3.b37.vcf "
            "-resourse:omni, known=false,training=true,truth=true,prior=12.0 {BUNDLE}/1000G_omni2.5.b37.vcf "
            "-resourse:mills, known=false,training=true,truth=true,prior=12.0 {BUNDLE}/Mills_and_1000G_gold_standard.indels.b37.vcf "
            "-resourse:{DBSNP}, known=true,training=false,truth=false,prior=2.0 {BUNDLE}/{DBSNP}_138.b37.vcf "
            "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff "
            "-mode SNP "
            "-recalFile {output} "
            "-tranchesFile run_{run}.sample_{sample}.snp.tranches "
            "-rscriptFile run_{run}.sample_{sample}.snp.plots.R"

rule applying_snpfilter:
    input: rules.genotype_gvcfs.output
    output: "run_{run}.sample_{sample}.snp.recal.vcf"
    shell:
        "Java -Xmx64g -jar -gatk -T ApplyRecalibration "
            "-R {reference} "
            "-input {input} "
            "-mode SNP "
            "-recalFile{ run_{run}.sample_{sample}.snp.recal }"
            "-tranchesFile run_{run}.sample_{sample}.snp.tranches "
            "-o  {output} "
            "-ts_filter_level {SNPVQSLOD}"

rule modelling_indels:
    input: rules. applying_snpfilter.output
    output: "run_{run}.sample_{sample}.indel.recal"
    shell:
        "Java -Xmx64g -jar gatk -T VariantRecalibrator "
            "-R {reference} "
            "-input {input} "
            "--maxGaussians 4 "
            "-resourse:mills, known=false,training=true,truth=true,prior=12.0 {BUNDLE}/Mills_and_1000G_gold_standard.indels.b37.vcf "
            "-resourse:{DBSNP}, known=true,training=false,truth=false,prior=2.0 {BUNDLE}/{DBSNP}_138.b37.vcf "
            "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff "
            "-mode INDEL "
            "-recalFile {output} "
            "-tranchesFile run_{run}.sample_{sample}.indel.tranches "
            "-rscriptFile run_{run}.sample_{sample}.indel.plots.R"

rule applying_indelfilter:
    input: rules. applying_snpfilter.output
    output: "run_{run}.sample_{sample}.recal.vcf"
    shell:
        "Java -Xmx64g -jar -gatk -T ApplyRecalibration "
            "-R {reference} "
            "-input {input} "
            "-mode INDEL "
            "-recalFile{ run_{run}.sample_{sample}.indel.recal }"
            "-tranchesFile run_{run}.sample_{sample}.indel.tranches "
            "-o  {output} "
            "-ts_filter_level {INDELVQSLOD}"

rule select_variants:
    input: rules.applying_indelfilter.output
    output: "run_{run}.sample_{sample}.retained.vcf"
    shell:
        "Java -Xmx64g -jar gatk -T SelectVariants -R {reference} "
            "-V {input} "
            "-se 'ls gvcfs/*.g.vcf | grep -oP 'sample_{sample}\K[^\.]+' - | sed 's/^/S' | paste -s -d '|'' "
            "-select 'vc.isNotFiltered()' --excludeNonVariants "
            "-o {output}"

rule mk_vqsr:
    input: rules.applying_indelfilter.output
    output: ".status"
    shell:
        "mkdir vqsr"
        "mv *.recal* vqsr"
        "mv *tranches* vqsr"
        "mv *plots.* vqsr"
        "mv *.raw* vqsr"

rule calculate_genotypeposteriors:
    input: rules.select_variants.output
    output: "run_{run}.sample_{sample}.GR.vcf"
    shell:
        "Java -Xmx64g -jar gatk -T CalculateGenotypePosteriors -R {reference} -V {input} "
            "--supporting {BUNDLE}/1000G_phase3_v4_20130502.sites.vcf -o {output}"

rule variant_filtration:
    input: rules.calculate_genotypeposteriors.output
    output: "run_{run}.sample_{sample}.GF.vcf"
    shell:
        "Java -Xmx64g -jar gatk -T VariantFiltration -R {reference} -V {input} -G_filter 'CQ < 20.0' -G_filterName lowCQ -o {output}"

rule mk_gt_refine:
    input: rules. calculate_genotypeposteriors.output
    output: ".status"
    shell:
        "mkdir gt_refine"
        "mv *.retained.* gt_refine"
        "mv *.GR.* gt_refine"

rule change_gf:
    input: rules.variant_filtration.output
    output: "run_{run}.sample_{sample}.GF1.vcf"
    shell:
        "sed -I 's/MT\t/M\t/' {input} > {output}"

rule breakmultitool:
    input: rules.change_gf.output
    output: "run_{run}.sample_{sample}.S1.vcf"
    shell:
        "$BREAKMULTITOOL {input} > {output}" 

rule DBSNP:
    input: rules.breakmultitool.output
    output: "run_{run}.sample_{sample}.S2.vcf"
    shell:
        "Java -Xmx64g -jar ${SNPSIFT} annotate -v -noInfo {DBSNP} {input} > {output}"

rule CLINVAR:
    input: rules.DBSNP.output
    output: "run_{run}.sample_{sample}.S3.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} annotate -v -noId -info PM,OM,CLNSIG,CLNDBN,CLNREVSTAT {CLINVAR} {input} > {output}"

rule KGENOMES:
    input: rules.CLINVAR.output
    output: "run_{run}.sample_{sample}.S4.vcf"
    shell:
        "Java -Xmx64g -jar ${SNPSIFT} annotate -v -noId -info X1000Gp3_AF {KGENOMES} {input} > {output}"

rule exac:
    input: rules.KGENOMES.output
    output: "run_{run}.sample_{sample}.S5.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} annotate -v -noId -info ExAC_AF EXAC {input} > {output}"

rule esp:
    input: rules.exac.output
    output: "run_{run}.sample_{sample}.S6.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} annotate -v -noId -info ESP6500_AF $ESP {input} > {output}"

rule LOCALAFDB:
    input: rules.esp.output
    output: "run_{run}.sample_{sample}.S7.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} annotate -v -noId -info LOCALAF {LOCALAFDB} {input} > {output}"

rule GWASCAT:
    input: rules.LOCALAFDB.output
    output: "run_{run}.sample_{sample}.S8.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} qwasCat -db {GWASCAT} {input} > {output}"

rule positive_train_site:
    input: rules.GWASCAT.output
    output: "run_{run}.sample_{sample}.S9.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} rmInfo {input} 'AC' 'AN' 'BaseQRankSum "
            "'ClippingRankSum' 'DP' 'ExcessHet' 'InbreedingCoeff' "
            "'MLEAC' 'MLEAF' 'MQRankSum' 'PG' 'ReadPosRankSum' 'SOR' "
            "'culprit' 'POSITIVE_TRAIN_SITE' > {output}"

rule DBNSFP:
    input: rules.positive_train_site.output
    output: "run_{run}.sample_{sample}.S10.vcf"
    shell:
        "Java -Xmx64g -jar {SNPSIFT} {DBNSFP} -db {DBNSFP} -f PROVEAN_pred,SIFT_pred,Polyphen2_HVAR_pred -v {input} > {output}"

rule SNPEFF:
    input: rules.DBNSFP.output
    output: "run_{run}.sample_{sample}.Final.vcf"
    shell:
        " Java -Xmx64g -jar {SNPEFF} -v hg19 {input} > {output}"

rule mk_annotation:
    input: rules.change_gf.output
    output: ".status"
    shell:
        "mkdir annotation"
        "mv *.S*.vcf* annotation/"
        "mv run_{run}.sample_{sample}.GF.vcf* annotation/"

rule splitting_vcf:
    input: rules.SNPEFF.output
    output: "run_{run}.sample_{sample}.{sample}.vcf"
    shell:
        "ls -lh ./gvcfs/*.g.vcf | grep -oP '[^\/]+sample_{sample}[^\.]+' | perl -pe 's|.*?sample_{sample}|S|g' | sort -u > samples.list"
        "for sample in $(cat samples.list)"
        "do" 
            "Java -Xmx4g -jar gatk -T SelectVariants -R {reference} -V {input} -sn $sample --excludeNonVariants -o {output}"
        "done"

