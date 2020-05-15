#!/usr/bin/env python3.6

import re
import sys

vcf = sys.argv[1]

csq_order = ["transcript_ablation",
"splice_acceptor_variant",
"splice_donor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost",
"start_lost",  # new in v81
"initiator_codon_variant",  # deprecated
"transcript_amplification",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"protein_altering_variant",  # new in v79
"splice_region_variant",
"incomplete_terminal_codon_variant",
"stop_retained_variant",
"synonymous_variant",
"coding_sequence_variant",
"mature_miRNA_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"non_coding_transcript_exon_variant",
"non_coding_exon_variant",  # deprecated
"intron_variant",
"NMD_transcript_variant",
"non_coding_transcript_variant",
"nc_transcript_variant",  # deprecated
"upstream_gene_variant",
"downstream_gene_variant",
"start_retained_variant",
"TFBS_ablation",
"TFBS_amplification",
"TF_binding_site_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
"feature_elongation",
"regulatory_region_variant",
"feature_truncation",
"intergenic_variant",
""]


with open(vcf, 'r') as vf_handle:
    for line in vf_handle:
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('##INFO=<ID=CSQ'):
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
            elif line.startswith('#CHROM'):
                print('##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomad af global">')
                print('##INFO=<ID=gnomAD_AFR_AF,Number=1,Type=Float,Description="gnomad af african">')
                print('##INFO=<ID=gnomAD_AMR_AF,Number=1,Type=Float,Description="gnomad af american">')
                print('##INFO=<ID=gnomAD_EAS_AF,Number=1,Type=Float,Description="gnomad af east asia">')
                print('##INFO=<ID=gnomAD_SAS_AF,Number=1,Type=Float,Description="gnomad af south asia">')
                print('##INFO=<ID=gnomAD_NFE_AF,Number=1,Type=Float,Description="gnomad af non fins">')
                print('##INFO=<ID=gnomAD_FIN_AF,Number=1,Type=Float,Description="gnomad af finnish">')
                print('##INFO=<ID=worst_csq,Number=1,Type=String,Description="worst variant consequence">')
                print('##INFO=<ID=gene_name,Number=1,Type=String,Description="gnomad af global">')
                print('##INFO=<ID=clinvar_significane,Number=1,Type=String,Description="Clinsig from clinvar">')
                print('##INFO=<ID=nonref_AF,Number=1,Type=Float,Description="total frequency of alt alleles">')
            print(line)
        else:
#            print(vep_field_names)
            nonref_af = str(sum([float(x) for x in re.findall('AF=([^;]+)', line)[0].split(',')]))
#            print(re.findall('CSQ=([^;^\t]+)', line)[0].split(','))
            annotations = [dict(zip(vep_field_names, x.split('|'))) for x in re.findall('CSQ=([^;^\t]+)', line)[0].split(',')]
            csqs = [x['Consequence'] for x in annotations]
            csq_order_dict = {csq:i for i,csq in enumerate(csq_order)}
            csqs_clean = []
            for csq in csqs:
                for elem in csq.split('&'):
                    csqs_clean.append(elem)
#            print(csqs_clean)
            worst_csq = csq_order[min([csq_order_dict[csq] for csq in csqs_clean])]
            gene_name = annotations[0]['SYMBOL'] if annotations[0]['SYMBOL'] != '' else '.'
            clnsig = annotations[0]['CLIN_SIG'] if annotations[0]['CLIN_SIG'] != '' else 'no_value_provided'
            gnomad_af = str(sum(set([float(x['gnomAD_AF']) for x in annotations if x['gnomAD_AF'] != ''])))
            gnomad_afr_af = str(sum(set([float(x['gnomAD_AFR_AF']) for x in annotations if x['gnomAD_AFR_AF'] != ''])))
            gnomad_amr_af = str(sum(set([float(x['gnomAD_AMR_AF']) for x in annotations if x['gnomAD_AMR_AF'] != ''])))
            gnomad_eas_af = str(sum(set([float(x['gnomAD_EAS_AF']) for x in annotations if x['gnomAD_EAS_AF'] != ''])))
            gnomad_sas_af = str(sum(set([float(x['gnomAD_SAS_AF']) for x in annotations if x['gnomAD_SAS_AF'] != ''])))
            gnomad_nfe_af = str(sum(set([float(x['gnomAD_NFE_AF']) for x in annotations if x['gnomAD_NFE_AF'] != ''])))
            gnomad_fin_af = str(sum(set([float(x['gnomAD_FIN_AF']) for x in annotations if x['gnomAD_FIN_AF'] != ''])))
            content = line.strip().split('\t')
            content[7] = f'{content[7]};nonref_AF={nonref_af};worst_csq={worst_csq};clinvar_significance={clnsig};gene_name={gene_name};gnomAD_AF={gnomad_af};gnomAD_AFR_
AF={gnomad_afr_af};gnomAD_AMR_AF={gnomad_amr_af};gnomAD_EAS_AF={gnomad_eas_af};gnomAD_SAS_AF={gnomad_sas_af};gnomAD_NFE_AF={gnomad_nfe_af};gnomAD_FIN_AF={gnomad_fin_af}'
            print('\t'.join(content))




