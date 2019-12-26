# a simple script to make sample GVCFs in a map file in python

import re

with open("/gatk/output_folder/vc_proj/all_outputs/joint_vcf/gvcfs.list") as file:

# This regex expression will get everything after the last backslash
# ([^\/]+$)

    for line in file:
        #tag = re.search(r'([^\/]+$)', line.strip('\n')).group(0).strip('.dedup.bam')
        print("-V {}".format(line))
