import sys
import pandas as pd
csv_file=open(sys.argv[1])
sample_name=sys.argv[2]

def read_df(csv_file):
	df_rc = pd.read_csv(csv_file,
				sep='\t',
				comment='#',
				names=['#CHROM','POS','REF','col4','col5','geno1','geno2'])
	
	return df_rc

def count_all_reads(chrompos):
	global rc_all
	this_df=rc_all[rc_all['chrompos']==chrompos]
	return this_df['reads'].sum()
		
df_rc = read_df(csv_file)
df_rc = df_rc.drop([0])
df_rc = df_rc.drop(columns=['col4','col5'])

df_rc_geno2 = df_rc[df_rc["geno2"].notnull()]
df_rc_geno2 = df_rc_geno2.drop(columns=['geno1'])

#df_rc_geno2.columns=['chrom','pos','ref','geno1']
df_rc_geno2.rename(columns={'geno2': 'geno1'}, inplace=True)
df_rc_geno2.head()


df_rc_geno1= df_rc.drop(columns = ['geno2'])

rc_all = pd.concat([df_rc_geno1, df_rc_geno2])
#quit()
rc_all = rc_all.sort_values(by=['#CHROM','POS'])
rc_all['chrompos']=rc_all['#CHROM'].map(str)+"_"+rc_all['POS'].map(str)

rc_all['ALT']=rc_all['geno1'].str.split(":").str[0]
rc_all['reads']=rc_all['geno1'].str.split(":").str[1]
rc_all['strands']=rc_all['geno1'].str.split(":").str[2]
rc_all['avg_qual']=rc_all['geno1'].str.split(":").str[3]
rc_all['map_qual']=rc_all['geno1'].str.split(":").str[4]
rc_all['plus_reads']=rc_all['geno1'].str.split(":").str[5]
rc_all['minus_reads']=rc_all['geno1'].str.split(":").str[6]
rc_all['DP']=rc_all['chrompos'].apply(count_all_reads)

rc_all= rc_all[rc_all['REF']!=rc_all['ALT']]

rc_all['ID']='.'
rc_all['FILTER']='NaN'
rc_all['INFO']='NaN'
rc_all['QUAL']='NaN'
rc_all['FORMAT']='GT:AD:DP:GQ:PL'
rc_all['DP']=rc_all['chrompos'].apply(count_all_reads)
rc_all['SAMPLE']= "0/1:"+rc_all['reads'].map(str)+':'+rc_all['DP'].map(str)+':NaN:NaN'

rc_vcf=rc_all[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']]

rc_vcf.rename(columns={'SAMPLE': sample_name}, inplace=True)
rc_vcf.to_csv(sample_name+".w_alt_alleles.vcf", sep='\t', index=False)

