import sys
import os
import argparse
import logging
import plot
import pandas as pd
import scipy.stats as stats
import seaborn as sns

# set basic pandas options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# parser command line
parser = argparse.ArgumentParser(description="Analysis correlation and different biomarkers for the PET/LIT study.")
parser.add_argument('-in', '--in_file', required=True, help='Path to in-file .')
parser.add_argument('-out', '--out_folder', required=True, help='Path to out-folder.')
parser.add_argument('-log', '--log', required=False, help='Path to log file.')
args = parser.parse_args()

# initiate logging
if args.log is not None:
    logging.basicConfig(filename=args.log, encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
else:
    logging.basicConfig(encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

# logging
logging.info('Running {}'.format(os.path.basename(__file__)))
logging.info('Arguments loaded with ArgumentParser:')
for arg,value in sorted(vars(args).items()):
    logging.info(' Argument {}: {}'.format(arg,value))

analysis_folder = args.out_folder
out_folder = analysis_folder
patients = pd.read_excel(args.in_file, sheet_name='Patients')
variants = pd.read_excel(args.in_file, sheet_name='LB-Variants')
diff_pet_single = pd.read_excel(args.in_file, sheet_name='LB_PET-Single')
diff_pet_pair = pd.read_excel(args.in_file, sheet_name='LB_PET-Pair')
diff_s100 = pd.read_excel(args.in_file, sheet_name='LB_S100')
diff_ldh = pd.read_excel(args.in_file, sheet_name='LB_LDH')

# set file paths
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

# set seaborn style
sns.set_theme(style='white')

# set basic parameters
exclude_patients_pet = []
exclude_patients_bio = []
max_offset_pet = 21
max_offset_bio = 7
minimum_tlg = -1
log = True

# dictionary to collect figure data for export
source_data = {}

######
# single time points
########
logging.info('Single PET-analysis')
# set filters for PET single
diff_pet_single['FILTER'] = diff_pet_single['FILTER'].fillna('')
diff_pet_single.loc[diff_pet_single['PATIENT-ID'].isin(exclude_patients_pet), 'FILTER'] += 'PATIENT-EXCLUDED'
diff_pet_single.loc[pd.isnull(diff_pet_single['PETCT-OFFSET']), 'FILTER'] += 'NO-PETCT'

tmp_single_final = diff_pet_single[diff_pet_single['FILTER']==''].copy()
if max_offset_pet>=0:
    diff_pet_single.loc[abs(diff_pet_single['PETCT-OFFSET'])>max_offset_pet,'FILTER'] += "MORE-THAN-" + str(max_offset_pet) + "-DAYS"
if minimum_tlg>=0:
    diff_pet_single.loc[diff_pet_single['PETCT-TLG']<minimum_tlg,'FILTER'] += "TLG-LESS-THAN-" + str(minimum_tlg)

# reduce dataset to palliative patients only
tmp_diff1 = diff_pet_single[(diff_pet_single['FILTER']=='') & (diff_pet_single['TREATMENT']=='Pall')].copy()

# correlate of TLG, MTV and cfDNA
column1 = 'PETCT-TLG'
column2 = 'PLASMA-MEAN-MULTI-AF'
column3 = 'PETCT-MTV'

# normaltest, H0 sample came from normal distribution, p-value < 0.05 H0 has to be rejected
p = stats.normaltest(tmp_diff1[column1])[1]
logging.info(' Normaltest {}: p={:.2f}.'.format(column1,p))
p = stats.normaltest(tmp_diff1[column2])[1]
logging.info(' Normaltest {}: p={:.2f}.'.format(column2,p))
p = stats.normaltest(tmp_diff1[column3])[1]
logging.info(' Normaltest {}: p={:.2f}.'.format(column3,p))

# calculate Spearman's correlation treatment start
rho1, pval1 = stats.spearmanr(tmp_diff1.loc[(tmp_diff1['PLASMA-T0-TREATMENT-START']=='yes'),column1],tmp_diff1.loc[(tmp_diff1['PLASMA-T0-TREATMENT-START']=='yes'),column2])
logging.info(' Number of LBs at T0: {}, Rho {}/{}: {}'.format(len(tmp_diff1[tmp_diff1['PLASMA-T0-TREATMENT-START']=='yes']),column1,column2,rho1))
rho2, pval2 = stats.spearmanr(tmp_diff1.loc[(tmp_diff1['PLASMA-T0-TREATMENT-START']=='yes'), column3],tmp_diff1.loc[(tmp_diff1['PLASMA-T0-TREATMENT-START']=='yes'),column2])
logging.info(' Number of LBs at T0: {}, Rho {}/{}: {}'.format(len(tmp_diff1[tmp_diff1['PLASMA-T0-TREATMENT-START']=='yes']),column3,column2,rho2))

# calculate Spearman's correlation all samples and scatter plot
n = len(tmp_diff1['PATIENT-ID'].unique())
logging.info(' Number of patients: {}'.format(n))
rho1, pval1 = stats.spearmanr(tmp_diff1[column1],tmp_diff1[column2])
logging.info(' Number of LBs all Ts: {}, Rho {}/{}: {}, p-value: {}'.format(len(tmp_diff1),column1,column2,rho1,pval1))
plot.plot_correlation(x=tmp_diff1[column1], y=tmp_diff1[column2], title='Correlation of liquid biopsies and PET/CT', xlabel='TLG PET/CT', ylabel='Multi-AF',rho=rho1, p=pval1, out=(out_folder + 'Figure_Correlation-TLG.svg'), log=log)
rho2, pval2 = stats.spearmanr(tmp_diff1[column3],tmp_diff1[column2])
logging.info(' Number of LBs all Ts: {}, Rho {}/{}: {}, p-value: {}'.format(len(tmp_diff1),column3,column2,rho2,pval2))
plot.plot_correlation(x=tmp_diff1[column3], y=tmp_diff1[column2], title='Correlation of liquid biopsies and PET/CT', xlabel='MTV PET/CT', ylabel='Multi-AF',rho=rho2, p=pval2, out=(out_folder + 'Figure_Correlation-MTV.svg'), log=log)

source_data.update({'Figure 3 C':tmp_diff1[['PATIENT-ID',column3,column2]]})
source_data.update({'Supplementary Fig. 5':tmp_diff1[['PATIENT-ID',column1,column2]]})

######
# S100
########
# set filters
diff_s100['FILTER'] = diff_s100['FILTER'].fillna("")
diff_s100.loc[diff_s100['PATIENT-ID'].isin(exclude_patients_bio), 'FILTER'] += 'PATIENT-EXCLUDED'
if max_offset_bio>=0:
    diff_s100.loc[abs(diff_s100['S100-OFFSET'])>max_offset_bio,'FILTER'] += "MORE-THAN-" + str(max_offset_bio) + "-DAYS"

#
logging.info('S100')
offset_s100 = 0.1

# reduce dataset to palliative patients only
tmp_diff1 = diff_s100[ (diff_s100['FILTER']=='') & (diff_s100['S100-OFFSET'].notna()) & (diff_s100['TREATMENT']=='Pall')].copy()

#
column1 = 'S100-CONC'
column2 = 'PLASMA-MEAN-MULTI-AF'

p = stats.normaltest(tmp_diff1[column1])[1]
logging.info(' Normaltest {}: p={:.2f}.'.format(column1,p))
p = stats.normaltest(tmp_diff1[column2])[1]
logging.info(' Normaltest {}: p={:.2f}.'.format(column2,p))

n = len(tmp_diff1['PATIENT-ID'].unique())
logging.info(' Number of patients: {}'.format(n))
rho1, pval1 = stats.spearmanr(tmp_diff1[column1],tmp_diff1[column2])
logging.info(' S100 {} Spearman\'s rho: {:.2f}, test p-value: {}'.format(len(tmp_diff1),rho1,pval1))
plot.plot_correlation(x=tmp_diff1[column1], y=tmp_diff1[column2], title='Correlation of liquid biopsies and S100', xlabel='S100', ylabel='Multi-AF',rho=rho1, p=pval1, out=(out_folder + 'Figure_Correlation-S100.svg'), log=log)

source_data.update({'Figure 3 A':tmp_diff1[['PATIENT-ID',column1,column2]]})

########
# LDH
###
logging.info('LDH')
diff_ldh['FILTER'] = diff_ldh['FILTER'].fillna("")
diff_ldh.loc[diff_ldh['PATIENT-ID'].isin(exclude_patients_bio), 'FILTER'] += 'PATIENT-EXCLUDED'
if max_offset_bio>=0:
    diff_ldh.loc[abs(diff_ldh['LDH-OFFSET'])>max_offset_bio,'FILTER'] += "MORE-THAN-" + str(max_offset_bio) + "-DAYS"

#
offset_ldh = 250
tmp_diff1 = diff_ldh[ (diff_ldh['FILTER']=='') & (diff_ldh['LDH-OFFSET'].notna()) & (diff_ldh['TREATMENT']=='Pall')].copy()

#
column1 = 'LDH-CONC'
column2 = 'PLASMA-MEAN-MULTI-AF'

n = len(tmp_diff1['PATIENT-ID'].unique())
logging.info(' Number of patients: {}'.format(n))
p = stats.normaltest(tmp_diff1[column1])
logging.info(' Normaltest {}: {}'.format(column1,p[1]))
p = stats.normaltest(tmp_diff1[column2])
logging.info(' Normaltest {}: {}'.format(column2,p[1]))

rho1, pval1 = stats.spearmanr(tmp_diff1[column1],tmp_diff1[column2])
logging.info(' LDH {} Spearman\'s rho: {:.2f}, p-value: {}'.format(len(tmp_diff1),rho1,pval1))
plot.plot_correlation(x=tmp_diff1[column1], y=tmp_diff1[column2], title='Correlation of liquid biopsies and LDH', xlabel='LDH', ylabel='Multi-AF',rho=rho1, p=pval1, out=(out_folder + 'Figure_Correlation-LDH.svg'), log=log)

source_data.update({'Figure 3 B':tmp_diff1[['PATIENT-ID',column1,column2]]})

# write excel list for figure data
logging.info('Writing Excel file...')
with pd.ExcelWriter(path=out_folder+'Data_Correlation.xlsx') as writer:
    workbook = writer.book
    for key in source_data.keys():
        sheet_name = key
        worksheet = workbook.add_worksheet(key)
        worksheet.write(0,0,key)
        writer.sheets[key] = worksheet
        source_data[key].to_excel(writer, sheet_name=sheet_name, index=False, startrow=2, startcol=2)

logging.info('Done.')