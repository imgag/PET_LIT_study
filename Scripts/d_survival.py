import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import plot
import lifelines
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import subprocess as sp
from lifelines.statistics import logrank_test
from scipy import stats

# set basic pandas options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# parser command line
parser = argparse.ArgumentParser(description="Analysis survival and different biomarkers for the PET/LIT study.")
parser.add_argument('-in', '--in_file', required=True, help='Path to in-file .')
parser.add_argument('-out', '--out_folder', required=True, help='Path to out-folder.')
parser.add_argument('-log', '--log', required=False, help='Path to log file.')
args = parser.parse_args()

# initiate logging
if args.log is not None:
    logging.basicConfig(filename=args.log, encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
else:
    logging.basicConfig(encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

#
logging.info('Running {}'.format(os.path.basename(__file__)))
logging.info('Arguments loaded with ArgumentParser:')
for arg,value in sorted(vars(args).items()):
    logging.info(' Argument {}: {}'.format(arg,value))

source_data = {}

# set file paths
out_folder = args.out_folder
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

# extract PFS and OS for all patients and add them to the current DF
patients = pd.read_excel(args.in_file, sheet_name='Patients')
samples = pd.read_excel(args.in_file, sheet_name='LB-Samples')
timepoints = pd.read_excel(args.in_file, sheet_name='Timepoints')
timepoints[['PFS-DURATION','PFS-EVENT','OS-DURATION','OS-EVENT']] = pd.DataFrame([[None,None,None,None]], index=timepoints.index).astype(float)
for i,r in patients.iterrows():
    timepoints.loc[timepoints['PATIENT-ID']==r['PATIENT-ID'], ['PFS-DURATION','PFS-EVENT','OS-DURATION','OS-EVENT']] = [r['PFS-DURATION'],r['PFS-EVENT'],r['OS-DURATION'],r['OS-EVENT']]

if (timepoints['PFS-DURATION'].isnull().values.any()) or (timepoints['PFS-EVENT'].isnull().values.any()):
    logging.error('Not all patients with PFS.')
    sys.exit(0)
if (timepoints['OS-DURATION'].isnull().values.any()) or (timepoints['OS-EVENT'].isnull().values.any()):
    logging.error('Not all patients with OS.')
    sys.exit(0)

logging.info('PFS and OS for all patients available.')

#######
## CTDNA
## Kaplan-Meier and logrank
#####

# change stdout
orig_stdout = sys.stdout
f = open(out_folder + 'd_Survival-Logrank.txt','w')
sys.stdout = f

##
col = 'CTDNA-T0-MRD-POS'
timepoints[col] = np.nan
timepoints.loc[timepoints['CTDNA-T0-MRD']<0.05,col] = 1
timepoints.loc[timepoints['CTDNA-T0-MRD']>=0.05,col] = 0
tmp1 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==1)]
tmp2 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
print('p = {:.10f}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-1.svg',title='Pall CTDNA T0 PFS (n={})'.format((len(tmp1)+len(tmp2))),label1='Tumour detected',label2='No tumour detected')

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-2.svg',title='Pall CTDNA T0 OS (n={})'.format((len(tmp1)+len(tmp2))),label1='Tumour detected',label2='No tumour detected')

##
col = 'CTDNA-T1-MRD-POS'
timepoints[col] = np.nan
timepoints.loc[timepoints['CTDNA-T1-MRD']<0.05,col] = 1
timepoints.loc[timepoints['CTDNA-T1-MRD']>=0.05,col] = 0
tmp1 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==1)]
tmp2 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-3.svg',title='Pall CTDNA T1 PFS (n={})'.format((len(tmp1)+len(tmp2))),label1='Tumour detected',label2='No tumour detected')

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-4.svg',title='Pall CTDNA T1 OS (n={})'.format((len(tmp1)+len(tmp2))),label1='Tumour detected',label2='No tumour detected')

##
col = 'CTDNA-T2-MRD-POS'
timepoints[col] = np.nan
timepoints.loc[timepoints['CTDNA-T2-MRD']<0.05,col] = 1
timepoints.loc[timepoints['CTDNA-T2-MRD']>=0.05,col] = 0
tmp1 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==1)]
tmp2 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-5.svg',title='Pall CTDNA T2 PFS (n={})'.format((len(tmp1)+len(tmp2))),label1='Tumour detected',label2='No tumour detected')

(group1_d, group1_e, group2_d, group2_e) = (tmp1['OS-DURATION'].tolist(), tmp1['OS-EVENT'].tolist(), tmp2['OS-DURATION'].tolist(), tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e, event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-6.svg',title='Pall CTDNA T2 OS (n={})'.format((len(tmp1)+len(tmp2))),label1='Tumour detected',label2='No tumour detected')

## Kaplan-Meier CTDNA T0-T1
col = 'CTDNA-T0-T1-ASC'
timepoints[col] = np.nan
for i,r in timepoints.iterrows():
    if (not pd.isnull(r['CTDNA-T0'])) and (not pd.isnull(r['CTDNA-T1'])):
        delta = (r['CTDNA-T1-AF'] - r['CTDNA-T0-AF'])
        if delta>=0:
            timepoints.loc[i,col] = 1
        elif delta<0:
            timepoints.loc[i,col] = 0
tmp1 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==1)]
tmp2 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-7.svg',title='Change of ctDNA level T0-T1 (PFS)\nin palliative combined ICI patients'.format((len(tmp1)+len(tmp2))),label1='AF asc',label2='AF desc')

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-8.svg',title='Change of ctDNA level T0-T1 (OS)\nin palliative combined ICI patients'.format((len(tmp1)+len(tmp2))),label1='AF asc',label2='AF desc')

source_data.update({'Figure 4 A+B':timepoints.loc[(timepoints['TREATMENT']=='Pall') & (timepoints[col].isin([0,1])),['PATIENT-ID','TREATMENT','CTDNA-T0-T1-ASC','PFS-DURATION','PFS-EVENT','OS-DURATION','OS-EVENT']]})

## Kaplan-Meier CTDNA T1-T2
col = 'CTDNA-T1-T2-ASC'
timepoints[col] = np.nan
for i,r in timepoints.iterrows():
    if (not pd.isnull(r['CTDNA-T1'])) and (not pd.isnull(r['CTDNA-T2'])):
        delta = (r['CTDNA-T2-AF'] - r['CTDNA-T1-AF'])
        if delta>=0:
            timepoints.loc[i,col] = 1
        elif delta<0:
            timepoints.loc[i,col] = 0
tmp1 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==1)]
tmp2 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-9.svg',title='Change of ctDNA level T1-T2 (PFS)\nin palliative combined ICI patients'.format((len(tmp1)+len(tmp2))),label1='AF asc',label2='AF desc')

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print(' p-value = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-10.svg',title='Change of ctDNA level T1-T2 (OS)\nin palliative combined ICI patients'.format((len(tmp1)+len(tmp2))),label1='AF asc',label2='AF desc')

source_data.update({'Figure 4 C+D':timepoints.loc[(timepoints['TREATMENT']=='Pall') & (timepoints[col].isin([0,1])),['PATIENT-ID','TREATMENT','CTDNA-T1-T2-ASC','PFS-DURATION','PFS-EVENT','OS-DURATION','OS-EVENT']]})

## TMB
col = 'TMB-high'
tmb_cutoff = 7
patients[col] = np.nan
patients.loc[patients['TMB']>tmb_cutoff,col] = 1
patients.loc[patients['TMB']<=tmb_cutoff,col] = 0
# n.b. one patient is external and has empty TMB
tmp1 = patients[pd.isnull(patients['FILTER']) & (patients[col]==1)]
tmp2 = patients[pd.isnull(patients['FILTER']) & (patients[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-11.svg',title='Pall TMB PFS (n={},cut-off={})'.format((len(tmp1)+len(tmp2)),tmb_cutoff),label1='TMB-high',label2='TMB-low')

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-12.svg',title='Pall TMB OS (n={},cut-off={})'.format((len(tmp1)+len(tmp2)),tmb_cutoff),label1='TMB-high',label2='TMB-low')

source_data.update({'Supplementary Fig. 7 A': patients.loc[pd.isnull(patients['FILTER']) & patients[col].isin([0,1]),['PATIENT-ID','TREATMENT',col,'PFS-DURATION','PFS-EVENT']]})

## LDH
# LDH increased > 250 U/l
cutoff = 250
col = 'LDH-T0-POS'
timepoints[col] = np.nan
timepoints.loc[(timepoints['TREATMENT']=='Pall') & (timepoints['LDH-T0'].notnull()) & (timepoints['LDH-T0-VALUE']>cutoff),col] = 1
timepoints.loc[(timepoints['TREATMENT']=='Pall') & (timepoints['LDH-T0'].notnull()) & (timepoints['LDH-T0-VALUE']<=cutoff),col] = 0
tmp1 = timepoints[(timepoints[col]==1)]
tmp2 = timepoints[(timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-13.svg',title='Pall LDH T0 PFS (n={})'.format((len(tmp1)+len(tmp2))),label1='LDH >{}U/l'.format(cutoff),label2='LDH <={}U/l'.format(cutoff))

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-14.svg',title='Pall LDH T0 OS (n={})'.format((len(tmp1)+len(tmp2))),label1='LDH >{}U/l'.format(cutoff),label2='LDH <={}U/l'.format(cutoff))

source_data.update({'Supplementary Fig. 7 B': timepoints.loc[timepoints[col].isin([0,1]),['PATIENT-ID','TREATMENT',col,'PFS-DURATION','PFS-EVENT']]})

## PET/CT
col = 'PET-T0-T2-ASC'
timepoints[col] = np.nan
for i,r in timepoints.iterrows():
    if (not pd.isnull(r['PET-T0-MTV'])) and (not pd.isnull(r['PET-T2-MTV'])):
        delta = (r['PET-T2-MTV'] - r['PET-T0-MTV'])
        if delta>=0:
            timepoints.loc[i,col] = 1
        elif delta < 0:
            timepoints.loc[i,col] = 0
tmp1 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==1)]
tmp2 = timepoints[(timepoints['TREATMENT']=='Pall') & (timepoints[col]==0)]

(group1_d,group1_e,group2_d,group2_e) = (tmp1['PFS-DURATION'].tolist(),tmp1['PFS-EVENT'].tolist(),tmp2['PFS-DURATION'].tolist(),tmp2['PFS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} PFS'.format(col))
print('p = {}'.format(log.p_value))
print('p = {:.10f}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-15.svg',title='Pall PET T0-T2 PFS (n={})'.format((len(tmp1)+len(tmp2))),label1='MTV_asc',label2='MTV_desc')

(group1_d,group1_e,group2_d,group2_e) = (tmp1['OS-DURATION'].tolist(),tmp1['OS-EVENT'].tolist(),tmp2['OS-DURATION'].tolist(),tmp2['OS-EVENT'].tolist())
log = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e)
print('###')
print('{} OS'.format(col))
print('p = {}'.format(log.p_value))
log.print_summary()
print('')
plot.plot_kaplan_meier(group1_d,group1_e,group2_d,group2_e,out_folder + r'Figure_Survival-16.svg',title='Pall PET T0-T2 OS (n={})'.format((len(tmp1)+len(tmp2))),label1='MTV_asc',label2='MTV_desc')

sys.stdout = orig_stdout
f.close()

#######
## CTDNA
## Cox
#####

# change stdout
orig_stdout = sys.stdout
f = open(out_folder + 'd_Survival-CoxRegr.txt','w')
sys.stdout = f

logging.info('COX Regression Pall')
print('T0-T1 ASC PFS and OS')
tmp_timepoints = timepoints[(timepoints['TREATMENT']=='Pall')  & (timepoints['CTDNA-T0-T1-ASC'].notnull())]
cph = lifelines.CoxPHFitter()
tmp_cox = tmp_timepoints[['PFS-DURATION','PFS-EVENT','CTDNA-T0-T1-ASC']]
cph.fit(tmp_cox, duration_col='PFS-DURATION', event_col='PFS-EVENT')
print('###')
cph.print_summary()
print('')
tmp_cox = tmp_timepoints[['OS-DURATION','OS-EVENT','CTDNA-T0-T1-ASC']]
cph.fit(tmp_cox, duration_col='OS-DURATION', event_col='OS-EVENT')
print('###')
cph.print_summary()
print('')

print('T1-T2 ASC PFS and OS')
tmp_timepoints = timepoints[(timepoints['TREATMENT']=='Pall')  & (timepoints['CTDNA-T1-T2-ASC'].notnull())]
cph = lifelines.CoxPHFitter()
tmp_cox = tmp_timepoints[['PFS-DURATION','PFS-EVENT','CTDNA-T1-T2-ASC']]
cph.fit(tmp_cox, duration_col='PFS-DURATION', event_col='PFS-EVENT')
print('###')
cph.print_summary()
print('')
tmp_cox = tmp_timepoints[['OS-DURATION','OS-EVENT','CTDNA-T1-T2-ASC']]
cph.fit(tmp_cox, duration_col='OS-DURATION', event_col='OS-EVENT')
print('hazard ratio: {}'.format(cph.hazard_ratios_))
print('confidence intervals: {}'.format(cph.confidence_intervals_))
print('###')
cph.print_summary()
print('')
sys.stdout = orig_stdout
f.close()

#######
## CTDNA
## swimmer plot adjuvant
#####

count = 0
tmp_patients = patients[((patients['TREATMENT']=='Adj') & (pd.isnull(patients['FILTER'])))].sort_values(by=['PFS-EVENT','PATIENT-ID'], ascending=[True,False])
print(tmp_patients)

mm = 0.1 / 2.54
y_pos = np.arange(len(tmp_patients))
min_x = tmp_patients['OS-DURATION'].min()
max_x = tmp_patients['OS-DURATION'].max()
fig = plt.figure(figsize=(2*160*mm,2*60*mm))
ax = fig.add_subplot(111)
ax.tick_params(axis='x', labelrotation=90)

tmp_mask = ((samples['PATIENT-ID'].isin(tmp_patients['PATIENT-ID'])) & (samples['TREATMENT']=='Adj') & (pd.isnull(samples['FILTER'])))
tmp = {'Patients':tmp_patients[['PATIENT-ID','TREATMENT','PFS-EVENT','PFS-DURATION','OS-EVENT','OS-DURATION']],'Samples':samples.loc[tmp_mask,['PATIENT-ID','TREATMENT','TREATMENT-DAY','COUNT-VARIANTS-DETECTED']]}
source_data.update({'Figure 4 E':tmp})

for i1,p in tmp_patients.iterrows():
    relapse = '#D9AFD7'
    no_relapse = '#69addb'

    color = no_relapse
    if p['PFS-EVENT'] == 1:
        color = relapse
    ax.barh(count,p['OS-DURATION']+25, left=-25, color=color, align='center')

    if p['PFS-EVENT'] == 1:
        ax.plot(p['PFS-DURATION'], count, marker='s', markerfacecolor='None', markeredgecolor=(0, 0, 0, 1), markersize=6.5, markeredgewidth=1.5, linestyle='', label="Relapse")
    if p['OS-EVENT'] == 1:
        ax.plot(p['OS-DURATION']-8, count, marker='s', markerfacecolor=(0, 0, 0, 1), markeredgecolor=(0, 0, 0, 1), markersize=6.5, markeredgewidth=0, linestyle='', label='Death')
    if p['OS-EVENT'] == 0:
        ax.plot(p['OS-DURATION']+8, count, marker='>', markerfacecolor=(0, 0, 0, 1), markeredgecolor=(0, 0, 0, 1), markersize=6.5, markeredgewidth=1, linestyle='', label='Ongoing')

    tmp_samples = samples[((samples['PATIENT-ID']==p['PATIENT-ID']) & (samples['TREATMENT']=='Adj') & (pd.isnull(samples['FILTER'])))]
    for i2,s in tmp_samples.iterrows():
        if s['TREATMENT-DAY']<min_x:
            min_x = s['TREATMENT-DAY']

        if s['COUNT-VARIANTS-DETECTED'] < 3:
            ax.plot(s['TREATMENT-DAY'], count, marker='o', markerfacecolor='#07D400', markeredgecolor='White', markersize=6.0, markeredgewidth=0, linestyle='', label='ctDNA negative')
        else:
            ax.plot(s['TREATMENT-DAY'], count, marker='o', markerfacecolor='#E63A1D', markeredgecolor='White', markersize=6.0, markeredgewidth=0, linestyle='', label='ctDNA positive')

    count += 1

handles,labels = ax.get_legend_handles_labels()

# remove duplicates
tmp_handles = []
tmp_labels = []
for (h,l) in zip(handles, labels):
    if not l in tmp_labels:
        tmp_handles.append(h)
        tmp_labels.append(l)

# sort legend
order = [0,4,2,1,3]
tmp_handles = [tmp_handles[i] for i in order]
tmp_labels = [tmp_labels[i] for i in order]
relapse = mlines.Line2D([], [], color=relapse, marker='None', linestyle='-', linewidth=5, label='relapse')
no_relapse = mlines.Line2D([], [], color=no_relapse, marker='None', linestyle='-', linewidth=5, label='no relapse')
tmp_handles.append(relapse)
tmp_handles.append(no_relapse)
tmp_labels.append('relapse')
tmp_labels.append('no relapse')
ax.legend(tmp_handles, tmp_labels)

ax.set_title('Relapse detection in adjuvant melanoma patients')
ax.set_yticks(y_pos,labels=tmp_patients['PATIENT-ID'])
ax.set_xlabel('Days post surgery')
ax.set_ylabel('Patients')
offset = (max_x-min_x)*0.05
plt.xlim(min_x-offset,max_x+offset+200)
fig.tight_layout()
fig.savefig(out_folder + r'Figure_Survival-Swimmer.svg', format='svg')
plt.close(fig)

# write excel list for figure data
logging.info('Writing Excel file...')
with pd.ExcelWriter(path=out_folder+'Data_Survival.xlsx') as writer:
    workbook = writer.book
    for key in source_data.keys():
        sheet_name = key
        worksheet = workbook.add_worksheet(key)
        worksheet.write(0,0,key)
        writer.sheets[key] = worksheet

        if key != 'Figure 4 E':
            source_data[key].to_excel(writer, sheet_name=sheet_name, index=False, startrow=2, startcol=2)

        if key == 'Figure 4 E':
            print(source_data[key])
            print(key)
            startrow = 2
            startcol = 2
            worksheet.write(startrow+1,startcol-1,'Patients')
            source_data[key]['Patients'].to_excel(writer, sheet_name=sheet_name, index=False, startrow=startrow, startcol=startcol)
            startrow = startrow + len(source_data[key]['Patients']) + 4
            worksheet.write(startrow+1,startcol-1,'Samples')
            source_data[key]['Samples'].to_excel(writer, sheet_name=sheet_name, index=False, startrow=startrow, startcol=startcol)

logging.info('Done.')