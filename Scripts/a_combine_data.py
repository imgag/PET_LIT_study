import sys
import os
import argparse
import logging
import functions
import pandas as pd

# set basic pandas options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# parser command line
parser = argparse.ArgumentParser(description="Combine PET/LIT data for further analysis.")
parser.add_argument('-in', '--in_file', required=True, help='Path to in excel file.')
parser.add_argument('-out', '--out_file', required=True, help='Path to out file.')
parser.add_argument('-log', '--log', required=False, help='Path to log file.')
args = parser.parse_args()

# create folder if not exists
tmp = os.path.dirname(args.out_file)
if not os.path.exists(tmp):
    os.makedirs(tmp)
    logging.info('Created output folder {}'.format(tmp))

# logging
if args.log is not None:
    logging.basicConfig(filename=args.log, encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
else:
    logging.basicConfig(encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

logging.info('Running {}'.format(os.path.basename(__file__)))
logging.info('Arguments loaded with ArgumentParser:')
for arg,value in sorted(vars(args).items()):
    logging.info(' Argument {}: {}'.format(arg,value))

# read excel file and prepare column formats
patients = pd.read_excel(args.in_file, sheet_name='Patients')
patients['FILTER'] = patients['FILTER'].fillna('')
lb_samples = pd.read_excel(args.in_file, sheet_name='LB-Samples')
lb_samples['FILTER'] = lb_samples['FILTER'].fillna('')
lb_variants = pd.read_excel(args.in_file, sheet_name='LB-Variants')
lb_variants['FILTER'] = lb_variants['FILTER'].fillna('')
lb_variants['FILTER-VARIANT-FILE'] = lb_variants['FILTER-VARIANT-FILE'].fillna('')
petct = pd.read_excel(args.in_file, sheet_name='PET-CT')
s100 = pd.read_excel(args.in_file, sheet_name='S100')
ldh = pd.read_excel(args.in_file, sheet_name='LDH')

# define allowed offset for different parameters
offset_7d = 7

# variants, set filter column
# TODO refactor to just loop the variants instead of all patients repeatedly
for i3,r3 in lb_variants.iterrows():
    tmp_filter = r3['FILTER-VARIANT-FILE']

    # set filter column which is supposed to be empty
    # Filter parameters considered: tumour AF filter, chrX filter, low conf filter, gnomad / in-house filter with in-house < 5%
    # Filter parameters tried: r3['INDEL']=='yes', if 'LOW_mALT_COUNT' in tmp, if 'STRAND_BIAS' in tmp, if r3['Multi_AF'] != 0, if r3['Multi_ALT'] < 3
    # nb for variant statistics only INDEL, HIGH-POP-AF, OFF-TARGET and HOMOPOLYMER count
    filter = []
    if 'INDEL' in tmp_filter:
        filter.append('INDEL')
    if 'HIGH_POP_AF' in tmp_filter:
        filter.append('HIGH-POP-AF')
    if 'off-target' in tmp_filter:
        filter.append('OFF-TARGET')
    if 'OUTLIER' in tmp_filter:
        filter.append('OUTLIER')
    if 'HOMOPOLYMER' in tmp_filter:
        filter.append('HOMOPOLYMER')
    if r3['MULTI-DEPTH'] < 1000:
        filter.append('LOW-DEPTH')

    # set filter column
    lb_variants.loc[i3, 'FILTER'] = functions.vcf_add_filter_list(lb_variants.loc[i3, 'FILTER'], filter)

# samples, set filter column
logging.info('Adding  summary columns for samples...')
mean_multi_af = []
mean_multi_depth = []
count_variants_per_sample = []
lb_samples[['T0-TREATMENT-START','T1-TREATMENT-D21','FIRST-SAMPLE','LAST-SAMPLE','VARIANTS-DETECTED','COUNT-VARIANTS-DETECTED','TIME-TO-PROGRESSION','PFS-EVENT','CLOSEST-TO-PROGRESSION']] = ['no','no','no','no','no',None,None,None,'no']
lb_samples[['MEAN-MULTI-AF','MEAN-MULTI-DEPTH','COUNT_VARIANTS']] = [None,None,None]
for i1,r1 in lb_samples.iterrows():
    # check variants for this sample
    tmp1 = lb_variants[(lb_variants['PLASMA-ID']==r1['PLASMA-ID']) & (lb_variants['FILTER']=='')]

    # min of 3 variants required for downstream analysis
    if len(tmp1)==0:
        lb_samples.loc[i1, 'FILTER'] = functions.vcf_add_filter(lb_samples.loc[i1, 'FILTER'], 'NO-VARIANTS')
    if len(tmp1)<3:
        lb_samples.loc[i1, 'FILTER'] = functions.vcf_add_filter(lb_samples.loc[i1, 'FILTER'], 'TOO-FEW-VARIANTS')
    # check how many variants were detected and if enough enough variants were detected to be ctDNA positive
    if len(tmp1.index) > 0:
        lb_samples.loc[i1, 'COUNT-VARIANTS-DETECTED'] = len(tmp1[tmp1['PVAL']<0.05].index)
        if lb_samples.loc[i1, 'COUNT-VARIANTS-DETECTED'] >= 3:
            lb_samples.loc[i1, 'VARIANTS-DETECTED'] = 'yes'

    # store additional variables for later use
    lb_samples.loc[i1,'MEAN-MULTI-AF'] = tmp1.loc[tmp1['FILTER']=='','MULTI-ALLELEFREQ'].mean()
    lb_samples.loc[i1,'MEAN-MULTI-DEPTH'] = tmp1.loc[tmp1['FILTER']=='','MULTI-DEPTH'].mean()
    lb_samples.loc[i1,'COUNT_VARIANTS'] = len(tmp1)

    # check if patient's first Sample
    first_sample = lb_samples[lb_samples['PATIENT-ID'] == r1['PATIENT-ID']].sort_values(by='TREATMENT-DAY',ascending=True).iloc[0]['PLASMA-ID']
    if r1['PLASMA-ID'] == first_sample:
        lb_samples.loc[i1, 'FIRST-SAMPLE'] = 'yes'

    # check if last sample
    last_sample = lb_samples[lb_samples['PATIENT-ID'] == r1['PATIENT-ID']].sort_values(by='TREATMENT-DAY',ascending=False).iloc[0]['PLASMA-ID']
    if r1['PLASMA-ID'] == last_sample:
        lb_samples.loc[i1, 'LAST-SAMPLE'] = 'yes'

    # check if treatment start sample, i.e. the closest sample T0 (T0 + 7d)
    tmp =  lb_samples[(lb_samples['PATIENT-ID'] == r1['PATIENT-ID']) & (lb_samples['TREATMENT-DAY'] <= offset_7d)].sort_values(by='TREATMENT-DAY',ascending=False)
    if len(tmp.index)>0:
        if r1['PLASMA-ID'] == tmp.iloc[0]['PLASMA-ID']:
            lb_samples.loc[i1, 'T0-TREATMENT-START'] = 'yes'

    # check if day 21 (PET/CT) sample, i.e. the closest sample to T21 (T21 +/- 7d), nb looks for min of abs values
    day21 = 21
    tmp = lb_samples[(lb_samples['PATIENT-ID'] == r1['PATIENT-ID']) & (abs(day21 - lb_samples['TREATMENT-DAY']) <= offset_7d)].sort_values(by='TREATMENT-DAY',ascending=False)
    if len(tmp.index)>0:
        if r1['PLASMA-ID'] == tmp.iloc[0]['PLASMA-ID']:
            lb_samples.loc[i1, 'T1-TREATMENT-D21'] = 'yes'

    # check if closest to progression
    count = 0
    time_to_progression = patients.loc[patients['PATIENT-ID']==r1['PATIENT-ID'],'PFS-DURATION'].iloc[0]
    lb_samples.loc[(lb_samples['PATIENT-ID']==r1['PATIENT-ID']), 'TIME-TO-PROGRESSION'] = abs(lb_samples.loc[lb_samples['PATIENT-ID']==r1['PATIENT-ID'],'TREATMENT-DAY']-time_to_progression)
    lb_samples.loc[(lb_samples['PATIENT-ID']==r1['PATIENT-ID']), 'PFS-EVENT'] = patients.loc[patients['PATIENT-ID']==r1['PATIENT-ID'],'PFS-EVENT'].iloc[0]
    sample_closest_to_progression = lb_samples[lb_samples['PATIENT-ID']==r1['PATIENT-ID']].sort_values(by='TIME-TO-PROGRESSION',ascending=True).iloc[0]['PLASMA-ID']
    if r1['PLASMA-ID'] == sample_closest_to_progression:
        lb_samples.loc[i1, 'CLOSEST-TO-PROGRESSION'] = 'yes'


# patients add summary columns
count_samples_per_patient = []
patients['COUNT-UNFILTERED-SAMPLES'] = None
for i,r in patients.iterrows():
    tmp1 = lb_samples[(lb_samples['PATIENT-ID']==r['PATIENT-ID']) & (lb_samples['FILTER']=='')]
    patients.loc[i,'COUNT-UNFILTERED-SAMPLES'] = len(tmp1.index)
    if len(tmp1.index)<1:
        patients.loc[i, 'FILTER'] = functions.vcf_add_filter(patients.loc[i, 'FILTER'], 'MISSING-LBS')


# S100 add summary columns
s100[['T0-TREATMENT-START','T1-TREATMENT-21D']] = ['no','no']
for i1, r1 in patients.iterrows():
    # T0 - Treatment start
    count = 0

    for i2,r2 in s100[(s100['PATIENT-ID']==r1['PATIENT-ID']) & (s100['TREATMENT-DAY']<=offset_7d)].sort_values(by='TREATMENT-DAY',ascending=False).iterrows():
        if count==0:
            s100.loc[i2, 'T0-TREATMENT-START'] = 'yes'
            count += 1

    # T1 - Treatment day 21
    day21 = 21
    count = 0
    for i2,r2 in s100[(s100['PATIENT-ID']==r1['PATIENT-ID']) & (abs(day21-s100['TREATMENT-DAY'])<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        if count==0:
            s100.loc[i2, 'T1-TREATMENT-21D'] = 'yes'
            count += 1


# LDH add summary columns
ldh[['T0-TREATMENT-START','T1-TREATMENT-21D']] = ['no','no']
for i1, r1 in patients.iterrows():
    # T0 - Treatment start
    count = 0
    for i2,r2 in ldh[(ldh['PATIENT-ID']==r1['PATIENT-ID']) & (ldh['TREATMENT-DAY']<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        if count==0:
            ldh.loc[i2, 'T0-TREATMENT-START'] = 'yes'
            count += 1

    # T1 - Treatment day 21
    day21 = 21
    count = 0
    for i2,r2 in ldh[(ldh['PATIENT-ID']==r1['PATIENT-ID']) & (abs(day21-ldh['TREATMENT-DAY'])<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        if count==0:
            ldh.loc[i2, 'T1-TREATMENT-21D'] = 'yes'
            count += 1


# combine PET/CT and cfDNA time points
logging.info('Combining data...')
diff_pet_single = pd.DataFrame(columns=['PATIENT-ID', 'TREATMENT', 'PLASMA-ID', 'PLASMA-T0-TREATMENT-START', 'PLASMA-MEAN-MULTI-AF', 'TUMOR-DETECTED-PVALUE', 'PET-ID', 'PETCT-TLG', 'PETCT-MTV', 'PETCT-OFFSET','FILTER'])
filter = []
missing_pet = []
for i1, r1 in patients[patients['FILTER']==''].iterrows():
    patient_id = r1['PATIENT-ID']
    treatment = r1['TREATMENT']

    tmp_petct = petct[(petct['PATIENT-ID'] == patient_id)].copy()
    if len(tmp_petct)<1:
        missing_pet.append(str(patient_id))
        continue

    for i2,r2 in lb_samples[(lb_samples['PATIENT-ID']==patient_id) & (lb_samples['FILTER']=='')].iterrows():
        plasma_id = r2['PLASMA-ID']
        plasma_treatment_start = r2['T0-TREATMENT-START']

        # find closest PET/CT
        [pet_id,petct_tlg,petct_mtv,petct_offset] = ['','','','']
        if len(tmp_petct)>0:
            tmp_petct['DIFF-PET'] = tmp_petct['TREATMENT-DAY'] - r2['TREATMENT-DAY']
            tmp_petct = tmp_petct.sort_values(by='DIFF-PET',key=abs).reset_index(drop=True)
            pet_id = tmp_petct.loc[0,'PET-ID']
            petct_tlg = tmp_petct.loc[0,'PET-TLG']
            petct_mtv = tmp_petct.loc[0,'PET-MTV']
            petct_offset = tmp_petct.loc[0,'DIFF-PET']

        row = [patient_id, treatment, plasma_id, plasma_treatment_start, r2['MEAN-MULTI-AF'], r2['TUMOR-DETECTED-PVALUE'], pet_id, petct_tlg, petct_mtv, petct_offset, ';'.join(filter)]
        diff_pet_single.loc[len(diff_pet_single)] = row

# clean up diff_pet_single and check that each PET is only associated with one liquid biopsy
petids = diff_pet_single['PET-ID'].unique()
for petid in petids:
    tmp_diff = diff_pet_single.loc[diff_pet_single['PET-ID']==petid]

    if len(tmp_diff)<2:
        continue

    # keep smallest offset
    # consider refactoring and replacing part with row labels
    tmp_diff = tmp_diff.sort_values(by='PETCT-OFFSET',key=abs)
    for i2,r2 in tmp_diff[1:].iterrows():
        diff_pet_single.loc[i2,'FILTER'] = 'MULTIPLE-USAGE-OF-PET'

logging.info(' PET/CT singles missing for patients {}.'.format(','.join(missing_pet)))


# diff PET/CT pairs
diff_pet_pair = pd.DataFrame(columns=['PATIENT-ID', 'TREATMENT', 'PET1-ID', 'PETCT-TLG1', 'PETCT-MTV1', 'PET2-ID', 'PETCT-TLG2', 'PETCT-MTV2', 'PETCT-TLG-DIFF', 'PETCT-MTV-DIFF', 'PLASMA-ID1', 'TUMOR-DETECTED-PVALUE1', 'PLASMA-MEAN-MULTI-AF1', 'PLASMA-DAYS-OFFSET1', 'PLASMA-ID2', 'TUMOR-DETECTED-PVALUE2','PLASMA-MEAN-MULTI-AF2', 'PLASMA-DAYS-OFFSET2', 'PLASMA-DIFF-MULTI-AF','FILTER'])
missing_pet = []
filter = []
for i1, r1 in patients[patients['FILTER']==""].iterrows():

    patient_id = r1['PATIENT-ID']
    treatment = r1['TREATMENT']

    # PET information
    tmp_petct = petct[(petct['PATIENT-ID'] == patient_id)].copy().sort_values(by='TREATMENT-DAY').reset_index(drop=True)

    if len(tmp_petct)<2:
        missing_pet.append(str(patient_id))
        continue

    petct_tlg1 = tmp_petct.iloc[0]['PET-TLG']
    petct_mtv1 = tmp_petct.iloc[0]['PET-MTV']
    pet1_id = tmp_petct.iloc[0]['PET-ID']
    petct_tlg2 = tmp_petct.iloc[1]['PET-TLG']
    petct_mtv2 = tmp_petct.iloc[1]['PET-MTV']
    pet2_id = tmp_petct.iloc[1]['PET-ID']

    # find closest samples to PET/CT
    tmp_samples = lb_samples[(lb_samples['PATIENT-ID']==patient_id) & (lb_samples['FILTER']=='')].copy()

    if len(tmp_samples)<1:
        logging.warning(' No sample for patient {}. Skipping.'.format(str(patient_id)))
        continue
    if tmp_samples['MEAN-MULTI-AF'].isnull().any():
        logging.warning(' MEAN-MULTI-AF for patient {} missing. Skipping.'.format(str(patient_id)))
        continue

    # first PET/CT
    tmp_samples['DIFF-PET1'] = tmp_samples['TREATMENT-DAY'] - tmp_petct.iloc[0]['TREATMENT-DAY']
    tmp_samples = tmp_samples.sort_values(by='DIFF-PET1', key=abs).reset_index(drop=True)
    plasma_af1 = tmp_samples.iloc[0]['MEAN-MULTI-AF']
    plasma_offset1 = tmp_samples.iloc[0]['DIFF-PET1']
    plasma_id1 = tmp_samples.iloc[0]['PLASMA-ID']
    plasma_mrd1 = tmp_samples.iloc[0]['TUMOR-DETECTED-PVALUE']

    # second PET/CT
    tmp_samples['DIFF-PET2'] = tmp_samples['TREATMENT-DAY']-tmp_petct.iloc[1]['TREATMENT-DAY']
    tmp_samples = tmp_samples.sort_values(by='DIFF-PET2', key=abs).reset_index(drop=True)
    plasma_af2 = tmp_samples.iloc[0]['MEAN-MULTI-AF']
    plasma_offset2 = tmp_samples.iloc[0]['DIFF-PET2']
    plasma_id2 = tmp_samples.iloc[0]['PLASMA-ID']
    plasma_mrd2 = tmp_samples.iloc[0]['TUMOR-DETECTED-PVALUE']

    diff_pet_pair.loc[len(diff_pet_pair)] = [patient_id, treatment, pet1_id, petct_tlg1, petct_mtv1, pet2_id, petct_tlg2, petct_mtv2, (petct_tlg2-petct_tlg1), (petct_mtv2-petct_mtv1), plasma_id1, plasma_mrd1, plasma_af1, plasma_offset1, plasma_id2, plasma_mrd2, plasma_af2, plasma_offset2, (plasma_af2-plasma_af1),""]

logging.info(' PET/CT pairs missing for patients {}.'.format(','.join(missing_pet)))


# diff LB biomarker and S100
diff_s100 = pd.DataFrame(columns=['PATIENT-ID', 'TREATMENT', 'PLASMA-ID', 'TUMOR-DETECTED-PVALUE', 'PLASMA-MEAN-MULTI-AF', 'S100-ID', 'S100-CONC', 'S100-OFFSET','FILTER'])
filter = []
for i1, r1 in patients[patients['FILTER']==''].iterrows():
    patient_id = r1['PATIENT-ID']
    treatment = r1['TREATMENT']

    for i2,r2 in lb_samples[(lb_samples['PATIENT-ID']==patient_id) & (lb_samples['FILTER']=='')].iterrows():

        # find closest S100
        tmp_s100 = s100[(s100['PATIENT-ID'] == patient_id)].copy()

        s100_conc = ''
        s100_offset = ''
        s100_id = ''
        # TODO consider .iloc instead and keep index
        if len(tmp_s100.index)>0:
            tmp_s100['DIFF-S100']  = tmp_s100['TREATMENT-DAY'] - r2['TREATMENT-DAY']
            tmp_s100 = tmp_s100.sort_values(by='DIFF-S100', key=abs).reset_index(drop=True)
            s100_id = tmp_s100.loc[0,'S100-ID']
            s100_conc = tmp_s100.loc[0,'VALUE']
            s100_offset = tmp_s100.loc[0,'DIFF-S100']

        row = [patient_id, treatment,r2['PLASMA-ID'],  r2['TUMOR-DETECTED-PVALUE'], r2['MEAN-MULTI-AF'], s100_id, s100_conc, s100_offset, ';'.join(filter)]
        diff_s100.loc[len(diff_s100)] = row

# clean up diff_s100 and check that each S100 measurement is only associated with one liquid biopsy
tmp_s100 = diff_s100['S100-ID'].unique()
for ts in tmp_s100:
    tmp_diff = diff_s100.loc[diff_s100['S100-ID']==ts]
    if len(tmp_diff)<2:
        continue
    # keep smallest offset
    # TODO consider refactoring without row labels
    tmp_diff = tmp_diff.sort_values(by='S100-OFFSET', key=abs)
    for i2,r2 in tmp_diff[1:].iterrows():
        diff_s100.loc[i2,'FILTER'] = 'MULTIPLE-USAGE-OF-S100'


# diff LB biomarker and LDH
diff_ldh = pd.DataFrame(columns=['PATIENT-ID', 'TREATMENT', 'PLASMA-ID', 'TUMOR-DETECTED-PVALUE', 'PLASMA-MEAN-MULTI-AF', 'LDH-ID', 'LDH-CONC','LDH-OFFSET','FILTER'])
filter = []
for i1, r1 in patients[patients['FILTER']==''].iterrows():
    patient_id = r1['PATIENT-ID']
    treatment = r1['TREATMENT']

    for i2,r2 in lb_samples[(lb_samples['PATIENT-ID']==patient_id) & (lb_samples['FILTER']=='')].iterrows():
        # find closest LDH
        tmp_ldh = ldh[(ldh['PATIENT-ID'] == patient_id)].copy()
        ldh_conc = ''
        ldh_offset = ''
        # TODO consider .iloc and keep index
        if len(tmp_ldh)>0:
            tmp_ldh['DIFF-LDH'] = tmp_ldh['TREATMENT-DAY'] - r2['TREATMENT-DAY']
            tmp_ldh = tmp_ldh.sort_values(by='DIFF-LDH', key=abs).reset_index(drop=True)
            ldh_id = tmp_ldh.loc[0,'LDH-ID']
            ldh_conc = tmp_ldh.loc[0,'VALUE']
            ldh_offset = tmp_ldh.loc[0,'DIFF-LDH']

        row = [patient_id, treatment, r2['PLASMA-ID'], r2['TUMOR-DETECTED-PVALUE'], r2['MEAN-MULTI-AF'], ldh_id, ldh_conc, ldh_offset, ';'.join(filter)]
        diff_ldh.loc[len(diff_ldh)] = row

# clean up diff_pet_single and check that each PET is only associated to one liquid biopsy
tmp_ldh = diff_ldh['LDH-ID'].unique()
for ts in tmp_ldh:
    tmp_diff = diff_ldh.loc[diff_ldh['LDH-ID']==ts]
    if len(tmp_diff)<2:
        continue
    # keep smalles offset
    tmp_diff = tmp_diff.sort_values(by='LDH-OFFSET', key=abs)
    for i2,r2 in tmp_diff[1:].iterrows():
        diff_ldh.loc[i2,'FILTER'] = 'MULTIPLE-USAGE-OF-LDH'


# add timepoints excel sheet, T0 - TREATMENT-START <=+7d, T1 - TREATMENT-DAY 21 +/-7d, T2 - PET +/-7d
day21 = 21
offset_7d = 7
offset_21d = 21
rows = []
# TODO consider pd.Series.abs for sorting and Time TO XYZ
# TODO consider to filter samples without empty filter column
for i1, r1 in patients[patients['FILTER']==''].iterrows():
    patient_id = r1['PATIENT-ID']
    treatment = r1['TREATMENT']

    #
    pet_t0 = [None, None]
    for i2,r2 in petct[(petct['PATIENT-ID']==r1['PATIENT-ID']) & (petct['TREATMENT-DAY']<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        pet_t0 = [r2['PET-ID'],r2['PET-MTV']]
        break
    #
    t2 = None
    pet_t2 = [None, None]
    tmp_petct = petct[(petct['PATIENT-ID']==r1['PATIENT-ID'])].sort_values(by='TREATMENT-DAY', ascending=False)
    if len(tmp_petct)>1:
        for i2,r2 in tmp_petct.iterrows():
            t2 = r2['TREATMENT-DAY']
            pet_t2 = [r2['PET-ID'], r2['PET-MTV']]
            break

    ctdna_t0 = [None, None, None]
    for i2,r2 in lb_samples[(lb_samples['PATIENT-ID']==r1['PATIENT-ID']) & (lb_samples['TREATMENT-DAY']<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        ctdna_t0 = [r2['PLASMA-ID'], r2['MEAN-MULTI-AF'], r2['TUMOR-DETECTED-PVALUE']]
        break

    ctdna_t1 = [None, None, None]
    for i2,r2 in lb_samples[(lb_samples['PATIENT-ID']==r1['PATIENT-ID']) & (abs(day21-lb_samples['TREATMENT-DAY'])<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        ctdna_t1 = [r2['PLASMA-ID'], r2['MEAN-MULTI-AF'], r2['TUMOR-DETECTED-PVALUE']]
        break

    ctdna_t2 = [None, None, None]
    if t2 is not None:
        for i2,r2 in lb_samples[(lb_samples['PATIENT-ID']==r1['PATIENT-ID']) & (abs(t2-lb_samples['TREATMENT-DAY'])<=offset_21d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
            ctdna_t2 = [r2['PLASMA-ID'], r2['MEAN-MULTI-AF'], r2['TUMOR-DETECTED-PVALUE']]
            break

    #
    s100_t0 = [None, None]
    for i2,r2 in s100[(s100['PATIENT-ID']==r1['PATIENT-ID']) & (s100['TREATMENT-DAY']<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        s100_t0 = [r2['S100-ID'], r2['VALUE']]
        break
    #
    s100_t1 = [None, None]
    for i2,r2 in s100[(s100['PATIENT-ID']==r1['PATIENT-ID']) & (abs(day21-s100['TREATMENT-DAY'])<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        s100_t1 = [r2['S100-ID'], r2['VALUE']]
        break
    #
    s100_t2 = [None, None]
    if t2 is not None:
        for i2,r2 in s100[(s100['PATIENT-ID']==r1['PATIENT-ID']) & (abs(t2-s100['TREATMENT-DAY'])<=offset_21d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
            s100_t2 = [r2['S100-ID'], r2['VALUE']]
            break

    #
    ldh_t0 = [None, None]
    for i2,r2 in ldh[(ldh['PATIENT-ID']==r1['PATIENT-ID']) & (ldh['TREATMENT-DAY']<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        ldh_t0 = [r2['LDH-ID'], r2['VALUE']]
        break
    #
    ldh_t1 = [None, None]
    for i2,r2 in ldh[(ldh['PATIENT-ID']==r1['PATIENT-ID']) & (abs(day21-ldh['TREATMENT-DAY'])<=offset_7d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
        ldh_t1 = [r2['LDH-ID'], r2['VALUE']]
        break
    #
    ldh_t2 = [None, None]
    if t2 is not None:
        for i2,r2 in ldh[(ldh['PATIENT-ID']==r1['PATIENT-ID']) & (abs(t2-ldh['TREATMENT-DAY'])<=offset_21d)].sort_values(by='TREATMENT-DAY', ascending=False).iterrows():
            ldh_t2 = [r2['LDH-ID'], r2['VALUE']]
            break

    rows.append([patient_id,treatment] + ctdna_t0 + ctdna_t1 + ctdna_t2 + pet_t0 + pet_t2 + s100_t0 + s100_t1 + s100_t2 + ldh_t0 + ldh_t1 + ldh_t2)

cols = ['PATIENT-ID','TREATMENT','CTDNA-T0','CTDNA-T0-AF','CTDNA-T0-MRD','CTDNA-T1','CTDNA-T1-AF','CTDNA-T1-MRD','CTDNA-T2','CTDNA-T2-AF','CTDNA-T2-MRD','PET-T0','PET-T0-MTV','PET-T2','PET-T2-MTV','S100-T0','S100-T0-VALUE','S100-T1','S100-T1-VALUE','S100-T2','S100-T2-VALUE','LDH-T0','LDH-T0-VALUE','LDH-T1','LDH-T1-VALUE','LDH-T2','LDH-T2-VALUE']
timepoints = pd.DataFrame(data=rows,columns=cols)


# write data files
logging.info('Writing Excel file...')
with pd.ExcelWriter(path=args.out_file) as writer:
    patients.to_excel(writer, sheet_name='Patients', index=False)
    lb_samples.to_excel(writer, sheet_name='LB-Samples', index=False)
    lb_variants.to_excel(writer, sheet_name='LB-Variants', index=False)
    petct.to_excel(writer, sheet_name='PET-CT', index=False)
    s100.to_excel(writer, sheet_name='S100', index=False)
    ldh.to_excel(writer, sheet_name='LDH', index=False)
    diff_pet_single.to_excel(writer, sheet_name='LB_PET-Single', index=False)
    diff_pet_pair.to_excel(writer, sheet_name='LB_PET-Pair', index=False)
    diff_s100.to_excel(writer, sheet_name='LB_S100', index=False)
    diff_ldh.to_excel(writer, sheet_name='LB_LDH', index=False)
    timepoints.to_excel(writer, sheet_name='Timepoints', index=False)

logging.info('Done.')