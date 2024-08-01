import os
import sys
import argparse
import logging
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import pandas as pd
from matplotlib.legend_handler import HandlerPatch
import subprocess as sp

# set basic pandas options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

logging.basicConfig(encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

logging.info('Running {}'.format(os.path.basename(__file__)))

# read in file
script_folder   = os.path.dirname(os.path.realpath(__file__)) + '/'
in_file = os.path.dirname(os.path.realpath(__file__)) + '/../Results/PetLit_Data_02.xlsx'
df_pat = pd.read_excel(in_file, sheet_name='Patients')
df_pat = df_pat[df_pat['TREATMENT']=='Pall']
df_ctdna = pd.read_excel(in_file, sheet_name='LB-Samples')
df_s100 = pd.read_excel(in_file, sheet_name='S100')
df_ldh = pd.read_excel(in_file, sheet_name='LDH')

#
rows = []
patients = df_pat.loc[pd.isnull(df_pat['FILTER']),'PATIENT-ID'].to_list()

column = 'COUNT-VARIANTS-DETECTED'
for pat in patients:
    ctdna = None
    filter = (df_ctdna['PATIENT-ID']==pat) & (df_ctdna['TREATMENT-DAY']<=7) & (df_ctdna['TREATMENT-DAY']>=-7)
    tmp = df_ctdna[filter].sort_values(by=['TREATMENT-DAY'],ascending=True)
    if len(tmp.index>0):
        ctdna = tmp.iloc[0][column]

    s100 = None
    filter = (df_s100['PATIENT-ID']==pat) & (df_s100['TREATMENT-DAY']<=7) & (df_s100['TREATMENT-DAY']>=-7)
    tmp = df_s100[filter].sort_values(by=['TREATMENT-DAY'],ascending=True)
    if len(tmp.index>0):
        s100 = tmp.iloc[0]['VALUE']

    ldh = None
    filter = (df_ldh['PATIENT-ID']==pat) & (df_ldh['TREATMENT-DAY']<=7) & (df_ldh['TREATMENT-DAY']>=-7)
    tmp = df_ldh[filter].sort_values(by=['TREATMENT-DAY'],ascending=True)
    if len(tmp.index>0):
        ldh = tmp.iloc[0]['VALUE']

    rows.append([pat,ctdna,s100,ldh])

df_venn = pd.DataFrame(data=rows,columns=['PATIENT-ID',column,'S100','LDH'])
df_venn = df_venn[(pd.notnull(df_venn['PATIENT-ID'])) & (pd.notnull(df_venn[column])) & (pd.notnull(df_venn['S100'])) & (pd.notnull(df_venn['LDH']))]

with pd.ExcelWriter(path=os.path.dirname(os.path.realpath(__file__)) + '/../Results/Data_Venn.xlsx') as writer:
    df_venn.to_excel(writer, sheet_name='Supplementary Fig. 6', index=False)
