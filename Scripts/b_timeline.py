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

#
export_patients = {'Figure 2 A1':'PAT-004','Figure 2 A2':'PAT-027','Figure 2 B1':'PAT-007','Figure 2 B2':'PAT-023','Figure 2 C1':'PAT-074','Figure 2 C2':'PAT-085','Supplementary Fig. 4':'PAT-017'}

# parser command line
parser = argparse.ArgumentParser(description="Create cfDNA Plots.")
parser.add_argument('-in', '--in_file', required=True, help='Path to in-file .')
parser.add_argument('-out', '--out_folder', required=True, help='Path to plot-folder.')
parser.add_argument('-log', '--log', required=False, help='Path to log file.')
args = parser.parse_args()

plot_folder = args.out_folder + '/Timelines'
prefix_plot = '{}{}'.format(plot_folder,'/Figure_')

# logging
if args.log is not None:
    logging.basicConfig(filename=args.log, encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
else:
    logging.basicConfig(encoding='utf-8', format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

logging.info('Running {}'.format(os.path.basename(__file__)))
logging.info('Arguments loaded with ArgumentParser:')
for arg,value in sorted(vars(args).items()):
    logging.info(' Argument {}: {}'.format(arg,value))

# set file paths
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

patients = pd.read_excel(args.in_file, sheet_name='Patients')
samples = pd.read_excel(args.in_file, sheet_name='LB-Samples')
variants = pd.read_excel(args.in_file, sheet_name='LB-Variants')
petct = pd.read_excel(args.in_file, sheet_name='PET-CT')
s100 = pd.read_excel(args.in_file, sheet_name='S100')
ldh = pd.read_excel(args.in_file, sheet_name='LDH')

# only keep unfiltered variants
variants = variants.fillna('')
variants = variants.loc[variants['FILTER']=='']

# plot variant data
logging.info('plotting')
ymin = 0
source_data = {}
for i,p in patients.iterrows():

    pid = p['PATIENT-ID']
    sid = 'na'
    treatment_start = 0

    if pid not in export_patients.values():
        continue

    # reduce dataset to relevant parameters for figure
    tmp_variants1 = variants.loc[(variants['PATIENT-ID']==pid), ['PATIENT-ID','PLASMA-ID','TREATMENT-DAY','VARIANT','MULTI-ALLELEFREQ']]
    tmp_samples1 = samples.loc[(samples['PATIENT-ID']==pid), ['PATIENT-ID','PLASMA-ID','TREATMENT-DAY','COUNT-VARIANTS-DETECTED']]
    tmp_petlit = petct.loc[(petct['PATIENT-ID']==pid), ['PATIENT-ID','TREATMENT-DAY','PET-TLG']]

    if pid in export_patients.values():
        key = list(export_patients.keys())[list(export_patients.values()).index(pid)]
        tmp_df = {'Samples':tmp_samples1,'Variants':tmp_variants1,'PET/CT':tmp_petlit}
        source_data.update({key:tmp_df})

    if (tmp_variants1.empty):
        logging.info(" No variants for patient {} ({})".format(pid,','.join(tmp_samples1['PLASMA-ID'])))
        continue

    fig = plt.figure(figsize=(10, 3 * 5))
    NUM_COLORS = len(tmp_variants1['VARIANT'].unique())
    cm = plt.get_cmap('summer')
    colors = [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)]

    ###############################
    # plot 1: all variants
    # format y-axis
    yaxis = 'MULTI-ALLELEFREQ'
    ax0 = plt.subplot(511)
    ax0.set_ylabel('Multi-AF')

    first_title = 'Liquid biopsies of {}\n'.format(pid)
    second_title = ''
    if p['TREATMENT'] == 'Pall':
        second_title = 'Response with combined ICI'
        if p['PROGRESS']==1:
            second_title = 'Progressive disease with combined ICI'
    elif p['TREATMENT'] == 'Adj':
        second_title = 'Response with adjuvant ICI'
        if p['PROGRESS'] == 1:
            second_title = 'Progressive disease with adjuvant ICI'
    # special title for patient 74
    if pid == 'PAT-074':
        second_title = 'Progressive disease with adjuvant ICI followed by response with combined ICI'
    ax0.set_title('{}{}'.format(first_title,second_title))

    tmp_max = max(tmp_variants1[yaxis])
    if tmp_max < ymin:
        tmp_max = ymin
    ymin1 = -tmp_max * 0.1
    ymax1 = tmp_max + (tmp_max-ymin1)*0.2
    if(ymax1==0):
        ymax1 = 0.05
    ax0.set_ylim(ymin1, ymax1)

    # format x-axis
    ax0.set_xlabel('Time (days)')
    tmp = tmp_variants1['TREATMENT-DAY'].tolist() + tmp_samples1['TREATMENT-DAY'].tolist() + [0]
    if p['TREATMENT'] == 'Pall':
        tmp += tmp_petlit['TREATMENT-DAY'].tolist()
    xmin1 = min(tmp)
    xmax1 = max(tmp)
    if (xmax1 - xmin1) < 30:
        xmax1 = 30
    offset = (xmax1 - xmin1) * 0.1
    xmin1 = xmin1 - offset
    xmax1 = xmax1 + offset
    # n.b. the following could remove PET/CTs from the plot for adjuvant patients
    ax0.set_xlim(xmin1, xmax1)

    # plot treatment start annotation
    treatment_start_delta = treatment_start - offset * 0.2
    if not pd.isnull(treatment_start):
        ax0.vlines(treatment_start, ymin=ymin1, ymax=ymax1, linestyles='dashed', colors='red')
    if pid == 'PAT-074':
        # start of combined IPI/NIVO on day 111
        ax0.vlines(111, ymin=ymin1, ymax=ymax1, linestyles='dashed', colors='steelblue')

    # plot PET annotation
    if p['TREATMENT'] == 'Pall':
        pet_delta = offset * 0.2
        if not tmp_petlit['PET-TLG'].empty:
            for i2,r2 in tmp_petlit.iterrows():
                blue_arrow = mpatches.FancyArrowPatch((r2['TREATMENT-DAY'],-(ymax1-ymin1)*0.22),(r2['TREATMENT-DAY'],-(ymax1-ymin1)*0.1),facecolor='steelblue',edgecolor='steelblue',label='PET/CT',arrowstyle='->',linewidth='2',mutation_scale=10,zorder=2)
                blue_arrow.set_clip_on(False)
                ax0.add_patch(blue_arrow)

    # plot stars
    for i2,r2 in tmp_samples1.iterrows():
        if(r2.loc['COUNT-VARIANTS-DETECTED'] >= 3):
            stars = '*'
            size = 6
            if (r2.loc['COUNT-VARIANTS-DETECTED'] >= 5):
                stars = '**'
                size = 10
                if (r2.loc['COUNT-VARIANTS-DETECTED'] >= 10):
                    stars = '***'
                    size = 18
            tmp_ymax1 = max(tmp_variants1.loc[tmp_variants1['PLASMA-ID']==r2['PLASMA-ID'],'MULTI-ALLELEFREQ'])
            ax0.plot(r2['TREATMENT-DAY'],tmp_ymax1+ymax1*0.1,marker='$%s$'%stars,label=stars,markersize=size,markeredgewidth=0.5,markeredgecolor='black',linewidth=1,markerfacecolor='black',zorder=3)

    # data points - variant per day
    variant = tmp_variants1['VARIANT'].unique()
    if len(variant)>len(colors):
        logging.warning('More variants than colors available for patient {} ({} vs {}).'.format(pid,len(variant),len(colors)))
    count_v = 0
    for v in variant:
        tmp_variants2 = tmp_variants1.loc[tmp_variants1['VARIANT']==v]
        tmp_variants2 = tmp_variants2.sort_values(by='TREATMENT-DAY')

        color = colors[count_v]
        count_v += 1

        ax0.plot(tmp_variants2['TREATMENT-DAY'], tmp_variants2[yaxis], label=v, marker='o', color=color)

    # add legend
    label = 'Surgery'
    if p['TREATMENT'] == 'Pall':
        label = 'Treatment start'

    legend_elements = []
    legend_elements.append(mlines.Line2D([], [], color='red', linestyle='--', label='Treatment start'))
    if p['TREATMENT'] == 'Pall':
        legend_elements.append(blue_arrow)
    if pid == 'PAT-074':
        legend_elements.append(mlines.Line2D([], [], color='steelblue', linestyle='--', label='Treatment change'))
    legend_elements.append(mlines.Line2D([],[],linestyle='None',label='\nctDNA detected'))
    legend_elements.append(mlines.Line2D([],[],marker=r'$%s$'%'*',linestyle='None',markeredgewidth=0.5,markeredgecolor='black',label='>= 3 signif. variants',markersize=6,markerfacecolor='black'))
    legend_elements.append(mlines.Line2D([],[],marker=r'$%s$'%'**',linestyle='None',markeredgewidth=0.5,markeredgecolor='black',label='>= 5 signif. variants',markersize=10,markerfacecolor='black'))
    legend_elements.append(mlines.Line2D([],[],marker=r'$%s$'%'***',linestyle='None',markeredgewidth=0.5,markeredgecolor='black',label='>= 10 signif. variants',markersize=18,markerfacecolor='black'))

    def make_legend_fancy_arrow_patch(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
        p = mpatches.FancyArrowPatch((0,0.5*height),(width,0.5*height),facecolor='steelblue',edgecolor='steelblue',label='PET/CT',arrowstyle='->',linewidth='2',mutation_scale=10)
        return p

    leg = ax0.legend(handles=legend_elements,handler_map={mpatches.FancyArrowPatch: HandlerPatch(patch_func=make_legend_fancy_arrow_patch), },bbox_to_anchor=(1.35, 1))
    for item,label in zip(leg.legendHandles, leg.texts):
        label.set_ha('left')
        if label._text in ['\nctDNA detected']:
            width = item.get_window_extent(fig.canvas.get_renderer()).width
            label.set_position((-1.4 * width, 0))

    ###############################
    # plot mean allele frequency of all variants
    tmp_samples1 = tmp_samples1.sort_values(by=['TREATMENT-DAY'])

    tmp_samples2 = samples.loc[samples['PATIENT-ID'] == pid]
    tmp_samples2 = tmp_samples2.sort_values(by=['TREATMENT-DAY'])

    ax1 = plt.subplot(512, sharex=ax0, sharey=ax0)
    ax1.set_xlabel('Time (Days)')
    yaxis = 'MEAN-MULTI-AF'
    ax1.set_ylabel(yaxis)
    ax1.set_title('Mean Allele Frequency')

    if not pd.isnull(treatment_start):
        ax1.vlines(treatment_start, ymin=ymin1, ymax=ymax1, linestyles='dashed', colors='red')

    color = colors[0]

    ax1.plot(tmp_samples2['TREATMENT-DAY'], tmp_samples2[yaxis], label='v', marker='o', color=color)

    for i2,ts1 in tmp_samples2.iterrows():
        tmp_color = color
        if (ts1.loc['TUMOR-DETECTED-PVALUE'] < 0.05):
            tmp_color = 'white'
        ax1.plot(ts1['TREATMENT-DAY'], ts1[yaxis], label='v', marker='o', markerfacecolor=tmp_color, markeredgecolor=color, )

    ###############################
    # cfDNA conc.
    tmp_samples2 = samples.loc[samples['PATIENT-ID'] == pid]
    tmp_samples2 = tmp_samples2.sort_values(by='TREATMENT-DAY')

    ax2 = plt.subplot(513, sharex=ax0)
    ax2.set_xlabel('Time (Days)')
    yaxis = 'CONC'
    ax2.set_ylabel(yaxis)
    ax2.set_title('Plasma conc. [ng/ml]')

    tmp_max = max(list(filter(None, tmp_samples2['CONC'])))
    if tmp_max < ymin:
        tmp_max = ymin
    if(pd.isna(tmp_max)):
        tmp_max = 0.05
    tmp_max += 0.01 * tmp_max
    ymin1 = -tmp_max * 0.1
    ymax1 = tmp_max * 1.1
    ax2.set_ylim(ymin1, ymax1)

    if not pd.isnull(treatment_start):
        ax2.vlines(treatment_start, ymin=ymin1, ymax=ymax1, linestyles='dashed', colors='red')
        ax2.text(treatment_start_delta, 0, 'start of treatment', rotation=90, color='red')

    ax2.plot(tmp_samples2['TREATMENT-DAY'], tmp_samples2[yaxis], label='v', marker='o', color='tab:gray')

    ###############################
    # plot biomarkers
    ax2 = plt.subplot(514, sharex=ax0)
    ax2.set_title('S100 and LDH')
    ax2.set_xlabel('Time (Days)')
    tmp_biomarker = s100.loc[s100['PATIENT-ID'] == pid]
    tmp_biomarker = tmp_biomarker.sort_values(by='TREATMENT-DAY')
    tmp_biomarker = tmp_biomarker[(tmp_biomarker['TREATMENT-DAY'] > xmin1) & (tmp_biomarker['TREATMENT-DAY'] < xmax1)]
    tmp_biomarker = tmp_biomarker.sort_values(by='TREATMENT-DAY')
    if not tmp_biomarker['VALUE'].empty:
        tmp_max = max(tmp_biomarker['VALUE'])
        tmp_max = tmp_max + 0.01 * tmp_max
        if tmp_max < 1:
            tmp_max = 1
        ymax1 = tmp_max * 1.1
        ymin1 = -tmp_max * 0.1

        ax2.set_ylim(ymin1, ymax1)
        ax2.set_ylabel('S100', color='tab:red')  # we already handled the x-label with ax1
        ax2.plot(tmp_biomarker['TREATMENT-DAY'], tmp_biomarker['VALUE'], label='v', marker='<', color='tab:red')

    if not pd.isnull(treatment_start):
        ax2.vlines(treatment_start, ymin=ymin1, ymax=ymax1, linestyles='dashed', colors='red')
        ax2.text(treatment_start_delta, 0, 'start of treatment', rotation=90, color='red')

    tmp_biomarker2 = ldh.loc[ldh['PATIENT-ID'] == pid]
    tmp_biomarker2 = tmp_biomarker2.sort_values(by='TREATMENT-DAY')
    tmp_biomarker2 = tmp_biomarker2[(tmp_biomarker2['TREATMENT-DAY'] > xmin1) & (tmp_biomarker2['TREATMENT-DAY'] < xmax1)]
    tmp_biomarker2 = tmp_biomarker2.sort_values(by='TREATMENT-DAY')
    if not tmp_biomarker2['VALUE'].empty:
        twin2 = ax2.twinx() # instantiate a second axes that shares the same x-axis
        tmp_max = max(tmp_biomarker2['VALUE'])
        tmp_max = tmp_max + 0.01 * tmp_max
        if tmp_max < 1:
            tmp_max = 1
        ymax1 = tmp_max * 1.1
        ymin1 = -tmp_max * 0.1
        twin2.set_ylim(ymin1, ymax1)
        twin2.set_ylabel('LDH', color='tab:gray')  # we already handled the x-label with ax1
        twin2.plot(tmp_biomarker2['TREATMENT-DAY'], tmp_biomarker2['VALUE'], label='v', marker='>', color='tab:gray')

    ###############################
    # plot PET/CT
    ax3 = plt.subplot(515, sharex=ax0)
    ax3.set_title('PET/CT TLG')
    ax3.set_xlabel('Time (Days)')

    if not tmp_petlit['PET-TLG'].empty:

        tmp_max = max(tmp_petlit['PET-TLG'])
        if tmp_max == 0:
            tmp_max = 0.01
        ymax1 = tmp_max * 1.1
        ymin1 = -tmp_max * 0.1
        ax3.set_ylim(ymin1, ymax1)

        ax3.plot(tmp_petlit['TREATMENT-DAY'], tmp_petlit['PET-TLG'], label='v', marker='o')

    if not pd.isnull(treatment_start):
        ax3.vlines(treatment_start, ymin=ymin1, ymax=ymax1, linestyles='dashed', colors='red')
        ax3.text(treatment_start_delta, 0, 'start of treatment', rotation=90, color='red')

    fig.tight_layout()
    fig.savefig("{}{}_{}.svg".format(prefix_plot,p['TREATMENT'],pid), format='svg')
    plt.close(fig=fig)

# write excel list for figure data
logging.info('Writing Excel file...')
with pd.ExcelWriter(path=args.out_folder+'Data_Timeline.xlsx') as writer:
    workbook = writer.book
    for key in source_data.keys():
        sheet_name = key
        worksheet = workbook.add_worksheet(key)
        worksheet.write(0,0,key)
        writer.sheets[key] = worksheet
        startrow = 2
        startcol = 2
        worksheet.write(startrow+1,startcol-1,'Variants')
        source_data[key]['Variants'].to_excel(writer, sheet_name=sheet_name, index=False, startrow=startrow, startcol=startcol)
        startrow = startrow + len(source_data[key]['Variants']) + 4
        worksheet.write(startrow+1,startcol-1,'Samples')
        source_data[key]['Samples'].to_excel(writer, sheet_name=sheet_name, index=False, startrow=startrow, startcol=startcol)
        startrow = startrow + len(source_data[key]['Samples']) + 4
        worksheet.write(startrow+1,startcol-1,'PET/CTs')
        source_data[key]['PET/CT'].to_excel(writer, sheet_name=sheet_name, index=False, startrow=startrow, startcol=startcol)

logging.info('Done.')