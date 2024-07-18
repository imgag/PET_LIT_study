import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import lifelines
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from matplotlib.offsetbox import AnchoredText

def plot_kaplan_meier(group1_d, group1_e, group2_d, group2_e, path, title = 'Survival', label1 = 'Group1', label2 = 'Group2'):

    kmf_group1 = lifelines.KaplanMeierFitter().fit(group1_d, event_observed=group1_e, label=label1)
    kmf_group2 = lifelines.KaplanMeierFitter().fit(group2_d, event_observed=group2_e, label=label2)
    p = logrank_test(durations_A=group1_d, durations_B=group2_d, event_observed_A=group1_e,event_observed_B=group2_e).p_value

    color_group1 = '#D9AFD7'
    color_group2 = '#69addb'

    # 'KAPLAN-DATE','KAPLAN-EVENT'
    mm = 0.1/2.54
    fig1,ax1 = plt.subplots(figsize=(2*90*mm,2*67*mm),dpi=300)
    ax1 = kmf_group1.plot_survival_function(ax=ax1,ci_show=False,show_censors=True,color=color_group1)
    ax1 = kmf_group2.plot_survival_function(ax=ax1,ci_show=False,show_censors=True,color=color_group2)
    ax1.set_ylim(0, 1.1)
    add_at_risk_counts(kmf_group1, kmf_group2, ax=ax1)
    ax1.set_title(title)
    ax1.set_ylabel('Survival')
    ax1.set_xlabel('Timeline (days)')

    p_string = 'p = {:.3f}'.format(p)
    if p<0.001:
        p_string = 'p < 0.001'
    ax1.add_artist(AnchoredText('{}'.format(p_string), loc=8, frameon=False))

    fig1 = ax1.get_figure()
    fig1.tight_layout()
    fig1.savefig(path, format='svg')
    plt.close(fig1)


def plot_correlation(x, y, title, xlabel, ylabel, rho, p, out, log=False):

    # check that series are of same length and sort values
    frame = {x.name : x, y.name : y}
    df = pd.DataFrame(frame)
    df = df.sort_values(x.name, ascending=True)
    x = df[x.name]
    y = df[y.name]
    offset_x1 = min(x)*0.95
    offset_x2 = max(x)*1.05
    offset_y1 = min(y)*0.95
    offset_y2 = max(y)*1.10
    y_max = max(y)

    p = round(p,3)
    p_string = 'p={}'.format(p)
    if p<0.001:
        p_string = 'p<0.001'

    # plot data points and print correlation
    fig = plt.figure(figsize=(6, 6))
    ax0 = fig.add_subplot(111)

    d = np.polyfit(x, y, 1)
    f = np.poly1d(d)
    reg_x = x
    reg_y = f(x)

    if log:
        ax0.set_xscale('log')
        ax0.set_yscale('log')

        x = x + x[x!=0].min()/2     # Add .min()/2 to zeros
        y = y + y[y!=0].min()/2     # Add .min()/2 to zeros
        x_log = np.log10(x)
        y_log = np.log10(y)

        tmp = abs((max(x_log)-min(x_log))*0.05)
        offset_x1 = 10**(min(x_log) - tmp)
        offset_x2 = 10**(max(x_log) + tmp)

        tmp = abs((max(y_log)-min(y_log))*0.05)
        offset_y1 = 10**(min(y_log) - tmp)
        offset_y2 = 10**(max(y_log) + tmp*2)

        d = np.polyfit(x_log, y_log, 1)
        f = np.poly1d(d)
        reg_x = x
        reg_y = 10**f(x_log)

    # limit regression values to max y value
    tmp_rx = []
    tmp_ry = []
    for rx,ry in zip(reg_x,reg_y):
        if ry<=y_max:
            tmp_rx.append(rx)
            tmp_ry.append(ry)
    reg_x = tmp_rx
    reg_y = tmp_ry

    ax0.set_title(title)
    ax0.set_xlabel(xlabel.replace('-', ' '))
    ax0.set_ylabel(ylabel.replace('-', ' '))
    ax0.set_xlim(offset_x1, offset_x2)
    ax0.set_ylim(offset_y1, offset_y2)

    ax0.scatter(x, y)
    ax0.plot(reg_x, reg_y, linestyle='dashed', linewidth=2.0, color='grey')
    text = "n = {}, Spearman\'s rho = {:.2f} ({})".format(len(y), round(rho,2), p_string)
    ax0.text(0.99, 0.96, text, transform=ax0.transAxes, size=12, horizontalalignment='right', alpha=1, color='k')

    plt.tight_layout()
    fig.savefig(out, format='svg')
    plt.close(fig)
