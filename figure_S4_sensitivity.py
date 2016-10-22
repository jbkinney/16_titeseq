#!/usr/bin/env python
import numpy as np
import pdb
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas
from scipy.stats import spearmanr
import os
import errno

from labeler import Labeler

def plot_rank_v_KD(ax, lib, yticks):
    #pA = lib['p_A']
    label_x = ['0', '10^-9.5', '10^-9', '10^-8.5', '10^-8', '10^-7.5', '10^-7', '10^-6.5', '10^-6', '10^-5.5', '10^-5']
    pA = []
    for kk in range(lib.shape[0]):
        prob = np.zeros((11,4))
        for ii in range(len(label_x)):
            for jj in range(4):
                prob[ii,jj] = float(lib['prob_fluorescein'+label_x[ii]+'bin'+str(jj)][kk])
        pA.append(prob)
    
    KD = np.array(lib['fit_KD'])
    lr = []
    for ii in range(11):
        lr.append([-np.log10(p[ii,2:4].sum()) + np.log10(p[ii,:2].sum()) for p in pA])
        
    lr = np.array(lr).T
    min_KD = max(-10,int(np.round(np.log10(np.nanmin(KD)))))
    max_KD = min(int(np.round(np.log10(np.nanmax(KD)))), -4)
    num_bins = int(np.round((max_KD-min_KD)*2+1))
    ranges = np.logspace(min_KD, max_KD, num_bins)
    ranges = np.hstack((0, ranges))
    cs = []
    for ii in range(ranges.shape[0]-2):
        usethis = (KD>ranges[ii]) & (KD<ranges[ii+2])
        subset = lr[usethis]
        subset[~np.isfinite(subset)]=0
        correlations = [spearmanr(subset[:,jj], np.log10(KD[usethis]))[0] for jj in range(11)]
        cs.append(correlations)
        
    im = ax.pcolor(np.array(cs).T, cmap=mpl.cm.coolwarm, vmin=-1, vmax=1)
    
    max_rho = max(np.array(cs).flatten())
    print 'max rho = %f'%max_rho
    
    num_conditions = 11
    ax.set_yticks(np.arange(num_conditions)+.5)

    # Option to turn off yticks
    if not yticks:
        ax.set_yticklabels([])
        ax.set_ylabel('')
    else:
        ax.set_ylabel('fluorescein [M]')
        ax.set_yticklabels(['$0$'] + ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)])
    
    ax.set_ylim([0,11])

    ticks = range(0,num_bins + 1, 2)
    ticks = [t+0.5 for t in ticks]
    ax.set_xticks(ticks)
    tick_labels = range(min_KD, max_KD+1)
    tick_labels = [r'$10^{'+ str(tick_label) +'}$' for \
        tick_label in tick_labels]
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel('$K_D$ [M]')
    ax.set_xlim([1,num_bins-1])

    return im, max_rho


rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')

# Needed for proper focusing
plt.ion()
plt.close('all')

# Create figure with subplots and specified spacing
figsize=(5.13,1.9)
rows = 1
cols = 3 
fig, axes = plt.subplots(rows,cols,figsize=figsize)
right = 0.85
bottom = 0.2
top = 0.9
plt.subplots_adjust(
    bottom=bottom,
    top=0.9,
    left=0.15,
    right=right,
    wspace=0.05,
    hspace=0.2)

# Plot results
im, max_rho = plot_rank_v_KD(axes[0],rep1, yticks=True)
title = 'replicate 1  ($\\rho_\mathrm{max}=%0.2f$)'%max_rho
axes[0].set_title(title,fontsize=mpl.rcParams['font.size'])

im, max_rho = plot_rank_v_KD(axes[1],rep2, yticks=False)
title = 'replicate 2  ($\\rho_\mathrm{max}=%0.2f$)'%max_rho
axes[1].set_title(title,fontsize=mpl.rcParams['font.size'])

im, max_rho = plot_rank_v_KD(axes[2],rep3, yticks=False)
title = 'replicate 2  ($\\rho_\mathrm{max}=%0.2f$)'%max_rho
axes[2].set_title(title,fontsize=mpl.rcParams['font.size'])

# Add colorbar (only need one im handel)
cax = fig.add_axes([right+.02, bottom, 0.03, top-bottom])
cbar = fig.colorbar(im, cax=cax, orientation='vertical', ticks=[-1, 0, 1])
cbar.set_label(r'$\rho$ ( $\log_{10} K_D$, enrichment)', labelpad=-1)
cbar.solids.set_rasterized(True)

# Save plot
plt.show()
plt.savefig('./pdfs/figure_S4_sensitivity.pdf')
#plt.close()
