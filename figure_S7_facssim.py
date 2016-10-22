import numpy as np
import pdb
import pandas
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib import gridspec
import sys
from make_simulated import *
from labeler import Labeler
import glob
import re
import os
import readfcs
from shapely.geometry import MultiPoint, Point
import math

#########################################################
# experimental data
#########################################################
# fluorescein concentrations
fl = np.hstack((np.zeros(1),np.logspace(-9.5, -5, 10)))

#total read counts. Should be good but uneven
x = ['0', '10^-9.5', '10^-9', '10^-8.5', '10^-8', '10^-7.5', '10^-7', '10^-6.5', '10^-6', '10^-5.5', '10^-5']

# bin boundaries
# April 15 gates
bin_vals_16_4_15 = np.array([[1.4775788014225415, 2.245106453825187],
[2.2459637507318497, 2.9683536647430646],
[2.96837636283699, 3.6663288893768895],
[3.6675170811229103, 5.2856967691125885]])

#April 19 gates
bin_vals_16_4_19 = np.array([[1.4775788014225415, 2.245106453825187],
[2.2459637507318497, 2.85027248745035],
[2.8650326346119392, 3.474446976276228],
[3.4849588153986915, 5.2856967691125885]])

#April 21 gates
bin_vals_16_4_21 = np.array([[1.4775788014225415,2.200826012340419],
[2.216443456408671,2.8355123402887603],
[2.8355123402887603,3.474446976276228],
[3.4849588153986915,5.2856967691125885]])
#########################################################################
bin_vals = bin_vals_16_4_15

boundaries = np.array([0, 10**bin_vals[0,1], 10**bin_vals[1,1], 10**bin_vals[2,1], 1e20])

########
# basal mean values taken from flow cytometry reports
##########
basal = [194, 243, 177]
basal = basal[0]

rep1 = pandas.read_csv('data/replicate_1.csv')

sort_counts = pandas.read_csv('data/sort_counts_16.4.15.txt', delimiter='\t')
sort_counts = np.array(sort_counts.iloc[1:,1:]).astype(float)

all_N = np.zeros((11,4))
label_x = ['0', '10^-9.5', '10^-9', '10^-8.5', '10^-8', '10^-7.5', '10^-7', '10^-6.5', '10^-6', '10^-5.5', '10^-5']
for ii in range(len(label_x)):
    for jj in range(4):
        all_N[ii,jj] = int(rep1['fluorescein'+label_x[ii]+'bin'+str(jj)].sum())
######################
# end of experimental data
######################

#######################
# Simulate cells
#######################

bin_range = bin_vals[3,1] - bin_vals[3,0]

amplitude = 10**bin_vals[3,0] + np.random.rand()*10**bin_range 
usethis = np.where((rep1['CDR1H'] ==  'ACTTTTAGTGACTACTGGATGAACTGGGTC') & (rep1['CDR3H'] ==	 'ACCCCAGTAGTCCATACCATAGTAAGAACC'))[0]
KD = float(rep1['fit_KD'][usethis])
fit_fraction = float(rep1['fit_fraction'][usethis])

R, truex, fc_bound = generate_counts(KD, amplitude, basal, 0.1, sort_counts, all_N, fl, boundaries, [np.sqrt(1), np.sqrt(0.5)])
sort_counts = sort_counts[1:]
for ii in range(truex.shape[0]):
    sort_counts[ii] = truex[ii] * sort_counts[ii].sum()

bins = np.linspace(1,5, 100)

# Needed for proper focusing
plt.ion()
plt.close('all')

# Create figure with subplots and specified spacing
figsize=(3.5,5.6)
rows=14
cols=1
col = 1
fig, axes = plt.subplots(figsize=figsize)
gs = gridspec.GridSpec(23, 2) 
plt.subplots_adjust(
    bottom=0.06,
    top=0.95,
    left=0.17,
    right=0.96,
    wspace=0.6,
    hspace=0.0)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.12,ypad=.01,fontsize=10)

# For CDR1H and CDR3H
conc_labels = ['$0$'] + \
              ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)]

file_labels = ['0M'] + \
              ['10^%1.1fM'%x for x in np.arange(-9.5,-4.5,0.5)]

  
# Make plots
for ii in range(11):
    ax = plt.subplot(gs[ii*2:(ii*2+2),0])
    curr = np.log10(np.array(fc_bound[ii]))
    for jj in range(4):
        b0 = min(bins, key=lambda x:abs(x-np.log10(boundaries[jj])))
        b1 = min(bins, key=lambda x:abs(x-np.log10(boundaries[jj+1])))
        temp = curr[(curr<b1)&(curr>b0)]
        ax.hist((temp), bins = bins, histtype='stepfilled', color=[jj/3.,0,0], lw=0, weights=np.ones(temp.shape)/curr.shape[0])
    ax.set_yticks([])
    ax.set_ylabel(conc_labels[ii])
    ax.set_ylim([0,0.07])
    if ii == 0:
        ax.set_title('simulated')

    if not(ii ==10):
        ax.set_xticks([])
    else:
        ax.set_xticks([1,3,5])
        ax.set_xticklabels([r'$10^1$',r'$10^3$',r'$10^5$'])
        ax.set_xlabel('fluorescence')
        
    ax.set_ylabel(conc_labels[ii],rotation=0,ha='right')

    conc_label = conc_labels[ii]
    ax.set_ylabel(conc_label,rotation=0,ha='right', fontsize=mpl.rcParams['font.size'])
    if ii!=10:
        ax.set_xlabel('')
        ax.set_xticks([])
    
    if ii == 5:
        ax.text(-0.35,0.5,'fluorescein [M]', rotation=90,ha='right', fontsize=mpl.rcParams['font.size'],transform=ax.transAxes)
    if ii == 0:
        labeler.label_subplot(ax,'A')

vmin = 1e0
vmax = 1E7
#panel B
ax = plt.subplot(gs[0:9,1])
im = ax.imshow(sort_counts , interpolation='nearest', cmap=mpl.cm.hot, norm=LogNorm(vmin=vmin, vmax=vmax))
ax.set_ylabel('fluorescein [M]',labelpad=10)
ticks = range(11)
labels = ['0'] + ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)]
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.set_xlabel('bin', labelpad=1)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.title('cells',fontsize=mpl.rcParams['font.size'])
labeler.label_subplot(ax,'B')

divider = make_axes_locatable(ax)
width = axes_size.AxesY(ax, aspect=1/30.)
pad = axes_size.Fraction(0.75, width)
cax = divider.append_axes("right", size=width, pad=pad)

cbar = fig.colorbar(im, cax=cax, orientation='vertical')
cbar.solids.set_rasterized(True)

# Panel C
ax = plt.subplot(gs[13:22,1])

im = ax.imshow(R, interpolation='nearest', cmap=mpl.cm.hot, norm=LogNorm(vmin=vmin, vmax=vmax))
ax.set_ylabel('fluorescein [M]',labelpad=10)
ticks = range(11)
labels = ['0'] + ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)]
ax.set_yticks(ticks)
ax.set_yticklabels(labels)

ax.set_xlabel('bin', labelpad=1)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.title('reads',fontsize=mpl.rcParams['font.size'])

labeler.label_subplot(ax,'C')
divider = make_axes_locatable(ax)
width = axes_size.AxesY(ax, aspect=1/30.)
pad = axes_size.Fraction(0.75, width)
cax = divider.append_axes("right", size=width, pad=pad)
cbar = fig.colorbar(im, cax=cax, orientation='vertical')
cbar.solids.set_rasterized(True)

plt.show()
plt.savefig('./pdfs/figure_S7_facssim.pdf')


