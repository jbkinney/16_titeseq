#!/usr/bin/env python
#from __future__ import division
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdb
from scipy.stats import pearsonr
from matplotlib.colors import LogNorm
import pandas
from helper import *

from labeler import Labeler

plt.ion()
plt.close('all')

class Results: pass;

# Needed for proper focusing
plt.ion()
plt.close('all')

# Define colors
red = [0.8,0,0]
blue = [0,0,0.8]
gray = [0.4,0.4,0.4]
lightgray = [0.8,0.8,0.8,0.9]
black = [0.,0.,0.]

zorder = 0

[seq_hash, seq, seq_cdr] = make_Sequence_Hash(
    'data/CDR_library_July_5_2013_sequences.txt')


rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')

cdr1_wtseq = 'TFSDYWMNWV'
cdr3_wtseq = 'GSYYGMDYWG'
wt_key = cdr1_wtseq+cdr3_wtseq

def plot_fraction_compare(rep1, rep2, ax, make_colorbar=False):
    # Use global zorder variable
    global zorder
    
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    
    pairs = dict()
    normalize1 = 0
    normalize2 = 0
    for k1, k3, fraction in zip(rep1['CDR1H'],rep1['CDR3H'], rep1['fit_fraction']):
        pairs[k1+k3] = [fraction, np.nan]
        if np.isfinite(fraction):
            normalize1+=fraction
            
    for k1, k3, fraction in zip(rep2['CDR1H'],rep2['CDR3H'], rep2['fit_fraction']):
        if np.isfinite(fraction):
            normalize2+=fraction
        if k1+k3 not in pairs:
            pairs[k1+k3] = [np.nan, fraction]
        else:
            pairs[k1+k3][1] = fraction
    
    frac1 = np.log10(np.array([v[0]/normalize1 for v in pairs.values()]))
    frac2 = np.log10(np.array([v[1]/normalize2 for v in pairs.values()]))
    usethis = np.isfinite(frac1) & np.isfinite(frac2) 
    #ax.errorbar(np.log10(KD1), np.log10(KD2), xerr=KD_err_1, yerr=KD_err_2, fmt='o')
    log_bounds = [-6, -2]
    nbins = 32
    H, xedges, yedges = np.histogram2d(
        frac1[usethis].flatten(), 
        frac2[usethis].flatten(),
        bins=[np.linspace(log_bounds[0], log_bounds[1], nbins),
             np.linspace(log_bounds[0], log_bounds[1], nbins)])

    xedges = xedges[0:-1] + (xedges[1] - xedges[0]) / 2
    yedges = yedges[0:-1] + (yedges[1] - yedges[0]) / 2
    [xx, yy] = np.meshgrid(xedges, yedges)

    lvls = np.logspace(start=-0.5,stop=3.,num=8,endpoint=True)
    lvls2 = np.logspace(start=0,stop=3,num=4,endpoint=True)
    H[H<0.1]=0.1
    im = ax.contourf(xx, yy, H.transpose(), \
        levels=lvls, norm = LogNorm(), \
        cmap = 'Greys', zorder= 2, linestyles=None)

    ax.set_aspect(1.)
    if make_colorbar:
        p3 = ax.get_position().get_points()
        x00, y0 = p3[0]
        x01, y1 = p3[1]

        # [left, bottom, width, height]
        position = ([x01+0.03, y0, 0.01, y1-y0])
        cbar = plt.colorbar(im, cax=plt.gcf().add_axes(position), orientation='vertical', ticks=lvls2)
        cbar.ax.set_yticklabels([r'$10^{%d}$'%np.log10(t) for t in lvls2])
        cbar.set_label(r'relative density',labelpad=2)



def zipf(rep1, ax):
    # Use global zorder variable
    global zorder
    
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    
    pairs = dict()
    normalize = 0.
    for k1,k3, fraction in zip(rep1['CDR1H'],rep1['CDR3H'], rep1['fit_fraction']):
        pairs[k1+k3] = [fraction]
        if np.isfinite(fraction):
            normalize += fraction
    
    fractions = np.array([v[0]/normalize for v in pairs.values()])

    ranks = np.arange(len(fractions))
    fracs = np.sort(fractions)[::-1]
    baseline = 1E-7*np.ones(ranks.shape)
    #pdb.set_trace()
    ax.fill_between(ranks,baseline,fracs, color=[0.5,0.5,0.5])
    plt.yscale('log', nonposy='clip')
    ax.set_ylim([1e-6, 1e-1])
    

# Create figure with subplots and specified spacing
figsize=(6,4.0)
rows=2
cols=3
#fig, axes = plt.subplots(rows,cols,figsize=figsize)

bottom=0.15
top=0.95
width = 0.7
height=top-bottom
pad = 0.1
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom=0.1,
    top=0.95,
    left=0.1,
    right=0.9,
    wspace=0.6,
    hspace=0.4)
    
labeler = Labeler(xpad=.08,ypad=.01,fontsize=10)

ax = axes[0,0]
labeler.label_subplot(ax,'A')

lims = [-6,-2]
plot_fraction_compare(rep1, rep3, ax)
ax.plot(lims, lims, '--', c='k', zorder=10)
ax.set_ylabel('fraction of population\nreplicate 1', labelpad=2)
ax.set_xlabel('fraction of population\nreplicate 3', labelpad=2)
ticks = range(lims[0], lims[1]+1)
tick_labels = [r'$10^{'+str(t)+'}$' for t in ticks]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(tick_labels)
ax.set_yticklabels(tick_labels)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax = axes[0,1]
lims = [-6,-2]
plot_fraction_compare(rep1, rep2, ax)
ax.plot(lims, lims, '--', c='k', zorder=10)
ax.set_ylabel('fraction of population\nreplicate 1', labelpad=2)
ax.set_xlabel('fraction of population\nreplicate 2', labelpad=2)
ticks = range(lims[0], lims[1]+1)
tick_labels = [r'$10^{'+str(t)+'}$' for t in ticks]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(tick_labels)
ax.set_yticklabels(tick_labels)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax = axes[0,2]
lims = [-6,-2]
plot_fraction_compare(rep3, rep2, ax, make_colorbar=True)
ax.plot(lims, lims, '--', c='k', zorder=10)
ax.set_ylabel('fraction of population\nreplicate 3', labelpad=2)
ax.set_xlabel('fraction of population\nreplicate 2', labelpad=2)
ticks = range(lims[0], lims[1]+1)
tick_labels = [r'$10^{'+str(t)+'}$' for t in ticks]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(tick_labels)
ax.set_yticklabels(tick_labels)
ax.set_xlim(lims)
ax.set_ylim(lims)


ax = axes[1,0]
labeler.label_subplot(ax,'B')

zipf(rep1, ax)
ax.set_xticks([0,1000,2000,3000,4000])
ax.set_ylabel('fraction of population', labelpad=2)
#ax.set_xlabel('rank', labelpad=2)
plt.title('replicate 1')

ax = axes[1,1]
zipf(rep2, ax)
ax.set_xticks([0,1000,2000,3000,4000])
ax.set_xlabel('rank', labelpad=2)
plt.title('replicate 2')

ax = axes[1,2]
zipf(rep3, ax)
ax.set_xticks([0,1000,2000,3000,4000])
#ax.set_xlabel('rank', labelpad=2)
plt.title('replicate 3')


# Show and save plot
plt.show()
plt.savefig('pdfs/figure_S11_composition.pdf')
#plt.close()



