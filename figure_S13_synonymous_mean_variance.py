#!/usr/bin/env python
#from __future__ import division
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdb
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LogNorm
import pandas
from helper import *

from labeler import Labeler

class Results: pass;

# Needed for proper focusing
plt.ion()
plt.close('all')

log_kd_low = -9.5
log_kd_hi = -5.0

# Define colors
red = [0.8,0,0]
blue = [0,0,0.8]
gray = [0.4,0.4,0.4]
lightgray = [0.8,0.8,0.8,0.9]
black = [0.,0.,0.]

[seq_hash, seq, seq_cdr] = make_Sequence_Hash(
    'data/CDR_library_July_5_2013_sequences.txt')


rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')

cdr1_wtseq = 'TFSDYWMNWV'
cdr3_wtseq = 'GSYYGMDYWG'
wt_key = cdr1_wtseq+cdr3_wtseq

def plot_KD_compare(lib, ax, plt_title='',make_colorbar=False):
    # Use global zorder variable
    #global zorder
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    log_bounds = [-9.5,-5]
    ab_pairs = dict()
    read_count_names = [k for k in lib.keys() if ('fluorescein' in k) and ('bin' in k) and ('prob' not in k)]
    for ID in lib.index:
        antibody = lib.loc[ID]
        k1 = antibody['CDR1H_AA']
        k3 = antibody['CDR3H_AA']
        KD = antibody['fit_KD']
        N = np.nansum(antibody[read_count_names])
        if k1+k3 in ab_pairs:
            ab_pairs[k1+k3].append([KD, N])
        else:
            ab_pairs[k1+k3] = [[KD, N]]
            
    
    KD_SSs = []
    mus = []
    for pairs in ab_pairs.values():
        if len(pairs)>1:
            mus.append(float(np.nanmean(np.log10(np.array(pairs)[:,0]))))
            KD_SSs.append(float(np.nanstd(np.log10(np.array(pairs)[:,0]))))

    mus = np.array(mus).flatten()
    KD_SSs = np.array(KD_SSs).flatten()
    ind = np.argsort(mus)
    mus = mus[ind]
    KD_SSs = KD_SSs[ind]
    
    nbins = 21
    freq, xedges, yedges = np.histogram2d(
        mus, 
        KD_SSs,
        bins=[np.linspace(log_bounds[0], log_bounds[1], nbins),
             np.linspace(0, log_bounds[1] - log_bounds[0], nbins)])

    xedges = xedges[0:-1] + (xedges[1] - xedges[0]) / 2
    yedges = yedges[0:-1] + (yedges[1] - yedges[0]) / 2
    [xx, yy] = np.meshgrid(xedges, yedges)

    lvls = np.logspace(start=-0.5,stop=2.,num=6,endpoint=True)
    lvls2 = np.logspace(start=0,stop=2,num=3,endpoint=True)
    freq[freq<0.1]=0.1
    freq[freq>1e2] = 1e2
    im = ax.contourf(xx, yy, freq.transpose(), \
        levels=lvls, norm = LogNorm(), \
        cmap = 'Greys', zorder= 2, linestyles=None)
    if len(plt_title) != 0 :
        plt.title(plt_title)
    
    ax.set_aspect(1)
    
    if make_colorbar:
        p3 = ax.get_position().get_points()
        x00, y0 = p3[0]
        x01, y1 = p3[1]

        # [left, bottom, width, height]
        position = ([x01-0.03, y0, 0.01, y1-y0])
        cbar = plt.colorbar(im, cax=plt.gcf().add_axes(position), orientation='vertical', ticks=lvls2)
        #cbar = plt.colorbar(im, orientation='vertical', ticks=lvls2)
        cbar.ax.set_yticklabels([r'$10^{%d}$'%np.log10(t) for t in lvls2])
        cbar.set_label(r'density',labelpad=2)
    
    ticks = [-10,-9, -8, -7, -6, -5]
    #tick_labels = [r'$10^{-10}$',r'$10^{-9}$',r'$10^{-8}$',r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$']
    ax.set_xticks(ticks)
    ax.set_yticks(range(5))
    #ax.set_xticklabels(tick_labels)
    #ax.set_yticklabels(tick_labels)
    ax.set_xlim(log_bounds)
    ax.set_ylim([0,log_bounds[1] - log_bounds[0]])
    

def plot_expression_compare(lib, ax, make_colorbar=False):
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    
    ab_pairs = dict()
    wt_ind = np.where((lib['CDR1H_AA'] == cdr1_wtseq) & (lib['CDR3H_AA'] == cdr3_wtseq))[0]
    wt_mE = np.nanmedian(lib['expression'][wt_ind])
    read_count_names = [k for k in lib.keys() if ('cmyc' in k) and ('prob' not in k)]
    for ID in lib.index:
        antibody = lib.loc[ID]
        k1 = antibody['CDR1H_AA']
        k3 = antibody['CDR3H_AA']
        expression = antibody['expression']/wt_mE
        N = np.nansum(antibody[read_count_names])
        if k1+k3 in ab_pairs:
            ab_pairs[k1+k3].append([expression, N])
        else:
            ab_pairs[k1+k3] = [[expression, N]]
            
    
    expression_SSs = []
    mus = []
    for pairs in ab_pairs.values():
        if len(pairs)>1:
            mus.append(float(np.nanmean((np.array(pairs)[:,0]))))
            expression_SSs.append(float(np.nanstd((np.array(pairs)[:,0]))))
    
    total_noise = np.nanvar((lib['expression']))
    expression_SSs = np.array(expression_SSs).flatten()
    mus = np.array(mus).flatten()
    
    ind = np.argsort(mus)
    mus = mus[ind]
    expression_SSs = expression_SSs[ind]
    
    nbins = 21
    freq, xedges, yedges = np.histogram2d(
        mus, 
        expression_SSs,
        bins=[np.linspace(0, 2, nbins),
             np.linspace(0, 2, nbins)])

    xedges = xedges[0:-1] + (xedges[1] - xedges[0]) / 2
    yedges = yedges[0:-1] + (yedges[1] - yedges[0]) / 2
    [xx, yy] = np.meshgrid(xedges, yedges)

    lvls = np.logspace(start=-0.5,stop=2.,num=6,endpoint=True)
    lvls2 = np.logspace(start=0,stop=2,num=3,endpoint=True)
    freq[freq<0.1]=0.1
    freq[freq>1e2] = 1e2
    im = ax.contourf(xx, yy, freq.transpose(), \
        levels=lvls, norm = LogNorm(), \
        cmap = 'Greys', zorder= 2, linestyles=None)
    ax.set_aspect(1)
    
    if make_colorbar:
        p3 = ax.get_position().get_points()
        x00, y0 = p3[0]
        x01, y1 = p3[1]

        # [left, bottom, width, height]
        position = ([x01-0.03, y0, 0.01, y1-y0])
        cbar = plt.colorbar(im, cax=plt.gcf().add_axes(position), orientation='vertical', ticks=lvls2)
        cbar.ax.set_yticklabels([r'$10^{%d}$'%np.log10(t) for t in lvls2])
        cbar.set_label(r'density',labelpad=2)
    
    ticks = [0, 1, 2]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    #ax.set_xticklabels(tick_labels)
    #ax.set_yticklabels(tick_labels)
    #ax.set_xlim(log_bounds)
    #ax.set_ylim(log_bounds)

figsize=(6,5*2./3)
rows=2
cols=3
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom=0.1,
    top=0.92,
    left=0.1,
    right=0.9,
    wspace=0.0,
    hspace=0.4)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.05,ypad=.02,fontsize=10)

# Panel A

print "A1"
ax = axes[0,0]
plot_KD_compare(rep1, ax, plt_title = 'replicate 3')
labeler.label_subplot(ax,'A')
ax.set_ylabel('standard deviation, log$_{10} K_D$', labelpad=2)
ax.set_xlabel(r'average log$_{10} K_D$', labelpad=2)

print "B1"
ax = axes[1,0]
plot_expression_compare(rep1, ax)
labeler.label_subplot(ax,'B')
ax.set_ylabel('standard deviation, $E$', labelpad=2)
ax.set_xlabel(r'average $E$', labelpad=2)

print "A2"
ax = axes[0,1]
plot_KD_compare(rep2, ax, plt_title = 'replicate 2')
ax.set_xlabel(r'average log$_{10} K_D$', labelpad=2)

print "B2"
ax = axes[1,1]
plot_expression_compare(rep2, ax)
ax.set_xlabel(r'average $E$', labelpad=2)

print "A3"
ax = axes[0,2]
plot_KD_compare(rep3, ax, plt_title = 'replicate 3', make_colorbar=True)
ax.set_xlabel(r'average log$_{10} K_D$', labelpad=2)

print "B3"
ax = axes[1,2]
plot_expression_compare(rep3, ax, make_colorbar=True)
ax.set_xlabel(r'average $E$', labelpad=2)

'''p3 = axes[2,1].get_position().get_points()
x00, y0 = p3[0]
x01, y1 = p3[1]

p3 = axes[1,2].get_position().get_points()
x10, junk = p3[0]
x11, junk = p3[1]

# [left, bottom, width, height]
ax.set_position([x10, y0, x11-x10, y1-y0])
'''
# Show and save plot
plt.show()
plt.savefig('pdfs/figure_synonymous_mean_variance.pdf')
#plt.close()
