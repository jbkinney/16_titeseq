#!/usr/bin/env python
#from __future__ import division
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas
from helper import *
from matplotlib import rc
from labeler import Labeler

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

def plot_KD(rep1, ax):
    # Use global zorder variable
    global zorder
    
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    log_bounds = [-9.5,-5]
    num_muts = lambda key: sum([not(k1==k2) for k1, k2 in zip(wt_key, key)])
    for ii in range(1, 4):
        keys = [s1+s3 for s1, s3 in zip(rep1['CDR1H_AA'], rep1['CDR3H_AA'])]
        ii_keys = [num_muts(key)==ii for key in keys]
        inds = np.where(ii_keys)[0]
        
        KDs = rep1['fit_KD'][inds]
        logKDs = np.log10(KDs)
        logKDs[logKDs<log_bounds[0]] = log_bounds[0]
        logKDs[logKDs>log_bounds[1]] = log_bounds[1]
        logKDs = logKDs[np.isfinite(logKDs)]
        freq, bins = np.histogram(logKDs, np.linspace(log_bounds[0], log_bounds[1], 20), normed=1)
        bins = np.linspace(log_bounds[0], log_bounds[1], 19)
        ax.semilogy(bins,freq, label=str(ii) +r' mutations$ $')
        plt.yscale('symlog', linthreshy=0.01, linscaley=0.5)
    
    ticks = [-9,-8,-7,-6,-5]
    tick_labels = [r'$10^{%i}$'%(t) for t in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    #ax.set_yticks([])
    fontProperties = {'size' : 8}
    rc('font',**fontProperties)

def plot_expression(rep1, ax):
    # Use global zorder variable
    global zorder
    
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    num_muts = lambda key: sum([not(k1==k2) for k1, k2 in zip(wt_key, key)])
    wt_ind = np.where((rep1['CDR1H_AA'] == cdr1_wtseq) & (rep1['CDR3H_AA'] == cdr3_wtseq))[0]
    wt_mE = np.nanmedian(rep1['expression'][wt_ind])
    for ii in range(1, 4):
        keys = [s1+s3 for s1, s3 in zip(rep1['CDR1H_AA'], rep1['CDR3H_AA'])]
        ii_keys = [num_muts(key)==ii for key in keys]
        inds = np.where(ii_keys)[0]
        
        mE = rep1['expression'][inds]
        mE = mE[np.isfinite(mE)]/wt_mE
        freq, bins = np.histogram(mE, np.linspace(0, 2, 20), normed=1)
        bins = np.linspace(0, 2, 19)
        ax.semilogy(bins, freq, label=str(ii) +' muts')
        plt.yscale('symlog', linthreshy=0.01, linscaley=0.5)
    

log_bounds = [-10,-4]

figsize=(5.13,4)
rows=2
cols=3
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom=0.1,
    top=0.95,
    left=0.15,
    right=0.95,
    wspace=0.02,
    hspace=0.4)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.08,ypad=.01,fontsize=10)

# Set colorbar limits
vmin = 1e3
vmax = 1E7

# Panel A
ax = axes[0,0]
plot_KD(rep1, ax)
plt.ylabel('relative probability', fontsize=8)
plt.title('replicate 1', fontsize=8)
labeler.label_subplot(ax,'A')


ax = axes[0,1]
plot_KD(rep2, ax)
ax.set_yticks([]) 
plt.title('replicate 2', fontsize=8)
plt.xlabel(r'$K_D$\ [M]')

  
ax = axes[0,2]
plot_KD(rep3, ax)
ax.set_yticks([])   
plt.title('replicate 3', fontsize=8)
plt.legend(loc='upper left', frameon=False, fontsize=8)

# Panel B
ax = axes[1,0]
plot_expression(rep1, ax)
ticks = [0, 0.5, 1, 1.5, 2]
ax.set_xticks(ticks)

labeler.label_subplot(ax,'B')
plt.ylabel('relative probability', fontsize=8)


ax = axes[1,1]
plot_expression(rep2, ax)
ax.set_yticks([])  
plt.xlabel(r'$E$', fontsize=8)
ticks = [ 0.5, 1, 1.5, 2]
ax.set_xticks(ticks)

ax = axes[1,2]
plot_expression(rep3, ax)
ax.set_yticks([])   
ticks = [ 0.5, 1, 1.5, 2]
ax.set_xticks(ticks)

# Show and save plot
plt.show()
plt.savefig('pdfs/figure_S12_multipoint.pdf')
#plt.close()



