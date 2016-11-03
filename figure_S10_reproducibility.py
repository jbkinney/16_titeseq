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

def plot_KD_z_compare(lib1, lib2, ax):
    
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    log_bounds = [-10,-4.5]
    
    pairs = dict()
    for k1,k3, KD, KD_sigma in zip(lib1['CDR1H'],lib1['CDR3H'], lib1['fit_KD'], lib1['fit_KD_sigma']):
        pairs[k1+k3] = [KD, np.nan, KD_sigma, np.nan]
            
    for k1,k3, KD, KD_sigma in zip(lib2['CDR1H'],lib2['CDR3H'], lib2['fit_KD'], lib2['fit_KD_sigma']):
        if k1+k3 not in pairs:
            pairs[k1+k3] = [np.nan, KD, np.nan, KD_sigma]
        else:
            pairs[k1+k3][1] = KD
            pairs[k1+k3][3] = KD_sigma

    KD1 = np.log10(np.array([v[0] for v in pairs.values()]))
    KD2 = np.log10(np.array([v[1] for v in pairs.values()]))

    KD_err_1 = np.array([v[2] for v in pairs.values()])
    KD_err_2 = np.array([v[3] for v in pairs.values()])
    z = (KD1 - KD2)*1./np.sqrt((KD_err_1+KD_err_2))
    
    usethis = np.isfinite(z)
    z = z[usethis]
    
    freq, bins = np.histogram(z, np.linspace(-6,6, 24))
    ax.plot(bins[:-1], freq/float(np.sum(freq)))

def plot_KD_compare(lib1, lib2, ax, make_colorbar=False):
    # Use global zorder variable
    #global zorder
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds

    log_bounds = [-10,-4.5]
    pairs = dict()
    for k1,k3, KD, KD_sigma in zip(lib1['CDR1H'],lib1['CDR3H'], lib1['fit_KD'], lib1['fit_KD_sigma']):
        pairs[k1+k3] = [KD, np.nan, KD_sigma, np.nan]
            
    for k1,k3, KD, KD_sigma in zip(lib2['CDR1H'],lib2['CDR3H'], lib2['fit_KD'], lib2['fit_KD_sigma']):
        if k1+k3 not in pairs:
            pairs[k1+k3] = [np.nan, KD, np.nan, KD_sigma]
        else:
            pairs[k1+k3][1] = KD
            pairs[k1+k3][3] = KD_sigma
                
    KD1 = np.clip(np.log10(np.array([v[0] for v in pairs.values()])),\
        log_kd_low, log_kd_hi)
    KD2 = np.clip(np.log10(np.array([v[1] for v in pairs.values()])),\
        log_kd_low, log_kd_hi)
    #KD1[KD1<log_bounds[0]]=log_bounds[0]
    #KD1[KD1>log_bounds[1]]=log_bounds[1]
    #KD2[KD2<log_bounds[0]]=log_bounds[0]
    #KD2[KD2>log_bounds[1]]=log_bounds[1]
    KD_err_1 = np.array([v[2] for v in pairs.values()])
    KD_err_2 = np.array([v[3] for v in pairs.values()])
    usethis = np.isfinite(KD1) & np.isfinite(KD2) & (np.sqrt(KD_err_1)<0.5) & (np.sqrt(KD_err_2)<0.5)
    #ax.errorbar(np.log10(KD1), np.log10(KD2), xerr=KD_err_1, yerr=KD_err_2, fmt='o')
    nbins = 32
    freq, xedges, yedges = np.histogram2d(
        KD1[usethis].flatten(), 
        KD2[usethis].flatten(),
        bins=[np.linspace(log_bounds[0], log_bounds[1], nbins),
             np.linspace(log_bounds[0], log_bounds[1], nbins)])

    xedges = xedges[0:-1] + (xedges[1] - xedges[0]) / 2
    yedges = yedges[0:-1] + (yedges[1] - yedges[0]) / 2
    [xx, yy] = np.meshgrid(xedges, yedges)

    lvls = np.logspace(start=-0.5,stop=3.,num=8,endpoint=True)
    lvls2 = np.logspace(start=0,stop=3,num=4,endpoint=True)
    freq[freq<0.1]=0.1
    freq[freq>1e3] = 1e3
    im = ax.contourf(xx, yy, freq.transpose(), \
        levels=lvls, norm = LogNorm(), \
        cmap = 'Greys', zorder= 2, linestyles=None)
    ax.set_aspect(1.)
    
    if make_colorbar:
        p3 = ax.get_position().get_points()
        x00, y0 = p3[0]
        x01, y1 = p3[1]

        # [left, bottom, width, height]
        position = ([x01, y0, 0.01, y1-y0])
        cbar = plt.colorbar(im, cax=plt.gcf().add_axes(position), orientation='vertical', ticks=lvls2)
        #cbar = plt.colorbar(im, orientation='vertical', ticks=lvls2)
        cbar.ax.set_yticklabels([r'$10^{%d}$'%np.log10(t) for t in lvls2])
        cbar.set_label(r'density',labelpad=2)
    
    #usethis = np.isfinite(KD1) & np.isfinite(KD2)
    [my_corr, pval] = pearsonr(KD1[usethis], KD2[usethis])
    
    
    ax.text(-7.25,-9, r'$r=%0.2f$'%(my_corr), zorder=20)
    ticks = [-10,-9, -8, -7, -6, -5]
    tick_labels = [r'$10^{-10}$',r'$10^{-9}$',r'$10^{-8}$',r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$']
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_yticklabels(tick_labels)
    ax.set_xlim(log_bounds)
    ax.set_ylim(log_bounds)
    ax.axvline(log_kd_low, linestyle=':', color='k')
    ax.axvline(log_kd_hi, linestyle=':', color='k')
    ax.axhline(log_kd_low, linestyle=':', color='k')
    ax.axhline(log_kd_hi, linestyle=':', color='k')
    ax.plot(log_bounds, log_bounds, '--', c='k')
    

def plot_expression_compare(lib1, lib2, ax, make_colorbar=False):

    
    # Always do this when plotting on a specified axis
    plt.sca(ax)
    # Affinity bounds
    lims = [0,2]
    
    pairs = dict()
    wt_ind1 = np.where((lib1['CDR1H_AA'] == cdr1_wtseq) & (lib1['CDR3H_AA'] == cdr3_wtseq))[0]
    wt_mE1 = np.nanmedian(lib1['expression'][wt_ind1])
    wt_ind2 = np.where((lib2['CDR1H_AA'] == cdr1_wtseq) & (lib2['CDR3H_AA'] == cdr3_wtseq))[0]
    wt_mE2 = np.nanmedian(lib2['expression'][wt_ind2])
    
    
    for k1,k3, mE, E in zip(lib1['CDR1H'],lib1['CDR3H'], lib1['expression'], np.array(lib1[['prob_cmyc0','prob_cmyc1','prob_cmyc2','prob_cmyc3']])):
        pairs[k1+k3] = [mE, np.nan, (E), np.ones(4)*np.nan ]
            
    for k1,k3, mE, E in zip(lib2['CDR1H'],lib2['CDR3H'], lib2['expression'], np.array(lib2[['prob_cmyc0','prob_cmyc1','prob_cmyc2','prob_cmyc3']])):
        if k1+k3 not in pairs:
            pairs[k1+k3] = [np.nan, mE, np.ones(4)*np.nan, (E)]
        else:
            pairs[k1+k3][1] = mE
            pairs[k1+k3][3] = (E)#-np.sum(E*np.log(E+(E==0)))

    mE1 = np.array([v[0] for v in pairs.values()])/wt_mE1
    mE2 = np.array([v[1] for v in pairs.values()])/wt_mE2
    E1 = np.array([v[2].tolist() for v in pairs.values()])
    E2 = np.array([v[3].tolist() for v in pairs.values()])
    
    mE1[mE1>lims[1]]=lims[1]
    mE2[mE2>lims[1]]=lims[1]
    usethis = np.isfinite(mE1) & np.isfinite(mE2) & (np.min(E1, axis=1)>0) & (np.min(E2, axis=1) >0)#& (mE1<2) & (mE2<2)
    nbins = 32
    [my_corr, pval] = pearsonr(mE1[usethis], mE2[usethis])
    H, xedges, yedges = np.histogram2d(
        mE1[usethis].flatten(), 
        mE2[usethis].flatten(),
        bins=[np.linspace(lims[0], lims[1], nbins),
             np.linspace(lims[0], lims[1], nbins)])

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
        #position = ax.get_position()
        #pdb.set_trace()
        #position[0] = position[0]+position[2]
        #position[2] = position[2]/10.
        #cbar = plt.colorbar(im, cax=plt.gcf().add_axes(position), orientation='vertical', ticks=lvls2)
        
        p3 = ax.get_position().get_points()
        x00, y0 = p3[0]
        x01, y1 = p3[1]

        # [left, bottom, width, height]
        position = ([x01, y0, 0.01, y1-y0])
        cbar = plt.colorbar(im, cax=plt.gcf().add_axes(position), orientation='vertical', ticks=lvls2)
        #cbar = plt.colorbar(im, orientation='vertical', ticks=lvls2)
        cbar.ax.set_yticklabels([r'$10^{%d}$'%np.log10(t) for t in lvls2])
        cbar.set_label(r'density',labelpad=2)
    #pdb.set_trace()
    #[my_corr, pval] = pearsonr(mE1[usethis], mE2[usethis])
    ax.text(0.1, 1.7, r'$r=$'+str(np.round(my_corr,2)), zorder=20) 
    
    #labeler.label_subplot(ax,'B')
    
    ticks = [0, 0.5, 1, 1.5, 2]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)  
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.plot(lims, lims, '--', c='k')

figsize=(6,5*2./3)
rows=2
cols=3
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom=0.1,
    top=0.95,
    left=0.1,
    right=0.9,
    wspace=0.4,
    hspace=0.4)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.05,ypad=.01,fontsize=10)

# Panel A
ax = axes[0,0]
plot_KD_compare(rep1, rep3, ax)
labeler.label_subplot(ax,'A')
ax.set_ylabel('$K_D$ [M], replicate 1', labelpad=2)
ax.set_xlabel('$K_D$ [M], replicate 3', labelpad=2)
#ax.set_aspect(1.)

ax = axes[1,0]
plot_expression_compare(rep1, rep3, ax)
labeler.label_subplot(ax,'B')
#ax.set_aspect(1.)
ax.set_ylabel('expression, replicate 1', labelpad=2)
ax.set_xlabel('expression, replicate 3', labelpad=2)

  
#ax = axes[2,0]
#plot_KD_z_compare(rep1, rep3, ax)
#ax.set_ylabel('probability', labelpad=2)
#ax.set_xlabel(r'$(K_{D,1} - K_{D,3})/\sqrt{\sigma_1^2+\sigma_3^2}$', labelpad=2)
#labeler.label_subplot(ax,'C')
#ax.set_aspect(12./0.7)

ax = axes[0,1]
plot_KD_compare(rep1, rep2, ax)
ax.set_ylabel('$K_D$ [M], replicate 1', labelpad=2)
ax.set_xlabel('$K_D$ [M], replicate 2', labelpad=2)

#ax.set_aspect(1.)


ax = axes[1,1]
plot_expression_compare(rep1, rep2, ax)
ax.set_ylabel('expression, replicate 1', labelpad=2)
ax.set_xlabel('expression, replicate 2', labelpad=2)
#ax.set_aspect(1.)

  
#ax = axes[2,1]
#plot_KD_z_compare(rep1, rep2, ax)
#ax.set_ylabel('probability', labelpad=2)
#ax.set_xlabel(r'$(K_{D,1} - K_{D,2})/\sqrt{\sigma_1^2+\sigma_2^2}$', labelpad=2)
#labeler.label_subplot(ax,'C')
#ax.set_aspect(12./0.6)


ax = axes[0,2]
plot_KD_compare(rep3, rep2, ax, make_colorbar=True)
ax.set_ylabel('$K_D$ [M], replicate 3', labelpad=2)
ax.set_xlabel('$K_D$ [M], replicate 2', labelpad=2)
#ax.set_aspect(1.)


ax = axes[1,2]
plot_expression_compare(rep3, rep2, ax, make_colorbar=True)
ax.set_ylabel('expression, replicate 3', labelpad=2)
ax.set_xlabel('expression, replicate 2', labelpad=2)


#ax = axes[2,2]

#plot_KD_z_compare(rep3, rep2, ax)
#ax.set_ylabel('probability', labelpad=2)
#ax.set_xlabel(r'$(K_{D,3} - K_{D,2})/\sqrt{\sigma_3^2+\sigma_2^2}$', labelpad=2)

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
plt.savefig('pdfs/figure_S10_reproducibility.pdf')
#plt.close()
