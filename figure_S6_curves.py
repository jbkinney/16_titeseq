#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas 
from helper import *
from labeler import Labeler
from figure_utils import plot_clone, get_clone_data, plot_titeseq

plt.ion()

xlim = [0, 1E-4]
xticks = [1E-9,1E-7,1E-5]

flow_ylim = [10**1.9,10**4.3]
flow_yticks = [1E2,1E3,1E4]

seq_ylim = [0, 3.99]
seq_yticks = [1E0,1E1,1E2]

kd_low = 10**-9.5
kd_hi = 10**-5.0

pad = 2

[seq_hash, seq, seq_cdr] = make_Sequence_Hash(
    'data/CDR_library_July_5_2013_sequences.txt')


rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')
all_reps = [rep1, rep2, rep3]

###
### Make plot
###

# Needed for proper focusing
plt.ion()
plt.close('all')

# Create figure with subplots and specified spacing
figsize=(6,7)
rows=10
cols=4
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    top=.98,
    bottom=.05,
    left=.05,
    right=.95,
    hspace=0,
    wspace=.5)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.03,ypad=-.01,fontsize=10)

# fluorescein grid

summary = get_clone_data()
clones = summary.keys()
inds = np.argsort([np.nanmean(np.log10(np.array(summary[k]['KD']))) for k in clones])

# Panel B
labeler.label_subplot(axes[0,0],'A')
fl = np.array([ 0,10**-9.5,10**-9, 10**-8.5, 10**-8,10**-7.5,10**-7,10**-6.5,10**-6,10**-5.5,10**-5])

plot_titeseq(axes[0,0],all_reps,clones[inds[0]])
plot_titeseq(axes[1,0],all_reps,clones[inds[1]])
plot_titeseq(axes[2,0],all_reps,clones[inds[2]])
plot_titeseq(axes[3,0],all_reps,clones[inds[3]])
plot_titeseq(axes[4,0],all_reps,clones[inds[4]],ylabel=True)
plot_titeseq(axes[5,0],all_reps,clones[inds[5]])
plot_titeseq(axes[6,0],all_reps,clones[inds[6]])
plot_titeseq(axes[7,0],all_reps,clones[inds[7]])
plot_titeseq(axes[8,0],all_reps,clones[inds[8]])
plot_titeseq(axes[9,0],all_reps,clones[inds[9]],xticklabels=True)

plot_titeseq(axes[0,1],all_reps,clones[inds[10]])
plot_titeseq(axes[1,1],all_reps,clones[inds[11]])
plot_titeseq(axes[2,1],all_reps,clones[inds[12]])
plot_titeseq(axes[3,1],all_reps,clones[inds[13]])
plot_titeseq(axes[4,1],all_reps,clones[inds[14]],ylabel=True)
plot_titeseq(axes[5,1],all_reps,clones[inds[15]])
plot_titeseq(axes[6,1],all_reps,clones[inds[16]])
plot_titeseq(axes[7,1],all_reps,clones[inds[17]])
plot_titeseq(axes[8,1],all_reps,clones[inds[18]],xticklabels=True)
axes[9,1].axis('off')

# Panel A
labeler.label_subplot(axes[0,2],'B')

plot_clone(fl,axes[0,2],clones[inds[0]])
plot_clone(fl,axes[1,2],clones[inds[1]])
plot_clone(fl,axes[2,2],clones[inds[2]])
plot_clone(fl,axes[3,2],clones[inds[3]])
plot_clone(fl,axes[4,2],clones[inds[4]], ylabel=True)
plot_clone(fl,axes[5,2],clones[inds[5]])
plot_clone(fl,axes[6,2],clones[inds[6]])
plot_clone(fl,axes[7,2],clones[inds[7]])
plot_clone(fl,axes[8,2],clones[inds[8]])
plot_clone(fl,axes[9,2],clones[inds[9]],xticklabels=True)

plot_clone(fl,axes[0,3],clones[inds[10]])
plot_clone(fl,axes[1,3],clones[inds[11]])
plot_clone(fl,axes[2,3],clones[inds[12]])
plot_clone(fl,axes[3,3],clones[inds[13]])
plot_clone(fl,axes[4,3],clones[inds[14]], ylabel=True)
plot_clone(fl,axes[5,3],clones[inds[15]])
plot_clone(fl,axes[6,3],clones[inds[16]])
plot_clone(fl,axes[7,3],clones[inds[17]])
plot_clone(fl,axes[8,3],clones[inds[18]],xticklabels=True)
axes[9,3].axis('off')

plt.show()
plt.savefig('pdfs/figure_S6_curves.pdf')
#plt.close('all')


