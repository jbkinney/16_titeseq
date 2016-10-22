#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb
import pandas
from helper import *
plt.close('all')
plt.ion()

from labeler import Labeler

[seq_hash, seq, seq_cdr] = make_Sequence_Hash(
    'data/CDR_library_July_5_2013_sequences.txt')

# Load data sets
rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')
all_reps = [rep1, rep2, rep3]
cdr1_wtseq = 'TFSDYWMNWV'
cdr3_wtseq = 'GSYYGMDYWG'
wt_seq = cdr1_wtseq+cdr3_wtseq

aff_fun = lambda x, ind: np.nanmean(np.log10(np.array(x['fit_KD'])[ind]))
exp_fun = lambda x, ind: np.nanmean((np.array(x['expression'])[ind]))
aas = np.array(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])
by_group = [5,0,17,7,9,10,4,19,18,15,16,11,13,1,12,6,8,14,2,3]
aa_map = {aas[int_val]:ii for ii, int_val in zip(range(20), by_group)}

def c_matrix(lib, map_fun):
    synonymous = np.where((lib['CDR1H_AA'] == cdr1_wtseq) & (lib['CDR3H_AA'] == cdr3_wtseq))
    wt_val = map_fun(lib, synonymous)
    usethis = []
    for ii in range(lib.shape[0]):
        temp = np.sum([c1 != wt1 for c1,wt1 in zip(lib['CDR1H_AA'][ii], cdr1_wtseq)])+np.sum([c3 != wt3 for c3,wt3 in zip(lib['CDR3H_AA'][ii], cdr3_wtseq)])
        usethis.append(temp==1)
    
    single = lib[usethis]
    
    out = []
    for ii in range(20):
        temp = []
        for jj in range(20):
            curr_seq = list(wt_seq)
            curr_seq[ii] = aas[jj]
            curr_seq = ''.join(curr_seq)
            usethis = [s1+s3 == curr_seq for s1, s3 in zip(lib['CDR1H_AA'], lib['CDR3H_AA'])]
            ind = np.where(usethis)[0]
            temp.append(map_fun(lib, ind))
        out.append(temp)
        
    out = np.array(out)
    practical_min = np.min(out[np.isfinite(out)])
    out[~np.isfinite(out)] = practical_min
    return out, wt_val

def heatmap_color(vals, zero, plottype):
    myrange = np.max(vals)-np.min(vals)
    zero = (zero - np.min(vals))/myrange
    if plottype == 'kd':
        cdict = {'red':   ((0.0,  0, 1),
                   (zero, 1, 1.0),
                   (1.0,  0, 0)),

                   'green': ((0.0, 0, 0),
                   (zero, 1, 1),
                   (1.0,  0, 0)),

                   'blue':  ((0.0,  0, 0),
                   (zero,  1, 1),
                    (1.0,  1.0, 1))}
    else:
        cdict = {'red':   ((0.0,  0, 0),
                   (zero, 1, 1.0),
                   (1.0,  1, 1)),

                   'green': ((0.0, 0, 0),
                   (zero, 1, 1),
                   (1.0,  0, 0)),

                   'blue':  ((0.0,  0, 1),
                   (zero,  1, 1),
                    (1.0,  0, 0))}        
                    
    return cdict



# This is the function that does all of the plotting
def plot_panel(ax, heatmap, zero, wtseq, pos, plottype, optseq_dict):
    
    if plottype=='me':
        vlim = [0.6,1.4]
    elif plottype=='kd':
        vlim = [-9.5, -5.0]

    plt.sca(ax)
    heatmap_array = np.array(heatmap.T[by_group])
    # Fix WT value to that passed as "zero"; this is set to an average
    # of the value for both CDR1 and CDR3 data sets
    for ii, aa in enumerate(wtseq):
         heatmap_array[aa_map[aa],ii] = zero

    temp_map = mpl.colors.LinearSegmentedColormap('my_colormap', \
        heatmap_color(vlim, zero, plottype),256)
    #cax=ax.imshow(heatmap_array, interpolation='nearest', \
    #              cmap=temp_map, vmin=vlim[0], vmax=vlim[1])
    cax = ax.pcolor(heatmap_array, cmap=temp_map, vmin=vlim[0], vmax=vlim[1])
    opt_pos = optseq_dict.keys()
    opt_aa = optseq_dict.values()
    for ii, aa in enumerate(wtseq):
        plt.scatter(ii+0.5, aa_map[aa]+0.5, \
            marker='o', c=[0.8,0.2,0.8], linewidths=0.5, s=10)

        # Plot OPT seq mutation if any occurs at position ii
        if pos[ii] in opt_pos:
            plt.scatter(ii+0.5,aa_map[optseq_dict[pos[ii]]]+0.5, \
                marker='o', c='Lime', linewidths=0.5, s=10)

    ax.set_yticks(np.linspace(0.5,19.5,20))
    ax.set_yticklabels([aas[ind] for ind in by_group], ha='left')
    ax.set_xlabel('VH position',labelpad=2)
    ax.set_xticks(np.linspace(0.5,9.5,4))
    ax.set_xticklabels([str(ii) for ii in pos[::3]])
    [tick.set_color(aa_colors2[ii]) for (ii,tick) in zip(aas[by_group],ax.yaxis.get_ticklabels())]
    ax.tick_params(axis='y', which='major', pad=10)

    fig = ax.get_figure()
    if plottype=='me':
        cbar = fig.colorbar(cax, orientation='vertical')
        ticks = np.linspace(start=vlim[0], stop=vlim[1],endpoint=True,num=5)
        ticklabels = [r'$%0.1f$'%t for t in ticks]
        if vlim[0]>0:
            ticklabels[0]= r'$\leq$' + ticklabels[0]
        else:
            ticklabels[0] = ticklabels[0]
        ticklabels[-1]= r'$\geq$' + ticklabels[-1]
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels(ticklabels)
        cbar.set_label(r'$E$',labelpad=0)

    elif plottype=='kd':
        vlim = [-9.5,-5.0]
        ticks = [-9.5, -9.0, -8.0, -7.0, -6.0, -5.0] 
        ticklabels = [r'$10^{%0.1f}$'%t for t in ticks]
        ticklabels[0]= r'$\leq$' + ticklabels[0]
        ticklabels[-1]= r'$\geq$' + ticklabels[-1]
        cbar = fig.colorbar(cax, orientation='vertical', ticks=ticks)
        cbar.ax.set_yticklabels(ticklabels)
        cbar.ax.tick_params()
        cbar.set_label(r'$K_D$ [M]',labelpad=0)

    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none') 

    plt.xlim([0,10])
    plt.ylim([0,20])
    cbar.solids.set_rasterized(True)

    # Return heatmap_array for further analysis
    return heatmap_array

# Needed for proper focusing
plt.ion()
plt.close('all')

# Define colors
red = [0.8,0,0]
blue = [0,0,0.8]
gray = [0.4,0.4,0.4]
lightgray = [0.8,0.8,0.8,0.9]
black = [0.,0.,0.]

# Create figure with subplots and specified spacing
figsize=(3.42,4.5)
rows = 2
cols = 2    
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom = 0.07,
    top = 0.95,
    left = 0.07,
    right = 0.88,
    hspace = 0.4,
    wspace = 0.6)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.035,ypad=.015,fontsize=10)

wtseq1 = 'TFSDYWMNWV'
seq1pos = np.arange(28,38)
optseq1_dict = {30:'G',31:'H'}
wtseq2 = 'GSYYGMDYWG'
seq2pos = np.arange(100,110)
optseq2_dict = {101:'A',102:'S',106:'E',108:'L'}

# Get affinity zero
A_heatmaps = []
A_wts = []
for rep in all_reps:
    temp_hm, wt_temp = c_matrix(rep, aff_fun)
    A_heatmaps.append(temp_hm)
    A_wts.append(wt_temp)

A_heatmap = np.zeros(A_heatmaps[0].shape)
for ii in range(A_heatmap.shape[0]):
    for jj in range(A_heatmap.shape[1]):
        A_heatmap[ii,jj] = np.nanmedian([heatmap[ii,jj] for heatmap in A_heatmaps])

A_wt = np.median(A_wts)

# Affinity plot, lib1
ax = axes[0,0]
labeler.label_subplot(ax,'A')
A_1h_map = plot_panel(ax, A_heatmap[:10], A_wt, wtseq1, seq1pos, 'kd', optseq1_dict)
ax.set_title('1H', fontsize=mpl.rcParams['font.size'])

# Affinity plot, lib2
ax = axes[0,1]
labeler.label_subplot(ax,'B')
A_3h_map = plot_panel(ax, A_heatmap[10:], A_wt, wtseq2, seq2pos, 'kd', optseq2_dict)
ax.set_title('3H', fontsize=mpl.rcParams['font.size'])

# Get expression zero
E_heatmaps = []
E_wts = []
for rep in all_reps:
    temp_hm, wt_temp = c_matrix(rep, exp_fun)
    E_heatmaps.append(temp_hm)
    E_wts.append(wt_temp)

E_heatmap = np.zeros(E_heatmaps[0].shape)
for ii in range(E_heatmap.shape[0]):
    for jj in range(E_heatmap.shape[1]):
        E_heatmap[ii,jj] = np.nanmedian([heatmap[ii,jj] for heatmap in E_heatmaps])
E_wt = np.median(E_wts)

# Expression plot, lib1
print '1H Results:'
ax = axes[1,0]
labeler.label_subplot(ax,'C')
E_1h_map = plot_panel(ax, E_heatmap[:10]/E_wt, 1, wtseq1, seq1pos, 'me', optseq1_dict)
ax.set_title('1H', fontsize=mpl.rcParams['font.size'])

# Expression plot, lib2
print '3H Results:'
ax = axes[1,1]
labeler.label_subplot(ax,'D')
E_3h_map = plot_panel(ax, E_heatmap[10:]/E_wt, 1, wtseq2, seq2pos, 'me', optseq2_dict)
ax.set_title('3H', fontsize=mpl.rcParams['font.size'])

plt.show()
plt.savefig('./pdfs/figure_5_landscape.pdf')
#plt.close()

A_1h_vals = A_1h_map[A_1h_map != A_wt]
A_3h_vals = A_3h_map[A_3h_map != A_wt]

E_1h_vals = E_1h_map[E_1h_map != 1]
E_3h_vals = E_3h_map[E_3h_map != 1]

