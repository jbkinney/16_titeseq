#!/usr/bin/env python
from __future__ import division
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb
import pandas
from helper import *
from scipy import stats
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


# Needed for proper focusing
plt.ion()
plt.close('all')

# Define colors
red = [0.8,0,0]
blue = [0,0,0.8]
gray = [0.4,0.4,0.4]
lightgray = [0.8,0.8,0.8,0.9]
black = [0.,0.,0.]


wtseq1 = 'TFSDYWMNWV'
seq1pos = np.arange(28,38)
optseq1_dict = {30:'G',31:'H'}
wtseq2 = 'GSYYGMDYWG'
seq2pos = np.arange(100,110)
optseq2_dict = {101:'A',102:'S',106:'E',108:'L'}

# Get affinity hetmaps
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
A_1h_map = np.array(A_heatmap[:10].T[by_group]) 
A_1h_map = np.clip(A_1h_map,-9.5,-5)
A_3h_map = np.array(A_heatmap[10:].T[by_group]) 
A_3h_map = np.clip(A_3h_map,-9.5,-5)

# Get non-wt values
A_1h_vals = A_1h_map[A_1h_map != A_wt]
A_3h_vals = A_3h_map[A_3h_map != A_wt]

# Get expression hetmaps
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
E_1h_map = np.array(E_heatmap[:10].T[by_group]/E_wt) 
E_3h_map = np.array(E_heatmap[10:].T[by_group]/E_wt) 

# Get non-wt values
E_1h_vals = E_1h_map[E_1h_map != 1]
E_3h_vals = E_3h_map[E_3h_map != 1]

# 
# This is the function that does the plotting
#
def histogram(ax, vals, wt, hist_type, nbins=19, title=''):
    fontsize=7

    # Compute histogram bins and xlim
    if hist_type=='A':
        bounds = [-9.5,-5.0]
        ylim = [0,120]
    elif hist_type=='E':
        bounds = [0.0,1.6]
        ylim = [0,50]
    dx = (bounds[1]-bounds[0])/nbins
    xlim = [bounds[0]-dx, bounds[1]+dx]
    bins = np.linspace(xlim[0],xlim[1],nbins+2)
    # Modify elements that hit boundaries
    vals[vals <= (bounds[0]+1e-2)] = bounds[0] - dx/2
    vals[vals >= (bounds[1]-1e-2)] = bounds[1] + dx/2

    gray=[.6,.6,.6]
    red=[1,.4,.4]
    blue=[.4,.4,1]

    # Histgoram values
    ax.plot([wt,wt],ylim,'--',color=black,linewidth=1, zorder=0)
    ax.hist(vals,bins=bins[1:-1],color=gray,linewidth=0, zorder=1)
    
    #ax.set_xticks(fontsize=fontsize)
    #ax.set_yticks(fontsize=fontsize)
    ax.set_ylim(ylim)
    ax.set_ylabel('substitutions',fontsize=fontsize)
    ax.set_title(title,fontsize=fontsize)

    if hist_type=='A':
        ticks = [-9.0, -8.0, -7.0, -6.0, -5.0] 
        #ticklabels = [r'$10^{%0.1f}$'%t for t in ticks]
        #ticklabels[0]= r'$\leq$' + ticklabels[0]
        #ticklabels[-1]= r'$\geq$' + ticklabels[-1]
        ax.hist(vals[vals <= bounds[0]],bins=bins[:2],color=red,linewidth=0, zorder=2)
        ax.hist(vals[vals >= bounds[1]],bins=bins[-2:],color=blue,linewidth=0, zorder=3)
        ax.plot([bounds[0],bounds[0]],ylim,':',color=black,linewidth=1, zorder=0)
        ax.plot([bounds[1],bounds[1]],ylim,':',color=black,linewidth=1, zorder=0)
        ax.set_xlim(xlim)
        ax.set_xticks(ticks)
        #ax.set_xticklabels(ticklabels)
        ax.set_xlabel('$\log_{10}\ K_D/ \mathrm{M}$',fontsize=fontsize)
    elif hist_type=='E':
        ticks = [0.0, 0.4, 0.8, 1.2, 1.6] 
        ax.set_xlim(bounds)
        ax.set_xticks(ticks)
        ax.set_xlabel('$E$',fontsize=fontsize)


# Perform simulations to assess the level of optimization
def sample_rand_seq_values(map,num_samples,num_char=20):
    
    # Verify map shape reflects the correct number of amino acids
    assert map.shape[0]==num_char, 'map.shape[0]==%d!, does not match num_char==%d'%(map.shape[0],num_char)

    # Get sequence length
    L = map.shape[1]

    # Initialize matrix to hold contributions to value
    chosen_vals = np.zeros([num_samples,L])

    # Choose random values for each position
    for i in range(L):
        chosen_vals[:,i] = np.random.choice(map[:,i],num_samples)

    # Return sum of values
    return chosen_vals.sum(axis=1)

# Perform simulations to assess the level of optimization
def sample_rand_opt_values(map,num_samples,num_char=20,num_muts=6):
    
    # Verify map shape reflects the correct number of amino acids
    assert map.shape[0]==num_char, 'map.shape[0]==%d!, does not match num_char==%d'%(map.shape[0],num_char)

    # Get sequence length
    L = map.shape[1]

    # Initialize matrix to hold contributions to value
    chosen_vals = np.zeros([num_samples,L])
    for i in range(L):
        chosen_vals[:,i] = np.random.choice(map[:,i],num_samples)

    # Choose random values for each position
    mask_template = np.zeros(L)
    mask_template[:num_muts] = 1
    mask = np.zeros([num_samples,L])
    for n in range(num_samples):
        mask[n,:] = np.random.permutation(mask_template)

    # Return sum of values
    return (chosen_vals*mask).sum(axis=1)

num_samples = int(1E7)
num_opt_samples = int(1E7)

#
# Discuss A
#

# Print statistics
print 'Mutations in 3H have a larger effect on $K_D$ than mutations in 1H.'
print 'Median log_10 K_D values. 1H:%0.2f, 3H:%0.2f'%\
    (np.median(A_1h_vals),np.median(A_3h_vals))
print 'One-sided Mann-Whitney U test for 1H K_Ds < 3H K_Ds: P=%0.1e'%\
    stats.mannwhitneyu(A_1h_vals,A_3h_vals,alternative='less').pvalue

num_muts = 190.
A_1h_kill_num = sum(A_1h_vals >= -5)
A_1h_kill_pct = A_1h_kill_num*100./num_muts

A_3h_kill_num = sum(A_3h_vals >= -5)
A_3h_kill_pct = A_3h_kill_num*100./num_muts

A_1h_great_num = sum(A_1h_vals <= -9.5)
A_1h_great_pct = A_1h_great_num*100./num_muts

A_3h_great_num = sum(A_3h_vals <= -9.5)
A_3h_great_pct = A_3h_great_num*100./num_muts

A_1h_weaken_num = sum(A_1h_vals >= A_wt)
A_1h_weaken_pct = A_1h_weaken_num*100./num_muts

A_3h_weaken_num = sum(A_3h_vals >= A_wt)
A_3h_weaken_pct = A_3h_weaken_num*100./num_muts

print '''
In both regions, mutations weakened binding, i.e., increased $K_D$:
%d (%d%%) for 1H and %d (%d%%) for 3H. 
'''%(A_1h_weaken_num,A_1h_weaken_pct,A_3h_weaken_num,A_3h_weaken_pct)

print '''
In both regions, a substantial fraction of mutations increased $K_D$ above
our detection limit of 1E-5 M: %d (%d%%) for 1H and %d (%d%%) for 3H. 
'''%(A_1h_kill_num,A_1h_kill_pct,A_3h_kill_num,A_3h_kill_pct)

print '''
Very few mutations reduced $K_D$ below our detection limit of 
1E-9.5 M: %d (%d%%) for 1H and %d (%d%%) for 3H. 
'''%(A_1h_great_num,A_1h_great_pct,A_3h_great_num,A_3h_great_pct)

# Simulate random sequences
A_1h_pvalue = sum(sample_rand_seq_values(A_1h_map - A_wt,num_samples) < 0)/num_samples
A_3h_pvalue = sum(sample_rand_seq_values(A_3h_map - A_wt,num_samples) < 0)/num_samples

print '''
To quantify the level of optimization in 1H and 3H, we computed the 
log_10 K_D values of %.1e randomly generated sequences, assuming additivity 
in $\Delta \log_10 K_D$ values of the mutations at each position. 
We find that both 1H and 3H sequences are highly optimized, with only %0.1e 
end %0.1e random sequences in each region having an energy predicted to be
as low or lower than WT. 
'''%(num_samples,A_1h_pvalue,A_3h_pvalue)

# 1H mutations in OPT
opt_1h_muts = [(0,2),(15,3)]
opt_3h_muts = [(1,1),(9,2),(19,6),(4,8)]

# Estimate OPT affinity
A_opt = sum([A_1h_map[c]-A_wt for c in opt_1h_muts]) + \
        sum([A_3h_map[c]-A_wt for c in opt_3h_muts]) + \
        A_wt

print '''
If we assume additivity in $\Delta \log_10 K_D$, the 1H and 3H mutations 
in OPT are predicted to reduce $K_D$ from the WT value of 1E%0.2fM to 1E%0.2fM.
This should not be taken as a quantitative prediction for OPT affinity. First, 
OPT contains other mutations outside of 1H and 3H, which might further
increase affinity. Second, on of the mutations on its own lowers K_D below our
detection threshold of 1E-9.5M; our $K_D$ calculation assumes a value of
1E-9.5, which underestimates the affinity-increasing effect of this mutation. 
'''%(A_wt,A_opt)

A_map = np.concatenate((A_1h_map,A_3h_map),axis=1)
rand_opt_vals = A_wt + \
    sample_rand_opt_values(A_map-A_wt,num_opt_samples,num_char=20,num_muts=6)
better_opt_frac = sum(rand_opt_vals <= A_opt)/num_opt_samples


print '''
Still, this result can be taken as an indication of the correctness of our
$K_D$ values. Indeed, OPT differs from WT by 6 mutations within the 1H and
3H variable regions. If 6 mutations are chosen at random within these 2
regions, the fraction of the time that they are expected to decrease as much as
or more than the calculation for OPT is %0.1e.
'''%(better_opt_frac)


#
# Discuss E
#

print 'Mutations in 3H and in 1H have similar effects on E on average.'
print 'Median E values. 1H:%0.3f, 3H:%0.3f'%\
    (np.median(E_1h_vals),np.median(E_3h_vals))
print 'two-sided Mann-Whitney U test for 1H Es != 3H Es: P=%0.1e'%\
    stats.mannwhitneyu(E_1h_vals,E_3h_vals,alternative='two-sided').pvalue

print '''
However, 3H mutations have a larger variance variance in mutational effect
than do mutations in 1H (P = %0.1e, Levene's test). Thus suggests that 1H
is more highly optimized for expression than 3H.
'''%\
stats.levene(E_1h_vals,E_3h_vals).pvalue

# Simulate random sequences
E_1h_pvalue = sum(sample_rand_seq_values(np.log(E_1h_map),num_samples) > 0)/num_samples
E_3h_pvalue = sum(sample_rand_seq_values(np.log(E_3h_map),num_samples) > 0)/num_samples

print '''
To quantify the level of optimization in 1H and 3H, we estimated the E
values of %.1e randomly generated sequences, assuming additivity 
in $\Delta \log_10 E$ values of the mutations at each position. 
We find that both 1H and 3H sequences are highly optimized. However, 1H is
substantially more optimized than 3H in this regard, with only a fraction
 %0.1e (for 1H) vs. %0.1e (for 3H) of random sequences having predicted
expression above WT.
'''%(num_samples,E_1h_pvalue,E_3h_pvalue)


# Activate if want to make plots
if True:

    # Create figure with subplots and specified spacing
    figsize=(3.42,3.42)
    rows = 2
    cols = 2    
    fig, axes = plt.subplots(rows,cols,figsize=figsize)
    plt.subplots_adjust(
        bottom = 0.12,
        top = 0.95,
        left = 0.12,
        right = 0.95,
        hspace = 0.5,
        wspace = 0.5)

    # Make a labler to add labels to subplots
    labeler = Labeler(xpad=.09,ypad=.0,fontsize=10)

    ax = axes[0,0]
    labeler.label_subplot(ax,'A')
    histogram(ax, A_1h_vals, wt=A_wt, hist_type='A', title='1H')

    ax = axes[0,1]
    labeler.label_subplot(ax,'B')
    histogram(ax, A_3h_vals, wt=A_wt, hist_type='A', title='3H')

    ax = axes[1,0]
    labeler.label_subplot(ax,'C')
    histogram(ax, E_1h_vals, wt=1.0, hist_type='E', title='1H')

    ax = axes[1,1]
    labeler.label_subplot(ax,'D')
    histogram(ax, E_3h_vals, wt=1.0, hist_type='E', title='3H')

    plt.show()
    plt.savefig('pdfs/figure_S9_histograms.pdf')


