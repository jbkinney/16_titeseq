#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas 
from helper import *
from labeler import Labeler
from scipy.stats import pearsonr
from matplotlib.colors import LogNorm

xlim = [0, 1E-4]
xticks = [0, 1E-9,1E-7,1E-5]

flow_ylim = [10**1.5,10**4.5]
flow_yticks = [1E2,1E3,1E4]

seq_ylim = [-0.25, 3.25]
seq_yticks = [1E0,1E1,1E2]

kd_low = 10**-9.5
kd_hi = 10**-5.0

pad = 2

[seq_hash, seq, seq_cdr] = make_Sequence_Hash(
    'data/CDR_library_July_5_2013_sequences.txt')

def KD_fun(x, y, y_std, s, lb, KD0):
    out = [lb, KD0, s]
    y[y<0] = np.nan
    error = np.nanmean((1./y_std**2 * (np.log(y) - np.log(hill(out, x))))**2)
    return out, error

def plot_clone(fl, ax, name, xticklabels=False, ylabel=False):
    plt.sca(ax)
    check_clones = get_clone_data()
    clone = check_clones[name]
    base_color = np.array(clone['color'])

    num_replicates = len(clone['KD'])
    mean_log_KD = np.mean(np.log10(clone['KD']))
    std_log_KD = np.std(np.log10(clone['KD']))/np.sqrt(len(clone['KD']))
    print('Clone %s: KD = %.2f +- %.2f (%d replicates)'%\
        (name, mean_log_KD, std_log_KD, num_replicates))

    # Iterate over replicates
    fade_color = np.array([0.,0.,0.])
    weights = np.linspace(0.4,1,len(clone['ave']))
    for replicate in range(len(clone['ave'])):
        color = base_color*weights[replicate] + \
            fade_color*(1-weights[replicate])
        y = clone['ave'][replicate]# rma 15.12.31
        #y -= y[0] #rma 15.12.31
        x = np.array(fl)
        max_error = np.inf #rma 15.12.31
        s0 = np.nanmax(y)# rma 15.12.31, take previously calculated K_D, and fit the amplitudes, given 0 basal
        for s in np.logspace(np.log10(s0)-2, np.log10(s0)+1, 100):# rma 15.12.31
            temp, error = KD_fun(x, y, 1, s, y[0], clone['KD'][replicate])# rma 15.12.31
            if error<max_error:# rma 15.12.31
                pfit = temp# rma 15.12.31
                max_error = error #rma 15.12.31
        xsample = np.logspace(np.log10(10**-12),np.log10(xlim[1]),200)
        ysample = hill(pfit, xsample)
        ax.loglog(xsample, ysample, lw=2, c=color,label=name)
        usethis = np.isfinite(y)
        ax.scatter(fl[usethis], y[usethis], s=20, c=color, zorder=6, edgecolor='k',lw=0.5)
        #pdb.set_trace()
        KD = np.clip(clone['KD'][replicate],kd_low,kd_hi)
        ax.axvline(KD, linestyle='-', c=color, lw=1, \
            alpha=0.3, zorder=0)

    # # Label clone
    # if 'pJK' in name:
    #     title = 'WT'
    # elif 'opt' in name:
    #     title = 'OPT'
    # elif 'delta' in name:
    #     title = r'$\Delta$'
    # else:
    #     title = 'C%s'%name
    # ax.set_title(title, va='top', ha='left', x=0.1, y=.8, \
    #     fontsize=mpl.rcParams['font.size'])
    # x=0.05
    # y=0.8
    # if 'pJK' in name:
    #     title = 'WT' #%num_syn
    # elif 'opt' in name:
    #     title = 'OPT'
    #     x=0.7
    #     y=0.15
    # elif 'delta' in name:
    #     title = r'$\Delta$'
    # elif '107' == name:
    #     x=0.7
    #     y=0.15
    #     title = 'C%s'%name #,num_syn)
        
    # else:
    #     title = 'C%s'%name #,num_syn)
    # ax.set_title(title, va='top', ha='left', x=x, y=y, \
    #     fontsize=mpl.rcParams['font.size'])
        
    # Specify clone names 
    if 'pJK' in name:
        title = 'WT'
    elif 'opt' in name:
        title = 'OPT'
    else:
        title = 'C%s'%name

    # Specify placement of clone names
    x = 0.5
    lower_clones = ['OPT','C107','WT']
    y = 0.15 if title in lower_clones else 0.8

    # Display clone names
    ax.set_title(title, va='top', ha='center', x=x, y=y, \
        fontsize=mpl.rcParams['font.size'])

    plt.xscale('symlog', linthreshx=10**-9.75, linscalex=0.5)
    ax.set_xlim(xlim)
    ax.set_ylim(flow_ylim)
    
    ax.set_xticks(xticks)
    ax.set_yticks(flow_yticks)
    ax.set_xlabel('fluorescein [M]',labelpad=pad)
    ax.set_ylabel('mean fluroescence [au]',labelpad=pad, ha='center')

    if not xticklabels:
        ax.set_xticklabels([])
        ax.set_xlabel('')

    if not ylabel:
        ax.set_ylabel('',labelpad=pad)

    # Plot K_D bounds
    plt.plot([kd_low,kd_low],flow_ylim,':k',zorder=-100)
    plt.plot([kd_hi,kd_hi],flow_ylim,':k',zorder=-100)


def hill(p, F):
    p = np.abs(p)
    amplitude = p[2]
    return p[0]+amplitude * F/(F+p[1])


def get_clone_data():
    summary = []
    summaryg = []
    curves = pandas.read_csv('data/titration_curves.csv')


    check_clones = {}
    keys = set(curves['name'])
    conc = ['0', '-9.5','-9', '-8.5','-8', '-7.5','-7', '-6.5','-6', '-5.5','-5']
    for k in keys:
        inds = np.where(k==curves['name'])[0]
        CDR1 = curves[' CDR1H'][inds[0]][1:]
        CDR3 = curves[' CDR3H'][inds[0]][1:]
        if (CDR1 in seq_hash) and (CDR3 in seq_hash) and ((CDR1 == seq[1]) or (CDR3==seq[0])) or (k == 'opt') or (k=='delta'):
            if k == 'opt':
                color = [0,0.8,0]
            elif k == 'delta':
                color = [0.3,0.3,0.3]
            elif k == 'pJK36':
                color = [0.8,0.2,0.8]
            elif CDR1 == seq[1]:
                color = [0.8,0,0]
            else:
                color = [0,0,0.8]
            
            temp = {}
            temp['CDR1H'] = CDR1
            temp['CDR3H'] = CDR3
            temp['CDR1HAA'] = curves[' CDR1HAA'][inds[0]][1:]
            temp['CDR3HAA'] = curves[' CDR3HAA'][inds[0]][1:]
            #temp['color'] = colors[k]
            temp['color'] = color
            inds = np.where(k==curves['name'])[0]
            PE_vals = [curves[' PE_log10_fluorescein='+fl][inds].tolist() for fl in conc]
            
            temp['ave'] = np.array(np.array(PE_vals).T)
            temp['exp'] = np.array(curves[' expression'][inds])
            temp['KD'] = np.array(curves[' KD'][inds])

            check_clones[k] = temp
    return check_clones


def plot_titeseq(ax, libs, clone_name, xticklabels=False, ylabel=False):
    plt.sca(ax)
    x = [0, 10**-9.5, 10**-9, 10**-8.5, 10**-8, 10**-7.5, 10**-7, 10**-6.5, 10**-6, 10**-5.5, 10**-5]
    check_clones = get_clone_data()
    clone = check_clones[clone_name]
    
    num_syn = 0
    activelib = []
    already_plotted = False
    for libnum, lib in enumerate(libs):
        lib_max = np.max(lib['fit_saturation'])
        ind = np.where((clone['CDR1H']==lib['CDR1H']) & (clone['CDR3H']==lib['CDR3H']))[0]
        
        plotx = np.array(x)
        amp = []
        Kd = []
        
        base_color = np.array(clone['color'])
        fade_color = np.array([0.,0.,0.])
        weights = [.4,.7,1.0]
        for n, jj in enumerate(ind):
            color = base_color*weights[libnum] + \
                fade_color*(1-weights[libnum])
            xs = np.logspace(np.log10(10**-12),np.log10(xlim[1]),200)
            ys = hill([lib['basal'][jj], lib['fit_KD'][jj], lib['fit_saturation'][jj]], xs)
            ys = 3 * np.log10(ys)/np.log10(lib_max)    
            ax.semilogx(xs, ys, c=color, zorder=10, lw=2)
            KD = np.clip(lib['fit_KD'][jj],kd_low,kd_hi)
            ax.axvline(KD, linestyle='-', c=(1+color)/2., lw=1,\
                zorder=0)

            amp.append(lib['fit_saturation'][jj])
            Kd.append(lib['fit_KD'][jj])
            num_syn += 1

    plt.xscale('symlog', linthreshx=10**-9.75, linscalex=0.5)
    ax.set_xticks(xticks)
    ax.set_xlim(xlim)

    # Specify clone names 
    if 'pJK' in clone_name:
        title = 'WT'
    elif 'opt' in clone_name:
        title = 'OPT'
    else:
        title = 'C%s'%clone_name

    # Specify placement of clone names
    x = 0.5
    lower_clones = ['OPT','C107','WT']
    y = 0.15 if title in lower_clones else 0.8

    # Display clone names
    ax.set_title(title, va='top', ha='center', x=x, y=y, \
        fontsize=mpl.rcParams['font.size'])

    ax.set_xlabel('fluorescein [M]',labelpad=pad)
    ax.set_ylabel('mean bin',labelpad=pad,ha='center')

    if not xticklabels:
        ax.set_xticklabels([])
        ax.set_xlabel('')

    if not ylabel:
        ax.set_ylabel('',labelpad=pad)
    ax.set_yticks(range(4))
    ax.set_ylim(seq_ylim)

    # Plot K_D bounds
    plt.plot([kd_low,kd_low],seq_ylim,':k',zorder=-100)
    plt.plot([kd_hi,kd_hi],seq_ylim,':k',zorder=-100)

###
### Stuff for panel A

def plot_combine_clones(libs, log_bounds, ax):
    # Use global zorder variable
    global zorder
    
    # Always do this when plotting on a specified axis
    plt.sca(ax)

    check_clones = get_clone_data()

    # Compute average flow variance. Use this to estimate uncertainty in 
    # Flow KD values. 
    flow_KD_vars = {}
    lines = []
    KDs = []
    flow_KDs = []
    seq_KDs = []
    for key in check_clones.keys():
        # Always do this when plotting on a specified axis
        plt.sca(ax)
        clone = check_clones[key]

        # Get flow cytometry information on clones
        clone_KD = np.clip(clone['KD'],kd_low,kd_hi)

        flow_KD_mean = np.nanmean(np.log10(clone_KD))
        KDs.append(flow_KD_mean)
        flow_KD_std = np.nanstd(np.log10(clone_KD))/np.sqrt(len(clone_KD))
        print 'clone %s: %d replicates'%(key, len(clone_KD))
        flow_KD_vars[key] = np.nanvar(np.log10(clone_KD))/len(clone_KD)
        
        markersize = 30
        # Get Tite-Seq information on clones
        titseqKD = []
        titseqKD_sigma = []
        for lib in libs:
            ind = np.where((clone['CDR1H']==lib['CDR1H']) & (clone['CDR3H']==lib['CDR3H']))[0]
            if ind.shape:
                lib_color = clone['color']
                # Compute Tite-Seq mean log10 KD
                titseqKD = titseqKD + [np.clip(lib['fit_KD'][jj],kd_low,kd_hi) for jj in ind]
                titseqKD_sigma = titseqKD_sigma + [lib['fit_KD_sigma'][jj] for jj in ind]
            
        seq_KD_mean = np.nanmean(np.log10(np.array(titseqKD)))
        # Compute Uncertainty in Tite-Seq KD. 
        seq_KD_std = KD_std = np.std(np.log10(titseqKD))/np.sqrt(len(titseqKD))
            
        print '%s seq: log10_KD = %f +- %f'%(key,seq_KD_mean,seq_KD_std)
        print '%s flow: log10_KD = %f +- %f'%(key,flow_KD_mean,flow_KD_std)
        if 'pJK' in key:
            title = 'WT' #%num_syn
        elif 'opt' in key:
            title = 'OPT'
        elif 'delta' in key:
            title = r'$\Delta$'
        else:
            title = 'C%s'%key #,num_syn)
        lines.append('%s & %s & %s & $10^{%.2f \pm %.2f}$& $10^{%.2f \pm %.2f}$\\\\ \\hline \n'%(title, clone['CDR1HAA'],clone['CDR3HAA'], flow_KD_mean, flow_KD_std, seq_KD_mean, seq_KD_std))
        
        # If the mean is out of boundary, ignore
        #if (seq_KD_mean <= log_bounds[0] or seq_KD_mean >= log_bounds[1] or flow_KD_mean <= log_bounds[0] or flow_KD_mean >= log_bounds[1]) and (not (key == 'delta')) and (not (key == 'opt')) :
        #    continue
        
        seq_KD_mean = np.max([log_bounds[0], seq_KD_mean])
        seq_KD_mean = np.min([log_bounds[1], seq_KD_mean])
        flow_KD_mean = np.max([log_bounds[0], flow_KD_mean])
        flow_KD_mean = np.min([log_bounds[1], flow_KD_mean])
        
        dot_color = lib_color
        error_color = lib_color

        # Plot data
        ax.errorbar(flow_KD_mean, seq_KD_mean, 
            xerr=flow_KD_std, yerr = KD_std, \
            c = error_color, lw=1, capsize=0, zorder=zorder)
        zorder += 1
        ax.scatter(flow_KD_mean, seq_KD_mean, \
            c = dot_color, s=markersize, lw=0.5, edgecolor='k', zorder=zorder)
        zorder += 1

        # Record KDs
        flow_KDs.append(flow_KD_mean)
        seq_KDs.append(seq_KD_mean)

    ind = np.argsort(KDs)
    for ii in ind:
        KD_table.write(lines[ii])

    # Plot KD bounds
    ax.axvline(np.log10(kd_low), linestyle=':', color='k', zorder=-100)
    ax.axvline(np.log10(kd_hi), linestyle=':', color='k', zorder=-100)
    ax.axhline(np.log10(kd_low), linestyle=':', color='k', zorder=-100)
    ax.axhline(np.log10(kd_hi), linestyle=':', color='k', zorder=-100)

    # Compute R and P, then return these values
    R, P = pearsonr(seq_KDs,flow_KDs)
    return R,P

###
###
zorder = 0
red = [0.8,0,0]
blue = [0,0,0.8]
gray = [0.4,0.4,0.4]
lightgray = [0.8,0.8,0.8,0.9]
black = [0.,0.,0.]

KD_table = open('./pdfs/tableS1_comparison.tex', 'w')
header = '''\\begin{center}
\\begin{tabular}{ | l | l | l | l | l |}
\\hline\n
Clone & CDR1H & CDR3H & Clonal flow $K_D$ (M) & Tite-seq $K_D$ (M) \\\\ \\hline \n'''
KD_table.write(header)


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
figsize=(3.42,6)
rows=5
cols=2
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    #top=.55,
    #bottom=.05,
    top=0.98,
    bottom=0.38,
    left=.12,
    right=.95,
    hspace=0,
    wspace=.5)

# Tite-Seq vs Flow

# Panel C

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.07,ypad=-.01,fontsize=10)

# Position panel
#bottom=0.62
#top=0.98
bottom=0.05
top=0.30
left=0.30
right=0.75
height=top-bottom
width=right-left
ax = fig.add_axes([left,bottom,width,height])
labeler.label_subplot(ax,'C',xpad_adjust=.05,ypad_adjust=0)

log_bounds = [-10,-4.5]
lims = log_bounds

R,P = plot_combine_clones(all_reps, log_bounds, ax)

ax.plot(lims, lims, '--', c=gray, zorder=-10)
ax.set_ylabel('$K_D$ [M], Tite-Seq', labelpad=2)
ax.set_xlabel('$K_D$ [M], flow', labelpad=2)
ticks = range(int(np.ceil(lims[0])), int(lims[1])+1)
tick_labels = [r'$10^{%i}$'%(t) for t in ticks]

ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(tick_labels)
ax.set_yticklabels(tick_labels)
ax.set_xlim(lims)
ax.set_ylim(lims)

# Compute and how R^2 value
#xs = np.array([l.get_data()[0][0] for l in ax.lines[:-1]])
#ys = np.array([l.get_data()[1][0] for l in ax.lines[:-1]])
#usethis = np.isfinite(xs) & np.isfinite(ys)

#R, P = pearsonr(xs[usethis],ys[usethis])
P_order = np.floor(np.log10(P))
digits = np.round(P/(10**P_order),2)
ax.text(x=-9,y=-5.5,s='$R^2 = %0.2f$\n$P = %1.1f \cdot 10^{%i}$'%(R**2,digits, P_order),va='top',ha='left')

# Tite-Seq

# Panel A
labeler.label_subplot(axes[0,0],'A')
plot_titeseq(axes[0,0], all_reps,'opt')
plot_titeseq(axes[1,0], all_reps,'107')
plot_titeseq(axes[2,0], all_reps,'pJK36',ylabel=True)
plot_titeseq(axes[3,0], all_reps,'5')
plot_titeseq(axes[4,0], all_reps,'45',xticklabels=True)


# Flow

# Panel B
labeler.label_subplot(axes[0,1],'B')
fl = np.array([ 0,10**-9.5,10**-9, 10**-8.5, 10**-8,10**-7.5,10**-7,10**-6.5,10**-6,10**-5.5,10**-5])
plot_clone(fl,axes[0,1],'opt')
plot_clone(fl,axes[1,1],'107')
plot_clone(fl,axes[2,1],'pJK36',ylabel=True)
plot_clone(fl,axes[3,1],'5')
plot_clone(fl,axes[4,1],'45',xticklabels=True)



# Save figure
plt.show()
plt.savefig('pdfs/figure_4_accuracy.pdf')

# Needed for proper focusing
#plt.ion()
#plt.close('all')
    

footer = '''    \\end{tabular} \n
\\end{center}'''
KD_table.write(footer)

rsq = []
for rep in all_reps:
    KD = rep['fit_KD']
    mE = rep['expression']
    usethis = np.isfinite(KD) & np.isfinite(mE)
    r, pval = pearsonr(KD[usethis], mE[usethis])
    rsq.append(r**2)
    
print 'KD, expression correlation: ' + str(np.mean(rsq))
