#!/usr/bin/env python
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib import gridspec
import math
import readfcs
import sys
import pandas
import pdb
from matplotlib import path
from labeler import Labeler


def count_reads(ax, lib, name, ylabel, vmax, vmin):
    # Always do this when working with a specific axes
    plt.sca(ax)

    all_N = np.zeros((12,4))
    label_x = ['0', '10^-9.5', '10^-9', '10^-8.5', '10^-8', '10^-7.5', '10^-7', '10^-6.5', '10^-6', '10^-5.5', '10^-5']
    for ii in range(len(label_x)):
        for jj in range(4):
            all_N[ii+1,jj] = int(lib['fluorescein'+label_x[ii]+'bin'+str(jj)].sum())
    
    for jj in range(4):
        all_N[0, jj] = int(lib['cmyc'+str(jj)].sum())
    
    
    #for NA, NE in zip(lib['A'], lib['E']):
    #    all_N[0] += NE
    #    all_N[1:] += NA

    
    im = ax.imshow(all_N, interpolation='nearest', cmap=mpl.cm.hot, norm=LogNorm(vmin=vmin, vmax=vmax))
    if ylabel:
        ax.set_ylabel('fluorescein [M]',labelpad=-10)
        ticks = range(12)
        labels = ['expression  ', '0'] + ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)]
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
    else:
        ax.set_yticks([])

    ax.set_xlabel('bin', labelpad=1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.title(name,fontsize=mpl.rcParams['font.size'])

    return im

def count_cells(ax, datafile, name, ylabel, vmax, vmin):
    # Always do this when working with a specific axes
    plt.sca(ax)

    df = pandas.read_csv(datafile, delimiter='\t')
    all_N = np.array(df.iloc[1:,1:]).astype(float)

    im = ax.imshow(all_N, interpolation='nearest', cmap=mpl.cm.hot, norm=LogNorm(vmin=vmin, vmax=vmax))
    if ylabel:
        ax.set_ylabel('fluorescein [M]',labelpad=-10)
        ticks = range(12)
        labels = ['expression  ', '0'] + ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)]
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
    else:
        ax.set_yticks([])

    ax.set_xlabel('bin', labelpad=1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.title(name,fontsize=mpl.rcParams['font.size'])

    return im

def in_gate(data, gate, name1, name2):
    x_data = data[name1]
    y_data = data[name2]
    x_data = [math.log10(max(x, 1e-10)) for x in x_data]
    y_data = [math.log10(max(y, 1e-10)) for y in y_data]
    poly = path.Path(gate)
    data = np.array([x_data, y_data])
    usethis = poly.contains_points(data.T)
    return usethis

def filter_fsc_ssc(data):
    sub_data = lambda indata, gate: indata[in_gate(indata, gate, 'FSC-A','SSC-A')]

    gate = [(4.33487278721722, 4.336993061421767),
    (4.004403501790057, 4.010254050207731),
    (3.5539963164652484, 3.324117555127123),
    (3.761167625032098, 3.10094159611305),
    (4.797024167866345, 4.285490917033903),
    (4.701406640835492, 4.474332113122737),
    (4.33487278721722, 4.336993061421767)]
    
    data = sub_data(data, gate)

    return data

def filter_ps(data):
    inv_data = lambda indata, gate: indata[[not ii for ii in in_gate(indata, gate, 'FITC-A', 'PE-A')]]
    sub_data = lambda indata, gate: indata[[ii for ii in in_gate(indata, gate, 'FITC-A', 'PE-A')]]
    
    gate = [(2.1356696655075864, 1.6073794088650155),
    (2.8528011182389887, 1.6588815532528791),
    (2.9802911542801267, 2.6717570595475224),
    (4.637661622814923, 4.491499494585357),
    (3.7133588615166713, 5.058023082851852),
    (1.7850720663944564, 5.040855701389231),
    (1.5141557398070378, 1.7790532234912266),
    (1.6575820303533182, 1.487207738626668),
    (2.1356696655075864, 1.6073794088650155)]
    
    inside = sub_data(data, gate)

    return inside

def plot_affinity(ax, data, gates):
    plt.sca(ax)
    colors = [(0,0,0),(1./3,0,0),(2./3,0,0),(1,0,0)]
    x = np.array(data['PE-A'])
    bins = np.logspace(gates[0][0],gates[3][1], np.round((gates[3][1]-gates[0][0])/0.01))
    n, bins= np.histogram(x, bins)
    norm = 1./np.max(n)
    for gate, color in zip(gates, colors):
        usethis = (x>10**gate[0])&(x<10**gate[1])
        if sum(usethis):
            bins = np.unique((x[usethis]))
            bins = np.logspace(gate[0],gate[1], np.round((gate[1]-gate[0])/0.01))
            N = np.sum(usethis)
            n, bins, patches = plt.hist(x[usethis], bins, histtype='stepfilled', color=color, lw=0, weights = np.ones(N)*float(norm))
                  
    ax.set_xlabel(r'PE signal [au]',labelpad=3)
    ax.set_yticks([])
    #ax.set_ylabel(name)
    ax.set_xscale('log') 
    ax.set_xlim([30,10**gates[-1,-1]])
    ax.set_ylim([0,1])

    # Print mean log x
    usethis = (x>10)&(x<1E5)

    #print '%s: %f'%(name,np.mean(np.log10(x[usethis])))


def plot_expression(ax, data, gates):
    plt.sca(ax)
    colors = [(0,0,0),(1./3,0,1./3),(2./3,0,2./3),(1,0,1)]
    x = np.array(data['Brilliant Violet 421-A'])
    bins = np.logspace(gates[0][0],gates[3][1], np.round((gates[3][1]-gates[0][0])/0.01))
    n, bins= np.histogram(x, bins)
    norm = 1./np.max(n)
    for gate, color in zip(gates, colors):
        usethis = (x>10**gate[0])&(x<10**gate[1])
        if sum(usethis):
            bins = np.unique((x[usethis]))
            bins = np.logspace(gate[0],gate[1], np.round((gate[1]-gate[0])/0.01))
            N = np.sum(usethis)
            n, bins, patches = plt.hist(x[usethis], bins, histtype='stepfilled', color=color, lw=0, weights = np.ones(N)*float(norm))
    
    ax.set_xlabel(r'BV signal [au]',labelpad=3)
    ax.set_yticks([])
    ax.set_xscale('log') 
    ax.set_xlim([30,10**gates[-1,-1]])
    ax.set_ylim([0,1])

def get_filenames(path):
    files = []
    for infile in glob.glob( os.path.join(path, '*.fcs') ):
        files.append(infile)
    return files

out_names = ['./pdfs/figure_3_procedure.pdf', './pdfs/figure_S2_rep2.pdf', './pdfs/figure_S3_rep3.pdf']
rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')
reps = [rep1, rep2, rep3]
sort_counts = ['data/sort_counts_16.4.15.txt', 'data/sort_counts_16.4.19.txt', 'data/sort_counts_16.4.21.txt']
directories = ['data/fcs1/','data/fcs2/','data/fcs3/']
bin_vals_16_4_15 = np.array([[1.4775788014225415, 2.245106453825187],
[2.2459637507318497, 2.9683536647430646],
[2.96837636283699, 3.6663288893768895],
[3.6675170811229103, np.log10(3e4)]])

#April 19 gates
bin_vals_16_4_19 = np.array([[1.4775788014225415, 2.245106453825187],
[2.2459637507318497, 2.85027248745035],
[2.8650326346119392, 3.474446976276228],
[3.4849588153986915, np.log10(1e5)]])

#April 21 gates
bin_vals_16_4_21 = np.array([[1.4775788014225415,2.200826012340419],
[2.216443456408671,2.8355123402887603],
[2.8355123402887603,3.474446976276228],
[3.4849588153986915,np.log10(3e4)]])

aff_gates = [bin_vals_16_4_15, bin_vals_16_4_19, bin_vals_16_4_21]

exp_gate_16_4_15 = np.array([[1.4917982838529373, 3.337952074987111],
[3.358239479285288, 3.80681652987833],
[3.8214685440936806, 4.272299750719853],
[4.275680984769549, 5.]])

exp_gate_16_4_19 = np.array([[1.4917982838529373, 2.5027872647121274],
[2.523074669010305, 3.352604089202462],
[3.3672561034178123, 4.242995722289152],
[4.246376956338848, 5.345278022490142]])

exp_gate_16_4_21 = np.array([[1.4917982838529373,3.411212146063864],
[3.416847536146691,3.80681652987833],
[3.8214685440936806,4.242995722289152],
[4.246376956338848,5]])

exp_gates = [exp_gate_16_4_15, exp_gate_16_4_19, exp_gate_16_4_21]

for out_name, rep, sort_name, rep_number, directory, aff_gate, exp_gate in zip(out_names, reps, sort_counts, range(1,4), directories, aff_gates, exp_gates):
    # Needed for proper focusing
    plt.ion()
    plt.close('all')

    # Create figure with subplots and specified spacing
    figsize=(3.5,5.6)
    rows=14
    cols=1
    col = 1
    fig, axes = plt.subplots(figsize=figsize)
    gs = gridspec.GridSpec(28, 2) 
    plt.subplots_adjust(
        bottom=0.06,
        top=0.95,
        left=0.17,
        right=0.96,
        wspace=0.6,
        hspace=0.0)

    # Make a labler to add labels to subplots
    labeler = Labeler(xpad=.13,ypad=.01,fontsize=10)

    # For CDR1H and CDR3H
    conc_labels = ['$0$'] + \
                  ['$10^{%1.1f}$'%x for x in np.arange(-9.5,-4.5,0.5)]

    file_labels = ['0M'] + \
                  ['10^%1.1fM'%x for x in np.arange(-9.5,-4.5,0.5)]

    #csv_name = 'out.csv'
    filenames = get_filenames(directory)
    names = [re.search('Sort (\d+)', ii) for ii in filenames]
    condition = [n.group(1) for n in names]
        
    # Make plots
    for [filename, well] in zip(filenames, condition):
        well = int(well)
        library_data = readfcs.readfcs(filename)
        data = library_data.rename(columns={'SSC-W': 'SSC-H', 'SSC-H': 'SSC-W', 'FSC-W': 'FSC-H', 'FSC-H': 'FSC-W'})

        data = filter_fsc_ssc(data) # filter by fsc and ssc
        maskeddata = filter_ps(data) # pre-sort filter based on FITC and PE signal
        # If expression bin, compute expression distribution
        if well == 1:
            ax = plt.subplot(gs[26:28,0])
            name = 'expression'
            plot_expression(ax, data, exp_gate)
            labeler.label_subplot(ax,'B')
            #plt.title('expression', \
            #    fontsize=mpl.rcParams['font.size'])
            #data.to_csv('./data/expression.csv')

        else:
            ind = (well-2)*2
            ax = plt.subplot(gs[ind:(ind+2),0])
            name = 'bin %d'%well
            plot_affinity(ax, maskeddata, aff_gate)

            conc_label = conc_labels[well-2]
            ax.set_ylabel(conc_label,rotation=0,ha='right', fontsize=mpl.rcParams['font.size'])
            #maskeddata.to_csv('./data/'+file_labels[well-2]+'.csv')
            if well!=12:
                ax.set_xlabel('')
                ax.set_xticks([])
            
            if well == 7:
                
                ax.text(-0.35,0.5,'fluorescein [M]', rotation=90,ha='right', fontsize=mpl.rcParams['font.size'],transform=ax.transAxes)
            if well==2:
                labeler.label_subplot(ax,'A')
                #plt.title('affinity', \
                #    fontsize=mpl.rcParams['font.size'])



        #plt.subplot(gs[12,0]).set_visible(False)
        #plt.subplot(gs[11,0]).set_visible(False)
    labeler = Labeler(xpad=.12,ypad=.01,fontsize=10)

    vmin = 1e3
    vmax = 1E7
    #panel C
    ax = plt.subplot(gs[0:12,1])
    im = count_cells(ax,sort_name,'cells', \
        ylabel=True, vmax=vmax, vmin=vmin)
    labeler.label_subplot(ax,'C')
    
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1/30.)
    pad = axes_size.Fraction(0.75, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    
    #cax = fig.add_axes([0.92, 0.5375, 0.025, 0.4125])
    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    #cbar.set_label(r'number of sorted cells')
    cbar.solids.set_rasterized(True)

    # Panel D
    ax = plt.subplot(gs[16:28,1])
    im = count_reads(ax, rep,'reads' ,ylabel=True, vmax=vmax, vmin=vmin)
    labeler.label_subplot(ax,'D')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1/30.)
    pad = axes_size.Fraction(0.75, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    
    #cax = fig.add_axes([0.92, 0.06, 0.025, 0.38])
    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    #cbar.set_label(r'number of sequence reads')
    cbar.solids.set_rasterized(True)
    plt.show()
    plt.savefig(out_name)
    