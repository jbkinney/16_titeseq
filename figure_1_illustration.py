#!/usr/bin/env python
from __future__ import division
import pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import scipy as sp

from labeler import Labeler

def log_hill_function(x,A,B,log10_K): 
    return np.log10(B+A*(x/(x+10.0**log10_K)))

# Hill function PDF
# noisiness = std of noise
def sample_hill_pdf(N,conc,amp=700,kd=1E-5,bg=10,noisiness=1):
    log_normal_noise = np.exp(noisiness*np.random.randn(N))
    return (bg + amp*conc/(conc + kd))*log_normal_noise

# Class to hold clone information
class Clone:

    def __init__(self, kd, amp, color):
        self.kd = kd
        self.amp = amp
        self.color = color
        self.noisiness = 1
        self.bg = 10

    def sample_pdf(self,N,conc):
        return sample_hill_pdf(N, conc, amp=self.amp,kd=self.kd, \
                               bg=self.bg, noisiness=self.noisiness)

    def median_signal(self,concs):
        return sample_hill_pdf(N=1, conc=concs, amp=self.amp, kd=self.kd, \
                               bg=self.bg, noisiness=0)

    def pdf_values(self,flouesence_grid,conc):
        mu = self.median_signal(conc)
        log_mu = np.log(mu)
        log_fls = np.log(flouesence_grid)
        pdf_vals = (1/np.sqrt(2*np.pi*self.noisiness**2)) * \
             np.exp(-((log_fls-log_mu)**2)/(2.0*self.noisiness**2))
        return pdf_vals

# Needed for proper focusing
plt.ion()
plt.close('all')

# Define colors
red = [0.9,0.3,0.3]
blue = [0.4,0.0,1]
gray = [0.4,0.4,0.4]
lightgray = [0.7,0.7,0.7]
black = [0.,0.,0.]
cyan = np.array([86., 201., 236.])/256.
orange = np.array([243., 156., 79.])/256.

print 'Simulating data...'

# Perform experiments on these clones
clone1 = Clone(kd=1.2E-9, amp=300, color=cyan)
clone2 = Clone(kd=1E-6, amp=1000, color=orange)
clones = [clone1, clone2]
num_clones = len(clones)

# Perform simulations for clones
concs = np.hstack((0, np.logspace(start=-9.5,stop=-5,num=10,base=10,endpoint=True)))
conc_grid = np.hstack((0,np.logspace(start=-10,stop=-2,num=1000,base=10,endpoint=True)))
conc_ticks = np.logspace(start=-9,stop=-5,num=3,base=10,endpoint=True)
#concs = 10.0**(np.arange(-8,-3.5,.5))
num_concs = len(concs)
conc_lim = [min(concs)/10, max(concs)*10]

# Sort into bins defined by these intervals
bin_edges = [0,30,300,np.Inf]
bins = [(bin_edges[i],bin_edges[i+1]) for i in range(len(bin_edges)-1)]
num_bins = len(bins)
fl_grid = 10.0**(np.arange(0,4.1,.01))

# Other relevant numbers
num_samples = 1000
num_replicates = 1

### data, binned_data ### 
# Dimension 0: samples
# Dimension 1: replicates
# Dimension 2: concentrations
# Dimension 3: clones
data = np.zeros([num_samples, num_replicates, num_concs, num_clones])
binned_data = np.zeros(data.shape)
#pdb.set_trace()

### hist_data ### 
# Dimension 0: bins
# Dimension 1: replicates
# Dimension 2: concentrations
# Dimension 3: clones
hist_data = np.zeros([num_bins, num_replicates, num_concs, num_clones])

### mean_bin_data ### 
# Dimension 0: replicates
# Dimension 1: concentrations
# Dimension 2: clones
mean_bin_data = np.zeros([num_replicates, num_concs, num_clones])


# Simulate all data
for rep_num in range(num_replicates):
    for clone_num, clone in enumerate(clones):
        for conc_num, conc in enumerate(concs):
            data[:, rep_num, conc_num, clone_num] = \
                clone.sample_pdf(num_samples,conc)

# Bin data
for bin_num, bin in enumerate(bins):
    indices = (data >= bin[0]) & (data < bin[1])
    binned_data[indices] = bin_num

# Compute histogram and mean bin numbers
for rep_num in range(num_replicates):
    for clone_num, clone in enumerate(clones):
        for conc_num, conc in enumerate(concs):

            # Histogram the data
            for bin_num in range(num_bins):
                hist_data[bin_num,rep_num,conc_num,clone_num] = \
                    sum(binned_data[:,rep_num,conc_num,clone_num] == bin_num)

            # Comptue mean bin number
            mean_bin_data[rep_num,conc_num,clone_num] = \
              sum(np.arange(num_bins)*hist_data[:,rep_num,conc_num,clone_num]) / sum(hist_data[:,rep_num,conc_num,clone_num])

print 'Done!'

row1_labels = 'BC'
row2_labels = 'DE'

# Create figure with subplots and specified spacing
figsize=(3.42,4.5)
rows = 3
cols = 2  
fig, axes = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom = 0.07,
    top = 0.95,
    left = 0.12,
    right = 0.98,
    hspace = 0.2,
    wspace = 0.5)

# Turn axes off on upper left corner plots
for i in range(2):
    for j in range(2):
        plt.setp(axes[i,j],frame_on=False, xticks=[], yticks=[])

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.07,ypad=0.0,fontsize=10)

# Label upper left corner
#ax = axes[0,0]
ax = plt.subplot(451)
labeler.label_subplot(ax,'A')
plt.setp(ax,frame_on=False, xticks=[], yticks=[])

ax = plt.subplot(453)
labeler.label_subplot(ax,'B')
plt.setp(ax,frame_on=False, xticks=[], yticks=[])

ax = plt.subplot(455)
labeler.label_subplot(ax,'C')
plt.setp(ax,frame_on=False, xticks=[], yticks=[])


#### panel D
#ax = plt.subplot2grid((rows,cols), (3,0), rowspan=2, colspan=1)
ax = axes[2,0]
labeler.label_subplot(ax,'D')

# Draw lines
for clone in clones:
    # Get data points
    ys = np.log10(clone.median_signal(concs))
    xs = concs

    # Fit a hill function
    A0 = clone.amp
    B0 = clone.bg
    log10_K0 = np.log10(clone.kd)
    p0 = [A0, B0, log10_K0]
    p, pconv =  sp.optimize.curve_fit(log_hill_function,xs,ys,p0)
    
    # Plot inferred KD
    KD = 10.0**p[2]
    ax.axvline(KD, linestyle='-', c=clone.color,\
            lw=1, alpha=0.5, zorder=0)

    # Plot fitted curve
    fit_ys = 10.0**log_hill_function(conc_grid,p[0],p[1],p[2])
    ax.loglog(conc_grid, fit_ys, '-', zorder=2, c=clone.color, lw=2)

    # Plot data points
    ax.loglog(xs, 10.**ys, '.', markerfacecolor=clone.color, \
        markeredgecolor='k', markeredgewidth=.5, markersize=10, zorder=2)


plt.xscale('symlog', linthreshx=1e-10, linscalex=0.5)
ax.set_xlabel(r'antigen [M]', labelpad=2)
ax.set_ylabel('fluorescence [au]', labelpad=2)
ax.set_xlim(conc_lim)
ax.set_xticks(conc_ticks)
ax.set_ylim([5,2000])

# #### panel E
#ax = plt.subplot2grid((rows,cols), (3,1), rowspan=2, colspan=1)
ax = axes[2,1]
labeler.label_subplot(ax,'E')

# Draw lines
for clone_num, clone in enumerate(clones):
    #ax.axvline(clone.kd, linestyle='-', c=clone.color, \
    #    lw=1, alpha=0.3, zorder=0)
    ys = mean_bin_data[0,:,clone_num]
    xs = concs
    #ax.semilogx(xs, ys, '.', markerfacecolor=clone.color, \
    #    markeredgecolor='k', markeredgewidth=.5, markersize=10, zorder=1)

    # Fit a hill function
    A0 = ys[-1] - ys[0]
    B0 = ys[0]
    log10_K0 = np.log10(clone.kd)
    p0 = [A0, B0, log10_K0]
    p, pconv =  sp.optimize.curve_fit(log_hill_function,xs,ys,p0)
    
    # Plot inferred KD
    KD = 10.0**p[2]
    ax.axvline(KD, linestyle='-', c=clone.color,\
            lw=1, alpha=0.5, zorder=0)

    # Plot fitted curve
    fit_ys = log_hill_function(conc_grid,p[0],p[1],p[2])
    indices = fit_ys > 0
    ax.semilogx(conc_grid[indices], fit_ys[indices], '-',\
        zorder=2, c=clone.color, lw=2)

    # Plot data points
    usethis = (ys > 0)
    ax.semilogx(xs, ys, '.', markerfacecolor=clone.color, \
        markeredgecolor='k', markeredgewidth=.5, markersize=10, zorder=2)
    
    plt.xscale('symlog', linthreshx=1e-10, linscalex=0.5)
        
plt.xscale('symlog', linthreshx=1e-10, linscalex=0.5)
ax.set_xlabel('antigen [M]', labelpad=2)
ax.set_ylabel('mean bin', labelpad=2)
ax.set_ylim([0,num_bins-1])
ax.set_yticks(range(num_bins))
ax.set_xlim(conc_lim)
ax.set_xticks(conc_ticks)

plt.show()
plt.savefig('pdfs/figure_1_illustration.pdf')
