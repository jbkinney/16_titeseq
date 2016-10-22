#!/usr/bin/env python
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import readfcs
from shapely.geometry import MultiPoint, Point
import pandas
import pdb
from scipy.optimize import minimize, curve_fit

from labeler import Labeler

wt_color = [0.8,0.2,0.8]
opt_color = [0,0.8,0]

# Needed for proper focusing
plt.ion()
plt.close('all')
fontsize=7

def in_gate(data, gate, name1, name2):
    x_data = data[name1]
    y_data = data[name2]
    x_data = [np.log10(max(x, 1e-10)) for x in x_data]
    y_data = [np.log10(max(y, 1e-10)) for y in y_data]
    poly = MultiPoint(gate).convex_hull
    usethis = [poly.contains(Point(x, y)) for x, y in zip(x_data, y_data)]
    return usethis

def gate_fsc(data):
    sub_data = lambda indata, gate: indata[in_gate(indata, gate, 'FSC-A','SSC-A')]
    data = sub_data(data,[(3.6,2.9),(3.6, 3.8),(4.4,4.4),(4.4,4)])
    return data

def plot_affinity(ax, data1, data2):
    sub_data = lambda indata, gate: indata[[not ii for ii in in_gate(indata, gate, 'PE-A', 'Alexa Fluor 405-A')]]
    PE1 = np.log10(np.array(data1['PE-A']))
    BV1 = np.log10(np.array(data1['Alexa Fluor 405-A']))
    PE2 = np.log10(np.array(data2['PE-A']))
    BV2 = np.log10(np.array(data2['Alexa Fluor 405-A']))
    gridsize=128
    bounds = np.array([3.5,5.5])
    bins = np.linspace(bounds[0],bounds[1],gridsize)
    
    H1, xedges, yedges = np.histogram2d(PE1.flatten(), \
        BV1.flatten(), bins=[bins, bins])
    H2, xedges, yedges = np.histogram2d(PE2.flatten(), \
        BV2.flatten(), bins=[bins, bins])
    H1/=float(H1.max())
    H2/=float(H2.max())

    def f(x,b):
        return x+b

    # Tighter bounds for regression
    bounds2 = np.array([4.5,5])

    # Do regression for WT: PE1 vs BV1 
    indices = (BV1 > bounds2[0]) & (BV1 < bounds2[1])
    p, pconv = curve_fit(f,BV1[indices],PE1[indices])
    b1 = p[0]
    grid_b1 = b1*gridsize/(bounds[1]-bounds[0])

    # Do regression of OPT: PE2 vs BV2
    indices = (BV2 > bounds2[0]) & (BV2 < bounds2[1])
    p, pconv = curve_fit(f,BV2[indices],PE2[indices])
    b2 = p[0]
    grid_b2 = b2*gridsize/(bounds[1]-bounds[0])

    # Determine ratio of OPT to WT
    opt_over_wt = 10**(b2-b1)
    print 'Functional expression: OPT/WT = %f'%opt_over_wt

    min_color_intensity = 1./2 # sets the minimum color displayed
    denominator = 1./(1-min_color_intensity)
    offset = denominator - 1
    im = np.zeros((H1.shape[0], H1.shape[1], 3), dtype=float)
    usethis = (H1==0)
    H1+=offset
    H1/=np.max(H1)
    H1[usethis]=0
    im[:,:,0] = H1 # Set pJK36 color here, 0 for red, 1 for green, 2 for blue
    im[:,:,2] = H1
    usethis = (H2==0)
    H2+=offset
    H2/=np.max(H2)
    H2[usethis]=0
    im[:,:,1] = H2 # Set pJK37 color here, 0 for red, 1 for green, 2 for blue
    usethis = (H1==0) & (H2==0)
    im[usethis,:]=1
    
    # If there is actually an axes to plot on
    if ax:
        ax.imshow(im, origin='lower', interpolation='nearest')
        ax.plot(-1,-1,'.',label='WT', ms=5, mec='none', c=[0.8,0,0.8]) # Set pJK36 legend color here
        ax.plot(-1,-1,'.',label='OPT', ms=5, mec='none', c=[0.0,0.8,0.]) # Set pJK37 legend color here

        grid_x = np.arange(gridsize)
        ax.plot(grid_x,grid_x + grid_b1,'-', lw=3, alpha=0.5, c=[1,0,1])
        ax.plot(grid_x,grid_x + grid_b2,'-', lw=3, alpha=0.5, c=[0.0,1,0.0])

        ax.set_xlabel(r'BV signal [au]', labelpad=2)
        ax.set_ylabel(r'PE signal [au]', labelpad=2)
        tick = [0.5*gridsize/2,1.5*gridsize/2]
        #tick=[gridsize/4,gridsize/4,2*gridsize/4,3*gridsize/4,gridsize]
        labels = ['$10^4$','$10^5$']
        ax.set_xlim([0,gridsize-1])
        ax.set_ylim([0,gridsize-1])
        #ax.set_xlim([0,])
        #ax.set_ylim([0.5*gridsize/2,1.5*gridsize/2])
        ax.set_xticks(tick)
        ax.set_yticks(tick)
        ax.set_xticklabels(labels)
        ax.set_yticklabels(labels)
        plt.legend(loc='upper left', borderaxespad=0., labelspacing=0.25, \
            handlelength=1,borderpad=0.2, handletextpad=0.1,\
            fontsize=mpl.rcParams['font.size'])

    return opt_over_wt

def KD_fun(x, y, y_std, s, lb, KD0):
    out = [lb, KD0, s]
    y[y<0] = np.nan
    error = np.nanmean((1./y_std**2 * (np.log(y) - np.log(hill(out, x))))**2)
    return out, error

def KD_optimize(params, x, y, y_std, sbounds, bbounds):
    params=impose_bounds(params, sbounds, bbounds)
    y[y<0] = np.nan
    error = np.nanmean((1./y_std**2 * (np.log(y) - np.log(hill(params, x))))**2)
    return error

def impose_bounds(params, sbounds, bbounds):
    params = [np.abs(ps) for ps in params]
    params[0] = np.max([bbounds[0], params[0]])
    params[0] = np.min([bbounds[1], params[0]])
    params[1] = np.max([1e-10, params[1]])
    params[1] = np.min([1, params[1]])
    params[2] = np.max([sbounds[0], params[2]])
    params[2] = np.min([sbounds[1], params[2]])
    return params

def get_KD(fl, y, y_std, a0, b0):
    lbs, amplitudes, KD_try = np.meshgrid(np.hstack((0, np.logspace(0, 2, 21))), np.logspace(-0, 1.0, 21), np.logspace(-10, -1, 20))
    lbs = lbs.flatten()
    amplitudes = amplitudes.flatten()
    KD_try = KD_try.flatten()
    
    y = np.array(y)
    x = np.array(fl)
    lowest_error = np.Inf
    out = [0, 1e-3, 100]
    for s, b, k in zip(amplitudes, lbs, KD_try):
        pfit, error = KD_fun(x, y, y_std, s*a0, b*b0/10., k)
        if error <= lowest_error:
            out = pfit
            lowest_error = error
    sbounds = [np.min(amplitudes*a0), np.max(amplitudes*a0)]
    bbounds = [np.min(lbs*b0/10.), np.max(lbs*b0/10.)]
    SSE = lambda input:KD_optimize(input, x, y, y_std, sbounds, bbounds )
    out = minimize(SSE, out, method='nelder-mead')
    out = np.abs(out['x']).tolist()
    out = impose_bounds(out, sbounds, bbounds)
    return out

def hill(p, F):
    p = np.abs(p)
    amplitude = p[2]
    return p[0]+amplitude * F/(F+p[1])
    
def plot_clone(fl, y, ax, color):
    x = np.array(fl)
    y[1]=np.nan
    pfit = get_KD(x, y, 1, np.nanmean(y[1:]), 0.01)
    
    xl = [1E-10, 3E-4]
    xticks = [1E-10,1E-8,1E-6,1E-4]

    #xl = [1E-9, 1E-4]
    yl = [3E-2,3E0]

    # Plot lines
    xsample = np.logspace(np.log10(xl[0]),np.log10(xl[1]),100)
    #xsample = np.hstack((np.linspace(0,1e-9,100),np.logspace(-9,-4, 100)))
    ysample = hill(pfit, xsample)-pfit[0]
    #ax.plot(xsample, ysample, lw=1, c=color)
    ax.loglog(xsample,ysample, lw=2, c=color)

    # Plot points
    usethis = np.isfinite(y) & (fl > 10**-8.4) & (fl < 10**-4.6)
    ax.scatter(fl[usethis], y[usethis]-y[0], s=20, c=color,\
        zorder=6, edgecolor='k', lw=0.5)
    yval = np.nanmax(pfit[2])/100.
    yprime = (yval)/pfit[2]
    xval = 1e-9
    yval = pfit[2] * xval/(xval + pfit[1])
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks(xticks)

# Create figure with subplots and specified spacing
figsize=(3.42,3.42)
rows = 1
cols = 1
fig, ax = plt.subplots(rows,cols,figsize=figsize)
plt.subplots_adjust(
    bottom=0.1,
    top=0.9,
    left=0.1,
    right=0.98,
    wspace=0.7,
    hspace=0.6)

# # Make a labler to add labels to subplots
# labeler = Labeler(xpad=.11,ypad=-.01,fontsize=10)

# #########################################################################################################################################
# # plots A-B, titration curves
# x = np.array([ 0, 10**-8.5, 10**-8,10**-7.5,10**-7,10**-6.5,10**-6,10**-5.5,10**-5,10**-4.5,10**-4])

# # pJK36 = np.array([[348,331,361,392,407,463,455,444,352,326,382],
# #                   [393,893,671,879,939,1064,1249,1363,1629,2500,5978],
# #                   [402,1010,571,906,912,997,1363,1398,1601,2009,2907]], dtype=float)

# # pJK36[0,-2:] = np.nan #badly washed
# # pJK36_delta = np.array([[366,361,365,357,347,384,336,340,373,371,512],
# #                   [366,383,370,401,416,406,410,419,461,627,2394],
# #                   [398,358,370,403,409,391,439,450,578,631,1375]])


# # #fig = plt.figure()
# # #ax = fig.add_axes([.25, .25, .63, .68])
# # ax = axes[0,0]
# # labeler.label_subplot(ax,'A')
# # base_color = np.array([0.8,0,0])
# # fade_color = np.array([0.,0.,0.])
# # weights = [.5,.75,1.0]
# # for ii in range(pJK36.shape[0]):
# #     color = base_color*weights[ii] + fade_color*(1-weights[ii])
# #     curr = pJK36[ii] - pJK36_delta[ii]
# #     plot_clone(x, curr/float(np.nanmax(curr)-curr[0]), ax, color)

# # ax.set_xlabel('fluorescein [M]', labelpad=2)
# # ax.set_ylabel('(adjusted)\nfluorescence [au]', labelpad=2)
# # ax.set_title('WT', va='top', ha='left', x=0.1, y=.8,\
# #     fontsize=mpl.rcParams['font.size'])

# pJK37 = np.array([[319.47742527754008, 826.26547885824493, 743.55841248944773, 779.25333538606674, 742.20815219318217, 728.62032946410682, 754.57203522585849, 773.93826329387537, 866.86071980777285, 932.54932247461556, 884.72637818895294],
#                   [314.30528211935155, 1253.2489674385663, 1466.9629186874058, 1507.3860751791995, 1469.2618510899479, 1562.1839438962918, 1439.1694621522881, 1484.701376825526, 1658.519762053427, 1840.1490493663687, 2113.7425671049664],
#                   [322.98310330996236, 1080.8712404442526, 922.17729763268744, 952.10362534509022, 902.19171700895231, 895.25051874491851, 928.07857928074475, 929.18094179255149, 997.14574399718367, 1072.5453567505513, 1744.1519278936009],
#                   [309.75845905141813, 2235.3923607258666, 1760.3574156444668, 1809.8532445309879, 1782.2513548115398, 1581.9271369558312, 1752.463864145973, 1700.9305376589064, 1778.1526025423364, 1947.7103984906362, 3557.9623751553777]], dtype=float)
# pJK37_delta = np.array([[312.85295507065098, 318.47298052479607, 321.90696847008485, 314.54888473365611, 298.62190146064432, 321.46160668059213, 440.51948014445003, 315.40842310595588, 366.52435098719275, 422.25391483991109, 969.7149504835586],
#                         [301.41066658471493, 272.76047344280926, 319.30191200918773, 303.52273375222177, 362.94043520335123, 382.63825780463424, 322.67615513628817, 328.21709503562994, 390.28971119340133, 404.17957457262844, 792.71926970647962],
#                         [325.56016905573267, 309.483914520116, 339.19976743868466, 318.47983094382442, 309.82401839459891, 301.99601868775176, 343.86564280083104, 340.38645229127371, 376.54534840893871, 459.56701157871157, 1115.7139118379484],
#                         [352.55010574630501, 310.78023191495021, 339.64355539686323, 364.34508081074318, 353.39250303015353, 379.70235455268687, 371.77230628686118, 400.72535615068625, 415.4202603711214, 1070.9795578764672, 2801.8366031922574]])

# ax = axes[0]
# labeler.label_subplot(ax,'A')
# base_color = np.array([0,0,.8])
# fade_color = np.array([0.,0.,0.])
# weights = [.4,.6,.8,1.0]
# for ii in range(pJK37.shape[0]):
#     color = base_color*weights[ii] + fade_color*(1-weights[ii])
#     curr = pJK37[ii] - pJK37_delta[ii]
#     plot_clone(x, curr/float(np.nanmax(curr)-curr[0]), ax, color)

# ax.set_xlabel('fluorescein [M]', labelpad=2)
# ax.set_ylabel('(adjusted)\nfluorescence [au]', labelpad=2)
# ax.set_title('OPT', va='top', ha='left', x=0.1, y=.8, \
#     fontsize=mpl.rcParams['font.size'])


#########################################################################################################################################

# Flow cytometry data sets
data_infos = [('./data/CDR1and3_Feb23_2015_pJK36_001.fcs',
             './data/CDR1and3_Feb23_2015_pJK37_001_003.fcs'),
             ('./data/CDR1and3_Feb27_2015_pJK36_001.fcs',
             './data/CDR1and3_Feb27_2015_pJK37_001_003.fcs'),
             ('./data/CDR1and3_Mar3_2015_pJK36_001.fcs',
             './data/CDR1and3_Mar3_2015_pJK37_001_003.fcs'),
             ('./data/CDR3_Feb12_2015_pJK36_001.fcs',
             './data/CDR3_Feb12_2015_pJK37_002.fcs')]


# Compute ratio OPT/WT ratio for all four experiments
ratios = []
for i, info in enumerate(data_infos):
    print i
    #ax = plt.subplot2grid((rows,cols), (1,0), colspan=2, rowspan=2)
    pJK36 = readfcs.readfcs(info[0])
    pJK37 = readfcs.readfcs(info[1])
    pJK36 = pJK36.rename(columns={'SSC-W': 'SSC-H', 'SSC-H': 'SSC-W', 'FSC-W': 'FSC-H', 'FSC-H': 'FSC-W'})
    pJK37 = pJK37.rename(columns={'SSC-W': 'SSC-H', 'SSC-H': 'SSC-W', 'FSC-W': 'FSC-H', 'FSC-H': 'FSC-W'})
    pJK36 = gate_fsc(pJK36)
    pJK37 = gate_fsc(pJK37)# First, gate cells by fsc and ssc
    ratio = plot_affinity(False, pJK36, pJK37) #then compare fluorescein and c-myc signals
    ratios.append(ratio)

ratios = np.array(ratios)
print 'OPT/WT functional expression ratio = %0.2f +- %0.2f'%(np.mean(ratios), np.std(ratios)/np.sqrt(len(ratios)))

# Only need to plot a single 2D FACS curve, since replicates show the same thing. 
info = data_infos[1]
#ax = plt.subplot(111) #axes[1]
pJK36 = readfcs.readfcs(info[0])
pJK37 = readfcs.readfcs(info[1])
pJK36 = pJK36.rename(columns={'SSC-W': 'SSC-H', 'SSC-H': 'SSC-W', 'FSC-W': 'FSC-H', 'FSC-H': 'FSC-W'})
pJK37 = pJK37.rename(columns={'SSC-W': 'SSC-H', 'SSC-H': 'SSC-W', 'FSC-W': 'FSC-H', 'FSC-H': 'FSC-W'})
pJK36 = gate_fsc(pJK36)
pJK37 = gate_fsc(pJK37)# First, gate cells by fsc and ssc

#labeler.label_subplot(ax,'CDEF'[ii])
#labeler.label_subplot(ax,'B',xpad_adjust=-0.04)

plot_affinity(ax, pJK36, pJK37) #then compare fluorescein and c-myc signals
ax.set_title('2 $\mu$M fluorescein',fontsize=mpl.rcParams['font.size'])

plt.show()
plt.savefig('./pdfs/figure_S5_facs.pdf')
#plt.close()
