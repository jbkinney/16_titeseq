#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas
import scipy.stats as stats
import numpy as np
from matplotlib.ticker import MaxNLocator
import pdb
from helper import *

from labeler import Labeler

# Needed for proper focusing
plt.ion()
plt.close('all')

[seq_hash, seq, seq_cdr] = make_Sequence_Hash(
    'data/CDR_library_July_5_2013_sequences.txt')

#CDR1 and 3 concatenated
CDRs = 'TFSDYWMNWVGSYYGMDYWG'

#Raw pdb files
CDR_atoms = pandas.read_csv('./data/CDR13_positions.txt',  delimiter=r"\s+", header=None)
all_atoms = pandas.read_csv('./data/all_aas.txt',  delimiter=r"\s+", header=None)

#Data.frame files of just the alpha carbons and fluorescein molecules
usethis = all_atoms[2] == 'CA'
allalphas = all_atoms[usethis]
usethis = CDR_atoms[2] == 'CA'
CDRalphas = CDR_atoms[usethis]
fluor = CDR_atoms[3] == 'FLU'
justfluor = CDR_atoms[fluor]

#Keys for matching CDR to all alpha carbon data set
all_alpha_keys = allalphas[1].tolist()
CDR_alpha_keys = CDRalphas[1].tolist()

#coordinates of sub-data sets
CDRalphacoords = np.array([CDRalphas[6].tolist(), CDRalphas[7].tolist(), CDRalphas[8].tolist()])
allCDRcoords = np.array([CDR_atoms[6].tolist(), CDR_atoms[7].tolist(), CDR_atoms[8].tolist()])
allcoords = np.array([all_atoms[6].tolist(), all_atoms[7].tolist(), all_atoms[8].tolist()])
allalphacoords = np.array([allalphas[6].tolist(), allalphas[7].tolist(), allalphas[8].tolist()])
fluor_coords = np.array([justfluor[5].tolist(), justfluor[6].tolist(), justfluor[7].tolist()])

#distance between CDR alpha carbons and fluorescein
fluor_dist = np.array([np.sqrt(((aa - fluor_coords.T)**2).sum(axis=1)).min() for aa in CDRalphacoords.T])

contact_distance = 4
fluor_contact = np.array([(np.sqrt(((aa - fluor_coords.T)**2).sum(axis=1))).min() for aa in allCDRcoords.T])
usethis = fluor_contact < contact_distance
contact = np.where(usethis)[0]
usethis = (CDR_atoms[5][usethis])
contact = np.unique([np.where(CDRalphas[5] == ind)[0][0] for ind in usethis])
fluor_contact =np.zeros(20)
fluor_contact[contact] = 1

#distances between CDR alpha carbons and all atoms
protein_dist = np.array([np.sqrt(((CDRalphacoords.T - aa)**2).sum(axis=1)).min() for aa in allcoords.T])

#get all connections to the heavy chain then light chain
def count_all_connections(connect_distance):
    connection = [[] for jj in range(20)]
    extra_connection = [[] for jj in range(20)]
    intra_connection = [[] for jj in range(20)]

    usethis = protein_dist < connect_distance
    usethis &= np.array(all_atoms[4] == 'H')
    original = np.where(usethis)[0]
    usethis = (all_atoms[5][usethis])
    adjacentH = [np.where((allalphas[5] == ind) & (allalphas[4] == 'H'))[0][0] for ind in usethis]

    for nearest, ii in zip(original, adjacentH):
        temp = np.sqrt(((allcoords.T[nearest] - CDRalphacoords.T)**2).sum(axis=1))
        temp = np.where(temp < connect_distance)[0]
        temp = np.unique(np.array(CDRalphas[5])[temp])
        temp = [np.where(CDRalphas[5] == ind)[0][0] for ind in temp]
        for jj in temp:
            if all_alpha_keys[ii] != CDR_alpha_keys[jj]:
                connection[jj].append(ii)

            subset = range(10 * (int(jj)/10), 10 * (int(jj)/10 + 1))
            subset_keys = [CDR_alpha_keys[thisCDR] for thisCDR in subset]
            if all_alpha_keys[ii] not in subset_keys:
                extra_connection[jj].append(ii)

            if all_alpha_keys[ii] in subset_keys and all_alpha_keys[ii] != CDR_alpha_keys[jj]:
                intra_connection[jj].append(ii)

    usethis = protein_dist < connect_distance
    usethis &= np.array(all_atoms[4] == 'L')
    original = np.where(usethis)[0]
    usethis = (all_atoms[5][usethis])
    adjacentL = [np.where((allalphas[5] == ind) & (allalphas[4] == 'L'))[0][0] for ind in usethis]

    for nearest, ii in zip(original, adjacentL):
        temp = np.sqrt(((allcoords.T[nearest] - CDRalphacoords.T)**2).sum(axis=1))
        temp = np.where(temp < connect_distance)[0]
        temp = np.unique(np.array(CDRalphas[5])[temp])
        temp = [np.where(CDRalphas[5] == ind)[0][0] for ind in temp]
        for jj in temp:
            if all_alpha_keys[ii] != CDR_alpha_keys[jj]:
                connection[jj].append(ii)

            subset = range(10 * (int(jj)/10), 10 * (int(jj)/10 + 1))
            subset_keys = [CDR_alpha_keys[thisCDR] for thisCDR in subset]
            if all_alpha_keys[ii] not in subset_keys:
                extra_connection[jj].append(ii)

            if all_alpha_keys[ii] in subset_keys and all_alpha_keys[ii] != CDR_alpha_keys[jj]:
                intra_connection[jj].append(ii)

    connection = [list(set(cs)) for cs in connection]
    extra_connection = [list(set(cs)) for cs in extra_connection]
    intra_connection = [list(set(cs)) for cs in intra_connection]
    return connection, extra_connection, intra_connection, adjacentH, adjacentL


connect_distance = 6
connection, extra_connection, intra_connection, adjacentH, adjacentL = count_all_connections(connect_distance)


def plotLetterAdjacency_manuscript(ax, x, blue, red, filename, xax_label, yax_label, connections):
    label_connections = r'Number contacts '
    numconnects = [len(connections[ii]) for ii in range(20)]
    
    green = np.zeros(x.shape)
    colors = np.array([red/red.max(), green, blue/blue.max()]).T
    ax.scatter(numconnects, x, c=colors, s=20, alpha=1, \
        edgecolors='k',lw=0.5)
    plt.xlabel(xax_label,labelpad=2)
    plt.ylabel(yax_label,labelpad=2)
    plt.xlim([0,15])
    ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(5, integer=True))

    # Label with correlation coefficient
    R, P = stats.pearsonr(numconnects,x)
    ax.set_title('$R^2=%0.2f$,  $P=%0.3f$'%(R**2,P), \
        fontsize=mpl.rcParams['font.size'])


def make_2D_colorbar(ax, x, y):
    [xx, yy] = np.meshgrid(np.linspace(0, 1, 101), np.linspace(0, 1, 101))
    zz = yy * 0
    R = np.zeros([xx.shape[0], xx.shape[1], 3], dtype="d")
    R[:,:,0] = xx
    R[:,:,1] = zz
    R[:,:,2] = yy
    ax.imshow(R)
    plt.gca().invert_yaxis()
    ax.set_xticks([0,101])
    ax.set_yticks([0,101])
    p=int(np.floor(np.log10(y.max())))
    xstr=r'$10^{'+str(np.round(float(y.max())/10**p,1))+'}$'

    ax.set_xticklabels(['low','hi'])
    ax.set_yticklabels(['low','hi'],rotation=0)
    ax.set_ylabel(r'$S_E$', labelpad=-2)
    ax.set_xlabel(r'$S_K$', labelpad=-2)
    ax.tick_params(axis=u'both', which=u'both',length=0)



aff_fun = lambda x, ind: np.nanmean(np.log10(np.array(x['fit_KD'])[ind]))
exp_fun = lambda x, ind: np.nanmean((np.array(x['expression'])[ind]))
aas = np.array(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])
by_group = [5,0,17,7,9,10,4,19,18,15,16,11,13,1,12,6,8,14,2,3]
aa_map = {aas[int_val]:ii for ii, int_val in zip(range(20), by_group)}
cdr1_wtseq = 'TFSDYWMNWV'
cdr3_wtseq = 'GSYYGMDYWG'
wt_seq = cdr1_wtseq+cdr3_wtseq

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
    

def plot_distance(ax, x, y, blue, red):
    green = np.zeros(x.shape)
    colors = np.array([red/red.max(), green, blue/blue.max()]).T
    ax.scatter(y, x, c=colors, s=20, zorder=21, alpha=1, lw=0.5)
    ax.set_xlabel(r'distance to antigen ($\mathrm{\AA}$)',labelpad=2)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(4, integer=True))

    # Label with correlation coefficient
    R, P = stats.pearsonr(x,y)
    ax.set_title('$R^2=%0.2f$,  $P=%0.3f$'%(R**2,P),\
        fontsize=mpl.rcParams['font.size'])



# Load data sets
rep1 = pandas.read_csv('data/replicate_1.csv')
rep2 = pandas.read_csv('data/replicate_2.csv')
rep3 = pandas.read_csv('data/replicate_3.csv')

all_reps = [rep1, rep2, rep3]

# Create figure with subplots and specified spacing
fontsize=7
figsize=(3.42,6)
rows = 4
cols = 2    
fig = plt.figure(figsize=figsize)
plt.subplots_adjust(
    bottom = 0.07,
    top = 0.95,
    left = 0.15,
    right = 0.95,
    hspace = 0.6,
    wspace = 0.8)

# Make a labler to add labels to subplots
labeler = Labeler(xpad=.07,ypad=0.01,fontsize=10)

SE = []
SA = []
for rep in all_reps:
    temp_A, wt_A = \
        c_matrix(rep, aff_fun)
    temp_A = np.sqrt(np.mean((temp_A-wt_A)**2, axis = 1))
    
    temp_E, wt_E = \
        c_matrix(rep, exp_fun)
    temp_E = np.sqrt(np.mean((temp_E-wt_E)**2, axis = 1))
    
    SE.append(temp_E)
    SA.append(temp_A)

SE = np.nanmedian(np.array(SE), axis=0)
SA = np.nanmedian(np.array(SA), axis=0)

#Figure 6A
ax = plt.subplot2grid((rows,cols), (0,0))
labeler.label_subplot(ax,'A')
plt.axis('off')

pymol = open('./structure/pymol_color.py','w')
pymol.write('from pymol import cmd, stored\n')
positions = range(28, 38) + range(100, 110)
for ii in range(len(positions)):
    red = str(SA[ii]/SA.max())
    blue = str(SE[ii]/SE.max())
    pymol.write('cmd.set_color(\'mycolor'+str(ii)+'\',['+red+',0,'+blue+'])\n')
    pymol.write('cmd.color(\'mycolor'+str(ii)+'\', \'chain \\\'H\\\' and i. '+str(positions[ii]) + '\')\n')

pymol.close()
ax = plt.subplot2grid((6,3), (2,2))
make_2D_colorbar(ax, SE, SA)

connect_distance = 6
connection, extra_connection, intra_connection, adjacentH, adjacentL = count_all_connections(connect_distance)


# Panel B: sigma_K vs. num contacts
ax = plt.subplot2grid((rows,cols), (2,0))
labeler.label_subplot(ax,'B')
plotLetterAdjacency_manuscript(ax, SA, SE, SA, '6D', 'number of contacts', r'$S_K$', connection)

# Panel C: sigma_K vs. distance to antigen
ax = plt.subplot2grid((rows,cols), (2,1))
labeler.label_subplot(ax,'C')
plot_distance(ax, x=SA, y=fluor_dist, blue=SE, red=SA)
ax.set_ylabel(r'$S_K$')

# Panel A: sigma_E vs. num contacts
ax = plt.subplot2grid((rows,cols), (3,0))
labeler.label_subplot(ax,'D')
plotLetterAdjacency_manuscript(ax, SE, SE, SA, '6D', 'number of contacts', r'$S_E$', connection)

# Panel B: sigma_E vs. distance to antigen
ax = plt.subplot2grid((rows,cols), (3,1))
labeler.label_subplot(ax,'E')
plot_distance(ax, x=SE, y=fluor_dist, blue=SE, red=SA)
ax.set_ylabel(r'$S_E$')


plt.show()
plt.savefig('pdfs/figure_6_structure.pdf')
#plt.close()


