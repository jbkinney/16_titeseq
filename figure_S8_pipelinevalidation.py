import numpy as np
from labeler import Labeler
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb 
import pandas

plt.ion()
plt.close('all')

simulation = dict()
simulation[10000] = pandas.read_csv('data/ideal_0_10000.0.csv')
simulation[1] = pandas.read_csv('data/ideal_0_1.0.csv')
simulation[0.001] = pandas.read_csv('data/ideal_0_0.001.csv')

plot_lims = [-10,-4.5]

figsize=(4,3)
fig = plt.figure(figsize = (4,4))
ax = fig.add_axes([0.2,0.2,0.7,0.7])
boundaries = np.linspace(plot_lims[0],plot_lims[1], 60)
colors = ['#00FF00', '#008F00','#000000']
for ii, sample in enumerate([0.001,1,10000]):
    KD = np.array(simulation[sample]['true_KD'])
    fit = np.array(simulation[sample]['fit_KD'])
    x = []
    y = []
    y_std = []
    for jj in range(boundaries.shape[0]-1):
        usethis = (KD<(10**boundaries[jj+1]))&(KD>(10**boundaries[jj]))
        x.append(np.mean(np.log10(KD[usethis])))
        y.append(np.mean(np.log10(fit[usethis])))
        y_std.append(np.max([0.05,np.std(np.log10(fit[usethis]))]))
    x = np.array(x)
    y = np.array(y)
    y_std = np.array(y_std)
    ax.fill_between(x, y-y_std, y+y_std, facecolor=colors[ii],linewidth=0)
    
ax.set_xticks(range(-10,-3))
ax.set_xticklabels([r'$10^{'+str(ii)+'}$' for ii in range(-10,-3)])
ax.set_yticks(range(-10,-3))
ax.set_yticklabels([r'$10^{'+str(ii)+'}$' for ii in range(-10,-3)])
plt.xlabel(r'$K_D$, true')
plt.ylabel(r'$K_D$, fit')
ax.set_aspect(1.)
ax.set_xlim(plot_lims)
ax.set_ylim(plot_lims)
plt.plot(plot_lims,plot_lims,'--w')
plt.axhline(-9.5,linestyle=':',color='k')
plt.axhline(-5.0,linestyle=':',color='k')
plt.show()
plt.savefig('pdfs/figure_S8_pipelinevalidation.pdf')
#plt.close()