from __future__ import division
import pdb
import pylab
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Various matplotlib parameters
mpl.rcParams['legend.scatterpoints'] = 1
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'cmGeneral'

# Class used to label subplots
class Labeler:
    def __init__(self,xpad=.1,ypad=.1,fontsize=10):
        self.xpad = xpad
        self.ypad = ypad
        self.fontsize = fontsize

    def label_subplot(self,ax,label):
        bounds = ax.get_position().bounds
        fig = ax.get_figure()
        y = bounds[1] + bounds[3] + self.ypad
        x = bounds[0] - self.xpad
        fig.text(x,y,label,fontsize=self.fontsize,ha='right',va='bottom')

# # Figure container class
# class Figure:

#     def __init__(self, width, height, rows, cols, 
#         bottom=0.1, 
#         top=0.9, 
#         left=0.1,
#         right=0.9,
#         hspace=0.2,
#         wspace=0.2,
#         fontsize=7,
#         label_fontsize=10,
#         units='in',
#         label_offset_x=0, 
#         label_offset_y=0
#         ):

#         # Set fontsize
#         pylab.rcParams['font.size'] = fontsize

#         # Create figure of specified size
#         if units=='in':
#             figsize = (width, height)
#         elif units=='mm':
#             figsize = (width*0.0393701, height*0.0393701)
#         elif units=='cm':
#             figsize = (width*0.393701, height*0.393701)
#         else:
#             assert False, 'Unknown length units for figure dimensions'

#         # Create figure of specified size
#         self.fig =  plt.figure(figsize=figsize)

#         # Specify placement params
#         self.rows = rows
#         self.cols = cols
#         self.bottom = bottom
#         self.top = top
#         self.left = left
#         self.right = right
#         self.hspace = hspace
#         self.wspace = wspace
#         self.fontsize = fontsize
#         self.label_offset_x = label_offset_x
#         self.label_offset_y = label_offset_y
#         self.label_fontsize = label_fontsize

#         self.figwidth = self.right - self.left
#         self.axwidth = self.figwidth/(self.cols + (self.cols-1)*self.wspace)
#         self.figheight = self.top - self.bottom
#         self.axheight = self.figheight/(self.rows + (self.rows-1)*self.hspace)

#         # Adjust subplot spacing
#         plt.subplots_adjust(
#             bottom=self.bottom,
#             top=self.top,
#             left=self.left,
#             right=self.right,
#             wspace=self.wspace,
#             hspace=self.hspace)

#     # Creates a panel in the given figure
#     def make_panel(self, label, row, col, rowspan=1, colspan=1, \
#         axonly=True, shape=None):
#         panel = Panel(self, label, row, col, rowspan=rowspan, 
#             colspan=colspan, shape=shape)
#         if axonly:
#             return panel.ax
#         else:
#             return panel

# # Panel container class
# class Panel:

#     def __init__(self, figure, label, row, col, 
#         rowspan=1, 
#         colspan=1,
#         shape=None
#         ): 

#         self.fig = figure
#         self.label = label
#         self.row = row
#         self.col = col
#         self.rowspan = rowspan
#         self.colspan = colspan

#         if shape:
#             assert len(shape)==2, \
#                 'Shape argument has length %d, not length 2'%len(shape)
#         else:
#             shape = (self.fig.rows, self.fig.cols)

#         self.label_pos_x = figure.left + \
#             self.col*figure.axwidth*(1+figure.wspace) + figure.label_offset_x
#         self.label_pos_y = figure.top - \
#             self.row*figure.axheight*(1+figure.hspace) + figure.label_offset_y
#         self.fontsize = self.fig.fontsize

#         plt.figtext(self.label_pos_x, self.label_pos_y, self.label,
#             fontsize=figure.label_fontsize)

#         loc = (self.row, self.col)
#         self.ax = plt.subplot2grid(shape, loc, \
#             rowspan=self.rowspan, colspan=self.colspan)
