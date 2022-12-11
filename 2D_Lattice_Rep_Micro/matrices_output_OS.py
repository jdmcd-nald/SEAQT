#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:04:08 2020
plt.imshow support
Interpolation:
Supported values are 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom',
 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'.
@author: alejomonbar
"""
import matplotlib.pyplot as plt



import numpy as np
import os
import sys
from matplotlib.colors import ListedColormap
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

viridis = cm.get_cmap('viridis', 10)
newcolors = viridis(np.linspace(0, 1, 10))

#pink = np.array([1, 1, 1, 1])
#pink = np.array([170/256, 170/256, 170/256, 1])
pink = np.array([.281412, .155834, .469201, 1])
newcolors[:1, :] = pink
newcmp = ListedColormap(newcolors)


def plot_examples(cms):
    """
    helper function to plot two colormaps
    """
    np.random.seed(19680801)
    data = np.random.randn(30, 30)

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    for [ax, cmap] in zip(axs, cms):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()

print(newcolors)

def colored(matrix):
    new_matrix = 1 * matrix
    elem = Counter(list(matrix.reshape(-1)))
    sort = sorted(elem.keys(), key = lambda x:elem[x], reverse = True)
    i = 9
    for s in sort:
        if s != 0:
            new_matrix[matrix == s] = i
    return new_matrix



import platform
is_windows = any(platform.win32_ver());
if is_windows:
    delim="\\"
else:
    delim="/"

path = os.getcwd() +delim;
files = os.listdir(path) #Save the data in the folder Data
# matrices = []
for file in files:
    if file[1]=="x":
	
        matrix = np.loadtxt(path + file)
        colored_matrix = colored(matrix)
        # matrices.append(matrices)
        fig, ax = plt.subplots(1,1)
        cmap = ListedColormap(['#AAAAAA','#481567' ,'#453781', '#39568c', '#2d708e', '#238a8d', '#20a387' ,'#3cbb75','#73d055','#B8DE29'])
        ccmap=ListedColormap(newcolors)
        cs = ax.matshow(colored_matrix, cmap=ccmap) # See above to produce other styles
        fig.colorbar(cs)
        fig.savefig(os.getcwd() +delim+"img"+file.replace(".", "_") + ".png")
        plt.close(fig)
    else:
        yy=0
