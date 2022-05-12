import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc
from scipy.stats import norm

# Fixing random state for reproducibility
#np.random.seed(19680801)
def  random_data():
    # some random data
    x = 2.0*np.random.rand(500)-1.0
    y = 2.0*np.random.rand(500)-1.0
    return x,y
def Sobol_1d_serial():
    #Sobol 1d random data
    sampler = qmc.Sobol(d=1, scramble=False)
    sample = sampler.random_base2(m=3)
    x = []
    y = []
    while (len(x)<1000):
        x.append(2.0*sampler.random()[0][0]-1.0)
        y.append(2.0*sampler.random()[0][0]-1.0)
    return x,y
def Sobol_1d():
    #Sobol 1d random data
    sampler = qmc.Sobol(d=1, scramble=False)
    sample = sampler.random_base2(m=3)
    #x = []
    #y = []
    x = 2.0*sampler.random(1000)[:,0]-1.0
    y = 2.0*sampler.random(1000)[:,0]-1.0
    return x,y
def Pure_Sobol_2d(size):
    #Pure Sobol 2d random data
    sampler = qmc.Sobol(d=2, scramble=False)
    sample = sampler.random_base2(m=size)
    #x = []
    #y = []
    data1 = sampler.random(1000)
    x = 2.0*data1[:,0]-1.0
    y = 2.0*data1[:,1]-1.0
    return x,y
def Sobol_2d():
    #Sobol 2d alt random data
    sampler = qmc.Sobol(d=2, scramble=False)
    sample = sampler.random_base2(m=3)
    #x = []
    #y = []
    x = 2.0*sampler.random(500)[:,0]-1.0
    y = 2.0*sampler.random(500)[:,1]-1.0
    return x,y
def random_sym_data():
    # random u, -u data
    #x = []
    #y = []
    x = 2.0*np.random.rand(250)-1.0
    y = 2.0*np.random.rand(250)-1.0
    x=np.append(x,-x)
    y=np.append(y,-y)
    return x,y
def Sobol_sym_data(size):
    # Sobol 2d alt u , -u data
    sampler = qmc.Sobol(d=2, scramble=False)
    sample = sampler.random_base2(m=size)
    #x = []
    #y = []
    x = 2.0*sampler.random(250)[:,0]-1.0
    y = 2.0*sampler.random(250)[:,1]-1.0
    x=np.append(x,-x)
    y=np.append(y,-y)
    return x,y
def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')

#x,y = Sobol_sym_data(4)
#x,y = Pure_Sobol_2d(3)
x,y = Sobol_2d()

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# use the previously defined function
scatter_hist(x, y, ax, ax_histx, ax_histy)

plt.show()
## start with a square Figure
#fig = plt.figure(figsize=(8, 8))
#
## Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
## the size of the marginal axes and the main axes in both directions.
## Also adjust the subplot parameters for a square plot.
#gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
#                      left=0.1, right=0.9, bottom=0.1, top=0.9,
#                      wspace=0.05, hspace=0.05)
#
#ax = fig.add_subplot(gs[1, 0])
#ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
#ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
#
## use the previously defined function
#scatter_hist(x, y, ax, ax_histx, ax_histy)
#
#plt.show()
