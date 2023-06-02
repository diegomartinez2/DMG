import matplotlib.pyplot as plt
import numpy as np

plt.style.use('_mpl-gallery-nogrid')
# read data
namefile = "xbw"
with open(namefile) as file:
     lines = [line.rsplit() for line in file]
Spec = np.zeros(255051)
for i in range(255051):
     Spec[i]=lines[i][6]
data = np.resize(Spec,(51,int(len(Spec)/51)))

fig, ax1 = plt.subplots(1,1)
cax = ax1.imshow(data.T,
#	vmin = 0.0 , vmax = 0.004,
#	vmin = 0.0 , vmax = 0.3,
	cmap=plt.colormaps['jet'], origin='lower',
	interpolation='gaussian', aspect='auto')
ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

cbar = fig.colorbar(cax)
plt.show()
