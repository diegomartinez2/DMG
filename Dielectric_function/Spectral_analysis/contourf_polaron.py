"""
=================
contourf(X, Y, Z)
=================

See `~matplotlib.axes.Axes.contourf`.
"""
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('_mpl-gallery-nogrid')
# read data
namefile = "xae"
with open(namefile) as file:
     lines = [line.rsplit() for line in file]
#Qx = np.zeros(5001)
Qy = np.zeros(51)
Frequency = np.zeros(5001)
Spec = np.zeros(255051)
for i in range(51):
    Qy[i]=lines[i][1]
for i in range(5001):
    Frequency[i]=lines[i][2]
for i in range(255051):
    Spec[i]=lines[i][6]

#for i in range(255051):
#     print (lines[i][2])
#  if (lines[i][0]==0.04346):
     #Qx[i]=lines[i][0]
     #Qy[i]=lines[i][1]
     #Frequency[i]=lines[i][2]
     #Spec[i]=lines[i][6]
# make data
#X, Y = np.meshgrid(np.linspace(-3, 3, 256), np.linspace(-3, 3, 256))
#Z = (1 - X/2 + X**5 + Y**3) * np.exp(-X**2 - Y**2)
X, Y = np.meshgrid(Qy, Frequency)
Z = np.reshape(Spec, (5001,51))
data = np.resize(Spec,(int(len(X)/51),51))
print (Z, Z.min(), Z.max())
levels = np.linspace(Z.min(), Z.max(), 70000)
#levels = np.linspace(0.0, 0.3, 7)
print (levels)
# plot
fig, ax = plt.subplots()

#ax.contourf(X, Y, Z, levels=levels)
#ax.contourf(Z , levels=levels)
ax.imshow(Z, aspect= 'auto',cmap ='Greens', 
		vmin = 0.0, vmax = 0.01,
#		vmin = Z.min(), vmax = Z.max(),
#                 extent =[X.min(), X.max(), Y.min(), Y.max()],
                    interpolation ='nearest', origin ='lower')

plt.show()
fig, ax1 = plt.subplots(1,1)
cax = ax1.imshow(data.T, 
	vmin = 0.0 , vmax = 0.004,
	cmap=plt.colormaps['jet'], origin='lower')
ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

cbar = fig.colorbar(cax)
plt.show()
