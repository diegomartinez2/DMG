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
Qy = np.zeros(5001)
Frequency = np.zeros(5001)
Spec = np.zeros(255051)
for i in range(5001):
    Qy[i]=lines[i][1]
    Frequency[i]=lines[i][2]
for j in range(255051):
    Spec[i]=lines[i][6]
Z    np.reshape(A, (-1, 2))
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
z = np.reshape(Spec, (5001,5001))
levels = np.linspace(Z.min(), Z.max(), 7)

# plot
fig, ax = plt.subplots()

ax.contourf(X, Y, Z, levels=levels)

plt.show()
