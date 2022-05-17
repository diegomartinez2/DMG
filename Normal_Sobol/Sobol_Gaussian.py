import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt

# #rand = np.random.normal(size = (size, n_modes))
# #sampler = qmc.Sobol(d=n_modes, scramble=True)
# size=500
# sampler = qmc.Sobol(d=1, scramble=True)
# size_sobol = int(np.log(size)/np.log(2))+1
# sample = sampler.random_base2(m=size_sobol)
#
# data=[]
# while (len(data)<size):
#  data1 = sampler.random()
#  data2 = sampler.random()
#  v1=2.0*data1-1.0
#  v2=2.0*data2-1.0
# # Riq=np.linalg.norm(v1)**2+np.linalg.norm(v2)**2
#  Riq=v1**2+v2**2
# # Riq_n=np.linalg.norm(Riq)
# # if(Riq_n <= 1 and Riq_n >= 0):
#  if(0 <= Riq <= 1):
#   data3=np.sqrt(-2*np.log(Riq)/Riq)
#   data.append(v1[0][0]*data3[0][0])
#
# plt.hist(data, bins=50)
# plt.show()

###########################################
size=100
sampler = qmc.Sobol(d=2, scramble=False)
size_sobol = int(np.log(size)/np.log(2))+1
print("size = ",size," ;size Sobol = ",size_sobol)
sample = sampler.random_base2(m=size_sobol)

data1=[]
while (len(data1)<size):
 data = sampler.random()
 v1=2.0*data[0][0]-1.0
 v2=2.0*data[0][1]-1.0
# Riq=np.linalg.norm(v1)**2+np.linalg.norm(v2)**2
 Riq=v1**2+v2**2
# Riq_n=np.linalg.norm(Riq)
# if(Riq_n <= 1 and Riq_n >= 0):
 if(0 <= Riq <= 1):
     #print("Riq=",Riq)
     data3=np.sqrt(-2*np.log(Riq)/Riq)
     #print("data=",v1*data3)
     data1.append(v1*data3)

plt.hist(data1, bins=50)
print("Sobol graphics")
plt.show()
plt.scatter(data1,range(len(data1)))
print("Sobol graphics")
plt.show()
###########################################################
data1=[]
while (len(data1)<size):
 v1=2.0*np.random.rand()-1.0
 v2=2.0*np.random.rand()-1.0
 Riq=v1**2+v2**2
 if(0 <= Riq <= 1):
     #print("Riq=",Riq)
     data3=np.sqrt(-2*np.log(Riq)/Riq)
     #print("data=",v1*data3)
     data1.append(v1*data3)

plt.hist(data1, bins=50)
print("random normal graphics")
plt.show()
plt.scatter(data1,range(len(data1)))
print("random normal graphics")
plt.show()

data1=[]
while (len(data1)<size):
     data1.append(np.random.normal())

plt.hist(data1, bins=50)
print("numpy normal graphics")
plt.show()
plt.scatter(data1,range(len(data1)))
print("numpy normal graphics")
plt.show()
