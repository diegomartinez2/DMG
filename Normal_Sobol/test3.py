import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
  
# Generate some data for this 
# demonstration.
data=[]
while (len(data)<250):
 data1 = np.random.random()
 data2 = np.random.random()
 v1=2.0*data1-1.0
 v2=2.0*data2-1.0
 Riq=v1**2+v2**2
 if(Riq <= 1 and Riq >= 0):
  data3=np.sqrt(-2*np.log(Riq)/Riq)
  data.append(v1*data3) 
# Fit a normal distribution to
# the data:
# mean and standard deviation
mu, std = norm.fit(data) 
  
# Plot the histogram.
plt.hist(data, bins=25, density=True, alpha=0.6, color='b')
  
# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
  
plt.plot(x, p, 'k', linewidth=2)
title = "Fit Values: {:.2f} and {:.2f}".format(mu, std)
plt.title(title)
  
plt.show()
