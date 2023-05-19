import numpy as np
import matplotlib.pyplot as plt
import sys, os
import scipy, scipy.optimize


# Load all the dynamical matrices and compute volume
DIRECTORY = "thermal_expansion"
FILE = os.path.join(DIRECTORY, "thermal_expansion.dat")

# Load the data from the final data file
temperatures, volumes = np.loadtxt(FILE, unpack = True)


# Prepare the figure and plot the V(T) from the sscha data
plt.figure(dpi = 150)
plt.scatter(temperatures, volumes, label = "SSCHA data")

# Fit the data to estimate the volumetric thermal expansion coefficient
def parabola(x, a, b, c):
    return a + b*x + c*x**2
def diff_parab(x, a, b, c):
    return b + 2*c*x

popt, pcov = scipy.optimize.curve_fit(parabola, temperatures, volumes,
                                      p0 = [0,0,0])

# Evaluate the volume thermal expansion
vol_thermal_expansion = diff_parab(300, *popt) / parabola(300, *popt)
plt.text(0.6, 0.2, r"$\alpha_v = "+"{:.1f}".format(vol_thermal_expansion*1e6)+r"\times 10^6 $ K$^{-1}$",
         transform = plt.gca().transAxes)


# Plot the fit
_t_ = np.linspace(np.min(temperatures), np.max(temperatures), 1000)
plt.plot(_t_, parabola(_t_, *popt), label = "Fit")

# Adjust the plot adding labels, legend, and saving in eps
plt.xlabel("Temperature [K]")
plt.ylabel(r"Volume [$\AA^3$]")
plt.legend()
plt.tight_layout()
plt.savefig("thermal_expansion.png")
plt.show()
