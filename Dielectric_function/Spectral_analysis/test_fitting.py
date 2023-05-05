import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
#%matplotlib inline
x_array = np.linspace(1,300,250)

amp1 = 50
cen1 = 100
wid1 = 5

amp2 = 100
cen2 = 150
wid2 = 10

amp3 = 50
cen3 = 200
wid3 = 5

y_array_3lorentz = (amp1*wid1**2/((x_array-cen1)**2+wid1**2)) + \
                    (amp2*wid2**2/((x_array-cen2)**2+wid2**2)) +\
                     (amp3*wid3**2/((x_array-cen3)**2+wid3**2))

# creating some noise to add the the y-axis data
y_noise_3lorentz = (((np.random.ranf(250))))*5
y_array_3lorentz += y_noise_3lorentz
fig = plt.figure(figsize=(4,3))
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0])

ax1.plot(x_array, y_array_3lorentz, "ro")

#ax1.set_xlim(-5,105)
#ax1.set_ylim(-0.5,5)

ax1.set_xlabel("x_array",family="serif",  fontsize=12)
ax1.set_ylabel("y_array",family="serif",  fontsize=12)

ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

ax1.tick_params(axis='both',which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
ax1.tick_params(axis='both',which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)

fig.tight_layout()
fig.savefig("raw_3Lorentz.png", format="png",dpi=1000)

def _1Lorentzian(x, amp, cen, wid):
    return amp*wid**2/((x-cen)**2+wid**2)

def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
            (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                (amp3*wid3**2/((x-cen3)**2+wid3**2))

popt_3lorentz, pcov_3lorentz = scipy.optimize.curve_fit(_3Lorentzian, x_array, y_array_3lorentz, p0=[amp1, cen1, wid1, \
                                                                                    amp2, cen2, wid2, amp3, cen3, wid3])

perr_3lorentz = np.sqrt(np.diag(pcov_3lorentz))

pars_1 = popt_3lorentz[0:3]
pars_2 = popt_3lorentz[3:6]
pars_3 = popt_3lorentz[6:9]
lorentz_peak_1 = _1Lorentzian(x_array, *pars_1)
lorentz_peak_2 = _1Lorentzian(x_array, *pars_2)
lorentz_peak_3 = _1Lorentzian(x_array, *pars_3)

# this cell prints the fitting parameters with their errors
print ("-------------Peak 1-------------")
print ("amplitude = %0.2f (+/-) %0.2f" % (pars_1[0], perr_3lorentz[0]))
print ("center = %0.2f (+/-) %0.2f" % (pars_1[1], perr_3lorentz[1]))
print ("width = %0.2f (+/-) %0.2f" % (pars_1[2], perr_3lorentz[2]))
print ("area = %0.2f" % np.trapz(lorentz_peak_1))
print ("--------------------------------")
print ("-------------Peak 2-------------")
print ("amplitude = %0.2f (+/-) %0.2f" % (pars_2[0], perr_3lorentz[3]))
print ("center = %0.2f (+/-) %0.2f" % (pars_2[1], perr_3lorentz[4]))
print ("width = %0.2f (+/-) %0.2f" % (pars_2[2], perr_3lorentz[5]))
print ("area = %0.2f" % np.trapz(lorentz_peak_2))
print ("--------------------------------")
print ("-------------Peak 3-------------")
print ("amplitude = %0.2f (+/-) %0.2f" % (pars_3[0], perr_3lorentz[6]))
print ("center = %0.2f (+/-) %0.2f" % (pars_3[1], perr_3lorentz[7]))
print ("width = %0.2f (+/-) %0.2f" % (pars_3[2], perr_3lorentz[8]))
print ("area = %0.2f" % np.trapz(lorentz_peak_3))
print ("--------------------------------")

residual_3lorentz = y_array_3lorentz - (_3Lorentzian(x_array, *popt_3lorentz))

ig = plt.figure(figsize=(4,4))
gs = gridspec.GridSpec(2,1, height_ratios=[1,0.25])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
gs.update(hspace=0)

ax1.plot(x_array, y_array_3lorentz, "ro")
ax1.plot(x_array, _3Lorentzian(x_array, *popt_3lorentz), 'k--')#,\
         #label="y= %0.2f$e^{%0.2fx}$ + %0.2f" % (popt_exponential[0], popt_exponential[1], popt_exponential[2]))

# peak 1
ax1.plot(x_array, lorentz_peak_1, "g")
ax1.fill_between(x_array, lorentz_peak_1.min(), lorentz_peak_1, facecolor="green", alpha=0.5)

# peak 2
ax1.plot(x_array, lorentz_peak_2, "y")
ax1.fill_between(x_array, lorentz_peak_2.min(), lorentz_peak_2, facecolor="yellow", alpha=0.5)

# peak 3
ax1.plot(x_array, lorentz_peak_3, "c")
ax1.fill_between(x_array, lorentz_peak_3.min(), lorentz_peak_3, facecolor="cyan", alpha=0.5)

# residual
ax2.plot(x_array, residual_3lorentz, "bo")

#ax1.set_xlim(-5,105)
#ax1.set_ylim(-0.5,8)

#ax2.set_xlim(-5,105)
#ax2.set_ylim(-0.5,0.75)

ax2.set_xlabel("x_array",family="serif",  fontsize=12)
ax1.set_ylabel("y_array",family="serif",  fontsize=12)
ax2.set_ylabel("Res.",family="serif",  fontsize=12)

ax1.legend(loc="best")

ax1.xaxis.set_major_locator(ticker.MultipleLocator(20))
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

ax1.xaxis.set_major_formatter(plt.NullFormatter())

ax1.tick_params(axis='x',which='major', direction="out", top="on", right="on", bottom="off", length=8, labelsize=8)
ax1.tick_params(axis='x',which='minor', direction="out", top="on", right="on", bottom="off", length=5, labelsize=8)
ax1.tick_params(axis='y',which='major', direction="out", top="on", right="on", bottom="off", length=8, labelsize=8)
ax1.tick_params(axis='y',which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)

ax2.tick_params(axis='x',which='major', direction="out", top="off", right="on", bottom="on", length=8, labelsize=8)
ax2.tick_params(axis='x',which='minor', direction="out", top="off", right="on", bottom="on", length=5, labelsize=8)
ax2.tick_params(axis='y',which='major', direction="out", top="off", right="on", bottom="on", length=8, labelsize=8)
ax2.tick_params(axis='y',which='minor', direction="out", top="off", right="on", bottom="on", length=5, labelsize=8)

fig.tight_layout()
fig.show()
fig.savefig("fit3Lorentzian_peaks_resid.png", format="png",dpi=1000)
