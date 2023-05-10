import matplotlib.pyplot as plt
import scipy
from scipy.signal import find_peaks
import numpy as np

#with open("xaa") as file:
with open("xbv") as file:
#with open("xzaez") as file:
#with open("xzabu") as file:
#with open("xtg") as file:
     lines = [line.rsplit() for line in file]
Frequency = np.zeros(5001)
Spec = np.zeros(5001)
for i in range(5001):
     print (lines[i][2])
     Frequency[i]=lines[i][2]
     Spec[i]=lines[i][6]

#Frequency = np.array([ 5.0, 6.3, 8.0, 10.0, 12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0]) #third octave band spectrum, 19 Values
#Spec = np.array([ 40, 45, 51, 42, 44, 56, 42, 55, 57, 58, 45, 40, 38, 36, 32, 30, 28, 30, 29]) #noise level, 19 Values

peaks, _ = find_peaks(Spec, prominence=1)
m = np.zeros(Frequency.shape, dtype=bool)
m[peaks] = True

#-------------------------------
amp1 = 0.03
cen1 = 0.3
wid1 = 0.1

amp2 = 0.14
cen2 = 0.7
wid2 = 0.3

amp3 = 0.1
cen3 = 0.9
wid3 = 0.1
x_array = Frequency
y_array_3lorentz = Spec
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
#------------------------------


#plt.plot(Frequency[peaks], Spec[peaks], "xr", label="prominence")
plt.plot(Frequency, Spec, label="original")
plt.plot(Frequency, lorentz_peak_1+lorentz_peak_2+lorentz_peak_3, "ro", label="ajuste")
plt.plot(Frequency, lorentz_peak_1, label="1")
plt.plot(Frequency, lorentz_peak_2, label="2")
plt.plot(Frequency, lorentz_peak_3, label="3")
plt.legend()

plt.tight_layout()
plt.show()
