import matplotlib.pyplot as plt
import scipy
from scipy.signal import find_peaks
import numpy as np

#with open("xaa") as file:
#with open("xbv") as file:
#with open("xzaez") as file:
#with open("xzabu") as file:
#with open("xtg") as file:
#with open("xjz") as file:

namefile="xju"
#namefile="xjv"
#namefile="xjw"
#namefile="xjx"
with open(namefile) as file:
     lines = [line.rsplit() for line in file]
Frequency = np.zeros(5001)
Spec = np.zeros(5001)
for i in range(5001):
     print (lines[i][2])
     Frequency[i]=lines[i][2]
     Spec[i]=lines[i][6]


peaks, _ = find_peaks(Spec, prominence=1)
m = np.zeros(Frequency.shape, dtype=bool)
m[peaks] = True

#-------------------------------
amp1 = 0.03 #03
cen1 = 0.06 #11
wid1 = 0.03 #1

amp2 = 0.014
cen2 = 0.07
wid2 = 0.03

amp3 = 0.01
cen3 = 0.1
wid3 = 0.01

amp4 = 0.002
cen4 = 0.07
wid4 = 0.06

x_array = Frequency
y_array_3lorentz = Spec
def _base(x):
    if (x < 0):
        x_out = 0
    elif (x <= 0.0271):
        x_out = -1.9985*x*x+0.0651067*x+0.000119668
    elif (x >= 0.0271):
        x_out = -1.99329*x*x+0.264442*x-0.00527247
    else:
        x_out = 0
    return x_out

def _1Lorentzian(x, amp, cen, wid):
    return amp*wid**2/((x-cen)**2+wid**2)+_base(x)

def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
            return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                    (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                        (amp3*wid3**2/((x-cen3)**2+wid3**2))

def _4Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3, amp4,cen4,wid4):
            return (amp1*wid1**2/((x-cen1)**2+wid1**2)) + (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                        (amp3*wid3**2/((x-cen3)**2+wid3**2)) + (amp4*wid4**2/((x-cen4)**2+wid4**2))

#--------------------------------------------
plt.plot(Frequency, Spec, label="original")
base = np.zeros(Frequency.shape)
for i in range(len(Frequency)):
   base[i]=_base(Frequency[i])
plt.plot(Frequency,base, label="base")
plt.xlim([0, 0.13])
plt.ylim([0,0.015])
plt.tight_layout()
plt.savefig("Original_{}".format(namefile))
plt.show()
#--------------------------------------------
popt_1lorentz, pcov_1lorentz = scipy.optimize.curve_fit(_1Lorentzian, x_array, y_array_3lorentz, p0=[amp1, cen1, wid1])
#popt_3lorentz, pcov_3lorentz = scipy.optimize.curve_fit(_3Lorentzian, x_array, y_array_3lorentz, p0=[amp1, cen1, wid1, amp2, cen2, wid2, amp3, cen3, wid3])
#popt_4lorentz, pcov_4lorentz = scipy.optimize.curve_fit(_4Lorentzian, x_array, y_array_3lorentz, p0=[amp1, cen1, wid1, amp2, cen2, wid2, amp3, cen3, wid3, amp4, cen4, wid4])

perr_1lorentz = np.sqrt(np.diag(pcov_1lorentz))
#perr_3lorentz = np.sqrt(np.diag(pcov_3lorentz))
#perr_4lorentz = np.sqrt(np.diag(pcov_4lorentz))

pars_1 = popt_1lorentz[0:3]
#pars_1 = popt_3lorentz[0:3]
#pars_2 = popt_3lorentz[3:6]
#pars_3 = popt_3lorentz[6:9]
#pars_1 = popt_4lorentz[0:3]
#pars_2 = popt_4lorentz[3:6]
#pars_3 = popt_4lorentz[6:9]
#pars_4 = popt_4lorentz[9:12]
lorentz_peak_1 = _1Lorentzian(x_array, *pars_1)
#lorentz_peak_2 = _1Lorentzian(x_array, *pars_2)
#lorentz_peak_3 = _1Lorentzian(x_array, *pars_3)
#lorentz_peak_4 = _1Lorentzian(x_array, *pars_4)

## this cell prints the fitting parameters with their errors
print ("-------------Peak 1-------------")
print ("amplitude = {} (+/-) {}".format(pars_1[0], perr_1lorentz[0]))
print ("center = {} (+/-) {}".format(pars_1[1], perr_1lorentz[1]))
print ("width = {} (+/-) {}".format(pars_1[2], perr_1lorentz[2]))
print ("area = {}".format(np.trapz(lorentz_peak_1)))
print ("--------------------------------")
print ("-------------Peak 2-------------")
#print ("amplitude = {} (+/-) {}".format(pars_2[0], perr_3lorentz[3]))
#print ("center = {} (+/-) {}".format(pars_2[1], perr_3lorentz[4]))
#print ("width = {} (+/-) {}".format(pars_2[2], perr_3lorentz[5]))
#print ("area = {}".format(np.trapz(lorentz_peak_2)))
print ("--------------------------------")
print ("-------------Peak 3-------------")
#print ("amplitude = {} (+/-) {}".format(pars_3[0], perr_3lorentz[6]))
#print ("center = {} (+/-) {}".format(pars_3[1], perr_3lorentz[7]))
#print ("width = {} (+/-) {}".format(pars_3[2], perr_3lorentz[8]))
#print ("area = {}".format(np.trapz(lorentz_peak_3)))
print ("--------------------------------")

#print ("-------------Peak 1-------------")
#print ("amplitude = %0.2f (+/-) %0.2f" % (pars_1[0], perr_4lorentz[0]))
#print ("center = %0.2f (+/-) %0.2f" % (pars_1[1], perr_4lorentz[1]))
#print ("width = %0.2f (+/-) %0.2f" % (pars_1[2], perr_4lorentz[2]))
#print ("area = %0.2f" % np.trapz(lorentz_peak_1))
#print ("--------------------------------")
#print ("-------------Peak 2-------------")
#print ("amplitude = %0.2f (+/-) %0.2f" % (pars_2[0], perr_4lorentz[3]))
#print ("center = %0.2f (+/-) %0.2f" % (pars_2[1], perr_4lorentz[4]))
#print ("width = %0.2f (+/-) %0.2f" % (pars_2[2], perr_4lorentz[5]))
#print ("area = %0.2f" % np.trapz(lorentz_peak_2))
#print ("--------------------------------")
#print ("-------------Peak 3-------------")
#print ("amplitude = %0.2f (+/-) %0.2f" % (pars_3[0], perr_4lorentz[6]))
#print ("center = %0.2f (+/-) %0.2f" % (pars_3[1], perr_4lorentz[7]))
#print ("width = %0.2f (+/-) %0.2f" % (pars_3[2], perr_4lorentz[8]))
#print ("area = %0.2f" % np.trapz(lorentz_peak_3))
#print ("--------------------------------")
#print ("-------------Peak 4-------------")
#print ("amplitude = %0.2f (+/-) %0.2f" % (pars_4[0], perr_4lorentz[9]))
#print ("center = %0.2f (+/-) %0.2f" % (pars_4[1], perr_4lorentz[10]))
#print ("width = %0.2f (+/-) %0.2f" % (pars_4[2], perr_4lorentz[11]))
#print ("area = %0.2f" % np.trapz(lorentz_peak_4))
#print ("--------------------------------")

residual_1lorentz = y_array_3lorentz - (_1Lorentzian(x_array, *popt_1lorentz))
#residual_3lorentz = y_array_3lorentz - (_3Lorentzian(x_array, *popt_3lorentz))
#residual_4lorentz = y_array_3lorentz - (_4Lorentzian(x_array, *popt_4lorentz))
#------------------------------


plt.plot(Frequency, Spec, label="original")
#plt.plot(Frequency, lorentz_peak_1+lorentz_peak_2+lorentz_peak_3, "r",ls=':', label="ajuste")
#plt.plot(Frequency, lorentz_peak_1+lorentz_peak_2+lorentz_peak_3+lorentz_peak_4, "ro", label="ajuste")
plt.plot(Frequency, lorentz_peak_1, label="1")
#plt.plot(Frequency, lorentz_peak_2, label="2")
#plt.plot(Frequency, lorentz_peak_3, label="3")
#plt.plot(Frequency, lorentz_peak_4, label="4")
plt.legend()
plt.xlim([0, 0.13])
plt.ylim([0,0.015])
plt.tight_layout()
plt.savefig("Ajuste_{}".format(namefile))
plt.show()
