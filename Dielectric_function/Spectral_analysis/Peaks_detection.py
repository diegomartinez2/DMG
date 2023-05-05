import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np

Frequency = np.array([ 5.0, 6.3, 8.0, 10.0, 12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0]) #third octave band spectrum, 19 Values
Spec = np.array([ 40, 45, 51, 42, 44, 56, 42, 55, 57, 58, 45, 40, 38, 36, 32, 30, 28, 30, 29]) #noise level, 19 Values

peaks, _ = find_peaks(Spec, prominence=1)
m = np.zeros(Frequency.shape, dtype=bool)
m[peaks] = True

plt.plot(Frequency[peaks], Spec[peaks], "xr", label="prominence")
plt.plot(Frequency, Spec, label="original")
plt.plot(Frequency[~m], Spec[~m], label="'smoothed'")
plt.legend()

plt.tight_layout()
plt.show()
