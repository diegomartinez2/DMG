# DATOS FINALES
grep -A 10 "V_eq" npt.log

# GRÁFICO CONVERGENCIA
grep Step npt.log | awk '{print $1,$5}' > p_vs_v.dat
python -c "
import matplotlib.pyplot as plt; import numpy as np
d=np.loadtxt('p_vs_v.dat'); plt.plot(d[:,1],d[:,0],'b-'); plt.xlabel('Vol (Å³)'); plt.ylabel('Press (bar)'); plt.axhline(0,c='r'); plt.title('V_eq chimesFF'); plt.show()
"
