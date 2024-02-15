#import Spectral_plasmons_analysis
#import Eliashberg
import Eliashberg_new2 as Eliashberg
import numpy as np
import os.path
import matplotlib.pyplot as plt
from scipy import integrate

def main(arg):

    out = Eliashberg.read_1_excel_file(filename="1DP_c") #filenames=('1DP','HPI','HPII');filenames=('1DP_c','HPI_c','HPII'_c)
    print(out) #test to see the format of out
    pass #remove later
    qx,qy,Omega,Gamma,Ratio = Eliashberg.Excel_data_parser(out)
    #qx,qy,Omega,Gamma,Ratio,Omega2,Gamma2,Ratio2,Omega3,Gamma3,Ratio3 = Eliashberg.Excel_data_parser(out)
    superconductor = Eliashberg.Eliashberg2(qx,qy,Omega,Gamma,Ratio)
    superconductor.read_Ne()
    print("Omega range:",np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),'::',np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))))
    #frequencies = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000) #test
    frequencies = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000) #test

    lambda_1 = superconductor.Lambda_new(frequencies)

    #np.savetxt('Lambda.txt', (lambda_1))
    print('Lambda_1=',lambda_1,'[Lambda calculated from Lambda_q]?') # Lambda calculated from Lambda_q
    #print('Lambda_1=',lambda_1*2*superconductor.Ne[np.where(superconductor.energy==0.0)[0][0]]) # Lambda calculated from Lambda_q#self.Ne[np.where(self.energy==0.0)[0][0]]
    print('Lambda_2=',superconductor.lambda_2,'[Lambda calculated fron Eliashberg function]?') #Lambda calculated fron Eliashberg function

    print('Test Gaussian->')
    gauss = np.zeros(100)
    for w in range(100): #Make this from -100 to 100??
        gauss[w] = superconductor.gaussian(w, 50.0,10.0)
    plt.figure(figsize=(10,6))
    plt.plot(gauss)
    plt.show()
    print(np.trapz(gauss))
    w = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000)
    #superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1]))
    #suma = integrate.simpson(superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1])), w)
    suma = integrate.simpson(superconductor.test_gaussian(w,superconductor.Omega,len(superconductor.Omega)), w)
    print(suma,'::',suma/len(superconductor.Omega))
    #np.savetxt('salida.txt',(q_x,q_y,superconductor.N_ef,superconductor.pars[:,1],superconductor.pars[:,2],superconductor.lambda_q_lista))
    #print(len(superconductor.pars[:,1]),len(superconductor.pars[:,2]),len(superconductor.lambda_q_lista))
    #a=np.arange(0,0.43458,0.00869)
    #b=np.repeat(a,50)
    #c=np.tile(a,50)
    #print('len=',len(b),'=',len(c))
    np.savetxt('salida.txt',np.c_[superconductor.qx,superconductor.qy,np.full(len(superconductor.Omega), superconductor.N_ef),superconductor.Omega,superconductor.Gamma,superconductor.lambda_q_lista],header='#---q_x---q_y---N_ef[States/Spin/eV/Unit Cell]---Omega[eV]---Gamma[eV]---Lambda_q')
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
