#import Spectral_plasmons_analysis
#import Eliashberg
import Eliashberg_new2 as Eliashberg
import numpy as np
import os.path
import matplotlib.pyplot as plt
from scipy import integrate

def Excel_data(filename):
    out = Eliashberg.read_1_excel_file(filename=filename) #filenames=('1DP','HPI','HPII');filenames=('1DP_c','HPI_c','HPII'_c)
    #qx,qy,Omega,Gamma,Ratio = Eliashberg.Excel_data_parser(out)
    qx=out[:,0]
    qy=out[:,1]
    Omega=out[:,2]
    Gamma=out[:,3]
    Ratio=out[:,4]
    if False: #make true for "1DP" or "1DP_c"
        qx=np.append(qx,out[:,0])
        qy=np.append(qy,out[:,1])
        Omega=np.append(Omega,out[:,5])
        Gamma=np.append(Gamma,out[:,6])
        Ratio=np.append(Ratio,out[:,7])
        qx=np.append(qx,out[:,0])
        qy=np.append(qy,out[:,1])
        Omaga=np.append(Omega,out[:,8])
        Gamma=np.append(Gamma,out[:,9])
        Ratio=np.append(Ratio,out[:,10])
    return qx,qy,Omega,Gamma,Ratio

def main(arg):

    # out = Eliashberg.read_1_excel_file(filename="HPI_c") #filenames=('1DP','HPI','HPII');filenames=('1DP_c','HPI_c','HPII'_c)
    # #qx,qy,Omega,Gamma,Ratio = Eliashberg.Excel_data_parser(out)
    # qx=out[:,0]
    # qy=out[:,1]
    # Omega=out[:,2]
    # Gamma=out[:,3]
    # Ratio=out[:,4]
    # if False: #make true for "1DP" or "1DP_c"
    #     qx=np.append(qx,out[:,0])
    #     qy=np.append(qy,out[:,1])
    #     Omega=np.append(Omega,out[:,5])
    #     Gamma=np.append(Gamma,out[:,6])
    #     Ratio=np.append(Ratio,out[:,7])
    #     qx=np.append(qx,out[:,0])
    #     qy=np.append(qy,out[:,1])
    #     Omaga=np.append(Omega,out[:,8])
    #     Gamma=np.append(Gamma,out[:,9])
    #     Ratio=np.append(Ratio,out[:,10])
    file_HP = "HPII"
    qx,qy,Omega,Gamma,Ratio = Excel_data(filename="{}_c".format(file_HP))
    superconductor = Eliashberg.Eliashberg2(qx,qy,Omega,Gamma,Ratio)
    superconductor.read_Ne()
    # qx,qy,Omega,Gamma,Ratio = Excel_data(filename="HPII_c")
    # superconductor2 = Eliashberg.Eliashberg2(qx,qy,Omega,Gamma,Ratio)
    # superconductor2.read_Ne()
    print("Omega range:",np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),'::',np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))))
    #frequencies = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000) #test
    frequencies = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000) #test
    #superconductor2.plot_contour_d(data=[],mask_value=0, diagonal = True)
    lambda_1 = superconductor.Lambda_new(frequencies)
    # lambda_1_2 = superconductor2.Lambda_new(frequencies)
    #np.savetxt('Lambda.txt', (lambda_1))
    print('Lambda_1[{}]='.format(file_HP),lambda_1,'[Lambda calculated from Lambda_q for {}]?'.format(file_HP)) # Lambda calculated from Lambda_q
    #print('Lambda_1=',lambda_1*2*superconductor.Ne[np.where(superconductor.energy==0.0)[0][0]]) # Lambda calculated from Lambda_q#self.Ne[np.where(self.energy==0.0)[0][0]]
    print('Lambda_2[{}]='.format(file_HP),superconductor.lambda_2,'[Lambda calculated fron Eliashberg function for {}]?'.format(file_HP)) #Lambda calculated fron Eliashberg function

    # print('Lambda_1[HPII]=',lambda_1_2,'[Lambda calculated from Lambda_q for HPII]?')
    # print('Lambda_2[HPII]=',superconductor2.lambda_2,'[Lambda calculated fron Eliashberg function for HPII]?')
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
#-----plot-start
    fig_lambda_q = plt.figure(figsize=(10,6))
    ax = fig_lambda_q.add_subplot(1, 1, 1)
    #ax.plot(superconductor.lambda_w_lista,frequencies[1:])
    #ax.plot(superconductor.lambda_w_lista,superconductor.w_0[1:])
    ax.plot(superconductor.lambda_w_lista,superconductor.w_0)
    ax.set_title('$\lambda$ vs. $\omega$')
    ax.set_xlabel('$\lambda(\omega)$')
    ax.set_ylabel('$\omega$')
    #ax.set_xticks([0,len(superconductor.w)])
    #ax.set_yticks([0,5001])
    #ax.set_xticklabels(["0","$\pi$"])
    #ax.set_xticklabels([superconductor.w[0],superconductor.w[-1]])
    #ax.set_yticklabels(["0","1"])
    plt.tight_layout()
    plt.show()
    fig_lambda_q.savefig("Ajuste_d_{}".format("lambda_w"))
    a2F_lista = []
    #print("frequencies=",frequencies)
    #for w in frequencies:
    #    a2F_lista.append(superconductor.a2F_new(w))
    a2F_lista = superconductor.a2F_new(frequencies)

    #print ("a2F(w)=",a2F_lista)
    fig_a2F = plt.figure(figsize=(10,6))
    ax = fig_a2F.add_subplot(1, 1, 1)
    ax.plot(a2F_lista,frequencies)
    ax.set_title('a2F vs. $\omega$')
    ax.set_ylabel('$\omega$ (meV)')
    ax.set_xlabel('a2F')

    plt.plot ()
    plt.tight_layout()
    plt.show()
    fig_a2F.savefig("Ajuste_d_{}".format("a2F"))
#---plot-end
    # superconductor.lambda_2 += float(superconductor2.lambda_2)
    np.savetxt("Lambda_{}".format(file_HP),[superconductor.lambda_2])
    np.savetxt("a2F_{}".format(file_HP),np.vstack((frequencies, a2F_lista)).T)
    lambda_HPI = np.loadtxt("Lambda_HPI_c")
    lambda_HPII = np.loadtxt("Lambda_HPII_c")
    frequencies_HPI,a2F_HPI = np.loadtxt("a2F_HPI")
    frequencies_HPII,a2F_HPII = np.loadtxt("a2F_HPII")
    # mask = a2F_HPII[:,  0] == a2F_HPI[:,  0]
    # print("mask=",mask)
    # np.add.at(a2F_HPII, (mask,  1), a2F_HPI[mask,  1])
    fig_a2Ftot = plt.figure(figsize=(10,6))
    ax = fig_a2Ftot.add_subplot(1, 1, 1)
    ax.plot(a2F_HPI,frequencies_HPI)
    ax.plot(a2F_HPII,frequencies_HPII)
    ax.set_title('a2F total vs. $\omega$')
    ax.set_ylabel('$\omega$ (meV)')
    ax.set_xlabel('a2F HPI+HPII')
    plt.plot ()
    plt.tight_layout()
    plt.show()

    T_c = superconductor.T_c(mu_par=0.1,lambda_t=(lambda_HPI+lambda_HPII)) #mu*=0.1 y mu*=0.15. Son los valores típicos.
    print("T_c=",T_c,"eV:: T_c=",T_c*11604,"K")
    np.savetxt("lambda_and_T_C.txt",(superconductor.lambda_2,T_c,T_c*11604), header='Lambda, T_C (eV), T_c (K)')
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
