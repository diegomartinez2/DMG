#import Spectral_plasmons_analysis
#import Eliashberg
import Eliashberg_new_clean as Eliashberg
import numpy as np
import os.path
import matplotlib.pyplot as plt
from scipy import integrate

def Excel_data(filename,flag_1DP = False):
    """
    Read the data for the plasmons in cuprate from the excel files.
    If the file is for the 1DP...
    """
    out = Eliashberg.read_1_excel_file(filename=filename) #filenames=('1DP','HPI','HPII');filenames=('1DP_c','HPI_c','HPII'_c)
    #qx,qy,Omega,Gamma,Ratio = Eliashberg.Excel_data_parser(out)
    qx=out[:,0]
    qy=out[:,1]
    Omega=out[:,2]
    Gamma=out[:,3]
    Ratio=out[:,4]
    if flag_1DP: #make true for "1DP" or "1DP_c"
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
    return qx,qy,Omega/1000,Gamma/1000,Ratio

def main(arg):
    calculate_all_hyperplasmons = False
    calculate_all_hyperplasmons_and_1DP = False
    plot_contour_flag = False

    if calculate_all_hyperplasmons:
        file_HP = "HPI"
        qx,qy,Omega,Gamma,Ratio = Excel_data(filename="{}_c2".format(file_HP))
        file_HP = "HPII"
        qx_temp,qy_temp,Omega_temp,Gamma_temp,Ratio_temp = Excel_data(filename="{}_c2".format(file_HP))
        qx = np.append(qx,qx_temp)
        qy = np.append(qy,qy_temp)
        Omega = np.append(Omega,Omega_temp)
        Gamma = np.append(Gamma,Gamma_temp)
        Ratio = np.append(Ratio,Ratio_temp)
        if calculate_all_hyperplasmons_and_1DP:
            file_HP = "1DP"
            qx_temp,qy_temp,Omega_temp,Gamma_temp,Ratio_temp = Excel_data(filename="{}_c3".format(file_HP),flag_1DP=False)
            qx = np.append(qx,qx_temp)
            qy = np.append(qy,qy_temp)
            Omega = np.append(Omega,Omega_temp)
            Gamma = np.append(Gamma,Gamma_temp)
            Ratio = np.append(Ratio,Ratio_temp)
    else:
        file_HP = "1DP"
        qx,qy,Omega,Gamma,Ratio = Excel_data(filename="{}_c3".format(file_HP),flag_1DP=False)
    superconductor = Eliashberg.Eliashberg2(qx,qy,Omega,Gamma,Ratio)
    superconductor.read_Ne()
    print("Omega range:",np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),'::',np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))))
    frequencies = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),10000) #test
    if plot_contour_flag:
        superconductor.plot_contour_d(data=[],mask_value=0, diagonal = True)
    lambda_1 = superconductor.Lambda_new(frequencies)
    #np.savetxt('Lambda.txt', (lambda_1))
    print('Lambda_1[{}]='.format(file_HP),lambda_1,'[Lambda calculated from Lambda_q for {}]?'.format(file_HP)) # Lambda calculated from Lambda_q
    print('Lambda_2[{}]='.format(file_HP),superconductor.lambda_2,'[Lambda calculated fron Eliashberg function for {}]?'.format(file_HP)) #Lambda calculated fron Eliashberg function

    print('Test Gaussian->')
    gauss = np.zeros(100)
    for w in range(100): #Make this from -100 to 100??
        gauss[w] = superconductor.gaussian(w, 50.0,10.0)
    plt.figure(figsize=(10,6))
    plt.plot(gauss)
    plt.tight_layout()
    plt.show()
    print(np.trapz(gauss))
    w = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000)
    #superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1]))
    #suma = integrate.simpson(superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1])), w)
    suma = integrate.simpson(superconductor.test_gaussian(w,superconductor.Omega,len(superconductor.Omega)), w)
    print(suma,'::',suma/len(superconductor.Omega))
    np.savetxt('salida.txt',np.c_[superconductor.qx,superconductor.qy,np.full(len(superconductor.Omega), superconductor.N_ef),superconductor.Omega,superconductor.Gamma,superconductor.lambda_q_lista],header='#---q_x---q_y---N_ef[States/Spin/eV/Unit Cell]---Omega[eV]---Gamma[eV]---Lambda_q')
#-----plot-start
    fig_lambda_q = plt.figure(figsize=(10,6))
    ax = fig_lambda_q.add_subplot(1, 1, 1)
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
    a2F_lista = superconductor.a2F_new(frequencies)

    fig_a2F = plt.figure(figsize=(10,6))
    ax = fig_a2F.add_subplot(1, 1, 1)
    ax.plot(a2F_lista,frequencies)
    ax.set_title('a2F vs. $\omega$')
    ax.set_ylabel('$\omega$ (eV)')
    ax.set_xlabel('a2F')

    plt.plot ()
    plt.tight_layout()
    plt.show()
    fig_a2F.savefig("Ajuste_d_{}".format("a2F"))
#---plot-end
    # # superconductor.lambda_2 += float(superconductor2.lambda_2)
    # np.savetxt("Lambda_{}".format(file_HP),[superconductor.lambda_2])
    # np.savetxt("a2F_{}".format(file_HP),np.vstack((frequencies, a2F_lista))) #.T)
    # lambda_HPI = np.loadtxt("Lambda_HPI_c")
    # lambda_HPII = np.loadtxt("Lambda_HPII_c")
    # frequencies_HPI,a2F_HPI = np.loadtxt("a2F_HPI")
    # frequencies_HPII,a2F_HPII = np.loadtxt("a2F_HPII")
    # # mask = a2F_HPII[:,  0] == a2F_HPI[:,  0]
    # # print("mask=",mask)
    # # np.add.at(a2F_HPII, (mask,  1), a2F_HPI[mask,  1])
    # fig_a2Ftot = plt.figure(figsize=(10,6))
    # ax = fig_a2Ftot.add_subplot(1, 1, 1)
    # ax.plot(a2F_HPI,frequencies_HPI)
    # ax.plot(a2F_HPII,frequencies_HPII)
    # ax.set_title('a2F total vs. $\omega$')
    # ax.set_ylabel('$\omega$ (meV)')
    # ax.set_xlabel('a2F HPI+HPII')
    # plt.plot ()
    # plt.tight_layout()
    # plt.show()
    print(superconductor.lambda_w_lista[-1])
    #T_c = superconductor.T_c(mu_par=0.1,lambda_t=superconductor.lambda_2) #(lambda_HPI+lambda_HPII)) #mu*=0.1 y mu*=0.15. Son los valores típicos.
    T_c = superconductor.T_c(mu_par=0.1,lambda_t=superconductor.lambda_w_lista[-1]) #(lambda_HPI+lambda_HPII)) #mu*=0.1 y mu*=0.15. Son los valores típicos.
    print("T_c=",T_c,"K")
    print("test T_C*K_b²=",T_c*8.617333262e-5*8.617333262e-5)
    np.savetxt("lambda_and_T_C.txt",(superconductor.lambda_2,T_c), header='Lambda, T_C (K)')
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
