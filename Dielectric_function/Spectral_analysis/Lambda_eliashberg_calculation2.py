#import Spectral_plasmons_analysis
#import Eliashberg
import Eliashberg_new as Eliashberg
import numpy as np
import os.path
import matplotlib.pyplot as plt
from scipy import integrate

def main(arg):
    #pars_create=True
    #file1 = './pars_data.txt'
    #file2 = './frequencies_data.txt'
    #if pars_create:
    #if not (os.path.isfile(file1) and os.path.isfile(file2)):
        #print('------------------------')
        #print('Creating parameters file')
        #print('------------------------')
        #filelist=('1DP','HPI','HPII')
        #'xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal',
        #'xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay',
        #'xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg','xbh','xbi','xbj','xbk','xbl',
        #'xbm','xbn','xbo','xbp','xbq','xbr','xbs','xbt','xbu','xbv','xbw','xbx')
        #for namefile in filelist:
            #plasmon = Eliashberg.Plasmon_analysis(arg,namefile)
            #data_all=Eliashberg.read_1_excel_file(file_directory=".",filename=namefile)
            #data, frequencies, qx= plasmon.load_data()
            #print (np.shape(data),"=(51,5001)?")
            #print (np.shape(data_all))
            #plasmon.fitting_Lorentz(frequencies,data)
            #if (namefile == filelist[0]):
                #pars = plasmon.pars
                #print('pars_1=',pars)
            #else:
                #pars = np.concatenate((pars,plasmon.pars), axis=0)
                #print('pars_n=',pars)
        #print (np.shape(pars),"=(2550,3)?")
        #np.savetxt('pars_txt.dat', pars) #no negative values
        #np.savetxt(file1, pars) #no negative values
        #np.savetxt(file2, frequencies) #no negative values
    #else:
        #np.loadtxt('pars_txt.dat', pars) #no negative values
        #pars = [] #np.zeros((2550,3))
        #pars = np.loadtxt(file1) #no negative values,
        #frequencies = np.loadtxt(file2) #no negative values
    #print(pars[:,0])
    #print('Min(pars[0])=',np.amin(pars[:,0]))
    out = Eliashberg.read_1_excel_file(filename="1DP")
    qx,qy,Omega,Gamma,Ratio = Eliashberg.Excel_data_parser(out)
    superconductor = Eliashberg.Eliashberg2(qx,qy,Omega,Gamma,Ratio)
    superconductor.read_Ne()
    print("Omega range:",np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),'::',np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))))
    frequencies = np.linspace(np.min(superconductor.Omega-np.abs(np.max(superconductor.Gamma))),np.max(superconductor.Omega+np.abs(np.max(superconductor.Gamma))),20000) #test
    #lambda_1, lambda_2 = superconductor.Lambda(frequencies)
#    lambda_1 = superconductor.Lambda(frequencies)

    lambda_1 = superconductor.Lambda_new(frequencies)

    #np.savetxt('Lambda.txt', (lambda_1))
    print('Lambda_1=',lambda_1,'[]?') # Lambda calculated from Lambda_q
    #print('Lambda_1=',lambda_1*superconductor.Ne[np.where(superconductor.energy==0.0)[0][0]]) # Lambda calculated from Lambda_q#self.Ne[np.where(self.energy==0.0)[0][0]]
    print('Lambda_2=',superconductor.lambda_2,'[]?') #Lambda calculated fron Eliashberg function
    # Frequncies = frequencies+1
    # Frequencies = np.append(frequencies,Frequncies, axis=0)
    #print(len(Frequencies[1:]),len(superconductor.lambda_2))
    print('Test Gaussian->')
    gauss = np.zeros(100)
    for w in range(100): #Make this from -100 to 100??
        gauss[w] = superconductor.gaussian(w, 50.0,10.0)
    plt.figure(figsize=(10,6))
    plt.plot(gauss)
    plt.show()
    print(np.trapz(gauss))
    w = np.linspace(-100,100,20000)
    #superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1]))
    suma = integrate.simpson(superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1])), w)
    print(suma,'::',suma/len(pars[:,1]))
    #np.savetxt('salida.txt',(q_x,q_y,superconductor.N_ef,superconductor.pars[:,1],superconductor.pars[:,2],superconductor.lambda_q_lista))
    #print(len(superconductor.pars[:,1]),len(superconductor.pars[:,2]),len(superconductor.lambda_q_lista))
    a=np.arange(0,0.43458,0.00869)
    b=np.repeat(a,50)
    c=np.tile(a,50)
    print('len=',len(b),'=',len(c))
    np.savetxt('salida.txt',np.c_[c,b,np.full(len(superconductor.pars[:,1]), superconductor.N_ef),superconductor.pars[:,1],superconductor.pars[:,2],superconductor.lambda_q_lista],header='#---q_x---q_y---N_ef[States/Spin/eV/Unit Cell]---Omega[eV]---Gamma[eV]---Lambda_q')
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
