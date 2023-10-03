import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from ase.units import create_units
#from scipy.integrate import quad
from scipy import integrate

import scipy as scipy
from scipy.signal import find_peaks

class Plasmon_analysis(object):
    """Code for the analysis of Silkin's plasmon data."""

    def __init__(self, arg, namefile):
        super(Plasmon_analysis, self).__init__()
        self.arg = arg
        self.namefile = namefile
        self.pars = np.zeros((51,3))
        self.perr_lorentz = np.zeros((51,3))

    def load_data(self):
        with open(self.namefile) as file:
             lines = [line.rsplit() for line in file]
        Spec = np.zeros(255051)
        for i in range(255051):
             Spec[i]=lines[i][6]
        #data = np.resize(Spec,(51,int(len(Spec)/51)))
        data = np.resize(Spec,(51,5001))
        Frequency = np.zeros(5001)
        for i in range(5001):
             Frequency[i]=lines[i][2]
        print ("q_x=",lines[0][0])
        return data, Frequency, lines[0][0]

    def load_big_file(self, index, filename="A7_EPS.dat", diagonal=False):
        q_x, q_y, w, epsilon = np.loadtxt(filename,usecols=(0,1,2,6), unpack=True)
        print (len(np. unique(q_x))) #np. unique(my_array, return_counts=True)
        print (len(np. unique(q_y)))
        print (len(np. unique(w)))
        #self.data = np.resize(epsilon,(51,51,5001))
        self.data = np.resize(epsilon,(len(np. unique(q_x)),
                                        len(np. unique(q_y)),
                                        len(np. unique(w))))
        #plot.contour(data[i])
        """
        direccion q_y -> data[i]
        direccion q_x=q_y=i -> data = data[i][i] for i in {0..51}
        direccion
        """
        data = np.zeros((len(np. unique(q_x)),len(np. unique(w))))
        print ("+++++")
        if (diagonal == True):
            for i in range(len(np. unique(q_x))):
                print ("[" + "=" * i +  " " * ((len(np. unique(q_x)) - i)) + "]" +  str(i) + "%",end='\r', flush=True)
                for j in range(len(np. unique(w))):
                    data[i,j] = self.data[i,i,j]

        print ("*************")
        print (data.shape)
        # print (np.diagonal(self.data,axis1=0, axis2=1).T)
        # print (np.diagonal(self.data,axis1=0, axis2=1).shape)
        print ("--------------")

        print (np.shape(self.data[int(index)]),"=(51,5001)?")
        if (diagonal == True):
            return data, w[:5001]
        else:
            return self.data[int(index)], w[:5001]

    def plot_contour(self,data):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1.imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

        cbar = fig.colorbar(cax)
        plt.show()
        pass

    def plot_contour2(self,data,data2):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(2,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1[0].imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.175,
        	vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        cax2 = ax1[1].imshow(data2.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        #ax1[1].set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        print (self.pars[:,1])
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = self.pars[:,1]*5000,c='k',marker='x',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]+self.pars[:,2])*5000,c='k',marker='1',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]-self.pars[:,2])*5000,c='k',marker='2',s=10)

        cbar = fig.colorbar(cax)
        cbar2 = fig.colorbar(cax2)
        if (len( sys.argv) == 3):
            plt.savefig("Ajuste_{}".format("original_vs_fitted"))
        else:
            plt.savefig("Ajuste_{}_{}".format(self.namefile,"original_vs_fitted"))
        plt.show()
        pass

    def locate_1Lorenztian(self, Frequency, Spec, index):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        print (peaks,Frequency[peaks])
        print (peaks,_["widths"])
        amp = 0.009
        cen = 0.009
        wid = 0.001
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
        pars = popt_lorentz[:]
        lorentz_peak = _1Lorentzian(Frequency, *pars)
        print ("-------------Peak ",index,"-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars[0], perr_lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars[1], perr_lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars[2], perr_lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak))
        print ("--------------------------------")

        # plt.plot(Frequency, Spec, label="original")
        # plt.plot(Frequency, lorentz_peak, "r",ls=':', label="ajuste")
        # #plt.plot(Frequency[peaks], Spec[peaks], "xr", label="peaks")
        # #plt.plot(Frequency[peaks], Spec[peaks], "1k", label="center")
        # plt.legend()
        # plt.xlim([0, Frequency.max()])
        # plt.ylim([0, Spec.max()])
        # plt.tight_layout()
        # plt.savefig("Ajuste_{}_{}".format(self.namefile,index))
        # plt.show()
        self.pars[index] = pars
        self.perr_lorentz[index] = perr_lorentz
        #np.savetxt(f, (pars[0], perr_lorentz[0],pars[1], perr_lorentz[1],pars[2], perr_lorentz[2]) , fmt='%1.3f', newline=", ")
        #f.write("\n")
        return lorentz_peak

    def fitting_Lorentz(self,frequencies,data):
        self.Fitted_data = np.zeros((51,5001))
        for i in range(51):
            self.Fitted_data[i] = self.locate_1Lorenztian(frequencies,data[i],i)
        pass

class Eliashberg(object):
    """docstring for Eliashberg.
    This object calculates the Lambda from two methods:
    From the Eliashberg function amd fron the lambda at the q points.
    """

    def __init__(self, pars):
        units = create_units('2014')   #new way of defining units
        #super(Eliashberg, self).__init__()
        self.pars = pars
        #self.from_cm1_to_Hartree = 4.55633e-6 # from cm-1 to Hartree
        #self.from_cm1_to_Hartree = 1/219474.6 # from cm-1 to Hartree
        self.from_cm1_to_Hartree = units.invcm/units.Hartree #0.0000045563352812122295
        self.from_GHz_to_Hartree = self.from_cm1_to_Hartree /29.9792458 # from GHz to Hartree
        #self.from_GHz_to_Hartree = 1/6579683.879634054
        #self.from_GHz_to_Hartree = / 6.5796839 × 10^9
        #self.from_GHz_to_Hartree = 1.5198298554969748e-7
        self.from_Ry_to_Hartree = units.Ry/units.Hartree #0.5 # from Ry to Hartree
        """
        cm^-1 to Hartree: Hartree = cm^-1 / 219474.6
        Hartree to cm^-1: cm^-1 = Hartree * 219474.6

        GHz to Hartree: Hartree = GHz / 6.5796839 × 10^9
        Hartree to GHz: GHz = Hartree * 6.5796839 × 10^9

        eV to Hartree: Hartree = eV*0,0367493
        """
        self.from_cm1_to_eV = units.invcm/units.eV #0.00012398425731484318
        self.from_GHz_to_eV = 0.000004135669661004052

    def read_Ne(self,filename="out_DOS.dat"):
        self.energy, self.Ne = np.loadtxt(filename,usecols=(0,1), unpack=True)
#        self.Ne_0[np.where(energy==0.0)]
        #numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0])+1],dx=(self.energy[9]-self.energy[0])/10)
        numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0]+1)],dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy)) # maybe this is better??
        #numero_de_elementos2 = integrate.simpson(self.Ne[:(np.where(self.energy==0.0)[0][0]+1)])
        total_states = np.trapz(self.Ne,dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
        print("Número de elementos:",numero_de_elementos,"|total=",total_states)
        #test the selection
        #print (self.Ne[:(np.where(self.energy==0.0)[0][0]+1)])
        #print (self.Ne[:(np.where(self.energy==0.0)[0][0])+1])

        #------^--------
        print ("Density of states at Fermi level per cell=",self.Ne[np.where(self.energy==0.0)])
        print ("Number of elements (for the factor)=",numero_de_elementos)
        pass

    def Lambda_q(self,gamma_q,omega_q,Nef):
        """
        Calculates the Lambda(q) functions
        ---input---
        gamma_q: widths of the lorentzian fitting of the plasmon
        omega_q: Frequencies of the plasmon, fitted to lorentzian
        Nef: Density of states at Fermi level per cell
        ---output---
        Lamb_q: Lambda(q)
        """
        Lamb_q=(1/(np.pi*Nef)) * (gamma_q/omega_q**2) #fix from omega to omega²
        return Lamb_q

    def a2F(self,x):
        """
        Calculates the Eliashberg spectral functions
        ---input---
        x: coordiate to calculate the Eliashberg function.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """
        #print ("shape = ",np.shape(self.pars))
        center = self.pars[:,1]
        #center = np.absolute(center) #not necessary and better if not
        width = self.pars[:,2]
        width = np.absolute(width)
        gauss_width = 0.004 #from 0.01 [0.1,0.05,0.01,0.005,0.001,0.0005,0.0001] 0.004 is the best option for the test
        summa = 0
        factor1 = 1 / (2*len(center)) #a2F(w)=1/2N Sum{Lambda*Omega*delta(w-Omega)}
        for i in range(len(center)):
            summa += (width[i]*center[i]) * self.gaussian(x,center[i],gauss_width)
        return factor1*summa

    def gaussian(self,x, center,gauss_width):
        """
        Gaussian distribution
        ---input---
        x: x position of the gaussian graph, in this case the frequency
        center: center of the gaussian distribution, in this case the frequency of the plasmon Lorenztian approximation
        gauss_width: width of the gaussian distribution, in this case the witdh of the plasmon lorentzian aproximation
        ---output---
        g: gaussian value in 'x'
        """
        #d = self.pars[:,0] #amplitude
        #center = self.pars[:,1]
        mu = center #self.pars[:1] #center
        sigma = gauss_width #self.pars[:,2] #width
        #d = 1/(2*np.sqrt(2*np.pi)*sigma)
        d = 1/(sigma*np.sqrt(2*np.pi))
        g = d*np.exp(-(x-mu)**2 / ( 2.0 * sigma**2 )  )
        #exp(-(w_aux-w(j,l))**2.0d0/(2.0d0*broad**2.0d0))/ (broad*sqrt(twopi)) #copy from qe-5.1.0_elph
        return g

    def test_gaussian(self,w,w_q,N_q):
        """
        Test for the gaussian, the result must be the DOS.
        """
        width = 5#*self.from_cm1_to_eV #self.from_cm1_to_Hartree*2 #1,5,10 cm-1
        test_dos = 1/N_q # is this equal to len(W_q)?
        for i in range(len(w_q)): #q=6x6x6 (???) for the test data. q=50x50=250 for the plasmon
            test_dos += self.gaussian(w,w_q[i],width) #This must be the DOS, integrate this to check the gaussian.
        print("Test of the gaussian: this must be the integral of the DOS=",test_dos)
        return test_dos


    def plot_lambda(self,x):
        plt.figure(figsize=(10,6))
        #plt.figure().set_figwidth(15)
        #plt.figure().set_figheight(2)
        plt.plot(x)
        plt.show()
        pass

    def Lambda(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        """
        center = self.pars[:,1]
        width = self.pars[:,2]
        print('len(Center)=',len(center))
        Nef = self.Ne[np.where(self.energy==0.0)[0][0]]  #test
        #Nef = self.Nef
        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):
            summa1 += self.Lambda_q(width[i],center[i],Nef)
        Lambda_1=summa1/len(center) #* 33 #misterious factor... joking, this is the number of nodes in the example.
        self.lambda_2=[]
#        Frequncies = Frequencies+1
#        Frequencies = np.append(Frequencies,Frequncies, axis=0)
        if (Frequencies[0] != 0):
             Frequencies = np.append(Frequencies,Frequencies[:len(Frequencies)//3]+Frequencies[-1])
        else:
             Frequencies = np.append(Frequencies,Frequencies[1:len(Frequencies)//3]+Frequencies[-1])
#        for w in Frequencies[0:int(2*len(Frequencies)/3)]:
        for w in Frequencies[0:int(len(Frequencies))]: #edit for test
             print('Frequencie(',w,')              ', end="\r", flush=True)
             if (w == 0):
                 a2F_x = []
             else:
                 a2F_x.append(self.a2F(w)/w)
#---test--v--multiprocessing*****
        #import multiprocessing as mp #multiprocessing calculation of the a2F()
        w = Frequencies[Frequencies != 0]
        # w = np.append(w,w[:len(w)//3]+w[-1]) # to expand the frequencies to cover the widths (and some extra) for the calculations
        with mp.Pool() as pool:
           res = pool.map(self.a2F,w)
        # a2F_x = res / w #or np.divide(res, w)
        a2F_x = np.divide(res, w)
#---test--^--multiprocessing*****
        #self.lambda_2 = 2*np.trapz(a2F_x,dx=(Frequencies[-1]-Frequencies[0])/len(Frequencies)) <-- this makes lambda1 20 times bigger than lambda2 (???)
        self.lambda_2 = 2*np.trapz(a2F_x) #test <-- this lamba2 is 250 times lambda1 (that corresponds to q_x*q_y ???)
        self.plot_lambda(a2F_x)
        return Lambda_1

    def Lambda_new(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        """
        center = self.pars[:,1]
        width = self.pars[:,2]
        print('len(Center)=',len(center),'(=50*50=2500? or =50*51=2550?')
        #method 1 ----------------------------
        Nef = self.Ne[np.where(self.energy==0.0)[0][0]]  #comment for test and uncomment line under

        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):  #summa in q (6x6x6 or q_x*q_y)
            summa1 += self.Lambda_q(width[i],center[i],Nef) #1/N_q Sum_q( Lamb_q )
        Lambda_1=summa1/len(center)
        #method 2 -------------------------------
        self.lambda_2=[]
        if (Frequencies[0] != 0):
            Frequencies = np.append(Frequencies,Frequencies[:len(Frequencies)//3]+Frequencies[-1])
        else:
            Frequencies = np.append(Frequencies,Frequencies[1:len(Frequencies)//3]+Frequencies[-1])
        w = Frequencies[Frequencies != 0]
        #w = np.linspace(0.0001,0.9999,20000) #test
        with mp.Pool() as pool:
            res = pool.map(self.a2F_new,w)
        a2F_x = np.divide(res, w)
        #a2F_x *= self.from_GHz_to_eV
        #a2F_x /= self.from_cm1_to_Hartree *0.5#self.from_cm1_to_Hartree /29.9792458
        #self.lambda_2 = 2*np.trapz(a2F_x,dx=(Frequencies[-1]-Frequencies[0])/len(Frequencies))
        #self.lambda_2_test = 2*integrate.simpson(a2F_x)
        #self.lambda_2_test2 = 2*integrate.simpson(np.divide(res,w)*self.from_cm1_to_Hartree /29.9792458,w)
        #print(self.lambda_2_test2)
        self.lambda_2_test2 = 2*integrate.simpson(a2F_x,w)
        #print(self.lambda_2_test2)
        #w = np.linspace(-100,100,20000)
        #self.lambda_2_test2 = 2*integrate.simpson(np.divide(self.a2F_new(w),w),w)
        #self.plot_lambda(a2F_x)
        self.plot_lambda(res)
        return Lambda_1
    def a2F_new(self,x):
        """
        Calculates the Eliashberg spectral functions
        ---input---
        x: coordiate to calculate the Eliashberg function.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """
        center = self.pars[:,1] #*put units correctly...
        width = self.pars[:,2] #*put units the same as center
        width = np.absolute(width)
        units = create_units('2014')
        gauss_width = 5*self.from_cm1_to_eV#(units.invcm/units.Hartree) #0.00002 # test the units of this... should be aprox. 5 cm-1 (1, 5 or 10)
        summa = 0
        factor1 = 1 / (2*len(center)) #a2F(w)=1/2N Sum{Lambda*Omega*delta(w-Omega)}
        for i in range(len(center)):
            summa += (width[i]*center[i]) * self.gaussian(x,center[i],gauss_width) #check the units of the gaussian...
        return factor1*summa
class Eliashberg_test(object):
    """docstring for Eliashberg.
    This object calculates the Lambda from two methods:
    From the Eliashberg function amd fron the lambda at the q points.
    """

    def __init__(self, pars):
        units = create_units('2014')   #new way of defining units
        #super(Eliashberg, self).__init__()
        self.pars = pars

        self.from_cm1_to_Hartree = units.invcm/units.Hartree #0.0000045563352812122295
        self.from_GHz_to_Hartree = self.from_cm1_to_Hartree /29.9792458 # from GHz to Hartree
        self.from_Ry_to_Hartree = units.Ry/units.Hartree #0.5 # from Ry to Hartree
        """
        cm^-1 to Hartree: Hartree = cm^-1 / 219474.6
        Hartree to cm^-1: cm^-1 = Hartree * 219474.6

        GHz to Hartree: Hartree = GHz / 6.5796839 × 10^9
        Hartree to GHz: GHz = Hartree * 6.5796839 × 10^9

        eV to Hartree: Hartree = eV*0,0367493
        """
        self.from_cm1_to_eV = units.invcm/units.eV #0.00012398425731484318
        self.from_GHz_to_eV = 0.000004135669661004052

    # def read_Ne(self,filename="out_DOS.dat"):
    #     self.energy, self.Ne = np.loadtxt(filename,usecols=(0,1), unpack=True)
    #     numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0]+1)],dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy)) # maybe this is better??
    #     total_states = np.trapz(self.Ne,dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
    #     print("Número de elementos:",numero_de_elementos,"|total=",total_states)
    #     #test the selection
    #     print (self.Ne[:(np.where(self.energy==0.0)[0][0]+1)])
    #     print (self.Ne[:(np.where(self.energy==0.0)[0][0])+1])
    #
    #     #------^--------
    #     print ("Density of states at Fermi level per cell=",self.Ne[np.where(self.energy==0.0)])
    #     print ("Number of elements (for the factor)=",numero_de_elementos)
    #     pass

    def Lambda_q(self,gamma_q,omega_q,Nef):
        """
        Calculates the Lambda(q) functions
        ---input---
        gamma_q: widths of the lorentzian fitting of the plasmon
        omega_q: Frequencies of the plasmon, fitted to lorentzian
        Nef: Density of states at Fermi level per cell
        ---output---
        Lamb_q: Lambda(q)
        """
        Lamb_q=(1/(np.pi*Nef)) * (gamma_q/omega_q**2) #fix from omega to omega²
        return Lamb_q


    def gaussian(self,x, center,gauss_width):
        """
        Gaussian distribution
        ---input---
        x: x position of the gaussian graph, in this case the frequency
        center: center of the gaussian distribution, in this case the frequency of the plasmon Lorenztian approximation
        gauss_width: width of the gaussian distribution, in this case the witdh of the plasmon lorentzian aproximation
        ---output---
        g: gaussian value in 'x'
        """
        mu = center #self.pars[:1] #center
        sigma = gauss_width #self.pars[:,2] #width
        d = 1/(sigma*np.sqrt(2*np.pi))
        g = d*np.exp(-(x-mu)**2 / ( 2.0 * sigma**2 )  )
        #exp(-(w_aux-w(j,l))**2.0d0/(2.0d0*broad**2.0d0))/ (broad*sqrt(twopi)) #copy from qe-5.1.0_elph
        return g

    def test_gaussian(self,w,w_q,N_q):
        """
        Test for the gaussian, the result must be the DOS.
        """
        width = 5#*self.from_cm1_to_eV #self.from_cm1_to_Hartree*2 #1,5,10 cm-1
        test_dos = 1/N_q # is this equal to len(W_q)?
        for i in range(len(w_q)): #q=6x6x6 (???) for the test data. q=50x50=250 for the plasmon
            test_dos += self.gaussian(w,w_q[i],width) #This must be the DOS, integrate this to check the gaussian.
        print("Test of the gaussian: this must be the integral of the DOS=",test_dos)
        return test_dos


    def plot_lambda(self,x):
        plt.figure(figsize=(10,6))
        #plt.figure().set_figwidth(15)
        #plt.figure().set_figheight(2)
        plt.plot(x)
        plt.show()
        pass


    def Lambda_new(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        """
        center = self.pars[:,1]
        width = self.pars[:,2]
        print('len(Center)=',len(center),'(=50*50=2500? or =50*51=2550?')
        #method 1 ----------------------------
        Nef = self.Ne[np.where(self.energy==0.0)[0][0]]  #comment for test and uncomment line under
        #Nef = self.Nef
        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):  #summa in q (6x6x6 or q_x*q_y)
            summa1 += self.Lambda_q(width[i],center[i],Nef) #1/N_q Sum_q( Lamb_q )
        Lambda_1=summa1/len(center) #* 33 #misterious factor... joking, this is the number of nodes in the example.
        #method 2 -------------------------------
        self.lambda_2=[]
        # if (Frequencies[0] != 0):
        #      Frequencies = np.append(Frequencies,Frequencies[:len(Frequencies)//3]+Frequencies[-1])
        # else:
        #      Frequencies = np.append(Frequencies,Frequencies[1:len(Frequencies)//3]+Frequencies[-1])
        w = Frequencies[Frequencies != 0]
        with mp.Pool() as pool:
           res = pool.map(self.a2F_new,w)
        a2F_x = np.divide(res, w)
        #a2F_x *= self.from_GHz_to_eV
        #a2F_x /= self.from_cm1_to_Hartree *0.5#self.from_cm1_to_Hartree /29.9792458
        self.lambda_2_test2 = 2*integrate.simpson(a2F_x,w)
        self.plot_lambda(res)
        return Lambda_1
    def a2F_new(self,x):
        """
        Calculates the Eliashberg spectral functions
        ---input---
        x: coordiate to calculate the Eliashberg function.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """
        center = self.pars[:,1] #*put units correctly...
        width = self.pars[:,2] #*put units the same as center
        width = np.absolute(width)
        #units = create_units('2014')
        gauss_width = 5*self.from_cm1_to_Hartree#*(units.invcm/units.Hartree) #0.00002 # test the units of this... should be aprox. 5 cm-1 (1, 5 or 10)
        summa = 0
        factor1 = 1 / (2*len(center)) #a2F(w)=1/2N Sum{Lambda*Omega*delta(w-Omega)}
        for i in range(len(center)):
            summa += (width[i]*center[i]) * self.gaussian(x,center[i],gauss_width) #check the units of the gaussian...
        return factor1*summa
