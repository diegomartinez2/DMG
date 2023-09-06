import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from ase.units import create_units
from scipy.integrate import quad

class Eliashberg(object):
    """docstring for Eliashberg."""

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
        total_states = np.trapz(self.Ne,dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
        print("Número de elementos:",numero_de_elementos,"|total=",total_states)
        #test the selection
        print (self.Ne[:(np.where(self.energy==0.0)[0][0]+1)])
        print (self.Ne[:(np.where(self.energy==0.0)[0][0])+1])

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
        width = 1 #1,5,10 cm-1
        test_dos = 1/N_q # is this equsl to len(W_q)?
        for i in range(len(w_q)): #q=6x6x6 (???) for the test data. q=50x50=250 for the plasmon
            test_dos += self.gaussian(w,w_q[i],width) #This must be the DOS, integrate this to check the gaussian.
        print("Test of the gaussian: this must be the integral of the DOS=",)


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
        Nef = self.Ne[np.where(self.energy==0.0)[0][0]]  #test
        #Nef = self.Nef
        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):  #summa in q (6x6x6 or q_x*q_y)
            summa1 += self.Lambda_q(width[i],center[i],Nef) #1/N_q Sum_q( Lamb_q )
        Lambda_1=summa1/len(center) #* 33 #misterious factor... joking, this is the number of nodes in the example.
        #method 2 -------------------------------
        self.lambda_2=[]
        if (Frequencies[0] != 0):
             Frequencies = np.append(Frequencies,Frequencies[:len(Frequencies)//3]+Frequencies[-1])
        else:
             Frequencies = np.append(Frequencies,Frequencies[1:len(Frequencies)//3]+Frequencies[-1])
        w = Frequencies[Frequencies != 0]
        with mp.Pool() as pool:
           res = pool.map(self.a2F_new,w)
        a2F_x = np.divide(res, w)
        self.lambda_2 = 2*np.trapz(a2F_x) #test <-- this lamba2 is 250 times lambda1 (that corresponds to q_x*q_y ???)
        self.plot_lambda(a2F_x)
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
        gauss_width = 0.004 # test the units of this... should be aprox. 5 cm-1 (1, 5 or 10)
        summa = 0
        factor1 = 1 / (2*len(center)) #a2F(w)=1/2N Sum{Lambda*Omega*delta(w-Omega)}
        for i in range(len(center)):
            summa += (width[i]*center[i]) * self.gaussian(x,center[i],gauss_width) #check the units of the gaussian...
        return factor1*summa
