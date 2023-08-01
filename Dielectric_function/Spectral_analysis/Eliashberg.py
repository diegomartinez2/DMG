import matplotlib.pyplot as plt
import numpy as np

class Eliashberg(object):
    """docstring for Eliashberg."""

    def __init__(self, pars):
        #super(Eliashberg, self).__init__()
        self.pars = pars
        self.from_cm1_to_Hartree = 4.55633e-6 # from cm-1 to Hartree
        self.from_GHz_to_Hartree = self.from_cm1_to_Hartree /29.9793 # from GHz to Hartree

    def read_Ne(self,filename="out_DOS.dat"):
        self.energy, self.Ne = np.loadtxt(filename,usecols=(0,1), unpack=True)
#        self.Ne_0[np.where(energy==0.0)]
        numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0])+1],dx=(self.energy[9]-self.energy[0])/10)
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
        center = np.absolute(center)
        width = self.pars[:,2]
        width = np.absolute(width)
        gauss_width = 0.01
        summa = 0
        factor1 = 1 / (2*len(center))
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
#        Nef = self.Ne[np.where(self.energy==0.0)[0][0]]  #test
        Nef = self.Nef
        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):
            summa1 += self.Lambda_q(width[i],center[i],Nef)
        Lambda_1=summa1/len(center)
        self.lambda_2=[]
#        Frequncies = Frequencies+1
#        Frequencies = np.append(Frequencies,Frequncies, axis=0)
#        for w in Frequencies[0:int(2*len(Frequencies)/3)]:
        for w in Frequencies[0:int(2*len(Frequencies))]:
            print('Frequencie(',w,')              ', end="\r", flush=True)
            if (w == 0):
                a2F_x = []
            else:
                a2F_x.append(self.a2F(w)/w)
        self.lambda_2 = 2*np.trapz(a2F_x,dx=(Frequencies[9]-Frequencies[0])/9)
        self.plot_lambda(a2F_x)
        return Lambda_1
