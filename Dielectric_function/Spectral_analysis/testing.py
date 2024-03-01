import numpy as np
import matplotlib.pyplot as plt
#import Eliashberg_new as Eliashberg
import Eliashberg_new_clean as Eliashberg
w_a, a2F_1DP, a2F_HPI, a2F_HPII, a2F_Total = np.loadtxt("a2F_grupo.txt", usecols=(0, 1, 2, 3, 4), unpack=True)
w, Lambda_1DP, Lambda_HPI, Lambda_HPII, Lambda_Total = np.loadtxt("Lambda_grupo.txt", usecols=(0, 1, 2, 3, 4), unpack=True)

for i in range(len(w)):
	print(w[i],(Lambda_1DP[i]+Lambda_HPI[i]+Lambda_HPII[i]), Lambda_Total[i],"::",(Lambda_1DP[i]+Lambda_HPI[i]+Lambda_HPII[i])-Lambda_Total[i])

plt.plot(w, (Lambda_1DP+Lambda_HPI+Lambda_HPII), label="Summa")
plt.plot(w, Lambda_Total, label="Total")
plt.legend()
plt.show()
print(len(a2F_Total),"::",len(w),":",len(w_a))
print(w[0],w_a[1])
mask = Lambda_Total <= 0
print(mask)
print(Lambda_Total[mask])
mask = w <= 0
print("w mask=",mask)
print(w[mask])
def T_c(w,mu_par,lambda_t,a2F):
        """
        Allen-Dynes formula (Advanced McMillan's equation and its application for the analysis of highly-compressed superconductors)
        Tc=W_log/1.2*exp(-1.04*(1+self.lambda_2)/(self.lambda_2-mu_*(1+0.62*self.lambda_2)))
        Note this is an approximation formula and is not accurate for some superconductors
        formula from DOI 10.1007/s10948-017-4295-y
        """
        AllenDynes = True
        out2=-1.04*(1+lambda_t)/(lambda_t-mu_par*(1+0.62*lambda_t))
        if (AllenDynes):
            print ("f1,f2=",f_1(mu_par,lambda_t),f_2(mu_par,w,lambda_t,a2F))
            out=f_1(mu_par,lambda_t)*f_2(mu_par,w,lambda_t,a2F)*w_log(w,lambda_t,a2F)/1.2 * np.exp(out2)
        else:
            out=w_log(w,lambda_t,a2F)/1.2 * np.exp(out2)

        return out/8.617333262e-5 #Boltzman constant in eV/K

def w_log(w,lambda_t,a2F):
        """
        w_log calculated from Eliashberg
        """
        #eV_to_K=11604
        #eV_to_K = 11604.5250061657

        #mask = self.w >=  0
        #mask = self.w >  0
        #w = self.w#[mask]#*eV_to_K
        #print("max(w)=",max(w))
        for i in range(len(w)):
          if w[i]<0:
            print("error en w<0",i)
            exit()
        a2F_w = a2F/w #np.divide(a2F, w)
        log_w = np.log(w)
        w_log = np.exp((2.0/lambda_t)*np.trapz(a2F_w*log_w,w))
            #integrate.simpson((np.divide(self.a2F_new(self.w), self.w)*np.log(self.w)),self.w))
            #np.trapz(a2F_w*np.log(w)),w)
        return w_log

def w_2(w,lambda_t,a2F):
        #eV_to_K=11604
        #w = w#*eV_to_K
        a2F_w = a2F*w
        w_2 = np.sqrt((2.0/lambda_t)*
            #integrate.simpson(
            np.trapz(a2F_w,w))
        return w_2

def f_1(mu_par,lambda_t):
        LAMBDA_temp = 2.46*(1+3.8*mu_par)
        return np.power(1+np.power(lambda_t/LAMBDA_temp,3/2),1/3)

def f_2(mu_par,w,lambda_t,a2F):
        LAMBDA_temp =1.82*(1+6.3*mu_par)*(w_2(w,lambda_t,a2F)/w_log(w,lambda_t,a2F))
        return 1 + (( w_2(w,lambda_t,a2F)/w_log(w,lambda_t,a2F) - 1) * lambda_t**2)/(
        (lambda_t**2) + (LAMBDA_temp**2))

a2F = a2F_Total[1:]
#Lambda = Lambda_1DP+Lambda_HPI+Lambda_HPII
Lambda = Lambda_Total
mask = Lambda > 0
print(mask)

Tc = T_c(w[mask],0.1,Lambda[-1],a2F[mask])
print("T_c=",Tc,"(K) con Lambda_total")
Lambda = Lambda_1DP+Lambda_HPI+Lambda_HPII
#Lambda = Lambda_Total
mask = Lambda > 0
print(mask)

Tc = T_c(w[mask],0.1,Lambda[-1],a2F[mask])
print("T_c=",Tc,"(K) sumando las Lambdas")
