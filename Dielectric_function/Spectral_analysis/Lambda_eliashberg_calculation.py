import Spectral_plasmons_analysis
import numpy as np


def main(arg):
    filelist=('xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal',
    'xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay',
    'xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg','xbh','xbi','xbj','xbk','xbl',
    'xbm','xbn','xbo','xbp','xbq','xbr','xbs','xbt','xbu','xbv','xbw','xbx')
    for namefile in filelist:
        plasmon = Spectral_plasmons_analysis.Plasmon_analysis(arg,namefile)
        data, frequencies, qx= plasmon.load_data()
        print (np.shape(data),"=(51,5001)?")
        plasmon.fitting_Lorentz(frequencies,data)
        if (namefile == filelist[0]):
            pars = plasmon.pars
            #print('pars_1=',pars)
        else:
            pars = np.concatenate((pars,plasmon.pars), axis=0)
            #print('pars_n=',pars)
    #np.savetxt('pars_txt.dat', pars) #no negative values
    superconductor = Spectral_plasmons_analysis.Eliashberg(pars)
    superconductor.read_Ne()
    lambda_1 = superconductor.Lambda(frequencies)
    np.savetxt('Lambda.txt', (lambda_1))
    print('Lambda_1=',np.sum(lambda_1)/len(lambda_1))
    #np.savetxt('Lambda_from_a2F.txt', np.array((frequencies[1:],superconductor.lambda_2)).T, header='frequencies,Lambda')
    print('Lambda_2=',np.sum(superconductor.lambda_2))
    print("len(freq[1:]),len(superconductor.lambda_2)")
    print(len(frequencies[1:]),len(superconductor.lambda_2))


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
