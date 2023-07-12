import Spectral_plasmons_analysis



def main():
    filelist=('xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal',
    'xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay',
    'xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg','xbh','xbi','xbj','xbk','xbl',
    'xbm','xbn','xbo','xbp','xbq','xbr','xbs','xbt','xbu','xbv','xbw','xbx')
    pars = []
    for i in filelist:
        data, frequencies, qx= polaron.load_data()
        print (np.shape(data),"=(51,5001)?")
        polaron.fitting_Lorentz(frequencies,data)
        np.append(pars,polaron.pars)

        

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
