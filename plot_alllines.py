'''
Written May 2016 by JTF

Plots normalized spectrum with normalized model
'''


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
from glob import glob

#Get list of directories in current directory. Only keep if it is a date by looking for ./2
directories = [x[0] for x in os.walk('./') if x[0][0:3]=='./2']

plots_pdf = PdfPages('all_plotfits.pdf')

for xdir in directories:
    os.chdir(xdir)
    wds = sorted(glob('norm_*txt'))
    models = sorted(glob('model_*txt'))
    #print xdir
    #print wds
    #print models
    #print '\n'
    for y in wds:
        wdlamb,wdinten, wdsigma = np.genfromtxt(y,unpack=True,skip_header=1)
        modlamb,modinten = np.genfromtxt(models[wds.index(y)],unpack=True)
        #print y, models[wds.index(y)]
        title = str(xdir) + ': ' + str(y[10:y.find('930')])

        #Break up spectrum into individual lines for plotting
        alphalow = 6413.
        alphahigh = 6713.
        betalow = 4721.
        betahigh = 5001.
        gammalow = 4220.
        gammahigh = 4460.
        deltalow = 4031. #4031
        deltahigh = 4191. #4191
        epsilonlow = 3925. #3925
        epsilonhigh = 4021. # 4021
        heightlow = 3859. #3859
        heighthigh = 3925. # 3925
        hninelow = 3815. #3815
        hninehigh = 3855. #3855
        htenlow = 3785. #3785
        htenhigh = 3815.
        
        if wdlamb.max() > 6000.:
            wdlambalpha, wdintenalpha = wdlamb[np.where((wdlamb > alphalow) & (wdlamb < alphahigh+1))], wdinten[np.where((wdlamb > alphalow) & (wdlamb < alphahigh+1.))]
            modlambalpha, modintenalpha = modlamb[np.where((modlamb > alphalow) & (modlamb < alphahigh+1))], modinten[np.where((modlamb > alphalow) & (modlamb < alphahigh+1.))]
        
        wdlambbeta, wdintenbeta = wdlamb[np.where((wdlamb > betalow) & (wdlamb < betahigh+1))], wdinten[np.where((wdlamb > betalow) & (wdlamb < betahigh+1.))]
        modlambbeta, modintenbeta = modlamb[np.where((modlamb > betalow) & (modlamb < betahigh+1))], modinten[np.where((modlamb > betalow) & (modlamb < betahigh+1.))]
        
        wdlambgamma, wdintengamma = wdlamb[np.where((wdlamb > gammalow) & (wdlamb < gammahigh+1))], wdinten[np.where((wdlamb > gammalow) & (wdlamb < gammahigh+1.))]##wdlamb[349:460], wdinten[349:460]
        modlambgamma, modintengamma = modlamb[np.where((modlamb > gammalow) & (modlamb < gammahigh+1))], modinten[np.where((modlamb > gammalow) & (modlamb < gammahigh+1.))]##modlamb[349:459], modinten[349:459]
        
        wdlambdelta, wdintendelta = wdlamb[np.where((wdlamb > deltalow) & (wdlamb < deltahigh+1))], wdinten[np.where((wdlamb > deltalow) & (wdlamb < deltahigh+1.))]##wdlamb[230:349], wdinten[230:349]
        modlambdelta, modintendelta = modlamb[np.where((modlamb > deltalow) & (modlamb < deltahigh+1))], modinten[np.where((modlamb > deltalow) & (modlamb < deltahigh+1.))]##odlamb[230:348], modinten[230:348]
        
        wdlambepsilon, wdintenepsilon = wdlamb[np.where((wdlamb > epsilonlow) & (wdlamb < epsilonhigh+1))], wdinten[np.where((wdlamb > epsilonlow) & (wdlamb < epsilonhigh+1.))]##wdlamb[150:230], wdinten[150:230]
        modlambepsilon, modintenepsilon = modlamb[np.where((modlamb > epsilonlow) & (modlamb < epsilonhigh+1))], modinten[np.where((modlamb > epsilonlow) & (modlamb < epsilonhigh+1.))]##modlamb[150:229], modinten[150:229]
        
        wdlamb8, wdinten8 = wdlamb[np.where((wdlamb > heightlow) & (wdlamb < heighthigh+1))], wdinten[np.where((wdlamb > heightlow) & (wdlamb < heighthigh+1.))]##wdlamb[85:149], wdinten[85:149]
        modlamb8, modinten8 = modlamb[np.where((modlamb > heightlow) & (modlamb < heighthigh+1))], modinten[np.where((modlamb > heightlow) & (modlamb < heighthigh+1.))]##modlamb[85:148], modinten[85:148]
        
        wdlamb9, wdinten9 = wdlamb[np.where((wdlamb > hninelow) & (wdlamb < hninehigh+1))], wdinten[np.where((wdlamb > hninelow) & (wdlamb < hninehigh+1.))]##wdlamb[37:85], wdinten[37:85]
        modlamb9, modinten9 = modlamb[np.where((modlamb > hninelow) & (modlamb < hninehigh+1))], modinten[np.where((modlamb > hninelow) & (modlamb < hninehigh+1.))]##modlamb[37:84], modinten[37:84]
        
        wdlamb10, wdinten10 = wdlamb[np.where((wdlamb > htenlow) & (wdlamb < htenhigh+1))], wdinten[np.where((wdlamb > htenlow) & (wdlamb < htenhigh+1.))]##wdlamb[:37], wdinten[:37]
        modlamb10, modinten10 = modlamb[np.where((modlamb > htenlow) & (modlamb < htenhigh+1))], modinten[np.where((modlamb > htenlow) & (modlamb < htenhigh+1.))]##modlamb[:37], modinten[:37]



        plt.clf()
        if wdlamb.max() > 6000.:
            plt.plot(wdlambalpha-6564.6,wdintenalpha-0.3,'k')
            plt.plot(modlambalpha-6564.6,modintenalpha-0.3,'r')
        plt.plot(wdlambbeta-4862.6,wdintenbeta,'k')
        plt.plot(modlambbeta-4862.6, modintenbeta,'r')
        plt.plot(wdlambgamma-4341.6, wdintengamma+0.3,'k')
        plt.plot(modlambgamma-4341.6, modintengamma+0.3,'r')
        plt.plot(wdlambdelta-4102.9, wdintendelta+0.6,'k')
        plt.plot(modlambdelta-4102.9, modintendelta+0.6,'r')
        plt.plot(wdlambepsilon-3971.1, wdintenepsilon+0.9,'k')
        plt.plot(modlambepsilon-3971.1, modintenepsilon+0.9,'r')
        plt.plot(wdlamb8-3890.1,wdinten8+1.2,'k')
        plt.plot(modlamb8-3890.1, modinten8+1.2,'r')
        plt.plot(wdlamb9-3836.4,wdinten9+1.5,'k')
        plt.plot(modlamb9-3836.4, modinten9+1.5,'r')
        plt.plot(wdlamb10-3798.9,wdinten10+1.8,'k')
        plt.plot(modlamb10-3798.9, modinten10+1.8,'r')
        plt.xlabel('Relative Wavelength (Ang.)')
        plt.ylabel('Relative Flux')
        plt.title(title)
        plt.xlim(-150,150)
        plt.ylim(0,3.5)
        #plt.savefig('VPHAS1813-2138_12000_625.png',format='png')#,dpi=12000)
        plots_pdf.savefig()
        #plt.show()
    os.chdir('../')


plots_pdf.close()
