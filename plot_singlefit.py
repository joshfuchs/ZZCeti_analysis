'''
Written May 2016 by JTF

Plots normalized spectrum with normalized model
'''


import matplotlib.pyplot as plt
import numpy as np
import os

#Read in normalized spectrum
#night1 = '/srv/two/jtfuchs/Interpolated_Models/Bergeron_new/bottom10000_700'
#os.chdir(night1)
spec1 = 'norm_PG1149p058toSDSSg_06-01_4.0.txt'
wdlamb,wdinten,wdsigma = np.genfromtxt(spec1,unpack=True,skip_header=1)
title = spec1[11:-9] #+ ' : ' + spec1[-14:-9]
filename = spec1[11:-9]#+'_'+spec1[-14:-9]
##os.chdir('../')
#Read in normalized model
#night2 = '/srv/two/jtfuchs/Interpolated_Models/Koester_ML2alpha08/bottom10000_700'
#os.chdir(night2)
spec2 = 'model_PG1149p058toSDSSg_06-01_4.0.txt'
#spec2 = 'norm_wcftb.WD1116p026_930_blue_flux_master_noext_03-14_3.91.txt'
#title = 'Teff: ' + spec2[2:7] + ' K, ' + 'log(g): ' + spec2[10] + '.' + spec2[11:13]
modlamb,modinten = np.genfromtxt(spec2,unpack=True)#,skip_header=1)

#plot a second spectrum if you wish
#spec3 = 'model_WD1422p095_930_blue_modelflux_07-13_6.65.txt'
#lamb3,inten3 = np.genfromtxt(spec3,unpack=True,skip_header=1)

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
plt.plot(wdlambbeta-4862.6,wdintenbeta,'k')#,label='Bergeron')
plt.plot(modlambbeta-4862.6, modintenbeta,'r')#,label='Koester')
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
plt.xlim(-150,150)
plt.ylim(0,3.0)
plt.legend()
#plt.title(title)
#plt.savefig(filename+'.png',format='png')#,dpi=12000)
plt.show()
