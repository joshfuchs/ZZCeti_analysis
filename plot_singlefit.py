'''
Written May 2016 by JTF

Plots normalized spectrum with normalized model
'''


import matplotlib.pyplot as plt
import numpy as np
import os

#Read in normalized spectrum

os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/soar/pipeline/2016/2016-10-11/GOODMAN/VPHAS1813-2138')
spec1 = 'norm_VPHAS1813-2138_930_blue_flux_10-13_3.35.txt'

#os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/soar/pipeline/2016/2016-08-20/GOODMAN/ZZCETIS')
#spec1 = 'norm_VPHAS1813-2138_930_blue.FIT_10-13_4.83.txt'

wdlamb,wdinten, wdsigma = np.genfromtxt(spec1,unpack=True)


#spec2 = 'norm_WD1422p095_930_blue_flux_CD32_1_07-08_4.75_narrow.txt'
#lamb2,inten2,sigma2 = np.genfromtxt(spec2,unpack=True)

#Go get the normalized, convoled models
#os.chdir('/srv/two/jtfuchs/Interpolated_Models/Koester_ML2alpha06/bottom11500_750')
#os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/clemens/students/group/modelfitting/Bergeron_New')
os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/clemens/students/group/modelfitting/Koester_06')

#spec2 = 'da7000.0_625.0VPHAS1813-2138_930_blue.FIT_4.83_norm.txt'
spec2 = 'da12000.0_625.0VPHAS1813-2138_930_blue_flux_3.35_norm.txt'
title = 'Teff: ' + spec2[2:7] + ' K, ' + 'log(g): ' + spec2[10] + '.' + spec2[11:13]
modlamb,modinten = np.genfromtxt(spec2,unpack=True,skip_header=1)

#plot a second spectrum if you wish
#spec3 = 'model_WD1422p095_930_blue_modelflux_07-13_6.65.txt'
#lamb3,inten3 = np.genfromtxt(spec3,unpack=True,skip_header=1)

#Break up spectrum into individual lines for plotting

betalow = 4812.
betahigh = 4911.
gammalow = 4281.
gammahigh = 4374.
deltalow = 4052. #4031
deltahigh = 4152. #4191
epsilonlow = 3940. #3925
epsilonhigh = 4007. # 4021
heightlow = 3865. #3859
heighthigh = 3920. # 3925
hninelow = 3815. #3815
hninehigh = 3855. #3855
htenlow = 3785. #3785
htenhigh = 3815.



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
#plt.savefig('VPHAS1813-2138_12000_625.png',format='png')#,dpi=12000)
plt.show()
