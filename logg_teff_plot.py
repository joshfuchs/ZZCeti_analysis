'''
Written by JT Fuchs, UNC - Chapel Hill

Read in fitting_solutions.txt file and plot all the solutions in the logg-Teff plane.

'''
import numpy as np
import matplotlib.pyplot as plt


'''
fitting_solutions.txt contains: blue filename, red filenames, best model, best Teff, best logg, FWHM, best chi-square, date-time of fit.
'''


'''
arr = np.genfromtxt('all_fit_solutions.txt',dtype=None,delimiter='\t')
blue_file, red_file, best_model, Teff, logg, fwhm, chi_square, date = [], [], [], np.zeros(len(arr)), np.zeros(len(arr)), np.zeros(len(arr)), np.zeros(len(arr)), []

for m in np.arange(len(arr)):
    blue_file.append(arr[m][0])
    red_file.append(arr[m][1])
    best_model.append(arr[m][2])
    Teff[m] = arr[m][3]
    logg[m] = arr[m][4]
    fwhm[m] = arr[m][5]
    chi_square[m] = arr[m][6]
    date.append(arr[m][7])

logg = logg / 1000.
'''
arr = np.genfromtxt('all_teff_logg.txt',dtype=None,delimiter='\t')
wd, date, Teff, Tefferr, logg, loggerr = [], [], np.zeros(len(arr)), np.zeros(len(arr)), np.zeros(len(arr)), np.zeros(len(arr))

for m in np.arange(len(arr)):
    wd.append(arr[m][0])
    date.append(arr[m][1])
    Teff[m] = arr[m][2]
    Tefferr[m] = arr[m][3]
    logg[m] = arr[m][4]
    loggerr[m] = arr[m][5]
    


#Print solutions that meet particular value
#for x in np.arange(len(Teff)):
#    if Teff[x] == 15000. or Teff[x] == 10000. or logg[x] == 9.5 or logg[x] == 7.0:
#        print Teff[x], logg[x], blue_file[x]
plt.clf()
plt.plot(Teff,logg,'bo')
plt.xlim(15500,9500)
plt.ylim(9.75,7.0)
#Plot Gianninas 2011 instability strip
plt.plot([14000.,11600.],[8.82,7.36],'r')
plt.plot([11800,10600],[8.88,7.35],'r')
#Plot 4 horsemen
#plt.plot(12740,8.0,'ro') #WD1150-153 from 2015-02-10
#plt.plot(12790,8.0,'ro') #GD 165 from 2015-04-26
#plt.plot(12780,8.05,'ro') #GD 133 from 2015-05-21
#plt.plot(12790,8.05,'ro') #L19-2 from 2015-05-22
#plot ensemble from Chris
#plt.plot(12301,8.063,'ro') #SDSS J1151+0525, EPIC 201802933
#plt.plot(14280,8.502,'ro') #SDSS J0837+1856, EPIC 211914185 low S/N?
plt.plot(12787,8.069,'ro') #L19-2 from 2015-05-22
plt.plot(12788,8.015,'ro') #GD 165 from 2015-04-26
plt.plot(12782,8.054,'ro')#GD 133
plt.plot(12722,8.009,'ro')#WD 1150

plt.plot(12615,8.037,'mo') #GD 133
plt.plot(12691,8.09,'mo') #L19-2
plt.plot(12629,8.098,'mo') #L19-2
plt.plot(12603,8.045,'mo') #L19-2
plt.plot(12516,8.018,'mo') #GD165
#outbursters
plt.plot(11142,8.112,'k^')
plt.plot(11315,8.068,'k^')
plt.plot(11260,8.128,'k^')
plt.plot(11505,8.125,'k^')
plt.plot(11187,8.127,'k^')

#plt.plot(12457,7.955,'ro') #R548 from 2015-08-23
#plt.plot(11832,7.944,'ro') #SDSS J0840+1303, EPIC 228682478 from 2016-01-07
#plt.plot(12609,8.019,'ro') #SDSS J0051+0339, EPIC 220347759
plt.show()

