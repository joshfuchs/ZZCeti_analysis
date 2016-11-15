'''
Written by JT Fuchs, UNC - Chapel Hill

Read in fitting_solutions.txt file and plot all the solutions in the logg-Teff plane.

'''
import numpy as np
import matplotlib.pyplot as plt


'''
Some code here to go find all the text files
'''

'''
fitting_solutions.txt contains: blue filename, red filenames, best model, best Teff, best logg, FWHM, best chi-square, date-time of fit.
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

plt.clf()
plt.plot(Teff,logg,'bo')
plt.xlim(15500,9500)
plt.ylim(9.75,7.0)
#Plot Gianninas 2011 instability strip
plt.plot([14000.,11600.],[8.82,7.36],'r')
plt.plot([11800,10600],[8.88,7.35],'r')
plt.show()

