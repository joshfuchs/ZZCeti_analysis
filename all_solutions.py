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

arr = np.genfromtxt('fitting_solutions.txt',dtype=None,delimiter='\t')
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
plt.xlim(13000,11700)
plt.ylim(8.2,7.9)
plt.show()

