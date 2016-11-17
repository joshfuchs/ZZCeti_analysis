'''
Written by JT Fuchs, UNC - Chapel Hill

Read in fitting_solutions.txt file and print requested information

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

task = raw_input('What would you like to see? (star)')

if task == 'star':
    star_name = raw_input('Name of the star? ')
    for x in blue_file:
        if x.lower().__contains__(star_name.lower()) == True:
            print x, Teff[blue_file.index(x)], logg[blue_file.index(x)]
