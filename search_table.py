'''
Written by JT Fuchs, UNC - Chapel Hill

Read in fitting_solutions.txt file and print requested information

'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

'''
Some code here to go find all the text files
'''

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
'''
arr = np.genfromtxt('all_teff_logg.txt',dtype=None,delimiter='\t')
blue_file, date, Teff, Tefferr, logg, loggerr = [], [], np.zeros(len(arr)), np.zeros(len(arr)), np.zeros(len(arr)), np.zeros(len(arr))

for m in np.arange(len(arr)):
    blue_file.append(arr[m][0])
    date.append(arr[m][1])
    Teff[m] = arr[m][2]
    Tefferr[m] = arr[m][3]
    logg[m] = arr[m][4]
    loggerr[m] = arr[m][5]
'''

catalog = Table.read('full_catalog.txt',format='ascii')  


task = raw_input('What would you like to see? (star, range, filter, duplicates) ')

if task == 'star':
    star_name = raw_input('Name of the star? ')
    for x in range(0,len(catalog['FILENAME'])):
        if catalog['FILENAME'][x].lower().__contains__(star_name.lower()) == True:
            print catalog['FILENAME'][x], catalog['DATE-OBS'][x], catalog['Teff'][x], catalog['Tefferr'][x],  catalog['logg'][x], catalog['loggerr'][x]

if task == 'range':
    trange = raw_input('Would you like to set a temperature range? (yes/no) ')
    if trange == 'yes':
        tlower = float(raw_input('Lower temperature limit: '))
        tupper = float(raw_input('Upper temperature limit: '))
    else:
        tlower = 1.
        tupper = 100000.
    grange = raw_input('Would you like to set a log(g) range? (yes/no) ')
    if grange == 'yes':
        glower = float(raw_input('Lower log(g) limit: '))
        gupper = float(raw_input('Upper log(g) limit: '))
    else:
        glower = 1.
        gupper = 20.
    for x in np.arange(len(catalog['Teff'])):
        if catalog['Teff'][x] <= tupper and catalog['Teff'][x] >= tlower and catalog['logg'][x] <= gupper and catalog['logg'][x] >= glower:
            print catalog['FILENAME'][x], catalog['DATE-OBS'][x], catalog['Teff'][x], catalog['Tefferr'][x], catalog['logg'][x], catalog['loggerr'][x]#,catalog['DAV'][x]

if task == 'filter':
    which_filter = raw_input('Name of filter (DAV,K2,outburst): ')
    yesno = raw_input('Yes or no?(0,1) ')
    for x in range(0,len(catalog['FILENAME'])):
        #print str(catalog[which_filter][x]), yesno, str(catalog[which_filter][x]) == yesno
        if str(catalog[which_filter][x]) == yesno:
            print catalog['FILENAME'][x], catalog['DATE-OBS'][x], catalog['Teff'][x], catalog['logg'][x], catalog[which_filter][x]


if task == 'duplicates':
    for x in range(0,len(catalog['FILENAME'])):
        for y in range(0,len(catalog['FILENAME'])):
            if np.abs(catalog['RA'][x]-catalog['RA'][y]) < 0.01 and catalog['DATE-OBS'][x] != catalog['DATE-OBS'][y]:
                print catalog['FILENAME'][x], catalog['DATE-OBS'][x], catalog['FILENAME'][y], catalog['DATE-OBS'][y]
