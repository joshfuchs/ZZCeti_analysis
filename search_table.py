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



catalog = Table.read('catalog_flux_clean.txt',format='ascii')  


task = raw_input('What would you like to see? (star, range, filter, duplicates, allstars, K2) ')

if task == 'star':
    star_name = raw_input('Name of the star? ')
    for x in range(0,len(catalog['FILENAME'])):
        if catalog['FILENAME'][x].lower().__contains__(star_name.lower()) == True:
            print catalog['FILENAME'][x], catalog['DATE-OBS'][x], catalog['b10teff'][x], catalog['b10tefferr'][x],  catalog['b10logg'][x], catalog['b10loggerr'][x], catalog['EXPTIME'][x], catalog['Magnitude'][x]
            #print 'b10: ', catalog['b10teff'][x], catalog['b10logg'][x]
            #print 'alpha: ', catalog['ateff'][x], catalog['alogg'][x]
            #print 'beta: ', catalog['bteff'][x], catalog['blogg'][x]
            #print 'gamma: ', catalog['gteff'][x], catalog['glogg'][x]
            #print 'delta: ', catalog['dteff'][x], catalog['dlogg'][x]
            #print 'epsilon: ', catalog['eteff'][x], catalog['elogg'][x]
            #print 'H8: ', catalog['H8teff'][x], catalog['H8logg'][x]
            #print 'H9: ', catalog['H9teff'][x], catalog['H9logg'][x]
            #print 'H10: ', catalog['H10teff'][x], catalog['H10logg'][x]
            #print catalog[x]

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
    for x in np.arange(len(catalog['ateff'])):
        if catalog['ateff'][x] <= tupper and catalog['ateff'][x] >= tlower and catalog['alogg'][x] <= gupper and catalog['alogg'][x] >= glower:
            print catalog['FILENAME'][x], catalog['DATE-OBS'][x], catalog['ateff'][x], catalog['alogg'][x]#,catalog['DAV'][x]

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

if task == 'allstars':
    for x in range(0,len(catalog['FILENAME'])):
        print catalog['WD'][x]

if task == 'K2':
    for x in range(0,len(catalog['FILENAME'])):
        if catalog['K2'][x] == 1:
            print catalog['WD'][x], catalog['DATE-OBS'][x], catalog['b10teff'][x], catalog['b10logg'][x], catalog['g10teff'][x], catalog['g10logg'][x]
