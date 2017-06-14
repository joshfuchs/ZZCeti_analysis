'''
Search through table and find duplicate observations. Use desired criteria to combine observations.
Remove duplicate rows.

The criterion to find duplicates is a difference in the decliation of less than 0.1 degrees. If there are targets with declinations closer than this to each other, you should change this.

Written by JT Fuchs, 7 April 2017

INPUT: Table of results, sorted by RA.

OUTPUT: Table of results with merged duplicates

'''

import numpy as np
from astropy.table import Table

catalog = Table.read('windycatalog_master.txt',format='ascii')
#print catalog.colnames
#exit()

duplicate_rows = []

x = 0
while x < len(catalog['WD']):
    #print ''
    if catalog['Duplicates'][x] == 1:
        #print catalog['WD'][x]
        duplicate_rows.append(x) #Include this if you want the rows of ALL duplicates
        if (len(catalog['WD']) - x) < 10:
            searchrange = (len(catalog['WD']) - x)
        else:
            searchrange = 11
        for m in range(1,searchrange):
            #print m
            #print catalog['DEC2'][x],catalog['DEC2'][x+m]
            if np.abs(catalog['DEC2'][x]-catalog['DEC2'][x+m]) < 0.1:
                #print 'found match'
                y = m
                #print x+y
                duplicate_rows.append(x+y)
        #print x, y, m
        #print catalog['WD'][x:x+y+1]
        '''
        #Combine the values for each teff and logg measurement. Replace the values from the first row and remove duplicate rows
        catalog['a10teff'][x] = np.mean(catalog['a10teff'][x:x+y+1])
        catalog['a10logg'][x] = np.mean(catalog['a10logg'][x:x+y+1])

        catalog['b10teff'][x] = np.mean(catalog['b10teff'][x:x+y+1])
        catalog['b10logg'][x] = np.mean(catalog['b10logg'][x:x+y+1])

        catalog['g10teff'][x] = np.mean(catalog['g10teff'][x:x+y+1])
        catalog['g10logg'][x] = np.mean(catalog['g10logg'][x:x+y+1])


        catalog['b9teff'][x] = np.mean(catalog['b9teff'][x:x+y+1])
        catalog['b9logg'][x] = np.mean(catalog['b9logg'][x:x+y+1])

        catalog['b8teff'][x] = np.mean(catalog['b8teff'][x:x+y+1])
        catalog['b8logg'][x] = np.mean(catalog['b8logg'][x:x+y+1])


        catalog['ateff'][x] = np.mean(catalog['ateff'][x:x+y+1])
        catalog['alogg'][x] = np.mean(catalog['alogg'][x:x+y+1])

        catalog['bteff'][x] = np.mean(catalog['bteff'][x:x+y+1])
        catalog['blogg'][x] = np.mean(catalog['blogg'][x:x+y+1])

        catalog['gteff'][x] = np.mean(catalog['gteff'][x:x+y+1])
        catalog['glogg'][x] = np.mean(catalog['glogg'][x:x+y+1])

        catalog['dteff'][x] = np.mean(catalog['dteff'][x:x+y+1])
        catalog['dlogg'][x] = np.mean(catalog['dlogg'][x:x+y+1])

        catalog['eteff'][x] = np.mean(catalog['eteff'][x:x+y+1])
        catalog['elogg'][x] = np.mean(catalog['elogg'][x:x+y+1])

        catalog['H8teff'][x] = np.mean(catalog['H8teff'][x:x+y+1])
        catalog['H8logg'][x] = np.mean(catalog['H8logg'][x:x+y+1])

        catalog['H9teff'][x] = np.mean(catalog['H9teff'][x:x+y+1])
        catalog['H9logg'][x] = np.mean(catalog['H9logg'][x:x+y+1])

        catalog['H10teff'][x] = np.mean(catalog['H10teff'][x:x+y+1])
        catalog['H10logg'][x] = np.mean(catalog['H10logg'][x:x+y+1])
        '''
        #if m == 0:
        #    catalog.remove_row(x+y)
        #else:
        #    catalog.remove_rows([x+1,x+y+1])




        x += y + 1
    else:
        x += 1

print duplicate_rows
#catalog.remove_rows(duplicate_rows)
#catalog.write('test_nodup.txt',format='ascii')


#Create table
#Save WD, Filename, Othername, DAte obs, Ra, DEC, Airmass, SNR, EXPTIME, Seeing, b10teff, b10logg, Magnitude

#long
#newtable = Table([catalog['WD'][duplicate_rows],catalog['OtherName'][duplicate_rows],catalog['DATE-OBS'][duplicate_rows],catalog['RA'][duplicate_rows],catalog['DEC'][duplicate_rows],catalog['Airmass'][duplicate_rows],catalog['SNR'][duplicate_rows],catalog['EXPTIME'][duplicate_rows],catalog['Seeing'][duplicate_rows],catalog['MOUNTAZ'][duplicate_rows],catalog['WINDAZ'][duplicate_rows],catalog['WINDSPEED'][duplicate_rows],catalog['b10teff'][duplicate_rows],catalog['b10logg'][duplicate_rows],catalog['Magnitude'][duplicate_rows]],names=['WD','OTHERNAME','DATE-OBS','RA','DEC','AIRMASS','SNR','EXPTIME','SEEING','MOUNTAZ','WINDAZ','WINDSPEED','b10teff','b10logg','MAGNITUDE'])

#short
newtable = Table([catalog['WD'][duplicate_rows],catalog['DATE-OBS'][duplicate_rows],catalog['Airmass'][duplicate_rows],catalog['SNR'][duplicate_rows],catalog['EXPTIME'][duplicate_rows],catalog['Seeing'][duplicate_rows],catalog['MOUNTAZ'][duplicate_rows],catalog['WINDAZ'][duplicate_rows],catalog['WINDSPEED'][duplicate_rows],catalog['b10teff'][duplicate_rows],catalog['b10logg'][duplicate_rows],catalog['Magnitude'][duplicate_rows]],names=['WD','DATE-OBS','AIRMASS','SNR','EXPTIME','SEEING','MOUNTAZ','WINDAZ','WINDSPEED','b10teff','b10logg','MAGNITUDE'])


#print newtable

#newtable.write('WINDY_duplicates.txt',format='ascii')
#newtable.write('WINDY_duplicates_latex.txt',format='latex')
