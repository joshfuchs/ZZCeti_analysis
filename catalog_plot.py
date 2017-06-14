import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table



catalog = Table.read('catalog_flux_clean_nodup.txt',format='ascii')
#print catalog.colnames

#print catalog
#print len(catalog)
#exit()

'''
for x in range(len(catalog['RA'])):
    if np.abs(catalog['RA'][x] - catalog['RA2'][x]) > 0.01:
        print 'RA: ', catalog['WD'][x],catalog['WD2'][x],catalog['RA'][x], catalog['RA2'][x], catalog['DATE-OBS'][x],catalog['DATE-OBS2'][x]
    if np.abs(catalog['DEC'][x] - catalog['DEC2'][x]) > 0.01:
        print 'DEC: ', catalog['WD'][x],catalog['WD2'][x],catalog['DEC'][x], catalog['DEC2'][x], catalog['DATE-OBS'][x],catalog['DATE-OBS2'][x]
'''

#catalog.show_in_browser(jsviewer=True)

#print np.median(catalog['SNR'])
#print 0.3*np.median(catalog['Seeing'])
#print 0.3*np.mean(catalog['Seeing'])
#print 0.3*np.std(catalog['Seeing'])
#print 0.3*np.min(catalog['Seeing'])
#print 0.3*np.max(catalog['Seeing'])

#np.savetxt('magnitudes.txt',np.transpose([catalog['Magnitude']]))

#for x in range(len(catalog['b10teff'])):
#    if (catalog['SNR'][x] > 400.):
#        print catalog['WD'][x], catalog['DATE-OBS'][x]


#exit()
plt.clf()

#plt.plot(catalog['b10teff'],catalog['b10logg'],'bo')
#plt.plot(catalog['RA'],catalog['RA2'],'bo')

plt.xlim(15500,9500) #15500, 9500
plt.ylim(9.75,7.0) #9.75, 7.0

teffchoice = 'b10teff'
loggchoice = 'b10logg'

for x in range(len(catalog['b10teff'])):
    #print catalog['WD'][x]
    '''
    if catalog['Magnitude'][x] < 17.0:
        plt.plot(catalog['b10teff'][x],catalog['b10logg'][x],'bo')
    if catalog['Magnitude'][x] >= 17.0:
        plt.plot(catalog['b10teff'][x],catalog['b10logg'][x],'ro')
    '''
    
    if catalog['DAV'][x] == 1:
       plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'bo',markersize=7.0)
       #if catalog['K2'][x] == 1:
       #    plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'go',markersize=7.0)
       #plt.annotate(catalog['WD'][x],(catalog[teffchoice][x],catalog[loggchoice][x]))
        #if catalog['Teff'][x] > 13000:
        #    print catalog['WD'][x],catalog['Teff'][x],catalog['logg'][x],catalog['DATE-OBS'][x]
    if catalog['DAV'][x] == 0:
        plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'r^',markersize=7.0)
        #if catalog['K2'][x] == 1:
        #   plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'g^',markersize=7.0)
        #plt.annotate(catalog['WD'][x],(catalog[teffchoice][x],catalog[loggchoice][x]))
    
    #    #print catalog['WD'][x],catalog['Teff'][x],catalog['logg'][x]
    #if catalog['SNR'][x] > 300:
    #    plt.plot(catalog['b10teff'][x],catalog['b10logg'][x],'g^',markersize=10.0)
    #if catalog['K2'][x] == 1:
    #    plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'g^',markersize=10.0)
    '''
    if catalog['Duplicates'][x] == 1:
        a10teffdup = []
        a10loggdup = []
        b10teffdup = []
        b10loggdup = []
        g10teffdup = []
        g10loggdup = []
        b9teffdup = []
        b9loggdup = []
        b8teffdup = []
        b8loggdup = []
        for m in range(10):
            if np.abs(catalog['DEC'][x]-catalog['DEC'][x+m]) < 0.1:
                plt.plot(catalog[teffchoice][x+m],catalog[loggchoice][x+m],'bo')
                a10teffdup.append(catalog['a10teff'][x+m])
                a10loggdup.append(catalog['a10logg'][x+m])
                b10teffdup.append(catalog['b10teff'][x+m])
                b10loggdup.append(catalog['b10logg'][x+m])
                g10teffdup.append(catalog['g10teff'][x+m])
                g10loggdup.append(catalog['g10logg'][x+m])
                b9teffdup.append(catalog['b9teff'][x+m])
                b9loggdup.append(catalog['b9logg'][x+m])
                b8teffdup.append(catalog['b8teff'][x+m])
                b8loggdup.append(catalog['b8logg'][x+m])


                print catalog['WD'][x+m], catalog['DATE-OBS'][x+m], catalog['Duplicates'][x+m]
        print 'Teff range: ', np.ptp(np.asarray(b10teffdup))
        print 'Teff rstd: ', np.std(np.asarray(b10teffdup))
        print 'logg range: ', np.ptp(np.asarray(b10loggdup))
        print 'logg std: ', np.std(np.asarray(b10loggdup))
        #print 'A-10'
        #print np.mean(np.asarray(a10teffdup)), np.mean(np.asarray(a10loggdup))
        #print 'B-10'
        #print np.mean(np.asarray(b10teffdup)), np.mean(np.asarray(b10loggdup))
        #print 'G-10'
        #print np.mean(np.asarray(g10teffdup)), np.mean(np.asarray(g10loggdup))
        #print 'B-9'
        #print np.mean(np.asarray(b9teffdup)), np.mean(np.asarray(b9loggdup))
        #print 'B-8'
        #print np.mean(np.asarray(b8teffdup)), np.mean(np.asarray(b8loggdup))

        plt.show()
        print ''
    '''
                    
    
    #if catalog['outburst'][x] == 1:
    #    plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'g^',markersize=10.0)
    #    #print catalog['WD'][x]
    if ('1425-811' in catalog['WD'][x]):
        plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'g^',markersize=10.0)
    if ('1422p095' in catalog['WD'][x]):
        plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'c^',markersize=10.0)
    if ('1116' in catalog['WD'][x]) :
        plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'m^',markersize=10.0)
    if ('-153' in catalog['WD'][x]):
        plt.plot(catalog[teffchoice][x],catalog[loggchoice][x],'k^',markersize=10.0)
    
    #for x in range(0,len(catalog['WD'])):
    '''
    if '0940' in catalog['WD'][x]:
        print catalog['WD'][x]
        plt.clf()
        plt.xlim(15500,9500) #15500, 9500
        plt.ylim(9.75,7.0) #9.75, 7.0
        plt.annotate('alpha',(catalog['ateff'][x],catalog['alogg'][x]))
        plt.annotate('beta',(catalog['bteff'][x],catalog['blogg'][x]))
        plt.annotate('gamma',(catalog['gteff'][x],catalog['glogg'][x]))
        plt.annotate('delta',(catalog['dteff'][x],catalog['dlogg'][x]))
        plt.annotate('epsil',(catalog['eteff'][x],catalog['elogg'][x]))
        plt.annotate('8',(catalog['H8teff'][x],catalog['H8logg'][x]))
        plt.annotate('9',(catalog['H9teff'][x],catalog['H9logg'][x]))
        plt.annotate('10',(catalog['H10teff'][x],catalog['H10logg'][x]))
        plt.annotate('A10',(catalog['a10teff'][x],catalog['a10logg'][x]))
        plt.annotate('B10',(catalog['b10teff'][x],catalog['b10logg'][x]))
        plt.show()
    '''
    

    
#Plot Gianninas 2011 instability strip
#plt.plot([14000.,11600.],[8.82,7.36],'r')
#plt.plot([11800,10600],[8.88,7.35],'r')
#plt.plot(0.,0.,'bo',label='DAV')
#plt.plot(0.,0.,'r^',label='NOV')
#plt.plot(0.,0.,'go',label='K2 DAV')
#plt.plot(0.,0.,'g^',label='K2 NOV')
#plt.legend(loc=3)
plt.xlabel('Temperature (K)' )
plt.ylabel('log (g)')
plt.show()


