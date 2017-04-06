'''
Combine coords_fits_names.txt with catalog_data.txt by coordinates

'''

from astropy.table import Table
import numpy as np

#Catalog matching tools from www.stsci.edu/~ferguson/software/pygoodsdist/pygoods/
def angsep(ra1deg,dec1deg,ra2deg,dec2deg):
    """ Determine separation in degrees between two celestial objects 
        arguments are RA and Dec in decimal degrees. 
    """
    ra1rad=ra1deg*np.pi/180
    dec1rad=dec1deg*np.pi/180
    ra2rad=ra2deg*np.pi/180
    dec2rad=dec2deg*np.pi/180
    
    # calculate scalar product for determination
    # of angular separation
    
    x=np.cos(ra1rad)*np.cos(dec1rad)*np.cos(ra2rad)*np.cos(dec2rad)
    y=np.sin(ra1rad)*np.cos(dec1rad)*np.sin(ra2rad)*np.cos(dec2rad)
    z=np.sin(dec1rad)*np.sin(dec2rad)
    
    rad=np.arccos(x+y+z) # Sometimes gives warnings when coords match
    
    # use Pythargoras approximation if rad < 1 arcsec
    sep = np.choose( rad<0.000004848 , (
        np.sqrt((np.cos(dec1rad)*(ra1rad-ra2rad))**2+(dec1rad-dec2rad)**2),rad))
        
    # Angular separation
    sep=sep*180/np.pi

    return sep


def matchsorted(ra,dec,ra1,dec1,tol):
    """ Find closest ra,dec within tol to a target in an ra-sorted list of ra,dec.
        Arguments:
          ra - Right Ascension decimal degrees (numpy sorted in ascending order)
          dec - Declination decimal degrees (numpy array)
          ra1 - RA to match (scalar, decimal degrees)
          ra1 - Dec to match (scalar, decimal degrees)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - index of the best match within tol; -1 if no match within tol
          sep - separation (defaults to tol if no match within tol)
    """
    i1 = np.searchsorted(ra,ra1-tol)-1
    i2 = np.searchsorted(ra,ra1+tol)+1
    if i1 < 0:
        i1 = 0
    sep = angsep(ra[i1:i2],dec[i1:i2],ra1,dec1)
    # print "tolerance ",tol
    indices=np.argsort(sep)
    # print sep
    if sep[indices[0]] > tol:
        return -1, tol
    ibest = indices[0] + i1
    return ibest, sep[indices[0]]
    

def matchpos(ra1,dec1,ra2,dec2,tol):
    """ Match two sets of ra,dec within a tolerance.
        Longer catalog should go first
        Arguments:
          ra1 - Right Ascension decimal degrees (numpy array)
          dec1 - Declination decimal degrees (numpy array)
          ra2 - Right Ascension decimal degrees (numpy array)
          dec2 - Declination decimal degrees (numpy array)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - indices of the best matches within tol; -1 if no match within tol
          sep - separations (defaults to tol if no match within tol)
    """
    indices = np.argsort(ra1)
    rasorted = ra1[indices]
    decsorted = dec1[indices]
    ibest = []
    sep = []
    for i in range(len(ra2)):
        j,s = matchsorted(rasorted,decsorted,ra2[i],dec2[i],tol)
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
        sep += [s]
    return np.array(ibest),np.array(sep)


observed_table = Table.read('model_results.txt',format='ascii')

print observed_table.colnames

catalog_data = Table.read('catalog_data.tsv',format='ascii')

print catalog_data.colnames

#use matchpos to find matching indices
RAobserved = np.asarray(observed_table['RA'])
DECobserved = np.asarray(observed_table['DEC'])


RAcat = np.asarray(catalog_data['RA'])
DECcat = np.asarray(catalog_data['DEC'])

tolerance = 0.1 #Tolerance in decimal degrees

ibest, sep = matchpos(RAcat,DECcat,RAobserved,DECobserved,tolerance)
#print ibest
#exit()
#print len(ibest)
#print len(sep)
#print len(RAobserved), len(RAcat)

#print ibest
new_order= np.sort(ibest)

#for x in np.range(len(RAobserved)):
catalog_RA =  np.array(catalog_data['RA'][ibest])
catalog_DEC =  np.array(catalog_data['DEC'][ibest])
catalog_WDname =  np.array(catalog_data['WD Name'][ibest])
catalog_othername = np.array(catalog_data['Other Name'][ibest])
catalog_dav =  np.array(catalog_data['DAV/NOV'][ibest])
catalog_k2 =  np.array(catalog_data['K2'][ibest])
catalog_outburst =  np.array(catalog_data['Outburster'][ibest])
catalog_mag =  np.array(catalog_data['Magnitude'][ibest])
catalog_novlim =  np.array(catalog_data['NOV Limit'][ibest])
catalog_p = np.array(catalog_data['<P>'][ibest])

#print catalog_WDname[0:5],catalog_RA[0:5],catalog_DEC[0:5]


#Take those matching indices and combine the tables

info = Table([observed_table['WD'],observed_table['FILENAME'],catalog_WDname,catalog_othername,observed_table['DATE-OBS'],observed_table['DATE-OBS2'],observed_table['RA'],observed_table['DEC'],observed_table['Airmass'],observed_table['SNR'],observed_table['EXPTIME'],observed_table['Seeing'],observed_table['a10teff'],observed_table['a10tefferr'],observed_table['a10logg'],observed_table['a10loggerr'],observed_table['b10teff'],observed_table['b10tefferr'],observed_table['b10logg'],observed_table['b10loggerr'],observed_table['g10teff'],observed_table['g10tefferr'],observed_table['g10logg'],observed_table['g10loggerr'],observed_table['b9teff'],observed_table['b9tefferr'],observed_table['b9logg'],observed_table['b9loggerr'],observed_table['b8teff'],observed_table['b8tefferr'],observed_table['b8logg'],observed_table['b8loggerr'],observed_table['ateff'],observed_table['atefferr'],observed_table['alogg'],observed_table['aloggerr'],observed_table['bteff'],observed_table['btefferr'],observed_table['blogg'],observed_table['bloggerr'],observed_table['gteff'],observed_table['gtefferr'],observed_table['glogg'],observed_table['gloggerr'],observed_table['dteff'],observed_table['dtefferr'],observed_table['dlogg'],observed_table['dloggerr'],observed_table['eteff'],observed_table['etefferr'],observed_table['elogg'],observed_table['eloggerr'],observed_table['H8teff'],observed_table['H8tefferr'],observed_table['H8logg'],observed_table['H8loggerr'],observed_table['H9teff'],observed_table['H9tefferr'],observed_table['H9logg'],observed_table['H9loggerr'],observed_table['H10teff'],observed_table['H10tefferr'],observed_table['H10logg'],observed_table['H10loggerr'],catalog_RA,catalog_DEC,catalog_dav,catalog_k2,catalog_outburst,catalog_mag,catalog_novlim,catalog_p],names=['WD','FILENAME','WDName','OtherName','DATE-OBS','DATE-OBS2','RA','DEC','Airmass','SNR','EXPTIME','Seeing','a10teff','a10tefferr','a10logg','a10loggerr','b10teff','b10tefferr','b10logg','b10loggerr','g10teff','g10tefferr','g10logg','g10loggerr','b9teff','b9tefferr','b9logg','b9loggerr','b8teff','b8tefferr','b8logg','b8loggerr','ateff','atefferr','alogg','aloggerr','bteff','btefferr','blogg','bloggerr','gteff','gtefferr','glogg','gloggerr','dteff','dtefferr','dlogg','dloggerr','eteff','etefferr','elogg','eloggerr','H8teff','H8tefferr','H8logg','H8loggerr','H9teff','H9tefferr','H9logg','H9loggerr','H10teff','H10tefferr','H10logg','H10loggerr','RA2','DEC2','DAV','K2','outburst','Magnitude','NOV_Lim','P'])

#Sort the table by RA
info.sort('RA')

#Add column for duplicates
duplicates = np.zeros(len(info['RA']))
for x in range(len(info['RA'])):
    if x == (len(info['RA'])-1):
        if np.abs(info['RA'][x] - info['RA'][x-1]) < 0.05 and np.abs(info['DEC'][x] - info['DEC'][x-1]) < 0.05:
            duplicates[x] = 1.
    elif np.abs(info['RA'][x] - info['RA'][x+1]) < 0.05 and np.abs(info['DEC'][x] - info['DEC'][x+1]) < 0.05:
        #print  info['WD'][x], info['RA'][x], info['DEC'][x],catalog['WD'][x+1], info['RA'][x+1], info['DEC'][x+1]
        duplicates[x] = 1.
    elif np.abs(info['RA'][x] - info['RA'][x-1]) < 0.05 and np.abs(info['DEC'][x] - info['DEC'][x-1]) < 0.05:
        duplicates[x] = 1.

    #print info['WD'][x], duplicates[x]

info['Duplicates'] = duplicates


print info.colnames
info.write('catalog_model.txt',format='ascii')
