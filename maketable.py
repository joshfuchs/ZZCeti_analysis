'''
Creates a single table for all the fitted white dwarfs.




'''






import os
from glob import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits
import analysis_tools as at


#===============================================
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

#===============================================
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
    
#===============================================
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




#===============================================
def decimal_dec(hdu_str):
    # Read header strings in "hh:mm:ss" or "dd:mm:ss" fromat 
    # and outputs the value as a decimal. 
    val_list = [float(n) for n in hdu_str.split(':')]
    if str(val_list[0])[0] == '-':
        sng = -1
        val_list[0] = sng*val_list[0]
    else:
        sng = 1
    val_deci =  sng*(val_list[0]+((val_list[1]+(val_list[2]/60.0))/60.0))
    return val_deci


#===============================================
def decimal_ra(hdu_str):
    # Read header strings in "hh:mm:ss" or "dd:mm:ss" fromat 
    # and outputs the value as a decimal. 
    val_list = [float(n) for n in hdu_str.split(':')]
    if val_list[0] < 0 :
        sng = -1
        val_list[0] = sng*val_list[0]
    else:
        sng = 1
    val_deci =  15*sng*(val_list[0]+((val_list[1]+(val_list[2]/60.0))/60.0))

    return val_deci

#===============================================
def read_files(a,b):
    alphafile = a + 'alpha' + b
    betafile = a + 'beta' + b
    gammafile = a + 'gamma' + b
    deltafile = a + 'delta' + b
    epsilonfile = a + 'epsilon' + b
    H8file = a + 'H8' + b
    H9file = a + 'H9' + b
    H10file = a + 'H10' + b

    try:
        alphachi = np.genfromtxt(alphafile,dtype='d')
    except:
        print 'No H-alpha file \n'
        pass
    betachi = np.genfromtxt(betafile,dtype='d')
    gammachi = np.genfromtxt(gammafile,dtype='d')
    deltachi = np.genfromtxt(deltafile,dtype='d')
    epsilonchi = np.genfromtxt(epsilonfile,dtype='d')
    H8chi = np.genfromtxt(H8file,dtype='d')
    H9chi = np.genfromtxt(H9file,dtype='d')
    H10chi = np.genfromtxt(H10file,dtype='d')


#===============================================
#===============================================

#Set file type (model_short or master)
filetype = 'master'

#Set up grid. This is saved in the header of the chi*txt file
bottomt = 10000.
topt = 15000.
stept = 10.
teff = np.linspace(bottomt,topt,(topt-bottomt)/stept+1.,endpoint=True)

bottomg = 7.0 
topg = 9.5
stepg = 0.05
logg = np.linspace(bottomg,topg,(topg-bottomg)/stepg+1.,endpoint=True)



#Get list of directories in current directory. Only keep if it is a date by looking for ./2
directories = [x[0] for x in os.walk('./') if x[0][0:3]=='./2']

filename = []
direc = [] 
a10teffs= []
a10tefferrs= []
a10loggs= []
a10loggerrs= []
b10teffs= []
b10tefferrs= []
b10loggs= []
b10loggerrs= []
g10teffs= []
g10tefferrs= []
g10loggs= []
g10loggerrs= []
b9teffs= []
b9tefferrs= []
b9loggs= []
b9loggerrs= []
b8teffs= []
b8tefferrs= []
b8loggs= []
b8loggerrs= []
ateffs= []
atefferrs= []
aloggs= []
aloggerrs= []
bteffs= []
btefferrs= []
bloggs= []
bloggerrs= []
gteffs= []
gtefferrs= []
gloggs= []
gloggerrs= []
dteffs= []
dtefferrs= []
dloggs= []
dloggerrs= []
eteffs= []
etefferrs= []
eloggs= []
eloggerrs= []
H8teffs= []
H8tefferrs= []
H8loggs= []
H8loggerrs= []
H9teffs= []
H9tefferrs= []
H9loggs= []
H9loggerrs= []
H10teffs= []
H10tefferrs= []
H10loggs= []
H10loggerrs= []




for xdir in directories:
    os.chdir(xdir)
    print xdir

    file_search = 'chi*' + filetype + '*beta*txt'
    file_list = sorted(glob(file_search))
    print file_list
    for new_file in file_list:
        print new_file
        first_part = new_file[0:new_file.find(filetype)] + filetype + '_'
        second_part = new_file[-15:]
        if second_part[0] == 'a':
            second_part = second_part[1:]
        #print first_part, second_part
        #read_files(first_part,second_part)

        alphafile = first_part + 'alpha' + second_part
        betafile = first_part + 'beta' + second_part
        gammafile = first_part + 'gamma' + second_part
        deltafile = first_part + 'delta' + second_part
        epsilonfile = first_part + 'epsilon' + second_part
        H8file = first_part + 'H8' + second_part
        H9file = first_part + 'H9' + second_part
        H10file = first_part + 'H10' + second_part
        
        #Read in grid size and spacing from header of beta file
        with open(betafile,'r') as f:
            first_line = f.readline()
        try:
            bottomg,stepg,topg,bottomt,stept,topt,numpoints = [float(x) for x in first_line[2:].split(",")]
        except:
            bottomg,stepg,topg,bottomt,stept,topt = [float(x) for x in first_line[2:].split(",")]
        teff = np.linspace(bottomt,topt,(topt-bottomt)/stept+1.,endpoint=True)
        logg = np.linspace(bottomg,topg,(topg-bottomg)/stepg+1.,endpoint=True)
        teffgrid, logggrid = np.meshgrid(teff,logg)


        try:
            alphachi = np.genfromtxt(alphafile,dtype='d')
        except:
            print 'No H-alpha file \n'
            alphachi = None
            pass
        betachi = np.genfromtxt(betafile,dtype='d')
        gammachi = np.genfromtxt(gammafile,dtype='d')
        deltachi = np.genfromtxt(deltafile,dtype='d')
        epsilonchi = np.genfromtxt(epsilonfile,dtype='d')
        H8chi = np.genfromtxt(H8file,dtype='d')
        H9chi = np.genfromtxt(H9file,dtype='d')
        H10chi = np.genfromtxt(H10file,dtype='d')


        #Combine different surfaces 
        try:
            a10chi = alphachi + betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
        except:
            print 'No a10chi surface'
            a10chi = None
            pass
        b10chi = betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
        b9chi = betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi
        b8chi = betachi + gammachi + deltachi + epsilonchi + H8chi
        g10chi = gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
        
       #Fit and find all solutions
        try:
            #print '\nAlpha:'
            ateff,atefferr,alogg,aloggerr = at.find_solution(alphachi,logg,teff)
        except:
            print '\nNo Alpha:'
            ateff,atefferr,alogg,aloggerr = 0.,0.,0.,0.
            pass
        
        try:
            #print '\nAlpha10:'
            a10teff,a10tefferr,a10logg,a10loggerr = at.find_solution(a10chi,logg,teff)
        except:
            print '\nNo Alpha10:'
            a10teff,a10tefferr,a10logg,a10loggerr =  0.,0.,0.,0.
            pass

        bteff,btefferr,blogg,bloggerr = at.find_solution(betachi,logg,teff)
        gteff,gtefferr,glogg,gloggerr = at.find_solution(gammachi,logg,teff)
        dteff,dtefferr,dlogg,dloggerr = at.find_solution(deltachi,logg,teff)
        eteff,etefferr,elogg,eloggerr = at.find_solution(epsilonchi,logg,teff)
        H8teff,H8tefferr,H8logg,H8loggerr = at.find_solution(H8chi,logg,teff)
        H9teff,H9tefferr,H9logg,H9loggerr = at.find_solution(H9chi,logg,teff)
        H10teff,H10tefferr,H10logg,H10loggerr = at.find_solution(H10chi,logg,teff)
 
        b10teff,b10tefferr,b10logg,b10loggerr = at.find_solution(b10chi,logg,teff)
        b9teff,b9tefferr,b9logg,b9loggerr = at.find_solution(b9chi,logg,teff)
        b8teff,b8tefferr,b8logg,b8loggerr = at.find_solution(b8chi,logg,teff)
        g10teff,g10tefferr,g10logg,g10loggerr = at.find_solution(g10chi,logg,teff)
        

        filename.append(new_file[9:new_file.find('_930')])
        direc.append(xdir)
        a10teffs.append(str(a10teff))
        a10tefferrs.append(str(a10tefferr))
        a10loggs.append(str(a10logg))
        a10loggerrs.append(str(a10loggerr))
        b10teffs.append(str(b10teff))
        b10tefferrs.append(str(b10tefferr))
        b10loggs.append(str(b10logg))
        b10loggerrs.append(str(b10loggerr))
        g10teffs.append(str(g10teff))
        g10tefferrs.append(str(g10tefferr))
        g10loggs.append(str(g10logg))
        g10loggerrs.append(str(g10loggerr))
        b9teffs.append(str(b9teff))
        b9tefferrs.append(str(b9tefferr))
        b9loggs.append(str(b9logg))
        b9loggerrs.append(str(b9loggerr))
        b8teffs.append(str(b8teff))
        b8tefferrs.append(str(b8tefferr))
        b8loggs.append(str(b8logg))
        b8loggerrs.append(str(b8loggerr))
        ateffs.append(str(ateff))
        atefferrs.append(str(atefferr))
        aloggs.append(str(alogg))
        aloggerrs.append(str(aloggerr))
        bteffs.append(str(bteff))
        btefferrs.append(str(btefferr))
        bloggs.append(str(blogg))
        bloggerrs.append(str(bloggerr))
        gteffs.append(str(gteff))
        gtefferrs.append(str(gtefferr))
        gloggs.append(str(glogg))
        gloggerrs.append(str(gloggerr))
        dteffs.append(str(dteff))
        dtefferrs.append(str(dtefferr))
        dloggs.append(str(dlogg))
        dloggerrs.append(str(dloggerr))
        eteffs.append(str(eteff))
        etefferrs.append(str(etefferr))
        eloggs.append(str(elogg))
        eloggerrs.append(str(eloggerr))
        H8teffs.append(str(H8teff))
        H8tefferrs.append(str(H8tefferr))
        H8loggs.append(str(H8logg))
        H8loggerrs.append(str(H8loggerr))
        H9teffs.append(str(H9teff))
        H9tefferrs.append(str(H9tefferr))
        H9loggs.append(str(H9logg))
        H9loggerrs.append(str(H9loggerr))
        H10teffs.append(str(H10teff))
        H10tefferrs.append(str(H10tefferr))
        H10loggs.append(str(H10logg))
        H10loggerrs.append(str(H10loggerr) )


    print '\n'
    os.chdir('../')


fit_table = Table([filename,direc,a10teffs,a10tefferrs,a10loggs,a10loggerrs,b10teffs,b10tefferrs,b10loggs,b10loggerrs,g10teffs,g10tefferrs,g10loggs,g10loggerrs,b9teffs,b9tefferrs,b9loggs,b9loggerrs,b8teffs,b8tefferrs,b8loggs,b8loggerrs,ateffs,atefferrs,aloggs,aloggerrs,bteffs,btefferrs,bloggs,bloggerrs,gteffs,gtefferrs,gloggs,gloggerrs,dteffs,dtefferrs,dloggs,dloggerrs,eteffs,etefferrs,eloggs,eloggerrs,H8teffs,H8tefferrs,H8loggs,H8loggerrs,H9teffs,H9tefferrs,H9loggs,H9loggerrs,H10teffs,H10tefferrs,H10loggs,H10loggerrs],names=['Filename','DATOBS','a10teff','a10tefferr','a10logg','a10loggerr','b10teff','b10tefferr','b10logg','b10loggerr','g10teff','g10tefferr','g10logg','g10loggerr','b9teff','b9tefferr','b9logg','b9loggerr','b8teff','b8tefferr','b8logg','b8loggerr','ateff','atefferr','alogg','aloggerr','bteff','btefferr','blogg','bloggerr','gteff','gtefferr','glogg','gloggerr','dteff','dtefferr','dlogg','dloggerr','eteff','etefferr','elogg','eloggerr','H8teff','H8tefferr','H8logg','H8loggerr','H9teff','H9tefferr','H9logg','H9loggerr','H10teff','H10tefferr','H10logg','H10loggerr'])
#fit_table.write('all_teff_logg_flux.txt',format='ascii')


##########################
#Search through blue fits headers and get spectral details
##########################


spec_files = []
date_obs = []
RAs = []
DECs = []
airmass = []
snr = []
nimages = []
exptime = []
seeing = []
mount_az = []
wind_az = []
wind_speed = []
for xdir in directories:
    os.chdir(xdir)
    #print xdir
    specfile_search = '*blue*' + filetype + '*fits'
    blue_files = sorted(glob(specfile_search))
    #print blue_files
    for spectrum in blue_files:
        hdulist = fits.open(spectrum)
        spec_files.append(spectrum)
        RAs.append(hdulist[0].header['RA'])
        DECs.append(hdulist[0].header['DEC'])
        airmass.append(hdulist[0].header['AIRMASS'])
        snr.append(hdulist[0].header['SNR'])
        nimages.append(hdulist[0].header['NCOMBINE'])
        exptime.append(hdulist[0].header['EXPTIME'])
        seeing.append(hdulist[0].header['SPECFWHM'])
        mount_az.append(hdulist[0].header['MOUNT_AZ'])
        wind_az.append(hdulist[0].header['ENVDIR'])
        wind_speed.append(hdulist[0].header['ENVWIN'])
        date_obs.append(xdir)
    


    os.chdir('../')

exp = np.asarray(nimages)*np.asarray(exptime)

spec_details_table = Table([spec_files,date_obs,RAs,DECs,airmass,snr,exp,seeing,mount_az,wind_az,wind_speed],names=['File Name','DATE-OBS','RA','DEC','Airmass','SNR','EXPTIME','SEEING','MOUNTAZ','WINDAZ','WINDSPEED'])
#spec_details_table.write('spec_details.txt',format = 'ascii')


##########################
#Combine spec_details_table and fits_table into single file
##########################


newwd = []
newfile = []
newdate = []
newdate2 = []
newra = []
newdec = []
newairmass = []
newsnr = []
newexp = []
newseeing = []
newmountaz  = []
newwindaz = []
newwindspeed = []
a10teffs= []
a10tefferrs= []
a10loggs= []
a10loggerrs= []
b10teffs= []
b10tefferrs= []
b10loggs= []
b10loggerrs= []
g10teffs= []
g10tefferrs= []
g10loggs= []
g10loggerrs= []
b9teffs= []
b9tefferrs= []
b9loggs= []
b9loggerrs= []
b8teffs= []
b8tefferrs= []
b8loggs= []
b8loggerrs= []
ateffs= []
atefferrs= []
aloggs= []
aloggerrs= []
bteffs= []
btefferrs= []
bloggs= []
bloggerrs= []
gteffs= []
gtefferrs= []
gloggs= []
gloggerrs= []
dteffs= []
dtefferrs= []
dloggs= []
dloggerrs= []
eteffs= []
etefferrs= []
eloggs= []
eloggerrs= []
H8teffs= []
H8tefferrs= []
H8loggs= []
H8loggerrs= []
H9teffs= []
H9tefferrs= []
H9loggs= []
H9loggerrs= []
H10teffs= []
H10tefferrs= []
H10loggs= []
H10loggerrs= []
for x in range(0,len(spec_details_table)):
    for y in range(0,len(fit_table)):
        if fit_table['Filename'][y] in spec_details_table['File Name'][x] and spec_details_table['DATE-OBS'][x] == fit_table['DATOBS'][y]:
            #print coord_table['File Name'][x], coord_table['DATE-OBS'][x], fit_table['WD'][y], fit_table['DATE-OBS'][y],fit_table['Teff'][y]
            newfile.append(spec_details_table['File Name'][x])
            newdate2.append(spec_details_table['DATE-OBS'][x])
            newra.append(spec_details_table['RA'][x])
            newdec.append(spec_details_table['DEC'][x])
            newairmass.append(spec_details_table['Airmass'][x])
            newsnr.append(spec_details_table['SNR'][x])
            newexp.append(spec_details_table['EXPTIME'][x])
            newseeing.append(spec_details_table['SEEING'][x])
            newmountaz.append(spec_details_table['MOUNTAZ'][x])
            newwindaz.append(spec_details_table['WINDAZ'][x])
            newwindspeed.append(spec_details_table['WINDSPEED'][x])

            newwd.append(fit_table['Filename'][y])
            newdate.append(fit_table['DATOBS'][y])
            a10teffs.append(fit_table['a10teff'][y])
            a10tefferrs.append(fit_table['a10tefferr'][y])
            a10loggs.append(fit_table['a10logg'][y])
            a10loggerrs.append(fit_table['a10loggerr'][y])
            b10teffs.append(fit_table['b10teff'][y])
            b10tefferrs.append(fit_table['b10tefferr'][y])
            b10loggs.append(fit_table['b10logg'][y])
            b10loggerrs.append(fit_table['b10loggerr'][y])
            
            g10teffs.append(fit_table['g10teff'][y])
            g10tefferrs.append(fit_table['g10tefferr'][y])
            g10loggs.append(fit_table['g10logg'][y])
            g10loggerrs.append(fit_table['g10loggerr'][y])

            b9teffs.append(fit_table['b9teff'][y])
            b9tefferrs.append(fit_table['b9tefferr'][y])
            b9loggs.append(fit_table['b9logg'][y])
            b9loggerrs.append(fit_table['b9loggerr'][y])
            b8teffs.append(fit_table['b8teff'][y])
            b8tefferrs.append(fit_table['b8tefferr'][y])
            b8loggs.append(fit_table['b8logg'][y])
            b8loggerrs.append(fit_table['b8loggerr'][y])

            ateffs.append(fit_table['ateff'][y])
            atefferrs.append(fit_table['atefferr'][y])
            aloggs.append(fit_table['alogg'][y])
            aloggerrs.append(fit_table['aloggerr'][y])
            bteffs.append(fit_table['bteff'][y])
            btefferrs.append(fit_table['btefferr'][y])
            bloggs.append(fit_table['blogg'][y])
            bloggerrs.append(fit_table['bloggerr'][y])
            gteffs.append(fit_table['gteff'][y])
            gtefferrs.append(fit_table['gtefferr'][y])
            gloggs.append(fit_table['glogg'][y])
            gloggerrs.append(fit_table['gloggerr'][y])
            dteffs.append(fit_table['dteff'][y])
            dtefferrs.append(fit_table['dtefferr'][y])
            dloggs.append(fit_table['dlogg'][y])
            dloggerrs.append(fit_table['dloggerr'][y])
            eteffs.append(fit_table['eteff'][y])
            etefferrs.append(fit_table['etefferr'][y])
            eloggs.append(fit_table['elogg'][y])
            eloggerrs.append(fit_table['eloggerr'][y])
            H8teffs.append(fit_table['H8teff'][y])
            H8tefferrs.append(fit_table['H8tefferr'][y])
            H8loggs.append(fit_table['H8logg'][y])
            H8loggerrs.append(fit_table['H8loggerr'][y])
            H9teffs.append(fit_table['H9teff'][y])
            H9tefferrs.append(fit_table['H9tefferr'][y])
            H9loggs.append(fit_table['H9logg'][y])
            H9loggerrs.append(fit_table['H9loggerr'][y])
            H10teffs.append(fit_table['H10teff'][y])
            H10tefferrs.append(fit_table['H10tefferr'][y])
            H10loggs.append(fit_table['H10logg'][y])
            H10loggerrs.append(fit_table['H10loggerr'][y])
            
            #if fit_table['DATE-OBS'][y] == './2016-07-14' or '.2016-08-06':
            #    print fit_table['WD'][y] ,fit_table['DATE-OBS'][y],coord_table['File Name'][x],coord_table['RA'][x],coord_table['DEC'][x]

#Transform RA and DEC to decimal notation
radec = np.zeros(len(newra))
decdec = np.zeros(len(newdec))
for x in range(0,len(newra)):
    radec[x] = decimal_ra(newra[x])
    decdec[x] = decimal_dec(newdec[x])





observed_table = Table([newwd,newfile,newdate,newdate2,radec,decdec,newairmass,newsnr,newexp,newseeing,newmountaz,newwindaz,newwindspeed,a10teffs,a10tefferrs,a10loggs,a10loggerrs,b10teffs,b10tefferrs,b10loggs,b10loggerrs,g10teffs,g10tefferrs,g10loggs,g10loggerrs,b9teffs,b9tefferrs,b9loggs,b9loggerrs,b8teffs,b8tefferrs,b8loggs,b8loggerrs,ateffs,atefferrs,aloggs,aloggerrs,bteffs,btefferrs,bloggs,bloggerrs,gteffs,gtefferrs,gloggs,gloggerrs,dteffs,dtefferrs,dloggs,dloggerrs,eteffs,etefferrs,eloggs,eloggerrs,H8teffs,H8tefferrs,H8loggs,H8loggerrs,H9teffs,H9tefferrs,H9loggs,H9loggerrs,H10teffs,H10tefferrs,H10loggs,H10loggerrs],names=['WD','FILENAME','DATE-OBS','DATE-OBS2','RA','DEC','Airmass','SNR','EXPTIME','Seeing','MOUNTAZ','WINDAZ','WINDSPEED','a10teff','a10tefferr','a10logg','a10loggerr','b10teff','b10tefferr','b10logg','b10loggerr','g10teff','g10tefferr','g10logg','g10loggerr','b9teff','b9tefferr','b9logg','b9loggerr','b8teff','b8tefferr','b8logg','b8loggerr','ateff','atefferr','alogg','aloggerr','bteff','btefferr','blogg','bloggerr','gteff','gtefferr','glogg','gloggerr','dteff','dtefferr','dlogg','dloggerr','eteff','etefferr','elogg','eloggerr','H8teff','H8tefferr','H8logg','H8loggerr','H9teff','H9tefferr','H9logg','H9loggerr','H10teff','H10tefferr','H10logg','H10loggerr'])

'''
observed_table.sort('RA')
observed_table.write('test_flux_results.txt',format='ascii')
'''


#################
#Match observed_table with details from catalog_data.tsv
#################

catalog_data = Table.read('catalog_data.tsv',format='ascii')

#use matchpos to find matching indices
RAobserved = np.asarray(observed_table['RA'])
DECobserved = np.asarray(observed_table['DEC'])


RAcat = np.asarray(catalog_data['RA'])
DECcat = np.asarray(catalog_data['DEC'])

tolerance = 0.1 #Tolerance in decimal degrees

ibest, sep = matchpos(RAcat,DECcat,RAobserved,DECobserved,tolerance)
new_order= np.sort(ibest)

#Sort the catalog data correctly 
catalog_name = np.array(catalog_data['Name'][ibest])
catalog_WDname =  np.array(catalog_data['WD Name'][ibest])
catalog_othername = np.array(catalog_data['Other Name'][ibest])
catalog_RA =  np.array(catalog_data['RA'][ibest])
catalog_DEC =  np.array(catalog_data['DEC'][ibest])
catalog_dav =  np.array(catalog_data['DAV/NOV'][ibest])
catalog_k2 =  np.array(catalog_data['K2'][ibest])
catalog_outburst =  np.array(catalog_data['Outburster'][ibest])
catalog_mag =  np.array(catalog_data['Magnitude'][ibest])
catalog_novlim =  np.array(catalog_data['NOV Limit'][ibest])
catalog_p = np.array(catalog_data['<P>'][ibest])
catalog_wmp = np.array(catalog_data['WMP'][ibest])
catalog_shortp = np.array(catalog_data['Shortest P'][ibest])
catalog_longp = np.array(catalog_data['Longest P'][ibest])






#Take those matching indices and combine the tables

info = Table([observed_table['WD'],observed_table['FILENAME'],catalog_WDname,catalog_othername,observed_table['DATE-OBS'],observed_table['DATE-OBS2'],observed_table['RA'],observed_table['DEC'],observed_table['Airmass'],observed_table['SNR'],observed_table['EXPTIME'],observed_table['Seeing'],observed_table['MOUNTAZ'],observed_table['WINDAZ'],observed_table['WINDSPEED'],observed_table['a10teff'],observed_table['a10tefferr'],observed_table['a10logg'],observed_table['a10loggerr'],observed_table['b10teff'],observed_table['b10tefferr'],observed_table['b10logg'],observed_table['b10loggerr'],observed_table['g10teff'],observed_table['g10tefferr'],observed_table['g10logg'],observed_table['g10loggerr'],observed_table['b9teff'],observed_table['b9tefferr'],observed_table['b9logg'],observed_table['b9loggerr'],observed_table['b8teff'],observed_table['b8tefferr'],observed_table['b8logg'],observed_table['b8loggerr'],observed_table['ateff'],observed_table['atefferr'],observed_table['alogg'],observed_table['aloggerr'],observed_table['bteff'],observed_table['btefferr'],observed_table['blogg'],observed_table['bloggerr'],observed_table['gteff'],observed_table['gtefferr'],observed_table['glogg'],observed_table['gloggerr'],observed_table['dteff'],observed_table['dtefferr'],observed_table['dlogg'],observed_table['dloggerr'],observed_table['eteff'],observed_table['etefferr'],observed_table['elogg'],observed_table['eloggerr'],observed_table['H8teff'],observed_table['H8tefferr'],observed_table['H8logg'],observed_table['H8loggerr'],observed_table['H9teff'],observed_table['H9tefferr'],observed_table['H9logg'],observed_table['H9loggerr'],observed_table['H10teff'],observed_table['H10tefferr'],observed_table['H10logg'],observed_table['H10loggerr'],catalog_RA,catalog_DEC,catalog_dav,catalog_k2,catalog_outburst,catalog_mag,catalog_novlim,catalog_p,catalog_wmp,catalog_shortp,catalog_longp],names=['WD','FILENAME','WDName','OtherName','DATE-OBS','DATE-OBS2','RA','DEC','Airmass','SNR','EXPTIME','Seeing','MOUNTAZ','WINDAZ','WINDSPEED','a10teff','a10tefferr','a10logg','a10loggerr','b10teff','b10tefferr','b10logg','b10loggerr','g10teff','g10tefferr','g10logg','g10loggerr','b9teff','b9tefferr','b9logg','b9loggerr','b8teff','b8tefferr','b8logg','b8loggerr','ateff','atefferr','alogg','aloggerr','bteff','btefferr','blogg','bloggerr','gteff','gtefferr','glogg','gloggerr','dteff','dtefferr','dlogg','dloggerr','eteff','etefferr','elogg','eloggerr','H8teff','H8tefferr','H8logg','H8loggerr','H9teff','H9tefferr','H9logg','H9loggerr','H10teff','H10tefferr','H10logg','H10loggerr','RA2','DEC2','DAV','K2','outburst','Magnitude','NOV_Lim','<P>','WMP','ShortestP','LongestP'])

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


#print info.colnames
newname = 'catalog_' + filetype + '.txt'
info.write(newname,format='ascii')
