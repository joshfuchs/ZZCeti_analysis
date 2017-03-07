'''
Written by JT Fuchs, November 2016

Fit chi-square surfaces 

'''

import os
from glob import glob
import numpy as np


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


def find_solution(combined,logg,teff):
    #======================================
    # Fit a parabola in the logg direction, then take the minimum from those fits, and fit a parabola in the Teff direction
    # Find minimum point, move + and - 10 and 10  K to fit 5 parabolas in logg, take those centers, fit parabola
    combinedindex = np.unravel_index(combined.argmin(),combined.shape)
    combinedlogg, combinedteff = logg[combinedindex[0]], teff[combinedindex[1]]
    rangeg = 3. #Number of grid points in logg space around lowest value to pick
    ranget = 3. #Number of grid points in Teff space around lowest value to pick
    
    #pick out region of grid with spacing of rangeg and ranget around the minimum
    if combinedindex[0]-rangeg < 0:
        ##loggsmall = logg[0:combinedindex[0]+rangeg+1]
        loggsmall = logg[0:2*rangeg+1]
        logglow = 0
        ##logghigh = combinedindex[0]+rangeg+1
        logghigh = 2*rangeg+1
    elif combinedindex[0]+rangeg >= len(logg):
        ##loggsmall = logg[combinedindex[0]-rangeg:-1]
        loggsmall = logg[-2*rangeg-1:]
        ##logglow = combinedindex[0]-rangeg
        logglow = -2*rangeg-1
        logghigh = -1
    else:
        loggsmall = logg[combinedindex[0]-rangeg:combinedindex[0]+rangeg+1]
        logglow = combinedindex[0]-rangeg
        logghigh = combinedindex[0]+rangeg+1
    if combinedindex[1]-ranget < 0:
        ##teffsmall = teff[0:combinedindex[1]+ranget+1]
        teffsmall = teff[0:2*ranget+1]
        tefflow = 0
        ##teffhigh = combinedindex[1]+ranget+1
        teffhigh = 2*ranget+1
    elif combinedindex[1]+ranget >= len(teff):
        ##teffsmall = teff[combinedindex[1]-ranget:-1]
        teffsmall = teff[-2*ranget-1:]
        ##tefflow = combinedindex[1]-ranget
        tefflow = -2*ranget-1
        teffhigh = -1
    else:
        teffsmall = teff[combinedindex[1]-ranget:combinedindex[1]+ranget+1]
        tefflow = combinedindex[1]-ranget
        teffhigh = combinedindex[1]+ranget+1
    
    
    
    #Get the low and high values for each
    lowg, highg = loggsmall[0], loggsmall[-1]
    lowt, hight = teffsmall[0], teffsmall[-1]
    teffsmallgrid, loggsmallgrid = np.meshgrid(teffsmall,loggsmall)
    if (logghigh == -1) and (teffhigh == -1):
        combinedsmall = combined[logglow:,tefflow:]
    elif logghigh == -1:
        combinedsmall = combined[logglow:,tefflow:teffhigh]
    elif teffhigh == -1:
        combinedsmall = combined[logglow:logghigh,tefflow:]
    else:
        combinedsmall = combined[logglow:logghigh,tefflow:teffhigh]

    #Create finer small grid with spacing of 1 K and 0.005 logg
    lenteffgrid = np.round(hight-lowt+1.) #Round to ensure we get the correct number of points. Otherwise, occasionally face a strange int/float issue.
    teffsmallfine = np.linspace(lowt,hight,lenteffgrid,endpoint=True)
    lenlogggrid = np.round((highg-lowg)*1000.+1.)
    loggsmallfine = np.linspace(lowg,highg,lenlogggrid,endpoint=True)
    teffsmallfinegrid, loggsmallfinegrid = np.meshgrid(teffsmallfine,loggsmallfine)

    #Fit a polynomial to different Teff values to find center of logg
    loggval = np.zeros(len(combinedsmall[:,0]))
    chival = np.zeros(len(combinedsmall[:,0]))
    for x in np.arange(len(combinedsmall[:,0])):
        pol = np.polyfit(loggsmall,combinedsmall[:,x],2)
        pc = np.poly1d(pol)
        if x == np.median(np.arange(len(combinedsmall[:,0]))):
            pckeep = np.poly1d(pol)
        loggval[x] = loggsmallfine[pc(loggsmallfine).argmin()]
        chival[x] = pc(loggval[x])
        #plt.clf()
        #plt.plot(loggsmall,combinedsmall[:,x],'b^')
        #plt.plot(loggsmallfine,pc(loggsmallfine))
        #plt.show()

    #Now take these values and fit a polynomial in the Teff direction
    tpol = np.polyfit(teffsmall,chival,2)
    tpp = np.poly1d(tpol)

    bestteff = teffsmallfine[tpp(teffsmallfine).argmin()]
    #print chival, chival.min()
    #print tpp(bestteff)
    #plt.clf()
    #plt.plot(teffsmall,chival,'g^')
    #plt.plot(teffsmallfine,tpp(teffsmallfine),'k')
    #plt.plot(bestteff,tpp(bestteff),'ms')
    #plt.show()

    #Take these solutions and find the errors
    deltateff = tpp(teffsmallfine)-tpp(bestteff)
    lowterr = teffsmallfine[np.where(deltateff < 2.3)].min()-1.
    highterr =  teffsmallfine[np.where(deltateff < 2.3)].max()+1.
    tefferr = ((bestteff-lowterr) + (highterr-bestteff)) / 2.
    
    deltalogg = pckeep(loggsmallfine) - pckeep(loggval.min())
    bestlogg = loggval.min()
    lowloggerr =  loggsmallfine[np.where(deltalogg < 2.3)].min()-0.001
    highloggerr =  loggsmallfine[np.where(deltalogg < 2.3)].max()+0.001
    loggerr = ((bestlogg-lowloggerr) + (highloggerr-bestlogg)) / 2.
    
    print 'Teff = ', bestteff , '+/-' , tefferr
    print 'logg = ', bestlogg,  '+/-' , loggerr
    
    #plt.clf()
    #plt.plot(loggsmallfine,deltalogg)
    #plt.show()
    return bestteff, tefferr, bestlogg, loggerr

#===============================================



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

f = open('all_teff_logg.txt','w')

for xdir in directories:
    os.chdir(xdir)
    print xdir

    file_list = glob('chi*beta*txt')
    #print file_list
    for new_file in file_list:
        print new_file
        first_part = new_file[0:new_file.find('model')] + 'model_'
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
        
        try:
            alphachi = np.genfromtxt(alphafile,dtype='d')
        except:
            #print 'No H-alpha file \n'
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
            a10 = alphachi + betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
        except:
            pass
        b10 = betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
        g10 = gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
        
       #Fit and find all solutions
       try:
           print '\nAlpha:'
           ateff,atefferr,alogg,aloggerr = find_solution(alphachi,logg,teff)
           a10teff,a10tefferr,a10logg,a10loggerr = find_solution(a10chi,logg,teff)
       except:
           ateff,atefferr,alogg,aloggerr = 0.,0.,0.,0.
           a10teff,a10tefferr,a10logg,a10loggerr =  0.,0.,0.,0.
           pass
        

        bteff,btefferr,blogg,bloggerr = find_solution(betachi,logg,teff)
        gteff,gtefferr,glogg,gloggerr = find_solution(gammachi,logg,teff)
        dteff,dtefferr,dlogg,dloggerr = find_solution(deltachi,logg,teff)
        eteff,etefferr,elogg,eloggerr = find_solution(epsilonchi,logg,teff)
        H8teff,H8tefferr,H8logg,H8loggerr = find_solution(H8chi,logg,teff)
        H9teff,H9tefferr,H9logg,H9loggerr = find_solution(H9chi,logg,teff)
        H10teff,H10tefferr,H10logg,H10loggerr = find_solution(H10chi,logg,teff)
 
        b10teff,b10tefferr,b10logg,b10loggerr = find_solution(b10chi,logg,teff)
        g10teff,g10tefferr,g10logg,g10loggerr = find_solution(g10chi,logg,teff)
        



        info = new_file[9:new_file.find('_930')]+ '\t' + xdir + '\t' + str(a10teff) + '\t' + str(a10tefferr) + '\t' + str(a10logg) + '\t' + str(a10loggerr) + '\t' + str(b10teff) + '\t' + str(b10tefferr) + '\t' + str(b10logg) + '\t' + str(b10loggerr) + '\t' + str(g10teff) + '\t' + str(g10tefferr) + '\t' + str(g10logg) + '\t' + str(g10loggerr)  + '\t' + str(ateff) + '\t' + str(atefferr) + '\t' + str(alogg) + '\t' + str(aloggerr) + '\t' + str(bteff) + '\t' + str(btefferr) + '\t' + str(blogg) + '\t' + str(bloggerr) + '\t' + str(gteff) + '\t' + str(gtefferr) + '\t' + str(glogg) + '\t' + str(gloggerr)  + '\t' + str(dteff) + '\t' + str(dtefferr) + '\t' + str(dlogg) + '\t' + str(dloggerr) + '\t' + str(eteff) + '\t' + str(etefferr) + '\t' + str(elogg) + '\t' + str(eloggerr) + '\t' + str(H8teff) + '\t' + str(H8tefferr) + '\t' + str(H8logg) + '\t' + str(H8loggerr) + '\t' + str(H9teff) + '\t' + str(H9tefferr) + '\t' + str(H9logg) + '\t' + str(H9loggerr)+ '\t' + str(H10teff) + '\t' + str(H10tefferr) + '\t' + str(H10logg) + '\t' + str(H10loggerr) 
        f.write(info + '\n')

    print '\n'
    os.chdir('../')


f.close()
