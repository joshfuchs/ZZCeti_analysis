import numpy as np



#===============================================
#===============================================
#===============================================

#Corrections from ML2/alpha = 0.8 to 3D from Tremblay et al. (2013)
def ml28_to_3d_teff(teff,logg):
    c = np.zeros(8)
    c[0] = 1.0947335E-03
    c[1] = -1.8716231E-01
    c[2] = 1.9350009E-02
    c[3] = 6.4821613E-01
    c[4] = -2.2863187E-01
    c[5] = 5.8699232E-01
    c[6] = -1.0729871E-01
    c[7] = 1.1009070E-01

    gx = logg - 8.0
    tx = (teff - 10000.)/1000.

    shift = c[0] + (c[1] + c[6]*tx + c[7]*gx) * np.exp(-(c[2]+c[4]*tx+c[5]*gx)**2. * (tx-c[3])**2.)
    
    teff3d = teff + (shift*1000.)
    
    return teff3d

#===============================================
#===============================================
#===============================================

def ml28_to_3d_logg(teff,logg):
    d = np.zeros(12)
    d[0] = 7.5209868E-04
    d[1] = -9.2086619E-01
    d[2] = 3.1253746E-01
    d[3] = -1.0348176E+01
    d[4] = 6.5854716E-01
    d[5] = 4.2849862E-01
    d[6] = -8.8982873E-02
    d[7] = 1.0199718E+01
    d[8] = 4.9277883E-02
    d[9] = -8.6543477E-01
    d[10] = 3.6232756E-03
    d[11] = -5.8729354E-02

    gx = logg - 8.0
    tx = (teff - 10000.)/1000.

    shift = (d[0] + d[4]*np.exp(-1.*d[5]*(tx-d[6])**2.)) + d[1]*np.exp(-1.*d[2]*(tx-(d[3]+d[7]*np.exp(-1.*(d[8]+d[10]*tx+d[11]*gx)**2.*(tx-d[9])**2.)))**2.)

    logg3d = logg + shift

    return logg3d


#===============================================
#===============================================
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

    #print combinedsmall.shape
    #print lowt, hight
    #print lowg, highg

    #Create finer small grid with spacing of 1 K and 0.005 logg
    lenteffgrid = np.round(hight-lowt+1.) #Round to ensure we get the correct number of points. Otherwise, occasionally face a strange int/float issue.
    teffsmallfine = np.linspace(lowt,hight,lenteffgrid,endpoint=True)
    lenlogggrid = np.round((highg-lowg)*1000.+1.)
    loggsmallfine = np.linspace(lowg,highg,lenlogggrid,endpoint=True)
    teffsmallfinegrid, loggsmallfinegrid = np.meshgrid(teffsmallfine,loggsmallfine)

    #####################

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
        #print teffsmall[x], loggval[x], chival[x]
        #plt.clf()
        #plt.plot(loggsmall,combinedsmall[:,x],'b^')
        #plt.plot(loggsmallfine,pc(loggsmallfine))
        #plt.show()

    #Now take these values and fit a polynomial in the Teff direction
    lowtanswer = lowt - 100.
    hightanswer = hight + 100.
    lenteffanswer = np.round(hightanswer-lowtanswer+1.)
    teffanswer = np.linspace(lowtanswer,hightanswer,lenteffanswer,endpoint=True)
    tpol = np.polyfit(teffsmall,chival,2)
    tpp = np.poly1d(tpol)

    bestteff = teffanswer[tpp(teffanswer).argmin()]
    #plt.clf()
    #plt.plot(teffsmallfine,tpp(teffsmallfine),'k')
    #plt.plot(teffanswer,tpp(teffanswer),'ro')
    #plt.plot(teffsmall,chival,'g^')
    #plt.plot(bestteff,tpp(bestteff),'ms')
    #plt.show()

    ########################
    #print 'Now logg'
    #Fit a polynomial to different logg values to find center of Teff
    teffval = np.zeros(len(combinedsmall[0,:]))
    chivalteff = np.zeros(len(combinedsmall[0,:]))
    for x in np.arange(len(combinedsmall[0,:])):
        polteff = np.polyfit(teffsmall,combinedsmall[x,:],2)
        pcteff = np.poly1d(polteff)
        teffval[x] = teffsmallfine[pc(teffsmallfine).argmin()]
        chivalteff[x] = pcteff(teffval[x])
        #print loggsmall[x], teffval[x], chivalteff[x]
        #plt.clf()
        #plt.plot(teffsmall,combinedsmall[x,:],'b^')
        #plt.plot(teffsmallfine,pcteff(teffsmallfine))
        #plt.show()

    #Now take these values and fit a polynomial in the Teff direction
    lowganswer = lowg - 0.25
    highganswer = highg + 0.25
    lenlogganswer = np.round((highganswer-lowganswer)*1000.+1.)
    logganswer = np.linspace(lowganswer,highganswer,lenlogganswer,endpoint=True)
    gpol = np.polyfit(loggsmall,chivalteff,2)
    gpp = np.poly1d(gpol)
    bestlogg = logganswer[gpp(logganswer).argmin()]
    #plt.clf()
    #plt.plot(loggsmallfine,gpp(loggsmallfine),'k')
    #plt.plot(logganswer,gpp(logganswer),'ro')
    #plt.plot(loggsmall,chivalteff,'g^')
    #plt.plot(bestlogg,gpp(bestlogg),'ms')
    #plt.show()
    
    ########################

    #Take these solutions and find the errors
    deltateff = tpp(teffanswer)-tpp(bestteff)
    lowterr = teffanswer[np.where(deltateff < 2.3)].min()-1.
    highterr =  teffanswer[np.where(deltateff < 2.3)].max()+1.
    tefferr = ((bestteff-lowterr) + (highterr-bestteff)) / 2.
    
    deltalogg = gpp(logganswer) - gpp(bestlogg)
    lowloggerr = logganswer[np.where(deltalogg < 2.3)].min()-0.001
    highloggerr = logganswer[np.where(deltalogg < 2.3)].max()+0.001
    loggerr = ((bestlogg-lowloggerr) + (highloggerr-bestlogg)) / 2.

    print 'Teff = ', bestteff , '+/-' , tefferr
    print 'logg = ', bestlogg,  '+/-' , loggerr
    #plt.clf()
    #plt.plot(loggsmallfine,deltalogg)
    #plt.show()
    return bestteff, tefferr, bestlogg, loggerr
