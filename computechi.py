'''
Written May 2016 by JTF

Reads in grids of chi-square values and computes minimum for each one.

Eventually will want to create surface plots.
Can use Axes3D.scatter to plot individual points

To Do:
- Determine actual minumun chi square value at lowest point for plotting. 
'''

import numpy as np
import os
import mpfit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RectBivariateSpline
import analysis_tools as at


def parabola(x,p):
    return p[0] + p[1]*x + p[2]*x**2.

def fitparabola(p,fjac=None,x=None,y=None,err=None):
    model = parabola(x,p)
    status = 0
    return([status,(y-model)/err])

def polynomial2(x,y,p):
    return p[0] + p[1]*x + p[2]*x**2. + p[3]*x*y + p[4]*y**2. + p[5]*y

def fitpolynomial2(p,fjac=None,x=None,y=None,z=None,err=None):
    model = polynomial2(x,y,p)
    status = 0
    return([status,(z-model)/err])

def polynomial3(x,y,p):
    return p[0] + p[1]*x + p[2]*x**2. + p[3]*x**3. + p[4]*x**2.*y + p[5]*x*y + p[6]*x*y**2. + p[7]*y**3. + p[8]*y**2. + p[9]*y

def fitpolynomial3(p,fjac=None,x=None,y=None,z=None,err=None):
    model = polynomial3(x,y,p)
    status = 0
    return([status,(z-model)/err])

def paraboloid(x,y,p):
    return p[0]*(((x-p[1])/p[2])**2. + ((y-p[3])/p[4])**2.) + p[5]

def fitparaboloid(p,fjac=None,x=None,y=None,z=None,err=None):
    model = paraboloid(x,y,p)
    status = 0
    return([status,(z-model)/err])


#os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/clemens/students/group/modelfitting/Koester_06/RESULTS')


wdname = 'wcftb.WD1422p095_930_blue_flux_master'
datedone = '06-20_4.28.txt'
cutteff = False
teff_limit = 14000.


#Set up filenames to read
allfile = 'chi_' + wdname + '_' + datedone
print 'File: ', allfile
alphafile = 'chi_' + wdname + '_alpha_' + datedone
betafile = 'chi_' + wdname + '_beta_' + datedone
gammafile = 'chi_' + wdname + '_gamma_' + datedone
deltafile = 'chi_' + wdname + '_delta_' + datedone
epsilonfile = 'chi_' + wdname + '_epsilon_' + datedone
H8file = 'chi_' + wdname + '_H8_' + datedone
H9file = 'chi_' + wdname + '_H9_' + datedone
H10file = 'chi_' + wdname + '_H10_' + datedone

#Read in first grid to determine range of spacing of variables

with open(allfile,'r') as f:
    first_line = f.readline()
try:
    bottomg,stepg,topg,bottomt,stept,topt,numpoints = [float(x) for x in first_line[2:].split(",")]
except:
    bottomg,stepg,topg,bottomt,stept,topt = [float(x) for x in first_line[2:].split(",")]
teff = np.linspace(bottomt,topt,(topt-bottomt)/stept+1.,endpoint=True)
logg = np.linspace(bottomg,topg,(topg-bottomg)/stepg+1.,endpoint=True)
teffgrid, logggrid = np.meshgrid(teff,logg)
'''
#Set up grid. This is saved in the header of the chi*txt file
bottomt = 10000.
topt = 15000.
stept = 10.
teff = np.linspace(bottomt,topt,(topt-bottomt)/stept+1.,endpoint=True)

bottomg = 7.0 
topg = 9.5
stepg = 0.05
logg = np.linspace(bottomg,topg,(topg-bottomg)/stepg+1.,endpoint=True)

teffgrid, logggrid = np.meshgrid(teff,logg)
'''




#Read in saved grids
allchi  = np.genfromtxt(allfile,dtype='d')
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
'''
#Convert to reduced chi-square
allchi /= numpoints
try:
    alphachi /= numpoints
except:
    pass
betachi /= numpoints
gammachi /= numpoints
deltachi /= numpoints
epsilonchi /= numpoints
H8chi /= numpoints
H9chi /= numpoints
H10chi /= numpoints
'''
#combine different lines
print 'Shape: ', allchi.shape
combined =  betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi#alphachi + betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
try:
    a10chi = alphachi + betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
except:
    pass
b10chi = betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
g10chi = gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi
g9chi = gammachi + deltachi + epsilonchi + H8chi + H9chi
g8chi = gammachi + deltachi + epsilonchi + H8chi
b9chi = betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi
b8chi = betachi + gammachi + deltachi + epsilonchi + H8chi


#specify a portion of the grid to extract
#lowg, highg = 7.75, 8.25
#lowgindex, highgindex = np.where(logg == lowg), np.where(logg == highg)
#loggsmall = logg[lowgindex[0]:highgindex[0]+1]

#lowt, hight = 12250., 12750.
#lowtindex, hightindex = np.where(teff == lowt), np.where(teff == hight)
#teffsmall = teff[lowtindex[0]:hightindex[0]+1]


#plot wireframe and scatter plot of chi-square values
fig = plt.figure(1)
ax = fig.gca(projection='3d')
sur = ax.plot_wireframe(teffgrid,logggrid,combined,rstride=1,cstride=1)
#ax.scatter(np.ravel(teffgrid),np.ravel(logggrid),np.ravel(combined),marker='o',s=30,c='r')
#plt.show()
#exit()

#Determine minimum values of each grid
allindex = np.unravel_index(allchi.argmin(),allchi.shape)
alllogg, allteff = logg[allindex[0]], teff[allindex[1]]
#print 'All: ' , alllogg, allteff

try:
    alphaindex = np.unravel_index(alphachi.argmin(),alphachi.shape)
    alphalogg, alphateff = logg[alphaindex[0]], teff[alphaindex[1]]
    #print 'Alpha: ' , alphalogg, alphateff
except:
    pass

betaindex = np.unravel_index(betachi.argmin(),betachi.shape)
betalogg, betateff = logg[betaindex[0]], teff[betaindex[1]]
#print 'Beta: ' , betalogg, betateff

gammaindex = np.unravel_index(gammachi.argmin(),gammachi.shape)
gammalogg, gammateff = logg[gammaindex[0]], teff[gammaindex[1]]
#print 'Gamma: ' , gammalogg, gammateff

deltaindex = np.unravel_index(deltachi.argmin(),deltachi.shape)
deltalogg, deltateff = logg[deltaindex[0]], teff[deltaindex[1]]
#print 'Delta: ' , deltalogg, deltateff

epsilonindex = np.unravel_index(epsilonchi.argmin(),epsilonchi.shape)
epsilonlogg, epsilonteff = logg[epsilonindex[0]], teff[epsilonindex[1]]
#print 'Epsilon: ' , epsilonlogg, epsilonteff

H8index = np.unravel_index(H8chi.argmin(),H8chi.shape)
H8logg, H8teff = logg[H8index[0]], teff[H8index[1]]
#print 'H8: ' , H8logg, H8teff

H9index = np.unravel_index(H9chi.argmin(),H9chi.shape)
H9logg, H9teff = logg[H9index[0]], teff[H9index[1]]
#print 'H9: ' , H9logg, H9teff

H10index = np.unravel_index(H10chi.argmin(),H10chi.shape)
H10logg, H10teff = logg[H10index[0]], teff[H10index[1]]
#print 'H10: ' , H10logg, H10teff

combinedindex = np.unravel_index(combined.argmin(),combined.shape)
combinedlogg, combinedteff = logg[combinedindex[0]], teff[combinedindex[1]]
#print 'Combined: ' , combinedlogg, combinedteff
#exit()

#Print the chi-square value of a particular grid at a particular point
#loggwant = np.abs(logg-8.25).argmin()
#teffwant = np.abs(teff-13000).argmin()
#print allchi[loggwant,teffwant]

#Print values along a particular row
#teffwant = np.where(teff == 13900)
#loggwant = np.where(logg == 7.95)
#print H10chi[loggwant,:]
#print H10chi[:,teffwant]
#plt.clf()
#plt.plot(teff,np.array(H10chi[loggwant,:][0][0]))
#plt.show()


#======================================
#Remove part of the chi-square grid in a secondary solution is being found.
if cutteff:
    ########
    #Upper limit on Teff
    ########
    #teff_limit = 14500.
    print combined.shape
    teffcut = np.abs(teff-teff_limit).argmin()
    print teffcut
    teff = teff[0:teffcut]
    #print len(teff_new)
    combined = combined[:,0:teffcut]
    betachi = betachi[:,0:teffcut]
    H9chi = H9chi[:,0:teffcut]
    H8chi = H8chi[:,0:teffcut]
    H10chi = H10chi[:,0:teffcut]
    b10chi = b10chi[:,0:teffcut]
    g10chi = g10chi[:,0:teffcut]
    g9chi = g9chi[:,0:teffcut]
    deltachi = deltachi[:,0:teffcut]
    epsilonchi = epsilonchi[:,0:teffcut]
    b9chi = b9chi[:,0:teffcut]
    b8chi = b8chi[:,0:teffcut]
    g8chi = g8chi[:,0:teffcut]
    try:
        alphachi = alphachi[:,0:teffcut]
        a10chi = a10chi[:,0:teffcut]
    except:
        pass

'''
#######
#Lower limit on Teff
teff_limit = 14700.
print combined.shape
teffcut = np.abs(teff-teff_limit).argmin()
print teffcut
teff = teff[teffcut:]
#print len(teff_new)
combined = combined[:,teffcut:]
H10chi = H10chi[:,teffcut:]
'''

#Find solution for whatever combinations you want
try:
    print '\nAlpha:'
    ateff,atefferr,alogg,aloggerr = at.find_solution(alphachi,logg,teff)
    print '\nAlpha-10:'
    a10teff,a10tefferr,a10logg,a10loggerr = at.find_solution(a10chi,logg,teff)
except:
    pass
#exit()

print '\nBeta:'
bteff,btefferr,blogg,bloggerr = at.find_solution(betachi,logg,teff)

print '\nGamma:'
gteff,gtefferr,glogg,gloggerr = at.find_solution(gammachi,logg,teff)
print '\nDelta:'
dteff,dtefferr,dlogg,dloggerr = at.find_solution(deltachi,logg,teff)
print '\nEpsilon:'
eteff,etefferr,elogg,eloggerr = at.find_solution(epsilonchi,logg,teff)

print '\nH8:'
H8teff,H8tefferr,H8logg,H8loggerr = at.find_solution(H8chi,logg,teff)


print '\nH9:'
H9teff,H9tefferr,H9logg,H9loggerr = at.find_solution(H9chi,logg,teff)


print '\nH10:'
H10teff,H10tefferr,H10logg,H10loggerr = at.find_solution(H10chi,logg,teff)

print '\nBeta - H10:'
b10teff,b10tefferr,b10logg,b10loggerr = at.find_solution(b10chi,logg,teff)

print '\nGamma - H10:'
g10teff,g10tefferr,g10logg,g10loggerr = at.find_solution(g10chi,logg,teff)

print '\nGamma - H9:'
g9teff,g9tefferr,g9logg,g9loggerr = at.find_solution(g9chi,logg,teff)

print '\nBeta - H9:'
b9teff,b9tefferr,b9logg,b9loggerr = at.find_solution(b9chi,logg,teff)
print '\nBeta - H8:'
b8teff,b8tefferr,b8logg,b8loggerr = at.find_solution(b8chi,logg,teff)

print '\nGamma - H8:'
g8teff,g8tefferr,g8logg,g8loggerr = at.find_solution(g8chi,logg,teff)


#print '\nCombined:'
#combinedteff,combinedtefferr,combinedlogg,combinedloggerr = at.find_solution(combined,logg,teff)

#exit()

#interpolation = RectBivariateSpline(loggsmall,teffsmall,combinedsmall,kx=3,ky=3,s=0) 
interpolation = RectBivariateSpline(logg,teff,combined,kx=3,ky=3,s=0) 
#lowchi = interpolation(loggval.min(),bestteff)
levels = [1,2,3,10,100,200,300,400,500,600,700] # range(0,1000,300)
#plot contour plot
plt.figure()
#CS = plt.contour(teff,loggsmall,combinedsmall-lowchi)#,levels=levels)
CS = plt.contourf(teff,logg,b9chi,100,cmap='jet')#,levels=levels)
plt.colorbar(CS)
plt.xlim(15000,10000)
plt.ylim(9.5,7.0)
#plt.plot(bestteff,loggval.min(),'^')
#plt.xlim(bestteff+250.,bestteff-250.)
#plt.ylim(loggval.min()+0.25,loggval.min()-0.25)
#plt.clabel(CS,inline=1,fontsize=10)
plt.show()

#Check out the following with a smaller grid
###cs = plt.pcolor(teffsmall,loggsmall,combinedsmall-tpp(bestteff))
###cb = plt.colorbar(cs)








exit() #Below this is some code to fit an elliptic paraboloid to the surface, as well as doing a cubic spline interpolation. These are just other options.





'''
#Fit a different order polynomial
guess = np.zeros(10)
xes = allchi[5,:8]
yes = allchi[:8,5]
pol = np.polyfit(xes,yes,3)
pol2 = np.polyfit(yes,xes,3)
guess[0] = allchi.min()
guess[1] = pol[2]
guess[2] = pol[1]
guess[3] = pol[0]
guess[4] = 1.
guess[5] = 1.
guess[6] = 5.
guess[7] = pol2[0]
guess[8] = pol2[1]
guess[9] = pol2[2]
fa = {'x':np.ravel(teffgrid),'y':np.ravel(logggrid),'z':np.ravel(allchi),'err':np.ravel(error)}
params = mpfit.mpfit(fitpolynomial3,guess,functkw=fa,quiet=True)
zz = polynomial3(teffgrid,logggrid,params.params)

#Fine minimum of fit from coarse grid
fitindex = np.unravel_index(zz.argmin(),zz.shape)
fitlogg, fitteff = logg[fitindex[0]],teff[fitindex[1]]
print 'Fit: ',  fitlogg, fitteff

zztest = polynomial3(tefftestgrid,loggtestgrid,params.params)
fitindextest = np.unravel_index(zztest.argmin(),zztest.shape)
fitloggtest, fittefftest = loggtest[fitindextest[0]],tefftest[fitindextest[1]]
print 'Fit: ', fitloggtest, fittefftest

#Plot all Chi square points and the fit
fig3 = plt.figure(3)
ax3 = fig3.gca(projection='3d')
surf3 = ax3.plot_surface(teffgrid,logggrid,zz,rstride=1,cstride=1,shade=False,cmap='jet')
plt.draw()
ax3.scatter(np.ravel(teffgrid),np.ravel(logggrid),np.ravel(allchi),marker='o',s=30)
surf3.set_edgecolors(surf3.to_rgba(surf3._A))
surf3.set_facecolors('white')
#plt.show()

#Calculate residuals and show those too
residuals = zz - allchi
fig4 = plt.figure(4)
ax4 = fig4.gca(projection='3d')
surf4 = ax4.plot_surface(teffgrid,logggrid,residuals,rstride=1,cstride=1,shade=False,cmap='jet')
#plt.draw() #use this if you don't want it filled in
surf4.set_edgecolors(surf4.to_rgba(surf4._A))
surf4.set_facecolors('white')
plt.show()
'''


#Try fitting a polynomial to the smaller subset
error = np.ones([len(loggsmall),len(teffsmall)])
#print loggsmall
#print teffsmall
xes = combinedsmall[len(loggsmall)//2,:]
#print xes
yes = combinedsmall[:,len(teffsmall)/2]
#print yes
pol = np.polyfit(teffsmall,xes,2)
pol2 = np.polyfit(loggsmall,yes,2)

#2-order in both directions
#guess = np.zeros(6)
#guess[0] = combinedsmall.min()
#guess[1] = pol[1]
#guess[2] = pol[0]
#guess[3] = 1.
#guess[4] = pol2[1]
#guess[5] = pol2[0]
'''
#3-order in both directions
guess = np.zeros(10)
guess[0] = combinedsmall.min()
guess[1] = pol[2]
guess[2] = pol[1]
guess[3] = pol[0]
guess[4] = 1.
guess[5] = 1.
guess[6] = 5.
guess[7] = pol2[0]
guess[8] = pol2[1]
guess[9] = pol2[2]
'''
#Elliptic paraboloid

guess = np.zeros(6)
guess[0] = 0.4
guess[1] = teff[combinedindex[1]]
guess[2] = (teffsmall[1] - teffsmall[0]) / (combinedsmall[0,1] - combinedsmall[0,0] ) #15.
guess[3] = logg[combinedindex[0]]
guess[4] = (loggsmall[1] - loggsmall[0]) / (combinedsmall[1,0] - combinedsmall[0,0] )#0.005
guess[5] = combinedsmall.min()
#print guess
fa = {'x':np.ravel(teffsmallgrid),'y':np.ravel(loggsmallgrid),'z':np.ravel(combinedsmall),'err':np.ravel(error)}
#params = mpfit.mpfit(fitpolynomial2,guess,functkw=fa,quiet=True)
params = mpfit.mpfit(fitparaboloid,guess,functkw=fa,quiet=True,maxiter=1000)
print params.status, params.niter, params.fnorm, params.dof


#zz = polynomial2(teffsmallgrid,loggsmallgrid,params.params)
zz = paraboloid(teffsmallgrid,loggsmallgrid,params.params)
guessfit = paraboloid(teffsmallgrid,loggsmallgrid,guess)
#Fine minimum of fit from coarse grid
fitindex = np.unravel_index(zz.argmin(),zz.shape)
fitlogg, fitteff = loggsmall[fitindex[0]],teffsmall[fitindex[1]]
#print 'Fit: ',  fitlogg, fitteff


zztest = paraboloid(teffsmallfinegrid,loggsmallfinegrid,params.params)
#zztest = polynomial2(teffsmallfinegrid,loggsmallfinegrid,params.params)
fitindextest = np.unravel_index(zztest.argmin(),zztest.shape)
fitloggtest, fittefftest = loggsmallfine[fitindextest[0]],teffsmallfine[fitindextest[1]]
print 'logg fit: ', fitloggtest
print 'Teff fit: ', fittefftest


#Plot all Chi square points and the fit
fig3 = plt.figure(3)
ax3 = fig3.gca(projection='3d')
surf3 = ax3.plot_surface(teffsmallgrid,loggsmallgrid,zz,rstride=1,cstride=1,shade=False,cmap='jet')
#surf3 = ax3.plot_surface(teffsmallgrid,loggsmallgrid,guessfit,rstride=1,cstride=1,shade=False,cmap='jet')
#surf3 = ax3.plot_surface(teffsmallfinegrid,loggsmallfinegrid,zztest,rstride=1,cstride=1,shade=False,cmap='jet')
plt.draw()
ax3.scatter(np.ravel(teffsmallgrid),np.ravel(loggsmallgrid),np.ravel(combinedsmall),marker='o',s=30)
ax3.scatter(fittefftest,fitloggtest,zztest.min(),marker='o',c='r',s=60)
surf3.set_edgecolors(surf3.to_rgba(surf3._A))
surf3.set_facecolors('white') 
#plt.show()

#Calculate residuals and show those too
residuals = zz - combinedsmall
fig4 = plt.figure(4)
ax4 = fig4.gca(projection='3d')
surf4 = ax4.plot_surface(teffsmallgrid,loggsmallgrid,residuals,rstride=1,cstride=1,shade=False,cmap='jet')
#plt.draw() #use this if you don't want it filled in
surf4.set_edgecolors(surf4.to_rgba(surf4._A))
surf4.set_facecolors('white')
#plt.show()


#Find delta chi square == 1 surface
deltazztest = zztest - zztest.min()
oldglist = []
oldtlist = []

n,m = 0,0
for j in teffsmallfine:
    m = 0
    for i in loggsmallfine:
        if deltazztest[m,n] <= 1.:
            oldglist.append(i)
            oldtlist.append(j)
        m += 1
    n += 1
#print np.amin(oldglist),np.amax(oldglist)
#print np.amin(oldtlist),np.amax(oldtlist)

print 'logg error: ',(np.amax(oldglist)-np.amin(oldglist))/2.
print 'Teff error: ',(np.amax(oldtlist)-np.amin(oldtlist))/2.

plt.show()



#===================
#cubic spline interpolation of grid
interpolation = RectBivariateSpline(logg,teff,combined,kx=3,ky=3,s=0) 


glist = []
tlist = []
newgrid = np.empty([len(loggsmallfine),len(teffsmallfine)])
print 'Reading off new values'
n,m = 0,0
for j in teffsmallfine:
    m = 0
    for i in loggsmallfine:
        newgrid[m,n] = interpolation(i,j)
        if interpolation(i,j) <= 1000: #number is the max delta chi square we want
            #print i,j,out(i,j)
            glist.append(i)
            tlist.append(j)
        m += 1
    n += 1
print 'Done reading off new values'

interpindex = np.unravel_index(newgrid.argmin(),newgrid.shape)
interplogg, interpteff = loggsmallfine[interpindex[0]], teffsmallfine[interpindex[1]]
print 'Interpolation: ' , interplogg, interpteff
