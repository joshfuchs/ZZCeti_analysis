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

os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/clemens/students/group/modelfitting/Koester_06/RESULTS')
#os.chdir('/srv/two/jtfuchs/Interpolated_Models/10teff05logg/RESULTS')
#os.chdir('/srv/two/jtfuchs/Interpolated_Models/Koester_ML2alpha06/RESULTS')

wdname = 'WD1422p095_930_blue_flux_cd32_1'
datedone = '06-06_4.75.txt'


#Set up grid
teff = np.linspace(11000,14000,13,endpoint=True)
logg = np.linspace(7.00,9.50,11,endpoint=True)

#Fine Grid
#teff = np.linspace(10000,13000,301,endpoint=True)
#logg = np.linspace(7.0,9.0,41,endpoint=True)

teffgrid, logggrid = np.meshgrid(teff,logg)

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

#Read in saved grids
allchi  = np.genfromtxt(allfile,dtype='d')
#alphachi = np.genfromtxt(alphafile,dtype='d')
betachi = np.genfromtxt(betafile,dtype='d')
gammachi = np.genfromtxt(gammafile,dtype='d')
deltachi = np.genfromtxt(deltafile,dtype='d')
epsilonchi = np.genfromtxt(epsilonfile,dtype='d')
H8chi = np.genfromtxt(H8file,dtype='d')
H9chi = np.genfromtxt(H9file,dtype='d')
H10chi = np.genfromtxt(H10file,dtype='d')

#combine different lines
print allchi.shape
combined =  betachi + gammachi + deltachi + epsilonchi + H8chi + H9chi + H10chi


#specify a portion of the grid to extract
#lowg, highg = 7.75, 8.25
#lowgindex, highgindex = np.where(logg == lowg), np.where(logg == highg)
#loggsmall = logg[lowgindex[0]:highgindex[0]+1]

#lowt, hight = 12250., 12750.
#lowtindex, hightindex = np.where(teff == lowt), np.where(teff == hight)
#teffsmall = teff[lowtindex[0]:hightindex[0]+1]


#plot wireframe and scatter plot of chi-square values
#fig = plt.figure(1)
#ax = fig.gca(projection='3d')
#sur = ax.plot_wireframe(teffgrid,logggrid,combined,rstride=1,cstride=1)
#ax.scatter(np.ravel(teffsmallgrid),np.ravel(loggsmallgrid),np.ravel(combinedsmall),marker='o',s=30,c='r')
#plt.show()
#exit()

#Determine minimum values of each grid
allindex = np.unravel_index(allchi.argmin(),allchi.shape)
alllogg, allteff = logg[allindex[0]], teff[allindex[1]]
print 'All: ' , alllogg, allteff

#alphaindex = np.unravel_index(alphachi.argmin(),alphachi.shape)
#alphalogg, alphateff = logg[alphaindex[0]], teff[alphaindex[1]]
#print 'Alpha: ' , alphalogg, alphateff

betaindex = np.unravel_index(betachi.argmin(),betachi.shape)
betalogg, betateff = logg[betaindex[0]], teff[betaindex[1]]
print 'Beta: ' , betalogg, betateff

gammaindex = np.unravel_index(gammachi.argmin(),gammachi.shape)
gammalogg, gammateff = logg[gammaindex[0]], teff[gammaindex[1]]
print 'Gamma: ' , gammalogg, gammateff

deltaindex = np.unravel_index(deltachi.argmin(),deltachi.shape)
deltalogg, deltateff = logg[deltaindex[0]], teff[deltaindex[1]]
print 'Delta: ' , deltalogg, deltateff

epsilonindex = np.unravel_index(epsilonchi.argmin(),epsilonchi.shape)
epsilonlogg, epsilonteff = logg[epsilonindex[0]], teff[epsilonindex[1]]
print 'Epsilon: ' , epsilonlogg, epsilonteff

H8index = np.unravel_index(H8chi.argmin(),H8chi.shape)
H8logg, H8teff = logg[H8index[0]], teff[H8index[1]]
print 'H8: ' , H8logg, H8teff

H9index = np.unravel_index(H9chi.argmin(),H9chi.shape)
H9logg, H9teff = logg[H9index[0]], teff[H9index[1]]
print 'H9: ' , H9logg, H9teff

H10index = np.unravel_index(H10chi.argmin(),H10chi.shape)
H10logg, H10teff = logg[H10index[0]], teff[H10index[1]]
print 'H10: ' , H10logg, H10teff

combinedindex = np.unravel_index(combined.argmin(),combined.shape)
combinedlogg, combinedteff = logg[combinedindex[0]], teff[combinedindex[1]]
print 'Combined: ' , combinedlogg, combinedteff
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
elif combinedindex[0]+rangeg > len(logg):
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
elif combinedindex[1]+ranget > len(teff):
    ##teffsmall = teff[combinedindex[1]-ranget:-1]
    teffsmall = teff[-2*rangeg-1:]
    ##tefflow = combinedindex[1]-ranget
    tefflow = -2*rangeg-1
    teffhigh = -1
else:
    teffsmall = teff[combinedindex[1]-ranget:combinedindex[1]+ranget+1]
    tefflow = combinedindex[1]-ranget
    teffhigh = combinedindex[1]+ranget+1



#Get the low and high values for each
lowg, highg = loggsmall[0], loggsmall[-1]
lowt, hight = teffsmall[0], teffsmall[-1]

teffsmallgrid, loggsmallgrid = np.meshgrid(teffsmall,loggsmall)
if logghigh == -1:
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
tpol = np.polyfit(teffsmall,chival,3)
tpp = np.poly1d(tpol)

bestteff = teffsmallfine[tpp(teffsmallfine).argmin()]
print chival, chival.min()
print tpp(bestteff)
#plt.clf()
#plt.plot(teffsmall,chival,'g^')
#plt.plot(teffsmallfine,tpp(teffsmallfine),'k')
#plt.plot(bestteff,tpp(bestteff),'ms')
#plt.show()

#Take these solutions and find the errors
deltateff = tpp(teffsmallfine)-tpp(bestteff)
lowterr = teffsmallfine[np.where(deltateff < 2.3)].min()-1.
highterr =  teffsmallfine[np.where(deltateff < 2.3)].max()+1.

deltalogg = pckeep(loggsmallfine) - pckeep(loggval.min())
lowloggerr =  loggsmallfine[np.where(deltalogg < 2.3)].min()-0.001
highloggerr =  loggsmallfine[np.where(deltalogg < 2.3)].max()+0.001

print 'Teff = ', bestteff , '+/-' , highterr-bestteff
print 'logg = ', loggval.min(),  '+/-' , highloggerr-loggval.min()

#plt.clf()
#plt.plot(loggsmallfine,deltalogg)
#plt.show()

interpolation = RectBivariateSpline(loggsmall,teffsmall,combinedsmall,kx=3,ky=3,s=0) 
lowchi = interpolation(loggval.min(),bestteff)
levels = range(0,1000,300)
#plot contour plot
plt.figure()
CS = plt.contourf(teffsmall,loggsmall,combinedsmall-lowchi,levels=levels)
plt.colorbar(CS)
plt.plot(bestteff,loggval.min(),'^')
plt.xlim(bestteff-250.,bestteff+250.)
plt.ylim(loggval.min()-0.25,loggval.min()+0.25)
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
