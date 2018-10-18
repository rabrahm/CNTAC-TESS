import pandas as pd
import matplotlib.pyplot as plt
import astropy
import astropy.coordinates as coord
import astropy.units as u
import sys
import numpy as np
from astropy import constants as const

G    = const.G.value
Msun = const.M_sun.value
Rsun = const.R_sun.value
Mjup = const.M_jup.value
Rjup = const.R_jup.value
AU   = const.au.value
Mearth = const.M_earth.value
Rearth = const.R_earth.value

def get_K(P,mp,ms=1.,e=0.,i=.5*np.pi):
	#return the velocity semiamplitude: P->[day], mp->[jupiter masses], ms->[solar masses]
	K = 203. * P**(-1./3.) * mp * np.sin(i) * (1- e**2)**(-1./2.) / (ms + 9.548e-4*mp)**(2./3.)
	return K

def get_mass(rp1,tq1):
	nd = np.array([5673.,1.036,1.083,14.89,0.478,0.1204,0.315,0.027,0.027,0.847,0.013,0.013,9.677])
	#Ms P ecc sma Mp Rp
	#d = np.loadtxt('tepcat.txt',skiprows=0,usecols=[1,7,10,19,20,23,26,29])
	d = np.loadtxt('tepcat2.txt',skiprows=0,usecols=[1,7,10,19,20,23,26,27,28,29,30,31,43])
	d2 = np.loadtxt('tepcat2.txt',skiprows=0,usecols=[0],dtype='str')

	d = np.vstack((d,nd))
	d2 = np.hstack((d2,np.array(['NewPlanet'])))
	#print gfd
	Mp = d[:,6]
	Mpl = d[:,7]
	Mpu = d[:,8]
	Mpe = 0.5*(Mpl+Mpu)
	Rp = d[:,9]
	Rpl = d[:,10]
	Rpu = d[:,11]
	Rpe = 0.5*(Rpl+Rpu)
	Ms = d[:,1]
	sma = d[:,5]
	per = d[:,3]
	ecc = d[:,4]
	teff = d[:,0]
	Rs = d[:,2]
	mJ = d[:,-1]

	J = np.where((Mp>0.02)&(Rp!=-1)&(Mpe/Mp<0.40)&(Rp!=-1)&(Rpe/Rp<0.40)&(mJ>0))


	Mp,Rp,Ms,sma,per,ecc,teff,Rs,names,mJ = Mp[J],Rp[J],Ms[J],sma[J],per[J],ecc[J],teff[J],Rs[J],d2[J],mJ[J]
	#print tfrd


	#I = np.where((Mp>0.054)&(Mp<0.18)&(Rp>0.))
	#print names[I],fds
	#print I,fds
	#I = np.where((ecc>0.6) & (per<30))[0]
	#print I,fd

	"""
	A = (63./4.)*np.sqrt(G*(Ms*Msun)**3)*((Rp*Rjup)**5/(Qp*Mp*Mjup))
	tau_circ = A * (sma*AU)**(-6.5)
	tau_circ = 1. / tau_circ
	tau_circ = tau_circ/(3600.*24*365*1e9)
	I = np.where(tau_circ>100)[0]
	tau_circ[I] = 100.
	I = np.where(tau_circ<0.001)[0]
	tau_circ[I] = 0.001

	print tau_circ
	#print ds
	print per[-1],ecc[-1],np.log10(Mp[-1]),np.log10(tau_circ)[-1]
	print np.log10(Mp)
	#print fds

	"""
	tq = teff * np.sqrt((Rs*Rsun)/(2*sma*AU))
	st = 0.2
	while True:
		if tq1 < 1200:
		    I = np.where((tq<1200) & (Rp>rp1-rp1*st) & (Rp<rp1+rp1*st))[0]
		else:
		    I = np.where((tq<tq1+tq1*st) & (tq>tq1-tq1*st) & (Rp>rp1-rp1*st) & (Rp<rp1+rp1*st))[0]
		if len(I)>0:
			break
		else:
			st += 0.1


	return Mp[I][np.random.randint(len(I))]


barclay = pd.read_csv("detected_planets.csv",header=42)

# cuts
# DEC: < 25
# no distance cut, the TLC is all nearby
# Radius: > 4 R_earth
# insolation: < 150 F_earth


flux_limit = 2e8 / 1360800
print(flux_limit)

barclay_s = barclay.query('(PlanetRadius > 4) & (DEdeg < 0) & (Vmag < 11) & (PlanetPeriod > 10) ')

refrad = np.arange(0.4,1.,0.01)

teqs =  np.array((barclay_s['Insol']*1361/(4*5.67e-8))**(1./4.))
rpls =  np.array(barclay_s['PlanetRadius'])/11.2
eccs =  np.array(barclay_s['Ecc'])
Ps =  np.array(barclay_s['PlanetPeriod'])
mst =  np.array(barclay_s['StarMass'])

mpls = []
ks = []
for i in range(len(teqs)):
	print rpls[i],teqs[i]
	mpl = get_mass(rpls[i],teqs[i])
	K = get_K(Ps[i],mpl,mst[i],eccs[i])
	mpls.append(mpl)
	ks.append(K)
mpls = np.array(mpls)
ks = np.array(ks)
print teqs
print rpls
print mpls
plt.hist(np.log10(ks))
plt.axvline(np.log10(30),color='black')
plt.xlabel(r'log$_{10}$(K) [m/s]')
plt.ylabel(r'N')
print np.log10(0.5)
I1 = np.where(ks>30)[0]
I2 = np.where(mpls>0.3)[0]

print len(ks),len(I1),len(I2)
plt.show()
print ks
print fds
