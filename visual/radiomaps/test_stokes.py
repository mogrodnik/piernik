#!/usr/bin/env python
from math import *
import numpy as np
import pylab as py
from stokes import stokes_params
import settings as stg

# This program is aiming to test the routine 'stokes_params along a single line of sight
B0=1.0  # indukcja pola magnetycznego
cJnu=1.0 # stala intensywnosci promieniowania synchrotronowego
x_ion = 1.0 #stopien jonizacji
pp=1.0   #indeks widmowy energii elektrycznej, (wykladnik potegi rozkladu elektronow)
Nu=1.0 # czestotliwosc promieniowania
Ecr0=1.0#kappa, gestosc energii promieniowania kosmicznego
dd0=1.0 #gestosc gazu
Zmin=0.0
lbd=1.0 #dlugosc promieniowania
Zmax=2.0*pi/(0.812*lbd**2*dd0*B0*x_ion)
n=100 # ilosc komorek
RM_suma=0.0
dz=(Zmax-Zmin)/n
rho=np.zeros(n) # gestosc materii
Z=np.linspace(Zmin,Zmax,n) #os wzdluz ktorej patrzymy
Ecr = np.zeros(n)
#INDUKCJA
Bp = np.zeros(n) #os x
Bq = np.zeros(n) #os y
Bn = np.zeros(n) #os z
PI=(pp+1.0)/(pp+7.0/3.0) #poczatkowy stopien polaryzacji
for ii in range (n):
	#if Z[ii]>=0.0 and Z[ii] <= 10.0:
	rho[ii]= dd0
	Ecr[ii] = Ecr0
	Bp[ii] = B0 # 1.0
	Bq[ii] = B0 # 0.0
	Bn[ii] = B0 # 0.0

wave_data = [[Nu], [lbd]]
stg.p 	= pp
stg.cJnu = cJnu
stg.N_nulbd = len([Nu])
stg.test = True
stg.print_TP, stg.print_PI = True, True
stg.x_ion = x_ion
plot_data_arrays=stokes_params(Bp,Bq,Bn,rho,Ecr,wave_data,dz,n)
					  #stokes_params(Bp,Bq,Bn,rho_ion,Ecrp,wave_data,ds,n3,Ecre=[],Ncre=[],ncre=0):
print('I= ',plot_data_arrays[0],'Q= ',plot_data_arrays[1],'U= ',plot_data_arrays[2])

Q=plot_data_arrays[1]
U=plot_data_arrays[2]
psi_obs=0.5*np.arctan2(U,Q)
print("PI / I =", np.sqrt(Q[0]**2 + U[0]**2) / plot_data_arrays[0][0])
py.figure(1)
py.subplot(211)
py.title('Q')
py.xlabel('Z')
py.ylabel(r'$\Delta Q$')
py.plot(plot_data_arrays[5])
py.subplot(212)
py.title("U")
py.ylabel(r'$\Delta U$')
py.xlabel('Z')
py.plot(plot_data_arrays[6])
py.show()
