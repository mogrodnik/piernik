#!/usr/bin/env python
import numpy as np

mode = 'simple'   # maps generated on the base of single CR protons distribution
#mode = 'spectral' # maps generated on the base of CR electrons spectrum

Hz = 1 ; MHz = Hz * 1.0e6 ; GHz = MHz * 1.0e3 ; clight = 3.e8 #m/s
metr = 1 ; cm = 0.01 * metr

const_synch = 1.e8
p_min_fix = 5.e0     #TODO !> if possible - read the below from h5file, if available
p_max_fix = 3.e5
ncre      = 14
maxcren   = ncre # 12
q_eps     = 0.01
q_big     = 30.
arr_dim   = 200      #TODO <! if possible - read the above from h5file, if available

# cJnu - a constant scaling synchrotron emissivity.
cJnu = 1.0
# x_ion - ionization degree of the thermal gas (a number in the range [0-1])
x_ion = 0.01 #1.0 #0.173
# spectral index of CR energy density - typically ~2.8
# In the spectral approach it is computed from the CR energy distribution
# In the simple approach it is set as a cosntant
p = 2.8
# Settings for the radiotelescope beam parameters
nbeam = 30
sigma = 1.0

xlo, xhi = 67, 317
zlo, zhi = 21, 171
klo, khi = -1, -1

#
# normalization paramer for vectors on the map (stg)
normalise_exponent_PI = 0.5     # 0.7
# vectors scaling for quiver method (stg)
scale_vec_PI = 0.01

print_vec = False
print_TP = False
print_PI = False
print_SI = False
print_RM = False
print_df = False
print_log = False

suffix = 'png', #'pdf',

def RMunit():
   # RM = 0,812 * B_{||} * rho_i * ds * (cm/pc)**3 * (mgs/muG) * (Msun/mH) * (102/136) [rad/m**2]
   mgs_muG = 2.8519
   cm_pc   = 3.2407793e-19
   Msun_g  = 1.9884e33
   mH_g    = 1.673559e-24
   RMunit = 0.812 * (cm_pc)**3 * (mgs_muG) * (Msun_g/mH_g) * 3./4.
   #print('RMunit = ',RMunit)
   return RMunit

def krange(n3,axis):
   if (axis ==2):
      klo, khi = zlo, zhi
   else:
      klo, khi = xlo, xhi
   kl, kh = klo, khi
   if klo < 0:
      kl = 0
   if khi < 0:
      kh = n3
   return kl, kh

# radio frequency
# nu = 0.150 ; 1.0   ; 1.4   ; 2.0   ; 3.0   ; 4.0    ; 5.0  ; 6.0  ; 7.0    ; 8.0     ; 9.0     ; 10.0 ; 12.0   ; 15.0 ; 20.0   ; 30.0
#lbd = 2 m   ; 30 cm ; 21 cm ; 15 cm ; 10 cm ; 7.5 cm ; 6 cm ; 5 cm ; 4.3 cm ; 3.75 cm ; 3.33 cm ; 3 cm ; 2.5 cm ; 2 cm ; 1.5 cm ; 1 cm
def set_nuandlbd(ns,ls,ss):
   if (ss == 1):
      nu  = 1.0*GHz
      lbd = 0.001*metr
   else:
      nu  = 30.0*GHz
      lbd = 0.01*metr
   if (ns >= 0.0) and (ls >= 0.0):
      print(ss, ': Do not set  nu and lambda simultaneously. I use nu')
   # radiation wavelength
   if (ns >= 0.0):
      nu = ns*GHz
      lbd = clight/nu
   elif (ls >= 0.0):
      lbd = ls
      nu = clight/lbd
   else:
      lbd = clight/nu
   print(ss, ': nu = ', nu/GHz, ' GHz, lambda = ', lbd/cm, ' cm')
   return nu, lbd

#labels:
lab_TP = 'TP'
lab_PI = 'PI'
lab_SI = 'SI'
lab_RM = 'RM'

# cbar_label - title of color bar
def ety(lab):
   etyk = lab
   if print_log:
      etyk = r'log$_{10}$('+lab+')'
   if (lab == lab_RM):
      etyk = etyk+r' $[rad / m^2]$'
   return etyk

# colormap - name of the colormap
def colormap(lab):
   if (lab == lab_TP):
      scmp = 'jet'
   elif (lab == lab_PI):
      scmp = 'jet'
   elif (lab == lab_SI):
      scmp = 'jet'
   elif (lab == lab_RM):
      scmp = 'RdBu_r'
   return scmp

# title - title of the plot
def title(lab,attr):
   time, l1, l2 = attr
   lbdcm = ' | '+str(l1/cm)+' cm'
   if lab == lab_SI:
      lbdcm = lbdcm+' & '+str(l2/cm)+' cm'
   if lab == lab_RM:
      lbdcm = ''
   return lab+lbdcm+' | t = '+str(round(time))+' Myr'

def etyfil(ax,file_name,nu):
   return "_"+ax+"_"+file_name.split('/')[-1]+'_nu'+str(int(nu/MHz)).zfill(3)+'MHz'

# dok - parameter determining spatial separation of vectors on the map
def dokv(ff,ax_set):
   if ff:
      if ax_set == 2:
         dok = 15
      else:
         dok = 15
   else:
      dok = sigma
   return dok

def dokvec(wp,wq,X,Y,axis,ff):
   dok = dokv(ff,axis)
   wp, wq, X, Y = np.array(wp), np.array(wq), np.array(X), np.array(Y)
   print(np.shape(wp))
   wp = wp[::dok,::dok]
   wq = wq[::dok,::dok]
   X  = X[::dok,::dok]
   Y  = Y[::dok,::dok]
   print('max vector ', np.amax(np.sqrt(wq**2 + wp**2)))
   return wp, wq, X, Y

# parameter determining normalization of colormaps
def norm(lab):
   if (lab == lab_TP):
      norml = 10.0
   elif (lab == lab_PI):
      norml = 10.0
   elif (lab == lab_SI):
      norml = 10.0
   elif (lab == lab_RM):
      norml = 0.003
   return norml

def fvmax(lab,ff,data):
   print, ff
   if print_df:
      if lab == 'RM':
         vmx = np.max(data)*norm(lab)
         return -vmx, vmx
      if ff:
         return 0.0, np.mean(data)*norm(lab)
      else:
         return 0.0, np.max(data)*norm(lab)
   else:
      if lab == 'TP':
         vmn, vmx = 0.0, 20.
      if lab == 'PI':
         vmn, vmx = 0.0 , 20000.
      if lab == 'SI':
         vmn, vmx = -2.0 , 0.0
      if lab == 'RM':
         vmn, vmx = -100., 100.
      if print_log:
         vmn, vmx = -6., -2.
   return vmn, vmx
