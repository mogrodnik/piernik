#!/usr/bin/env python
from math import *
import numpy as np
from draw_map import plotext
#import sys
from sys import stdout as sys_stdout
from stokes import stokes_params
import settings as stg
from settings import ncre
import yt

# wczytac dane przez yt

# Stokes przetwarza dane piksel po pikselu, wgłąb wybranej osi.

# przed podaniem do stokes_params dane nalezy uszeregowac poprzez
# sorted = np.argsort(ray_object["t"])
# i pozniej podawac poprzeez ray_object[field][sorted]

# Stokes powinien dostac zamiast ds wektor dlugosci komorek, byc moze calkowac trzeba bedax_map[z]e po polowie komorek.

ax_map = {'x': 0, 'y': 1, 'z': 2}
high_verbosity = False

def data_h5_yt(filename, ax_set, wave_data, imresw, imdepth):
# function reading data from an HDF5 file (.h5) and using it to compute Stokes parameters.
# filename   - filename.h5
# ax_set     - parameter storing the user's choice of aax_map[x]s along which the
#             projection is executed
# ax_map     - helper dictionary
# x_ion      - degree of ionization of the ionized gas (stg)
# ds         - cell size vector along the line of sight
# imres      - image resolution; number of rays taken in image axes
# imdepth    - range of integration along ax_set
   N_nl = stg.N_nulbd

   h5ds = yt.load(filename, unit_system="code")
   I, Q, U, RM, SI = [], [], [], [], []

   time = h5ds.current_time.v[0]
   DLedge = h5ds.domain_left_edge.v
   DRedge = h5ds.domain_right_edge.v

   #WARNING in read_data this are divided by 1000.
   # Default values of map limits are assumed here
   xmin = DLedge[ax_map["x"]]
   xmax = DRedge[ax_map["x"]]
   ymin = DLedge[ax_map["y"]]
   ymax = DRedge[ax_map["y"]]
   zmin = DLedge[ax_map["z"]]
   zmax = DRedge[ax_map["z"]]

   print("Domain dimensions: ", 'xmin, xmax =(', xmin,',', xmax, '), ymin, ymax =(', ymin,',', ymax, '), zmin, zmax =(', zmin,',', zmax,')')

   # Labels are defined below
   rho_lbl = "density"
   ecrp_lbl = "cr01"
   ecre_lbl = []
   ncre_lbl = []
   bx_lbl = "mag_field_x"
   by_lbl = "mag_field_y"
   bz_lbl = "mag_field_z"
   sorT_lbl = "t"
   cre_iter = list(range(ncre)) # Do not use range generator each time (saves time if stg.spectral_mode)
   for ic in cre_iter:
      ecre_lbl.append('cree'+str(ic+1).zfill(2))
      ncre_lbl.append('cren'+str(ic+1).zfill(2))

   [wlo, whi, hlo, hhi] = plotext(ax_set)
   [wlo, whi, hlo, hhi] = [wlo * 1000., whi * 1000., hlo * 1000., hhi * 1000.] # WARNING these are PROBABLY in kpc, data will be in pc, thus multiplying them by 1000

   # Prepare limits and pre-organize B arrays, according to the settings and ax_set
   if ax_set==0:
      bset = [by_lbl, bz_lbl, bx_lbl]
      if (imdepth) != "max":
         xmin, xmax = -imdepth, imdepth # else: already defined
      ymin, ymax = wlo, whi
      zmin, zmax = hlo, hhi
   if ax_set==1:
      bset = [bx_lbl, bz_lbl, by_lbl]
      if (imdepth) != "max":
         ymin, ymax = -imdepth, imdepth # else: already defined
      xmin, xmax = wlo, whi
      zmin, zmax = hlo, hhi
   if ax_set==2:
      bset = [bx_lbl, by_lbl, bz_lbl]
      if (imdepth) != "max":
         zmin, zmax = -imdepth, imdepth # else: already defined
      xmin, xmax = wlo, whi
      ymin, ymax = hlo, hhi

   # Effective integration limits
   rbeg = [xmin, ymin, zmin]
   rend = [xmax, ymax, zmax]
   if (imdepth == "max"):
      print("(WARNING) Parameter 'yt_depth' (depth of integration along 'ax_set') was not provided, assuming default maximal range: -%10.2f:%10.2f" %(rbeg[ax_set], rend[ax_set]))

   # Prepare indices pointing width-ax and height-ax of the demanded image
   i_w = ax_map["x"] if (ax_set == ax_map["y"] or ax_set == ax_map["z"]) else ax_map["y"]
   i_h = ax_map["z"] if (ax_set == ax_map["x"] or ax_set == ax_map["y"]) else ax_map["y"]

   # Check & prepare resolution of the image
   if (imresw == 0):
      imres = [0, 0]
      imres[0] = 64
      imres[1] = ceil( (rend[i_h] - rbeg[i_h] ) / (rend[i_w] - rbeg[i_w]) * imres[0])
      print("(WARNING) Parameter 'imres' (image resolution in width) was not provided, assuming default: %ix%i" %(imres[0], imres[1]) )
   else: # scale height-resolution to aspect ratio determined in settings.plotext
      imres = [0, 0]
      imres[0] = imresw
      imres[1] = ceil( (rend[i_h] - rbeg[i_h] ) / (rend[i_w] - rbeg[i_w]) * imresw)

   # Calculate pixel distances
   dw = ( rend[i_w] - rbeg[i_w] ) / imres[0]
   dh = ( rend[i_h] - rbeg[i_h] ) / imres[1]

   # Prepare arrays of synthetic observable quantities
   I  = np.zeros((N_nl, imres[0], imres[1]))
   Q  = np.zeros((N_nl, imres[0], imres[1]))
   U  = np.zeros((N_nl, imres[0], imres[1]))
   RM = np.zeros((N_nl, imres[0], imres[1]))
   SI = np.zeros((N_nl, imres[0], imres[1]))
   Ecrp = []

   #print(np.shape(I), N_nl)
   #print( I[0][16,16] )
   # Prepare the limits of plotted area to be returned, NOTICE formally these are the same as provided by settings.plotext
   x = np.linspace(rbeg[i_w] / 1000., rend[i_w] / 1000., imres[0])
   y = np.linspace(rbeg[i_h] / 1000., rend[i_h] / 1000., imres[1])
   figext = [(rbeg[i_w])/ 1000., (rend[i_w])/ 1000., (rbeg[i_h] - 0.5 * dw)/ 1000., (rend[i_h] - 0.5 * dw)/ 1000.]

   # Prepare variable coordinates for iteration in domain
   rbegv = [0, 0, 0]
   rendv = [0, 0, 0]
   rbegv[:] = rbeg[:]
   rendv[:] = rbeg[:]            # NOTICE only coordinate at ax_set is supposed to vary here between rbegv and rendv
   rendv[ax_set] = rend[ax_set]  # NOTICE remaining coordinates will start to vary during iteration

   rbegv[i_w] = rbegv[i_w] + 0.5 * dw
   rbegv[i_h] = rbegv[i_h] + 0.5 * dh
   print('Sendig data to map computation \n')
   # This part should works regardless of chosen ax_set
   for i in range(imres[0]):
      # Cursor up one line
      sys_stdout.write("\033[F")

      rbegv[i_w] = rbeg[i_w] + i * dw
      rendv[i_w] = rbegv[i_w]

      for j in range(imres[1]):
         rbegv[i_h] = rbeg[i_h] + j * dh
         rendv[i_h] = rbegv[i_h]

         R = h5ds.ray(rbegv, rendv)
         iRsort = np.argsort(R[sorT_lbl])
         lends    =  R.shape[0]

         if (lends == 0): # Countermeasure againt empty vector from yt - probably bug due to coordinates pointing right in between cells.
            rbegv[i_h] = rbeg[i_h] + j * dh + 1e-6 # WARNING: if pointing right in between cells I'm assuming it's supposed to be the next one
            rendv[i_h] = rbegv[i_h]
            R = h5ds.ray(rbegv, rendv)
            iRsort = np.argsort(R[sorT_lbl])
            lends    =  R.shape[0]

         # WARNING yt.ray object is unsorted by default, sorting adds ~20% to total execution time
         ds       =  R.fwidth[iRsort].v[:, ax_set]
         rho_ion  =  R[rho_lbl][iRsort].v * stg.x_ion
         Bp       =  R[bset[0]][iRsort].v
         Bq       =  R[bset[1]][iRsort].v
         Bn       =  R[bset[2]][iRsort].v

         if (stg.spectral_mode):
            Ecre = np.array([R[ecre_lbl[l]][iRsort].v for l in cre_iter])
            Ncre = np.array([R[ncre_lbl[l]][iRsort].v for l in cre_iter])
            plot_data_arrays = stokes_params(Bp, Bq, Bn, rho_ion, Ecrp, wave_data, ds, lends, Ecre=Ecre, Ncre=Ncre, ncre=ncre) # Does not contain khi-klo! - irrelevant for ray (supplying length of ds)
         else:
            Ecrp =  R[ecrp_lbl][iRsort].v
            plot_data_arrays = stokes_params(Bp, Bq, Bn, rho_ion, Ecrp, wave_data, ds, lends) # Does not contain khi-klo! - irrelevant for ray (supplying length of ds)

         if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
            I[:,i,j] = plot_data_arrays[0][:]
         if stg.print_PI or stg.print_SI or stg.print_vec:
            Q[:,i,j], U[:,i,j] = plot_data_arrays[1:3][:]
         if stg.print_RM:
            RM[:,i,j] = plot_data_arrays[3][:]
         if stg.print_SI:
            SI[:,i,j] = plot_data_arrays[4][:]

         if (high_verbosity): print("[%10.2f %10.2f %10.2f] , [%10.2f %10.2f %10.2f] | %3i / %3i | >>> %e" %(rbegv[0], rbegv[1], rbegv[2], rendv[0], rendv[1], rendv[2], i, j, I[i,j]) )
      print('Done: %5.1f%% (%i out of %i rows)' %(100 * float(i+1) / float(imres[0]), i+1, imres[0] ) )

   return I, Q, U, RM, SI, x, y, figext, time
