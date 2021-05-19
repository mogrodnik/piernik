#!/usr/bin/env python

from math import *
import numpy as np
import h5py as h5
import sys
from stokes import stokes_params
from analytical_data import data_raw
import settings as stg
from settings import ncre
#from __future__ import print_function
#import electrons


def data_a(ax_set,wave_data):

# The function is computing Stokes parameters for the disk constructed in analytical_data.py routine
#nx, ny, nz - resolution of the domain in x, y and z directions.
#ax_set     - parameter storing the user's choice of axis along which the
#             projection is executed
#x_ion      - degree of ionization of the ionized gas (stg)
#nu         - observation frequency of synchrotron radiation
#ldb        - observation  wavelength of synchrotron radiation
#cJnu       - a constant scaling synchrotron emissivity.
#ds         - cell size along the line of sight
   I, Q, U, RM, SI = [], [], [], [], []
   # - current cell size along the line of sight for analytical data
   ds = 0.01
   # number of cells along x,y and z directions for analytiacl data
   nx, ny, nz = 50, 50, 50

# Initialization of the arrays storing physical data
   shape = (nx,ny,nz)
   Bp = np.zeros(shape)
   Bq = np.zeros(shape)
   Bn = np.zeros(shape)
   Ecrp = np.zeros(shape)
   Ecre = np.zeros((stg.ncre,nx,ny,nz))
   rho_ion = np.zeros(shape)

# filling the arrays with analytical data
   data = data_raw(nx,ny,nz)
   rho_ion = data[3]
   Ecrp = data[4]
   Ecre = data[5]

   if ax_set == 0:
      idr = [2,1,0]
      ndr = [nz,ny,nx]
   elif ax_set == 1:
      idr = [0,2,1]
      ndr = [nx,nz,ny]
   elif ax_set == 2:
      idr = [0,1,2]
      ndr = [nx,ny,nz]
   if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
      I  = np.zeros((ndr[0],ndr[1]))
   if stg.print_PI or stg.print_SI or stg.print_vec:
      Q  = np.zeros((ndr[0],ndr[1]))
      U  = np.zeros((ndr[0],ndr[1]))
   if stg.print_RM:
      RM = np.zeros((ndr[0],ndr[1]))
   if stg.print_SI:
      SI = np.zeros((ndr[0],ndr[1]))
   Bp =  data[idr[2]]
   Bq =  data[idr[1]]
   Bn = -data[idr[0]]

# projection along x-axis
   if ax_set == 0:
      for i in range (ny):
         for j in range (nz):
            plot_data_arrays = stokes_params(Bp[:,i,j], Bq[:,i,j], Bn[:,i,j], rho_ion[:,i,j], Ecrp[:,i,j], Ecre[:,:,i,j], wave_data, ds, ncre, ndr[2])
            if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
               I[j,i] = plot_data_arrays[0]
            if stg.print_PI or stg.print_SI or stg.print_vec:
               Q[j,i], U[j,i] = plot_data_arrays[1:3]
            if stg.print_RM:
               RM[j,i] = plot_data_arrays[3]
            if stg.print_SI:
               SI[j,i] = plot_data_arrays[4]
# projection along y-axis
   elif ax_set == 1:
      for k in range (nz):
         for i in range (nx):
            plot_data_arrays = stokes_params(Bp[i,:,k], Bq[i,:,k], Bn[i,:,k], rho_ion[i,:,k], Ecrp[i,:,j], Ecre[:,i,:,k], wave_data, ds, ncre, ndr[2])
            if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
               I[i,k] = plot_data_arrays[0]
            if stg.print_PI or stg.print_SI or stg.print_vec:
               Q[i,k], U[i,k] = plot_data_arrays[1:3]
            if stg.print_RM:
               RM[i,k] = plot_data_arrays[3]
            if stg.print_SI:
               SI[i,k] = plot_data_arrays[4]
# projection along z-axis
   elif ax_set == 2:
      for i in range (nx):
         for j in range (ny):
            plot_data_arrays = stokes_params(Bp[i,j,:], Bq[i,j,:], Bn[i,j,:], rho_ion[i,j,:], Ecrp[i,j,:], Ecre[:,i,j,:], wave_data, ds, ncre, ndr[2])
            if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
               I[i,j] = plot_data_arrays[0]
            if stg.print_PI or stg.print_SI or stg.print_vec:
               Q[i,j], U[i,j] = plot_data_arrays[1:3]
            if stg.print_RM:
               RM[i,j] = plot_data_arrays[3]
            if stg.print_SI:
               SI[i,j] = plot_data_arrays[4]
   x, y = np.arange(ndr[1]), np.arange(ndr[0])
   figext = [0, ndr[0], 0, ndr[1]]
   time = 0.0
   return I, Q, U, RM, SI, x, y, figext, time

def data_h5(plik,ax_set,wave_data):
# function computing Stokes parameters for simulation data stored in an HDF5 file (.h5)
#plik   - filename.h5
#ax_set     - parameter storing the user's choice of axis along which the
#             projection is executed
#x_ion      - degree of ionization of the ionized gas (stg)
#nu         - observation frequency of synchrotron radiation
#ldb        - observation  wavelength of synchrotron radiation
#cJnu       - a constant scaling synchrotron emissivity.
#ds         - cell size along the line of sight

   I, Q, U, RM, SI = [], [], [], [], []
   h5f = h5.File(plik, 'r')
   time = h5f.attrs['time'][0]
   attrs = h5f['domains']['base'].attrs
   xmin = attrs['x-edge_position'][0]/1000.
   xmax = attrs['x-edge_position'][1]/1000.
   ymin = attrs['y-edge_position'][0]/1000.
   ymax = attrs['y-edge_position'][1]/1000.
   zmin = attrs['z-edge_position'][0]/1000.
   zmax = attrs['z-edge_position'][1]/1000.
   print("Domain dimensions: ", 'xmin, xmax =(', xmin,',', xmax, '), ymin, ymax =(', ymin,',', ymax, '), zmin, zmax =(', zmin,',', zmax,')')
   nxd, nyd, nzd = attrs['n_d'][0:3]
   Ecrp     = np.zeros((nxd,nyd,nzd))
   Ecre     = np.zeros((ncre,nxd,nyd,nzd))
   rho     = np.zeros((nxd,nyd,nzd))
   Bp      = np.zeros((nxd,nyd,nzd))
   Bq      = np.zeros((nxd,nyd,nzd))
   Bn      = np.zeros((nxd,nyd,nzd))
   if ax_set==0:
      bset = ['mag_field_y', 'mag_field_z', 'mag_field_x']
      x = np.linspace(ymin, ymax, nyd)
      y = np.linspace(zmin, zmax, nzd)
      figext = [ymin, ymax, zmin, zmax]
      ds = (xmax - xmin) / nxd * 1000.
      n1, n2, n3 = nyd, nzd, nxd
   if ax_set==1:
      bset = ['mag_field_x', 'mag_field_z', 'mag_field_y']
      zmin = zmin + 1.2
      zmax = zmax - 1.2
      x = np.linspace(xmin, xmax, nxd)
      y = np.linspace(zmin, zmax, nzd)
      figext = [xmin, xmax, zmin, zmax]
      ds = (ymax - ymin) / nyd * 1000.
      n1, n2, n3 = nxd, nzd, nyd
   if ax_set==2:
      bset = ['mag_field_x', 'mag_field_y', 'mag_field_z']
      #print(np.shape(Bp), Bp, np.shape(h5yt.arr['mag_field_x']))
      x = np.linspace(xmin, xmax, nxd)
      y = np.linspace(ymin, ymax, nyd)
      figext = [xmin, xmax, ymin, ymax]
      ds = (zmax - zmin) / nzd * 1000.
      n1, n2, n3 = nxd, nyd, nzd

   print('Reconstructing domain from cg parts')
   grid = h5f['grid_dimensions']
   for ig in range(grid.shape[0]):
      h5g = h5f['data']['grid_0000000'+str(ig).zfill(3)]
      off = h5g.attrs['off']
      ngb = h5g.attrs['n_b']
      n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
      ce  = n_b+off
      rho[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] =  h5g['density'][:,:,:].swapaxes(0,2)
      Bp[ off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] =  h5g[bset[0]][:,:,:].swapaxes(0,2)
      Bq[ off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] =  h5g[bset[1]][:,:,:].swapaxes(0,2)
      Bn[ off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] = -h5g[bset[2]][:,:,:].swapaxes(0,2)
      if stg.mode == 'simple':
         Ecrp[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] =  h5g['cr01'][:,:,:].swapaxes(0,2)
      elif stg.mode == 'spectral':
         for ic in range(ncre):
            Ecre[ic,off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] =  h5g['cren'+str(ic+1).zfill(2)][:,:,:].swapaxes(0,2)

   h5f.close()
   rho_ion = np.zeros_like(rho)
   rho_ion=rho*stg.x_ion
   nx, ny, nz = s = Bp.shape
   I  = np.zeros((n1,n2))
   Q  = np.zeros((n1,n2))
   U  = np.zeros((n1,n2))
   RM = np.zeros((n1,n2))
   if stg.print_SI:
      SI = np.zeros((n1,n2))

   klo, khi = stg.krange(n3,ax_set)
   print('Sendig data to map computation')
   print('')
# projection along x-axis
   if ax_set == 0:
      for i in range (ny):
         sys.stdout.write("\033[F") # Cursor up one line
         print('Layer: ', i, ny)
         for j in range (nz):
            plot_data_arrays = stokes_params(Bp[klo:khi,i,j], Bq[klo:khi,i,j], Bn[klo:khi,i,j], rho_ion[klo:khi,i,j], Ecrp[klo:khi,i,j], Ecre[:,klo:khi,i,j], wave_data, ds, ncre, khi-klo)
            if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
               I[i,j] = plot_data_arrays[0]
            if stg.print_PI or stg.print_SI or stg.print_vec:
               Q[i,j], U[i,j] = plot_data_arrays[1:3]
            if stg.print_RM:
               RM[i,j] = plot_data_arrays[3]
            if stg.print_SI:
               SI[i,j] = plot_data_arrays[4]
      print('B mean, max x = ', np.mean(np.abs(Bn[klo:khi,:,:])), np.amax(np.abs(Bn[klo:khi,:,:])))
# projection along y-axis
   elif ax_set == 1:
      for k in range (nz):
         sys.stdout.write("\033[F") # Cursor up one line
         print('Layer: ', k, nz)
         for i in range (nx):
            plot_data_arrays = stokes_params(Bp[i,klo:khi,k], Bq[i,klo:khi,k], Bn[i,klo:khi,k], rho_ion[i,klo:khi,k], Ecrp[i,klo:khi,k], Ecre[:,i,klo:khi,k], wave_data, ds, ncre, khi-klo)
            if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
               I[i,k] = plot_data_arrays[0]
            if stg.print_PI or stg.print_SI or stg.print_vec:
               Q[i,k], U[i,k] = plot_data_arrays[1:3]
            if stg.print_RM:
               RM[i,k] = plot_data_arrays[3]
            if stg.print_SI:
               SI[i,k] = plot_data_arrays[4]
# # projection along z-axis
   elif ax_set == 2:
      for i in range (nx):
         sys.stdout.write("\033[F") # Cursor up one line
         print('Layer: ', i, nx)
         for j in range (ny):
            plot_data_arrays = stokes_params(Bp[i,j,klo:khi], Bq[i,j,klo:khi], Bn[i,j,klo:khi], rho_ion[i,j,klo:khi], Ecrp[i,j,klo:khi], Ecre[:,i,j,klo:khi], wave_data, ds, ncre, khi-klo)
            if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
               I[i,j] = plot_data_arrays[0]
            if stg.print_PI or stg.print_SI or stg.print_vec:
               Q[i,j], U[i,j] = plot_data_arrays[1:3]
            if stg.print_RM:
               RM[i,j] = plot_data_arrays[3]
            if stg.print_SI:
               SI[i,j] = plot_data_arrays[4]
      print('B mean, max z = ', np.mean(np.abs(Bn[:,:,klo:khi])), np.amax(np.abs(Bn[:,:,klo:khi])))
   return I, Q, U, RM, SI, x, y, figext, time
