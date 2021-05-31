#!/usr/bin/env python
from math import *
import numpy as np
import pylab as py
import scipy as sc
import sys
import os
import getopt
import h5py as h5
import settings as stg

from stokes import vector_direction, polarized
from read_data import *
from convolution import data_beam_convolve
from gauss_beam import gauss_beam
from draw_map import draw_map
from electrons import initialize_crspectrum_tools

file_name = ""
from_file = False
ax = 'x'
ax_set= 0
convol= False
handSI, handPI, handRM, handTP, handVC, pretVC = False, False, False, False, False, False
lbd_set  = -1.0
lbd2_set = -1.0
nu_set   = -1.0
nu2_set  = -1.0

def cli_params(argv):
   # The function serves for reading and interpretation of the comand line input parameters
   try:
      opts,args=getopt.getopt(argv,"adhf:ik:l:m:n:R:prtsuvxyz",["help","file","convolve","log"])
      #print opts,"op",args,"arg"
   except getopt.GetoptError:
      print("Error: unknown parameter")
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h", "--help"):
         print("-f [-file] filename.h5 generates maps from an hdf5 file. Running the script without the -f parameter generates a map based on analytical data (it is currently broken)  \n -l sets the wavelengths \n -k sets the second wavelength for sectral index maps. \n -n sets the frequency \n -m sets the second frequency for spectral index maps \n -s convolves the resulting data with a 2D Gauss function representing the angular characteristic of the radiotelescope beam \n -x generates projection parallel to x-axis (default option)  \n -y along y-axis \n -z along z-axis \n -i generates the map of Spectral Index (SI) \n -p the map of Polarized Intensity (PI) \n -t produces the map of Total Power (TP) \n -v to add vectors and \n -u not to add vectors \n --log to drow the map in logarythmic scale \n -R integer reads and plots only one refinement level")
         sys.exit()
      elif opt == '-x':
         global ax_set, ax
         ax_set = 0
         ax = 'x'

      elif opt == '-y':
         ax_set = 1
         ax = 'y'

      elif opt == '-z':
         ax_set = 2
         ax = 'z'

      elif opt == "-l":
         global lbd_set
         lbd_set = float(arg)

      elif opt == "-k":
         global lbd2_set
         lbd2_set = float(arg)

      elif opt == "-m":
         global nu2_set
         nu2_set = float(arg)

      elif opt == "-n":
         global nu_set
         nu_set = float(arg)

      elif opt == "-i":
         global handSI
         handSI = True

      elif opt == "-p":
         global handPI
         handPI = True

      elif opt == "-r":
         global handRM
         handRM = True

      elif opt == "-t":
         global handTP
         handTP = True

      elif opt == "-v":
         global handVC, pretVC
         pretVC = True
         handVC = True

      elif opt == "-u":
         pretVC = True
         handVC = False

      elif opt == "-d":
         stg.print_df = True

      elif opt in ("--log"):
         stg.print_log = True

      elif opt in ("-c", "--convolve"):
         global convol
         convol=True

      elif opt in ("-R"):
         stg.lvl_only = int(arg)
         stg.one_level = True

      elif opt in ("-f", "--file"):
         global file_name
         file_name = str(arg)
         try:
            with open(file_name):
               global from_file
               from_file=True
         except IOError:
            print("No such file")
            sys.exit()

cli_params(sys.argv[1:])

if handPI or handRM or handSI or handTP:
   stg.print_PI, stg.print_RM, stg.print_SI, stg.print_TP = handPI, handRM, handSI, handTP
if pretVC:
   stg.print_vec = handVC

nu,   lbd   = stg.set_nuandlbd(nu_set,  lbd_set, 1)
nu_2, lbd_2 = stg.set_nuandlbd(nu2_set, lbd2_set,2)
wave_data = [nu, lbd, nu_2, lbd_2]
#sys.exit()
initialize_crspectrum_tools(ncre) # TODO FIXME probably should be initialized ONLY IF stg.mode == "spectral"

if not from_file:
   # Analytical data might be used,for testing of the maping routines, to generate 3D arrays of CR energy density, gas density and magnetic fild.
   plot_data_arrays = data_a(ax_set, wave_data)
else:
   # To visualize simulation results the 3D arrays are read from hdf5 files
   plot_data_arrays = data_h5(file_name, ax_set, wave_data)

I, Q, U, RM, SI, x, y, figext, time = plot_data_arrays

if convol:
   # Generation of the Gaussian profile of the radiotelescop beam
   beam = gauss_beam(nbeam, sigma)
   # We convolve the resulting Stokes parameters tables with the beam function
   if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
      I  = data_beam_convolve(I,  beam, nbeam)
   if stg.print_PI or stg.print_SI or stg.print_vec:
      Q  = data_beam_convolve(Q,  beam, nbeam)
      U  = data_beam_convolve(U,  beam, nbeam)
   if stg.print_RM:
      RM = data_beam_convolve(RM, beam, nbeam)

if stg.print_PI or stg.print_vec:
   PI, PI_obs = polarized(I, Q, U)
if stg.print_vec:
   wp, wq = vector_direction(PI_obs, Q, U)
   # Creation of a mesh of points to place vectors
   X, Y = np.meshgrid(x,y)
   vecs = stg.dokvec(wp, wq, X, Y, ax_set, from_file)
else:
   vecs = []

if from_file and convol:
   if stg.print_TP:
      I  = I**stg.normalise_exponent_PI
   if stg.print_PI:
      PI = PI**stg.normalise_exponent_PI

print("From_file: ", from_file)

etyfil = stg.etyfil(ax,file_name,nu)
attr = [time,lbd,lbd_2]
if stg.print_TP:
   # Drawing Total Power (TP) map
   draw_map( I.T, vecs, figext, ax_set, attr, etyfil, 'TP', from_file)
if stg.print_PI:
   # Drawing Polarized Intensity (PI) map
   draw_map(PI.T, vecs, figext, ax_set, attr, etyfil, 'PI', from_file)
if stg.print_SI:
   draw_map(SI.T, vecs, figext, ax_set, attr, etyfil, 'SI', from_file)
if stg.print_RM:
   if np.max(RM) != 1.0 or np.min(RM) != 0.0:
      # We draw the Faraday rotation - Rotation measue (RM) only if RM != 0
      draw_map(RM.T/stg.norm('RM'), vecs, figext, ax_set, attr, etyfil, 'RM', from_file)
   else:
      print('RM = 0; I do not create the map.')

#py.show()
