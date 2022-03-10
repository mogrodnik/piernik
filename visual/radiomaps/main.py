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
from draw_map import draw_map, dump_data
from electrons import initialize_crspectrum_tools
from profiles  import plot_profile


save_data = False
file_name = ""
from_file = False
ax = 'x'
ax_set= 0
convol= False
new_cp= False
handSI, handPI, handRM, handTP, handVC, pretVC = False, False, False, False, False, False
lbd_set  = [-1.0]
lbd2_set = -1.0   # DEPRECATED
nu_set   = [-1.0]
nu2_set  = -1.0   # DEPRECATED
yt_imres = 0
yt_depth = "max"

def cli_params(argv):
   # The function serves for reading and interpretation of the comand line input parameters
   try:

      opts,args=getopt.getopt(argv,"adhf:i:k:l:m:n:R:cpPrtSuvxyz",["dump","help","file=","convolve","contours","lin=","log","spectral","suffix=","yt=","tz=","pz=","iz=","rz=","pr=","vp=","cp="])
      #print (opts,"op",args,"arg")
   except getopt.GetoptError:
      print("Error: unknown parameter")
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h", "--help"):
         print("-f [-file] filename.h5 generates maps from an hdf5 file. Running the script without the -f parameter generates a map based on analytical data (it is currently broken)  \n -l sets the wavelengths \n -n sets the frequency \n -c convolves the resulting data with a 2D Gauss function representing the angular characteristic of the radiotelescope beam \n -x generates projection parallel to x-axis (default option)  \n -y along y-axis \n -z along z-axis \n -i pair1,pair2,.. generates the map of Spectral Index (SI), requires additional input: list indexes of of pairs of wavelengths or frequencies provided in -l or -n (e.g., 01,02)  \n -p the map of Polarized Intensity (PI) \n -t produces the map of Total Power (TP) \n -r produces the map of Rotation measue (RM) \n -v to add vectors and \n -u not to add vectors \n --log to drow the map in logarythmic scale \n -S --spectral to generate synchrotron radiation maps using electron number and energy density data \n -R <integer> reads and plots only one refinement level\n -P to additionally plot profiles of S/P/T intensity \n --log to drow the map in logarythmic scale \n --lin {SI,TP,PI,RM} to override logarythmic scaling for chosen components \n --yt RESX,DEPTH or RESX,DEPTHL,DEPTHH \t use yt package for reading data and construct image of resolution RESX:(RESX / <aspect ratio>); type 'max' for maximal practical resolution, with depth working as -DEPTH:DEPTH with one DEPTH argument ('max' for maximal depth) or DEPTHL:DEPTHH for layer with 2 DEPTH arguments. Slower method, but works well with all levels of AMR data\n --{tz|pz|rz|iz} vmin,vmax \t to set plot limits for the given field (TP, PI, SI or RM, respectively). \n \t \t \t \t For TP and PI it must match the number of wavelengths provided, for SI the number of pairs, for RM one set of limits (wavelength independent quantity) \n --pr width,height to set the absolute limits of plotted physical space \n --vp DVEC to set vector coverage \n --cp nbeam1,sigma1:nbeam2,sigma2:... alter convolution parameters \n \n --contours draw map contours (Default: 10 levels)\n --suffix <text> to have output files with given suffix before extension\n --dump saves data to ASCII file (e.g., to plot it using other script)")
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
         lbd_set = [ float(item) for item in arg.split(",")]

      elif opt == "-n":
         global nu_set
         nu_set = [ float(item) for item in arg.split(",")]

      elif opt == "-i":
         global handSI
         handSI = True
         stg.SI_set = [ (int(item[0]), int(item[1])) for item in arg.split(",")]

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

      elif opt in ("--lin"):
         aux = arg.split(",")
         for item in aux: stg.apply_linear.append(item.upper())

      elif opt == "-P":
         stg.print_prof = True

      elif opt in ("-c", "--convolve"):
         global convol
         convol = True

      elif opt in ("-R"):
         stg.lvl_only = int(arg)
         stg.one_level = True

      elif opt in ("-S", "--spectral"):
         stg.spectral_mode = True
         stg.mode = "spectral"   # DEPRECATED
         print("Proceeding with spectral method")

      elif opt in ("--yt"):
         global yt_imres, yt_depth
         aux        = arg.split(",")
         if (len(aux) < 2): sys.exit("Problem in reading parameters RESX,DEPTH from --yt arguments; got: %s" %str(aux))
         yt_imres = aux[0] if aux[0] == "max" else int(aux[0])
         try:
            if len(aux[0:-1]) == 2:
               yt_depth = [ float(aux[1]), float(aux[2])]
            else:
               yt_depth = [-float(aux[1]), float(aux[1])]
            print("Modules data_h5_yt and yt will be imported, range %10.1f:%10.1f, resolution: %s" %(yt_depth[0], yt_depth[1], str(yt_imres)))
         except:
            yt_depth   = "max"
         stg.use_yt = True

      elif opt in ("--tz"):
         stg.Tvmn_user, stg.Tvmx_user = [], []
         for auxi in arg.split(":"):
            stg.Tvmn_user.append(float(auxi.split(",")[0]))
            stg.Tvmx_user.append(float(auxi.split(",")[1]))

      elif opt in ("--pz"):
         stg.Pvmn_user, stg.Pvmx_user = [], []
         for auxi in arg.split(":"):
            stg.Pvmn_user.append(float(auxi.split(",")[0]))
            stg.Pvmx_user.append(float(auxi.split(",")[1]))

      elif opt in ("--rz"):   # RM is independent on wavelength
         stg.Rvmn_user, stg.Rvmx_user = -100., 100.
         stg.Rvmn_user = float(arg.split(",")[0])
         stg.Rvmx_user = float(arg.split(",")[1])

      elif opt in ("--iz"):
         stg.Svmn_user, stg.Svmx_user = [], []
         for auxi in arg.split(":"):
            stg.Svmn_user.append(float(auxi.split(",")[0]))
            stg.Svmx_user.append(float(auxi.split(",")[1]))

      elif opt in ("--pr"):
         aux = arg.split(",")
         stg.px_user, stg.pz_user = float(aux[0]), float(aux[1])

      elif opt in ("--vp"):
         stg.dokv_user = int(arg)

      elif opt in ("--cp"):
         global new_cp
         stg.nbeam = [] ; stg.sigma = []
         new_cp = True
         for auxi in arg.split(":"):
            stg.nbeam.append(int(auxi.split(",")[0]))
            stg.sigma.append(float(auxi.split(",")[1]))

      elif opt in ("--contours"):
         stg.print_cnt = True

      elif opt in ("--dump"):
         global save_data
         save_data = True

      elif opt in ("--suffix"):
         stg.suffix = str(arg)

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

if (len(nu_set) > len(lbd_set)):
   lbd_set = [-1. for i in range(len(nu_set))]
else:
   nu_set =  [-1. for i in range(len(lbd_set))]

nu  = [ 0. for i in range(len(lbd_set)) ]
lbd = [ 0. for i in range(len(lbd_set)) ]

stg.N_nulbd = len(lbd_set)

#If user-set limits are to applied, make sure not to accept improper number of limits: stop, if sizes do not match
if ( stg.N_nulbd > 1 ):
   if (stg.print_TP and (stg.Tvmn_user != False and stg.Tvmx_user != False  )):
      if (stg.N_nulbd != np.shape(stg.Tvmn_user)[0] or stg.N_nulbd != np.shape(stg.Tvmx_user)[0]): sys.exit("(TP range) when setting ranges, they must be set for all wavelengths (--tz option).")
   if (stg.print_PI and (stg.Pvmn_user != False and stg.Pvmx_user != False  )):
      if (stg.N_nulbd != np.shape(stg.Pvmn_user)[0] or stg.N_nulbd != np.shape(stg.Pvmx_user)[0]): sys.exit("(PI range) when setting ranges, they must be set for all wavelengths (--pz option).")
   if (stg.print_SI and (stg.Svmn_user != False and stg.Svmx_user != False  )):
      if (np.shape(stg.SI_set)[0] != np.shape(stg.Svmn_user)[0] or np.shape(stg.SI_set)[0] != np.shape(stg.Svmx_user)[0]): sys.exit("(SI range) when setting ranges, they must be set for all pairs of wavelengths (--iz option).")
   if (new_cp == True and (len(stg.nbeam) != stg.N_nulbd or  len(stg.sigma) != stg.N_nulbd)): sys.exit("(convolution parameters) when altering nbeam and sigma, they must be set for all wavelengths (got %i, should be %i, check --cp option)." %(len(stg.nbeam), stg.N_nulbd ))

for i in range(stg.N_nulbd):
   nu[i],   lbd[i]   = stg.set_nuandlbd(nu_set[i],  lbd_set[i], i)

print("Wavelengths (m):  ", lbd)
print("Frequencies (MHz):", [ round(item / 1.e6, 2) for item in nu])

wave_data = [nu, lbd]
if (stg.spectral_mode):
   read_CRESP_params(file_name)
   initialize_crspectrum_tools(stg.ncre, stg.p_min_fix, stg.p_max_fix, nu)

if not from_file:
   # Analytical data might be used,for testing of the maping routines, to generate 3D arrays of CR energy density, gas density and magnetic fild.
   plot_data_arrays = data_a(ax_set, wave_data)
else:
   # To visualize simulation results the 3D arrays are read from hdf5 files
   if (stg.use_yt):  # yt-projects' yt module can be used if it is available (useful for AMR data, but slow)
      from yt_read_data import data_h5_yt
      plot_data_arrays = data_h5_yt(file_name, ax_set, wave_data, yt_imres, yt_depth)
   else:
      plot_data_arrays = data_h5(file_name, ax_set, wave_data)

I, Q, U, RM, SI, x, y, figext, time = plot_data_arrays

for i_nl in range(stg.N_nulbd):
   if convol:
      # Generation of the Gaussian profile of the radiotelescop beam
      sigma = stg.sigma[i_nl] if new_cp else stg.sigma[0]
      nbeam = stg.nbeam[i_nl] if new_cp else stg.nbeam[0]
      beam = gauss_beam(nbeam, sigma)
      # We convolve the resulting Stokes parameters tables with the beam function
      if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
         I[i_nl]  = data_beam_convolve(I[i_nl],  beam, nbeam)
      if stg.print_PI or stg.print_SI or stg.print_vec:
         Q[i_nl]  = data_beam_convolve(Q[i_nl],  beam, nbeam)
         U[i_nl]  = data_beam_convolve(U[i_nl],  beam, nbeam)
      if stg.print_RM:
         RM = data_beam_convolve(RM, beam, nbeam)

   if stg.print_PI or stg.print_vec:
      PI, PI_obs = polarized(I[i_nl], Q[i_nl], U[i_nl])
   if stg.print_vec:
      wp, wq = vector_direction(PI_obs, Q[i_nl], U[i_nl])
      # Creation of a mesh of points to place vectors
      X, Y = np.meshgrid(x,y)
      vecs = stg.dokvec(wp, wq, X, Y, ax_set, from_file)
   else:
      vecs = []

   if from_file and convol:
      if stg.print_TP:
         I[i_nl]  = I[i_nl]**stg.normalise_exponent_PI
      if stg.print_PI:
         PI = PI**stg.normalise_exponent_PI

   print("From_file: ", from_file)

   etyfil = stg.etyfil(ax,file_name,nu[i_nl])
   attr = [time,lbd[i_nl], lbd[0]]
   etyyt = etyspec = ""

   if (stg.use_yt):
      if (yt_depth != "max"):
         etyyt = "_yt_res"+str(yt_imres)+"range"+str(round(yt_depth[0]))+":"+str(round(yt_depth[1]))
      else:
         etyyt = "_yt_res"+str(yt_imres)+"range"+yt_depth.upper()
   if (stg.spectral_mode): etyspec = "_spectral"

   etyfil = etyfil + etyyt + etyspec + stg.suffix  # Full name + optional suffix

   if stg.print_TP:
      # Drawing Total Power (TP) map
      draw_map( I[i_nl].T, vecs, figext, ax_set, attr, etyfil, 'TP', from_file, i_nl)
      if (stg.print_prof): plot_profile(I[i_nl].T, figext, ax_set, etyfil, 'TP', attr)
      if (save_data): dump_data(I[i_nl].T, lbd[i_nl], nu[i_nl], figext, time, etyfil, 'TP')
   if stg.print_PI:
      # Drawing Polarized Intensity (PI) map
      draw_map(PI.T, vecs, figext, ax_set, attr, etyfil, 'PI', from_file, i_nl)
      if (stg.print_prof): plot_profile(PI.T, figext, ax_set, etyfil, 'PI', attr)
      if (save_data):
         dump_data(PI.T, lbd[i_nl], nu[i_nl], figext, time, etyfil, 'PI')
         dump_data(Q[i_nl].T, lbd[i_nl], nu[i_nl], figext, time, etyfil, 'Q')
         dump_data(U[i_nl].T, lbd[i_nl], nu[i_nl], figext, time, etyfil, 'U')
for i_pair in stg.SI_set:
   attr = [time,lbd[i_pair[0]],lbd[i_pair[1]]]
   if stg.print_SI:
      etyfil = stg.etyfil(ax,file_name,nu[i_pair[0]]) + "vs"+str(int(nu[i_pair[1]]/stg.MHz)).zfill(3)+'MHz' # Provides proper wavelengths / frequencies to filename
      etyfil = etyfil + etyyt + etyspec + stg.suffix # Full name + optional suffix

      draw_map(SI[stg.SI_set.index(i_pair)].T, vecs, figext, ax_set, attr, etyfil, 'SI', from_file, stg.SI_set.index(i_pair))
      if (stg.print_prof): plot_profile(SI[stg.SI_set.index(i_pair)].T, figext, ax_set, etyfil, 'SI', attr)
      if (save_data): dump_data(SI[stg.SI_set.index(i_pair)].T, [lbd[i_pair[0]], lbd[i_pair[1]]], [nu[i_pair[0]], nu[i_pair[1]]], figext, time, etyfil, 'SI')

if stg.print_RM:
   if np.max(RM) != 1.0 or np.min(RM) != 0.0:
      # We draw the Faraday rotation - Rotation measue (RM) only if RM != 0
      draw_map(RM.T/stg.norm('RM'), vecs, figext, ax_set, attr, etyfil, 'RM', from_file, i_nl)
      if (save_data): dump_data(RM.T/stg.norm('RM'), lbd[i_nl], nu[i_nl], figext, time, etyfil, 'RM')
   else:
      print('RM = 0; I do not create the map.')
