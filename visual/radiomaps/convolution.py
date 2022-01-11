from math import *
import numpy as np
from scipy import signal

def data_beam_convolve(data,beam,nbeam):
# function convolving the synchrotron image of simulation data with the radiotelescope beam.
# data - array of input data containg the raw synchrotron image of simulation data
# beam - array of radiotelescope beam
# nbeam - nbeam x nbeam is the grid size of radiotelescope beam

   d_shape=np.shape(data)
   d_len=len(d_shape)
   npm=d_shape[0]
   nqm=d_shape[1]

   data_ext=np.zeros((npm+2*nbeam,nqm+2*nbeam))
   data_ext[nbeam:npm+nbeam,nbeam:nqm+nbeam] = data

   image_ext=signal.convolve(data,beam)
   image =image_ext[int(nbeam/2):int(npm+nbeam/2),int(nbeam/2):int(nqm+nbeam/2)]

   return image
