from math import *
import numpy as np
import h5py as h5
import sys
from analytical_data import data_raw
import settings as stg
import matplotlib.pyplot as plt

#przeznaczeniem tego modulu jest przygotowanie usrednionych w plaszczyznie X profili emisji i innych obserwabli
#w zalozeniu profil rzutuje na os X emisje w funkcji odleglosci od dysku (z lub h). Os Y prezentuje wartosc obserwabli.
#nxboxes - liczba zakresow
#wrange  - calkowity zakres szerokosci, ktory zostanie podzielony na nxboxes kawalkow i w ktorym obserwabla bedzie usredniana
#hrange  - wysokosc w kpc, rzutuje na zakres w osi X profili
#irange  - zakres wartosci obserwabli

fsize   = 14

wranged = [-7.37, 7.37] # ranges from disk center in kpc (as in Mulcachy et al 2018, fig. 14)
hranged = [-3.22, 3.22] # ranges from disk plane  in kpc (as in Mulcachy et al 2018, fig. 14)
iranged = [0.1,  30.]   # range in observable            (as in Mulcachy et al 2018, fig. 14)
nxboxesd= 9             # number of boxes (as in Mulcachy et al 2018, fig. 14)
kpcasecd= 46./1000.     # 1'' equivalent to 46 pc as default
labelnam= {"SI":"spectral index", "TP":"total intensity", "PI":"polarized intensity"}

def plot_profile(data, figext_tot, ax_set, ety_file, label, **kwargs):

   if (ax_set == 2):
      print("\033[93m Axis set is 'z', omitting profile plotting. \033[0m")
      return # do nothing - no profiles to plot, if disk is shown face-on

   data_shape = np.shape(data) # expecting 2-dim data, WARNING - obtained data is reversed
   data[data == -np.inf] = 0.     # remove infs

   datawshape = max(data_shape)
   datahshape = min(data_shape)

   if stg.print_log:
      data = np.log10(data)

   wrange = kwargs.get("w", wranged)
   hrange = kwargs.get("h", hranged)
   irange = kwargs.get("z", iranged)
   nxboxes= kwargs.get("n", nxboxesd)
   kpcasec= kwargs.get("s", kpcasecd)

   figw_tot = figext_tot[1] - figext_tot[0]
   figh_tot = figext_tot[3] - figext_tot[2]
   figw = [figext_tot[0], figext_tot[1]]
   figh = [figext_tot[2], figext_tot[3]]

   wcell = figw_tot / data_shape[0]
   hcell = figh_tot / data_shape[1]

   if (nxboxes is int(nxboxes**0.5)**2):
      n_r = int(nxboxes**0.5)
      n_c = int(nxboxes**0.5)
   else:
      n_r = int(nxboxes**0.5)
      n_c = math.ceil(float(nxboxes) / float(n_r))

   wwidth   = (wrange[1] - wrange[0]) / nxboxes
   wxrng    = [wrange[0]]
   for i in range(nxboxes):
      wxrng.append(wxrng[0] + (i+1) * wwidth)

   for iprof in range(nxboxes):
      fig = plt.figure(figsize=(7,5), dpi=150)

      ind_h, h_adj_edges = get_encompassed_range(datahshape, figh, hrange)

      mid_coord = (wxrng[iprof+1] + wxrng[iprof])*0.5
      ihb = max(ind_h[0]-1, 0)   ;  ihe = min(ind_h[-1]+1, datahshape)

      means = []
      for i in range(ihb, ihe, 1):
         means.append(avg_vec_in_range(data[i][:], datawshape, [figw[0], figw[1]], [wxrng[iprof], wxrng[iprof+1]] ) )

      hdata = np.linspace(hrange[0], hrange[1], ihe-ihb)

      plt.axes().tick_params(axis='y', which='minor', left=True)
      plt.xlabel("Distance from the disk plane (kpc)", fontsize=fsize)
      plt.ylabel("Averaged %s at %5.1f kpc" %(labelnam[label], mid_coord), fontsize=fsize)

      if (stg.print_log and (label == "TP" or label == "PI")):
         means = np.power(10., means)
         plt.yscale("log")
      plt.scatter(hdata, means, color="xkcd:red", marker = "X")
      plt.plot(   hdata, means, color="xkcd:blue", linewidth=0.75)

      plt.savefig(label + "_profile_" + str(iprof+1) + ety_file + ".png")
      plt.savefig(label + "_profile_" + str(iprof+1) + ety_file + ".pdf")

def avg_vec_in_range(vec_data, sec_dim, sec_lims, sec_where):
   twidth = (sec_lims[1] - sec_lims[0])
   dwidth = twidth / sec_dim              # WARNING Assumes uniform distances
   indices_in, incrs = get_encompassed_range(sec_dim, sec_lims, sec_where)

   avg_whole = np.mean(vec_data[indices_in])
   i_beg = max(indices_in[0],  0)
   i_end = min(indices_in[-1], sec_dim)

   avg_section = avg_whole
   num_avg     = float(len(indices_in))
   avg_section = update_average(avg_section, num_avg, vec_data[max(indices_in[0]-1, 0)], incrs[0])
   num_avg     = num_avg + incrs[0]
   avg_section = update_average(avg_section, num_avg, vec_data[max(indices_in[-1]+1, 0)], incrs[1])
   return avg_section

def get_encompassed_range(sec_dim, sec_lims, sec_where):
   indices_whole = []   ;     i = 0
   bnd_factors   = [0., 0.]
   dwidth = (sec_lims[1] - sec_lims[0]) / sec_dim
   w_edges= [ (sec_lims[0] +  dwidth * i) for i in range(sec_dim + 1)]

   for item in w_edges:
      if (item > sec_where[0] and item < sec_where[1]):
         indices_whole.append(i)
      i = i+1

   bnd_factors[0] = abs(sec_where[0] - w_edges[indices_whole[0]]) / dwidth # gives remaining fraction of cell, checked
   bnd_factors[1] = abs(sec_where[1] - w_edges[indices_whole[-1]]) / dwidth # gives remaining fraction of cell, checked
   return indices_whole, bnd_factors

def update_average(avg_val, avg_num, avg_update, incr):
   denom_new = avg_num + incr
   avg_val_new = avg_val + (avg_update - avg_val) / (avg_num + incr) # This allows to update cell-average by factor of cell
   return avg_val_new
