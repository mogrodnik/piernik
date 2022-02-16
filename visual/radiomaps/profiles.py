from math import *
import numpy as np
import h5py as h5
import sys
from analytical_data import data_raw
import settings as stg
import matplotlib.pyplot as plt

#przeznaczeniem tego modulu jest przygotowanie usrednionych w plaszczyznie X profili emisji i innych obserwabli
#w zalozeniu profil rzutuje na os X emisje w funkcji odleglosci od dysku (z lub h). Os Y prezentuje wartosc obserwabli.
#nwboxes - liczba zakresow
#wrange  - calkowity zakres szerokosci, ktory zostanie podzielony na nwboxes kawalkow i w ktorym obserwabla bedzie usredniana
#hrange  - wysokosc w kpc, rzutuje na zakres w osi X profili
#irange  - zakres wartosci obserwabli

fsize   = 14

wranged = [-7.37, 7.37] # ranges from disk center in kpc (as in Mulcachy et al 2018, fig. 14)
hranged = [-3.22, 3.22] # ranges from disk plane  in kpc (as in Mulcachy et al 2018, fig. 14)
iranged = [0.1,  30.]   # range in observable            (as in Mulcachy et al 2018, fig. 14)
nwboxesd= 9             # number of boxes (as in Mulcachy et al 2018, fig. 14)
kpcasecd= 46./1000.     # 1'' equivalent to 46 pc as default
labelnam= {"SI":"Spectral index", "TP":"Total intensity", "PI":"Polarized intensity"}
ylimsdef= {"SI":(-3.,0.5), "TP":(1.e-2, 5.e0), "PI":(1.e-3, 5.e-1)}
plt_grds= [True, True]
use_def_ranges = True
plot_errorbars = True
use_def_ylims  = True
plot_one_fig = True
save_data    = True  # in case the profiles are to be plotted in a different way
return_axes  = False # return all plotted objects as list of axes
inches_ax = 2.
aspect_ax = 6.
dpi_plot  = 150.

# main plotting function: processess data and produces plot files
def plot_profile(data, figext_tot, ax_set, ety_file, label, attributes, **kwargs):

# establish 'width' (ordinate) axis
   if (ax_set == 2):
      print("\033[93m Axis set is 'z', omitting profile plotting. \033[0m")
      return # do nothing - no profiles to plot, if disk is shown face-on
   elif (ax_set == 0):
      ax_w = "x"
   elif (ax_set == 1):
      ax_w = "y"
# Unpack time and wavelength info
   time = attributes[0]
   lbd1 = attributes[1]
   lbd2 = attributes[2]

# get shape of the data
   data_shape = np.shape(data) # expecting 2-dim data, WARNING - obtained data is reversed
   data[data == -np.inf] = 0.     # remove infs

# assume larger data length in width
   datawshape = max(data_shape)
   datahshape = min(data_shape)

   if stg.print_log:
      data = np.log10(data)

# get keyword arguments
   wrange = kwargs.get("w", wranged)
   hrange = kwargs.get("h", hranged)
   irange = kwargs.get("z", iranged)
   nwboxes= kwargs.get("n", nwboxesd)
   kpcasec= kwargs.get("s", kpcasecd)

# Get physical size of the data
   figw_tot = figext_tot[1] - figext_tot[0]
   figh_tot = figext_tot[3] - figext_tot[2]
   figw     = [figext_tot[0], figext_tot[1]]
   figh     = [figext_tot[2], figext_tot[3]]

   if (not use_def_ranges):
      wrange = [figext_tot[0] , figext_tot[1] ]
      hrange = [figext_tot[0] , figext_tot[1] ]

# wcell and hcell are array cell sizes:
   wcell = figw_tot / data_shape[0]
   hcell = figh_tot / data_shape[1]

# Establish width ranges for the sections
   wwidth   = (wrange[1] - wrange[0]) / nwboxes
   wxrng    = [wrange[0]]
   for i in range(nwboxes):
      wxrng.append(wxrng[0] + (i+1) * wwidth)

   mfig, axs = plt.subplots(nrows=nwboxes, ncols=1, figsize=(inches_ax * aspect_ax, inches_ax * nwboxes), dpi=150)
   means_save = []
   stds_save  = []
# prepare and save profile(s)
   for iprof in range(nwboxes):
      ind_h, h_adj_edges = get_encompassed_range(datahshape, figh, hrange)

      mid_coord = (wxrng[iprof+1] + wxrng[iprof])*0.5
      ihb = max(ind_h[0]-1, 0)   ;  ihe = min(ind_h[-1]+1, datahshape-1)

      means = []
      stds  = []
      for i in range(ihb, ihe, 1):
         avg, std = avg_std_section(data[i][:], datawshape, [figw[0], figw[1]], [wxrng[iprof], wxrng[iprof+1]])
         means.append(avg)
         stds.append(std)

      if (stg.print_log and (label == "TP" or label == "PI")):
         means = np.power(10., means)
         scaletype = "log"
      else:
         scaletype = "linear"

      hdata = np.linspace(hrange[0], hrange[1], ihe-ihb)
      plot_labelsxy = [ not return_axes, not plot_one_fig ]

      xlab = "Distance from the disk plane (kpc)"
      ylab = "Averaged %s at %5.1f kpc" %(labelnam[label], mid_coord)
      title = "Time = %10.2f Myr" %time   # WARNING assuming units
      if (return_axes): all_axes = []

      if (not plot_one_fig):
         fig_one, ax = plt.subplots(figsize=(7,5), dpi=150)
         ax = plot_one_profile(ax, hdata, means, stds, xlab, ylab, title, label, scaletype, ylims, text_overplot, plot_labelsxy, plot_title=False, ticks_visiblexy = [ (return_axes or iprof != nwboxes -1), False ])

         if (return_axes):
            all_axes.append(ax)
         else:
            plt.savefig(label + "_profile_" + str(iprof+1) + ety_file + ".png")
            plt.savefig(label + "_profile_" + str(iprof+1) + ety_file + ".pdf")

         ax.clear()
         plt.close(fig_one)

      axm = axs[iprof]
      xlab = "z (kpc) at %s = %5.1f kpc" %(ax_w, mid_coord)
      axm = plot_one_profile(axm, hdata, means, stds, xlab, ylab, title, label, scaletype, ylims, text_overplot, plot_labelsxy, plot_title=False, ticks_visiblexy = [ iprof != nwboxes -1, False ])
      means_save.append(means)
      stds_save.append(stds)

   mfig.suptitle(title, fontsize = fsize*1.35)
   mfig.supylabel("%s (averaged over %s) at $\lambda=$ %8.2f m" %(labelnam[label], ax_w, lbd1), fontsize = fsize*1.5)
   #plt.tight_layout()
   plt.tight_layout(rect=[0., 0., 1., 0.99]) # tight_layout may clip suptitle

   if (not return_axes):
      plt.savefig(label + "_profiles" + ety_file + ".png")
      plt.savefig(label + "_profiles" + ety_file + ".pdf")
      plt.close(mfig)

   if (save_data):
      fP = open(label+ety_file+"_profile_data.dat", "w+")
      fP.write("#     dist\t")
      for wi in range(nwboxes):
         fP.write(" %8s(p %2i)\t  %8s(p %2i)\t" %("mean",wi,"std",wi))
      fP.write("\n")
      for i in range(len(hdata)):
         fP.write("%12.5f\t" %hdata[i])
         for wi in range(nwboxes):
            fP.write(" %14.8e\t  %14.8e\t " %(means_save[wi][i], stds_save[wi][i]))
         fP.write("\n") # should add new line
      fP.close()

   if (return_axes): return mfig, all_axes

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
   avg_section = update_average(avg_section, num_avg, vec_data[max(indices_in[-1], 0)], incrs[1])
   return avg_section

def avg_std_section(vec_data, sec_dim, sec_lims, sec_where):
   twidth = (sec_lims[1] - sec_lims[0])
   dwidth = twidth / sec_dim              # WARNING Assumes uniform distances
   indices_in, incrs = get_encompassed_range(sec_dim, sec_lims, sec_where)
   indices_in.insert(0, max(indices_in[0],  0))
   indices_in.insert(0, min(indices_in[-1], sec_dim))

   avg_vec = np.mean(vec_data[indices_in])
   std_vec = np.std( vec_data[indices_in])

   return avg_vec, std_vec

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

def plot_one_profile(ax, xdata, ydata, yerrdata, xlabel, ylabel, title, label, scaletype, ylims, coords_text, plot_labels, plot_title, ticks_visiblexy):
   if plot_title:
      ax.set_title(title, fontsize=fsize)

   ax.set_yscale(scaletype)
   if use_def_ylims: # if true, using global 'ylimsdef' via passed 'label' variable
      ax.set_ylim(ylimsdef[label][0], ylimsdef[label][1])
   else:
      ax.set_ylim(ylims)

   ax.grid(plt_grds[0], 'major', 'x', ls='--', lw=.5, c='k', alpha=.3)
   ax.grid(plt_grds[1], 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)

   if (ticks_visiblexy[0]): ax.xaxis.set_ticklabels([])
   if (ticks_visiblexy[1]): ax.yaxis.set_ticklabels([])

   if (plot_errorbars):
      ax.errorbar(xdata, ydata, yerr=yerrdata, color="xkcd:red",  marker = "+")
   else:
      ax.scatter(xdata, ydata, color="xkcd:red",  marker = "+")

   ax.plot(   xdata, ydata, color="xkcd:blue", linewidth=0.75)

   if plot_labels[0]: ax.set_xlabel(xlabel, fontsize=fsize)
   if plot_labels[1]: ax.set_ylabel(ylabel, fontsize=fsize)

   ax.tick_params(labelsize = fsize * 0.95)
   return ax
