import pylab as py
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.gridspec as gridspec
import settings as stg
import os
from math import ceil

fsize = 24

def figsz(axis):
   if axis == 0:
      fs=(11.5,7.5)
   if axis == 1:
      fs=(11.5,5.5)
   if axis == 2:
      fs=(11.5,10.5)
   return fs

def ledge(axis):
   if axis == 2:
      le = 0.1
   else:
      le = 0.15
   return le

def ledgecb(axis):
   if axis == 2:
      le1, le2 = 0.015, 0.92
   else:
      le1, le2 = -0.005, 0.955
   return le1, le2

def chosxlab(axis):
   labs = ['y [kpc]', 'x [kpc]', 'x [kpc]']
   return labs[axis]

def chosylab(axis):
   labs = ['z [kpc]', 'z [kpc]', 'y [kpc]']
   return labs[axis]

def axticks(ax,axis,fsize):
   ax.set_xticks([-30, -20, -10, 0, 10, 20, 30])
   ax.set_xticklabels(['-30', '-20', '-10', '0', '10', '20', '30'], fontsize=fsize)
   if axis == 2:
      ax.set_yticks([-30, -20, -10, 0, 10, 20, 30])#, fontsize=fsize)
      ax.set_yticklabels(['-30', '-20', '-10', '0', '10', '20', '30'], fontsize=fsize)
   else:
      ax.set_yticks([-20, -10, 0, 10, 20])
      ax.set_yticklabels(['-20', '-10', '0', '10', '20'], fontsize=fsize)
   return ax

def scalevec(axis):
   if axis == 2:
      return 0.10
   else:
      return 0.08

def plotext(axis):
   if (stg.px_user != False):
      px = stg.px_user
   else:
      px = 25.0
   if (stg.pz_user != False):
      pz = stg.pz_user
   else:
      pz = 19.2
   if axis == 2:
      return [-px,px,-px,px]
   else:
      return [-px,px,-pz,pz]

def prepare_draw(lab,attr,axis):
   le = ledge(axis)
   if lab == 'RM':
      lef = 0.085
   else:
      lef = 0.09

   gs00 = gridspec.GridSpec(8, 8, left=0.0)
   gs00.update(left=lef, right=0.835, bottom=le, top=0.95, wspace=0.1, hspace=0.4)
   ax = py.subplot(gs00[:,:])

   ax.set_title(stg.title(lab,attr),fontsize=fsize)

   ax.set_xlabel(chosxlab(axis),fontsize=fsize)
   ax.set_ylabel(chosylab(axis),fontsize=fsize)
   ax = axticks(ax,axis,fsize)
   return ax

def draw_cb(lab,axis,cax):
   le = ledge(axis)
   le1, le2 = ledgecb(axis)
   cbar2ax = py.axes([0.85, le+le1, 0.025, le2-le])
   cb = py.colorbar(cax, cax=cbar2ax)
   cb.set_label(stg.ety(lab),fontsize=fsize)
   cb.set_label(stg.ety(lab),fontsize=fsize)

   if lab == 'RM':
      caxmin, caxmax = cax.get_clim()
      if (not stg.print_log or lab in stg.apply_linear):
         ticks = [ item for item in np.linspace(caxmin, caxmax, min(abs(caxmax - caxmin), 10)) ]
         tickl = [str(round(item,2)) for item in ticks]
         cb.set_ticks(ticks)
         cb.set_ticklabels(tickl)
      else: # Forces symlog-like scale. Not great, not terrible, but works. TODO improve ticks coverage of cbar, arbitrary parameters applied here.
         powers = [item for item in np.linspace(np.log10(stg.lin_threshold), caxmax, min(abs(ceil(caxmax) - ceil(caxmin)), 10))]
         ticks  = [item for item in powers]
         powers = [item for item in np.linspace(np.sign(caxmin) * np.log10(stg.lin_threshold), caxmin, min(ceil(caxmax) - ceil(caxmin), 10))]
         for item in powers: ticks.append(item)
         tickl  = [str(round(np.sign(item) * 10**abs(item),1)) for item in ticks]
         ticks.append(0);  tickl.append('0.')
         cb.set_ticks(ticks)
         cb.set_ticklabels(tickl)

   for t in cb.ax.get_yticklabels():
      t.set_fontsize(fsize)
   return cb

def draw_vecs(ax,vecs,axis):
# wp, wq - components of vectors to be displayed
# X, Y - a mesh of points to locate vectors

   wp, wq, X, Y = vecs

   cellsize = (X[-1][-1] - X[0][0]) / (max(np.shape(X)) + 1) # Assumes square-like cell, serves to scale vectors to cell-diagonal * dokvec
   scale_u = 2. * np.sqrt(2) * stg.scale_vec_maxPD / ( ax.figure.dpi * cellsize ) if stg.use_vec_scaling else scalevec(ax)

   r = 1 #np.sqrt(wq**2 + wp**2)
   Q1 = ax.quiver(X, Y, +0.5*wq/r, +0.5*wp/r, headwidth=0, minlength=0, color='black', width=0.002, scale_units = "dots" if stg.use_vec_scaling else "xy", scale = scale_u)
   Q2 = ax.quiver(X, Y, -0.5*wq/r, -0.5*wp/r, headwidth=0, minlength=0, color='black', width=0.002, scale_units = "dots" if stg.use_vec_scaling else "xy", scale = scale_u)

   qk = ax.quiverkey(Q1, 0.83, 0.04, 0.6, 'p = 60%',labelpos='E', coordinates='figure', color='black', fontproperties={'weight': 'bold', 'size': '24'})
   return ax

def draw_map(data,vecs,figext,axis,attr,plot_file,lab,ff,i_nl):
# data - table containing data to be displayed
# vmin_, vmax_ -  minimum i maksimum of the color scale

   if stg.print_log:
      if (lab not in stg.apply_linear):
         if (lab == "TP" or lab == "PI"):
            data = np.log10(data)
         else:
            data = np.where(data >  stg.lin_threshold,  np.log10(data),  data)
            data = np.where(data < -stg.lin_threshold, -np.log10(-data), data)

   vmin_, vmax_ = stg.fvmax(lab,ff,data,i_nl)

   img = py.figure(figsize=figsz(axis))
   ax = prepare_draw(lab,attr,axis)

   cax = ax.imshow(data, origin='lower',vmin=vmin_,vmax=vmax_,extent=figext,cmap=stg.colormap(lab))
   ax.axis(plotext(axis))

   cb = draw_cb(lab,axis,cax)

   if stg.print_vec:
      ax = draw_vecs(ax,vecs,axis)

   py.draw()
   if not os.path.exists('./radiomaps/'):
      os.makedirs('./radiomaps/')
   for ss in stg.suffix:
      py.savefig('./radiomaps/'+lab+plot_file+'.'+ss)
      print("Image storred in file: ",'./radiomaps/'+lab+plot_file+'.'+ss)
