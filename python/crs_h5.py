#!/usr/bin/python
from pylab import zeros, sqrt, size
import matplotlib.pyplot as plt
from numpy import log10, log, pi, asfarray, array, linspace, sign
import h5py
import os
import sys

#------ default values of parameters --------
e_small = 1.0e-6
eps     = 1.0e-15
ncre      = 45
p_min_fix = 0.4e0
p_max_fix = 1.65e4
cre_eff = 0.01
q_big   = 30.

arr_dim = 200
helper_arr_dim = int(arr_dim / 4)

c = 1.0 # 0.3066067E+06  # PSM -> 0.3066067E+06, SI -> 0.2997925E+09

first_run = True
fixed_width = True
got_q_tabs = False
q_explicit = True

def nq2f(n,q,p_l,p_r):
      if p_r> 0.0 and p_l > 0 :
        pr_by_pl = p_r / p_l
        nq2f = n / (4* pi * p_l**3)
        if abs(q-3.0) > eps:
            nq2f = nq2f *(3.0 - q) / (( pr_by_pl )**( 3.0 - q ) - 1.0)
        else:
            nq2f = nq2f / log(p_r/p_l)
      return nq2f
# 1D Newton-Raphson algorithm (to find q):
#
def nr_get_q(q_start, alpha, p_ratio,exit_code):
  iter_limit = 30
  tol_f      = 1.0e-9
  x = q_start
  df = 1.0
  for i in range(iter_limit):
        if abs(x) >= q_big:
           x = (x/abs(x)) * q_big
           exit_code = True
           break
        dx = min(x*1e-3,10e-2)
        dx = sign(dx) * max(abs(dx), 1.0e-10)
        df = 0.5*(fun(x+dx, alpha, p_ratio)-fun(x-dx, alpha, p_ratio))/dx
        delta = -fun(x, alpha, p_ratio)/df
        if abs(delta) <= tol_f:
            exit_code = False
            return x, exit_code
        else:
            x = x + delta
  nr_get_q = x

  return nr_get_q, exit_code

# function used to find q: ----------------------
def fun(x, alpha, p_ratio):
      if abs(x - 3.0) < 1.0e-3:
         fun = -alpha + (-1.0 + p_ratio)/log(p_ratio)
      elif abs(x - 4.0) < 1.0e-3:
         fun = -alpha + p_ratio*log(p_ratio)/(p_ratio - 1.0)
      else:
         fun = -alpha + ((3.0-x)/(4.0-x))*((p_ratio**(4.0-x)-1.0)/(p_ratio**(3.0-x)-1.0))
      return fun

def prepare_q_tabs():
   global alpha_tab_q, q_grid, q_big, q_space
   q_grid = zeros(arr_dim) # for later interpolation
   q_space= zeros(helper_arr_dim) # for start values

   q_grid[:]      = q_big
   q_grid[int(arr_dim/2):] = -q_big

   def ln_eval_array_val(i, arr_min, arr_max, min_i, max_i):
      b = (log(float(max_i)) -log(float(min_i)))/ (arr_max - arr_min)
      ln_eval_array_val = (arr_min-log(float(min_i))/b ) + log(float(i)) / b
      return ln_eval_array_val

   for i in range(1,int(0.5*helper_arr_dim)):
      q_space[i-1] = ln_eval_array_val(i, q_big, float(0.05), 1 , int(0.5*helper_arr_dim -1))

   for i in range(0, int(0.5*helper_arr_dim)+1):
      q_space[int(0.5 * helper_arr_dim)+i-1] = -q_space[int(0.5 * helper_arr_dim)-i]

   a_max_q     = 10. #* p_fix_ratio# / clight
   a_min_q     = 1.00000005
   alpha_tab_q = zeros(arr_dim)
   alpha_tab_q[:] = a_min_q

   j = arr_dim - int(arr_dim/(arr_dim/10.))
   while (q_grid[j] <= (-q_big) and (q_grid[arr_dim-1] <= (-q_big)) ):
      a_max_q = a_max_q * 0.95
      for i in range(0,arr_dim):
         alpha_tab_q[i]  = a_min_q * 10.0**((log10(a_max_q/a_min_q))/float(arr_dim)*float(i))
      fill_q_grid() # computing q_grid takes so little time, that saving the grid is not necessary.
   return

def fill_q_grid():
   global q_grid, alpha_tab_q, p_fix_ratio, arr_dim, helper_arr_dim, q_space
   previous_solution = q_grid[int(len(q_grid)/2)]
   exit_code         = True
   x                 = previous_solution
   for i in range(1,arr_dim,1):
      x,exit_code = nr_get_q(previous_solution, alpha_tab_q[i],p_fix_ratio, exit_code)
      if exit_code == True: # == True
         for j in range(1,helper_arr_dim,1):
            x = q_space[j]
            x,exit_code = nr_get_q(x,alpha_tab_q[i],p_fix_ratio, exit_code)
            if exit_code == False:
               q_grid[i] = x
               prev_solution = x
      else: # exit_code == false
         q_grid[i]         = x
         previous_solution = x
   return

def interpolate_q(alpha):
   global arr_dim, alpha_tab_q, q_grid
   index = int((log10(alpha/alpha_tab_q[0])/log10(alpha_tab_q[-1]/alpha_tab_q[0])) * (arr_dim - 1))# + 1
   if (index < 0 or index > arr_dim-1):
      index = max(0, min(arr_dim-1, index))
      q_out = q_grid[index]
   else:
      index2 = index+1
      q_out  = q_grid[index] + (alpha - alpha_tab_q[index]) * ( q_grid[index] - q_grid[index2]) / (alpha_tab_q[index] - alpha_tab_q[index2])

   return q_out
# plot data ------------------------------------
def plot_data(plot_var, pl, pr, fl, fr, q, time, location, i_lo_cut, i_up_cut):
   global first_run, e_small
   f_lo_cut = fl[0] ;      f_up_cut = fr[-1]
   p_lo_cut = pl[0] ;   p_up_cut = pr[-1]

   if plot_var  == 'f' :
      plot_var_l      = fl
      plot_var_r      = fr
      plot_var_lo_cut  = f_lo_cut
      plot_var_up_cut = f_up_cut
   elif plot_var == 'n' :
      plot_var_l      = 4*pi*fl*pl**2
      plot_var_r      = 4*pi*fr*pr**2
      plot_var_lo_cut = 4*pi*f_lo_cut*p_lo_cut**2
      plot_var_up_cut = 4*pi*f_up_cut*p_up_cut**2
   elif plot_var == 'e' :
      plot_var_l      = 4*pi*c**2*fl*pl**3
      plot_var_r      = 4*pi*c**2*fr*pr**3
      plot_var_lo_cut = 4*pi*c**2*f_lo_cut*p_lo_cut**3
      plot_var_up_cut = 4*pi*c**2*f_up_cut*p_up_cut**3

   s = plt.subplot(122)
   plt.cla()
   s.set_xscale('log')
   s.set_yscale('log')
   plt.xlabel('p')
   plt.ylabel(plot_var)
   global plot_p_min, plot_p_max, plot_var_min, plot_var_max
   if first_run :
      plot_p_min    =  p_lo_cut
      plot_p_max    =  p_up_cut
      plot_var_min = 0.1*e_small
      first_run = False
      if (plot_var == "e"):
        plot_var_min = 1.0e-4 * e_small
      elif (plot_var == "f" ):
        plot_var_min = e_small / (4*pi * (c ** 2)  * p_max_fix **3) /10.
      elif (plot_var == "n" ):
        plot_var_min = 0.1 / (4*pi * p_max_fix ** 2)

   if fixed_width == True:
      plt.xlim (0.1 * plot_p_min   ,  10.*plot_p_max )
      plt.xlim

   plt.ylim (plot_var_min , plot_ymax)
   plt.grid()

#plot floor value
   p_range = linspace(s.get_xlim()[0],s.get_xlim()[1])
   e_smalls = zeros(len(p_range))
   e_smalls[:] = e_small
   if (plot_var == "e"):
      plt.plot(p_range, e_smalls, color="green", label="$e_{small}$")
   elif(plot_var == "n"):
      plt.plot(p_range, e_small/(c*p_range), color="green",label="$n_{small}$")

   for i in range(0, size(fr)) :
      plt.plot([pl[i],pr[i]],[plot_var_l[i],plot_var_r[i]],color='r')
      plt.plot([pl[i],pl[i]], [plot_var_min,plot_var_l[i]],color='r')
      plt.plot([pr[i],pr[i]], [plot_var_min,plot_var_r[i]],color='r')
   s.set_facecolor('white')
   plt.title("Spectrum of %s(p) \n Time = %7.3f | location: %7.2f %7.2f %7.2f " % (plot_var, time, location[0],location[1],location[2]) )

   return s

def simple_plot_data(plot_var, p, var_array, time, location, i_lo_cut, i_up_cut):
   global first_run, plot_p_min, plot_p_max, fixed_width
   p_lo_cut = p[0] ;   p_up_cut = p[-1]
   var_array = var_array / p # does this do the trick with correct q?
   s = plt.subplot(122)
   plt.cla()
   s.set_xscale('log')
   s.set_yscale('log')
   plt.xlabel('p')
   plt.ylabel(plot_var)

   if first_run :
      plot_p_min    =  p_lo_cut
      plot_p_max    =  p_up_cut
      plot_var_min = 0.1*e_small
      first_run = False
      if (plot_var == "e"):
        plot_var_min = e_small
      elif (plot_var == "f" ):
        plot_var_min = e_small / (4*pi * (c ** 2)  * p_max_fix **3) /10.
      elif (plot_var == "n" ):
        plot_var_min = 0.1 / (4*pi * p_max_fix ** 2)

   if fixed_width == True:
      plt.xlim (0.25 * plot_p_min   ,  5.*plot_p_max )
      plt.xlim

   plt.ylim (plot_var_min , plot_ymax)
   plt.grid()

   plt.plot(p,var_array,color='r')
   for i in range(i_lo_cut,i_up_cut):
       plt.plot([p[i],p[i]],[plot_var_min, var_array[i]],color="r")

#plot floor values
   p_range = linspace(s.get_xlim()[0],s.get_xlim()[1])
   e_smalls = zeros(len(p_range))
   e_smalls[:] = e_small

   s.set_facecolor('white')
   plt.title(" %s(p) \n Time = %7.3f | location: %7.2f %7.2f %7.2f " % (plot_var, time, location[0],location[1],location[2]) )

   return s
#-----------------------------------------------------------------

def crs_plot_main(parameter_names, parameter_values, plot_var, ncrs, ecrs, field_max, time, location, use_simple):
    global first_run, got_q_tabs

    try:
        for i in range(len(parameter_names)):
            exec("%s = %s" %(parameter_names[i], parameter_values[i]))
    except:
        print "Exiting: len(names) not equal len(values)"
        sys.exit()

# TODO -------- do it under *args TODO
    fixed_width = True
# TODO -------- do it under *args TODO

    first_run = True
# -------------------
    global plot_ymax, p_fix_ratio
    plot_ymax = field_max * cre_eff
    edges = []
    p_fix = []
    edges[0:ncre] = range(0,ncre+1, 1)
    p_fix[0:ncre] = zeros(ncre+1)

# WARNING !!! cutoff momenta are not precisely described here !!! WARNING
    log_width   = (log10(p_max_fix/p_min_fix))/(ncre-2.0)
# TODO: do it in the first run / calling script TODO
    if first_run:
        for i in range(0,ncre-1):
            p_fix[i+1]  = p_min_fix * 10.0**( log_width * edges[i])
            p_fix_ratio = 10.0 ** log_width
            p_fix[0]    = ( sqrt(p_fix[1] * p_fix[2]) ) / p_fix_ratio
            p_fix[ncre]    = ( sqrt(p_fix[ncre-2]* p_fix[ncre-1]) ) * p_fix_ratio
            p_fix = asfarray(p_fix)

    p_mid_fix = zeros(ncre)
    p_mid_fix[1:ncre-1] = sqrt(p_fix[1:ncre-1]*p_fix[2:ncre])
    p_mid_fix[0]        = p_mid_fix[1] / p_fix_ratio
    p_mid_fix[ncre-1]   = p_mid_fix[ncre-2] * p_fix_ratio
    p_mid_fix = asfarray(p_mid_fix)
    i_lo = 0
    i_up = ncre
    empty_cell = True

    e_small = e_small / 10.
#------------ locate cutoff indices
    for i in range(1,ncre):
        i_lo = i
        if ecrs[i] >= e_small:
            empty_cell = False # not plotting if empty
            break
    for i in range(ncre,1,-1):
        i_up = i
        if ecrs[i-1] >= e_small: break

    sys.stdout.write('Time = %6.2f |  i_lo = %2d, i_up = %2d, %11s.'%(time, i_lo if not empty_cell else 0, i_up if not empty_cell else 0, '(empty cell)' if empty_cell else ' '))
    pln = p_fix[i_lo:i_up]
    prn = p_fix[i_lo+1:i_up+1]
    pln = array(pln)
    prn = array(prn)
    ncrs_act = ncrs[i_lo-2:i_up-2]
    ecrs_act = ecrs[i_lo-2:i_up-2]
    q_nr = [] ; fln = [] ; frn = []

    if (not got_q_tabs):
      prepare_q_tabs()
      got_q_tabs = True
    if (not q_explicit):
       print "Spectral indices q will be interpolated"
    else:
       print "Spectral indices q will be obtained explicitly"

    print "n = ", ncrs
    print "e = ", ecrs

    for i in range(0,i_up - i_lo):
         if (q_explicit == True):
            q_tmp = 3.5 ; exit_code = False
            q_tmp, exit_code = nr_get_q(q_tmp, ecrs[i+i_lo]/(ncrs[i+i_lo]*c*pln[i]), prn[i]/pln[i], exit_code)
         else:
            q_tmp = interpolate_q(ecrs[i+i_lo]/(ncrs[i+i_lo]*c*pln[i]))
         q_nr.append(q_tmp)
         fln.append(nq2f(ncrs[i+i_lo], q_nr[-1], pln[i], prn[i]))

    print "q = ", q_nr

    q_nr = array(q_nr)
    fln  = array(fln)
    frn  = array(fln)
    frn  = frn * (prn/pln) ** (-q_nr)
    plot = False

    if empty_cell==False:
         if (use_simple):
            var_array = zeros(ncre)
            if (plot_var == "e"):
               var_array = array(ecrs)
            elif (plot_var == "n"):
               var_array = array(ncrs)
            else: # plotvar == "f"
               var_array = array(fln) # TODO not fully prepared
            plot = simple_plot_data(plot_var, p_mid_fix, var_array, time, location, i_lo, i_up)
         else:
            plot = plot_data(plot_var, pln, prn, fln, frn, q_nr, time, location, i_lo, i_up)

    return plot, empty_cell
