from numpy import log, log10, zeros, asfarray, sqrt, sign, pi
from settings import MHz, const_synch, maxcren, arr_dim_q, q_big, q_eps # WARNING these are hard-coded in settings

# General constants and arrays here
# Optimize (minimize duplicated executions)
p_fix               = []           # initialized below (in initialize_crspectrum_tools)
p_fix3              = []           # initialized below (in initialize_crspectrum_tools), p_fix in 3rd power
_16p1_MHz           = 16.1 * MHz   # needed by nu_B_to_p
_fourXpi            = 4. * pi      # needed by nq2f and crenppfq
_inittab_dim        = 0            # initialized below (in prepare_q_grid)
_log10_pmax_by_pmin = 1.0          # initialized below (in initialize_crspectrum_tools)
_log10_enpc_ratio   = 1.0          # initialized below (in prepare_q_grid)
_ncre_m2            = 0.           # initialized below (in initialize_crspectrum_tools)
_const_1956_x_sqrt_nu_by_16p1MHz      = []
_smallB             = 1.e-24
_one                = 1.0
_zero               = 0.0
# FIXME / BEWARE for some reason ncre, p_min_fix and p_max_fix retain their default, hard-coded values from settings,
# even if 'read_CRESP_params' is run at 'main' before 'initialize_crspectrum_tools'.
# The below three serves as in-module storage for the right values, provided at initialization
_p_min_fix          = 0.0
_p_max_fix          = 0.0
_ncre               = 0

def initialize_crspectrum_tools(ncre, p_min_fix, p_max_fix, nu):
   global p_fix, p_fix3, _ncre_m2, _log10_pmax_by_pmin, _const_1956_x_sqrt_nu_by_16p1MHz, _p_min_fix, _p_max_fix, _ncre
# var_names now should be initialized with values from problem.par@hdf file
# Initialize quantities dependent on read parameters and remaining constant throughout computation
   _p_min_fix           = p_min_fix
   _p_max_fix           = p_max_fix
   _ncre                = ncre
   _ncre_m2             = float(ncre - 2.)
   _log10_pmax_by_pmin  = log10(_p_max_fix / _p_min_fix)
# Reconstruct momenta
   edges = []
   edges[0:ncre+1] = range(0,ncre+1, 1)
   p_fix[0:ncre+1] = zeros(ncre+1)

   log_width   = (_log10_pmax_by_pmin / _ncre_m2)
   p_fix_ratio = 10.0 ** log_width

   for i in range(0,ncre-1):
      p_fix[i+1] = p_min_fix * p_fix_ratio**float(i)
   p_fix[0]    = ( sqrt(p_fix[1] * p_fix[2]) ) / p_fix_ratio               #let us rid of zero values
   p_fix[ncre] = ( sqrt(p_fix[ncre-2]* p_fix[ncre-1]) ) * p_fix_ratio
   p_fix = asfarray(p_fix)
   p_fix3 = asfarray( p_fix ** 3)

   prepare_q_grid(p_fix_ratio)
   initialize_freqs(nu)

   return
#=============================================================================================================================
def initialize_freqs(nu):
   global _const_1956_x_sqrt_nu_by_16p1MHz
   for i in range(len(nu)):
      _const_1956_x_sqrt_nu_by_16p1MHz.append( 1956.8521 * sqrt( nu[i] / _16p1_MHz) )

#=============================================================================================================================
def nu_B_to_p(nu, B_perp): # Based on approximation by Mulcahy et al. (2018) (eqn. 2), https://arxiv.org/abs/1804.00752
   # WARNING !!! cutoff momenta are not precisely described here !!! WARNING
   p = _const_1956 * sqrt( (nu / _16p1_MHz ) / ( B_perp + _smallB ) ) # [B_perp] = mGs - assumed at input, [nu] = 16.1 MHz, but already present
   return p
#=============================================================================================================================
def nu_all_B_to_p(B_perp, nui): # Based on approximation by Mulcahy et al. (2018) (eqn. 2), https://arxiv.org/abs/1804.00752
   # As in function nu_B_to_p, but omitting some multiplications and powers
   p = _const_1956_x_sqrt_nu_by_16p1MHz[nui] / ( B_perp + _smallB )**0.5 # [B_perp] = mGs - assumed at input, [nu] = 16.1 MHz, but already present
   return p
#=============================================================================================================================
def nu_to_ind(nu_ind, B_perp, ncre):
   p_nu = nu_all_B_to_p(B_perp, nu_ind)
   nu_to_ind = int( (log10(p_nu/_p_min_fix)/_log10_pmax_by_pmin) * _ncre_m2 + _one)

   if nu_to_ind > ncre: nu_to_ind = ncre
   if nu_to_ind < 0: nu_to_ind = 0
   if ( not p_nu < p_fix[nu_to_ind] and p_nu > p_fix[nu_to_ind-1]): # WARNING temporary fix
      nu_to_ind = min(nu_to_ind +1,ncre)

   return nu_to_ind, p_nu
#=============================================================================================================================
def crenpp(nu_s, ncre, bperp, ecr):
   p_ind, p_nu = nu_to_ind(nu_s, bperp, ncre)
   n_ind = min(p_ind - 1, maxcren - 1)
   cren_i = ecr[n_ind]
   p_i    = p_fix[p_ind]
   if (p_nu < p_fix[p_ind] and p_nu > p_fix[p_ind-1]):
      delta_p = p_i - p_fix[p_ind-1]
   else:
      delta_p = _p_max_fix - p_i  # WARNING temporary fix
   elfq = const_synch * cren_i / delta_p
   return elfq

#=============================================================================================================================
def crenppfq(nu_ind, ncre, bperp, ecr, ncr):      # recovers spectral index q from electron energy and number density
   p_ind, p_nu = nu_to_ind(nu_ind, bperp, ncre)   # and calculates emissivity
   i_bin  = min(p_ind - 1, maxcren-1)
   cren_i = ncr[i_bin]
   cree_i = ecr[i_bin]
   p_i    = p_fix[p_ind]
   p_im1  = p_fix[max(p_ind - 1, 0)]
   p_im13 = p_fix3[max(p_ind - 1, 0)]

   q1     = 4.1
   npq_nu = _zero
   if (cree_i > _zero and cren_i > _zero):
      enpc_bin = cree_i / (cren_i * _one * p_im1)
      q_bin     = interpolate_q(enpc_bin, q1)  # WARNING TODO MAGIC NUMER; rather read and use q_init value from 'problem.par'
      # slv_error = False
      # q1, slv_error = ne_to_qNR(q1, enpc_bin, p_fix_ratio, slv_error) # possible, but computationally expensive
      f_binXfourXpi  = nq2fXfourXpi(cren_i, q_bin, p_im1, p_im13, p_i) # zwraca wartosc f(p,q) na lewej scianie wybranego binu
      npq_nu         = f_binXfourXpi * ((p_nu / p_im1)**(- q_bin)) * (p_nu)**2   # recovered N(p_nu)

   elfq = const_synch * npq_nu
   return elfq

#=============================================================================================================================
def prepare_q_grid(p_fix_ratio):   # fills grid with values of spectral indices 'q' in range of given e / (n p c), for further q interpolation
   global enpc_tab_q, q_grid, _inittab_dim, _log10_enpc_ratio
   _inittab_dim = int(arr_dim_q / 4)

   q_grid    = zeros(arr_dim_q)    # will contain data for later interpolation
   q_inittab = zeros(_inittab_dim) # contains initial guesses for root finding procedure

   q_grid[:]               = q_big
   q_grid[int(arr_dim_q/2):] = -q_big

   for i in range(1, int(0.5 * _inittab_dim)):
      q_inittab[i-1] = ln_eval_array_val(i, q_big, float(0.05), 1, int(0.5 * _inittab_dim -1))
   for i in range(0, int(0.5 * _inittab_dim)+1):
      q_inittab[int(0.5 * _inittab_dim) + i - 1] = -q_inittab[int(0.5 * _inittab_dim) - i]

   enpc_max       = 1.1 * p_fix_ratio# / clight
   enpc_min       = 1.0 + 5e-8      # WARNING Magic number
   enpc_tab_q     = zeros(arr_dim_q)
   enpc_tab_q[:]  = enpc_min

   j = min(arr_dim_q - int(arr_dim_q / (arr_dim_q / 100.)), arr_dim_q - 1)
   while (q_grid[j] <= (-q_big) and (q_grid[arr_dim_q - 1] <= (-q_big))):
      enpc_max = enpc_max - enpc_max * 0.005    # WARNING Magic number
      for i in range(0, arr_dim_q):
         enpc_tab_q[i]  = enpc_min * 10.0**((log10(enpc_max / enpc_min)) / float(arr_dim_q) * float(i))
      fill_q_grid(p_fix_ratio, q_inittab)

   _log10_enpc_ratio = log10(enpc_tab_q[-1] / enpc_tab_q[0])
   return
#=============================================================================================================================
def ln_eval_array_val(i, arr_min, arr_max, min_index, max_index): #DEPRECATED here
      b = (log(float(max_index)) -log(float(min_index))) / (arr_max - arr_min)
      ln_eval_array_val = (arr_min - log(float(min_index)) / b ) + log(float(i)) / b
      return ln_eval_array_val
#=============================================================================================================================
def fill_q_grid(p_fix_ratio, q_init):
   global q_grid, enpc_tab_q
   previous_solution = q_grid[int(len(q_grid)/2)]
   solving_error     = True
   qx                = previous_solution

   for i in range(1, arr_dim_q, 1):
      qx, solving_error = ne_to_qNR(previous_solution, enpc_tab_q[i], p_fix_ratio, solving_error)
      if (solving_error):
         for j in range(1, _inittab_dim, 1):
            qx = q_init[j]
            qx, solving_error = ne_to_qNR(qx, enpc_tab_q[i], p_fix_ratio, solving_error)
            if (not solving_error):
               q_grid[i]     = qx
               prev_solution = qx
      else:
         q_grid[i]         = qx
         previous_solution = qx
   return
#=============================================================================================================================
def ne_to_qNR(x_init_guess, alpha, p_ratio, slv_error):   # Newton-Raphson-type 1-D root finding algorithm
  iter_limit = 30
  tol_f      = 1.0e-9
  min_dx     = 1.e-10
  x = x_init_guess
  df = 1.0

  for i in range(iter_limit):
        if abs(x) >= q_big:
           x = (x/abs(x)) * q_big
           slv_error = True
           break
        dx = min(x * 1e-3, 10.e-2)
        dx = max(abs(dx),  min_dx)  # restrict dx to positive values
        df = 0.5 * (alpha_to_q(x+dx, alpha, p_ratio)-alpha_to_q(x-dx, alpha, p_ratio))/dx
        delta = -alpha_to_q(x, alpha, p_ratio)/df
        if abs(delta) <= tol_f:
            slv_error = False
            return x, slv_error
        else:
            x = x + delta
  x_result = x

  return x_result, slv_error
#=============================================================================================================================
def alpha_to_q(x, alpha, p_ratio_4_q): # Using this function (Miniati, 2001, eq. 29; present also in cresp_NR_method.F90)
      q_in3 = 3. - x                   # we nunerically find spectral index 'q'.
      q_in4 = 1. + q_in3
      if (abs(q_in3) < q_eps):
         alpha_to_q = (p_ratio_4_q**q_in4 - 1.)/log(p_ratio_4_q)
      elif (abs(q_in4) < q_eps):
         alpha_to_q = q_in3 * log(p_ratio_4_q)/(p_ratio_4_q**q_in3 - 1.)
      else:
         alpha_to_q = q_in3/q_in4 * (p_ratio_4_q**q_in4 - 1.)/(p_ratio_4_q**q_in3 - 1.)
      alpha_to_q = alpha_to_q - alpha

      return alpha_to_q
#=============================================================================================================================
def interpolate_q(alpha, q_error):  # Finds value of spectral index 'q' by provided 'alpha' = e/npc; faster than running ne_to_qNR
   index = int((log10(alpha / enpc_tab_q[0]) / _log10_enpc_ratio) * (arr_dim_q - 1) ) # + 1
   if (index < 0 or index > arr_dim_q - 2):
      q_out = q_error # NOTICE assuming more natural value, instead of catastrophically high / low q_big
   else:
      index2 = index + 1      # Prepare and linearly interpolate
      q_out  = q_grid[index] + (alpha - enpc_tab_q[index]) * ( q_grid[index] - q_grid[index2]) / (enpc_tab_q[index] - enpc_tab_q[index2])

   return q_out
#==============================================================================================================================

def nq2fXfourXpi(n, q, p_l, p_l3, p_r):
      qm3       = q - 3.
      three_m_q = 3. - q
      if p_r > _zero and p_l > _zero :
        pr_by_pl = p_r / p_l
        #nq2f = n / (_fourXpi * p_l3) # NOTICE if we want f
        nq2fXfourXpi = n / (p_l3)     # NOTICE but we want f * 4 * pi
        if ( abs(qm3) > q_eps ):
            nq2fXfourXpi = nq2fXfourXpi * three_m_q / (( pr_by_pl )**three_m_q - _one)
        else:
            nq2fXfourXpi = nq2fXfourXpi / log(pr_by_pl)
      return nq2fXfourXpi
