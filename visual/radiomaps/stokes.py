import math as mth
import pylab as py
import numpy as np
import electrons
import settings as stg

def stokes_params(Bp,Bq,Bn,rho_ion,Ecrp,wave_data,ds,n3,**kwargs):
# This is the core procedure to compute Stokes parameters I, Q, U
# and the rotation measur RM for polarized radio waves
# emitted along the line of sight.
#
# Procedure arguments:
#
# Bp,Bq - components of magnetic fild in the plane of sky
# which is perpendicular to the line of sight
# Bn - component of magnetic field parallel to the line of sight
# rho_ion - density of ionized gas, which is responsible for Faraday rotation,
# formally it is the number density of thermal electrons n_e.
# Ecr - energy density of cosmic rays. In the spectral version of the code it is the spectral energy density.
# p - speatral index of CR electron energy spectrum (stg) - typically ~2.8
# nu_s - frequency of the observed synchrotron radiation
# cJnu - a constant scaling synchrotron emissivity.
# cJnu= 2.344*(1.60219)^(p-2)*10^(8-24p)*a(p) (stg)

# Loaded if stg.mode == "spectral"
   if (stg.mode == "spectral"):
      Ecre = kwargs.get("Ecre", [])
      Ncre = kwargs.get("Ncre", [])
      ncre = kwargs.get("ncre", 0)

   I_sum, Q_sum, U_sum, RM_sum, SI, Q, U = [], [], [], [], [], [], []
   nu_s, lambda_s, nu_2, lambda_2 = wave_data
   n = len(rho_ion)
# modulus of magnetic field vector
   #B_tot = np.sqrt((Bp**2+Bq**2+Bn**2))
# modulus of the perpendiculat (to the line of sight) component of B
   B_perp = np.sqrt(Bp**2 + Bq**2)
# magnetic field component parallet to the line of sight
   B_paral = Bn
# Faraday rotation for a single cell
   dRM = B_paral*rho_ion*ds
# TODO: The constant 0.808 is missing here
# RM is measured in rad m^{}−2}, the line-of-sight magnetic field B_{reg,\paral} in \muG, the thermal electron density ne in cm^{−3} and the line of sight l in pc.

# Integral of the rotation measure along the line of sight - in the discrete representation this is cumulative sum starting from the cell on the line of sight, which is closest to the observer. This quantity is needed for computations of the rotation of the polarization plane of the radiation emitted along the line of sight.
   RM = dRM.cumsum()*stg.RMunit()
   if stg.print_RM:
      # Total rotation measure along the line of sight
      RM_sum=RM[n-1]

   if stg.print_PI or stg.print_SI or stg.print_vec:
      # The change of polarization angle across the whole emitting region.
      delta_phi = RM*lambda_s**2

      # sin and cos delta_phi will be needed to compute Stokes parameters for the radiation received by an observer.
      sin_delta_phi = np.sin(delta_phi)
      cos_delta_phi = np.cos(delta_phi)
      # intrinsic polarization angle for the readiation emitted in each cell
      #phi_int is defined by  sin(phi_int) and cos(phi_int)
      sin_phi_int = np.zeros_like(B_perp)
      cos_phi_int = np.zeros_like(B_perp)
      bpn0 = np.where(B_perp != 0.0)
      sin_phi_int[bpn0] = Bq[bpn0] / B_perp[bpn0]
      cos_phi_int[bpn0] = Bp[bpn0] / B_perp[bpn0]
      # polarization angle of the radiation received by the observer is equal the sum of the intrinsic polarization angle and the rotation angle on the way fro the cell to the observer.
      sin_psi = sin_phi_int * cos_delta_phi + cos_phi_int * sin_delta_phi
      cos_psi = cos_phi_int * cos_delta_phi - sin_phi_int * sin_delta_phi
      # We need a double polarization angle to compute Stokes parameters
      sin_2psi = 2.0*sin_psi*cos_psi
      cos_2psi = cos_psi**2 - sin_psi**2

   if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
      # The total intensity of synchrotron radiation emitted form each cell of lengths ds along the line of sight
      I = np.zeros_like(B_perp)
      if stg.mode == 'simple':
         for i3 in range(n3):
            I[i3] = stg.cJnu*B_perp[i3]**((stg.p+1.0)/2.0) * (1.0/nu_s)**((stg.p-1.0)/2.0) * Ecrp[i3] * ds

      elif stg.mode == 'spectral':
         for i3 in range(n3):
            #elfq = electrons.crenpp(nu_s, ncre, B_perp[i3], Ecre[:,i3])            # DEPRECATED
            elfq = electrons.crenppfq(nu_s, ncre, B_perp[i3], Ecre[:,i3], Ncre[:,i3])
            I[i3] = np.sqrt(nu_s*B_perp[i3]) * elfq

         if stg.print_SI:
            I2 = np.zeros_like(B_perp)
            for i3 in range(n3):
               #elfq2 = electrons.crenpp(nu_2, ncre, B_perp[i3], Ecre[:,i3])        # DEPRECATED
               elfq2 = electrons.crenppfq(nu_2, ncre, B_perp[i3], Ecre[:,i3], Ncre[:,i3])
               I2[i3] = np.sqrt(nu_2*B_perp[i3]) * elfq2

   if stg.print_PI or stg.print_SI or stg.print_vec:
      # degree of polarization of synchrotron radiation emitted in a single cell
      PI = (stg.p + 1.0) / (stg.p + 7.0/3.0)   # BEWARE: stg.p should be replkaced by the spectral index of CR energy density derived  from the spectrum
      # Stokes parameters of each cell, corrected for Faraday rotation on the Way to the observer.
      Q = PI*I*cos_2psi
      U = PI*I*sin_2psi

# Stokes parameters received by the observer are summed up for all cells along the line of sight
   if stg.print_PI or stg.print_SI or stg.print_vec or stg.print_TP:
      I_sum = I.sum()
   if stg.print_PI or stg.print_SI or stg.print_vec:
      Q_sum = Q.sum()
      U_sum = U.sum()
   if stg.print_SI:
      I2_sum = I2.sum()
      SI = np.log10(I2_sum/(I_sum+1.e-240)+1.e-200) / np.log10(nu_2/nu_s)

# The function returns the summed Stokes parameters I, Q and U, and rotation measure RM. The spectral index should be removed from here and computed in the plot_maps routine.
   return I_sum, Q_sum, U_sum, RM_sum, SI

def polarized(I,Q,U):
# I - total intensity of the synchrotron radiation (Total Power - TP).
# Q,U - Stokes parameters.
   PI = np.sqrt(Q**2 + U**2)
   PI_obs = (PI / (I+1e-6))
   print('Mean PI/I:             ', np.mean(PI_obs))
   return PI, PI_obs


def vector_direction(PI_obs,Q,U):
# PI_obs - intensity of the polarized radio emission.
# Q,U - Stokes parameters.
   psi_obs = Wp_pi = Wq_pi = np.zeros_like(PI_obs)
   psi_obs = 0.5*np.arctan2(U,-Q) + 0.5*np.pi
   Wp_pi = PI_obs * np.sin(psi_obs)
   Wq_pi = PI_obs * np.cos(psi_obs)
   return Wp_pi.T, -Wq_pi.T
