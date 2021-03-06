!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
module fluids_pub
! pulled by ANY
   implicit none
   public    ! QA_WARN nothing to hide here

#ifdef IONIZED
   logical, protected                         :: has_ion = .true.
#else /* !IONIZED */
   logical, protected                         :: has_ion = .false.
#endif /* !IONIZED */

#ifdef NEUTRAL
   logical, protected                         :: has_neu = .true.
#else /* !NEUTRAL */
   logical, protected                         :: has_neu = .false.
#endif /* !NEUTRAL */

#ifdef DUST
   logical, protected                         :: has_dst = .true.
#else /* !DUST */
   logical, protected                         :: has_dst = .false.
#endif /* !DUST */

   real :: cs2_max ! maximum of isothermal sound speeds over all fluids

end module fluids_pub
