module cresp_variables ! & constants
!pulled by COSM_RAY_ELECTRONS
!    use initcosmicrays, only: ncre
  implicit none
! indices to operate within cg%u structure
  
!    type cre_args
       integer(kind=4),allocatable, dimension(:)  :: cren ! < indices pointing cosmic ray electrons number per bin in cg%u(iarr_cre,...) in the range of (1:ncre)
       integer(kind=4),allocatable, dimension(:)  :: cree ! < indices pointing cosmic ray electrons energy per bin in cg%u(iarr_cre,...) in the range of (ncre+1:2*ncre)
       integer(kind=4)                            :: crepu! < index pointing cosmic ray electrons upper momentum cut in in cg%u(iarr_cre,...) at (2*ncre+1)
       integer(kind=4)                            :: crepl! < index pointing cosmic ray electrons upper momentum cut in in cg%u(iarr_cre,...) at (2*ncre+2)
!    end type cre_args

! type (cre_args) cre_table             ! cosmic ray electrons
       integer(kind=4),allocatable, dimension(:)  :: cre_table
  
  
  integer          , parameter :: order = 3
!   integer          , parameter :: nbin = 5
! 
!   real(kind=8)     , parameter :: t_max = 6.0d0 !10d0
  
  real(kind=8)     , parameter :: f_init  = 1.0d0
  real(kind=8)     , parameter :: q_init  = 5.0
  real(kind=8)                 :: p_min_fix = 1.0d4
  real(kind=8)                 :: p_max_fix = 1.0d5
  real(kind=8)                 :: p_lo = 1.05d4 !p_min_fix  !/10. !* 20.0
  real(kind=8)                 :: p_up = 1.05d5 !p_max_fix ! / 2.0    !9.9d0
  real(kind=8)                 :: w
  
!   real(kind=8)     , parameter :: dt_ini = 0.1

  real(kind=8)     , parameter :: q_big = 10d0
!   real(kind=8)     , parameter      :: cfl_cr  = 0.1d0 ! cfl factor for CR

  ! these will most probably be in types and will be modificated by the driver (piernik)
!   real(kind=8)                 :: u_d0 = 7.5d-1 ! 5.0d-1
!   real(kind=8)                 :: u_d
!   real(kind=8)                 :: u_b = 1d-7 !0d0 !5d-7!
!   real(kind=8)                 :: div_v = 0d0
!   real(kind=8)                 :: omega_d = 0.5d0 !0.1d0    ! frequency of div(v) oscilations



  integer                      :: c2nd, c3rd
  real(kind=8)                 :: p_lo_next, p_up_next  ! momemntum for spectrum cut-offs
  real(kind=8)                 :: n_tot, n_tot0, e_tot, e_tot0


end module cresp_variables
