!>\file wamphys_swdef.F90
!! This file contains aalocatable array/scalar definitions for WAM physics.

!>\ingroup mod_GFS_phys_time_vary
!! This module defines arrays/scalars in SW drivers of WAM physics

module wamphys_swdef
!> \section arg_table_wamphys_swdef
!! \htmlinclude wamphys_swdef.html
!!
! swio-data     
   use machine , only : kind_phys
   implicit none  
      character(len=255)  :: swio_file='swio_waminput.nc'
      integer :: sw_ntd
      real(kind=kind_phys), dimension(:), allocatable   ::   sw_f107, sw_f107d,  sw_kp, sw_kpa
      real(kind=kind_phys), dimension(:), allocatable   ::   sw_nhp,   sw_nhpi,  sw_shp,sw_shpi   
      real(kind=kind_phys), dimension(:), allocatable   ::   sw_den, sw_ang,  sw_bz, sw_bt, sw_vel
      real(kind=kind_phys), dimension(:), allocatable   ::   sw_ap, sw_apa     
      real(kind=kind_phys), dimension(:), allocatable   ::   sw_time
      real(kind=kind_phys), dimension(:), allocatable   ::   sw_days_since                 
!  
! curent step values of  swio-data
!      
!      real(kind=kind_phys) :: csw_f107, csw_f107d,  csw_kp, csw_kpa, csw_ap, csw_apa  
!      real(kind=kind_phys) :: csw_nhp,  csw_nhpi,  csw_shp,   csw_shpi, csw_den, csw_ang 
!      real(kind=kind_phys) :: csw_bz,   csw_bt,    csw_vel,   csw_time     
!      integer              :: sw_ktprev
!        
!wamphys_def_ipewam_cpl
!
!
     integer, parameter                  :: lowipe_lev150 = 80
!    integer, parameter                  :: lonr = 2
!    integer, parameter                  :: latr = 2   
     integer, parameter                  :: levsr = 150 
     integer                             :: imx           ! horzontal dimension of current PE)               
!
!  GSM allocation on "lonr, lats_node, k8o:nlevs"  => FV3WAM-style (im,   k8o:nlevs)
!  ALLOCATE(ZMT(lonr, lats_node_r_max, lowst_ipe_level:levs))     
!     ../gsm/phys/gfs_physics_initialize_mod.f:        lowst_ipe_level = 80
!
!      real(kind=kind_phys), allocatable ::     ZMT_3d(:,:,:), MMT_3d(:,:,:)
!      real(kind=kind_phys), allocatable ::     JHR_3d(:,:,:), SHR_3d(:,:,:) 
!      real(kind=kind_phys), allocatable ::     O2DR_3d(:,:,:) 
      
!      real(kind=kind_phys), allocatable ::     zmt(:,:), mmt(:,:)
!      real(kind=kind_phys), allocatable ::     jhr(:,:), shr(:,:) 
!      real(kind=kind_phys), allocatable ::     o2dr(:,:) 
!               
!      fv3-allocation:    mmt(im, lowst_ipe_lev150: levs)
!      vs   allocate(mmt_3d(lonr, lats_node_r_max, lowst_ipe_level:levs)) 
!
!../../framework/scripts/ccpp_fortran_to_metadata.py wamphys_swdef.F90
end module wamphys_swdef      


  

