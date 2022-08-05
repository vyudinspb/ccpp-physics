module wamphys_common
!
     use machine,  only : kind_phys
                         
     implicit none
     
      real(kind=kind_phys)   ::  pi, pi2, pih, rad_to_deg, deg_to_rad
      real(kind=kind_phys)   ::  pid12, pi_24hr, pi2_365d, pid3, pid9, pid6
      real(kind=kind_phys)   ::  arad, p0s     
      real(kind=kind_phys)   ::  grav, grav2, rgrav, rgrav2, g0r2e
      real(kind=kind_phys)   ::  cpd,  rd, rv, fv            
      real(kind=kind_phys)   ::  rdi,  rcpd, rcpd2  
      
      real(kind=kind_phys)   ::  gor,  gr2,  grcp, gocp, rcpdl, grav2cpd
      real(kind=kind_phys)   ::  bnv2min, bnv2max
      real(kind=kind_phys)   ::  dw2min,  velmin, minvel        
      real(kind=kind_phys)   ::  omega1, omega2,   omega3   
      real(kind=kind_phys)   ::  hpscale, rhp, rhp2, rh4, rhp4, khp
      real(kind=kind_phys)   ::  mkzmin, mkz2min,  mkzmax, mkz2max, cdmin    
      real(kind=kind_phys)   ::  rcpdt
      real(kind=kind_phys)   ::  hpskm  != 7.00           
!
end module wamphys_common

module wamphys_const                
!
  use machine, only: kind_phys
  use physcons
  public
  
!SIK_VAY 2022> \name WAM/idea physics related constants (idea_composition for radiation)

  real(kind=kind_phys),parameter:: con_amo   =0.5_kind_phys*con_amo2               !< molecular wght of O (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amn  =14.067_kind_phys                      !< molecular wght of N (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amn2  =2.0_kind_phys*con_amN                !< molecular wght of N2 (\f$g/mol\f$) 
  real(kind=kind_phys),parameter:: con_amno  =con_amo+con_amn                      !< molecular wght of NO (\f$g/mol\f$) 
  real(kind=kind_phys),parameter:: con_amH  = 1.00784_kind_phys                    !< molecular wght of H (\f$g/mol\f$) 
  real(kind=kind_phys),parameter:: con_amHe  = 4.0026_kind_phys                    !< molecular wght of He (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amH2o  =  con_amo +  con_amH+con_amH  
  real(kind=kind_phys),parameter:: con_pi2     = 2.0d0*con_pi                      !< 2*pi     
  real(kind=kind_phys),parameter:: con_dtr     = con_pi/180.0d0                    !< pi/180.
  real(kind=kind_phys),parameter:: con_rtd     = 1.0d0/con_dtr                     !< 180./pi        
!........................................!

!> \name Secondary constants

     
!........................................!
  real(kind=kind_phys) , parameter :: wam_avgd = con_avgd*1000.0d0          !wam_avgd=6.022e26 ! avogadro const 1/kmol
  real(kind=kind_phys) , parameter :: rbz   = 1.d0/ con_boltz   
  real(kind=kind_phys) , parameter :: rmo   = wam_avgd/con_amo
  real(kind=kind_phys) , parameter :: rmo3   = wam_avgd/con_amo/3.  
  real(kind=kind_phys) , parameter :: rmo2  = wam_avgd/con_amo2
  real(kind=kind_phys) , parameter :: rmh2o = wam_avgd/con_amh2o
  real(kind=kind_phys) , parameter :: rmn2  = wam_avgd/con_amn2 
  real(kind=kind_phys) , parameter :: ELCH  = 1.602e-19    
  real(kind=kind_phys) , parameter :: pi = con_pi
  real(kind=kind_phys) , parameter :: pi2 =  pi + pi
  real(kind=kind_phys) , parameter :: Pid6 =PI/6., Pid2 =PI/2., Pid4 =PI/4.
  real(kind=kind_phys) , parameter :: Pid3 =PI/3., Pid9 =PI/9., Pid12 =PI/12.
  real(kind=kind_phys) , parameter :: Pid18 =PI/18., pi2_365d = pi2/365.
  real(kind=kind_phys) , parameter :: pi_24hr = pi2/24.0_kind_phys
  
  real(kind=kind_phys) , parameter :: DTR =PI/180.0, R_2_D =1.0/DTR, fac_lst = R_2_D/15.0
  
  real(kind=kind_phys) , parameter ::  rad_to_deg=r_2_d, deg_to_rad=dtr
  real(kind=kind_phys) , parameter :: REARTH=6.370E06
  real(kind=kind_phys) , parameter :: YDAYS = 365.0, RYDAYS =1./YDAYS
  
  real(kind=kind_phys) , parameter ::  vmr_nzero = 1.e-36, mmr_nzero = vmr_nzero ! min values for mixing ratios 
  real(kind=kind_phys) , parameter ::  con_nzero = vmr_nzero*1.e19               ! 1/m3  
  real(kind=kind_phys) , parameter ::  mmr_min = 1.e-32, mmr_max=0.9999999
              
  end module wamphys_const     
!
!
