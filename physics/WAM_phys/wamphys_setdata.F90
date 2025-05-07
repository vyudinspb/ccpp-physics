!======================================= 
module wamphys_setfix_tracers
   use machine , only : kind_phys

   implicit none  
      integer, parameter   :: np=65                          ! number of pressure levels of orig
      integer, parameter   :: nvmr=15
      integer, parameter   :: npj=17       
      real(kind=kind_phys) :: ohi(np),ho2i(np)               ! units of conc-n 1/cm3 or 1/m3
      data ohi/9.5e12,1.3e13,1.6e13,1.8e13,2.0e13,  &
        2.1e13,2.2e13,1.9e13,1.3e13,6.0e12,2.1e12,  &
        1.0e12,3.0e11,1.0e11,4.0e10,1.6e10,7.0e9,   &
        3.2e9,1.2e9,5.0e8,2.0e8,44*0.0/
      data ho2i/7.0e12,9.0e12,1.2e13,1.3e13,1.6e13, &
        1.7e13,1.7e13,1.3e13,7.0e12,2.5e12,7.0e11,  &
        1.5e11,3.5e10,8.0e9,2.5e9,9.0e8,3.0e8,      &
        1.4e8,6.0e7,2.5e7,1.0e7,44*0.0/
	
! old J17	
      real(kind=kind_phys) :: JI(npj),FHT(npj),CI(npj)
      real(kind=kind_phys) :: J17(npj)      
      DATA CI/8*0.900,0.680,0.43,0.18,6*-0.066/
      DATA JI/.4e-8,.78e-8,1.5e-8,3.e-8,6.8e-8,.15e-6,.34e-6,.77e-6,  &  
       1.07e-6,1.35e-6,1.6e-6,1.81e-6,2.05e-6,2.23e-6,2.36e-6,2.5e-6, 2.57e-6/        
      DATA FHT/8*1.2,1.85,2.50,3.150,6*3.8/

!                                         see IDEA_TRACERS_INPUT
   
end module wamphys_setfix_tracers
 
module wamphys_set_merge_rad

   use machine , only : kind_phys
   implicit none  
! merging scheme for heating rates

      real(kind=kind_phys), parameter     :: xb=7.0+0.6, xt=8.5+1.5      ! for Hp = 7 km: 52.5 km < Z_logp < 59.5 km 
      real(kind=kind_phys), parameter     :: xbl= 0.99*xb        ! xbl < xb
      real(kind=kind_phys), parameter     :: rdx=1./(xt-xb)
      real(kind=kind_phys), parameter     :: xlogps = 11.5129    ! alog(1.e5=Ps_in_Pa)
      real(kind=kind_phys), parameter     :: prdot02 = 1.e5*exp(-xbl)      ! mb because pr = pr_idea in (mb)
      integer                             :: npsrad              ! layer where Pressure < 0.02   (2Pa)      
      ! nps-pressure index to start WAM-solar/phot <= 52.5 km 
 end module wamphys_set_merge_rad    
 
 module wamphys_set_data_solar
      use machine , only : kind_phys
      implicit none  
     
      integer,  parameter :: nwaves = 37
      integer,  parameter :: nwaves_euv = 22
      integer,  parameter :: nwaves_src =nwaves-nwaves
      integer,  parameter :: lyman_a_num  = nwaves+1-12   
      integer,  parameter :: nsp_euv = nwaves
      
      real(kind=kind_phys), parameter   ::  PCC=1.985E-25    !  wavelength/energy conversion SI E=hc/lamda
      real(kind=kind_phys), parameter   ::  eccentric=1.     !  eccentricity of Earth's orbit      

      real(kind=kind_phys)                    :: euv37(nsp_euv)
      real(kind=kind_phys), dimension(nwaves) :: sfmin, afac,rlmeuv
      real(kind=kind_phys), dimension(nwaves) :: csao, csao2, csan2, csao3
      real(kind=kind_phys), dimension(nwaves) :: csio, csio2, csin2
      real(kind=kind_phys), dimension(nwaves) :: csdo2, csdeo2
      real(kind=kind_phys), dimension(nwaves) :: rwpcc
      
      real(kind=kind_phys)                    :: unit_conv = 1.e-18                     ! xsections into cm^2

      integer, parameter                      :: nz_euv = 17
      real(kind=kind_phys), dimension(nz_euv) :: effeuv17,effuv17, p17, z17
      
! o2 scale factor for 15 pressure level 
      integer, parameter                      ::  nz_jo2sf = 15
      real(kind=kind_phys)                    ::  jj_scale_factor(nz_jo2sf)
      real(kind=kind_phys)                    ::  z15(nz_jo2sf)
! tiegcm srb-efficiency factors on 63-pressures
      integer, parameter                      ::  nz_63 = 63
      real(kind=kind_phys)                    ::  pres63(nz_63), srbeff63(nz_63)
      real(kind=kind_phys)                    ::  z63(nz_63)
      
      data jj_scale_factor/ 1.0, 4.465536, 4.365480, 3.904985,  3.367959, &       
           3.202786, 2.378429, 1.636311, 1.423021, 1.452178, 1.588099,    &     
           1.714328,  1.811639, 1.907779, 1.987971/

      data effuv17/0.59, 0.59, 0.58, 0.57, 0.56, 0.52, 0.48, 0.43,        &
                   .4,.4,.4,.39,.34,.26,.19,.17,.16/

      data effeuv17/8*1.0,.75,.6,.62,.54,.49,.41,.33,.30,.30/
                         
      data srbeff63/1.000,1.000,1.000,1.000,1.000,      &                  
                  1.000,1.000,1.000,1.000,1.000,1.000,  &                             
                  1.000,1.000,1.000,.980,.983,.982,     &                               
                  .970,.945,.913,.880,.852,.832,.820,   &                              
                  .817,.816,.810,.795,.777,.765,.764,   &  
                  .759,.730,.664,.579,.500,.446,.416,   &  
                  .400,.393,.390,.390,.391,.391,.390,   &  
                  .388,.384,.380,.375,.366,.350,.324,   &  
                  .291,.260,.234,.214,.200,.190,.184,   &  
                  .180,.176,.173,.170/
	
      data pres63/6.90775528,    6.57442194,                & 
                  6.24108859,  5.90775525,  5.57442191,     &                           
                  5.24108856,  4.90775522,  4.57442188,     &                           
                  4.24108853,  3.90775519,  3.57442185,     &                           
                  3.2410885,   2.90775516,  2.57442182,     &                           
                  2.24108847,  1.90775513,  1.57442179,     &                           
                  1.24108844,  0.9077551,   0.574421757,    &                          
                  0.241088414, -0.0922449296,-0.425578273,  &                       
                  -0.758911616,-1.09224496,  -1.4255783,    &                          
                  -1.75891165, -2.09224499,  -2.42557833,   & 
                  -2.75891168, -3.09224502,  -3.42557836,   & 
                  -3.75891171, -4.09224505,  -4.42557839,   & 
                  -4.75891174, -5.09224508,  -5.42557842,   & 
                  -5.75891177, -6.09224511,  -6.42557845,   & 
                  -6.75891179, -7.09224514,  -7.42557848,   & 
                  -7.75891182, -8.09224517,  -8.42557851,   & 
                  -8.75891185, -9.0922452,   -9.42557854,   & 
                  -9.75891188, -10.0922452,  -10.4255786,   & 
                  -10.7589119, -11.0922453,  -11.4255786,   & 
                  -11.7589119, -12.0922453,  -12.4255786,   & 
                  -12.758912,  -13.0922453,  -13.4255787,   & 
                  -13.758912/
	         
!     TIEGCM data

      real                    :: dsfmin(nwaves), dafac(nwaves)
      real                    :: drlmeuv(nwaves) ! wavelengths (cm)
!                                           ! absoption cross sections (x1e18cm^2)
      real, dimension(nwaves) :: sigeuv_o, sigeuv_o2, sigeuv_N2, sigeuv_o3
!                                           !ionization branching ratios (off absorption)
      real, dimension(nwaves) ::  BphotonI_O, BphotonI_O2,BphotonI_N2
!                                           ! O2 dissociation branching ratios
      real, dimension(nwaves) ::  bro2DPh(nwaves)
!                                           !O2 Photo-e dissociation branching ratios
      real, dimension(nwaves) ::  bro2Del(nwaves)
!      
! from subroutine idea_solar_init data statements
!
      DATA drlmeuv/1.725e-05, 1.675e-05, 1.625e-05, 1.575e-05, &
		   1.525e-05, 1.475e-05, 1.425e-05, 1.375e-05, &
                   1.325e-05, 1.275e-05, 1.225e-05, 1.216e-05, &
                   1.175e-05, 1.125e-05, 1.075e-05, 1.038e-05, &
                   1.007e-05, 9.810e-06, 9.440e-06, 9.440e-06, &
                   9.440e-06, 8.555e-06, 8.555e-06, 8.555e-06, &
                   7.240e-06, 7.240e-06, 5.950e-06, 4.300e-06, &
                   3.050e-06, 2.570e-06, 1.895e-06, 1.125e-06, &
                   5.100e-07, 2.500e-07, 1.300e-07, 6.000e-08, &
                   2.250e-08/

! O2 absorption coefficient: (1,2,3,4)= (O2,O, N2,O3) 
! Cross Section cm^2 after scaling by 10-18 (unit_conv)
      DATA  sigeuv_O2 / &
		5.00e-01, 1.50e+00, 3.40e+00, 6.00e+00, 1.00e+01, &
                1.30e+01, 1.50e+01, 1.20e+01, 2.20e+00, 4.00e-01, &
                1.30e+01, 1.00e-02, 1.40e+00, 4.00e-01, 1.00e+00, &
                1.15e+00, 1.63e+00, 1.87e+01, 3.25e+01, 1.44e+01, &
                1.34e+01, 1.33e+01, 1.09e+01, 1.05e+01, 2.49e+01, &
                2.36e+01, 2.70e+01, 2.03e+01, 1.68e+01, 1.32e+01, &
                7.63e+00, 2.63e+00, 6.46e-01, 2.10e-01, 2.25e-01, &
                3.40e-02, 4.54e-03/

! O absorption coefficient:
      DATA  sigeuv_O /  & 
		0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 3.79e+00, 4.10e+00, 3.00e+00, 4.79e+00, &
                8.52e+00, 1.31e+01, 1.07e+01, 7.72e+00, 6.02e+00, &
                3.78e+00, 1.32e+00, 3.25e-01, 1.05e-01, 1.13e-01, &
                1.70e-02, 2.27e-03/

! N2 absorption coefficient:
      DATA  sigeuv_N2 /   &                                              
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 2.55e+00, 1.15e+02, 1.44e+01, &
                2.18e+00, 7.17e+01, 1.31e+01, 2.14e+00, 5.45e+01, &
                2.30e+01, 2.31e+01, 1.97e+01, 1.17e+01, 9.94e+00, &
                5.09e+00, 1.53e+00, 3.46e-01, 1.14e+00, 1.41e-01, &
                2.01e-02, 2.53e-03/
!
! O3 absorption cross section
     DATA sigeuv_o3/     &
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
               0.00e+00, 0.00e+00, 0.00e+00, 1.25e+01, 9.20e+00, &
               9.20e+00, 9.20e+00, 9.16e+00, 9.50e+00, 9.50e+00, &
               9.50e+00, 1.47e+01, 1.47e+01, 1.47e+01, 2.74e+01, &
               2.02e+01, 3.33e+01, 7.74e+00, 0.00e+00, 0.00e+00, &
               0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
               0.00e+00, 0.00e+00/
! 
! The three major species' ionization branching ratio (off absorption):
! O2
      DATA  BPhotonI_O2 / &                                              
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 6.13e-01, 8.30e-01, 6.20e-01, 7.86e-01, &
                7.56e-01, 5.34e-01, 5.74e-01, 5.49e-01, 4.76e-01, &
                6.73e-01, 9.83e-01, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00/

      DATA  BPhotonI_O /  &                                              
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00/
! N2
      DATA  BPhotonI_N2 /  &                                             
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 4.29e-01, &
                6.80e-01, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00/
! O2 photon dissociation branching ratio
      DATA  bro2DPh  /       &                                           
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 3.87e-01, 1.70e-01, 3.80e-01, 2.14e-01, &
                2.44e-01, 4.66e-01, 4.26e-01, 4.51e-01, 5.24e-01, &
                3.27e-01, 1.74e-02, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00/
! O2 photoelectron dissociation branching ratio
      DATA  bro2DEl  /      &                                            
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 1.10e-02, 6.53e-01, 7.62e-01, 9.96e-01, &
                1.27e+00, 2.04e+00, 4.11e+00, 5.70e+01, 1.78e+01, &
                2.03e+01, 8.79e+01/
     
! Solar spectrum based on EUVAC and glow for wave length less than 1050 A
! and Woods for wavelength greater than 1050 A
! solar minimum flux (when P_index=80, unit:photon cm^-2 S^-1)
   DATA  dsfmin /3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10, &
                 5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10, &
                 2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11, &
                 8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09, &
                 4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09, &
                 3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09, &
                 1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09, &
                 8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09, &
                 5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04, &
                 5.010e+01/
! scaling factor A as defined in EUVAC model
    DATA dafac /5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03, &
                1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03, &
                2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03, &
                2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03, &
                5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03, &
                4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03, &
                3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02, &
                7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03, &
                1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01, &
                6.240e-01/
				
 end module wamphys_set_data_solar    
 
 
!=======================================    
module wamphys_set_data_ion
   use machine , only : kind_phys
   implicit none
!=====================================================================     
!  input data on the magnetic grid >>>  iondata_file    = 'iondata_tjr.nc'
!=====================================================================        
      integer, parameter     :: nymag = 91, jcen =(nymag-1)/2+1    ! dimension of glat
      integer, parameter     :: nxmag = 20                         ! dimension of glon
      real (kind=kind_phys)  :: glat(nymag)
      real (kind=kind_phys)  :: glon(nxmag)
      
      real(kind=kind_phys), dimension(nxmag, nymag) :: cormag, btot, dipang 
      real(kind=kind_phys), parameter  ::   pi = 3.1459268, pi2 = 2.*pi        
      real(kind=kind_phys), parameter  ::   ddlat= pi/(nymag-1)   !3.4906585033333331e-002
      real(kind=kind_phys), parameter  ::   ddlon= pi2/nxmag 
      
      
!=====================================================================            
! input data  >>>  tirosdata_file    = 'tiros_tjr.nc'   
       
      integer, parameter :: nt_21 = 21            ! tiros-21  1st dimension
      integer, parameter :: nt_20 = 20   ! tiros-20  2nd dimension
      integer, parameter :: nt_7  =  7   ! tiros-7   3rd dimension
!
      real(kind=kind_phys), dimension(nt_21, nt_20, nt_7) :: emaps, cmaps     ! tiros

      integer, parameter     :: n_flx    = 15               ! tiros
      integer, parameter     :: n_bnd    = 21               ! tiros
      real(kind=kind_phys)   :: djspectra(n_flx, n_bnd)     ! tiros
!=====================================================================      
      integer, parameter :: jmaxwell = 6
      real(kind=kind_phys), parameter :: width_maxwell = 0.050

      real(kind=kind_phys)       :: en_maxwell(jmaxwell)               
      real(kind=kind_phys)       :: width(n_flx), en(n_flx)            ! defined by data-stat
      real(kind=kind_phys)       :: te11(n_bnd),te15(n_bnd)            ! precomp_iondata_fixed
!
      real ::  ratio(21)                                              !precomp_iondata_fixed
      real ::  ionchr(21), rlam(21)  
      real ::  ion_recomb(8),lognpres(8)                              ! defined by data-stat
!
      data en    /.37,.6,.92,1.37,2.01,2.91,4.19,6.,8.56,12.18, &
                  17.3,24.49,36.66,54.77,81.82/
      data width/.158,.315,.315,.63,.631,1.261,1.26,2.522,   &
                 2.522,5.043,5.043,10.,14.81,22.13,33.06/
		 
      data rlam/1.49,1.52,1.51,1.48,1.43,1.37,1.30,1.22, 1.12,1.01,0.895, &
                0.785,0.650,0.540,0.415,0.320,0.225, 0.14,0.08,0.04,0.0/

      data ionchr/.378 , .458 , .616 , .773 , .913 , 1.088 , 1.403 , &
             1.718 , 2.033 , 2.349 , 2.979 , 3.610 , 4.250 , 4.780 , &
             6.130 , 7.392 , 8.653 , 9.914 , 12.436 , 14.957 , 17.479/

       
      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,-10.507/

      data ion_recomb/3.20e-13,3.20e-13,2.75e-13,1.45e-13,1.13e-13, &
                      8.30e-14,3.70e-14,2.00e-14/   


end module wamphys_set_data_ion
 
module wamphys_set_data_vg150


! only : lev_wam,  h2ora150, o3ra150, prlog_70_150, prlog150

   use machine , only : kind_phys
   implicit none
   integer, parameter     :: lev_wam = 150 , lev80 = 80
   
   real(kind=kind_phys), dimension(lev80) ::  h2ora150, o3ra150, prlog_70_150  
   real(kind=kind_phys), dimension(lev_wam) ::  prlog150
!
! can be written as the data-file: h2ora150, o3ra150, prlog_70_150, prlog150
!     
! 71-150 in levs=150
! o3(71-150)
      data o3ra150/4.10541952E-06,3.47100766E-06,2.87068966E-06, 2.35683753E-06,& 
                   1.96476323E-06,1.68001584E-06,1.46059012E-06,1.28086944E-06, &      
                   1.12287103E-06,9.73440677E-07,8.31057093E-07,6.96823493E-07, &  
                   5.70485075E-07,4.54900920E-07,3.51380290E-07,2.59055385E-07, &  
                   1.83987938E-07,1.33985182E-07,9.93050813E-08,8.12517455E-08, &  
                   1.04879335E-07,1.96984693E-07,3.40876799E-07,5.63920720E-07, &  
                   8.83452184E-07,1.23309195E-06,1.61560931E-06,1.90510281E-06, &  
                   2.00312741E-06,1.98334669E-06,1.75853471E-06,1.44161553E-06, &  
                   1.11576928E-06,7.89776361E-07,5.25719302E-07,3.33307290E-07, &  
                   1.90201852E-07,9.50490959E-08,4.25181927E-08,1.71517381E-08, &  
                   6.31168787E-09,2.32353325E-09,2.00874504E-09,1.66279638E-09, &  
                   1.36930561E-09,1.12419760E-09,9.19829659E-10,7.49814512E-10, &  
                   6.08729657E-10,4.91976560E-10,3.95658450E-10,3.16475520E-10, &  
                   2.51635148E-10,1.98775150E-10,1.55898314E-10,1.21316667E-10, &  
                   9.36041570E-11,7.15566049E-11,5.41579312E-11,4.05517867E-11, &  
                   3.00178145E-11,2.19518467E-11,1.58493829E-11,1.12917456E-11, &  
                   7.93434889E-12,5.49657616E-12,3.75284443E-12,2.52454277E-12, &  
                   1.67265129E-12,1.09096239E-12,6.99914181E-13,4.41092526E-13, &  
                   2.72463151E-13,1.64366989E-13,9.62680762E-14,5.41996795E-14, &  
                   2.88221148E-14,1.39894852E-14,5.72118432E-15,4.70438733E-16/   
! h2o 71-150 in levs=150
      data h2ora150/4.15074772E-06,4.13699000E-06,4.11797890E-06, 4.09487986E-06,    &   
                    4.06858733E-06, 4.03597828E-06, 3.99688515E-06, 3.95067808E-06,  &      
                    3.89717454E-06, 3.83486354E-06, 3.76154928E-06, 3.67776509E-06,  &   
                    3.57952092E-06, 3.45696758E-06, 3.30616948E-06, 3.13086436E-06,  &    
                    2.91936568E-06, 2.64976784E-06, 2.33136751E-06, 1.97812350E-06,  &    
                    1.56715103E-06, 1.18281856E-06, 8.41511396E-07, 5.69260876E-07,  &    
                    3.88780697E-07, 2.50438515E-07,1.54300660E-07, 1.02009581E-07,   &    
                    6.65450034E-08, 4.17382808E-08, 2.82805186E-08, 2.01512556E-08,  &    
                    1.41564448E-08, 1.02806445E-08, 7.94408149E-09, 6.32637731E-09,  &    
                    5.12551203E-09, 4.27811892E-09, 3.70565449E-09, 3.31366890E-09,  &    
                    3.03512593E-09, 2.86004858E-09, 3.14079315E-09, 3.43411317E-09,  &    
                    3.75162719E-09, 4.09541203E-09, 4.46698364E-09, 4.86779007E-09,  &   
                    5.29913960E-09,5.76212751E-09, 6.25754557E-09, 6.78577268E-09,   &    
                    7.34664587E-09, 7.93931252E-09, 8.56206704E-09, 9.21217910E-09,  &    
                    9.88572497E-09, 1.05774394E-08, 1.12806111E-08, 1.19870504E-08,  &   
                    1.26871579E-08, 1.33701220E-08, 1.40242598E-08, 1.46374975E-08,  &    
                    1.51979607E-08, 1.56946171E-08, 1.61178886E-08, 1.64601425E-08,  &    
                    1.67159770E-08, 1.68822374E-08, 1.69577319E-08,1.69426375E-08,   &     
                    1.68375826E-08, 1.66423366E-08, 1.63538713E-08, 1.59631314E-08,  &    
                    1.54486221E-08, 1.47606491E-08, 1.37697900E-08, 6.83803988E-09/  
!
                   data prlog150/                                                   &   
                   -.010495013621173093,-.0047796645053569788,.0017317939011674947, &  
                     .0091445549523354423,.017575964483718530,.027156409259219756,  &     
                       .038029776798164390,.050354098813263921,.064301975566456532, &       
                        .080060604331725002,.097831430661753094,.11782928094771801, &        
                          .1402811398792534,.16542406580956714,.19350235130401269,  &          
                          .22476418567991183,.25945740055444533,.29782445285326098, &          
                          .34009706723024435,.38649056296119455,.43719785635654262, &          
                          .49238384410129205,.55218040859471984,.61668213110539583, &         
                          .68594338857491133,.75997679704427235,.83875307597091064, &          
                          .92220240072958426,1.0102172332720634,1.1026563113455261, &          
                          1.1993493125965171,1.3001020446943199,1.4047022170280321, &          
                          1.5129293855161694,1.6245648325694693,1.7393953508635152, &          
                          1.8572172092097674,1.9778421006194884,2.1011027392385673, &          
                          2.2268591758946630,2.3550060174351386,2.4854802051244675, &          
                          2.6182709767607482,2.7534197208872264,2.8909811818384683, &          
                          3.0309947994264324,3.1734952229335223,3.3185141840359855, &          
                          3.4660806887940865,3.6162208681691572,3.7689579941351727, &          
                          3.9243126000728270,4.0823023952198705,4.2429421890712113, &          
                          4.4062440367083493,4.5722171074490712,4.7408678113152334, &          
                          4.9121997598293738,5.0862136813095233,5.2629074813405223, &          
                          5.4422762630064341,5.6243123514382578,5.8090052129405310, &          
                          5.9963414712632250,6.1863048973029393,6.3788765414713984, &          
                          6.5740346053768723,6.7717544276443897,6.9720085218050194, & 
	                  7.1747665967259380,7.3799955384809381,7.5876594664001917, &        
                          7.7977196240960511,8.0101343615235372,8.2248592986929712, &          
                          8.4418472650864640,8.6610481622066615,8.8824091352082579, &          
                          9.1058744786331118,9.3313856279708016,9.5588812882688448, &         
                          9.7882972917212427,10.019566564224711,10.252619320659370, &          
                          10.487382934517129,10.723781889667228,10.961737871454984, &          
                          11.201169705433793,11.441993457513469,11.684122340664617, &          
                          11.927466717621803,12.171934093143912,12.417429143544323, &          
                          12.663853787366328,12.911107059366048,13.159085134752051, &          
                          13.407681405572704,13.656786420959016,13.906287867607862, &          
                          14.156070558238182,14.406016545469104,14.656016536601980, &         
                          14.906016520565826,15.156016549153117,15.406016561713120, &         
                          15.656016540294088,15.906016534936020,16.156016568807406, &        
                          16.406016534338949,16.656016529668964,16.906016577452895, &         
                          17.156016561083050,17.406016546123016,17.656016537621376, &         
                          17.906016506464592,18.156016518186942,18.406016541248878, &         
                          18.656016561606052,18.906016537947128,19.156016514271940, &       
                          19.406016513968229,19.656016519436953,19.906016540106812, &          
                          20.156016505341348,20.406016517771079,20.656016529319292, &          
                          20.906016501516998,21.156016508903946,21.406016525111486, &          
                          21.656016521424860,21.906016510407873,22.156016517462234, &          
                          22.406016507059238,22.656016506270586,22.906016547648477, &          
                          23.156016547534154,23.406016516157738,23.656016508155140, &          
                          23.906016513141179,24.156016519666711,24.406016497318937, &          
                          24.656016492779738,24.906016503952380,25.156016485200560, &          
                          25.406016484453687,25.656016491410533,25.906016498929535, &          
                          26.156016525642059,26.406016533196180,27.231955945328760/      
!          	               
	data prlog_70_150/7.3799955384809381,7.5876594664001917, &  
                          7.7977196240960511,8.0101343615235372,8.2248592986929712, &          
                          8.4418472650864640,8.6610481622066615,8.8824091352082579, &          
                          9.1058744786331118,9.3313856279708016,9.5588812882688448, &         
                          9.7882972917212427,10.019566564224711,10.252619320659370, &          
                          10.487382934517129,10.723781889667228,10.961737871454984, &          
                          11.201169705433793,11.441993457513469,11.684122340664617, &          
                          11.927466717621803,12.171934093143912,12.417429143544323, &          
                          12.663853787366328,12.911107059366048,13.159085134752051, &          
                          13.407681405572704,13.656786420959016,13.906287867607862, &          
                          14.156070558238182,14.406016545469104,14.656016536601980, &         
                          14.906016520565826,15.156016549153117,15.406016561713120, &         
                          15.656016540294088,15.906016534936020,16.156016568807406, &        
                          16.406016534338949,16.656016529668964,16.906016577452895, &         
                          17.156016561083050,17.406016546123016,17.656016537621376, &         
                          17.906016506464592,18.156016518186942,18.406016541248878, &         
                          18.656016561606052,18.906016537947128,19.156016514271940, &       
                          19.406016513968229,19.656016519436953,19.906016540106812, &          
                          20.156016505341348,20.406016517771079,20.656016529319292, &          
                          20.906016501516998,21.156016508903946,21.406016525111486, &          
                          21.656016521424860,21.906016510407873,22.156016517462234, &          
                          22.406016507059238,22.656016506270586,22.906016547648477, &          
                          23.156016547534154,23.406016516157738,23.656016508155140, &          
                          23.906016513141179,24.156016519666711,24.406016497318937, &          
                          24.656016492779738,24.906016503952380,25.156016485200560, &          
                          25.406016484453687,25.656016491410533,25.906016498929535, &          
                          26.156016525642059,26.406016533196180,27.231955945328760/
			  


end module wamphys_set_data_vg150

module wam_efieldw05_read_data
 
   use machine , only : kind_phys
      implicit none 
      integer, parameter  :: lu =99  
      integer , parameter ::   iulog = 6
      character(len=12), parameter :: model='epot'

      integer,  parameter :: csize=28, d1_pot=15, d2_pot=18
      integer             :: ab(csize), ls(csize), ms(csize)
      integer             :: maxl_pot,maxm_pot 
  
! Data read from SCHAtable.dat
      integer :: maxk_scha,maxm_scha

      integer,parameter      :: d1_scha=19, d2_scha=7, d3_scha=68

      real(kind=kind_phys)   :: allnkm(d1_scha,d2_scha,d3_scha)

      real(kind=kind_phys)   :: th0s(d3_scha) 
      real(kind=kind_phys)   :: alschfits(d2_pot,csize), schfits(d1_pot,csize), ex_pot(2)        
      ! Data read from W05scBndy.dat

      integer,parameter    :: na=6, nb=7
      real(kind=kind_phys) :: bndya(na),bndyb(nb),ex_bndy(2)
           
      !real(kind=kind_phys) :: rad2deg,deg2rad           ! set by SetModel_new
      real(kind=kind_phys) :: bndyfitr                  ! calculated by setboundary
      real(kind=kind_phys) :: esphc(csize),bsphc(csize) ! calculated by SetModel_new
      real(kind=kind_phys) :: tmat(3,3),ttmat(3,3)      ! from setboundary
       
      integer,parameter    :: mxtablesize=500 
      real(kind=kind_phys) :: plmtable(mxtablesize,csize),colattable(mxtablesize)
      real(kind=kind_phys) :: nlms(csize)

     contains
     
     subroutine read_potential(infile)
! Read ascii data file W05scEpot.dat or W05scBpot.dat, written by 
!   pro write_potential (write_data.pro)
            implicit none
! Args:
            character(len=*),intent(in) :: infile
! Local:
            character(len=16) :: fname
            integer :: i
            integer, parameter :: iulog = 6     
            integer :: csize_rd,d1_rd,d2_rd
      !
            open(lu,file=trim(infile),status='old', ACCESS ='SEQUENTIAL')
            read(lu,"(a)") fname
            read(lu,"(28i3)") ab
            read(lu,"(3i3)") csize_rd,d1_rd,d2_rd
            if (csize_rd /= csize) then
                  write(iulog,"('>>> read_potential: file ',a,': incompatable csize: ', &
                  'csize_rd=',i4,' csize=',i4)") fname,csize_rd,csize
                  stop 'csize'
            endif
            if (d1_rd /= d1_pot) then
                  write(iulog,"('>>> read_potential: file ',a,': incompatable d1: ', &
                  'd1_rd=',i4,' d1_pot=',i4)") fname,d1_rd,d1_pot
                  stop 'd1'
            endif
            if (d2_rd /= d2_pot) then
                  write(iulog,"('>>> read_potential: file ',a,': incompatable d2: ', &
                  'd2_rd=',i4,' d2_pot=',i4)") fname,d2_rd,d2_pot
                  stop 'd2'
            endif
            do i=1,csize
                  read(lu,"(6e20.9)") alschfits(:,i)
            enddo
            read(lu,"(2f10.3)") ex_pot
            read(lu,"(28i3)") ls
            read(lu,"(2i3)") maxl_pot,maxm_pot
            read(lu,"(28i3)") ms

            do i=1,csize
                  read(lu,"(6e20.9)") schfits(:,i)
            enddo
            close(lu)
      end subroutine read_potential
!-----------------------------------------------------------------------
! Read ascii data file SCHAtable.dat, 
      subroutine read_schatable(infile)
 
            implicit none
! Args:
            character(len=*),intent(in) :: infile
! Local:
            character(len=16) :: fname
            integer :: i,j
 
            open(lu,file=trim(infile),status='old', ACCESS = 'SEQUENTIAL')
            read(lu,"(a)") fname
            read(lu,"(2i3)") maxk_scha,maxm_scha
            do i=1,d3_scha
            do j=1,d2_scha
                  read(lu,"(6e20.9)") allnkm(:,j,i)
            enddo
            enddo
            read(lu,"(8f10.4)") th0s

!           print *, 'read_schatable', th0s

            close(lu)
      end subroutine read_schatable
!-----------------------------------------------------------------------
! Read ascii data file W05scBndy.dat, written by pro write_bndy
!   (write_data.pro)     
      subroutine read_bndy(infile)

            implicit none
! Args:
            character(len=*),intent(in) :: infile
! Local:
            character(len=16) :: fname
            integer :: rd_na,rd_nb
            integer, parameter :: iulog = 6 

            open(lu,file=trim(infile),status='old', ACCESS = 'SEQUENTIAL')
            read(lu,"(a)") fname
            read(lu,"(2i3)") rd_na,rd_nb
            if (rd_na /= na) then
                  write(iulog,"('>>> read_potential: file ',a,': incompatable na: ', &
                  'rd_na=',i4,' na=',i4)") fname,rd_na,na
                  stop 'na'
            endif
            if (rd_nb /= nb) then
                  write(iulog,"('>>> read_potential: file ',a,': incompatable nb: ', &
                  'rd_nb=',i4,' nb=',i4)") fname,rd_nb,nb
                  stop 'nb'
            endif
            read(lu,"(8e20.9)") bndya
            read(lu,"(8e20.9)") bndyb
            read(lu,"(8e20.9)") ex_bndy

!           print *, 'read_bndy', ex_bndy
            close(lu)     
      end subroutine read_bndy
                  
 end module wam_efieldw05_read_data
! 
 module wam_efield_setdef_data
 
   use machine , only : kind_phys
   use physcons , only : con_pi
   implicit none

!---------------------------------------------------------------------- 
!  wam_efield_init: index for factors f_m(mlt),f_l(UT),f_-k(d)
!----------------------------------------------------------------------
      real(kind=kind_phys), parameter :: r_e  =  6.371e6
      real(kind=kind_phys), parameter :: h_r  =  130.0e3
      real(kind=kind_phys), parameter :: dy2yr= 365.24
      real(kind=kind_phys), parameter :: dy2mo= 30.6001  
      real(kind=kind_phys), parameter :: pi = con_pi                     
      real(kind=kind_phys), parameter :: ef_max  = 0.015  ! max e-field for high latitude boundary location [V/m]
      real(kind=kind_phys), parameter :: lat_sft = 54.	  ! shift of highlat_bnd to 54 deg 
      
      real(kind=kind_phys), parameter :: EPOCH=1980.,TH0=11.19,PH0=-70.76, DIPOLE=.30574     
      real(kind=kind_phys), parameter :: bnd_wei = 44.         ! colat. [deg]
           
      integer,              parameter :: nmax_sin = 2 ! max. wave number to be represented      
             
      integer, parameter :: ni = 1091                         ! for n=12 m=-18:18

      integer,dimension(0:ni) :: kf,lf, mf, nf, jf
      real(kind=kind_phys)    :: ft(1:3,0:2)                     ! used for f_-k(season,k)

      real(kind=kind_phys) ::  a_klnm(0:ni)                   !  A_klm
      real(kind=kind_phys) ::  a_lf(0:ni)                     !  A_klmn^lf for minimum  
      real(kind=kind_phys) ::  a_hf(0:ni)                     !  A_klmn^hf for maximum
         
      integer, parameter   ::  nm   = 19,     mm   = 18    					
      integer, parameter   ::  nmp  = nm + 1, mmp  = mm + 1
      
      real(kind=kind_phys) :: r(0:nm,0:mm)                     ! R_n^m
      real(kind=kind_phys) :: pmopmmo(0:mm)                    ! sqrt(1+1/2m)
      
! fixed magnetic grids

      integer, parameter   ::   nmlon = 180,nmlonh= nmlon/2, nmlonp1 = nmlon+1      ! mlon 
      integer, parameter   ::   nmlat = 90, nmlath= nmlat/2, nmlatp1 = nmlat+1      ! mlat
      integer, parameter   ::   iseasav = 0  ! flag for season 
        
      real(kind=kind_phys) :: ylatm(0:nmlat), ylonm(0:nmlon), dlonm,  dlatm   
      real(kind=kind_phys) :: sinIm_mag(0:nmlat)
      real(kind=kind_phys) :: potent(0:nmlon,0:nmlat)
      real(kind=kind_phys) :: ed1(0:nmlon,0:nmlat), ed2(0:nmlon,0:nmlat)
      
      real(kind=kind_phys) ::  date, day, ut                                       ! iyear+iday+ut iday+ut  
      integer              ::  iday, iyear, imo, iday_m                           ! day of month 
      integer              ::  imax,  ilat_sft, jmin, jmax, nmlat_wei
      
      real(kind=kind_phys) :: ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP               ! defined in readcoef
      real(kind=kind_phys) :: CX(9),ST(6),CT(6),AM(3,3,11)  
      real(kind=kind_phys) :: rtd, dtr, sqr2, hr2rd, dy2rd,deg2mlt, mlt2deg         
!---------------------------------------------------------------------
      contains
!subroutine magnetic_grids        
!subroutine set_readcoef
!subroutine efread_acoef 
!subroutine read_acoef_efield
!subroutine pnm   
!subroutine prep_fk
!subroutine index_quiet
!subroutine ff
!subroutine prep_pnm

 subroutine magnetic_grids
      implicit none
      
      integer              :: i,j
      real(kind=kind_phys) :: fac, yy

      rtd     = 180./con_pi 	        ! radians -> deg
      dtr     = con_pi/180.	        ! deg -> radians
      sqr2    = sqrt(2.e0)
      hr2rd   = con_pi/12.	        ! pi/12 hrs
      dy2rd   = 2.*con_pi/dy2yr         ! 2*pi/365.24  average year
      deg2mlt = 24./360.                ! convert degrees to MLT hours
      mlt2deg = 360./24.                ! for mlt to mlon       

!-------------------------------------------------------------------
! Set grid deltas:
!-------------------------------------------------------------------
      dlatm = 180./nmlat
      dlonm = 360./nmlon

!-------------------------------------------------------------------
! Set magnetic latitude array 
!-------------------------------------------------------------------
      do j = 0,nmlat
            ylatm(j) = j*dlatm
            yy = (ylatm(j) - 90.)*dtr
            fac = cos(yy)               ! sinIm = 2*sin(lam_m)/sqrt[4-3*cos^2(lam_m)]
            fac = 4. - 3.*fac*fac
            fac = 2./sqrt( fac )
            sinIm_mag(j) = fac*sin( yy )
      end do 

!------------------------------------------------------------------
! Set magnetic longitude array
!------------------------------------------------------------------
      do i = 0,nmlon
        ylonm(i) = i*dlonm
      end do
      
!-----------------------------------------------------------------
! find boundary index for weimer
!------------------------------------------------------------------
      do j = 0,nmlat
        nmlat_wei = j
        if( bnd_wei <= ylatm(j) ) then
           exit
        end if
      end do 

!-------------------------------------------------------------------
! find latitudinal shift
!-------------------------------------------------------------------
      do j = 0,nmlat
        ilat_sft = j
        if( lat_sft <= ylatm(j) ) then
           exit
        end if
      end do 

!------------------------------------------------------------------
! find index for linear interpolation of ed2 at mag.equator 
! use 12 deg - same as in TIEGCM      
!------------------------------------------------------------------
      do j = 0,nmlat
        yy = ylatm(j) - 90.
        if( yy <= -12. ) then
	  jmin = j
        else if( yy > 12. ) then
	  jmax = j
	  exit
       end if
      end do      
 end subroutine magnetic_grids 
 
 subroutine set_readcoef
      implicit none
      real(kind=kind_phys) :: step, stpr, r2d
      r2d = 180./con_pi
      STEP = 10.
      STPR = STEP/6671.
      STPD = STPR * R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMN = 45.
      ALAMX = 90. - STPD
      ALAMR = ALAMN / r2d
      
 end subroutine set_readcoef
    
! prep_fk, prep_pnm,  index_quiet 
 subroutine prep_pnm    
! define array  R    real(kind=kind_phys) :: r(0:19,0:18)                     ! R_n^m  
      implicit none                

      integer  :: mp, m, n
      real(kind=kind_phys) :: xms, xns, den

      do mp = 1, mmp                      ! m+1 = 1,mm+1                                     
	m = mp - 1                                               
	xms = m*m                                                
	if( mp /= 1 ) then
           pmopmmo(m) = sqrt( 1. + .5/M )
	end if
	do n = m,nm                        ! n = m,N                                     
	  xns    = n*n                                       
	  den    = max(4.*xns - 1.,1.)
	  r(n,m) = sqrt( (xns  - xms)/den )
	end do                 
      end do 
      
 end subroutine prep_pnm
      
 subroutine ff( ph, mt, f )                                                    
!-----------------------------------------------------------------------
!Purpose: calculate F for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:  f_m(phi) = sqrt(2) sin(m phi) m > 0
!                   = 1                  m = 0
!                   = sqrt(2) cos(m phi) m < 0
!
! Author: A. Maute Nov 2003  am 11/18/03
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      integer,intent(in)   :: mt
      real(kind=kind_phys) ,intent(in)      :: ph	! geo. longitude of 0SLT (ut*15)
      real(kind=kind_phys) ,intent(out)     :: f(-mt:mt)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer                  :: m, i, j, mmo
      real(kind=kind_phys)     :: sp, cp    

      sp   = sin( ph/rtd )
      cp   = cos( ph/rtd )
      f(0) = 1.e0
                                                                
      f(-1) = sqr2*cp
      f(1)  = sqr2*sp      								 
      do m = 2,mt
        mmo   = m - 1  
        f(m)  = f(-mmo)*sp + cp*f(mmo)
        f(-m) = f(-mmo)*cp - sp*f(mmo)
      end do      

 end subroutine ff       
      
 subroutine index_quiet
!-----------------------------------------------------------------
!Arrays:  kf lf nf mf jf
!
! Purpose: set up index for factors f_m(mlt),f_l(UT),f_-k(d) to
!    describe the electric potential Phi for the empirical model   
!
! Method:
!    Phi = sum_k sum_l sum_m sum_n [ A_klmn * P_n^m *f_m(mlt)*f_l(UT)*f_-k(d)]
!    - since the electric potential is symmetric about the equator
!      n+m odd terms are set zero resp. not used
!    - in the summation for calculation Phi the index have the following
!      range n=1,12 and m=-n,n, k=0,2 l=-2,2
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------       

      implicit none

!----------------------------------------------------------------      
!	... local variables
!----------------------------------------------------------------                                                                   
      integer :: i, j, k, l, n, m

      i = 0 	! initialize
      j = 1 
      do k = 2,0,-1
        do l = -2,2
          if( k == 2 .and. abs(l) == 2 ) then
             cycle
          end if
          do n = 1,12
            do m = -18,18 
              if( abs(m) <= n ) then		    !  |m| < n
                if( (((n-m)/2)*2) == (n-m) ) then   ! only n+m even
             	  if( n-abs(m) <= 9 ) then	    ! n-|m| <= 9 why?
             	    kf(i) = 2-k
             	    lf(i) = l
             	    nf(i) = n
             	    mf(i) = m
             	    jf(i) = j
             	    i	  = i + 1	 ! counter
                  end if
                end if
              end if
            end do ! m
          end do ! n
        end do ! l
      end do ! k

      imax = i - 1  
      if(imax /= ni ) then    ! check if imax == ni 
!       write(iulog,'(a19,i5,a18,i5)') 'index_quiet: imax= ',imax,  &
!         ' not equal to ni =',ni 
        stop ' in index_quiet of wamphys_fix_hindcast ' 
      end if							
!     if(debug) write(*,*) ' index_quiet, imax=',imax

 end subroutine index_quiet 
      
 subroutine prep_fk
      
!-------------------------------------------------------------------
! Purpose: set up constants factors for f_-k(day) used for empirical model
!     to calculate the electric potential Author: A. Maute Nov 2003  am 11/19/03
!-------------------------------------------------------------------

      ft(1,0) = .75*sqrt( 6.e0 )/con_pi			
      ft(1,1) = 2.e0*ft(1,0)					      
      ft(1,2) = 1.e0						      
      ft(2,0) = ft(1,0) 					      
      ft(2,1) = -ft(1,1)					      
      ft(2,2) = 1.e0						      
      ft(3,0) = ft(2,1) 					      
      ft(3,1) = 0.						      
      ft(3,2) = 1.e0							   

 end subroutine prep_fk
      
 subroutine pnm( ct, p )
!----------------------------------------------------------------      
! Purpose: normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
! Method:
!   P_m^m    = sqrt(1+1/2m)*si*P_m-1^m-1                  m>0
!   P_n^m    = [cos*P_n-1^m - R_n-1^m*P_n-2^m ]/R_n^m     n>m>=0
!   dP/d phi = n*cos*P_n^m/sin-(2*n+1)*R_n^m*P_n-1^m/sin  n>=m>=0
!   R_n^m    = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!--------------------------------------------------------------------                                                                   

      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      real(kind=kind_phys), intent(inout) :: ct              ! cos(colat)                 
      real(kind=kind_phys), intent(inout) :: p(0:nm,0:mm)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n, np
      real(kind=kind_phys) :: pm2, st

!      ct = min( ct,.99 )		! cos(colat)
      st = sqrt( 1. - ct*ct ) 	        ! sin(colat)

      p(0,0) = 1.  
      do mp = 1,mmp  ! m+1=1,mm+1
        m = mp - 1
	if( m >= 1 ) then
           p(m,m) = pmopmmo(m)*p(m-1,m-1)*st 			
	end if
	pm2 = 0.                                                                  
	do n = mp,nm                    ! n=m+1,N
	   np     = n + 1
	   p(n,m) = (ct*p(n-1,m) - r(n-1,m)*pm2)/r(n,m)
	   pm2    = p(n-1,m)
        end do
      end do

 end subroutine pnm
      
 subroutine read_acoef_efield(efield_lflux_file, efield_hflux_file)
!----------------------------------------------------------------     
! Purpose:
!-----------------------------------------------------------------------
! low/midlatitude potential from Scherliess model
!-----------------------------------------------------------------------
!             global_idea_coeff_hflux.dat 
!             global_idea_coeff_lflux.dat
!    1. read in coefficients A_klmn^lf for solar cycle minimum and
!      A_klmn^hf for maximum 
!    2. adjust S_a (f107d) such that if S_a<80 or S_a > 220 it has reasonable numbers
!      S_aM = [atan{(S_a-65)^2/90^2}-a90]/[a180-a90]
!      a90  = atan [(90-65)/90]^2
!      a180 = atan [(180-65)/90]^2
!    3. inter/extrapolation of the coefficient to the actual flux which is
!      given by the user
!      A_klmn = S_aM [A_klmn^hf-A_klmn^lf]/90. + 2*A_klmn^lf-A_klmn^hf
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!---------------------------------------------------------------
      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file

      integer  :: i,ios,unit,istat

      integer, parameter  :: lu = 99
      character (len=256):: locfn

!------------------------------------------------------------------    
!  get coefficients file for solar minimum: 
!-----------------------------------------------------------------                                                                   
      unit     = lu


      locfn = efield_lflux_file

!------------------------------------------------------------------    
! open datafile with coefficients A_klnm
!------------------------------------------------------------------     


      open(unit=unit,file=trim(efield_lflux_file), &
          status = 'old',iostat = ios)
!     if(ios.gt.0) then
!     write(iulog,*) 
!    &'read_acoef: error in opening coeff_lf file',
!    &' unit ',unit
!       call endrun
!     end if

!----------------------------------------------------------------------------                                                                   
! read datafile with coefficients A_klnm
!--------------------------------------------------------------------   
!     write(iulog,*) 'read_acoef: read file ',trim(locfn),' unit ',unit

      read(unit,*,iostat = ios) a_lf
      
      close(unit)
!--------------------------------------------------------------------  
!  get coefficients file for solar maximum: 
!--------------------------------------------------------------------

      unit = lu
      locfn= efield_hflux_file

!-------------------------------------------------------------------
! open datafile with coefficients A_klnm
!------------------------------------------------------------------
!     write(iulog,*) 'read_acoef: open file ',trim(efield_hflux_file),' unit ',unit
      open(unit=unit,file=trim(efield_hflux_file), status = 'old',iostat = ios)
!     if(ios.gt.0) then
!      write(iulog,*) 
!    &'read_acoef: error in opening coeff_hf file',' unit ',unit
!       call endrun
!     end if

!-----------------------------------------------------------------
! read datafile with coefficients A_klnm
!----------------------------------------------------------------
!     write(iulog,*) 'read_acoef: read file ',trim(locfn)

      read(unit,*,iostat = ios) a_hf
      
!     if(ios.gt.0) then
!      write(iulog,*) 
!    &'read_acoef: error in reading coeff_hf file',' unit ',unit
!       call endrun
!     end if


      close(unit)

 end subroutine read_acoef_efield
  
 subroutine efread_acoef (efield_lflux_file, efield_hflux_file)
      
      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file

      integer  :: i,ios,unit,istat      
      character (len=256):: locfn
       unit= 10
       locfn=efield_lflux_file
        open(unit=unit,file=trim(locfn), status = 'old',iostat = ios)
        read(unit,*,iostat = ios) a_lf
       close(unit)
       
       locfn= efield_hflux_file
        open(unit=unit,file=trim(locfn), status = 'old',iostat = ios)
        read(unit,*,iostat = ios) a_hf
       close(unit)       
 end subroutine efread_acoef     
            
	                		    
end module wam_efield_setdef_data












!
! goto: idea_ion_empirmodels.f idea_solar_heating.f idea_merge_ipe_to_wam.f ideaca.f efield*f
! 10 subs       h2oc.f h2ohdc.f
!
! total 33 subs
!tracer_const_h.f
!namelist_wamphysics_def.f
!idea_cal_advance.f
!idea_co2.f
!idea_composition.f
!idea_dissipation.f
!idea_getno_snoe.f
!idea_h2o.f
!idea_imf_input.f
!idea_interpol_datab.f
!idea_io_units.f
!idea_ion.f
!idea_ion_empirmodels.f
!idea_ion_input.f
!idea_merge_ipe_to_wam.f
!idea_mpi_def.f
!idea_o2_o3.f
!idea_phys.f
!idea_solar_heating.f
!idea_solar_init.f
!idea_solar_input.f
!idea_tracer.f
!idea_tracers_input.f
!idea_trad_geopgrav.f
! ideaca.f
! h2oc.f  
! h2ohdc.f
! co2hc.f
!
! wamphys_init_module.F90
! wamphys_setdata.F90  !efield.f w05sc_efield_merge.f
!                      module wam_efieldw05_read_data
!                      module wam_efield_setdef_data
!                      sub read_acoef_efield sub efread_acoef
!
!wamphys_weimer2005.F90
!       w05sc_efield_merge.f
!       wamphys_weimer2005.F90:  use wam_efieldw05_read_data
!      wamphys_weimer2005.F90:  use wam_efield_setdef_data

! wamphys_math_interp.F90
! efield.f
!
!wamphys.F90:!    use efield_wam, only          :  iday,iyear,iday_m,imo 
