
module  wamphys_init_module
!   
!....................................................................................
!  order => radiation => wam_physics => gfs_physics
!...................................................................................
!
    use  machine,                only : kind_phys
    use  wamphys_common,         only : arad, pi, pi2, hpscale, hpskm 
    
    use  wamphys_set_data_vg150, only : lev_wam,  h2ora150, o3ra150, prlog_70_150, prlog150
    
    use  wamphys_set_data_solar, only : nwaves, sigeuv_o, sigeuv_o2, sigeuv_N2, sigeuv_o3
    use  wamphys_set_data_solar, only : BphotonI_O, BphotonI_O2,BphotonI_N2, bro2DPh, bro2Del    
    use  wamphys_set_data_solar, only : dsfmin, dafac,  drlmeuv, unit_conv, pcc, rwpcc 
    use  wamphys_set_data_solar, only : csao, csao2, csan2, csao3, csio, csio2, csin2, csdo2, csdeo2 

 
    use wamphys_swdef

    implicit none
    logical              :: module_is_initialized
    integer, parameter   :: iulog = 6
    
!    integer, parameter   :: lowipe_lev150 = 90
    
    character(len=255)   :: wamphys_nml2 ='wamphys_nml'
    character(len=255)   :: wamphys_swio='swio_waminput.nc'   
    
    logical            :: do_init_tides = .false.          ! initialize tidal diag (24-12 and 8-hr tides)
    logical            :: do_init_swio  = .true.           ! initialize read swio-input
    logical            :: do_init_wamtracers  = .true.     ! initialize read of fixed/input constiuents
    logical            :: do_init_wamion      = .true.     ! initialize read of the ionosphere input
    logical            :: do_init_wamrad      = .true.     ! initialize read of the wam-radiationinput  
    
    character(len=255)     :: wamrad  
! co2-para [vmrs .....150-720 ppmv]
! co2hc.f   :        open(30,file='global_idea_coeff_lte.150',status = 'OLD')subroutine co2cin
    character(len=255)     ::  file_lte150='global_idea_coeff_lte.150'
    character(len=255)     ::  file_lte360='global_idea_coeff_lte.360'  
    character(len=255)     ::  file_lte540='global_idea_coeff_lte.540'  
    character(len=255)     ::  file_lte720='global_idea_coeff_lte.720'  
          
!h2oc.f:        OPEN(11,FILE='global_idea_ggww_in4.par',STATUS='OLD')
!h2oc.f:        OPEN(71,FILE='global_idea_h2ort_kg7t.par',STATUS='OLD')
!h2oc.f:        OPEN(11,FILE='global_idea_ggww_in1.par',STATUS='OLD')
!h2oc.f:        OPEN(71,FILE='global_idea_h2ovb_kg7t.par',STATUS='OLD')

    character(len=255)     ::  file_ggww_in4 = 'global_idea_ggww_in4.par'
    character(len=255)     ::  file_ggww_in1 = 'global_idea_ggww_in1.par'
    character(len=255)     ::  file_h2ort = 'global_idea_h2ort_kg7t.par'
    character(len=255)     ::  file_h2ovb = 'global_idea_h2ovb_kg7t.par'        
  
    character(len=255)     :: efield_lflux_file = 'global_idea_coeff_lflux.dat'  
    character(len=255)     :: efield_hflux_file = 'global_idea_coeff_hflux.dat'   
    character(len=255)     :: W05SC_HAtbl_file  = 'global_idea_coeff_W05SCHAtable.dat'
    character(len=255)     :: W05SC_Bndy_file = 'global_idea_coeff_W05scBndy.dat'
    character(len=255)     :: W05SC_Bpot_file = 'global_idea_coeff_W05scBpot.dat'
    character(len=255)     :: W05SC_Epot_file = 'global_idea_coeff_W05scEpot.dat'          
    character(len=255)     :: iondata_file    = 'iondata_tjr.nc'
    character(len=255)     :: tirosdata_file  = 'tiros_tjr.nc' 
     
    character(len=255)     ::  efield_model   = 'epot'          ! or 'bpot'
    character(len=255)     ::  SPW_DRIVERS    = 'swpc_fst'
    character(len=255)     ::  SWIN_DRIVERS   = 'swin_wam'
!
!   tracers and NO-model
!
    character(len=255)     :: file_glob_tracers = 'wam_rm3_globcomp.nc' 
    character(len=255)     :: file_no_snoe = 'snoe_eof.nc'  
! solar-annual
!   wam_nems_imf_2012.nc   wasolar_dan_20161019.nc
!             

    real(kind=kind_phys) :: soro_gg = 0.      ! sea level for gravity in upper layers    
    real(kind=kind_phys) :: wam_JH0 =  1.75
    real(kind=kind_phys) :: wam_JH_tanh = 0.5 
    real(kind=kind_phys) :: wam_JH_semiann = 0.5 
    real(kind=kind_phys) :: wam_JH_ann =  0.0
    real(kind=kind_phys) :: wam_JH_sto =  25000.
    real(kind=kind_phys) :: wam_JH_st1 =  5000.
    real(kind=kind_phys) :: wam_skeddy0 =  140.
    real(kind=kind_phys) :: wam_skeddy_semiann = 60. 
    real(kind=kind_phys) :: wam_skeddy_ann = 0. 
    real(kind=kind_phys) :: wam_tkeddy_ann  =  0.
    real(kind=kind_phys) :: wam_tkeddy0 =  280.
    real(kind=kind_phys) :: wam_tkeddy_semiann = 0. 
    real(kind=kind_phys) :: wam_f107_fix   = 120.
    real(kind=kind_phys) :: wam_f107a_fix  = 120.
    real(kind=kind_phys) :: wam_kp_fix     = 3.0   
    real(kind=kind_phys) :: wam_gW_fix     = 30. 
    real(kind=kind_phys) :: wam_tiros_activity_fix =7. 
    
    
     real(kind=kind_phys) :: JH0 =  1.75
     real(kind=kind_phys) :: JH_tanh = 0.5 
     real(kind=kind_phys) :: JH_semiann = 0.5 
     real(kind=kind_phys) :: JH_ann =  0.0
     real(kind=kind_phys) :: JH_st0 =  25000.
     real(kind=kind_phys) :: JH_st1 =  5000.
     real(kind=kind_phys) :: skeddy0 =  140.
     real(kind=kind_phys) :: skeddy_semiann = 60. 
     real(kind=kind_phys) :: skeddy_ann = 0. 
     real(kind=kind_phys) :: tkeddy_ann  =  0.
     real(kind=kind_phys) :: tkeddy0 =  280.
     real(kind=kind_phys) :: tkeddy_semiann = 0. 
     real(kind=kind_phys) :: f107_fix   = 120.
     real(kind=kind_phys) :: f107a_fix  = 120.
     real(kind=kind_phys) :: kp_fix     = 3.0   
     real(kind=kind_phys) :: kpa_fix     = 3.0 
     real(kind=kind_phys) :: gW_fix     = 30. 
     real(kind=kind_phys) :: tiros_activity_fix =7. 
    
    integer  :: tiros_switch  = 0
    namelist /wamphys_nml/ JH0,      JH_tanh, JH_semiann, JH_ann, JH_st0, JH_st1,      &
                           skeddy0, skeddy_semiann, skeddy_ann,                        &
                           tkeddy0, tkeddy_semiann, tkeddy_ann,                        &
                           f107_fix, f107a_fix, kp_fix, kpa_fix,                       &
			   gw_fix, tiros_activity_fix, tiros_switch,                   & 
			   spw_drivers, efield_model, swin_drivers,do_init_tides,      &
			   do_init_swio, do_init_wamtracers, do_init_wamion, do_init_wamrad
!
!SPW_DRIVERS, SWIN_DRIVERS
! allocatable arrays, initilized during "cires_ugwp_init" &
!                     released   during "cires_ugwp_finalize"
!
   real(kind=kind_phys), allocatable :: kvg(:), ktg(:), krad(:), kion(:)
   real(kind=kind_phys), allocatable :: zkm(:), pmb(:), zlog(:)
   real(kind=kind_phys), allocatable :: rfdis(:), rfdist(:)
!
! wamphys-grid and global profiles defined on (1:levs, surf to top).
!
   real(kind=kind_phys), allocatable ::  pr_idea(:), prlog(:)
   real(kind=kind_phys), allocatable ::  h2ora(:),o3ra(:)
   real(kind=kind_phys), allocatable ::  gglob(:), muglob(:), zglob(:), tglob(:)
   
! wamphys-solar and radiation, global eff-profiles  

   real(kind=kind_phys), allocatable :: eff_hart(:), eff_hugg(:), eff_chap(:), eff_herz(:)   
   real(kind=kind_phys), allocatable :: eff_srb(:), eff_src(:), eff_lya(:) 
   real(kind=kind_phys), allocatable :: o2_scale_factor(:),  srbeff(:)   
   real(kind=kind_phys), allocatable ::  effeuv(:), effuv(:)
   
    integer            ::  nps       ! start WAM radiation level 71 for 150L-grid
! wamphys NO-SNOE model variables
   
   integer, parameter  ::  no_ny33=33
   integer, parameter  ::  no_nz16=16
   integer, parameter  ::  no_neofs=7
   
   real(kind=kind_phys), allocatable   ::  no_eof(:, : ,:)
   real(kind=kind_phys) , allocatable  ::  no_m(:,:)
   real(kind=kind_phys), allocatable   ::  no_zkm(:), no_mlat(:) 
!
! co2/h2o
!   
    real(kind=kind_phys), allocatable   :: co2my(:)
    real(kind=kind_phys),allocatable,dimension(:) ::  & 
         gh2ort,gh2ovb,dg1rt,dg2rt, dg1vb,dg2vb,gdp,xx,wvmmrc,coeff      
!=============================================================
! Below is the specific GSMWAM-150L  vertical grid markers in km
! specified on wam_150L levels
!=============================================================  

   integer              ::   nlev_h2o,nlevc_h2o,nlev_co2
   integer              ::   k41,k71,k110,k105,k100,k43
   integer              ::   k91,k47,k64,k81,k87
			  
!=====================================================================     
! external global tracer, temp, geop, and molweight vertical profiles
! char names(var) ;names:long_name = "names: h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm
			  
   real (kind=kind_phys), allocatable :: vmr_wam(:, :), wam_hya(:), wam_hyb(:)
   
   integer, parameter  :: ind_zg =15, ind_mu = ind_zg-1, ind_tg = ind_mu -1

   real (kind=kind_phys), allocatable :: wam_oh(:), wam_ho2(:)    ! simplified chemistry
         
   real (kind=kind_phys)  :: bz, rbz, amo,amn2, amo2, amo3, amh2o, amno
   real (kind=kind_phys)  :: rmo, rmo2, rmn2, rmh2o,rmo3
   		      
      real (kind=kind_phys), parameter:: muo =3.9e-7    ! viscosity coefficient
!                                                          of O (kg/m/s) 
      real (kind=kind_phys), parameter:: muo2=4.03e-7   ! viscosity coefficient
!                                                          of O2 (kg/m/s) 
      real (kind=kind_phys), parameter:: mun2=3.43e-7   ! viscosity coefficient
!                                                         of N2 (kg/m/s) 
      real (kind=kind_phys), parameter:: lao =75.9e-5   ! thermal conductivity
!                                                      coefficient of O (W/m/K)
      real (kind=kind_phys), parameter:: lao2=56.e-5    ! thermal conductivity
!                                                      coefficient of O2(W/m/K)
      real (kind=kind_phys), parameter:: lan2=56.e-5    ! thermal conductivity
!                                                      coefficient of N2(W/m/K)
      real (kind=kind_phys), parameter:: cpo =2.5       !specific heats of o
      real (kind=kind_phys), parameter:: cpo2=3.5       !specific heats of o2
      real (kind=kind_phys), parameter:: cpn2=3.5       !specific heats of n2
      
      real (kind=kind_phys), parameter::  pref = 1.e5, rpref =1./pref  
                
   contains
   
!   call wamphys_init_all(me, master, nlunit, logunit, jdat,input_nml_file, fn_nml2,   &
!              lonr, latr, levs, ak, bk, dtp, errmsg, errflg) 
! 
!    
!-----------------------------------------------------------------------------------

  subroutine  wamphys_init_all(me, master, nlunit, logunit, jdat_gfs, &
       fn_nml23, fn_nml2, lonr, latr, levs, ak, bk, dtp, errmsg, errflg)
!
!
    use  netcdf   
    use  wamphys_set_merge_rad,  only : npsrad,  prdot02       != 1.e5*exp(-xbl)   hpa  

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr
    integer, intent (in) :: jdat_gfs(8)
    real(kind=kind_phys),    intent (in) :: ak(levs+1), bk(levs+1)
    real(kind=kind_phys),    intent (in) :: dtp   
    
    character(len=64), intent  (in) :: fn_nml23   
    character(len=64), intent (in) :: fn_nml2
    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg  
    
    integer :: ios
    logical :: exists
    
    integer :: ncid,  iernc, vid, dimid, status
    integer :: k
    
    integer :: ddd_wam,    curday_wam
!    real(kind=kind_phys) ::  pref
    real(kind=kind_phys), allocatable ::  prpa(:)  
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
!    pref = 1.e5
    


    if (me == master) print *, ' READ input.nml fnml2 '
    if (me == master) print *, trim (fn_nml2), ' namelist in wamphys_init'
    inquire (file =trim (fn_nml2) , exist = exists)
!
    if (.not. exists) then
       if (me == master) write(6,'(3a)') 'wamphys_init: namelist file: ', trim (fn_nml2), ' does not exist'
        errflg = 1
       if (me == master) print *, ' wamphys_nml errflg=', errflg, trim (fn_nml2)
       stop ' not. exists fn_nml2 for wamphys_nml '
    else
        open (unit = nlunit, file = trim(fn_nml2), action = 'read', status = 'old', iostat = ios)
    endif
    
    rewind (nlunit)
    read   (nlunit, nml = wamphys_nml)
    close  (nlunit)

     
    curday_wam = jdat_gfs(1)*10000 + jdat_gfs(2)*100 +jdat_gfs(3)    
    call wamphys_day_of_year(jdat_gfs(1), jdat_gfs(2), jdat_gfs(3), ddd_wam)
      
! write version number and namelist to log file
    if (me == master) then

	
        write (6, *) " ================================================================== "
        write (6, *) "CCPP wamphys_namelist from ", trim (fn_nml2)
        write (6, nml = wamphys_nml)		
        write (6, *) " ================================================================== "	
        write (6, *) "calendar_ugwp ddd_ugwp=", ddd_wam
        write (6, *) "calendar_ugwp curday_ugwp=", curday_wam
        write (6, *) " ================================================================== "	
        write (6, *) 	 ddd_wam, ' jdat_gfs ddd of year '	
    endif
!
! Allocate fixed WAM vertical profiles on the mid-layer grid
!         
    allocate( pr_idea(levs),   prlog(levs)  )
    allocate( h2ora(levs),     o3ra(levs)   )
    allocate( zkm(levs),   pmb(levs), zlog(levs))
!    allocate( prsilvl(levs) )
    allocate(gglob(levs), muglob(levs), zglob(levs), tglob(levs))
!  
123 format(i4,4(2x, e10.3))  
   do k=1, levs
!       pmb(k)   = ak(k)*100. + pref*bk(k)         ! Pa -unit  Pref = 1.e5 pa, pmb = Pa
       
       pmb(k)   = (ak(k)+ak(k+1))*.5 + pref*(bk(k)+bk(k+1))*.5
       
       pr_idea(k) = 0.01* pmb(k)             ! in mb or hPa 
       zkm(k)   = -hpskm*alog(pmb(k)*rpref)
       zlog(k)  = -alog(pmb(k))          ! unknown Pressure-units of IDEA ????? for eff euv/uv ??
       prlog(k) = alog(1000./pr_idea(k))
!       write(6,123)  k , ak(k), ak(k)*100,  pmb(k), bk(k)
       
!===========================       
! use wamphys_set_merge_rad, only : npsrad define the 1-st WAM-radiation level 
!                            prdot02 =22 Pa or 0.2 mb zlopgp7 = 58.9 km Hp = 0.99*8.5
!===========================
    enddo   
   
      do k=1,levs
        if(pmb(k).le.prdot02) then
          npsrad=k
          exit
        endif
      enddo  

    
 if (me == master) then 
       print *, ' npsrad ', npsrad,   pmb(npsrad), prdot02 
        print *  
       print *, 'pref, rpref ',    pref, rpref
   do k=1, levs   
   
   write(6,121) k, ak(k), pmb(k), bk(k), pr_idea(k), prlog(k)
   enddo   
   print *
   print *, ' zlog ', zlog(1), zlog(levs)
   print *, ' prlog ', prlog(1), prlog(43), prlog(levs)     
 endif       
121 format(i4, 5(2x, e10.3))


      call wam_composition_init150(levs, me, master)  ! defines k71, k43 etc.. 
      if (me == master) print *, ' wamphys_init_module after wam_composition_init150 '
      
      call get_global_grav_mu_zgeo (levs, me,master)
      
       if (me == master) print *, ' wamphys_init_module after get_global_grav_mu_zgeo'  
         
!   call gravco2(levs, philco2, soro_gg, gg) 
! check dimension of pmb  ..... sounds like Pa
!
!wamphys_radsolar_init.F90:! efficiencies suggested by Mlynchack and Solomon for UV-heating rates 40-110 km
!
! init solar-radiation efficiencies
!   
 if (do_init_wamrad) then 
 
   allocate(co2my(levs))    
   call wam_co2cin(prlog(k43),pmb(k43), muglob(k43),gglob(k43),co2my(k43), &
                    nlev_co2,me,master) 
   co2my(1:k43-1) = co2my(k43)		      

   if (.not.allocated(prpa)) allocate(prpa(nlevc_h2o))  
   prpa(1:nlevc_h2o) = pmb(k71:levs)
   
!nlevc_h2o=levs-k71+1   

   allocate(gh2ort(levs), gh2ovb(levs),dg1rt(levs),dg2rt(levs)  )
   allocate( dg1vb(levs),dg2vb(levs),gdp(levs))
   allocate( xx(levs),wvmmrc(levs),coeff(levs))
   
   call wam_h2ocin(prpa, nlevc_h2o,  me, master, &
        gh2ort(k71),gh2ovb(k71),dg1rt(k71),dg2rt(k71), &
        dg1vb(k71),dg2vb(k71),gdp(k71),xx(k71),wvmmrc(k71),coeff(k71) )
	
!    if(me == master) then
!      print*, ' wam_h2ocin dims ', levs, k71, nlevc_h2o
!      do k=k71, levs  
!       print *, ' h2ocin ', k, xx(k), gh2ort(k), dg1rt(k)
!      enddo
!      print*, ' h2ocin ++++++++ ', me
!!    endif  
   deallocate (prpa) 
   
!      if(me == master) stop ' wam_h2ocin '
   allocate (effeuv(levs), effuv(levs))  
    allocate (eff_hart(levs), eff_hugg(levs), eff_chap(levs), eff_herz(levs))
    allocate (eff_srb(levs), eff_src(levs), eff_lya(levs))
    allocate (o2_scale_factor(levs),  srbeff(levs) )
    allocate (wam_oh(levs),  wam_ho2(levs) )         
!    
! old-WAM use single "ef-array" for O2-O3 UV-radiation
!   
  
      eff_srb(1:levs) =1.00
      eff_lya(1:levs) =0.95
      eff_src(1:levs) =0.85   ! should be modified to introduce vertical profile with min at ~90km
      eff_hugg(1:levs)=1.00
      eff_chap(1:levs)=1.00
      eff_herz(1:levs)=1.00
      eff_hart(1:levs)=1.00   ! should be modified to suppress ozone mesospheric heating 
      srbeff(1:levs)  =1.00  
           
      call heat_uv_eff(levs, eff_hart, eff_src, me, master)                        !from idea_o2_o3.f  
      
      call get_eff_zlog( levs, zlog, effuv, effeuv,  o2_scale_factor, srbeff, me, master) 
!
      call wam_readno_snoewx(file_no_snoe, me, master)
      
      call wam_solar_radphoto_init(me, master) 
      
 endif    ! init wam radiation     
       
    if (do_init_wamtracers) call wam_tracer_init(levs, me, master)

 if (do_init_wamion) then 
 
      call wam_ion_init(levs, me, master) 
      efield_model = 'epot'
      call wam_efield_init(efield_model, me, master)
      
 endif    ! init wam-ionosphere          
!
! FV3 has ca-adjustments  
!           call ideaca_init(prsilvl,levs+1)
   
!
!======================
    module_is_initialized = .true.
    
    if (me == master) print *, ' wamphys_init is initialized ', module_is_initialized

    end subroutine wamphys_init_all
!
!    
    subroutine get_global_grav_mu_zgeo (levs, me, master)   ! defined in module gglob, muglob, zglob, soro_gg 
     use  netcdf   
     use  physcons , only : re => con_rerth, grav => con_g
     use wamphys_math_interp, only : idea_comp150_interp
     implicit none
     integer, intent(in)             :: levs, me, master    
     
!=======================================================================================================    
!/scratch1/NCEPDEV/swpc/Valery.Yudin/save/BASE_SVN/BASE_WAM_DATA/WAM_COMP/wam_rm3_globcomp.nc  
!        lev = 150 ;
!        var = 15 ;
!variables:char names(var) ;names:long_name = "names: h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm
!        float whyam(lev) ; whyam:long_name = "hybrid A Coefficient at layer midpoints" ;
!        float whybm(lev) ; whybm:long_name = "hybrid B Coefficient at layer midpoints" ;
!        char names(var) ;names:long_name = "names: h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm" ;
!        float vmr_wam(var, lev) ;vmr_wam:long_name = " tracers global vertical profiles in 1/m3 "
!=======================================================================================================  
     
    integer :: ncid,  iernc, vid, dimid, status         
    integer :: k, nlev, nvar
!
    real(kind=kind_phys), allocatable :: pdata(:), tdata(:), zdata(:), mudata(:) 
    real(kind=kind_phys) :: g0re2       
    iernc=NF90_OPEN(trim(file_glob_tracers), nf90_nowrite, ncid)
    g0re2 = grav*re*re
    if(iernc.ne.0) then         			    
       print *, 'cannot file_glob_tracers file=',trim(file_glob_tracers)			    
          return
        else
	
       status = nf90_inq_dimid(ncid, "lev", DimID)
!      if (status /= nf90_noerr) call handle_err(status)
!
       status = nf90_inquire_dimension(ncid, DimID,  len = nlev)      
       status = nf90_inq_dimid(ncid, "var", DimID)
       status = nf90_inquire_dimension(ncid, DimID,  len =nvar )
       
           if (me == master)  print *, nlev, nvar, ' nlev-nvar dimd of file_glob_tracers'
	   if (nlev .le. 0 .or. nvar .le. 0) then 
	       print *, 'tracer-file=',    trim(file_glob_tracers)	   
	       print *, ' nlev =', nlev, 'nvar=',nvar
	       stop
	   endif
  	   
        if (.not.allocated(vmr_wam))  allocate (vmr_wam(nlev, nvar))
        if (.not.allocated(wam_hya))  allocate (wam_hya(nlev)) 
        if (.not.allocated(wam_hyb))  allocate (wam_hyb(nlev)) 	        

        iernc=nf90_inq_varid( ncid, 'vmr_wam', vid )
        iernc= nf90_get_var( ncid, vid, vmr_wam)
        iernc=nf90_inq_varid( ncid, 'whyam', vid )
        iernc= nf90_get_var( ncid, vid, wam_hya)
        iernc=nf90_inq_varid( ncid, 'whybm', vid )
        iernc= nf90_get_var( ncid, vid, wam_hyb)		
  	iernc=nf90_close(ncid)
	
	endif 
	if (levs == nlev ) then
	 zglob(:) = vmr_wam(:, ind_zg)*1.e3
	 tglob(:) = vmr_wam(:, ind_tg)
	 muglob(:) = vmr_wam(:, ind_mu)
	 gglob(:)  = g0re2/((re+zglob(:))*(re+zglob(:)))
	 if (me ==  master) then
	 print *,  ' ggglob ', minval(gglob), maxval(gglob)
	 print *,  ' zglob ', minval(zglob), maxval(zglob) 
	 print *,  ' tglob ', minval(tglob), maxval(tglob) 
	 endif	 
	else
!	 print *,  ' levs =/= 150 ', levs, nlev	
	 allocate(pdata(nlev), tdata(nlev), zdata(nlev), mudata(nlev)) 
	 pdata(:) = -alog(wam_hyb +wam_hya)
	 zdata(:) = vmr_wam(:, ind_zg)
	 tdata(:) = vmr_wam(:, ind_tg) 
	 mudata(:) = vmr_wam(:, ind_mu)
	 call idea_comp150_interp(zdata, pdata, nlev, prlog, zglob, levs)
	 call idea_comp150_interp(tdata, pdata, nlev, prlog, tglob, levs) 
	 call idea_comp150_interp(mudata, pdata, nlev, prlog, muglob, levs) 	
	 gglob(:)  = g0re2/((re+zglob(:)*1.e3)*(re+zglob(:)*1.e3)) 	  
!	
! different vertical model grids
!	interpolate zdata => zglob
!     
	endif   
	 if (me ==  master) then
	 print *,  ' gglob ', minval(gglob), maxval(gglob)
	 print *,  ' zglob ', minval(zglob), maxval(zglob) 
	 print *,  ' tglob ', minval(tglob), maxval(tglob) 
	 print *,  ' muglob ', minval(muglob), maxval(muglob)
	 endif	
	 
	 if (allocated (pdata)) then 
	    deallocate (pdata)
	    deallocate (tdata)
	    deallocate (zdata)	
	    deallocate (mudata)	
	 endif   
	 
	 if(allocated(wam_hya)) then
	   deallocate( vmr_wam)
	   deallocate(wam_hya)
	   deallocate(wam_hyb)
	 endif 
	        	    
	return  
    end  subroutine get_global_grav_mu_zgeo    
    
  subroutine wam_solar_radphoto_init(me, master)    
!
     implicit none
     integer, intent(in)                 :: me, master   
     integer :: j, jinv
     
        do j = 1, nwaves  
        jinv = nwaves+1-j
!
         csao(j)   = sigeuv_o(jinv)*unit_conv
         csao2(j)  = sigeuv_o2(jinv)*unit_conv
         csao3(j)  = sigeuv_o3(jinv)*unit_conv
         csan2(j)  = sigeuv_n2(jinv)*unit_conv
	 
         csio(j)   = csao(j) *bphotoni_o(jinv)
         csio2(j)  = csao2(j)*bphotoni_o2(jinv)
         csin2(j)  = csan2(j)*bphotoni_n2(jinv)
         csdo2(j)  = csao2(j)*bro2dph(jinv)
         csdeo2(j) = csio2(j)*bro2del(jinv)
	 
         rwpcc(j)  = pcc/(drlmeuv(jinv)*1.e-2)      
         enddo
!
      RETURN  
  end subroutine wam_solar_radphoto_init
  
        
    subroutine wam_composition_init150(levs, me, master)
    
    use wamphys_math_interp, only : idea_comp150_interp    
    
    implicit none
    integer, intent(in)                 :: levs             ! number of pressure levels
    integer, intent(in)                 :: me, master    
!   

    integer             :: k, npn, np        
    integer             :: mpn, mp	
        if(levs.eq.150) then
!          prlog150 = prlog(1:levs)	 ! -log(P150_Pa)
          k41=41
          k110=110
          k71=71
          k105=105
          k100=100
! co2
          k43=43
! ion
          k91=91
! merge
          k47=47
          k64=64
          k81=81
          k87=87
	  
        else
!
! different V-grid from "levs=150 GSMWAM grid"
!       prlog(k) = alog(1000./pr_idea(k))
!	
          k71=levs
          k81=levs
          k87=levs
          k91=levs
          k100=levs
          k105=levs
          k110=levs
	  k41=levs
	  k64 = levs
	  k43 =levs
          do k=1,levs-10
          if(prlog(k).ge.prlog150(41).and.prlog(k-1).lt.prlog150(41))   k41=k
          if(prlog(k).ge.prlog150(71).and.prlog(k-1).lt.prlog150(71))   k71=k
          if(prlog(k).le.prlog150(110).and.prlog(k+1).gt.prlog150(110)) k110=k
          if(prlog(k).ge.prlog150(100).and.prlog(k-1).lt.prlog150(100)) k100=k
          if(prlog(k).le.prlog150(105).and.prlog(k+1).gt.prlog150(105)) k105=k
          if(prlog(k).ge.prlog150(43).and.prlog(k-1).lt.prlog150(43))   k43=k
          if(prlog(k).ge.prlog150(91).and.prlog(k-1).lt.prlog150(91))   k91=k
          if(prlog(k).ge.prlog150(47).and.prlog(k-1).lt.prlog150(47))   k47=k
          if(prlog(k).ge.prlog150(64).and.prlog(k-1).lt.prlog150(64))   k64=k
          if(prlog(k).ge.prlog150(81).and.prlog(k-1).lt.prlog150(81))   k81=k
          if(prlog(k).ge.prlog150(87).and.prlog(k-1).lt.prlog150(87))   k87=k
          enddo
        endif
	
          nlev_h2o=k110-k41+1
          nlevc_h2o=levs-k71+1
          nlev_co2=levs-k43+1
	  nps = k71
	  
      if(levs.eq.150) then
      
          h2ora(k71:levs)=h2ora150
          o3ra(k71:levs)=o3ra150
	  
      else
!
! needs extra-work and check
!      
           
	   npn =150 
	   np = npn-71+1
	   mpn = levs
	   mp = mpn-k71+1
          call idea_comp150_interp(h2ora150,prlog150(71:npn), np, prlog(k71:levs), h2ora(k71:levs), mp)
          call idea_comp150_interp(o3ra150, prlog150(71:npn), np, prlog(k71:levs), o3ra(k71:levs), mp)
	  
      endif	
	  h2ora(1:k71-1)=0. ! or constant values  h2ora(k71) > 
	   o3ra(1:k71-1)=0. !  o3ra(k71) replace them to the global WX-profiles of H2O and O3    
!	  
! better to use the Global H2O and O3 profiles (WX  with SABER/MLS compilation	  	  
!          call idea_interp(o3ra150, 71,150,80, o3ra, levs)
         
      if ( me == master ) then
         print *, ' wam_composition_init150 - performed '
!   integer              ::   nlev_h2o,nlevc_h2o,nlev_co2
!   integer              ::   k41,k71,k110,k105,k100,k43
!   integer              ::   k91,k47,k64,k81,k87

        print *, ' wam_comp k41,k71,k110,k105,k100,k43	 ', &
	        k41,k71,k110,k105,k100,k43
	print *	
        print *, ' wam_comp k91,k47,k64,k81,k87 ', k91,k47,k64,k81,k87
	print *		
	print *,' wam_comp nlev_h2o/co2', nlev_h2o,nlevc_h2o,nlev_co2
	print *, 'wamphys-prlog', prlog(1), prlog(43), prlog(levs)
	do k=1,levs
	 write(6,333) k, prlog(k), o3ra(k)*1.e6, h2ora(k)*1.e6
	enddo 
	print *, 'wamphys-prlog', prlog(1), prlog(43), prlog(levs)		
!	pause	' wam_composition_init150 '	
      endif
 333  format(i4, 3(2x, E10.3))     
    end subroutine wam_composition_init150
    
!=============================================
  subroutine wam_tracer_init(levs, me, master)
  
  use wamphys_math_interp, only :  z65toz
  use wamphys_setfix_tracers, only : np, ohi, ho2i
  use physcons, only : avgd => con_avgd
    implicit none
    integer, intent(in)                 :: levs                   ! number of model levels
    integer, intent(in)                 :: me, master 
    real(kind=kind_phys)    ::   rkgavgd, mo, mo2, mo3, mn2, mh2o
!
! initialize wam_oh(levs) & wam_ho2(levs) from some-global profiles (z65)
!      
      call z65toz(np, levs,ohi, zlog, wam_oh, 0.)
      call z65toz(np, levs,ho2i,zlog, wam_ho2,0.)
!
! assign set of costants employed in idea_tracer
!     use idea_composition, only   : bz,  amo,amn2,  amo2, amo3, amh2o
!     use idea_composition, only   : rbz, rmo, rmo2, rmn2, rmh2o,rmo3
      bz=1.3806505e-23
      rbz = 1./bz
      
      amo  =15.9994
      amo2 = 2.* amo
      amo3 = 3.* amo
      amh2o  = 18.0154
      amn2=28.013 
      amno=30.0061        
       
       rkgavgd =1.e-3/avgd
       mo=amo *rkgavgd
       mo2=amo2 *rkgavgd
       mn2=amn2 *rkgavgd
       mh2o=amh2o *rkgavgd
       mo3=amo3 *rkgavgd

       rmn2  = 1./mn2
       rmo   = 1./mo
       rmo3  = 1./mo3
       rmo2  = 1./mo2
       rmh2o = 1./mh2o
    
  end subroutine wam_tracer_init
  
!=============================================  
  subroutine wam_ion_init(levs, me, master)
!
! read and prepare time-invariant input ion-data
!  
    use netcdf
    use  wamphys_set_data_ion, only : emaps, cmaps, djspectra
    use  wamphys_set_data_ion, only : en, width, ratio, te15, te11,en_maxwell,jmaxwell, width_maxwell
    use  wamphys_set_data_ion, only : cormag, btot, dipang, glon, glat     
    implicit none
    integer, intent(in)                 :: levs                   ! number of model levels
    integer, intent(in)                 :: me, master   
    
    integer ::  ierr
    integer ::  ncid, vid, ierNC
    integer ::  astat
    integer ::  m, iband, j, iflux
!
! ASCII-files/arrays of WAMGSM-ion emaps1/cmaps1/djspectra1
!    
     real(kind=kind_phys)  :: emaps1(21,20,7),cmaps1(21,20,7),djspectra1(15,21)
     
     character (len=100) :: string_dum
     integer :: istat, i
     integer, parameter :: unit3=13
     integer, parameter :: unit7=17
!
! read emaps/cmaps/djspectra data and prepare arrays      subroutine precomp_iondata_fixed
!    
       ierNC=NF90_OPEN(trim(tirosdata_file), nf90_nowrite, ncid)   
        if (iernc /=0) write(6,*) ncid, 'ncid ', iernc, ' iernc '
!
        iernc=nf90_inq_varid( ncid, 'emaps', vid )
        iernc= nf90_get_var( ncid, vid, emaps)
        iernc=nf90_inq_varid( ncid, 'cmaps', vid )
        iernc= nf90_get_var( ncid, vid, cmaps)
        iernc=nf90_inq_varid( ncid, 'djspectra', vid )
        iernc= nf90_get_var( ncid, vid, djspectra)
	
       iernc=nf90_close(ncid) 
       
!emaps cmaps djspectra
! ascii-files for 
!
      close(unit=unit3)
      open(unit=unit3,file='ionprof',status='old',form='formatted', iostat=istat)
         read (unit3,99001) emaps1
         read (unit3,99001) cmaps1
       close(unit=unit3)
99001 format (1x,6e13.6)
      close(unit=unit7)
      open(unit=unit7,file='tiros_spectra',status='old', form='formatted', iostat=istat)

         read(unit=unit7,fmt=*) string_dum
         read(unit=unit7,fmt=*) string_dum
         read(unit7,fmt=*)
         do iband=1,21
         read(unit7,fmt=*) string_dum
         read(unit7,fmt=*)
         read(unit=unit7,fmt='(1x,5e10.4)')(djspectra1(iflux,iband),iflux=1,15)
         read(unit7,fmt=*)
         enddo
       close(unit=unit7)
       
!emaps cmaps djspectra
! ascii-files for 
!
	
      do iband=1,21
        ratio(iband) = (iband-1)*0.05     
        te15(iband)=0.0
        te11(iband)=0.0
       do m=1,15
         te15(iband)=te15(iband)+djspectra(m,iband)*en(m)*width(m)*1.6e-06
       enddo
       do  m=1,11
         te11(iband)=te11(iband)+djspectra(m,iband)*en(m)*width(m)*1.6e-06
       enddo         
      enddo                   ! iband=1,21
      do j = 1,jmaxwell
         en_maxwell(j) = j*width_maxwell - 0.025
      enddo  	
!
!
!
       iernc=nf90_open(trim(iondata_file), nf90_nowrite, ncid)   
       if (iernc /=0) write(6,*) ncid, 'ncid ', iernc, ' iernc '
!
        iernc=nf90_inq_varid( ncid, 'glon', vid )
        iernc= nf90_get_var( ncid, vid, glon)
        iernc=nf90_inq_varid( ncid, 'glat', vid )
        iernc= nf90_get_var( ncid, vid, glat)
        iernc=nf90_inq_varid( ncid, 'btot', vid )
        iernc= nf90_get_var( ncid, vid, btot)
        iernc=nf90_inq_varid( ncid, 'cormag', vid )
        iernc= nf90_get_var( ncid, vid, cormag)
        iernc=nf90_inq_varid( ncid, 'dipang', vid )
        iernc= nf90_get_var( ncid, vid, dipang)
	
       iernc=nf90_close(ncid) 

    
  end subroutine wam_ion_init
 
  subroutine wam_efield_init(efmodel, me, master)
  
  use wam_efield_setdef_data, only : magnetic_grids, set_readCoef, prep_fk, prep_pnm,  index_quiet
  
  use wam_efieldw05_read_data, only : read_potential, read_schatable, read_bndy
  use wam_efield_setdef_data, only : read_acoef_efield
    implicit none  
    character(len=*)    :: efmodel
    integer, intent(in) :: me, master
    
    call magnetic_grids
    
    call read_acoef_efield (efield_lflux_file, efield_hflux_file)
    
      call index_quiet                ! set up index for f_m(mlt),f_l(UT),f_-k(d)
      call prep_fk	              ! set up the constant factors for f_k
      call prep_pnm	              ! set up the constant factors for P_n^m & dP/d phi  
      
      
      call set_readcoef
!      
! defines      use efield_wam, only: ALAMN =>ALAMN,ALAMX=>ALAMX,ALAMR=>ALAMR,
!     &STPD=>STPD,STP2=>STP2,CSTP=>CSTP,SSTP=>SSTP
      
      if (trim(efmodel) == 'epot') then     
        call read_potential(W05SC_Epot_file)
      else
        call read_potential(W05SC_Bpot_file)
      endif

      call read_schatable(W05SC_HAtbl_file)
      call       read_bndy(W05SC_Bndy_file)

      
   end subroutine  wam_efield_init 
   
  subroutine wamphys_dealloc
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate them"
! before "end" of the FV3GFS
!
    implicit none
!
!   deallocate arrays employed in:
!     wamphys_init
!
   if (allocated (pr_idea)) deallocate (pr_idea)
   if (allocated (prlog))   deallocate (prlog)
   if (allocated (h2ora))   deallocate (h2ora)   
   if (allocated (o3ra))    deallocate (o3ra)   
!   if (allocated (gg))      deallocate (gg) 
   deallocate(gglob, muglob, zglob, tglob )
     
   if (allocated (wam_oh)) deallocate (wam_oh)
   if (allocated (wam_ho2)) deallocate (wam_ho2)       
!
   deallocate( zkm,   pmb, zlog) 
   deallocate (effeuv, effuv)  
   deallocate (eff_hart, eff_hugg, eff_chap, eff_herz)
   deallocate (eff_srb, eff_src, eff_lya)
   deallocate (o2_scale_factor,  srbeff )    
   
   
  end subroutine wamphys_dealloc


 end module wamphys_init_module



! INIT Rad-efficiences
!
!  
     subroutine  heat_uv_eff(levs, eff_hart, eff_src, me, master)
!
      use wamphys_init_module, only : pr_idea
!
      use machine,           only: kind_phys 
      implicit none
      integer,intent(in) :: levs             ! number of pressure levels
      integer, intent(in) :: me, master          
      real(kind=kind_phys) ,  intent(out) :: eff_hart(levs)  ! heat-eff by ML-93 for Hartley
      real(kind=kind_phys) ,  intent(out) :: eff_src(levs)  !

      integer i
      real(kind=kind_phys) c0(3),c1(3),c2(3),c3(3),lp,x
      real(kind=kind_phys) EFF_SRC17(17), EFF_EUV17(17)
!=====================================
!MS-93
      data c0/0.66965,0.92621,    0.75349/
      data c1/-0.009682,0.13396,  0.0036/
      data c2/0.033093,-0.076863, 0.0595/
      data c3/0.017938,0.006897, -0.0228/
!
! pr_idea in mb
!
!     DATA EFF_SRC17/5*.28,.29,.32,.38,.4,.4,.4,.39,.34,.26,.19,.17,.16/         
      do i=1,levs
        lp=log10(pr_idea(i))       
 
         if(lp.ge.0.) then
            eff_hart(i)=c0(2)+c1(2)+c2(2)+c3(2)      ! vay should be ~1 or 0.99
          elseif(lp.ge.-2) then
            x=1.+lp
            eff_hart(i)=c0(2)+c1(2)*x+c2(2)*x*x+c3(2)*x*x*x  
          elseif(lp.ge.-4) then
            x=3.+lp
            eff_hart(i)=c0(1)+c1(1)*x+c2(1)*x*x+c3(1)*x*x*x  
            eff_src(i)=c0(3)+c1(3)*x+c2(3)*x*x+c3(3)*x*x*x  
          else
            eff_hart(i)=c0(1)-c1(1)+c2(1)-c3(1) 
            eff_src(i)=c0(3)-c1(3)+c2(3)-c3(3)  
          endif
       enddo
      return
      end  subroutine heat_uv_eff  
      
      subroutine get_eff_zlog( levs, zlog, effuv, effeuv,  o2_scale_factor, srbeff, me, master)

      
      use machine,           only: kind_phys
      use wamphys_math_interp, only : interpol_wamz,interpol_wamz_down
            
      implicit none
      integer, intent(in) :: levs
      integer, intent(in) :: me, master        
      real(kind=kind_phys) , dimension(levs),  intent(in)  ::  zlog    
      real(kind=kind_phys) , dimension(levs),  intent(out) :: effuv, effeuv,  o2_scale_factor, srbeff
!
! fixed data
!
      integer, parameter        :: nz_euv = 17
      real(kind=kind_phys), dimension(nz_euv)   :: effeuv17,effuv17, p17, z17
! O2 scale factor for 15 pressure level 
      integer, parameter        ::  nz_Jo2sf = 15
      real(kind=kind_phys)                ::  JJ_scale_factor(nz_Jo2sf)
      real(kind=kind_phys)                ::  z15(nz_Jo2sf)  
! TIEGCM SRB-efficiency factors on 63-pressures
      integer, parameter  ::  nz_63 = 63
      real(kind=kind_phys)                ::  pres63(nz_63), SRBEFF63(nz_63)
      real(kind=kind_phys)                ::  z63(nz_63) 
    
       DATA EFFEUV17/8*1.0,.75,.6,.62,.54,.49,.41,.33,.30,.30/   
       DATA EFFUV17/0.59, 0.59, 0.58, 0.57, 0.56, 0.52, 0.48, 0.43, &
                   .4,.4,.4,.39,.34,.26,.19,.17,.16/      
                
       DATA JJ_scale_factor/ 1.0, 4.465536, 4.365480, 3.904985, &        
                         3.367959, 3.202786, 2.378429, 1.636311, &      
                         1.423021, 1.452178, 1.588099, 1.714328, &     
                         1.811639, 1.907779, 1.987971/
			 
       DATA pres63/6.90775528,  6.57442194,     &                        
        6.24108859,  5.90775525,  5.57442191,   &                      
        5.24108856,  4.90775522,  4.57442188,   &                       
        4.24108853,  3.90775519,  3.57442185,   &                       
        3.2410885,   2.90775516,  2.57442182,   &                       
        2.24108847,  1.90775513,  1.57442179,   &                       
        1.24108844,  0.9077551,   0.574421757,  &                       
        0.241088414, -0.0922449296,-0.425578273,&                      
        -0.758911616,-1.09224496,  -1.4255783,  &                      
        -1.75891165, -2.09224499,  -2.42557833, &
        -2.75891168, -3.09224502,  -3.42557836, &
        -3.75891171, -4.09224505,  -4.42557839, &
        -4.75891174, -5.09224508,  -5.42557842, &
        -5.75891177, -6.09224511,  -6.42557845, &
        -6.75891179, -7.09224514,  -7.42557848, &
        -7.75891182, -8.09224517,  -8.42557851, &
        -8.75891185, -9.0922452,   -9.42557854, &
        -9.75891188, -10.0922452,  -10.4255786, &
        -10.7589119, -11.0922453,  -11.4255786, &
        -11.7589119, -12.0922453,  -12.4255786, &
        -12.758912,  -13.0922453,  -13.4255787, &
        -13.758912/
!
! local
!
      integer :: k, kup, kdw
 			 
      do k=1, nz_euv
        p17(k)=5.2285*exp(1.-k)
        z17(k)=-log(p17(k))
      enddo	
      
      do k=1, nz_Jo2sf      
        z15(k)=-log(1.0376*exp(1.-k))
      enddo
      
      do k=1, nz_63
        z63(k) = -pres63(k)
      enddo      
      		       
      CALL  interpol_wamz( nz_euv, z17, effuv17,  levs, zlog, effuv,  kup,kdw)
      CALL  interpol_wamz( nz_euv, z17, effeuv17, levs, zlog, effeuv, kup,kdw) 
      
      CALL  interpol_wamz_down(nz_63, z63, SRBEFF63, levs, zlog, SRBEFF, 1.0) 
      
      CALL  interpol_wamz_down(nz_Jo2sf, z15,JJ_scale_factor , levs, zlog, &
                     o2_scale_factor , 1.0) 

     
      end  subroutine  get_eff_zlog
             
      subroutine wam_readno_snoewx(file, me, master)
      
       use wamphys_init_module, only : iulog, no_nz16,no_ny33,no_neofs
       use wamphys_init_module, only : no_mlat, no_zkm, no_m,no_eof
       use netcdf 
              
       implicit none
       character(len=*), intent(in) :: file
       integer, intent(in) ::  me, master  
!locals
       integer :: nz, ny, neofs
       integer :: istat, ierr, astat
       integer :: dim_id, var_id 
       integer :: ncid, vid, iernc
       character(len=256) :: locfn
       integer, dimension(nf90_max_var_dims) :: dimidT  
!----------------------------------------------------------------------
!	... open the netcdf file
!----------------------------------------------------------------------
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
       if (iernc /=0) write(6,*) ncid, 'ncid ', iernc, ' iernc '     
!----------------------------------------------------------------------
!	... read the snoe dimensions
!----------------------------------------------------------------------
      iernc=nf90_inq_varid( ncid, 'EOF', vid )
       ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
       iernc = nf90_inquire_dimension(ncid, dimidT(3), len=neofs)
       iernc = nf90_inquire_dimension(ncid, dimidT(2), len=nz)
       iernc = nf90_inquire_dimension(ncid, dimidT(1), len=ny)
       if(me == master) then
        if (nz .ne. no_nz16 .or. ny.ne.no_ny33 .or. neofs.ne.no_neofs) then
         write(iulog,*)'snoe_rdeof: failed to read expected neofs=nz=ny'
        endif
       endif
         
!----------------------------------------------------------------------
!	... allocate snoe variables
!----------------------------------------------------------------------
 
       allocate( no_mlat(ny),  no_zkm(nz),stat=astat )  
       allocate( no_m(ny, nz),stat=astat )
       allocate( no_eof(ny, nz, neofs),stat=astat )
       
       if( astat /= 0 ) then
       write(iulog,*) ' alloc_err in read_no_snoe no_eof ' 
       write(iulog,*) 'solar_readno_snoew: failed to allocate eofs; error = ',astat
       end if 

!----------------------------------------------------------------------
!	... read the snoe variables
!----------------------------------------------------------------------
        iernc=nf90_inq_varid( ncid, 'lat', vid )
        iernc= nf90_get_var( ncid, vid, no_mlat)
        iernc=nf90_inq_varid( ncid, 'z', vid )
        iernc= nf90_get_var( ncid, vid, no_zkm)

        iernc = nf90_inq_varid( ncid, 'NO', var_id )
        iernc = nf90_get_var( ncid, var_id, no_m )

        iernc = nf90_inq_varid( ncid, 'EOF', var_id )
        iernc = nf90_get_var( ncid, var_id, no_eof )  !(/1,1,1/), (/ny, nz, neofs/), no_eof )

        iernc=nf90_close(ncid)
	if ( me == master ) then
         print *, ' solar_readno_snoewx - performed '
        endif    
      end subroutine wam_readno_snoewx
      
    subroutine wamphys_idate_calendar(idate, fhour, ddd, fddd) 
    
    use machine, only: kind_phys    		 
    implicit none  
! input     
    integer, intent(in)                 :: idate(4)
    real(kind=kind_phys), intent(in)   :: fhour
!out    
    integer, intent(out)                :: ddd    
    real(kind=kind_phys), intent(out)  :: fddd  
!
!locals
!
      real(kind=kind_phys) :: rinc(5), rjday
      integer              :: jdow, jdoy, jday
      real(4)              :: rinc4(5)
      integer              :: w3kindreal, w3kindint
      
      integer ::  iw3jdn
      integer :: jd1, jddd
      
      integer  idat(8),jdat(8)  
          
       
      idat(1:8)    = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc(1:5)    = 0.
      rinc(2) = fhour
!    
!      call w3kind(w3kindreal,w3kindint)
      w3kindreal=8
      if(w3kindreal==4) then
        rinc4 = rinc
        call w3movdat(rinc4, idat,jdat)
      else
        call w3movdat(rinc,  idat,jdat)
      endif           
!     jdate(8)- date and time (yr, mo, day, [tz], hr, min, sec)
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow, ddd, jday)
      fddd = float(ddd) + jdat(5) / 24.   
end  subroutine wamphys_idate_calendar  
           
    subroutine wamphys_day_of_year(yr, mm, dd, ddd)
!
! computes day of year to get tau_limb forcing written with 1-day precision
!    
    implicit none
    integer, intent(in) :: yr, mm, dd
    integer, intent(out) :: ddd
    
    integer ::  iw3jdn
    integer :: jd1, jddd
    jd1     = iw3jdn(yr,1,1)
    jddd    = iw3jdn(yr,mm,dd)
    ddd = jddd-jd1+1
    
    end subroutine wamphys_day_of_year
 
		       
    subroutine wamphys_presolar(im, solhr,slag,sdec,cdec, dayno,     &                              
                         xlon,xlat, sinlat, coslat,xmu,utsec,sda,    &              
                         maglat,maglon,wbtot,wdipang,essa)
!------------------------------------------------------------------------
! calculate solar zenith angle
!------------------------------------------------------------------------
      use machine     , only : kind_phys
      use physcons, pi => con_pi
      
      use wamphys_math_interp,  only :  interp2_ionfield, spole_ion
      use wamphys_set_data_ion, only :  cormag, btot, dipang, glon, glat

      implicit none
!argument
! input
      integer              :: im
      real(kind=kind_phys) :: sdec,slag,solhr,cdec
      real(kind=kind_phys) :: sinlat(im),coslat(im),xlon(im),xlat(im)
      integer, intent(in)  :: dayno                                    ! day of year     
! output
      real(kind=kind_phys), intent(out) ::     xmu(im)       !cos solar zenith angle
      
! output magnetic and electric parameters 


      real(kind=kind_phys), intent(out) :: maglon(im)  !magnetic longitude (rad)
      real(kind=kind_phys), intent(out) :: maglat(im)  !magnetic latitude (rad)
      real(kind=kind_phys), intent(out) :: wbtot(im)    !mapgnetic field strength
      real(kind=kind_phys), intent(out) :: wdipang(im)  !dip angle (degree)
      real(kind=kind_phys), intent(out) :: essa(im)    !magnetic local time
      real(kind=kind_phys), intent(out) :: sda         ! solar declination angle (rad)
      real(kind=kind_phys), intent(out) :: utsec       !universal time      
!local   
      real(kind=kind_phys) :: wcormag(im),   wcmorg(im)
   
      real(kind=kind_phys) :: cns, ss, cc,  ch, ty
      real(kind=kind_phys) ::  pid2, dtr
      integer :: i
!  compute cosine of solar zenith angle for both hemispheres.
      
      pid2 = .5*pi
      dtr=pi/180.
      
      utsec=solhr*3600.
      cns = pi*(solhr-12.)/12.+slag
      do i=1,im
        ss     = sinlat(i) * sdec
        cc     = coslat(i) * cdec
        ch     = cc * cos(xlon(i)+cns)
        xmu(i) = ch + ss

      enddo
!
! get solar declination angle
      ty = (dayno+15.5)*12./365.
      if ( ty > 12.0 ) ty = ty - 12.0
      sda = atan(0.434*sin(pi/6.0*(ty-3.17)))
!   
!     print*,'www8',sda,asin(sdec)
!
!           get maglat maglon btot from  "getmag" in idea_ion.f
!   

      call interp2_ionfield(im, xlat, xlon, wcormag, wbtot, wdipang, cormag, btot, dipang, glon, glat)
      
      call spole_ion(im,xlat,xlon,utsec, sda, maglon, essa, wcmorg)
!      
!old     call getmag(im,utsec,xlat,xlon,sda,btot,dipang,maglon,maglat,essa)
!
! inside interp_ionfield

      do i=1,im
         maglat(i)= pid2-wcormag(i)*DTR
      enddo	       
      return
      
      end subroutine wamphys_presolar
      

