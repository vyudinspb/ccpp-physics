!>  \file wamphys.F90
!>\defgroup wamphys_run WAM Physics General Algorithm

module wamphys

!!!   wamphys_swdef.F90   wamphys_swdef.meta  wamphys_setdata.F90
!!!   wamphys_init.F90 inside of  wamphys.F90
!
!Model%do_wamphys_diag, Model%do_wamphys  Model%do_wamgfs_rad Model%do_wamipe
!wamphys_init_module.F90
!wamphys_co2.F90
!wamphys_co2hc.F90
!wamphys_diag2d.F90
!wamphys_efield.F90
!wamphys_extra_constants.F90
!wamphys_get_pzgeo.F90
!wamphys_h2o.F90
!wamphys_h2oc.F90
!wamphys_ion.F90
!wamphys_ion_empirmodels.F90
!wamphys_math_interp.F90
!wamphys_merge_ipe2wam.F90
!wamphys_molec_dissipation.F90
!wamphys_multigases.F90
!wamphys_rad_o2_o3.F90
!wamphys_setdata.F90
!wamphys_sheat_jrates.F90
!wamphys_tracer_run.F90
!wamphys_weimer2005.F90
! ccpp cannot handle global-scalar time-dep input
!!!wamphys_swdef.F90
!!!wamphys_swio_data.F90
!
! wamphys_tides.F90: online diagnostics for daily-mean state and diurnal/subdiurnal modes
!
!==========================

    use machine,                 only: kind_phys

   
    implicit none

    private

    public wamphys_init, wamphys_run, wamphys_finalize

    logical :: is_initialized = .False.

contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES WAM physics
! ------------------------------------------------------------------------
!>@brief The subroutine initializes the wam physics
!> \section arg_table_wamphys_init Argument Table
!! \htmlinclude wamphys_init.html
!!
! -----------------------------------------------------------------------
!
    subroutine wamphys_init(do_wamipe, do_wamgfs_rad, do_wamphys_diag,         &
                me, master, nlunit, input_nml_file, logunit,                   &
                fn_nml2, jdat, lonr, latr, levs, ak, bk, dtp,                  &
		ntrac, nto1, nto2, nto3, ntqv,                                 &
                con_pi, con_rerth, con_p0,                                     &
                con_g, con_omega,  con_cp, con_rd, con_rv,con_fvirt,           &  
		errmsg, errflg)
!		
!cd /scratch1/NCEPDEV/swpc/Svetlana.Karol/save/BASE_SVN/WAM_STANDALONE/WAMPHYS_STD
!
    use wamphys_common
    use wamphys_init_module,     only: wamphys_init_all
    use wamphys_multigases,      only: wamphys_set_major_tracers
!----  initialization of unified_ugwp
    implicit none
    logical,              intent (in) :: do_wamphys_diag, do_wamgfs_rad,  do_wamipe
    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    character(len=256),   intent (in) :: input_nml_file(:)
    integer,              intent (in) :: logunit
    character(len=64),    intent (inout) :: fn_nml2    
    integer,              intent (in) :: jdat(8)
    integer,              intent (in) :: lonr
    integer,              intent (in) :: latr   
    integer,              intent (in) :: levs
    real(kind=kind_phys), intent (in) :: ak(:)
    real(kind=kind_phys), intent (in) :: bk(:)
    real(kind=kind_phys), intent (in) :: dtp
    integer,              intent (in) :: ntrac, nto1, nto2, nto3, ntqv 
    
          
    real(kind=kind_phys), intent (in) :: con_pi, con_rerth, con_p0
    real(kind=kind_phys), intent (in) :: con_g, con_omega, con_cp, con_rd, con_rv, con_fvirt
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg					

    integer :: ios
    logical :: exists

    integer :: k
    character(len=64) :: fn_nml23 
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    fn_nml2 = 'input.nml'
!
!
! ============
!
! initialize wamphys_common
!   con_pi, con_rerth, con_p0,  con_g, con_omega,  con_cp, con_rd, con_rv, con_fvirt
!
!==========================
!       call wamphys_common_init
!==========================
    pi    = con_pi
    arad  = con_rerth
    p0s   = con_p0
    grav  = con_g
    omega1= con_omega
    cpd   = con_cp
    rd    = con_rd
    rv    = con_rv
    fv    = con_fvirt
    
    grav2  = grav + grav; rgrav  = 1.0/grav ; rgrav2 = rgrav*rgrav
    rdi    = 1.0 / rd ; rcpd = 1./cpd;  rcpd2 = 0.5/cpd
    gor    = grav/rd
    gr2    = grav*gor
    grcp   = grav*rcpd
    gocp   = grcp
    rcpdl  = cpd*rgrav
    grav2cpd = grav*grcp

    pi2      = 2.*pi ;  pih = .5*pi
    pi_24hr = pi2/24. 
    pi2_365d = pi2/365. 
    Pid12 =PI/12.
    Pid3 =PI/3.
    Pid9 =PI/9.
    rad_to_deg=180.0/pi
    deg_to_rad=pi/180.0

    bnv2min = (pi2/1800.)*(pi2/1800.)
    bnv2max = (pi2/30.)*(pi2/30.)
    dw2min  = 1.0
    velmin  = sqrt(dw2min)
    minvel  = 0.5

    omega2  = 2.*omega1
    omega3  = 3.*omega1

    hpscale = 7000. ; hpskm = hpscale*1.e-3
    rhp     = 1./hpscale
    rhp2 = 0.5*rhp; rh4 = 0.25*rhp
    rhp4 = rhp2 * rhp2
    khp  = rhp* rd/cpd
    mkzmin  = pi2/80.0e3
    mkz2min = mkzmin*mkzmin
    mkzmax  = pi2/500.
    mkz2max = mkzmax*mkzmax
    cdmin   = 2.e-2/mkzmax

    rcpdt  = rcpd/dtp
    
    if (me == master ) then
      print *, ' input wamphys fn_nml2 ', trim(fn_nml2)
      print *, ' input input_nml_file ', input_nml_file      
    endif
    

  
    call wamphys_set_major_tracers(me, master, nlunit, fn_nml2, ntrac, nto1, nto2, nto3, ntqv)
        

 
    fn_nml23 = fn_nml2
    call wamphys_init_all(me, master, nlunit, logunit, jdat, fn_nml23, fn_nml2,   &
              lonr, latr, levs, ak, bk, dtp, errmsg, errflg) 

			       			       
       if (errflg/=0) return

    if (me == master) then
       print *,  ' ccpp: wam_init_cuacires   '

       print *,  ' ccpp wamphys_initwam_init  wam_diag '   , do_wamphys_diag
       print *,  ' ccpp wamphys_initwam_init  wam_ipe_cpl ', do_wamipe
       print *,  ' ccpp wamphys_init  wamgfs_rad '         , do_wamgfs_rad
 

       print *, ' ccpp: wam_init_cuacires   '
       
    endif



    is_initialized = .true.


    end subroutine wamphys_init


! -----------------------------------------------------------------------
! finalize of wamphys   (_finalize)
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the wamphysics
!> \section arg_table_wamphys_finalize Argument Table
!! \htmlinclude wamphys_finalize.html
!!

    subroutine wamphys_finalize(errmsg, errflg)

    implicit none
!
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

!!!    call wamphys_dealloc.... use deallocation as needed

    is_initialized = .false.

    end subroutine wamphys_finalize
!
! 
! -----------------------------------------------------------------------
!   driver is called before the lower atmosphere physics
!     Updates of Temperature-Winds-Tracers are performed inside
!     tendency are only for diagnostics          
! -----------------------------------------------------------------------------------------
!<group name="physics"   <subcycle loop="1">
!  order = >get_prs_fv3 >>GFS_suite_interstitial_1 >> wamphys-wamphys_post .....>dcyc2t3
!   
! -----------------------------------------------------------------------------------------
!>@brief These subroutines and modules execute the CUA-CIRES WAM physics
!> \section arg_table_wamphys_run Argument Table
!! \htmlinclude wamphys_run.html
!!
!> \section gen_wamphys_run CIRES WAM Physics Scheme General Algorithm
!! @{
     subroutine wamphys_run(do_wamipe, do_wamgfs_rad, do_wamphys_diag,         &
          me, master, im,  levs, ntrac, lonr,                                  &
	  dtp, fhzero, kdt, jdat, nto1, nto2, nto3, ntqv,                      & 
	  dx, xlon, xlat, xlon_d,xlat_d, sinlat, coslat, area,                 & 
	  oro, solhr,slag,sdec,cdec, coszen, htrsw, htrlw,                     &
	  f107, f107d, kp, kpa,                                                &
	  nhp, nhpi, shp, shpi, swbt, swang, swvel, swbz, swden,               &	 
          ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk, phii, phil,               &
          dudt, dvdt, dtdt, dudt_iwamph, dvdt_iwamph, dtdt_iwamph,             &
	  do1dt_iwamph,  do2dt_iwamph,  dqdt_iwamph,                           &	 
          gzmt, gmmt, gjhr, gshr, go2dr, errmsg, errflg)
	  
!
!########################################################################
!  Attention New Arrays and Names must be ADDED inside
! new location: and rearrange of types May 2022 from the "core-ccpp" team YYYX
! a) ccpp/data/GFS_typedefs.meta
! b) ccpp/data/GFS_typedefs.F90
! c) ccpp/data/CCPP_typedefs.F90/meta "diag-cs is not tested"
!########################################################################
           
!    use wamphys_def_ipewam_cpl, only :    levipe => lowipe_lev150

    use wamphys_init_module, only : nlev_h2o,nlev_co2, nlevc_h2o 
    
    use wamphys_init_module, only : co2my
    use wamphys_init_module, only : gh2ort,gh2ovb,dg1rt,dg2rt, dg1vb,dg2vb
    use wamphys_init_module, only : gdp,  xx, wvmmrc, coeff   
    use wamphys_init_module, only :         
    use wamphys_init_module, only : spw_drivers, swin_drivers
    
    use wamphys_init_module, only : f107_fix, f107a_fix, kp_fix, kpa_fix

    use wamphys_ion,         only : wam_ion_run
!    use efield_wam, only          :  iday,iyear,iday_m,imo    
    implicit none
    logical,                intent (in) :: do_wamipe    
    logical,                intent (in) :: do_wamgfs_rad  
    logical,                intent (in) :: do_wamphys_diag
    
    integer,                 intent(in) :: me, master, im, levs, ntrac,lonr
    integer,                 intent(in) :: nto1, nto2, nto3, ntqv    
    
    real(kind=kind_phys),    intent(in) :: dtp, fhzero
    integer,                 intent(in) :: kdt, jdat(8)

!         
    real(kind=kind_phys), intent(in), dimension(im) :: dx, xlon, xlat, xlon_d,xlat_d, sinlat, coslat
	 
    real(kind=kind_phys),intent(in), dimension(im) :: area, oro, coszen  
	  
    real(kind=kind_phys), intent(in)  :: solhr, slag, sdec, cdec
  
    real(kind=kind_phys) , intent(in) :: nhp, nhpi, shp, shpi, swbt, swang, swvel, swbz, swden 	 
    real(kind=kind_phys),  intent(inout)  ::  f107, f107d, kp, kpa 
    
    real(kind=kind_phys),  intent(in), dimension(im, levs) ::	htrsw, htrlw  
     
    real(kind=kind_phys),intent(inout), dimension(im, levs)        ::  ugrs,vgrs, tgrs   	
    real(kind=kind_phys),intent(inout), dimension(im, levs, ntrac) :: 	 qgrs    
          
    real(kind=kind_phys) , intent(in)                     ::   prsi(im,levs+1)
    real(kind=kind_phys) , intent(in)                     ::   phii(im, levs+1)        
	
    real(kind=kind_phys), intent(in),   dimension(im, levs) :: 	prslk, phil,  prsl	 

    
	 
    real(kind=kind_phys), intent(inout), dimension(im, levs) ::	 gzmt, gmmt, gjhr, gshr, go2dr	 

!    real, intent(inout), dimension(im,lowst_ipe_level:levs)  ::        gzmt, gmmt, gjhr, gshr, go2dr
    	 
    real(kind=kind_phys), intent(inout),dimension(im, levs)   ::	dudt, dvdt, dtdt
    
    
    real(kind=kind_phys)                      ::        rdtp
    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg
    	
         
! WAM-diagnostics (out)

         real(kind=kind_phys), intent(out),dimension(im, levs) :: dudt_iwamph, dvdt_iwamph,   dtdt_iwamph
         real(kind=kind_phys), intent(out), dimension(im, levs) :: do1dt_iwamph, do2dt_iwamph, dqdt_iwamph
	 
!=========================================================
!    
!  Standalone diagnostics     
!         real(kind=kind_phys), dimension(im, levs) ::   Deddy
!         real(kind=kind_phys), dimension(im, levs) ::	ion_jh, ion_drx, ion_dry
!         real(kind=kind_phys), dimension(im, levs) ::	no_snoe2d,el2d, sped2d, shal2d  
!         real(kind=kind_phys), dimension(im, levs) ::  dt_qsrc, dt_qsrb, dt_qlya,amu2d, jo2_out
	 
!==========================================================
! local variables
!==========================================================
    real(kind=kind_phys)  :: bnhp, bnhpi, bshp, bshpi, bswbt, bswang, bswvel, bswbz, bswden
     	   
    real(kind=kind_phys), dimension(im, levs) ::	wtot
      
    integer, parameter  :: ntrac_i=2                  ! number of 2 WAM chem. tracers (O-O2)
    integer, parameter  :: lowst_ipe_level = 90
    integer, parameter  :: ipr = 10    
    real(kind=kind_phys), parameter  :: sdf = 86400.   
    
    real(kind=kind_phys), dimension(im)    :: xpk_low, xpk_high
 
    
    real(kind=kind_phys), dimension(im)    :: maglat,maglon,btot, dipang, essa
     
    integer, dimension(im) :: plow, phigh    
    integer                :: dayno 
    integer                :: i, k, j1,j2
    
    
    real(kind=kind_phys), dimension(im)       ::  cospass, xmu  
    real(kind=kind_phys), dimension(im, levs) :: nair, rho, am, am29 
    real(kind=kind_phys), dimension(im, levs) ::  o3_n, o3_ng  
    real(kind=kind_phys), dimension(im, levs) ::  o_n,  o2_n, n2_n, dtincw, dtincr
    
!      real(kind=kind_phys), dimension(im, levs) ::	dtno,dteuv  
				
    real(kind=kind_phys), dimension(im, levs) ::   dtRAD,  dtco2c, dtco2h, dth2oh, dth2oc, dto3
    real(kind=kind_phys), dimension(im, levs) ::   cp, grav, zgeo
    real(kind=kind_phys), dimension(im, levs) ::   dudt_ion, dvdt_ion, dtdt_ion  
       
       

       real(kind=kind_phys)  ::  rdelp(im, levs), del(im, levs)     ! rdelp(i,k) = 1./(PRSI(i,k) - PRSI(i,k+1))
                                                  
       real(kind=kind_phys)  ::  kappa(im, levs)                 ! mid-lev values k=R/cp for eddy heat cond
       real(kind=kind_phys)  ::  exner_i(im, levs+1), exner(im,levs) !   (1.e5/P)^(R/Cp)
       real(kind=kind_phys)  ::   amu2d(im, levs),   jo2_out(im,levs)
!        
! geo-solar-related vars
       real(kind=kind_phys)  :: utsec,   sda
       integer               ::  ddd_year       
       integer               ::  nth2o
                 

       
       real(kind=kind_phys)   :: wn1, wsum                
       logical               :: wamNAN

    ! Initialize CCPP error handling variables
       real(kind=kind_phys)  :: heatmax, sdtincw, sdtincr
    errmsg = ''
    errflg = 0
    nth2o  = ntqv          ! standard first index for H2O, ntqv - number of water-based species
    
    heatmax = 0.002          !  100 K/day
!---------------------------------------------
! two ways to setup from namelist/SW-arrays
! ---------------------------------------------   
    f107  = f107_fix
    f107d = f107a_fix
    kp    = kp_fix
    kpa   = kpa_fix	  
!=========================================================================	  
! compute    zg, grav, exner, exner_i, kappa_i cp =>   :prsik, prslk
!=========================================================================
   if (.not.do_wamphys_diag) return 
	dudt = 0.
	dvdt = 0.
        dtdt = 0.


    if (do_wamphys_diag) then
!
       do1dt_iwamph(:,:) =qgrs(:,:, nto1)
       do2dt_iwamph(:,:) =qgrs(:,:, nto2)
       dqdt_iwamph (:,:) =qgrs(:,:, nto3)
       dudt_iwamph(:,:) = ugrs
       dvdt_iwamph(:,:) = vgrs
       dtdt_iwamph(:,:) = tgrs   
    endif   
  
        if ( me == master ) then
	    print *, ' tgrs-W1: ', maxval(tgrs), minval(tgrs) 
        endif

      
      call wamphys_day_of_year(jdat(1), jdat(2), jdat(3), ddd_year) 
      
!      iday = ddd_year
!      iyear = jdat(1)
!      iday_m = jdat(3)
!      imo = jdat(3)
      
      do k=1, levs 
         del(:,k) = prsi(:,k) - prsi(:, k+1)
      enddo 
!      
! wamphys_get_pzgeo.F90: wam_get_rdmulti(im, levs,ntrac,q,xr)     
!   call wam_get_cpr(im,levs,ntrac,qhrs, cp, rdmulti)   >>cv = cp -rdmulti <<
!      
      call wamphys_zgrav(im, levs, ntrac, tgrs, qgrs,          &
                prsl,  prsi, phii, phil, del, oro,             & 
                zgeo, grav, exner, exner_i, kappa, cp, rdelp)   
		
      call wamphys_presolar(im,solhr,slag,sdec,cdec,ddd_year, xlon,xlat,      &                  
           sinlat, coslat, cospass, utsec, sda, maglat,maglon,btot,dipang,essa)

         if ( me == master ) then 
!	   print *, ' bf wamphys_tracer_run ind-ces ', ntrac,ntrac_i,nto1, nto2, nto3
!	   print *,  ' nto1 ', nto1, maxval(qgrs(:,:, nto1)), minval(qgrs(:,:, nto1))
!	   print *,  ' nto2 ', nto2, maxval(qgrs(:,:, nto2)), minval(qgrs(:,:, nto2))
!	   print *,  ' nto3 ', nto3, maxval(qgrs(:,:, nto3)), minval(qgrs(:,:, nto3))
         endif   

!       IF(do_wamipe) call get_vertical_parameters_for_merge      &
!           gzmt, prsl, im, lowst_ipe_level, levs, plow, phigh, xpk_low, xpk_high)
	
	 
! define arrays	 
!  go2dr, plow, phigh, xpk_low, xpk_high	
!    
       go2dr(:,:)=0.
       gzmt(:,:) =0.
       gmmt(:,:) =0. 
       gjhr(:,:) =0. 
       gshr(:,:) =0.
!       
       do i=1, im
         plow(i)  = 110
	 phigh(i) = 110+2
	 xpk_low(i)  =    prsl(i,plow(i))
	 xpk_high(i) =    prsl(i,phigh(i))
       enddo
      	 
       
      call wamphys_tracer_run(me, master,im, levs,ntrac,ntrac_i,nto1, nto2, nto3,   &
                   ntqv, ddd_year, dtp, grav,prsi,prsl, tgrs, qgrs,                 &
                   o_n,o2_n, o3_n, n2_n,nair,rho,am, amu2d, cospass, zgeo, jo2_out, &
		   f107,f107d, go2dr, plow, phigh, xpk_low, xpk_high)
		   
         if ( me == master ) then 
!	   print *, ' Af wamphys_tracer_run ind-ces ', ntrac,ntrac_i,nto1, nto2, nto3
!	   print *,  ' nto1 ', nto1, maxval(qgrs(:,:, nto1)), minval(qgrs(:,:, nto1))
!	   print *,  ' nto2 ', nto2, maxval(qgrs(:,:, nto2)), minval(qgrs(:,:, nto2))
!	   print *,  ' nto3 ', nto3, maxval(qgrs(:,:, nto3)), minval(qgrs(:,:, nto3))
         endif 		   
!      
!      call wamphys_getcp (im, levs, ntrac, qgrs, cp, nto1, nto2, nto3, ntqv)  
!		   	  
       call wam_get_cp(im, levs,  ntrac, qgrs,  cp)     
!=========================

      call wamphys_molec_dissipation(im, levs,grav,prsi,prsl,  &
                   ugrs, vgrs, tgrs,o_n,o2_n,n2_n,dtp,cp, rho)	
     	   
       call wamrad_o3prof(im, levs,ntrac, nto3,  qgrs ,am, nair, o3_ng)
       
!   
       call wamphys_heat_uveuv(im, levs, tgrs, cospass,o_n,o2_n,o3_n,n2_n, rho, cp, &
            ddd_year, prsl, zgeo, grav,am, maglat, f107, f107d, kpa, dtrad)


!diag-standalone  dteuv, dt_qsrc, dt_qsrb, dt_qlya, dtno, no_snoe2d) 

	    
       call wamrad_co2(im,  levs,nlev_co2,ntrac,nto1, nto2, nto3, ntqv, &
                       co2my, grav,cp, qgrs,tgrs, dtco2c,cospass,dtco2h)
		       
       	
       call wamrad_h2o(im,  levs,nlev_h2o,nlevc_h2o,ntrac,nto1, nto2, nto3, ntqv, &
                       grav,cp, qgrs, tgrs, dth2oh, cospass, dth2oc,              & 
		gh2ort,gh2ovb,dg1rt,dg2rt, dg1vb,dg2vb,gdp,xx,wvmmrc,coeff)

       
       call wamrad_o2_o3(im,levs,cospass, tgrs ,o2_n, o3_ng, rho, cp, zgeo, grav, dto3) 

            
      if ( me == master .and. kdt == 1) print *,  ' after wamrad_o2_o3 '
    
 
       	    
       IF(do_wamipe) then 
          do k = lowst_ipe_level, levs
          do i = 1, im
            gjhr(i, k) = gjhr(i, k) / cp(i, k)
            gshr(i, k) = gshr(i, k) / cp(i, k)	  
	  enddo
	  enddo
!	   call idea_merge_ipe_to_wam(GSHR, dtrad, im,, levs, lowst_ipe_level, prsl, &
!	    plow, phigh, xpk_low, xpk_high)
        endif	
	    		       
      if(do_wamgfs_rad) then 
 
      do i=1,im
!coszen (im)  - real, avg of cosz over daytime sw time-call-interval
! daytime, XCOSZ=cosine of solar zenith angle at current time
! if ( xcosz(i) > f_eps .and. coszen(i) > f_eps )xmu(i) = xcosz(i) / coszen(i) ~1.
! for day take instant values radiation every timestep
!
        if(cospass(i) > 1.e-4 .and. coszen(i) > 1.e-4) then
             xmu(i) = cospass(i)/coszen(i)
        else
!night
          xmu(i) = 0.
        endif
      enddo

	 call wamphys_rad_merge(me, master, im ,levs, xmu, prsl, htrlw, htrsw, wtot, &
                             dtrad, dtco2c,dtco2h,dth2oh,dth2oc,dto3) 
			     		
         tgrs = tgrs + wtot * dtp  


       endif
!
! fixed solar/geo parameters
!      
    bnhp=27.42402
    bnhpi= 7   
    bshp =33.016
    bshpi =7
    bswbt =6.670438
    bswang =213.0435
    bswvel =450.7061
    bswbz =-3.442665
    bswden =6.949437      
! 
   
	    
      call wam_ion_run(im, levs, jdat, prsl, solhr,cospass,zgeo, grav, o_n,o2_n,n2_n, &
            cp, ugrs, vgrs, tgrs,dudt_ion,dvdt_ion,dtdt_ion,rho,xlat,xlon,            &
            ddd_year,utsec,sda,maglon,maglat,btot,dipang,essa, f107, f107d, kp,       &
            bnhp, bnhpi, bshp, bshpi, bswbt, bswang, bswvel, bswbz, bswden) 
!	    swin_drivers, spw_drivers)

  
!       IF(do_wamipe) then 
!          call idea_merge_ipe_to_wam(GZMT, dudt_ion,&
!              im, levs, lowst_ipe_level, prsl, plow, phigh, xpk_low, xpk_high)
!          call idea_merge_ipe_to_wam(GMMT, dvdt_ion,&
!              im, levs, lowst_ipe_level, prsl, plow, phigh, xpk_low, xpk_high)
!	   call idea_merge_ipe_to_wam(GJHR, dtdt_ion,&
!              im, levs, lowst_ipe_level, prsl, plow, phigh, xpk_low, xpk_high)
!       endif	

        ugrs = ugrs +dtp*dudt_ion
	vgrs = vgrs +dtp*dvdt_ion
	tgrs = tgrs +dtp*dtdt_ion
	
!	call wamphys_dadj(prsi, tgrs, im, levs, me, master)
!	call wamphys_dadj_or(prsi, tgrs, im, levs+1, me, master)	
 345    continue 
 
 			      
        if ( me == master ) then
	    print *, ' tgrs-W2: ', maxval(tgrs), minval(tgrs) 
!	    print *, ' ugrs: ', maxval(ugrs), minval(ugrs) 
!	    print *, ' vgrs: ', maxval(vgrs), minval(vgrs) 	    	    
	endif	
	
	if (do_wamphys_diag) then
       rdtp = 1./dtp
       do1dt_iwamph(:,:) = (qgrs(:,:, nto1)  -do1dt_iwamph(:,:))*rdtp
       do2dt_iwamph(:,:) = (qgrs(:,:, nto2) -do2dt_iwamph(:,:))*rdtp
       dqdt_iwamph (:,:) = (qgrs(:,:, nto3) -dqdt_iwamph(:,:))*rdtp
       dudt_iwamph(:,:) = (ugrs(:,:) - dudt_iwamph(:,:))*rdtp
       dvdt_iwamph(:,:) = (vgrs(:,:) - dvdt_iwamph(:,:))*rdtp
       dtdt_iwamph(:,:) = (tgrs(:,:) - dtdt_iwamph(:,:))*rdtp 	
	endif 
	return
	
!!! tests: rt.sh -l wamphys_rt.conf -wk 

    end subroutine wamphys_run
!! @}
end module wamphys
