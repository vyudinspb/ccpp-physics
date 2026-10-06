!>\file cua_ugwpv1_ngw_ufswam.F90
!!

module cua_ugwpv1_ngw_ufswam

      use machine,          only : kind_phys
      use cires_ugwpv1_module,only : amo, amo2, amn2, amhe, invbz
      use cires_ugwpv1_module,only : Runiv, avgdbz
      use cires_ugwpv1_module,only : muo, muo2, mun2, muhe
      use cires_ugwpv1_module,only : lao, lao2, lan2, lahe       
      use cires_ugwpv1_module,only : cpo, cpo2, cpn2, cphe  
            
      use ugwp_common ,     only : rgrav,  grav,  cpd,    rd,  rv, rcpdl, grav2cpd,    &
                                   omega2,  rcpd,   rcpd2,  pi,    pi2, fv,            &
                                   rad_to_deg, deg_to_rad,                             &
                                   rdi,        gor,    grcp,   gocp,                   &
                                   bnv2min,  bnv2max,  dw2min, velmin, gr2,            &
                                   hpscale, rhp, rh4, grav2, rgrav2, mkzmin, mkz2min
				   
      use ugwp_common ,     only : arad		
contains


!----------------------------------------------------------------
!  Broad spectrum FVS-1993, mkz^nSlope with nSlope = 0, 1,2
!  dissipative solver with NonHyd/ROTational-effects;
!  "Reflected" GWs are treated as waves with "zero-upward" flux,
!  they are "out" of given vertical column
!  No "lateral" prop-n of GWs
!
!  Jan-May 2026, V. Yudin CUA new solver for NGWs adapted for UFSWAM
!                moderate horizontal resolutions C96/C192
!
!   cires_ugwpv1_ngw_solv_ufswam 
!----------------------------------------------------------------

      subroutine ugwpv1_ngw_solv_ufswam(mpi_id, master, im, levs, ntrac, kdt, dtp, &
                 tau_ngw, tm , um, vm, q3m, qm, prsl, prsi, zmet,  zmeti, delp,  &
                 xlatd, sinlat, coslat,                                     &
                 pdudt, pdvdt, pdtdt, dked, zngw)
!
!--------------------------------------------------------------------------------
!      nov 2015 alternative gw-solver for nggps-wam
!      nov 2017 nh/rotational gw-modes for nh-fv3gfs
!      oct 2019 adding empirical satellite-based
!          source function and *F90 CIRES-style of the code
!      oct 2020  Diagnostics of "tauabs, wrms, trms" is taken out
!      May 2025 V. Yudin: New 3 knobs for UFSWAM Space Weather Model
! --------------------------------------------------------------------------------
! 
      use cires_ugwpv1_module,only : knob_ugwp_ufswam, knob_ugwp_wam_med, knob_ugwp_msp_ind 
      
      use cires_ugwpv1_module,only :  krad, kvg, kion, ktg, iPr_ktgw, Pr_kdis, Pr_kvkt

      use cires_ugwpv1_module,only : knob_ugwp_doheat, knob_ugwp_dokdis, idebug_gwrms
      
      use cires_ugwpv1_module,only : knob_ugwp_ufswam, knob_ugwp_wam_med
      
      use cires_ugwpv1_module,only : msp_ind => knob_ugwp_msp_ind, msp_dim  => knob_ugwp_msp_dim     
      
      use cires_ugwpv1_module,only :  psrc => knob_ugwp_palaunch, zsp_gw => knob_ugwp_sponge

      use cires_ugwpv1_module,only : maxdudt, maxdtdt, max_eps, dked_min, dked_max


      use ugwp_wmsdis_init, only : v_kxw,  rv_kxw,   v_kxw2, tamp_mpa, tau_min, ucrit, &
                                   gw_eff,                                             &
                                   nslope,  ilaunch, zms,                              &
                                   zci,     zdci,    zci4, zci3, zci2,                 &
                                   zaz_fct, zcosang, zsinang,  nwav,    nazd,          &
                                   zcimin, zcimax, rimin, sc2, sc2u, ric, lturb
!
      implicit none
!
!      real(kind=kind_phys), parameter   :: zsp_gw  = 106.5e3  ! sponge for GWs above the model top namelist
!
      real(kind=kind_phys), parameter   :: linsat2 = 1.0
      real(kind=kind_phys), parameter   :: dturb_max = 150.0   !!!250.0  ! max val  of eddy-coeff for kdis=2
      real(kind=kind_phys), parameter   :: fg_cool  = 0., fg_ktur = 0.      
      integer,              parameter   :: ener_norm =0
      integer,              parameter   :: ener_lsat=0
      
      integer,              parameter   :: nstdif = 1          ! 3-point smmoother application 1-onetime  
      logical,              parameter   :: hb93_kt=.false.     ! HB-93 CAM free atm-re eddy diff-n 
      logical,              parameter   :: debug_conv = .false.          

      integer, intent(in)  :: levs                            ! vertical level
      integer, intent(in)  :: im                              ! horiz tiles
      integer, intent(in)  :: ntrac                           ! number of tracers that include O-O2-N2      
      integer, intent(in)  :: mpi_id, master, kdt

      real(kind=kind_phys) ,intent(in)   :: dtp               ! model time step
      real(kind=kind_phys) ,intent(in)   :: tau_ngw(im)

      real(kind=kind_phys) ,intent(in)   :: vm(im,levs)       ! meridional wind
      real(kind=kind_phys) ,intent(in)   :: um(im,levs)       ! zonal wind
      real(kind=kind_phys) ,intent(in)   :: qm(im,levs)       ! spec. humidity
      real(kind=kind_phys) ,intent(in)  :: q3m(im,levs,ntrac) ! all tracers    
!
! Upgrades for UFSWAM-2025 C96/C192 need all species to compute molec. dissipation     
!           
      
      real(kind=kind_phys) ,intent(in)   :: tm(im,levs)       ! kinetic temperature

      real(kind=kind_phys) ,intent(in)   :: prsl(im,levs)     ! mid-layer pressure
      real(kind=kind_phys) ,intent(in)   :: delp(im,levs)     ! mid-layer thickness d(prsi)
      real(kind=kind_phys) ,intent(in)   :: zmet(im,levs)     ! meters now !!!!!       phil =phil=grav*Z
      real(kind=kind_phys) ,intent(in)   :: prsi(im,levs+1)   !  interface pressure
      real(kind=kind_phys) ,intent(in)   :: zmeti(im,levs+1)  !  interface geopi/meters
      real(kind=kind_phys) ,intent(in)   :: xlatd(im)         ! xlat_d in degrees
      real(kind=kind_phys) ,intent(in)   :: sinlat(im)
      real(kind=kind_phys) ,intent(in)   :: coslat(im)
!
! out-gw effects on u-v-t and for UFSWAM 
!                   "consider" next diffusion of major species O-O2-N2-He with coupled equations Mdif(i,j)+ Deddy
!
!
      real(kind=kind_phys) ,intent(out) :: pdudt(im,levs)     ! zonal momentum tendency
      real(kind=kind_phys) ,intent(out) :: pdvdt(im,levs)     ! meridional momentum tendency
      real(kind=kind_phys) ,intent(out) :: pdtdt(im,levs)     ! gw-heating (u*ax+v*ay)/cp and cooling
      real(kind=kind_phys) ,intent(out) :: dked(im,levs)      ! gw-eddy diffusion
      real(kind=kind_phys) ,intent(out) :: zngw(im)           ! launch height
      
!     real(kind=kind_phys) ,intent(out) :: pdqdt(im,levs, ntrac)  only [msp_ind:msp_ind+2]
!
!
! local ===========================================================================================

      real(kind=kind_phys)              :: tauabs(im,levs)    !
      real(kind=kind_phys)              :: wrms(im,levs)      !
      real(kind=kind_phys)              :: trms(im,levs)      !
                                          
      real(kind=kind_phys)              :: zwrms(nwav,nazd), wrk1(levs), wrk2(levs)
      real(kind=kind_phys)              :: atrms(nazd, levs),awrms(nazd, levs), akzw(nwav,nazd, levs+1)
!
! local ===========================================================================================

      real(kind=kind_phys)              :: no2(levs), no1(levs), nn2(levs), nhe(levs), nair(levs)  
      real(kind=kind_phys)              :: vumol(levs+1), ktmol(levs+1), grav_vari(levs+1)  
      real(kind=kind_phys)              :: cp_msp(levs), mu_msp(levs),   grav_var(levs)
                  
      real(kind=kind_phys)              :: taux(levs+1)         ! EW component of vertical momentum flux (pa)
      real(kind=kind_phys)              :: tauy(levs+1)         ! NS component of vertical momentum flux (pa)
      real(kind=kind_phys)              :: fpu(nazd, levs+1)    ! az-momentum flux
      real(kind=kind_phys)              :: ui(nazd, levs+1)     ! azimuthal wind

      real(kind=kind_phys)              :: fden_bn(levs+1)         ! density/brent
      real(kind=kind_phys)              :: flux  (nwav, nazd) ,   flux_m  (nwav, nazd)                
!   
       real(kind=kind_phys)  :: bn(levs+1)                    ! interface BV-frequency
       real(kind=kind_phys)  :: bn2(levs+1)                   ! interface BV*BV-frequency
       real(kind=kind_phys)  :: rhoint(levs+1)                ! interface density
       real(kind=kind_phys)  :: uint(levs+1)                  ! interface zonal wind
       real(kind=kind_phys)  :: vint(levs+1)                  ! meridional wind
       real(kind=kind_phys)  :: tint(levs+1)                  ! temp-re
       real(kind=kind_phys)  :: kappat(levs)                  ! R/cp
       real(kind=kind_phys)  :: irhodz_mid(levs)
       real(kind=kind_phys)  :: suprf(levs+1)                 ! RF-super linear dissipation
       real(kind=kind_phys)  :: cstar(levs+1) ,cstar2(levs+1)
       real(kind=kind_phys)  :: v_zmet(levs+1)
       real(kind=kind_phys)  :: vueff(levs+1)
       real(kind=kind_phys)  :: dfdz_v(nazd, levs), dfdz_heat(nazd, levs)    ! axj = -df*rho/dz  directional Ax

       real(kind=kind_phys), dimension(levs)   ::  atm , aum, avm, aqm, aprsl, azmet, dz_met
       real(kind=kind_phys), dimension(levs+1) ::  aprsi, azmeti, dz_meti, arcpd

       real(kind=kind_phys), dimension(levs)   :: wrk3
       real(kind=kind_phys), dimension(levs)   :: uold, vold, told, unew, vnew, tnew, ptold, pkold
       real(kind=kind_phys), dimension(levs)   :: rho, rhomid, adif, cdif, acdif
       real(kind=kind_phys), dimension(levs)   :: Qmid, AKT
       real(kind=kind_phys), dimension(levs+1) :: dktur, Ktint, Kvint
       real(kind=kind_phys), dimension(levs+1) :: fden_lsat, fden_bnen

       integer,              dimension(levs)   :: Anstab

       real(kind=kind_phys)  :: sig_u2az(nazd),  sig_u2az_m(nazd)
       real(kind=kind_phys)  :: wave_dis(nwav, nazd), wave_disaz(nazd)
       real(kind=kind_phys)  :: rdci(nwav),  rci(nwav)
       real(kind=kind_phys)  :: wave_act(nwav, nazd)           ! active waves at given vert-level
       real(kind=kind_phys)  :: ul(nazd)                       ! velocity in azimuthal direction at launch level
!
! wamphys
!       
       real(kind=kind_phys)  :: rdmult                         ! rdi => 287. onto Rd/mu_multi       
!
! scalars
!
       real(kind=kind_phys)  :: bvi, bvi2, bvi3, bvi4, rcms       ! BV at launch level
       real(kind=kind_phys)  :: c2f2, cf1, wave_distot


       real(kind=kind_phys)  ::  flux_norm                         ! norm-factor
       real(kind=kind_phys)  ::  taub_src, rho_src, zcool, vmdiff
!
       real(kind=kind_phys)  :: zthm, dtau, cgz, ucrit_maxdc
       real(kind=kind_phys)  :: vm_zflx_mode, vc_zflx_mode
       real(kind=kind_phys)  :: kzw2, kzw3, kdsat, cdf2, cdf1, wdop2,v_cdp2
       real(kind=kind_phys)  :: ucrit_max
       real(kind=kind_phys)  :: pwrms, ptrms
       real(kind=kind_phys)  :: zu, zcin, zcin2, zcin3, zcin4, zcinc
       real(kind=kind_phys)  :: zatmp, fluxs, zdep,  ze1, ze2
       real(kind=kind_phys)  :: rdi_var, rdi_vari, cpint
!
       real(kind=kind_phys)  :: zdelp, zdelm, taud_min
       real(kind=kind_phys)  :: tvc,  tvm, ptc, ptm
       real(kind=kind_phys)  :: umfp, umfm, umfc, ucrit3
       real(kind=kind_phys)  :: fmode, expdis, fdis
       real(kind=kind_phys)  :: v_kzi, v_kzw, v_cdp, v_wdp, tx1, fcorsat, dzcrit
       real(kind=kind_phys)  :: v_wdi, v_wdpc
       real(kind=kind_phys)  :: ugw, vgw, ek1, ek2, rdtp, rdtp2, rhp_wam

       integer :: j, jj, k, kk, inc, jk, jkp, jl, iaz
       integer ::  ksrc, km2, km1, kp1, ktop
!
! Kturb-part
!
      real(kind=kind_phys)     :: uz, vz, shr2 , ritur, ktur, zlturb, hpwam

      real(kind=kind_phys)     :: kamp, zmetk, zgrow
      real(kind=kind_phys)     :: stab, stab_dt, dtstab
      real(kind=kind_phys)     :: nslope3
!
      integer                  :: nstab, ist
      real(kind=kind_phys)     :: w1, w2, w3, dtdif

      real(kind=kind_phys)     ::  dzmetm, dzmetp, dzmetf, bdif, bt_dif, apc, kturp
      real(kind=kind_phys)     ::  rstar, rstar2

      real(kind=kind_phys) :: snorm_ener, sigu2, flux_2_sig, ekin_norm
      real(kind=kind_phys) :: taub_ch, sigu2_ch
      real(kind=kind_phys) :: Pr_kdis_eff,  mf_diss_heat, iPr_max
      real(kind=kind_phys) :: exp_sponge, mi_sponge, gipr, pkappa

!--------------------------------------------------------------------------
!
        nslope3  = nslope + 3.0
        Pr_kdis_eff = gw_eff*pr_kdis
	pkappa =rd/cpd
        iPr_max = max(1.0,  iPr_ktgw)
        gipr  = grav* Ipr_ktgw
!
! test for input fields
!        if (mpi_id == master .and. kdt < 2) then
!              print *, im, levs, dtp, kdt,    ' vay-solv2-v1'
!              print *,  minval(tm), maxval(tm), ' min-max-tm '
!          print *,  minval(vm), maxval(vm), ' min-max-vm '
!          print *,  minval(um), maxval(um), ' min-max-um '
!          print *,  minval(qm), maxval(qm), ' min-max-qm '
!          print *,  minval(prsl), maxval(prsl), ' min-max-Pmid '
!          print *,  minval(prsi), maxval(prsi), ' min-max-Pint '
!          print *,  minval(zmet), maxval(zmet), ' min-max-Zmid '
!          print *,  minval(zmeti), maxval(zmeti), ' min-max-Zint ' 
!          print *,  minval(prslk), maxval(prslk), ' min-max-Exner '
!          print *,  minval(tau_ngw), maxval(tau_ngw), ' min-max-taungw '
!    endif
!

    if (idebug_gwrms == 1) then
      tauabs=0.0; wrms =0.0 ; trms =0.0
    endif

       rci(:) = 1./zci(:)
       rdci(:) = 1./zdci(:)

       rdtp = 1./dtp
       rdtp2 = 0.5*rdtp

       ksrc= max(ilaunch, 3)
       km2 = ksrc - 2
       km1 = ksrc - 1
       kp1 = ksrc + 1
       ktop= levs+1

        suprf(ktop) = kion(levs)
#ifdef MULTI_GASES
!        print *, ' recompute zmet/zmeti, variable grav(z), kappat(z) '
!       kappat(1:levs) = rd/cpd
!      cal cires_dry_adjust(levs, tm, kappa, prsl, prsi, delp)	
#endif
       do k=1,levs
          suprf(k) = kion(k)               ! approximate 1-st order damping with Fast super-RF of FV3
            pdvdt(:,k) = 0.0
            pdudt(:,k) = 0.0
            pdtdt(:,k) = 0.0
            dked(: ,k) = 0.0
        enddo
	
       if (knob_ugwp_wam_med == 1) return ! the Mesosccale configuration of UFSWAM-MED (C384-C768) does not use NGW
	                                  ! it returns "zero" tendencies and dked
!-----------------------------------------------------------
! column-based j=1,im pjysics with 1D-arrays
!
! for FV3WAM need to add computation of molecular visc/cond-ty
!
!-----------------------------------------------------------
       DO j=1, im
          jl =j
        tx1           = omega2 * sinlat(j) *rv_kxw
        cf1 = abs(tx1)
        c2f2      = tx1 * tx1
     ucrit_max = max(ucrit, cf1)
     ucrit3 = ucrit_max*ucrit_max*ucrit_max
!
! ngw-fluxes at all gridpoints (with tau_min at least)
!          
     aprsl(1:levs) = prsl(jl,1:levs)
!
! ksrc-define "aprsi(1:levs+1)   redefine "ilaunch"
!
        do k=1, levs
             if (aprsl(k) .lt. psrc ) exit
           enddo
           ilaunch = max(k-1, 3)
           ksrc= max(ilaunch, 3)

           zngw(j) = zmet(j, ksrc)/grav

        km2 = ksrc - 2
        km1 = ksrc - 1
        kp1 = ksrc + 1

!=====ksrc

     aum(1:levs)      = um(jl,1:levs)
     avm(1:levs)      = vm(jl,1:levs)
     atm(1:levs)      = tm(jl,1:levs)
     aqm(1:levs)      = qm(jl,1:levs)
     
     
     aprsi(1:levs+1)  = prsi(jl,1:levs+1)
      
     azmet(1:levs)    = zmet(jl,1:levs)     
     azmeti(1:levs+1) = zmeti(jl,1:levs+1)
     
     nair = aprsl/atm*invbz
     rho_src = aprsl(ksrc)*rdi/atm(ksrc)
     azmeti(levs+1) = azmeti(levs+1)/(grav -azmeti(levs+1)/arad)
     do k=1, levs
       azmet(k) = azmet(k)/(grav -azmet(k)/arad)
       azmeti(k) = azmeti(k)/(grav -azmeti(k)/arad)
       no1(k) = q3m(j, k, msp_ind)
       no2(k) = q3m(j, k, msp_ind+1)
       nn2(k) = 1. - no1(k)-no2(k)
     
     if  (msp_dim == 3) then
       nhe(k) = q3m(j, k, msp_ind+2)
       nn2(k) = 1. -nhe(k)
       else
       nhe(1:levs) = 0. 
      endif 
     enddo 
         
     call get_ufswam_med(levs, msp_dim, no1, no2, nn2, nhe, nair, &
	  atm, aprsl, azmet, azmeti, cp_msp, ktmol, vumol, mu_msp, grav_var, grav_vari) 
!
! recompute Zmet, Zmeti for variable gravity and Rdgas = Runiv/mu_msp (287 => 500 oxygen-helium)
!	  
	 
         taub_ch = max(tau_ngw(jl), tau_min)
         taub_src = taub_ch            


       sigu2  = taub_src/rho_src/v_kxw * zms
       sig_u2az(1:nazd) = sigu2
!
! compute diffusion-based arrays km2:levs
!
       do jk = km2, levs
          dz_meti(jk) = azmeti(jk+1)-azmeti(jk)
          dz_met(jk)  =  azmet(jk)-azmeti(jk-1)
       enddo
!       ---------------------------------------------
!   interface mean flow parameters launch -> levs+1
!       ---------------------------------------------
       do jk= km1,levs
           tvc = atm(jk)          !* (1. +fv*aqm(jk))
           tvm = atm(jk-1)        !* (1. +fv*aqm(jk-1))
	   
           zthm          = 2.0/(tvc+tvm)
	   
!interface
           uint(jk)   = 0.5*(aum(jk-1)+aum(jk))
           vint(jk)   = 0.5*(avm(jk-1)+avm(jk))
           tint(jk)   = 0.5*(tvc+tvm)
	   rdi_var    = Runiv/mu_msp(jk)
	   rdi_vari   = 2*Runiv/(mu_msp(jk)+mu_msp(jk-1))
	   rhp_wam = 0.25*zthm/rdi_vari*grav_vari(jk)
           rhomid(jk) = aprsl(jk)*rdi_var/atm(jk)
           rhoint(jk) = aprsi(jk)*rdi_vari*zthm             !  rho = p/(RTv)
	   
           zdelp      = dz_meti(jk)                         !  >0 ...... dz-meters
           v_zmet(jk)  = 2.*zdelp                           ! 2*kzi*[Z_int(k+1)-Z_int(k)]
           zdelm          = 1./dz_met(jk)                   ! 1/dz  ...... 1/meters
!
!          bvf2 = grav2*zdelm*(ptc-ptm)/(ptc + ptm) ! N2=[g/PT]*(dPT/dz)

           wrk3(jk)= 1./zdelp/rhomid(jk)                     ! 1/rho_mid(k)/[Z_int(k+1)-Z_int(k)]
           irhodz_mid(jk) = rdtp*zdelp*rhomid(jk)/rho_src
! 
!cires_ugwpv1_initialize.F90:      parameter :: ulturb=150.,         sc2u = ulturb* ulturb
!cires_ugwpv1_initialize.F90:      parameter :: lturb = 30. (m),     sc2  = lturb*lturb   
! cires_ugwpv1_sporo.F90:          parameter :: lsc2 = lturb*lturb,  usc2 = uturb*uturb
!
!	
!          rcpdl = cpd/grav  grcp   = grav*rcpd grav2cpd = grav*grcp =grav*grav/Cp
!WAM corrections:    grav(z) and cp_multi(z) interface gradients and averaged values
!cp_msp, ktmol, vumol, mu_msp, grav_var, grav_vari

! UFS          bn2(jk)    = grav2cpd*zthm*(1.0+rcpdl*(tvc-tvm)*zdelm)

! UFSWAM
           cpint = .5*(cp_msp(jk)+cp_msp(jk-1))
	   bn2(jk)    = grav_vari(jk)*zthm*(grav_vari(jk)/cpint+ (tvc-tvm)*zdelm)
	   
           uz = aum(jk) - aum(jk-1)
           vz = avm(jk) - avm(jk-1)
           shr2 = (max(uz*uz+vz*vz, dw2min)) * zdelm *zdelm
	   
	   if (bn2(jk) < 0. ) then 
	     if (debug_conv) then	   
	       print *, 'GWPV1-INP CINSTAB', tvc, tvm, nint(azmet(jk)*1.e-3)
	       print *, 'GWPV1-INP CIN_RIC', bn2(jk)/shr2, ' lat ', nint(xlatd(j))
	     endif
	   endif
	   
           bn2(jk)    = max(min(bn2(jk), bnv2max), bnv2min)
           bn(jk)     = sqrt(bn2(jk))
	   ritur = bn2(jk)/shr2	   
           zmetk  =  azmet(jk)* rhp_wam                     ! rh4-ufs mid-layer height k_int => k_int+1
           zgrow = exp(zmetk)	   
	   zlturb = min(lturb*zmetk, zdelp*.5)
           kamp = sqrt(shr2)* zlturb*zlturb
           w1 = 1./(1. + 5.*ritur)
	   w1 = w1*w1
	   if (hb93_kt) then
	     w1 = 1./(1. + 10.*ritur)/(1.+8.*ritur)
	   endif
!	   
!HB-1993 Compute_eddy_diff	   
!ritur < 0: fofri = sqrt(max(1. - 18.*ritur,0.)) else 	
!  fofri = 1.0_r8/(1.0_r8 + 10.0_r8*ri(i,k)*(1.0_r8 + 8.0_r8*ri(i,k)))
!   gipr  = grav* Ipr_ktgw 
!	   
             ktur= min(max(kamp * w1, dked_min), dked_max)
             zmetk =  azmet(jk)* rhp
	     
             vueff(jk)  = vumol(jk) + ktur*fg_ktur + kvg(jk)       !kvg-MLT see "cires_ugwpv1_initialize.F90"

             akt(jk) = gipr/tvc
          enddo

        if (idebug_gwrms == 1) then
         do jk= km1,levs
            wrk1(jk) = rv_kxw/rhoint(jk)
            wrk2(jk)=  rgrav2*zthm*zthm*bn2(jk)    ! dimension [K*K]*(c2/m2)
          enddo
        endif

!
! extrapolating values for ktop = levs+1 (lev-interface for prsi(levs+1) =/= 0)
!
           jk = levs

	   rdi_vari   = 2*Runiv/(mu_msp(k)+mu_msp(k-1))
           rhoint(ktop) = 0.5*aprsi(levs)*rdi_vari/atm(jk)
           tint(ktop)  = atm(jk)       !*(1. +fv*aqm(jk))
	   
           uint(ktop)  = aum(jk)
           vint(ktop)  = avm(jk)

           v_zmet(ktop) =  v_zmet(jk)
           vueff(ktop)  = vueff(jk)
           bn2(ktop)    = bn2(jk)
           bn(ktop)    = bn(jk)
!
! akt_mid *KT = -g*(1/H + 1/T*dT/dz)*KT     ... grav/tvc     for eddy heat conductivity
!  GW-eddy cooling:  [m2/s3] 
!  gor =grav / rd=rdi_vari

     do jk=km1, levs
           rdi_vari   = .5*grav_vari(jk)*(mu_msp(jk)+mu_msp(jk-1))/Runiv
           akt(jk) = -akt(jk)*(rdi_vari + (tint(jk+1)-tint(jk))/dz_meti(jk) )
     enddo


           bvi = bn(ksrc);  bvi2 = bvi * bvi;
       bvi3 = bvi2*bvi; bvi4 = bvi2 * bvi2; rcms = zms/bvi
!
! project winds at ksrc
!               
        do iaz=1, nazd
           ul(iaz) = zcosang(iaz) *uint(ksrc) + zsinang(iaz) *vint(ksrc)
        enddo
!

          do jk=ksrc, ktop    
             cstar(jk) = bn(jk)/zms
             cstar2(jk) = cstar(jk)*cstar(jk)

             fden_lsat(jk) = rhoint(jk)/bn(jk)*v_kxw*Linsat2
       
           do iaz=1, nazd
               zu = zcosang(iaz)*uint(jk) + zsinang(iaz)*vint(jk)
               ui(iaz, jk) =  zu      !!!- ul(iaz)*0. shift of winds
             enddo
          enddo

          rstar = 1./cstar(ksrc)
          rstar2 = rstar*rstar
!      -----------------------------------------
!       set launch momentum flux spectral density
!       -----------------------------------------

          fpu(1:nazd, km2:ktop) =0.

        do inc=1,nwav

           zcin  = zci(inc)*rstar

!
! integrate (flux(cin) x dcin  ) old tau-flux and normalization
!
           flux(inc,1) =  rstar*(zcin*zcin)/(1.+ zcin**nslope3)
!
! fsat  = rstar*(zcin*zcin) * taub_src / SN * [rho/rho_src *N_src/N]
!
            fpu(1,ksrc) = fpu(1,ksrc) + flux(inc,1)*zdci(inc)    ! dc/cstar = dim-less
             
           do iaz=1,nazd   
             akzw(inc, iaz, ksrc) = bvi*rci(inc)
           enddo

         enddo
!
! adjust rho/bn vertical factors for saturated fluxes (E(m) ~m^-3)

         flux_norm  = taub_src / fpu(1, ksrc)             ! [Pa * dc/cstar *dim_less]
         ze1 =  flux_norm * bvi/rhoint(ksrc) *rstar *rstar2
      do jk=ksrc, ktop
          fden_bn(jk) = ze1* rhoint(jk) / bn(jk)          ! [Pa]/[m/s] * rstar2
      enddo
!
      do inc=1, nwav
          flux(inc,1) = flux_norm*flux(inc,1)
      enddo


  if (ener_norm == 1) then
     snorm_ener =  0.
       do inc=1,nwav
     zcin  = zci(inc)*rstar
                
     ze2 = zcin /(1.+ zcin**nslope3)

     snorm_ener = snorm_ener + ze2*zdci(inc)*rstar  !dim-less
         flux(inc,1) = ze2 * zcin           
       enddo

     ekin_norm = 1./snorm_ener
     
!   taub_src = sigu2 * rho_src * [v_kxw / zms ]
!   sigu2  = taub_src*zms/(rho_src/v_kxw)
!   ze1 = sigu2*ks*dens/Ns = taub*zms/Ns

        ze1 = taub_src*zms/bvi * ekin_norm
            taub_src = 0.
            
      do inc=1,nwav
       flux(inc,1) = ze1* flux(inc,1)
       taub_src = taub_src + flux(inc,1)*zdci(inc)
      enddo
      ze1 = ekin_norm * v_kxw * rstar2
      do jk=ksrc, ktop
         fden_bnen(jk) = rhoint(jk) / bn(jk) *ze1    ! mult on => sigu2(z)*cdf2 => flux_sat
      enddo
             
   endif
!
      do iaz=1,nazd
         fpu(iaz, ksrc) = taub_src
         fpu(iaz, km1)  = taub_src
      enddo

!     copy flux-1 into other azimuths
!     --------------------------------


      do iaz=2, nazd
        do inc=1,nwav
             flux(inc,iaz) = flux(inc,1)
          enddo
      enddo

!       if (mpi_id  == master .and. ener_norm == 1) then
!          print *
!          print *, 'vay_norm: ', taub_src,  taub_ch,  sigu2,  flux_norm, ekin_norm
!          print *
!       endif

       if (idebug_gwrms == 1) then
      pwrms =0.
      ptrms =0.      
      tx1 = real(nazd)/rhoint(ksrc)*rv_kxw
      ze2 = wrk2(ksrc)    !  (bvi*atm(ksrc)*rgrav)**2
     do inc=1, nwav
      v_kzw = bvi*rci(inc)
      ze1 = flux(inc,1)*zdci(inc)*tx1*v_kzw
      pwrms = pwrms + ze1
      ptrms = ptrms + ze1 * ze2
     enddo
       wrms(jl, ksrc) = pwrms
       trms(jl, ksrc) = ptrms 
        endif

!     --------------------------------
    wave_act(:,:) = 1.0
!                                        vertical do-loop
        do jk=ksrc, levs

       jkp = jk+1
!                                        azimuth do-loop
         do iaz=1, nazd

         sig_u2az_m(iaz) = sig_u2az(iaz)

          umfp = ui(iaz, jkp)
             umfm = ui(iaz, jk)
             umfc = .5*(umfm + umfp)
!                                        wave-cin loop
             dfdz_v(iaz, jk)  = 0.0
             dfdz_heat(iaz, jk)  = 0.0
             fpu(iaz, jkp)    = 0.0
         sig_u2az(iaz) =0.0
!
!             wave_dis(iaz, :) = vueff(jk)
          do inc=1, nwav
            flux_m(inc, iaz) = flux(inc, iaz)

            zcin  = zci(inc)          ! zcin =/0  by definition
            zcinc = rci(inc)

             if(wave_act(inc,iaz) == 1.0) then
!=======================================================================
! discrete mode
! saturated limit    wfit = kzw*kzw*kt; wfdt = wfit/(kxw*cx)*betat
! & dissipative      kzi = 2.*kzw*(wfdm+wfdt)*dzpi(k)
!=======================================================================

               v_cdp =  zcin - umfp
           v_cdp2=v_cdp*v_cdp
           cdf2  = v_cdp2 - c2f2
         if (v_cdp .le. ucrit_max .or. cdf2 .le. 0.0) then
!
! between layer [k-1,k or jk-jkp]  (Chi - Uk) -> ucrit_max, wave's absorption
!
            wave_act(inc,iaz) =0.
                akzw(inc, iaz, jkp) = pi/dz_meti(jk)   ! pi2/dzmet
        fluxs = 0.0                           !max(0., rhobnk(jkp)*ucrit3)*rdci(inc)
        flux(inc,iaz)   = fluxs

         else

              v_wdp = v_kxw*v_cdp
              wdop2 = v_wdp* v_wdp

!
! rotational cut-off
!
              kzw2 = (bn2(jkp)-wdop2)/Cdf2
!
!cires_ugwp_initialize.F90:      real, parameter :: mkzmin = pi2/80.0e3
! Lzmax = 80 km for mkzmin =6.28/Lzmax
!
          if ( kzw2 > mkz2min ) then
                v_kzw = sqrt(kzw2)
                akzw(inc, iaz, jkp) = v_kzw
!
!linsatdis:  kzw2, kzw3, kdsat, c2f2,  cdf2, cdf1
!
!kzw2 = (bn2(k)-wdop2)/Cdf2  - rhp4 - v_kx2w  ! full lin DS-NGW (N2-wd2)*k2=(m2+k2+[1/2H]^2)*(wd2-f2)
!              Kds_sat = kxw*Cdf1*rhp2/kzw3
!krad, kvg, kion, ktg
                v_cdp  = sqrt( cdf2 )
                v_wdp  = v_kxw *  v_cdp
        v_wdi = kzw2*vueff(jk) + kion(jk)                  ! UFSWAM -Ion drag coeff-nt is needed
        v_wdpc = sqrt(v_wdp*v_wdp +v_wdi*v_wdi)
        v_kzi  = v_kzw*v_wdi/v_wdpc

!
                ze1 = v_kzi*v_zmet(jk)

        if (ze1 .ge. 1.e-2) then
            expdis = max(exp(-ze1), 0.01)  ! ten-times attenuation of the flux/energy between 2-layers
        else
            expdis = 1.- ze1 + .5*ze1*ze1
        endif

!
                wave_act(inc,iaz) = 1.0
                fmode =  flux(inc,iaz)

               flux_2_sig = v_kzw/v_kxw/rhoint(jkp)
               w1 = v_wdpc/kzw2/v_kzw/v_zmet(jk)
           else                                           ! kzw2 <= mkz2min large "Lz"-reflection

                expdis = 1.0
                v_kzw  = mkzmin

                v_cdp  = 0.                               ! no effects of reflected waves with Lz > Lzmax =80 km
                wave_act(inc,iaz) = 0.0
                akzw(inc, iaz, jkp) = v_kzw
		
                fmode = 0.
                w1 =0.
             endif

              fdis  =  fmode*expdis*wave_act(inc,iaz)
!==============================================================================
!
! Saturated Fluxes and Energy: Spectral and Dicrete Modes
!
! S2003             fluxs= fden_bn(jk)*(zcin-ui(jk,iaz))**2/zcin
! WM2001            fluxs= fden_bn(jk)*(zcin-ui(jk,iaz))
! saturated flux + wave dissipation - Keddy_gwsat in UGWP-V1
!  linsatdis = 1.0 , here:   u'^2 ~ linsatdis* [v_cdp*v_cdp]
!
! old-sat                    fluxs= fden_bn(jkp)*cdf2*zcinc*wave_act(inc,iaz)
!                            fluxs= fden_bn(jkp)*cdf2*zcinc*wave_act(inc,iaz)
! new sat                    fluxs= fden_bn(jkp)*sqrt(cdf2)*wave_act(inc,iaz)
!
!                    fluxs= fden_bn(jkp)*sqrt(cdf2)*wave_act(inc,iaz)

!
!
! old spectral  sat-limit with "mapping to source-level"  sp_tau(cd) = fden_bn(jkp)*sqrt(cdf2)
! new spectral  sat-limit with "mapping to source-level"  sp_tau(cd) = fden_bn(jkp)*cdf2*rstar2
! [fden_bn(jkp)] = Pa/dc
! fsat  = rstar*(zcin*zcin) * [taub_src / SN * [ rstar3*rho/rho_src *N_src/N] = fden_bn ]

           if (ener_norm == 0) fluxs= fden_bn(jkp)*cdf2*wave_act(inc,iaz)     ! dim-n:  Pa/[m/s]
!
! single mode saturation limit: [rho(z)/bn(z)*kx *linsat2* cd^3] /dc
!
           if (ener_lsat == 1) fluxs= fden_Lsat(jkp)*cdf2*sqrt(cdf2)*rdci(inc)*wave_act(inc,iaz)

       if (ener_norm == 1) then

! spectral saturation limit
      
              if (ener_lsat == 0)  fluxs= fden_bnen(jk)*cdf2*wave_act(inc,iaz)*sig_u2az_m(iaz)

! single mode saturation limit: [rho(z)/bn(z)*kx *linsat2* cd^3] /dc

              if (ener_lsat == 1) fluxs= fden_Lsat(jkp)*cdf2*sqrt(cdf2)*rdci(inc)*wave_act(inc,iaz)
!
            endif
!----------------------------------------------------------------------------
!    dicrete mode saturation fden_sat(jkp) = rhoint(jkp)/bn(jkp)*v_kxw
!             fluxs = fden_sat(jkp)*cdf2*sqrt(cdf2)/zdci(inc)*L2sat
!             fluxs_src = fden_sat(ksrc)*cdf2*sqrt(cdf2)/zdci(inc)*L2sat
!----------------------------------------------------------------------------
              zdep = fdis-fluxs           ! dimension [Pa/dc] *dc = Pa
	      
              if(zdep > 0.0 ) then
! subs on sat-limit
                 ze1 = flux(inc,iaz)
                 flux(inc,iaz)    = fluxs
         ze2 = log(ze1/fluxs)*w1     ! Kdsat-compute damping of mode =>df = f-fluxs
                                     ! here we can add extra-dissip for the next layer
              else
! assign dis-ve flux
                 flux(inc,iaz) =   fdis
              endif

         dtau = flux_m(inc,iaz)-flux(inc,iaz)
         if (dtau .lt. 0) flux(inc,iaz)   = flux_m(inc,iaz)
         
!
! GW-sponge domain: saturate all "GW"-modes above "zsp_gw"
!
         if ( azmeti(jkp) .ge. zsp_gw) then
            mi_sponge = .5/dz_meti(jk)
            ze2 = v_wdp /v_kzw * mi_sponge                     ! Ksat*v_kzw2 = [mi_sat*wdp/kzw]
            v_wdi = ze2 + v_wdi*0.25                           ! diss-sat GW-sponge
            v_wdpc = sqrt(v_wdp*v_wdp +v_wdi*v_wdi)
            v_kzi  = v_kzw*v_wdi/v_wdpc
            ze1 = v_kzi*v_zmet(jk)
            exp_sponge = exp(-ze1)
!
! additional sponge near the model top
!
             flux(inc,iaz) =  flux(inc,iaz) *exp_sponge
          endif

         endif  !  coriolis or CL condition-checkif  => (v_cdp .le. ucrit_max) then
	  
        endif   !  only for waves w/o CL-absorption/refelections   wave_act=1


!
! sum for given (jk, iaz)  all active "wave-mode" contributions
!
        if (wave_act(inc,iaz) == 1) then

          zcinc =zdci(inc)
              vc_zflx_mode      = flux(inc,iaz)
              vmdiff = max(0., flux_m(inc,iaz)-vc_zflx_mode)
              if (vmdiff <= 0. ) vc_zflx_mode = flux_m(inc,iaz)
              ze1 = vc_zflx_mode*zcinc
              fpu(iaz, jkp) = fpu(iaz,jkp)  + ze1                ! flux (pa)   at
            sig_u2az(iaz) = sig_u2az(iaz) + ze1*flux_2_sig       ! ekin(m2/s2) at z+dz

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! (heat deposition integration over spectral mode for each azimuth
!      later sum over selected azimuths as "non-negative" scalars)
!      cdf1 = sqrt( (zci(inc)-umfc)**2-c2f2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             zdelp = wrk3(jk)*cdf1 *zcinc

             zdelp = wrk3(jk)* v_cdp *zcinc * vmdiff


!          zcool = 1.                           !  COOL=(-3.5 + Pr)/Pr
!          zcool = [Kv/Pr]*N2*(Pr-Cp/R)/cp
!             edis  = (c-u)*ax/cp = Kv_dis*N2/cp
!             cool  = -Kt*N2/R
!  add heat-conduction "bulk" impact: 1/Pr*(g*g*rho)* d [rho*Kv(dT/dp- R/Cp *T/p)]
!
              dfdz_v(iaz, jk) = dfdz_v(iaz,jk)  +  zdelp          ! +cool  !heating & simple cooling  < 0
              dfdz_heat(iaz, jk) = dfdz_heat(iaz,jk) + zdelp      ! heating -only      > 0
	      
            endif !wave_act(inc,iaz) == 1) 
!
          enddo      ! wave-summation   loop

          ze1 =fpu(iaz, jk)
          if (fpu(iaz, jkp) > ze1 ) fpu(iaz, jkp) = ze1
!
! compute wind and temp-re rms
!
    if (idebug_gwrms == 1) then
       pwrms =0.
       ptrms =0.                
      do inc=1, nwav
      if (wave_act(inc,iaz) > 0.) then
           v_kzw =akzw(inc, iaz, jk)
           ze1 = flux(inc,iaz)*v_kzw*zdci(inc)*wrk1(jk)
       pwrms = pwrms + ze1
       ptrms = ptrms + ze1*wrk2(jk)
      endif 
     enddo
           Awrms(iaz, jk) = pwrms   ! wind-rms
           Atrms(iaz, jk) = ptrms   ! temp-rms
    endif

! --------------
        enddo                              ! end  Azimuth do-loop

!
! eddy wave dissipation  to  limit GW-rms from current "dz" =>"dz_up", next layer
!
           tx1 = sum(abs(dfdz_heat(1:nazd, jk)))/bn2(jk)
           ze1=max(dked_min, tx1)
           ze2=min(dked_max, ze1)
           vueff(jkp) = ze2 + vueff(jkp)
!
       enddo                               ! end Vertical do-loop
!
! top-layers constant interface-fluxes and zero-heat
! we allow non-zero momentum fluxes and thermal effects
!            fpu(1:nazd,levs+1) = fpu(1:nazd, levs)
!            dfdz_v(1:nazd, levs) = 0.0

! ------------------------------------------------------------
!   sum contributions for total zonal and meridional fluxes +
!           energy dissipation
! ------------------------------------------------------------
!
!========================================================================
! at the source level and below taux = 0 (taux_E=-taux_W by assumption)
!========================================================================

         do jk=ksrc,  levs
            taux(jk) = 0.0
            tauy(jk) = 0.0
           do iaz=1,nazd
             taux(jk)  = taux(jk)  + fpu(iaz,jk)*zcosang(iaz)
             tauy(jk)  = tauy(jk)  + fpu(iaz,jk)*zsinang(iaz)
             pdtdt(jl,jk) = pdtdt(jl,jk) + dfdz_v(iaz,jk)
             dked(jl,jk) =  dked(jl,jk)  + dfdz_heat(iaz,jk)
          enddo
         enddo
             jk = ktop; taux(jk)=0.; tauy(jk)=0.
          do iaz=1,nazd
             taux(jk)  = taux(jk)  + fpu(iaz,jk)*zcosang(iaz)
             tauy(jk)  = tauy(jk)  + fpu(iaz,jk)*zsinang(iaz)
          enddo

       if (idebug_gwrms == 1) then
         do jk=kp1,  levs
           do iaz=1,nazd
            wrms(jl,jk)  =wrms(jl,jk)   + Awrms(iaz,jk)
            trms(jl,jk)  =trms(jl,jk)   + Atrms(iaz,jk)
            tauabs(jl,jk)=tauabs(jl,jk) + fpu(iaz,jk)
           enddo
         enddo
       endif
!

       do jk=ksrc+1,levs
          jkp = jk + 1
           zdelp = wrk3(jk)*gw_eff
           ze1   = (taux(jkp)-taux(jk))* zdelp
           ze2   = (tauy(jkp)-tauy(jk))* zdelp

           if (abs(ze1) >= maxdudt ) then
             ze1 = sign(maxdudt, ze1)
           endif
           if (abs(ze2) >= maxdudt ) then
             ze2 = sign(maxdudt, ze2)
           endif

           pdudt(jl,jk) = -ze1
           pdvdt(jl,jk) = -ze2
!
! Cx =0 based Cx=/= 0. above
!
!
        if (knob_ugwp_doheat == 1) then
!
!maxdtdt=  dked_max * bnfix2
!
           pdtdt(jl,jk) = pdtdt(jl,jk)*gw_eff
           ze2 = pdtdt(jl,jk)
           if (abs(ze2) >= max_eps ) pdtdt(jl,jk) = sign(max_eps, ze2)

           dked(jl,jk) =  dked(jl,jk)/bn2(jk)   ! "dked-rhs" contains only GW-heating, we don't use K-Lindzen diss-n   
           ze1  = max(dked_min, dked(jl,jk))
           dked(jl,jk)  = min(dked_max, ze1)
           qmid(jk) = pdtdt(j,jk)
        endif
	
       enddo
!----------------------------------------------------------------------------------
! Update heat = ek_diss/cp and aply 1-2-1 smoother for "dked" => dktur
!               here with "u_new = u +dtp*dudt ; vnew = v + v +dtp*dvdt
!               can check "stability" in the column and "add" ktur-estimation
!               to suppress instability as needed so dked = dked_gw + ktur_ric
!----------------------------------------------------------------------------------
   
             dktur(1:levs) = dked(jl,1:levs)
!
        do ist= 1, nstdif
             do jk=ksrc,levs-1
              adif(jk)  =.25*(dktur(jk-1)+ dktur(jk+1)) + .5*dktur(jk)
             enddo
              dktur(ksrc:levs-1) = adif(ksrc:levs-1)
        enddo
       
         dktur(levs) = .5*( dked(jl,levs)+ dked(jl,levs-1))
         dktur(levs+1) = dktur(levs)

       do jk=ksrc,levs+1
          ze1 = .5*( dktur(jk) +dktur(jk-1) )
          kvint(jk) = ze1
          ktint(jk) = ze1*iPr_ktgw
        enddo

!
! Thermal budget qmid = qheat + qcool
!
       do jk=ksrc+1,levs
           ze2  = qmid(jk) + (dktur(jk)*Akt(jk) + grav_vari(jk)*(ktint(jk+1)-ktint(jk))/dz_meti(jk))*fg_cool
           qmid(jk) = ze2
         if (abs(ze2) >= max_eps ) qmid(jk) = sign(max_eps, ze2)
	   arcpd(jk) = 2./(cp_msp(jk)+cp_msp(jk-1))
           pdtdt(jl,jk) = qmid(jk) * arcpd(jk)
           dked(jl, jk) = dktur(jk)
         enddo
!
! perform explicit eddy "diffusive" 3-point smoothing of "u-v-t"
! from the surface/launch-gw to the "top"
!
!
! update by source function X(t+dt) = X(t) + dtp * dXdt
!
      uold(km2:levs) = aum(km2:levs)+pdudt(jl,km2:levs)*dtp
      vold(km2:levs) = avm(km2:levs)+pdvdt(jl,km2:levs)*dtp
      told(km2:levs) = atm(km2:levs)+pdtdt(jl,km2:levs)*dtp
      pkold(km2:levs) = exp(-pkappa*alog(aprsl(km2:levs)))
      ptold(km2:levs) = told(km2:levs)*pkold(km2:levs)
!
! diagnose turb-profile using "stability-check" relying on the free-atm diffusion
!          sc2 = 30m x 30m
!
           dktur(km2:levs) = dked_min

      do jk=km1,levs
            uz = uold(jk) - uold(jk-1)
            vz = vold(jk) - vold(jk-1)
            ze1    = dz_met(jk)
            zdelm  = 1./ze1

            tvc = told(jk)     !* (1. +fv*aqm(jk))
            tvm = told(jk-1)   !* (1. +fv*aqm(jk-1))
            zthm          = 2.0 / (tvc+tvm)
            shr2 = (max(uz*uz+vz*vz, dw2min)) * zdelm *zdelm

!UFS            bn2(jk)    = grav2cpd*zthm  * (1.0+rcpdl*(tvc-tvm)*zdelm)

             cpint = .5*(cp_msp(jk)+cp_msp(jk-1))
	     bn2(jk)    = grav_vari(jk)*zthm*(grav_vari(jk)/cpint+ (tvc-tvm)*zdelm)
	     	    
	     if (bn2(jk) < 0. ) then 
	       if (debug_conv) then	     
	         print *, 'GWPV1-OUT CINSTAB', tvc, tvm, nint(azmet(jk)*1.e-3)
	         print *, 'GWPV1-OUT CIN_RIC', bn2(jk)/shr2, ' lat ', nint(xlatd(j))
	       endif
	   endif

	    
            bn2(jk)    = max(min(bn2(jk), bnv2max), bnv2min)
	    
           rdi_vari   = 2*Runiv/(mu_msp(jk)+mu_msp(jk-1))
	   rhp_wam = 0.25*zthm/rdi_vari*grav_vari(jk)
           zmetk  =  azmet(jk)* rhp_wam                     ! mid-layer height k_int => k_int+1
           zgrow = exp(zmetk)
           ritur = bn2(jk)/shr2
	   zlturb = min(lturb*zmetk, zdelp*.5)
           kamp = sqrt(shr2)* zlturb*zlturb
           w1 = 1./(1. + 5.*ritur)
	   w1 = w1*w1
	   if (hb93_kt) then
	     w1 = 1./(1. + 10.*ritur)/(1.+8.*ritur)
	   endif
! Smag-type of eddy diffusion K_smag = Sqrt(Deformation - N2/Pr)* L2 *const
!
            kamp = kamp * w1 * w1
            ktur= min(max(kamp, dked_min), dked_max)
	    
            dktur(jk) = ktur*fg_ktur    ! eddy mixing associated with "new" U-V-T profiles
!
! update of dked = dked_gw  + k_turb_mf
!                
         dked(jl, jk) = dked(jl, jk) +dktur(jk)

      enddo

	     
!
! apply eddy effects due to GWs:  explicit scheme Kzz*dt/dz2 < 0.5 stability
!
      if (knob_ugwp_dokdis == 2) then

        do jk=ksrc,levs
          ze1 = min(.5*(dktur(jk) +dktur(jk-1)), dturb_max)
          kvint(jk) = kvint(jk) + ze1                    ! add dktur-background eddy mixing for "GW-adjusted" profiles
        enddo 
          kvint(km1) =  kvint(ksrc)
          kvint(ktop) = kvint(levs)

        dzmetm =  1./dz_met(km1)
        Adif(km1:levs)  =  0.
        Cdif(km1:levs)  =  0.
          do jk=km1,levs-1

           dzmetp = 1./dz_met(jk+1)
           dzmetf = 1./(dz_meti(jk)*rhomid(jk))


        ktur = kvint(jk)  *rhoint(jk)   * dzmetf
        kturp =Kvint(jk+1)*rhoint(jk+1) * dzmetf
                
        Adif(jk) = ktur  * dzmetm
        Cdif(jk) = kturp * dzmetp
        ApC = adif(jk)+cdif(jk)
        ACdif(jk) = ApC

        w1 =  ApC*iPr_max
        if (rdtp  <  w1 ) then
          Anstab(jk) = floor(w1*dtp) + 1
        else
          Anstab(jk) = 1
        endif
          dzmetm = dzmetp
       enddo

        nstab = maxval( Anstab(ksrc:levs-1))

!        if (nstab .ge. 3) print *, 'nstab ', nstab
!
! Attention "k" instead "jk" vertical loop
!
        dtdif = dtp/real(nstab)
          ze1 = 1./dtdif

        do ist= 1, nstab
          do k=ksrc,levs-1  
            Bdif   = ze1 - ACdif(k)
            Bt_dif = ze1 - ACdif(k)* iPr_ktgw                                  ! ipr_Ktgw = 1./Pr <1
            unew(k) = uold(k)*Bdif  + uold(k-1)*Adif(k) + uold(k+1)*Cdif(k)
            vnew(k) = vold(k)*Bdif  + vold(k-1)*Adif(k) + vold(k+1)*Cdif(k)
!            tnew(k) = told(k)*Bt_dif+(told(k-1)*Adif(k) + told(k+1)*Cdif(k))*iPr_ktgw
            tnew(k) = ptold(k)*Bt_dif+(ptold(k-1)*Adif(k) + ptold(k+1)*Cdif(k))*iPr_ktgw	    
         enddo

           uold(ksrc:levs-1) = unew(ksrc:levs-1)*dtdif    ! value du/dtp *dtp = du
           vold(ksrc:levs-1) = vnew(ksrc:levs-1)*dtdif
           told(ksrc:levs-1) = tnew(ksrc:levs-1)*dtdif/pkold(ksrc:levs-1)
!
! smoothing the boundary points: "k-1" = ksrc-1 and  "k+1" = levs
!
           uold(levs) = uold(levs-1)
           vold(levs) = vold(levs-1)
           told(levs) = told(levs-1)
            enddo
!
! compute "smoothed" tendencies by molecular + GW-eddy diffusions
!
             do k=ksrc,levs-1
!         
! final updates of tendencies and diffusion
!
          ze2 = rdtp*(uold(k) - aum(k))
          ze1 = rdtp*(vold(k) - avm(k))
	  
           if (abs(ze1) >= maxdudt ) then
             ze1 = sign(maxdudt, ze1)
           endif
           if (abs(ze2) >= maxdudt ) then
             ze2 = sign(maxdudt, ze2)
           endif
!
           pdudt(jl, k) = ze2
           pdvdt(jl, k) = ze1
!--------	   
! pdtdt	   
               uz = uold(k+1) - uold(k-1)
               vz = vold(k+1) - vold(k-1)
           ze2 = 1./(dz_met(k+1)+dz_met(k) )
	   
           mf_diss_heat = arcpd(k)*kvint(k)*(uz*uz +vz*vz)*ze2*ze2  ! vert grad local heat GW-deddy
	   
          ze1= rdtp*( told(k) - atm(k) )		    
          ze2= ze1 + mf_diss_heat              ! extra heat due to eddy viscosity
          if (abs(ze2) >= maxdtdt ) ze2 = sign(maxdtdt, ze2)	  
!          pdtdt(jl,k)= ze2                                         !nocool noheat
!----
         enddo


      ENDIF    !  dissipative IF-loop for vertical eddy difusion u-v-t

       enddo        ! J-loop
!
       RETURN

!=================================   diag print after "return" ======================
       if (kdt ==1 .and. mpi_id == master) then
!
          print *, ' ugwpv1: nazd-nw-ilaunch=', nazd, nwav,ilaunch, maxval(kvg), ' kvg '
          print *,  'ugwpv1: zdci(inc)=' ,  maxval(zdci), minval(zdci)
          print *,  'ugwpv1: zcimax=' , maxval(zci) ,' zcimin=' , minval(zci)
!         print *,  'ugwpv1: tau_ngw=' , maxval(taub_src)*1.e3,  minval(taub_src)*1.e3, tau_min

          print *

        endif

       if (kdt == 1 .and. mpi_id == master) then
         print *, 'vgw done nstab ', nstab
!
         print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw ax  ugwp'
         print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw ay  ugwp'
         print *, maxval(dked)*1.,  minval(dked)*1,  'vgw keddy m2/sec ugwp'
         print *, maxval(pdtdt)*86400.,  minval(pdtdt)*86400,'vgw eps  ugwp'
!
!        print *, ' ugwp -heating rates '
       endif
!=================================
       return       
     end subroutine ugwpv1_ngw_solv_ufswam
!     
     subroutine cires_dry_adjust(levs, t, kappa, prsl, prsi, delp)
!     
! vay oct 2023 for conv-instability in the Mesosphere and Thermosphere 
! make dry adjustment of tkin(levs) +> convectively stable configuration
! "warm the conv-unstable layers": relax unstable T-re profile => T-stab
!  Apply GW-scheme only for in the conv-stable atmosphere
!  
      use machine , only : kind_phys

      implicit none
      integer             , parameter ::    k_trop = 40       ! SMT-domain
      integer             , parameter ::    niter =  15                
      integer, intent(in)  :: levs
      real(kind=kind_phys),intent(in)    :: kappa(levs), prsl(levs), delp(levs)
      real(kind=kind_phys),intent(in)    :: prsi(levs+1)      
      real(kind=kind_phys),intent(inout) ::  t(levs)  
! local
!
!
      real(kind=kind_phys)     ::    gammad, dtdp , zeps, zepsdp, zgamma, rdenom    
      logical                  ::     dodad, ilconv 
      integer                  ::   k, jiter
      real(kind=kind_phys), dimension(levs)     ::   c1dad, c2dad, c3dad, c4dad
      real(kind=kind_phys)     :: dtm, dtins
! from exobase
! to k_trop define the "instability" in the column
!        

       zeps = 2.0e-5     ! set convergence criteria     
       dtdp = (t(levs-1)-t(levs))/(prsl(levs-1)-prsl(levs))
       gammad = .25*(kappa(levs-1)+kappa(levs))/(t(levs-1)+t(levs))/prsi(levs)  
       dodad = (dtdp + zeps) .gt. gammad
       
   do k=levs-1,k_trop, -1   
         gammad = .25*(kappa(k-1)+kappa(k))/(t(k-1)+t(k))/prsi(k)
         dtdp = (t(k-1) - t(k))/(prsl(k-1) - prsl(k))
         dodad = dodad .or. (dtdp + zeps).gt.gammad
   end do
   
      if( .not. dodad) return
!
! dodat=.T. => make dry adjustment can be "multiple" layers in FV3       
         zeps = 2.0e-5
         do k=levs-1,k_trop, -1
            c1dad(k) = 0.25*(kappa(k-1)+kappa(k))*(prsl(k-1)-prsl(k))/prsi(k)
            c2dad(k) = (1. - c1dad(k))/(1. + c1dad(k))
            rdenom = 1./(delp(k)*c2dad(k) + delp(k-1))
            c3dad(k) = rdenom*delp(k)
            c4dad(k) = rdenom*delp(k-1)
         end do
	 
50       do jiter=1,niter
            ilconv = .true.
            do k=levs-1,k_trop, -1
               zepsdp = zeps*(prsl(k-1) - prsl(k))
               zgamma = c1dad(k)*(t(k) + t(k-1))
	       dtins = zgamma + zepsdp
	       dtm = t(k-1)-t(k)
               if ( dtm >= dtins) then
                  ilconv = .false.
                  t(k-1) = t(k)*c3dad(k) + t(k-1)*c4dad(k)
                  t(k) = c2dad(k)*t(k-1)
               endif
            end do
            if (ilconv) go to 80 ! convergence => next longitude
         end do
!
! Crude (double) convergence criterion if no convergence in niter iterations
!
          zeps = zeps + zeps
	 
         if (zeps > 1.e-4) then
            write(6,*)'DADADJ: No convergence in dry adiabatic adjustment'
           else
            write(6,810) zeps
            go to 50
         endif
      
80    continue
!
!
!800   format(' lat,lon = ',2i5,', zeps= ',e9.4)
810   format(//,'UGWPP-V1 DADADJ: Convergence doubled to EPS=',E9.4)
	     
     return
     end subroutine cires_dry_adjust
     
     subroutine get_ufswam_med(levs, msp_dim, no1, no2, nn2, nhe, nair, tkin, pmid,  &
	                      zm, zmi, cp_msp, ktmol, vumol, mu_msp, grav_var, grav_vari)
			      
! Banks & Kockarts, 1973: ktmol, vumol	

     implicit none	 
      integer, intent(in):: levs, msp_dim
      real(kind=kind_phys), intent(in)      :: tkin(levs), pmid(levs), zm(levs) , zmi(levs+1)      
      real(kind=kind_phys), intent(in)      :: no2(levs), no1(levs), nn2(levs), nhe(levs), nair(levs)  
      
      real(kind=kind_phys), intent(out)     :: vumol(levs+1), ktmol(levs+1), grav_vari(levs+1)  
      real(kind=kind_phys), intent(out)     :: cp_msp(levs), mu_msp(levs),   grav_var(levs)
!
!  local
!    
     integer :: k, kk      
     real(kind=kind_phys) :: amu, muk, lak, rhok, t69k, o_n, o2_n, n2_n, he_n      
     real(kind=kind_phys) :: rarad
!
! bottom layer in physics k=1    avgdbz= avgd * bz
!  
     kk = 1
     k =  1
     rarad  = 1./arad  
     
       grav_vari(kk) = grav/(1.+zmi(kk)*rarad)/(1.+zmi(kk)*rarad) 
       cp_msp(k) = cpo*no1(k) + cpo2*no2(k) + cpn2*nn2(k)+ cphe*nhe(k)
       amu = no1(k)/amo + no2(k)/amo2 + nn2(k)/amn2 + nhe(k)/amhe
       mu_msp(k) = 1./amu
       
       o_n = no1(k)*mu_msp(k)/amo
       o2_n= no2(k)*mu_msp(k)/amo2
       n2_n= nn2(k)*mu_msp(k)/amn2
       he_n= nhe(k)*mu_msp(k)/amhe
       
        muk = o_n*muo + o2_n*muo2 + n2_n*mun2
	lak = (o_n*lao + o2_n*lao2 + n2_n*lan2)/cp_msp(k)
	t69k = tkin(k) ** 1.69
	
	rhok =avgdbz/pmid(k)*amu  ! surface value
	  
        vumol(kk) = muk*t69k*rhok
        ktmol(kk) = lak*t69k*rhok
	
		  
     do k=2, levs
        cp_msp(k) = cpo*no1(k) + cpo2*no2(k) + cpn2*nn2(k)+ cphe*nhe(k)
	amu = no1(k)/amo + no2(k)/amo2 + nn2(k)/amn2 + nhe(k)/amhe
	mu_msp(k) = 1./amu
!vmr from mmr	

       o_n = no1(k)*mu_msp(k)/amo
       o2_n= no2(k)*mu_msp(k)/amo2
       n2_n= nn2(k)*mu_msp(k)/amn2
       he_n= nhe(k)*mu_msp(k)/amhe
       
        muk = o_n*muo + o2_n*muo2 + n2_n*mun2
	lak = (o_n*lao + o2_n*lao2 + n2_n*lan2)/cp_msp(k)
	t69k = tkin(k) ** 1.69
	
	rhok =avgdbz/pmid(k)*amu
	  
        vumol(k) = .5*(vumol(k-1)+muk*t69k*rhok)
        ktmol(k) = .5*(ktmol(k-1)+lak*t69k*rhok)
	
	grav_var(k)  = grav/(1.+zm(k)*rarad)/(1.+zm(k)*rarad)
	grav_vari(k) = grav/(1.+zmi(k)*rarad)/(1.+zmi(k)*rarad)
     enddo 
     
       k= levs+1
       grav_vari(k) = grav/(1.+zmi(k)*rarad)/(1.+zmi(k)*rarad) 
       vumol(k) = vumol(k-1)
       ktmol(k) = ktmol(k-1)
           	 
     end subroutine get_ufswam_med
     
     
     	 
end module cua_ugwpv1_ngw_ufswam
