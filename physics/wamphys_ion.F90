!--------------------------------------------------------------------
!      module wamphys_ion
!====================================================================
! apr 06 2012 henry juang, initial implement for nems
! dec    2012    jun wang, move init out of column physics
! nov    2015   vay  sub-ne idea_ion_init(levs) & idea_ion_input.f 
! nov    2016   vay  upgrades related to  idea_ion_input.f
!                                          idea_ion_empirmodels.f   
! jul    2017   zhuxiao li, rashid and tim add the joule heating factor(jh_fac) to 
!               include the seasonal variation and semiannual variation of the
!               joule heating.
! mar    2018   zhuxiao li and tzu-wei fang add the new features to put
! the driving parameters into the idea_geteb before call get_efield. 
!
! mar 2022      svetlana karol ccpp-version of idea_ion.f
!                wamphys_ion.f90 (without mpi-statemments and kind_phys precision)
!                separates "idea_ion_init" from "idea_ion" and ion_empir_models
!===================
module wamphys_ion
      
      use  machine,             only : kind_phys
      use  wamphys_const,       only : rad_to_deg, deg_to_rad, pid12, pi2_365d 
      use  wamphys_const,       only : pi, pi2, pi_24hr 
         
!      use wamphys_diag2d, only : el_diag, shal_diag, sped_diag 
      implicit none
        
      contains
!
 subroutine wam_ion_run(im,levs, jdat,                     &
 pres,solhr,cospass,zg,grav,o_n,o2_n,n2_n,cp,              &    
 adu,adv,adt,dudt,dvdt,dtdt,rho,rlat,rlon,                 & 
 dayno,utsec,sda,maglon,maglat,btot,dipang,essa,           & 
 f107, f107d, kp, nhp, nhpi, shp, shpi,                    & 
 swbt, swang, swvel, swbz, swden )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! driver      dtdt(i,k)=jh(i,k)/cp(i,k), dudt dvdt
!              ion darge and joule heating
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      use wamphys_init_module    , only :jh0, jh_semiann,   &
                                         jh_ann, jh_tanh,   & 
                                         jh_st0, jh_st1
      use wamphys_init_module    , only : swin_drivers, spw_drivers					 
      use wamphys_init_module    , only : k91 
      use wamphys_const    ,       only : fac_lst      
      implicit none
      

      integer, intent(in)     :: im !number of logitude
      integer, intent(in)     :: levs ! number of pres grid
      integer, intent(in)     :: dayno !calender day
      integer, intent(in)     ::  jdat(8)
      real(kind=kind_phys), intent(in)     :: pres(im,levs) ! pressure, pa

      real(kind=kind_phys), intent(in)     :: o_n(im,levs)  ! number density o (/m3)
      real(kind=kind_phys), intent(in)     :: o2_n(im,levs)
      real(kind=kind_phys), intent(in)     :: n2_n(im,levs)
      real(kind=kind_phys), intent(in)     :: rho(im,levs)  ! mass density (kg/m3)
      real(kind=kind_phys), intent(in)     :: zg(im,levs)   !  height (m)
      real(kind=kind_phys), intent(in)     :: grav(im,levs) !  gravity (m/s**2)
      real(kind=kind_phys), intent(in)     :: cp(im,levs)   !  (j/kg/k)
      real(kind=kind_phys), intent(in)     :: cospass(im)! cos solar zenith angle (rad) 
      real(kind=kind_phys), intent(in)     :: rlat(im) ! latitude (rad)
      real(kind=kind_phys), intent(in)     :: rlon(im) ! longitude (rad)
      real(kind=kind_phys), intent(in)     :: solhr     ! universal time (h)
  
! state

      real(kind=kind_phys), intent(in)     :: adt(im,levs)  ! temperature (k)
      real(kind=kind_phys), intent(in)     :: adu(im,levs)  ! zonal wind (m/s)
      real(kind=kind_phys), intent(in)     :: adv(im,levs)  ! meridional wind (m/s)

! input magnetic and electric parameters 

      real(kind=kind_phys), intent(in) :: maglon(im)  !magnetic longitude (rad)
      real(kind=kind_phys), intent(in) :: maglat(im)  !magnetic latitude (rad)
      real(kind=kind_phys), intent(in) :: btot(im)    !mapgnetic field strength
      real(kind=kind_phys), intent(in) :: dipang(im)  !dip angle (degree)
      real(kind=kind_phys), intent(in) :: essa(im)    !magnetic local time
      real(kind=kind_phys), intent(in) :: sda         ! solar declination angle (rad)
      real(kind=kind_phys), intent(in) :: utsec       !universal time
      
      real(kind=kind_phys)  , intent(in)   :: f107, f107d, kp  
					       
      real(kind=kind_phys)  , intent(in)   :: nhp, nhpi, shp, shpi            ! solar-geo inputs from wam_f107_kp.txt
      real(kind=kind_phys)  , intent(in)   :: swbz, swvel, swbt, swang, swden ! solar-geo inputs from wam_f107_kp.txt
   
! output
      real(kind=kind_phys), intent(out)     :: dtdt(im,levs)  ! temperature change (k/s)
      real(kind=kind_phys), intent(out)     :: dudt(im,levs)  ! zonal wind change (m/s2)
      real(kind=kind_phys), intent(out)     :: dvdt(im,levs)  ! meridional change wind (m/s2)

! local
! define jh_fac, storm_fac zhuxiao li 

      real(kind=kind_phys) :: rlt(im),sza(im),jh(im,levs),rinc(5),jh_fac, st_fac, vbz
      real(kind=kind_phys) :: rrho  
      
!diags                            
      real(kind=kind_phys), dimension(im,levs)   ::  el2d, sped2d, shal2d
      integer   i,k

! get vbz swvel*swbz

      vbz = swvel*swbz

! get sza in rad
      sza=acos(cospass)
      rlt =(solhr +rlon*fac_lst)*pi_24hr  
          
! get local solar time in rad:  rlt=(rlon/(15.*pi/180.)+solhr)/24.*2.*pi
!                               rlt=( rlon*r_2_d/15. + solhr )*(2.*pi/24.)
! get ion_drag & joule heating


       
      call wgetionparams(im,levs,k91,jdat, pres,                        &
        dayno,utsec,f107,f107d,kp,nhp,nhpi,spw_drivers,sda,sza,rlat,zg, &
        grav,o_n, o2_n, n2_n,adu,adv,adt,rho,rlt,rlon,                  &    
        btot,dipang,maglon,maglat,essa,swbt,swang,swvel,swbz,swden,     &
        dudt,dvdt,jh) 

! update to ........ k/sec
      do i=1,im

!   joule heating factor to consider the seasonal variation and
!   semiannual variation, zhuxiao.li

!  jh0_5
!           jh_fac = 1.75+0.5*tanh(2.*rlat(i))*cos((dayno+9.)*pi2_365d)    &
!                       +0.5*(cos(2.*pi2_365d*(dayno-80.)))

        jh_fac = jh0+jh_tanh*tanh(2.*rlat(i))*cos((dayno+9.)*pi2_365d)     &
                   +jh_semiann*cos(2.*pi2_365d*(dayno-80.))

! vbz adjustment

        if (abs(vbz).le.jh_st1) then
           st_fac = 1.
        else
           st_fac = (jh_st0+jh_st1)/(jh_st0+ abs(vbz))
        endif
	  
!       st_fac = 1.
	  jh_fac=jh_fac*st_fac
        do k=1,levs
	     rrho = 1./rho(i,k)
           dudt(i,k) = rrho*dudt(i,k)
           dvdt(i,k) = rrho*dvdt(i,k)	   
           dtdt(i,k) = rrho*jh(i,k)/cp(i,k)*jh_fac
        enddo
      enddo

      return
  end subroutine wam_ion_run

!=================================

  subroutine wgetionparams(im,levs,lev1, jdat, pres,                   &
    dayno,utsec,f107,f107d,kp,hp,hpi,spw_drivers,sda,sza,rlat,ht,grav, &
    o_n, o2_n, n2_n,adu,adv,adt,rho,rlt,rlon,                          &
    btot,dipang,maglon,maglat,essa,swbt,swang,swvel,swbz,swden,        & 
    dudt,dvdt,jh)
	   
      use  wamphys_set_data_ion,   only : emaps, cmaps, djspectra
      use  wam_efield_setdef_data, only : nmlon, nmlat, ylonm
      
      use  wam_efield_setdef_data, only : ed1, ed2
      use  wamphys_init_module, only  : gw_fix, tiros_activity_fix, tiros_switch 
      use wamphys_const,        only  : elch

!      use idea_ion_input, only : emaps1, cmaps1, djspectra1
!      use idea_ion_input, only : tiros_switch

      use wamphys_ion_empirmodels, only : wam_earth_chiu_model, wam_tiros_ionize
      use wamphys_ion_empirmodels, only : wam_tiros_ionize_data      

      implicit none
      integer,    intent(in)     :: dayno  !day 
      integer,    intent(in)     :: jdat(8) 
      integer,    intent(in)     :: im !number of logitude
      integer,    intent(in)     :: levs ! number of pres grid
      integer,    intent(in)     :: lev1 ! lowest pres level to start k91 from idea_comp...
!
!   
      real(kind=kind_phys), intent(in)     :: pres(im,levs) ! pressure, pa
      real(kind=kind_phys), intent(in)     :: o_n(im,levs)  ! number density o (/m3)
      real(kind=kind_phys), intent(in)     :: o2_n(im,levs)
      real(kind=kind_phys), intent(in)     :: n2_n(im,levs)
      real(kind=kind_phys), intent(in)     :: adt(im,levs)  ! temperature (k)
      real(kind=kind_phys), intent(in)     :: adu(im,levs)  ! zonal wind (m/s)
      real(kind=kind_phys), intent(in)     :: adv(im,levs)  ! meridional wind (m/s)
      real(kind=kind_phys), intent(in)     :: rho(im,levs)  ! mass density (kg/m3)
      real(kind=kind_phys), intent(in)     :: ht(im,levs)   ! geopotential height (m)
      real(kind=kind_phys), intent(in)     :: grav(im,levs) !  gravity (m/s**2)
      real(kind=kind_phys), intent(in)     :: f107, f107d, kp 
      real(kind=kind_phys), intent(in)     :: hp, hpi
      
      character,intent(in) :: spw_drivers
      
      real(kind=kind_phys), intent(in)     :: sda      ! solar declination angle (rad)
      real(kind=kind_phys), intent(in)     :: sza(im)  ! solar zenith angle (rad) 
      real(kind=kind_phys), intent(in)     :: rlt(im)  ! local time (rad) 
      real(kind=kind_phys), intent(in)     :: rlat(im) ! latitude (rad)
      real(kind=kind_phys), intent(in)     :: rlon(im) ! longitude (rad)
      real(kind=kind_phys), intent(in)     :: utsec    ! universal time (s)
!     real(kind=kind_phys), intent(in) :: elx(im)
!     real(kind=kind_phys), intent(in) :: ely(im)     !electric field
      real(kind=kind_phys), intent(in) :: maglon(im)  !magnetic longitude (rad)
      real(kind=kind_phys), intent(in) :: maglat(im)  !magnetic latitude (rad)
      real(kind=kind_phys), intent(in) :: btot(im)    !mapgnetic field strength
      real(kind=kind_phys), intent(in) :: dipang(im)  !dip angle (degree)
      real(kind=kind_phys), intent(in) :: essa(im)    !magnetic local time
      real(kind=kind_phys), intent(in) :: swbt, swang, swvel, swbz, swden


      real(kind=kind_phys), intent(out) :: dvdt(im,levs) !(m/s2)
      real(kind=kind_phys), intent(out) :: dudt(im,levs) !(m/s2)
      real(kind=kind_phys), intent(out) :: jh(im,levs)   ! (j/kg/s)
! local

      real(kind=kind_phys) :: ht1(levs),v1(levs),nden(levs),o2n(levs),on(levs),  &         
         n2n(levs),elx(im),ely(im),ssa,elz(im),ee1(im),                          & 
         ee2(im),cosdif,sindif,sdip,cdip,btheta,bphi,elecx,                      &
         elecy,dif,dlat,dlon
	 
!     real(kind=kind_phys) :: ut, bz, by, ed1, ed2, 
      integer              :: iday
      	 
      integer  ::  k,i,n,   tiros_activity_level 
      integer  :: ii  
!     ion drag variables :
! 
!     teff(levs)      1d local array of temperature
!     pion1(levs)     number density o+
!     pion2(levs)     number desntiy no+
!     pion3(levs)     number density o2+
!     r                
!     sigped           pedersen conductivity
!     sighall          hall conductivity
!     jphi(levs)     eastward curreil
!     jth(levs)      southward
!     rvin(levs)     ion/neutral collision  frequency param    
!     ramin(levs)    mean ion mass
!     a5               meridional ion drag term
!     b5               zonal ion drag term
!     c7               joule heating term
!     
      real(kind=kind_phys)      :: teff(levs), pion1(levs), pion2(levs)
      real(kind=kind_phys)      :: sigped, sighal, pion3(levs)
      real(kind=kind_phys)      :: rvin(levs), ramin(levs)
      real(kind=kind_phys)      :: r, brad, bth, dip
      real(kind=kind_phys)      :: a5, b5, c7
      real(kind=kind_phys)      :: jth,jrad  
      real(kind=kind_phys)      :: jphi,   gw 
      real(kind=kind_phys)      :: eden(im,levs)  !electron density
      real(kind=kind_phys)      :: eden_chiu(im,levs)  !electron density from chiu
      real(kind=kind_phys)      :: eden_aurora(im,levs)  !electron density from tiros
!
      real(kind=kind_phys),dimension(levs) :: pr1,rho1,grav1, adt1
      real(kind=kind_phys),dimension(levs) :: o3p_1,o2_1,n2_1,aur_1
!   
      real(kind=kind_phys),dimension(im, levs) :: el_diag, shal_diag, sped_diag  
!===================================================================
!         calculate electric field and magnetic field         
!===================================================================
! vay-2016: (im, ix) => (im)
!
      call wamphys_geteb(im, dayno,utsec,jdat, f107,f107d,kp,maglat,maglon,   &       
          essa,swbt,swang,swvel,swbz,swden,ee1,ee2)
	  
!     print *, maxval(ee1), maxval(ee2)
!	pause ' wamphys_geteb '	  


! ===================================================================
! =                   calculate electron density                    =
! ===================================================================
! 
! set the power index 1 to 10 eventually read in
! set the number of gigawatts of auroral power (used if activity level = 10)
!
      if (trim(spw_drivers)=='swpc_fst') then
         tiros_activity_level = hpi
         gw = hp
      else          
         tiros_activity_level = tiros_activity_fix
         gw = gw_fix
      endif	 
 
      eden(:, :)=0.      ! or min-value
      do i = 1,im 

         ii = i
         do k=1,levs
           ht1(k)=ht(i,k)
         enddo

         call wam_earth_chiu_model(sda,sza(i),maglat(i),   &              
              maglon(i),rlt(i), rlat(i), f107,             &              
              dipang(i)*deg_to_rad, dayno, ht1, eden_chiu, &            
              ii,lev1, levs,im)

! tiros_ionize returns ionization rates for o, o2, and n2 for a given
! geomagnetic latitude gl and magnetic local time mlt based on
! tiros/noaa statistical maps of energy influx and characteristic energy
! the ionization rates assumes an observed spectrum based on the
! characteristic energy ch at the location, the tiros spectrum is
! between 300ev - 100kev, below 300ev assumes a maxwellian where
! the average energy of the distribution is ch
! the model atmosphere should be provided by calling program (e.g., ipe)
! a sample model profile is read in from ipe_nh1_data and ipe_nh2_data
!
! input
! geomagnetic latitude mlat in radians
! magnetic hour angle from noon essa in degrees
!      essa = essa_in
! output eden_aurora units number/m3   ion density from aurora
  
!
! netcdf input files
!
         if(tiros_switch == 0) then	     
           do k=1, levs
              pr1(k) =pres(i,k)
              rho1(k)=rho(i,k)
              grav1(k) = grav(i,k)
              adt1(k) = adt(i,k)
              o3p_1(k) =o_n(i,k)
              o2_1(k) =o2_n(i,k)
              n2_1(k)=n2_n(i,k)
           enddo
   
           call wam_tiros_ionize(lev1,levs, pr1, rho1, ht1,   &
             grav1,o3p_1,o2_1, n2_1,adt1,maglat(i),essa(i),   &
             tiros_activity_level, gw, aur_1)
       
           do k=1, levs
              eden_aurora(i,k) = aur_1(k)
           enddo
         endif  		
   
   ! ascii input files
   
         if (tiros_switch ==1) then
         
            call wam_tiros_ionize_data(pres(i,:), lev1,levs,ht1,emaps,cmaps,  &
                 djspectra, grav(i,:),o_n(i,:),o2_n(i,:),                     & 
                 n2_n(i,:),adt(i,:),maglat(i),essa(i),tiros_activity_level,   &
                 gw, eden_aurora(i,:) )                                          !don't use, eflux, ch)
   	   
         endif

         eden(i,:)=sqrt(eden_chiu(i,:)**2   + eden_aurora(i,:)**2)
  
      enddo
      
!     print *, maxval(eden_aurora), minval(eden_aurora)
!	pause ' eden_aurora wamphys '

      el_diag = eden
       
!=================================================================
!         calculate ion drag, joule heating and particle        =
!                     precipitation terms                       =
!=================================================================
      do i = 1, im 
!======================================horizontal dep-nt loops
  
          dip     = dipang(i)* deg_to_rad
          sdip    =sin(dip)                              !new
          cdip    =cos(dip)                              !new

          elecx   =ee2(i)*sdip
          elecy   =ee1(i)

          ssa=rlon(i)+(utsec/3600.-12.)*pid12              ! radians  ????? for south pole
          dif     =essa(i)*deg_to_rad -   ssa                     !check unit
          cosdif  =cos(dif)
          sindif  =sin(dif)

          elz(i)  =-1.*ee2(i)*cdip
          if(sdip.ge.0.) then
            elx(i)  =elecx*cosdif-elecy*sindif
            ely(i)  =elecx*sindif+elecy*cosdif
          else
            elx(i)  =elecx*cosdif+elecy*sindif
            ely(i)  =-1.*elecx*sindif+elecy*cosdif
          endif

          btheta  =-1.*btot(i)*cdip*cosdif               ! new
          bphi    =-1.*btot(i)*cdip*sindif               ! new
          bth     = btot(i)                              ! in teslas, so no *1.e-9
          brad    = -1.*bth*sdip
!====================================== vertical dep-nt loops
          do k = lev1,levs
            nden(k)=o_n(i,k)+o2_n(i,k)+n2_n(i,k)
            teff(k) = adt(i,k)
  	  
            v1(k)=           -adv(i,k)      ! v1 positive south
  	  
            on(k)=o_n(i,k)
            o2n(k)=o2_n(i,k)
            n2n(k)=n2_n(i,k)
          enddo
	
          do k=lev1,levs    
            pion1(k) = on(k)
            pion2(k) = 0.5*(o2n(k)+n2n(k)) 
            pion3(k) = 0.5*(o2n(k)+n2n(k))
          enddo
! get ion neutral collision frequency

         call wam_ionneut(on,o2n,n2n,pion1,pion2,pion3,      &
                     teff,rvin,ramin,levs,lev1)
		  
! calculate ion drag and electron deposition
! jth                     - n/s electrical conductivity
! jphi                    - e/w electrical conductivity

        do k=lev1,levs 
          r       = (ramin(k)*rvin(k))/(elch*bth)
          sigped  = (eden(i,k)*elch*r)/(bth*(1.0+r**2))
          sighal  = sigped*r
	  
          sped_diag(i,k) = sigped
          shal_diag(i,k) = sighal
	  
          jphi = sigped*(ely(i)-v1(k)*brad)              &               
              - sighal*(elx(i) + adu(i,k)*brad)*sdip         !new fix of 201612=> /sdip => *sdip
	      
          jth  = sigped*(elx(i) + adu(i,k)*brad) +       &               
              sighal*(ely(i)-v1(k)*brad)*sdip                !new
	      
          jrad = sigped*(elz(i) - adu(i,k)*btheta+v1(k)*bphi)  &         
             -sighal*(ely(i)-v1(k)*brad)*cdip                !new
	     
          a5 =    jphi*brad-jrad*bphi                         !new w/o rho-factors
          b5 = -(jth*brad-jrad*btheta)                        !new
	     
          c7 =jth*(elx(i)+adu(i,k)*brad)+        &             
              jphi*(ely(i)-v1(k)*brad)+jrad*elz(i)            !new
! calculation of ion drag terms end
          dvdt(i,k)   = -a5
          dudt(i,k)   =  b5
          jh(i,k)     = c7
        enddo
        do k=1,lev1-1 
          dvdt(i,k) = 0.
          dudt(i,k) = 0.
          jh(i,k) = 0.
        enddo
      enddo                          ! i-hor.-pixel loop

!        if (mpi_id == 0) then
!           print *, 'vay-ion getionparams'
!        endif

      return
  end subroutine wgetionparams

  subroutine wam_ionneut(p1,p2,p3, pion1,pion2,pion3, t, vin, a_min, nmax, n0)
         
      implicit none 
      
      integer, intent(in)                 :: nmax, n0   !levs-150, levs1-90
      real(kind=kind_phys), intent(in)    ::  p1(nmax) , p2(nmax) , p3(nmax) 
      real(kind=kind_phys), intent(in)    ::   t(nmax) 
      real(kind=kind_phys), intent(in)    :: pion1(nmax) , pion2(nmax), pion3(nmax)      
      
      real(kind=kind_phys), intent(out)    ::         vin(nmax) ,  a_min(nmax)
 
!                                                               
      real(kind=kind_phys)    :: sum_vay ,  v1 , v2 
      integer :: n
      real(kind=kind_phys)    ::     a(3) , b(3) 
      
      real(kind=kind_phys), parameter :: amu=1.66e-27  
      real(kind=kind_phys), parameter :: factor=1.0
      real(kind=kind_phys), parameter :: ftiny=1.e-36
      real(kind=kind_phys)    :: mi1 , mi2, mi3, summol
      real(kind=kind_phys)    ::  logtn    
      data mi1 , mi2, mi3/16. , 30., 32./
!********************************************************************
! the following a,b, are cooeficients used to caculate ion-neutral
! collision frequency. tim's thesis  3.5a,3.5b.  mjh 1.9.97
!********************************************************************
      
      data a/3.42e-11 , 6.66e-10 , 6.82e-10/
      data b/2.44e-10 , 4.28e-10 , 4.34e-10/
      
      do  n = n0 , nmax
         summol = pion2(n) + pion3(n)
         if(summol.lt.ftiny) summol=0.0
         sum_vay =pion2(n) + pion1(n)  + pion3(n)
         v2 = b(1)*p1(n) + b(2)*p2(n) + b(3)*p3(n)
         logtn =log10(t(n))
         v1 = a(3)*p3(n) + a(2)*p2(n) + a(1)*p1(n)*factor*sqrt(t(n))    &
            *(1.08 - 0.139*logtn+4.51e-03*logtn*logtn)
         if(v1.lt.ftiny) v1=0.0    ! v1 = max(v1, ftiny)
         if(v2.lt.ftiny) v2=0.0    ! v2 = max(v2, ftiny)

         vin(n) = (v1*pion1(n)+v2*summol)*1.e-06/sum_vay
         a_min(n) = (pion1(n)*mi1+pion2(n)*mi2+pion3(n)*mi3)*amu/sum_vay
	 
      enddo

      return
  end subroutine wam_ionneut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! tiros empirical model
! idea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wamphys_geteb(im, dayno,utsec,jdat, f107,f107a,kp,maglat,maglon,   & 
        essa,swbt,swang,swvel,swbz,swden,ee1,ee2)
	
      use  wamphys_efield              !only  : wamget_efield    
!     use efield_wam                  !  iday,iyear,iday_m,imo,f107d,by,bz,ut,v_sw     
!     use date_def                   ! needs extra work to pass idate
      use wam_efield_setdef_data
      implicit none
      
      integer, intent(in) :: im     ! number of data points in efield 
      integer, intent(in) :: dayno  !  day of year
      integer, intent(in) :: jdat(8)
              
      real(kind=kind_phys), intent(in)    :: utsec  ! second
      real(kind=kind_phys), intent(in)    :: f107, f107a  ! 
      real(kind=kind_phys), intent(in)    :: kp  ! 
      real(kind=kind_phys), intent(in)    :: maglat(im)  ! magnetic latitude (rad)
      real(kind=kind_phys), intent(in)    :: maglon(im)  ! magnetic longitude (rad)
      real(kind=kind_phys), intent(in)    :: essa(im)    ! degree
      real(kind=kind_phys), intent(in)    :: swbt, swang, swvel, swbz, swden
      real(kind=kind_phys), intent(out)   :: ee1(im)     ! electric field x direction mv/m
      real(kind=kind_phys), intent(out)   :: ee2(im)     ! electric field y direction mv/m

      integer ::   i,k,iref,jref
!     real(kind=kind_phys)  :: bz, by
      real(kind=kind_phys)  :: utsec_last
      real(kind=kind_phys)  :: dx,dy,aa,bb,maglond,maglatd               
      real(kind=kind_phys)  :: ed11(0:nmlon,0:nmlat),ed22(0:nmlon,0:nmlat)
      real(kind=kind_phys)  :: v_sw

      data utsec_last/-1./
      save utsec_last,ed11,ed22

! initiate
! calculate efield only if diff time step

      if(utsec.ne.utsec_last) then
        utsec_last=utsec	
!     use f107d directly from observation sheet, revised by tzu-wei and zhuxiao 
!     efield.f:   tilt = get_tilt( iyear, imo, iday_m, ut)
!     use efield         !  iday,iyear,iday_m,imo, ut, f107d, by, bz, v_sw 

 
        f107d = f107a         ! taken by efield by common-use
	
        ut = utsec/3600.	! taken by efield by common-use
        iday = dayno          ! day of year
        imo  =  jdat(2)
        iday_m =  jdat(3)
	  iyear =  jdat(1)
	
        bz = .433726 - kp*(.0849999*kp + .0810363)    &                  
             + f107a*(.00793738 - .00219316*kp)
        by = 0.
        v_sw = 450.00
!================================================   call from "module efield"
        call  wamget_efield(swbt, swang, swvel, swbz, swden)

!       print*,'www'
!       print'(8f10.4)',potent(0:180,68)
        ed11=ed1
        ed22=ed2
!       print*,'ed2',ed2(149,65)
      endif
!
      do k=1,im
        maglatd=maglat(k)*rad_to_deg

        jref=0
        dy=0.0
	
!efield.f:     &         ylonm, ylatm   ! magnetic longitudes/latitudes (degc .....ylatm(0:nmlat)ylonm(0:nmlon)

        do i=0,nmlat-1
!         if(maglatd.ge.ylatm1(i)-90..and.maglatd.le.ylatm1(i+1)-90.)   
          if(maglatd.ge.ylatm (i)-90..and.maglatd.le.ylatm (i+1)-90.) then   
            jref=i
!hmhj       dy=(maglatd-ylatm1(i)+90.)/(ylatm1(i+1)-ylatm1(i))
            dy=(maglatd-ylatm (i)+90.)/(ylatm (i+1)-ylatm (i))
          endif
        enddo 

!       maglond=maglon(k)/pi*180.
        maglond=essa(k)+180.
        if(maglond.lt.0.)   maglond=maglond+360.
        if(maglond.gt.360.) maglond=maglond-360.
! ylonm = [0,360]
        iref=0
        dx=0.0
        do i=0,nmlon-1
          if(maglond.ge.ylonm (i).and.maglond.le.ylonm (i+1)) then
            iref=i
            dx=(maglond-ylonm (i))/(ylonm (i+1)-ylonm (i))
          endif
        enddo 

        aa=(1.-dx)*ed11(iref,jref)+dx*ed11(iref+1,jref)
        bb=(1.-dx)*ed11(iref,jref+1)+dx*ed11(iref+1,jref+1)
        ee1(k)=(1.-dy)*aa+dy*bb
        aa=(1.-dx)*ed22(iref,jref)+dx*ed22(iref+1,jref)
        bb=(1.-dx)*ed22(iref,jref+1)+dx*ed22(iref+1,jref+1)
        ee2(k)=(1.-dy)*aa+dy*bb

      enddo
      return
  end subroutine wamphys_geteb
!idea_geteb
end module wamphys_ion
