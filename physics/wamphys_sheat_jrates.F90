!
! ccpp-version of UV/EUV solar radiation/jrates in column physics
!
! 2017/1018: VAY jo3_rate is also caluclated for ox-chemistry       o3dissociation_rate)
!
! vay-2015   solar_heat_dissociation_tiegcm   1d_height subroutine
! 2022/06   S.I. Karol ccpp-updates of WAMPHYS vay-2015
!
!	   
   subroutine wamphys_sheat_jrates(np,nps,o,o2,o3,n2, ho, ho2, hn2,         &
         f107, f107d, cospass, dayno, height, sheat,sh1, shsrc,shsrb,shlya, &
	 jo2_rate, jo3_rate)
       
       use machine,                      only : kind_phys 
       use  wamphys_const,               only : pi, pi2, pid2, r0 => rearth       
      
       use wamphys_set_data_solar,       only  : euv37, nsp_euv            ! euv-flux(nsp_euv=37)
       use wamphys_init_module   ,       only  : o2_scale_factor,  srbeff  ! (levs)-arrays
       use wamphys_init_module,          only  : effuv, effeuv             ! ctip-orig
       use wamphys_init_module,          only  : eff_src, eff_srb, eff_lya ! (levs)-arrays
       use wamphys_set_data_solar,       only  : rwpcc                     ! pcc/wavelength
       use wamphys_set_data_solar,       only  : csao, csao2, csan2, csao3
       use wamphys_set_data_solar,       only  : csio, csio2, csin2
       use wamphys_set_data_solar,       only  : csdo2, csdeo2
       use wamphys_set_data_solar,       only  : dsfmin, dafac              !sfmin, afac
       use wamphys_set_data_solar,       only  : nwaves, nwaves_euv
       use wamphys_set_data_solar,       only  : lyman_a_num,nwaves_src
       
!-------------------------------------------------------------------------
! calculate solar heating and dissociation rates
! 
!  **
!  calculates solar heating and dissociation rates 
!     from euv and uv(src+lya+srb)
!  
!  input:
!  o atomic oxygen number density profile m-3
!  o2 molecular oxygen number density profile m-3
!  n2 molecular nitrogen number density profile m-3
!  ho atomic oxygen scale height profile m
!  ho2 molecular oxygen scale height profile m
!  hn2 molecular nitrogen scale height profile m
!  f107 solar flux 
!  f107d 81 day mean f10.7
!  cospass cosine of solar zenith angle
!  dayno day of year ( 1 to 366)
!
!  output:
!  sheat, sh1(euv), sh2(uv) heating rates j/m3
!  o2dissociation_rate dissociation rates s-1
!---------------------------------------------------------------------------
      implicit none
      
!      integer, parameter  :: idea_solar_fix=1 
!       logical, parameter  :: para_schemes=.false.
!
! input
!
      integer, intent(in) :: np                                        ! number of pressure levels
      integer, intent(in) :: nps                                       ! pressure index to start 2Pa/~75.7 km
                                                                       ! now from the merging domain 8.5*7 =59.5km/52.5km
								       !
      integer, intent(in) :: dayno                                     ! day of year 
           
      real(kind=kind_phys), intent(in)    :: o(np),o2(np),n2(np)       ! number density 1/m3
      real(kind=kind_phys), intent(in)    :: o3(np)                    ! ozone density  1/m3
      real(kind=kind_phys), intent(in)    :: ho(np),ho2(np),hn2(np)    ! scale height(m)

      real(kind=kind_phys), intent(in)    :: f107       ! f10.7cm
      real(kind=kind_phys), intent(in)    :: f107d      ! 81 day mean f10.7
      real(kind=kind_phys), intent(in)    :: cospass    ! cos zenith angle
      real(kind=kind_phys), intent(in)    :: height(np) !layer height (m)

! output

      real(kind=kind_phys), intent(out)   :: sheat(np),sh1(np)  ! w/m3 heating rate
      real(kind=kind_phys), intent(out)   :: shsrc(np), shsrb(np), shlya(np)   
       
      real(kind=kind_phys), dimension(np) :: sh2

      real(kind=kind_phys), intent(out)   :: jo2_rate(np)               ! down to xb =7.5 km 
      real(kind=kind_phys), intent(out)   :: jo3_rate(np)               ! diag-now down to xb =7.5 km 

! locals 

      real(kind=kind_phys) :: spaeuv              ! vay scalar instead of 2d-array paeuv
      real(kind=kind_phys) :: z, coschi
      real(kind=kind_phys) :: seco,seco2, seco3, secn2
      real(kind=kind_phys) :: wo,wo2, wo3, wn2
      real(kind=kind_phys) :: tau,tauo,tauo2,tauo3,taun2 
      real(kind=kind_phys) :: nightfac
      real(kind=kind_phys) :: attenuation, local_flux
      real(kind=kind_phys) :: rnight_o,rnight_o2,rnight_o3,rnight_n2

! solar variability factor for srb dissociation calculation
      real(kind=kind_phys) :: flux(nwaves)
      real(kind=kind_phys) :: fmxfmn
      real(kind=kind_phys) :: pind
      real(kind=kind_phys), dimension(np) :: sco2, sco, scn2,sco3 ! slant columns o2/o/n2
      
! local  o2 dissociation rates

      real(kind=kind_phys) :: pre_loc, paeuv_loc,  jsrc_loc                   
      real(kind=kind_phys) :: jlya_loc, jsrb_loc,  jeuv_loc
	  
! o3- local
      real(kind=kind_phys) :: jo3_euvloc,  jo3_srcloc, jo3_lyaloc, pdnolya
      real(kind=kind_phys) :: jo3_harloc, jo3_chhloc
      
! lyman-a and src, srb and euv o2 dissociation rate coefficients

      real(kind=kind_phys) :: jlya(np), jsrc(np), jsrb(np),jeuv(np)
! srb
      real(kind=kind_phys) :: srband, srband_fac,srbheat(np)
      


! vayudin-2015
      real(kind=kind_phys) :: vay1, vay2, vay3, vay4, rnight
!   
      real(kind=kind_phys) :: vay_fmxfmn, vay_srband_fac, vay_srb3
      real(kind=kind_phys) :: sfeps
      real(kind=kind_phys) :: temp1, sco3t, smfac
      
!  


      integer i,j,jinv

! raa: solar 10-cm flux normalized to 1 au

      real(kind=kind_phys) :: f107_1au, f107d_1au

       nightfac=1.e-9                     ! ratio of nightime ionisation to sec=1.
       smfac = 1.e-4                      !conversion

! raa: reintroduce the sed factor, normalize solar flux to 1 au
        sfeps   = (1.0+0.0167*cos(pi2*(dayno-3.)/366.0))**2
        f107_1au =  f107/sfeps
        f107d_1au = f107d/sfeps

! vay-2017 fmxfmn=(f107-71.)/(220-71.) isn't right for f107~70.

        fmxfmn=(f107_1au-65.)/165.   
        srband_fac=0.784591675
        vay_fmxfmn=(1.+0.11*fmxfmn)
        vay_srband_fac=(1.+srband_fac)*.1

! raa: replace eccentric with sfeps
        vay_srb3=vay_srband_fac * vay_fmxfmn * sfeps
        pind=0.5*(f107_1au+f107d_1au) -80.
!
      do j = 1, nwaves  
        jinv = nwaves+1-j
        flux(j)=dsfmin(jinv)*max(0.8, (1.0+dafac(jinv)*pind) )*sfeps
        if (flux(j) .lt. 0.0) flux(j) = 0.0	 
      enddo
!      
          sheat(1:np)=0.
          sh1(1:np)=0.
          sh2(1:np)=0.
          shsrc(:) =0.
          shsrb(:) =0.
          shlya(:) =0.
          jo2_rate(1:np)=0. 
          jo3_rate(1:np)=0. 
!
! compute slant columns in 1/cm2
!
!      call wam_slantcolumns(np, nps, r0, height, cospass,    &
!          ho, ho2, hn2, o, o2, o3, n2, sco, sco2, sco3, scn2, nightfac, smfac)
!
! add extra: (seco, rnight_o)-arrays to avoid 2-nd calls of sub_chapman
! below idea-style for slant columns : wo2=o2(i)*ho2(i)*seco2*1.e-4
!   
      do i=nps,np
!                                   set total height in m rad+zmeters
       z = r0 + height(i)

       jlya(i)   = 0.
       jsrc(i)   = 0.
       jsrb(i)   = 0.
       jeuv(i)   = 0.
       srband    = 0.
       srbheat(i)= 0.

       pre_loc     = 0.
       jsrc_loc    = 0.
       jsrb_loc    = 0.
       jlya_loc    = 0.
       jeuv_loc    = 0.
       jo3_euvloc    = 0.
       jo3_lyaloc    = 0.
       jo3_srcloc    = 0.
       tau=0.
!  
! calculate sec(za), incorporating chapmann grazing incidence function
      call wamsub_chapman(cospass, ho(i), z, nightfac, seco, rnight_o) 
      call wamsub_chapman(cospass, ho2(i),z, nightfac, seco2,rnight_o2)
      call wamsub_chapman(cospass, ho(i)*.333,z,nightfac,seco3,rnight_o3)
      call wamsub_chapman(cospass, hn2(i),z, nightfac, secn2,rnight_n2) 

! for tau convert column abundance from m^-2 to cm^-2 (10^-4)

        wo= o(i) *ho(i) *seco*  smfac
        wo2=o2(i)*ho2(i)*seco2* smfac
        wn2=n2(i)*hn2(i)*secn2* smfac
	
        wo3 =sco3(i)

!  loop over all 37 bands
!       sh1(i)=0.
!       sh2(i)=0.
!       jo2_rate(i)=0.
!       jo3_rate(i)=0.

      sh1(i)=0.
      do j = 1, nwaves
         tauo =csao(j)*wo
         tauo2=csao2(j)*wo2
         taun2=csan2(j)*wn2
         tauo3=csao3(j)*wo3
         tau=tauo+tauo2+taun2

! use parameters: solar_tuny = 1.e-36, min_flux =1.e-20, tau_max=200.
      if (tau .ge.  1.e-36 .and. tau <= 200.) then 
           attenuation = exp(-tau)
       else  if (tau < 1.e-36 ) then  
            attenuation = 1.0
       else  
            attenuation=0.
      endif

! raa: remove eccentric, flux is already normalized to sed
          local_flux=flux(j)*attenuation

          if(local_flux < 1.e-20 ) local_flux=0.

! euv heating calculation(w/m3)
     if(j <= nwaves_euv) then
 
         spaeuv=local_flux*rwpcc(j)*(csao(j)*o(i)*rnight_o +   &            
            csao2(j)*o2(i)*rnight_o2+csan2(j)*n2(i)*rnight_n2)
	 
         sh1(i)    =sh1(i)  + spaeuv

! euv jo2
          jeuv_loc = jeuv_loc+local_flux*(csdo2(j)+csdeo2(j))

! euv jo3
          jo3_euvloc = jo3_euvloc +local_flux*csao3(j)
       else

!  uv src/lya channels 
         spaeuv=local_flux*csao2(j)*o2(i)*rwpcc(j)*rnight_o2

! calculate jsrc, excluding jlya
           if(j /= lyman_a_num) then 
              jsrc_loc=jsrc_loc+local_flux*csao2(j)
              shsrc(i)  =shsrc(i)  +spaeuv
              jo3_srcloc = jo3_srcloc +local_flux*csao3(j)
           else                                                

! calculate jlya 
              jlya_loc=local_flux*csao2(j)
              shlya(i)  =  spaeuv                                        !* effuv(i) 
              pdnolya = (0.68431  *exp(-8.22114e-21*sco2(i))+    &
                  0.229841 *exp(-1.77556e-20*sco2(i))+           &
                  0.0865412*exp(-8.22112e-21*sco2(i)))*          &
                  flux(lyman_a_num)
              jo3_lyaloc =2.27e-17*pdnolya
           endif                                               ! src/lya
      endif         ! uv/euv
    enddo           ! end of wavelength loop - j-index

!==========================================================
!  qtotal(k,i,lat) = qtotal(k,i,lat)+ho2src(k,i)+ho2srb(k,i) 
!  calculate o2 schumunn runge band heating  see "o2srbc.f"
!==========================================================
      if (wo2 < 1.e18) then
         srband=o2(i)*2.43e-25               !!1.e-19*1.e-6  strobel-78 exp(21) 20% accuracy
      else
         srband=o2(i)*1.e-6/(0.66*wo2+3.44e9*sqrt(wo2))
      endif
      
      srband =srband * vay_srb3                ! (vay_srband_fac*vay_fmxfmn*eccentric)
      if(srband.lt.0.)   srband  = 0.

!      sum srb+src+lya heating into the total uv heat
!=====================================================
!      apply heat-efficiency factors & rnight_o2

      shsrb(i) = srband  *rnight_o2
      sh1(i)   = sh1(i)  *effeuv(i)                
      shsrc(i) = shsrc(i)*effuv(i)                 
      shlya(i) = shlya(i)*effuv(i)                
      sh2(i)   = shsrc(i) + shlya(i) +shsrb(i)     
! total
      sheat(i) = sh1(i) +sh2(i) ! euv + uv

!===============================================
! calculate jo2 due to srb
      if (wo2.lt.1.e19) then
! raa: add sed factor
        jsrb_loc = vay_fmxfmn*1.1e-7*exp(-1.97e-10*(wo2**0.522))*sfeps         
      else
        jsrb_loc =vay_fmxfmn*1.45e8*(wo2**(-0.83))*sfeps
      endif

! sum dissociation rates, eccentric is inclued in the local_flux and srb
 
      jsrc(i)   = jsrc_loc
      jsrb(i)   = jsrb_loc          
      jlya(i)   = jlya_loc
      jeuv(i)   = jeuv_loc
!     
      jo2_rate(i) =o2_scale_factor(i)*rnight_o2*      &
          (jsrc_loc+jlya_loc+jsrb_loc+jeuv_loc)
!
! Jo3_rate: jo3_euvloc+jo3_srcloc+jo3_lyaloc+jo3_harloc+jo3_chhloc
! jo3_harloc & jo3_chhloc
!         
      sco3t = sco3(i)                        
      if (sco3t < 1.e+5) sco3t = 1.e+5
      temp1 = 1.0e-3*exp(-1.5577e-13*sco3t**0.6932)
      if (sco3t < 1.6e+20)                              &
          temp1 = 1.04e-2*exp(-1.0217e-6*sco3t**0.3587)
! hartley bands of o3:
!          jo3_harloc = sfeps*0.68*                          &
!           (temp1 + 1.085*exp(-1.4912e-14*sco2(i)**0.5298)*  &
!           ( 4.053e-4  *exp(-8.1381e-16*sco3t**0.8856) +     &
!             4.700e-6  *exp(-1.5871e-14*sco3t**0.7665)   )*  &
!            *exp(-1.4655e-25*sco2(i)**1.0743) )
      jo3_harloc =                                          &
          (temp1*0.68+exp(-1.4912e-14*sco2(i)**0.5298)*    &
          (4.053e-4  *exp(-8.1381e-16*sco3t**0.8856)+      & 
          4.7e-6    *exp(-1.5871e-14*sco3t**0.7665))*      &
          1.085     *exp(-1.4655e-25*sco2(i)**1.0743)*0.68)*sfeps   
!
! chappius and huggins bands:
!
      jo3_chhloc = sfeps*                                 &
          (    4.5e-4*exp(-3.4786e-11*sco2(i)**0.3366     &
          -1.0061e-20*sco3t  **0.9719)                    & 
          + ( 7.5e-4 *exp(-2.7663e-19*sco3t**1.0801 )     &
          +2.5e-4/(1.+1.5772e-18  *sco3t**0.9516))        &
          *exp(-1.0719e-10*sco2(i)**0.3172 )  )
!
        jo3_rate(i) = rnight_o3 * (jo3_euvloc+jo3_srcloc+jo3_lyaloc+jo3_harloc+jo3_chhloc)
      enddo    ! height loop
      
      return
!
! below nps of tiegcm (~76 km)
!
      do i=1,nps-1
         sco3t =  exp(-abs(float(i-nps))*0.2) 
!
! jo3 (o3p+o1d) below ~50 km is ~ constant with small decrease 
! to keep o2 -constant below 80 km
!
         jo2_rate(i) = jo2_rate(nps)*sco3t 
         jo3_rate(i) = jo3_rate(nps)*sco3t  
      enddo
!
      return
      end subroutine wamphys_sheat_jrates
!      
!
      subroutine wamphys_heat_uveuv(im, levs, te, cospass, o_n,o2_n,o3_n,n2_n, rho, cp, &
            dayno, prsl, zg, grav,am, maglat, f107, f107d, kpa, dtrad)
!	    
! diagnostics seuv_diag, shsrc_diag, shsrb_diag, shlya_diag, qno_diag, no_diag)	
!    
! S. I. Karol May 2022 ccpp-version of idea_solar_heating.f
!
! interface to compute heating rates after updates of O-O2 by chemistry and diffusion + NO-cooling
!	    
! 
	 
      use machine, only                : kind_phys
      use wamphys_init_module, only    : amo2, amn2, amo
      use wamphys_set_merge_rad, only  : npsrad 
      
      use physcons,  only              : rgas=>con_rgas, avgd => con_avgd
      
      implicit none   
      
      integer, intent(in) :: im       !number of data piont in te
 
      integer, intent(in) :: levs     !number of press level
      integer, intent(in) :: dayno    ! calender day
!
!VAY-2016
      real(kind=kind_phys), intent(in)    :: f107, f107d, kpa ! solar-geo drivers
!
      real(kind=kind_phys), intent(in)    :: te(im,levs)  !temperature
      real(kind=kind_phys), intent(in)    :: cospass(im)  ! cos zenith angle
      real(kind=kind_phys), intent(in)    :: maglat(im)   ! in radians for NO-snoe
      real(kind=kind_phys), intent(in)    :: cp(im,levs)  ! 
      real(kind=kind_phys), intent(in)    :: o_n(im,levs) !number density of O(/cm3) 
      real(kind=kind_phys), intent(in)    :: o2_n(im,levs)!number density of O2
      real(kind=kind_phys), intent(in)    :: o3_n(im,levs)!number density of O2
      real(kind=kind_phys), intent(in)    :: n2_n(im,levs)!number density of N2
      real(kind=kind_phys), intent(in)    :: am(im,levs)  !mass of mix (kg)
      real(kind=kind_phys), intent(in)    :: prsl(im,levs)!layer press (Pa)
      real(kind=kind_phys), intent(in)    :: zg(im,levs)  !layer height (m)
      real(kind=kind_phys), intent(in)    :: grav(im,levs)! (m/s2)
      real(kind=kind_phys), intent(in)    :: rho(im,levs)  ! density (kg/m3) 
       
      real(kind=kind_phys), dimension (im, levs), intent(out) :: dtrad   ! (K/s) solar heating QS-QNO
            
      real(kind=kind_phys)  :: jo2_rate(levs), jo3_rate(levs)     
! Locals

      integer  i,k
      real(kind=kind_phys) :: t(levs),n2(levs), o(levs),o2(levs),o3(levs)
      real(kind=kind_phys) :: ho(levs), ho2(levs),hn2(levs)
      real(kind=kind_phys) :: sheat(levs),qno(levs),no_new(levs)
      real(kind=kind_phys) :: amm(levs),prr(levs),alt(levs),nn(levs),sh1(levs),sh2(levs)
      real(kind=kind_phys), dimension(levs) :: shsrc,shsrb,shlya
      real(kind=kind_phys) :: ht(levs)
      
!out-diag      
!   

      real(kind=kind_phys), dimension (im, levs)  :: dtdt_euv, dtdt_uvr, dtdt_qno        ! (K/s)     

      real(kind=kind_phys) :: vay_rgas_o, vay_rgas_o2,vay_rgas_n2 

      real(kind=kind_phys) :: vay_avgd, tg_vay, vay_cprho, rcpro, ron2

!==================================================================

      vay_rgas_o =  1.e3*rgas/amo
      vay_rgas_o2 = 1.e3*rgas/amo2
      vay_rgas_n2 = 1.e3*rgas/amn2
      vay_avgd    = 1.e3*avgd 
      ron2        = amo/amn2          

      do i=1,im
      
        do k=1,levs
          o(k)=o_n(i,k)                            !/m3 in idea_phys
          o2(k)=o2_n(i,k)   
          o3(k)=o3_n(i,k)        
          n2(k)=n2_n(i,k)     
          t(k) = te(i,k)                 
          tg_vay = te(i,k)/grav(i,k)
          ho(k) =vay_rgas_o*tg_vay                  !m
          ho2(k)=.5*ho(k) 
          hn2(k)=ron2*ho(k)  
          ht(k) =  zg(i,k)
          alt(k) = 1.e-3*ht(k)
          prr(k)=prsl(i,k)
          nn(k)=o(k)+o2(k)+o3(k)+n2(k)
          amm(k)=am(i,k)*vay_avgd 
        enddo
       call wam_getno1d(levs,f107,kpa,maglat(i),dayno,alt,prr,nn,amm,no_new)	
!        no_diag(i,:) = no_new(:)  
! no-cooling
      call wamphys_coolno1(levs, npsrad, t, o, o2, no_new, qno) 
     
! get heating
!   
      call wamphys_sheat_jrates(levs,npsrad,o,o2,o3,n2, ho, ho2, hn2, &    
          f107,f107d, cospass(i), dayno, ht, sheat,sh1, shsrc,shsrb,shlya, jo2_rate,jo3_rate)


        do k=npsrad,levs
!    
          rcpro = 1./(cp(i,k)*rho(i,k))
!	  
! new diag-arrays	  
!
          dtdt_euv(i,k)=sh1(k)*rcpro
          dtdt_uvr(i,k)=(shsrc(k)+shsrb(k)+shlya(k))*rcpro
	  dtdt_qno(i,k) =qno(k)*rcpro
!net q
          dtrad(i,k)=(sheat(k)-qno(k))*rcpro
! 
        enddo ! k
	        
          do k=1, npsrad-1
	   dtdt_euv(i,k)= 0.
	   dtdt_uvr(i,k)= 0.
	   dtdt_qno(i,k) =0.
	   dtrad(i,k)=0.
	  enddo
        enddo !i 	  
	return      
      end  subroutine wamphys_heat_uveuv
      subroutine wamsub_chapman(coschi,scale_ht, ht, nfac, seco, rnight)
!=============================================================
! vay oct 26/2016 review and corrections for undefined "pi"
!     ... etc..function chapmann(coschi,scale_ht,ht)
!         function is correponds to sec(zenith angle) 1./coschi
! may 2022 updates for FV3WAM/ccpp    kind_phys
!=============================================================
      use machine,     only: kind_phys
      use wamphys_const, only : pi, pi2, pid2, r0 => rearth, r_2_d
      
      
      implicit none     ! adding this .... major bugs are on compilation
!
      real(kind=kind_phys), intent(in)  :: scale_ht          
      real(kind=kind_phys), intent(in)  :: ht
      real(kind=kind_phys), intent(in)  :: coschi
      real(kind=kind_phys), intent(in)  :: nfac    
      real(kind=kind_phys), intent(out) :: seco
      real(kind=kind_phys), intent(out) :: rnight
      
!      real, parameter :: r0    = 6.370e06    ! radius of earth (m)

      real(kind=kind_phys) y_err, rad_to_z, chi, chid, erfcL, schi

      chi  = acos(coschi)    
      SCHI = sqrt(1. - coschi*coschi)
      chid = chi*r_2_d
!
! daytime
!
      rnight=1.0 
      if(chid.le.75.) then 
          seco = 1./coschi
          return
       endif   
        
! calculate chapman bit for angles over 75-105 degrees  
! spherical lighting
!     
      if(chid.gt.75. .and. chid .lt. 105.) then

          rad_to_z = ht/scale_ht

! step1: calculate error function

          y_err = sqrt(0.5*rad_to_z) * abs(coschi)
          
!bug-diag  if (y_err.gt.100.0.and.ht.lt.(2*r0*1.e2))   
!       if (y_err.gt.100.0.and.ht.lt.(2*r0))            
!     &    write(6,*) 'warning problem in chapman function',  y_err,ht

        if (y_err.le.8.0) then
          erfcL = (1.0606963 + 0.55643831* y_err) /(1.0619898 + 1.7245609* y_err + y_err *y_err)
        else
            erfcL = 0.56498823 / (0.06651874 + y_err)
        endif
      
! step2. calculate chapmann
! for solar zenith angles 75 < theta <= 90
        
        if(chid.le.90.)then
            seco = sqrt(pid2 * rad_to_z) * erfcL
        else
! for solar zenith angles > 90 (equation 15)  pi2 or pid2
        seco = sqrt(pid2 * rad_to_z)*    &                  
             ( sqrt(schi)*exp(rad_to_z*(1-schi)) -0.5*erfcL)
!         write(6,*)'chapman over 90', chid,chapmann
         endif

      endif
!
! full-night
!       
       if(seco < 0.) then
       
          rnight=nfac        ! nightfac=1.e-9
          seco=1./0.07       ! should be big-number but in tw-code =1. deviated
                             ! exp(-1./0.07) =6.24875e-07
       endif


      end subroutine wamsub_chapman      
      
      subroutine wam_slantcolumns(levs, lev1, r0, zgi, cospass,   &
        ho, ho2, hn2, xo, xo2, xo3, xn2, sco, sco2, sco3, scn2, nightfac, smfac)
       
       use machine,     only: kind_phys
       implicit none
!
! compute slant columns for 3-major species
!

      real(kind=kind_phys),  intent(in)   :: nightfac, smfac 
      
      integer,intent(in) :: levs,lev1
      real(kind=kind_phys), intent(in)   :: r0, zgi(levs)
      real(kind=kind_phys), intent(in)   :: cospass
      real(kind=kind_phys), intent(in), dimension(levs)   :: ho, ho2, hn2            ! m
      real(kind=kind_phys), intent(in), dimension(levs)   :: xo, xo2, xo3, xn2       ! 1/m3
      real(kind=kind_phys), intent(out), dimension(levs)  :: sco, sco2, sco3,scn2    ! 1/cm2
! locals      
      real(kind=kind_phys) :: z ! meters
      
      real(kind=kind_phys) :: seco,  rnight_o
      real(kind=kind_phys) :: seco2, rnight_o2
      real(kind=kind_phys) :: seco3, rnight_o3, ho3      
      real(kind=kind_phys) :: secn2, rnight_n2
      real(kind=kind_phys), dimension(levs)  :: vco, vco2,vco3, vcn2               ! 1/cm2
      real(kind=kind_phys) :: dzm
      integer :: i,k,j
!
! from top 2 bottom 
!
      k=levs
      vco(levs)  = xo(levs)*smfac*ho(levs)
      vco2(levs) = xo2(levs)*smfac*ho2(levs)
      vcn2(levs) = xn2(levs)*smfac*hn2(levs) 
      ho3 =ho(levs)*.3333
      vco3(levs) = xo3(levs)*smfac*ho3
         z = r0 +zgi(k)
        call wamsub_chapman(cospass, ho(k), z,  nightfac,  seco, rnight_o ) 
        call wamsub_chapman(cospass, ho2(k), z, nightfac, seco2, rnight_o2)
        call wamsub_chapman(cospass, ho3,    z, nightfac, seco3, rnight_o3) 
        call wamsub_chapman(cospass, hn2(k), z, nightfac, secn2, rnight_n2) 


        sco(k)= vco(k)*seco*smfac
        sco2(k)= vco2(k)*seco2*smfac
        sco3(k)= vco3(k)*seco3*smfac
        scn2(k)= vco2(k)*seco2*smfac

      do k=levs-1, lev1, -1
      dzm = .5*(zgi(k+1)-zgi(k))
      vco(k)  = vco(k+1)  + (xo(k+1)+xo(k))*dzm 
      vco2(k) = vco2(k+1) + (xo2(k+1)+xo2(k))*dzm
      vco3(k) = vco3(k+1) + (xo3(k+1)+xo3(k))*dzm      
      vcn2(k) = vcn2(k+1) + (xn2(k+1)+xn2(k))*dzm 
        
         z = r0 +zgi(k)
         ho3 =ho(k)*0.3333
         call wamsub_chapman(cospass, ho(k), z,  nightfac,  seco, rnight_o ) 
         call wamsub_chapman(cospass, ho2(k), z, nightfac, seco2, rnight_o2)
         call wamsub_chapman(cospass, ho3,    z, nightfac, seco3, rnight_o3)     
         call wamsub_chapman(cospass, hn2(k), z, nightfac, secn2, rnight_n2) 
         
!
!  transform vcol  => scol
! 

        sco(k)=  vco(k) *seco*smfac
        sco2(k)= vco2(k)*seco2*smfac
        sco3(k)= vco3(k)*seco3*smfac
        scn2(k)= vcn2(k)*secn2*smfac

!
!        wn2=n2(i)*hn2(i)*secn2*1.e-4
!        n2(i)*hn2(i) = vcn2(k)
!
      enddo
      end subroutine wam_slantcolumns        
!      
      subroutine wam_getno1d(levs,f107in,kpain,mlatrad,doy,alt,pr, nair, am, no)
!      
! S. I. Karol May 2022 ccpp-version of idea_getno_snoe.f
!
! Oct 2016 VAY: new interface  getno1d ... mozaic of bugs in 2012-14 versions of WAM
!               (a) rad/deg (b) pr(cb) but expected in Pa (c) no-NO for lat > 80deg
!               (d) dangerous interp-n and DATA for eofs etc...tropical NO ????
!
! Mar 2018 Zhuxiao Li and Tzu-Wei Fang,read in the 24hr ave kp (kpa)
!       from driving parameters file instead of reading kp.  ??? REASONS
!
       use wamphys_math_interp, only : interpol_wamz
       use wamphys_init_module, only : eof => no_eof, nom => no_m, z16 => no_zkm, lat33 => no_mlat
       use wamphys_init_module, only : no_ny33, no_nz16, amno       
       use  wamphys_const,      only : pi, pi2, con_nzero, d2r => dtr, r2d => r_2_d
       use machine,           only : kind_phys         
!
      implicit none
! input
      integer, intent(in)  :: levs         !number of pressure level
      integer, intent(in)  :: doy          ! day of year from 1 to 365   

      real(kind_phys),    intent(in)  :: f107in       !f10.7 index instant 
      real(kind_phys),    intent(in)  :: kpain        ! 24hr average kp index
      real(kind_phys),    intent(in)  :: mlatrad      ! magnetic latitude in radians

      real(kind_phys),    intent(in)  :: alt(levs)    ! in km wam
      real(kind_phys),    intent(in)  :: pr(levs)     ! in pa
      real(kind_phys),    intent(in)  :: am(levs)     ! avg mass g/mol
      real(kind_phys),    intent(in)  :: nair(levs)      !/m3 number density
! out
      real(kind_phys),    intent(out) :: no(levs)     ! number density of no (/m3) 
!
! locals
!  we have         33, 16,7 now...
!      real(kind_phys) :: eof(33,16,3),nom(33,16),z16(16), lat33(33)

      real(kind_phys) ::  mlat          ! degrees  in snoe data
      real(kind_phys) ::  zm(no_nz16)
      real(kind_phys) ::  dz(levs)
      real(kind_phys) ::  f107, kpa
      real(kind_phys) ::  dx, dl,m1,m2,m3,theta0,dec
!

      integer :: iref,kref(levs)
      integer :: i,k,il,k1,k2
      integer :: kup, kdw      ! [NO]-domain in WAM-layers
!
       kpa  = kpain
       f107 = f107in

       mlat = mlatrad *r2d

!       print *, kpa, f107, ' getno1d-vay-mlat-deg ', mlat

       if(kpa .lt. 0.7)   kpa=0.7
       if(f107.lt.70.0)   f107 = 70.0
!
! find interp latitude lat33 in degrees...
!     
! vay-2016: mlat-bug
        do i=1,no_ny33-1
          if(mlat.gt.lat33(i).and.mlat.lt.lat33(i+1)) then
            iref=i
            dl=(mlat-lat33(i))/(lat33(i+1)-lat33(i))
            exit
          endif
        enddo
!vay oct-2016 add abs(mlat) > 80.
         if (mlat.le.lat33(1)) then 
            dl = 0.
            iref =1
         endif
         if (mlat.ge.lat33(no_ny33)) then 
            dl = 1.
            iref =no_ny33-1
         endif
!
!  snoe no interpolated to model grid (molecules/cm^3)
!... eof1 - kpa
!     m1 =  kpa * 0.689254 - 1.53366
      m1 =  kpa * 0.785760 - 1.94262           !waccmx-2015 eof1
      
!... eof2 - declination
      theta0 = pi2/365.*float(doy - 1)
      dec = 0.006918 - 0.399912 * cos(theta0)   + 0.070257 * sin(theta0)   &      
          - 0.006758 * cos(2.*theta0) + 0.000907 * sin(2.*theta0)          &
          - 0.002697 * cos(3.*theta0) + 0.001480 * sin(3.*theta0)
      dec = dec * r2d  
      
      m2 = -.319782 + dec*(.0973109 + dec*(.00048981 - dec*.000103608))
      m3 =  alog10(f107) * 6.44069 - 13.9832                 !waccmx-2015 eof3  
!   
!
!... zonal mean distrib. is sum of mean and eofs
! (dl)*no(i+1) + (1-dl)*no(i) interpolation to wam "mlat"
!  wx-code     do k = 1,nlev
!  zm(:,k) = no_mean(:,k) - m1 * eofs(:,k,1) + m2 * eofs(:,k,2) - m3 * eofs(:,k,3) 
!              end do
    
   do k=1,no_nz16
     zm(k) =  &
       dl*(nom(iref+1,k)- m1*eof(iref+1,k,1)+ m2*eof(iref+1,k,2)- m3*eof(iref+1,k,3))+ &                              (1.-dl)*(nom(iref,k)- m1*eof(iref,k,1)+ m2*eof(iref,k,2)- m3*eof(iref,k,3))
       
     if (zm(k).le. 0.0) zm(k)=con_nzero
    enddo
    zm = zm*1.e6          ! zm transform fom cm-3 to m-3  ok... due to data units
!
! vertical interp, from k1 to k2-1, extend k2 to levs, keep 
! cons 1 to k1-1
!

!interpolate
 
!        print *, ' getno-zmno ', maxval(zm), minval(zm)
 
        no(1:levs)=con_nzero
        call  interpol_wamz( no_nz16, z16, zm, levs, alt, no, kup, kdw ) 

!        print *, ' getno-zmnoi ', maxval(no), minval(no)
       
!
!extrapolate  above kup
!
        do k=kup+1,levs
           dx =log(pr(k-1))-log(pr(k))
           no(k)=no(k-1)*nair(k)/nair(k-1)*               &                    
              exp(dx*(1.-.5*amno*(1./am(k-1)+1./am(k))))
        enddo
	
!	
!attenuation below kdw
!	
        do k=kdw-1,1,-1
!          
           no(k)=no(kdw)*exp(-(kdw-k)*(kdw-k)*.1)
        enddo
	
!
! extra-check for positive no
!
        do k=1,levs
          no(k)=max(no(k), con_nzero)
        enddo
!
! incease [no] by a factor when kpa gt 5.0, same with ctipe, based on champ
! neutral density obervations 
!  
        do k=1,levs
          if (kpa.gt.5.0.and.kpa.le.6.0) then
            no(k) = no(k)*1.5
          elseif (kpa.gt.6.0.and.kpa.le.7.0) then
            no(k) = no(k)*2.5
          elseif (kpa.gt.7.0.and.kpa.le.8.0) then
            no(k) = no(k)*3.5
          elseif (kpa.gt.8.0.and.kpa.le.9.0) then
            no(k) = no(k)*4.5
          endif
        enddo
!
!
! check no
!       print *, ' getno-zmnoi2 ', maxval(no), minval(no)
!        print *, 'kup getno, iref', kup, iref, con_nzero
      return
      end subroutine wam_getno1d    
!           
     subroutine wamphys_coolno1(np,nps,t,o_n, o2_n, no_n, qno)   
!     
! sik may   2022 ccpp-version 
!vay oct 1  2015 clean-up
!    oct 28 2016 for new trunk
!      use   physcons, only : bz =>  con_boltz             
!-------------------------------------------------------------------------
! calculate no cooling   Kockarts, G. (1980). Nitric oxide cooling
!-------------------------------------------------------------------------
!  **
!  input:
!  t temperature profile k
!  o atomic oxygen number density profile m-3
!  no nitric oxide number density profile m-3
!  output:
!  qno: no cooling rate j/m3
!  ** 
      use machine,           only : kind_phys 
      implicit none
      integer, intent(in):: np                          ! numer of pressure levels
      integer, intent(in):: nps                         ! pressure index to start
      real(kind_phys), intent(in)   :: o_n(np),no_n(np) ! number density/m3   
      real(kind_phys), intent(in)   :: o2_n(np)           
      real(kind_phys), intent(in)   :: t(np)        ! temp (k)   
!out  
      real(kind_phys), intent(out)  :: qno(np)
      real(kind_phys) :: qnox(np)        
!locals     
      real(kind_phys) :: k10,hv, a10, g
      real(kind_phys) :: a1,a2,a3,om1,om     
      real(kind_phys) :: a10ghv, hvbz, a23   
      real(kind_phys) :: bz
      real(kind_phys) ::  no_deact, o1_rate, o2_rate,phot_e    
      integer i, k
!
!vay-2015/16
!                           ! hv=phot_e  = 3.726e-13_r8   at 5.3 mkm (erg)
!
       o1_rate = 2.7e-11*1.e-6      !3.3-3.6(-12)
       o2_rate = 2.4e-14*1.e-6
       phot_e  = 3.726e-13
       
       a10=13.3             !  trans_prob = 13.3_r8       
       bz=1.38e-23          ! boltzman
       k10=2.7e-17          ! k10*[o3p] vs 2.7 in waccm  k10=3.6e-17 
                            ! decreased because cm3/s = 10(-11)
                            !                   m3/s    10(-17)                                       
       hv=3.726e-20         ! in joules   1 erg is equal to 1.0e-7 joule.
       hvbz = hv/bz                                                                                                    
       g=1.0                                                             
!                                                     
       a2=5.4e-6*(1./(exp(hvbz/5800.)-1.))    
       a3=0.5*exp(-hvbz/247.5)
       a23 = a2+a3
       a10=13.3     
       a10ghv = a10*g*hv
!-------------------------------------------------
! k10/a10, hbvz, a10ghv, a23 => idear_solar_init
!=================================================
      qno(1:nps-1) =0.
      
      do i=nps,np                                                   
        om1=k10*o_n(i)          
        om=om1/(om1+a10)
        a1=exp(-hvbz/t(i)) 
        qno(i)=a10ghv*no_n(i)*om*(a1-a23)      ! should be in "j/kg/s" ~"j/s*[1/m3-no]"
	if (qno(i) < 0) qno(i) = 0.
      enddo
!
! Calculate NO cooling (ref: Kockarts, GRL, vol. 7, pp 137-140, 1980)
! WACCM
!      qnox(1:nps-1) =0.
!      do k=nps,np  
!      
!       no_deact =  o1_rate * o_n(k) + o2_rate * o2_n(k)      
!       qnox(k) = 1.e-7 * phot_e * A10 * no_n(k)* exp(-2700./t(k)) &
!                       * (no_deact / (no_deact + A10)) 
!      enddo      
!      qno = qnox
!====================================
! deactivation term trans_prob = 13.3_r8  hvbz=2700
!
!   nocool(i,k) = -1.e-4 * phot_e * trans_prob* exp(-2700./t(i,k)) &                     
!                  * no_conc * (no_deact / (no_deact + trans_prob)) 
! [no] 1/m3
! [om] dimensionless
! [a10] -dimensionless
!      a1*a10*hv*[no]/cp/rho = [K/sec]
!      a1 -dimens
! units:     hv*[no]/cp
!====================================
      return                                                            
      end  subroutine wamphys_coolno1
      
!
! merge solar radiation of WAM and GFS-lower atmosphere
!      
   subroutine wamphys_rad_merge(me, master, im,levs,xmu, prsl, hlw, swh, wtot,    & 
              dtrad, dtco2c,dtco2h,dth2oh,dth2oc,dto3)
	      
! 	 call wamphys_rad_merge(me, master, im ,levs, xmu, prsl, htrlw, htrsw, wtot, &
!                          dtrad, dtco2c,dtco2h,dth2oh,dth2oc,dto3)  
  
     use machine ,              only : kind_phys
     use wamphys_init_module,   only : rpref
     implicit none
!
      integer, intent(in) :: im  ! number of data points in hlw,dt..(first dim)
      integer, intent(in) :: me, master  
      integer, intent(in) :: levs             ! number of pressure levels
      
      real(kind=kind_phys), parameter     :: xb = 7.5, xt = 8.5, rdx=1./(xt-xb)
      
      
      real(kind=kind_phys), intent(in)    :: hlw(im,levs)     ! GFS lw rad (K/s)
      real(kind=kind_phys), intent(in)    :: swh(im,levs)     ! GFS sw rad (K/s)

      real(kind=kind_phys), intent(in)    :: prsl(im,levs)    ! pressure
      real(kind=kind_phys), intent(in)    :: xmu(im)          ! mormalized cos zenith angle
      
      real(kind=kind_phys), intent(in)    :: dtco2c(im,levs)  ! wam co2 cooling(K/s)
      real(kind=kind_phys), intent(in)    :: dtco2h(im,levs)  ! wam co2 heating(K/s)
      real(kind=kind_phys), intent(in)    :: dth2oc(im,levs)  ! wam h2o cooling(K/s)
      real(kind=kind_phys), intent(in)    :: dth2oh(im,levs)  ! wam h2o heating(K/s)
      real(kind=kind_phys), intent(in)    :: dto3(im,levs)    ! wam o3 heating(K/s)
      real(kind=kind_phys), intent(in)    :: dtrad(im,levs)   ! wam euv/srb heating(K/s) heating(K/s)
      real(kind=kind_phys), intent(out)   :: wtot(im,levs)    ! merged total radiation
      
!     local
      real(kind=kind_phys) :: rad_low(im, levs), rad_wam(im, levs)
      real(kind=kind_phys) :: xk,wl,wh, ah
      integer i,k,j
!
      do k=1,levs
        do i=1,im
	
          xk =-log(prsl(i,k) * rpref)
	  
          wh = dtco2c(i,k)+dth2oc(i,k)+dtco2h(i,k)+dth2oh(i,k)+dto3(i,k)+dtrad(i,k)
          wl = hlw(i,k) + swh(i,k)*xmu(i)
!	  rad_low(i,k) = wl
!	  rad_wam(i,k) = wh
          if(xk < xb) then
             wtot(i,k) = wl
          elseif(xk >= xb .and. xk <= xt) then
	     ah = (xk-xb)*rdx
             wtot(i,k) = wl*(1.-ah) + wh*ah
          else
             wtot(i,k) = wh
          endif
!	  if ( me == master .and. i == 15) then
!	    write(6,111) k, xk*7., wh*86400., wl*86400., wtot(i,k)*86400., xb
!	  endif
        enddo
      enddo
      
      return
      
  111   format(i4, 5(2x, F12.4), ' merge z-wh-wl-wm')    
          if ( me == master ) then
	     print *
	     print *, ' wamphys_rad_merge wtot ', maxval(wtot)*86400., minval(wtot)*86400.
	     print *, ' wamphys_rad_merge wlow ', maxval(rad_low)*86400., minval(rad_low)*86400.
	     print *, ' wamphys_rad_merge wupp ', maxval(rad_wam)*86400., minval(rad_wam)*86400.
	     print *, ' wamphys_rad_merge xmu ',  maxval(xmu), minval(xmu)
	     print *	     
	  endif  
      return
      end subroutine wamphys_rad_merge		
