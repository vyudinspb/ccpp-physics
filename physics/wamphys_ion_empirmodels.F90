!**  $id: chiu_model.f90,v 1.1.1.1 2006/06/04 18:19:13 cwplot exp $
!r
!r   chui_model.f90       chiu ionosphere, to return electron density.
!
!r   tiros_ionize_data  added the scaling of qiont with large gw and high tiros_activity_level
!

! wamphys_const, only : pi, elch, r_2_d, dtr
! wamphys_const, only : pi, dtr, pid12, pid6, rydays, pid2, pid3, pid4, pid9, pid18
! wamphys_set_data_ion, only  :nt_21, nt_20, nt_7, n_flx, n_bnd 
! wamphys_set_data_ion, only : emaps, cmaps, djspectra  
! wamphys_set_data_ion, only : jmaxwell, e0, width_maxwell, en_maxwell
! wamphys_set_data_ion, only : ratio, rlam, ion_recomb,lognpres
! wamphys_set_data_ion, only : width, en, te11,te15, ionchr 
!=================================================================
      module wamphys_ion_empirmodels
!
       use  machine,              only : kind_phys

       use  wamphys_const,       only :  pi, pi2, pid12, pid6, rydays, &
                                         pid2, pid3, pid4, pid9, pid18, dtr,  r_2_d 
				  
      
       use wamphys_const,         only :   elch, ydays, rydays
       
       contains
      
!=========================================================================
      subroutine wam_earth_chiu_model(sda,sza,thmag,phimr,rlt,rlat,    &     
              f107, dip, nday,ht1d,eden3d,ilon,lev1,ht_dim,lon_dim)
!
!vay-2016, cosmetic corrections for implicit none (ty and ilon)
! jun  2022 S. Karol ccpp-version
!
      
      implicit none
!
      integer,  intent(in)   :: nday,ht_dim,lon_dim
      integer,  intent(in)   :: lev1                    ! first level of ion-re K91

      real(kind=kind_phys),     intent(in)   :: sda, sza, thmag, phimr, rlt
      real(kind=kind_phys),     intent(in)   :: rlat, f107, dip

      real(kind=kind_phys),     intent(in)   :: ht1d(ht_dim)
! out  
      real(kind=kind_phys),     intent(out)  :: eden3d(lon_dim,ht_dim)  ! should 1d-array no neeeds for ilon-index



!locals

      integer                :: ilon, i, n
!
!*** start of declarations inserted by spag

      real(kind=kind_phys) :: abstmg , beta , cbp , cosrlt , costmg , cosza
      real(kind=kind_phys) :: g , g5 , g6 , g7 , g8 , gel , gel1 ,   gsm  
	 
      real(kind=kind_phys) ::   qel , rgamma , rho , rk  ,  sap ,sintmg 
      real(kind=kind_phys) ::   ty, ty1 , ty2  , w , wr , x , y
 
!*** end of declarations inserted by spag
  
!- define parameters

  
      real(kind=kind_phys) :: f(3) , pb(3) , s(3) , rd(3) , rl(3) , rt(3) , e(3) ,        &   
        u(3) , v(3) , p(3) , flong(3) , dipf(3), alp(3), a(3) ,rh(3), rr(3), fz(3), fn(3)
  
      real(kind=kind_phys) :: z(ht_dim)
  
!  
! absolutely no idea what these are. imported from tucan.f ????
! danger
!
      data alp/.5 , .5 , 1./
      data p/110. , 180. , 0./
      data a/1.36 , 2.44 , 0.66/
      data rh/10. , 34. , 0./
!
! data works only during the first call
!      
      save alp, p, a, rh
      
      rho = (f107-50.)*0.01
      ty = (nday+15.5)*12.*rydays       !/365.
      if ( ty > 12.0 ) ty = ty - 12.0
  
      abstmg = abs(thmag)
      cosza =  cos(sza)
      sintmg = sin(thmag)
      costmg = cos(thmag)
      cosrlt = cos(rlt)


      ty1 = sin(pid12*ty)
      ty2 = cos(pid6 *ty)
      p(1) = 110.
      p(2) = 180.
      f(1) = 0.0
      f(2) = 0.0
      pb(1) = 1.0
      pb(2) = 1.0
      s(1) = sqrt(1.0+1.15*rho)
      s(2) = sqrt(1.0+1.24*rho+0.25*rho*rho)
      rl(1) = 1.0
      rl(2) = 1.0
      e(1) = 1.0
      e(2) = 1.0
      flong(1) = 1.0
      flong(2) = 1.0
      dipf(1) = 1.0
      dipf(2) = 1.0
      g5 = sin(phimr)
      g6 = sin(phimr*.5)
      g7 = sqrt(abs(g5))
      g8 = cos(.5*(phimr-pi))
      sap = sin(sda)*sintmg
      f(3) = exp(-(2.92*sin(pid2-abstmg))**6)
      if ( thmag <= 0.0 ) then
         cbp = 0.0
         if ( g7 /= 0.0 ) cbp = ty1*(0.5*g6-0.5*g5-g6**8)-(1.0+ty1)   &  
         *ty2*g5/g7*exp(-4.0*g6*g6)
         pb(3) = (2.5+2.0*rho+ty2*(0.5+(1.3+0.5*rho)*g8**4)           &  
         +(1.3+0.5*rho)*cos(rlt-pi*(1.0+cbp)))                        &  
         *(1.0+0.4*ty1*ty1*exp(-ty1*g8**4))
      else
        wr = exp(-1.2*(cos(thmag-dtr*23.5*cosrlt)-costmg))
        pb(3) = (2.0+1.0*rho)*wr*(1.0+0.3*ty1)
      endif
      s(3) = (1.0+rho+0.204*rho**2+0.05*rho**3)
      if ( rho > 1.1 ) s(3) = 2.41 + 1.53*(s(3)-2.41)*(sintmg)**2
      p(3) = 240 + 75.0*rho + 83.0*rho*sap*costmg +                 &    
      30.0*cos(rlt-4.5*abs(thmag)-pi)                               &    
      + 10.0*costmg*cos(pid3*(ty-4.5))
      
      rd(1) = exp(2.0*(cosza/abs(cosza)*sqrt(abs(cosza))-1.0))
      rd(2) = exp((1.0+0.5*log(1.0+30.0*rho))*(cosza/abs(cosza)*sqrt(abs(cosza))-1.0))
      rd(3) = (0.9+0.32*sap)*(1.0+sap*(cos(rlt+pid4))**2)*exp(-1.1*(1.0+cos(rlt-0.873)))
      
      qel = 1.0 - 0.15*exp(-sqrt((12.0*thmag+1.05)**2+(ty/2.0-3.0)**2))
      rl(3) = (1.2-0.5*(costmg)**2)                                &     
        *(1.0+0.05*rho*(sintmg)**3*cos(pid6*ty))                   &   
        *(exp(3.0*cos(0.5*thmag*(sin(rlt)-1.0))))*qel
      w = cos(rlat+sda*cosrlt) - cos(rlat)
      rt(1) = exp(-0.4*w)
      rt(2) = exp(-0.25*w)
      beta = 1.3 + 0.139*rho**2 + 0.009*rho**3
      
      rk = 1.0 + 0.085*(cos(thmag-pid6)*(cos(pid12*(ty-2.0)))**3    &   
         +cos(thmag+pid4)*(cos(pid12*(ty-8.0)))**2)
      
      x = 0.7*(rk+0.178*rho**2/s(3)*cos(pid3*(ty-4.3)))          &     
        *exp(-beta*(cos(thmag+sda*cosrlt)-costmg))
      
      y = 0.2*(1.0-sin(abstmg-0.524))*(1.0+0.6*cos(pid3*(ty-3.94)))    & 
        *cos(pid6*(ty-1.0)) + (0.13-0.06*sin(abs(abstmg-pid9)))        & 
        *cos(pid3*(ty-4.5)) - (0.15+0.3*sin(abstmg))*(1.-cosrlt)**0.25 &
	*(cos(thmag+sda))**3
      rt(3) = x + y/s(3)
      g = (1.0+0.6*sqrt(rho)-0.2*rho)*exp(0.25*(1.0+cos(rlt-4.01)))
      gel = (costmg)**8*(cos(abstmg-0.262))**12
      gel1 = 1.0 + 0.05*(0.5-cos(pid3*ty)+cos(pid6*ty))
      e(3) = (1.0-0.4*(costmg)**10)                            &         
      *(1.0+0.6*(costmg)**10*(cos(rlt+pid4))**2)*(1.0+g*gel)   &      
      *gel1
      rgamma = 1.0 + 0.03*(0.5-cos(pid3*ty)+cos(pid6*ty))
      gsm=0.15-(1.0+rho)*(sin(thmag/2.0))**2*exp(-0.33*(ty-6.0)**2)
      flong(3) = 1.0 + 0.1*(costmg)**3*cos(2.0*(phimr-7.0*pid18))
      dipf(3)=rgamma*(1.0+gsm*exp(-18.0*(abs(dip)-2.0*pid9)**2))
      do i = 1 , 3
         u(i) = s(i)*rd(i)*rl(i)*rt(i)*e(i)*flong(i)*dipf(i)
         v(i) = f(i)*pb(i) + (1.0-f(i))*u(i)
      enddo
!
! vertica loop
!
      do n = lev1 , ht_dim
         z(n) = ht1d(n)/1000.
         if ( z(n) <= p(3) ) rh(3) = 2.0*(20.0+0.1*z(n))
         if ( z(n) > p(3)  ) rh(3) = 2.0*(20.0+0.1*p(3))
      do i = 1 , 3
        rr(i) = (z(n)-p(i))/rh(i)
        fz(i) = exp(alp(i)*(1.0-rr(i)-exp(-rr(i))))
        fn(i) = a(i)*fz(i)*v(i)
      enddo
       eden3d(ilon,n) = (fn(1)+fn(2)+fn(3))*1.e11
      enddo
      do n=1,lev1-1
            eden3d(ilon,n)=0.
      enddo

      return
      end subroutine wam_earth_chiu_model
!
       subroutine wam_tiros_ionize_data                               &
       (pres, lev1,levs,z,emaps,cmaps,djspectra,                  &
       grav,on,o2n, n2n,tn,gm_lat,essa1,tiros_activity_level,gw,  &
       eden_aurora1d)
!
!vay-2015: pass den = rho from the top-level program take-out comput/arrays ntot & meanmass
!          version with data/xxx/-statements.......eden_aurora1d   3-density due to aurora
!
!         output: of tiros_ionize_dat  CIRES-version 2015 with every time step in-n and data-stat
 
      implicit none
!      
      integer, intent(in) :: levs, lev1
      integer, parameter  :: jmaxwell = 6
      integer, intent(in) ::  tiros_activity_level
      real(kind=kind_phys), intent(in) :: pres(levs), z(levs),grav(levs),on(levs),    &
             o2n(levs),n2n(levs), tn(levs), gm_lat
      
      real(kind=kind_phys)  :: emaps(21,20,7),cmaps(21,20,7),djspectra(15,21)      
      real(kind=kind_phys) ::  essa1 
      real(kind=kind_phys) ::  gw
              
      real(kind=kind_phys) :: pres1(levs)

      real(kind=kind_phys) ::  gl, mlt

      real(kind=kind_phys), intent(out) :: eden_aurora1d(levs)     
     	     
      real(kind=kind_phys) :: bz, gscon, amu, e0
!

      real(kind=kind_phys) ::  ntot(levs),meanmass(levs)
!
      real(kind=kind_phys) ::  qiont(levs),ratio(21),rlam(21)   &
      ,den(levs),dl_lower,dl_upper,qiont_lower,qiont_upper      &
      ,mo,mo2,mn2,alpha, rno,range_en,pr,ratioz,rlamz,mh,mhe,q  & 
      ,qiont_o(levs),qiont_o2(levs),qiont_n2(levs)
      
      real(kind=kind_phys) :: width(15),en(15),te11(21),te15(21),width_maxwell
      
      real(kind=kind_phys) :: ionchr(21),ratio_ch,en_maxwell(jmaxwell),dl(jmaxwell)   &
      ,qion_maxwell(levs),lognpres(8),   ion_recomb(8),logpres,rr
!
      real(kind=kind_phys) :: ch , chi , dfac , diff , dprof , ed ,      &
          eflux , qdmsp , qt(levs),                              &
          ri , rj , th , swbz , offset, thmagd
	  
      integer i1 , i2 , j1 , j2 , k , kk , ld , n , nn , jj , jjj
      
      integer :: j, i, m, l, iband
            
      data en/.37,.6,.92,1.37,2.01,2.91,4.19,6.,8.56,12.18,           &
      17.3,24.49,36.66,54.77,81.82/
      data width/.158,.315,.315,.63,.631,1.261,1.26,2.522,            &
      2.522,5.043,5.043,10.,14.81,22.13,33.06/
      data rlam/1.49,1.52,1.51,1.48,1.43,1.37,1.30,1.22,              &
      11.12,1.01,0.895,0.785,0.650,0.540,0.415,0.320,0.225,           & 
      20.14,0.08,0.04,0.0/
!
      data ionchr/.378 , .458 , .616 , .773 , .913 , 1.088 , 1.403 ,    &
          1.718 , 2.033 , 2.349 , 2.979 , 3.610 , 4.250 , 4.780 ,      & 
          6.130 , 7.392 , 8.653 , 9.914 , 12.436 , 14.957 , 17.479/
	  
      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,  &
      -10.507/
      
      data ion_recomb/3.20e-13,3.20e-13,2.75e-13,1.45e-13,1.13e-13,    & 
       8.30e-14,3.70e-14,2.00e-14/
       
!      data dprof/4.23e19 , 5.62e19 , 5.77e19 , 5.70e19 , 1.04e19 ,    &
!         1.03e20 , 1.22e20 , 1.23e20 , 0.00e19 , 8.36e19 , 2.37e20 ,  &
!         2.61e20 , 0.00e19 , 0.00e18 , 3.07e20 , 5.26e20 , 0.00e19 ,  &
!         0.00e18 , 0.00e18 , 8.57e20/
!      data qdmsp/15*0.0/
    
      do 16 m=1,21
   16 ratio(m) = (m-1)*0.05
   
      do iband=1,21
        te15(iband)=0.0
        te11(iband)=0.0
! the ionization rates will need to be normalized to te11, which is
! the energy flux between 300ev and 20kev, which is provided by the
! tiros energy influx maps emaps, rather than the energy from
! 300ev to 100kev, which is what the spectra were normalized to
!
! check the energy influx is normalized to 1 erg/cm2/s

      do 17 m=1,15
   17 te15(iband)=te15(iband)+djspectra(m,iband)*en(m)*width(m)*1.6e-06
   
! normalize with the energy influx 300ev to 20kev

      do 18 m=1,11
   18 te11(iband)=te11(iband)+djspectra(m,iband)*en(m)*width(m)*1.6e-06
   
!      print *, iband, te11(iband), te15(iband)
      enddo
      bz = 1.38e-23
      gscon = 8.314e3
      mo = 16.
      mo2 = 32.
      mn2 = 28.
      mh = 1.
      mhe = 4.
      amu = 1.661e-27
      e0=0.035
      width_maxwell=0.050
      do j = 1,jmaxwell
      en_maxwell(j) = j*0.05 - 0.025
      enddo
! initialize qiont, etc
      do i=1,levs
      qiont(i) = 0.0
      qion_maxwell(i) = 0.0
      qiont_o(i) = 0.0
      qiont_o2(i) = 0.0
      qiont_n2(i) = 0.0
      eden_aurora1d(i) = 0.0
      enddo
! convert magnetic latitude from radians to degrees
      thmagd = gm_lat * r_2_d
!     print *, 'gm_lat   thmagd  essa1', gm_lat, thmagd, essa1
      th = abs(thmagd) - 50.
      if(abs(thmagd).le.50.) goto 200
! calculate magnetic hour angle from noon in gregrees
!      essa1 = (mlt + 12.)*15.
! now passed essa1 directly
      if ( essa1.ge.360.0 ) then
          essa1 = essa1 - 360.0
      elseif ( essa1.lt.0.0 ) then
          essa1 = essa1 + 360.
      endif
!!  **
      l = tiros_activity_level - 2
      if ( l.lt.1 ) l = 1
      if ( l.gt.7 ) l = 7
! define dfac to scale qiont later with large gw and tiros_activity_level
! added by zhuxiao.li

      dfac = 1.0
      if (tiros_activity_level.gt.9 .and. gw.gt.96.0)dfac = gw/96.0

   
      ri = essa1/18.0 + 11.
      i1 = ri                   ! i1 =int(ri) ?
      ri = ri - i1
      if ( i1.gt.20 ) i1 = i1 - 20
      i2 = i1 + 1
      if ( i2.gt.20 ) i2 = i2 - 20
      rj = th/2. + 1.
      j1 = rj                   ! j1 =int(rj) ?
      rj = rj - j1
      j2 = j1 + 1
!
      eflux = rj*ri*emaps(j2,i2,l) + (1.-rj)*ri*emaps(j1,i2,l)              &
              + rj*(1.-ri)*emaps(j2,i1,l) + (1.-rj)*(1.-ri)*emaps(j1,i1,l)
      eflux = 10.**(eflux)*1.e-3
!      print *, 'eflux   ', eflux
!
      ch = rj*ri*cmaps(j2,i2,l) + (1.-rj)*ri*cmaps(j1,i2,l) + rj*(1.-ri) &
           * cmaps(j2,i1,l) + (1.-rj)*(1.-ri)*cmaps(j1,i1,l)
!      print *, 'ch   ', ch
! validation tests:
! to compare with figure 4 or 5 f-r and evans 1987
! set ch to 5 different mean energies and
! set eflux to 1.0 mw/m2
!       ch = 2.98
!       eflux=1.0
! a useful thing to compare is the ionization rate profile for ch=2.98
! with the equivalent output assuming a maxwellian spectrum using the
! other code ionize_ipe_3 with the same mean energy of 2.98
! in this case the profiles are similar, at other other values of ch
! they can be quite different.
!
!      print *, 'set for test ch eflux', ch, eflux
!
      if ( ch.lt.0.378 ) ch = 0.379
!      if ( ch.gt.17.479 ) write (6,99001) ch
!
      do 300 kk = 2 , 21
         if ( ch.le.ionchr(kk) ) then
            k = kk - 1
            goto 400
         endif
 300  continue
 400  chi = ch - ionchr(k)
      diff = ionchr(kk) - ionchr(k)
      ratio_ch = chi/diff
!      if(ratio_ch.gt.1.) print *, 'ratio_ch out of bounds', ratio_ch

99001 format ('  ch value outof bound in tiros',f10.6)
! loop through ipe height levels


      do 2000 i=lev1,levs
! stop at 1000km altitude

      if(z(i)*1.e-3 .le. 1000.) then                ! zwam < 1000, vay    instead goto    

      ntot(i) = on(i)+o2n(i)+n2n(i)
      meanmass(i) = (on(i)*mo+o2n(i)*mo2+n2n(i)*mn2)/ntot(i)
      den(i) = pres(i)*meanmass(i)/(gscon*tn(i))
!
! calculate ion recombination rate
!      lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,-10.507/
      logpres = log(pres(i))

      if (logpres.ge.-4.835) rr=3.20e-13

      
      if (logpres.le.-10.507) rr=2.00e-14

         if (logpres.gt.-10.507.and.logpres.lt.-4.835)  then    
          do jjj=3,8
             if(lognpres(jjj).le.logpres)then
                  jj=jjj-1
                  exit
           endif
          enddo
      

!err      rr = ion_recomb(jj)-logpres*(ion_recomb(jjj)-ion_recomb(jj))/
!err     &(lognpres(jjj)-lognpres(jj))
!update vay-2016/10
!nems.x             0000000000adf0a4  tiros_ionize_data         449/450  idea_ion_empirmodels.f
!ion_recomb(jjj =>  ion_recomb(jj)
!
      rr = ion_recomb(jjj)-(lognpres(jjj)-logpres)*                  &
      (ion_recomb(jjj)-ion_recomb(jj))/(lognpres(jjj)-lognpres(jj))
      
     endif
     
!
!

 3000 format(1x,i5,2f10.2,1p3e9.2)
! loop through all energy bands


      do 10 l=1,15
!      dl(l)=rno*en(l)*exp(-en(l)/alpha)
      dl_lower = djspectra(l,k)
      dl_upper = djspectra(l,kk)
      range_en=4.57e-05*en(l)**1.75
      pr=range_en*grav(i)
      ratioz=pres(i)/pr
      
      if(ratioz.gt.1.0) goto 20
      
      do 12 m=1,21
      if(ratioz.gt.ratio(m)) goto 12
      rlamz=rlam(m-1)+(ratioz-ratio(m-1))*(rlam(m)-rlam(m-1))/(ratio(m)-ratio(m-1))
      goto 13
   12 continue
   13 continue
      goto 21
      
   20 rlamz=0.0
   21 continue
   
      qiont_lower=den(i)*en(l)*rlamz*dl_lower*width(l)*1.e7/range_en/e0
      qiont_upper=den(i)*en(l)*rlamz*dl_upper*width(l)*1.e7/range_en/e0
      qiont(i)=qiont(i)+(ratio_ch*qiont_upper+(1.-ratio_ch)*qiont_lower)* &
               eflux/(ratio_ch*te11(kk)+(1.-ratio_ch)*te11(k))
   11 continue
   10 continue
! add 0 - 300ev as maxwellian with ch, and eflux
      alpha = ch/2.
      rno=eflux*6.24e12/2./alpha**3
      do 110 l=1,jmaxwell
      dl(l)=rno*en_maxwell(l)*exp(-en_maxwell(l)/alpha)
      range_en=4.57e-05*en_maxwell(l)**1.75
      pr=range_en*grav(i)
      ratioz=pres(i)/pr
      
      if(ratioz.gt.1.0) goto 120
      
      do 112 m=1,21
      if(ratioz.gt.ratio(m)) goto 112
      rlamz=rlam(m-1)+(ratioz-ratio(m-1))*(rlam(m)-rlam(m-1))/(ratio(m)-ratio(m-1))
       goto 113
  112 continue
  
  113 continue
  
      goto 121
      
  120 rlamz=0.0
  
  121 continue
  
      qion_maxwell(i)=qion_maxwell(i)+den(i)*en_maxwell(l)*rlamz*dl(l)   &
          *width_maxwell/range_en/e0
	  
  110 continue
      qiont(i) = qiont(i) + qion_maxwell(i)

!  the qiont scaled by dfac with large gw and tiros_activity level
!  added by zhuxiao
       qiont(i) = qiont(i)*dfac
!
!vay-2016, extra security >0 and non-zero
!
      if (rr.gt.0.) then 
         eden_aurora1d(i)=sqrt(qiont(i)/rr)
      else
         eden_aurora1d(i)=0.
      endif
!
      q=qiont(i)/(0.92*n2n(i)+1.5*o2n(i)+0.56*on(i))
      qiont_o(i)=(0.5*o2n(i)+0.56*on(i))*q
      qiont_o2(i)=o2n(i)*q

      endif
 2000 continue

  200 continue

    
!
      return
      end subroutine wam_tiros_ionize_data
      
!=======================================================================================
! clean version with wamphys_set_data_ion
!=======================================================================================
      subroutine wam_tiros_ionize(lev1,levs, pres, den,        &
       z, grav,on,o2n,n2n,tn,                              &  
       gm_lat,essa1,tiros_activity_level,gw, eden_aurora1d)
!
! version with idea_ion_input
!
!oct 2016 vay: a) err-interp logpress
!              b) range (special fort function) => range_en
!              c) clarify if ( ch.lt.0.378 ) ch = 0.379             ! ???????? 378-379
!              d)        if ( ch.gt.17.479 ) write (6,99001) ch     ! out of bounds ???
!
!              e) p-wam = constabt, suggestion rlamz(15, levs) -----constant with time
!
      use wamphys_set_data_ion, only : nt_21, nt_20, nt_7, n_flx, n_bnd 
      use wamphys_set_data_ion, only : emaps, cmaps, djspectra
!
      use wamphys_set_data_ion, only : jmaxwell     ! 6-21-15
      use wamphys_set_data_ion, only : ratio, rlam, ion_recomb,lognpres, &
                                width, en, te11,te15, ionchr,     &
                                width_maxwell, en_maxwell
!
!
      implicit none
!
      integer, intent(in) ::  levs, lev1
      integer, intent(in) ::  tiros_activity_level
      real(kind=kind_phys) ::  essa1 
      real(kind=kind_phys) ::  gw  
      real(kind=kind_phys) ::  z(levs), grav(levs), pres(levs), den(levs)
      real(kind=kind_phys) ::  tn(levs)
      real(kind=kind_phys) ::  on(levs),o2n(levs),n2n(levs)
      real(kind=kind_phys) ::  gl, mlt,  gm_lat
!out
      real(kind=kind_phys), intent(out) ::   eden_aurora1d(levs) 
!
!
!local
!
!
      real(kind=kind_phys) ::  dl (jmaxwell)
      real(kind=kind_phys) ::  alpha
      real(kind=kind_phys) ::  qiont(levs)
      real(kind=kind_phys) ::  dl_lower,dl_upper,qiont_lower,     &
                               range_en,pr,ratioz,rlamz,q, ratio_ch
			       
      real(kind=kind_phys) :: qmid
      real(kind=kind_phys) :: qiont_o(levs),qiont_o2(levs),qiont_n2(levs), qt(levs)
      
      real(kind=kind_phys) :: qion_maxwell(levs)
      real(kind=kind_phys) :: qiont_upper
      real(kind=kind_phys) :: rno

      real(kind=kind_phys) :: ch , chi , dfac , diff , dprof , ed ,logpres,rr
      real(kind=kind_phys) :: eflux, qdmsp,  ri , rj , th , swbz , offset, thmagd
      
      real(kind=kind_phys) :: e0      
          
	  

      integer i1 , i2 , j1 , j2 , k , kk , ld , n , nn , jj , jjj
      integer :: j, i, m, l, iband 
! initialize qiont, etc

      e0 =0.035
      do i=1,levs
      qiont(i) = 0.0
      qion_maxwell(i) = 0.0
      qiont_o(i) = 0.0
      qiont_o2(i) = 0.0
      qiont_n2(i) = 0.0
!...................................... zero output
      eden_aurora1d(i) = 0.0
      enddo
    
! convert magnetic latitude from radians to degrees
      thmagd = gm_lat * r_2_d
!!!!      print *, 'gm_lat   thmagd  essa1', gm_lat, thmagd, essa1
      th = abs(thmagd) - 50.
      if(abs(thmagd).le.50.) goto 200
! calculate magnetic hour angle from noon in gregrees
!      essa1 = (mlt + 12.)*15.
! now passed essa1 directly
      if ( essa1.ge.360.0 ) then
          essa1 = essa1 - 360.0
      elseif ( essa1.lt.0.0 ) then
          essa1 = essa1 + 360.
      endif
!
      l = tiros_activity_level - 2
!
! wam-limits from 1 to 7
!
      if ( l.lt.1 ) l = 1
      if ( l.gt.7 ) l = 7
! define dfac to scale qiont later with large gw and tiros_activity_level
! added by zhuxiao.li
      dfac = 1.0
      if (tiros_activity_level.gt.9 .and. gw.gt.96.0)  &
        dfac = gw/96.0
!
! what's this vay
!
      ri = essa1/18.0 + 11.
      i1 = int(ri)                   ! i1 = ri
      ri = ri - i1                   ! ri = ri -int(ri) ???? now
!
      if ( i1.gt.20 ) i1 = i1 - 20
      i2 = i1 + 1
      if ( i2.gt.20 ) i2 = i2 - 20
      rj = th/2. + 1.
      j1 = int(rj)                   ! j1 =rj
      rj = rj - j1
      j2 = j1 + 1
! both rj <1 & ri < 1
      eflux = rj*ri*emaps(j2,i2,l) + (1.-rj)*ri*emaps(j1,i2,l)    &
             + rj*(1.-ri)*emaps(j2,i1,l) + (1.-rj)*(1.-ri)        &
             *emaps(j1,i1,l)
!
      eflux = 10.**(eflux) * 1.e-3

!
!
      ch = rj*ri*cmaps(j2,i2,l) + (1.-rj)*ri*cmaps(j1,i2,l) + rj*(1.-ri)     &
          *cmaps(j2,i1,l) + (1.-rj)*(1.-ri)*cmaps(j1,i1,l)
!
!      if (mpi_id == 0) print *, 'ion-ch   ', ch
! validation tests:
! to compare with figure 4 or 5 f-r and evans 1987
! set ch to 5 different mean energies and
! set eflux to 1.0 mw/m2
!       ch = 2.98
!       eflux=1.0
! a useful thing to compare is the ionization rate profile for ch=2.98
! with the equivalent output assuming a maxwellian spectrum using the
! other code ionize_ipe_3 with the same mean energy of 2.98
! in this case the profiles are similar, at other other values of ch
! they can be quite different.
!       if (mpi_id == 0) then
!       print *, 'set for test ch eflux', ch, eflux
!       endif
      if ( ch.lt.0.378 ) ch = 0.379             ! ???????? 378-379
      if ( ch.gt.17.479 ) write (6,99001) ch
!
      do 300 kk = 2 , 21
         if ( ch.le.ionchr(kk) ) then
            k = kk - 1
            goto 400
         endif
 300  continue
 400  chi = ch - ionchr(k)
      diff = ionchr(kk) - ionchr(k)
      ratio_ch = chi/diff
!      if(ratio_ch.gt.1.) print *, 'ratio_ch out of bounds', ratio_ch
!
99001 format ('  ch value outof bound in tiros',f10.6)

! loop through ipe height levels or wam-levels

      do 2000 i=lev1,levs
!

      if(z(i)*1.e-3 .gt. 1000.) exit   ! do we need this z-wam < 650 km


!
! new z-interpolation   for rr > 0 (apparently)
!
      logpres = log(pres(i))
      if (logpres.ge.-4.835)then
      rr=3.20e-13
      goto 450
      endif
      if (logpres.le.-10.507)then
      rr=2.00e-14
      goto 450
      endif
      do jjj=3,8
         if(lognpres(jjj).le.logpres)then
         jj=jjj-1
         goto 451
         endif
      enddo
  451 continue

!fix vay oct-2016
!
      rr = ion_recomb(jjj)-(lognpres(jjj)-logpres)*     &
      (ion_recomb(jjj)-ion_recomb(jj))/(lognpres(jjj)-lognpres(jj))


  450 continue

 
!
! loop through all energy bands
      do 10 l=1,15
!      dl(l)=rno*en(l)*exp(-en(l)/alpha)
      dl_lower = djspectra(l,k)
      dl_upper = djspectra(l,kk)

      range_en=4.57e-05*en(l)**1.75
      pr=range_en*grav(i)
      ratioz=pres(i)/pr
      if(ratioz.gt.1.0) goto 20
      do 12 m=1,21
      if(ratioz.gt.ratio(m)) goto 12
      rlamz=rlam(m-1)+(ratioz-ratio(m-1))*(rlam(m)-rlam(m-1))/     &
           (ratio(m)-ratio(m-1))
      goto 13
   12 continue
   13 continue
      goto 21
   20 rlamz=0.0
   21 continue
! ...................... again interpolation in ratio,,,ratioz
! similar suggestion rlamz(15, levs) -----constant with time
!
      qmid = den(i)*en(l)*rlamz*width(l)*1.e7/range_en/e0
      qiont_lower= qmid*dl_lower
      qiont_upper= qmid*dl_upper
      qiont(i)=qiont(i)+(ratio_ch*qiont_upper+(1.-ratio_ch)*qiont_lower)    &
         *eflux/(ratio_ch*te11(kk)+(1.-ratio_ch)*te11(k))
   11 continue

   10 continue     ! loop through all energy bands

! add 0 - 300ev as maxwellian with ch, and eflux

      if (ch.le.0.) then
!         print *, ' tiros: ch=0 ', ch
         ch =2.98
      endif
      alpha = .5*ch

      rno=eflux*6.24e12/2./alpha**3
!
      do 110 l=1,jmaxwell
      dl(l)=rno*en_maxwell(l)*exp(-en_maxwell(l)/alpha)
      range_en=4.57e-05*en_maxwell(l)**1.75
      pr=range_en*grav(i)
      ratioz=pres(i)/pr
      if(ratioz.gt.1.0) goto 120
      do 112 m=1,21
      if(ratioz.gt.ratio(m)) goto 112
      rlamz=rlam(m-1)+(ratioz-ratio(m-1))*(rlam(m)-rlam(m-1))/     &
           (ratio(m)-ratio(m-1))
      goto 113
  112 continue
  113 continue
      goto 121
  120 rlamz=0.0
  121 continue
!......................similar rlamz_jmx(levs, 6)
!                       
      qion_maxwell(i)=qion_maxwell(i)       &
       +den(i)*en_maxwell(l)*rlamz*dl(l)*width_maxwell/range_en/e0
!
  110 continue   !  loop for300ev as maxwellian with ch, and eflux
!
      qiont(i) = qiont(i) + qion_maxwell(i)
!  the qiont scaled by dfac with large gw and tiros_activity level
!  added by zhuxiao
       qiont(i) = qiont(i)*dfac

!vay-2016
      if (rr.gt.0.) then 
         eden_aurora1d(i)=sqrt(qiont(i)/rr)
      else
         eden_aurora1d(i)=0.
      endif
!
! proxies for charge from neutrals ????
!
      q=qiont(i)/(0.92*n2n(i)+1.5*o2n(i)+0.56*on(i))
      qiont_o(i)=(0.5*o2n(i)+0.56*on(i))*q
      qiont_o2(i)=o2n(i)*q
      qiont_n2(i)=0.92*n2n(i)*q
      
 2000 continue
 
  200 continue 
!
      return
      end subroutine wam_tiros_ionize
!
      end module wamphys_ion_empirmodels
