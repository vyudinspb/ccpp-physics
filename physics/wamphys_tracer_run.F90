! Apr 06 2012    Henry Juang, initial implement for nems
! Oct 20 2015    Weiyu Yang,  add f10.7 inputted data.
! Oct 15 2016    VAY, put f10.7 as a parameter
! Nov    2016    Correction of JO2 and oxygen chemistry
! September 2017 Weiyu Yang, add IPE back coipling WAM code.
! October 2017   Rashid Akmaev, eddy mixing, corrections and clean-up
! May 2022       Svetlana Karol and CUA
!-----------------------------------------------------------------------
  subroutine wamphys_tracer_run(me, master, im,levs,ntrac,ntrac_i, nto, nto2, nto3, ntqv,   & 
       dayno, dtp, grav,prsi,prsl, adt, q,   n1, n2, ozn, n3, nair, rho, am, am29,          &
       cospass, zg, jrates_O2, f107,f107d, go2dr, plow, phigh, xpk_low, xpk_high) 
      
      use machine,  only             : kind_phys 
      use physcons, only             : avgd => con_avgd             
      use wamphys_init_module, only  : bz,amo,amn2, amo2, amo3, amh2o
      use wamphys_init_module, only  : rbz, rmo, rmo2, rmn2, rmh2o,rmo3
      use wamphys_init_module, only  : lowst_ipe_level => lowipe_lev150     
      
!      use module_IPE_to_WAM, only: low_ipe_level   !  wamphys_merge_ipe2wam
      
      implicit none
! Argument
      integer, intent(in) :: me, master              ! current PE
      integer, intent(in) :: im                      ! number of horiz data points
      
      integer, intent(in) :: nto, nto2, nto3, ntqv   ! more tracers
      
      integer, intent(in) :: levs            ! number of pressure levels
      integer, intent(in) :: ntrac           ! number of tracer (total)
      integer, intent(in) :: ntrac_i         ! number of tracer add by IDEA
      integer, intent(in) :: dayno           ! calendar day
      
      real(kind=kind_phys), intent(in)    :: dtp             ! time step in second      
      real(kind=kind_phys), intent(in)    :: cospass(im)     ! cos zenith angle
      real(kind=kind_phys), intent(in)    :: zg(im,levs)     !layer height (m)
      real(kind=kind_phys), intent(in)    :: f107, f107d     ! variable F107 to recomput JJ
      real(kind=kind_phys), intent(in)    :: prsi(im,levs+1) ! interface pressure in Pa
      real(kind=kind_phys), intent(in)    :: prsl(im,levs)   ! layer pressure in Pa
      real(kind=kind_phys), intent(in)    :: grav(im,levs)   ! (m/s2)
      real(kind=kind_phys), intent(in)    :: adt(im,levs)    ! input  temp 
      
      real(kind=kind_phys), intent(in), dimension(im)    :: xpk_low, xpk_high    
      integer,              intent(in), dimension(im)    :: plow, phigh    
!      real(kind=kind_phys), intent(in)    :: go2dr(im, lowst_ipe_level:levs)
      real(kind=kind_phys), intent(in)    :: go2dr(im,levs)
            
      real(kind=kind_phys), intent(inout) :: q(im,levs,ntrac)   ! input output tracer
      real(kind=kind_phys), intent(out)   :: n1(im,levs)     ! number density of o (/m3)
      real(kind=kind_phys), intent(out)   :: n2(im,levs)     ! number density of o2 (/m3)
      real(kind=kind_phys), intent(out)   :: ozn(im,levs)    ! number density of o2 (/m3)
      real(kind=kind_phys), intent(out)   :: n3(im,levs)     ! number density of n2 (/m3)
      real(kind=kind_phys), intent(out)   :: nair(im,levs)   ! total number density (/m3)
      real(kind=kind_phys), intent(out)   :: rho(im,levs)    ! density of  (kg/m3)
      real(kind=kind_phys), intent(out)   :: am(im,levs)     ! avg mass of mix  (kg)
      
      real, intent(out)   :: am29(im,levs)   ! avg mass of air 28.84 => 16.
! local argument
      real(kind=kind_phys) ::  dq1(im,levs,ntrac_i),dq2(im,levs,ntrac_i),  &
                               dq3(im,levs,ntrac_i)
	      
      real(kind=kind_phys) ::  mh2o,mo3, mo,mo2,mn2        
      real(kind=kind_phys) ::  qin(im,levs,ntrac_i)
      real(kind=kind_phys) ::  qn2
      integer i,k,in, itr, i2
      real(kind=kind_phys), intent(out)  :: Jrates_O2(im, levs)            ! O2 dissociation rate
      
      real(kind=kind_phys)               :: Jrates_O3(im, levs)
            
      logical, parameter :: debug =.false.
      
!     Two WAM tracers: O and O2   indices: nto, nto2, nto3
      itr = nto
      if (nto > nto2) itr =nto2     
        do in=1,ntrac_i  
            i2 =  itr-1+in    
          do i=1,im
            do k=1,levs
              qin(i,k,in)=q(i, k, i2)   ! take two last tracers O and O2
  !           qin(i,k,1)=q(i, k, nto) 
  !	        qin(i,k,2)=q(i, k, nto2)  
            enddo
          enddo
        enddo
  
! mean mass, mass and number densities
!     here n,n1,n2 in /m3 , rho in kg/m3
      do i=1,im
        do k=1,levs
	
           qn2=1.-q(i,k,ntqv)-q(i,k,nto3)-qin(i,k,1)-qin(i,k,2)

! mean molecular mass of gaseous tracers
           am(i,k)=1./(qin(i,k,1)*rmo+qin(i,k,2)*rmo2+q(i,k,ntqv)*rmh2o+     &  
               q(i,k,nto3)*rmo3+qn2*rmn2)

! total number density and mass density
           nair(i,k)=rbz*prsl(i,k)/adt(i,k)
           rho(i,k)=am(i,k)*nair(i,k)

! partial number densities for radiation and chemistry,
!     make sure non-negative
           n1(i,k)=max(qin(i,k,1)*rho(i,k)*rmo, 0.)
           n2(i,k)=max(qin(i,k,2)*rho(i,k)*rmo2,0.)
           n3(i,k) =max(qn2*rho(i,k)*rmn2,0.)
 	  
           ozn(i,k)= max(q(i,k,nto3)*rho(i,k)*rmo3,0.)
	  
        enddo
      enddo
      
      if (debug) then
         print *, ' qin1 ', maxval(qin(:,:,1)), minval(qin(:,:,1))   
         print *, ' qin2 ', maxval(qin(:,:,2)), minval(qin(:,:,2)) 
                    
         print *, ' nair ', maxval(nair), minval(nair)  
         print *, ' rho ', maxval(rho), minval(rho) 
         print *, ' amu ', maxval(am), minval(am)     
         print *, ' rmo-s ', rmo, rmn2, rmo2, rmo3, rbz
         print *, ' zg ', maxval(zg), minval(zg)
         print *, ' adt ', maxval(adt), minval(adt)       
         print *, ' grav ', maxval(grav), minval(grav)
         print *, ' o_n ',  maxval(n1), minval(n1)          
         print *, ' o2_n ', maxval(n2), minval(n2)
         print *, ' n2_n ', maxval(n3), minval(n3) 
         print *, ' o3_n ', maxval(ozn), minval(ozn)     
      endif
! Eddy mixing

      call wam_tracer_eddy(im,levs,ntrac_i,grav,prsi,prsl,rho,dtp,  &
          qin,dayno,dq3)

! Mutual molecular diffusion of major thermospheric species O, O2, N2

      call wam_tracer_m(me, master,im,levs,ntrac_i,grav,prsi,prsl,adt,dtp,     &
          qin,am,dq1)
 
		
! O2 dissociation rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!no radiation!!!
      call wamphys_dissociation(im, levs, adt, cospass, n1, n2,  ozn,  n3,  &
               dayno, zg, grav,  f107,  f107d, Jrates_O2, Jrates_O3)

! Merge the  IPE back coupling WAM variable arrays into WAM.  ????   GO2DR => jrates_O2
!
!      IF (ipe_to_wam_coupling) THEN
!        call wamphys_merge_ipe2wam(go2dr, jrates_O2,         &
!         im, levs, lowst_ipe_level, prsl, plow, phigh,  &
!         xpk_low, xpk_high)
!      END IF

! Oxygen chemistry + ad hoc HOx sinks of O
      call wam_tracer_c(im,levs,ntrac_i,adt,dtp,Jrates_O2,  &
          n1,n2,nair,rho, qin,dq2)

! Update mmr of O and O2 due to mixing, diffusion, and chemistry
      do in=1,ntrac_i
        i2 =  itr-1+in    
        do i=1,im
          do k=1,levs
!           q(i,k,nto) = q(i,k,nto)+dq1(i,k,1)+dq2(i,k,1)+dq3(i,k,1)
!	      q(i,k,nto2)= q(i,k,nto2)+dq1(i,k,2)+dq2(i,k,2)+dq3(i,k,2)
            q(i,k, i2)= q(i,k,i2) +dq1(i,k,in)+dq2(i,k,in)+dq3(i,k,in)
          enddo
        enddo
      enddo

! Update number densities and rho (nair= p/(kT) is conserved) and mean mass
      do i=1,im
         do k=1,levs
            qn2=1.-q(i,k,ntqv)-q(i,k,nto3)-q(i,k,nto)-q(i,k,nto2)
            am(i,k)=1./(q(i,k,nto)*rmo+q(i,k,nto2)*rmo2+q(i,k,ntqv)*rmh2o+      & 
                q(i,k,nto3)*rmo3+qn2*rmn2)
	      rho(i,k) =am(i,k)*nair(i,k)
            n1(i,k)=max(q(i,k,nto)*rho(i,k)*rmo,0.)
            n2(i,k)=max(q(i,k,nto2)*rho(i,k)*rmo2,0.)
            n3(i,k) =max(qn2*rho(i,k)*rmn2,0.)
            ozn(i,k)= max(q(i,k,nto3)*rho(i,k)*rmo3,0.)
         enddo
      enddo
      
      am29 = am * 1.e3 * avgd

      return
  end subroutine wamphys_tracer_run

  subroutine wam_tracer_m(me, master,im,levs,ntrac_i,grav,prsi,prsl,adt,  & 
      dtp,qin,am,dq)
!-----------------------------------------------------------------------
! Calculate tracer changes by mutual molecular diffusion of
!     major thermospheric species O, O2, N2
! 2006 Rashid Akmaev and Fei Wu
!-----------------------------------------------------------------------
      use machine, only   : kind_phys
      use physcons,  only :rgas=>con_rgas, avgd => con_avgd
      
      use wamphys_init_module, only:  amo, amo2, amn2, bz, rbz
      
      implicit none
! Argument
      integer, intent(in) :: im    ! number of long data points in fields
      integer, intent(in) :: me, master 
      integer, intent(in) :: levs  ! number of pressure levels
      integer, intent(in) :: ntrac_i ! number of tracer add by IDEA
      real(kind=kind_phys),    intent(in) :: dtp   ! time step in second
      real(kind=kind_phys), intent(in)    :: prsi(im,levs+1) ! interface pressure in KPa
      real(kind=kind_phys), intent(in)    :: prsl(im,levs)   ! layer pressure in KPa
      real(kind=kind_phys), intent(in)    :: grav(im,levs)   ! (m/s2)
      real(kind=kind_phys), intent(in) :: adt(im,levs)   ! input  temp at dt=0
      real(kind=kind_phys), intent(in) :: qin(im,levs,ntrac_i)   ! input tracer
      real(kind=kind_phys), intent(in)   :: am(im,levs)   ! avg mass of mix  (kg)
      real(kind=kind_phys), intent(out):: dq(im,levs,ntrac_i) ! output tracer changes
!local  variables

      real(kind=kind_phys) :: n1_i(levs+1),n2_i(levs+1),n3_i(levs+1),n_i(levs+1)
      real(kind=kind_phys) :: t_i(levs+1),am_i(levs+1),qout(im,levs,ntrac_i)
      real(kind=kind_phys) :: beta(2,2,levs+1),a(2,2,levs),b(2,2,levs),c(2,2,levs)
      real(kind=kind_phys) :: ggg(2,2),ee(2,2,levs+1),f(2,levs+1),   &                      
          d12,d13,d23,a12,a13,a23,s12,s13,s23,mo,mo2,mn2,            &  
          dp1(levs),dp1_i(levs+1)
      real(kind=kind_phys) :: partb_i(levs+1),parta(levs),hold1,dtp1,hold2
      integer k,i,kk,kk1,in
! change unit from g/mol to kg
!     if ( me == master ) then
!	    print *, ' grav_molec_m ', maxval(grav(:,147)), minval(grav(:,147)) 
!	    print *, ' am_molec_m ', maxval(am), minval(am) 
!     endif

      mo=amo*1.e-3/avgd
      mo2=amo2*1.e-3/avgd
      mn2=amn2*1.e-3/avgd
! some constants
      a12=9.69e18
      a13=9.69e18
      a23=8.3e18

      s12=0.774
      s13=0.774
      s23=0.724
! set boundary
      beta(1:2,1:2,1)=0.
      beta(1:2,1:2,levs+1)=0.
      a(1:2,1:2,1)=0.
      c(1:2,1:2,levs)=0.
      ee(1:2,1:2,levs+1)=0.
      f(1:2,levs+1)=0.
!
      dtp1=1./dtp
      t_i=0.
      am_i=0.
      n_i=0.
      n1_i=0.
      n2_i=0.
      n3_i=0.
!
! for each longitude
!
      do i=1,im
! calculate temp in interface pressure levels
! get compositions at interface pressure levels
        do k=2,levs
          t_i(k)=(adt(i,k-1)+adt(i,k))*.5
          am_i(k)=.5*(am(i,k-1)+am(i,k))
          n_i(k)=prsi(i,k)/bz/t_i(k) 
          n1_i(k)=max(0.,                &
              .5*(qin(i,k,1)+qin(i,k-1,1))*am_i(k)*n_i(k)/mo)
          n2_i(k)=max(0.,                & 
              .5*(qin(i,k,2)+qin(i,k-1,2))*am_i(k)*n_i(k)/mo2)
          n3_i(k)=max(0.,n_i(k)-n1_i(k)-n2_i(k))
        enddo

! calculate beta at interface pressure
        do k=2,levs
          d12=a12*t_i(k)**(s12)
          d13=a13*t_i(k)**(s13)
          d23=a23*t_i(k)**(s23)
          hold1=1./(n1_i(k)*d23+n2_i(k)*d13+n3_i(k)*d12)
          beta(1,1,k)=hold1*d13*mo*(n1_i(k)*mn2*d23+    &                
                 (n2_i(k)*mo2+n3_i(k)*mn2)*d12)
          beta(2,2,k)=hold1*d23*mo2*(n2_i(k)*mn2*d13+   &               
                 (n1_i(k)*mo+n3_i(k)*mn2)*d12)
          beta(1,2,k)=hold1*d23*mo*n1_i(k)*(mn2*d13-mo2*d12)
          beta(2,1,k)=hold1*d13*mo2*n2_i(k)*(mn2*d23-mo*d12)
        enddo

! solve tridiagonal problem
        do k=1,levs
          dp1(k)=1./(prsi(i,k)-prsi(i,k+1))
          parta(k)=dtp*grav(i,k)*dp1(k)/bz
        enddo
        do k=2,levs
          dp1_i(k)=1./(prsl(i,k-1)-prsl(i,k))
          partb_i(k)=.5*(grav(i,k)+grav(i,k-1))/t_i(k)
        enddo
        do k=2,levs
          hold1=parta(k)*partb_i(k)
          hold2=am(i,k-1)*prsl(i,k-1)*dp1_i(k)
          a(1,1,k)=hold1*beta(1,1,k)*(hold2/mo-.5)
          a(1,2,k)=hold1*beta(1,2,k)*(hold2/mo2-.5)
          a(2,1,k)=hold1*beta(2,1,k)*(hold2/mo-.5)
          a(2,2,k)=hold1*beta(2,2,k)*(hold2/mo2-.5)
        enddo
        do k=1,levs-1
          hold1=parta(k)*partb_i(k+1)
          hold2=am(i,k+1)*prsl(i,k+1)*dp1_i(k+1)
          c(1,1,k)=hold1*beta(1,1,k+1)*(hold2/mo+.5)
          c(1,2,k)=hold1*beta(1,2,k+1)*(hold2/mo2+.5)
          c(2,1,k)=hold1*beta(2,1,k+1)*(hold2/mo+.5)
          c(2,2,k)=hold1*beta(2,2,k+1)*(hold2/mo2+.5)
        enddo
        do k=2,levs-1
          hold1=am(i,k)*prsl(i,k)*dp1_i(k+1)
          hold2=am(i,k)*prsl(i,k)*dp1_i(k)
          b(1,1,k)=1.+parta(k)*(partb_i(k+1)*beta(1,1,k+1)*(hold1/mo-.5)    &
                             +partb_i(k)*beta(1,1,k)*(hold2/mo+.5))
          b(2,2,k)=1.+parta(k)*(partb_i(k+1)*beta(2,2,k+1)*(hold1/mo2-.5)   &
                             +partb_i(k)*beta(2,2,k)*(hold2/mo2+.5))
          b(1,2,k)=parta(k)*(partb_i(k+1)*beta(1,2,k+1)*(hold1/mo2-.5)      &
                             +partb_i(k)*beta(1,2,k)*(hold2/mo2+.5))
          b(2,1,k)=parta(k)*(partb_i(k+1)*beta(2,1,k+1)*(hold1/mo-.5)       &
                             +partb_i(k)*beta(2,1,k)*(hold2/mo+.5))
        enddo
        hold1=am(i,1)*prsl(i,1)*dp1_i(2)
        b(1,1,1)=1.+parta(1)*partb_i(2)*beta(1,1,2)*(hold1/mo-.5)
        b(2,2,1)=1.+parta(1)*partb_i(2)*beta(2,2,2)*(hold1/mo2-.5)
        b(1,2,1)=parta(1)*partb_i(2)*beta(1,2,2)*(hold1/mo2-.5)
        b(2,1,1)=parta(1)*partb_i(2)*beta(2,1,2)*(hold1/mo-.5)
        hold2=am(i,levs)*prsl(i,levs)*dp1_i(levs)
        b(1,1,levs)=1.+parta(levs)*partb_i(levs)*beta(1,1,levs)*    &      
        (hold2/mo+.5)
        b(2,2,levs)=1.+parta(levs)*partb_i(levs)*beta(2,2,levs)*    &      
        (hold2/mo2+.5)
        b(1,2,levs)=parta(levs)*partb_i(levs)*beta(1,2,levs)*       &       
        (hold2/mo2+.5)
        b(2,1,levs)=parta(levs)*partb_i(levs)*beta(2,1,levs)*       &       
        (hold2/mo+.5)
        do k=levs,1,-1
          ggg(1,1)=b(2,2,k)-c(2,1,k)*ee(1,2,k+1)-c(2,2,k)*ee(2,2,k+1)
          ggg(2,2)=b(1,1,k)-c(1,1,k)*ee(1,1,k+1)-c(1,2,k)*ee(2,1,k+1)
          ggg(1,2)=-1.*b(1,2,k)+c(1,1,k)*ee(1,2,k+1)+c(1,2,k)*ee(2,2,k+1)
          ggg(2,1)=-1.*b(2,1,k)+c(2,1,k)*ee(1,1,k+1)+c(2,2,k)*ee(2,1,k+1)
          hold1=1./(ggg(1,1)*ggg(2,2)-ggg(1,2)*ggg(2,1))
          ggg=ggg*hold1
          ee(1,1,k)=ggg(1,1)*a(1,1,k)+ggg(1,2)*a(2,1,k)       
          ee(1,2,k)=ggg(1,1)*a(1,2,k)+ggg(1,2)*a(2,2,k)       
          ee(2,1,k)=ggg(2,1)*a(1,1,k)+ggg(2,2)*a(2,1,k)       
          ee(2,2,k)=ggg(2,1)*a(1,2,k)+ggg(2,2)*a(2,2,k)       
          f(1,k)=ggg(1,1)*(qin(i,k,1)+c(1,1,k)*f(1,k+1)   +             &               
           c(1,2,k)*f(2,k+1))+ggg(1,2)*(qin(i,k,2)+c(2,1,k)*f(1,k+1) +  &       
           c(2,2,k)*f(2,k+1))
          
          f(2,k)=ggg(2,1)*(qin(i,k,1)+c(1,1,k)*f(1,k+1)  +             &              
           c(1,2,k)*f(2,k+1))+ggg(2,2)*(qin(i,k,2)+c(2,1,k)*f(1,k+1)+  &       
           c(2,2,k)*f(2,k+1))
        enddo
        do in=1,ntrac_i
          qout(i,1,in)=f(in,1)
          dq(i,1,in)=qout(i,1,in)-qin(i,1,in)
        enddo
        do k=2,levs
          qout(i,k,1)=ee(1,1,k)*qout(i,k-1,1)+ee(1,2,k)*qout(i,k-1,2)+f(1,k) 
                  
          qout(i,k,2)=ee(2,1,k)*qout(i,k-1,1)+ee(2,2,k)*qout(i,k-1,2)+f(2,k)
                  
          do in=1,ntrac_i
            dq(i,k,in)=qout(i,k,in)-qin(i,k,in)
          enddo
        enddo
      enddo !i
!      if ( me == master ) then
!	    print *, ' dq,lev=147 ', maxval(dq(:,147,:)), minval(dq(:,147,:)) 
!      endif
      return
      end subroutine wam_tracer_m

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
  subroutine wam_tracer_c(im,levs,ntrac_i,adt,dtp,jo2,n1,n2,    & 
             n,rho,qin,dq)
!-----------------------------------------------------------------------
! calculate tracer changes caused by chemistry reaction
! 2006 Rashid Akmaev and Fei Wu
!      Tim Fuller-Rowell HOx oxygen sinks
!-----------------------------------------------------------------------
      use physcons, only : rgas=>con_rgas
      use physcons, only : avgd => con_avgd
      use machine,  only : kind_phys
      use wamphys_init_module,  only : amo, amn2, amo2
      use wamphys_init_module,  only : oh => wam_oh, ho2 => wam_ho2, vmr_wam

      implicit none

! Argument
      integer, intent(in) :: im          ! number of long data points
 
      integer, intent(in) :: levs        ! number of pressure levels
      integer, intent(in) :: ntrac_i     ! number of tracer add by IDEA
      real(kind=kind_phys),    intent(in) :: dtp         ! time step in second
      real(kind=kind_phys), intent(in) :: adt(im,levs)   ! input  temp at dt=0

      real(kind=kind_phys), intent(in) :: qin(im,levs,ntrac_i)   ! input tracer

      real(kind=kind_phys), intent(in) :: jo2(im, levs)   ! input photo diss rate
      real(kind=kind_phys), intent(in) :: n1(im,levs)     ! number density of o
      real(kind=kind_phys), intent(in) :: n2(im,levs)     ! number density of o2
      real(kind=kind_phys), intent(in) :: n(im,levs)      ! number density of mixture
      real(kind=kind_phys), intent(in) :: rho(im,levs)    ! density of mixture
      real(kind=kind_phys), intent(out):: dq(im,levs,ntrac_i) ! output
! Local variables
      real(kind=kind_phys), dimension(levs) :: noh, nh, nho2, natom, ndens
      real(kind=kind_phys) :: k1,k2,p1,p2,L1,L2
      real(kind=kind_phys) :: k3, k4, k5
      real(kind=kind_phys) :: mo,mo2,mn2
      integer k,i
      real(kind=kind_phys) :: q1new, q2new

      mo=amo*1.e-3/avgd
      mo2=amo2*1.e-3/avgd
      mn2=amn2*1.e-3/avgd

! O2 +hv =JJ => O + O
! O + O+ M  =k1 => O2+M 
! O + O2    =k2=>  O3
! O +OH  =k3=> O2 +H 
! O + O     =k4=>  O2 ?
! O +HO2  =k5=> O2 +OH 
! coefficents
        k3=4.2e-17
        k4=1.0e-26
        k5=3.5e-17

      do k=1,levs
         do i=1,im
            k1=4.7e-45*(300./adt(i,k))**2
            k2=6.e-46*(300./adt(i,k))**(2.4)

            p1=2.*Jo2(i,k)*n2(i,k) * mo/rho(i,k) !O-production mmr
            L1=2.*k1*n1(i,k)*n(i,k)+k2*n2(i,k)*n(i,k)      & ! O-loss
                +k3*oh(k)+k5*ho2(k)+2.*k4*n1(i,k)
            q1new=(qin(i,k,1)+p1*dtp)/(1.+L1*dtp)

            dq(i,k,1)=q1new-qin(i,k,1)
! conserve oxygen
            dq(i,k,2)= -dq(i,k,1)

        enddo
      enddo

      return
  end subroutine wam_tracer_c

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
  subroutine wam_tracer_eddy(im,levs,ntrac_i,grav,prsi,prsl,   &
          rho,dtp,qin,dayno,dq)

! Calculate major species changes by eddy mixing O, O2, and
!     (indirectly) N2
! October 2017 Rashid Akmaev
!
      use machine,     only: kind_phys
      
      use wamphys_init_module, only : skeddy0, skeddy_semiann, skeddy_ann
      
					 
      implicit none
! Arguments
      integer, intent(in) :: im    ! number of long data points in fields
 
      integer, intent(in) :: levs  ! number of pressure levels
      integer, intent(in) :: ntrac_i ! number of addtl tracers in IDEA
      integer, intent(in) :: dayno ! for semiannual variation
      real(kind=kind_phys), intent(in) :: dtp      ! time step in second
      real(kind=kind_phys), intent(in) :: prsi(im,levs+1) ! interface pressure in Pa
      real(kind=kind_phys), intent(in) :: prsl(im,levs)   ! layer pressure in Pa
      real(kind=kind_phys), intent(in) :: grav(im,levs)   ! (m/s**2)
      real(kind=kind_phys), intent(in) :: rho(im,levs)   ! mass density (kg/m**3)
      real(kind=kind_phys), intent(in) :: qin(im,levs,ntrac_i)   ! input tracers
      real(kind=kind_phys), intent(out):: dq(im,levs,ntrac_i) ! output tracer changes
! Locals
      real(kind=kind_phys) :: alpha(levs+1),beta(levs),qout(ntrac_i)
      real(kind=kind_phys) :: a(levs),b(levs),c(levs)
      real(kind=kind_phys) :: e(levs),d(levs,ntrac_i), temp
      integer i,k

!-----------------------------------------------------------------------
! Calculate Keddy (m**2/s) (move this to init subroutine/module later)
! Keddy parameters: mean, width in scale heights, height of max

      real(kind=kind_phys), parameter:: pi = 3.141592653
! semiannual amp


      real(kind=kind_phys), parameter:: dkeddy = 0.
      real(kind=kind_phys), parameter:: dx = 2.    
      real(kind=kind_phys), parameter:: xmax = 15.
      real(kind=kind_phys) :: keddy(levs+1),x, kedmax
      
!     skeddy0=140., skeddy_semiann=60., skeddy_ann=0.,
!     tkeddy0=280., tkeddy_semiann=0., tkeddy_ann=0., 

      if(dkeddy <= 1e-10) then

! Add semiannual variation
!          keddy(:) = skeddy0 +  skeddy_semiann*(cos(4.*pi*(dayno+9.)/365.))   ! WAM-GSM

        do k=1,levs+1
          x = alog(1e5/prsi(1,k))
	    kedmax =skeddy0 +  skeddy_semiann*(cos(4.*pi*(dayno+9.)/365.)) 
          keddy(k)= kedmax*exp(-((x-xmax)/dx)**2) +.5            
        enddo                    
	 		    
      else
         do k=1,levs+1
! height in scale heights
            x = alog(1e5/prsi(1,k))
            keddy(k)= skeddy0*exp(-((x-xmax)/dx)**2) +.5 
         enddo
      endif
!-----------------------------------------------------------------------
! Boundary conditions
      a(1) = 0.
      c(levs) = 0.

      do i = 1,im
! Auxiliary arrays
! at interfaces
         do k = 2,levs
            alpha(k) = keddy(k)*(.5*(rho(i,k-1)+rho(i,k)))**2*    &
                (.5*(grav(i,k-1)+grav(i,k)))/(prsl(i,k-1) -       & 
                prsl(i,k))
         enddo

! in layers
         do k = 1,levs
            beta(k) = dtp*grav(i,k)/(prsi(i,k) - prsi(i,k+1))
         enddo

! Coefficients a(k), c(k) and b(k)
         do k = 2,levs
            a(k) = beta(k)*alpha(k)
         enddo
         do k = 1,levs-1
            c(k) = beta(k)*alpha(k+1)
         enddo
         do k = 1,levs
            b(k) = 1. + a(k) + c(k)
         enddo

! Solve tridiagonal problem for each tracer
! boundary conditions
         e(levs) = a(levs)/b(levs)
         d(levs,:) = qin(i,levs,:)/b(levs)

! go down, find e(k) and d(k)
         do k = levs-1,1,-1
!vay-2022 
            temp =1./(b(k) - c(k)*e(k+1))	 
            e(k) = a(k)*temp       !/(b(k) - c(k)*e(k+1))
            d(k,:) = (c(k)*d(k+1,:) + qin(i,k,:)) *temp
         enddo
         
! go up, find solution
         qout(:) = d(1,:)
         dq(i,1,:) = qout(:) - qin(i,1,:)
         do k = 2,levs
            qout(:) = e(k)*qout(:) + d(k,:)
            dq(i,k,:) = qout(:) - qin(i,k,:)
         enddo

      enddo
  end subroutine wam_tracer_eddy
!
!      
  subroutine wamphys_dissociation(im,  levs, te, cospass, o_n, o2_n, o3_n, n2_n, &
                dayno, zg, grav,  f107, f107d, jo2rad, jo3rad)
!----------------------------------------------------------------------------
! calculete solar dissociation
!----------------------------------------------------------------------------
      use wamphys_set_merge_rad, only : npsrad 
      use wamphys_init_module,   only : effuv, effeuv

      use wamphys_init_module,   only : amo, amn2, amno, amo2
      
      use physcons, only         : rgas => con_rgas, avgd => con_avgd
      
      use machine, only    : kind_phys
      implicit none

      integer, intent(in) :: im           !number of data piont in te

      integer, intent(in) :: levs         !number of press level
      integer, intent(in) :: dayno        ! calender day
      
      real(kind=kind_phys), intent(in)    :: f107, f107d     
      real(kind=kind_phys), intent(in)    :: te(im,levs)     !temperature
      real(kind=kind_phys), intent(in)    :: cospass(im)     ! cos zenith angle
      real(kind=kind_phys), intent(in)    :: o_n(im,levs)    !number density of O(1/m3)
      real(kind=kind_phys), intent(in)    :: o2_n(im,levs)   !number density of O2
      real(kind=kind_phys), intent(in)    :: o3_n(im,levs)   !number density of O2      
      real(kind=kind_phys), intent(in)    :: n2_n(im,levs)   !number density of N2
      real(kind=kind_phys), intent(in)    :: zg(im,levs)     !layer height (m)
      real(kind=kind_phys), intent(in)    :: grav(im,levs)   ! (m/s2)

! VAY out dissociation_rate2d

      real(kind=kind_phys), intent(out)   :: jo2rad(im, levs)
      real(kind=kind_phys), intent(out)   :: jo3rad(im, levs)
!locals

      real(kind=kind_phys)    :: n2(levs), o(levs),o2(levs), o3(levs)
      real(kind=kind_phys)    :: ho(levs), ho2(levs),hn2(levs)
      real(kind=kind_phys)    :: sheat(levs),sh1(levs),sh2(levs)

      real(kind=kind_phys)    :: shsrc(levs),shsrb(levs),shlya(levs)

      real(kind=kind_phys)    ::  ht(levs)
      real(kind=kind_phys)    ::  jo2_wrk(levs),jo3_wrk(levs)

      real(kind=kind_phys) :: vay_rgas_o   !, vay_rgas_o2,vay_rgas_n2 

      real(kind=kind_phys) :: tg_vay, rodn2
      integer :: i,k
!nullify
      
      vay_rgas_o =  1.e3*rgas/amo

      rodn2 = amo/amn2

      do i=1,im
	  jo2rad(i,1:npsrad-1) =0.0 
	  jo3rad(i,1:npsrad-1) =0.0 	       
        do k=1,levs

          o(k)=o_n(i,k)                             !/m3
          o2(k)=o2_n(i,k)
          o3(k)=o3_n(i,k)	  
          n2(k)=n2_n(i,k)
      
          ht(k)  = zg(i,k)                          !layer height (m)
          tg_vay = te(i,k)/grav(i,k)
          ho(k) =vay_rgas_o*tg_vay                  !m
          ho2(k)=.5*ho(k) 
          hn2(k)= rodn2*ho(k)  
!         ho(k)=1.e3*rgas*t(k)/(amo*grav(i,k))      !m      
        enddo
	
! vay-next
!       call wamphys_sheat_jrates(levs,npsrad,o,o2,o3, n2, ho,ho2,hn2,         &
!           f107,f107d,cospass(i),dayno,ht,  sheat, sh1, shsrc, shsrb, shlya,  &
!	      jo2_wrk, jo3_wrk)

        call wamphys_jrates_shuveuv(levs,npsrad,o,o2,o3, n2, ho,ho2,hn2,         &
           f107,f107d,cospass(i),dayno,ht,  sheat, sh1, sh2,                     &
	     jo2_wrk, jo3_wrk)
	   
	  jo2rad(i,npsrad:levs) = jo2_wrk(npsrad:levs) 
        jo3rad(i,npsrad:levs) = jo3_wrk(npsrad:levs)

            
      enddo        !i-hor index
      
!      print *, 'wamphys_Jo2', maxval(jo3rad(1:im,npsrad:levs)), minval(jo3rad(1:im,npsrad:levs)), npsrad

      return
  end subroutine wamphys_dissociation

