      subroutine wamrad_h2o(im, levs, nlev, nlevc,ntrac,nto, nto2, nto3, nth2o, &
                            grav,cp,adr, adt,dth,cosz,dtc, &
		gh2ort,gh2ovb,dg1rt,dg2rt, dg1vb,dg2vb,gdp,xx,wvmmrc,coeff)			    
			    
!       call wamrad_h2o(im,  levs,nlev_h2o,nlevc_h2o,ntrac,nto1, nto2, nto3, ntqv, &
!                       grav,cp, qgrs, tgrs, dth2oh, cospass, dth2oc) 
! Apr 06 2012  Henry Juang, initial implement for nems
! Dec    2012    Jun Wang,  move init step out of column physics
!
      use machine,    only : kind_phys
!      use physcons,  only : amo3=>con_amo3, amh2o=>con_amw
      use wamphys_const,       only : mmr_min, vmr_nzero
      use wamphys_init_module, only : pmb,prlog, h2ora, k41,k110,k71,k100,k105 
      use wamphys_init_module, only : amo,amo2, amo3, amn2,amh2o 

      implicit none
! Argument
      integer, intent(in) :: im  ! number of data points in adt (first dim)
 
      
      integer, intent(in) :: levs   ! number of pressure levels in GFS
      integer, intent(in) :: nlev    ! number of pressure levels in heating
      integer, intent(in) :: nlevc   ! number of pressure levels in cooling
      integer, intent(in) :: ntrac  ! number of tracer
      
      integer, intent(in) :: nto, nto2, nto3, nth2o    
      
      real(kind=kind_phys), intent(in)    :: adr(im,levs,ntrac) ! tracer
      real(kind=kind_phys), intent(in)    :: adt(im,levs) ! temp (k)
      real(kind=kind_phys), intent(in)    :: grav(im,levs)    ! (m/s2)
      real(kind=kind_phys), intent(in)    :: cp(im,levs)    ! J/kg/k
      real(kind=kind_phys), intent(in)    :: cosz(im)        ! cos zenith angle
      
      
      real(kind=kind_phys),dimension(levs) ::  gh2ort,gh2ovb,dg1rt,dg2rt    
      real(kind=kind_phys),dimension(levs) ::  dg1vb,dg2vb,gdp,xx,wvmmrc,coeff 
      
      real(kind=kind_phys), intent(out)   :: dtc(im,levs)    ! cooling rate k/s
      real(kind=kind_phys), intent(out)   :: dth(im,levs)    ! heating rate k/s
!
      real(kind=kind_phys) :: pmodi(nlev),ggg(nlev), h2ommr(nlev),mu(nlev),rcp(nlev),dthi(nlev)                       
      real(kind=kind_phys) :: adrn2,rate,dx
      
      real(kind=kind_phys) ::  h2ommrc(nlevc),temp(nlevc),qr(nlevc),qv(nlevc),prpa(nlevc)
      integer i,k,k1
      logical :: qvNAN, qrNaN
!
!     print*,'www1',nlev_h2o,nlevc_h2o,k41,k110,k71,k100,k105
!     print*,'www1',h2ora(71),h2ora(150)
!
      dtc(:,:)=0.
      dth(:,:)=0.
      
!
! VAY-2016: don't need to zero "dtc" and "dth" by "special" loops below
!
!
      do k=1,nlev
        pmodi(k)=pmb(k41-1+k)
      enddo
!
! VAY-2016 with NRL H2O_phys no needs to scale H2O > k71  check it later
!
      do i=1,im
           rate=adr(i,k71,nth2o)/h2ora(k71)
          do k=1,nlev
            k1=k41-1+k
              if(k1.le.k71-1) then
                h2ommr(k)=adr(i,k1,nth2o)
              else
                h2ommr(k)=rate*h2ora(k1)     ! h2ora .... global Profile of H2O ???
              endif
	     h2ommr(k)=max(h2ommr(k), mmr_min)  
            adrn2=1.-adr(i,k1,nto)-adr(i,k1,nto2)-adr(i,k1,nto3)-h2ommr(k)
            ggg(k)=grav(i,k1)
            mu(k)=1./(adr(i,k1,nto)/amo+adr(i,k1,nto2)/amo2+    &             
                   h2ommr(k)/amh2o+adr(i,k1,nto3)/amo3+adrn2/amn2)	    
          enddo
         dthi(1:nlev)=0.
! get heating
        call wam_h2ohdc(cosz(i), pmodi, h2ommr, ggg, mu, dthi, nlev)
! k1=k41-1+k
        do k=k41,k110
          dth(i,k)= dthi(k-k41+1)/cp(i,k)
        enddo
!vay don't need          dth(i,1:k41-1)=0.
      enddo
!
! merge to 0. on top in 5-layers by LINEAR extrap-n????
!
      dx=prlog(k105)-prlog(k100)
      do i=1,im
        do k=k100+1,k105-1
           dth(i,k)=dth(i,k)*(prlog(k105)-prlog(k))/dx
        enddo
       do k=k105, levs
           dth(i,k)=0.
       enddo
      enddo
      
! cooling
! cooling idea pressure level 71-150 upward
!

      prpa(1:nlevc)=pmb(k71:levs)
      
      do i=1,im
        rate=adr(i,k71,nth2o)/h2ora(k71)
          do k=1,nlevc
             h2ommrc(k)=rate*h2ora(k+k71-1)
             h2ommrc(k)=max(h2ommrc(k),mmr_min)
             temp(k)=adt(i,k+k71-1)
          enddo
	  qr(1:nlevc) = 0.
	  qv(1:nlevc) = 0.
        call wam_h2occ(temp,pmb(k71:levs),h2ommrc,qr,qv,nlevc,     &
	 gh2ort(k71),gh2ovb(k71),dg1rt(k71),dg2rt(k71), dg1vb(k71),&
	 dg2vb(k71), gdp(k71),xx(k71),wvmmrc(k71),coeff(k71)  )	
	 
	do k=k71, k105	 
	  dtc(i,k)= qr(k+1-k71)+qv(k+1-k71) 
	enddo	

      enddo
      return
      end subroutine wamrad_h2o
