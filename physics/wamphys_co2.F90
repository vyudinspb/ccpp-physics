      subroutine wamrad_co2(im,levs,nlev,ntrac, nto, nto2, nto3, nth2o, & 
                 co2my, grav,cp,adr,adt, dtdt,cosz,dtdth)
!
!
! Apr 06 2012   Henry Juang, initial implement for nems
! Dec 13 2012   Jun Wang     move init step out of column physics
! Feb 13 2012   Jun Wang     move gravity array gg to idea_compistion module
! May 10 2022   Svetlana Karol & CUA FV3WAM-ccpp version and review

!      use wam_co2pro_mod, only: co2my
      use machine,    only: kind_phys
      
      use physcons, only : amo2=>con_amo2, amo3=>con_amo3, amh2o=>con_amw
      
      use wamphys_const, only : amo=>con_amo, amn2=>con_amn2
      use wamphys_init_module, only : prlog, k43
      implicit none
! Argument
      integer, intent(in) :: im  ! number of data points in adt (first dim)

      integer, intent(in) :: levs   ! number of pressure levels
      integer, intent(in) :: nlev   ! number of pressure levels in calculation
      integer, intent(in) :: ntrac  ! number of tracer
      integer, intent(in) :: nto, nto2, nto3, nth2o    
       
      real(kind=kind_phys), intent(in)    :: adr(im,levs,ntrac) ! tracer
      real(kind=kind_phys), intent(in)    :: adt(im,levs)    ! temperature
      real(kind=kind_phys), intent(in)    :: cp(im,levs)    ! J/kg/k
      real(kind=kind_phys), intent(in)    :: grav(im,levs)    ! g (m/s2)
      real(kind=kind_phys), intent(in)    :: co2my(levs)     
      real(kind=kind_phys), intent(in)    :: cosz(im)    !cos solar zenith angle 

      real(kind=kind_phys), intent(out)   :: dtdt(im,levs)    ! cooling rate k/s
      real(kind=kind_phys), intent(out)   :: dtdth(im,levs)    ! heating rate k/s
!
      real(kind=kind_phys) :: pmod(levs),q_n2(im,nlev),ma(im,nlev)      &                   
      ,q_o(im,nlev),q_o2(im,nlev),hold(levs)
      integer i,k,kk
!
! precalling
      dtdth(:,:)=0.
      dtdt(:,:) =0.
!
      do i=1,im
        do k=k43,levs
          kk=k-k43+1
          q_n2(i,kk)=1.-adr(i,k,nto)-adr(i,k,nto2)-adr(i,k,nto3)-adr(i,k,nth2o)
	  
          ma(i,kk)=1./(adr(i,k,nto)/amo+adr(i,k,nto2)/amo2+adr(i,k,nth2o)/amh2o+     &
                  adr(i,k,nto3)/amo3+q_n2(i,kk)/amn2)
		  
          q_o(i,kk)=adr(i,k,nto)*ma(i,kk)/amo
          q_o2(i,kk)=adr(i,k,nto2)*ma(i,kk)/amo2
          q_n2(i,kk)=q_n2(i,kk)*ma(i,kk)/amn2
	  
        enddo
      enddo
!     print*,'www2',im,q_o(1:im1,nlev)
! CO2 cooling
      call wam_co2cc(im, prlog,adt,levs,prlog(k43),     &                  
                dtdt(1,k43),nlev,ma,q_o,q_o2,q_n2, co2my(k43))
! J/kg/s to k/s
      do i=1,im
        do k=k43,levs
          dtdt(i,k)=dtdt(i,k)/cp(i,k)
        enddo
!vay-16          dtdt(i,1:k43-1)=0.
      enddo
! CO2 heating
      do i=1,im
        call wam_qnirc(cosz(i),prlog(k43),co2my,hold(k43),nlev)
        do k=k43,levs
!       dtdth(i,k)=hold(k-k43+1)
        dtdth(i,k)=hold(k)
        enddo
!vay-16        dtdth(i,1:k43-1)=0.
      enddo
      return
      end
