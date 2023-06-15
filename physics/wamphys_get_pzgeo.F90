subroutine wam_get_prsi(levs, im, ak5, bk5, psurf, prsi)
!
! hyb2press_gc.f:          prsi(i,k)  = ak5(k)+bk5(k)*pgr(i)
! 
      use machine ,              only : kind_phys
      implicit none
      integer :: levs, im, ix
      real(kind=kind_phys) :: psurf(im)  , ak5(levs+1), bk5(levs+1)    
      real(kind=kind_phys) :: prsi(im,levs+1)     
      integer i, k, n
        do i=1,im
          do k=1, levs
            prsi(i,k) = ak5(k)*100.+bk5(k)*psurf(i)  ! ak5 in mb
          enddo
          prsi(i,levs+1) = 0.
        enddo
end subroutine wam_get_prsi

subroutine wamphys_get_przgeo(im,levs,ntrac,t,q, &
                   prsi,prki,prsl,prkl,phii,phil,del)
!
      use machine ,              only : kind_phys
      use physcons ,             only : cp => con_cp, nu => con_fvirt
      use physcons ,             only : rd => con_rd, rkap => con_rocp
      use  wamphys_multigases,   only : ri , cpi
      implicit none
!
      integer im, levs, ntrac      
      real(kind=kind_phys) :: prsi(im,levs+1), prki(im,levs+1), phii(im,levs+1)
      real(kind=kind_phys) :: prsl(im,levs), phil(im,levs), prkl(im,levs) 
      real(kind=kind_phys) :: del(im,levs), t(im,levs), q(im,levs,ntrac)     
      real(kind=kind_phys) :: xcp(im,levs), xr(im,levs), kappa(im,levs)
      real(kind=kind_phys) :: tem, dphib, dphit, dphi
      real(kind=kind_phys), parameter :: zero=0.0, p00i=1.0e-5
      real(kind=kind_phys), parameter :: rkapi=1.0/rkap, rkapp1=1.0+rkap
      integer i, k, n
!
      do k=1,levs
        do i=1,im
          del(i,k) = prsi(i,k) - prsi(i,k+1)
        enddo
      enddo

      call wam_get_cpr(im,levs,ntrac, q, xcp, xr)	      	    
!
      do k=1,levs
        do i=1,im
          kappa(i,k) = xr(i,k)/xcp(i,k)
          prsl(i,k)  = (prsi(i,k) + prsi(i,k+1))*0.5
          prkl(i,k)  = (prsl(i,k)*p00i) ** kappa(i,k)
        enddo
      enddo
      do k=2,levs
        do i=1,im
          tem = 0.5 * (kappa(i,k) + kappa(i,k-1))
          prki(i,k-1) = (prsi(i,k)*p00i) ** tem
        enddo
      enddo
      do i=1,im
        prki(i,1) = (prsi(i,1)*p00i) ** kappa(i,1)
      enddo
      k = levs + 1
      if (prsi(1,k) .gt. 0.0) then
        do i=1,im
          prki(i,k) = (prsi(i,k)*p00i) ** kappa(i,levs)
        enddo
      endif
!
      do i=1,im
        phii(i,1)   = 0.0           ! ignoring topography height here
      enddo
      do k=1,levs
        do i=1,im
          tem         = xr(i,k) * t(i,k)
          dphi        = (prsi(i,k) - prsi(i,k+1)) * tem/(prsi(i,k) + prsi(i,k+1))
          phil(i,k)   = phii(i,k) + dphi
          phii(i,k+1) = phil(i,k) + dphi
        enddo
      enddo
      return
end subroutine wamphys_get_przgeo
           
subroutine wam_get_cpr(im,levs,ntrac,q,xcp,xr)
!
      use machine ,      only : kind_phys
      use wamphys_multigases, only : ri , cpi
      implicit none
!
      real(kind=kind_phys), parameter :: zero=0.0
      integer :: im,  levs, ntrac
      real(kind=kind_phys) :: q(im,levs,ntrac)
      real(kind=kind_phys) :: xcp(im,levs),xr(im,levs),sumq(im,levs)
      integer i, k, n
!
      sumq = zero
      xr   = zero
      xcp  = zero
      do n=1,ntrac
        if( ri(n) > 0.0 ) then
          do k=1,levs
            do i=1,im
              xr(i,k)   = xr(i,k)   + q(i,k,n) * ri(n)
              xcp(i,k)  = xcp(i,k)  + q(i,k,n) * cpi(n)
              sumq(i,k) = sumq(i,k) + q(i,k,n)
            enddo
          enddo
        endif
      enddo
      do k=1,levs
        do i=1,im
          xr(i,k)    = (1.-sumq(i,k))*ri(0)  + xr(i,k)
          xcp(i,k)   = (1.-sumq(i,k))*cpi(0) + xcp(i,k)
        enddo
      enddo
!
      return
end subroutine wam_get_cpr
      
subroutine wam_get_rdmulti(im, levs,ntrac,q,xr)
!
      use machine ,            only : kind_phys
      use  wamphys_multigases, only : ri     
      implicit none
!
      real (kind=kind_phys), parameter :: zero=0.0
      integer im,  levs, ntrac
      real(kind=kind_phys) q(im,levs,ntrac)
      real(kind=kind_phys) xr(im,levs),sumq(im,levs)
      integer i, k, n
!
      sumq = zero
      xr   = zero
      do n=1,ntrac
        if( ri(n) > 0.0 ) then
          do k=1,levs
            do i=1,im
              xr(i,k)   = xr(i,k)   + q(i,k,n) * ri(n)
              sumq(i,k) = sumq(i,k) + q(i,k,n)
            enddo
          enddo
        endif
      enddo
      do k=1,levs
        do i=1,im
          xr(i,k)    = (1.-sumq(i,k))*ri(0)  + xr(i,k)
        enddo
      enddo
!
      return
end subroutine wam_get_rdmulti
      
subroutine wam_get_cp(im, levs,ntrac,q, xcp)
!
      use machine ,      only : kind_phys
      use  wamphys_multigases, only : cpi
      implicit none
!
      real (kind=kind_phys), parameter :: zero=0.0
      integer im, levs, ntrac
      real(kind=kind_phys) q(im,levs,ntrac)
      real(kind=kind_phys) xcp(im,levs),sumq(im,levs)
      integer i, k, n
!
      sumq = zero
      xcp  = zero
      do n=1,ntrac
        if( cpi(n) > 0.0 ) then
          do k=1,levs
            do i=1,im
              xcp(i,k)  = xcp(i,k)  + q(i,k,n) * cpi(n)
              sumq(i,k) = sumq(i,k) + q(i,k,n)
            enddo
          enddo
        endif
      enddo
      do k=1,levs
        do i=1,im
          xcp(i,k)   = (1.-sumq(i,k))*cpi(0) + xcp(i,k)
        enddo
      enddo
!
      return
end subroutine wam_get_cp
      
subroutine get_phi_fv3wam(im, levs, xrd, gt0, del_gz, phii, phil)
     use machine ,              only : kind_phys
     implicit none

     ! Interface variables
     integer, intent(in) :: im, levs
  
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: gt0
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: xrd
     real(kind=kind_phys), dimension(:,:),     intent(inout) :: del_gz
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: phii
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: phil
 
     ! Local variables
     integer :: i, k

! SJL: Adjust the height hydrostatically in a way consistent with FV3 discretization
     do i=1,im
        phii(i,1) = 0.0
     enddo
     do k=1,levs
       do i=1,im
         del_gz(i,k) = del_gz(i,k)*gt0(i,k)           !* (1.0 + con_fvirt*max(zero,gq01(i,k)))
         phii(i,k+1) = phii(i,k) + del_gz(i,k)
         phil(i,k)   = 0.5*(phii(i,k) + phii(i,k+1))
       enddo
     enddo

end subroutine get_phi_fv3wam
    
subroutine get_prs_fv3wam(im, levs, phii, prsi, tgrs, xrd, del, del_gz)
     use physcons,       only : con_rd
     use machine ,       only : kind_phys
     implicit none

     ! Interface variables
     integer, intent(in) :: im, levs
 
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: xrd    
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: phii
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: prsi
     real(kind=kind_phys), dimension(:,:),     intent(in)    :: tgrs
     
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: del
     real(kind=kind_phys), dimension(:,:),     intent(out)   :: del_gz
     real(kind=kind_phys) :: rdry
     
     ! Local variables
     integer :: i, k
      
     rdry = 1./con_rd

! SJL: Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
! del_gz is a temp array recording the old info before (t,q) are adjusted
!  dz = dZgeo/T  => Rv T rho = P => dZgeo/dP = -RT/P => dZgeo/T = -RdP/P  
     do k=1,levs
       do i=1,im
         del(i,k) = prsi(i,k) - prsi(i,k+1)
         del_gz(i,k) = (phii(i,k+1) - phii(i,k)) /tgrs(i,k)/(1.0 + xrd(i,k)*rdry)
       enddo
     enddo

end subroutine get_prs_fv3wam
!==================================================================================
!
! the program define the pressure-vertical grid for the major WAM-species (O-O2-N2)
!                        and the gravity accelearation
!            exner, exner_i, kappa_i - arrays are based on Cp/R after the dycore
!              
!==================================================================================   
subroutine wamphys_zgrav(im, levs, ntrac, tgrs, qgrs,             &
                prsl,  prsi, phii, phil, del, oro,                & 
                zgeo, grav, exner, exner_i, kappa, xcp, rdel ) 
		
     use machine ,       only : kind_phys
     use physcons,       only : con_rerth, con_g  !6.3712e+6 meters
     implicit none

     ! Interface variables
     
     integer, intent(in) :: im, levs, ntrac
 
     real(kind=kind_phys), dimension(im),          intent(in)    :: oro    
     real(kind=kind_phys), dimension(im, levs+1),  intent(in)    :: phii
     real(kind=kind_phys), dimension(im, levs+1),  intent(in)    :: prsi
     real(kind=kind_phys), dimension(im,levs),     intent(in)    :: phil
     real(kind=kind_phys), dimension(im,levs),     intent(in)    :: prsl     
     real(kind=kind_phys), dimension(im,levs),     intent(in)    :: del  
        
     real(kind=kind_phys), dimension(im,levs),     intent(in)    :: tgrs
     real(kind=kind_phys), dimension(im,levs,ntrac),intent(in)    :: qgrs 
     
     real(kind=kind_phys), dimension(im,levs),     intent(out)   :: zgeo 
     real(kind=kind_phys), dimension(im,levs),     intent(out)   :: grav
     real(kind=kind_phys), dimension(im,levs),     intent(out)   :: xcp     
     real(kind=kind_phys), dimension(im,levs),     intent(out)   :: exner  
     real(kind=kind_phys), dimension(im,levs+1),   intent(out)   :: exner_i  
     real(kind=kind_phys), dimension(im,levs),     intent(out)   :: rdel  
     real(kind=kind_phys), dimension(im,levs),     intent(out)   :: kappa   
!
!local
!                
     integer                  ::  i, k    
     real(kind=kind_phys)     ::  rzrad, rer, gphil, p00i
     real(kind=kind_phys)     ::  goro(im) 
     real(kind=kind_phys)     ::  xr(im, levs)    
     real(kind=kind_phys)     ::  inv_exner, tem
         
     rer = 1./con_rerth
	   p00i =1.e-5
     do i = 1,im   
	     goro(i) = con_g * oro(i) /(1. +oro(i)*rer)
	   enddo   
     
     do k = 1,levs
       do i = 1,im
      	   gphil = goro(i)+phil(i,k)
      	   zgeo(i,k) = gphil/(con_g - gphil*rer)	   
      	   rzrad = 1./(1.+zgeo(i,k)*rer)	    
           grav(i,k) = con_g *rzrad *rzrad
	         rdel(i,k) = 1./del(i,k)
       enddo
     enddo
	 
!exner, exner_i, kappa_i

     call wam_get_cpr(im,levs, ntrac, qgrs, xcp, xr)	      	    

     do k=1,levs
       do i=1,im
          kappa(i,k) = xr(i,k)/xcp(i,k)
		      inv_exner =(prsl(i,k)*p00i) ** kappa(i,k)
          exner(i,k)  = 1./inv_exner
       enddo
     enddo

     do k=2,levs
       do i=1,im
          tem = 0.5 * (kappa(i,k) + kappa(i,k-1))
		      inv_exner =(prsi(i,k)*p00i) ** tem
          exner_i(i,k) = 1./inv_exner
       enddo
     enddo
	   k=1
     do i=1,im
       exner_i(i,k) = 1./((prsi(i,k)*p00i) ** kappa(i,k))
     enddo
     k = levs + 1
     if (prsi(1,k) .gt. 0.0) then
       do i=1,im
         exner_i(i,k) = 1./((prsi(i,k)*p00i) ** kappa(i,levs))
       enddo
     endif
   end subroutine wamphys_zgrav	
   	
