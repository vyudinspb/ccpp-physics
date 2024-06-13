module wamphys_math_interp
       use machine , only : kind_phys
!
! subroutine interpol_quad        ionospheric-data
! subroutine idea_comp150_interp  H2o/O3 profiles above 70 km => FV3WAM-grid 
! subroutine interp_field
!
      integer, parameter :: iulog = 6
      contains
!============================== ALL-diverse WAM-interpolation subroutines are below
! linear, quad, spline and special with data-grids and values
!=======================================================      
  SUBROUTINE WAM_POLINT(XA,YA,N,X,Y,DY)
!  Polynomial interpolation from "Numerical Recipes" in wamphys_h2oc.F90
      use machine,                only: kind_phys
      implicit none
      integer :: N
      integer, PARAMETER :: NMAX=10 
      real(kind=kind_phys) :: XA(N),YA(N),C(NMAX),D(NMAX)
      real(kind=kind_phys) :: DY, DIF, DIFT , DEN, W, HP, HO, X, Y
      integer :: NS, M , I          
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
!!!compiler warning          IF(DEN.EQ.0.) PAUSE
          IF(DEN.EQ.0.) STOP 'DEN.EQ.0. in POLINT'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE

      END subroutine wam_polint

!***********************************************************************
  subroutine wam_lint1(x1,y1,n1,yleft,yright,x,y)
! wam_co2hc.F90
! A very simple linear interpolation of y1(x1) into y(x)
! ***x1 is assumed in ascending order***
!
      implicit none
! Array dimension
      integer,intent(in):: n1
      real(kind=kind_phys),intent(in):: x1(n1),y1(n1),yleft,yright,x
      real(kind=kind_phys),intent(out):: y
! Work memory
      integer:: i
      real(kind=kind_phys):: dx

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(x.lt.x1(1)) then
         y=yleft
      elseif(x.gt.x1(n1)) then
         y=yright
      else
         do i=2,n1
            dx=x-x1(i)
            if(dx.le.0.) then
               y=y1(i)+(y1(i)-y1(i-1))*dx/(x1(i)-x1(i-1))
               return
            endif
         enddo
      endif
  end subroutine wam_lint1    

!***********************************************************************

  subroutine wam_splin1(x1,y1,x2,y2,n1,n2)
! A simple routine to interpolate y1[x1(n1)] to y2[x2(n2)] using cubic
! spline.
! Both x1 and x2 are assumed to be ordered in THE SAME, ascending or
! descending, order.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Oct 2006: Rashid Akmaev

      use machine,                only: kind_phys
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine arguments
! INPUT


      integer,intent(in):: n1,n2
      real(kind=kind_phys),intent(in):: x1(n1),y1(n1),x2(n2)
! OUTPUT
      real(kind=kind_phys),intent(out):: y2(n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Internal parameters, work space
!
      real(kind=kind_phys),parameter:: one_third=1./3.
      integer:: i,k,l,nvs
      real(kind=kind_phys):: a(n1),dx,dxmh,dy(n1),e(n1),f(n1),g,h(n1),wx1(n1),wx2(n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Initialize output
!
      y2(:)=0.
!
! Simple check of argument order
!
      if(x1(1).lt.x1(n1)) then
         wx1(:)=x1(:)
         wx2(:)=x2(:)
      else
!
! Reverse x1 and x2 (changing sign seems easier than changing order)
!
         wx1(:)=-x1(:)
         wx2(:)=-x2(:)
      endif
!
! Prepare spline coefficients (Note: they depend on y and so have to
! be recalculated every time)

      nvs=n1-1
      do k=1,nvs
         h(k)=wx1(k+1)-wx1(k)
         dy(k)=(y1(k+1)-y1(k))/h(k)
      enddo
      a(1)=0.
      a(n1)=0.
      e(n1)=0.
      f(n1)=0.
      do k=nvs,2,-1
         g=1./(h(k)*e(k+1)+2.*(h(k-1)+h(k)))
         e(k)=-g*h(k-1)
         f(k)=g*(3.*(dy(k)-dy(k-1))-h(k)*f(k+1))
      enddo
      do k=2,nvs
         a(k)=e(k)*a(k-1)+f(k)
      enddo
!
! Calculate spline values
!
      l=1
      do i=1,n2
         do k=l,nvs
            dx=wx2(i)-wx1(k)
            dxmh=dx-h(k)
            l=k
            if(dxmh.le.0.) exit
         enddo
         g=dx/h(l)
         y2(i)=y1(l)+dx*(dy(l)+one_third*dxmh*(a(l)*(2.-g)+    &         
             a(l+1)*(1.+g)))
      enddo

  end subroutine wam_splin1

!***********************************************************************

  subroutine wam_splin2(x1,y1,x2,y2,n1,n2,jm,km)

! A simple routine to interpolate km arrays y1[x1(n1)], specified on
! the same grid x1(n1), to km arrays y2[x2(n2)] on the same grid x2(n2)
! using cubic spline, where km<=jm and jm is the first dimension of 
! arrays y1 and y2.
! Both grids x1 and x2 are assumed to be ordered in THE SAME, ascending
! or descending, order.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Dec 2006: Rashid Akmaev
! Made from splin1

      use machine,                only: kind_phys
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine arguments
! INPUT
      integer,intent(in):: jm,km,n1,n2
      real(kind=kind_phys),intent(in):: x1(n1),y1(jm,n1),x2(n2)
! OUTPUT
      real(kind=kind_phys),intent(out):: y2(jm,n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Internal parameters, work space
!
      real(kind=kind_phys),parameter:: one_third=1./3.
      integer:: i,j,k,l,nvs
      real(kind=kind_phys):: a(km,n1),dx,dxmh,dy(km,n1),e(n1),f(km,n1),g,g2,h(n1), &
            wx1(n1),wx2(n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Initialize output
!
      y2(:,:)=0.
!
! Simple check of argument order
!
      if(x1(1).lt.x1(n1)) then
         wx1(:)=x1(:)
         wx2(:)=x2(:)
      else
!
! Reverse x1 and x2 (changing sign seems easier than changing order)
!
         wx1(:)=-x1(:)
         wx2(:)=-x2(:)
      endif
!
! Prepare spline coefficients
!
      nvs=n1-1
      do k=1,nvs
         h(k)=wx1(k+1)-wx1(k)
         do j=1,km
            dy(j,k)=(y1(j,k+1)-y1(j,k))/h(k)
         enddo
      enddo
      e(n1)=0.
      do j=1,km
         a(j,1)=0.
         a(j,n1)=0.
         f(j,n1)=0.
      enddo
! 
! Calculate e and f coefficients
!
      do k=nvs,2,-1
         g=1./(h(k)*e(k+1)+2.*(h(k-1)+h(k)))
         e(k)=-g*h(k-1)
         do j=1,km
            f(j,k)=g*(3.*(dy(j,k)-dy(j,k-1))-h(k)*f(j,k+1))
         enddo
      enddo
! 
! Calculate a coefficients
!
      do k=2,nvs
         do j=1,km
            a(j,k)=e(k)*a(j,k-1)+f(j,k)
         enddo
      enddo
!
! Calculate spline values
!
      l=1
      do i=1,n2
         do k=l,nvs
            dx=wx2(i)-wx1(k)
            dxmh=dx-h(k)
            l=k
            if(dxmh.le.0.) exit
         enddo
         dxmh=one_third*dxmh
         g=1.+dx/h(l)
         g2=3.-g
         do j=1,km
            y2(j,i)=y1(j,l)+dx*(dy(j,l)+dxmh*(a(j,l)*g2+a(j,l+1)*g))
         enddo
      enddo

  end subroutine wam_splin2      
!===============================================================         
  subroutine interpol_quad(v,x,u,p)
!  WEMER-2005
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
        use machine , only : kind_phys
        implicit none
! Args:
        real(kind=kind_phys),intent(in) :: v(:),x(:),u(:)
        real(kind=kind_phys),intent(out) :: p(:)
! Local:
        integer :: nv,nx,nu,i,ix
        real(kind=kind_phys) :: x0,x1,x2

        nv = size(v)
        nx = size(x)
        nu = size(u)
        if (nx /= nv) then
          write(iulog,"('>>> interpol_quad: nx /= nv: nx=',i4,' nv=',i4)") nx,nv
          p(:) = 0.
          return
        endif
        do i=1,nu
          ix = value_locate(x,u(i))
          if (ix <= 1.or.ix >= nx) then
!           write(iulog,"('>>> interpol_quad: ix out of range: nu=',i4,' ix=',i4)") nu,ix
            p(i) = 0.
            cycle
          endif
          x1 = x(ix)
          x0 = x(ix-1)
          x2 = x(ix+1)
          if(x0.eq.0..and.x1.eq.0.)then
           p(i) =0.0
          else
           p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2) / ((x0-x1) * (x0-x2)) + &
                  v(ix)   * (u(i)-x0) * (u(i)-x2) / ((x1-x0) * (x1-x2)) + &
                  v(ix+1) * (u(i)-x0) * (u(i)-x1) / ((x2-x0) * (x2-x1))
          endif 
        enddo
!       write(iulog,"('interpol_quad: nu=',i4,' p=',/,(1pe12.4)") nu,p

  end subroutine interpol_quad
       
  subroutine idea_comp150_interp(ain,zin, np, zout, aout, mp)
!      
! call idea_comp150_interp(h2ora150,prlog150(71:npn), np, prlog(k71:levs), h2ora(k71:levs), mp)
! old idea_composition  
      use machine , only : kind_phys
      implicit none
      real(kind=kind_phys) :: ain(np), zin(np)
      real(kind=kind_phys) :: aout(np), zout(np)
      integer              :: np, mp     
       
      real(kind=kind_phys) :: zk,     dz      
     
      integer              :: kref,k,i
      
      aout(1)  = ain(1)
      aout(mp) = ain(np)
      
      do k=1,mp
         kref=0
	   zk = zout(k)
         do i=1,np-1
           if(zk.ge.zin(i).and.zk.le.zin(i+1)) then
             kref=i
             dz=(zk-zin(i))/(zin(i+1)-zin(i))
            endif
         enddo
         if(kref.ne.0) aout(k)=dz*ain(kref+1)+(1.-dz)*ain(kref)
         if(zk > zin(np)) aout(k) = ain(np)
         if(zk < zin(1))  aout(k) = ain(1)	 
      enddo
      return
  end subroutine idea_comp150_interp
       
  subroutine z15toz(ain,levs, Z, aout,down)
      
      use machine , only : kind_phys      
! interpolate 15 pressure levels (from Tim's grid) to
! idea pressure grid Z(levs)
      implicit none
      integer, parameter  :: np=15                     !number of pressure levels of input
      integer, intent(in) :: levs                      !number of pressure levels of output 
      real(kind=kind_phys),    intent(in) :: Z(levs)   !model -log(pressure) grid
      real(kind=kind_phys),    intent(in) :: ain(np)   !input field in 15 pressure grid
      real,    intent(in) :: down      !field value below 1.0376Pa => surface
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real(kind=kind_phys) :: p15(np),z15(np),dz
      integer kref,k,i

      do k=1,np
        p15(k)=1.0376*exp(1.-k)
        z15(k)=-log(p15(k))
      enddo
!
! grids
! interpolation
!
      do k=1,levs
        do i=1,np-1
          if(z(k).ge.z15(i).and.z(k).le.z15(i+1)) then
             dz=(z(k)-z15(i))/(z15(i+1)-z15(i))*(ain(i+1)-ain(i))
             aout(k)=ain(i) +dz
          endif
        enddo   
        if(z(k).lt.z15(1))  aout(k)=down           ! zero
        if(z(k).gt.z15(15)) aout(k)=ain(15)        ! constant-value    
      enddo
      return
  end subroutine z15toz
      
  subroutine z63toz(ain, levs, Z, aout,down)
! interpolate 63 pressure levels (from Tim's grid) to
! idea pressure grid Z(levs)
      use machine , only : kind_phys
 
      implicit none
      integer, parameter  :: np=63     !number of pressure levels of input
      integer, intent(in) :: levs      !number of pressure levels of output 
      real(kind=kind_phys),    intent(in) :: Z(levs)   ! model grid
      real(kind=kind_phys),    intent(in) :: ain(np)   !input field in 63 pressure grid
      real(kind=kind_phys),    intent(in) :: down      !field value under 6.9Pa
      real(kind=kind_phys),    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real(kind=kind_phys) p63(np),z63(np),dz
      integer kref,k,i
!
      DATA p63/6.90775528,  6.57442194,            &
          6.24108859,  5.90775525,  5.57442191,    &
          5.24108856,  4.90775522,  4.57442188,    &
          4.24108853,  3.90775519,  3.57442185,    &
          3.2410885,   2.90775516,  2.57442182,    &
          2.24108847,  1.90775513,  1.57442179,    &
          1.24108844,  0.9077551,   0.574421757,   &
          0.241088414, -0.0922449296,-0.425578273, &
        -0.758911616,-1.09224496,  -1.4255783,     &
        -1.75891165, -2.09224499,  -2.42557833,    &
        -2.75891168, -3.09224502,  -3.42557836,    &
        -3.75891171, -4.09224505,  -4.42557839,    &
        -4.75891174, -5.09224508,  -5.42557842,    &
        -5.75891177, -6.09224511,  -6.42557845,    &
        -6.75891179, -7.09224514,  -7.42557848,    &
        -7.75891182, -8.09224517,  -8.42557851,    &
        -8.75891185, -9.0922452,   -9.42557854,    &
        -9.75891188, -10.0922452,  -10.4255786,    &
        -10.7589119, -11.0922453,  -11.4255786,    &
        -11.7589119, -12.0922453,  -12.4255786,    &
        -12.758912,  -13.0922453,  -13.4255787,    &
        -13.758912/
      do k=1,np
        z63(k)=-p63(k)
      enddo
      do k=1,levs
        do i=1,np-1
          if(z(k).ge.z63(i).and.z(k).le.z63(i+1)) then
              dz=(z(k)-z63(i))/(z63(i+1)-z63(i))*(ain(i+1)-ain(i))
              aout(k)=dz+ain(i)
          endif
        enddo
        if(z(k).lt.z63(1))  aout(k)=down
        if(z(k).gt.z63(63)) aout(k)=ain(63)     
      enddo
      return
  end subroutine z63toz
      
  subroutine interpol_wamz( nin, xin, yin, nout, xout, yout, kup, kdw )
      use machine , only : kind_phys
      implicit none
!-----------------------------------------------------------------------
!       ... linear interpolation in vertical
!           does not extrapolate, but repeats edge values      implicit none
!-----------------------------------------------------------------------
      integer,  intent(in) :: nin, nout
      real(kind=kind_phys), intent(in)    :: xin(nin)   ! reverse 150 => 100
      real(kind=kind_phys), intent(in)    :: yin(nin)
      real(kind=kind_phys), intent(in)    :: xout(nout)
      real(kind=kind_phys), intent(out)   :: yout(nout)
      integer,intent(out) ::  kup      ! zout(kout) > zin(1) = 150 km. zin(nin)=100.
      integer,intent(out) ::  kdw         
      real(kind=kind_phys):: dxin      ! top => bot dxin < 0; bot => top dxin > 0
!-----------------------------------------------------------------------
!       ... local variables 
!-----------------------------------------------------------------------
      integer :: i, j
      kup = nout                       ! for extrap-n up
      kdw = 1
      dxin  = xin(2)-xin(1)
      IF (dxin < 0) THEN               ! snoe-no
       do j = 1,nout
        if( xout(j) .ge. xin(1) ) then
           yout(j) = yin(1)
           kup = min(j, kup)
        else if( xout(j) .le. xin(nin) ) then 
           yout(j) = yin(nin)
 	    kdw = max(j, kdw)
        else
          do i = 1, nin-1
             if ((xout(j) >= xin(i+1)) .and. (xout(j) < xin(i)) ) &
             yout(j) =  yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / (xin(i+1) - xin(i))
          end do
        end if
       end do
! dxin < 0
      ELSE
! dxin > 0     
       do j = 1,nout
        if( xout(j) .le. xin(1) )   then
           yout(j) = yin(1)
	     kdw = max(j, kdw)
        else if( xout(j) .ge. xin(nin) ) then
            yout(j) = yin(nin)
            kup = min(j, kup)
        else
         do i = 1, nin-1
            if ((xout(j) .gt. xin(i)) .and. (xout(j) .lt. xin(i+1) )) &
             yout(j) = yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / (xin(i+1) - xin(i))       
         end do
        end if
       end do
      ENDIF
  end subroutine interpol_wamz

!            interpol_wamz_down(nz_63, z63, SRBEFF63, levs, zlog, SRBEFF, 1.0)      
  subroutine interpol_wamz_down( nin,  xin, yin,      nout, xout, yout,  down )
!-----------------------------------------------------------------------
!       ... linear interpolation in vertical
!           does not extrapolate, but repeats edge values
!-----------------------------------------------------------------------
      use machine , only : kind_phys
      implicit none
!-----------------------------------------------------------------------
      integer,  intent(in) :: nin, nout
      real(kind=kind_phys), intent(in)    :: xin(nin)   ! reverse 150 => 100
      real(kind=kind_phys), intent(in)    :: yin(nin)
      real(kind=kind_phys), intent(in)    :: xout(nout)
      real(kind=kind_phys), intent(in)    :: down      ! zeroes below
      real(kind=kind_phys), intent(out)   :: yout(nout)
      integer             ::  kup      ! zout(kout) > zin(1) = 150 km. zin(nin)=100.
      real(kind=kind_phys):: dxin      ! top => bot dxin < 0; bot => top dxin > 0
!-----------------------------------------------------------------------
!       ... local variables 
!-----------------------------------------------------------------------
      integer :: i, j
      kup = nout                       ! for extrap-n up

      dxin  = xin(2)-xin(1)
      IF (dxin < 0) THEN               ! snoe-no
      do j = 1,nout
       if( xout(j) .gt. xin(1) ) then
          yout(j) = yin(1)
          kup = min(j, kup)
       else   if( xout(j) .lt. xin(nin) ) then 
          yout(j) = down
       else
         do i = 1, nin-1
           if ((xout(j) >= xin(i+1)) .and. (xout(j) <= xin(i)) ) &
           yout(j) =  yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / (xin(i+1) - xin(i))
         end do
       end if
      end do
! dxin < 0
      ELSE
! dxin > 0
       do j = 1,nout
        if( xout(j) .lt. xin(1) )   then
           yout(j) = down
        else  if( xout(j) .gt. xin(nin) ) then
            yout(j) = yin(nin)
            kup = min(j, kup)
        else
         do i = 1, nin-1
           if ((xout(j) .ge. xin(i)) .and. (xout(j) .le. xin(i+1) )) &
           yout(j) = yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / (xin(i+1) - xin(i))        
         end do
        end if
       end do
      ENDIF
  end subroutine interpol_wamz_down

  subroutine z65toz(np, levs, ain, prin, aout, down)
      
! interpolate 65 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)

      use machine , only : kind_phys
      implicit none
      
      integer, intent(in) :: np  !number of pressure levels of input
      integer, intent(in) :: levs  !number of pressure levels of output 
      real(kind=kind_phys),    intent(in) :: ain(np)   !input field in 65 pressure grid
      real(kind=kind_phys),    intent(in) :: down      !field value under 5.2285Pa
      real(kind=kind_phys),    intent(in) :: prin(levs)!   
         
      real(kind=kind_phys),    intent(out):: aout(levs)!output in levs pressure grid
!local variable

      real p65(np),z65(np),z(levs),dz
      integer kref,k,i

      do k=1,np
        p65(k)=5.2285*exp((1.-k)*.25)
        z65(k)=-log(p65(k))
      enddo

      do k=1,levs
        z(k)= prin(k)  !-log(pr(k)*100.)  ! mb or Prlog in Pa
      enddo

      do k=1,levs
        kref=0
        do i=1,np-1
          if(z(k).ge.z65(i).and.z(k).le.z65(i+1)) then
            kref=i
            dz=(z(k)-z65(i))/(z65(i+1)-z65(i))
          endif
        enddo
        if(kref.ne.0) then
          aout(k)=dz*ain(kref+1)+(1.-dz)*ain(kref)
        elseif(z(k).lt.z65(1)) then
          aout(k)=down
        elseif(z(k).gt.z65(65)) then
          aout(k)=ain(np)
        endif
      enddo
      return
  end subroutine z65toz

!---------------------------------------------------------------------
!***********************************************************************
! File splin.f
!***********************************************************************
! December 2006: created by Rashid Akmaev

! Simple spline interpolation subroutines optimized for various numbers
!     of input/output arrays.

! Contains
!      subroutine splin1(x1,y1,x2,y2,n1,n2)
! for y1[x1(n1)] -> y2[x2(n2)]
!      subroutine splin2(x1,y1,x2,y2,n1,n2,jm,km)
! for km transforms y1[km,x1(n1)] -> y2[km,x2(n2)] on the same grids
! x1 and x2, and jm - first dimension of y1 and y2

!***********************************************************************

  subroutine splin1(x1,y1,x2,y2,n1,n2)
      use machine , only : kind_phys
! A simple routine to interpolate y1[x1(n1)] to y2[x2(n2)] using cubic
! spline.
! Both x1 and x2 are assumed to be ordered in THE SAME, ascending or
! descending, order.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Oct 2006: Rashid Akmaev

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine arguments
! INPUT
      integer,intent(in):: n1,n2
      real(kind=kind_phys),intent(in):: x1(n1),y1(n1),x2(n2)
! OUTPUT
      real(kind=kind_phys),intent(out):: y2(n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Internal parameters, work space
!
      real(kind=kind_phys),parameter:: one_third=1./3.
      integer:: i,k,l,nvs
      real(kind=kind_phys):: a(n1),dx,dxmh,dy(n1),e(n1),f(n1),g,h(n1),wx1(n1),wx2(n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Initialize output
!
      y2(:)=0.
!
! Simple check of argument order
!
      if(x1(1).lt.x1(n1)) then
         wx1(:)=x1(:)
         wx2(:)=x2(:)
      else
      ! Reverse x1 and x2 (changing sign seems easier than changing order)
!
         wx1(:)=-x1(:)
         wx2(:)=-x2(:)
      endif
!
! Prepare spline coefficients (Note: they depend on y and so have to
! be recalculated every time)

      nvs=n1-1
      do k=1,nvs
         h(k)=wx1(k+1)-wx1(k)
         dy(k)=(y1(k+1)-y1(k))/h(k)
      enddo
      a(1)=0.
      a(n1)=0.
      e(n1)=0.
      f(n1)=0.
      do k=nvs,2,-1
         g=1./(h(k)*e(k+1)+2.*(h(k-1)+h(k)))
         e(k)=-g*h(k-1)
         f(k)=g*(3.*(dy(k)-dy(k-1))-h(k)*f(k+1))
      enddo
      do k=2,nvs
         a(k)=e(k)*a(k-1)+f(k)
      enddo
!
! Calculate spline values
!
      l=1
      do i=1,n2
         do k=l,nvs
            dx=wx2(i)-wx1(k)
            dxmh=dx-h(k)
            l=k
            if(dxmh.le.0.) exit
         enddo
         g=dx/h(l)
         y2(i)=y1(l)+dx*(dy(l)+one_third*dxmh*(a(l)*(2.-g)+a(l+1)*(1.+g)))     
      enddo

  end subroutine splin1
!***********************************************************************
  subroutine splin2(x1,y1,x2,y2,n1,n2,jm,km)
        use machine , only : kind_phys
! A simple routine to interpolate km arrays y1[x1(n1)], specified on
! the same grid x1(n1), to km arrays y2[x2(n2)] on the same grid x2(n2)
! using cubic spline, where km<=jm and jm is the first dimension of 
! arrays y1 and y2.
! Both grids x1 and x2 are assumed to be ordered in THE SAME, ascending
! or descending, order.
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine arguments
! INPUT
      integer,intent(in):: jm,km,n1,n2
      real(kind=kind_phys),intent(in):: x1(n1),y1(jm,n1),x2(n2)
! OUTPUT
      real(kind=kind_phys),intent(out):: y2(jm,n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Internal parameters, work space
!
      real(kind=kind_phys),parameter:: one_third=1./3.
      integer:: i,j,k,l,nvs
      real(kind=kind_phys):: a(km,n1),dx,dxmh,dy(km,n1),e(n1),f(km,n1),g,g2,h(n1), &
                             wx1(n1),wx2(n2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Initialize output
!
      y2(:,:)=0.
!
! Simple check of argument order
!
      if(x1(1).lt.x1(n1)) then
         wx1(:)=x1(:)
         wx2(:)=x2(:)
      else
!
! Reverse x1 and x2 (changing sign seems easier than changing order)
!
         wx1(:)=-x1(:)
         wx2(:)=-x2(:)
      endif
!
! Prepare spline coefficients
!
      nvs=n1-1
      do k=1,nvs
         h(k)=wx1(k+1)-wx1(k)
         do j=1,km
            dy(j,k)=(y1(j,k+1)-y1(j,k))/h(k)
         enddo
	 
      enddo
      e(n1)=0.
      do j=1,km
         a(j,1)=0.
         a(j,n1)=0.
         f(j,n1)=0.
      enddo
! 
! Calculate e and f coefficients
!
      do k=nvs,2,-1
         g=1./(h(k)*e(k+1)+2.*(h(k-1)+h(k)))
         e(k)=-g*h(k-1)
         do j=1,km
            f(j,k)=g*(3.*(dy(j,k)-dy(j,k-1))-h(k)*f(j,k+1))
         enddo
      enddo
! 
! Calculate a coefficients
!
      do k=2,nvs
         do j=1,km
            a(j,k)=e(k)*a(j,k-1)+f(j,k)
         enddo
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Calculate spline values
!
      l=1
      do i=1,n2
         do k=l,nvs
            dx=wx2(i)-wx1(k)
            dxmh=dx-h(k)
            l=k
            if(dxmh.le.0.) exit
         enddo
         dxmh=one_third*dxmh
         g=1.+dx/h(l)
         g2=3.-g
         do j=1,km
            y2(j,i)=y1(j,l)+dx*(dy(j,l)+dxmh*(a(j,l)*g2+a(j,l+1)*g))
         enddo
      enddo

      end subroutine splin2

! special interpolation for pre-specified [20,91] data:cormag, btot, dipang, glat, glon, nxmag,nymag     
  subroutine interp2_ionfield(im, rlat,rlon,cormago,btoto,dipango,cormag, btot, dipang, glon, glat)
! VAY DANGER !!!!!!special interpolation 
! interp works only for given [20,91]  with 
!         j=46 center + fixed ddlat=180/90? and ddlon=360/20?
!
!
!      USE IDEA_ION_INPUT, only :
!     & cormag, btot, dipang, glat, glon, nxmag,nymag
!
      use machine , only : kind_phys
      implicit none
      
      integer,intent(in)  :: im            ! number of longitude
      real(kind=kind_phys),   intent(in)  :: rlat(im)      ! latitude (rad)
      real(kind=kind_phys),   intent(in)  :: rlon(im)      ! longitude (rad)
      
      real(kind=kind_phys),   intent(in)  :: cormag(20,91),btot(20,91),dipang(20,91),glat(91),glon(20)
      
      real(kind=kind_phys),   intent(out) :: cormago(im),btoto(im),dipango(im) ! field value 

      real(kind=kind_phys) :: dll,dl,ddlat,ddlon,a1,a2,b1,b2,aa,bb
      integer i,iref, iref1, jref,jref1
      integer ::  jcen, ixdim, jydim 
!
! lat lon interval  << needs to  be rewritten >>
!
      ddlat= 3.4906585033333331E-002
      ddlon= 0.3141592653000000
      jcen = 46 
      
      ixdim = 20
      jydim = 91    
! 
      do i=1,im
! get latitude index
        iref=int(rlat(i)/ddlat)+  jcen
        dl=(rlat(i)-glat(iref))/ddlat
	  iref1 = iref+1
	  if(iref1.gt.jydim) iref1=iref
! get longitude index
        jref=int(rlon(i)/ddlon)+1
        jref1=jref+1
        if(jref1.gt.ixdim) jref1=jref1 - ixdim
        dll=(rlon(i)-glon(jref))/ddlon
!	
        a1=cormag(jref,iref)
        a2=cormag(jref1,iref)
        b1=cormag(jref,iref1)
        b2=cormag(jref1,iref1)
	
        aa=(1.-dll)*a1+dll*a2
        bb=(1.-dll)*b1+dll*b2
	
        cormago(i)=(1.-dl)*aa+dl*bb
!
        a1=btot(jref,iref)
        a2=btot(jref1,iref)
        b1=btot(jref,iref1)
        b2=btot(jref1,iref1)
        aa=(1.-dll)*a1+dll*a2
        bb=(1.-dll)*b1+dll*b2
	
        btoto(i)=(1.-dl)*aa+dl*bb
!
        a1=dipang(jref,iref)
        a2=dipang(jref1,iref)
        b1=dipang(jref,iref1)
        b2=dipang(jref1,iref1)
        aa=(1.-dll)*a1+dll*a2
        bb=(1.-dll)*b1+dll*b2
	
        dipango(i)=(1.-dl)*aa+dl*bb
      enddo
      
      btoto =btoto*1.e-9
      
      return
  end subroutine interp2_ionfield
     
  subroutine spole_ion(im, rlat, phir,utsec, sda, phimr, essa, cmorg)
      use machine     , only : kind_phys
      implicit none
      real(kind=kind_phys), parameter    :: pi=3.141592653,dtr=pi/180., pid2 = .5*pi
      
      integer,intent(in)                 :: im          ! number of longitude 
      real(kind=kind_phys),   intent(in) :: rlat(im)    !geo latitude (rad)
      real(kind=kind_phys),   intent(in) :: phir(im)    !geo longitude (rad)
      real(kind=kind_phys),   intent(in) :: utsec       !ut second
      real(kind=kind_phys),   intent(in) :: sda         !solar declination angle (rad)
      real(kind=kind_phys),   intent(out):: phimr(im)   !maglongitude (rad)        
      real(kind=kind_phys),   intent(out):: essa(im)    !magnetic local time    
      real(kind=kind_phys),   intent(out):: cmorg(im)          
! local variables
      real(kind=kind_phys) :: th,th1,phi1,sinth,sinth1,costh1,sinph1,cosph1,ac1,bc1,cc1   
      real(kind=kind_phys) :: ac2,bc2,cc2,phim,ssp,sspr,csda,as1,bs1,cs1,as2,bs2,cs2,gml, cmag
      integer i

      do i=1,im
      th=pid2-rlat(i)
!
! set pole coord. for each hemis.
!
      if (rlat(i).ge.0.0) then
       th1=9.25*dtr
       phi1=-78.0*dtr
      else   
       th1=16.32*dtr
       phi1=-54.0*dtr
      end if
!
      sinth=sin(th)
      sinth1=sin(th1)
      costh1=cos(th1)
      sinph1=sin(phi1)
      cosph1=cos(phi1)

!     do i=1,im
      ac1=sinth*cos(phir(i))
      bc1=sinth*sin(phir(i))
      cc1=cos(th)
      ac2=ac1*costh1*cosph1+bc1*costh1*sinph1-cc1*sinth1
      if((abs(ac2)).lt.0.001)ac2=0.001
      bc2=-ac1*sinph1+bc1*cosph1
      cc2=ac1*sinth1*cosph1+bc1*sinth1*sinph1+cc1*costh1
      cmorg(i)=acos(cc2)
      phimr(i)=atan2(bc2,ac2)
      phim=phimr(i)/dtr
!     ssp=360.-utsec/240.
      ssp=180.-utsec/240.
      sspr=ssp*dtr
      csda=pid2-sda
      as1=cos(sspr)*sin(csda)
      bs1=sin(sspr)*sin(csda)
      cs1=cos(csda)
      as2=as1*costh1*cosph1+bs1*costh1*sinph1-cs1*sinth1
      if((abs(as2)).lt.0.001)as2=0.001
      bs2=-as1*sinph1+bs1*cosph1
      cs2=as1*sinth1*cosph1+bs1*sinth1*sinph1+cs1*costh1
      gml=atan2(bs2,as2)/dtr
      essa(i)=phim-gml
      enddo
      
      return
  end subroutine spole_ion
           
     
!=====================================================================      
!      efield.f: JULDAY_WAM,  SVD + Spherical harmonics + Adjust
!=====================================================================
  SUBROUTINE ADJUST(ANGLE)
        use machine , only : kind_phys
!-----------------------------------------------------------------------
!	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
!
        implicit none 
!
!------------------------------Arguments--------------------------------
!
        real(kind=kind_phys) ::  angle

        IF(ANGLE.LT.0.)  ANGLE=ANGLE+360.
        IF(ANGLE.GE.360.)ANGLE=ANGLE-360.
        
	RETURN
  END SUBROUTINE ADJUST

  INTEGER FUNCTION JULDAY_WAM(MM,ID,IYYY)
!
!-----------------------------------------------------------------------
!
!     use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!------------------------------Arguments--------------------------------
!
      integer mm, id, iyyy
!
!-----------------------------Parameters------------------------------
!
      integer igreg
      PARAMETER (IGREG=15+31*(10+12*1582))
!
!---------------------------Local variables-----------------------------
!
      integer ja, jm, jy
!
!-----------------------------------------------------------------------
!
!!!compiler warning      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.EQ.0) STOP 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY_WAM=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY_WAM=JULDAY_WAM+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
  END FUNCTION JULDAY_WAM
      
  subroutine svdcmp( a, m, n, mp, np, w, v )
!------------------------------------------------------------------------- 
! purpose: singular value decomposition
!
! method:
! given a matrix a(1:m,1:n), with physical dimensions mp by np,
! this routine computes its singular value decomposition,
! the matrix u replaces a on output. the
! diagonal matrix of singular values w is output as a vector
! w(1:n). the matrix v (not the transpose v^t) is output as
! v(1:n,1:n).
!
! author: a. maute dec 2003      
! (* copyright (c) 1985 numerical recipes software -- svdcmp *!
! from numerical recipes 1986 pp. 60 or can be find on web-sites
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      integer, intent(in)     :: m
      integer, intent(in)     :: n
      integer, intent(in)     :: mp
      integer, intent(in)     :: np
      real, intent(inout) :: a(mp,np)
      real, intent(out)   :: v(np,np)
      real, intent(out)   :: w(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, its, j, k, l, nm
      real :: anorm
      real  :: c
      real  :: f
      real  :: g
      real  :: h
      real  :: s
      real  :: scale
      real  :: x, y, z
      real  :: rv1(nmax)
      logical  :: cnd1
      logical  :: cnd2

      g     = 0.0
      scale = 0.0
      anorm = 0.0

      do i = 1,n  !loop1
        l = i + 1
        rv1(i) = scale*g
        g     = 0.0
        s     = 0.0
        scale = 0.0
        if( i <= m ) then
          do k = i,m
            scale = scale + abs(a(k,i))
          end do
          if( scale /= 0.0 ) then
            do k = i,m
              a(k,i) = a(k,i)/scale
              s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,i) = f - g
            if( i /= n ) then
              do j = l,n
                s = 0.0
                do k = i,m
                  s = s + a(k,i)*a(k,j)
                end do
                f = s/h
                do k = i,m
                  a(k,j) = a(k,j) + f*a(k,i)
                end do
              end do
            end if
            do k = i,m
              a(k,i) = scale*a(k,i)
            end do
          endif
        endif
        w(i) = scale *g
        g     = 0.0
        s     = 0.0
        scale = 0.0
        if( i <= m .and. i /= n ) then
          do k = l,n
            scale = scale + abs(a(i,k))
          end do
          if( scale /= 0.0 ) then
            do k = l,n
              a(i,k) = a(i,k)/scale
              s      = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k = l,n
              rv1(k) = a(i,k)/h
            end do
            if( i /= m ) then
              do j = l,m
                s = 0.0
                do k = l,n
                  s = s + a(j,k)*a(i,k)
                end do
                do k = l,n
                  a(j,k) = a(j,k) + s*rv1(k)
                end do
              end do
            end if
            do k = l,n
              a(i,k) = scale*a(i,k)
            end do
          end if
        end if
        anorm = max( anorm,(abs(w(i)) + abs(rv1(i))) )
      enddo !loop1

      do i = n,1,-1
        if( i < n ) then
          if( g /= 0.0 ) then
            do j = l,n
              v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l,n
              s = 0.0
              do k = l,n
                s = s + a(i,k)*v(k,j)
              end do
              do k = l,n
                v(k,j) = v(k,j) + s*v(k,i)
              end do
            end do
          end if
          do j = l,n
            v(i,j) = 0.0
            v(j,i) = 0.0
          end do
        end if
        v(i,i) = 1.0
        g = rv1(i)
        l = i
      end do

      do i = n,1,-1
        l = i + 1
        g = w(i)
        if( i < n ) then
          do j = l,n
            a(i,j) = 0.0
          end do
        end if
        if( g /= 0.0  ) then
          g = 1./g
          if( i /= n ) then
            do j = l,n
              s = 0.0
              do k = l,m
                s = s + a(k,i)*a(k,j)
              end do
              f = (s/a(i,i))*g
              do k = i,m
                a(k,j) = a(k,j) + f*a(k,i)
              end do
            end do
          end if
          do j = i,m
            a(j,i) = a(j,i)*g
          end do
        else
          do j = i,m
            a(j,i) = 0.0
          end do
        end if
        a(i,i) = a(i,i) + 1.0
      end do

      do k = n,1,-1
        do its = 1,30 !loop2
          do l = k,1,-1
            nm = l - 1
            cnd1 = abs( rv1(l) ) + anorm == anorm
            if( cnd1 ) then
              cnd2 = .false.
              exit
            end if
            cnd2 = abs( w(nm) ) + anorm == anorm
            if( cnd2 ) then
              cnd1 = .true.
              exit
            else if( l == 1 ) then
              cnd1 = .true.
              cnd2 = .true.
            end if
          end do

          if( cnd2 ) then
            c = 0.0
            s = 1.0
            do i = l,k
              f = s*rv1(i)
              if( (abs(f) + anorm) /= anorm ) then
                g = w(i)
                h = sqrt(f*f + g*g)
                w(i) = h
                h = 1.0/h
                c = (g*h)
                s = -(f*h)
                do j = 1,m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
                end do
              end if
            end do
          end if

          if( cnd1 ) then
            z = w(k)
            if( l == k ) then
              if( z < 0.0 ) then
                w(k) = -z
                do j = 1,n
                  v(j,k) = -v(j,k)
                end do
              end if
!             exit loop2
              go to 20
            end if
          end if

          x = w(l)
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y)
          g = sqrt( f*f + 1.0 )
          f = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x
          c = 1.0
          s = 1.0
          do j = l,nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = sqrt( f*f + h*h )
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do nm = 1,n
              x = v(nm,j)
              z = v(nm,i)
              v(nm,j) = (x*c)+(z*s)
              v(nm,i) = -(x*s)+(z*c)
            end do
            z = sqrt( f*f + h*h )
            w(j) = z
            if( z /= 0.0 ) then
              z = 1.0/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do nm = 1,m
              y = a(nm,j)
              z = a(nm,i)
              a(nm,j) = (y*c)+(z*s)
              a(nm,i) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0
          rv1(k) = f
          w(k)   = x
        end do  !loop2
   20 continue
      end do
      
  end subroutine svdcmp

!-------------------------------------------------------------------------      
! purpose: solves a*x = b
!
! method:     
! solves a*x = b for a vector x, where a is specified by the arrays
! u,w,v as returned by svdcmp. m and n
! are the logical dimensions of a, and will be equal for square matrices.
! mp and np are the physical dimensions of a. b(1:m) is the input right-hand 
! side. x(1:n) is the output solution vector. no input quantities are 
! destroyed, so the routine may be called sequentially with different b
!
! author:  a. maute dec 2002   
! (* copyright (c) 1985 numerical recipes software -- svbksb *!
! from numerical recipes 1986 pp. 57 or can be find on web-sites
!-------------------------------------------------------------------------      

  subroutine svbksb( u, w, v, m, n, mp, np, b, x )
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      integer, intent(in)   :: mp
      integer, intent(in)   :: np
      real , intent(in)  :: u(mp,np)
      real , intent(in)  :: w(np)
      real , intent(in)  :: v(np,np)
      real , intent(in)  :: b(mp)
      real , intent(out) :: x(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, j, jj
      real :: s
      real :: tmp(nmax)

      do j = 1,n
        s = 0. 
        if( w(j) /= 0. ) then
          do i = 1,m
            s = s + u(i,j)*b(i)
          end do
          s = s/w(j)
        endif
        tmp(j) = s
      end do

      do j = 1,n
        s = 0. 
        do jj = 1,n
          s = s + v(j,jj)*tmp(jj)
        end do
        x(j) = s
      end do

  end subroutine svbksb
!
! Efield/Weimer Functions
!      
  integer function value_locate(vec,val)
!
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val >= vec(i+1)
! Input vec must be monotonically increasing
!
        implicit none
! Args:
        real,intent(in) :: vec(:),val
!
! Local:
        integer :: n,i
!
        value_locate = 0
        n = size(vec)
        if (val < vec(1)) return
        if (val > vec(n)) then
          value_locate = n
          return
        endif
        do i=1,n-1
          if (val >= vec(i) .and. val <= vec(i+1)) then
            value_locate = i
            return
          endif
        enddo

  end function value_locate
!-----------------------------------------------------------------------
  real function lngamma(xx)
!
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
!
      implicit none
      real,intent(in) :: xx
      real :: x,y,tmp,ser
      real :: cof(6) = (/76.18009172947146, -86.50532032941677, 24.01409824083091, &
                        -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5/)
      integer :: j

      y = xx
      x = xx
      tmp = x+5.5
      tmp = tmp-(x+0.5)*log(tmp)
      ser = 1.000000000190015
      do j=1,5
        y = y+1
        ser = ser+cof(j)/y
      enddo
      lngamma = -tmp+log(2.5066282746310005*ser/x)
  end function lngamma
!-----------------------------------------------------------------------
  real function factorial(n)
      implicit none
      integer,intent(in) :: n
      integer :: m
      if (n <= 0) then
        write(iulog,"('>>> factorial: n must be positive: n=',i4)") n
        factorial = 0.
        return
      endif
      if (n == 1) then
        factorial = 1.
        return
      endif
      factorial = float(n)
      do m = n-1,1,-1
        factorial = factorial * float(m)
      enddo
  end function factorial
!-----------------------------------------------------------------------     
end module wamphys_math_interp
