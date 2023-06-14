module wamphys_weimer2005
!================================
!
!================================
  use machine,  only: kind_phys
  use physcons, only: con_pi 
  use wam_efieldw05_read_data
  use wam_efield_setdef_data
  use wamphys_const, only : rad2deg=>rad_to_deg, deg2rad=>deg_to_rad
  use wamphys_math_interp, only :interpol_quad
  implicit none
!      real(kind_phys) :: rad2deg,deg2rad           ! set by SetModel_new
!      real(kind_phys) :: bndyfitr                  ! calculated by setboundary
!      real(kind_phys) :: esphc(csize),bsphc(csize) ! calculated by SetModel_new
!      real(kind_phys) :: tmat(3,3),ttmat(3,3)      ! from setboundary

!      integer,parameter :: mxtablesize=500 
      
!      real(kind_phys) :: plmtable(mxtablesize,csize),colattable(mxtablesize)
!      real(kind_phys) :: nlms(csize)
  contains
!-----------------------------------------------------------------------
  subroutine weimer05(epoto)
    implicit none
! Args:
!  real(kind_phys), intent(in)  :: angle, bt, swvel, swden, tilt

      real(kind_phys), parameter :: angle= 0.0, bt= 4.3368, swvel=343.66
      real(kind_phys), parameter :: swden= 5.0, tilt  = 0.0
      real(kind_phys), intent(out) :: epoto(2,22,180)
! Local:
      real(kind_phys)              :: sangle, stilt

      epoto = 0.0

! for northern hemisphere:  tilt,  angle
      call SetModel_new(angle,bt,tilt,swvel,swden)
      call get_elec_field(epoto(1,1:22,1:180))

! for southern hemisphere:  -tilt, BY will be 360-BY(north)
      sangle = mod(360.0-angle, 360.0)
      stilt  = -tilt
      
      call SetModel_new(sangle,bt,stilt,swvel,swden)
      call get_elec_field(epoto(2,1:22,1:180))

  end subroutine weimer05
!-----------------------------------------------------------------------
  subroutine get_elec_field(epot)
    implicit none
! Args:
! real(kind_phys), intent(out) :: ex(22,180), ey(22,180)
      real(kind_phys), intent(out) :: epot(22,180)
! Local: 
      real(kind_phys), parameter   :: fill = 1.0e36
      real(kind_phys), parameter   :: radi = 6371.e3
!  real(kind_phys), parameter   :: dellat = 1.0, dellt = 0.10, delmlt = 1.20
      real(kind_phys), parameter   :: dellat = 1.0, dellt = 0.10, delmlt = 0.133333333333
      integer           :: l, m

      real(kind_phys)              :: gmlt, gmlte, gmltw, dely
      real(kind_phys)              :: glat, glatu, glatd, delx
      real(kind_phys)              :: epotu, epotd, epote, epotw
      real(kind_phys)              :: epx, epy, phir, cphi, sphi
      real(kind_phys)              :: xerz
      
!     delx = radi * 2.0 * dellat * rad2deg 
      xerz = radi*2.0*dellt*15.0*deg2rad
      epot = 0.

 mlt_loop: do l = 1, 180 
 
        gmlt  = real(mod((l-1)*delmlt+12.0,24.0d0))
        gmlte = real(mod(gmlt+dellt,24.d0))
        gmltw = real(mod(gmlt-dellt,24.d0))
	
 lat_loop: do m = 2, 22 
        glat  = 90.0 - (m-1.0)*2.0
        glatu = glat + dellat
        glatd = glat - dellat
        dely  = xerz*cos(glat*deg2rad)

        call EpotVal_new(glat,  gmlt, epot(m,l))

    enddo lat_loop
    enddo mlt_loop

  end subroutine get_elec_field
!-----------------------------------------------------------------------
  subroutine SetModel_new(angle,bt,tilt,swvel,swden)
    implicit none
! Args:
! file_path: directory in which to find data file (must have "/" at end)
      real(kind_phys),intent(in) :: angle,bt,tilt,swvel,swden
! Local:
      integer :: i,j
      real(kind_phys) :: stilt,stilt2,sw,swp,swe,c0,rang,cosa,sina,cos2a,sin2a
      real(kind_phys) :: cfits(d1_pot,csize),a(d1_pot)
!
      if (trim(model) /= 'epot'.and.trim(model) /= 'bpot') then
        write(iulog,"('>>> model=',a)") trim(model)
        write(iulog,"('>>> SetModel_new: model must be either epot or bpot')")
        stop 'SetModel_new' 
      endif
! Read data:	
      call setboundary(angle,bt,tilt,swvel,swden)

      stilt = sin(tilt*deg2rad)
      stilt2 = stilt**2
      sw = bt*swvel/1000.
      swe = (1.-exp(-sw*ex_pot(2)))*sw**ex_pot(1)
      c0 = 1.
      swp = swvel**2 * swden*1.6726e-6
      rang = angle*deg2rad
      cosa = cos(rang)
      sina = sin(rang)
      cos2a = cos(2.*rang)
      sin2a = sin(2.*rang)
      if (bt < 1.) then ! remove angle dependency for IMF under 1 nT
        cosa = -1.+bt*(cosa+1.)
        cos2a = 1.+bt*(cos2a-1.)
        sina = bt*sina
        sin2a = bt*sin2a
      endif
      cfits = schfits ! schfits(d1_pot,csize) is in module w05read_data
      a = (/c0      , swe       , stilt      , stilt2     , swp, &
            swe*cosa, stilt*cosa, stilt2*cosa, swp*cosa, &
            swe*sina, stilt*sina, stilt2*sina, swp*sina, &
            swe*cos2a,swe*sin2a/)
      if (trim(model) == 'epot') then
        esphc(:) = 0.
        do j=1,csize
          do i=1,int(d1_pot)
            esphc(j) = esphc(j)+cfits(i,j)*a(i)
          enddo
        enddo
!     write(iulog,"('SetModel_new: esphc=',/,(6e12.4))") esphc
      else
        bsphc(:) = 0.
        do j=1,csize
          do i=1,int(d1_pot)
            bsphc(j) = bsphc(j)+cfits(i,j)*a(i)
          enddo
        enddo
!      write(iulog,"('SetModel_new: bsphc=',/,(6e12.4))") bsphc
      endif
  end subroutine SetModel_new
!-----------------------------------------------------------------------
  subroutine setboundary(angle,bt,tilt,swvel,swden)        ! Zhuxiao
    implicit none
! Args:
! file_path: directory in which to find data file (must have "/" at end)
!     character(len=*),intent(in) :: file_path   ! by Zhuxiao
      real(kind_phys),intent(in) :: angle,bt,tilt,swvel,swden
! Local:
      integer :: i
      real(kind_phys) :: swp,xc,theta,ct,st,tilt2,cosa,btx,x(na),c(na)
      real(kind_phys), parameter :: num_0 = 0., num_1 = 1. 
! Calculate the transformation matrix to the coordinate system
! of the offset pole.
      xc = 4.2
      theta = xc*(deg2rad)
      ct = cos(theta)
      st = sin(theta)

      tmat(1,:) = (/ ct, num_0, st/)   ! avoid conflict
      tmat(2,:) = (/ 0., 1., 0./) 
      tmat(3,:) = (/-st, num_0, ct/)   ! avoid conflict

      ttmat(1,:) = (/ct, num_0,-st/)
      ttmat(2,:) = (/ 0.,1., 0./)
      ttmat(3,:) = (/st, num_0, ct/)   ! avoid conflict

      swp = swden*swvel*swvel*1.6726e-6   ! pressure
      tilt2 = tilt**2
      cosa = cos(angle*deg2rad)
      btx = 1.-exp(-bt*ex_bndy(1))
      if (bt > 1.) then
        btx = btx*bt**ex_bndy(2)
      else
        cosa = 1.+bt*(cosa-1.) ! remove angle dependency for IMF under 1 nT
      endif
      x = (/num_1, cosa, btx, btx*cosa, swvel, swp/)
      c = bndya
      bndyfitr = 0.
      do i=1,na
        bndyfitr = bndyfitr+x(i)*c(i)
      enddo
!     write(iulog,"('setboundary: cosa=',f8.3,' btx=',f8.3)") cosa,btx
!     write(iulog,"('setboundary: bndyfitr=',e12.4)") bndyfitr
  end subroutine setboundary
!-----------------------------------------------------------------------
  subroutine EpotVal_new(lat,mlt,epot)
    implicit none
! Args:
      real(kind_phys),intent(in) :: lat,mlt
      real(kind_phys),intent(out) :: epot
! Local:
      integer :: inside,j,m,mm,skip
      real(kind_phys) :: z,phir,plm,colat,nlm
      real(kind_phys) :: phim(2),cospm(2),sinpm(2)
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
      call checkinputs(lat,mlt,inside,phir,colat)

      if (inside == 0) then
        epot = 0.
        return
      endif

! IDL code: 
!   phim=phir # replicate(1,maxm) * ((indgen(maxm)+1) ## replicate(1,n_elements(phir)))
!   where the '#' operator multiplies columns of first array by rows of second array,
!   and the '##' operator multiplies rows of first array by columns of second array.
!   Here, maxm == maxm_pot == 2 (from w05read_data module), and phir is a scalar. The 
!   above IDL statement then becomes: phim = ([phir] # [1,1]) * ([1,2] ## [phir]) where
!   phim will be dimensioned [1,2]
!
      phim(1) = phir
      phim(2) = phir*2.
      cospm(:) = cos(phim(:))
      sinpm(:) = sin(phim(:))

      z = 0.
      skip=0
      do j=1,csize
        if (skip == 1) then
          skip = 0
          cycle
        endif
        m = ms(j)
        if (ab(j)==1) then

          plm = scplm(j,colat,nlm) ! scplm function is in this module

          skip = 0
          if (m == 0) then
            z = z+plm*esphc(j)
          else
            z = z+plm*(esphc(j)*cospm(m)+esphc(j+1)*sinpm(m))
            skip = 1
          endif

        endif ! ab(j)
      enddo
      epot = z 
  end subroutine EpotVal_new
!-----------------------------------------------------------------------
  subroutine mpfac(lat,mlt,mpmpfac)
    implicit none
! Args:
      real(kind_phys),intent(in) :: lat,mlt
      real(kind_phys),intent(out) :: mpmpfac
! Local:
      integer :: j,m,inside,skip
      real(kind_phys) :: phim(2),cospm(2),sinpm(2),cfactor
      real(kind_phys) :: re,z,phir,plm,colat,nlm, ppi

      re = 6371.2 + 110. ! km radius (allow default ht=110)

!  checkinputs returns inside=1 if lat is inside model boundary,
!  inside=0 otherwise. Phir and colat are also returned by checkinputs.
      call checkinputs(lat,mlt,inside,phir,colat)
      if (inside == 0) then
        mpmpfac = 0.
        return
      endif

      phim(1) = phir
      phim(2) = phir*2.
      cospm(:) = cos(phim(:))
      sinpm(:) = sin(phim(:))

      z = 0.
      skip=0
      jloop: do j=1,csize
        if (skip == 1) then
          skip = 0
          cycle
        endif
        if (ls(j) >= 11) exit jloop
        m = ms(j)
        if (ab(j) == 1) then
          plm = scplm(j,colat,nlm) ! colat and nlm are returned (both reals)
          plm = plm*(nlm*(nlm+1.))
  ! bsphc was calculated in SetModel_new (when SetModel_new called with 'bpot')
          if (m==0) then
            z = z-plm*bsphc(j)
          else
            z = z-(plm*(bsphc(j)*cospm(m)+bsphc(j+1)*sinpm(m)))
            skip = 1
          endif
        endif
      enddo jloop ! j=1,csize
!     ppi = 4.*atan(1.)
      cfactor = -1.e5/(4.*con_pi*re**2) ! convert to uA/m2
!     cfactor = -1.e5/(4.*ppi*re**2) 
      z = z*cfactor
      mpmpfac = z
!   write(iulog,"('mpfac: lat=',f8.3,' mlt=',f8.3,' mpmpfac=',1pe12.4)") lat,mlt,mpmpfac
  end subroutine mpfac
!-----------------------------------------------------------------------
  real function scplm(indx,colat,nlm)
    implicit none
! Args:
      integer,intent(in) :: indx
      real(kind_phys),intent(in) :: colat
      real(kind_phys),intent(out) :: nlm
! Local:
      integer,save :: tablesize
      integer :: istat,i,j,l,m,skip
      real(kind_phys) :: th0,output(1),colata(1),plm1
      real(kind_phys) :: cth(mxtablesize)
      real(kind_phys),save :: prevth0=1.e36

      scplm = 0.
      th0 = bndyfitr
      if (prevth0 /= th0) then
        tablesize = 3*nint(th0)
        if (tablesize > mxtablesize) then 
          write(iulog,"('>>> tablesize > mxtablesize: tablesize=',i5,      &
              ' mxtablesize=',i5,' tn0=',e12.4)") tablesize,mxtablesize,th0
          stop 'tablesize'
        endif
!       write(iulog,"('scplm: indx=',i3,' colat=',f8.3,' th0=',e12.4, &
!       ' tablesize=',i3)") indx,colat,th0,tablesize
        do i=1,tablesize
          colattable(i) = float(i-1)*(th0/float(tablesize-1))
          cth(i) = cos(colattable(i)*deg2rad)
        enddo

!       write(iulog,"('scplm: tablesize=',i4,' colattable=',/,(6f8.3))") &
!       tablesize,colattable(1:tablesize)
        prevth0 = th0
        nlms = 0. ! whole array init 
        skip=0
        do j=1,csize
          if (skip == 1) then
            skip = 0
            cycle
          endif
          l = ls(j)
          m = ms(j)
          nlms(j) = nkmlookup(l,m,th0) ! nkmlookup in this module
!         write(iulog,"('scplm after nkmlookup: j=',i3,' l=',i3,' m=',i3, &
!         ' nlms(j)=',f8.4)") j,l,m,nlms(j)
          call pm_n(m,nlms(j),cth,plmtable(1:tablesize,j),tablesize)
          skip = 0
          if (m /= 0 .and. ab(j) > 0) then
            plmtable(1,j+1) = plmtable(1,j)
            nlms(j+1) = nlms(j)
            skip = 1
          endif
        enddo ! j=1,csize
      endif ! prevth0
      nlm = nlms(indx)
      colata(1) = colat

      call interpol_quad(plmtable(1:tablesize,indx),colattable(1:tablesize), &
                        colata,output)
      scplm = output(1)
  end function scplm
!-----------------------------------------------------------------------
  subroutine pm_n(m,r,cth,plmtable,tablesize)
    implicit none
! Args:
      integer,intent(in) :: m,tablesize
      real(kind_phys),intent(in) :: r
      real(kind_phys),intent(in) :: cth(tablesize)
      real(kind_phys),intent(out) :: plmtable(tablesize)
! Local:
      integer :: i,k,ii
      real(kind_phys) :: rm,rk,div,ans,xn
      real(kind_phys),dimension(tablesize) :: a,x,tmp,table

      if (m == 0) then 
        a = 1. ! whole array op
      else
        do i=1,tablesize
          a(i) = sqrt(1.-cth(i)**2)**m
        enddo
      endif
      xn = r*(r+1.)
      x(:) = (1.-cth(:))/2.

      table = a ! whole array init

      k = 1
      pmn_loop: do         ! repeat-until loop in idl code
        do i=1,tablesize
          rm = float(m)
          rk = float(k)
          a(i) = a(i)*(x(i)*((rk+rm-1.)*(rk+rm)-xn)/(rk*(rk+rm)))
          table(i) = table(i)+a(i) ! "result" in idl code
        enddo
  !     write(iulog,"('pm_n: k=',i3,' a=',/,(6(1pe12.4)))") k,a
  !     write(iulog,"('pm_n: k=',i3,' table=',/,(6(1pe12.4)))") k,table
        k = k+1
        do i=1,tablesize
          div = abs(table(i))
          if (div <= 1.e-6) div = 1.e-6
          tmp(i) = abs(a(i)) / div
        enddo
  !     write(iulog,"('pm_n: k=',i3,' abs(a)=',/,(6(1pe12.4)))") k,abs(a)
  !     write(iulog,"('pm_n: k=',i3,' abs(table)=',/,(6(1pe12.4)))") k,abs(table)
        if (maxval(tmp) < 1.e-6) exit pmn_loop
      enddo pmn_loop
      ans = km_n(m,r)

      plmtable(:) = table(:)*ans

      end subroutine pm_n
!-----------------------------------------------------------------------
  real function km_n(m,rn)
    implicit none
! Args:
      integer,intent(in) :: m
      real(kind_phys),intent(in) :: rn
! Local:
      integer :: i,n
      real(kind_phys) :: rm

      if (m == 0) then 
        km_n = 1.
        return
      endif
      
      rm = float(m)
      km_n = sqrt(2.*exp(lngamma(rn+rm+1.)-lngamma(rn-rm+1.))) / (2.**m*factorial(m))

  end function km_n
!-----------------------------------------------------------------------
  real function nkmlookup(k,m,th0)
    implicit none
! Args:
      integer,intent(in) :: k,m
      real(kind_phys),intent(in) :: th0
! Local:
      integer :: kk,mm
      real(kind_phys) :: th0a(1),out(1)

      if (th0 == 90.) then
        nkmlookup = float(k)
        return
      endif
      th0a(1) = th0
      kk = k+1
      mm = m+1
      if (kk > maxk_scha) then
!       write(iulog,"('>>> nkmlookup: kk > maxk: kk=',i4,' maxk=',i4)") kk,maxk_scha
        call interpol_quad(allnkm(maxk_scha,mm,:),th0s,th0a,out)
      endif
      if (mm > maxm_scha) then
!       write(iulog,"('>>> nkmlookup: mm > maxm: kk=',i4,' maxm=',i4)") kk,maxm_scha
        call interpol_quad(allnkm(kk,maxm_scha,:),th0s,th0a,out)
      endif
      if (th0 < th0s(1)) then
!       write(iulog,"('>>> nkmlookup: th0 < th0s(1): th0=',e12.4,' th0s(1)=',e12.4)")
!           th0,th0s(1)
      endif
!     write(iulog,"('nkmlookup call interpol: kk=',i3,' mm=',i3,' th0=',e12.4,&
!     &' allnkm=',/,(6(1pe12.4)))") kk,mm,th0a,allnkm(kk,mm,:)

      call interpol_quad(allnkm(kk,mm,:),th0s,th0a,out)

      nkmlookup = out(1)

  end function nkmlookup
!-----------------------------------------------------------------------
  subroutine checkinputs(lat,mlt,inside,phir,colat)
    implicit none
! Args:
      real(kind_phys),intent(in) :: lat,mlt
      integer,intent(out) :: inside
      real(kind_phys),intent(out) :: phir,colat
! Local:
      real(kind_phys) :: lon,tlat,tlon,radii

      lon = mlt*15.
      call dorotation(lat,lon,tlat,tlon)
      radii = 90.-tlat
      inside = 0
      if (radii <= bndyfitr) inside = 1 ! bndyfitr from setboundary
      phir = tlon*deg2rad
      colat = radii

  end subroutine checkinputs
!-----------------------------------------------------------------------
  subroutine dorotation(latin,lonin,latout,lonout)
    implicit none
! Args:
      real(kind_phys),intent(in) :: latin,lonin
      real(kind_phys),intent(out) :: latout,lonout
! Local:
      real(kind_phys) :: latr,lonr,stc,ctc,sf,cf,a,b,pos(3)
      integer :: i

      latr = latin*deg2rad
      lonr = lonin*deg2rad
      stc = sin(latr)
      ctc = cos(latr)
      sf = sin(lonr)
      cf = cos(lonr)
      a = ctc*cf
      b = ctc*sf
! IDL code: Pos= TM ## [[A],[B],[STC]]
! The ## operator multiplies rows of first array by columns of second array.
! Currently, TM(3,3) = Tmat (or TTmat if "reversed" was set)
! If called w/ single lat,lon, then a,b,stc are dimensioned (1), and
! Pos is then (1,3)
      do i=1,3
        pos(i) = tmat(1,i)*a + tmat(2,i)*b + tmat(3,i)*stc
      enddo
    
      latout = asin(pos(3))*rad2deg
      lonout = atan2(pos(2),pos(1))*rad2deg
      
  end subroutine dorotation
!-----------------------------------------------------------------------    
!-----------------------------------------------------------------------
  integer function value_locate(vec,val)
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val >= vec(i+1)
! Input vec must be monotonically increasing
    implicit none
! Args:
      real(kind_phys),intent(in) :: vec(:),val
  ! Local:
      integer :: n,i

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
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
    implicit none
      real(kind_phys),intent(in) :: xx
      real(kind_phys) :: x,y,tmp,ser
      real(kind_phys) :: cof(6) = (/76.18009172947146, -86.50532032941677, 24.01409824083091, &
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
end module wamphys_weimer2005
