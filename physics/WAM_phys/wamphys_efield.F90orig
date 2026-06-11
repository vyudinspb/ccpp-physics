
module  wamphys_efield
      
!--------------------------------------------------------------------- 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
! - low/midlatitudes electric potential is from an empirical model from
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! - the transition zone is smoothed
! - output is the horizontal global electric field in magnetic coordinates direction
!  at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real:: ut,           ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real ::               
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!  
! Author: A. Maute Dec 2003  am 12/30/03
!
!   Apr 06 2012 Henry Juang, initial implement for nems
!   Nov 20 2014 Jun   Wang,  change JULDAY to JULDAY_WAM
!   Mar 18 2017 Zhuxiao Li and Tzu-Wei, add option to read the solar
!   wind related driving parameteres from outside parameter file.
!   May 2018 Zhuxiao Li, add the code to call Weimer2005 (w05sc) .
!------------------------------------------------------------------------------ 
     
      use wam_efieldw05_read_data, only: read_bndy, &
                                      read_potential, read_schatable, model
				      
      use wamphys_weimer2005, only:EpotVal_new, SetModel_new				      
      use wam_efield_setdef_data
      use wamphys_math_interp, only : adjust
!subroutine magnetic_grids        
!subroutine set_readcoef
!subroutine efread_acoef 
!subroutine read_acoef_efield
!subroutine pnm   
!subroutine prep_fk
!subroutine index_quiet
!subroutine ff
!subroutine prep_pnm    
   
      implicit none

!---------------------------------------------------------------------- 
! solar parameters
!---------------------------------------------------------------------- 
      real(kind=kind_phys) ::   f107d           ! 10.7 cm solar flux
      real(kind=kind_phys) ::   by              ! By component of IMF [nT]
      real(kind=kind_phys) ::   bz              ! Bz component of IMF [nT]
      
      integer, parameter   :: iulog =6
      logical, parameter :: iutav=.false.      ! .true.  means UT-averaging          

      logical, parameter :: debug =.false.
      
      
      
!========================================================================
      contains      

  subroutine wamget_efield(swbt, swang, swvel, swbz, swden)
!-----------------------------------------------------------------------
! Purpose: calculates the global electric potential field on the
!          geomagnetic grid (MLT in deg) and derives the electric field 
!
! Method:
!
! Author: A. Maute Dec 2003  am 12/17/03    
!-----------------------------------------------------------------------
      real, intent(in) :: swbt, swang, swvel, swbz, swden
      integer :: idum1, idum2, tod ! time of day [s] 
      real    :: kp

      tod=ut*3600.
!     ut = tod/3600.                   ! UT of day [sec]


      call adj_S_a
!-----------------------------------------------------------------------
! calculate global electric potential    
!-----------------------------------------------------------------------
      call GlobalElPotential(swbt, swang, swvel, swbz, swden)
!     print*,'pot_efield',potent(149,66),potent(149,64)

!-----------------------------------------------------------------------
! calculate derivative of global electric potential 
!-----------------------------------------------------------------------
      call DerivPotential
!     print*,'ed2_efield',ed2(149,65),potent(149,66),potent(149,64)

  end subroutine wamget_efield

  subroutine GlobalElPotential(swbt, swang, swvel, swbz, swden)
!-----------------------------------------------------------------------
! Purpose: calculates the global electric potential field on the
!          geomagnetic grid (MLT in deg) 
!
! Method: rewritten code from Luedger Scherliess (11/20/99 LS)
!     routine to calculate the global electric potential in magnetic
!     Apex coordinates (Latitude and MLT).
!     High Latitude Model is Weimer 1996.
!     Midlatitude model is Scherliess 1999.
!     Interpolation in a transition region at about 60 degree 
!     magnetic apex lat
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-----------------------------------------------------------------------
      real, intent(in) :: swbt, swang, swvel, swbz, swden
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: ilon, ilat, idlat
      integer  :: ihlat_bnd(0:nmlon)     ! high latitude boundary
      integer  :: itrans_width(0:nmlon)  ! width of transition zone
      
      
      real :: mlt, mlon, mlat, mlat_90, pot
      real :: pot_midlat(0:nmlon,0:nmlat) ! potential from L. Scherliess model
      real :: pot_highlat(0:nmlon,0:nmlat) ! potential from Weimer model
      real :: pot_highlats(0:nmlon,0:nmlat)! smoothed potential from Weimer model

!-----------------------------------------------------------------------
! convert to date and day	
!-----------------------------------------------------------------------
      day  = iday + ut/24.
      date = iyear + day/dy2yr

      do ilat = 0,nmlath                        ! Calculate only for one magn. hemisphere
	  mlat = ylatm(ilat)                      ! mag. latitude
        do ilon = 0,nmlon	                ! lon. loop
          call efield_mid( mlat, ylonm(ilon), pot )
	    pot_midlat(ilon,ilat+nmlath) = pot	! SH/NH symmetry 
	    pot_midlat(ilon,nmlath-ilat) = pot
        end do
      end do

!-----------------------------------------------------------------------
! hight latitude potential from Weimer model
! at the poles Weimer potential is not longitudinal dependent
!-----------------------------------------------------------------------
      call prep_weimer(swbt, swang, swvel, swbz, swden) ! calculate IMF angle & magnitude, tilt

!$omp parallel do private(ilat, ilon, mlat_90, pot)
      do ilat = 0,nmlat_wei  ! Calculate only for one magn. hemisphere
        mlat_90 = 90. - ylatm(ilat)  ! mag. latitude
        do ilon = 0,nmlon

          call EpotVal_new(mlat_90, ylonm(ilon)*deg2mlt, pot )
          pot = 1000.*pot
!-----------------------------------------------------------------------
! NH/SH symmetry
!-----------------------------------------------------------------------
          pot_highlat(ilon,ilat)        = pot
          pot_highlat(ilon,nmlat-ilat)  = pot
          pot_highlats(ilon,ilat)       = pot
          pot_highlats(ilon,nmlat-ilat) = pot
        end do
      end do     
!     print*,'www2','highlat',ut,by,bz,pot_highlat(0:180,68),nmlat_wei

!-----------------------------------------------------------------------
! weighted smoothing of high latitude potential
!-----------------------------------------------------------------------
      idlat = 2              ! smooth over -2:2 = 5 grid points
      call pot_latsmo( pot_highlats, idlat )
!     print*,'www2','highlat',ut,pot_highlat(0:180,45)
!-----------------------------------------------------------------------
! calculate the height latitude bounday ihl_bnd
! 1. calculate E field from weimar model
!    boundary is set where the total electric field exceeds
!    0.015 V/m (corresp. approx. to 300 m/s)
! 2. moved halfways to 54 deg 
! output : index 0-pole nmlath-equator
!-----------------------------------------------------------------------
      call highlat_getbnd( ihlat_bnd )
!-----------------------------------------------------------------------
! 3. adjust high latitude boundary sinusoidally
!    calculate width of transition zone
!-----------------------------------------------------------------------
      call bnd_sinus( ihlat_bnd, itrans_width ) 
!-----------------------------------------------------------------------
! 4. ajust high latitude potential to low latitude potential      
!-----------------------------------------------------------------------
!     print*,'www30',ihlat_bnd
      call highlat_adjust( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd )
!     print*,'www3','highlat',ut,pot_highlat(145:153,68)
!     print*,'www3','midlat',ut,pot_midlat(145:153,68)
!-----------------------------------------------------------------------
! interpolation of high and low/midlatitude potential in the
! transition zone and put it into global potent array
!-----------------------------------------------------------------------
      call interp_poten( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd, itrans_width) 
!     print*,'www4','potent',ut,by,bz,potent(0:181,68)
!-----------------------------------------------------------------------
! potential weighted smoothing in latitude
!-----------------------------------------------------------------------
      idlat = 2                 ! smooth over -2:2 = 5 grid points
      call pot_latsmo2( potent, idlat )
!     print*,'www5','pot_efield',potent(149,68)
!-----------------------------------------------------------------------
! potential smoothing in longitude
!-----------------------------------------------------------------------
      idlat = nmlon/48          ! smooth over -idlat:idlat grid points
      call pot_lonsmo( potent, idlat )
!     print*,'www6','pot_efield',ut,by,bz,potent(0:180,68)
!-----------------------------------------------------------------------
! output
!-----------------------------------------------------------------------
! output ( change later to netcdf file)
!      do ilat=0,nmlat
!       do ilon=0,nmlon
!         write(iulog,'(4(x,f12.5))') ylatm(ilat),ylonm(ilon), &
!           potent(ilon,ilat),potent(ilon,nmlat-ilat)
!         write(iulog,'(4(x,f12.5))') ylatm(ilat),ylonm(ilon), &
!           potent(ilon,ilat),potent(ilon,nmlat-ilat)
!	write(iulog,'(f10.3)') potent(ilon,ilat)
!       end do
!      end do

  end subroutine GlobalElPotential


                                                                           
  subroutine adj_S_a
!------------------------------------------------------------------
! adjust S_a -> S_aM   eqn.8-11 Scherliess draft
!------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: i
      real :: x2, y2, a90, a180, S_aM

      x2 = 90.*90.
      y2 = (90. - 65.)
      y2 = y2*y2
      a90  = atan2(y2,x2)
      y2 = (180. - 65.)
      y2 = y2*y2
      a180 = atan2(y2,x2)
!     y2 = (S_a-65.)
      y2 = (f107d - 65.)
      y2 = y2*y2
      S_aM = (atan2(y2,x2) - a90)/(a180 - a90) 
      S_aM = 90.*(1. + S_aM)
!     if(debug) write(iulog,*) 'f107d=',f107d,' S_aM =',S_aM
!     if(debug) write(iulog,*) 'By=',by

!-----------------------------------------------------------------
! inter/extrapolate to S_a (f107d)
!----------------------------------------------------------------
      do i = 0,ni                       ! eqn.8 Scherliess draft
      
        a_klnm(i) = S_aM*(a_hf(i)-a_lf(i))/90.+ 2.*a_lf(i)- a_hf(i)
	
! for testing like in original code
!        a_klnm(i)=S_a*(a_hf(i)-a_lf(i))/90.+2.*a_lf(i)-a_hf(i)
!        a_klnm(i)=f107d*(a_hf(i)-a_lf(i))/90.+2.*a_lf(i)-a_hf(i)
      end do

  end subroutine adj_S_a

  subroutine set_fkflfs( fk, fl, fs )
!------------------------------------------------------------------
! Purpose:  set f_-k(day) depending on seasonal flag used for empirical model
!     to calculate the electric potential
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------
      real, intent(out) ::  &
         fk(0:2),  fl(-2:2), fs(2)    ! f_-k(day) f_l(ut) f_s(f10.7) 

!------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: lp
      real :: ang
      real :: lon_ut

!------------------------------------------------------------------
! f_-k(day) 
! use factors for iseasav == 0 - Scherliess had iseasav as an input parameter
!------------------------------------------------------------------
      lp = iseasav
      if( iseasav == 0 ) then
        ang   = (day + 9.)*dy2rd
        fk(0) = sqr2*cos( 2.*ang )
        fk(1) = sqr2*cos( ang )
        fk(2) = 1.
      else if( iseasav >= 1 .and. iseasav <= 3 ) then
        fk(0) = ft(lp,0)
        fk(1) = ft(lp,1)
        fk(2) = ft(lp,2)
      else if( iseasav == 4 ) then
        fk(0) =0.
        fk(1) =0.
        fk(2) =1.
      end if

!-----------------------------------------------------------------
! f_l(ut) 
!-----------------------------------------------------------------
      lon_ut = 15.*ut        ! 15.*mlt - xmlon + 69. 
      call ff( lon_ut, 2, fl )                                                 
      if( iutav ) then  	! UT-averaging    
	  ang   = fl(0)
        fl(:) = 0.
        fl(0) = ang	
      end if

!-----------------------------------------------------------------
! f_s(f10.7)  only fs(1) used  	
!-----------------------------------------------------------------
      fs(1) = 1.
!     fs(2) = S_a			  
      fs(2) = f107d			  

  end subroutine set_fkflfs

  subroutine efield_mid( mlat, mlon, pot )
!------------------------------------------------------------------
! Purpose: calculate the electric potential for low and 
!      midlatitudes from an empirical model (Scherliess 1999)
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------

!------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      real, intent(in)  :: mlat, mlon
      real, intent(out) :: pot               ! electric potential (V)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: i, mp, np, nn
      real :: mod_mlat, ct, x
      real :: fk(0:2)      	    ! f_-k(day) 
      real :: fl(-2:2)          ! f_l(ut)  
      real :: fs(2)	            ! f_s(f10.7) 
      real :: f(-18:18)
      real :: p(0:nm,0:mm)      ! P_n^m	 spherical harmonics

      pot = 0. ! initialize                                        

      mod_mlat = mlat
      if( abs(mlat) <= 0.5 ) then
         mod_mlat = 0.5                     ! avoid geomag.equator
      end if

!------------------------------------------------------------------
! set f_-k, f_l, f_s depending on seasonal flag
!------------------------------------------------------------------
      call set_fkflfs( fk, fl, fs ) 
      
!------------------------------------------------------------------
! spherical harmonics 
!------------------------------------------------------------------
      ct = cos( (90. - mod_mlat)*dtr )  ! magnetic colatitude 
      call pnm( ct, p )	                   ! calculate P_n^m
      call ff( mlon, 18, f )               ! calculate f_m (phi) why 18 if N=12                              

      do i = 0,imax  
        mp  = mf(i)                                                      
        np  = nf(i)
        nn  = abs(mp)                      !   P_n^m = P_n^-m  
        x   = a_klnm(i)* fl(lf(i)) * fk(kf(i)) * fs(jf(i))
	  pot = pot + x*f(mp)*p(np,nn) 
      end do 
      
  end subroutine efield_mid                                              

  subroutine prep_weimer(swbt, swang, swvel, swbz, swden)

!-----------------------------------------------------------------
! Purpose:  for Weimer model calculate IMF angle, IMF magnitude
!  tilt of earth
!
! Method: using functions and subroutines from Weimer Model 1996
!     output:  angle, &  ! IMF angle
!     	       bt,    &  ! IMF magnitude
!     	       tilt      ! tilt of earth
!
! Author: A. Maute Nov 2003  am 11/20/03
!    more output: v_sw, swden

!-----------------------------------------------------------------

      use wamphys_init_module, only : SPW_DRIVERS, SWIN_DRIVERS

      real, intent(in) :: swbt, swang, swvel, swbz, swden

!-----------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------
      real ::   &
       angle,   &   ! IMF angle
       bt,      &   ! IMF magnitude
       tilt,    &   ! tilt of earth
       v_sw,    &   ! sw velocity
       den          ! sw density
!-----------------------------------------------------------------
! function declarations
!-----------------------------------------------------------------
      real :: get_tilt	 ! in wei96.f

!------by Zhuxiao-----
       if (trim(SWIN_DRIVERS)/='swin_wam') then
         den = 5.

         if( by == 0. .and. bz == 0.) then
           angle = 0.
         else
           angle = atan2( by,bz )
         end if

         angle = angle*rtd
         call adjust( angle )
         bt = sqrt( by*by + bz*bz )
      else
         angle = swang
         bt    = swbt
         den   = swden
         v_sw  = swvel
      end if
      if(debug) then
       write(iulog,"(/,'efield prep_weimer:')")
       write(iulog,"(/,'by code:')")
       write(iulog,*)  '  Bz   =',bz
       write(iulog,*)  '  By   =',by
       write(iulog,*)  '  Bt   =',bt
       write(iulog,*)  '  angle=',angle
       write(iulog,*)  '  VSW  =',v_sw
       write(iulog,*)  '  tilt =',tilt
       write(iulog,*)  '  swden =',swden
      end if

!-------------------------------------------------------------------
! use month and day of month - calculated with average no.of days per month
! as in Weimer
!-------------------------------------------------------------------
!     if(debug) write(iulog,*) 'prep_weimer: day->day of month',
!    &iday,imo,iday_m,ut

      call wam_get_tilt( iyear, imo, iday_m, ut, tilt )

      call SetModel_new(angle,bt,tilt,v_sw,den)

      if(debug) then
         write(iulog,"(/,'efield prep_weimer:')")
         write(iulog,"(/,'after SetModel_new:')")
         write(iulog,*)  '  Bz   =',bz
         write(iulog,*)  '  By   =',by
         write(iulog,*)  '  Bt   =',bt
         write(iulog,*)  '  angle=',angle
         write(iulog,*)  '  VSW  =',v_sw
         write(iulog,*)  '  tilt =',tilt
      end if

  end subroutine prep_weimer

  subroutine pot_latsmo( pot, idlat )  ! pots == pot_highlats
!--------------------------------------------------------------------
! Purpose: smoothing in latitude of  potential
!
! Method: weighted smoothing in latitude 
! assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(in)     :: idlat
      real, intent(inout) :: pot(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: ilon, ilat, id
      real :: wgt, del
      real :: w(-idlat:idlat)
!     real :: pot_smo(0:nmlat) ! temp array for smooth. potential
      real :: pot_smo(0:nmlon,0:nmlat_wei) ! temp array for smooth. potential

!------------------------------------------------------------------
! weighting factors (regular grid spacing) 
!------------------------------------------------------------------
      wgt = 0. 
      do id = -idlat,idlat
        del   = abs(id)*dlatm	! delta lat_m
        w(id) = 1./(del + 1.)
        wgt   = wgt + w(id)
      end do
      wgt = 1./wgt

      do ilat = idlat,nmlat_wei-idlat
         pot_smo(:,ilat) = matmul( pot(:,ilat-idlat:ilat+idlat),w )*wgt
      end do

      do ilat = idlat,nmlat_wei-idlat
         pot(:,ilat)       = pot_smo(:,ilat)
         pot(:,nmlat-ilat) = pot_smo(:,ilat)
      end do

  end subroutine pot_latsmo

  subroutine pot_latsmo2( pot, idlat ) 
!------------------------------------------------------------------
! Purpose:  smoothing in latitude of  potential
!
! Method: weighted smoothing in latitude 
!         assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!------------------------------------------------------------------

!------------------------------------------------------------------
!	... dummy arguments
!------------------------------------------------------------------
      integer, intent(in)     :: idlat
      real, intent(inout) :: pot(0:nmlon,0:nmlat)

!------------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: ilon, ilat, id
      real :: wgt, del
      real :: w(-idlat:idlat)
!     real :: pot_smo(0:nmlat) ! temp array for smooth. potential
      real :: pot_smo(0:nmlon,0:nmlath) ! temp array for smooth. potential

!-------------------------------------------------------------------
! weighting factors (regular grid spacing)  
!-------------------------------------------------------------------
      wgt = 0.
      do id = -idlat,idlat
	  del   = abs(id)*dlatm	! delta lat_m
	  w(id) = 1./(del + 1.)
        wgt   = wgt + w(id)
      end do
      wgt = 1./wgt

      do ilat = idlat,nmlath-idlat
         pot_smo(:,ilat) = matmul( pot(:,ilat-idlat:ilat+idlat),w )*wgt
      end do

      do ilat = idlat,nmlath-idlat
         pot(:,ilat) = pot_smo(:,ilat)
      end do

  end subroutine pot_latsmo2

  subroutine pot_lonsmo( pot, idlon ) 
!-------------------------------------------------------------------
! Purpose: smoothing in longitude of potential
!
! Method:  weighted smoothing in longitude
!          assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(in)     :: idlon
      real, intent(inout) :: pot(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: ilon, ilat, id, iabs
      real :: wgt, del
      real :: w(-idlon:idlon)
      real :: pot_smo(0:nmlath) ! temp array for smooth. potential
      real :: tmp(-idlon:nmlon+idlon) ! temp array for smooth. potential

!-------------------------------------------------------------------
! weighting factors (regular grid spacing) 
!-------------------------------------------------------------------
      wgt = 0.
      do id = -idlon,idlon
        del   = abs(id)*dlonm	! delta lon_m
        w(id) = 1./(del + 1.)
        wgt   = wgt + w(id)
      end do
      wgt = 1./wgt
	
      do ilat = 0,nmlath
          tmp(0:nmlon)   = pot(0:nmlon,ilat)
          tmp(-idlon:-1) = pot(nmlon-idlon:nmlon-1,ilat)
          tmp(nmlon+1:nmlon+idlon) = pot(1:idlon,ilat)
          do ilon = 0,nmlon
             pot(ilon,ilat)=dot_product(tmp(ilon-idlon:ilon+idlon),w)*wgt
             pot(ilon,nmlat-ilat) = pot(ilon,ilat)
          end do   
      end do
!     print*,'www8','pot_efield',pot(149,66)
      
  end subroutine pot_lonsmo

  subroutine highlat_getbnd( ihlat_bnd ) 
!------------------------------------------------------------------
! Purpose: calculate the height latitude bounday index ihl_bnd
!
! Method:
! 1. calculate E field from weimar model
!    boundary is set where the total electric field exceeds
!    0.015 V/m (corresp. approx. to 300 m/s)
! 2. moved halfways to 54 deg not necessarily equatorwards as in the
!    original comment from L. Scherliess- or?
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(out) :: ihlat_bnd(0:nmlon)

!------------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: ilon, ilat, ilat_sft_rvs
      real :: mlat, mlt, es, ez, e_tot

      ilat_sft_rvs = nmlath - ilat_sft          ! pole =0, equ=90

      do ilon = 0,nmlon                         ! long.
	  ihlat_bnd(ilon) = 0
        mlt  = ylonm(ilon)*deg2mlt              ! mag.local time ?
        do ilat = nmlat_wei+1,0,-1              ! lat. loop moving torwards pole
	    mlat = 90. - ylatm(ilat)           ! mag. latitude pole = 90 equator = 0
	  
          call gecmp( mlat, mlt, es, ez )	! get electric field
	  
          e_tot = sqrt( es**2 + ez**2 )
          if( abs(e_tot) >= ef_max ) then                        ! e-filed > limit -> boundary
            ihlat_bnd(ilon) = ilat - (ilat - ilat_sft_rvs)/2     ! shift boundary to lat_sft (54deg)
            exit
          end if
        end do
      end do     

!     write(iulog,"('highlat_getbnd: ihlat_bnd=',/,(12i6))") ihlat_bnd

  end subroutine highlat_getbnd

  subroutine bnd_sinus( ihlat_bnd, itrans_width )  
!------------------------------------------------------------------
! Purpose: 
!   1. adjust high latitude boundary (ihlat_bnd) sinusoidally
!   2. width of transition zone from midlatitude potential to high latitude
!      potential (itrans_width)
!
! Method:
! 1.adjust boundary sinusoidally
!   max. wave number to be represented nmax_sin
!   RHS(mi) = Sum_phi Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi)*hlat_bnd(phi) 
!   U(mi,mk)   = Sum_phi Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi) *
!                Sum_(mk=-nmax_sin)^_(mk=nmax_sin) f_mk(phi)
!   single values decomposition of U
!   solving U*LSG = RHS 
!   calculating hlat_bnd:
!   hlat_bnd = Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi)*LSG(mi)
!
! 2. width of transition zone from midlatitude potential to high latitude
!    potential
!    trans_width(phi)=8.-2.*cos(phi) 
!
! Author: A. Maute Nov 2003  am 11/20/03
!------------------------------------------------------------------

      use wamphys_math_interp, only : svdcmp, svbksb
      
!     use sv_decomp, only : svdcmp, svbksb
     
!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(inout) :: ihlat_bnd(0:nmlon)    ! loaction of boundary
      integer, intent(out)   :: itrans_width(0:nmlon) ! width of transition zone

!-----------------------------------------------------------------
! local variables
!-----------------------------------------------------------------
      integer, parameter :: nmax_a = 2*nmax_sin+1 ! absolute array length
      integer, parameter :: ishf   = nmax_sin+1
      integer  :: ilon, i, i1, j, bnd
      real :: sum, mlon
      real :: rhs(nmax_a)
      real :: lsg(nmax_a)
      real :: u(nmax_a,nmax_a)
      real :: v(nmax_a,nmax_a)
      real :: w(nmax_a,nmax_a)
      real :: f(-nmax_sin:nmax_sin,0:nmlon)

!------------------------------------------------------------------
!    Sinusoidal Boundary calculation
!------------------------------------------------------------------
      rhs(:) = 0.
      lsg(:) = 0.
      u(:,:) = 0.
      v(:,:) = 0.
      w(:,:) = 0.

      do ilon = 0,nmlon                  ! long.
        bnd  = nmlath - ihlat_bnd(ilon) ! switch from pole=0 to pole =90
        call ff( ylonm(ilon), nmax_sin, f(-nmax_sin,ilon) )
        do i = -nmax_sin,nmax_sin
	    i1 = i + ishf
          rhs(i1) = rhs(i1) + f(i,ilon) * bnd
!	  write(iulog,*) 'rhs ',ilon,i1,bnd,f(i,ilon),rhs(i+ishf)
          do j = -nmax_sin,nmax_sin 
            u(i1,j+ishf) = u(i1,j+ishf) + f(i,ilon)*f(j,ilon)
!	    write(iulog,*) 'u ',ilon,i1,j+ishf,u(i+ishf,j+ishf)
          end do
        end do
      end do

!     if (debug) write(iulog,*) ' Single Value Decomposition'
      call svdcmp( u, nmax_a, nmax_a, nmax_a, nmax_a, w, v )

!     if (debug) write(iulog,*) ' Solving'
      call svbksb( u, w, v, nmax_a, nmax_a, nmax_a, nmax_a, rhs, lsg )
!      
      do ilon = 0,nmlon  ! long.
!     sum = 0.
	   sum = dot_product( lsg(-nmax_sin+ishf:nmax_sin+ishf),f(-nmax_sin:nmax_sin,ilon) )  
!       do i = -nmax_sin,nmax_sin
!         sum = sum + lsg(i+ishf)*f(i,ilon)  
!       end do
           ihlat_bnd(ilon)    = nmlath - int( sum + .5 )                                ! closest point
           itrans_width(ilon) = &
           int( 8. - 2.*cos( ylonm(ilon)*dtr ) + .5 )/dlatm  ! 6 to 10 deg.
      end do      
!     write(iulog,"('bnd_sinus: ihlat_bnd=',/,(12i6))") ihlat_bnd
!     write(iulog,"('bnd_sinus: itrans_width=',/,(12i6))") itrans_width

  end subroutine bnd_sinus

  subroutine highlat_adjust( pot_highlats, pot_highlat, pot_midlat, ihlat_bnd )
!------------------------------------------------------------------
! Purpose: Adjust mid/low latitude electric potential and high latitude
!          potential such that there are continous across the mid to high 
!          latitude boundary
!
! Method:
! 1. integrate Phi_low/mid(phi,bnd) along the boundary mid to high latitude
! 2. integrate Phi_high(phi,bnd) along the boundary mid to high latitude
! 3. adjust Phi_high by delta =
!    Int_phi Phi_high(phi,bnd) d phi/360. - Int_phi Phi_low/mid(phi,bnd) d phi/360.
!
! Author: A. Maute Nov 2003  am 11/21/03
!------------------------------------------------------------------

!------------------------------------------------------------------
!	... dummy arguments
!------------------------------------------------------------------
      implicit none
      integer, intent(in)     :: ihlat_bnd(0:nmlon) ! boundary mid to high latitude
      real, intent(in)    :: pot_midlat(0:nmlon,0:nmlat) ! low/mid latitude potentail
      real, intent(inout) :: pot_highlat(0:nmlon,0:nmlat)! high_lat potential
      real, intent(inout) :: pot_highlats(0:nmlon,0:nmlat)! high_lat potential! smoothed high_lat potential

!------------------------------------------------------------------
! local:     
!------------------------------------------------------------------
      integer  :: bnd, ilon, ilat, ilatS, ibnd60, ibnd_hl
      real :: pot60, pot_hl, del

!-------------------------------------------------------------------
! 1. integrate Phi_low/mid(phi,bnd) along the boundary mid to high latitude
! 2. integrate Phi_high(phi,bnd) along the boundary mid to high latitude
!-------------------------------------------------------------------
      pot60  = 0.
      pot_hl = 0.
      do ilon = 1,nmlon  ! long.            ! bnd -> eq to pole 0:90
        ibnd60  = nmlat - ihlat_bnd(ilon)   ! 0:180 pole to pole
        ibnd_hl = ihlat_bnd(ilon)           ! colatitude
        pot60   = pot60 + pot_midlat(ilon,ibnd60)
        pot_hl  = pot_hl + pot_highlats(ilon,ibnd_hl)
      end do
      pot60  = pot60/(nmlon)
      pot_hl = pot_hl/(nmlon)
!     print*,'www300',pot60,pot_hl,nmlat_wei,nmlon
      
!     if (debug) write(iulog,*) 'Mid-Latitude Boundary Potential =',
!    &pot60
!     if (debug) write(iulog,*) 'High-Latitude Boundary Potential=',
!    &pot_hl

!-------------------------------------------------------------------
! 3. adjust Phi_high by delta =
!    Int_phi Phi_high(phi,bnd) d phi/360. - Int_phi Phi_low/mid(phi,bnd) d phi/360.
!-------------------------------------------------------------------
      del = pot_hl - pot60

      do ilat = 0,nmlat_wei      ! colatitude
        ilats = nmlat - ilat
        do ilon = 0,nmlon
            pot_highlat(ilon,ilat)   = pot_highlat(ilon,ilat)   - del
            pot_highlat(ilon,ilats)  = pot_highlat(ilon,ilats)  - del
            pot_highlats(ilon,ilat)  = pot_highlats(ilon,ilat)  - del
            pot_highlats(ilon,ilats) = pot_highlats(ilon,ilats) - del
        end do
      end do

  end subroutine highlat_adjust

  subroutine interp_poten( pot_highlats, pot_highlat, pot_midlat, &
        ihlat_bnd, itrans_width ) 
!-------------------------------------------------------------------
! Purpose: construct a smooth global electric potential field 
!
! Method: construct one global potential field
! 1. low/mid latitude: |lam| < bnd-trans_width
!   Phi(phi,lam) = Phi_low(phi,lam)
! 2. high latitude: |lam| > bnd+trans_width
!   Phi(phi,lam) = Phi_hl(phi,lam)
! 3. transition zone: bnd-trans_width <= lam <= bnd+trans_width 
! a. interpolate between high and low/midlatitude potential
!   Phi*(phi,lam) = 1/15*[ 5/(2*trans_width) * {Phi_low(phi,bnd-trans_width)*
!   [-lam+bnd+trans_width] + Phi_hl(phi,bnd+trans_width)*
!   [lam-bnd+trans_width]} + 10/(2*trans_width) {Phi_low(phi,lam)*
!   [-lam+bnd+trans_width] + Phi_hl(phi,lam)*
!   [lam-bnd+trans_width]}]
! b.  Interpolate between just calculated Potential and the high latitude
!    potential in a 3 degree zone poleward of the boundary:
!    bnd+trans_width < lam <= bnd+trans_width+ 3 deg 
!   Phi(phi,lam) = 1/3 { [3-(lam-bnd-trans_width)]* Phi*(phi,lam) +
!   [lam-bnd-trans_width)]* Phi_hl*(phi,lam) }
!
! Author: A. Maute Nov 2003  am 11/21/03      
!------------------------------------------------------------------

!------------------------------------------------------------------
!	... dummy arguments
!------------------------------------------------------------------
      integer, intent(in)  :: ihlat_bnd(0:nmlon)
      integer, intent(in)  :: itrans_width(0:nmlon)
      real, intent(in) :: pot_highlats(0:nmlon,0:nmlat)
      real, intent(in) :: pot_highlat(0:nmlon,0:nmlat)
      real, intent(in) :: pot_midlat(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      real, parameter :: fac = 1./3.
      integer  :: ilon, ilat
      integer  :: ibnd, tw, hb1, hb2, lat_ind
      integer  :: min_ilat 
      integer  :: j1, j2
      real :: a, b, lat, b1, b2
      real :: wrk1, wrk2

!$omp parallel do private(ilat,ilon,ibnd,tw)
      do ilon = 0,nmlon
        ibnd = ihlat_bnd(ilon)     ! high latitude boundary index
	  tw   = itrans_width(ilon)  ! width of transition zone (index)
!-------------------------------------------------------------------
! 1. low/mid latitude: |lam| < bnd-trans_width
!   Phi(phi,lam) = Phi_low(phi,lam)
!-------------------------------------------------------------------
        do ilat = 0,nmlath-(ibnd+tw+1)
          potent(ilon,nmlath+ilat) = pot_midlat(ilon,nmlath+ilat)
          potent(ilon,nmlath-ilat) = pot_midlat(ilon,nmlath+ilat)
        end do
!------------------------------------------------------------------
! 2. high latitude: |lam| > bnd+trans_width
!   Phi(phi,lam) = Phi_hl(phi,lam)
!------------------------------------------------------------------
        do ilat = 0,ibnd-tw-1
          potent(ilon,ilat)       = pot_highlats(ilon,nmlat-ilat)
          potent(ilon,nmlat-ilat) = pot_highlats(ilon,nmlat-ilat)
        end do
      end do
!------------------------------------------------------------------
! 3. transition zone: bnd-trans_width <= lam <= bnd+trans_width 
!------------------------------------------------------------------
! a. interpolate between high and low/midlatitude potential
! update only southern hemisphere (northern hemisphere is copied
! after smoothing)
!------------------------------------------------------------------
!!$omp parallel do private(ilat,ilon,ibnd,tw,a,b,b1,b2,hb1,hb2,
!    &lat_ind,j1,j2,wrk1,wrk2)
      do ilon = 0,nmlon
        ibnd = ihlat_bnd(ilon)          ! high latitude boundary index
	  tw   = itrans_width(ilon)       ! width of transition zone (index)
        a    = 1./(2.*tw)
        b1   = (nmlath - ibnd + tw)*a
        b2   = (nmlath - ibnd - tw)*a
        hb1  = nmlath - (ibnd + tw)
        j1   = nmlath - hb1
        hb2  = nmlath - (ibnd - tw)
        j2   = nmlath - hb2
        if (j2 < 0) j2 = 0              ! Tomoko's fix - j2 >= 0
        wrk1 = pot_midlat(ilon,j1)
        wrk2 = pot_highlats(ilon,j2)
!       write(iulog,*) 'pot_all ',ilon,hb1,hb2,nmlath -ibnd,tw
        min_ilat = ibnd-tw
        if (min_ilat < 0) min_ilat = 0  ! Tomoko's fix
        do ilat = min_ilat,ibnd+tw      ! do ilat = ibnd-tw,ibnd+tw
	    lat_ind = nmlath - ilat
          potent(ilon,ilat) = fac*((wrk1 + 2.*pot_midlat(ilon,ilat))*(b1 - a*lat_ind)  &
         	  + (wrk2 + 2.*pot_highlats(ilon,ilat))*(a*lat_ind - b2))
          potent(ilon,nmlat-ilat) = potent(ilon,ilat)
        end do
!------------------------------------------------------------------
! b.  Interpolate between just calculated Potential and the high latitude
!    potential in a 3 degree zone poleward of the boundary
!------------------------------------------------------------------
	  do ilat = hb2+1,nmlath
	    a = max( 3./dlatm - (ilat - hb2 - 1),0. )
	    b = 3./dlatm - a
          potent(ilon,nmlath-ilat) = (a*potent(ilon,nmlath-ilat) &  
           + b*pot_highlat(ilon,nmlath-ilat))/3.*dlatm
          potent(ilon,nmlath+ilat) = potent(ilon,nmlath-ilat)
        end do
      end do      

  end subroutine interp_poten

  subroutine DerivPotential
!-----------------------------------------------------------------
! Purpose: calulates the electric field [V/m] from the electric potential
!
! Method:  Richmond [1995] eqn 5.9-5.10
! ed1(:,:) = Ed1 = - 1/[R cos lam_m] d PHI/d phi_m
! ed2(:,:) = Ed2 = 1/R d PHI/d lam_m /sin I_m
! R = R_e + h_r we assume a reference height of 130 km which is also
! used in the TIEGCM code
!
! Author: A. Maute Dec 2003  am 12/16/03
!-----------------------------------------------------------------

      integer  :: i, j, ip1f, ip2f, ip3f
      real :: coslm, r, fac, wrk
      real :: wrk1d(0:nmlon)

      r = r_e + h_r  ! earth radius + reference height [m]
!-----------------------------------------------------------------
! ed2= Ed2 is the equatorward/downward component of the electric field, at all 
! geomagnetic grid points (central differencing)
!-----------------------------------------------------------------
      fac = .5/(dlatm*dtr*r)
!$omp parallel do private(j, i, wrk )
      do j = 1,nmlath-1		! southern hemisphere
! idea
        wrk = fac/sinIm_mag(j)
!       wrk = fac
        do i = 0,nmlon
          ed2(i,j) = (potent(i,j+1) - potent(i,j-1))*wrk
        end do
      end do

!$omp parallel do private(j, i, wrk )
      do j = nmlath+1,nmlat-1	! northern hemisphere
        wrk = fac/sinIm_mag(j)
        do i = 0,nmlon
          ed2(i,j) = (potent(i,j+1) - potent(i,j-1))*wrk
        end do
      end do

!-----------------------------------------------------------------------
! Interpolate of ed2 between between -12 <= lam_m <= 12 degrees:
!-----------------------------------------------------------------------
      wrk1d(:) = ed2(:,jmax) - ed2(:,jmin)
      do j = jmin+1,jmax-1
        fac = (ylatm(j) - ylatm(jmin))/(ylatm(jmax) - ylatm(jmin))
        do i = 0,nmlon
	    ed2(i,j) = ed2(i,jmin) + fac*wrk1d(i)
	end do
      end do

!-----------------------------------------------------------------------
! ed1= Ed1 is the zonal component of the electric field, at all 
! geomagnetic grid points (central differencing)
!-----------------------------------------------------------------------
      fac = .5/(dlonm*dtr*r)

      do j = 1,nmlat-1
        coslm = ylatm(j) - 90.
        coslm = cos( coslm*dtr )
        wrk = fac/coslm
        do i = 1,nmlon-1
          ed1(i,j) = -(potent(i+1,j) - potent(i-1,j))*wrk
        end do
        i = 0
        ed1(i,j)  = -(potent(i+1,j) - potent(nmlon-1,j))*wrk
        ed1(nmlon,j) = ed1(i,j)
      end do

!-----------------------------------------------------------------------
! Poles:
!-----------------------------------------------------------------------
      do i = 0,nmlon
        ip1f = i + nmlon/4
        if( ip1f > nmlon ) then
           ip1f = ip1f - nmlon
        end if
        ip2f = i + nmlon/2
        if( ip2f > nmlon ) then
           ip2f = ip2f - nmlon
        end if
        ip3f = i + 3*nmlon/4
        if( ip3f > nmlon ) then
           ip3f = ip3f - nmlon
        end if
        ed1(i,0)=.25*(ed1(i,1)-ed1(ip2f,1)+ed2(ip1f,1)-ed2(ip3f,1))
        ed1(i,nmlat) = .25*(ed1(i,nmlat-1) - ed1(ip2f,nmlat-1)  &
             + ed2(ip1f,nmlat-1) - ed2(ip3f,nmlat-1))
	     
        ed2(i,0)=.25*(ed2(i,1)-ed2(ip2f,1)-ed1(ip1f,1)+ed1(ip3f,1))
	
        ed2(i,nmlat) = .25*(ed2(i,nmlat-1) - ed2(ip2f,nmlat-1)  &
             - ed1(ip1f,nmlat-1) + ed1(ip3f,nmlat-1))
      end do

  end subroutine DerivPotential

!      end module wamphys_efield
 

!================================================================================================

!*********************** Copyright 1996, Dan Weimer/MRC ***********************

!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following routines (translib.for) were added to return the dipole tilt. C
!  GET_TILT was initially a procedure (TRANS), here it has been changed into   C
!  a function which returns the dipole tilt.                                   C
! Barbara Emery (emery@ncar.ucar.edu) and William Golesorkhi, HAO/NCAR (3/96)  C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! COORDINATE TRANSFORMATION UTILITIES
!**********************************************************************        
  subroutine WAM_GET_TILT(YEAR,MONTH,DAY,HOUR, get_tilt)
!
!-----------------------------------------------------------------------
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following line initially was:                                           C
!       SUBR TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)                            C
!  It has been changed to return the dipole tilt from this function call.      C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         
!      THIS SUBR DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
!      TRANSFORMATIONS, IDENTIFIED BY K.
!          K=1 TRANSFORMS GSE to GEO
!          K=2     "      GEO to MAG
!          K=3     "      GSE to MAG
!          K=4     "      GSE to GSM
!          K=5     "      GEO to GSM
!          K=6     "      GSM to MAG
!          K=7     "      GSE to GEI
!          K=8     "      GEI to GEO
!          K=9     "      GSM to SM 
!	   K=10    "      GEO to SM 
!	   K=11    "      MAG to SM 
!
!      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
!      FILE UNIT=IDBUG
!
!      The formal names of the coordinate systems are:
!	GSE - Geocentric Solar Ecliptic
!	GEO - Geographic
!	MAG - Geomagnetic
!	GSM - Geocentric Solar Magnetospheric
!	SM  - Solar Magnetic
!	
!      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
!      ST(I) AND CT(I) ARE SINES & COSINES.       
!
!      Program author:  D. R. Weimer
!
!      Some of this code has been copied from subs which had been
!      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
!      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
!
!      The formulas for the calculation of Greenwich mean sidereal time (GMST)
!      and the sun's location are from "Almanac for Computers 1990",
!      U.S. Naval Observatory.
!
!-----------------------------------------------------------------


        use wam_efield_setdef_data, only: CX,ST,CT,AM,EPOCH, TH0,PH0,DIPOLE

        implicit none 
!
!-----------------------------Return Value--------------------------
!
        real :: get_tilt


        INTEGER :: YEAR, MONTH, DAY
        REAL    :: HOUR
!
!-----------------------------Parameters------------------------------
        INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
        INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
        INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
        INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
        INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG
        
        PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
        PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
        PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
        PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
        PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
        PARAMETER (MAGSM =11,SMMAG =-11)
!---------------------------Local variables-----------------------------
!
        integer IDBUG
        integer j, k, jd, iyr, i, mjd

        REAL UT, T0, GMSTD, GMSTH, ECLIP, MA, LAMD, SUNLON, pi
        real b32, b33, b3
!
!-------------------------External Functions----------------------------
!
!        integer julday_wam
!        external julday_wam
!
!-----------------------------------------------------------------------
!
!       EPOCH=1980.
!       TH0=11.19
!       PH0=-70.76
!       DIPOLE=.30574
        pi=3.141592653
!       pi=2.*ASIN(1.)
!------ NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! IDBUG=0 to prevent printing data to the screen or writing data to a file.    C
        IDBUG = 0
!------ NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	
	JD=JULDAY_WAM(MONTH,DAY,IYR)
	
	MJD=JD-2400001
!     T0=(real(MJD,r8)-51544.5)/36525.0
	T0=(float(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 +36000.770*T0 +3.87933E-4*T0*T0 +15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
      MA=357.528 + 35999.050*T0 + 0.041066678*UT
      CALL ADJUST(MA)
      LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
      CALL ADJUST(LAMD)
      SUNLON=LAMD + (1.915-0.0048*T0)*SIN(MA*pi/180.) + 0.020*SIN(2.*MA*pi/180.)
      CALL ADJUST(SUNLON) 

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
! Derived later:
!       CX(6) = Dipole tilt angle  
!       CX(7) = Angle between sun and magnetic pole
!       CX(8) = Subsolar point latitude
!       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIN(CX(I)*pi/180.)
	  CT(I) = COS(CX(I)*pi/180.)
	ENDDO
!         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
!         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) + &
        AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
!         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
!         
      DO I=1,3   
      DO J=1,3   
         AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) + &
         AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
!         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3  
!         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) + &
          AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) + &
           AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO 
!
	ST(6) = AM(3,1,GSEMAG)        
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIN(ST(6)*pi/180.)  

      AM(1,1,GSMSM) = CT(6)
      AM(1,2,GSMSM) = 0.
      AM(1,3,GSMSM) = -ST(6)
      AM(2,1,GSMSM) = 0.
      AM(2,2,GSMSM) = 1.
      AM(2,3,GSMSM) = 0.
      AM(3,1,GSMSM) = ST(6)
      AM(3,2,GSMSM) = 0.
      AM(3,3,GSMSM) = CT(6)  
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +  &
           AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +  &
          AM(I,2,GSMSM)*AM(J,2,GSMMAG) +AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
      
!
      CX(7)=ATAN2( AM(2,1,11) , AM(1,1,11) )
      
      CX(7)=CX(7)*180./pi
      CX(8)=ASIN( AM(3,1,1)*pi/180. )
      CX(9)=ATAN2( AM(2,1,1) , AM(1,1,1) )
      CX(9)=CX(9)*180./pi

      IF(IDBUG.NE.0)THEN
! 
        DO K=1,11
!        WRITE(IDBUG,1001) K
         DO I=1,3
!          WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The next line was added to return the dipole tilt from this function call.  C

      GET_TILT = CX(6)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
  END SUBROUTINE WAM_GET_TILT

!======================================================================

  INTEGER FUNCTION JULDAY_WAM(MM,ID,IYYY)

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


!================================================================================================

  SUBROUTINE GECMP (AMLA,RMLT,ET,EP)
!
!-----------------------------------------------------------------------
!          Get Electric field components for the Weimer electrostatic
!          potential model.  Before use, first load coefficients (CALL
!          READCOEF) and initialize model conditions (CALL SETMODEL).
!
!          INPUTS:
!            AMLA = Absolute value of magnetic latitude (deg)
!            RMLT = Magnetic local time (hours).
!          RETURNS:
!            ET = Etheta (magnetic equatorward*) E field component (V/m)
!            EP = Ephi   (magnetic eastward)     E field component (V/m)
!
!          * ET direction is along the magnetic meridian away from the
!            current hemisphere; i.e., when ET > 0, the direction is
!              southward when RMLA > 0
!              northward when RMLA < 0
!
!          NCAR addition (Jan 97).  R.Barnes
!-----------------------------------------------------------------------


      use wam_efield_setdef_data, only: ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP

      use wamphys_weimer2005, only:EpotVal_new 
      

      implicit none 
!
!
      real amla, rmlt, et, ep
!
!-----------------------------Parameters------------------------------
!
      real d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847)
      PARAMETER (R2D = 57.2957795130823208767981548147)
!
!---------------------------Local variables-----------------------------
!
      real p1, p2
      real xmlt, xmlt1, kpol, dphi, amla1
!
!-----------------------------------------------------------------------
!
      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

!          Calculate -(latitude gradient) by stepping 10 km along the
!          meridian in each direction (flipping coordinates when going
!          over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
        AMLA1 = 180. - AMLA1
        XMLT1 = XMLT1 + 12.
      ENDIF
      call EpotVal_new(AMLA1    , XMLT1, P1 )
      call EpotVal_new(AMLA-STPD, XMLT,  P2 )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

!          Calculate -(lon gradient).  For most latitudes, step along a
!          great circle.  However, limit minimum latitude to the model
!          minimum (distorting the path onto a latitude line).  Also,
!          avoid a divide by zero at the pole avoid by using Art's trick
!          where Ephi(90,lon) = Etheta(90,lon+90)
      IF (AMLA .LT. ALAMX) THEN
        AMLA1 = MAX (ASIN(SIN(AMLA*D2R)*CSTP) , ALAMR)
        DPHI  = ASIN (SSTP/SIN(AMLA1))*R2D
        AMLA1 = AMLA1*R2D
        call EpotVal_new(AMLA1, XMLT+DPHI, P1 )
        call EpotVal_new(AMLA1, XMLT-DPHI, P2 )
      ELSE
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

!          Below model minimum lat, the potential is value at min lat
      IF (AMLA .LT. ALAMN) THEN
        ET = 0.
        EP = EP * COS(ALAMR)/COS(AMLA*D2R)
      ENDIF

  100 RETURN
  END SUBROUTINE GECMP

end module wamphys_efield
