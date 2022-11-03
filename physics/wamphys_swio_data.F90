module wamphys_swio_data

  use machine, only: kind_phys
  
  use wamphys_swdef     , only :  swio_file, sw_ntd, sw_time, sw_days_since
  
  use wamphys_swdef     , only :  sw_f107, sw_f107d,  sw_kp, sw_kpa,           &
                                  sw_nhp,   sw_nhpi,  sw_shp,sw_shpi,          & 
			          sw_den, sw_ang,  sw_bz, sw_bt, sw_vel, sw_ap,sw_apa
			   
!  use wamphys_swdef     , only :  csw_f107, csw_f107d,  csw_kp, csw_kpa, csw_ap, csw_apa, &   
!                           csw_nhp,   csw_nhpi,  csw_shp,   csw_shpi, csw_den, csw_ang, &  
!                           csw_bz,    csw_bt, csw_vel, csw_time
			   
			     
!
implicit none
   integer           :: ncal
   integer           :: it1, it2  

   
   public :: read_swio_data, swio_time_interp 

contains
  
   subroutine read_swio_data(me, master, errmsg, errflg)
  
    use  netcdf
    integer, intent(in) ::  me, master   
    integer :: ncid,  iernc, vid, dimid, status         
    integer :: k
    
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg    
!
      
      iernc=NF90_OPEN(trim(swio_file), nf90_nowrite, ncid)
     
       if(iernc.ne.0) then         
          write(errmsg,'(*(a))') "read_swio_data: cannot open swio_file data-file ",  &
                                    trim(swio_file)
	    print *, 'cannot open swio_file=',trim(swio_file)			    
          errflg = 1
          return
        else


       status = nf90_inq_dimid(ncid, "time", DimID)
!      if (status /= nf90_noerr) call handle_err(status)
!
       status = nf90_inquire_dimension(ncid, DimID,  len = sw_ntd)
       
       status = nf90_inq_dimid(ncid, "calformat", DimID)
       status = nf90_inquire_dimension(ncid, DimID,  len =ncal )
       
           if (me == master)  print *, sw_ntd, ncal, ' dimd of swio_file '
	   if (sw_ntd .le. 0 .or. ncal .le. 0) then 
	       print *, 'swio-file=',    trim(swio_file)	   
	       print *, 'sw_ntd=',sw_ntd, 'ncal=',ncal
	       stop
	   endif
	   	   
        if (.not.allocated(sw_f107))  allocate (sw_f107(sw_ntd ))
        if (.not.allocated(sw_time))  allocate (sw_time(sw_ntd ))         
        if (.not.allocated(sw_days_since))    allocate (sw_days_since(ncal))
	if (.not.allocated(sw_f107d))  allocate (sw_f107d(sw_ntd ))
        if (.not.allocated(sw_kp))  allocate (sw_kp(sw_ntd))  
	if (.not.allocated(sw_kpa))  allocate (sw_kpa(sw_ntd))
        if (.not.allocated(sw_nhp))  allocate (sw_nhp(sw_ntd ))  
	if (.not.allocated(sw_nhpi))  allocate (sw_nhpi(sw_ntd ))
        if (.not.allocated(sw_shp))  allocate (sw_shp(sw_ntd ))  
	if (.not.allocated(sw_shpi))  allocate (sw_shpi(sw_ntd ))
        if (.not.allocated(sw_den))  allocate (sw_den(sw_ntd))  
	if (.not.allocated(sw_ang))  allocate (sw_ang(sw_ntd ))
        if (.not.allocated(sw_bz))  allocate (sw_bz(sw_ntd ))  
	if (.not.allocated(sw_bt))  allocate (sw_bt(sw_ntd))  
	if (.not.allocated(sw_vel))  allocate (sw_vel(sw_ntd ))
        if (.not.allocated(sw_ap))  allocate (sw_ap(sw_ntd))  
	if (.not.allocated(sw_apa))  allocate (sw_apa(sw_ntd ))
	
	
	
	
	
	      
	iernc=nf90_inq_varid( ncid, 'time', vid )
        iernc= nf90_get_var( ncid, vid, sw_time)
	
	iernc=nf90_inq_varid( ncid, 'y4d3h2m2s2', vid )
        iernc= nf90_get_var( ncid, vid, sw_days_since)
		
	iernc=nf90_inq_varid( ncid, 'f107', vid )
        iernc= nf90_get_var( ncid, vid, sw_f107)
	iernc=nf90_inq_varid( ncid, 'kp', vid )
        iernc= nf90_get_var( ncid, vid, sw_kp)
	iernc=nf90_inq_varid( ncid, 'f107d', vid )
        iernc= nf90_get_var( ncid, vid, sw_f107d)
	iernc=nf90_inq_varid( ncid, 'kpa', vid )
        iernc= nf90_get_var( ncid, vid, sw_kpa)
	iernc=nf90_inq_varid( ncid, 'nhp', vid )
        iernc= nf90_get_var( ncid, vid, sw_nhp)
	iernc=nf90_inq_varid( ncid, 'nhpi', vid )
        iernc= nf90_get_var( ncid, vid, sw_nhpi)
	iernc=nf90_inq_varid( ncid, 'shp', vid )
        iernc= nf90_get_var( ncid, vid, sw_shp)
	iernc=nf90_inq_varid( ncid, 'shpi', vid )
        iernc= nf90_get_var( ncid, vid, sw_shpi)
	iernc=nf90_inq_varid( ncid, 'swbt', vid )
        iernc= nf90_get_var( ncid, vid, sw_bt)
	iernc=nf90_inq_varid( ncid, 'swang', vid )
        iernc= nf90_get_var( ncid, vid, sw_ang)
	iernc=nf90_inq_varid( ncid, 'swvel', vid )
        iernc= nf90_get_var( ncid, vid, sw_vel)
	iernc=nf90_inq_varid( ncid, 'swbz', vid )
        iernc= nf90_get_var( ncid, vid, sw_bz)
	iernc=nf90_inq_varid( ncid, 'swden', vid )
        iernc= nf90_get_var( ncid, vid, sw_den)
	iernc=nf90_inq_varid( ncid, 'ap', vid )
        iernc= nf90_get_var( ncid, vid, sw_ap)
	iernc=nf90_inq_varid( ncid, 'apa', vid )
        iernc= nf90_get_var( ncid, vid, sw_apa)
	
	
	
	
			
	iernc=nf90_close(ncid)
	
	endif    
	
  end  subroutine read_swio_data 
      
    subroutine swio_time_interp(me, master, im, idate, fhour,             &
	      csw_f107, csw_f107d,  csw_kp, csw_kpa, csw_ap, csw_apa,     &
	      csw_nhp,  csw_nhpi,  csw_shp,   csw_shpi, csw_den, csw_ang, &
	      csw_bz,   csw_bt,    csw_vel,   csw_time, sw_ktprev) 
	         
    use machine, only: kind_phys	           
    implicit none
    
!input    
    integer, intent(in)               :: me, master
    integer, intent(in)               :: im, idate(4)
    real(kind=kind_phys), intent(in)  :: fhour
      
    real(kind=kind_phys),intent(out)  :: csw_f107, csw_f107d,  csw_kp, csw_kpa, csw_ap, csw_apa  
         real(kind=kind_phys),intent(out)  :: csw_nhp,  csw_nhpi,  csw_shp,   csw_shpi, csw_den, csw_ang 
         real(kind=kind_phys),intent(out)  :: csw_bz,   csw_bt,    csw_vel,   csw_time  	 
         integer,              intent(out) ::	sw_ktprev     
!locals

    integer :: i, j1, j2, it1, it2 , kmin
    integer :: ddd, istsw   
    real(kind=kind_phys)  :: tx1, tx2, w1, w2, fmin 
!
! define day of year ddd ..... from the old-fashioned "GFS-style"
! sw_days_since(ncal): y4 ddd h2 m2 s2
!
         call swio_idate_calendar(idate, fhour, ddd, fmin)  
    	 csw_time = fmin
	 
            it1 = 2
	   istsw = min(1, sw_ntd)   !istsw = min(sw_ktprev, sw_ntd)
         do kmin=istsw, sw_ntd
	    if (fmin .lt. sw_time(kmin) ) then
	    it2 = kmin
	    exit
	    endif
	 enddo
	 
	 it2 = min(it2, sw_ntd)	 
	 it1 = max(it2-1,1)
	 if (it2 > sw_ntd ) then
	  print *, ' Error in time-interpolation for sw_interp '	 
	  print *, ' it1, it2, sw_ntd ', it1, it2, sw_ntd
	  print *, ' Error in time-interpolation see wamphys_swio_data.F90 '	  
	  stop
	 endif
!	 
	
	 w2 = (csw_time-sw_time(it1))/(sw_time(it2)-sw_time(it1))
	 w1 = 1.0-w2     
	 
       csw_f107 = sw_f107(it1)*w1 +sw_f107(it2)*w2
       csw_f107d = sw_f107d(it1)*w1 +sw_f107d(it2)*w2 
       csw_kpa = sw_kpa(it1)*w1 +sw_kpa(it2)*w2        
       csw_kp = sw_kp(it1)*w1 +sw_kp(it2)*w2    
       csw_nhp = sw_nhp(it1)*w1 +sw_nhp(it2)*w2 
       csw_nhpi= sw_nhpi(it1)*w1 +sw_nhpi(it2)*w2 
       csw_shp = sw_shp(it1)*w1 +sw_shp(it2)*w2 
       csw_shpi= sw_shpi(it1)*w1 +sw_shpi(it2)*w2 
       csw_den = sw_den(it1)*w1 +sw_den(it2)*w2 
       csw_ang = sw_ang(it1)*w1 +sw_ang(it2)*w2 
       csw_bz = sw_bz(it1)*w1 +sw_bz(it2)*w2 
       csw_bt = sw_bt(it1)*w1 +sw_bt(it2)*w2 
       csw_vel = sw_vel(it1)*w1 +sw_vel(it2)*w2 
       
       csw_ap = sw_ap(it1)*w1 +sw_ap(it2)*w2 
       csw_apa = sw_apa(it1)*w1 +sw_apa(it2)*w2 
       
       sw_ktprev  =it1
   
    end subroutine swio_time_interp  
    
    subroutine swio_idate_calendar(idate, fhour, ddd, fmin) 
    
    use machine, only: kind_phys    		 
    implicit none  
! input     
    integer, intent(in)                 :: idate(4)
    real(kind=kind_phys), intent(in)   :: fhour
!out    
    integer, intent(out)                :: ddd    
    real(kind=kind_phys),  intent(out)  :: fmin 
!    
!locals
!
      real(kind=kind_phys) :: rinc(5), rjday
      integer              :: jdow, jdoy, jday
      real(4)              :: rinc4(5)
      integer              :: w3kindreal, w3kindint
      
      integer ::  iw3jdn
      integer :: jd1, jddd
      
      integer  idat(8),jdat(8)  
          
       
      idat(1:8)    = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc(1:5)    = 0.
      rinc(2) = fhour
!    
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        rinc4 = rinc
        call w3movdat(rinc4, idat,jdat)
      else
        call w3movdat(rinc,  idat,jdat)
      endif           
!     jdate(8)- date and time (yr, mo, day, [tz], hr, min, sec)
!     sw_days_since                 syr, sddd, shr, smin, ssec
!
      jdow = 0
      jdoy = 0
      jday = 0
      
      call w3doxdat(jdat,jdow, ddd, jday)
!      fddd = float(ddd) + jdat(5) / 24. 
       fmin = fhour*60.      
    end  subroutine swio_idate_calendar   
! w3doxdat(jdat,jdow, ddd, jday)    
!input:jdat       integer (8) ncep absolute date and time
!   (year, month, day, time zone,hour, minute, second, millisecond)
!output:
!      jdow       integer day of week (1-7, where 1 is sunday)
!      jdoy       integer day of year (1-366, where 1 is january 1)
!      jday       integer julian day (day number from jan. 1,4713 b.c.)   
end  module wamphys_swio_data
