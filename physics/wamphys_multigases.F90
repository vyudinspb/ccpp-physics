!***********************************************************************
!>@brief The module 'multi_gases' peforms multi-gas WAM constitutents.
!>@author H.-M. H. Juang, NOAA/NWS/NCEP/EMC - dycore subs 2020
!>@author Valery. A. Yudin, NASA/GSFC and CUA - wam physics subs 2022
!>@ToDo check 'multi-gases' and water-related substances sphum+1:nwat
!>      make sure the order of "gases" passed to WAM-physics
! the current set and order of species is a) H2O; b) CLW/SNOW/GRAuP; 2-nwat ;c) nwat+1 (O), nwat+2(02) 

module wamphys_multigases

! 

      use machine, only: kind_phys
!      
! below => nonsense of CCPP * not-TODO: MAKE THIS INPUT ARGUMENTS con_rd, con_cp
!          it is a danger to screw up codes
! check also "module gfs_phy_tracer_config" => defines tracer-names spo2, spo, spo3 ????
!                                             "bad" analog for sphum
!!!!#ifdef MULTI_GASES
!        print *,' ++++ ntoz nto nto2 ',ntoz,nto,nto2
!        if(ntoz  > 0) gfs_phy_tracer%vname(ntoz)  = 'spo3'   'o3mr'
!        if(nto   > 0) gfs_phy_tracer%vname(nto)   = 'spo'    'o1mr'
!        if(nto2  > 0) gfs_phy_tracer%vname(nto2)  = 'spo2'   'o2mr'
!!!#else

      use physcons, only : rdgas => con_rd, cp_air => con_cp
!    
!< nwat number of hydrometeors in dcyore (including water vapor)
!    GFS_init_type%nwat; tracer_names(:)
!
      implicit none
      integer              ::  ntrac_wam      
      integer              ::  num_gas
      integer              ::  ind_gas, num_wat
      integer, parameter   ::  ind_h2o = 1
      integer, parameter   ::  ind_n2  = 0 
             
      integer              ::  ind_o1, ind_o2, ind_o3   
          
      real(kind=kind_phys), allocatable ::  ri(:),cpi(:)
      integer, allocatable              ::  ind_wamtr(:)

      real(kind_phys), allocatable :: vir(:)
      real(kind_phys), allocatable :: vicp(:)
      real(kind_phys), allocatable :: vicv(:)


      CONTAINS
      
    subroutine wamphys_set_major_tracers (me, master, nlunit, fn_nml2, ntrac, nto1, nto2, nto3, ntqv) 
      implicit none 
      integer, intent(in) ::    master,  me, nlunit
      character(len=64), intent (inout) :: fn_nml2
      integer, intent(in) ::   ntrac, nto1, nto2, nto3, ntqv     
      integer             ::    ntrac_nml, ntrac_wat      
      integer             ::    i
      integer :: ios
      logical :: exists
!ncnst =ntrac nwat-last-water-based  
    
      namelist /wamphys_tracer_cpi/ ri,cpi, ntrac_wat, ntrac_nml  
!      write(6,*) 'inmultigases  0000'
	 
!         write(6,*)   'fn_nml2 in multigases 00000', trim(fn_nml2)      
    inquire (file =trim (fn_nml2) , exist = exists)
    if (.not. exists) then
        write(6,'(3a)') 'wamphys_init: namelist file: ', trim (fn_nml2), ' does not exist'
        
	if (me  ==  master) print *, ' wamphys_set_major_tracers wamphys_nml error=', trim (fn_nml2)
!        return      
    endif 
    
      if (.not. allocated(ri)) then
        allocate( ri(0:ntrac))
        allocate(cpi(0:ntrac))  
        allocate(ind_wamtr(0:ntrac))   
      endif
      
 
		     
!=========================================
! default cpi/ri for: N2, H2O, O3, O, O2
! open(nlun_ion, file=trim(nml_ion), status='old' )
!=========================================      
! 
       if (me  ==  master)   write(6,*) 'inmultigases'
	 
!         write(6,*)   'fn_nml2 in multigases 1111', trim(fn_nml2) 
	open (unit = nlunit, file = 'input.nml', action = 'read', status = 'old', iostat = ios)
!        write(6,*)   'fn_nml2 in multigases', trim(fn_nml2) 
	rewind (nlunit)
        read(nlunit, wamphys_tracer_cpi)   
        close  (nlunit)
	
      	if ( ntrac_nml .ne. ntrac ) then
	  write(6,*)  ' wamphys_set_major_tracers dim-n of ntrac =/=  ntrac_nml '
	  write(6,*)  '	ntrac_nml = ', ntrac_nml
	  write(6,*)  '	ntrac_wam = ', ntrac  
	  if ( me == master ) then 
	    do i=0,ntrac_nml 
	       if (ri(i) > 0) write(6,*)  i, ' ri =', ri(i), ' cpi= ', cpi(i)  
	    enddo
	  endif
	endif
	
	  if ( me == master ) then 
	  write(6,*)  ' wamphys_set_major_tracers dim-n of ntrac_nml/nwat ', ntrac_nml, ntrac_wat 	  
	    do i=0,ntrac_nml 
	       write(6,*)  i, ' ri =', ri(i), ' cpi= ', cpi(i)  
	    enddo
	  endif	
	  
        ind_wamtr(ind_n2)  = 0
        ind_wamtr(ind_h2o) = ntqv
        ind_o3 = ntrac_wat+1
        ind_wamtr(ind_o3) = ind_o3 
        ind_o1 =  ind_wamtr(ind_o3)+1
        ind_o2 =  ind_o1+1      
        ind_wamtr(ind_o1) = ind_o1
        ind_wamtr(ind_o2) = ind_o2
	num_wat = ntrac_wat	  
      end  subroutine wamphys_set_major_tracers    
! --------------------------------------------------------
      subroutine multi_gases_init(ngas, nwat, ri, cpi, is_master)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Ouput: vir(i): ri/rdgas - r0/rdgas
!        vir(0): r0/rdgas
!        vicp(i): cpi/cp_air - cp0/cp_air
!        vicp(0): cp0/cp_air
!        cv_air = cp_air - rdgas
!        vicv(i): cvi/cv_air - cv0/cv_air
!        vicv(0): cv0/cv_air
!--------------------------------------------
      integer, intent(in):: ngas, nwat
      real(kind=kind_phys), intent(in):: ri(0:ngas)
      real(kind=kind_phys), intent(in):: cpi(0:ngas)
      logical, intent(in):: is_master
! Local:
      integer :: n, sphum, sphump1
      real (kind=kind_phys) ::  cvi(0:ngas)
      real(kind=kind_phys)  ::  cv_air
      logical :: default_gas=.false.

      sphum = 1
      sphump1 = sphum+1
!     
      ind_gas = nwat+1
      
      do n=0,ngas
        if( ri(n).ne.0.0 .or. cpi(n).ne.0.0 ) num_gas=n
      enddo
      
      if ( num_gas.eq.1 ) default_gas=.true.
      
      allocate( vir (0:num_gas) )
      allocate( vicp(0:num_gas) )
      allocate( vicv(0:num_gas) )

      cv_air = cp_air - rdgas
      do n=0,num_gas
        cvi(n) = cpi(n) - ri(n)
      enddo

      vir (0) =  ri(0)/rdgas
      vicp(0) = cpi(0)/cp_air
      vicv(0) = cvi(0)/cv_air
      if( default_gas ) then
        vir (0) = 1.0
        vicp(0) = 1.0
        vicv(0) = 1.0
      endif
      do n=1,num_gas
        vir(n) = 0.0
        if(  ri(n).gt.0.0 ) vir (n) =  ri(n)/rdgas -  vir (0)
        vicp(n) = 0.0
        if( cpi(n).gt.0.0 ) vicp(n) = cpi(n)/cp_air - vicp(0)
        vicv(n) = 0.0
        if( cvi(n).gt.0.0 ) vicv(n) = cvi(n)/cv_air - vicv(0)
      enddo

      if( is_master ) then
        write(*,*) ' ccpp multi_gases_init with ind_gas=',ind_gas
        write(*,*) ' ccpp multi_gases_init with num_gas=',num_gas
        write(*,*) ' ccpp multi_gases_init with vir =',vir
        write(*,*) ' ccpp multi_gases_init with vicp=',vicp
        write(*,*) ' ccpp multi_gases_init with vicv=',vicv
      endif

      return
      end subroutine multi_gases_init
! ----------------------------------------------------------------

! ----------------------------------------------------------------
      subroutine wamphys_multigases_finalize()

      if(allocated(ri )) deallocate(ri )
      if(allocated(cpi)) deallocate(cpi)
      if(allocated(ind_wamtr)) deallocate(ind_wamtr)

      return
      end subroutine wamphys_multigases_finalize
! ----------------------------------------------------------------

end module wamphys_multigases
