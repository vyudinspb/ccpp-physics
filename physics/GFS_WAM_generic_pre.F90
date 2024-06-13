!> \file GFS_WAM_generic_pre.F90 
!! This file contains the CCPP-compliant WAM-PHYS pre
!! This file can be removed since in "GFS_suite_interstitial_1" d?dt = 0;
!! interstitial codes.
module GFS_WAM_generic_pre

contains

!> \section arg_table_GFS_WAM_generic_pre_run Argument Table
!! \htmlinclude GFS_WAM_generic_pre_run.html
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!>  @{
      subroutine GFS_WAM_generic_pre_run(lssav, ldiag3d, qdiag3d, dtf,  dudt, dvdt, dtdt, dqdt,  &
      & flag_for_wam_generic_tend, dtend, dtidx, index_of_temperature, index_of_x_wind,  &
      & index_of_y_wind, index_of_process_wamphys, ntoz,nto1,nto2, do_wamphys_diag,errmsg, errflg)

      use machine, only : kind_phys
      implicit none
      
      logical, intent(in) :: lssav, ldiag3d, qdiag3d, flag_for_wam_generic_tend
      logical, intent(in) :: do_wamphys_diag     !< flag for wamphys Diagnostics

      real(kind=kind_phys), intent(in) :: dudt(:,:), dvdt(:,:), dtdt(:,:)
      real(kind=kind_phys), intent(in) :: dqdt(:,:,:)
      real(kind=kind_phys), intent(in) :: dtf
      

      ! dtend only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &  index_of_x_wind, index_of_y_wind, index_of_process_wamphys
      integer, intent(in) :: ntoz,nto1,nto2
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: idtend

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lssav) then
     
        if (ldiag3d .and. flag_for_wam_generic_tend.and. .not.do_wamphys_diag) then
          idtend = dtidx(index_of_temperature, index_of_process_wamphys)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) - dtdt*dtf
          endif

          idtend = dtidx(index_of_x_wind, index_of_process_wamphys)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) - dudt*dtf
          endif

          idtend = dtidx(index_of_y_wind, index_of_process_wamphys)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) - dvdt*dtf
          endif
        endif
        if (qdiag3d .and. flag_for_wam_generic_tend) then
          idtend = dtidx(100+ntoz, index_of_process_wamphys)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) - dqdt(:,:,ntoz)*dtf
          endif

          idtend = dtidx(100+nto1, index_of_process_wamphys)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) - dqdt(:,:,nto1)*dtf
          endif

          idtend = dtidx(100+nto2, index_of_process_wamphys)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dqdt(:,:,nto2)*dtf
          endif
        endif
      endif

    end subroutine GFS_WAM_generic_pre_run
!> @}

end module GFS_WAM_generic_pre
