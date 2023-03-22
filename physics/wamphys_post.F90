!>  \file wamphys_post.F90
!! This file contains
module wamphys_post

contains

!>\defgroup ugwpv1_gsldrag_post ugwpv1_gsldrag Scheme Post
!! @{

      subroutine wamphys_post_init ()
      end subroutine wamphys_post_init

!> \section arg_table_wamphys_post_run Argument Table
!! \htmlinclude wamphys_post_run.html
!!
!     subroutine wamphys_post_run ( im, levs, do_wamphys_diag, dtf,                &
!         dudt_wamph,  dvdt_wamph,  dtdt_wamph,  do1dt_wamph,  do2dt_wamph,   &
!         dudt_iwamph, dvdt_iwamph, dtdt_iwamph, do1dt_iwamph, do2dt_iwamph,    &	 
!	 dtdt, dudt, dvdt, errmsg, errflg)

     subroutine wamphys_post_run ( im, levs, do_wamphys_diag, dtf,           &
         dudt_wamph,  dvdt_wamph,  dtdt_wamph,  do1dt_wamph,  do2dt_wamph,   &
         dudt_iwamph, dvdt_iwamph, dtdt_iwamph, do1dt_iwamph, do2dt_iwamph,  &
         dudt_wammd,dvdt_wammd,dtdt_wammd,                                    &   !Add more diag vars wchen
         dudt_wamion,dvdt_wamion,dtdt_wamion,                                 &
         dtdt_wamrad,                                                         &	 	 
         du3dt_wammd,dv3dt_wammd,dt3dt_wammd,                                 &   
         du3dt_wamion,dv3dt_wamion,dt3dt_wamion,                              &   
         dt3dt_wamrad,lssav, ldiag3d, qdiag3d,dtend,                          &
         flag_for_wam_generic_tend, dtidx, index_of_temperature,              &
         index_of_x_wind,index_of_y_wind, index_of_process_wamphys,           &
         index_of_process_wammd,index_of_process_wamrad,                      &
         index_of_process_wamion,ntoz,nto1,nto2,                              & 	 !Add more diag vars wchen
	       dtdt, dudt, dvdt,dqdt, errmsg, errflg)



        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: do_wamphys_diag     !< flag for wamphys Diagnostics
        !Add more diag vars wchen
	      logical, intent(in) :: lssav, ldiag3d, qdiag3d, flag_for_wam_generic_tend
        real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
        integer, intent(in) :: dtidx(:,:)
        integer, intent(in) :: index_of_temperature,index_of_x_wind, index_of_y_wind
        integer, intent(in) :: index_of_process_wamphys
        integer, intent(in) :: index_of_process_wammd,index_of_process_wamrad,index_of_process_wamion
        integer, intent(in) :: ntoz,nto1,nto2
        real(kind=kind_phys), intent(in),    dimension(:,:) ::dudt_wammd,dvdt_wammd,dtdt_wammd
        real(kind=kind_phys), intent(in),    dimension(:,:) ::dudt_wamion,dvdt_wamion,dtdt_wamion
        real(kind=kind_phys), intent(in),    dimension(:,:) ::dtdt_wamrad
        real(kind=kind_phys), intent(inout),    dimension(:,:) ::du3dt_wammd,dv3dt_wammd,dt3dt_wammd
        real(kind=kind_phys), intent(inout),    dimension(:,:) ::du3dt_wamion,dv3dt_wamion,dt3dt_wamion
        real(kind=kind_phys), intent(inout),    dimension(:,:) ::dt3dt_wamrad

        real(kind=kind_phys), intent(in),    dimension(:,:) :: dtdt_iwamph, dvdt_iwamph
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dudt_iwamph
        real(kind=kind_phys), intent(in),    dimension(:,:) :: do1dt_iwamph, do2dt_iwamph 	
 
        real(kind=kind_phys), intent(inout),    dimension(:,:) :: dtdt_wamph, dvdt_wamph
        real(kind=kind_phys), intent(inout),    dimension(:,:) :: dudt_wamph
        real(kind=kind_phys), intent(inout),    dimension(:,:) :: do1dt_wamph, do2dt_wamph 	

        real(kind=kind_phys), intent(inout), dimension(:,:)    :: dtdt, dudt, dvdt
        real(kind=kind_phys), intent(inout), dimension(:,:,:)  :: dqdt
        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg
	      integer :: idtend         !Add more diag vars wchen


        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0
!
! post creates the "time-averaged" diagnostics"
!

    if (do_wamphys_diag) then   
          dtdt_wamph  = dtdt_wamph  + dtf *dtdt_iwamph
          dudt_wamph  = dudt_wamph  + dtf *dudt_iwamph
          dvdt_wamph  = dvdt_wamph  + dtf *dvdt_iwamph
          do1dt_wamph = do1dt_wamph + dtf *do1dt_iwamph
          do2dt_wamph = do2dt_wamph + dtf *do2dt_iwamph
          !Add more diag vars wchen
          du3dt_wammd = du3dt_wammd + dtf*dudt_wammd
          dv3dt_wammd = dv3dt_wammd + dtf*dvdt_wammd
          dt3dt_wammd = dt3dt_wammd + dtf*dtdt_wammd
          du3dt_wamion = du3dt_wamion + dtf*dudt_wamion
          dv3dt_wamion = dv3dt_wamion + dtf*dvdt_wamion
          dt3dt_wamion = dt3dt_wamion + dtf*dtdt_wamion
          dt3dt_wamrad = dt3dt_wamrad + dtf*dtdt_wamrad
        
        if (lssav) then
     
          if (ldiag3d .and. flag_for_wam_generic_tend) then
            idtend = dtidx(index_of_temperature, index_of_process_wamphys)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dtdt_iwamph*dtf
            endif
  
            idtend = dtidx(index_of_x_wind, index_of_process_wamphys)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dudt_iwamph*dtf
            endif
  
            idtend = dtidx(index_of_y_wind, index_of_process_wamphys)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dvdt_iwamph*dtf
            endif

            idtend = dtidx(index_of_x_wind, index_of_process_wammd)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dudt_wammd*dtf
            endif
            idtend = dtidx(index_of_y_wind, index_of_process_wammd)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dvdt_wammd*dtf
            endif
            idtend = dtidx(index_of_temperature, index_of_process_wammd)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dtdt_wammd*dtf
            endif

            idtend = dtidx(index_of_x_wind, index_of_process_wamion)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dudt_wamion*dtf
            endif
            idtend = dtidx(index_of_y_wind, index_of_process_wamion)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dvdt_wamion*dtf
            endif
            idtend = dtidx(index_of_temperature, index_of_process_wamion)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dtdt_wamion*dtf
            endif

            idtend = dtidx(index_of_temperature, index_of_process_wamrad)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + dtdt_wamrad*dtf
            endif
          endif
          if (qdiag3d .and. flag_for_wam_generic_tend) then
           
            idtend = dtidx(100+nto1, index_of_process_wamphys)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + do1dt_iwamph*dtf
            endif
  
            idtend = dtidx(100+nto2, index_of_process_wamphys)
            if(idtend>=1) then
              dtend(:,:,idtend) = dtend(:,:,idtend) + do2dt_iwamph*dtf
            endif
          endif
        endif
      endif
!=====================================================================
! Updates inside wamphys_run
!
!        dtdt = dtdt_iwamph+dtdt
!        dudt = dudt_iwamph+dudt
!        dvdt = dvdt_iwamph+dvdt
! or 
! zero tendencies for before the next chains of physics
!        dtdt =0.
!        dudt =0.
!        dvdt =0.
!=====================================================================
      end subroutine wamphys_post_run      

      subroutine wamphys_post_finalize ()
      end subroutine wamphys_post_finalize

!! @}
end module wamphys_post
