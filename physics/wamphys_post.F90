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

     subroutine wamphys_post_run ( im, levs, do_wamphys_diag, dtf,                &
         dudt_wamph,  dvdt_wamph,  dtdt_wamph,  do1dt_wamph,  do2dt_wamph,   &
         dudt_iwamph, dvdt_iwamph, dtdt_iwamph, do1dt_iwamph, do2dt_iwamph,    &	 
	 dtdt, dudt, dvdt,errmsg, errflg)



        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: do_wamphys_diag     !< flag for wamphys Diagnostics
	
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dtdt_iwamph, dvdt_iwamph
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dudt_iwamph
        real(kind=kind_phys), intent(in),    dimension(:,:) :: do1dt_iwamph, do2dt_iwamph 	
 
        real(kind=kind_phys), intent(inout),    dimension(:,:) :: dtdt_wamph, dvdt_wamph
        real(kind=kind_phys), intent(inout),    dimension(:,:) :: dudt_wamph
        real(kind=kind_phys), intent(inout),    dimension(:,:) :: do1dt_wamph, do2dt_wamph 	
	
        real(kind=kind_phys), intent(inout), dimension(:,:) :: dtdt, dudt, dvdt

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg

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
