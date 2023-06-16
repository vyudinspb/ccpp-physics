!--------------------------------------
! needs to be rewritten during 'wamphys_init' due to identical PRESSURE-grid for
! IPE-WAM domain:    xpk_low, xpk_high = constant & lev_low = constant
!                    plow, phigh = constant =/=f(im)
!--------------------------------------
! someone Weiyu/2017 => WAM/src/share
!      MODULE module_IPE_to_WAM
!      INTEGER :: lowst_ipe_level
!      LOGICAL :: ipe_to_wam_coupling
!      REAL, parameter :: dx=1.0, rdx=1.0/dx
!      REAL, DIMENSION(:, :, :), POINTER :: ZMT, MMT, JHR, SHR, O2DR
!======================================
module wamphys_module_IPE_to_WAM
      
      use machine,                 only: kind_phys
      IMPLICIT none
      
      real(kind=kind_phys), parameter   :: dx =0.25,  rdx = 1./dx
      integer                           :: low_ipe_level
      
end module wamphys_module_IPE_to_WAM
      
SUBROUTINE wamphys_par_for_ipemerge(gzmt,im, levs, lowst_ipe_level, prsl, plow, phigh, &
                                          xpk_low, xpk_high)
	   
       use machine,                 only: kind_phys
! Get all necessray vertical related parameters to merge the IPE 
! back coupling arrays into the WAM related arrays.
!---------------------------------------------------------------
        USE wamphys_module_IPE_to_WAM, only: dx, rdx, low_ipe_level
        IMPLICIT none
        INTEGER              :: im, lowst_ipe_level, levs
        real(kind=kind_phys) :: gzmt(im, lowst_ipe_level:levs)
        real(kind=kind_phys) :: prsl(im, levs) 
        INTEGER              :: plow(im), phigh(im)
        real(kind=kind_phys) :: xpk_low(im), xpk_high(im)	
	
        real(kind=kind_phys)   :: alowipe
        INTEGER                :: i, k

! Dtermine the number of the lowest level in IPE coupling arrays.
! Pre-set up the default values in the IPE arrays are 1E-50.
!----------------------------------------------------------------
        do i =1, im
          low_ipe_level = levs - 2
          do k = lowst_ipe_level, LEVS
            alowipe = max(-1.E30, gzmt(i, k))
            IF(alowipe > -1.E29) THEN
              low_ipe_level = k
              exit
            END IF
          end do
1002      CONTINUE

! Get the merging lower and higher levels.
!-----------------------------------------
          plow(i)    = low_ipe_level-1           ! pressure level above which merging is done
          xpk_low(i) = alog(prsl(i, plow(i)))    ! its log-pressure

! log-pressure below which merging is done, not necessarily a model level,
! note also that log-pressures are < 0
!-------------------------------------------------------------------------
          xpk_high(i) = xpk_low(i) - dx
!
          do k = low_ipe_level, levs
             if (prsl(i,k) < exp(xpk_high(i))) then
                phigh(i) = k - 1
                exit
             endif
          enddo
1003      CONTINUE
        enddo
END SUBROUTINE wamphys_par_for_ipemerge

SUBROUTINE wamphys_merge_ipe2wam(array_ipe, array_wam,               &
           im, levs, lev_low, prsl, plow, phigh, xpk_low, xpk_high)

! This subroutine is to take care the merging calcualtions
! for all six IPE back coupling variable arrays to merge into
! the corresponding WAM model arrays.
      use machine,           only: kind_phys
      USE wamphys_module_IPE_to_WAM, only: rdx
      IMPLICIT NONE

! September 2017 Weiyu Yang coding.
! September 2017 Rashid Akmaev, some changes suggested
!-----------------------------------------------------
      real(kind=kind_phys), DIMENSION(im, lev_low:levs) :: array_ipe
      real(kind=kind_phys), DIMENSION(im, levs)         :: array_wam
      real(kind=kind_phys), DIMENSION(im, levs)         :: prsl
      real(kind=kind_phys), DIMENSION(im)               :: xpk_low, xpk_high
      REAL :: xk

      INTEGER, DIMENSION(im) :: plow, phigh
      INTEGER :: i, k, lev_low, levs, ix, im

      do i=1,im
!
! dx-domain of merging, actually it is constant for all horizotal points
!      
        do k = plow(i)+1, phigh(i)
          xk  = alog(prsl(i,k))
          array_wam(i, k) = (array_ipe(i,k) * (xpk_low(i) - xk) + &
            array_wam(i,k) * (xk - xpk_high(i))) * rdx
        enddo
	
        do k = phigh(i) + 1, levs
           array_wam(i, k) = array_ipe(i, k)
        end do
	
      enddo

END SUBROUTINE wamphys_merge_ipe2wam
