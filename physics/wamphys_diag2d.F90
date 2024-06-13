   module wamphys_diag2d
      use machine,  only             :  kind_phys
      real(kind_phys), allocatable   :: no_diag(:,:)
      real(kind_phys), allocatable   :: el_diag(:,:) 
      real(kind_phys), allocatable   :: sped_diag(:,:)
      real(kind_phys), allocatable   :: shal_diag(:,:)
      real(kind_phys), allocatable   :: qno_diag(:,:) 
      real(kind_phys), allocatable   :: qeuv_diag(:,:)   
      real(kind_phys), allocatable   :: srb1_diag(:,:)  
      real(kind_phys), allocatable   :: srb2_diag(:,:)  
      real(kind_phys), allocatable   :: qlya_diag(:,:) 
                   
   contains
   subroutine wamphys_diag2d_init(im, levs)
   integer :: im, levs
   
     allocate(no_diag(im, levs))
     allocate(el_diag(im, levs))
     allocate(sped_diag(im, levs))
     allocate(shal_diag(im, levs))
     allocate(qno_diag(im, levs)) 
     allocate(qeuv_diag(im, levs)) 
     allocate(srb1_diag(im, levs))
     allocate(srb2_diag(im, levs)) 
     allocate(qlya_diag(im, levs)) 
              
   return
   end subroutine wamphys_diag2d_init
end module wamphys_diag2d
