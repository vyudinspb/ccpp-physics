
      subroutine wamrad_o2_o3(im,levs,cosz,adt,o2_n,o3_n,rho,cp,  &     
      zg,grav,dth)
!
! Apr 06 2012  Henry Juang, initial implement for nems
! Jan 02 2013  Jun Wang,    move o3ini out of column physics
! Nov 29 2016  VAY review/suggestions and SRB/SRC-out
! May 10 2022  Svetala Karol and CUA FV3WAM-ccpp version

      use physcons, only :    avgd => con_avgd  !, amo3=> con_amo3 
      use physcons, only :    bz => con_boltz        
      use wamphys_init_module, only : eff_hugg, eff_chap, eff_herz, eff_hart
      use wamphys_init_module, only : amo3, amo2, rmo3    
        
      use machine,     only : kind_phys 
!
      implicit none
! Argument
      integer, intent(in) :: im               ! number of data points in adt (first dim)
 
      integer, intent(in) :: levs             ! number of pressure levels
      real(kind=kind_phys), intent(in)    :: cosz(im)         ! cos zenith angle
      real(kind=kind_phys), intent(in)    :: adt(im,levs)     !temp(k) 
      real(kind=kind_phys), intent(in)    :: o2_n(im,levs)    ! /m3
      real(kind=kind_phys), intent(in)    :: o3_n(im,levs)    ! /m3
      real(kind=kind_phys), intent(in)    :: rho(im,levs)     ! kg/m3
      real(kind=kind_phys), intent(in)    :: cp(im,levs)      ! J/kg/k
      real(kind=kind_phys), intent(in)    :: zg(im,levs)      ! height (m)
      real(kind=kind_phys), intent(in)    :: grav(im,levs)    ! (m/s2)
      
      real(kind=kind_phys), intent(out)   :: dth(im,levs)     ! heating rate k/s
!
      real(kind=kind_phys) :: hc,fc,dc,hha,fha,dha,hhu,i1,i2,m,dhu,lams,laml,  &              
                              hhz,fhz,dhzo2,dhzo3,hsrb,fsrb,dsrb,ysrb,h1,rodfac
      real(kind=kind_phys) clmo2(levs),clmo3(levs) 
      real(kind=kind_phys)  :: rdzg
      real(kind=kind_phys)  :: tg_vay
      real(kind=kind_phys)  :: fo2_vay, fo3_vay
      real(kind=kind_phys)  :: fc_dc, fha_dha
      real(kind=kind_phys)  :: hu_exp1, hu_exp2, rm, i2_i1
      integer i,k
!
! very ....."dangerous" games with constants !!!!
!
      fc=370.     !J/m2/s
      dc=2.85E-25 !m2
      fc_dc = fc*dc
      fha=5.13  !J/m2/s
      dha=8.7E-22 !m2
      fha_dha = fha*dha
      i1=0.07   !J/m2/s/A
      i2=0.05
      m=0.01273   !/A
      rm = 1./m

      lams=2805.
      laml=3015.
      dhu=1.15e-6    !m2
      fhz=1.5        !J/m2/s
      dhzo2=6.e-28   !m2
      dhzo3=4.e-22   !m2
      fsrb=0.0128    !J/m2/s
      dsrb=2.07e-24  !m2
      ysrb=0.0152

       hu_exp1 = exp(-m*laml)*dhu
       hu_exp2 = exp(-m*lams)*dhu
       i2_i1 = -0.02
!
! rewrite and compute constants in the init_solar (?)
!   exp(-1.*m*laml) in loops
!   mmr*dp => <n>dZ
      dth(:,:)=0.
      fo3_vay = 1000.*avgd/amo3*bz
      fo2_vay = 1000.*avgd/amo2*bz
      do i=1,im
        if(cosz(i).ge.0.) then
          rodfac=35./sqrt(1224.*cosz(i)**2+1.)
          tg_vay = adt(i,levs)/grav(i,levs)
!         clmo2(levs)=1.e3*o2_n(i,levs)*bz*adt(i,levs)*avgd/(grav(i,levs)*amo2)
          clmo2(levs)=o2_n(i,levs)*tg_vay*fo2_vay
          clmo3(levs)=o3_n(i,levs)*tg_vay*fo3_vay
          do k=levs-1,1,-1
           rdzg =  .5*(zg(i,k+1)-zg(i,k))
          clmo2(k)=clmo2(k+1)+(o2_n(i,k+1)+o2_n(i,k))*rdzg                                     
          clmo3(k)=clmo3(k+1)+(o3_n(i,k+1)+o3_n(i,k))*rdzg              
!
          enddo
          clmo2=clmo2*rodfac   !rad path
          clmo3=clmo3*rodfac
!
!very ....."dangerous" games with constants !!!!exp(-1.*m*laml)  i2/i1 etc....  
!            can be restricted to MLT
!
          do k=1,levs
!
!Chappius                                   acc-cy 2% according Strobel 1978
!
            hc=fc_dc*exp(-dc*clmo3(k))
! Hartley
            hha=fha_dha*exp(-dha*clmo3(k))
! Huggins
           
            hhu=rm*(i1+i2_i1*exp(-clmo3(k)*hu_exp1)     &  
            -i2*exp(-clmo3(k)*hu_exp2)) /clmo3(k)
! Herzberg
            hhz=fhz*(dhzo2*o2_n(i,k)+dhzo3*o3_n(i,k))*   &      
              exp(-dhzo2*clmo2(k)-dhzo3*clmo3(k))
	      
!            hhu=(i1+(i2-i1)*exp(-dhu*clmo3(k)*exp(-m*laml))   &    
!            -i2*exp(-1.*dhu*clmo3(k)*exp(-m*lams)))/(m*clmo3(k))
!
! VAY-2016: hsrb 2 times above P > 0.02 orchetrating with SOLAR in idea_solar_heating
!           moved to idea_solar_heating.f
! Comments for Zhu's updates
!vay            h1=sqrt(1.+4.*dsrb*clmo2(k)/(pi*ysrb))
!vay            hsrb=fsrb*dsrb*o2_n(i,k)*exp(-.5*pi*ysrb*(h1-1.))/h1
!           dth(i,k)=((hc+hha+hhu)*o3_n(i,k)+hhz+hsrb)/   &              
!                 (cp(i,k)*rho(i,k))

            dth(i,k)=((hc+hha*eff_hart(k)+hhu)*o3_n(i,k)+hhz)/  &         
                  (cp(i,k)*rho(i,k))
          enddo
        else
          dth(i,1:levs)=0.
        endif
      enddo
      return
      end  subroutine wamrad_o2_o3
      
      subroutine wamrad_o3prof(im,levs,ntrac, nto3, q, am, nair,o3_n)
!      
! extend GFS-O3 in the thermosphere for WAM radiation
!        will be replaced if O3 is a part of model chemistry
!

      use wamphys_init_module, only : k71, o3ra, rmo3
      use machine,             only : kind_phys 
      

      implicit none
! IN/OUT
      integer, intent(in) :: im                    ! number of data points in PE

      integer, intent(in) :: levs                                     ! number of pressure levels
      integer, intent(in) :: ntrac                                    ! number of tracer
      integer, intent(in) :: nto3                                     !   # in q(:,:, nto3) 
        
      real(kind=kind_phys), intent(in)    :: q(im,levs,ntrac)         ! tracers in "mmr"
      real(kind=kind_phys), intent(in)    :: am(im,levs)              ! mixture mol weight kg
      real(kind=kind_phys), intent(in)    :: nair(im,levs)            ! number density  /m3
      real(kind=kind_phys), intent(out)   :: o3_n(im,levs)            ! /m3
! locals
      real(kind=kind_phys):: rate
      integer i,k

      do i=1,im
        do k=1,k71-1
          o3_n(i,k)=q(i,k,nto3)*am(i,k)*nair(i,k)*rmo3
        enddo
          rate=q(i,k71,nto3)/o3ra(k71)
        do k=k71,levs
          o3_n(i,k)=o3ra(k)*rate*am(i,k)*nair(i,k)*rmo3
        enddo
      enddo
      return
      
      print *, ' nto3 ', nto3
      i = 45
      do k=1, levs
        print*, k, o3_n(i,k), q(i,k,nto3)*am(i,k)*nair(i,k)*rmo3
      enddo
      pause 'wamrad_o3prof '
      return
      end subroutine wamrad_o3prof
