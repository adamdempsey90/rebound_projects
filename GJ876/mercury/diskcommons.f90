      module disk_commons
       use disk_params
       implicit none

       real(pre)::pi
       real(pre),allocatable,dimension(:)::surface_den,sound_speed,legendre_P,dg_r,rh,theta
       real(pre),allocatable,dimension(:,:)::dg_rho,qlm_less,qlm_more

      end module disk_commons
