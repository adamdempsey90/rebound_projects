!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MJP: Trivial alterations to prevent name conflicts with Mercury
! au       -->> dg_au
! msun     -->> dg_msun
! dr       -->> annulus
! rmax     -->> dg_rmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module disk_params
       implicit none
       integer,parameter::pre=8

       ! disk grid size

       integer::nr=300
       integer::ntheta=257
       integer,parameter::nrmax=1000
       integer,parameter::nthetamax = 1000

       ! disk structure

       real(pre)::annulus=0.01 ! AU
       real(pre)::angle=90. ! this should be given in degrees.
       real(pre)::power_sig=1.75
       real(pre)::toomre_q=1.5
       real(pre)::reference_r=1 ! AU
       real(pre)::reference_sig=2400 ! cgs
       real(pre)::mstar=1.!msun
       real(pre)::mmw=2.33
       real(pre)::rmin=0.01
       real(pre)::dg_rmax=3
       real(pre)::small_rho=1d-20

       ! gravity params

       real(pre)::pitch=12d0 ! pitch angle in deg
       real(pre)::rcorot=50.
       real(pre)::phase_drift=8d0 ! 2pi/phase_drift
       real(pre)::power_fac=1d0

       integer::mmax=3
       integer::mmin=2

       integer::nlegen=60
       integer,parameter::nlegenmax=60
       integer::use_nlegen=42
       logical::use_grid=.true.

       ! constants

       real(pre)::rgas=8.254e7
       real(pre)::dg_msun=2d33
       real(pre)::dg_au=1.5d13
       real(pre)::Gcgs=6.67d-8

       real(pre),parameter::one=1d0,two=2d0,three=3d0,four=4d0, &
     &  five=5d0,six=6d0,seven=7d0,eight=8d0,nine=9d0,zero=0d0, &
     &  ten=10d0,half=0.5d0

      end module disk_params

