!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MJP: Trivial alterations to prevent name conflicts with Mercury
! au       -->> dg_au
! msun     -->> dg_msun
! dr       -->> annulus
! rmax     -->> dg_rmax
! pi       -->> dg_pi
! r        -->> dg_r
! rho      -->> dg_rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MJP: Minor alterations for passing variables & minimizing calculations
! - Requires rcyl, rsph , z, sin(p) & cos(p) to be explicitly passed in 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module disk_grav
       use disk_params
       use derived_types
       implicit none

       integer::nr1
       real(pre)::dg_pi,dg_2pi,dg_halfpi,dg_2pipitch,dtheta,omega_pattern
       real(pre),dimension(nrmax)::surface_den,sound_speed,dg_r,rh
       real(pre),dimension(nthetamax)::theta
       real(pre),dimension(nlegenmax)::legendre_P,legendre_P_prime
!
!       real(pre),allocatable,dimension(:,:)::dg_rho,qlm_less,qlm_more
!       real(pre),allocatable,dimension(:,:)::phigrid
!       type(forcegrid),allocatable,dimension(:,:)::fgrid
!
       real(pre),dimension(nthetamax,nrmax)::dg_rho
       real(pre),dimension(0:nlegenmax-1,nrmax)::qlm_less,qlm_more
       real(pre),dimension(nthetamax,nrmax)::phigrid
       type(forcegrid),dimension(nthetamax,nrmax)::fgrid


       real(pre) :: timerstart=0,timerend=0,timertotal=0
       real(pre), dimension(2) :: dgelapsed
       real :: dgtotaltime
       contains

       subroutine initialize_gravity()
        
        nr1=nr+1
        dg_pi=acos(-one)
        dg_2pi=2*dg_pi
        dg_2pipitch=dg_2pi*pitch
        dg_halfpi=0.5*dg_pi
        angle=angle*dg_pi/180d0
        pitch=pitch*dg_pi/180d0
        dtheta=angle/dble(ntheta)
        omega_pattern=sqrt(mstar/rcorot**3)
        !allocate(surface_den(nr))
        !allocate(sound_speed(nr))
        !allocate(dg_r(nr1))
        !allocate(rh(nr))
        !allocate(theta(ntheta))
        !allocate(qlm_less(0:nlegen-1,nr1))
        !allocate(qlm_more(0:nlegen-1,nr1))
        !allocate(dg_rho(ntheta,nr))
        if(use_nlegen>60)then
          print *,"nlegen>60. The current program is not intended"
          print *,"to handle such large orders. Recall that you  "
          print *,"are requesting r^60. Killing the run..."
          stop
        endif
!        allocate(legendre_P(0:nlegen-1))
!        allocate(legendre_P_prime(0:nlegen-1))
!
!        if(use_grid)then
!          allocate(phigrid(ntheta,nr)) 
!          allocate(fgrid(ntheta,nr))
!        endif

       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine get_legendre(x)
         real(pre),intent(in)::x
         integer::ip
    
         legendre_P(0)=one
         legendre_P(1)=x
         legendre_P(2)=half*(three*x*x-one)

         do ip=3,nlegen-1
           legendre_P(ip)=two*x*legendre_P(ip-1)-legendre_P(ip-2) &
     &        -(x*legendre_P(ip-1)-legendre_P(ip-2))/dble(ip)
         enddo

       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine get_legendre_prime(x,y)
         real(pre),intent(in)::x,y
         integer::ip
    
         legendre_P_prime(0)=zero

         do ip=1,nlegen-1
          legendre_P_prime(ip)=dble(ip)*legendre_P(ip-1)+x &
              * legendre_P_prime(ip-1)
         enddo
         do ip=1,nlegen-1
           legendre_P_prime(ip)=-legendre_P_prime(ip)*y
         enddo

       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function powerSpectrum(m)
         real(pre)::m
         !powerSpectrum=0.5d0/(1d2+m**2)**1.3d0*power_fac
         powerSpectrum=0.5d0/(1d2)**1.3d0*power_fac/(m-one)**2
       end function
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function phase_wave(rad,r0,m)
         real(pre)::rad,r0,m
         phase_wave=-log(rad/r0)/pitch*m
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function phase_wave_prime(rad,m)
         real(pre)::rad,m
         phase_wave_prime=-one/(pitch*rad)*m
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function phase_angle(omega,m,time)
         real(pre)::omega,time,m
         phase_angle=cos(dg_2pi*omega*time*m)*dg_2pi/phase_drift&
     &              +dg_2pi*omega_pattern*time
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function wave_number_r(rad,m)
         real(pre)::rad,m
         wave_number_r=dg_2pi/(rad*(exp(dg_2pipitch/m))-one)
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function wave_number_r_prime(rad,m)
         real(pre)::rad,m
         wave_number_r_prime=-dg_2pi/(rad*rad*(exp(dg_2pipitch/m))-one)
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine make_disk_model()
         integer::ir,it
         real(pre)::h,omega,s,rho0,mass

         mass=zero
         do it=1,ntheta
           theta(it)=(dble(it)-half)*dtheta-angle*half
         enddo
         dg_r(nr1)=dble(nr1-1)*annulus
         do ir=1,nr
           dg_r(ir)=dble(ir-1)*annulus
           rh(ir)=dg_r(ir)+annulus*half
           if(rh(ir)>rmin.and.rh(ir)<dg_rmax)then
            surface_den(ir)=reference_sig &
     &                     *(reference_r/rh(ir))**power_sig
           else
            surface_den(ir)=zero
           endif
           omega=sqrt(mstar*dg_msun*Gcgs/(rh(ir)*dg_au)**3)
           sound_speed(ir)=toomre_q*surface_den(ir)*dg_pi*Gcgs &
     &                    /omega
           h=sound_speed(ir)/omega
           if(h>zero)then
            rho0=surface_den(ir)/(sqrt(dg_pi)*h)
            do it=1,ntheta
              s=rh(ir)*theta(it)
              dg_rho(it,ir)=max(rho0*exp(-(s*dg_au/h)**2),small_rho)
              mass=mass+dg_rho(it,ir)*rh(ir)**2*dtheta*dg_2pi*annulus*dg_au**3/dg_msun
            enddo
           else
            do it=1,ntheta
               dg_rho(it,ir)=zero
            enddo
           endif
         enddo
         ! create some basic sanity diagnostics

         open(unit=666,file="profiles_1d.dat")
         write(666,"(A)")"#R(AU),SIGMA(CGS),SOUND_SPEED(CGS)"
         write(666,"(A,1X,1pe15.8)")"#DISKMASS: ",mass
         do ir=1,nr
           write(666,"(4(1pe15.8,1X))")rh(ir),surface_den(ir), &
     &                             sound_speed(ir),dg_rho(ntheta/2+1,ir)
         enddo
         close(666)
         open(666,file="density.gnuplot")
         do ir=1,nr
           do it=1,ntheta
             write(666,"(3(1pe15.8,1X))a")rh(ir),theta(it),dg_rho(it,ir)
           enddo
         enddo
         close(666)


       end subroutine make_disk_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine get_spherical_moments()
        integer::ir,it,l
        real(pre)::x

        do ir=1,nr1
          do l=0,nlegen-1
           qlm_less(l,ir)=zero
           qlm_more(l,ir)=zero 
          enddo
        enddo
        do ir=2,nr1
         do l=0,nlegen-1
          qlm_less(l,ir)=qlm_less(l,ir-1)
          qlm_more(l,nr1-ir+1)=qlm_more(l,nr1-ir+2)
!          if (l.eq.0.and.ir.lt.12) then
!             open (unit=1,file='user.out',access='append')
!             write (1,*) 'In GetSphMom First loop, l=',l,' ir=',ir,'LESS=',qlm_less(l,ir)
!             close(unit=1)
!          end if
         enddo
         do it=1,ntheta
          x=cos(dg_halfpi-theta(it)) ! x is colatitude
          call get_legendre(x)
          do l=0,nlegen-1
            qlm_less(l,ir)=qlm_less(l,ir)+legendre_P(l)&
     &                  *rh(ir-1)**dble(l) &
     &                  *rh(ir-1)**2  &
     &                  *annulus*dg_au**3/dg_msun*dg_2pi*dtheta*dg_rho(it,ir-1)
            qlm_more(l,nr1-ir+1)=qlm_more(l,nr1-ir+1)+legendre_P(l)&
     &                  /rh(nr1-ir+1)**(dble(l+1))&
     &                  *rh(nr1-ir+1)**2*dg_au**3/dg_msun*dg_2pi*dtheta&
     &                  *dg_rho(it,nr1-ir+1)*annulus
!            if (l.eq.0.and.ir.lt.12) then
!               open (unit=1,file='user.out',access='append')
!               write (1,*) 'In GetSphMom, l=',l,' ir=',ir,' it=',it,'rh(ir-1)=',rh(ir-1)
!               write (1,*) 'rh(ir-1)**dble(l)=',rh(ir-1)**dble(l),' LP(l)',legendre_P(l)
!               write (1,*) 'dtheta & dg_rho(it,ir-1)=',dtheta,dg_rho(it,ir-1)
!               write (1,*) 'LESS=',qlm_less(l,ir),' MORE',qlm_more(l,nr1-ir+1),'nr1-ir+1=',nr1-ir+1
!               close(unit=1)
!            endif

          enddo
         enddo
        enddo 
        if (use_grid)call create_gravity_grid()
       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function potential2d(rad,t)
        real(pre),intent(in)::rad,t
        real(pre)::phi0,phi1,x,LP,R0,R1,RP0,RP1
        integer::ir,l

        ! determine ir
        x=cos(dg_halfpi-t)
        call get_legendre(x)
        phi0=zero;phi1=zero

        ir=int(rad/annulus)+1
        if(ir>=nr1)then
          ir=nr1
          do l=0,use_nlegen-1
            phi0=phi0-legendre_P(l)*qlm_less(l,ir)/rad**(l+1)
          enddo
          potential2d=phi0
        else
          do l=0,use_nlegen-1
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate some temporary variables to speed up calculations
             LP = legendre_P(l)
             R0 = dg_r(ir)**l ; R1 = R0*dg_r(ir) ; RP0 = dg_r(ir+1)**l ; RP1 = RP0*dg_r(ir+1) 
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate the potential
             if(ir>1)then
                phi0=phi0-LP*qlm_less(l,ir)/R1
             endif
             phi0=phi0-LP*qlm_more(l,ir)*R0
             phi1=phi1-LP*qlm_less(l,ir+1)/RP1
             phi1=phi1-LP*qlm_more(l,ir+1)*RP0
          enddo
          potential2d=phi0+(phi1-phi0)*(rad-dg_r(ir))/annulus
       endif
 
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function potential_pert(rad,t,p,time)
        real(pre),intent(in)::rad,t,p,time
        real(pre)::phi0,phi1,z,t0,omega,am,rcyl,m
        integer::im


        if (rad==zero)then
         potential_pert=zero
         return
        endif

        t0=zero 
        rcyl=rad*cos(t)
        z=rad*sin(t)
        if(use_grid) then
           phi0=smooth_potential_from_grid(rcyl,t0)
        else
           phi0=potential2d(rcyl,t0)
        endif

        phi1=zero
        omega=sqrt(mstar/rcyl**3)

        do im=mmin,mmax
          m=dble(im)
          am=phi0*powerSpectrum(m)
          phi1=phi1+am*(cos(m*(p-phase_wave(rcyl,rcorot,m))&
     &        +phase_angle(omega,m,time))) &
     &        *exp(-(wave_number_r(rcyl,m)*z)**2)
        enddo
        potential_pert=phi1
 
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(pre) function potential_all(rad,t,p,time)
        real(pre),intent(in)::rad,t,p,time
        real(pre)::phi0


        if (rad==zero)then
         potential_all=zero
         return
        endif

        if(use_grid) then
          phi0=smooth_potential_from_grid(rad,t)
        else
          phi0=potential2d(rad,t)
        endif
        potential_all=phi0+potential_pert(rad,t,p,time)
 
       end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine force_and_potn_smooth(rad,t,sint,cost,sinp,cosp,fx,fy,fz,potential2d)
        real(pre),intent(in)::rad,t,sint,cost,sinp,cosp
        real(pre),intent(out)::fx,fy,fz,potential2d
        real(pre)::fr0,fr1,ft0,ft1,x,y,frout,ftout
        real(pre)::QLM0,QLM1,LP,LPP,LP_QLM0,LP_QLM1,LPP_QLM1
        real(pre)::phi0,phi1,RM,R0,R1,R2,RPM,RP0,RP1,RP2
        integer::ir,l

        ! determine ir
        x=cos(dg_halfpi-t)
        y=sin(dg_halfpi-t)
        call get_legendre(x)
        call get_legendre_prime(x,y)
        fr0=zero;fr1=zero
        ft0=zero;ft1=zero
!
        phi0=zero;phi1=zero


        ir=int(rad/annulus)+1
        if(ir>=nr1)then
          ir=nr1
          do l=0,use_nlegen-1
             QLM0 = qlm_less(l,ir)
             LP = legendre_P(l)
             LP_QLM0 = LP*QLM0
             R1 = rad**(l+1)
!
             phi0=phi0-LP_QLM0/R1
!
             fr0=fr0-LP_QLM0/(R1*rad)*dble(l+1)
             ft0=ft0+legendre_P_prime(l)*QLM0/R1
          enddo
          potential2d=phi0
          frout=fr0;ftout=ft0/rad
        else
          do l=0,use_nlegen-1
             !!!!!!!!!! Calculate temporary variables used to cut-down on repeat calculations 
             QLM0 = qlm_less(l,ir) ; QLM1 = qlm_less(l,ir+1)
             LP = legendre_P(l) ; LPP = legendre_P_prime(l)
             LP_QLM0 = LP*QLM0 ; LP_QLM1 = LP*QLM1 ; LPP_QLM1 = LPP*QLM1
             RM  = dg_r(ir)**dble(l-1)   ; R0 = RM*dg_r(ir)    ; R1  = R0*dg_r(ir)   ; R2 =R1*dg_r(ir);
             RPM = dg_r(ir+1)**dble(l-1) ; RP0 = RPM*dg_r(ir+1); RP1 = RP0*dg_r(ir+1); RP2=RP1*dg_r(ir+1);
             !
             !!!!!!!!!! Now calculate the main quantities
             if(ir>1)then
                phi0=phi0-LP_QLM0/R1
                !
                fr0=fr0-LP_QLM0/R2*dble(l+1)
                fr0=fr0+LP_QLM0*RM*dble(l)
                ft0=ft0+LPP*QLM0/R1
             endif
             phi0=phi0-LP_QLM0*R0
             phi1=phi1-LP_QLM1/RP1
             phi1=phi1-LP_QLM1*RP0
             !
             ft0=ft0+LPP*QLM0*R0
             fr1=fr1-LP_QLM1/RP2*dble(l+1)
             fr1=fr1+LP_QLM1*RPM*dble(l)
             ft1=ft1+LPP_QLM1/RP1
             ft1=ft1+LPP_QLM1*RP0
          enddo
          potential2d=phi0+(phi1-phi0)*(rad-dg_r(ir))/annulus
          frout=fr0+(fr1-fr0)*(rad-dg_r(ir))/annulus
          ftout=zero
          if(rad>zero)ftout=(ft0+(ft1-ft0)*(rad-dg_r(ir))/annulus)/rad
       endif
       fx=frout*cosp*cost+ftout*sint*cosp
       fy=frout*sinp*cost+ftout*sint*sinp
       fz=frout*sint-ftout*cost
     end subroutine force_and_potn_smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine force_and_potn_pert(rcyl,z,p,sinp,cosp,time,fx,fy,fz,phi_pert)
        real(pre),intent(in)::rcyl,z,p,sinp,cosp,time
        real(pre),intent(out)::fx,fy,fz,phi_pert
        real(pre)::phi0,t0,omega,am,m,fr,fp
        real(pre)::wnr=0,wnr2=0,wnr2z2=0,pw=0,pa=0
        real(pre)::angle,am_ewnr2z2,s_am_e,c_am_e,phi1
        integer::im

        t0=zero 

        if(use_grid) then
           phi0=smooth_potential_from_grid(rcyl,t0)
        else
           phi0=potential2d(rcyl,t0)
        endif
        omega=omega_pattern              !sqrt(mstar/rcyl**3)
        phi1=zero
        fx=zero;fy=zero;fz=zero
        fr=zero;fp=zero

        do im=mmin,mmax

          m=dble(im)
          am=phi0*powerSpectrum(m)

          !!!!!!!!!! Calculate the temporary quantities used to speed-up the calculations
          wnr    = wave_number_r(rcyl,m)
          wnr2   = wnr*wnr
          wnr2z2 = wnr2*z*z
          angle  = m*(p - phase_wave(rcyl,rcorot,m)) + phase_angle(omega,m,time)
        
          am_ewnr2z2 = am*exp(-wnr2z2)
          s_am_e = sin(angle) * am_ewnr2z2
          c_am_e = cos(angle) * am_ewnr2z2

          !!!!!!!!!! Now get the potential and the forces.
          phi1=phi1+c_am_e
          fr=fr - s_am_e*m*phase_wave_prime(rcyl,m) &
     &      +c_am_e*z**2*wave_number_r_prime(rcyl,m)*wnr*two
          fp=fp + s_am_e*m/rcyl
          fz=fz + c_am_e*two*z*wnr2z2

          !!!!!!!!! Old version of the force code.
!          fr=fr-am*sin(m*(p-phase_wave(rcyl,rcorot,m)) &
!     &      +phase_angle(omega,m,time))*m&
!     &      *phase_wave_prime(rcyl,m)&
!     &      *exp(-(wave_number_r(rcyl,m)*z)**2) &
!     &      +am*cos(m*(p-phase_wave(rcyl,rcorot,m)) &
!     &      +phase_angle(omega,m,time))*z**2&
!     &      *wave_number_r_prime(rcyl,m)*wave_number_r(rcyl,m)*two &
!     &      *exp(-(wave_number_r(rcyl,m)*z)**2)
!          fp=fp+am*sin(m*(p-phase_wave(rcyl,rcorot,m))&
!     &        +phase_angle(omega,m,time))*exp(-(wave_number_r(rcyl,m)&
!     &        *z)**2)*m/rcyl
!          fz=fz+am*cos(m*(p-phase_wave(rcyl,rcorot,m)) &
!     &      +phase_angle(omega,m,time))*two*z &
!     &      *wave_number_r(rcyl,m)**2*exp(-(wave_number_r(rcyl,m)*z)**2)
         ! if (rcyl.gt.and.rcyl.lt.and.) then 
!          open (unit=1,file='user.out',access='append')
!          write (1,*)  'rcyl,z,p',rcyl,z,p
!          write (1,*) 'phi0,powerSpectrum(m)',phi0,powerSpectrum(m)
!          write (1,*)  'm,am,fr,fp,fz',m,am,fr,fp,fz
!          close(unit=1)
         ! end if
        enddo

        phi_pert=phi1

        fx=fr*cosp-fp*sinp
        fy=fr*sinp+fp*cosp
        
!        open (unit=1,file='user.out',access='append')
!        write (1,*)  'f(x,y,z,r,p)',fx,fy,fz,fr,fp
!        close(unit=1)
        
       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine force_smooth_2(rad,t,frout,ftout)
        real(pre),intent(in)::rad,t
        real(pre),intent(out)::frout,ftout
        real(pre)::fr0,fr1,ft0,ft1,x,y,p
        integer::ir,l

        ! determine ir
        x=cos(dg_halfpi-t)
        y=sin(dg_halfpi-t)
        p=zero
        call get_legendre(x)
        call get_legendre_prime(x,y)
        fr0=zero;fr1=zero
        ft0=zero;ft1=zero

        ir=int(rad/annulus)+1

        !open (unit=1,file='user.out',access='append')
        !write (1,*)  ' rad=',rad,' t=',t,'ir=',ir,' nr1=',nr1
        !close(unit=1)

        if(ir>=nr1)then
          ir=nr1
          do l=0,use_nlegen-1
           fr0=fr0-legendre_P(l)*qlm_less(l,ir)/rad**dble(l+2)*dble(l+1)
           ft0=ft0+legendre_P_prime(l)*qlm_less(l,ir)/rad**dble(l+1)
          enddo
          frout=fr0;ftout=ft0/rad
        else
          do l=0,use_nlegen-1
            if(ir>1)then
             fr0=fr0-legendre_P(l)*qlm_less(l,ir)/dg_r(ir)**dble(l+2) &
     &          *dble(l+1)
!             fr0=fr0+legendre_P(l)*qlm_more(l,ir)*dg_r(ir)**dble(l-1) &
!     &         *dble(l)
             ft0=ft0+legendre_P_prime(l)*qlm_less(l,ir)/dg_r(ir)**dble(l+1)
            endif
!        if (rad.gt.10.49.and.rad.lt.10.51.and.t.gt.-0.0001.and.t.lt.0.0001) then
!           open (unit=1,file='user.out',access='append')
!           write (1,*) 'l=',l,' ir=',ir,' fr0=',fr0,'  ft0=',ft0
!           write (1,*) 'LP(l)=',legendre_P(l),' QLMLess=',qlm_less(l,ir),' dg_r(ir)=',dg_r(ir)
!           close(unit=1)
!        endif
        fr0=fr0+legendre_P(l)*qlm_more(l,ir)*dg_r(ir)**dble(l-1) &
     &         *dble(l)
            ft0=ft0+legendre_P_prime(l)*qlm_more(l,ir)*dg_r(ir)**dble(l)
            fr1=fr1-legendre_P(l)*qlm_less(l,ir+1)/dg_r(ir+1)**dble(l+2)&
     &         *dble(l+1)
            fr1=fr1+legendre_P(l)*qlm_more(l,ir+1)*dg_r(ir+1)**dble(l-1)&
     &         *dble(l)
            ft1=ft1+legendre_P_prime(l)*qlm_less(l,ir+1)&
     &         /dg_r(ir+1)**dble(l+1)
            ft1=ft1+legendre_P_prime(l)*qlm_more(l,ir+1) &
     &         *dg_r(ir+1)**dble(l)
          enddo
          frout=fr0+(fr1-fr0)*(rad-dg_r(ir))/annulus
          ftout=zero
          if(rad>zero)ftout=(ft0+(ft1-ft0)*(rad-dg_r(ir))/annulus)/rad

          
        if (rad.gt.10.49.and.rad.lt.10.51.and.t.gt.-0.0001.and.t.lt.0.0001) then
           open (unit=1,file='user.out',access='append')
           write (1,*) 'ir=',ir,' nr1=',nr1,' frout=',frout,' ftout=',ftout
           close(unit=1)
        endif


        endif
 
       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine create_gravity_grid()
        integer::it,ir
        real(pre)::fr,ft

        print *,"# creating gravity grid"
        do ir=1,nr
         do it=1,ntheta
           phigrid(it,ir)=potential2d(rh(ir),theta(it))
           call force_smooth_2(rh(ir),theta(it),fr,ft)
           fgrid(it,ir)%fr=fr;
           fgrid(it,ir)%ft=ft
         enddo
        enddo
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(pre) function catmull_rom_spline(y0,y1,y2,y3,y4,dx,k)
        real(pre),intent(in)::y0,y1,y2,y3,y4,dx,k
        real(pre)::h00,h10,h01,h11,q,m1,m2,m3,h

        h=k
        if(k>dx)h=h-dx
        h=h/dx
        m2=(y3-y1)/(two*dx)
        h00=(one+two*h)*(one-h)**2
        h10=h*(one-h)**2
        h01=h*h*(three-two*h)
        h11=h*h*(h-one)

        if(h<dx)then
          m1=(y2-y0)/(two*dx)
          q=y1*h00+h10*m1*dx+h01*y2+h11*m2*dx
        else
          m3=(y4-y2)/(two*dx)
          q=y2*h00+h10*m2*dx+h01*y3+h11*m3*dx
        endif
         
        catmull_rom_spline=q

      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      real(pre) function natural_cubic_spline_fixeddx(y0,y1,y2,dx,h)
        real(pre),intent(in)::y0,y1,y2,dx,h
        real(pre)::a1,b1,a2,b2,t,q,k1,k2,k3,delb
        real(pre),dimension(3)::b
        real(pre),dimension(3,3)::a

        a(1,1)=two/(dx) ! use actual difference in case nonequal spacing
        a(1,2)=one/(dx)
        a(1,3)=zero
        a(2,1)=one/(dx)
        a(2,2)=four/dx
        a(2,3)=one/(dx)
        a(3,1)=zero
        a(3,2)=one/(dx)
        a(3,3)=two/(dx)

        b(1)=three*(y1-y0)/(dx)**2
        b(2)=three*((y1-y0)/(dx)**2+(y2-y1)/(dx)**2)
        b(3)=three*(y2-y0)/(dx)**2

        delb=b(1)-b(3)
        k2=(b(3)-(b(2)-delb*half))/(a(3,2)-two*a(3,3))
        k3=(b(3)-a(3,2)*k2)/a(3,3)
        k1=(b(1)-a(1,2)*k2)/a(1,1)

        a1=k1*dx-(y1-y0)
        b1=-k2*dx+y1-y0
        a2=k2*dx-(y2-y1)
        b2=-k3*dx+y2-y1 

        if(h<dx)then
          t=h/dx
          q=(one-t)*y0+t*y1+t*(one-t)*(a1*(one-t)+b1*t)
        else
          t=(h-dx)
          q=(one-t)*y1+t*y2+t*(one-t)*(a2*(one-t)+b2*t)
        endif
        natural_cubic_spline_fixeddx=q


      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(pre) function get_quad(y0,y1,y2,dx,h)
        real(pre),intent(in)::y0,y1,y2,dx,h
        real(pre)::a=0,b=0,c=0

        a=y0
        b=(three*(y2-a)-four*(y2-y1))/dx
        c=((y2-a)/dx-b)/dx
        get_quad=a+h*(b+c*h)

      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(pre) function smooth_potential_from_grid(rad,t)
        integer::it,ir
        real(pre),intent(in)::rad,t
        real(pre)::h
        real(pre)::p0,p1,p2,y0,y1,y2

        it=int((t+angle*half)/(dtheta))+1                   ! Original code
        !it=int((t - dg_halfpi + angle*half)/(dtheta))+1      ! I think this is what's required to make z=0 the half-way point
        ir=int(rad/annulus)+1
        if(it<2.or.it>ntheta-1.or.ir<2.or.ir>nr-1)then
           smooth_potential_from_grid=potential2d(rad,t)
           return
        endif

        ! use quadratics for interpolation
        ! interpolate first in the theta direction
        h=t-theta(it-1)
        y0=phigrid(it-1,ir-1)
        y1=phigrid(it  ,ir-1)
        y2=phigrid(it+1,ir-1)
        p0=get_quad(y0,y1,y2,dtheta,h)

        y0=phigrid(it-1,ir)
        y1=phigrid(it  ,ir)
        y2=phigrid(it+1,ir)
        p1=get_quad(y0,y1,y2,dtheta,h)

        y0=phigrid(it-1,ir+1)
        y1=phigrid(it  ,ir+1)
        y2=phigrid(it+1,ir+1)
        p2=get_quad(y0,y1,y2,dtheta,h)

        ! now once in the r direction
        h=rad-rh(ir-1)
        y0=p0
        y1=p1
        y2=p2
        smooth_potential_from_grid=get_quad(y0,y1,y2,annulus,h)


      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 

      subroutine smooth_force_from_spline(rad,t,sint,cost,sinp,cosp,fx,fy,fz,phi_from_grid) !catmull-rom
        integer::it,ir
        real(pre),intent(in)::rad,t,sint,cost,sinp,cosp
        real(pre),intent(out)::fx,fy,fz,phi_from_grid
        real(pre)::fr,ft,h
        real(pre)::p0,p1,p2,p3,p4,y0,y1,y2,y3,y4

        it=int((t+angle*half)/(dtheta))+1                   ! Original code
        !it=int((t - dg_halfpi + angle*half)/(dtheta))+1      ! I think this is what's required to make z=0 the half-way point
        ir=int(rad/annulus)+1
        if(it<3.or.it>ntheta-2.or.ir<3.or.ir>nr-2)then
           call force_and_potn_smooth(rad,t,sint,cost,sinp,cosp,fx,fy,fz,phi_from_grid)
           return
        endif

        ! use quadratics for interpolation
        ! interpolate first in the theta direction
        h=t-theta(it-1)
        y0=fgrid(it-2,ir-2)%ft
        y1=fgrid(it-1,ir-2)%ft
        y2=fgrid(it  ,ir-2)%ft
        y2=fgrid(it+1,ir-2)%ft
        y2=fgrid(it+2,ir-2)%ft
        p0=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir-1)%ft
        y1=fgrid(it-1,ir-1)%ft
        y2=fgrid(it  ,ir-1)%ft
        y2=fgrid(it+1,ir-1)%ft
        y2=fgrid(it+2,ir-1)%ft
        p1=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir)%ft
        y1=fgrid(it-1,ir)%ft
        y2=fgrid(it  ,ir)%ft
        y2=fgrid(it+1,ir)%ft
        y2=fgrid(it+2,ir)%ft
        p2=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir+1)%ft
        y1=fgrid(it-1,ir+1)%ft
        y2=fgrid(it  ,ir+1)%ft
        y2=fgrid(it+1,ir+1)%ft
        y2=fgrid(it+2,ir+1)%ft
        p3=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir+2)%ft
        y1=fgrid(it-1,ir+2)%ft
        y2=fgrid(it  ,ir+2)%ft
        y2=fgrid(it+1,ir+2)%ft
        y2=fgrid(it+2,ir+2)%ft
        p4=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        ! now once in the r direction
        h=rad-rh(ir-1)
        y0=p0
        y1=p1
        y2=p2
        y3=p3
        y4=p4
        ft=catmull_rom_spline(y0,y1,y2,y3,y4,annulus,h)

        ! now for fr

        h=t-theta(it-1)
        y0=fgrid(it-2,ir-2)%fr
        y1=fgrid(it-1,ir-2)%fr
        y2=fgrid(it  ,ir-2)%fr
        y2=fgrid(it+1,ir-2)%fr
        y2=fgrid(it+2,ir-2)%fr
        p0=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir-1)%fr
        y1=fgrid(it-1,ir-1)%fr
        y2=fgrid(it  ,ir-1)%fr
        y2=fgrid(it+1,ir-1)%fr
        y2=fgrid(it+2,ir-1)%fr
        p1=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir)%fr
        y1=fgrid(it-1,ir)%fr
        y2=fgrid(it  ,ir)%fr
        y2=fgrid(it+1,ir)%fr
        y2=fgrid(it+2,ir)%fr
        p2=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir+1)%fr
        y1=fgrid(it-1,ir+1)%fr
        y2=fgrid(it  ,ir+1)%fr
        y2=fgrid(it+1,ir+1)%fr
        y2=fgrid(it+2,ir+1)%fr
        p3=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        y0=fgrid(it-2,ir+2)%fr
        y1=fgrid(it-1,ir+2)%fr
        y2=fgrid(it  ,ir+2)%fr
        y2=fgrid(it+1,ir+2)%fr
        y2=fgrid(it+2,ir+2)%fr
        p4=catmull_rom_spline(y0,y1,y2,y3,y4,dtheta,h)

        ! now once in the r direction
        h=rad-rh(ir-1)
        y0=p0
        y1=p1
        y2=p2
        y3=p3
        y4=p4
        fr=catmull_rom_spline(y0,y1,y2,y3,y4,annulus,h)

        ! that was ft, now get fr

        fx=fr*cosp*cost+ft*sint*cosp
        fy=fr*sinp*cost+ft*sint*sinp
        fz=fr*sint-ft*cost

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_force_and_potn_from_grid(rad,t,sint,cost,sinp,cosp,fx,fy,fz,phi_from_grid)
        integer::it,ir
        real(pre),intent(in)::rad,t,sint,cost,sinp,cosp
        real(pre),intent(out)::fx,fy,fz,phi_from_grid
        real(pre)::fr,ft,h
        real(pre)::p0,p1,p2,y0,y1,y2,b,c
        real(pre)::phip0,phip1,phip2,phiy0,phiy1,phiy2,phib


        !!!!!!!!!!!!!!! Get the grid cell/element that the particle resides within.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        it=int((t+angle*half)/(dtheta))+1                   ! Original code
        !it=int((t - dg_halfpi + angle*half)/(dtheta))+1      ! I think this is what's required to make z=0 the half-way point
        ir=int(rad/annulus)+1



        !!!!!!!!!!!!!!!  If the position of the particle is outside the defined range of the grid, then do the explicit calculation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(it<2.or.it>ntheta-1.or.ir<2.or.ir>nr-1)then
           call force_and_potn_smooth(rad,t,sint,cost,sinp,cosp,fx,fy,fz,phi_from_grid)
           return
        endif



        !!!!!!!!!!!!!!! If the particle is within the confines of the grid, then interpolate to fine the force & potn at that point
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! use quadratics for interpolation
        ! interpolate first in the theta direction
        h=t-theta(it-1)
        y0=fgrid(it-1,ir-1)%ft   ;      phiy0=phigrid(it-1,ir-1)
        y1=fgrid(it  ,ir-1)%ft   ;      phiy1=phigrid(it  ,ir-1)
        y2=fgrid(it+1,ir-1)%ft   ;      phiy2=phigrid(it+1,ir-1)
!        p0=get_quad(y0,y1,y2,dtheta,h)   ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        b=(three*(y2-y0)-four*(y2-y1))/dtheta      ;    phib  = (three*(phiy2-phiy0)-four*(phiy2-phiy1))/dtheta
        p0=((y2-y0)/dtheta-b)/dtheta               ;    phip0 = ((phiy2-phiy0)/dtheta-phib)/dtheta
        p0 = y0+h*(b+p0*h)                         ;    phip0 = phiy0+h*(phib+phip0*h)


        y0=fgrid(it-1,ir)%ft     ;      phiy0=phigrid(it-1,ir)
        y1=fgrid(it  ,ir)%ft     ;      phiy1=phigrid(it  ,ir)
        y2=fgrid(it+1,ir)%ft     ;      phiy2=phigrid(it+1,ir)
        !p1=get_quad(y0,y1,y2,dtheta,h)   ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        b=(three*(y2-y0)-four*(y2-y1))/dtheta      ;    phib  = (three*(phiy2-phiy0)-four*(phiy2-phiy1))/dtheta
        p1=((y2-y0)/dtheta-b)/dtheta               ;    phip1 = ((phiy2-phiy0)/dtheta-phib)/dtheta
        p1 = y0+h*(b+p1*h)                         ;    phip1 = phiy0+h*(phib+phip1*h)

 
        y0=fgrid(it-1,ir+1)%ft     ;      phiy0=phigrid(it-1,ir+1)
        y1=fgrid(it  ,ir+1)%ft     ;      phiy1=phigrid(it  ,ir+1)
        y2=fgrid(it+1,ir+1)%ft     ;      phiy2=phigrid(it+1,ir+1)
        !p2=get_quad(y0,y1,y2,dtheta,h)   ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        b=(three*(y2-y0)-four*(y2-y1))/dtheta      ;    phib  = (three*(phiy2-phiy0)-four*(phiy2-phiy1))/dtheta
        p2=((y2-y0)/dtheta-b)/dtheta               ;    phip2 = ((phiy2-phiy0)/dtheta-phib)/dtheta
        p2 = y0+h*(b+p2*h)                         ;    phip2 = phiy0+h*(phib+phip2*h)



        ! now once in the r direction
        h=rad-rh(ir-1)
        y0=p0                      ;      phiy0 = phip0
        y1=p1                      ;      phiy1 = phip1
        y2=p2                      ;      phiy2 = phip2
        !ft=get_quad(y0,y1,y2,annulus,h)   ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        b=(three*(y2-y0)-four*(y2-y1))/annulus     ;    phib  = (three*(phiy2-phiy0)-four*(phiy2-phiy1))/annulus
        ft=h*(b + h*((y2-y0)/annulus-b)/annulus)   ;    phi_from_grid = h*(phib + h*((phiy2-phiy0)/annulus-phib)/annulus) 
        ft = y0  + ft                              ;    phi_from_grid = phiy0 + phi_from_grid



        ! that was ft, now get fr
        h=t-theta(it-1)

        y0=fgrid(it-1,ir-1)%fr
        y1=fgrid(it  ,ir-1)%fr
        y2=fgrid(it+1,ir-1)%fr
!        p0=get_quad(y0,y1,y2,dtheta,h)
        b=(three*(y2-y0)-four*(y2-y1))/dtheta   ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        c=((y2-y0)/dtheta-b)/dtheta
        p0 = y0+h*(b+c*h)

      
        y0=fgrid(it-1,ir)%fr
        y1=fgrid(it  ,ir)%fr
        y2=fgrid(it+1,ir)%fr
!        p1=get_quad(y0,y1,y2,dtheta,h)
        b=(three*(y2-y0)-four*(y2-y1))/dtheta   ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        p1=((y2-y0)/dtheta-b)/dtheta
        p1 = y0+h*(b+p1*h)
 
        y0=fgrid(it-1,ir+1)%fr
        y1=fgrid(it  ,ir+1)%fr
        y2=fgrid(it+1,ir+1)%fr
!        p2=get_quad(y0,y1,y2,dtheta,h) ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        b=(three*(y2-y0)-four*(y2-y1))/dtheta
        p2=((y2-y0)/dtheta-b)/dtheta
        p2 = y0+h*(b+p2*h)


        ! now once in the r direction
        h=rad-rh(ir-1)
        y0=p0
        y1=p1
        y2=p2
        !fr=get_quad(y0,y1,y2,annulus,h) ! ... Turning this function call off & doing the explicit calculation is marginally quicker
        b=(three*(y2-y0)-four*(y2-y1))/annulus
        fr=((y2-y0)/annulus-b)/annulus
        fr = y0+h*(b+fr*h)

        fx=fr*cosp*cost+ft*sint*cosp  
        fy=fr*sinp*cost+ft*sint*sinp  
        fz=fr*sint-ft*cost            
 
 
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_force_from_grid_lin(rad,t,sint,cost,sinp,cosp,fx,fy,fz,phi_from_grid)
        integer::it,ir
        real(pre),intent(in)::rad,t,sint,cost,sinp,cosp
        real(pre),intent(out)::fx,fy,fz,phi_from_grid
        real(pre)::fr,ft
        real(pre)::y0,y1,delt,delr

        it=int((t+angle*half)/(dtheta))+1                   ! Original code
        !it=int((t - dg_halfpi + angle*half)/(dtheta))+1      ! I think this is what's required to make z=0 the half-way point
        ir=int(rad/annulus)+1
        if(it<2.or.it>ntheta-1.or.ir<2.or.ir>nr-1)then
           call force_and_potn_smooth(rad,t,sint,cost,sinp,cosp,fx,fy,fz,phi_from_grid)
           return
        endif

        ! first find quadrant
        delr=rad-rh(ir)
        delt=t-theta(it)
         
        if(delr>=zero.and.delt>=zero)then ! first quadrant, use ir,ir+1, it, it+1
          y0=fgrid(it,ir)%ft+delt*(fgrid(it+1,ir)%ft)/dtheta     
          y1=fgrid(it,ir+1)%ft+delt*(fgrid(it+1,ir+1)%ft)/dtheta
          ft=y0+delr*(y1-y0)/annulus

          y0=fgrid(it,ir)%fr+delt*(fgrid(it+1,ir)%fr)/dtheta
          y1=fgrid(it,ir+1)%fr+delt*(fgrid(it+1,ir+1)%fr)/dtheta
          fr=y0+delr*(y1-y0)/annulus
        elseif(delr<=zero.and.delt>=zero)then ! second, use ir-1,ir, it, it+1
          y0=fgrid(it,ir-1)%ft+delt*(fgrid(it+1,ir-1)%ft)/dtheta
          y1=fgrid(it,ir)%ft+delt*(fgrid(it+1,ir)%ft)/dtheta
          ft=y0+delr*(y0-y1)/annulus

          y0=fgrid(it,ir-1)%fr+delt*(fgrid(it+1,ir-1)%fr)/dtheta
          y1=fgrid(it,ir)%fr+delt*(fgrid(it+1,ir)%fr)/dtheta
          fr=y0+delr*(y0-y1)/annulus
        elseif(delr<=zero.and.delt<=zero)then ! third, use ir-1,ir, it-1, it
          y0=fgrid(it-1,ir-1)%ft+delt*(fgrid(it,ir-1)%ft)/dtheta
          y1=fgrid(it-1,ir)%ft+delt*(fgrid(it,ir)%ft)/dtheta
          ft=y0+delr*(y0-y1)/annulus

          y0=fgrid(it-1,ir-1)%fr+delt*(fgrid(it,ir-1)%fr)/dtheta
          y1=fgrid(it-1,ir)%fr+delt*(fgrid(it,ir)%fr)/dtheta
          fr=y0+delr*(y0-y1)/annulus
        else ! use ir,ir+1,it-1,it
          y0=fgrid(it-1,ir)%ft+delt*(fgrid(it,ir)%ft)/dtheta
          y1=fgrid(it-1,ir+1)%ft+delt*(fgrid(it,ir+1)%ft)/dtheta
          ft=y0+delr*(y1-y0)/annulus

          y0=fgrid(it-1,ir)%fr+delt*(fgrid(it,ir)%fr)/dtheta
          y1=fgrid(it-1,ir+1)%fr+delt*(fgrid(it,ir+1)%fr)/dtheta
          fr=y0+delr*(y1-y0)/annulus
        endif
 

        fx=fr*cosp*cost+ft*sint*cosp
        fy=fr*sinp*cost+ft*sint*sinp
        fz=fr*sint-ft*cost

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine force_all(rsph,rcyl,z,t,sint,cost,p,sinp,cosp,time,fx,fy,fz,phi)
        real(pre),intent(in)::rsph,rcyl,z,t,sint,cost,p,sinp,cosp,time
        real(pre),intent(out)::fx,fy,fz,phi
        real(pre)::fxp=0,fyp=0,fzp=0,phipert

        phi=zero
        fx=zero;fy=zero;fz=zero
        fxp=zero;fyp=zero;fzp=zero

        if(use_grid)then
!           call smooth_force_and_potn_from_grid(rsph,t,sint,cost,sinp,cosp,fx,fy,fz,phi)
        else
!           call force_and_potn_smooth(rsph,t,sint,cost,sinp,cosp,fx,fy,fz,phi)
        endif        
!        call force_and_potn_pert(rcyl,z,p,sinp,cosp,time,fxp,fyp,fzp,phipert)    


!        open (unit=1,file='user.out',access='append')
!        write(1,*) ' fx,fy,fz & phipert & fxp,fyp,fzp= ',fx,' ',fy,' ',fz,' & ',phi,' & ',fxp,' ',fyp,' ',fzp
!        close(unit=1)  

        phi = phi + phipert
        fx=fx+fxp
        fy=fy+fyp
        fz=fz+fzp

 
       end subroutine force_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




       subroutine cleanup()

!         deallocate(surface_den,sound_speed,dg_r,rh,qlm_less,qlm_more, &
!     &     theta,dg_rho,legendre_P,legendre_P_prime)
!         if(use_grid)deallocate(fgrid,phigrid)

       end subroutine cleanup

      end module
