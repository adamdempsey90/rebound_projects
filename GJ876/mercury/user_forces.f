c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      USER_FORCES.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matthew John Payne
c
c Optional userdefined forces / routines that can be grafted onto Mercury
c -->> Non-conservative forces, etc, etc
c
c     (1) Drift Check (Force Early Finish)
c     (2) Forcible Migration 
c     (3) Gas Drag
c     (4) Forcible Mass Growth
c     (5) Integrator switches
c     (6) Relativistic Forces
c     (7) Tidal Drag
c     (8) Disk Potential
c     (9) Nagasawa Tidal Drag Model(s)
c
c------------------------------------------------------------------------------
c
c Start the subroutine
      subroutine user_forces (time,nbod,nbig,m,x,v,a,params,
     %     lmem,mem,h0,algor,initial,s,rphys)
c Parameter input.
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c Input/Output
      integer nbod, nbig, algor, initial(4)
      real*8 time,m(nbod),x(3,nbod),v(3,nbod),a(3,nbod),s(3,nbod)
      real*8 rphys(nbod)
      real*8 params(NPARAMS,nbod), h0
      integer lmem(NMESS)
      character*80 mem(NMESS)
      reaL*8 gm,semi,r,v2
c Local
      integer i,j,k
c
c
c------------------------------------------------------------------------------
c
c
c Allow the code to finish early if drifting occurs
      if (DRIFTSWITCH.eq.1) then
         call drift_check (time,nbig,m,x,v,params,lmem,mem,initial)
      end if 
c 
c Allow individual bodies to be forcibly migrated
      if (MIGRATIONSWITCH.eq.1) then 
         call migrate (time,nbig,m,x,v,a,params,h0)
      end if
c
c Allow the effects of gas drag to be modelled
      if (GASDRAGSWITCH.gt.0) then 
         if (GASDRAGSWITCH.eq.1) then 
            call simple_gas_drag (time,nbod,m,x,v,a,params,algor)
c            call old_gas_drag (time,nbod,m,x,v,a,params)
         else 
            if (GASDRAGSWITCH.eq.2) then 
               call eccentric_gas_drag (time,nbod,m,x,v,a,params)
               if (x(3,2).gt.1e-5) then ! If the binary companion is not in the z=0 plane, then things will fuck up
                  open (unit=1,file='user.out',access='append')
                  write(1,*) 'STOP: binary not in z=0 plane: x(3,2)=',
     %                 x(3,2)
                  close(unit=1)
                  stop
               end if
            end if
         end if
      end if
c
c Allow individual bodies to be forcibly increased in mass
      if (MASSGROWTHSWITCH.eq.1) then 
         call mass_growth (m,nbig,params,time,h0)
      end if
c
c Allow the integrator to be switched forcibly between HYBRID <--> BS modes
      if (BSTOHYBRIDSWITCH.eq.1) then 
         if (algor.eq.2) then
            call bs_to_hybrid(time,algor)
         end if
      end if
      if (HYBRIDTOBSSWITCH.eq.1) then 
         if (algor.eq.10) then
            call hybrid_to_bs(time,nbig,m,x,v,params,algor)
         end if
      end if
c
c Allow a simple relativistic forcing term to be included
      if (RELATIVITYSWITCH.eq.1) then 
         call relativity (time,nbod,m,x,v,a,params)
      end if
c
c
c
c Allow tidal forces between large bodies to be included
      if (TIDALSWITCH.eq.1) then
         call tidal (time,nbig,m,x,v,a,params,s,rphys,h0)
      end if
c For the Nagasawa tidal damping, all that happens in the user_forces
c routine is that the algorithm is forced to switch hybrid -> bs
c The actual tidal "kick" is implemented in the bs routine in
c the main body of the code.
      if(T4NAGASAWATIDALSWITCH.eq.1.or.T4NAGASAWABINARYSWITCH.eq.1) then
c         call t4nagasawatidal (time,nbig,m,x,v,params,h0)
         if (algor.eq.10) then
            call hybrid_to_bs(time,nbig,m,x,v,params,algor)
         end if
      end if
c     
c
c
c Allow the effect from a disk potential to be included
      if (DISKPOTENTIALSWITCH.eq.1) then
         call diskpotn (time,nbig,m,x,v,a,params)
      end if
c
c
c Allow the gravity from an extended disk to be included
c (This uses the routine(s) developed by Aaron Boley)
      if (DISKGRAVITYSWITCH.eq.1) then
         call diskgravity (time,nbod,m,x,v,a,params,h0)
      end if
c
c
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c
c
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      DRIFT_CHECK
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c If the DRIFTSWITCH switch = 1, then this routine will be called
c This routine 
c
c------------------------------------------------------------------------------
c
      subroutine drift_check (time,nbig,m,x,v,params,lmem,mem,initial)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig, initial(4)
      real*8 time,m(nbig),x(3,nbig),v(3,nbig),a(3,nbig)
      real*8 params(NPARAMS,nbig)
      integer lmem(NMESS)
      character*80 mem(NMESS)
c
c Local
      integer j,k
      real*8 gm, semi(nbig), r, v2
      integer offswitch
c
c------------------------------------------------------------------------------
c
      offswitch =0
      if (nbig.ne.initial(1)) then 
         offswitch = 1         
      end if
c
c     For each of the planets/massive bodies...
      do j = 2, nbig
c     
c     Get the r,v & a==semi (note use of Chambers' routine mco_x2a).
            gm = m(1) + m(j)
            call mco_x2a (gm,x(1,j),x(2,j),x(3,j),
     %           v(1,j),v(2,j),v(3,j),semi(j),r,v2)
c     
c     If the semi-major axis of the large body has changed significantly,
c     or a collision has occured, then do data dump & end the simulation            
            if (semi(j) .lt. (params(1,j)*(1-DRIFTTOL))) then 
               offswitch = 1
            end if
c
            if (offswitch.eq.1) then
c     ----------Do a final data dump (just to user.out)
            open (unit=1, file='user.out', access='append')
            do k = 2, nbig
               write(1,*) 't(yr)= ',time/365.25,' k=',k,
     %         ' m(/Sun)=',m(k)/K2,
     %         ' semi= ',semi(k), ' initial semi= ',params(1,k),
     %         ' x,y,z= ',x(1,k),x(2,k),x(3,k),
     %         ' vx,vy,vz= ',v(1,k),v(2,k),v(3,k)
            end do
            close(unit=1)
c     ----------Call mio_err (see E.g. L6319 & L5666)
            call mio_err (6,mem(81),lmem(81),
     %        mem(108),lmem(108),' ',1,mem(109),lmem(109))
            end if
            if (semi(j) .gt. (params(1,j)*(1+DRIFTTOL))) then
c     ----------Do a final data dump (just to user.out)
            open (unit=1, file='user.out', access='append')
            do k = 2, nbig
               write(1,*) 't(yr)= ',time/365.25,' k=',k,
     %         ' m(/Sun)=',m(k)/K2,
     %         ' semi= ',semi(k), ' initial semi= ',params(1,k),
     %         ' x,y,z= ',x(1,k),x(2,k),x(3,k),
     %         ' vx,vy,vz= ',v(1,k),v(2,k),v(3,k)
            end do
            close(unit=1)
c     ----------Call mio_err (see E.g. L6319 & L5666)
            call mio_err (6,mem(81),lmem(81),
     %        mem(108),lmem(108),' ',1,mem(109),lmem(109))
            end if
c
         end do
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIGRATE.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne
c 
c Date: October 2009  &  Oct 2010
c
c Want to be able to force a migration.
c
c Ultimately want this module to be capable of forcing migration
c using a variety of different models
c
c Of prime importance at the moment is reproducing the Lee & Peale '02 model
c The appendix gives Lea & Peale's expression for incrementing x_i & v_i in terms
c    of the elements a, e, i, ...
c However, in order to get these one needs to solve Kepler's eqn anyway, so
c    it would seem simpler just to directly increment a & i
c Will try both
c
c NB
c     gm = grav const * (central + secondary mass)
c     q = perihelion distance
c     e = eccentricity
c     i = inclination                 )
c     p = longitude of perihelion !!! )   in            -- \varpi = \omega + \Omega
c     n = longitude of ascending node ) radians         -- \Omega
c     l = mean anomaly                )
c     g = argument of perihelion                        -- \omega
c     x,y,z = Cartesian positions  ( units the same as a )
c     u,v,w =     "     velocities ( units the same as sqrt(gm/a) )
c
c
c------------------------------------------------------------------------------
c
      subroutine migrate (time,nbig,m,x,v,a,params,h0)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig
      real*8 time,m(nbig),x(3,nbig),v(3,nbig),a(3,nbig)
      real*8 params(NPARAMS,nbig),h0
c
c Local
      integer j,k
      integer ngflag_dummy,ngf_dummy,opt_dummy
      real*8 semi,r,q,e,i,p,n,l,gm,r_disk,oldsemi,temp,factor,v2
      real*8 xj(3,nbig),vj(3,nbig)
      real*8 jcen_dummy(3)
      real*8 qinner,einner,iinner,pinner,ninner,linner,period_o,period_i
      real*8 semiinner
c
c------------------------------------------------------------------------------
c     
c     If is within the correct time range
      if (MIGRATIONSTART.lt.(time/365.25)) then 
         if (MIGRATIONEND.gt.(time/365.25)) then 
c
c
c     If the migration mechanism = LP03, then 
c     ... directly damp a & e
c     ... using an exponential model
c     This is based on the model present in Lee & Peale 2003, in which 
c     a rate of change for a is specified (frac{\dot{a}}{a}), and then 
c     the eccentricity damping is scaled to this (frac{\dot{e}}{e} = K frac{\dot{a}}{a})
c     Typical values might have frac{\dot{a}}{a} ~ 10^{-5} yr^{-1}  and  K ~ 100
            if (MIGRATIONMECHANISM.eq.'LP03') then
               do j = 2, nbig
c
c     Only damp if rate above zero
                     if (params(2,j).gt.0.0) then
c
c     Calculate the Jacobian coordinates of the particles (doing as per Lee & Peale )
                        ngflag_dummy = 0
                        ngf_dummy = 0
                        opt_dummy = 0
                        call mco_h2j (time,jcen_dummy,nbig,nbig,
     %                       h0,
     %                       m,x,v,xj,vj,
     %                       ngf_dummy,ngflag_dummy,opt_dummy)
c     Convert the Jacobian coordinates to elements
                        gm = m(1) + m(j)
                        call mco_x2el (gm,xj(1,j),xj(2,j),xj(3,j),
     %                       vj(1,j),vj(2,j),vj(3,j),q,e,i,p,n,l)     
c     Find the semi
                        semi = q / (1.d0 - e)
c
c     Only damp if above the minimum semi-major axis
                        if (semi.gt.params(4,j)) then
c
c     Damp the quantities
                           semi  = semi * (1 - (h0 * params(2,j))) ! here params(2,j) =  frac{\dot{a}}{a}, i.e. input yr^-1, a rate
                                                                   !     THIS IS DIFFERENT TO THE USAGE IN THE WYATT03 MODEL BELOW
                           e     = e    * (1 - (h0 * params(3,j)))
                           i     = i    * (1 - (h0 * params(3,j))) ! Added in inclination damping as well. This is NOT in LP03, as coplanar
c     Recalc q
                           q     = semi * (1.d0 - e)      
c     Convert back to Jacobian/Cartesian coordinates
                           call mco_el2x (gm,q,e,i,p,n,l,
     %                          xj(1,j),xj(2,j),xj(3,j),
     %                          vj(1,j),vj(2,j),vj(3,j))
c
                           call mco_j2h(time,jcen_dummy,nbig,nbig,
     %                          h0,
     %                          m,xj,vj,x,v,
     %                          ngf_dummy,ngflag_dummy,opt_dummy)
                        end if
                    end if
                 end do
              end if
c     
c
c
c
c
c
c
c
c THIS IS A SLIGHTLY MODIFIED VERSION OF THE LEE & PEALE MIGRATION SCHEME ABOVE.
c IT IS IMPLEMENTED SUCH THAT THE PLANETS CAN BE FORCED TO STOP JUST OUTSIDE RESONANCE.
c
            if (MIGRATIONMECHANISM.eq.'LP03STALL') then
               do j = 2, nbig
c
c     Only damp if rate above zero
                     if (params(2,j).gt.0.0) then
c
c     Calculate the Jacobian coordinates of the particles (doing as per Lee & Peale )
                        ngflag_dummy = 0
                        ngf_dummy = 0
                        opt_dummy = 0
                        call mco_h2j (time,jcen_dummy,nbig,nbig,
     %                       h0,
     %                       m,x,v,xj,vj,
     %                       ngf_dummy,ngflag_dummy,opt_dummy)
c     Convert the Jacobian coordinates to elements
                        gm = m(1) + m(j)
                        call mco_x2el (gm,xj(1,j),xj(2,j),xj(3,j),
     %                       vj(1,j),vj(2,j),vj(3,j),q,e,i,p,n,l)     
c     Find the semi
                        semi = q / (1.d0 - e)
c
c     Only damp if above the minimum semi-major axis
                        if (semi.gt.params(4,j)) then
c
c     Only damp if the period ratio between the bojects is larger than the set value
c     (for halting just outside resonance)
c     Assume that only the outer body is being driven ...
                           if (j.eq.nbig) then
                             call mco_x2el (gm,xj(1,j),xj(2,j),xj(3,j),
     %                       vj(1,j),vj(2,j),vj(3,j),
     %                       qinner,einner,iinner,pinner,ninner,linner)     
                             semiinner = qinner/(1.d0 - einner)
                             period_o = semi**1.5
                             period_i = semiinner**1.5
                           if(period_o/period_i.gt.MIGRATIONRATIO) then  

c     Damp the quantities
                                semi  = semi * (1 - (h0 * params(2,j))) ! here params(2,j) =  frac{\dot{a}}{a}, i.e. input yr^-1, a rate
                                                                        !     THIS IS DIFFERENT TO THE USAGE IN THE WYATT03 MODEL BELOW
                                e     = e    * (1 - (h0 * params(3,j)))
                                i     = i    * (1 - (h0 * params(3,j))) ! Added in inclination damping as well. This is NOT in LP03, as coplanar
c     Recalc q
                                q     = semi * (1.d0 - e)      
c     Convert back to Jacobian/Cartesian coordinates
                                call mco_el2x (gm,q,e,i,p,n,l,
     %                               xj(1,j),xj(2,j),xj(3,j),
     %                               vj(1,j),vj(2,j),vj(3,j))
c
                                call mco_j2h(time,jcen_dummy,nbig,nbig,
     %                               h0,
     %                               m,xj,vj,x,v,
     %                               ngf_dummy,ngflag_dummy,opt_dummy)
                             end if
                          end if
                       end if  
                    end if
                 end do
              end if
c
c
c
c
c
c
c
c

c
c
c     If the migration mechanism = WYATT03, then
c     ... directly damp a using the model from Wyatt 2003
c     The acceleration is parallel to the velocity and is given by
c     {\bf accn} = 0.5*kappa*(G*Mstar/a^3)^(1/2)*{\bf v}/v
c     The planets thus experience /accn/ = 0.5*kappa*(G*Mstar/a^3)^(1/2)
c     This results in a constant rate of change of semi-major axis, \dot{a}=kappa
            if (MIGRATIONMECHANISM.eq.'WYATT03') then
               do j = 2, nbig
c     
c     Only damp if rate above zero
                     if (params(2,j).gt.0.0) then
c
c     Get the r,v & a==semi (note use of Chambers' routine mco_x2a).
                        gm = m(1) + m(j)
                        call mco_x2a (gm,x(1,j),x(2,j),x(3,j),
     %                       v(1,j),v(2,j),v(3,j),semi,r,v2)
c     Find the semi
                        semi = q / (1.d0 - e)
c     Only damp if above the minimum semi-major axis
                        if (semi.gt.params(4,j)) then

c     Increment the planet's accelerations.
                           factor=sqrt(gm/(semi**3))
                           factor=0.5*params(2,j)*factor/sqrt(v2)  ! params(2,j) = Migration rate, in AU/day (input as AU / year)
c                                                                  ! THIS IS A DIFFERENT USAGE TO THE LEE & PEALE STYLE MODEL ABOVE
                           a(1,j)=a(1,j)-factor*v(1,j)
                           a(2,j)=a(2,j)-factor*v(2,j)
                           a(3,j)=a(3,j)-factor*v(3,j)
     
                     end if     
                  end if
               end do
            end if 
c
c
c
c
c     If the migration mechanism = TL2003, then
c     ... directly damp a & e using the model from Thommes & Lissauer 2003
            if (MIGRATIONMECHANISM.eq.'TL2003') then
c
c     In this model we only apply damping to the outer planet
               j = nbig
               gm = m(1) + m(j)
               call mco_x2el (gm,x(1,j),x(2,j),x(3,j),
     %                 v(1,j),v(2,j),v(3,j),q,e,i,p,n,l)
c
c     Find the semi
               semi = q / (1.d0 - e)
c
c     Where is the inner edge of the disk (r_disk)? Does it matter?!?!?!
               r_disk = 0.1       ! Temporarily making this a constant - COME BACK TO THIS LATER.
c               r_disk = params(5,j)
c               ddt_r_disk = params(2,j)
c
c     Only damp if semi > r_disk
               if (semi.gt.r_disk) then
c     Damp semi & e directly (current implementation is ignoring any a_2 - r_disk dependence)
                  oldsemi = semi
                  semi = semi - params(2,j)*h0
                  e    = e    - params(3,j)*h0
c     Convert back to Cartesian coordinates
                  q    = semi*(1.d0 - e)
                  call mco_el2x (gm,q,e,i,p,n,l,
     %                 x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j))
c
c
c               temp = (time/365.25)/100.0
c               if (  semi .lt. 0.2) then 
c                  open (unit=1,file='user.out',access='append')
c                  write(1,*) 'Debug AFTER... j= ',j,',  ',
c     %                 'x= ',x(1,j),
c     %                 'y= ',x(2,j),
c     %                 'z= ',x(3,j),
c     %                 'vx= ',v(1,j),
c     %                 'vy= ',v(2,j),
c     %                 'vz= ',v(3,j)
c                  close(unit=1)
c               end if
c
c
            end if 
         end if
c     
c
c
c
         endif
      end if
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      SIMPLE_GAS_DRAG.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c Apply gas drag as described in Adachi '76, Rafikov '04 & Mendall '07
c
c Calculate the gas density rho(R,z) at a given radius & height.
c Calculate the resulting pressure supported gas velocity.
c Calculate the relative velocity between the body & the gas.
c Calculate the resultant gas drag on the body.
c
c /data/mpayne/mathematica_models/Gas_Drag_04.nb has some theory & calcs
c
c------------------------------------------------------------------------------
c
      subroutine simple_gas_drag (time,nbod,m,x,v,a,params,algor)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbod,algor
      real*8 time,m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
      real*8 params(NPARAMS,nbod)
c
c Local
      integer j
      real*8 R_GDP,R2_GDP,cosf,sinf
      real*8 rho_g,rho_g_time,time_exp
      real*8 Z02,z2_GDP_over_Z02
      real*8 v_factor, vx_KEP, vy_KEP, vg_over_vk, vx_GAS, vy_GAS
      real*8 vx_rel,vy_rel,vz_rel
      real*8 gm,semi,q,e,i,p,n,l
c
c------------------------------------------------------------------------------
c
c     If the simulation is within the correct time range...
      if (GASDRAGSTART.lt.(time/365.25)) then 
         if (GASDRAGEND.gt.(time/365.25)) then 
c
c     We can bring some of the disk-structure calculations out of the loop...
            if (DISSIPATION_OFFSET.gt.(time/365.25)) then 
               time_exp = 1
            else
               time_exp = ((time/365.25)-DISSIPATION_OFFSET ) / TAU_DISC
               time_exp = exp(-time_exp)
            end if 
c
c     It's probably pointless doing gas drag for extremely tenuous disks
            if(time_exp.gt.GASDRAGCUTOFF) then 
c
c     Calculate the normalisation of the gas density at this time 
               rho_g_time = RHO_DISC*time_exp
c
c     Apply the drag to all bodies
               do j = 2, nbod
c      
c     Only perform the calculation if the associated mass/density/radius > 0
                  if (params(8,j) .ne. 0.d0) then    ! Have already combined rho_p * r_p into param(8,)
c     
c     (1) Distances in the Disc Plane
                     R2_GDP  = x(1,j)*x(1,j) + x(2,j)*x(2,j)
                     R_GDP   = sqrt(R2_GDP)
                     
                     cosf    = x(1,j) / R_GDP
                     sinf    = x(2,j) / R_GDP
c
c
c     (2) Gas Disk Density at the location of the particle ------------------------------------------
c         - N.B. The time-dep has been taken out of the loop
c         - N.M.B. The model used here is effectively defined in the user_parameter file...
c                   z0 = z00 r^{\beta} = \sqrt{2} h r^{\beta}, i.e. I am defining z00 = \sqrt{2} h = \sqrt{2} c_s / \omega_k
c                   \rho \propto r^{-(\alpha + \gamma)} \propto r^{-\epsilon}
c                   \Sigma \propto z(r) * r^{-\epsilon} \propto r^{-\alpha}
c
c     (Squared) Vertical scale height / structure at the location of the particle.       N.B.      Z00_DISC = \sqrt{2}*H_DISC  as above
                     Z02 = Z002_DISC * (R2_GDP**GAMMA_DISC)
c     Gas density at the radial location of the particle (radial dependence in the plane gas-disk plane)
                     rho_g = rho_g_time*(R_GDP**(-EPSILON_DISC))
c     Gas density at the vertical location of the particle (exponential decay above the plane gas-disk plane)
                     z2_GDP_over_Z02 = (x(3,j) * x(3,j))/Z02               ! Gets used again below
                     rho_g = rho_g*exp(-z2_GDP_over_Z02)
c
c
c     (3) Relative velocities
c
c     Keplerian Velocity Vector of the (Circular) gas at the location of the particle
c     Look at Section 2.3 in Murray & Dermott: v_x = -\frac{na}{1-e^2}\sin{f} = -\sqrt{\frac{G(m_1+m_2)}{a(1-e^2)}} \sin{f}, etc, etc
c     Just using m(1) as (a) it includes K2, & (b) Just assume a test particle orbit (m=0)
                     v_factor = dsqrt(m(1) / R_GDP)         ! The eccentricity is zero
                     vx_KEP = - v_factor * sinf           
                     vy_KEP =   v_factor * cosf    
c 
c     Gas Velocity Vector ( scaling Vgas/Vkepeler)
c     -->> See eqn Takeuchi & Lin '02, Eqn 7
c     -->> Note the very slight difference in factors here because of the difference between h & z0 = qsrt{2} h
                     vg_over_vk = BETA_DISC*(1+z2_GDP_over_Z02)
                     vg_over_vk = (-EPSILON_DISC-vg_over_vk) 
                     vg_over_vk = 1.0+0.25*vg_over_vk*(Z02/R2_GDP)

                     vx_GAS = vg_over_vk * vx_KEP
                     vy_GAS = vg_over_vk * vy_KEP
c     
c     Relative velocities
                     vx_rel = v(1,j) - vx_GAS
                     vy_rel = v(2,j) - vy_GAS
                     vz_rel = v(3,j)
c
c
c     Call the routine which calculates the acceleration given the relative velocities  -------------------
                     call drag_factor (j,nbod,rho_g,a,
     %                    vx_rel,vy_rel,vz_rel,R_GDP,params,time)
c
c
                  end if
               end do
            end if
         end if
      end if
c
c
c------------------------------------------------------------x------------------
c
      return
      end
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASS_GROWTH.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c Apply a given mass growth rate to increment the mass of the planet
c
c------------------------------------------------------------------------------
c
      subroutine mass_growth (m,nbig,params,time,h0)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig
      real*8 m(nbig),params(NPARAMS,nbig),time,h0
c
c Local
      integer j
c
c------------------------------------------------------------------------------
c
c     If is within the correct time range
      if (MASSGROWTHSTART.lt.(time/365.25)) then 
         if (MASSGROWTHEND.gt.(time/365.25)) then 
c
c     Get the growth rates & do the appropriate calcs for each large body
            do j = 2, nbig
c     
c     If the mass of the plaent is ABOVE, the desired limit, do nothing
               if (m(j).lt.params(10,j)) then
c
c     Increment the mass of each planet / large body
c     -->> Interpret as exponential growth timescale (years)
                  if (MASSGROWTHMODEL.eq.'EXPONENTIAL') then 
                     m(j) = m(j)*(1 + params(9,j) * h0)
                  end if
c     -->> Interpret as an absolute linear growth rate (Msun / year)
                  if (MASSGROWTHMODEL.eq.'LINEAR') then 
                     m(j) = m(j) + params(9,j) * h0 * K2
                  end if
c
               end if     
            end do
         endif
      end if
c
c------------------------------------------------------------------------------
c     
      return
      end
c
c 
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      HYBRID_TO_BS
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c If HYBRIDTOBSSWITCH = 1, then force the routine to use the B-S
c integrator IF ... some criterion is fulfilled
c
c------------------------------------------------------------------------------
c
      subroutine hybrid_to_bs (time,nbig,m,x,v,params,algor)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig,algor
      real*8 time,m(nbig),x(3,nbig),v(3,nbig)
      real*8 params(NPARAMS,nbig)
c
c Local
      integer j
      real*8 gm, semi, r, v2
c
c------------------------------------------------------------------------------
c
c
      do j = 2, nbig
c     
c     Get the r,v & a==semi (note use of Chambers' routine mco_x2a).
         gm = m(1) + m(j)
         call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),
     %        semi,r,v2)
c     
c     If the semi-major axis of the large body is < HYBRIDTOBSDIST AU, 
c     switch to BS integrator
         if (HYBRIDTOBSCRITERION.eq.'DISTANCE') then 
            if (r .lt. HYBRIDTOBSDIST) then
               algor = 2
            end if
         end if
c
      end do

c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      BS_TO_HYBRID
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c If BSTOHYBRIDSWITCH = 1, then force the routine to use the B-S
c integrator IF ... some criterion is fulfilled
c
c------------------------------------------------------------------------------
c
      subroutine bs_to_hybrid (time,algor)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer algor
      real*8 time
c
c Local
      integer j
      real*8 gm, semi, r, v2
c
c------------------------------------------------------------------------------
c      
c     If time is > BSTOHYBRIDTIME yrs, 
c     switch to BS integrator

c      open (unit=1, file='user.out', access='append')
c      write(1,*) 'In bs-to-hybrid, before, t(yr)= ',time/365.25,
c     %     ' algor= ',algor
c      close(unit=1)

         if (BSTOHYBRIDCRITERION.eq.'DISTANCE') then 
            if ( (time/365.25) .gt. BSTOHYBRIDTIME) then
               algor = 10
            end if
         end if

c      open (unit=1, file='user.out', access='append')
c      write(1,*) 'In bs-to-hybrid, after      , t(yr)= ',time/365.25,
c     %     ' algor= ',algor
c      close(unit=1)
c
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      RELATIVISTIC PRECESSION
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c If RELATIVITYSWITCH = 1, then calculate an extra relativistic precession term
c
c Until I find something better, the equations used come from
c http://www.projectpluto.com/relativi.htm
c http://books.google.com/books?printsec=frontcover&vid=ISBN0935702687&vid=LCCN91065331#v=onepage&q=&f=false
c
c------------------------------------------------------------------------------
c
      subroutine relativity (time,nbod,m,x,v,a,params)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbod
      real*8 time,m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
      real*8 params(NPARAMS,nbod)
c Local
      integer j
      real*8 r,v2,f0,f1,f2
c
c------------------------------------------------------------------------------
c     
      do j = 2, nbod
         if (params(11,j).eq.1) then 
            
c Additional (post-Newtonian) acceleration term...
c Look at Benitez & Gallardo and references there-in...http://www.springerlink.com/content/d8r26w6610389256/fulltext.pdf
c                     -GM
c Additional accel =  ---- {(4GM / r - v^2) r + 4(v.r)v}
c                     r^3c^2
c
c                     -K2
c ................ =  ---- {(4K2 / r - v^2) r + 4(v.r)v}
c                     r^3c^2
c
c c2, the speed of light (squared) in (AU / day)^2 is defined in user_parameters
c
c Position & velocity (sqrd)
            r = sqrt ( x(1,j)*x(1,j) + x(2,j)*x(2,j) + x(3,j)*x(3,j) )
            v2 = v(1,j)*v(1,j) + v(2,j)*v(2,j) + v(3,j)*v(3,j)
c Calculational Steps
            f0 = K2 / (r*r*r*c2)
            f1 = f0*((4*K2 / r) - v2)
            f2 = f0*4*(v(1,j)*x(1,j) + v(2,j)*x(2,j) + v(3,j)*x(3,j))

c Acceleration Components
            a(1,j)=a(1,j)-f1*x(1,j)-f2*v(1,j)
            a(2,j)=a(2,j)-f1*x(2,j)-f2*v(2,j)
            a(3,j)=a(3,j)-f1*x(3,j)-f2*v(3,j)
c
c      open (unit=1, file='user.out', access='append')
c      write(1,*) 'Relativity switch      , t(yr)= ',time/365.25,
c     %     ' j=',j, '  f0=',f0,
c     %     ' a...=',a(1,j),' ',a(1,j),' ',a(1,j)
c      close(unit=1)
c
         end if
      end do

c
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      ECCENTRIC_GAS_DRAG.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Ji-Wei Xei & Matt Payne 
c
c Apply updated gas drag model
c Alterations by Ji-Wei allow basic modelling of eccentric discs
c
c New gas drag model can include gas-disk's
c 1)eccentricity 
c 2)inclination
c 3)precession of periastron
c 4)precession of ascend node
c
c Calculate the gas density rho(R,z) at a given radius & height.
c Calculate the resulting pressure supported gas velocity.
c Calculate the relative velocity between the body & the gas.
c Calculate the resultant gas drag on the body.
c
c------------------------------------------------------------------------------
c
      subroutine eccentric_gas_drag (time,nbod,m,x,v,a,params)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbod
      real*8 time,m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
      real*8 params(NPARAMS,nbod)
c
c Local
      integer j,k
      real*8 rho_g,time_exp
      real*8 r_disc,wk,pr_rate
      real*8 dN,dW,cosdN,sindN,cosdW,sindW
      real*8 x_temp,x_GDP,y_GDP,z_GDP,z2_GDP,R_GDP,R2_GDP
      real*8 cosf,sinf
      real*8 Z02,z2_GDP_over_Z02
      real*8 ad,ed,cb,cc,cn,cr,cta
      real*8 v_factor,vx_KEP,vy_KEP,vg_over_vk
      real*8 vx_temp,vy_temp,vx_GAS,vy_GAS,vz_GAS
      real*8 vx_rel, vy_rel, vz_rel
c
c------------------------------------------------------------------------------
c
c
c     (a) - Topline Switches  & Time Dependencies -------------------------------------------------------------
      if (GASDRAGSTART.lt.(time/365.25)) then  ! If we are within the time slot...  
         if (GASDRAGEND.gt.(time/365.25)) then ! ...where we want to model gas drag, then do calculation
c            
c     May not want to have a dissipating gas disc, if so, skip dissipation calculations
            if (TIMEDEPSWITCH.eq.0) then 
               time_exp=1.0
            else
               if (DISSIPATION_OFFSET.gt.(time/365.25)) then ! If we are still in the constant surf-density zone
                  time_exp=1.0
               else            
                  time_exp =((time/365.25)-DISSIPATION_OFFSET)/TAU_DISC
                  time_exp =exp(-time_exp)
               end if 
            end if
c      Time-Dependent Density Normalisation in the Mid-plane of the Gas Disk
            rho_g = RHO_DISC*time_exp
c
c     It's probably pointless doing gas drag for extremely tenuous disks
            if(time_exp.gt.GASDRAGCUTOFF) then 
c
c
c     (b) - Precession Angles / Rates --------------------------------------------------------------------------
               r_disc = dsqrt(x(1,j)*x(1,j) + x(2,j)*x(2,j))/COS_IBG              ! Radial distance in the gas disc. 
c           Precession Rate of Ascending Node
               if(I_BIN_GAS.lt.MIN_INC) then
                  dN=0.0
               else
                  if (NODE_PRECESS_MODEL.eq.0) then 
                     dN=dN_ARBITRARY
                  else
                     if (NODE_PRECESS_MODEL.eq.1) then 
                  !Precession rate is from Lawood et al. 1996 (Eqn 22)
                  !The constants are all absorbed into the paramter LARWOOD -->> See the user_parameter.inc file
c                 pr_rate=-0.75*tpi*cMb/cMa*(r_disc/aB)**3*SQRT_MASS_PRIM/(r_disc**1.5)  ! Thebault et al. 2006  
                        pr_rate = -LARWOOD/(r_disc**1.5)                                       
                        dN=dN_ARBITRARY*pr_rate*time+dN_INITIAL                                        
                     endif
                  end if
               end if
               cosdN = dcos(dN)
               sindN = dsin(dN)

c           Precession Rate of Periastron
               if (PERIASTRON_PRECESS_MODEL.eq.0) then 
                  dW = dW_ARBITRARY*pr_rate*time+dW_INITIAL
               else
                  if (PERIASTRON_PRECESS_MODEL.eq.1) then ! This is a fit to the Paardekooper '08 paper
                     if (r_disc.le.1)             dW = 0.0                               
                     if (r_disc.gt.1.and.r_disc.lt.4) then 
                        dW = -0.6*(r_disc-1.0)+5.0                                       
                     endif
                     if (r_disc.gt.4)             dW = PI                                
                  else
                     if (PERIASTRON_PRECESS_MODEL.eq.2) then 
                        dW=-(dW_ARBITRARY/PRECESS_TIMESCALE)
     %                       *TWOPI*time/365.25
                     end if
                  end if
               end if
               cosdW = dcos(dW)
               sindW = dsin(dW)
c
c
c     ________________________________________________________________________________________________________
c     Start the loop section of the calculation
c     --------------------------------------------------------------------------------------------------------
c     
               do j = 2, nbod                     ! Apply the drag to all bodies
c     
                  if (params(8,j) .ne. 0.d0) then ! Only perform the calculation if the associated mass/density/radius > 0
                                                  ! Have already combined rho_p * r_p into param(8,j)
c          
c
c     (1) - Coordinate Conversions -----------------------------------------------------------------------------
c           Want to calculate the coords of the solid body relative to the gas-disk plane (x_GDP, etc etc)
                     if (cosdN.ne.1.0.OR.COS_IBG.ne.1.0) then 
                        x_GDP = (x(1,j)*cosdN+x(2,j)*sindN)
                        y_GDP = (-x(1,j)*sindN+x(2,j)*cosdN)*COS_IBG 
     1                       + x(3,j)*SIN_IBG
                        z_GDP = (-x(1,j)*sindN+x(2,j)*cosdN)*(-SIN_IBG) 
     1                       + x(3,j)*COS_IBG
                     else 
                        x_GDP = x(1,j)
                        y_GDP = x(2,j)
                        z_GDP = x(3,j)
                     end if
                     z2_GDP=z_GDP*z_GDP    ! This is the quantity that is used below, never z_GDP
                     
                     if (cosdW.ne.1.0) then 
                        x_temp= x_GDP*cosdW + y_GDP*sindW
                        y_GDP =-x_GDP*sindW + y_GDP*cosdW
                        x_GDP = x_temp
                     end if
                     
                     R2_GDP  = x_GDP**x_GDP + y_GDP**y_GDP
                     R_GDP   = sqrt(R2_GDP)
                     cosf    = x_GDP / R_GDP
                     sinf    = y_GDP / R_GDP
c
c
c     (2) - Gas Disk Density at the location of the particle ------------------------------------------
c         - N.B. The time-dep has been taken out of the loop
c         - N.M.B. The model used here is effectively defined in the user_parameter file...
c                   z0 = z00 r^{\beta} = \sqrt{2} h r^{\beta}, i.e. I am defining z00 = \sqrt{2} h = \sqrt{2} c_s / \omega_k
c                   \rho \propto r^{-(\alpha + \gamma)} \propto r^{-\epsilon}
c                   \Sigma \propto z(r) * r^{-\epsilon} \propto r^{-\alpha}
c
c     (Squared) Vertical scale height / structure at the location of the particle.       N.B.      Z00_DISC = \sqrt{2}*H_DISC  as above
                     Z02 = Z002_DISC * (R2_GDP**GAMMA_DISC)
c     Gas density at the radial location of the particle (radial dependence in the plane gas-disk plane)
                     rho_g = rho_g*(R_GDP**(-EPSILON_DISC))
c     Gas density at the vertical location of the particle (exponential decay above the plane gas-disk plane)
                     z2_GDP_over_Z02 = z2_GDP/Z02               ! Gets used again below
                     rho_g = rho_g*exp(-z2_GDP_over_Z02)
c
c
c
c     (3) - Different Models for Disk Eccentricity -------------------------------------------------------------
c
                     if (GAS_ECCENTRICITY_MODEL.eq.-2) then ! Superbee, Paardekooper '08
                        ad = R_GDP
                        if(ad.le.0.5) ed=0.0
                        if(ad.gt.0.5.and.ad.le.3.0) ed=0.09*(ad-0.5)
                        if(ad.gt.3.0) ed=0.225-0.03*(ad-3.0)
                        ad = R_GDP * (1+ed*cosf)/(1-ed*ed)
                     elseif (GAS_ECCENTRICITY_MODEL.eq.-1) then ! Minmod, Paardekooper '08
                        ad = R_GDP
                        if(ad.le.1.0) ed=0.0
                        if(ad.gt.1.0) ed=0.007*(ad-1)
                        ad = R_GDP * (1+ed*cosf)/(1-ed*ed)
                     elseif (GAS_ECCENTRICITY_MODEL.eq.0) then ! Circular Gas Model
                        ed = 0.0
                        ad = R_GDP  !  R_GDP * (1-ed*cosf) / (1 - ed*ed)        
                     end if
c
c
c     (4) - Relative Velocities --------------------------------------------------------------------------------
c     
c     Keplerian Velocity Vector of the (Eccentric) gas at the location of the particle
c     Look at Section 2.3 in Murray & Dermott: v_x = -\frac{na}{1-e^2}\sin{f} = -\sqrt{\frac{G(m_1+m_2)}{a(1-e^2)}} \sin{f}, etc, etc
                     v_factor = dsqrt(m(1) / (ad * (1.0 - ed * ed)))
                     vx_KEP = - v_factor * sinf          
                     vy_KEP =   v_factor * (ed + cosf)   
c     
c     Gas Velocity Vector ( scaling Vgas/Vkepeler)
c     -->> See eqn Takeuchi & Lin '02, Eqn 7
c     -->> Note the very slight difference in factors here because of the difference between h & z0 = qsrt{2} h
                     vg_over_vk = BETA_DISC*(1+z2_GDP_over_Z02)
                     vg_over_vk = (-EPSILON_DISC-vg_over_vk) 
                     vg_over_vk = 1.0+0.25*vg_over_vk*(Z02/R2_GDP)
c
c     Convert velocity (of the gas) in gas disk coordinates to Binary Plane Coords (BPC).
c     N.B. v_KEP * vg_over_vk gives the gas velocity in the GDP, then the cosdW .. & cosdN converts the gas velocity to the BPC
                     vx_GAS = (vx_KEP*cosdW - vy_KEP*sindW) * vg_over_vk
                     vy_GAS = (vx_KEP*sindW + vy_KEP*cosdW) * vg_over_vk
                     vz_GAS = 0.d0
               
                     vx_temp = vx_GAS*cosdN                              ! Temporary step required to get to vx_GAS below
     %                         -(vy_GAS*COS_IBG-vz_GAS*SIN_IBG)*sindN                     
                     vy_temp = vx_GAS*sindN                     
     %                         + (vy_GAS*COS_IBG+vz_GAS*SIN_IBG)*cosdN

                     vz_GAS = vy_GAS*SIN_IBG+vz_GAS*COS_IBG
                     vy_GAS = vy_temp
                     vx_GAS = vx_temp

c     Relative Velocity of the particle w.r.t. the gas
                     vx_rel = v(1,j) - vx_GAS
                     vy_rel = v(2,j) - vy_GAS
                     vz_rel = v(3,j) - vz_GAS

c     (5) - Call the routine which calculates the acceleration given the relative velocities  -------------------
                     call drag_factor (j,nbod,rho_g,a,
     %                    vx_rel,vy_rel,vz_rel,R_GDP,params,time)
c
c
c
c     
                  end if
               end do 
            end if
         end if
      end if
c
c
c------------------------------------------------------------x------------------
c
      return
      end
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      DRAG_FACTOR.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c Calculates the acceleration on a particle due to gas drag.
c
c Uses the drag coefficient method of Woitke '03, Paardekooper '07 & Lyra '09
c  -->> This means that we don't have abrupt regime change, just have smooth extrapolation
c
c Requires as an input the velocities (vx_rel,vy_rel,vz_rel) of the particle relative to 
c the gas at the location of the particle
c
c------------------------------------------------------------------------------
c
      subroutine DRAG_FACTOR (j,nbod,rho_g,a,vx_rel,vy_rel,vz_rel,
     %     R,params,time)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer j,nbod
      real*8 rho_g
      real*8 a(3,nbod)
      real*8 vx_rel,vy_rel,vz_rel,R,params(NPARAMS,nbod)
      real*8 time
c
c Local
      real*8 lambda,Kn3,c_s,Ma,Re,C_D_E,C_D_S,C_D,factor
      real*8 v_rel,v2_rel
c
c------------------------------------------------------------------------------

c     (1) Magnitude of the relative velocities
      v2_rel = vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel
      v_rel  = dsqrt(v2_rel)

c     (2) Set up the various quantities
      lambda = MU_OVER_SIGMA / rho_g                        ! MU_OVER_SIGMA is in [\Msun / AU^2] from user_parameters
      Kn3 = 1.5*lambda / params(7,nbod)                     ! Knudsen Number, Kn = lambda / (2 * radius), again everything is in AU, etc, etc
      c_s = C_S0 * R                                        ! Sound Speed (C_S0 is in AU / day as required)
      Ma = v_rel / c_s                                      ! Mach Number,    Ma
      Re = SIXSQRTPIOVER8 * Ma * (params(7,nbod)/lambda)
      C_D_E = 2*sqrt(1+(ONE28_OVER_9PI/(Ma*Ma)))
      if (Re.lt.500) then 
         C_D_S = (24 / Re) + 3.6*(Re**(-0.313))
      else if (Re.lt.1500) then 
         C_D_S = 0.000095*(Re**1.397)
      else
         C_D_S = 2.61
      end if
      C_D = (Kn3*Kn3*C_D_E + C_D_S) /((Kn3 +1)*(Kn3 +1))                 

c     (3) - Gas Drag Acceleration Factor -----------------------------------------------------------------------
      factor = -0.375*rho_g*C_D*v_rel/params(8,j)        ! Have already combined rho_p * r_p into param(8,)

c     (4) - Implement the acceleration Component ---------------------------------------------------------------
      a(1,j)=a(1,j)+factor*vx_rel
      a(2,j)=a(2,j)+factor*vy_rel
      a(3,j)=a(3,j)+factor*vz_rel

c------------------------------------------------------------x------------------
c
      return
      end
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      OLD_GAS_DRAG.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c Apply gas drag as described in Adachi '76, Rafikov '04 & Mendall '07
c
c Calculate the gas density rho(R,z) at a given radius & height.
c Calculate the resulting pressure supported gas velocity.
c Calculate the relative velocity between the body & the gas.
c Calculate the resultant gas drag on the body.
c
c /data/mpayne/mathematica_models/Gas_Drag_04.nb has some theory & calcs
c
c------------------------------------------------------------------------------
c
      subroutine old_gas_drag (time,nbod,m,x,v,a,params)
c
      implicit none
      include 'mercury.inc'
      include 'user_parameters.inc'
c
c Input/Output
      integer nbod
      real*8 time,jcen(3),m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
      real*8 params(NPARAMS,nbod)
c
c Local
      integer j,k
      real*8 m_p,rhor
      real*8 RCYL,RCYL2,rsph,z2,Z2_0,rho_g
      real*8 v_G,v_Gx,v_Gy,v_relx, v_rely, v_relz,v_rel,factor
      real*8 vp,v2,one_over_a,semi,hx,hy,hz,h2,e
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,gm
      real*8 t13,t14,t15,t16,t17,t18,t19,t20,q,i,p,n,l
c
c------------------------------------------------------------------------------
c
c     Apply the drag to all bodies
      do j = 2, nbod
c
c     
c     Have already combined rho_p * r_p into param(4,)
         rhor = params(8,j)
c     
c     Only perform the calculation if the associated mass/density/radius > 0
         if (rhor .ne. 0.d0) then
c
c     Get the distances (Polar & Spherical) & velocities
            RCYL2 = x(1,j)**2+x(2,j)**2
            RCYL = sqrt(RCYL2)
            z2 = x(3,j)**2
            rsph = sqrt(RCYL2+z2)
c     
c     Shouldn't go hyperbolic, but might as well be careful
            one_over_a=2.d0/rsph-v2/m(1)
            if (one_over_a .lt. 0.d0) then
               one_over_a=-1.d0*one_over_a
            end if
c     
c     Vertical scale height / structure (squared)
            Z2_0 = Z002_DISC*(RCYL2**GAMMA_DISC)
c     
c     General density structure
            rho_g = RHO_DISC*(RCYL**(-EPSILON_DISC))*exp(-z2/z2_0)
c     Time dep density 
            rho_g=rho_g*exp(-(DISSIPATION_OFFSET+time/365.25)/TAU_DISC)
c     
c     Altered Keplerian Gas Velocity - only horizontal component of solar gravity
            gm = m(1) + m(j)
            v_G = sqrt(gm*RCYL2*(rsph**(-3.0)))
c     Pressure supported gas velocity (Sigma~R^-alpha, T~R^-beta)
c     alpha = epsilon - gamma
c     alpha = 1.5 & beta = 0.5 in Ida & Lin '04
c     v_G=v_KEP(1-eta), eta = ETA_DISC*(z2_0/RCYL2)
c            v_G = v_G*(1- ETA_DISC*(z2_0/RCYL2))
c     Assume velocity structure of gas is cylindrical (+/-ve signs!)
            v_Gx = v_G*(-x(2,j)/RCYL)
            v_Gy = v_G*(x(1,j)/RCYL)
c     
c     Relative velocities
            v_relx = v(1,j) - v_Gx
            v_rely = v(2,j) - v_Gy
            v_relz = v(3,j)
            v_rel = sqrt(v_relx**2+v_rely**2+v_relz**2)
c     
c     Characteristic acceleration (/v(i,j))
            factor = -0.2*rho_g*v_rel/(rhor)
c
c     Acceleration Components
c            a(1,j)=0
c            a(2,j)=0
c            a(3,j)=0
            a(1,j)=a(1,j)+factor*v_relx
            a(2,j)=a(2,j)+factor*v_rely
            a(3,j)=a(3,j)+factor*v_relz
c
c     Debugging the code
c                     gm = m(1) + m(j)
c                     call mco_x2el (gm,x(1,j),x(2,j),x(3,j),
c     %                    v(1,j),v(2,j),v(3,j),q,e,i,p,n,l)
c                     semi = q / (1-e)
c      if ((time/365.25).lt.100) then 
c      open (unit=1,file='user.out',access='append')
c                     open (unit=1,file='user.out',access='append')
c                     write(1,*) 'TIME =  ',time/365.25,' j =  ',j
c                     write(1,*) 'Old Drag Model, '
c                     write(1,*) 'XYZ, ',
c     %                    ' x(1,j)= ',x(1,j),
c     %                    ' x(2,j)= ',x(2,j),
c     %                    ' x(3,j)= ',x(3,j)
c                     write(1,*) 'UVW, ',
c     %                    ' v(1,j)= ',v(1,j),
c     %                    ' v(2,j)= ',v(2,j),
c     %                    ' v(3,j)= ',v(3,j)
c                     write(1,*) 'Elemenets, ',
c     %                    ' a=',semi,' e=',e,' i=',i
c                     write(1,*) 'Relative Velocity, ',
c     %                    ' v_relx= ',v_relx,
c     %                    ' v_rely= ',v_rely,
c     %                    ' v_relz= ',v_relz
c                     write(1,*) ' factor= ',factor
c                     write(1,*) 'Acceleration, ',
c     %                    ' a(1,j)= ',a(1,j),
c     %                    ' a(2,j)= ',a(2,j),
c     %                    ' a(3,j)= ',a(3,j)
c                     close(unit=1)
c                  end if                     

         end if
      end do
c
c
c------------------------------------------------------------x------------------
c
      return
      end
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      TIDAL DAMPING
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne 
c
c If TIDALSWITCH = 1, then calculate tidal intereactions between planets
c
c Until I find something better, the equations used come Mardling & Lin 2002
c
c------------------------------------------------------------------------------
c
      subroutine tidal (time,nbig,m,x,v,a,params,s,rphys,h0)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig
      real*8 time,m(nbig),x(3,nbig),v(3,nbig),a(3,nbig),s(3,nbig)
      real*8 params(NPARAMS,nbig),rphys(nbig),h0
c Local
      integer j,k
      real*8 gm,r,v2,semi,r_3,r_4,r_8,S1_5,Sj_5
      real*8 massratio1,massratioj,massratio2,semi_3,minus_six_n
      real*8 runit(3),rdotOmega1,rdotOmega1_2,rdotOmegaj,rdotOmegaj_2
      real*8 Omega1_2,Omegaj_2,rdotv,rcvmrO1(3),rcvmrOj(3)
      real*8 rcvmrOcr1(3),rcvmrOcrj(3)
      real*8 temp1,tempj,factor1,factorj,f_1(3),f_j(3)
      real*8 xrel(3),vrel(3),h0_mu_over_MoI1,h0_mu_over_MoI2
      real*8 h0_mu_over_MoIj,mu1j,mu2j,rcf1(3),rcfj(3)
c
c------------------------------------------------------------------------------
c     
c      open (unit=1,file='user.out',access='append')
c      write(1,*) '                     '
c      close(unit=1)
c      open (unit=1,file='user.out',access='append')
c      write(1,*) 'Debug tidal damping routine (start),',
c     %       'a(1,2)= ',a(1,2),
c     %       'a(2,2)= ',a(2,2),
c     %       'a(3,2)= ',a(3,2)
c      close(unit=1)
c
c     DO CALCULATIONS BETWEEN ALL BODIES AND THE PRIMARY
      do j = 2, nbig
c
c     Get the r,v & a==semi (note use of Chambers' routine mco_x2a).
         gm = m(1) + m(j)
         call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),
     %        semi,r,v2)
c
c     CHECK WHETHER THE SEPARATION IS LARGER THAN THE CUT-OFF DISTANCE
c         open (unit=1,file='user.out',access='append')
c         write(1,*) 'Debug tidal damping (start primary) j= ',j,',  ',
c     %        'r= ',r,
c     %        'TIDALCUT= ',TIDALCUT,
c     %       'a(1,2)= ',a(1,2),
c     %       'a(2,2)= ',a(2,2),
c     %       'a(3,2)= ',a(3,2),
c     %        's(1,j)= ',s(1,j),
c     %        's(2,j)= ',s(2,j),
c     %        's(3,j)= ',s(3,j)
c         close(unit=1)
c
c
c
         if (r.lt.TIDALCUT) then 
c
c     USEFUL SCALAR QUANTITIES
            r_3 = r*r*r
            r_4 = r_3*r
            r_8 = r_4 * r_4
            S1_5 = rphys(1)*rphys(1)*rphys(1)*rphys(1)*rphys(1)
c            S1_5 = params(7,1)*params(7,1)*params(7,1)*params(7,1)
c            S1_5 = S1_5 * params(7,1)
            Sj_5 = rphys(j)*rphys(j)*rphys(j)*rphys(j)*rphys(j)
c            Sj_5 = params(7,j)*params(7,j)*params(7,j)*params(7,j)
c            Sj_5 = Sj_5 * params(7,j)
            massratioj = m(j) / m(1)
            massratio1 = m(1) / m(j)
            mu1j = m(1) + m(j) / (m(1) + m(j))
            semi_3 = semi*semi*semi
            minus_six_n = -6.0*sqrt((m(1)+m(j))/semi_3)
c
c     USEFUL VECTOR QUANTITIES
c
c     \frac{r}{/r/}
            runit(1) = x(1,j) / r
            runit(2) = x(2,j) / r
            runit(3) = x(3,j) / r
c     r \dot{} \Omega
c         open (unit=1,file='user.out',access='append')
c         write(1,*) '         rdotOmega1 = ',rdotOmega1
c         close(unit=1)
            rdotOmega1 = s(1,1)*runit(1) + s(2,1)*runit(2) 
c         open (unit=1,file='user.out',access='append')
c         write(1,*) '         rdotOmega1 = ',rdotOmega1
c         close(unit=1)
            rdotOmega1 = rdotOmega1 + s(3,1)*runit(3)
c         open (unit=1,file='user.out',access='append')
c         write(1,*) '         rdotOmega1 = ',rdotOmega1
c         close(unit=1)
            rdotOmega1_2 = rdotOmega1 * rdotOmega1
            rdotOmegaj = s(1,j)*runit(1) + s(2,j)*runit(2) 
            rdotOmegaj = rdotOmegaj + s(3,j)*runit(3)
            rdotOmegaj_2 = rdotOmegaj * rdotOmegaj
c     \Omega^2
            Omega1_2 = s(1,1)*s(1,1) + s(2,1)*s(2,1) + s(3,1)*s(3,1)
            Omegaj_2 = s(1,j)*s(1,j) + s(2,j)*s(2,j) + s(3,j)*s(3,j)
c            ***** NEED TO DEFINE THE SPIN VECTORS / GET THEM FROM THE INPUT *****************************************
c            ************** ALSO NEED TO CHECK THE UNItS THAT THE SPIN VECTORS ARe IN **************************8
c
c     r \dot{} v
            rdotv = runit(1)*v(1,j) + runit(2)*v(2,j) + runit(3)*v(3,j)
c     r \cross v - r\Omega
            rcvmrO1(1) = runit(2)*v(3,j) - runit(3)*v(2,j) - r*s(1,1)
            rcvmrO1(2) = runit(3)*v(1,j) - runit(1)*v(3,j) - r*s(2,1)
            rcvmrO1(3) = runit(1)*v(2,j) - runit(2)*v(1,j) - r*s(3,1)
            rcvmrOj(1) = runit(2)*v(3,j) - runit(3)*v(2,j) - r*s(1,j)
            rcvmrOj(2) = runit(3)*v(1,j) - runit(1)*v(3,j) - r*s(2,j)
            rcvmrOj(3) = runit(1)*v(2,j) - runit(2)*v(1,j) - r*s(3,j)
c     (r \cross v - r\Omega) \cross r
            rcvmrOcr1(1) = rcvmrO1(2)*runit(3) - rcvmrO1(3)*runit(2)
            rcvmrOcr1(2) = rcvmrO1(3)*runit(1) - rcvmrO1(1)*runit(3)
            rcvmrOcr1(3) = rcvmrO1(1)*runit(2) - rcvmrO1(2)*runit(1)
            rcvmrOcrj(1) = rcvmrOj(2)*runit(3) - rcvmrOj(3)*runit(2)
            rcvmrOcrj(2) = rcvmrOj(3)*runit(1) - rcvmrOj(1)*runit(3)
            rcvmrOcrj(3) = rcvmrOj(1)*runit(2) - rcvmrOj(2)*runit(1)
c
c     CALCULATE THE FORCES (Acceleration requires both the star-planet and the planet-star terms)
c     Quadrupole Terms -----------------------------------------------------------------------------------------------------------------
            temp1 = S1_5*(1 + massratioj)*params(12,1) / r_4
            tempj = Sj_5*(1 + massratio1)*params(12,j) / r_4
c            ***** NEED TO DEFINE LOVE NUMBER IN THE params() ARRAY AND SET THAT UP IN THE INITIALIZATION ROUTINE ****
c
            factor1 = 5*rdotOmega1_2 - Omega1_2 - (12 * m(j) / (r_3))
            factorj = 5*rdotOmegaj_2 - Omegaj_2 - (12 * m(1) / (r_3))
c
c     This is the factor f_QD from the Mardling & Lin, 2002, Eqn3
            f_1(1) = temp1*(factor1*runit(1) -2*rdotOmega1*s(1,1))
            f_1(2) = temp1*(factor1*runit(2) -2*rdotOmega1*s(2,1))
            f_1(3) = temp1*(factor1*runit(3) -2*rdotOmega1*s(3,1))
            f_j(1) = tempj*(factorj*runit(1) -2*rdotOmegaj*s(1,j))
            f_j(2) = tempj*(factorj*runit(2) -2*rdotOmegaj*s(2,j))
            f_j(3) = tempj*(factorj*runit(3) -2*rdotOmegaj*s(3,j))
            
c
c
c     Damping Terms ---------------------------------------------------------- ---------------------------------------------------------
            temp1 = minus_six_n*params(13,1)*massratioj*S1_5*semi_3/r_8
            tempj = minus_six_n*params(13,j)*massratio1*Sj_5*semi_3/r_8    
c
c     Adding in the f_TF term from Mardling & Lin 2002, Eqn 4
            f_1(1) = f_1(1) + temp1*(3*rdotv*runit(1) + rcvmrOcr1(1))
            f_1(2) = f_1(2) + temp1*(3*rdotv*runit(2) + rcvmrOcr1(2))
            f_1(3) = f_1(3) + temp1*(3*rdotv*runit(3) + rcvmrOcr1(3))
            f_j(1) = f_j(1) + temp1*(3*rdotv*runit(1) + rcvmrOcrj(1))
            f_j(2) = f_j(2) + temp1*(3*rdotv*runit(2) + rcvmrOcrj(2))
            f_j(3) = f_j(3) + temp1*(3*rdotv*runit(3) + rcvmrOcrj(3))

            a(1,j) = a(1,j) +  f_1(1)
            a(2,j) = a(2,j) +  f_1(2)
            a(3,j) = a(3,j) +  f_1(3)
            a(1,j) = a(1,j) +  f_j(1)
            a(2,j) = a(2,j) +  f_j(2)
            a(3,j) = a(3,j) +  f_j(3)
c
c
c         open (unit=1,file='user.out',access='append')
c         write(1,*) 'Debug tidal damping routine (main)... j= ',j,',  ',
c     %        'a(1,j)= ',a(1,j),
c     %        'a(2,j)= ',a(2,j),
c     %        'a(3,j)= ',a(3,j),
c     %        's(1,j)= ',s(1,j),
c     %        's(2,j)= ',s(2,j),
c     %        's(3,j)= ',s(3,j),
c     %        '                     s(1,1)= ',s(1,1),
c     %        's(2,1)= ',s(2,1),
c     %        's(3,1)= ',s(3,1),
c     %        'temp1 = ',temp1,
c     %        'tempj = ',tempj,
c     %        'S1_5  = ',S1_5,
c     %        'rdotv = ',rdotv,
c     %        '         rdotOmega1 = ',rdotOmega1,
c     %        'rdotOmega1_2 = ',rdotOmega1_2,
c     %        'rdotOmegaj_2 = ',rdotOmegaj_2,
c     %        ' r = ',r,'        runit = ',runit
c         close(unit=1)
c
c     Evolve the spin...
c           r \cross f (See eqn 15, Mardling & Lin 2002)
            rcf1(1) = x(2,j)*f_1(3) - x(3,j)*f_1(2)
            rcf1(2) = x(3,j)*f_1(1) - x(1,j)*f_1(3)
            rcf1(3) = x(1,j)*f_1(2) - x(2,j)*f_1(1)
            rcfj(1) = x(2,j)*f_j(3) - x(3,j)*f_j(2)
            rcfj(2) = x(3,j)*f_j(1) - x(1,j)*f_j(3)
            rcfj(3) = x(1,j)*f_j(2) - x(2,j)*f_j(1)

            h0_mu_over_MoI1 = h0*mu1j / params(14,1) 
            h0_mu_over_MoIj = h0*mu1j / params(14,j) 

            s(1,1) = s(1,1)  - h0_mu_over_MoI1 * rcf1(1)
            s(2,1) = s(2,1)  - h0_mu_over_MoI1 * rcf1(2)
            s(3,1) = s(3,1)  - h0_mu_over_MoI1 * rcf1(3)
            s(1,j) = s(1,j)  - h0_mu_over_MoIj * rcfj(1)
            s(2,j) = s(2,j)  - h0_mu_over_MoIj * rcfj(2)
            s(3,j) = s(3,j)  - h0_mu_over_MoIj * rcfj(3)
c
         end if
      end do
c
c
c
c         open (unit=1,file='user.out',access='append')
c         write(1,*) 'Debug tidal damping routine (middle) j= ',j,',  ',
c     %        'a(1,2)= ',a(1,2),
c     %        'a(2,2)= ',a(2,2),
c     %        'a(3,2)= ',a(3,2),
c     %        '   nbig = ',nbig
c         close(unit=1)
c
c
c     Check whether there is a binary star in the simulation
c     If so, calculate the tidal interactions between this and all other bodies (except primary)
      if (TIDALBINARYSWITCH.eq.1) then
c
         do j = 3,nbig
c     
            xrel(1) = x(1,j) - x(1,2)
            xrel(2) = x(2,j) - x(2,2)
            xrel(3) = x(3,j) - x(3,2)
            vrel(1) = v(1,j) - v(1,2)
            vrel(2) = v(2,j) - v(2,2)
            vrel(3) = v(3,j) - v(3,2)
            gm = m(2) + m(j)
            call mco_x2a (gm,xrel(1),xrel(2),xrel(3),
     %           vrel(1),vrel(2),vrel(3),semi,r,v2)
c     
            if (r.lt.TIDALCUT) then 
c     
c     USEFUL SCALAR QUANTITIES
               r_3 = r*r*r
               r_4 = r_3*r
               r_8 = r_4 * r_4
               S1_5 = rphys(1)*rphys(1)*rphys(1)*rphys(1)*rphys(1)
c     S1_5 = params(7,1)*params(7,1)*params(7,1)*params(7,1)
c     S1_5 = S1_5 * params(7,1)
               Sj_5 = rphys(j)*rphys(j)*rphys(j)*rphys(j)*rphys(j)
c     Sj_5 = params(7,j)*params(7,j)*params(7,j)*params(7,j)
c     Sj_5 = Sj_5 * params(7,j)
               massratioj = m(j) / m(2)
               massratio2 = m(2) / m(j)
               mu2j = m(2) + m(j) / (m(2) + m(j))
               semi_3 = semi*semi*semi
               minus_six_n = -6.0*sqrt((m(2)+m(j))/semi_3)
c     
c     USEFUL VECTOR QUANTITIES
c     
c     \frac{r}{/r/}
               runit(1) = xrel(1) / r
               runit(2) = xrel(2) / r
               runit(3) = xrel(3) / r
c     r \dot{} \Omega
               rdotOmega1 = s(1,2)*runit(1) + s(2,2)*runit(2) 
               rdotOmega1 = rdotOmega1 + s(3,2)*runit(3)
               rdotOmega1_2 = rdotOmega1 * rdotOmega1
               rdotOmegaj = s(1,j)*runit(1) + s(2,j)*runit(2) 
               rdotOmegaj = rdotOmegaj + s(3,j)*runit(3)
               rdotOmegaj_2 = rdotOmegaj * rdotOmegaj
c     \Omega^2
               Omega1_2 = s(1,2)*s(1,2) + s(2,2)*s(2,2) + s(3,2)*s(3,2)
               Omegaj_2 = s(1,j)*s(1,j) + s(2,j)*s(2,j) + s(3,j)*s(3,j)
c     ***** NEED TO DEFINE THE SPIN VECTORS / GET THEM FROM THE INPUT *****************************************
c     r \dot{} v
               rdotv =runit(1)*v(1,j) +runit(2)*v(2,j) +runit(3)*v(3,j)
c     r \cross v - r\Omega
               rcvmrO1(1)= runit(2)*vrel(3) -runit(3)*vrel(2) -r*s(1,2)
               rcvmrO1(2)= runit(3)*vrel(1) -runit(1)*vrel(3) -r*s(2,2)
               rcvmrO1(3)= runit(1)*vrel(2) -runit(2)*vrel(1) -r*s(3,2)
               rcvmrOj(1)= runit(2)*vrel(3) -runit(3)*vrel(2) -r*s(1,j)
               rcvmrOj(2)= runit(3)*vrel(1) -runit(1)*vrel(3) -r*s(2,j)
               rcvmrOj(3)= runit(1)*vrel(2) -runit(2)*vrel(1) -r*s(3,j)
c     (r \cross v - r\Omega) \cross r
               rcvmrOcr1(1) = rcvmrO1(2)*runit(3) - rcvmrO1(3)*runit(2)
               rcvmrOcr1(2) = rcvmrO1(3)*runit(1) - rcvmrO1(1)*runit(3)
               rcvmrOcr1(3) = rcvmrO1(1)*runit(2) - rcvmrO1(2)*runit(1)
               rcvmrOcrj(1) = rcvmrOj(2)*runit(3) - rcvmrOj(3)*runit(2)
               rcvmrOcrj(2) = rcvmrOj(3)*runit(1) - rcvmrOj(1)*runit(3)
               rcvmrOcrj(3) = rcvmrOj(1)*runit(2) - rcvmrOj(2)*runit(1)
c     
c     CALCULATE THE FORCES (Acceleration requires both the star-planet and the planet-star terms)
c     Quadrupole Terms -----------------------------------------------------------------------------------------------------------------
               temp1 = S1_5*(1 + massratioj)*params(12,1) / r_4
               tempj = Sj_5*(1 + massratio2)*params(12,j) / r_4
c     ***** NEED TO DEFINE LOVE NUMBER IN THE params() ARRAY AND SET THAT UP IN THE INITIALIZATION ROUTINE ****
c     
               factor1 = 5*rdotOmega1_2 - Omega1_2 - (12 * m(j) / (r_3))
               factorj = 5*rdotOmegaj_2 - Omegaj_2 - (12 * m(2) / (r_3))
c     CHECK THAT THE G TErM HAS BEEN PROPErLY ACOUNTED FOR ...
c
c     This is the factor f_QD from the Mardling & Lin, 2002, Eqn3
               f_1(1) = temp1*(factor1*runit(1) -2*rdotOmega1*s(1,2))
               f_1(2) = temp1*(factor1*runit(2) -2*rdotOmega1*s(2,2))
               f_1(3) = temp1*(factor1*runit(3) -2*rdotOmega1*s(3,2))
               f_j(1) = tempj*(factorj*runit(1) -2*rdotOmegaj*s(1,j))
               f_j(2) = tempj*(factorj*runit(2) -2*rdotOmegaj*s(2,j))
               f_j(3) = tempj*(factorj*runit(3) -2*rdotOmegaj*s(3,j))
c    
c     
c     Damping Terms ---------------------------------------------------------- ----------------------------------------------------------   
               temp1=minus_six_n*params(13,1)*massratioj*S1_5*semi_3/r_8
               tempj=minus_six_n*params(13,j)*massratio2*Sj_5*semi_3/r_8    
c     NEED TO DEFINE THE Q-FACTORS ...............
c
c     Adding in the f_TF term from Mardling & Lin 2002, Eqn 4
               f_1(1) = f_1(1) + temp1*(3*rdotv*runit(1) + rcvmrOcr1(1))
               f_1(2) = f_1(2) + temp1*(3*rdotv*runit(2) + rcvmrOcr1(2))
               f_1(3) = f_1(3) + temp1*(3*rdotv*runit(3) + rcvmrOcr1(3))
               f_j(1) = f_j(1) + temp1*(3*rdotv*runit(1) + rcvmrOcrj(1))
               f_j(2) = f_j(2) + temp1*(3*rdotv*runit(2) + rcvmrOcrj(2))
               f_j(3) = f_j(3) + temp1*(3*rdotv*runit(3) + rcvmrOcrj(3))

               a(1,j) = a(1,j) +  f_1(1)
               a(2,j) = a(2,j) +  f_1(2)
               a(3,j) = a(3,j) +  f_1(3)
               a(1,j) = a(1,j) +  f_j(1)
               a(2,j) = a(2,j) +  f_j(2)
               a(3,j) = a(3,j) +  f_j(3)
c
c
c     Evolve the spin...
c           r \cross f (See eqn 15, Mardling & Lin 2002)
               rcf1(1) = xrel(2)*f_1(3) - xrel(3)*f_1(2)
               rcf1(2) = xrel(3)*f_1(1) - xrel(1)*f_1(3)
               rcf1(3) = xrel(1)*f_1(2) - xrel(2)*f_1(1)
               rcfj(1) = xrel(2)*f_j(3) - xrel(3)*f_j(2)
               rcfj(2) = xrel(3)*f_j(1) - xrel(1)*f_j(3)
               rcfj(3) = xrel(1)*f_j(2) - xrel(2)*f_j(1)
           
               h0_mu_over_MoI2 = h0*mu2j / params(14,2) 
               h0_mu_over_MoIj = h0*mu2j / params(14,j) 

               s(1,2) = s(1,2)  - h0_mu_over_MoI2 * rcf1(1)
               s(2,2) = s(2,2)  - h0_mu_over_MoI2 * rcf1(2)
               s(3,2) = s(3,2)  - h0_mu_over_MoI2 * rcf1(3)
               s(1,j) = s(1,j)  - h0_mu_over_MoIj * rcfj(1)
               s(2,j) = s(2,j)  - h0_mu_over_MoIj * rcfj(2)
               s(3,j) = s(3,j)  - h0_mu_over_MoIj * rcfj(3)
c
c         open (unit=1,file='user.out',access='append')
c         write(1,*) 'Debug tidal damping routine (end)... j= ',j,',  ',
c     %        'a(1,j)= ',a(1,j),
c     %        'a(2,j)= ',a(2,j),
c     %        'a(3,j)= ',a(3,j)
c         close(unit=1)
c
            end if
         end do     
      end if
c     
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      DISKPOTN.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne
c 
c Date: October 2010
c
c Want to be able to include the effects of an (external?) disk potential
c
c Use the model from Thommes & Lissauer 2003
c
c------------------------------------------------------------------------------
c
      subroutine diskpotn (time,nbig,m,x,v,a,params,h0)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig
      real*8 time,m(nbig),x(3,nbig),v(3,nbig),a(3,nbig)
      real*8 params(NPARAMS,nbig),h0
c
c Local
      integer j
      real*8 gm,r2D2,r2D,r2D3,z2,z3,rd,ri,ro,ri52,ro52
      real*8 diff52,diff72,diff92,const,fr,fz,semi,r,v2
c
c------------------------------------------------------------------------------
c     
c     DO CALCULATIONS ONLY FOR THE OUTER PLANET...
      do j = nbig, nbig
c
c     Get the r,v & a==semi (note use of Chambers' routine mco_x2a).
         gm = m(1) + m(j)
         call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),
     %        semi,r,v2)
c
c     Calculate the cylindrical positions
         r2D2 = x(1,j)*x(1,j) + x(2,j)*x(2,j)
         r2D  = sqrt(r2D2)
         r2D3 = r2D2*r2D
         z2   = x(3,j)*x(3,j)
         z3   = z2*x(3,j)
c
c     Calculate the inner (ri) and outer (ro) edges of the disk
         rd = 0.1
         ri = rd + 20* semi * (0.333333333333*m(j)/ m(1))**(0.333333333)
         ro = 100
         if (ro.le.ri) ro = 5*ri
c
c     Calculate useful powers of ri & ro 
         ri52 = ri*ri*sqrt(ri)
         ro52 = ro*ro*sqrt(ro)
         diff52 = (1/ri52) - (1/ri52)
         diff72 = (1/(ri52*ri)) - (1/(ro52*ro))
         diff92 = (1/(ri52*ri*ri)) - (1/(ro52*ro*ro))
c
c     Useful constant (\Sigma_0 G \pi )
c     JAN 2011 - PRESUMABLY THIS NEEDS ATO CHANGE!?!?!?!?!?
         const =  3.14159
c
c     Construct the power series...
c     WHY IS EVERYTHING NEGATIVE ?!?!?!?!?
         fr = 0.4*r2D*diff52 - 0.25*r2D3*diff92 + r2D*z2*diff72
         fr = const*fr
         fz = -0.8*x(3,j)*diff52 -r2D2*x(3,j)*diff72 +0.666667*z3*diff92
         fz = const*fz
c
c     Resolve the accelerations
         a(1,j) = a(1,j) - fr*x(1,j)/r2D
         a(2,j) = a(2,j) - fr*x(2,j)/r2D
         a(3,j) = a(3,j) + fz
c
c
      end do      
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      T4NAGASAWATIDAL.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne
c 
c Date: Jan 2011
c
c Want to replicate the models of tidal damping used in Nagasawa 2008 & 20011
c - It seems from their papers that they are generally using 
c        v -> v_new = Factor_v * v      i.e. it's just a scaling of the velocity
c - In some of their models (e.g. 2011, and T3 & T4 in 2008) they then add in 
c   the additional constraint that Delta L =0
c - The easiest way to satisfy this is just say that 
c        x -> x_new = Factor_x * x
c   with 
c        Factor_x = 1 / Factor_v
c
c Until I directly contact Nagasawa I don't know if the above is definitely 
c her approach, but it seems likely.
c
c------------------------------------------------------------------------------
c
      subroutine t4nagasawatidal (time,nbig,m,x,v,params,h0,j)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig,j
      real*8 time,m(nbig),x(3,nbig),v(3,nbig)
      real*8 params(NPARAMS,nbig),h0
c
c Local
      integer i
      real*8 gm,semi,r,v2,h2,p,temp,e,q,E_Anom_Now,M_Anom_Now,t_to_peri
      real*8 M_Anom_Peri,q3_2
      real*8 DeltaE,Scaling_V,hx,hy,hz
c
c------------------------------------------------------------------------------
c     
c     NEED TO GET PERICENTER ...
c      - Using logic similar to that of the mce_cent routine ...
c      - (i) We need to find the eccentricity / pericenter  -->> Use h = r x v
       hx = x(2,j)*v(3,j)-x(3,j)*v(2,j)
       hy = x(3,j)*v(1,j)-x(1,j)*v(3,j)
       hz = x(1,j)*v(2,j)-x(2,j)*v(1,j)
       h2 = hx*hx + hy*hy + hz*hz

       gm = m(1) + m(j)
       p = h2 / gm

       r = sqrt(x(1,j)*x(1,j) +x(2,j)*x(2,j) +x(3,j)*x(3,j))
       v2 = v(1,j)*v(1,j) +v(2,j)*v(2,j) +v(3,j)*v(3,j)
       temp = 1.d0 + p*(v2/gm - 2.d0/r)
       e = sqrt( max(temp,0.d0) )
       q = p / (1.d0 + e)
       semi = q / (1.d0 - e )
c
c     Calculate DeltaE
c     - Use the formulae from Eqn 6, Nagasawa 2011
       q3_2 = q**(1.5)
       DeltaE = (params(15,j) / q3_2 ) * exp(params(16,j) * q3_2)
       DeltaE = (DeltaE + params(17,j) / ( q**9 ))/m(j)
       Scaling_V = sqrt((2*DeltaE/v2) + 1)
c
c
c     Turn OFF the tidal damping if the planet has completely circularized
c     - Use the condition on p10 of Nagasawa 2011
c     N.B. For now I am NOT turning it off permanently, hence this ordering
c          (I want to allow for future scattering / evolution)
c
       if( (gm*gm).gt.(2*DeltaE + v2 -(2.d0*gm/r))*(r*r*v2) ) then 
c
c     Calculate x_new & v_new
c     - N.B. The v-scaling is as described in Nagasawa 2008/2011, however the 
c            x-scaling is my own reading-between-the-lines in order to keep 
c            Delta L = 0
c     - One could presumable implement an asymmetrical scaling to get the T1 
c       and T2 cases where Delta L != 0
          do i = 1,3
             v(i,j) = v(i,j)*Scaling_V
             x(i,j) = x(i,j)/Scaling_V
          end do
c     
c       else
       end if
c
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      T4NAGASAWABINARY.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne
c 
c Date: Jan 2011
c
c Uses the dynamical tide model of Nagasawa, but applies it to 
c the case of a planet interacting with binary stars
c
c------------------------------------------------------------------------------
c
      subroutine t4nagasawabinary (time,nbig,m,x,v,params,h0,j)
c
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c
c Input/Output
      integer nbig,j
      real*8 time,m(nbig),x(3,nbig),v(3,nbig)
      real*8 params(NPARAMS,nbig),h0
c
c Local
      integer i
      real*8 gm,semi,r,v2,h2,p,temp,e,q,E_Anom_Now,M_Anom_Now,t_to_peri
      real*8 M_Anom_Peri,q3_2
      real*8 DeltaE,Scaling_V,hx,hy,hz,dx(3,nbig),dv(3,nbig)
c
c------------------------------------------------------------------------------
c     
c     NEED TO GET PERICENTER ...
c      - Using logic similar to that of the mce_cent routine ...
c      - (i) We need to find the eccentricity / pericenter  -->> Use h = r x v
      do i = 1,3
         dx(i,j) = x(i,j) - x(i,2)
         dv(i,j) = v(i,j) - v(i,2)
      end do
      hx = dx(2,j)*dv(3,j)-dx(3,j)*dv(2,j)
      hy = dx(3,j)*dv(1,j)-dx(1,j)*dv(3,j)
      hz = dx(1,j)*dv(2,j)-dx(2,j)*dv(1,j)
      h2 = hx*hx + hy*hy + hz*hz

      gm = m(2) + m(j)
      p = h2 / gm

      r = sqrt(dx(1,j)*dx(1,j)+dx(2,j)*dx(2,j)+dx(3,j)*dx(3,j))
      v2 = dv(1,j)*dv(1,j)+dv(2,j)*dv(2,j)+dv(3,j)*dv(3,j)
      temp = 1.d0 + p*(v2/gm - 2.d0/r)
      e = sqrt( max(temp,0.d0) )
      q = p / (1.d0 + e)
      semi = q / (1.d0 - e )
c
c     Calculate DeltaE
c     - Use the formulae from Eqn 6, Nagasawa 2011
      q3_2 = q**(1.5)
      DeltaE = (params(15,j) / q3_2 ) * exp(params(16,j) * q3_2)
      DeltaE = (DeltaE + params(17,j) / ( q**9 ))/m(j)
      Scaling_V = sqrt((2*DeltaE/v2) + 1)
c
c
c     Turn OFF the tidal damping if the planet has completely circularized
c     - Use the condition on p10 of Nagasawa 2011
c     N.B. For now I am NOT turning it off permanently, hence this ordering
c          (I want to allow for future scattering / evolution)
c
      if( (gm*gm).gt.(2*DeltaE + v2 -(2.d0*gm/r))*(r*r*v2) ) then 
c
c     Calculate x_new & v_new
c     - N.B. The v-scaling is as described in Nagasawa 2008/2011, however the 
c            x-scaling is my own reading-between-the-lines in order to keep 
c            Delta L = 0
c     - One could presumable implement an asymmetrical scaling to get the T1 
c       and T2 cases where Delta L != 0
         do i = 1,3
            v(i,j) = v(i,2) + dv(i,j)*Scaling_V
            x(i,j) = x(i,2) + dx(i,j)/Scaling_V
         end do
c     
c       else
      end if
c
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      T4NAGASAWABINARY.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matt Payne
c 
c Date: Jan 2011
c
c Use a fitted function to give the additional force which whould be excerted on
c a body due to the presence of an extended mass distribution
c - The method should generalize, but the case considered when developing the
c   routine was that of an extended, time-varying and probably gravitatonally 
c   unstable disk within which planets and planetesimals are orbiting a young 
c   star. How does the additional, extended potential due to the massive gas disk 
c   change the orbital dynamics?

c
c------------------------------------------------------------------------------
c
      subroutine diskgravity (time,nbod,m,x,v,a,params,h0)
c
      use disk_params
      use disk_grav
c      use input
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'

c     
c Input/Output
      integer nbod
      real*8 time,m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
      real*8 params(NPARAMS,nbod),h0
c
c Local
      integer i,j
      real*8 rcyl,rsph,p,cosp,sinp,t,sint,cost,fxp,fyp,fzp,gm,potnpert
      real*8 mjp_x,mjp_y,mjp_z
c
c------------------------------------------------------------------------------
c
c We now use the defined functions/look-up table in conjunction with the particle 
c positions to calculate the additional force/accn on each particle.
c
c      open (unit=1,file='user.out',access='append')
c      write(1,*) 'Here about to do disk-grav force...'
c      close(unit=1)
c
      do j = 2, nbod
c
c Quickly calculate the various radii & angles required by the force routine
         rcyl = x(1,j)*x(1,j) + x(2,j)*x(2,j)
         rsph = rcyl + x(3,j)*x(3,j)
         rcyl = sqrt(rcyl)
         rsph = sqrt(rsph)
         p= atan2(x(2,j),x(1,j)); cosp= x(1,j)/rcyl; sinp= x(2,j)/rcyl
         t = datan2(x(3,j),rcyl);  sint= x(3,j)/rsph; cost= rcyl/rsph   ! Defining t=0 at the z=0 plane
         
c Call the force routine contained in the external diskgrav.f90 file
         call force_all(rcyl,rsph,x(3,j),t,sint,cost,
     %        p,sinp,cosp,time,
     %        fxp,fyp,fzp,potnpert)                         

c Update the forces using the additional quantity calculated above
c ---- The factors of K2 are to correct the units
c ---- In Aaron's code, G=Msun=1  ==> a_Earth = 1*1/1^2 = 1, v=1~30km/sec 
         gm = m(1) + m(j)
         a(1,j) = a(1,j) + gm*fxp
         a(2,j) = a(2,j) + gm*fyp
         a(3,j) = a(3,j) + gm*fzp

         !!!!! NEED TO CHECK THAT THESE UNITS ARE OK !!!!!
         params(18,j) = potnpert * gm
!
!
c
      end do
c     
c------------------------------------------------------------------------------
c     
      return
      end
c
c
