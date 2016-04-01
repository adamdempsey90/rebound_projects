c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      USER_INITIALISATION.F
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Matthew John Payne
c
c Set-up the params(i,nbod) array which used to record various parameters 
c     associated with each of the particles
c
c Use of params within the code
c     params(1) -->> Records the initial semi-major axes
c     params(2) -->> Forced migration / damping parameter for semi-major axes
c                    (input as \dot{a}/{a} = XXX / year, internally converted to / day)
c     params(3) -->> Forced migration / damping parameter for eccentricity
c                    (YYY x params(2), i.e. K>1 means it damps quicker)
c     params(4) -->> Minimum allowed semi (where to stop forced migration -- useful if I do a linear model)
c     params(5) -->> Undefined
c     params(6) -->> Associated Mass for gas damping
c     params(7) -->> Calculated radius of body
c                   (Doesn't actually require input in here, but needs density input elsewhere)
c     params(8) -->> Density x Radius
c     params(9) -->> Mass growth rate
c                    (input as \dot{M}/{M} = XXX / year, internally converted to / day)
c     params(10) -->> Absolute upper mass at which to stop growth (units of Msun)
c     params(11) -->> Relativistic Forces (On / Off for individual bodies)
c     params(12) -->> Tidal Love Number
c     params(13) -->> Tidal Q-factor
c     params(14) -->> Moment of Inertia Factor
c

c------------------------------------------------------------------------------
c
c Start the subroutine
      subroutine user_initialisation (nbod,nbig,m,x,v,rho,params,
     %                                initial,s,rphys)
c Parameter input.
      use disk_params
      use disk_grav
c      use input
      implicit none
      include 'mercury.inc' 
      include 'user_parameters.inc'
c Input/Output
      integer nbod,nbig,initial(4)
      real*8 m(nbod),x(3,nbod),v(3,nbod),rho(nbod)
      real*8 params(NPARAMS,nbod),s(3,nbod),rphys(nbod)   
c Local
      integer i
      real*8 gm,semi,r,v2,mag2, rhocgs, RJ, E, omega0, Q0, temp
c
c
c------------------------------------------------------------------------------
c
c Set some of the basic parameters for the central star
c Don't think *K2 is required, as am CALCULATING it using m(1)->m(1)*K2
      mag2=TIDALSTARSPINX**2 +TIDALSTARSPINY**2 +TIDALSTARSPINZ**2
      s(1,1) = TIDALSTARSPINX * TIDALSTARSPINMAG / sqrt(mag2)
      s(2,1) = TIDALSTARSPINY * TIDALSTARSPINMAG / sqrt(mag2)
      s(3,1) = TIDALSTARSPINZ * TIDALSTARSPINMAG / sqrt(mag2)

      params(1,1) = 0.0
      params(2,1) = 0.0
      params(3,1) = 0.0
      params(4,1) = 0.0
      params(5,1) = 0.0
      params(6,1) = m(1)
  
      params(7,1) = 3*params(6,1)
      params(7,1) = params(7,1)/(4*PI*rho(1))
      params(7,1) = params(7,1)**(1.0/3.0)

      params(8,1) = params(7,1) * rho(1)

      params(9,1) = 0.0
      params(10,1) = 0.0
      params(11,1) = 0.0

      if (TIDALSWITCH.eq.1) then 
         params(12,1) = TIDALSTARLOVE
         params(13,1) = TIDALSTARQ
         params(14,1) = TIDALSTARMOI
         params(14,1) = params(14,1)*(m(1)/K2)*params(7,1)*params(7,1)
      else
         params(12,1) = 0.0
         params(13,1) = 0.0
         params(14,1) = 0.0
      end if

      params(15,1) = 0.0
      params(16,1) = 0.0
      params(17,1) = 0.0

      params(18,1) = 0.0

c
c
c
c Loop over all of the non-central bodies...
      do i = 2, nbod
c
c (1) Record the initial semi-major axes & initial number of bodies
         gm = m(1) + m(i)
         call mco_x2a (gm,x(1,i),x(2,i),x(3,i),
     %        v(1,i),v(2,i),v(3,i),semi,r,v2)
         params(1,i) = semi
         
         initial(1) = nbig
         initial(2) = nbod
         initial(3) = 0
         initial(4) = 0
c
c (2) Forced migration / damping parameter (a)  
c     -->> (Convert AU/year to AU/day)
         params(2,i) = params(2,i) / (365.25d0)
c
c (3) Forced migration / damping parameter (e)  
c     -->> (K* rate of damping to the semi-major axis, i.e. K>1 => damps quicker than a)
         params(3,i) = params(3,i)*params(2,i)
c
c (4) Forced migration / damping minimum semi-major axis (AU)  
         params(4,i) = params(4,i)
c
c (5) -->> Undefined
         params(5,i) = params(5,i) 
c
c (6) Associated Mass for gas damping (only really matters for test particles)
         params(6,i) = params(6,i) * K2
c     
c (7) Calculated radius of body (if Ass Mass=0 => r=0) 
         if (params(6,i).eq.0.d0) then
            params(7,i) = 0d0
         else
            params(7,i) = 3*params(6,i)
            params(7,i) = params(7,i)/(4*PI*rho(i))
            params(7,i) = params(7,i)**(1.0/3.0)
c (8) Useful parameter is often density * radius
            params(8,i) = params(7,i) * rho(i)
            open (unit=1,file='user.out',access='append')
            write(1,*) 'Initialisation, params(6,i)= ',params(6,i),
     %           'params(7,i)= ',params(7,i),
     %           'params(8,i)= ',params(8,i)
            close(unit=1)
         end if
c (9) Mass growth rate
c     -->> Convert \dot{M}/{M} = XXX / year --->>>  XXX / day
         params(9,i) = params(9,i) / (365.25d0)
c (10) -->> Absolute upper mass
         params(10,i) = params(10,i) * K2
c (11) Relativistic Forces (Individual switch for the bodies)
         params(11,i) = params(11,i) 
c (12) Tidal damping parameter: Love Number
         params(12,i) = params(12,i) 
c (13) Tidal damping parameter: Q-factor, convert to love/Q, as that is what is used in routine.
         if (params(13,i).ne.0) then
            params(13,i) = params(12,i) / params(13,i)
         end if
c (14) Moment of Inertia parameter, \alpha. Convert to MoI
         params(14,i) = params(14,i)*m(i)*params(7,i)*params(7,i)
c
         open (unit=1,file='user.out',access='append')
         write(1,*) 'params(11,i)= ',params(11,i),
     %        'params(12,i)= ',params(12,i),
     %        'params(13,i)= ',params(13,i),
c     %        'Moment of inertia, params(14,i)= ',params(14,i),
c     %        'Need to divide by K2, to get Modot au^2 / day',
     %        'm(i) = ',m(i),
     %        'K2 = ',K2
c     %        'rho(i) = ',rho(i)
         close(unit=1)
c
c (15-17) Nagasawa 2008 / 2011 Tidal Damping Model(s)
c     Note that the energy kick can be expressed as ...
c     Delta E = ((params(15) / q^(3/2)) * exp(-params(2)* q^(3/2) ) + params(17) / (q**9))/m(i)
c     where q = a(1-e) = pericenter
c     And then v = sqrt(2 Delta E + v**2) * v / (/v/)
         rhocgs = AU * AU * AU * K2 / MSUN
         RJ = (3*0.001*K2 / (4*PI*(1.326 * rhocgs)))**(1.0/3.0)     ! Jupiter has a density of 1.326 g / cm^3
         E = m(i)*m(i) / params(7,i)              
         omega0 =  0.53*(params(7,i) / RJ) + 0.68
         Q0     = -0.12*(params(7,i) / RJ) + 0.68
         temp = sqrt((m(i))/(m(1) * (params(7,i)**3.0) ))
         params(15,i)=((-PI*PI)/(5*2.0**0.5))*omega0*Q0*Q0*E/temp
         params(16,i)=(-4*sqrt(2.0) / 3)*omega0*temp
         params(17,i)=-0.0036 * E / (temp**6)
         
c
c         open (unit=1,file='user.out',access='append')
c         write(1,*) 'params(15,i)= ',params(15,i),
c     %        'params(16,i)= ',params(16,i),
c     %        'params(17,i)= ',params(17,i),
c     %        ' r= params(7,i)= ',params(7,i),
c     %        ' Qo=',Q0,
c     %        ' omega0 = ',omega0,
c     %        ' temp=',temp,
c     %        ' RJ=',RJ,
c     %        ' rhocgs=',rhocgs,
c     %        ' E = ',E
c         close(unit=1)
c
c
      end do
c
c     
c
c
c
c     Do the initialization for the disk-gravity routine provided by Aaron
       params(18,i)=0.0;
c      call getarg(1,namelist_file)
c      call read_params(namelist_file)
c      call initialize_gravity()
c      call make_disk_model()
c      call get_spherical_moments()
c
c------------------------------------------------------------------------------
c     
      return
      end
c
c
c
c
