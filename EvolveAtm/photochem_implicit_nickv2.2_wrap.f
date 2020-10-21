      program test
        parameter(Nspecies=10,Ntime=500000)
        real*8 :: N_in(Nspecies)
        real*8 :: fH2O

C        real*8, allocatable :: Time(:)
C        real*8, allocatable :: T_surf(:)
C        real*8, allocatable :: N(:,:)
C        real*8, allocatable :: dNdt(:,:)
C        real*8, allocatable :: P(:,:)

        real*8 :: Time(Ntime)
        real*8 :: T_surf(Ntime)
        real*8 :: pressure(Ntime)
        real*8 :: mu(Ntime)
        real*8 :: N_tot(Ntime)
        real*8 :: f_tot(Ntime,Nspecies)
        real*8 :: dNdt(Ntime,Nspecies)


C        input
        fH2O = 1.e-6
        N_in(1) = 1.68675180e+28
        N_in(2) = 1.53771633e+26
        N_in(3) = 9.21714507e+22
        N_in(4) = 6.30403649e+25
        N_in(5) = 6.55363284e+24
        N_in(6) = 2.18545381e+25
        N_in(7) = 9.36588326e+22
        N_in(8) = 1.00000000e+10
        N_in(9) = 1.00000000e+10
        N_in(10) = 1.00000000e+10

C        call function
        call photochem(N_in,fH2O,Time,T_surf,pressure,mu,
     1                 N_tot,f_tot,dNdt)

C        print *,N

      end program



      subroutine photochem(N_in,fH2O,Time_out,T_surf,
     1                    pressure,mu,N_tot,f_tot,dNdt)
C     inputs:
C     N_in = initial column abundance of atmosphere (molecules/cm2)
C     fH2O = H2O mixing ratio in statosphere
C
C     outputs:
C     Time(time) = time (s)
C     T_surf(time) = surface temperature vs time (K)
C     N(time,species) = column abundance vs time (molecules/cm2),
C     dNdt(time,species) = production rate vs time (molecules/cm2/s)
C     P(time,species) = partial pressure vs time (bar)

c  this version has an implicit time integration which works much better
c
c idea is to evolve the atmosphere photochemically
c  this version simplifies by treating stratospheric
c  water vapor mixing ratio as a free parameter fH2O

c Created by Kevin Zahnle on 28 Nov 2019
c Copyright 2019 NASA.


c  these need to be dry atmospheres!

      parameter(Ntime=500000,Nspecies=10,Nsp=7,Nsun=7)
      implicit real*8(A-H,O-Z)
      real*8 p_tot, g, kB, mH  ! total input pressure, gravity etc
      real*8 Area, SSX   ! of Earth cm^2
      real*8 mass(Nspecies)   ! molar mass [g], Columns moles cm^-2
      real*8 NCO, NCO2, NCH4, NH2, NH2O, NNH3, NN2 ! inputs [bars] .ne. partial pressures
c      real*8 p_i(Nsp)  ! inputs actual partial pressures, atmospheres

      real*8 g_old(Nspecies)
      real*8 N_old(Nspecies)
      real*8 N_new(Nspecies)
      real*8 N_difference(Nspecies)
      real*8 N_avo  ! avogadros number

c we store this stuff
      real*8 P(Ntime,Nspecies)
      real*8 N(Ntime,Nspecies)
      real*8 dt(Ntime), kyr      ! 3.15e10 seconds  is unit of time
      real*8 N_j(Nspecies)!, N_tot(Ntime)       ! output columns in Moles/cm2
C      real*8 f_tot(Ntime,Nsp)     ! total mixing ratio
      real*8 N_strat(Ntime)       ! stratospheric column
      real*8 f_strat(Ntime,Nspecies)   ! stratosphere mixing ratio

C      input
      real*8, intent(in) :: fH2O
      real*8, intent(in) :: N_in(Nspecies) ! molecules/cm2

      real*8, intent(out) :: Time_out(Ntime)
      real*8, intent(out) :: T_surf(Ntime)
      real*8, intent(out) :: pressure(Ntime)
      real*8, intent(out) :: mu(Ntime)
      real*8, intent(out) :: N_tot(Ntime)
      real*8, intent(out) :: f_tot(Ntime,Nspecies)
      real*8, intent(out) :: dNdt(Ntime,Nspecies)


      real*8 :: Time(Ntime)

c time rate of change budget variables
      real*8 dCH4dt_ox(Ntime)
      real*8 total_C(Ntime), total_N(Ntime)
      real*8 total_H2loss(Ntime)
      real*8 N_i(Nspecies), dN_idt(Nspecies)

c stratospheric water vapor and surface temperature variables
      real*8 N_H2O_strat, NH2Ox, Nx
      real*8 T_trop,p_trop,p_surf,p_H2O_trop,p_H2O_surf
      real*8 p_H2O_tr, fH2O_strat(Ntime)
      real*8 p_dry(Ntime), N_dry(Ntime)
      real*8 T_s, mu_dry(Ntime)
      real*8 p_s,p_H2O_s, p_H2O(Ntime), N_H2O(Ntime)
      real*8 Mi

c photolysis variables
      real*8 tau_uv(Ntime), tau_vis(Ntime)
      real*8 tau_uv_i, Wolf_Toon
      real*8 tau(Ntime,Nsun), tau_i(Nsun)   ! the "optical depths"
      real*8 sigma(Nsun,Nsp), Flux(Nsun), Sun(Nsun)
      real*8 Phi(Ntime,Nsp), Phi_ion(Ntime,Nsp)
      real*8 Phi_i(Nsp), Phi_ion_i(Nsp)

c H escape variables
      real*8 N_t  ! mean molar mass in atm
      real*8 A_escape, B_escape, BBB, mubar, N_bar, C_escape
      real*8 H2_esc(Ntime), Cum_H2_esc(Ntime)

c chemistry variables
      real*8 T_strat
      real*8 N_HCN(Ntime), N_HCN_i
      real*8 NH3_HCN(Ntime), NH3_HCN_i  !
      real*8 O1D_OH(Ntime), O1D_H2(Ntime), O1D_CH4(Ntime)
      real*8 OH_CO(Ntime), OH_H2(Ntime), OH_CH4(Ntime)
      real*8 O_CO(Ntime), O_H2(Ntime), O_CH4(Ntime)
      real*8 O1D_OH_i, O1D_H2_i, O1D_CH4_i
      real*8 OH_CO_i, OH_H2_i, OH_CH4_i
      real*8 O_CO_i, O_H2_i, O_CH4_i
      real*8 k(8,8)
      real*8 age   ! age of Earth Myrs
      real*8 Mplanet, Rplanet  ! Mass, radius multiples of Earth
      real*8 Aplanet   ! EUV irradiation multiples of Earth

c  geology variables
      real*8 Phi_geo(Ntime,Nsp), Phi_geo_i(Nsp)

      real*8 Nsave(Nspecies), Ndotsave(Nspecies)
      real*8 Jac(Nspecies,Nspecies)
      real*8 KK(Nspecies,Nspecies), eps
      integer ipvt(Nspecies)

      integer buffer   ! flag.  buffer=0 if no buffer
                       !        buffer=10 QFM
                       !        buffer=9  QFM-1
                       !        buffer=20 IW
                       !        buffer=30 QFI
                       !        buffer = -1  does not use IW input

      CHARACTER*8 ISPEC(Nspecies), JSPEC(Nspecies), pspec(Nsp)
      CHARACTER*8 pispec(Nsp), gspec(Nsp), phispec(Nsp)
      character*12 planet

      Common /bblok/ mass, g, Area, mH, kB, kyr
      Common /Lblok/LH2,LCH4,LH2O,LCO2,LCO,LN2,LNH3,
     $   LHCN, LC2Hn, Lhaze
      Common /kblok/k, others3, others4  ! used in chemistry
      Common /Phblok/ sigma, Flux, Sun, SSX
      Common /dblok/ A_escape, B_escape  ! used for H escape

c     Data g/982./  ! cm/s^2
      Data mH/1.66e-24/
      Data kB/1.38e-16/
      Data mass/18.,2.,28.,44.,16.,28.,17.,27.,26.,99./  ! per mole
c     Data Area/5.1e18/  ! cm^2
      Data kyr/3.15e10/  ! 1000 years in seconds
      Data N_Avo/6.02e23/
      Data T_strat/300./
      data Wolf_Toon/3.e10/   ! haze production rate C atoms/cm2/s

      data T_trop, p_trop/180,0.1/  ! tropopause conditions
c     data fH2O/1.e-6/  ! free parameter

c some labels
      DATA LH2O,LH2,LCO,LCO2,LCH4,LN2,LNH3, LHCN, LC2Hn,
     $  Lhaze  /1,2,3,4,5,6,7,8,9,10/

      data eps /1.e-3/  ! used in comouting Jacobian

C      open(38, file='photochem.input',status='UNKNOWN')
C      open(39, file='IW.parameters',status='UNKNOWN')
C      open(76, file='Evolve.p',status='UNKNOWN')  ! output
C      open(77, file='Evolve.out',status='UNKNOWN')  ! output
C      open(78, file='Evolve.mix',status='UNKNOWN')  ! output

c   others  ?  dicey - the logic for CO + O is weak
c   what happens to the others?  makes water!  Good.  We don't track water.
      others3 = 1e6*1e4   !  scale height * density
      others4 = 1e6*1e4   !  scale height * density

c the actual skin temperature in runaway is
c     T_trp = (2.8e5/2/5.67e-5)**0.25 ! = 223 K
c     T_trop = 223
c     print *, T_trp
c

c ISPEC is for printing output
C  Define names of species
      ISPEC(LH2)   =   'H2'
      ISPEC(LCH4)  =   'CH4'
      ISPEC(LH2O)  =   'H2O'
      ISPEC(LCO2)  =   'CO2'
      ISPEC(LCO)   =   'CO'
      ISPEC(LN2)   =   'N2'
      ISPEC(LNH3)  =   'NH3'
      ISPEC(LHCN)  =   'HCN'
      ISPEC(LC2Hn) =   'C2Hn'
      ISPEC(Lhaze) =   'haze'

C  Define names of rates
      JSPEC(LH2)   =   'dH2dt'
      JSPEC(LCH4)  =   'dCH4dt'
      jSPEC(LH2O)  =   'dH2Odt'
      jspec(LCO2)  =   'dCO2dt'
      jspec(LCO)   =   'dCOdt'
      jspec(LN2)   =   'dN2dt'
      jspec(LNH3)  =   'dNH3dt'
      jspec(LHCN)  =   'dHCNdt'
      jspec(LC2Hn) =   'dC2Hndt'
      jspec(Lhaze) =   'dhazedt'

c  pressure
      pspec(LH2)   =   'pH2'
      pspec(LCH4)  =   'pCH4'
      pspec(LH2O)  =   'pH2O'
      pspec(LCO2)  =   'pCO2'
      pspec(LCO)   =   'pCO'
      pspec(LN2)   =   'pN2'
      pspec(LNH3)  =   'pNH3'

c  photolysis
      phispec(LH2)   =   'PhiH2'
      phispec(LCH4)  =   'PhiCH4'
      phispec(LH2O)  =   'PhiH2O'
      phispec(LCO2)  =   'PhiCO2'
      phispec(LCO)   =   'PhiCO'
      phispec(LN2)   =   'PhiN2'
      phispec(LNH3)  =   'PhiNH3'
c ions
      pispec(LH2)   =   'PH2i'
      pispec(LCH4)  =   'PCH4i'
      pispec(LH2O)  =   'PH2Oi'
      pispec(LCO2)  =   'PCO2i'
      pispec(LCO)   =   'PCOi'
      pispec(LN2)   =   'PN2i'
      pispec(LNH3)  =   'PNH3i'

c geological fluxes
      gspec(LH2)   =   'geoH2'
      gspec(LCH4)  =   'geoCH4'
      gspec(LH2O)  =   'geoH2O'
      gspec(LCO2)  =   'geoCO2'
      gspec(LCO)   =   'geoCO'
      gspec(LN2)   =   'geoN2'
      gspec(LNH3)  =   'geoNH3'

c assign N.  These are Moles/cm2.  This is where to start.

C      read (38,*) fH2O, buffer, planet, Mi
      planet = 'Earth'
C      Mi not used
C      buffer not used
C      fix planet above to Earth

C      read (39,*) NCO, NCO2, NH2, NH2O, NCH4, NN2, NNH3
C      print *, fH2O, NCO, NCO2, NH2, NH2O, NCH4, NN2, NNH3

      if (planet.eq.'Earth') then   ! fixed parameters
        area = 5.1e18
        g = 982.
        Aplanet = 1.    ! these can be read in through unit 38
        Mplanet = 1.
        Age     = 300.
        SSX = 1.
        B_escape = 0.006
        x_escape = 1.8*SSX
        C_escape = 2*exp(-1/x_escape)/(1+exp(-1/x_escape))
        A_escape = 2.e12*C_escape
      else
        stop 666   ! terra incognita
      endif

c  the program will work in molecules/cm2

c assign these too
C      N_in(LCO)  = NCO
C      N_in(LCO2) = NCO2
C      N_in(LH2)  = NH2
C      N_in(LCH4) = NCH4
C      N_in(LH2O) = NH2O   ! this will be revised?  No - leave it be
C      N_in(LN2)  = NN2
C      N_in(LNH3) = NNH3

c  the NH2O read in is *NOT* to be used, save for tests of whether all the water is evaporated

c assign these as initial columns
      do j =1,Nsp
        N(1,j)  = N_in(j)*1   ! convert to molecules
      enddo
      do j=Nsp+1,Nspecies
        N(1,j)  = 1.e10   ! arbitrary.  for comouting Jacobian, must not have zeroes in these
      enddo
      total_C(1) = N(1,LCO) + N(1,LCO2) + N(1,LCH4)    ! C conservation check
      total_N(1) = 2*N(1,LN2) + N(1,LHCN) + N(1,LNH3)  ! N conservation check
      total_H2loss(1) = 0.0

c  TTTTTTTTTTTTTT  tropopause   TTTTTTTTTTTTT   stratosphere  TTTTTTTTTTTTTT
c  its enough of a mess that it deserves a subroutine

c obtain dry atmosphere parameters
      mu_dry(1) = - mass(LH2O)*N(1,LH2O)
      N_dry(1) = - N(1,LH2O)    ! this will be in molecules/cm2
      do j=1,Nsp
          mu_dry(1) = mu_dry(1) + mass(j)*N(1,j)
          N_dry(1) = N_dry(1) + N(1,j)  ! the dry atmosphere
      enddo
      mu_dry(1) = mu_dry(1)/N_dry(1)
C      print *,'dry', N_dry(1), mu_dry(1)  ! N_dry is in molecules


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c the following segment is used to estimate the amount of water in
c  the atmosphere, needed for determining pressures.
c  it also determines how much of the atmosphere is above the troposphere

c  when the dry atmosphere is very big and almost entirely H2, it doesn't work with eta=0.14
c  to get it to work requires a smaller eta.
c      gamma = 1/(1-eta)
c      eta   = (gamma-1)/gamma
c  eta=0.14 is what gives your Earth.
c     gamma(eta=0.14) = 1.16
c     eta(gamma=1.15) = 0.13
c  My guess is that eta<0.14 for a wetter atmosphere
c  however, I can't make eta a function of water vapor because that would make the equations ridiculously complicated

c first guess
      p_surf  = mu_dry(1)*mH*g*N_dry(1)/1.e6     ! the dry weight

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  the little fix below does the job for the parameter space I'm in
c  the power eta should get smaller for wetter atmospheres,
c  but to be honest all of these eta models are pretty bogus
c     eta = 0.14*(1/p_surf)**0.03  ! the little fix
       eta = 0.14
C      print *, 'eta', eta
c   Roxana tells me that no one has been doing these moist adiabats all that well


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      T_s  = T_trop*(p_surf/p_trop)**(eta)      ! moist adiabat
      p_H2O_surf  = exp(5000./373.-5000./T_s )   ! water vapor pressure, bars
C      print *, 'surf', p_surf, T_s , p_H2O_surf

c  this is a newton solution for self-consistent total pressure p with eta
c  it includes an eta adiabat and the water vapor pressure
c   it treats all the dry species together as a single component
c  key is to put all the units into cgs to avoid stupids

c the following is unnecessary and untrustworthy and built in
c I don't see why I should keep it
c I need to determine if there is always a solution.
c  if there isn't, then what?
c  a simpler troposphere?
      p_s = 1.e6*(p_surf+p_H2O_surf)  ! dynes/cm2
      p_H2O_s = p_H2O_surf*1.e6       ! dynes/cm2
      do iloop=1,7
        FF = ((p_s-p_H2O_s)*mu_dry(1)*mH + p_H2O_s*mass(LH2O)*mH)
     $   * (N_dry(1)*(mass(LH2O)-mu_dry(1))*mH*g + p_s)
     $    - p_s**2*mass(LH2O)*mH
        dFFdp = (mu_dry(1)*mH + (mass(LH2O)-mu_dry(1))*mH
     $    *p_H2O_s*5000/T_trop*eta/p_s*(1e6*p_trop/p_s)**eta)
     $    * (N_dry(1)*(mass(LH2O)-mu_dry(1))*mH*g + p_s)
     $    + (p_s-p_H2O_s)*mu_dry(1)*mH + p_H2O_s*mass(LH2O)*mH
     $    - 2*p_s*mass(LH2O)*mH
C        print *, iloop, FF, p_s/1e6, T_s
        p_s = p_s - FF/dFFdp
        T_s = T_trop*(p_s/1e6/p_trop)**(eta)
        p_H2O_s = 1.e6*exp(5000./373.-5000./T_s)
C        print *, iloop, ff, p_s, T_s
      enddo

      T_surf(1) = T_s
      pressure(1) = p_s/1.e6
      p_H2O(1) = p_H2O_s/1.e6  ! convert to bars
      p_dry(1) = pressure(1) - p_H2O(1) ! partial pressure of dry gases
      mu(1) = (p_dry(1)*mu_dry(1) + p_H2O(1)*mass(LH2O))/pressure(1)
      N_H2O(1) = pressure(1)*1e6/mH/mu(1)/g - N_dry(1)   ! wet atmosphere

C      print *, 'T_s', T_s, pressure(1), p_H2O(1), p_dry(1)

c this statement changes for fH2O fixed
c determine N_H2O_trop for stratospheric mixing ratio
      p_H2O_tr    = exp(5000./373.-5000./T_trop)   ! bars
      fH2O_strat(1)  = p_H2O_tr/p_trop
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c use fixed fH2O
      p_H2O_trop  = fH2O*p_trop  ! bars

c the fraction of dry species above z_trop is
      fraction = (p_trop - p_H2O_trop)/p_dry(1)
c the stratospheric mixing ratio is  p_H2O_trop/p_trop
      NX = N_dry(1)*fraction  ! total dry column above tropopause, in molecules
      N_H2O_strat = p_H2O_trop/p_trop*NX   ! stratosphere column of H2O, in molecules
      N_strat(1)  = NX + N_H2O_strat       ! total stratospheric column
c the scaled H2O is (scaling to whole column is for chemistry and photolysis)
      NH2Ox = p_H2O_trop/p_trop*N_dry(1)
      N(1,LH2O)  = NH2Ox    ! this is the scaled N_in(LH2O) at the mixing ratio of stratosphere
                            ! used in photolyis and chemistry
c  TTTTTTTTTTTTTT   tropopause   TTTTTTTTTTTTT   stratosphere   TTTTTTTTTTTTTT

C      print *, p_H2O_trop, fraction, NX,  N_H2O_strat, NH2Ox

c   partial pressures and mixing ratios.
c   N_tot   includes the actual surface pressure of H2O
      N_t = 0.
      N_tot(1) = N_H2O(1) - N(1,LH2O)   ! correct for scaled N(H2O)
      do j=1,NSP
           N_t      = N_t + N(1,j)
           N_tot(1) = N_tot(1) + N(1,j)
      enddo
      do j=1,NSpecies
          p(1,j)       = pressure(1)*N(1,j)/N_tot(1)
      enddo
      p(1,LH2O) = pressure(1)*N_H2O(1)/N_tot(1)
c mixing ratios
      do j=1,NSpecies
          f_tot(1,j)   = p(1,j)/pressure(1)  ! total mixing ratio
          f_strat(1,j) = N(1,j)/N_t          ! stratospheric mixing ratio
      enddo

c     print *, f_strat(1,LH2O), fH2O, f_tot(1,LH2O)

      call k_rates(T_strat)
      call Photo_rates(T_strat,age)

c  geology
        do j=1,Nsp  ! initialize variables
          Phi_geo(1,j) = 0.0
          Phi_geo_i(j) = 0.0
        enddo

c  tau_uv.  Getting this wrong at the first step leads to some
c   ringing at the beginning.  So I've decided that the first
c   few time steps will be very short
      tau_uv(1)  = 100.   ! rayleigh if naught else
      tau_vis(1) = 1.  ! set this to something

c  keep time in kyrs
      time(1) = 0.0
      time_out = 0.0

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      dt(1) = 1.e-6  ! this will be important
      DO it=1,Ntime-1  ! outer loop
        time(it+1) = time(it) + dt(it)
        time_out(it+1) = time(it+1)*kyr
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c set up stuff for next time step
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c  geology
        do j=1,Nsp  ! initialize variables
          Phi_geo(it,j) = 0.0
          Phi_geo_i(j) = 0.0
        enddo

        tau_uv_i = tau_uv(it)
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        do j=1,Nspecies
          N_i(j) = N(it,j)
          N_new(j) = N(it,j)
          N_old(j) = N(it,j)
        enddo

C        print *, (dt(it)*kyr)
        DO ll=1,10
        N_old = N_new

C        print *, ll, N_new


        call Ndot(it,N_old, Phi_geo_i, tau_uv_i, NX,
     $  dN_idt, H2_esc_i, Phi_ion_i, Phi_i,
     $ tau_i, dCH4dt_ox_i, OH_CO_i, OH_H2_i, OH_CH4_i,
     $  O1D_OH_i, O1D_H2_i, O1D_CH4_i, O_CO_i,
     $  O_H2_i, O_CH4_i, N_HCN_i, NH3_HCN_i)  ! this just fills up the reaction rates
c          all the other stuff is just for storing


c save unperturbed values
        do j=1,Nspecies
          Nsave(j) = N_old(j)
          Ndotsave(j) = dN_idt(j)
          dNdt(it,j) = dN_idt(j)    ! this is saved for output,, not used in explicit
          g_old(j) = N_old(j)-N(it,j)-(dt(it)*kyr)*dN_idt(j)
        enddo
c  initialize Jacobian
        do j=1,Nspecies
          do i=1,Nspecies
            Jac(i,j) = 0.0
          enddo
        enddo
        do j =1,Nspecies
          do jj=1,Nspecies
            N_i(jj) = Nsave(jj)
          enddo
          N_i(j) = Nsave(j) *(1.+eps)   ! differential in N_i(j)
c call Ndot to determine dN_idt of N_i
          call Ndot(it,N_i, Phi_geo_i, tau_uv_i, NX,
     $      dN_idt, H2_esc_i, Phi_ion_i, Phi_i,
     $      tau_i, dCH4dt_ox_i, OH_CO_i, OH_H2_i, OH_CH4_i,
     $      O1D_OH_i, O1D_H2_i, O1D_CH4_i, O_CO_i,
     $      O_H2_i, O_CH4_i, N_HCN_i, NH3_HCN_i)  ! this just fills up the reaction rates
          do i=1,Nspecies
            Jac(i,j) = (dN_idt(i) - Ndotsave(i))/(eps*Nsave(j))
          enddo
        enddo   ! end of j loop

c       print *, Jac
c       stop
c fill KK
        do j=1,Nspecies
          do i = 1,Nspecies
            KK(i,j) = -(dt(it)*kyr)*Jac(i,j)
          enddo
        enddo
        do j=1,Nspecies
          KK(j,j) = -(dt(it)*kyr)*Jac(j,j) + 1.
        enddo

        CALL SGEFA(KK,10,10,IPVT,INFO)
        if(info.NE.0) then
          print 303, info
 303      FORMAT(1x,'info =',I3)
        else
          CALL SGESL(KK,10,10,IPVT,g_old,0)
        endif

        do j=1,Nspecies
          N_new(j) = N_old(j) - g_old(j)
        enddo

        do j=1,Nspecies
          N_difference(j) = (N_new(j)-N_old(j))/N_new(j)
        enddo

        if (maxval(N_difference) .le. 1e-10) then
          exit
        endif

        ENDDO
C       end newtons loop



C        exit


C        if (it.eq.100) then
C          exit
C        endif


c implicit integration
        do j=1,Nspecies
          N(it+1,j) = max(N_new(j),1.)
        enddo


c *********** H escape stuff  ***********
c integrate H2 loss
        H2_esc(it) = H2_esc_i
        total_H2loss(it+1) = total_H2loss(it) - H2_esc(it)*dt(it)*kyr
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c  save optical depth
        do L=1,Nsun
          tau(it,L) = tau_i(L)
        enddo

c ion production
        Do j=1,Nsp
          Phi_ion(it,j) = Phi_ion_i(j)
c       print *, ispec(j), Phi_ion(i,j)
        Enddo

c ******************** Total  ********************  Photolysis *********************
        Phi(it,LCO2) = Phi_i(LCO2)
        Phi(it,LCH4) = Phi_i(LCH4)
        Phi(it,LH2O) = Phi_i(LH2O)
        Phi(it,LN2)  = Phi_i(LN2)
        Phi(it,LNH3) = Phi_i(LNH3)
        Phi(it,LH2)  = 0.
        Phi(it,LCO)  = 0.

c   Store *********************** Chemistry  ***********************
c O1D
        O1D_H2(it)  = O1D_H2_i
        O1D_CH4(it) = O1D_CH4_i
        O1D_OH(it)  =  O1D_OH_i
c OH
        OH_CH4(it) = OH_CH4_i
        OH_H2(it)  = OH_H2_i
        OH_CO(it)  = OH_CO_i
c N > HCN
        N_HCN(it) = N_HCN_i
c  O + CO > CO2
        O_CH4(it) = O_CH4_i
        O_H2(it)  = O_H2_i
        O_CO(it)  = O_CO_i
c  NH3 > HCN
        NH3_HCN(it) =  NH3_HCN_i
c  oxidation of CH4 to make...
        dCH4dt_ox(it) = dCH4dt_ox_i


c $$$$$$$$$$$$$   Budgets   $$$$$$$$$$$$$   Budgets   $$$$$$$$$$$$$$$


c keep track of totals
        total_C(it+1) =  N(it+1,LCO) +N(it+1,LCO2) +N(it+1,LCH4)
     $       + N(it+1,Lhaze) + N(it+1,LHCN)
        total_N(it+1) = 2*N(it+1,LN2) + N(it+1,LHCN) + N(it+1,LNH3)
c       do j =1,NSpecies
c       print *, i, ispec(j), N(i,j), N(i+1,j)
c       enddo

c  tau_uv depends on haze production rate
        tau_uv(it+1) = 1. +
     $      10.*((dNdt(it,LHaze)+dNdt(it,LHCN))/Wolf_Toon)**0.8
        tau_vis(it+1)= 0.01 +
     $      0.5*((dNdt(it,LHaze)+dNdt(it,LHCN))/Wolf_Toon)**0.7
c       print *, i, tau_uv(i+1), tau_vis(i+1)

c set up tropopause conditions for next step

c obtain dry atmosphere parameters
        mu_dry(it+1) = - mass(LH2O)*N(it+1,LH2O)
        N_dry(it+1) = - N(it+1,LH2O)    ! this will be in molecules/cm2
        do j=1,Nsp
          mu_dry(it+1) = mu_dry(it+1) + mass(j)*N(it+1,j)
          N_dry(it+1) = N_dry(it+1) + N(it+1,j)  ! the dry atmosphere
        enddo
        mu_dry(it+1) = mu_dry(it+1)/N_dry(it+1)

c first guess use previous T_surf
        T_s = T_surf(it)
        p_H2O_surf  = exp(5000./373.-5000./T_s )   ! water vapor pressure, bars
c       print *, T_s, 'p_H2O_surf',p_H2O_surf

c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
        p_s = 1.e6*pressure(it)          ! dynes/cm2
        p_H2O_s = p_H2O_surf*1.e6       ! dynes/cm2
        do iloop=1,4
          FF = ((p_s-p_H2O_s)*mu_dry(it+1)*mH + p_H2O_s*mass(LH2O)*mH)
     $     * (N_dry(it+1)*(mass(LH2O)-mu_dry(it+1))*mH*g + p_s)
     $      - p_s**2*mass(LH2O)*mH
          dFFdp = (mu_dry(it+1)*mH + (mass(LH2O)-mu_dry(it+1))*mH
     $      *p_H2O_s*5000/T_trop*eta/p_s*(1e6*p_trop/p_s)**eta)
     $      * (N_dry(it+1)*(mass(LH2O)-mu_dry(it+1))*mH*g + p_s)
     $      + (p_s-p_H2O_s)*mu_dry(it+1)*mH + p_H2O_s*mass(LH2O)*mH
     $      - 2*p_s*mass(LH2O)*mH
c         print *, iloop, FF, p_s, T_s
          p_s = p_s - FF/dFFdp
          T_s = T_trop*(p_s/1e6/p_trop)**(eta)
          p_H2O_s = 1.e6*exp(5000./373.-5000./T_s)
c         print *, iloop, ff, p_s, T_s
        enddo

c columns are good.
c pressures are just so weird and counter-intuitive.
        T_surf(it+1) = T_s
        pressure(it+1) = p_s/1.e6
        p_H2O(it+1) = p_H2O_s/1.e6  ! convert to bars
        p_dry(it+1) = pressure(it+1) - p_H2O(it+1)
        mu(it+1) = (p_dry(it+1)*mu_dry(it+1)
     $    + p_H2O(it+1)*mass(LH2O))/pressure(it+1)
        N_H2O(it+1) = pressure(it+1)*1e6/mH/mu(it+1)/g - N_dry(it+1)

c       print *, N_H2O(i+1), N_dry(i+1), mu(i+1)

c determine N_H2O_trop for stratospheric mixing ratio

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c this statement changes for fH2O fixed
c determine N_H2O_trop for stratospheric mixing ratio
        p_H2O_tr         = exp(5000./373.-5000./T_trop)   ! bars
c store estimated fH2O
        fH2O_strat(it+1)  = p_H2O_tr/p_trop
c   but use fH2O fixed
        p_H2O_trop       = fH2O*p_trop  ! bars
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c the fraction of dry species above z_trop is
        fraction = (p_trop - p_H2O_trop)/p_dry(it+1)
c the stratospheric mixing ratio is  p_H2O_trop/p_trop
        NX = N_dry(it+1)*fraction  ! total dry column above tropopause, in molecules
        N_H2O_strat = p_H2O_trop/p_trop*NX   ! stratosphere column of H2O, in molecules
        N_strat(it+1)  = NX + N_H2O_strat     ! total stratospheric column
c the scaled H2O is (scaling to whole column is for chemistry and photolysis)
        NH2Ox = p_H2O_trop/p_trop*N_dry(it+1)
        N(it+1,LH2O)  = NH2Ox    ! this is the scaled N_i(LH2O) at the mixing ratio of stratosphere
                            ! used in photolyis and chemistry
c  TTTTTTTTTTTTTT   tropopause   TTTTTTTTTTTTT   stratosphere   TTTTTTTTTTTTTT

c   partial pressures.  *Include* surface pressure of H2O
        N_t        = 0.0
        N_tot(it+1) = N_H2O(it+1) - N(it+1,LH2O)
        do j=1,NSP
           N_t        = N_t  + N(it+1,j)
           N_tot(it+1) = N_tot(it+1) + N(it+1,j)
        enddo
        do j=1,Nspecies
          p(it+1,j) = pressure(it+1)*N(it+1,j)/N_tot(it+1)
        enddo
        p(it+1,LH2O) = pressure(it+1)*N_H2O(it+1)/N_tot(it+1)
c mixing ratios
        do j=1,NSpecies
          f_tot(it+1,j)   = p(it+1,j)/pressure(it+1)  ! total mixing ratio
          f_strat(it+1,j) = N(it+1,j)/N_t       ! stratospheric mixing ratio
        enddo

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c time step control
        If (it.gt. 10) Then  ! once the initial ringing has died out
          dtime = 1000.*kyr
          If ( abs(N(it+1,LCO)) .gt. 1.e10) Then
            dtime = min(dtime,abs(N(it+1,LCO)/dNdt(it,LCO)))
          Endif
          If ( abs(N(it+1,LCO2)) .gt. 1.e10) Then
            dtime = min(dtime,abs(N(it+1,LCO2)/dNdt(it,LCO2)))
          Endif
          If ( abs(N(it+1,LH2)) .gt. 1.e10) Then
            dtime = min(dtime,abs(N(it+1,LH2)/dNdt(it,LH2)))
          Endif
          If ( abs(N(it+1,LCH4)) .gt. 1.e10) Then
            dtime = min(dtime,abs(N(it+1,LCH4)/dNdt(it,LCH4)))
          Endif
          If ( abs(N(it+1,LNH3)) .gt. 1.e10) Then
            dtime = min(dtime,abs(N(it+1,LNH3)/dNdt(it,LNH3)))
          Endif
          dt(it+1) = 0.02*dtime/kyr  ! this will be important
        Else  ! start with a few small time steps to get tau_uv settled down
          dt(it+1) = dt(it)
        Endif
C        print *, it, dt(it+1)

c  decide if we are done
c  I will stop the integration when both N(i+1,LH2) and N(i+1,LCH4)
c   are smaller than 1.e10  molecules/cm2/s
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
        If ((f_tot(it,LCH4).lt.1.e-6).AND.(f_tot(it,LH2).lt.1.e-4)) Then
          imax = it
          goto 999
        Else
          imax = Ntime-1
        Endif
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      ENDDO  ! end outer Time loop

 999  Continue
      end
C      end main subroutine

      subroutine k_rates(T)  ! this just fills up the reaction rates
      implicit real*8(A-H,O-Z)
      real*8 k(8,8), T
      Common /kblok/k, others3, others4
      Common /Lblok/LH2,LCH4,LH2O,LCO2,LCO,LN2,LNH3,
     $   LHCN, LC2Hn, Lhaze

      do i=1,8
        do j=1,8
          k(i,j) = 0.0
        enddo
      enddo

c N2D reactions
      k(1,LH2) = 2.2e-12
      k(1,LCH4) = 4e-12
      k(1,LH2O) = 5e-11
      k(1,LCO2) = 3.6e-13
      k(1,LCO) = 1.9e-12
      k(1,LN2) = 1.7e-14

c O1D reactions
      k(2,LH2) = 1.1e-10
      k(2,LCH4) = 1.5e-10
      k(2,LH2O) = 2.2e-10
      k(2,LCO2) = 7.4e-11
      k(2,LCO) = 7e-12
      k(2,LN2) = 1.8e-11

c OH reactions
      k(3,LH2) = 6e-15
      k(3,LCH4) = 6e-15
      k(3,LCO) = 1.2e-13
      k(3,6) = 1.e-11    ! OH + others3
      k(3,7) = 1.2e-10   ! OH + CH2
      k(3,8) = 6e-11     ! OH + CH3

c O reactions
      k(4,LH2) = 1e-17
      k(4,LCH4) = 7e-18
      k(4,LCO) = 4e-17
      k(4,6)   = 1.e-11  !  O + others4
      k(4,7) = 1.2e-10   ! O + CH2
      k(4,8) = 1.2e-10   ! O + CH3

c N reactions
      k(5,7) = 1.2e-10   ! N + CH2
      k(5,8) = 1.1e-10

      k(7,8) = 7e-11     ! CHn + CHn

c H2 ion reactions
      k(8,LH2)   =  2.e-9
      k(8,LCH4)  =  3.8e-9
      k(8,LH2O)  =  7.e-9
      k(8,LCO2)  =  2.4e-9
      k(8,LCO)   =  3.e-9
      k(8,LN2)   =  2.e-9

      return
      end

      subroutine photo_rates(T,Age)  ! this just fills up the photolysis
      parameter(Nsun=7,Nsp=7)
      implicit real*8(A-H,O-Z)
      real*8 Sun(Nsun), Flux0(Nsun), SSX
      real*8 Sunpower(Nsun)
      real*8 Age, Sun_age    ! Age of Earth, Sun, Myrs
      real*8 sigma(Nsun,Nsp), Flux(Nsun)
      real*8 Lyalpha, Lygamma, T
      Common /Lblok/LH2,LCH4,LH2O,LCO2,LCO,LN2,LNH3,
     $   LHCN, LC2Hn, Lhaze
      Common /Phblok/sigma, Flux, Sun, SSX

c     data Sun/30,20,10,10,10, 4,2/
      data Sun_age/4567/     ! Myrs
c     data Sunpower/1.25,1.1,0.85,0.85,0.85,0.6,0.4/   ! probably better
      data Sunpower/1.25,1.1,0.85,0.85,0.85,0.5,0.3/   ! reproduces Sun as used in paper with Age=300

      do i=1,Nsun
        Sun(i) = (Sun_age/Age)**sunpower(i)
      enddo

      Lyalpha = 3.6e11
      Lygamma = 3.6e9
c this is unique to Sun
      Flux0(1) = 3e10     ! FXUV
      Flux0(2) = 2e10     ! FEUV
      Flux0(3) = 3.6e9    ! Lyman gamma +
      Flux0(4) = 3.6e11   ! Lyman alpha
      Flux0(5) = 4e11     ! FUV_CO2
      Flux0(6) = 2.5e12   ! FUV_H2O
      Flux0(7) = 6e13     ! FUV_NH3

c scale up to young Sun
      do i=1,Nsun
        Flux(i) = Flux0(i)*Sun(i)*SSX
      enddo

c cross sections
      do i=1,Nsun
        do j=1,Nsp
          sigma(i,j) = 0.0
        enddo
      enddo
      sigma(1,LH2)  = 2.e-18
      sigma(1,LCH4) = 2.e-17
      sigma(1,LH2O) = 1.3e-17
      sigma(1,LCO2) = 2.e-17
      sigma(1,LCO)  = 1.3e-17
      sigma(1,LN2)  = 1.3e-17

      sigma(2,LH2)  = 1.e-19
      sigma(2,LCH4) = 8.e-18
      sigma(2,LH2O) = 8.e-18
      sigma(2,LCO2) = 4.e-17
      sigma(2,LCO)  = 1.e-19
      sigma(2,LN2)  = 1.e-19

      sigma(3,LH2)   = 1.e-19
      sigma(3,LCH4)  = 8.e-18
      sigma(3,LH2O)  = 8.e-18
      sigma(3,LCO2)  = 4.e-17
      sigma(3,LCO)   = 1.e-19
      sigma(3,LN2)   = 3.e-16

c Lyman alpha
      sigma(4,LCH4)  = 8.e-18
      sigma(4,LH2O)  = 8.e-18
      sigma(4,LCO2)  = 5.e-20
      sigma(4,LNH3)  = 8.e-18   ! made up

      sigma(5,LCH4)  = 0.e-18   ! missing
      sigma(5,LH2O)  = 8.e-19
      sigma(5,LCO2)  = 8.e-18
      sigma(5,LNH3)  = 8.e-19   ! made up

      sigma(6,LH2O)  = 2.e-18
      sigma(6,LCO2)  = 3.e-20
      sigma(6,LNH3)  = 8.e-19   ! made up

      sigma(7,LNH3)  = 3.e-18

      return
      end

      subroutine sgefa(a,lda,n,ipvt,info)
      implicit real*8(A-H,O-Z)
      integer lda,n,ipvt(1),info
      real*8 a(lda,1)
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      real*8 t
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end

      subroutine sgesl(a,lda,n,ipvt,b,job)
      implicit real*8(A-H,O-Z)
      integer lda,n,ipvt(1),job
      real*8 a(lda,1),b(1)
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c
c     internal variables
c
      real*8 sdot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end


      subroutine sscal(n,sa,sx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit real*8(A-H,O-Z)
      real*8 sa,sx(1)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end

      real*8 function sdot(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit real*8(A-H,O-Z)
      real*8 sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end

      integer function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit real*8(A-H,O-Z)
      real*8 sx(1),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end

      subroutine saxpy(n,sa,sx,incx,sy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit real*8(A-H,O-Z)
      real*8 sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end

c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      subroutine Ndot(i, N_i, Phi_geo, tau_uv, NX,
     $  dN_idt, H2_esc, Phi_ion, Phi,
     $ tau, dCH4dt_ox, OH_CO, OH_H2, OH_CH4, O1D_OH, O1D_H2, O1D_CH4,
     $ O_CO, O_H2, O_CH4, N_HCN, NH3_HCN)  ! this just fills up the reaction rates

c  input: N_i(Nspecies)
c         Phi_geo(Nsp)
c         tau_uv  !  this is hard to compute implicitly
c         NX
c

c  output: dN_idt(Nspecies)
c          H2_esc
c          Phi_ion(Nsp), Phi(Nsp)
c          tau(Nsun)
c          dCH4dt_ox
c          OH_CO, OH_H2, OH_CH4
c          O1D_OH, O1D_H2, O1D_CH4
c          O_CO, O_H2, O_CH4
c          N_HCN, NH3_HCN


      parameter(Nspecies=10,Nsp=7,Nsun=7)
      implicit real*8(A-H,O-Z)
      real*8 N_i(Nspecies), dN_idt(Nspecies)
      real*8 x(Nsp),F(Nsp)
      real*8 k(8,8)
      real*8 mass(Nspecies)   ! molar mass [g], Columns moles cm^-2
      real*8 mubar, N_t, N_bar, NX
      real*8 A_escape, B_escape, BBB
      real*8 H2_esc
c time rate of change budget variables
      real*8 dCH4dt_ox
c  geology variables
      real*8 Phi_geo(Nsp)
c photolysis variables
      real*8 tau_uv, tau_vis
      real*8 tau(Nsun)   ! the "optical depths"
      real*8 sigma(Nsun,Nsp), Flux(Nsun), Sun(Nsun)
      real*8 Phi(Nsp), Phi_ion(Nsp)
      real*8 sum_ions, CH4_ion_breakup
c chemistry variables
      real*8 sumN2D, N_HCN, NH3_HCN
      real*8 sumO1D, O1D_OH, O1D_H2, O1D_CH4
      real*8 sumOH, OH_CO, OH_H2, OH_CH4
      real*8 sumO, O_CO, O_H2, O_CH4
      real*8 g, area, mH, kB, kyr
      real*8 fH2


      Common /Phblok/ sigma, Flux, Sun, SSX
      Common /kblok/k, others3, others4
      Common /Lblok/LH2,LCH4,LH2O,LCO2,LCO,LN2,LNH3,
     $   LHCN, LC2Hn, Lhaze
      Common /dblok/A_escape, B_escape
      Common /bblok/ mass, g, Area, mH, kB, kyr


c **************  hydrogen   **************  escape  **************
c determine mu of evolving stratosphere for H escape
        mubar = -N_i(LH2)*mass(LH2)  ! mu without H2
        N_t = 0.  ! this will be in molecules/cm2  its only used for mixing ratios
        N_bar = -N_i(LH2)
        do j=1,Nsp
          mubar = mubar + N_i(j)*mass(j)
          N_t = N_t + N_i(j)
          N_bar = N_bar + N_i(j)
        enddo
        mubar = mubar/N_bar

c correct for the diffusion limit determined by difference of inverse scale heights
c check this vs spreadsheet - not the same
        fH2 = N_i(LH2)/N_t
        BBB = ((mass(LCO2)-mass(LH2))/(mubar-mass(LH2)))**2
        H2_esc = -A_escape*Sun(2)*SSX*(fH2/(1+3.3*fH2))
     $     /sqrt(1 + B_escape*(Sun(2)*SSX)**2*BBB)
c  this can be generalized for differenf insolations and different gravities

c **************   Photolysis   ***************   stuff   ****************
c compute the "tau"s that weight the absorbers in each wavelength bin
c in spreadsheet tau is weighted by pressure
c  here tau is an actual tau for the entire column
c  But tau's should be weighted by N in the stratosphere
c  rather than do that, I weight the tau_uv to the entire column
c   as a practical matter the tau_uv may always be too small to matter
        do j=1,Nsp  ! initialize variables
          Phi_ion(j) = 0.0
          Phi(j)     = 0.0
        enddo

        do L=1,Nsun
          tau(L) = tau_uv*N_t/NX  ! weighted
        enddo

        do L=1,Nsun
          do j=1,Nsp
            tau(L) = tau(L) + sigma(L,j)*N_i(j)
          enddo
C        print *, i, L, tau(L)
        enddo
c       print *, 'uv', i, tau_uv*N_t/NX

c ion production
        Do j=1,Nsp
          Phi_ion(j) = N_i(j)*sigma(1,j)*flux(1)/4/tau(1)
        Enddo

c CO2  ******************** Total  ********************  Photolysis *********************
        Phi(LCO2) = 0.
        do L=2,6
          Phi(LCO2) = Phi(LCO2)
     $      + Flux(L)/4*sigma(L,LCO2)*N_i(LCO2)/tau(L)
        enddo
        Phi(LCO2) = Phi(LCO2) + 0.2*Phi_ion(LCO2)
c       print *, i, 'Phi(LCO2)', Phi(LCO2)/1e11
c       print *, i, 'Phi_ion(LCO2)', Phi_ion(LCO2)/1e11


c CH4 ********************  Total ********************   Photolysis ********************
        Phi(LCH4) = 0.
        do L=1,4
          Phi(LCH4) = Phi(LCH4)
     $      + Flux(L)/4*sigma(L,LCH4)*N_i(LCH4)/tau(L)
        enddo
c CH4 breakup by ions.  this is a mess
        sum_ions = 0.
        do j=1,6
          sum_ions = sum_ions + N_i(j)  !  *k(8,j)  not done right in sheet
        enddo
c       CH4_ion_breakup = Phi_H2_ions*k(8,2)*N(i,LCH4)/sum_ions
        CH4_ion_breakup = Phi_ion(LH2)/sum_ions
     $    + Phi_ion(LN2)*(N_i(LCH4)/(sum_ions - N_i(LN2))
     $    + N_i(LCO)/(sum_ions -N_i(LN2))*N_i(LCH4)    ! line 1  open parens
     $       /(sum_ions -N_i(LN2) -N_i(LCO))           ! line 2
     $    + N_i(LCO)/(sum_ions -N_i(LN2))*N_i(LCO2)    ! line 3
     $       /(sum_ions -N_i(LN2) -N_i(LCO))           ! line 5
     $    * N_i(LCH4)/(N_i(LH2)+N_i(LCH4)+N_i(LH2O))   ! line 5
     $    + N_i(LH2O)/(sum_ions - N_i(LN2))            ! line 6
     $    * N_i(LCH4)/(N_i(LCH4)+N_i(LH2O))  )         ! line 7 close parens
     $    + Phi_ion(LCO)*( N_i(LCH4)/(N_i(LCO2)+N_i(LCH4)+N_i(LH2O))   ! open paren
     $    + N_i(LCO2)/(N_i(LCO2)+N_i(LCH4)+N_i(LH2O))
     $    * N_i(LCH4)/(N_i(LH2)+N_i(LCH4)+N_i(LH2O)) )  ! line 3 close  parens
     $    + Phi_ion(LCO2)*N_i(LCH4)
     $      /(N_i(LH2)+N_i(LCH4)+N_i(LH2O))
     $    + Phi_ion(LH2O)*N_i(LCH4)/(N_i(LH2)+N_i(LCH4))

        Phi(LCH4) = Phi(LCH4) + CH4_ion_breakup

c H2O  ******************** Total  ********************  Photolysis ********************
        Phi(LH2O) = 0.
        do L=2,6
          Phi(LH2O) = Phi(LH2O)
     $      + Flux(L)/4*sigma(L,LH2O)*N_i(LH2O)/tau(L)
        enddo
c        print *, i, 'Phi(LH2O)', Phi(LH2O)/1e11

c N2  ******************** Total  ********************  Photolysis ********************
        Phi(LN2) = 0.
        do L=2,3
          Phi(LN2) = Phi(LN2)
     $      + Flux(L)/4*sigma(L,LN2)*N_i(LN2)/tau(L)
        enddo
c       print *, i, 'Phi(LN2)', Phi(LN2)

c NH3  ******************** Total  ********************  Photolysis ********************
        Phi(LNH3) = 0.
        do L=1,7
          Phi(LNH3) = Phi(LNH3)
     $      + Flux(L)/4*sigma(L,LNH3)*N_i(LNH3)/tau(L)
        enddo
c       print *, i, 'Phi(LN2)', Phi(LN2)

c  End *********************** Photolyses   *************************

c  Begin  ************* Chemistry **************************

c O1D
        sumO1D = 0.
        do j=1,6
          sumO1D = sumO1D + k(2,j)*N_i(j)
        enddo
        O1D_H2  = k(2,LH2)*N_i(LH2)/sumO1D
        O1D_CH4 = k(2,LCH4)*N_i(LCH4)/sumO1D
c O1D > OH
        O1D_OH = ( k(2,LH2)*N_i(LH2) + k(2,LCH4)*N_i(LCH4)
     $     + 2*k(2,LH2O)*N_i(LH2O) )/sumO1D
c OH
        sumOH = k(3,6)*others3   ! others3 is a background column
        do j=1,5
          sumOH = sumOH + k(3,j)*N_i(j)
        enddo
        OH_CH4 = k(3,LCH4)*N_i(LCH4)/sumOH
        OH_H2  = k(3,LH2)*N_i(LH2)/sumOH
        OH_CO  = k(3,LCO)*N_i(LCO)/sumOH
c       print *, i, 'OH_H2', OH_H2
c       print *, i, 'OH_CH4', OH_CH4
c       print *, i, 'OH_CO', OH_CO
c N > HCN
        sumN2D = 0.
        do j=1,6
          sumN2D = sumN2D + k(1,j)*N_i(j)
        enddo
        N_HCN = (k(1,LCH4)*N_i(LCH4) + k(1,LCO)*N_i(LCO)
     $      + k(1,LN2)*N_i(LN2) )/sumN2D
c       print *, i, 'N_HCN', N_HCN

c  O + CO > CO2
        sumO = k(4,6)*others4   ! others3 is a background column
        do j=1,5
          sumO = sumO + k(4,j)*N_i(j)
        enddo
        O_CH4 = k(4,LCH4)*N_i(LCH4)/sumO
        O_H2  = k(4,LH2)*N_i(LH2)/sumO
        O_CO  = k(4,LCO)*N_i(LCO)/sumO

c $$$$$$$$$$$$$   Budgets   $$$$$$$$$$$$$   Budgets   $$$$$$$$$$$$$$$

c  oxidation of CH4 to make...
        dCH4dt_ox = -Phi(LCO2)*O1D_CH4 -Phi(LH2O)*OH_CH4
c  NH3 > HCN
        NH3_HCN = (Phi(LCH4)-dCH4dt_ox)
     $   / (Phi(LCH4)-dCH4dt_ox + Phi(LH2O) +Phi(LCO2) )
c       print *, i, 'NH3_HCN', NH3_HCN

c  HCN
        dN_idt(LHCN) = 2*Phi(LN2)*N_HCN*
     $     (Phi(LCH4)-dCH4dt_ox)
     $      /( Phi(LCH4)-dCH4dt_ox + 2*Phi(LN2)*N_HCN )
     $    + Phi(LNH3)*NH3_HCN
c        print *, i,'dNdt(LHCN)', dN_idt(LHCN)

c C2Hn
        dN_idt(LC2Hn) = 0.5*(Phi(LCH4)-dCH4dt_ox)**3
     $   / (Phi(LCH4)-dCH4dt_ox + Phi(LH2O) +Phi(LCO2) )**2
c        print *, 'dNdt(LC2Hn)', dN_idt(LC2Hn)

c Haze
        dN_idt(LHaze) = (Phi(LCH4)-dCH4dt_ox) *
     $   ( (Phi(LCH4)-dCH4dt_ox)
     $   /(Phi(LCH4)-dCH4dt_ox + Phi(LH2O) +Phi(LCO2)) )**5
c        print *, 'dNdt(LHaze)', dN_idt(LHaze)

c  CH4
        dN_idt(LCH4) = -Phi(LCH4) + dCH4dt_ox + Phi_geo(LCH4)
     $    - Phi(LNH3)*NH3_HCN

c  CO2
        dN_idt(LCO2) = -Phi(LCO2) + Phi(LCO2)*O1D_OH*OH_CO
     $    + Phi(LH2O)*OH_CO + Phi(LCO2)*(1-O1D_OH)*O_CO
     $    - Phi_geo(LCO2)
c      print *, 'dNdt(LCO2)', dN_idt(LCO2)
c  CO
        dN_idt(LCO)  = Phi(LCO2)-(dCH4dt_ox-Phi(LCH4))
     $   - dN_idt(LHaze) - Phi(LCO2)*O1D_OH*OH_CO
     $   - Phi(LH2O)*OH_CO - Phi(LCO2)*(1-O1D_OH)*O_CO
     $   - dN_idt(LHCN) + Phi(LNH3)*NH3_HCN   ! correction
c       print *, 'dNdt(i,LCO)', dN_idt(LCO)
c  H2
        dN_idt(LH2) = H2_esc - 2*(dCH4dt_ox-Phi(LCH4))
     $   + 2*dN_idt(LCO2) + dN_idt(LCO) - dN_idt(LHaze) - dN_idt(LHCN)
     $   + Phi_geo(LH2) + 3*Phi(LNH3)*NH3_HCN
c  H2O
        dN_idt(LH2O) = 0.0     ! assumed constant
c  NH3
        dN_idt(LNH3) = -Phi(LNH3)     ! assumed lost
c N2  !
        dN_idt(LN2) = -0.5*dN_idt(LHCN) +0.5*Phi(LNH3)  ! note that *(1-NH3_HCN) is wrong here
c       print *,i,'dNdt(i,LN2)',dN_idt(LN2)
c       stop

      return
      end
