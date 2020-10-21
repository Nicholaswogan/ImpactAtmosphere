      Program Temperature_loop
c idea is to loop over temperature, call subroutine at each T
c  send initial conditions, return partial pressures

      parameter(Nu=77,Nsp=7)
      implicit real*8(A-H,O-Z)
      real*8 p_tot, g  ! total pressure, gravity
c  pressure must be in atmospheres
      real*8 Area   ! of Earth cm^2
      real*8 mass(Nsp)   ! molar mass [g], Columns moles cm^-2
      real*8 m_i, Moles_Fe, N_Fe, mi    ! impactor mass, total moles Fe, column Fe moles cm^-2
      real*8 pCO, pCO2, pCH4, pH2, pH2O, pNH3, pN2  ! inputs [bars] .ne. partial pressures
c input pressures are reservoirs, same as bars if single component atmosphere
      real*8 p_i(Nsp), p_j(Nsp)  ! first guess, actual partial pressures, atmospheres
      real*8 TQCH4(Nu), TQNH3(Nu)
      real*8 fraction, p_q(Nsp), N_q(Nsp)
      real*8 N_i(Nsp), N_tt, N_ttt         ! columns of unprocessed atm moles per cm2
      real*8 p_ttt
      real*8 N_j(Nsp), N_tot   ! columns moles per cm2, returned from Impact_gases
      real*8 E_i, v_i
      real*8 Q_w, C_H2O, C_CO2, C_N2, C_CO, C_w, M_w
      real*8 M_N2, M_CO2, M_H2O, M_CO, M_A, M_AC_A, mu_A, N_A
      real*8 minpH2O   ! minimum H2O for global event input in bars
      real*8 M_min_H2O ! minimum H2O for global event converted to mass
      real*8 mu_v, p_v, x, mH, M_steam, maxheat, minheat
      real*8 F_ir, T_w, T_0, T_h, T_c
      real*8 dT, dx, det, F1, F2, dF1dx,dF2dx, dF1dT, dF2dT
      real*8 eta      ! impact > water vapor coupling  0<eta<1
      
      real*8 NN_tot,NM_tot,pp_tot  ! used for pressures of post-impact dry atmosphere

c we store this stuff
      real*8 Temperature(Nu), pressure(Nu), p(Nu,Nsp)
      real*8 N(Nu,Nsp), N_total(Nu)
      integer test_N_O   ! not used right now

      integer buffer     ! flag.  buffer=1 if no buffer
                         !        buffer=2 QFM
                         !        buffer=3  QFM-1
                         !        buffer=4  QFM-2
                         !        buffer=5 IW
                         !        buffer=6 QFI
      Character*8 boof   !   for output
      character*8 planet !   Earth and Mars are only ones.  Can do Venus
      

      CHARACTER*8 ISPEC(Nsp)
      Common /ablok/ ISPEC
      Common /bblok/ mass, g, Area
      Common /Lblok/LH2O,LH2,LCO,LCO2,LCH4,LN2,LNH3

c     Data g/982./  ! cm/s^2
      Data mH/1.66e-24/  ! grams
      Data mass/18.,2.,28.,44.,16.,28.,17./  ! per mole
c     Data Area/5.1e18/  ! cm^2
c     data v_i/17.e5/  ! impact velocity cm/s - update for planet
      data F_ir/1.5e5/  ! radiative cooling ergs/cm2/s
      data C_H2O,C_w, C_CO2, C_N2,C_CO/2e7,2e7, 8e6,1.1e7,1.1e7/ ! heat capacities of gases, ergs/g/K
      data Q_w/2.5e10/  ! latent heat vaporization water, ergs/g
      
      data Fe_fraction/ 0.34 /

      data buffer/1/
      data T_h/1600/  ! hot steam, atm by impact
      data T_c/650/   ! cooled T of steam atm
      data T_w/5000/  ! water vapor saturation
      data T_0/300/   ! pre-impact T
      data T_qc/800/  ! CH4 quench T approx
      data T_qn/1100/ ! NH3 quench T approx

c some labels
      DATA LH2O,LH2,LCO,LCO2,LCH4,LN2,LNH3/1,2,3,4,5,6,7/

      open(45, file='IW.parameters',status='UNKNOWN')  ! output
      open(54, file='IW.pressure',status='UNKNOWN')  ! output
      open(56, file='IW.drypressure',status='UNKNOWN')  ! output
      open(55, file='IW.column',status='UNKNOWN')  ! output
      open(59, file='IW.out',status='UNKNOWN')  ! output
      open(69, file='IW.outcolumns',status='UNKNOWN')  ! output
      open(39, file='IW.input',status='OLD')  ! output

C  Define names of species
      ISPEC(1) =   'H2O'
      ISPEC(2) =   'H2'
      ISPEC(3) =   'CO'
      ISPEC(4) =   'CO2'
      ISPEC(5) =   'CH4'
      ISPEC(6) =   'N2'
      ISPEC(7) =   'NH3'

c  assign initial atm
c  these will be read in from the shell, eventually
c  these are reservoirs of single component atmospheres!
c  they need to be converted to moles
      read (39,*) buffer, mi, pCO, pCO2, pH2, pH2O, pCH4, pN2, pNH3
     $    , planet, minpH2O, eta

      m_i = 10.**mi
      print *, pCO, pCO2, pH2, pH2O, pCH4, pN2, pNH3
      print *, m_i, planet, minpH2O
      
c pick a planet
      if (planet.eq.'Earth') then
        area = 5.1e18
        g = 982.
        v_i = 17.e5
      elseif (planet .eq. 'Mars') then
        area = 1.4e18
        g = 375.
        v_i = 10.e5
      else 
        stop 666
      endif
      print *, area, g, v_i/1e5
      

      if (buffer .eq. 1) then
        boof = 'nobuffer'
      elseif (buffer .eq. 2) then
        boof = 'QFM'
      elseif (buffer .eq. 5) then
        boof = 'IW'
      elseif (buffer .eq. 6) then
        boof = 'QFI'
      elseif (buffer .eq. 3) then
        boof = 'QFM-1'
      elseif (buffer .eq. 4) then 
        boof = 'QFM-2'
      else
        boof = 'oops'
        stop 777
      endif  
      
      GOTO 98   ! the following stuff is  for code testing
        pCO     = 0.
        pCO2    = 100.
        pH2     = 0.0
        pH2O    = 500.
        pCH4    = 0.
        pN2     = 1.0
        pNH3    = 0.0
        minpH2O = 10.
 98   CONTINUE

c inventories, masses
      M_CO2     = pCO2*Area*1e6/g
      M_CO      = pCO*Area*1e6/g
      M_N2      = pN2*Area*1e6/g
      M_H2O     = pH2O*Area*1e6/g
      M_min_H2O = minpH2O*Area*1e6/g   !  30 bars is minimum for a proper runaway greenhouse

c input impact mass
c this will be read in from the shell
c these quantities are used if buffer=1
c     m_i = 3.1e23 ! grams
      Moles_Fe = Fe_fraction*m_i/56.
      N_Fe     = Moles_Fe/Area   !  Moles cm^-2

c cooling time *** this needs updating
      E_i = 0.5*m_i*v_i**2

c solve for mass of water evaporated
      M_A = M_CO2 + M_N2 + M_CO   ! total mass dry atm
      M_AC_A = M_CO2*C_CO2 + M_N2*C_N2 + M_CO*C_CO  ! total heat capacity dry atm

c the next two are only needed for Case 3
      mu_A = mH*M_A/(M_CO2/mass(LCO2) + M_N2/mass(LN2) + M_CO/mass(LCO))  ! mean molecular weight dry atm
      N_A = M_A/mu_A/Area  ! column number density dry atm molecules per cm2

c  assume that the heated atm/steam is heated to T_h
      maxheat = M_H2O*(Q_w + (T_h-T_0)*C_w) + M_AC_A*(T_h-T_0)
      minheat = M_min_H2O*(Q_w + (T_h-T_0)*C_w) + M_AC_A*(T_h-T_0)
      print *, m_i,eta*E_i, maxheat, minheat , M_AC_A*(T_h-T_0)

c  there are 3 cases - all the water evaporates, or not all evaporates
      If (eta*E_i .gt. maxheat) Then  ! all the water is evaporated is heated to T_h
        M_steam = M_H2O
        scale = 1
        print *, 'maxheat','scale', scale
      Elseif (eta*E_i .gt. minheat) then ! then only some of the water evaporates
        M_steam = (eta*E_i - M_AC_A*(T_h-T_0))/(Q_w + (T_h-T_0)*C_w)
        scale = 1
        print *, 'minheat', 'scale', scale
      Else ! only a fraction of the gases are heated.  This still doesn't work right
c I don't know how to do this yet
c  its probably all just invalid.
        M_steam = M_min_H2O 
        scale=eta*E_i/(M_min_H2O*(Q_w+(T_h-T_0)*C_w)+M_AC_A*(T_h-T_0))
c       scale=1
c  ! don't do anything.  just stuff values into output
        print *, '888', m_i, M_steam, M_min_H2O, scale
c       stop 888
      Endif
c     print *, m_i, M_steam, M_min_H2O, scale

c  Case 4 -  the water does not all evaporate, and eventually           ******
c   everything equilibrates.  This needs to be solve for total p and T  ******
c        two variables  x = total pressure cgs                          ******
c                       T = surface T                                   ******
c  This is really cool but entirely beside the point.  Its all debugged.******
c   Maybe I'll use it someday.                                          ******
      GOTO 99   ! bur skip it for now                                   ******
        x = (pCO2 + pN2 + pCO + 10. )*1.e6   !   case 3 first guesses   ******
        T = 600.                             !   case 3 first guesses   ******
        mu_v = mass(LH2O)*mH                 !                          ******
        do iloop =1,9  ! its  working but I don't understand why        ******
          p_v = 1.e6*exp(T_w/373-T_w/T)                             !   ******
          F1 =(x*mu_A + p_v*(mu_v-mu_A))*(N_A*(mu_v-mu_A)*g + x)    !   ******
     $       -x**2*mu_v                                             !   ******
          F2 =(Area*x/g - M_A)*(Q_w + (T-T_0)*C_w) + M_AC_A*(T-T_0) !   ******
     $       -eta*E_i                                               !   ******
          dF1dx = mu_A*(N_A*(mu_v-mu_A)*g + x) + x*mu_A             !   ******
     $      + p_v*(mu_v-mu_A) - 2*x*mu_v                            !   ******
          dF1dT = p_v*T_w/T**2*(mu_v-mu_A)*(N_A*(mu_v-mu_A)*g + x)  !   ******
          dF2dx = Area/g*(Q_w +(T-T_0)*C_w)                         !   ******
          dF2dT = (Area*x/g -M_A)*C_w + M_AC_A                      !   ******
          det = dF1dx*dF2dT - dF1dT*dF2dx                           !   ******
          dx  = (dF1dT*F2 - dF2dT*F1)/det                           !   ******
          dT  = (dF2dx*F1 - dF1dx*F2)/det                           !   ******
          x = x+dx                                                  !   ******
          T = T+dT                                                  !   ******
        enddo                                                       !   ******
        p = x  ! total pressure                                         ******
        M_steam = Area*x/g -M_A                                     !   ******
  99  CONTINUE                                                      !   ******
c  End of Case 3 -  the water does not all evaporate                    ******


      time_CH4 = (C_w*(T_h-T_qc)*M_steam+M_AC_A*(T_h-T_qc))/(F_ir*Area)
      time_NH3 = (C_w*(T_h-T_qn)*M_steam+M_AC_A*(T_h-T_qn))/(F_ir*Area)
      time_CH4 = time_CH4 *scale
      time_NH3 = time_NH3 *scale
      print *, time_CH4/1e10, time_NH3/1e10

c assign these as first guess pressures 
      p_i(LCO)  = pCO 
      p_i(LCO2) = pCO2 
      p_i(LH2)  = pH2 
      p_i(LCH4) = pCH4 
      p_i(LH2O) = pH2O*M_steam/M_H2O 
      p_i(LNH3) = pNH3 
      p_i(LN2)  = pN2 
              
c  convert initial "p_i(L)" into initial Columns (Moles cm^-2)
c  these 
c  this is the atmosphere pre-impact
      Do i=1,Nsp
        N_i(i) = p_i(i)*1.e6/g/mass(i)   ! mass(i) is of species in amu
      Enddo
      N_tt = 0.
      Do i=1,Nsp
        N_tt = N_tt + N_i(i)    ! mass(i) is of species in amu
      Enddo
      N_tot = N_tt  ! put something here
      
      print *, 'N_tt',N_tt
      
      DO iloop=1,Nu  ! outer loop

        T = 1420. - float(iloop)*10.
        Temperature(iloop) = T
       
c       print *, 'p_i(LH2O)', iloop, T, p_i(LH2O), p_i(LCO2)

c   I can put a decision here for buffer
c       print *, N_tot

        call Impact_gases(buffer,N_Fe,Moles_Fe,p_i,T,p_j,p_tot,
     $   N_j, test_N_O, N_tot)

c the returned values refer to conditions in the fraction "scale"
c   of the atmosphere.  I haven't got this sussed out.
c   pressure should use all the atmosphere, but that would require a different
c   algorithm for shocking, I think.   mixing/cooling is so fast,
c   I suspect one should just set Tquench to say 1500 and live with
c   that, while mixing unshocked and shocked
c   the problem is that the two end-member approximations will not meet
c   in the middle. 

c    
c  store returned values 
        do i=1,7
c         p(iloop,i) = p_j(i)   ! this is only correct if scale = 1
                                ! if scale<1, then pressure needs to recomputed
c         N(iloop,i) = N_j(i)   !  this is correct if scale=1    
          N(iloop,i) = N_j(i)*scale + (1-scale)*N_i(i)  ! this should be correct if scale<1      
        enddo
c       pressure(iloop) = p_tot  ! correct only if scale=1  

c when scale<1     These pressures and numbers include H2O
        P_ttt = 0.
        N_ttt = scale*N_tot + (1-scale)*N_tt
        N_total(iloop) = N_ttt
        do i=1,7
          P_ttt = p_ttt + N(iloop,i)*g*mass(i)/1e6   ! P_ttt is total bars
        enddo
        pressure(iloop) = P_ttt      ! this should be correct when scale < 1 
        do i=1,7
          p(iloop,i) = N(iloop,i)/N_ttt*p_ttt
        enddo
        print *, iloop, 'pressure=', pressure(iloop), p_tot, P_ttt

      ENDDO  ! end outer Temperature loop


c **********  QUENCHING   ************  QUENCHING   ************
c quench T - check these
      Do j=1,Nu
        TQCH4_1 = 42000/LOG(time_CH4*pressure(j)/0.000003)
        TQCH4_2 = 25000/LOG(time_CH4*pressure(j)**2/40)
        TQCH4(j) = min(TQCH4_1,TQCH4_2,1400.)

        TQNH3(j) = 52000/LOG(time_NH3*pressure(j)/0.0000001)
        TQNH3(j) = min(1400., TQNH3(j)) !  The quench T has to be in the range of T(Nu) 
      Enddo

c **********  QUENCHING   ************  QUENCHING   ************
c find H2, H2O, CO, CO2, CH4
      Do j=1,Nu
        If (TQCH4(j).gt.Temperature(j)) Then
          jCH4 = j   ! T(jCH4) < TQCH4 < T(jCH4-1)
          goto 20
        ENDIF
      Enddo
      stop 20  ! should not get here
 20   Continue

c linear interpolation
      fraction = (TQCH4(jCH4)-Temperature(jCH4))
     $  /(Temperature(jCH4-1)-Temperature(jCH4))
      Do i=1,5
        p_q(i) = p(jCH4,i) + fraction*(p(jCH4-1,i)-p(jCH4,i))  ! valid only if scale=1
        N_q(i) = N(jCH4,i) + fraction*(N(jCH4-1,i)-N(jCH4,i))
      Enddo
c     print *, jCH4, fraction, (p_q(i),i=1,5)

c **********  QUENCHING   ************  QUENCHING   ************
c find N2, NH3
      Do j=1,Nu
        If (TQNH3(j).gt.Temperature(j)) Then
          jNH3 = j   ! T(jNH3) < TQNH3 < T(jNH3-1)
          goto 30
        ENDIF
      Enddo
      stop 30  ! should not get here
 30   Continue

c linear interpolation (assumes trace species for N2, NH3)
      fraction = (TQNH3(jNH3)-Temperature(jNH3))
     $  /(Temperature(jNH3-1)-Temperature(jNH3))
      Do i=6,7
        p_q(i) = p(jNH3,i) + fraction*(p(jNH3-1,i)-p(jNH3,i))
        N_q(i) = N(jNH3,i) + fraction*(N(jNH3-1,i)-N(jNH3,i))
      Enddo
c     print *, jNH3, fraction, (p_q(i),i=1,Nsp)

c **********  OUTPUT   ************  OUTPUT   ************  OUTPUT   ************

c p_q is pressure including the water vapor
c  should also compute pressures after the water has condensed
c  use p_i(j) for this
c  
      NN_tot = -N_q(LH2O)
      NM_tot = -N_q(LH2O)*mass(LH2O)
      do j=1,Nsp
        NN_tot = NN_tot + N_q(j)        
        NM_tot = NM_tot + N_q(j) * mass(j)       
      enddo
      pp_tot = NM_tot*g/1e6          ! total pressure bars

c dry pressure  - this is not jiving with the 
      do j=1,Nsp
        p_i(j) = N_q(j)*pp_tot/NN_tot        
      enddo
      p_i(LH2O) = 0.
 
c **********  OUTPUT   ************  OUTPUT   ************  OUTPUT   ************

c  buffer is printed both as an integer buffer "N_O"  (atavism)
c  and as text name "boof"
c  the integer is for plotting with datagraph, for which I've yet to
c  figure out how to plot things against text.  Its annoying.
c  T_q are quench temperatures
c  pCO2 is input pCO2 inventory
c  other input inventories - minpH2O, pH2O, pCO, pN2 should be printed also in a header

 101  format(1x,'planet=',A8,1x,'buffer=',A8,1x,'pCO2=',F10.3,
     $  'pH2O=',F10.3,'minpH2O=',F10.3,'pN2=',F10.3,'pCO=',F10.3)

 103  Format(3x,'i',3x,'N_O',3x,'buffer',4x,'planet',7x,
     $   'PCO2',4x,'T_q(CH4)',2x,'T_q(NH3)',
     $   3x,'T[K]',4x,'P[bars]',6x,7(A8,3x))
 104  Format(1x,2I4,4x,A8,2x,A8,4F10.2,1P8E11.3)
 105  Format(12x,'M_i',4x,'N_O',3x,'buffer',3x,'planet',7x,
     $    'PCO2',5x,'T_q(CH4)',
     $  2x,'T_q(NH3)',13x,7(A8,3x))
 106  Format(6x,1PE11.2,I4,3x,A8,2x,A8,1x,0P3F10.2,10x,1P7E11.3)
 107  Format(9x,'M_i',6x,'N_O',4x,'buffer',4x,'planet',7x,
     $   'pCO2',4x,'T_q(CH4)',2x,'T_q(NH3)',
     $  4x,7(A8,3x))
 108  Format(6x,1PE11.3,I4,3x,A8,2x,A8,1x,0P3F10.2,1x,1P7E11.3)

      write(59,105) (ISPEC(i),i=1,7)
      write(59,106) M_i, buffer, boof, planet, pCO2,
     $   TQCH4(jCH4),TQNH3(jNH3), (N_q(i),i=1,7)
      write(59,106) M_i, buffer, boof, planet, pCO2,
     $   TQCH4(jCH4),TQNH3(jNH3), (p_q(i),i=1,7)
      write(59,106) M_i, buffer, boof, planet, pCO2,
     $   TQCH4(jCH4),TQNH3(jNH3), (p_i(i),i=1,7)

c **********  OUTPUT   ************  OUTPUT   ************  OUTPUT   ************
      write(59,101) planet,boof,pCO2,pH2O,minpH2O,pN2,pCO  ! header
      write(59,103) (ISPEC(i),i=1,7)
      Do j=1,Nu
        write(59,104) j, buffer, boof, planet, pCO2, TQCH4(j),
     $   TQNH3(j), Temperature(j), pressure(j), (p(j,i),i=1,7)
      Enddo

      write(69,101) planet,boof,pCO2,pH2O,minpH2O,pN2,pCO  ! header
      write(69,103) (ISPEC(i),i=1,7)
      Do j=1,Nu
        write(69,104) j, buffer, boof, planet, pCO2, TQCH4(j),
     $   TQNH3(j), Temperature(j), N_total(j), (N(j,i),i=1,7)
      Enddo

c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      write(45,*) N_q(LCO), N_q(LCO2), N_q(LH2), N_q(LH2O),
     $  N_q(LCH4), N_q(LN2), N_q(LNH3)

c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      write(54,107) (ISPEC(i),i=1,7)
      write(54,108) M_i, buffer, boof, planet, pCO2,
     $   TQCH4(jCH4),TQNH3(jNH3), (p_q(i),i=1,7)

c   dry pressure
      write(56,107) (ISPEC(i),i=1,7)
      write(56,108) M_i, buffer, boof, planet, pCO2,
     $   TQCH4(jCH4),TQNH3(jNH3), (p_i(i),i=1,7)

      write(55,107) (ISPEC(i),i=1,7)
      write(55,108) M_i, buffer, boof, planet, pCO2,
     $   TQCH4(jCH4),TQNH3(jNH3), (N_q(i),i=1,7)

      stop
      end

c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      Subroutine Impact_gases(buffer,N_Fe,Moles_Fe,p_i,T,p_j,p_tot,
     $ N, test_N_O, N_tot )
c input:
c     buffer to decide what buffer to use
c     p_i(Nsp)  input pressures = reservoirs
c     T    temperature
c     Moles_Fe, N_Fe
c     scale is not-yet-implemented attempt to do subglobal impacts
c     the idea is to restrict shocking to a fraction of the atmosphere
c     and then mixing the products with unshocked atmosphere
c     I did this in the spreadsheet, but I'm making some kind of elementary
c     mistake here

c return
c     p_j(Nsp) output pressures  [atm]
c     p_tot   [atm]
c     N(Nsp)   moles/cm2
c     test_N_O  =  1   the removal of O by Fe can work
c     test_N_O  =  -1   the removal of O by Fe fails
c        for these cases, set N_O to a small value and return a flag


      parameter(Nsp=7)
      implicit real*8(A-H,O-Z)
      real*8 fO2, xO2   ! oxygen fugacity, auxilliary variable
      real*8 p_tot, p_t, g  ! total pressure, gravity
      real*8 Mass_tot, Mass_t   ! total mass of atm
c  pressure must be in atmospheres
      real*8 Area   ! of Earth cm^2
      real*8 N_tot, N_C, N_H, N_O, N_t, N_N  ! total column densities, moles cm^-2
      real*8 mass(Nsp), N(Nsp)   ! molar mass [g], Columns moles cm^-2
      real*8 m_i, Moles_Fe, N_Fe    ! impactor mass, total moles Fe, column Fe moles cm^-2
      real*8 pCO, pCO2, pCH4, pH2, pH2O, pNH3, pN2  ! inputs [bars] .ne. partial pressures
c input pressures are reservoirs, same as bars if single component atmosphere
      real*8 p_i(Nsp)  ! input partial pressures, atmospheres
      real*8 p(Nsp)    ! output partial pressures, atmospheres
      real*8 p_j(Nsp)  ! output partial pressures, atmospheres
      real*8 K1, K2, K3, K4, K44, K5  ! equilibrium constants
      real*8 AK1(3), AK2(2), AK3(2), AK4(2), AK5(2)  ! for equilibrium constants
      real*8 AA(4),BB1(3),BB2(3),BB3(3),CC(3)  ! for mineral buffer fugacities
      real*8 IW, QFM, QFI   ! mineral buffer fO2
      real*8 Rgas           ! one of the gas constants

c used for solving
      real*8 oldx, oldy, oldF1, oldF2, oldNt, eps, dx, dy
      real*8 dF1dx, dF1dy, dF2dx, dF2dy, dF3dx, dF3dp, dF1dp
      real*8 deltax, deltay, Det, newx, newy
      real*8 oldp, oldF3, F3, dp, newp, deltap
      
      real*8 aq,bq,cq,y1,y2
      integer test_N_O

      integer buffer   ! flag.  buffer=1 if no buffer
                       !        buffer=2 QFM
                       !        buffer=3  QFM-1
                       !        buffer=4  QFM-2
                       !        buffer=5 IW
                       !        buffer=6 QFI

      CHARACTER*8 ISPEC(Nsp)
      Common /ablok/ ISPEC
      Common /bblok/ mass, g, Area
      Common /Lblok/LH2O,LH2,LCO,LCO2,LCH4,LN2,LNH3

      Data Rgas/8.314/  ! J/K/mole
      Data eps/1.e-3/

c compute Oxygen fugacity as f(T)
c QFM 900 - 1420 K
      DATA AA /0.09271, -203.3164, 1584.427, -587474/
      DATA BB1 /22.446, -33.182, -542941/  ! QFI 900-1042
      DATA BB2 /5.4771,	103.384, -562377/  ! QFI 1042-1184
      DATA BB3 /-27.3443, 369.74, -602739/ ! QFI 1184-1420
      DATA CC/ -8.474, 115.56, -244118/    ! IW
      
c coefficients of 4 equils
c   H2O*CO = K1*H2*CO2
c   (H2)^3*CO*K2 = CH4*H2O
c   (H2O/H2)^2 = fO2*K3
c   (NH3)^2 = K4*(H2)^3*(N2)
      DATA AK1/18.28, -2375.6, -8.69e5/
      DATA AK2/5.239e-14, 27285./
      DATA AK5/1.16e-6, 59911./

      DATA AK4/5.9e-13, 13207./

 101  Format(5x,'T',10x,7(A8,3x))
 102  Format(1x,F10.2,1P7E11.3)

      test_N_O = 1   ! if N_O < 0, this algorithm fails.

      Do i=1,Nsp
        p(i) = p_i(i)       ! set p(i) = input pressures = total reservoirs
      Enddo
c this next statement is true for the initial conditions, even though the p(i) need to be renormalized
      p_tot = 0.
      Do i=1,Nsp
        p_tot = p_tot + p(i)       ! total pressure of atmosphere 
      Enddo
      Mass_t = p_tot*1.e6/g*Area    ! total mass of atm
      
c  convert initial p(L) to initial Columns (Moles cm^-2)
c  THIS IS USED
      Do i=1,Nsp
        N(i) = p(i)*1.e6/g/mass(i)   ! mass(i) is of species in amu
      Enddo

c  compute composition as function of Temperature

c  *********************  ********************  *********************

c buffer should be read in
        IF (buffer .eq. 1) THEN  ! no buffer.  Do the conserved O case

c these are the reservoirs BEFORE the reaction with Fe
          N_tot = 0.0
          Do i=1,Nsp
            N_tot = N_tot + N(i)    ! Moles cm^-2
          Enddo

c compute partial pressures before decrementing FeO
          Do i=1,Nsp
            p(i) = N(i)/N_tot*p_tot
          Enddo

c  with buffer=1,  Oxygen is removed by Fe and pressure is subsequently constant
c   compute constant columns, taking into account loss of O to FeO
          N_C = N(LCO) + N(LCO2) + N(LCH4)
          N_H = N(LH2) + N(LH2O) + 2*N(LCH4)
          N_O = N(LCO) + 2*N(LCO2) + N(LH2O) - N_Fe
          N_N = 2*N(LN2) + N(LNH3)
          print *, 'N_O', N_O
          print *, N_H
          print *, N_C
          print *, N_N

c     It is not allowed for N_O < 0  Insert a test
          If (N_O .le. 0.) Then
            print *, 'N_O', N_O   ! need an error return
            test_N_O = -1
            N_O = 100.   ! something small as a placeholder
c           return   ! this will be a case where buffer has to be used
          Endif

c    If the Fe is exhausted, we use the case where O is scavenged from atmosphere
c    the total mass of the atmosphere is less massive
          Mass_tot = Mass_t - Moles_Fe*16   ! remove mass of O
c    the total pressure of the atmosphere is smaller
          p_tot = Mass_tot*g/Area/1.e6   ! bars

          K1 = AK1(1)*exp(AK1(2)/T+AK1(3)/T**2)  ! = CO*H2O/(CO2*H2)
          K2 = AK2(1)*exp(AK2(2)/T)*p_tot**2     ! K2*p_tot^2

c  there are 2 variables.  N_t is needed as another variable to put
c  the various species into partial pressures.
c  x = H2 and y = H2O  (moles/cm2)
c first guess
c  this needs to be done properly for an automagic numerical scheme
c  so I will make the guess by assuming that C is "negligible"
c
          y = N_O - 2*N_C   ! low estimate  I need a better first guess when C is not neglible
          x = N_H - y       ! this is can be too high 
          print *, x,y

c  a better 1st approx includes CO and CO2 but omits CH4.
c  the result is a quadratic for y
          aq = 1-K1
          bq = 2*N_C - N_O + K1*(N_O+N_H-N_C)
          cq = K1*N_H*(N_O-N_C)
          y1 = -0.5*bq/aq + 0.5*sqrt(bq**2 +4*aq*cq)/aq
          y2 = -0.5*bq/aq - 0.5*sqrt(bq**2 +4*aq*cq)/aq
c         print *, y1, y2
c y1 looks like its the right one
          y = y1
          x = N_H - y 

          do j=1,12

            N_t = x + y + N_C + N(LN2)  ! this ignores NH3
            F1 = (y-N_O)*(1+y/x/K1 +K2*x**3/y/N_t**2) +N_C*(1+2*y/x/K1)  ! = 0
            F2 = (x+y-N_H)*(1+y/x/K1) +K2*x**3/y/N_t**2*(x+y-N_H+2*N_c)

c    Jacobean terms, numerically evaluated
c       store old values
            oldx = x
            oldy = y
            oldF1 = F1
            oldF2 = F2
            oldNt = N_t

c       numerically evaluate dF/dx terms
            x = oldx*(1+eps)
            dx = oldx*eps
            N_t = x + y + N_C + N(LN2)  ! this ignores NH3
            F1 = (y-N_O)*(1+y/x/K1+K2*x**3/y/N_t**2) +N_C*(1+2*y/x/K1)
            F2 = (x+y-N_H)*(1+y/x/K1)+K2*x**3/y/N_t**2*(x+y-N_H+2*N_c)
            dF1dx = (F1-oldF1)/dx
            dF2dx = (F2-oldF2)/dx

c       dF/dy terms
            x = oldx   ! restore old x
            y = oldy*(1+eps)
            dy = oldy*eps
            N_t = x + y + N_C + N(LN2)  ! this ignores NH3
            F1 = (y-N_O)*(1+y/x/K1+K2*x**3/y/N_t**2) +N_C*(1+2*y/x/K1)
            F2 = (x+y-N_H)*(1+y/x/K1)+K2*x**3/y/N_t**2*(x+y-N_H+2*N_c)
            dF1dy = (F1-oldF1)/dy
            dF2dy = (F2-oldF2)/dy

c       determinant
            Det = dF1dx*dF2dy - dF1dy*dF2dx

C       restore old values
            y = oldy
            N_t = oldNt
            F1 = oldF1
            F2 = oldF2
c       estimate new x, y
            deltax = (dF1dy*F2 - dF2dy*F1)/Det
            deltay = (dF2dx*F1 - dF1dx*F2)/Det
            newx = x + deltax
            newy = y + deltay

c           print *, j, x, newx
c           print *, j, y, newy

            x = newx
            y = newy
          enddo

c       assign columns moles/cm2
          N(LH2)  = x
          N(LH2O) = y
          N(LCO)  = N_C/(1 + y/x/K1 + K2*x**3/y/N_t**2)
          N(LCO2) = N(LCO)*y/x/K1
          N(LCH4) = N_C - N(LCO) - N(LCO2)
          N_tot = N_t
          print *, N_t

c   nitrogen as minor species
          K4 = AK4(1)*exp(AK4(2)/T)
          K44 = K4 * (p_tot/N_t)**2 *x**3
          z = sqrt(K44**2/4 + 2*K44*N_N)/2 - K44/4
          N(LNH3) = z
          N(LN2) = 0.5*(N_N - z)
          print *, T, 'N(LNH3)',N(LNH3)

c    assign partial pressures
          do i=1,7
            p(i) = p_tot*N(i)/N_tot
          enddo

c         write(*,101) (ISPEC(i),i=1,7)
c         write(*,102) T,(p(i),i=1,7)

          do i=1,7
            p_j(i) = p(i)  ! return pressures
          enddo

c  *********************  ********************  *********************
        ELSE  ! the buffer case

c  bugs!  need to check this

c set up
c  in this case there is only the one variable x = N_H2
c  H2O is given by K5 with fO2 from the buffer
c  N_t and p_t are both variables, N_C and N_H are both conserved, and N_O is ignored

          If (buffer.eq.5) Then  ! IW
            xO2 = CC(3) +CC(2)*T +CC(1)*T*Log(T)
            IW = EXP(2*xO2/(Rgas*T))
            fO2 = IW   ! pressure of O2 in bars
          ElseIf (buffer.eq.4) Then   ! QFM-2
            xO2 = AA(4) +AA(3)*T +AA(2)*T*Log(T) +AA(1)*T**2
            QFM = EXP(xO2/(Rgas*T))
            fO2 = QFM/10.**2   ! pressure of O2 in bars
          ElseIf (buffer.eq.3) Then   ! QFM-1
            xO2 = AA(4) +AA(3)*T +AA(2)*T*Log(T) +AA(1)*T**2
            QFM = EXP(xO2/(Rgas*T))
            fO2 = QFM/10.   ! pressure of O2 in bars
          ElseIf (buffer.eq.2) Then  ! QFM
            xO2 = AA(4) +AA(3)*T +AA(2)*T*Log(T) +AA(1)*T**2
            QFM = EXP(xO2/(Rgas*T))
            fO2 = QFM   ! pressure of O2 in bars
          ElseIf (buffer.eq.6) Then  ! QFI buffer not working 22 Nov 
            If (T .lt. 1042.) then
              xO2 = BB1(3) +BB1(2)*T +BB1(1)*T*Log(T)
            elseif (T. lt. 1184.) then
              xO2 = BB2(3) +BB2(2)*T +BB2(1)*T*Log(T)
            else
              xO2 = BB3(3) +BB3(2)*T +BB3(1)*T*Log(T)
            endif
            QFI = EXP(xO2/(Rgas*T))
            fO2 = QFI   ! pressure of O2 in bars
          Else
            print *, 'unknown buffer'
            return
          Endif

c   the three equilibria
          K1 = AK1(1)*exp(AK1(2)/T+AK1(3)/T**2)
          K2 = AK2(1)*exp(AK2(2)/T)  ! this is different than in no buffer case because p varies
          K5 = sqrt(AK5(1))*exp(AK5(2)/2./T)*sqrt(fO2)  ! include fO2 in K5

c   these are the constant columns
          N_C = N(LCO) + N(LCO2) + N(LCH4)
          N_H = N(LH2) + N(LH2O) + 2*N(LCH4)
          N_N = 2*N(LN2) + N(LNH3)  ! implemented only for nitrogen

c   linear first guess.  this should be upgraded to include CH4 if possible
c   QFI fails for this reason I think.
          x = N_H/(1 + K5)
          y = x*K5
c   need a first guess for p_t   
          p_t = 0.7*p_tot   ! I do have a p_tot from the input atmosphere

c         print *, T

c  
c  p = sum_i N_i *mass_i *g

          do j=1,15
            N_t = N_C + N(LN2) + x*(1+K5)  ! = sum of Ns
            F1 = (1+K5/K1 + K2*p_t**2*x**2/(K5*N_t**2))*(x*(1+K5)-N_H)
     $       + 2*N_C*K2*p_t**2*x**2/(K5*N_t**2)
            F3 = (28 + 44*K5/K1 + 16*K2*p_t**2*x**2/(K5*N_t**2))
     $          /(1+K5/K1 + K2*p_t**2*x**2/(K5*N_t**2))
     $        + x*(2 + 18*K5) + 28*N(LN2) - p_t*1e6/g

c    Jacobean terms, numerically evaluated
c       store old values
            oldx = x
            oldp = p_t
            oldF1 = F1
            oldF3 = F3
            oldNt = N_t

c       numerically evaluate dF/dx terms
            x = oldx*(1+eps)
            dx = oldx*eps
            N_t = N_C + N(LN2) + x*(1+K5)
            F1 = (1+K5/K1 + K2*p_t**2*x**2/(K5*N_t**2))*(x*(1+K5)-N_H)
     $       + 2*N_C*K2*p_t**2*x**2/(K5*N_t**2)
            F3 = (28 + 44*K5/K1 + 16*K2*p_t**2*x**2/(K5*N_t**2))
     $         /(1+K5/K1 + K2*p_t**2*x**2/(K5*N_t**2))
     $        + x*(2 + 18*K5) + 28*N(LN2) - p_t*1e6/g
            dF1dx = (F1-oldF1)/dx
            dF3dx = (F3-oldF3)/dx

c       evaluate dF/dp terms
            x = oldx   ! restore old x
            p_t = oldp*(1+eps)
            dp = oldp*eps
            N_t = N_C + N(LN2) + x*(1+K5)
            F1 = (1+K5/K1 + K2*p_t**2*x**2/(K5*N_t**2))*(x*(1+K5)-N_H)
     $       + 2*N_C*K2*p_t**2*x**2/(K5*N_t**2)
            F3 = (28 + 44*K5/K1 + 16*K2*p_t**2*x**2/(K5*N_t**2))
     $          /(1+K5/K1 + K2*p_t**2*x**2/(K5*N_t**2))
     $        + x*(2 + 18*K5) + 28*N(LN2) - p_t*1e6/g
            dF1dp = (F1-oldF1)/dp
            dF3dp = (F3-oldF3)/dp

c       determinant
            Det = dF1dx*dF3dp - dF1dp*dF3dx

C       restore old values
            p_t = oldp
            N_t = oldNt
            F1 = oldF1
            F3 = oldF3
c       estimate new x, y  - ahh an error
            deltax = (dF1dp*F3 - dF3dp*F1)/Det
            deltap = (dF3dx*F1 - dF1dx*F3)/Det
            newx = x + deltax
            newp = p_t + deltap

c           print *, j, 'x', x, newx
c           print *, j, 'p', p_t, newp
            
            x = newx
            p_t = newp
          enddo

c      assign columns moles/cm2
          N(LH2)  = x
          N(LH2O) = x*K5
          N(LCO)  = N_C/(1 + K5/K1 + K2*x**2*p_t**2/K5/N_t**2)
          N(LCO2) = N(LCO)*K5/K1
          N(LCH4) = N_C - N(LCO) - N(LCO2)
          N_tot = N_t
          p_tot = p_t  ! update p_tot

c   nitrogen as minor species
          K4 = AK4(1)*exp(AK4(2)/T)
          K44 = K4 * (p_t/N_t)**2 *x**3
          z = sqrt(K44**2/4 + 2*K44*N_N)/2 - K44/4
          N(LNH3) = z
          N(LN2) = 0.5*(N_N - z)

c         assign partial pressures
          do i=1,7
            p(i) = p_tot*N(i)/N_tot
          enddo

c         write(*,101) (ISPEC(i),i=1,7)
c         write(*,102) T, (p(i),i=1,7)

          do i=1,7
            p_j(i) = p(i)  ! return pressures
          enddo
        ENDIF
c  *********************  ********************  *********************

      return
      end
