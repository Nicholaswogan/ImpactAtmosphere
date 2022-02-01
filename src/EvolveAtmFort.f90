module atmos
  implicit none

  double precision :: area = 5.1d18
  double precision :: g = 982.d0
  integer, parameter :: LH2O = 1
  integer, parameter :: LH2 = 2
  integer, parameter :: LCO = 3
  integer, parameter :: LCO2 = 4
  integer, parameter :: LCH4 = 5
  integer, parameter :: LN2 = 6
  integer, parameter :: LNH3 = 7
  integer, parameter :: LHCN = 8
  integer, parameter :: LC2Hn = 9
  integer, parameter :: LHaze = 10
  double precision, dimension(10) :: mass = (/ 0.018d3,0.002d3,0.028d3,0.044d3,0.016d3, &
                                     & 0.028d3,0.017d3,0.027d3,0.026d3,0.099d3 /)
  double precision :: n_avo = 6.022d23
  integer, parameter :: nsun = 7
  integer, parameter  :: nsp = 7
  integer, parameter  :: nspecies = 10
  double precision :: A_escape = 0.723d0*2d12
  double precision :: B_escape = 0.006d0
  double precision :: Wolf_Toon = 3.0d10
  double precision :: eta = 0.14d0
  double precision :: p_trop = 0.1d0
  double precision :: T_trop = 180.d0
  double precision :: fH2O = 1.d-6
  double precision, dimension(7) :: sun = (/ 30.d0 ,20.d0 ,10.d0 ,10.d0 ,10.d0 ,4.d0 ,2.d0 /)
  double precision :: others3 = 1.d6*1.d4   !  scale height * density
  double precision :: others4 = 1.d6*1.d4
  double precision :: tau_uv_init = 100.d0
  ! defined in setup
  double precision, dimension(8,8) :: k = 0.d0
  double precision, dimension(7,7) :: sigma = 0.d0
  double precision, dimension(7) :: Flux = 0.d0

  ! other
  double precision :: N_dry, mu_dry
  double precision, dimension(Nsp) :: Phi_geo = 0.d0
  double precision, dimension(Nspecies) :: NNN

  ! for errors
  logical :: ierr = .false.

contains

  include "setup.f90"

  subroutine rhs_all(y,tau_uv,dNdt,pressure,T_surf,N_H2O,mu)
    implicit none

    ! input
    double precision, dimension(10), intent(in) :: y
    double precision, dimension(10) :: N
    double precision, intent(in) :: tau_uv

    ! output
    double precision, dimension(10), intent(out) :: dNdt
    double precision, intent(out) :: pressure
    double precision, intent(out) :: T_surf
    double precision, intent(out) :: N_H2O
    double precision, intent(out) :: mu


    ! other
    integer L, i, j
    double precision :: p_surf, T_s, p_H2O_surf
    ! for non-linear solve
    integer, parameter :: nn = 3
    integer, parameter :: mm = 3
    double precision, dimension(nn) :: x_initial
    double precision, dimension(mm) :: fvec
    double precision :: tol = 1.49012d-8
    integer info
    integer, dimension(nn) :: iwa
    integer, parameter :: lwa = 54
    double precision, dimension(lwa) :: wa
    ! other
    double precision :: p_H2O, mubar, p_dry, p_H2O_trop
    double precision :: fraction, NX, N_H2O_strat, N_strat, NH2Ox
    double precision :: N_bar, N_t, BBB, fH2, H2_esc
    double precision, dimension(Nsun) :: tau
    double precision, dimension(Nsp) :: Phi_ion, Phi
    double precision :: sum_ions, CH4_ion_breakup, sumO1D, sumOH
    double precision :: sumN2D, sumO
    double precision :: dCH4dt_ox
    double precision :: N_HCN, NH3_HCN, O1D_CH4, O1D_H2, O1D_OH, O_CH4
    double precision :: O_CO, O_H2, OH_CH4, OH_CO, OH_H2

    dNdt = 0.0d0


    N = y

    !##### water vapor #######
    ! first thing to do is to calculate how much water is in
    ! the atmosphere... sorta correct some stuff
    N_dry = sum(N) - N(LH2O)
    ! average g/mol in dry atmosphere
    mu_dry = (sum(N*mass)-N(LH2O)*mass(LH2O))/N_dry

    ! This little routine self consistently calculates
    ! surface pressure, temperature, and the H2O pressure
    ! at the surface.

    ! initial guess
    p_surf = (mu_dry*N_dry*g/N_avo)
    T_s = T_trop*(p_surf/(p_trop*1.d06))**(eta)
    p_H2O_surf  = 1.d6*dexp(5000.d0/373.d0-5000.d0/T_s )
    x_initial = (/ T_s, P_surf,P_H2O_surf /)

    ! solve non-linear system with MINPACK
    call lmdif1(fcn,mm,nn,x_initial,fvec,tol,info,iwa,wa,lwa,100000)
    if ((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) then
      print*,'Non-linear solver failed in subroutine rhs_all ',info
      ierr = .true.
    endif

    T_surf = x_initial(1)
    pressure = x_initial(2)
    p_H2O = x_initial(3)
    pressure = pressure/1.d6
    p_H2O = dexp(5000.d0/373.d0-5000.d0/T_surf)
    p_dry = pressure-p_H2O
    mu = (p_dry*mu_dry + p_H2O*mass(LH2O))/pressure
    N_H2O = pressure*N_avo*1.d6/mu/g - N_dry

    ! water partial pressure in tropopause
    p_H2O_trop  = fH2O*p_trop

    ! fraction of dry atmosphre above tropopause
    fraction = (p_trop-p_H2O_trop)/p_dry
    NX = N_dry*fraction
    ! calculate statospheric H2O
    N_H2O_strat = p_H2O_trop/p_trop*NX

    ! total stratospheric column??? Not sure of purpose of this
    N_strat  = NX + N_H2O_strat

    ! statospheric H2O for photolysis purposes
    NH2Ox = p_H2O_trop/p_trop*N_dry
    N(LH2O) = NH2Ox

    !##### end water vapor #######

    !##### hydrogen escape
    N_bar = sum(N)-N(LH2)
    mubar = (sum(N(1:Nsp)*mass(1:Nsp))-N(LH2)*mass(LH2))/N_bar
    N_t = sum(N)

    ! escape bit
    BBB = ((mass(LCO2)-mass(LH2))/(mubar-mass(LH2)))**2.d0

    ! units of molecules/s
    fH2 = N(LH2)/N_t
    H2_esc = -A_escape*Sun(2)*1.d0*(fH2/(1.d0+2.3d0*fH2)) &
    /dsqrt(1.d0 + B_escape*(Sun(2)*1.d0)**2.d0*BBB)

    !#### Photolysis #####
    ! this is the optical thickness from the haze
    do L=1,Nsun
      tau(L) = tau_uv*N_t/NX  ! weighted
    enddo
    do L=1,Nsun
      do j=1,Nsp
        tau(L) = tau(L) + sigma(L,j)*N(j)
      enddo
    enddo



    ! ion production
    Phi_ion = 0.d0
    Do j=1,Nsp
      Phi_ion(j) = N(j)*sigma(1,j)*flux(1)/4.d0/tau(1)
    Enddo


    Phi = 0.d0
    !!! CO2 photolysis
    do L=2,6
      Phi(LCO2) = Phi(LCO2) &
                  + Flux(L)/4.d0*sigma(L,LCO2)*N(LCO2)/tau(L)
    enddo
    Phi(LCO2) = Phi(LCO2) + 0.2d0*Phi_ion(LCO2)


    !!! CH4 photolysis
    do L=1,4
      Phi(LCH4) = Phi(LCH4) &
                  + Flux(L)/4.d0*sigma(L,LCH4)*N(LCH4)/tau(L)
    enddo

    ! ion production
    sum_ions = 0.d0
    do j=1,6
      sum_ions = sum_ions + N(j)
    enddo

    CH4_ion_breakup = Phi_ion(LH2)*(N(LCH4)/sum_ions) &
                      + Phi_ion(LN2)*(N(LCH4)/(sum_ions - N(LN2)) &
                      + N(LCO)/(sum_ions -N(LN2))*N(LCH4) &
                         /(sum_ions -N(LN2) -N(LCO)) &
                      + N(LCO)/(sum_ions -N(LN2))*N(LCO2) &
                         /(sum_ions -N(LN2) -N(LCO)) &
                      * N(LCH4)/(N(LH2)+N(LCH4)+N(LH2O)) &
                      + N(LH2O)/(sum_ions - N(LN2)) &
                      * N(LCH4)/(N(LCH4)+N(LH2O))  ) &
                      + Phi_ion(LCO)*( N(LCH4)/(N(LCO2)+N(LCH4)+N(LH2O)) &
                      + N(LCO2)/(N(LCO2)+N(LCH4)+N(LH2O)) &
                      * N(LCH4)/(N(LH2)+N(LCH4)+N(LH2O)) ) &
                      + Phi_ion(LCO2)*N(LCH4) &
                        /(N(LH2)+N(LCH4)+N(LH2O)) &
                      + Phi_ion(LH2O)*N(LCH4)/(N(LH2)+N(LCH4))


    Phi(LCH4) = Phi(LCH4) + CH4_ion_breakup

    !!! H2O photolysis
    do L=2,6
      Phi(LH2O) = Phi(LH2O) &
       + Flux(L)/4.d0*sigma(L,LH2O)*N(LH2O)/tau(L)
    enddo

    !!! N2 photolysis
    do L=2,3
      Phi(LN2) = Phi(LN2) &
       + Flux(L)/4.d0*sigma(L,LN2)*N(LN2)/tau(L)
    enddo

    !!! NH3 photolysis
    do L=1,7
      Phi(LNH3) = Phi(LNH3) &
       + Flux(L)/4.d0*sigma(L,LNH3)*N(LNH3)/tau(L)
    enddo

    !#### end Photolysis #####

    !#### Photochemistry #####
    !!! O1D
    sumO1D = 0.d0
    do j=1,6
      sumO1D = sumO1D + k(2,j)*N(j)
    enddo
    O1D_H2  = k(2,LH2)*N(LH2)/sumO1D
    O1D_CH4 = k(2,LCH4)*N(LCH4)/sumO1D

    !!! O1D > OH
    O1D_OH = ( k(2,LH2)*N(LH2) + k(2,LCH4)*N(LCH4) &
              + 2.d0*k(2,LH2O)*N(LH2O) )/sumO1D

    !!! OH
    sumOH = k(3,6)*others3   ! others3 is a background column
    do j=1,5
      sumOH = sumOH + k(3,j)*N(j)
    enddo
    OH_CH4 = k(3,LCH4)*N(LCH4)/sumOH
    OH_H2  = k(3,LH2)*N(LH2)/sumOH
    OH_CO  = k(3,LCO)*N(LCO)/sumOH

    !!! N > HCN
    sumN2D = 0.d0
    do j=1,6
      sumN2D = sumN2D + k(1,j)*N(j)
    enddo
    N_HCN = (k(1,LCH4)*N(LCH4) + k(1,LCO)*N(LCO) &
            + k(1,LN2)*N(LN2) )/sumN2D

    !!! O + CO > CO2
    sumO = k(4,6)*others4   ! others3 is a background column
    do j=1,5
      sumO = sumO + k(4,j)*N(j)
    enddo
    O_CH4 = k(4,LCH4)*N(LCH4)/sumO
    O_H2  = k(4,LH2)*N(LH2)/sumO
    O_CO  = k(4,LCO)*N(LCO)/sumO

    !#### end Photochemistry #####

    !#### Budgets #####
    ! CH4 oxidiation
    dCH4dt_ox = -Phi(LCO2)*O1D_CH4 -Phi(LH2O)*OH_CH4

    ! NH3 > HCN
    NH3_HCN = (Phi(LCH4)-dCH4dt_ox) &
              / (Phi(LCH4)-dCH4dt_ox + Phi(LH2O) +Phi(LCO2) )

    ! CH4
    dNdt(LCH4) = -Phi(LCH4) + dCH4dt_ox + Phi_geo(LCH4) &
                    - Phi(LNH3)*NH3_HCN


    ! HCN
    dNdt(LHCN) = 2.d0*Phi(LN2)*N_HCN* &
                  (Phi(LCH4)-dCH4dt_ox) &
                  /( Phi(LCH4)-dCH4dt_ox + 2.d0*Phi(LN2)*N_HCN ) &
                  + Phi(LNH3)*NH3_HCN

    ! C2Hn
    dNdt(LC2Hn) = 0.5d0*(Phi(LCH4)-dCH4dt_ox)**3.d0 &
                  / (Phi(LCH4)-dCH4dt_ox + Phi(LH2O) +Phi(LCO2) )**2.d0

    ! Haze
    dNdt(LHaze) = (Phi(LCH4)-dCH4dt_ox) * &
                    ( (Phi(LCH4)-dCH4dt_ox) &
                    /(Phi(LCH4)-dCH4dt_ox + Phi(LH2O) +Phi(LCO2)) )**5.d0

    ! CO2
    dNdt(LCO2) = -Phi(LCO2) + Phi(LCO2)*O1D_OH*OH_CO &
                   + Phi(LH2O)*OH_CO + Phi(LCO2)*(1.d0-O1D_OH)*O_CO &
                   - Phi_geo(LCO2)

    ! CO
    dNdt(LCO)  = Phi(LCO2)-(dCH4dt_ox-Phi(LCH4)) &
                   - dNdt(LHaze) - Phi(LCO2)*O1D_OH*OH_CO &
                   - Phi(LH2O)*OH_CO - Phi(LCO2)*(1.d0-O1D_OH)*O_CO &
                   - dNdt(LHCN) + Phi(LNH3)*NH3_HCN
    ! H2
    dNdt(LH2) = H2_esc - 2.d0*(dCH4dt_ox-Phi(LCH4)) &
                  + 2.d0*dNdt(LCO2) + dNdt(LCO) - dNdt(LHaze) - dNdt(LHCN) &
                  + Phi_geo(LH2) + 3.d0*Phi(LNH3)*NH3_HCN

    ! H2O
    dNdt(LH2O) = 0.0d0 ! assumed constant

    ! NH3
    dNdt(LNH3) = -Phi(LNH3) ! assumed lost

    ! N2
    dNdt(LN2) = - 0.5d0*dNdt(LHCN) + 0.5d0*Phi(LNH3) !*(1.d0-NH3_HCN)

  end subroutine

  subroutine rhs(t,y,dNdt)
    implicit none

    ! input
    double precision, dimension(10), intent(in) :: y
    double precision, dimension(10), target :: N
    double precision, intent(in) :: t

    ! output
    double precision, dimension(10), intent(out) :: dNdt

    ! other
    integer i
    double precision pressure,T_surf,N_H2O,mu
    ! for nonlinear solve
    double precision, dimension(1) :: x, fvec
    double precision :: tol = 1.49012d-8
    integer info
    integer, dimension(1) :: iwa
    integer, parameter :: lwa = 10
    double precision, dimension(lwa) :: wa

    do i=1,nspecies
      N(i) = max(y(i),0.d0)
    enddo
    NNN = N ! global NNN
    x(1) = dlog(tau_uv_init) ! intial condtions
    
    ! solve nonlinear system
    call lmdif2(fcn2,1,1,x,fvec,tol,info,iwa,wa,lwa,10000)
    if ((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) then
      print*,'Non-linear solver failed in subroutine rhs ',info
      ierr = .true.
    endif
    tau_uv_init = dexp(x(1))
    
    call rhs_all(N,dexp(x(1)),dNdt,pressure,T_surf,N_H2O,mu)

  end subroutine

  subroutine rhs_verbose(t,y,dNdt,pressure, T_surf, tau_uv, N_H2O, mu)
    implicit none

    ! input
    double precision, dimension(10), intent(in) :: y
    double precision, dimension(10) :: N
    double precision, intent(in) :: t

    ! output
    double precision, dimension(10), intent(out) :: dNdt
    double precision, intent(out) :: pressure, T_surf, N_H2O, tau_uv, mu

    ! other
    integer i
    ! for nonlinear solve
    double precision, dimension(1) :: x, fvec
    double precision :: tol = 1.49012d-8
    integer info
    integer, dimension(1) :: iwa
    integer, parameter :: lwa = 10
    double precision, dimension(lwa) :: wa

    do i=1,nspecies
      N(i) = dmax1(y(i),0.d0)
    enddo

    NNN = N ! global NNN
    x(1) = dlog(tau_uv_init) ! intial condtions

    ! solve nonlinear system
    call lmdif2(fcn2,1,1,x,fvec,tol,info,iwa,wa,lwa,10000)
    if ((info.ne.1) .and. (info.ne.2) .and. (info.ne.3) .and. (info.ne.4)) then
      print*,'Non-linear solver failed in subroutine rhs_verbose ',info
      ierr = .true.
    endif
    tau_uv_init = dexp(x(1))
    tau_uv = dexp(x(1))

    call rhs_all(N,dexp(x(1)),dNdt,pressure,T_surf,N_H2O,mu)

  end subroutine

  subroutine fcn(m,n,x,fvec,iflag)
    implicit none
    integer m,n,iflag
    double precision, dimension(n) :: x
    double precision, dimension(m) :: fvec
    double precision :: T_s,P_s,P_H2O

    T_s = x(1)
    P_s = x(2)
    P_H2O = x(3)

    fvec(1) = T_s-T_trop*(P_s/(p_trop*1.d6))**eta
    fvec(2) = P_s-((N_dry + P_H2O*N_avo/(((P_s-P_H2O)*mu_dry+P_H2O*mass(LH2O))/P_s)/g) &
              *(((P_s-P_H2O)*mu_dry+P_H2O*mass(LH2O))/P_s)*g/N_avo)
    fvec(3) = P_H2O-1.d6*dexp(5000.d0/373.d0-5000.d0/T_s)
  end subroutine

  subroutine fcn2(m,n,x,fvec,iflag)
    implicit none
    integer m,n,iflag
    double precision, dimension(n) :: x
    double precision, dimension(m) :: fvec
    double precision :: tau_uv
    double precision, dimension(10) :: dNdt

    double precision pressure,T_surf,N_H2O,mu

    ! x is log tau_uv
    tau_uv = dexp(x(1))

    call rhs_all(NNN,tau_uv,dNdt,pressure,T_surf,N_H2O,mu)

    fvec(1) = tau_uv-10.d0*((dNdt(LHaze)+dNdt(LHCN))/Wolf_Toon)**0.8d0

  end subroutine

end module
