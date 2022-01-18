module diffusion
  implicit none

  ! free parameters
  double precision :: tope = 60.d5 ! top of atmosphere (cm)
  double precision :: T_trop = 180.d0 ! tropopuase T (K)
  double precision :: P_trop = 0.1d0 ! tropause P (bar)
  double precision :: LLL = 1.d0 ! g H2O/m3 of clouds
  double precision :: FFF = 0.05d0 ! fraction of the time it rains
  double precision :: gamma = 4.d5 ! average time of storm cycle (s)
  double precision :: rain_scaling = 1.d0 ! scale rain rate
  double precision :: vo = 1.2d-5 ! The turnover velocity of the ocean (cm/s)
  double precision :: zs = 100.d2 ! The depth of the surface ocean (cm)
  double precision :: zd = 4000.d2 ! The depth of the deep ocean (cm)

  ! planet specific
  double precision :: g = 981.d0 ! cm/s2

  ! constants
  double precision, parameter :: n_avo = 6.022d23 ! avo's number
  double precision, parameter :: k_boltz = 1.3807d-16 ! boltzmann constant (cgs units)
  double precision, parameter :: c1 = 1.d-6 !dynes/bar
  double precision, parameter :: c2 = 1.d-9 ![m3/cm3][L H2O/g H2O]
  double precision, parameter :: MH2O = 18.d0 ! g H2O/mol H2O
  double precision, parameter :: rho_H2O = 1000.d0 ! g H2O/L H2O
  double precision, parameter :: CC = n_avo/rho_H2O ! (molecules L)/(cm3 mol)
  double precision, parameter :: vp_HCN = 0.005d0 ! HCN piston velocity cm/s

  ! for making inputs global to module
  double precision :: TT_surf
  double precision :: ppressure
  double precision :: mmubar
  ! other global
  double precision :: ztrop
  double precision :: m_slope

  logical :: rainout_on = .true. ! if true then let it rain
  logical :: ierr
  double precision :: rtol = 1.d-5 ! Require HCN conservation to within this factor

contains
  subroutine hcn_transport(PhiHCN, Ts, Ps, mubar, Kzz, ocean_pH, nz, &
                           zm, nm, Tm, fHCN, mHCN_s, mHCN_d, WHCN, wH2O, rainout, ocean2atmos)
    implicit none

    ! input
    double precision, intent(in) :: PhiHCN ! HCN production rate (molecules/cm2/s)
    double precision, intent(in) :: Ts ! surface T (K)
    double precision, intent(in) :: Ps ! surface P (bar)
    double precision, intent(in) :: mubar ! mean molecular weight (g/mol)
    double precision, intent(in) :: Kzz ! diffusion (cm2/s)
    double precision, intent(in) :: ocean_pH ! ocean pH
    integer, intent(in) :: nz ! number of vertical layers

    ! output
    double precision, dimension(nz), intent(out) :: zm, nm, Tm
    double precision, dimension(nz),intent(out) :: fHCN
    double precision, dimension(nz),intent(out) :: WHCN, wH2O
    double precision, intent(out) :: mHCN_s, mHCN_d, rainout, ocean2atmos

    ! local
    integer i
    double precision :: bote = 0.d0
    double precision, dimension(nz+1) :: ze
    double precision, dimension(nz+1) :: Pe
    double precision, dimension(nz+1) :: Te
    double precision, dimension(nz+1) :: ne
    double precision, dimension(nz) :: Pm
    ! double precision, dimension(nz) :: Tm
    ! double precision, dimension(nz) :: nm
    double precision :: dz ! grid spacing (cm)
    ! matrix stuff
    double precision, dimension(nz+2,nz+2) :: AAA
    double precision, dimension(nz+2) :: bbbb
    integer, dimension(nz+2) :: ipiv
    ! for rainout stuff
    double precision, dimension(nz) :: k_star
    double precision, dimension(nz) :: f_H2Om
    double precision :: f_H2O_below, f_H2O_above, phiB_wh2o, phiT_wh2o
    double precision :: k_bar, Q_i, lnkH, H_HCN
    integer :: ind_trop, info
    ! for ocean
    double precision :: alpha_HCN, ktot, HCN_destroyed

    ! start without errors
    ierr = .false.

    ! make input vars global
    TT_surf = Ts
    ppressure = Ps
    mmubar = mubar

    !!!!!!! set up the atmosphere !!!!!!!

    ! tropopause altitude (cm)
    ztrop = -k_boltz/(mubar/N_avo*g)* &
              (T_trop-TT_surf)*(1.d0/log(T_trop/TT_surf)) &
              *log(p_trop/ppressure)

    ! slope of T(z) from ground to tropopause
    m_slope = ((T_trop-TT_surf)/ztrop)

    dz = tope/nz
    ! edge of grid
    ze(1) = 0.d0
    do i=2,nz+1
      ze(i) = ze(i-1)+dz
    enddo
    do i=1,nz+1
      call pressur(ze(i),Pe(i))
      call temper(ze(i),Te(i))
      call densty(ze(i),ne(i))
    enddo

    ! middle of grid
    zm(1) = 0.5d0*dz
    do i=2,nz
      zm(i) = zm(i-1)+dz
    enddo
    do i=1,nz
      call pressur(zm(i),Pm(i))
      call temper(zm(i),Tm(i))
      call densty(zm(i),nm(i))
    enddo
    ! water vapor in troposphere
    ind_trop = minloc(abs(zm - ztrop), dim=1)-1 ! index of troposphere
    if (zm(ind_trop)>ztrop) ind_trop = ind_trop-1 ! make it so that ind_trop is just below ztrop

    f_H2Om = 0.d0
    do i=1,ind_trop
      call f_H2O_saturation(zm(i),f_H2Om(i))
    enddo

    !!!!!!! calculate raining rate !!!!!!!

    ! middle of atmosphere
    wH2O = 0.d0
    do i=2,ind_trop-1
      wH2O(i) = (Kzz*ne(i+1)/dz**2.d0) * f_H2Om(i+1) &
              - (Kzz*ne(i+1)/dz**2.d0 + Kzz*ne(i)/dz**2.d0) * f_H2Om(i) &
              + (Kzz*ne(i)/dz**2.d0) * f_H2Om(i-1)
    enddo
    ! lower boundary
    call f_H2O_saturation(zm(1)-dz,f_H2O_below)
    ! f_H2O_below = ((f_H2Om(2)-f_H2Om(1))/dz)*(zm(1)-dz) &
                  ! + (f_H2Om(1)-((f_H2Om(2)-f_H2Om(1))/dz)*zm(1))
    phiB_wh2o = - Kzz*ne(1)*(f_H2Om(1) - f_H2O_below)/dz
    wH2O(1) = + (Kzz*ne(2)/dz**2.d0) * f_H2Om(2) &
              - (Kzz*ne(2)/dz**2.d0) * f_H2Om(1) &
              + phiB_wh2o/dz
    ! upper boundary. use linear extrapolation of f_H2Om
    ! call f_H2O_saturation(zm(ind_trop)+dz,f_H2O_above)
    f_H2O_above = ((f_H2Om(ind_trop)-f_H2Om(ind_trop-1))/dz)*(zm(ind_trop)+dz) &
                  + (f_H2Om(ind_trop)-((f_H2Om(ind_trop)-f_H2Om(ind_trop-1))/dz)*zm(ind_trop))
    phiT_wh2o = - Kzz*ne(ind_trop)*(f_H2O_above-f_H2Om(ind_trop))/dz
    wH2O(ind_trop) = - (Kzz*ne(ind_trop-1)/dz**2.d0)*f_H2Om(ind_trop) &
                     + (Kzz*ne(ind_trop-1)/dz**2.d0)*f_H2Om(ind_trop-1) &
                     - phiT_wh2o/dz
    wH2O = wH2O * rain_scaling
    ! k star
    k_star = 0.d0
    do i=1,ind_trop
      lnkH = 8205.7d0/max(Tm(i),273.d0)-25.323d0 ! in M/atm
      H_HCN = exp(lnkH)/1.013d0 ! M/bar
      k_bar = (C1*k_boltz*Tm(i)*H_HCN/(1.d0+C1*C2*n_avo*LLL*k_boltz*Tm(i)*H_HCN)) &
               * (WH2O(i)*MH2O/rho_H2O)
      Q_i = (1.d0-FFF) + (FFF/(gamma*k_bar))*(1.d0 - exp(-k_bar*gamma))
      k_star(i) = (1.d0/(gamma*Q_i)) * (1.d0 - exp(-k_bar*gamma))
    enddo
    if (.not. rainout_on) then
      k_star = 0.d0 ! turns off rainout
    endif

    !!!!!!! Set up the linear system of equations !!!!!!!

    AAA = 0.d0 ! the matrix

    ! The deep ocean
    call HCN_hydrolysis_rate(Te(1), ocean_pH, ktot)
    AAA(1,1) = -ktot - vo/zd
    AAA(1,2) = vo/zd

    ! The surface ocean
    lnkH = 8205.7d0/max(Ts,273.d0)-25.323d0 ! in M/atm
    alpha_HCN = exp(lnkH)/1.013d0 ! M/bar
    AAA(2,1) = vo/zs
    AAA(2,2) = -vp_HCN/zs - ktot - vo/zs
    AAA(2,3) = vp_HCN*alpha_HCN*Ps/zs + k_star(1)*nm(1)*dz/(CC*zs)
    do i = 2,ind_trop
      AAA(2,i+2) = k_star(i)*nm(i)*dz/(CC*zs) ! rain falling in ocean
    enddo

    ! First layer of atmosphere
    AAA(3,2) = vp_HCN*CC/(nm(1)*dz)
    AAA(3,3) = -Kzz*ne(2)/(nm(1)*dz**2.d0) - k_star(1) - vp_HCN*alpha_HCN*Ps*CC/(nm(1)*dz)
    AAA(3,4) = Kzz*ne(2)/(nm(1)*dz**2.d0)

    ! middle of atmosphere
    do i=2,nz-1
      AAA(i+2,i+3) = Kzz*ne(i+1)/dz**2.d0/nm(i) ! upper diagonal
      AAA(i+2,i+2) = (-Kzz*ne(i+1)/dz**2.d0 - Kzz*ne(i)/dz**2.d0)*(1.d0/nm(i)) - k_star(i) ! diagonal
      AAA(i+2,i+1) = Kzz*ne(i)/dz**2.d0/nm(i)! lower diagonal
    enddo

    ! upper boundary
    AAA(nz+2,nz+1) = Kzz*ne(nz)/(nm(nz)*dz**2.d0)
    AAA(nz+2,nz+2) = - Kzz*ne(nz)/(nm(nz)*dz**2.d0) - k_star(nz)

    ! setup the b matrix
    bbbb = 0.d0
    bbbb(nz+2) = - PhiHCN/nm(nz)/dz

    !!!!!!! solve the linear system of equations !!!!!!!
    call dgesv (nz+2, 1, AAA, nz+2, ipiv, bbbb, nz+2, info)
    if (info .ne. 0) then
      print*,'Linear solve failed in subroutine hcn_transport'
      ierr = .true.
    endif

    mHCN_d = bbbb(1)
    mHCN_s = bbbb(2)
    fHCN = bbbb(3:)
    WHCN = fHCN*nm*k_star ! HCN rainout rate vs altitude
    rainout = -sum(WHCN)*(zm(2)-zm(1)) ! HCN flux into the ocean from rain
    ocean2atmos = -vp_HCN*(alpha_HCN*Ps*fHCN(1)-mHCN_s)*CC ! HCN flux over surface boundary layer of ocean

    ! check for conservation of molecules
    HCN_destroyed = ktot*mHCN_d*CC*zd + ktot*mHCN_s*CC*zs
    if ((HCN_destroyed > PhiHCN + PhiHCN*rtol) .or. (HCN_destroyed < PhiHCN - PhiHCN*rtol))  then
      print*,'Warning: Poor conservation of HCN molecules.'
      print*,'HCN input = ',PhiHCN,' (molecules/cm2/s)'
      print*,'HCN destruction = ',HCN_destroyed,'(molecules/cm2/s)'
      ierr = .true.
    endif

  end subroutine

  subroutine temper(z,TTT)
    implicit none
    double precision, intent(in) :: z
    double precision, intent(out) :: TTT
    if (z <= ztrop) then
      TTT = m_slope*z+TT_surf
    elseif (z > ztrop) then
      TTT = T_trop
    endif
  end subroutine

  subroutine pressur(z,PPP)
    implicit none
    double precision, intent(in) :: z
    double precision, intent(out) :: PPP
    if (z <= ztrop) then
      PPP = 1.d6*ppressure*exp(-mmubar/N_avo*g/k_boltz &
              *(1.d0/m_slope)*log((m_slope*z+TT_surf)/TT_surf))
    elseif (z > ztrop) then
      PPP = 1.d6*ppressure*exp(-(mmubar/N_avo*g/k_boltz &
                    *(1.d0/m_slope)*log(T_trop/TT_surf) &
                    +mmubar/N_avo*g/(k_boltz*T_trop)*(z-ztrop)))
    endif
  end subroutine

  subroutine densty(z,dsty)
    implicit none
    double precision, intent(in) :: z
    double precision, intent(out) :: dsty
    double precision pres, temp
    call pressur(z,pres)
    call temper(z,temp)
    dsty = pres/(k_boltz*temp)
  end subroutine

  ! dynes from Zahnle (water saturation)
  subroutine p_H2O_saturation(z,pH2O)
    implicit none
    double precision, intent(in) :: z
    double precision, intent(out) :: pH2O
    double precision :: TTT
    call temper(z,TTT)
    pH2O = 1.0d6*exp(5000.0d0/373.0d0-5000.0d0/TTT)
  end subroutine

  subroutine f_H2O_saturation(z,fH2O)
    implicit none
    double precision, intent(in) :: z
    double precision, intent(out) :: fH2O
    double precision :: pH2O, PPP
    call p_H2O_saturation(z,pH2O)
    call pressur(z,PPP)
    fH2O = pH2O/PPP
  end subroutine

end module
