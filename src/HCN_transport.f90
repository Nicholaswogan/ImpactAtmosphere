module diffusion

  ! free parameters
  double precision :: LLL = 1.d0 ! g H2O/m3 of clouds
  double precision :: FFF = 0.05d0 ! fraction of the time it rains
  double precision :: gamma = 4.d5 ! average time of storm cycle (s)

  ! planet specific
  double precision :: g = 981.d0 ! cm/s2

  ! constants
  double precision, parameter :: n_avo = 6.022d23 ! avo's number
  double precision, parameter :: k_boltz = 1.3807d-16 ! boltzmann constant (cgs units)
  double precision, parameter :: c1 = 1.d-6 !dynes/bar
  double precision, parameter :: c2 = 1.d-9 ![m3/cm3][L H2O/g H2O]
  double precision, parameter :: MH2O = 18.d0 ! g H2O/mol H2O
  double precision, parameter :: rho_H2O = 1000.d0 ! g H2O/L H2O

  ! for making inputs global to module
  double precision :: TT_surf
  double precision :: ppressure
  double precision :: mmubar
  double precision :: ztrop
  double precision :: m_slope
  double precision :: T_trop = 180.d0 ! tropopuase T (K)
  double precision :: P_trop = 0.1d0 ! tropause P (bar)

  logical :: rainout_on = .true. ! if true then

contains
  subroutine hcn_transport(PhiHCN, Ts, Ps, mubar, vd_HCN, Kzz, tope, nz, &
                           T_trop1, P_trop1, zm, WHCN, fHCN)
    implicit none

    ! input
    double precision, intent(in) :: PhiHCN ! HCN production rate (molecules/cm2/s)
    double precision, intent(in) :: Ts ! surface T (K)
    double precision, intent(in) :: Ps ! surface P (bar)
    double precision, intent(in) :: mubar ! mean molecular weight (g/mol)
    double precision, intent(in) :: vd_HCN ! HCN deposition velocity (cm/s)
    double precision, intent(in) :: Kzz ! diffusion (cm2/s)
    double precision, intent(in) :: tope ! top of atmosphere (cm)
    integer, intent(in) :: nz ! number of vertical layers
    double precision, intent(in) :: T_trop1 ! top of atmosphere (cm)
    double precision, intent(in) :: P_trop1 ! top of atmosphere (cm)
    ! double precision, intent(in) :: HCN_alt ! altitude HCN is injected into atmosphere (cm)

    ! output
    double precision, dimension(nz),intent(out) :: fHCN
    double precision, dimension(nz),intent(out) :: WHCN
    double precision, dimension(nz), intent(out) :: zm

    ! other
    integer i
    double precision :: bote = 0.d0
    double precision, dimension(nz+1) :: ze
    double precision, dimension(nz+1) :: Pe
    double precision, dimension(nz+1) :: Te
    double precision, dimension(nz+1) :: ne

    double precision, dimension(nz) :: Pm
    double precision, dimension(nz) :: Tm
    double precision, dimension(nz) :: nm

    double precision, dimension(nz,nz) :: AAA
    double precision, dimension(nz) :: bbbb, DD, xc
    double precision, dimension(nz-1) :: DDL, DDU
    double precision :: D,D0,D1
    double precision :: E,E1,E2
    double precision :: dz ! grid spacing (cm)

    ! for rainout stuff
    double precision, dimension(nz) :: wH2O, k_star
    double precision, dimension(nz) :: f_H2Om
    double precision :: f_H2O_below, f_H2O_above, phiB_wh2o, phiT_wh2o
    double precision :: k_bar, Q_i, lnkH, H_HCN
    integer :: ind_trop

    ! make input vars global
    TT_surf = Ts
    ppressure = Ps
    mmubar = mubar
    T_trop = T_trop1
    P_trop = P_trop1

    !!!!!!! set up the atmosphere !!!!!!!

    ! troposphere altitude (cm)
    ztrop = -k_boltz/(mubar/N_avo*g)* &
              (T_trop-TT_surf)*(1.d0/dlog(T_trop/TT_surf)) &
              *dlog(p_trop/ppressure)

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
    ind_trop = minloc(dabs(zm - ztrop), dim=1)-1 ! index of troposphere
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
    ! k star
    k_star = 0.d0
    do i=1,ind_trop
      lnkH = 8205.7d0/max(Tm(i),273.d0)-25.323d0 ! in M/atm
      H_HCN = dexp(lnkH)/1.013d0 ! M/bar
      k_bar = (C1*k_boltz*Tm(i)*H_HCN/(1.d0+C1*C2*n_avo*LLL*k_boltz*Tm(i)*H_HCN)) &
               * (WH2O(i)*MH2O/rho_H2O)
      Q_i = (1.d0-FFF) + (FFF/(gamma*k_bar))*(1.d0 - dexp(-k_bar*gamma))
      k_star(i) = (1.d0/(gamma*Q_i)) * (1.d0 - dexp(-k_bar*gamma))
    enddo
    if (.not. rainout_on) then
      k_star = 0.d0 ! turns off rainout
    endif

    !!!!!!! Set up the linear system of equations !!!!!!!

    ! middle of atmosphere
    DDL = 0.d0
    DD = 0.d0
    DDU = 0.d0
    do i=2,nz-1
      DDU(i) = Kzz*ne(i+1)/dz**2.d0/nm(i)
      DD(i) = (-Kzz*ne(i+1)/dz**2.d0 - Kzz*ne(i)/dz**2.d0)*(1.d0/nm(i)) - k_star(i)
      DDL(i-1) = Kzz*ne(i)/dz**2.d0/nm(i)
    enddo

    ! lower boundary
    D = Kzz*ne(2)/(nm(1)*dz**2.d0)
    D0 = - D - k_star(1) - vd_HCN/dz
    D1 = D
    DD(1) = D0
    DDU(1) = D1

    ! upper boundary
    E = Kzz*ne(nz)/(nm(nz)*dz**2.d0)
    E1 = E
    E2 = - E - k_star(nz)
    DDL(nz-1) = E1
    DD(nz) = E2

    ! setup the b matrix
    bbbb = 0.d0
    bbbb(nz) = - PhiHCN/nm(nz)/dz

    !!!!!!! solve the linear system of equations !!!!!!!
    call tridiag(nz,DDL,DD,DDU,bbbb,xc)
    fHCN = xc
    WHCN = fHCN*nm*k_star

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
      PPP = 1.d6*ppressure*dexp(-mmubar/N_avo*g/k_boltz &
              *(1.d0/m_slope)*dlog((m_slope*z+TT_surf)/TT_surf))
    elseif (z > ztrop) then
      PPP = 1.d6*ppressure*dexp(-(mmubar/N_avo*g/k_boltz &
                    *(1.d0/m_slope)*dlog(T_trop/TT_surf) &
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
    pH2O = 1.0d6*dexp(5000.0d0/373.0d0-5000.0d0/TTT)
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

subroutine tridiag(l,a,b,c,d,xc)
  implicit none
  integer,intent(in) :: l
  double precision, dimension(l), intent(in) :: b,d
  double precision, dimension(l-1), intent(in) :: a,c
  double precision, dimension(l), intent(out) :: xc
  double precision, dimension(l) :: bc,dc
  double precision, dimension(l-1) :: ac,cc
  double precision :: mc
  integer i, ii
  ac = a
  bc = b
  cc = c
  dc = d
  do i=2,l
    mc = ac(i-1)/bc(i-1)
    bc(i) = bc(i) - mc*cc(i-1)
    dc(i) = dc(i) - mc*dc(i-1)
  enddo
  xc = bc
  xc(l) = dc(l)/bc(l)
  do i=1,l-1
    ii = l-i
    xc(ii) = (dc(ii)-cc(ii)*xc(ii+1))/bc(ii)
  enddo
end subroutine
