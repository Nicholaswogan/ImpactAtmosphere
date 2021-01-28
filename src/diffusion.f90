module diffusion

  double precision :: T_trop = 180.d0 ! tropopuase T (K)
  double precision :: P_trop = 0.1d0 ! tropause P (bar)

  double precision :: g = 981.d0 ! cm/s2
  double precision :: n_avo = 6.022d23 ! avo's number
  double precision :: k_boltz = 1.3807d-16 ! boltzmann constant (cgs units)

  double precision :: TT_surf
  double precision :: ppressure
  double precision :: mmubar
  double precision :: ztrop
  double precision :: m_slope

  logical :: ierr = .false.


contains
  subroutine diffuse(PhiHCN, Ts, Ps, mubar, vd_HCN, Kzz, tope, nz,fHCN)
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
    ! double precision, intent(in) :: HCN_alt ! altitude HCN is injected into atmosphere (cm)

    ! output
    double precision, dimension(nz),intent(out) :: fHCN

    ! other
    integer i
    double precision :: bote = 0.d0
    double precision, dimension(nz+1) :: ze
    double precision, dimension(nz+1) :: Pe
    double precision, dimension(nz+1) :: Te
    double precision, dimension(nz+1) :: ne

    double precision, dimension(nz) :: zm
    double precision, dimension(nz) :: Pm
    double precision, dimension(nz) :: Tm
    double precision, dimension(nz) :: nm

    double precision, dimension(nz,nz) :: AAA
    double precision, dimension(nz) :: bbbb, DD
    double precision, dimension(nz-1) :: DDL, DDU
    double precision :: A,B,C
    double precision :: D,D0,D1
    double precision :: E,E1,E2
    integer, dimension(nz) :: IPIV
    integer info
    double precision :: dz ! grid spacing (cm)

    ! make input vars global
    TT_surf = Ts
    ppressure = Ps
    mmubar = mubar

    ztrop = -k_boltz/(mubar/N_avo*g)* &
              (T_trop-TT_surf)*(1.d0/dlog(T_trop/TT_surf)) &
              *dlog(p_trop/ppressure)

    m_slope = ((T_trop-TT_surf)/ztrop)

    dz = tope/nz
    ! edge
    ze(1) = 0.d0
    do i=2,nz+1
      ze(i) = ze(i-1)+dz
    enddo
    do i=1,nz+1
      call pressur(ze(i),Pe(i))
      call temper(ze(i),Te(i))
      call densty(ze(i),ne(i))
    enddo

    ! middle
    zm(1) = 0.5d0*dz
    do i=2,nz
      zm(i) = zm(i-1)+dz
    enddo
    do i=1,nz
      call pressur(zm(i),Pm(i))
      call temper(zm(i),Tm(i))
      call densty(zm(i),nm(i))
    enddo

    DDL = 0.d0
    DD = 0.d0
    DDU = 0.d0

    ! middle of atmosphere
    do i=2,nz-1
      A = Kzz*ne(i+1)/dz**2.d0/nm(i+1)
      B = (-Kzz*ne(i+1)/dz**2.d0-Kzz*ne(i)/dz**2.d0)*(1.d0/nm(i))
      C = Kzz*ne(i)/dz**2.d0/nm(i-1)
      DDL(i-1) = C
      DD(i) = B
      DDU(i) = A
    enddo

    ! lower boundary
    D = Kzz*ne(2)/dz**2.d0
    D0 = -(D/nm(1))-vd_HCN/dz
    D1 = (D/nm(2))
    ! RHS[0] = D0*n[0] + D1*n[1]
    DD(1) = D0
    DDU(1) = D1

    ! upper boundary
    E = Kzz*ne(nz)/dz**2.d0
    E1 = (E/nm(nz-1))
    E2 = -(E/nm(nz))
    ! RHS[-1] = E1*n[-2] + E2*n[-1] + PhiT/dz
    DDL(nz-1) = E1
    DD(nz) = E2

    ! setup the b matrix
    bbbb = 0.d0
    bbbb(nz) = -PhiHCN/dz

    ! solve
    call dgtsv (nz, 1, DDL, DD, DDU, bbbb, nz, INFO)
    if (info.ne.0) then
      print*,'Linear solve failed in subroutine HCN_diffusion'
      ierr = .true.
    endif

    fHCN = bbbb/nm

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

end module
