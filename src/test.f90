program test
  use atmos
  use diffusion
  implicit none
  double precision, dimension(10) :: N
  double precision :: tau_uv

  double precision, dimension(10) :: dNdt
  double precision :: pressure
  double precision :: T_surf
  double precision :: N_H2O

  double precision PhiHCN, Ts, Ps, mubar, vd_HCN, Kzz, dz, tope
  integer, parameter :: nz = 60
  double precision,dimension(nz) :: fHCN, zm
  double precision :: T_trop1, P_trop1, HCN_alt


  dNdt = 0.d0
  pressure = 0.d0
  T_surf = 0.d0
  N_H2O = 0.d0

  N = (/ 1.68675180d+28, 1.53771633d+26, 9.21714507d+22, 6.30403649d+25, &
         6.55363284d+24, 2.18545381d+25, 9.36588326d+22, 1.00000000d+10, &
         1.00000000d+10, 1.00000000d+10 /)
  tau_uv = 100.d0

  call setup

  call rhs(0.d0,N,dNdt)

  ! print*,dNdt

  PhiHCN = 1.d9
  Ts = 280.d0
  Ps = 1.d0
  mubar = 28.d0
  vd_HCN = 5.d-3
  Kzz = 1.d5
  tope = 60.d5
  T_trop1 = 180.d0
  P_trop1 = 0.1d0

  call diffuse(PhiHCN, Ts, Ps, mubar, vd_HCN, Kzz, tope, nz, &
                T_trop1, P_trop1, zm, fHCN)


end program
