program test
  use atmos
  use diffusion
  implicit none

  ! for photochemistry
  double precision, dimension(10) :: N
  double precision :: tau_uv
  double precision, dimension(10) :: dNdt
  double precision :: pressure
  double precision :: T_surf
  double precision :: N_H2O

  ! for HCN transport
  double precision PhiHCN, Ts, Ps, mubar, vd_HCN, Kzz, ocean_pH
  integer, parameter :: nz = 200
  double precision,dimension(nz) :: fHCN, zm,nm, Tm,WHCN,wh2o
  double precision rainout, ocean2atmos, mHCN_s, mHCN_d


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


  PhiHCN = 1.d8
  Ts = 298.d0
  Ps = 1.d0
  mubar = 28.d0
  Kzz = 1.d5
  ocean_pH = 7.d0

  call hcn_transport(PhiHCN, Ts, Ps, mubar, Kzz, ocean_pH,nz, &
                     zm, nm, Tm, fHCN, mHCN_s, mHCN_d, WHCN,wh2o, rainout, ocean2atmos)

  ! print*,fhcn

end program
