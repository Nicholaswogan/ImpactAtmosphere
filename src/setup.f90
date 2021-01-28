subroutine setup
  implicit none
  double precision, dimension(7) :: flux0

  !!!!! reaction rates
  k = 0.d0
  ! N2D reactions
  k(1,LH2) = 2.2d-12
  k(1,LCH4) = 4d-12
  k(1,LH2O) = 5d-11
  k(1,LCO2) = 3.6d-13
  k(1,LCO) = 1.9d-12
  k(1,LN2) = 1.7d-14

  ! O1D reactions
  k(2,LH2) = 1.1d-10
  k(2,LCH4) = 1.5d-10
  k(2,LH2O) = 2.2d-10
  k(2,LCO2) = 7.4d-11
  k(2,LCO) = 7d-12
  k(2,LN2) = 1.8d-11

  ! OH reactions
  k(3,LH2) = 6d-15
  k(3,LCH4) = 6d-15
  k(3,LCO) = 1.2d-13
  k(3,6) = 1.d-11    ! OH + others3
  k(3,7) = 1.2d-10   ! OH + CH2
  k(3,8) = 6d-11     ! OH + CH3

  ! O reactions
  k(4,LH2) = 1d-17
  k(4,LCH4) = 7d-18
  k(4,LCO) = 4d-17
  k(4,6)   = 1.d-11  !  O + others4
  k(4,7) = 1.2d-10   ! O + CH2
  k(4,8) = 1.2d-10   ! O + CH3

  ! N reactions
  k(5,7) = 1.2d-10   ! N + CH2
  k(5,8) = 1.1d-10
  k(7,8) = 7d-11     ! CHn + CHn

  ! H2 ion reactions
  k(8,LH2)   =  2.d-9
  k(8,LCH4)  =  3.8d-9
  k(8,LH2O)  =  7.d-9
  k(8,LCO2)  =  2.4d-9
  k(8,LCO)   =  3.d-9
  k(8,LN2)   =  2.d-9

  !!!!! cross sections
  sigma = 0.d0
  sigma(1,LH2)  = 2.d-18
  sigma(1,LCH4) = 2.d-17
  sigma(1,LH2O) = 1.3d-17
  sigma(1,LCO2) = 2.d-17
  sigma(1,LCO)  = 1.3d-17
  sigma(1,LN2)  = 1.3d-17

  sigma(2,LH2)  = 1.d-19
  sigma(2,LCH4) = 8.d-18
  sigma(2,LH2O) = 8.d-18
  sigma(2,LCO2) = 4.d-17
  sigma(2,LCO)  = 1.d-19
  sigma(2,LN2)  = 1.d-19

  sigma(3,LH2)   = 1.d-19
  sigma(3,LCH4)  = 8.d-18
  sigma(3,LH2O)  = 8.d-18
  sigma(3,LCO2)  = 4.d-17
  sigma(3,LCO)   = 1.d-19
  sigma(3,LN2)   = 3.d-16

  ! Lyman alpha
  sigma(4,LCH4)  = 8.d-18
  sigma(4,LH2O)  = 8.d-18
  sigma(4,LCO2)  = 5.d-20
  sigma(4,LNH3)  = 8.d-18   ! made up

  sigma(5,LCH4)  = 0.d-18   ! missing
  sigma(5,LH2O)  = 8.d-19
  sigma(5,LCO2)  = 8.d-18
  sigma(5,LNH3)  = 8.d-19   ! made up

  sigma(6,LH2O)  = 2.d-18
  sigma(6,LCO2)  = 3.d-20
  sigma(6,LNH3)  = 8.d-19   ! made up

  sigma(7,LNH3)  = 3.d-18

  !!!!! solar flux
  Flux0 = 0.0d0
  ! This is the modern sun
  Flux0(1) = 3.d10     ! FXUV
  Flux0(2) = 2.d10     ! FEUV
  Flux0(3) = 3.6d9    ! Lyman gamma +
  Flux0(4) = 3.6d11   ! Lyman alpha
  Flux0(5) = 4.d11     ! FUV_CO2
  Flux0(6) = 2.5d12   ! FUV_H2O
  Flux0(7) = 6.d13     ! FUV_NH3

  flux = flux0*sun ! converts to ancient sun

end subroutine
