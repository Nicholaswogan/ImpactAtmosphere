

subroutine HCN_hydrolysis_rate(T, pH, ktot)
  implicit none
  double precision, intent(in) :: T
  double precision, intent(in) :: pH
  double precision, intent(out) :: ktot

  double precision :: H, OH, logk1H, k1H, logk1OH, k1OH
  double precision :: pKw, pka_HCN, kw, ka_HCN

  H=10.d0**(-1.d0*pH)
  OH=1.0d-14/H

  ! acid catalyzed: (in M^-1 s^-1)
  logk1H=-4950.d0/T + 8.43d0
  k1H=10.d0**(logk1H)
  ! base catalyzed: (in M^-1 s^-1)
  logk1OH=-4240.d0/T + 11.1d0
  k1OH=10.d0**(logk1OH)

  pKw=-6.0846d0+4471.33d0/T+.017053d0*T  !from Stribling and Miller 1987, from Miyakawa's original sources (Robinson and Stokes 1959 and Schlesinger and Miller 1973, resp.)
  pKa_HCN=-8.85d0+3802.d0/T+.01786d0*T
  Kw=10.d0**(-1.d0*pKw)
  Ka_HCN=10.d0**(-1.d0*pKa_HCN)

  ktot=k1H*H+(k1OH*Kw/(H+Ka_HCN))  !in s^-1

end subroutine
