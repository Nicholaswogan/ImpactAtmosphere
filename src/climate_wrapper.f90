module climate_wrapper
  implicit none
             
contains
  
  subroutine pahlevan_H2_clima(P_surf, mubar, grav, P_top, nz, &
                               zout, Pout, Tout, ztrop, err)
    use climate, only: pahlevan_H2_clima_ => pahlevan_H2_clima
  
    real(8), intent(in) :: P_surf
    real(8), intent(in) :: mubar, grav, P_top
    integer, intent(in) :: nz
  
    real(8), intent(out) :: zout(nz), Pout(nz), Tout(nz), ztrop
  
    character(len=1000), intent(out) :: err
  
    call pahlevan_H2_clima_(P_surf, mubar, grav, P_top, nz, &
                            zout, Pout, Tout, ztrop, err)
  
  end subroutine
  
  
end module
