module qn2re_thompson08
  ! wrote according to Thompson et al. (2008)
  ! DOI: 10.1175/2008MWR2387.1
  ! 

  use qn2re_wsm6, only: PI, PI_OVER_6
  real, parameter, private :: rho_w = 1000 ! kg/m3
  REAL, PARAMETER, PRIVATE :: rho_s = 100.0
  REAL, PARAMETER, PRIVATE :: rho_g = 500.0
  REAL, PARAMETER, PRIVATE :: rho_i = 890.0
contains
  elemental function qn2re_Thompson08_cloud(qw) result (Re)
    real, intent(in) :: qw ! cloud water mixing ratio (kg/kg)
    real :: Re
  end function qn2re_Thompson08_cloud
  ! ======= Rain ========= !
  ! Two moment rain scheme
  ! Assume mu=0 (exp. dist.)
  elemental function qn2re_Thompson08_rain(qr, den, Nr) result (Re)
    real, intent(in) :: qr ! rain mixing ratio (kg/kg)
    real, intent(in) :: den ! air density (kg/m3)
    real, intent(in) :: Nr ! rain number concentration
    real :: Re

    real, parameter :: N1 = 9.0e9 ! m^-4
    real, parameter :: N2 = 2.0e6 ! m^-4
    real, parameter :: qr0 = 1.0e-4 ! kg/kg

    real, parameter :: a = rho_w * PI_OVER_6
    real, parameter :: b = 3.0

    real :: N0r 
    real :: lmda
    
    real, parameter :: mu = 0
    
    ! Eq. A6 in Thompson et al. (2008)
    ! N0r = (N1-N2)*0.5 * tanh( (qr0 - qr) / (4*qr0) ) + (N1+N2)*0.5
    lmda = ( a * nr / (qr*den) * gamma(mu+b+1) / gamma(mu+1) ) ** (1.0/b)
    Re = gamma(3+mu+1) / gamma(2+mu+1) * (1/lmda) * 0.5

  end function qn2re_Thompson08_rain
  ! ======= Graupel ========= !
  elemental function qn2re_Thompson08_graup(qg, den) result (Re)
    real, intent(in) :: qg
    real, intent(in) :: den
    real :: Re

    real :: N0g
    real :: lmda

    N0g = max(1.0e4, min( 200/qg, 5.0e6) )
    lmda = exp (0.25 * log (pi * rho_g * n0g / (den*qg) ))
    re = 3. / (2*lmda)
  end function qn2re_Thompson08_graup

  ! ======= Snow ========= !
  ! Assume m(D) = 0.069 * D**2
  elemental function qn2re_Thompson08_snow(qr, den, T) result(Re)
    real, intent(in) :: qr
    real, intent(in) :: den
    real, intent(in) :: T ! K
    real :: Re

    real :: M2, M3
    real :: Tc

    Tc = T - 273.15
    
    M2 = qr * den / 0.069
    M3 = calc_Mn(m2, 3., Tc)
    Re = (M3/M2)/2
  end function qn2re_Thompson08_snow
  ! helper function
  ! Calculate nth moment of PSD
  ! according to Thompson et al. (2008) Eq. C9-C11
  elemental function calc_Mn(M2, n, Tc) result (Mn)
    real, intent(in) :: M2
    real, intent(in) :: n
    real, intent(in) :: Tc ! deg C
    real Mn

    real :: a, b


!      REAL, DIMENSION(10), PARAMETER, PRIVATE:: &
!      sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
!              0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
!      REAL, DIMENSION(10), PARAMETER, PRIVATE:: &
!      sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
!              0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

!         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
!     &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
!     &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
!     &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
!     &         + sa(10)*cse(3)*cse(3)*cse(3)
!         a_ = 10.0**loga_
!         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
!     &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
!     &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
!     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
!         smoz(k) = a_ * smo2(k)**b_
    a = 10.0**(5.065339         - 0.062659*Tc      - 3.032362*n &
             + 0.029469*Tc*n    - 0.000285*Tc**2   &
             + 0.312550*n**2    + 0.000204*Tc**2*n &
             + 0.003199*Tc*n**2 + 0.000000*Tc**3   &
             - 0.015952*n**3 )
    
    b =        0.476221         - 0.015896*Tc      + 0.165977*n &
             + 0.007468*Tc*n    - 0.000141*Tc**2   &
             + 0.060366*n**2    + 0.000079*Tc**2*n &
             + 0.000594*Tc*n**2 + 0.000000*Tc**3   &
             - 0.003577*n**3
    
    Mn = a * M2**b

  end function calc_Mn






end module
