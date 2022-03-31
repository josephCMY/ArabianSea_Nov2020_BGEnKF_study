module qn2re_wsm6
  ! functions to convert hydrometer 'density' ( q * rho_air) and number 
  ! concentration (N) from microphysical scheme to effective radius 
  ! Effective radius Re= E[D^3]/E[D^2]
  !
  ! Function signiture:
  ! Re = qn2re_XXXX_yyyy( rhoaq, N, T )
  ! where 
  !   XXXX: microphysical scheme name, e.g., WSM6
  !   yyyy: hydrometer type, e.g., cloud/rain/ice/snow/graup
  !   inputs: N or T may not be needed for some schemes/hydrometers 
  !     rhoaq = q * rho_air, q=mixing ratio (kg/kg), 
  !                          rho_air=air density (kg/m^3)
  !     N: number concentration of hydrometer (1/m^3)
  !     T: temperature (K)
  !   outputs:
  !     Re: effective radius
  !   All inputs and outputs are scalar and all subroutines are elemental
  implicit none
  
  real, parameter :: PI = 4.D0*DATAN(1.D0)
  real, parameter :: PI_OVER_6 = PI / 6.
contains
  !=================== WSM6 ===================!
  ! ------------------------------------------ !
  ! WSM6 cloud
  elemental function qn2re_WSM6_cloud(q, rhoa) result (Re)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real :: Re
    real :: rhoaq
    real, parameter :: RHO_CLOUD = 1000. ! kg/m^3
    real, parameter :: PAR1 = 6. / (PI*RHO_CLOUD*3.e8) ! 6 / (PI * rho_c * N_c)
    real, parameter :: PAR2 = 1./3.
    rhoaq = q * rhoa
    Re = (PAR1 * rhoaq) ** PAR2 * 0.5 ! convert diamater to radius
  end function qn2re_WSM6_cloud

  ! ------------------------------------------ !
  ! WSM6 ice: calculate 'packed' Re
  ! When Di is set to 500e-6, particles may have different density
  ! Need to pack information in Di and rhoi as such:
  ! Re(Di>500e-6) = Di/2 + rhoi*1e-6
  ! so that information on Di and rhoi can be easily retrieved
  ! 
  elemental function qn2re_WSM6_ice(q, rhoa) result (Re)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real :: Re
    real :: Ni, Mi, Di, rhoi
    real :: rhoaq
    real, parameter :: N_MAX=1.e6
    real, parameter :: N_MIN=1.e3
    real, parameter :: D_MAX=500.e-6
    real, parameter :: RHO_MAX = 920 ! kg/m^3
    rhoaq = q * rhoa 
    Ni = 5.38e7 * rhoaq**0.75
    if (Ni > N_MAX) then
      Ni = N_MAX
    elseif (Ni < N_MIN) then
      Ni = N_MIN
    end if
    Mi = rhoaq / Ni
    Di = 11.9 * sqrt(Mi) 
    if (Di>D_MAX) Di = D_MAX
    rhoi = Mi / (PI_OVER_6 * Di**3)
    if (rhoi > RHO_MAX) rhoi=RHO_MAX
    ! When Di is set to 500e-6, particles may have different density
    ! Need to pack information in Di and rhoi as such:
    ! Re(Di>500e-6) = Di/2 + rhoi*1e-6
    ! so that information on Di and rhoi can be easily retrieved
    if (Di >= D_MAX) then
      Re = Di * 0.5 + rhoi * 1e-6
    else
      Re = Di * 0.5
    end if
  end function qn2re_WSM6_ice
  ! ------------------------------------------ !
  ! WSM6 ice: calculate rho_ice used in scattering calculation
  ! 
  elemental function qn2re_WSM6_ice_rhoi(q, rhoa) result (rhoi)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real :: rhoi
    real :: Ni, Mi, Di
    real :: rhoaq
    real, parameter :: N_MAX=1.e6
    real, parameter :: N_MIN=1.e3
    real, parameter :: D_MAX=500.e-6
    real, parameter :: RHO_MAX = 920 ! kg/m^3
    rhoaq = q * rhoa
    Ni = 5.38e7 * rhoaq**0.75
    if (Ni > N_MAX) then
      Ni = N_MAX
    elseif (Ni < N_MIN) then
      Ni = N_MIN
    end if
    Mi = rhoaq / Ni
    Di = 11.9 * sqrt(Mi) 
    if (Di>D_MAX) Di = D_MAX
    rhoi = Mi / (PI_OVER_6 * Di**3)
    if (rhoi > RHO_MAX) rhoi=RHO_MAX
  end function qn2re_WSM6_ice_rhoi
  ! ------------------------------------------ !
  ! WSM6 ice: calculate UNPACKED Re used for validation
  ! 
  elemental function qn2re_WSM6_ice_unpacked(q, rhoa) result (Re)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real :: Re 
    real :: Ni, Mi, Di
    real :: rhoaq
    real, parameter :: N_MAX=1.e6
    real, parameter :: N_MIN=1.e3
    real, parameter :: D_MAX=500.e-6
    rhoaq = q * rhoa
    Ni = 5.38e7 * rhoaq**0.75
    if (Ni > N_MAX) then
      Ni = N_MAX
    elseif (Ni < N_MIN) then
      Ni = N_MIN
    end if
    Mi = rhoaq / Ni
    Di = 11.9 * sqrt(Mi) 
    if (Di>D_MAX) Di = D_MAX
    Re = Di * 0.5
  end function qn2re_WSM6_ice_unpacked
  ! ------------------------------------------ !
  ! WSM6 rain
  elemental function qn2re_WSM6_rain(q, rhoa) result (Re)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real :: Re
    real :: lmda
    real :: rhoaq
    real, parameter :: RHO_RAIN=1000.
    real, parameter :: PAR1 = PI * RHO_RAIN * 8e6 ! pi * rho_r * N0_r
    real, parameter :: LMDA_MAX=8e4
    rhoaq = q * rhoa
    lmda = (PAR1 / rhoaq)**0.25
    if (lmda>LMDA_MAX) lmda=LMDA_MAX
    Re = 3. / (2 * lmda ) ! Re = 3/(2*lmda)
  end function qn2re_WSM6_rain
  ! ------------------------------------------ !
  ! WSM6 snow
  elemental function qn2re_WSM6_snow(q, rhoa, T) result (Re)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real, intent(in)  :: T
    real :: Re
    real :: lmda, N0
    real :: rhoaq
    real, parameter :: RHO_SNOW=100.
    real, parameter :: PAR1 = PI * RHO_SNOW ! pi * rho_s (where rho_s=100kg/m^3)
    real, parameter :: N_MAX=1.e11
    real, parameter :: LMDA_MAX=1.e5
    rhoaq = q * rhoa
    N0 = 2.e6 * exp( 0.12*(273.15-T) ) 
    if (N0>N_MAX) N0=N_MAX
    lmda = (PAR1 * N0 / rhoaq)**0.25
    if (lmda>LMDA_MAX) lmda=LMDA_MAX
    Re = 3. / (2 * lmda ) ! Re = 3/(2*lmda)
  end function qn2re_WSM6_snow
  ! ------------------------------------------ !
  ! WSM6 graupel
  elemental function qn2re_WSM6_graup(q, rhoa) result (Re)
    real, intent(in)  :: q
    real, intent(in)  :: rhoa
    real :: Re
    real :: lmda
    real :: rhoaq
    real, parameter :: RHO_GRAUP=500.
    real, parameter :: PAR1 = PI * RHO_GRAUP * 4.e6 ! pi * rho_g * N0_g
    real, parameter :: LMDA_MAX=6e4
    rhoaq = q * rhoa
    lmda = (PAR1 / rhoaq)**0.25
    if (lmda > LMDA_MAX) lmda=LMDA_MAX
    Re = 3. / (2 * lmda ) ! Re = 3/(2*lmda)
  end function qn2re_WSM6_graup
end module qn2re_wsm6





