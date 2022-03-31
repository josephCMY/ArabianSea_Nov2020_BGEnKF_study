module qn2re_gfdlfv3
  ! wrote according to cloud_diagnosis() in gfdl_cloud_microphys.F90
  !
  ! Note: 
  !   The definition of effective radius used here for rain/snow/graupel 
  ! is not the same as those used in cloud_diagnosis(). Details are stated below. 
  !
  ! References:
  ! Heymsfield, A. J., and G. M. McFarquhar, 1996: High Albedos of Cirrus in the
  !   Tropical Pacific Warm Pool: Microphysical Interpretations from CEPEX and from
  !   Kwajalein, Marshall Islands. J. Atmos. Sci., 53, 2424–2451,
  !   doi:10.1175/1520-0469(1996)053<2424:HAOCIT>2.0.CO;2.
  ! Lin, Y.-L., R. D. Farley, and H. D. Orville, 1983: Bulk Parameterization of
  !   the Snow Field in a Cloud Model. J. Clim. Appl. Meteorol., 22, 1065–1092,
  !   doi:10.1175/1520-0450(1983)022<1065:BPOTSF>2.0.CO;2.

  use qn2re_wsm6, only: PI, PI_OVER_6

  real, parameter :: tice = 273.16

  real, parameter :: rhow = 1.0e3 ! density of water
  real, parameter :: rhor = 1.0e3 !         of rain
  real, parameter :: rhos = 0.1e3 !         of snow
  real, parameter :: rhog = 0.4e3 !         of graupel
  ! n0 in exponential PSD  
  real, parameter :: n0r = 8.0e6 ! for rain
  real, parameter :: n0s = 3.0e6 ! for snow
  real, parameter :: n0g = 4.0e6 ! for graupel 
  
  real, parameter :: qmin = 1.0e-7 
  real, parameter :: ccn  = 1.0e8 ! # concentration cloud water
  real, parameter :: beta = 1.22

  ! The following values are used by Lin et al. (1983) to match the fall
  ! velocity of a particle with radius equal to R_eff of the distribution, to
  ! the mean fall velocity of the distribution. We do not use this definition of
  ! R_eff. These values are not used.
  ! real, parameter :: alphar = 0.8     
  ! real, parameter :: alphas = 0.25 
  ! real, parameter :: alphag = 0.5
  ! real, parameter :: gammar = 17.837789 
  ! real, parameter :: gammas = 8.2850630 
  ! real, parameter :: gammag = 11.631769
  
  ! Min and max effective radius Re given in meters
  real, parameter :: rewmin = 5.0e-6,  rewmax = 10.0e-6    ! cld water
  real, parameter :: reimin = 10.0e-6, reimax = 150.0e-6   ! cld ice
  real, parameter :: rermin = 0.0,     rermax = 10000.0e-6 ! rain
  real, parameter :: resmin = 0.0,     resmax = 10000.0e-6 ! snow
  real, parameter :: regmin = 0.0,     regmax = 10000.0e-6 ! graupel
   
contains
  !======== GFDL MP in FV3 ====================!
  ! ------------------------------------------ !
  ! cloud water (martin et al., 1994)
  ! cloud water is assumed to be mono-dispersed
  ! cloud number concentration ccn is set to be 1.0e8
  ! ------------------------------------------ !
  elemental function qn2re_GFDLFV3_cloud(qw, den) result (Rew)
    real, intent(in)  :: qw  ! cloud water mixing ratio (kg/kg)
    real, intent(in)  :: den ! air density (kg/m^3)
    real :: Rew ! effective radius (m)

    real :: qcw
    if (qw .gt. qmin) then
      qcw = den * qw
      rew = exp (1.0 / 3.0 * log ((3 * den * qw ) / (4 * pi * rhow * ccn))) 
      rew = max (rewmin, min (rewmax, rew ))
    else
      qcw = 0.0
      rew = rewmin
    endif
  end function qn2re_GFDLFV3_cloud
  !======== GFDL MP in FV3 ====================!
  ! ------------------------------------------ !
  ! cloud ice (heymsfield and mcfarquhar, 1996), Eq. 5
  ! ------------------------------------------ !
  ! The following equation for R_eff is derived from their Eq. (5) and the 
  ! relationshp stated  ! in the text after their Eq. (5): A = 1.22*IWC/R_eff,
  ! where IWC in g m^(-3), A in mm^2 L^(-1), R_eff in mm,
  ! and the value 1.22 has unit of 10^(-6) g^(-1) m^3
  !
  ! It turns out that the ice particles have fixed density of 614.8 kg m^(-3)
  elemental function qn2re_GFDLFV3_ice(qi, den, t) result (Rew)
    real, intent(in)  :: qi  ! cloud ice mixing ratio (kg/kg)
    real, intent(in)  :: den ! air density (kg/m^3)
    real, intent(in)  :: t   ! temperature (K)
    real :: Rew ! effective radius (m)

    real :: qci    
    if (qi  .gt. qmin) then
      qci  = den  * qi 
      if (t  - tice .lt. - 50) then
        rei  = beta / 9.917 * exp ((1 - 0.891) * log (1.0e3 * qci )) * 1.0e-3
      elseif (t  - tice .lt. - 40) then
        rei  = beta / 9.337 * exp ((1 - 0.920) * log (1.0e3 * qci )) * 1.0e-3
      elseif (t  - tice .lt. - 30) then
        rei  = beta / 9.208 * exp ((1 - 0.945) * log (1.0e3 * qci )) * 1.0e-3
      else
        rei  = beta / 9.387 * exp ((1 - 0.969) * log (1.0e3 * qci )) * 1.0e-3
      endif
      rei  = max (reimin, min (reimax, rei ))
    else
      qci  = 0.0
      rei  = reimin
    endif
  end function qn2re_GFDLFV3_ice
  !======== GFDL MP in FV3 ====================!
  ! ------------------------------------------ !
  ! rain (lin et al., 1983)
  ! ------------------------------------------ !
  ! Note: in subroutine ‘cloud_diagnosis’ the relationship between R_eff and λ
  ! is different. It matches the fall velocity calculated using their Eq. (7) of
  ! a particle with radius R_eff to the mean fall velocity of the distribution
  ! calculated using their Eq. (11). That is not what we want here.
  elemental function qn2re_GFDLFV3_rain(qr, den) result (Rer)
    real, intent(in)  :: qr  ! rain ice mixing ratio (kg/kg)
    real, intent(in)  :: den ! air density (kg/m^3)
    real :: Rer ! effective radius (m)

    real :: qcr    
            
    if (qr  .gt. qmin) then
      qcr  = den  * qr 
      lambdar = exp (0.25 * log (pi * rhor * n0r / qcr ))
      !rer  = 0.5 * exp (log (gammar / 6) / alphar) / lambdar  ! Lin's approach
      rer = 3./(2*lambdar)
      rer  = max (rermin, min (rermax, rer ))
    else
      qcr  = 0.0
      rer  = rermin
    endif
  end function qn2re_GFDLFV3_rain
            
  !======== GFDL MP in FV3 ====================!
  ! ------------------------------------------ !
  ! snow (lin et al., 1983)
  ! ------------------------------------------ !
  ! Note: in subroutine ‘cloud_diagnosis’ the relationship between R_eff and λ
  ! is different. It matches the fall velocity calculated using their Eq. (8) of
  ! a particle with radius R_eff to the mean fall velocity of the distribution
  ! calculated using their Eq. (12). That is not what we want here.
  elemental function qn2re_GFDLFV3_snow(qs, den) result (Res)
    real, intent(in)  :: qs  ! snow ice mixing ratio (kg/kg)
    real, intent(in)  :: den ! air density (kg/m^3)
    real :: Res ! effective radius (m)

    real :: qcs    
            
    if (qs  .gt. qmin) then
      qcs  = den  * qs 
      lambdas = exp (0.25 * log (pi * rhos * n0s / qcs ))
      !res  = 0.5 * exp (log (gammas / 6) / alphas) / lambdas  ! Lin's approach
      res = 3./(2*lambdas)
      res  = max (resmin, min (resmax, res ))
    else
      qcs  = 0.0
      res  = resmin
    endif
  end function qn2re_GFDLFV3_snow
            
  !======== GFDL MP in FV3 ====================!
  ! ------------------------------------------ !
  ! graupel (lin et al., 1983)
  ! ------------------------------------------ !
  ! Note: in subroutine ‘cloud_diagnosis’ the relationship between R_eff and λ
  ! is different. It matches the fall velocity calculated using their Eq. (9) of
  ! a particle with radius R_eff to the mean fall velocity of the distribution
  ! calculated using their Eq. (13). That is not what we want here.
  elemental function qn2re_GFDLFV3_graup(qg, den) result (Reg)
    real, intent(in)  :: qg  ! graupel ice mixing ratio (kg/kg)
    real, intent(in)  :: den ! air density (kg/m^3)
    real :: Reg ! effective radius (m)

    real :: qcg    
            
    if (qg  .gt. qmin) then
      qcg  = den  * qg 
      lambdag = exp (0.25 * log (pi * rhog * n0g / qcg ))
      !reg  = 0.5 * exp (log (gammag / 6) / alphag) / lambdag  ! Lin's approach
      reg = 3./(2*lambdag)
      reg  = max (regmin, min (regmax, reg ))
    else
      qcg  = 0.0
      reg  = regmin
    endif
  end function qn2re_GFDLFV3_graup
            
end module qn2re_gfdlfv3
