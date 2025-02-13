  SUBROUTINE Load_CRTM_Structures_MW(atm, sfc, x1,y1,zmax,debug_atm_data)

   USE CRTM_Module
   use qn2re_wsm6
   use qn2re_gfdlfv3
   use qn2re_thompson08
   
   type(CRTM_Atmosphere_type), intent(in out), dimension(1) :: atm
   type(CRTM_Surface_type), intent(in out), dimension(1) :: sfc

   integer, intent(in) :: x1, y1,  zmax

   ! number of layers of ice to modify
   !integer, intent(in) :: N_CLOUDS
   logical, intent(in) :: debug_atm_data
   real :: T_air, rho_air

!*******************************************************************************
! satellite information
!*******************************************************************************

!  sat_dis = sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
!  SCAN_ANGLE = 180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))**2))
!  ZENITH_ANGLE = SCAN_ANGLE+180.0/3.14159*acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
!write(*,*) xlong(x,y),lon(x,y)
!write(*,*) lon(x,y)-sat_lon
!write(*,*) cos(lon(x,y)-sat_lon)
!write(*,*) acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
!write(*,*) 'SCAN_ANGLE',SCAN_ANGLE
!write(*,*) 'ZENITH_ANGLE',ZENITH_ANGLE

!*******************************************************************************
! load WRF data into CRTM structures
!*******************************************************************************


!--- calcurating delz
  do z=1,kx
   if(z.eq.1) then
    delz(z) = (PH(x,y,z+1) + PHB(x,y,z+1)) / 9.806 - hgt(x,y)
   else
    ! delz(z) = ((PH(x,y,z+1) + PHB(x,y,z+1))-(PH(x,y,z) + PHB(x,y,z)))/2/9.806  ! corrected 9 July 2020 after group meeting
    delz(z) = ((PH(x,y,z+1) + PHB(x,y,z+1))-(PH(x,y,z) + PHB(x,y,z)))/9.806
   endif
  enddo
  if (delz(1) <= 0.) delz(1) = delz(2)

!---Writing Atmospheric Profile
   atm(1)%Climatology         = TROPICAL
   atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
   atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
   atm(1)%Level_Pressure(0) = (pres(x,y,zmax)*3.0/2.0 - pres(x,y,zmax-1)/2.0) / 100.0  ! convert from Pa to hPA
!   atm(1)%Level_Pressure(0) = 50.0

   do z=zmax,1,-1
     if(z.eq.1) then
       atm(1)%Level_Pressure(zmax-z+1) = max(psfc(x,y), pres(x,y,1)*3.0/2.0 - pres(x,y,2)/2.0) / 100.0 
     else
       atm(1)%Level_Pressure(zmax-z+1) = ((pres(x,y,z-1) + pres(x,y,z)) / 2.0) / 100.0  ! convert from Pa to hPA
     endif
     atm(1)%Pressure(zmax-z+1)       = pres(x,y,z) / 100.0
     atm(1)%Temperature(zmax-z+1)    = tk(x,y,z)
     atm(1)%Absorber(zmax-z+1,1)     = qvapor(x,y,z)*1000.0
   enddo
! -- handle negative top-most level
    if (atm(1)%Level_Pressure(1) < 0) then
      atm(1)%Level_Pressure(1) = atm(1)%Level_Pressure(2)/2
    end if
    if (atm(1)%Level_Pressure(0) < 0) then
      atm(1)%Level_Pressure(0) = atm(1)%Level_Pressure(1)/2
    end if
    if (atm(1)%Pressure(1) < 0) then
      atm(1)%Pressure(1) = ( atm(1)%Level_Pressure(0) + atm(1)%Level_Pressure(1) )/2
    end if

!    atm(1)%Absorber(:,2) = &
!    (/2.87E+00,2.10E+00,1.45E+00,1.09E+00,1.01E+00,1.07E+00,1.11E+00,8.60E-01,4.45E-01,2.53E-01, &
!      9.39E-02,4.57E-02,6.40E-02,8.40E-02,8.27E-02,7.96E-02,6.70E-02,5.78E-02,4.92E-02,4.48E-02, &
!      3.77E-02,3.01E-02,9.08E-03,2.98E-03,5.12E-03,5.12E-03,1.16E-02,1.16E-02,1.16E-02,1.16E-02, &
!      1.16E-02,1.16E-02,1.16E-02,1.16E-02/)
   atm(1)%Absorber(:,2) = 0.0 ! 5.0E-02!&
!    (/1.26E+00, 5.55E-01, 3.24E-01, 1.07E-01, 7.03E-02, 5.87E-02, 6.15E-02, 6.43E-02, 6.99E-02, 7.17E-02,&
!      7.27E-02, 7.35E-02, 7.38E-02, 7.41E-02, 7.42E-02, 7.41E-02, 7.35E-02, 7.31E-02, 7.27E-02, 7.27E-02,&
!      7.27E-02, 7.26E-02, 7.17E-02, 7.05E-02, 6.80E-02, 6.73E-02, 6.73E-02, 6.76E-02, 6.72E-02, 6.62E-02,&
!      6.51E-02, 6.45E-02, 6.44E-02, 6.46E-02, 6.48E-02, 6.49E-02, 6.46E-02, 6.42E-02, 6.38E-02, 6.38E-02,&
!      6.42E-02, 6.48E-02, 6.56E-02, 6.64E-02, 6.64E-02, 6.72E-02, 6.84E-02, 6.84E-02, 6.84E-02, 6.94E-02,&
!      6.94E-02, 6.72E-02, 6.72E-02, 6.72E-02, 6.05E-02, 6.05E-02, 6.05E-02, 4.12E-02, 4.12E-02, 4.12E-02/)


!--- writing Cloud data
!


do z=1,N_CLOUDS
atm(1)%Cloud(z)%Type = 0
atm(1)%Cloud(z)%Water_Content = 0.0
atm(1)%Cloud(z)%Effective_Radius = 0.0
enddo

   icl = 0
!---calculating # of clouds (cloud and rain)
   IF ( atm(1)%n_Clouds > 0 ) THEN
!---Data for cloud

   icl = 1
   atm(1)%Cloud(icl)%Type = WATER_CLOUD
   atm(1)%Cloud(icl)%Effective_Radius(:) = 0
   atm(1)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qcloud(x,y,z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))
!       if (l_use_default_re) then
         atm(1)%Cloud(icl)%Effective_Radius(k1) =  16.8_fp   ! water cloud effective radius does not matter to the CRTM
!       else
!         atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_cloud(qcloud(x,y,z), rho_air) *1e6 !  16.8_fp * radius_ratio
!       end if
       atm(1)%Cloud(icl)%Water_Content(k1)    = &
           qcloud(x,y,z)* rho_air *delz(z)
     endif
   enddo

   icl = 2
   atm(1)%Cloud(icl)%Type = RAIN_CLOUD
   atm(1)%Cloud(icl)%Effective_Radius(:) = 0
   atm(1)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qrain(x,y,z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))
       if (l_use_default_re) then
         atm(1)%Cloud(icl)%Effective_Radius(k1) =  1000.0_fp 
       else
         if (mp_scheme == 0) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_rain(qrain(x,y,z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 1) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_GFDLFV3_rain(qrain(x,y,z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 2) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_Thompson08_rain(qrain(x,y,z), rho_air, nrain(x,y,z)) *1e6 ! 1000.0_fp * radius_ratio
         end if
       end if
       atm(1)%Cloud(icl)%Water_Content(k1)    = &
           qrain(x,y,z)* rho_air *delz(z)

       !write (*,*) x, y, z, atm(1)%Cloud(icl)%Effective_Radius(k1)
     endif
   enddo

   icl = 3
   atm(1)%Cloud(icl)%Type = ICE_CLOUD
   atm(1)%Cloud(icl)%Effective_Radius(:) = 0
   atm(1)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qice(x,y,z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))
!       if (l_use_default_re) then
         atm(1)%Cloud(icl)%Effective_Radius(k1) =  25.0_fp   ! cloud ice effective radius does not matter to the CRTM
!       else
!         atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_ice(qice(x,y,z), rho_air) *1e6 ! 125.0_fp * radius_ratio
!       end if
       
       atm(1)%Cloud(icl)%Water_Content(k1)    = & 
           qice(x,y,z)* rho_air *delz(z)
     endif
   enddo

   icl = 4
   atm(1)%Cloud(icl)%Type = SNOW_CLOUD
   atm(1)%Cloud(icl)%Effective_Radius(:) = 0
   atm(1)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qsnow(x,y,z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))
       if (l_use_default_re) then
         atm(1)%Cloud(icl)%Effective_Radius(k1) = 500.0_fp 
       else
         T_air = atm(1)%Temperature(k1)
         if (mp_scheme == 0) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_snow(qsnow(x,y,z), rho_air, T_air)  *1e6!500.0_fp * radius_ratio
         else if (mp_scheme == 1) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_GFDLFV3_snow(qsnow(x,y,z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 2) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_Thompson08_snow(qsnow(x,y,z), rho_air, T_air) *1e6 ! 1000.0_fp * radius_ratio
         end if
       end if
       atm(1)%Cloud(icl)%Water_Content(k1)    = & 
           qsnow(x,y,z)* rho_air *delz(z)
     endif
   enddo

   icl = 5
   atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
   atm(1)%Cloud(icl)%Effective_Radius(:) = 0
   atm(1)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qgraup(x,y,z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))
       if (l_use_default_re) then
         atm(1)%Cloud(icl)%Effective_Radius(k1) = 1500.0_fp 
       else
         if (mp_scheme == 0) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_graup(qgraup(x,y,z), rho_air)  *1e6!500.0_fp * radius_ratio
         else if (mp_scheme == 1) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_GFDLFV3_graup(qgraup(x,y,z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 2) then
           atm(1)%Cloud(icl)%Effective_Radius(k1) = qn2re_Thompson08_graup(qgraup(x,y,z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         end if
       end if
       atm(1)%Cloud(icl)%Water_Content(k1)    = & 
           qgraup(x,y,z)* rho_air *delz(z)
     endif
   enddo
   ENDIF

   ! remove negative water content
   do icl = 1, 5
     where (atm(1)%Cloud(icl)%Water_Content < 0) atm(1)%Cloud(icl)%Water_Content = 0.
   end do
!   print *, maxval(atm(1)%Cloud(1)%Effective_Radius(:)) &
!          , maxval(atm(1)%Cloud(2)%Effective_Radius(:)) &
!          , maxval(atm(1)%Cloud(3)%Effective_Radius(:)) &
!          , maxval(atm(1)%Cloud(4)%Effective_Radius(:)) &
!          , maxval(atm(1)%Cloud(5)%Effective_Radius(:)) 
!*******************************************************************************
! Surface data
!*************
   sfc(1)%Wind_speed     = windspeed(x,y)
   sfc(1)%Wind_direction = winddir(x,y)
   if(landmask(x,y).eq.1.0) then
    sfc(1)%Water_Coverage = 0.0_fp
    sfc(1)%Land_Coverage = 1.0_fp
    sfc(1)%Land_Temperature = tsk(x,y)
    sfc(1)%Soil_Temperature = tsk(x,y)
   else
    sfc(1)%Water_Coverage = 1.0_fp
    sfc(1)%Land_Coverage = 0.0_fp
    sfc(1)%Water_Type = 1  ! Sea water
    sfc(1)%Water_Temperature = tsk(x,y)
   endif


  END SUBROUTINE load_crtm_structures_mw

