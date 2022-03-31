program test_nc
  use lut_nc2crtm
  use netcdf  
  use CloudCoeff_Define, only: CloudCoeff_type
  
  character(len=*), parameter :: LUT_DIR = &
                   '/home/meteo/yxl232/wave/tmp/LUT/Kuo_fit/'
!                 '/home/meteo/yxl232/wave/tmp/LUT_compare_MieSphere/'
  integer :: ncid
  integer :: n_freq, n_temp, n_Re, n_Legendre
  integer :: status

  type(CloudCoeff_Type) :: cloudcoeff 

  character(len=256) :: filename_part
  character(len=256) :: filename_surfix

  call GETARG(1, filename_part)
  call GETARG(2, filename_surfix)
!  call nc2cloudcoeff( LUT_DIR // 'WSM6_GraupelLUT_09z_Liu08.nc', LUT_PHASE_LIQUID, cloudcoeff)
!
!  print *,  cloudcoeff%n_MW_Frequencies, &
!            cloudcoeff%n_Temperatures,   &
!            cloudcoeff%n_MW_Radii, &
!            cloudcoeff%n_Legendre_Terms, &
!            cloudcoeff%n_Densities, &
!            cloudcoeff%n_Phase_Elements


!  call nc2binary( LUT_DIR // 'WSM6_CloudIceLUT_MieSphere.nc', LUT_PHASE_SOLID, &
!                  LUT_DIR // 'WSM6_CloudIceLUT_MieSphere.bin' )

!  call nc2binary( LUT_DIR // 'WSM6_CloudWaterLUT_MieSphere.nc', LUT_PHASE_LIQUID, &
!                  LUT_DIR // 'WSM6_CloudWaterLUT_MieSphere.bin' )

!  call nc2binary( LUT_DIR // 'WSM6_GraupelLUT_09z_Liu08.nc', LUT_PHASE_SOLID, &
!                  LUT_DIR // 'WSM6_GraupelLUT_09z_Liu08.bin' )

!  call nc2binary( LUT_DIR // 'WSM6_SnowLUT_09z_Liu08.nc', LUT_PHASE_SOLID, &
!                  LUT_DIR // 'WSM6_SnowLUT_09z_Liu08.bin' )

!  call nc2binary( LUT_DIR // 'WSM6_RainLUT_09z_Liu08.nc', LUT_PHASE_LIQUID, &
!                  LUT_DIR // 'WSM6_RainLUT_09z_Liu08.bin' )


!  call nc2binary( LUT_DIR // 'WSM6_GraupelLUT_MieSphere_dRe10.nc', LUT_PHASE_SOLID, &
!                  LUT_DIR // 'WSM6_GraupelLUT_MieSphere_dRe10.bin' )

!  call nc2binary( LUT_DIR // 'WSM6_SnowLUT_MieSphere_dRe10.nc', LUT_PHASE_SOLID, &
!                  LUT_DIR // 'WSM6_SnowLUT_MieSphere_dRe10.bin' )

!  call nc2binary( LUT_DIR // 'WSM6_RainLUT_MieSphere_dRe10.nc', LUT_PHASE_LIQUID, &
!                  LUT_DIR // 'WSM6_RainLUT_MieSphere_dRe10.bin' )

  call nc2binary( trim(filename_part) // '_SnowLUT_' // trim(filename_surfix) // '.nc', LUT_PHASE_SOLID, &
                  trim(filename_part) // '_SnowLUT_' // trim(filename_surfix) // '.bin' )
  call nc2binary( trim(filename_part) // '_GraupelLUT_' // trim(filename_surfix) // '.nc', LUT_PHASE_SOLID, &
                  trim(filename_part) // '_GraupelLUT_' // trim(filename_surfix) // '.bin' )
  call nc2binary( trim(filename_part) // '_RainLUT_' // trim(filename_surfix) // '.nc', LUT_PHASE_LIQUID, &
                  trim(filename_part) // '_RainLUT_' // trim(filename_surfix) // '.bin' )


end program test_nc



