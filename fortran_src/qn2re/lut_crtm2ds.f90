program lut_crtm2ds
  use CloudCoeff_Binary_IO, only: CloudCoeff_Binary_ReadFile, &
                                  CloudCoeff_Binary_WriteFile
  use CloudCoeff_Define,    only: CloudCoeff_type, &
                                  CloudCoeff_Associated, &
                                  CloudCoeff_Destroy, &
                                  CloudCoeff_Create, &
                                  CloudCoeff_ValidRelease, &
                                  CloudCoeff_Info

  implicit none



  character(len=*), parameter :: filename_crtmlut = &
        '/home/meteo/yxl232/wave/tmp/LUT_crtm2ds/CloudCoeff.bin'
  character(len=*), parameter :: dir_DSlut = &
        '/home/meteo/yxl232/wave/tmp/LUT_crtm2ds/'

  type(CloudCoeff_type) :: cloudcoeff
  type(CloudCoeff_type) :: cloudcoeff_rain
  type(CloudCoeff_type) :: cloudcoeff_snow
  type(CloudCoeff_type) :: cloudcoeff_graup
  type(CloudCoeff_type) :: cloudcoeff_hail
  integer :: status

  

  ! read CloudCoeff from binary CRTM LUT file
  status = CloudCoeff_Binary_ReadFile( filename_crtmlut, cloudcoeff )

  ! create structures for new LUTs of rain, snow, and graupel
  call CloudCoeff_Create( cloudcoeff_rain, &
                          cloudcoeff%n_MW_Frequencies, &
                          cloudcoeff%n_MW_Radii, &
                          1, &        ! n_IR_freq
                          1, &        ! n_IR_Radii
                          cloudcoeff%n_Temperatures, &
                          1, &        ! n_density
                          cloudcoeff%n_Legendre_Terms, &
                          cloudcoeff%n_Phase_Elements)          ! n_phase_elements

  call CloudCoeff_Create( cloudcoeff_snow, &
                          cloudcoeff%n_MW_Frequencies, &
                          cloudcoeff%n_MW_Radii, &
                          1, & ! n_IR_freq
                          1, & ! n_IR_Radii
                          1, & ! n_temperature. Not needed for solid hydro
                          1, & ! n_density
                          cloudcoeff%n_Legendre_Terms, &
                          cloudcoeff%n_Phase_Elements)          ! n_phase_elements

  call CloudCoeff_Create( cloudcoeff_graup, &
                          cloudcoeff%n_MW_Frequencies, &
                          cloudcoeff%n_MW_Radii, &
                          1, & ! n_IR_freq
                          1, & ! n_IR_Radii
                          1, & ! n_temperature. Not needed for solid hydro
                          1, & ! n_density
                          cloudcoeff%n_Legendre_Terms, &
                          cloudcoeff%n_Phase_Elements)          ! n_phase_elements
  
  call CloudCoeff_Create( cloudcoeff_hail, &
                          cloudcoeff%n_MW_Frequencies, &
                          cloudcoeff%n_MW_Radii, &
                          1, & ! n_IR_freq
                          1, & ! n_IR_Radii
                          1, & ! n_temperature. Not needed for solid hydro
                          1, & ! n_density
                          cloudcoeff%n_Legendre_Terms, &
                          cloudcoeff%n_Phase_Elements)          ! n_phase_elements
  

  ! dimension vectors
  CloudCoeff_rain%Frequency_MW = CloudCoeff%Frequency_MW
  CloudCoeff_rain%Reff_MW      = CloudCoeff%Reff_MW
  CloudCoeff_rain%Temperature  = CloudCoeff%Temperature
  CloudCoeff_rain%Density      = 1   ! dummy
  CloudCoeff_rain%Frequency_IR = 200 ! dummy
  CloudCoeff_rain%Reff_IR      = 200 ! dummy

  CloudCoeff_snow%Frequency_MW = CloudCoeff%Frequency_MW
  CloudCoeff_snow%Reff_MW      = CloudCoeff%Reff_MW
  CloudCoeff_snow%Density      = 1   ! dummy
  CloudCoeff_snow%Frequency_IR = 200 ! dummy
  CloudCoeff_snow%Reff_IR      = 200 ! dummy

  CloudCoeff_graup%Frequency_MW = CloudCoeff%Frequency_MW
  CloudCoeff_graup%Reff_MW      = CloudCoeff%Reff_MW
  CloudCoeff_graup%Density      = 1   ! dummy
  CloudCoeff_graup%Frequency_IR = 200 ! dummy
  CloudCoeff_graup%Reff_IR      = 200 ! dummy

  CloudCoeff_hail%Frequency_MW = CloudCoeff%Frequency_MW
  CloudCoeff_hail%Reff_MW      = CloudCoeff%Reff_MW
  CloudCoeff_hail%Density      = 1   ! dummy
  CloudCoeff_hail%Frequency_IR = 200 ! dummy
  CloudCoeff_hail%Reff_IR      = 200 ! dummy

  ! scattering properties
  CloudCoeff_rain%ke_L_MW     = CloudCoeff%ke_L_MW
  CloudCoeff_rain%w_L_MW      = CloudCoeff%w_L_MW
  CloudCoeff_rain%g_L_MW      = CloudCoeff%g_L_MW
  CloudCoeff_rain%pcoeff_L_MW = CloudCoeff%pcoeff_L_MW
 
  ! snow density index in CRTM LUT is 1 
  CloudCoeff_snow%ke_S_MW(:,:,1)         = CloudCoeff%ke_S_MW(:,:,1)
  CloudCoeff_snow%w_S_MW (:,:,1)         = CloudCoeff%w_S_MW (:,:,1)
  CloudCoeff_snow%g_S_MW (:,:,1)         = CloudCoeff%g_S_MW (:,:,1)
  CloudCoeff_snow%pcoeff_S_MW(:,:,1,:,:) = CloudCoeff%pcoeff_S_MW(:,:,1,:,:)
  
  ! graupel density index in CRTM LUT is 2
  CloudCoeff_graup%ke_S_MW(:,:,1)         = CloudCoeff%ke_S_MW(:,:,2)
  CloudCoeff_graup%w_S_MW (:,:,1)         = CloudCoeff%w_S_MW (:,:,2)
  CloudCoeff_graup%g_S_MW (:,:,1)         = CloudCoeff%g_S_MW (:,:,2)
  CloudCoeff_graup%pcoeff_S_MW(:,:,1,:,:) = CloudCoeff%pcoeff_S_MW(:,:,2,:,:)
  
  ! hail density index in CRTM LUT is 3
  CloudCoeff_hail%ke_S_MW(:,:,1)         = CloudCoeff%ke_S_MW(:,:,3)
  CloudCoeff_hail%w_S_MW (:,:,1)         = CloudCoeff%w_S_MW (:,:,3)
  CloudCoeff_hail%g_S_MW (:,:,1)         = CloudCoeff%g_S_MW (:,:,3)
  CloudCoeff_hail%pcoeff_S_MW(:,:,1,:,:) = CloudCoeff%pcoeff_S_MW(:,:,3,:,:)
  
  ! output those new LUTs
  status = CloudCoeff_Binary_WriteFile( dir_DSlut // 'cloudcoeff_rain.bin', cloudcoeff_rain)
  status = CloudCoeff_Binary_WriteFile( dir_DSlut // 'cloudcoeff_snow.bin', cloudcoeff_snow)
  status = CloudCoeff_Binary_WriteFile( dir_DSlut // 'cloudcoeff_graup.bin', cloudcoeff_graup)
  status = CloudCoeff_Binary_WriteFile( dir_DSlut // 'cloudcoeff_hail.bin', cloudcoeff_hail)

  call CloudCoeff_Destroy( cloudcoeff )
  call CloudCoeff_Destroy( cloudcoeff_rain )
  call CloudCoeff_Destroy( cloudcoeff_snow )
  call CloudCoeff_Destroy( cloudcoeff_graup )
  call CloudCoeff_Destroy( cloudcoeff_hail )


end program lut_crtm2ds



