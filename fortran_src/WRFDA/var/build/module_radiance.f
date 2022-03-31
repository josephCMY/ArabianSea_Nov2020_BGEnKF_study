












module module_radiance

   
   
   

   use da_control, only : pi, use_landem, t_landem, t_kelvin
   use da_reporting, only : da_error,message


  
   
   
   

  
   USE CRTM_Module, only : graupel_cloud, rain_cloud, snow_cloud,crtm_adjoint, &
      crtm_atmosphere_create, crtm_surface_create, &
      crtm_atmosphere_destroy, crtm_surface_destroy, &
      crtm_forward,crtm_init,crtm_k_matrix, &
      crtm_tangent_linear, h2o_id,hail_cloud,ice_cloud, &
      o3_id, water_cloud, crtm_rtsolution_type, crtm_channelinfo_type, &
      crtm_atmosphere_type, crtm_surface_type, crtm_geometry_type, &
      crtm_surface_zero, crtm_atmosphere_zero, crtm_destroy, &
      climatology_model_name, &
      crtm_options_type, crtm_options_create, crtm_options_destroy, &
      crtm_rtsolution_create, crtm_rtsolution_destroy, crtm_rtsolution_associated, &
      crtm_irlandcoeff_classification
   USE CRTM_Atmosphere_Define, only: crtm_atmosphere_associated, &
      MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS
   USE CRTM_Surface_Define, only: crtm_surface_associated
   USE CRTM_Options_Define, only: crtm_options_associated

   USE CRTM_SensorInfo
   USE CRTM_Planck_Functions, only : CRTM_Planck_Temperature, &
      CRTM_Planck_Radiance, CRTM_Planck_Temperature_TL, &
      CRTM_Planck_Temperature_AD

   use gsi_kinds      ,  only : r_kind,r_double,i_kind,r_single
   use gsi_constants  ,  only : deg2rad, rad2deg,       &
                            init_constants_derived, &
                            one, three, zero, half, &
                            one_tenth, two, four

   
   implicit none
   
   real, parameter             :: q2ppmv = 1.60771704e+6   

  
  
  
  Character (len=8), Parameter :: rttov_platform_name(1:35) =          &
     & (/ 'noaa    ', 'dmsp    ', 'meteosat', 'goes    ', 'gms     ',  &
        & 'fy2     ', 'trmm    ', 'ers     ', 'eos     ', 'metop   ',  &
        & 'envisat ', 'msg     ', 'fy1     ', 'adeos   ', 'mtsat   ',  &
        & 'coriolis', 'jpss    ', 'gifts   ', 'tiros   ', 'meghatr ',  &
        & 'kalpana ', 'reserved', 'fy3     ', 'coms    ', 'meteor-m',  &
        & 'gosat   ', 'calipso ', 'reserved', 'gcom-w  ', 'nimbus  ',  &
        & 'himawari', 'mtg     ', 'saral   ', 'metop-ng', 'landsat '/)

  
  
  Character (len=8), Dimension(0:65) :: rttov_inst_name  =             &
     & (/ 'hirs    ', 'msu     ', 'ssu     ', 'amsua   ', 'amsub   ',  &
        & 'avhrr   ', 'ssmi    ', 'vtpr1   ', 'spare   ', 'tmi     ',  &
        & 'ssmis   ', 'airs    ', 'hsb     ', 'modis   ', 'atsr    ',  &
        & 'mhs     ', 'iasi    ', 'amsre   ', 'imager  ', 'atms    ',  &
        & 'mviri   ', 'seviri  ', 'imager  ', 'sounder ', 'imager  ',  &
        & 'vissr   ', 'mvisr   ', 'cris    ', 'spare   ', 'viirs   ',  &
        & 'windsat ', 'gifts   ', 'ssmt1   ', 'ssmt2   ', 'saphir  ',  &
        & 'madras  ', 'spare   ', 'imager  ', 'reserved', 'reserved',  &
        & 'mwts    ', 'mwhs    ', 'iras    ', 'mwri    ', 'abi     ',  &
        & 'mi      ', 'msumr   ', 'reserved', 'iir     ', 'mwr     ',  &
        & 'reserved', 'reserved', 'reserved', 'reserved', 'scams   ',  &
        & 'smmr    ', 'ahi     ', 'irs     ', 'altika  ', 'iasing  ',  &
        & 'tm      ', 'fci     ', 'amsr1   ', 'amsr2   ', 'vissr   ',  &
        & 'slstr   '/)

  
  
  
  
  
  
  Character (len=8), Parameter :: crtm_platform_name(1:35) =           &
     & (/ 'n       ', 'f       ', 'm       ', 'g       ', 'gms     ',  &
        & 'xxxxxxxx', 'trmm    ', 'ers     ', 'eos     ', 'metop   ',  &
        & 'envisat ', 'msg     ', 'xxxxxxxx', 'xxxxxxxx', 'mt      ',  &
        & 'coriolis', 'npp     ', 'gifts   ', 'tiros   ', 'meghat  ',  &
        & 'kalpana ', 'tiros   ', 'fy3     ', 'coms    ', 'xxxxxxxx',  &
        & 'xxxxxxxx', 'xxxxxxxx', 'reserved', 'gcom-w  ', 'xxxxxxxx',  &
        & 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx'/)

  
  
  
  
  
  
  Character (len=8), Dimension(0:65) :: crtm_sensor_name  =            &
     & (/ 'hirs    ', 'msu     ', 'ssu     ', 'amsua   ', 'amsub   ',  &
        & 'avhrr   ', 'ssmi    ', 'xxxxxxxx', 'spare   ', 'tmi     ',  &
        & 'ssmis   ', 'airs    ', 'hsb     ', 'modis   ', 'atsr    ',  &
        & 'mhs     ', 'iasi    ', 'amsre   ', 'imgr    ', 'atms    ',  &
        & 'mviri   ', 'seviri  ', 'imgr    ', 'sndr    ', 'imgr    ',  &
        & 'vissr   ', 'xxxxxxxx', 'cris    ', 'spare   ', 'viirs   ',  &
        & 'windsat ', 'xxxxxxxx', 'ssmt1   ', 'ssmt2   ', 'saphir  ',  &
        & 'madras  ', 'spare   ', 'imgr    ', 'reserved', 'reserved',  &
        & 'mwts    ', 'mwhs    ', 'iras    ', 'mwri    ', 'abi     ',  &
        & 'xxxxxxxx', 'xxxxxxxx', 'reserved', 'xxxxxxxx', 'xxxxxxxx',  &
        & 'reserved', 'reserved', 'reserved', 'reserved', 'xxxxxxxx',  &
        & 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx',  &
        & 'xxxxxxxx', 'xxxxxxxx', 'xxxxxxxx', 'amsr2   ', 'vissr   ',  &
        & 'xxxxxxxx'/)


   type satinfo_type
      integer, pointer   :: ichan(:)      
      integer, pointer   :: iuse (:)      
      real   , pointer   :: error(:)      
      real   , pointer   :: polar(:)      
      real   , pointer   :: error_factor(:) 
     
      real   , pointer   :: scanbias(:,:) 
      real   , pointer   :: scanbias_b(:,:,:) 
      real   , pointer   :: bcoef(:,:)   
      real   , pointer   :: bcoef0(:)    
      real   , pointer   :: error_std(:) 
   end type satinfo_type

   type (satinfo_type), pointer :: satinfo(:)

   CHARACTER( 80 ), allocatable, save :: Sensor_Descriptor(:)

end module module_radiance

