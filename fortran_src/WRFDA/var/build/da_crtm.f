












module da_crtm

   
   
   


   use module_domain, only : x_type, domain
   use da_define_structures, only : y_type, iv_type

   use module_radiance, only : CRTM_RTSolution_type,CRTM_ChannelInfo_type, &
      CRTM_Atmosphere_type, CRTM_Surface_type,CRTM_Geometry_type, &
      CRTM_Adjoint,CRTM_Forward,CRTM_Tangent_Linear, &
      CRTM_K_Matrix, CRTM_Planck_Temperature, CRTM_Planck_Temperature_TL, &
      CRTM_Planck_Temperature_AD, CRTM_Planck_Radiance, &
      CRTM_Atmosphere_Create,H2O_ID,GRAUPEL_CLOUD,ICE_CLOUD,HAIL_CLOUD, &
      rain_cloud,snow_cloud,O3_ID, &
      WATER_CLOUD, Sensor_Descriptor, MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS, &
      crtm_atmosphere_destroy, crtm_sensor_name, &
      crtm_surface_create,crtm_surface_destroy, &
      crtm_surface_zero, CRTM_Atmosphere_zero, satinfo, &
      crtm_platform_name, crtm_init, &
      rttov_inst_name,rttov_platform_name, climatology_model_name, &
      crtm_options_type, crtm_options_create, crtm_options_destroy, &
      crtm_atmosphere_associated, crtm_surface_associated, &
      crtm_options_associated, &
      crtm_rtsolution_create, crtm_rtsolution_destroy, crtm_rtsolution_associated, &
      crtm_irlandcoeff_classification

   use da_control, only : trace_use, crtm_cloud, gravity,stdout, biascorr, &
      biasprep, qc_rad,missing_r,rtminit_sensor,rtminit_nsensor, filename_len, &
      use_error_factor_rad,read_biascoef, analysis_date,time_window_max, &
      time_window_min,num_fgat_time,rtminit_platform, print_detail_rad, &
      rtminit_satid, global,kms,kme,ims,ime,jms,jme,kts,kte,use_clddet_mmr, &
      use_crtm_kmatrix, use_varbc, freeze_varbc, use_pseudo_rad, &
      use_antcorr, time_slots, use_satcv, use_simulated_rad, simulated_rad_io, &
      simulated_rad_ngrid, interp_option, use_mspps_emis, use_mspps_ts, calc_weightfunc, &
      use_clddet_ecmwf,its,ite,jts,jte, &
      crtm_coef_path, crtm_irwater_coef, crtm_mwwater_coef, crtm_irland_coef, crtm_visland_coef
   use da_interpolation, only : da_interp_lin_2d_partial,da_interp_lin_2d_adj_partial, &
      da_interp_2d_partial
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_reals
   use da_radiance1, only : da_biasprep,da_detsurtyp,da_biascorr, &
       da_biasprep,da_cld_eff_radius, da_mspps_emis, da_mspps_ts

   use da_reporting, only : da_error, message, da_warning, da_message
   use da_tools_serial, only : da_free_unit, da_get_unit
   use da_tools, only: da_get_time_slots, da_eof_decomposition
   use da_tracing, only : da_trace_entry, da_trace_exit

   TYPE (CRTM_ChannelInfo_type), allocatable, save :: ChannelInfo(:)

   
   
   
   integer, parameter :: n_soil_type = 16  
   integer, parameter :: USGS_n_type = 24
   integer, parameter :: IGBP_n_type = 20
   
   
   
   integer, parameter :: wrf_to_crtm_soil(n_soil_type) = &
      (/ 1, 1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9 /)
   
   
   integer, parameter :: usgs_to_crtm_mw(USGS_n_type) = &
      (/  7, 12, 12, 12, 12, 12,  7,  9,  8,  6, &
          2,  5,  1,  4,  3,  0,  8,  8, 11, 10, &
         10, 10, 11, 13 /)
   integer, parameter :: igbp_to_crtm_mw(IGBP_n_type) = &
      (/  4,  1,  5,  2,  3,  8,  9,  6,  6,  7, &
          8, 12,  7, 12, 13, 11,  0, 10, 10, 11 /)

contains

subroutine da_transform_xtoy_crtm (cv_size, cv, grid, iv, y )

   !---------------------------------------------------------------------------
   !  PURPOSE: transform from analysis increment to 
   !                          pertubation radiance.
   !
   !  METHOD:  delta_y = H delta_x
   !           1. input reference state of CRTM_TL
   !           2. interpolate analysis increment to obs location
   !           3. Call CRTM_TL
   !
   !  HISTORY: 11/16/2006 - Creation                        Zhiquan Liu
   !           11/25/2008 - Zero Jacobian for top level     Tom Auligne
   !
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)            :: cv_size         ! Size of cv array.
   real, intent(in)               :: cv(1:cv_size)   ! control variables.
   type (domain),  intent(in)     :: grid
   type (y_type),  intent(inout)  :: y        ! H' delta_x
   type (iv_type), intent(in)     :: iv       ! O-B structure.

   integer, parameter             :: AIRS_Max_Channels = 281

   integer                        :: k, l  ! Index dimension.

   integer           :: inst, num_rad, nchanl, n, icld
   integer           :: ipred, npred, gammapred, id
   real, allocatable :: temperature(:,:)
   real, allocatable :: absorber(:,:), xb_q(:,:)
   real, allocatable :: psfc(:)
!! for crtm_cloud
   real, allocatable :: qcw(:,:),qci(:,:),qrn(:,:),qsn(:,:),qgr(:,:)

   ! 1 local variables and types
   integer :: Allocate_Status
   integer :: n_layers, n_absorbers, n_clouds, n_aerosols
   type( CRTM_RTSolution_type ), ALLOCATABLE :: RTSolution(:,:),RTSolution_TL(:,:)
   type( CRTM_Atmosphere_type ), allocatable :: Atmosphere(:), Atmosphere_TL(:)
   type( CRTM_Surface_type ),    allocatable :: Surface(:), Surface_TL(:)
   type( CRTM_Geometry_type ),   allocatable :: GeometryInfo(:)
   type( CRTM_Options_type ), allocatable    :: Options(:)
   
   integer                        :: ts_index
   integer                        :: nclouds, ncv
   real, allocatable              :: cc_tl(:)
   real*8                         :: rad_clr, rad_cld   ! RT clear/cloudy radiances
   real, allocatable              :: rad_ovc(:)         ! RT overcast radiances
   real*8                         :: rad_tl, tb_tl

   ! Initializations for AIRS (MMR) Cloud Detection
   integer                        :: Band_Size(5), Bands(AIRS_Max_Channels,5) 
  
      Band_Size(1:5) = (/86, 0, 0, 16, 0 /)
      Bands(:,:)     = 0  
      Bands(1:Band_Size(1),1) = &
&    (/                                                 &              !&      1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
&                                                       &              !&     22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
&                                                       &              !&     52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
&                                                       &              !&     72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
&                     92,  93,  98,  99, 101, 104, 105, &              !&     83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
&     108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
&     138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
&     165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
&     180, 182, 185, 186, 190, 192,      198, 201, 204, &              !&     180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
&     207, 210,      215, 216,      221,      226, 227, &              !&     207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
&     232,                     252, 253, 256, 257, 261, &              !&     232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
&     262, 267, 272, 295, 299,      305,           310, &              !&     262, 267, 272, 295, 299, 300, 305, 308, 309, 310, &
&          321, 325, 333, 338, 355, 362, 375, 453, 475, &              !&     318, 321, 325, 333, 338, 355, 362, 375, 453, 475, &
&     484, 497, 528, 587, 672, 787, 791, 843, 870, 914, &
&     950 /)

!      Bands(1:Band_Size(2),2) = &
!&    (/ 1003, 1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082,  &
!&       1083, 1088, 1090, 1092, 1095, 1104, 1111, 1115, 1116, 1119,  &
!&       1120, 1123, 1130, 1138, 1142, 1178, 1199, 1206, 1221, 1237,  &
!&       1252, 1260, 1263, 1266, 1278, 1285 /)

!      Bands(1:Band_Size(3),3) = &
!&    (/       1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  !&    1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  
!&       1466,       1477,             1500, 1519,       1538, 1545, &  !&    1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &  
!&       1565, 1574, 1583, 1593,       1627, 1636,       1652, 1669, &  !&    1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, & 
!&                   1694, 1708,       1723, 1740, 1748,       1756, &  !&    1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
!&             1766, 1771, 1777,       1783, 1794, 1800,       1806, &  !&    1763, 1766, 1771, 1777, 1780, 1783, 1794, 1800, 1803, 1806, &
!&             1826, 1843  /)                                           !&    1812, 1826, 1843  /)

      Bands(1:Band_Size(4),4) = &
&    (/ 1852, 1865, 1866,       1868, 1869, 1872, 1873,       1876, &  !&    1852, 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, 
&             1881, 1882, 1883,                   1911, 1917, 1918, &  !&    1877, 1881, 1882, 1883, 1884, 1897, 1901, 1911, 1917, 1918, &
&                   1924, 1928        /)                               !&    1921, 1923, 1924, 1928, 1937  /)   

!      Bands(1:Band_Size(5),5) = &
!&    (/ 1938, 1939, 1941, 1946, 1947, 1948, 1958, 1971, 1973, 1988, &
!&       1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, 2104, &
!&       2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, &
!&       2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, 2134, &
!&       2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, 2234, &
!&       2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, 2355, &
!&       2363, 2370, 2371, 2377  /)  

!---------------------------------------------------------

   if ( iv%num_inst < 1 ) return

   if (trace_use) call da_trace_entry("da_transform_xtoy_crtm")

   sensor_loop : do inst = 1, iv%num_inst                 ! loop for sensor

      num_rad = iv%instid(inst)%info%n2 - iv%instid(inst)%info%n1 + 1
      if ( num_rad < 1 ) cycle

      allocate (Atmosphere   (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                Atmosphere_TL(iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",119, &
          (/"Error in allocatting Atmosphere"/))
      end if
      allocate (Surface      (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                Surface_TL   (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",126, &
          (/"Error in allocatting Surface"/))
      end if
      allocate (GeometryInfo (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",132, &
          (/"Error in allocatting GeometryInfo"/))
      end if
      allocate (Options      (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",138, &
          (/"Error in allocatting Options"/))
      end if

!----------------------------------------------------------------------------
! 1 allocation
!
! Atmosphere structure

      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         n_layers    = (kte-kts)+1   ! number of vertical levels
         n_absorbers = 2
         n_aerosols  = 0
         n_clouds    = 0
         if ( crtm_cloud ) n_clouds = 6

         call CRTM_Atmosphere_Create( Atmosphere(n), &
                                      n_layers,      &
                                      n_absorbers,   &
                                      n_clouds,      &
                                      n_aerosols )
         IF ( .NOT. CRTM_Atmosphere_Associated( Atmosphere(n) ) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",161, &
                         (/"Error in allocatting CRTM Atmosphere Structure"/))
         end if

         Atmosphere(n)%Absorber_ID(1)=H2O_ID
         Atmosphere(n)%Absorber_ID(2)=O3_ID
         Atmosphere(n)%Absorber_Units(1) = MASS_MIXING_RATIO_UNITS
         Atmosphere(n)%Absorber_Units(2) = VOLUME_MIXING_RATIO_UNITS
         Atmosphere(n)%Climatology=iv%instid(inst)%crtm_climat(n)

         if (crtm_cloud) then
            Atmosphere(n)%Cloud(1)%Type=WATER_CLOUD
            Atmosphere(n)%Cloud(2)%Type=ICE_CLOUD
            Atmosphere(n)%Cloud(3)%Type=RAIN_CLOUD
            Atmosphere(n)%Cloud(4)%Type=SNOW_CLOUD
            Atmosphere(n)%Cloud(5)%Type=GRAUPEL_CLOUD
            Atmosphere(n)%Cloud(6)%Type=HAIL_CLOUD
         end if

      end do

!-------------------------------------------------------------------------------

      nchanl    = ChannelInfo(inst)%n_channels
                                        
  ! Allocate forward model solution RTSolution array to number of channels
      allocate( RTSolution( ChannelInfo(inst)%n_Channels, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2 ), &
                RTSolution_TL( ChannelInfo(inst)%n_Channels, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2 ), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",191, &
                      (/"Error in allocatting RTSolution"/))
      END IF

      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         call CRTM_Surface_Create( Surface(n),  & ! Output
                                   nchanl )       ! Input
         IF ( .NOT. CRTM_Surface_Associated( Surface(n) ) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",200, &
                         (/"Error in allocatting CRTM Surface Structure"/))
         end if

  ! 1 Options structure
         call CRTM_Options_Create( Options(n),  &  ! Output
                                   nchanl )        ! Input
         IF ( .NOT. CRTM_Options_Associated( Options(n) ) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",208, &
                         (/"Error in allocatting CRTM Options Structure"/))
         endif
         if ( use_antcorr(inst) ) Options(n)%Use_Antenna_Correction = .true.

      end do

      allocate (temperature(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
      allocate (absorber(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
      allocate (xb_q(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
      allocate (psfc(iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))

      temperature(:,:) = 0.0
      absorber(:,:)    = 0.0
      xb_q(:,:)        = 0.0
      psfc(:)          = 0.0

      if (crtm_cloud) then

         allocate (qcw(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qci(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qrn(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qsn(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qgr(Atmosphere(iv%instid(inst)%info%n1)%n_layers, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         qcw(:,:) = 0.0
         qci(:,:) = 0.0
         qrn(:,:) = 0.0
         qsn(:,:) = 0.0
         qgr(:,:) = 0.0

      end if

      do k=kts,kte ! from bottom to top
         call da_interp_2d_partial (grid%xa%t(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            temperature(kte-k+1,:))
         call da_interp_2d_partial (grid%xa%q(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            absorber(kte-k+1,:))

         if (crtm_cloud) then

         call da_interp_2d_partial (grid%xa%qcw(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            qcw(kte-k+1,:))
         call da_interp_2d_partial (grid%xa%qci(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            qci(kte-k+1,:))
         call da_interp_2d_partial (grid%xa%qrn(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            qrn(kte-k+1,:))
         call da_interp_2d_partial (grid%xa%qsn(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            qsn(kte-k+1,:))
         call da_interp_2d_partial (grid%xa%qgr(:,:,k), iv%instid(inst)%info, k, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, &
            qgr(kte-k+1,:))

         end if

      end do

      call da_interp_2d_partial (grid%xa%psfc, iv%instid(inst)%info, 1, iv%instid(inst)%info%n1, iv%instid(inst)%info%n2, psfc(:))

      ! Gamma correction from VarBC
      !----------------------------

      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
            
      ! [1.1] Get horizontal interpolation weights:

      ! [1.3] Extract base state Atmosphere variables 
         Atmosphere(n)%level_pressure(0) = iv%instid(inst)%pf(0,n)
         do k=1,Atmosphere(n)%n_layers
            Atmosphere(n)%pressure(k)       = iv%instid(inst)%pm(k,n)
            Atmosphere(n)%level_pressure(k) = iv%instid(inst)%pf(k,n)
            Atmosphere(n)%temperature(k)    = iv%instid(inst)%tm(k,n)
            Atmosphere(n)%absorber(k,1)     = iv%instid(inst)%qm(k,n)
            xb_q(k,n) = 0.001*iv%instid(inst)%qm(k,n)/(1.0+iv%instid(inst)%qm(k,n)) ! specific humidity
         end do
         if (crtm_cloud) then
            do k=1,Atmosphere(n)%n_layers
              Atmosphere(n)%cloud(1)%water_content(k)=iv%instid(inst)%qcw(k,n)
              Atmosphere(n)%cloud(2)%water_content(k)=iv%instid(inst)%qci(k,n)
              Atmosphere(n)%cloud(3)%water_content(k)=iv%instid(inst)%qrn(k,n)
              Atmosphere(n)%cloud(4)%water_content(k)=iv%instid(inst)%qsn(k,n)
              Atmosphere(n)%cloud(5)%water_content(k)=iv%instid(inst)%qgr(k,n)
              Atmosphere(n)%cloud(6)%water_content(k)=iv%instid(inst)%qhl(k,n)
              Atmosphere(n)%cloud(1)%effective_radius(k)=iv%instid(inst)%rcw(k,n)
              Atmosphere(n)%cloud(2)%effective_radius(k)=iv%instid(inst)%rci(k,n)
              Atmosphere(n)%cloud(3)%effective_radius(k)=iv%instid(inst)%rrn(k,n)
              Atmosphere(n)%cloud(4)%effective_radius(k)=iv%instid(inst)%rsn(k,n)
              Atmosphere(n)%cloud(5)%effective_radius(k)=iv%instid(inst)%rgr(k,n)
              Atmosphere(n)%cloud(6)%effective_radius(k)=iv%instid(inst)%rhl(k,n)
           end do
         end if

      ! [1.4] User-supplied emissivity
           !Options%emissivity_switch = 1
           !Options%emissivity(1:Options%n_channels) = &
           !    iv%instid(inst)%emiss(1:Options%n_channels,n)

      ! [1.4] 1 Surface parameter data
         Surface(n)%Land_Coverage=iv%instid(inst)%land_coverage(n) 
         Surface(n)%Water_Coverage=iv%instid(inst)%water_coverage(n) 
         Surface(n)%Snow_Coverage=iv%instid(inst)%snow_coverage(n)
         Surface(n)%Ice_Coverage=iv%instid(inst)%ice_coverage(n)

         if (Surface(n)%Land_Coverage > 0.0) then
            Surface(n)%Land_Type=11 !GRASS_SOIL
            Surface(n)%Land_Temperature=iv%instid(inst)%ts(n)      ! K
            Surface(n)%Soil_Moisture_Content=iv%instid(inst)%smois(n)  !0.05    ! volumetric water content (g/cm**3)
            !Surface(n)%Canopy_Water_Content=0.05      ! gravimetric water content
            Surface(n)%Vegetation_Fraction=iv%instid(inst)%vegfra(n)
            Surface(n)%Soil_Temperature=iv%instid(inst)%tslb(n)
         end if
         if (Surface(n)%Water_Coverage > 0.0) then
            !Surface(n)%Water_Type=SEA_WATER          ! (Currently NOT used)
            Surface(n)%Water_Temperature=iv%instid(inst)%ts(n)     ! K
            Surface(n)%Wind_Speed=sqrt((iv%instid(inst)%u10(n))**2+ &
                                      (iv%instid(inst)%v10(n))**2)  ! m/sec
            !surface(n)%Wind_Direction=0.0            ! NOT used
            Surface(n)%Salinity=33.                   ! ppmv
         end if
         if (Surface(n)%Snow_Coverage > 0.0) then
            Surface(n)%Snow_Type=2 !NEW_SNOW
            Surface(n)%Snow_Temperature=iv%instid(inst)%ts(n)      ! K
            Surface(n)%Snow_Depth=iv%instid(inst)%snowh(n)         ! mm
            !Surface(n)%Snow_Density=0.2               ! g/cm**3
            !Surface(n)%Snow_Grain_Size=2.0            ! mm
         end if
         if (Surface(n)%Ice_Coverage > 0.0) then
            !Surface(n)%Ice_Type=FRESH_ICE             ! NO Table offered, single egrid%xample is FRESH_ICE
            Surface(n)%Ice_Temperature=iv%instid(inst)%ts(n)       ! K
            Surface(n)%Ice_Thickness=10.0              ! mm
            !Surface(n)%Ice_Density=0.9                ! g/cm**3
            !Surface(n)%Ice_Roughness=0.0               ! NO Table offered, single egrid%xample is ZERO
         end if
         Surface(n)%SensorData%n_channels        = nchanl
         Surface(n)%SensorData%Sensor_Channel    = ChannelInfo(inst)%Sensor_Channel
         Surface(n)%SensorData%Sensor_Id         = ChannelInfo(inst)%Sensor_Id
         Surface(n)%SensorData%WMO_Satellite_Id  = ChannelInfo(inst)%WMO_Satellite_Id
         Surface(n)%SensorData%WMO_Sensor_Id     = ChannelInfo(inst)%WMO_Sensor_Id
         Surface(n)%SensorData%Tb(1:nchanl) = iv%instid(inst)%tb_inv(1:nchanl,n) + &
                                              iv%instid(inst)%tb_xb(1:nchanl,n)

      ! -- Copy the TL atmosphere structure
         Atmosphere_TL(n) = Atmosphere(n)

      ! -- Copy the TL surface structure
         Surface_TL(n) = Surface(n)

    ! -- Zero the TL inputs
    ! Important: adjoint variables must be initialized
        call CRTM_Atmosphere_Zero( Atmosphere_TL(n) )
        call CRTM_Surface_Zero( Surface_TL(n) )

        ! Humidity in g/kg. Zero Jacobian for top level(s) (above 75hPa)
        do k = kts, kte
          if ( iv%instid(inst)%pm(k,n) < 75.0 ) absorber(k,n) = 0.0
        end do
        Atmosphere_TL(n)%Absorber(kts+1:kte,1)  = 1000.0 / (1.0-xb_q(kts+1:kte,n))**2  &
                                                  * absorber(kts+1:kte,n)

        Atmosphere_TL(n)%Temperature(kts+1:kte) = temperature(kts+1:kte,n)       ! Zero Jacobian for top level
        Atmosphere_TL(n)%Level_Pressure(Atmosphere_TL(n)%n_Layers) = 0.01 * psfc(n)

        if (crtm_cloud) then
           Atmosphere_TL(n)%Cloud(1)%Water_Content(kts:kte)             = qcw(kts:kte,n)
           Atmosphere_TL(n)%Cloud(2)%Water_Content(kts:kte)             = qci(kts:kte,n)
           Atmosphere_TL(n)%Cloud(3)%Water_Content(kts:kte)             = qrn(kts:kte,n)
           Atmosphere_TL(n)%Cloud(4)%Water_Content(kts:kte)             = qsn(kts:kte,n)
           Atmosphere_TL(n)%Cloud(5)%Water_Content(kts:kte)             = qgr(kts:kte,n)
           Atmosphere_TL(n)%Cloud(6)%Water_Content(kts:kte)             = 0.
        ! convert cloud content unit from kg/kg to kg/m^2
           do k=kts,kte
              do icld=1,Atmosphere(n)%n_Clouds

                 Atmosphere_TL(n)%Cloud(icld)%Water_Content(k)=  Atmosphere_TL(n)%Cloud(icld)%Water_Content(k)* &
                                (Atmosphere(n)%Level_Pressure(k)- Atmosphere(n)%Level_Pressure(k-1))*100./gravity
 
              enddo
           enddo
        end if

      ! Skin Temperature
      !-----------------
        if (use_satcv(1)) then
           ts_index = iv%instid(inst)%cv_index(n)%ts 	 
  	   if (Surface(n)%Land_Coverage  > 0.0) Surface_TL(n)%Land_Temperature  = cv(ts_index)
           if (Surface(n)%Water_Coverage > 0.0) Surface_TL(n)%Water_Temperature = cv(ts_index)
           if (Surface(n)%Snow_Coverage  > 0.0) Surface_TL(n)%Snow_Temperature  = cv(ts_index)
           if (Surface(n)%Ice_Coverage   > 0.0) Surface_TL(n)%Ice_Temperature   = cv(ts_index)	 
        end if 
	 
      ! Cloud cover(s)
      !---------------
        if (use_satcv(2)) then
	   nclouds = iv%instid(inst)%cv_index(n)%nclouds
	   ncv     = iv%instid(inst)%cv_index(n)%ncv
	   allocate (rad_ovc (nclouds))
	   allocate (cc_tl   (nclouds))

           cc_tl(:) = cv(iv%instid(inst)%cv_index(n)%cc(:))
!	    !---------------------------------------------------------------
!           ! Change of variable (preconditioning) 
!           !---------------------------------------------------------------
!	   do icld = 1, nclouds
!    	      cc_tl(icld) = SUM( cv(iv%instid(inst)%cv_index(n)%cc) * &
!	                            iv%instid(inst)%cv_index(n)%vtox(icld,1:ncv) )
!	   end do
	    	
	   do k = 1, nchanl
              if (ALL(iv%instid(inst)%ichan(k) /=  Bands(:,1))) cycle   ! Only Channels in Band 1 	       
              rad_clr    = iv%instid(inst)%rad_xb(k,n)             
              rad_ovc(:) = iv%instid(inst)%rad_ovc(k,kte-nclouds+1:kte,n)
              rad_cld    = SUM(cc_tl*rad_ovc) + (1.0 - SUM(cc_tl))*rad_clr
	      rad_tl     = SUM(cc_tl*(rad_ovc-rad_clr))
 	      call CRTM_Planck_Temperature_TL(inst,k,rad_clr,rad_tl,tb_tl)
	      y%instid(inst)%tb(k,n) = y%instid(inst)%tb(k,n) + tb_tl	          
   	   end do			     	 

	   deallocate(cc_tl, rad_ovc)
	else
	   y%instid(inst)%tb(:,n) = 0.0  
        end if 

      ! [1.5] 1 GeometryInfo Structure
        GeometryInfo(n)%Sensor_Zenith_Angle=iv%instid(inst)%satzen(n)
        GeometryInfo(n)%Source_Zenith_Angle=iv%instid(inst)%solzen(n)
        GeometryInfo(n)%iFOV=iv%instid(inst)%scanpos(n)
   !    GeometryInfo(n)%Satellite_Height=830.0
   !    GeometryInfo(n)%Sensor_Scan_Angle=
   !    GeometryInfo(n)%Sensor_Zenith_Angle=
   !    GeometryInfo(n)%Sensor_Scan_Angle=
   !    GeometryInfo(n)%Source_Zenith_Angle=

      end do


      ! [1.6] Call CRTM_TL model

      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( n )
      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
         if (.not. use_crtm_kmatrix) then
            call da_crtm_tl (1, nchanl, 1, Atmosphere(n),   &
                               Surface(n),      &
                               Atmosphere_TL(n),&
                               Surface_TL(n),   &
                               GeometryInfo(n), &
                               ChannelInfo(inst:inst),  &
                               RTSolution(:,n),   &
                               RTSolution_TL(:,n),&
                               Options(n))
         else
            RTSolution_TL(:,n)%brightness_temperature = 0.
            do l = 1, ChannelInfo(inst)%n_Channels
               do k = kts , kte
                  RTSolution_TL(l,n)%brightness_temperature = RTSolution_TL(l,n)%brightness_temperature + &
                                       iv%instid(inst)%t_jacobian(l,k,n) * Atmosphere_TL(n)%Temperature(k) +  &
                                       iv%instid(inst)%q_jacobian(l,k,n) * Atmosphere_TL(n)%absorber(k,1)
                  if (crtm_cloud) then
                     RTSolution_TL(l,n)%brightness_temperature = RTSolution_TL(l,n)%brightness_temperature + &
                        iv%instid(inst)%water_jacobian(l,k,n)   * Atmosphere_TL(n)%Cloud(1)%Water_Content(k) + &
                        iv%instid(inst)%ice_jacobian(l,k,n)     * Atmosphere_TL(n)%Cloud(2)%Water_Content(k) + &
                        iv%instid(inst)%rain_jacobian(l,k,n)    * Atmosphere_TL(n)%Cloud(3)%Water_Content(k) + &
                        iv%instid(inst)%snow_jacobian(l,k,n)    * Atmosphere_TL(n)%Cloud(4)%Water_Content(k) + &
                        iv%instid(inst)%graupel_jacobian(l,k,n) * Atmosphere_TL(n)%Cloud(5)%Water_Content(k) + &
                        iv%instid(inst)%hail_jacobian(l,k,n)    * Atmosphere_TL(n)%Cloud(6)%Water_Content(k)
                  end if
               end do
               RTSolution_TL(l,n)%brightness_temperature = RTSolution_TL(l,n)%brightness_temperature + &
                                       iv%instid(inst)%ps_jacobian(l,n) * Atmosphere_TL(n)%Level_Pressure(Atmosphere_TL(n)%n_Layers) 
            end do
         endif
      end do
     !$OMP END PARALLEL DO

      !-------------------------------------------------------------------
      ! [1.7] assign Hdx :
      !-------------------------------------------------------------------
      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         y%instid(inst)%tb(:,n) = y%instid(inst)%tb(:,n) + RTSolution_TL(:,n)%brightness_temperature
	 
         call CRTM_Atmosphere_Destroy( Atmosphere_TL(n) )
         IF ( CRTM_Atmosphere_Associated(Atmosphere_TL(n)) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",506, &
               (/"Error in deallocatting CRTM Atmosphere_TL Structure"/))
         end if

         call CRTM_Surface_Destroy(Surface_TL(n))
         IF ( CRTM_Surface_Associated(Surface_TL(n)) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",512, &
               (/"Error in deallocatting CRTM Surface_TL Structure"/))
         end if

      end do  ! end loop for pixels 

      deallocate (temperature)
      deallocate (absorber)
      deallocate (xb_q)
      deallocate (psfc)     

      if (crtm_cloud) then
         deallocate (qcw)
         deallocate (qci)
         deallocate (qrn)
         deallocate (qsn)
         deallocate (qgr)
      end if

      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
      !-------------------------------------------------------------------
      ! [2.0] Deallocating 1 structures
      ! -------------------------------------------------------------------
         call CRTM_Options_Destroy(Options(n))
         IF ( CRTM_Options_Associated(Options(n)) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",537, &
                         (/"Error in deallocatting CRTM Options Structure"/))
         end if

         call CRTM_Surface_Destroy(Surface(n))
         IF ( CRTM_Surface_Associated(Surface(n)) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",543, &
                         (/"Error in deallocatting CRTM Surface Structure"/))
         end if
      !-------------------------------------------------------------------
      ! [3.0] Deallocating 1 Atmosphere structures
      !-------------------------------------------------------------------
         call CRTM_Atmosphere_Destroy( Atmosphere(n) )
         IF ( CRTM_Atmosphere_Associated(Atmosphere(n)) ) THEN
            call da_error("da_transform_xtoy_crtm.inc",551, &
                         (/"Error in deallocatting CRTM Atmosphere Structure"/))
         end if

      end do 

      deallocate( RTSolution, RTSolution_TL, STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",559, &
                      (/"Error in deallocatting RTSolution"/))
      END IF

      deallocate( Atmosphere, Atmosphere_TL, STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",565, &
                      (/"Error in deallocatting Atmosphere"/))
      END IF

      deallocate( Surface, Surface_TL, STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",571, &
                      (/"Error in deallocatting Surface"/))
      END IF

      deallocate( Options, STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",577, &
                      (/"Error in deallocatting Options"/))
      END IF

      deallocate( GeometryInfo, STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm.inc",583, &
                      (/"Error in deallocatting GeometryInfo"/))
      END IF

   end do sensor_loop     ! end loop for sensor

   if (trace_use) call da_trace_exit("da_transform_xtoy_crtm")

end subroutine da_transform_xtoy_crtm

subroutine da_transform_xtoy_crtm_adj ( cv_size, cv, iv, jo_grad_y, jo_grad_x )

   !---------------------------------------------------------------------------
   ! PURPOSE: transform gradient from obs space to model grid space.
   !
   ! METHOD:  jo_grad_x = H^T jo_grad_y =  - H^T R^-1 ( d - H delta_x )
   !           1. input gradient in obs space and reference state of RTTOV
   !           2. call adjoint of RTM
   !           3. adjoint of interpolation from model grid to obs loc
   !
   !  HISTORY: 11/19/2006 - Creation                       Zhiquan Liu
   !           01/11/2008 - Add crtm_cloud                 Xiaoyan Zhang
   !           11/25/2008 - Zero Jacobian for top level    Tom Auligne
   !
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)            :: cv_size         ! Size of cv array.
   real, intent(inout)            :: cv(1:cv_size)   ! control variables.
   type (x_type), intent(inout)   :: jo_grad_x ! 
   type (y_type),  intent(in)     :: jo_grad_y ! H' delta_x
   type (iv_type), intent(in)     :: iv        ! O-B structure.

   integer, parameter             :: AIRS_Max_Channels = 281
   

   integer                        :: icld, jcld, i
   integer                        :: k  ! Index dimension.
   integer                        :: num_rad  ! Number of radiance obs
   integer                        :: inst, nchanl, n
   integer                        :: ipred, npred, gammapred, id
   real                           :: cv_local(1:cv_size)
   real, allocatable              :: q_ad(:,:)
   real, allocatable              :: t_ad(:,:)
   real, allocatable              :: p_ad(:)
   real, allocatable              :: xb_q(:,:)

!! for crtm_cloud
   real, allocatable              :: qcw_ad(:,:)
   real, allocatable              :: qci_ad(:,:)
   real, allocatable              :: qrn_ad(:,:)
   real, allocatable              :: qsn_ad(:,:)
   real, allocatable              :: qgr_ad(:,:)

!   type(infa_type), pointer       :: info                
!   integer                        :: i,j
!   real                           :: dx,dy,dxm,dym


   ! 1 local varaibles and types
   integer :: Allocate_Status
   integer :: n_layers, n_absorbers, n_clouds, n_aerosols
   type (CRTM_RTSolution_type ), allocatable :: RTSolution(:,:),RTSolution_AD(:,:)
   type (CRTM_Atmosphere_type ), allocatable :: Atmosphere(:), Atmosphere_AD(:)
   type (CRTM_Surface_type ),    allocatable :: Surface(:), Surface_AD(:)
   type (CRTM_Geometry_type ),   allocatable :: GeometryInfo(:)
   type (CRTM_Options_type ) , allocatable   :: Options(:)
   
   integer                        :: ts_index
   integer                        :: nclouds, ncv
   real, allocatable              :: cc_ad(:)
   real*8                         :: rad_clr            ! RT clear/cloudy radiances
   real, allocatable              :: rad_ovc(:)         ! RT overcast radiances
   real*8                         :: rad_ad, tb_ad

   ! Initializations for AIRS (MMR) Cloud Detection
   integer                        :: Band_Size(5), Bands(AIRS_Max_Channels,5) 
  
      Band_Size(1:5) = (/86, 0, 0, 16, 0 /)
      Bands(:,:)     = 0  
      Bands(1:Band_Size(1),1) = &
&    (/                                                 &              !&      1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
&                                                       &              !&     22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
&                                                       &              !&     52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
&                                                       &              !&     72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
&                     92,  93,  98,  99, 101, 104, 105, &              !&     83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
&     108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
&     138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
&     165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
&     180, 182, 185, 186, 190, 192,      198, 201, 204, &              !&     180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
&     207, 210,      215, 216,      221,      226, 227, &              !&     207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
&     232,                     252, 253, 256, 257, 261, &              !&     232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
&     262, 267, 272, 295, 299,      305,           310, &              !&     262, 267, 272, 295, 299, 300, 305, 308, 309, 310, &
&          321, 325, 333, 338, 355, 362, 375, 453, 475, &              !&     318, 321, 325, 333, 338, 355, 362, 375, 453, 475, &
&     484, 497, 528, 587, 672, 787, 791, 843, 870, 914, &
&     950 /)

!      Bands(1:Band_Size(2),2) = &
!&    (/ 1003, 1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082,  &
!&       1083, 1088, 1090, 1092, 1095, 1104, 1111, 1115, 1116, 1119,  &
!&       1120, 1123, 1130, 1138, 1142, 1178, 1199, 1206, 1221, 1237,  &
!&       1252, 1260, 1263, 1266, 1278, 1285 /)

!      Bands(1:Band_Size(3),3) = &
!&    (/       1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  !&    1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  
!&       1466,       1477,             1500, 1519,       1538, 1545, &  !&    1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &  
!&       1565, 1574, 1583, 1593,       1627, 1636,       1652, 1669, &  !&    1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, & 
!&                   1694, 1708,       1723, 1740, 1748,       1756, &  !&    1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
!&             1766, 1771, 1777,       1783, 1794, 1800,       1806, &  !&    1763, 1766, 1771, 1777, 1780, 1783, 1794, 1800, 1803, 1806, &
!&             1826, 1843  /)                                           !&    1812, 1826, 1843  /)

      Bands(1:Band_Size(4),4) = &
&    (/ 1852, 1865, 1866,       1868, 1869, 1872, 1873,       1876, &  !&    1852, 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, 
&             1881, 1882, 1883,                   1911, 1917, 1918, &  !&    1877, 1881, 1882, 1883, 1884, 1897, 1901, 1911, 1917, 1918, &
&                   1924, 1928        /)                               !&    1921, 1923, 1924, 1928, 1937  /)   

!      Bands(1:Band_Size(5),5) = &
!&    (/ 1938, 1939, 1941, 1946, 1947, 1948, 1958, 1971, 1973, 1988, &
!&       1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, 2104, &
!&       2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, &
!&       2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, 2134, &
!&       2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, 2234, &
!&       2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, 2355, &
!&       2363, 2370, 2371, 2377  /)  


!---------------------------------------------------------

   if ( iv%num_inst < 1 ) return

   if (trace_use) call da_trace_entry("da_transform_xtoy_crtm_adj")
   
   cv_local(:) = 0.0
   
   sensors_loop: do inst = 1, iv%num_inst                 ! loop for sensor

      num_rad = iv%instid(inst)%info%n2 - iv%instid(inst)%info%n1 + 1
      if ( num_rad < 1 ) cycle

      allocate (Atmosphere   (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                Atmosphere_AD(iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm_adj.inc",136, &
          (/"Error in allocatting Atmosphere"/))
      end if
      allocate (Surface      (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                Surface_AD   (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm_adj.inc",143, &
          (/"Error in allocatting Surface"/))
      end if
      allocate (GeometryInfo (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm_adj.inc",149, &
          (/"Error in allocatting GeometryInfo"/))
      end if
      allocate (Options      (iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
                STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm_adj.inc",155, &
          (/"Error in allocatting Options"/))
      end if

!----------------------------------------------------------------------------
! 1 allocation
!
! Atmosphere structure
!
      do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         n_layers    = (kte-kts)+1   ! number of vertical levels
         n_absorbers = 2
         n_aerosols  = 0
         n_clouds    = 0
         if ( crtm_cloud ) n_clouds = 6

         call CRTM_Atmosphere_Create( Atmosphere(n), &
                                      n_layers,      &
                                      n_absorbers,   &
                                      n_clouds,      &
                                      n_aerosols )
         IF ( .NOT. CRTM_Atmosphere_Associated( Atmosphere(n) ) ) THEN
            call da_error("da_transform_xtoy_crtm_adj.inc",178, &
              (/"Error in allocatting CRTM Atmosphere Structure"/))
         end if
         Atmosphere(n)%Absorber_ID(1)=H2O_ID
         Atmosphere(n)%Absorber_ID(2)=O3_ID
         Atmosphere(n)%Absorber_Units(1) = MASS_MIXING_RATIO_UNITS
         Atmosphere(n)%Absorber_Units(2) = VOLUME_MIXING_RATIO_UNITS
         Atmosphere(n)%Climatology=iv%instid(inst)%crtm_climat(n)

         if (crtm_cloud) then
            Atmosphere(n)%Cloud(1)%Type=WATER_CLOUD
            Atmosphere(n)%Cloud(2)%Type=ICE_CLOUD
            Atmosphere(n)%Cloud(3)%Type=RAIN_CLOUD
            Atmosphere(n)%Cloud(4)%Type=SNOW_CLOUD
            Atmosphere(n)%Cloud(5)%Type=GRAUPEL_CLOUD
            Atmosphere(n)%Cloud(6)%Type=HAIL_CLOUD
         end if
      end do

!-------------------------------------------------------------------------------

      allocate (q_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
      allocate (t_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
      allocate (p_ad(iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
      allocate (xb_q(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))

      q_ad = 0.0
      t_ad = 0.0
      p_ad = 0.0
      xb_q = 0.0

      if (crtm_cloud) then
         allocate (qcw_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qci_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qrn_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qsn_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         allocate (qgr_ad(kts:kte,iv%instid(inst)%info%n1:iv%instid(inst)%info%n2))
         qcw_ad = 0.0
         qci_ad = 0.0
         qrn_ad = 0.0
         qsn_ad = 0.0
         qgr_ad = 0.0
      end if

!      info => iv%instid(inst)%info

      nchanl    = ChannelInfo(inst)%n_channels
 
  ! Allocate forward model solution RTSolution array to number of channels
      allocate( RTSolution( ChannelInfo(inst)%n_Channels, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2 )   , &
                RTSolution_AD( ChannelInfo(inst)%n_Channels, iv%instid(inst)%info%n1:iv%instid(inst)%info%n2), &
               STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_transform_xtoy_crtm_adj.inc",231, &
                      (/"Error in allocatting RTSolution"/))
      END IF
 
      do n =  iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         call CRTM_Surface_Create( Surface(n),  &  ! Output
                                   nchanl )        ! Input
         IF ( .NOT. CRTM_Surface_Associated( Surface(n) ) ) THEN
           call da_error("da_transform_xtoy_crtm_adj.inc",240, &
             (/"Error in allocatting CRTM Surface Structure"/))
         end if
 
         ! 1 Options structure
         call CRTM_Options_Create( Options(n),  & ! Output
                                   nchanl )       ! Input
         IF ( .NOT. CRTM_Options_Associated( Options(n) ) ) THEN
           call da_error("da_transform_xtoy_crtm_adj.inc",248, &
             (/"Error in allocatting CRTM Options Structure"/))
         endif
         if ( use_antcorr(inst) ) Options(n)%Use_Antenna_Correction = .true.

         ! Gamma correction from VarBC
         !----------------------------
      end do

      do n =  iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         ! [1.0] Extract base state Atmosphere variables
         Atmosphere(n)%level_pressure(0) = iv%instid(inst)%pf(0,n)
         do k=1,Atmosphere(n)%n_layers
            Atmosphere(n)%pressure(k) = iv%instid(inst)%pm(k,n)
            Atmosphere(n)%level_pressure(k) = iv%instid(inst)%pf(k,n)
            Atmosphere(n)%temperature(k) = iv%instid(inst)%tm(k,n)
            Atmosphere(n)%absorber(k,1) = iv%instid(inst)%qm(k,n)
            xb_q(k,n) = 0.001*iv%instid(inst)%qm(k,n)/(1.0+iv%instid(inst)%qm(k,n)) ! specific humidity
         end do

         if (crtm_cloud) then
            do k=1,Atmosphere(n)%n_layers
               Atmosphere(n)%cloud(1)%water_content(k)=iv%instid(inst)%qcw(k,n)
               Atmosphere(n)%cloud(2)%water_content(k)=iv%instid(inst)%qci(k,n)
               Atmosphere(n)%cloud(3)%water_content(k)=iv%instid(inst)%qrn(k,n)
               Atmosphere(n)%cloud(4)%water_content(k)=iv%instid(inst)%qsn(k,n)
               Atmosphere(n)%cloud(5)%water_content(k)=iv%instid(inst)%qgr(k,n)
               Atmosphere(n)%cloud(6)%water_content(k)=iv%instid(inst)%qhl(k,n)
               Atmosphere(n)%cloud(1)%effective_radius(k)=iv%instid(inst)%rcw(k,n)
               Atmosphere(n)%cloud(2)%effective_radius(k)=iv%instid(inst)%rci(k,n)
               Atmosphere(n)%cloud(3)%effective_radius(k)=iv%instid(inst)%rrn(k,n)
               Atmosphere(n)%cloud(4)%effective_radius(k)=iv%instid(inst)%rsn(k,n)
               Atmosphere(n)%cloud(5)%effective_radius(k)=iv%instid(inst)%rgr(k,n)
               Atmosphere(n)%cloud(6)%effective_radius(k)=iv%instid(inst)%rhl(k,n)
            end do
         end if

         ! [1.1] User-supplied emissivity
         ! Options%emissivity_switch = 1
         ! Options%emissivity(1:Options%n_channels) = &
         !     iv%instid(inst)%emiss(1:Options%n_channels,n)

         ! [1.1] 1 Surface parameter data
         Surface(n)%Land_Coverage=iv%instid(inst)%land_coverage(n)
         Surface(n)%Water_Coverage=iv%instid(inst)%water_coverage(n)
         Surface(n)%Snow_Coverage=iv%instid(inst)%snow_coverage(n)
         Surface(n)%Ice_Coverage=iv%instid(inst)%ice_coverage(n)

         if (Surface(n)%Land_Coverage > 0.0) then
            Surface(n)%Land_Type=11 !GRASS_SOIL
            Surface(n)%Land_Temperature=iv%instid(inst)%ts(n)      ! K
            Surface(n)%Soil_Moisture_Content=iv%instid(inst)%smois(n)  !0.05    ! volumetric water content (g/cm**3)
            !Surface(n)%Canopy_Water_Content=0.05      ! gravimetric water content
            Surface(n)%Vegetation_Fraction=iv%instid(inst)%vegfra(n)
            Surface(n)%Soil_Temperature=iv%instid(inst)%tslb(n)
         end if
         if (Surface(n)%Water_Coverage > 0.0) then
            !Surface(n)%Water_Type=SEA_WATER          ! (Currently NOT used)
            Surface(n)%Water_Temperature=iv%instid(inst)%ts(n)     ! K
            Surface(n)%Wind_Speed=sqrt((iv%instid(inst)%u10(n))**2+ &
                               (iv%instid(inst)%v10(n))**2)  ! m/sec
            !surface(n)%Wind_Direction=0.0            ! NOT used
            Surface(n)%Salinity=33.0                  ! ppmv
         end if
         if (Surface(n)%Snow_Coverage > 0.0) then
            Surface(n)%Snow_Type=2 !NEW_SNOW
            Surface(n)%Snow_Temperature=iv%instid(inst)%ts(n)      ! K
            Surface(n)%Snow_Depth=iv%instid(inst)%snowh(n)         ! mm
            !Surface(n)%Snow_Density=0.2               ! g/cm**3
            !Surface(n)%Snow_Grain_Size=2.0            ! mm
         end if
         if (Surface(n)%Ice_Coverage > 0.0) then
            !Surface(n)%Ice_Type=FRESH_ICE             ! NO Table offered, single example is FRESH_ICE
            Surface(n)%Ice_Temperature=iv%instid(inst)%ts(n)       ! K
            Surface(n)%Ice_Thickness=10.0              ! mm
            !Surface(n)%Ice_Density=0.9                ! g/cm**3
            !Surface(n)%Ice_Roughness=0.0               ! NO Table offered, single example is ZERO
         end if
         Surface(n)%SensorData%n_channels        = nchanl
         Surface(n)%SensorData%Sensor_Channel    = ChannelInfo(inst)%Sensor_Channel
         Surface(n)%SensorData%Sensor_Id         = ChannelInfo(inst)%Sensor_Id
         Surface(n)%SensorData%WMO_Satellite_Id  = ChannelInfo(inst)%WMO_Satellite_Id
         Surface(n)%SensorData%WMO_Sensor_Id     = ChannelInfo(inst)%WMO_Sensor_Id
         Surface(n)%SensorData%Tb(1:nchanl) = iv%instid(inst)%tb_inv(1:nchanl,n) + &
                                              iv%instid(inst)%tb_xb(1:nchanl,n)
     
         ! -- Copy the adjoint atmosphere structure
         Atmosphere_AD(n) = Atmosphere(n)

         ! -- Copy the adjoint surface structure
         Surface_AD(n) = Surface(n)

         ! -- Zero the Adjoint outputs 
         ! Important: adjoint variables must be initialized
         call CRTM_Atmosphere_Zero( Atmosphere_AD(n) )
         call CRTM_Surface_Zero( Surface_AD(n) )

         ! [1.2] 1 GeometryInfo Structure
         GeometryInfo(n)%Sensor_Zenith_Angle=iv%instid(inst)%satzen(n)
         GeometryInfo(n)%Source_Zenith_Angle=iv%instid(inst)%solzen(n)
         GeometryInfo(n)%iFOV=iv%instid(inst)%scanpos(n)
   !     GeometryInfo(n)%Satellite_Height=830.0
   !     GeometryInfo(n)%Sensor_Scan_Angle=
   !     GeometryInfo(n)%Sensor_Zenith_Angle=
   !     GeometryInfo(n)%Sensor_Scan_Angle=
   !     GeometryInfo(n)%Source_Zenith_Angle=
	 
         ! [1.3] assign tb = R^-1 Re :

         do i = 1, ChannelInfo(inst)%n_Channels
            RTSolution_AD(i,n)%brightness_temperature = jo_grad_y%instid(inst)%tb(i,n)
            RTSolution_AD(i,n)%radiance = 0.0 ! must assign zero, since each call of AD model will return non-zero value
         end do

      end do

      ! [1.4] Call CRTM_AD model
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( n )
      do n =  iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
         if (.not. use_crtm_kmatrix) then
            call da_crtm_ad (1, nchanl, 1, Atmosphere(n),   &
                               Surface(n),      &
                               RTSolution_AD(:,n),&
                               GeometryInfo(n), &
                               ChannelInfo(inst:inst),  &
                               Atmosphere_AD(n),&
                               Surface_AD(n),   &
                               RTSolution(:,n),   &
                               Options(n))
         else
            do i = 1, ChannelInfo(inst)%n_Channels
               Atmosphere_AD(n)%Level_Pressure(Atmosphere(n)%n_layers) = &
                       Atmosphere_AD(n)%Level_Pressure(Atmosphere(n)%n_layers) + &
                       iv%instid(inst)%ps_jacobian(i,n) * RTSolution_AD(i,n)%brightness_temperature   
            end do

            do k = kts , kte
               do i = 1, ChannelInfo(inst)%n_Channels
                  Atmosphere_AD(n)%Temperature(k) = Atmosphere_AD(n)%Temperature(k) + &
                                iv%instid(inst)%t_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                  Atmosphere_AD(n)%absorber(k,1) = Atmosphere_AD(n)%absorber(k,1) + &
                                iv%instid(inst)%q_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                  if (crtm_cloud) then
                     Atmosphere_AD(n)%Cloud(1)%Water_Content(k) = Atmosphere_AD(n)%Cloud(1)%Water_Content(k) + &
                        iv%instid(inst)%water_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                     Atmosphere_AD(n)%Cloud(2)%Water_Content(k) = Atmosphere_AD(n)%Cloud(2)%Water_Content(k) + &
                        iv%instid(inst)%ice_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                     Atmosphere_AD(n)%Cloud(3)%Water_Content(k) = Atmosphere_AD(n)%Cloud(3)%Water_Content(k) + &
                        iv%instid(inst)%rain_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                     Atmosphere_AD(n)%Cloud(4)%Water_Content(k) = Atmosphere_AD(n)%Cloud(4)%Water_Content(k) + &
                        iv%instid(inst)%snow_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                     Atmosphere_AD(n)%Cloud(5)%Water_Content(k) = Atmosphere_AD(n)%Cloud(5)%Water_Content(k) + &
                        iv%instid(inst)%graupel_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                     Atmosphere_AD(n)%Cloud(6)%Water_Content(k) = Atmosphere_AD(n)%Cloud(6)%Water_Content(k) + &
                        iv%instid(inst)%hail_jacobian(i,k,n) * RTSolution_AD(i,n)%brightness_temperature
                  end if
               end do
            end do
         endif
      end do
      !$OMP END PARALLEL DO

      do n =  iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

      ! Skin Temperature
      !-----------------
	 if (use_satcv(1)) then
            ts_index = iv%instid(inst)%cv_index(n)%ts 	       
  	    if (Surface(n)%Land_Coverage  > 0.0) cv(ts_index) = Surface_AD(n)%Land_Temperature
            if (Surface(n)%Water_Coverage > 0.0) cv(ts_index) = Surface_AD(n)%Water_Temperature
            if (Surface(n)%Snow_Coverage  > 0.0) cv(ts_index) = Surface_AD(n)%Snow_Temperature
            if (Surface(n)%Ice_Coverage   > 0.0) cv(ts_index) = Surface_AD(n)%Ice_Temperature
	 end if 

      ! Cloud cover(s)
      !---------------
	 if (use_satcv(2)) then
	    nclouds = iv%instid(inst)%cv_index(n)%nclouds
	    ncv     = iv%instid(inst)%cv_index(n)%ncv
	    allocate (rad_ovc (nclouds))
	    allocate (cc_ad   (ncv))
	    do k = 1, nchanl
               if (ALL(iv%instid(inst)%ichan(k) /= Bands(:,1))) cycle   ! Only Channels in Band 1 	       
               rad_clr      = iv%instid(inst)%rad_xb(k,n)             
               do icld = kte-nclouds+1, kte
                  rad_ovc(icld)   = iv%instid(inst)%rad_ovc(k,icld,n)
               end do
	       rad_ad       = 0.0
	       tb_ad        = jo_grad_y%instid(inst)%tb(k,n)
	       call CRTM_Planck_Temperature_AD(inst,k,rad_clr,tb_ad,rad_ad)

 	       do icld = 1, nclouds
	          cc_ad(icld)        = rad_ad * (rad_ovc(icld)-rad_clr)
               end do
!	       !---------------------------------------------------------------
!              ! Change of variable (preconditioning) 
!              !---------------------------------------------------------------
!	       do icld = 1, ncv
!	          cc_ad(icld) = SUM(rad_ad * (rad_ovc(:)-rad_clr) * &
!		                iv%instid(inst)%cv_index(n)%vtox(icld,:) )
!	       end do

 	       do icld = 1, ncv
	          cv(iv%instid(inst)%cv_index(n)%cc(icld)) = cv(iv%instid(inst)%cv_index(n)%cc(icld)) + cc_ad(icld)
               end do
	    end do	
	    deallocate(rad_ovc, cc_ad)
	 end if 

      ! [1.5] Scale transformation and fill zero for no-control variable

        ! Convert cloud content unit from kg/kg to kg/m^2
         if (crtm_cloud) then
            do k=kts,kte
               do icld=1,Atmosphere(n)%n_Clouds
                  Atmosphere_AD(n)%Cloud(icld)%Water_Content(k) = &
		        Atmosphere_AD(n)%Cloud(icld)%Water_Content(k) * &
		        (Atmosphere(n)%Level_Pressure(k)- Atmosphere(n)%Level_Pressure(k-1))*100.0/gravity
               enddo
            enddo
         end if

       ! [1.6] Adjoint of Interpolate horizontally from ob to grid:

         if (crtm_cloud) then
            do k=kts,kte ! from bottom to top
               qcw_ad(k,n)=Atmosphere_AD(n)%Cloud(1)%Water_Content(kte-k+1)
               qci_ad(k,n)=Atmosphere_AD(n)%Cloud(2)%Water_Content(kte-k+1)
               qrn_ad(k,n)=Atmosphere_AD(n)%Cloud(3)%Water_Content(kte-k+1)
               qsn_ad(k,n)=Atmosphere_AD(n)%Cloud(4)%Water_Content(kte-k+1)
               qgr_ad(k,n)=Atmosphere_AD(n)%Cloud(5)%Water_Content(kte-k+1)
            end do
         end if

         do k=kts,kte-1 ! from bottom to top.  Zero Jacobian for top level
            if (atmosphere(n)%pressure(kte-k+1) >= 75.0) &        ! Zero Jacobian for top level(s)
               q_ad(k,n) = Atmosphere_AD(n)%Absorber(kte-k+1,1) * 1000.0 /  &
                           (1.0-xb_q(kte-k+1,n))**2               ! in g/kg
               t_ad(k,n) = Atmosphere_AD(n)%Temperature(kte-k+1)
         end do

         p_ad(n) = Atmosphere_AD(n)%Level_Pressure(atmosphere(n)%n_layers) * 0.01   ! in hPa

         call CRTM_Atmosphere_Destroy( Atmosphere_AD(n) )
         IF ( CRTM_Atmosphere_Associated(Atmosphere_AD(n)) ) THEN
            call da_error("da_transform_xtoy_crtm_adj.inc",511, &
                         (/"Error in deallocatting CRTM Atmosphere_AD Structure"/))
         end if

         call CRTM_Surface_Destroy(Surface_AD(n))
         IF ( CRTM_Surface_Associated(Surface_AD(n)) ) THEN
            call da_error("da_transform_xtoy_crtm_adj.inc",517, &
                         (/"Error in deallocatting CRTM Surface_AD Structure"/))
         end if

      end do       !  end loop for pixels
	 
         ! [1.7] Gamma correction to VarBC
      
        if (crtm_cloud) then
           call da_interp_lin_2d_adj_partial(jo_grad_x%qcw(:,:,kts:kte),iv%instid(inst)%info, kts,kte, qcw_ad)
           call da_interp_lin_2d_adj_partial(jo_grad_x%qci(:,:,kts:kte),iv%instid(inst)%info, kts,kte, qci_ad)
           call da_interp_lin_2d_adj_partial(jo_grad_x%qrn(:,:,kts:kte),iv%instid(inst)%info, kts,kte, qrn_ad)
           call da_interp_lin_2d_adj_partial(jo_grad_x%qsn(:,:,kts:kte),iv%instid(inst)%info, kts,kte, qsn_ad)
           call da_interp_lin_2d_adj_partial(jo_grad_x%qgr(:,:,kts:kte),iv%instid(inst)%info, kts,kte, qgr_ad)
        endif

         call da_interp_lin_2d_adj_partial(jo_grad_x%t(:,:,kts:kte),    iv%instid(inst)%info, kts,kte, t_ad)
         call da_interp_lin_2d_adj_partial(jo_grad_x%q(:,:,kts:kte),    iv%instid(inst)%info, kts,kte, q_ad)
         call da_interp_lin_2d_adj_partial(jo_grad_x%psfc, iv%instid(inst)%info, 1,1, p_ad)

         deallocate (q_ad)
         deallocate (t_ad)
         deallocate (p_ad)
         deallocate (xb_q)

         if (crtm_cloud) then
            deallocate (qcw_ad)
            deallocate (qci_ad)
            deallocate (qrn_ad)
            deallocate (qsn_ad)
            deallocate (qgr_ad)
         end if

      !-------------------------------------------------------------------
      ! [2.0] Deallocating 1 structures
      !-------------------------------------------------------------------
         deallocate( RTSolution, RTSolution_AD, STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) then
            call da_error("da_transform_xtoy_crtm_adj.inc",575, &
              (/"Error in deallocatting RTSolution"/))
         END IF

      do n =  iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel

         call CRTM_Options_Destroy(Options(n))
         IF ( CRTM_Options_Associated(Options(n)) ) THEN
            call da_error("da_transform_xtoy_crtm_adj.inc",583, &
               (/"Error in deallocatting CRTM Options Structure"/))
         end if

         call CRTM_Surface_Destroy(Surface(n))
         IF ( CRTM_Surface_Associated(Surface(n)) ) THEN
            call da_error("da_transform_xtoy_crtm_adj.inc",589, &
               (/"Error in deallocatting CRTM Surface Structure"/))
         end if

      !-------------------------------------------------------------------
      ! [3.0] Deallocating 1 Atmosphere structures
      !-------------------------------------------------------------------
         call CRTM_Atmosphere_Destroy( Atmosphere(n) )
         IF ( CRTM_Atmosphere_Associated(Atmosphere(n)) ) THEN
            call da_error("da_transform_xtoy_crtm_adj.inc",598, &
                         (/"Error in deallocatting CRTM Atmosphere Structure"/))
         end if

      end do

         deallocate( Atmosphere, Atmosphere_AD, STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) then
            call da_error("da_transform_xtoy_crtm_adj.inc",606, &
              (/"Error in deallocatting Atmosphere"/))
         END IF

         deallocate( Surface, Surface_AD, STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) then
            call da_error("da_transform_xtoy_crtm_adj.inc",612, &
              (/"Error in deallocatting Surface"/))
         END IF

         deallocate( Options, STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) then
            call da_error("da_transform_xtoy_crtm_adj.inc",618, &
              (/"Error in deallocatting Options"/))
         END IF

         deallocate( GeometryInfo, STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) then
            call da_error("da_transform_xtoy_crtm_adj.inc",624, &
              (/"Error in deallocatting GeometryInfo"/))
         END IF

   end do sensors_loop       ! end loop for sensor

       
   if (trace_use) call da_trace_exit("da_transform_xtoy_crtm_adj")
 
end subroutine da_transform_xtoy_crtm_adj

subroutine da_get_innov_vector_crtm ( it, grid, ob, iv )

   !---------------------------------------------------------------------------
   !  PURPOSE: Calculate innovation vector for radiance data.
   !
   !  METHOD:  d = y - H(x)
   !       1. interpolate grid%xb to obs location
   !       2. call foreward RTM to get simulated bright temperature 
   !       3. obs BT - simulated BT
   !  HISTORY: 12/15/2008 effective radius unit is micron     Zhiquan Liu
   !---------------------------------------------------------------------------

   implicit none
   
   integer,           intent(in)    :: it       ! External iteration.
   type (domain),     intent(in)    :: grid     ! first guess state.
   type (y_type),     intent(inout) :: ob       ! Observation structure.
   type (iv_type),    intent(inout) :: iv       ! O-B structure.

   integer, parameter             :: AIRS_Max_Channels = 281

   integer :: n, icld  ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: l        ! Index dimension.
   integer :: num_levs ! Number of obs levels.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   integer :: alloc_status(40)

   real, allocatable :: model_u10(:)
   real, allocatable :: model_v10(:)
   real, allocatable :: model_psfc(:)
   real              :: model_ptop
   real, allocatable :: model_ts(:)
   real, allocatable :: model_elv(:)
   real, allocatable :: model_smois(:)
   real, allocatable :: model_tslb(:)
   real, allocatable :: model_snowh(:)
   real, allocatable :: model_vegfra(:)
   real    :: model_isltyp, model_ivgtyp
   integer :: model_isflg

   real    :: model_qcw(kms:kme)
   real    :: model_rho(kms:kme)
   real, allocatable :: model_snow(:)  ! snow water equivalent, different from model_snowh,
                                       ! used in calculating reff_water
   real, allocatable :: em_mspps(:)    ! emissivity caluclated using MSPPS algorithm
   real              :: ts_mspps       ! surface temperature calcualted using MSPPS algorithm
   character(len=19) :: date_char

   integer :: month, crtm_climat

   ! real    :: model_tm(kms:kme)

   integer :: inst, nchanl, n1,n2
   integer :: ipred, npred, gammapred

   ! variables for computing clwp
   real    :: clw(kms:kme), dpf(kms:kme)
   real    :: clwp

   real    :: total_od

   ! 1 local varaibles and types
   integer :: Allocate_Status
   integer :: n_layers, n_absorbers, n_clouds, n_aerosols
!! for ecmwf cloud
   real    :: t_tropp(ims:ime,jms:jme),temp1
   integer :: k_tropp(ims:ime,jms:jme)
   type (CRTM_RTSolution_type), allocatable :: RTSolution(:,:)
   type (CRTM_Atmosphere_type)   :: Atmosphere(1)
   type (CRTM_Surface_type)      :: Surface(1)
   type (CRTM_Geometry_type)     :: GeometryInfo(1)
   type (CRTM_Options_type)      :: Options(1)
   type (CRTM_RTSolution_type), allocatable :: RTSolution_K(:,:)
   type (CRTM_Atmosphere_type), allocatable :: Atmosphere_K(:,:)
   type (CRTM_Surface_type),    allocatable :: Surface_K(:,:)

   real :: t(1), a(1), p(1)

   integer, allocatable :: wrf_to_crtm_mw(:)
   integer :: n_vegtype

!! for crtm cloud
   real :: qcw(1),qrn(1), qci(1),qsn(1),qgr(1)
   
   ! Initializations for AIRS (MMR) Cloud Detection
   integer           :: ilev, jlev, nclouds
   real, allocatable :: hessian(:,:)
   real*8, allocatable :: eignvec(:,:), eignval(:)
   real              :: rad_clr, rad_ovc_ilev, rad_ovc_jlev
   
   integer           :: Band_Size(5), Bands(AIRS_Max_Channels,5) 
  
      Band_Size(1:5) = (/86, 0, 0, 16, 0 /)
      Bands(:,:)     = 0  
      Bands(1:Band_Size(1),1) = &
&    (/                                                 &              !&      1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
&                                                       &              !&     22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
&                                                       &              !&     52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
&                                                       &              !&     72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
&                     92,  93,  98,  99, 101, 104, 105, &              !&     83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
&     108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
&     138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
&     165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
&     180, 182, 185, 186, 190, 192,      198, 201, 204, &              !&     180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
&     207, 210,      215, 216,      221,      226, 227, &              !&     207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
&     232,                     252, 253, 256, 257, 261, &              !&     232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
&     262, 267, 272, 295, 299,      305,           310, &              !&     262, 267, 272, 295, 299, 300, 305, 308, 309, 310, &
&          321, 325, 333, 338, 355, 362, 375, 453, 475, &              !&     318, 321, 325, 333, 338, 355, 362, 375, 453, 475, &
&     484, 497, 528, 587, 672, 787, 791, 843, 870, 914, &
&     950 /)

!      Bands(1:Band_Size(2),2) = &
!&    (/ 1003, 1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082,  &
!&       1083, 1088, 1090, 1092, 1095, 1104, 1111, 1115, 1116, 1119,  &
!&       1120, 1123, 1130, 1138, 1142, 1178, 1199, 1206, 1221, 1237,  &
!&       1252, 1260, 1263, 1266, 1278, 1285 /)

!      Bands(1:Band_Size(3),3) = &
!&    (/       1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  !&    1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &  
!&       1466,       1477,             1500, 1519,       1538, 1545, &  !&    1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &  
!&       1565, 1574, 1583, 1593,       1627, 1636,       1652, 1669, &  !&    1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, & 
!&                   1694, 1708,       1723, 1740, 1748,       1756, &  !&    1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
!&             1766, 1771, 1777,       1783, 1794, 1800,       1806, &  !&    1763, 1766, 1771, 1777, 1780, 1783, 1794, 1800, 1803, 1806, &
!&             1826, 1843  /)                                           !&    1812, 1826, 1843  /)

      Bands(1:Band_Size(4),4) = &
&    (/ 1852, 1865, 1866,       1868, 1869, 1872, 1873,       1876, &  !&    1852, 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, 
&             1881, 1882, 1883,                   1911, 1917, 1918, &  !&    1877, 1881, 1882, 1883, 1884, 1897, 1901, 1911, 1917, 1918, &
&                   1924, 1928        /)                               !&    1921, 1923, 1924, 1928, 1937  /)   

!      Bands(1:Band_Size(5),5) = &
!&    (/ 1938, 1939, 1941, 1946, 1947, 1948, 1958, 1971, 1973, 1988, &
!&       1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, 2104, &
!&       2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, &
!&       2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, 2134, &
!&       2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, 2234, &
!&       2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, 2355, &
!&       2363, 2370, 2371, 2377  /)  

   ! WHY? use argument it
   if (it==0) then; write(unit=stdout,fmt='(A)') "WHY? have argument it to"//"da_get_innov_vector_crtm.inc"; end if

   alloc_status (:) = 0

   if (trace_use) call da_trace_entry("da_get_innov_vector_crtm")

   ! 1 allocation

   ! 1 Atmosphere
   n_layers    = (kte-kts)+1   ! number of vertical levels
   n_absorbers = 2
   n_aerosols  = 0
   n_clouds    = 0
   if ( crtm_cloud ) n_clouds = 6

   call CRTM_Atmosphere_Create ( Atmosphere(1), &
                                 n_layers,      &
                                 n_absorbers,   &
                                 n_clouds,      &
                                 n_aerosols )
   IF ( .NOT. CRTM_Atmosphere_Associated( Atmosphere(1) ) ) THEN 
       call da_error("da_get_innov_vector_crtm.inc",164, &
         (/"Error in allocating CRTM Atmosphere Structure"/))
   END IF 

   Atmosphere(1)%Absorber_ID(1)=H2O_ID
   Atmosphere(1)%Absorber_ID(2)=O3_ID
   Atmosphere(1)%Absorber_Units(1) = MASS_MIXING_RATIO_UNITS
   Atmosphere(1)%Absorber_Units(2) = VOLUME_MIXING_RATIO_UNITS

   if (crtm_cloud) then
      Atmosphere(1)%Cloud(1)%Type=WATER_CLOUD
      Atmosphere(1)%Cloud(2)%Type=ICE_CLOUD
      Atmosphere(1)%Cloud(3)%Type=RAIN_CLOUD
      Atmosphere(1)%Cloud(4)%Type=SNOW_CLOUD
      Atmosphere(1)%Cloud(5)%Type=GRAUPEL_CLOUD
      Atmosphere(1)%Cloud(6)%Type=HAIL_CLOUD
   end if

   select case (trim(grid%mminlu))
      case ('USGS')
         n_vegtype = USGS_n_type
      case ('MODIFIED_IGBP_MODIS_NOAH')
         n_vegtype = IGBP_n_type
      case default
         call da_error("da_get_innov_vector_crtm.inc",188,(/"Unknown land use type"/))
   end select
   allocate(wrf_to_crtm_mw(n_vegtype))
   if ( n_vegtype == USGS_n_type ) then
         wrf_to_crtm_mw = usgs_to_crtm_mw
   else if ( n_vegtype == IGBP_n_type ) then
         wrf_to_crtm_mw = igbp_to_crtm_mw
   end if

   !------------------------------------------------------
   ! [1.0] calculate the background bright temperature
   !-------------------------------------------------------
   do inst = 1, iv%num_inst                 ! loop for sensor
      ! if ( iv%instid(inst)%num_rad < 1 ) cycle
      if (iv%instid(inst)%info%n2 < iv%instid(inst)%info%n1) cycle
      num_levs  = kte-kts+1 
      nchanl    = ChannelInfo(inst)%n_channels

      if ( ChannelInfo(inst)%Sensor_Type == 2 ) then  ! INFRARED_SENSOR=2
         if ( index(grid%mminlu, trim(CRTM_IRlandCoeff_Classification())) <= 0 ) then
            call da_error("da_get_innov_vector_crtm.inc",208, &
              (/"modify crtm_irland_coef to use the coef file that is consistent with your model land use type"/))
         end if
      end if

      ! Allocate forward model solution RTSolution array to number of channels
      allocate (RTSolution(ChannelInfo(inst)%n_Channels,1), STAT = Allocate_Status )
      if ( Allocate_Status /= 0 ) then
         call da_error("da_get_innov_vector_crtm.inc",216, (/"Error in allocating RTSolution"/))
      end if

      call CRTM_RTSolution_Create( RTSolution , &                  !Output
                                   Atmosphere(1)%n_layers )        !Input
      if ( .NOT.(ANY(crtm_rtsolution_associated(RTSolution))) ) then
         call da_error("da_get_innov_vector_crtm.inc",222, &
            (/"Error in allocating CRTM RTSolution Stucture"/))
      end if
            
      if (use_crtm_kmatrix) then
         allocate (RTSolution_K(ChannelInfo(inst)%n_Channels,1), &
                   Atmosphere_K(ChannelInfo(inst)%n_Channels,1), &
                   Surface_K(ChannelInfo(inst)%n_Channels,1),    &
                   STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) then
            call da_error("da_get_innov_vector_crtm.inc",232, (/"Error in allocating RTSolution_K"/))
         end if
         call CRTM_RTSolution_Create( RTSolution_K , &                !Output
                                      Atmosphere(1)%n_layers )        !Input
         if ( .NOT.(ANY(crtm_rtsolution_associated(RTSolution_K))) ) then
            call da_error("da_get_innov_vector_crtm.inc",237, &
               (/"Error in allocating CRTM RTSolution_K Stucture"/))
         end if
      end if
      
      call CRTM_Surface_Create ( Surface(1), & ! Output
                                 nchanl )      ! Input
      IF ( .NOT. CRTM_Surface_Associated( Surface(1) ) ) THEN 
         call da_error("da_get_innov_vector_crtm.inc",245, &
            (/"Error in allocating CRTM Surface Structure"/))
      end if

      ! 1 Options structure
      call CRTM_Options_Create( Options(1), &  ! Output
                                nchanl )       ! Input
      IF ( .NOT. CRTM_Options_Associated( Options(1) ) ) THEN 
        call da_error("da_get_innov_vector_crtm.inc",253, &
          (/"Error in allocatting CRTM Options Structure"/))
      endif
      if ( use_antcorr(inst) ) Options(1)%Use_Antenna_Correction = .true.

      ! do n= 1, iv%instid(inst)%num_rad           ! loop for pixel

      n1 = iv%instid(inst)%info%n1
      n2 = iv%instid(inst)%info%n2

      allocate (model_u10(n1:n2))
      allocate (model_v10(n1:n2))
      allocate (model_psfc(n1:n2))
      allocate (model_ts(n1:n2))
      allocate (model_elv(n1:n2))
      allocate (model_smois(n1:n2))
      allocate (model_tslb(n1:n2))
      allocate (model_snowh(n1:n2))
      allocate (model_snow(n1:n2))
      allocate (model_vegfra(n1:n2))

      model_u10(:)    = 0.0
      model_v10(:)    = 0.0
      model_psfc(:)   = 0.0
      model_ts(:)     = 0.0
      model_elv(:)    = 0.0
      model_smois(:)  = 0.0
      model_tslb(:)   = 0.0
      model_snowh(:)  = 0.0
      model_snow(:)   = 0.0
      model_vegfra(:) = 0.0

      ! Gamma correction from VarBC
      !----------------------------

!      ! Allocate Overcast Radiances for  Cloud Detection (MMR)
      !----------------------------------------------------------------
!#ifdef CRTM_MODIF
!      if (use_clddet_mmr .or. use_clddet_ecmwf ) then
!         do i = 1, nchanl
!            allocate(RTSolution(i,1)%Overcast(Atmosphere(1)%n_Layers))
!         end do		         
!      end if   
!#endif

      if ( use_mspps_emis(inst) ) allocate(em_mspps(nchanl))

      pixel_loop: do n=n1,n2
         ! if ( n > iv%instid(inst)%num_rad ) exit

         ! [1.1] Get horizontal interpolation weights:

         i = iv%instid(inst)%info%i(1,n)
         j = iv%instid(inst)%info%j(1,n)
         dx = iv%instid(inst)%info%dx(1,n)
         dy = iv%instid(inst)%info%dy(1,n)
         dxm = iv%instid(inst)%info%dxm(1,n)
         dym = iv%instid(inst)%info%dym(1,n)

         ! determine 1 climatology based on latitude and month
         date_char=iv%instid(inst)%info%date_char(n)
         read(unit=date_char(6:7), fmt='(i2.2)') month
         call da_det_crtm_climat(iv%instid(inst)%info%lat(1,n), month, crtm_climat)
         iv%instid(inst)%crtm_climat(n) = crtm_climat
         Atmosphere(1)%Climatology = iv%instid(inst)%crtm_climat(n)

         ! determine surface type of obs location
         !-----------------------------------------
         call da_detsurtyp ( grid%xb%snow, grid%xb%xice, grid%xb%landmask,  &
            grid%xb%ivgtyp, grid%xb%isltyp, &
            ims, ime, jms, jme, &
            i, j, dx, dy, dxm, dym, &
            model_isflg,model_ivgtyp, model_isltyp, &
            Surface(1)%Water_Coverage, Surface(1)%Ice_Coverage, &
            Surface(1)%Land_Coverage, Surface(1)%Snow_Coverage )

         call da_interp_2d_partial (grid%xb%snow, iv%instid(inst)%info, 1, n, n, model_snow(n:n))

         ! [1.2] Interpolate horizontally to ob:
         do k=kts,kte ! from bottom to top
            call da_interp_2d_partial (grid%xb%p(:,:,k), iv%instid(inst)%info,k,n,n,p)
            call da_interp_2d_partial (grid%xb%t(:,:,k), iv%instid(inst)%info,k,n,n,t)
            call da_interp_2d_partial (grid%xb%q(:,:,k), iv%instid(inst)%info,k,n,n,a) 

            Atmosphere(1)%Pressure(kte-k+1)    = 0.01   * p(1)  ! convert Pa to hPa
            Atmosphere(1)%Temperature(kte-k+1) = t(1)
            Atmosphere(1)%Absorber(kte-k+1,1)  = 1000.0 * a(1)/(1.0-a(1))  ! in g/kg mixing ratio

            ! NOTE: WRF high-level q values seems too big, replaced by constants
            if (p(1)*0.01 < 75.0) Atmosphere(1)%Absorber(kte-k+1,1) = 0.001

            call da_interp_2d_partial (grid%xb%qcw(:,:,k), iv%instid(inst)%info,k,n,n, model_qcw(kte-k+1:kte-k+1))
            
            if (crtm_cloud) then

               call da_interp_2d_partial (grid%xb%qci(:,:,k), iv%instid(inst)%info,k,n,n,qci)

               call da_interp_2d_partial (grid%xb%qrn(:,:,k), iv%instid(inst)%info,k,n,n,qrn)
 
               call da_interp_2d_partial (grid%xb%qsn(:,:,k), iv%instid(inst)%info, k,n,n,qsn)
 
               call da_interp_2d_partial (grid%xb%qgr(:,:,k), iv%instid(inst)%info, k,n,n,qgr)

               Atmosphere(1)%Cloud(1)%Water_Content(kte-k+1)=model_qcw(kte-k+1)
               Atmosphere(1)%Cloud(2)%Water_Content(kte-k+1)=qci(1)
               Atmosphere(1)%Cloud(3)%Water_Content(kte-k+1)=qrn(1)
               Atmosphere(1)%Cloud(4)%Water_Content(kte-k+1)=qsn(1)
               Atmosphere(1)%Cloud(5)%Water_Content(kte-k+1)=qgr(1)
               Atmosphere(1)%Cloud(6)%Water_Content(kte-k+1)=0.0

               call da_interp_2d_partial (grid%xb%rho(:,:,k), iv%instid(inst)%info, k,n,n, &
                  model_rho(k:k) )

               call da_cld_eff_radius(Atmosphere(1)%Temperature(kte-k+1),model_rho(k),&
                                      Atmosphere(1)%Cloud(2)%Water_Content(kte-k+1),  &  !qci
                                      Atmosphere(1)%Cloud(3)%Water_Content(kte-k+1),  &  !qrn
                                      Atmosphere(1)%Cloud(4)%Water_Content(kte-k+1),  &  !qsn
                                      Atmosphere(1)%Cloud(5)%Water_Content(kte-k+1),  &  !qgr
                                      model_snow(n),                                  &
                                      Surface(1)%Ice_Coverage, Surface(1)%Land_Coverage, 1, &
                                      Atmosphere(1)%Cloud(1)%Effective_Radius(kte-k+1), &
                                      Atmosphere(1)%Cloud(2)%Effective_Radius(kte-k+1), &
                                      Atmosphere(1)%Cloud(3)%Effective_Radius(kte-k+1), &
                                      Atmosphere(1)%Cloud(4)%Effective_Radius(kte-k+1), &
                                      Atmosphere(1)%Cloud(5)%Effective_Radius(kte-k+1) )

               ! reset the da_cld_eff_radius calcualted effective radius to constants if desired

               Atmosphere(1)%Cloud(1)%Effective_Radius(kte-k+1)=10.0  ! in micron
               Atmosphere(1)%Cloud(2)%Effective_Radius(kte-k+1)=30.0  ! in micron
               !Atmosphere(1)%Cloud(3)%Effective_Radius(kte-k+1)=300
               !Atmosphere(1)%Cloud(4)%Effective_Radius(kte-k+1)=600
               !Atmosphere(1)%Cloud(5)%Effective_Radius(kte-k+1)=600
               Atmosphere(1)%Cloud(6)%Effective_Radius(kte-k+1)=600

            end if
         end do

         call da_interp_2d_partial (grid%xb%u10,  iv%instid(inst)%info, 1, n, n, model_u10(n:n))
         call da_interp_2d_partial (grid%xb%v10,  iv%instid(inst)%info, 1, n, n, model_v10(n:n))
         call da_interp_2d_partial (grid%xb%psfc, iv%instid(inst)%info, 1, n, n, model_psfc(n:n))

         model_psfc(n) = 0.01*model_psfc(n)           ! convert to hPa
         model_ptop    = 0.01*grid%xb%ptop

         ! get 1 levels (0.005hPa at top) /model full level
         Atmosphere(1)%Level_Pressure(0)                      = model_ptop    ! to sigma level 51-->sigmaf=0
         Atmosphere(1)%Level_Pressure(Atmosphere(1)%n_Layers) = model_psfc(n) ! to sigma level 1->sigmaf=1

         ! calculating full-level P from half-level P using
         ! linear interpolation: (pf-pf1)/(sigmaf-sigmaf1)=(ph2-ph1)/(sigmah2-sigmah1)
         ! note: sigma is from bottem to top; pressure is from top to bottem
         do k=kts,kte-1 ! from top to bottem
            Atmosphere(1)%Level_Pressure(k)= (Atmosphere(1)%Pressure(k+1) - Atmosphere(1)%Pressure(k))/ &
                                             (grid%xb%sigmah(kte-k) - grid%xb%sigmah(kte-k+1)) * &
                                             (grid%xb%sigmaf(kte-k+1) - grid%xb%sigmah(kte-k+1)) + &
                                              Atmosphere(1)%Pressure(k)
         end do
 
         ! check for pressure monotonicity
         do k = kte, kts+1, -1  ! from bottom to top
            if ( Atmosphere(1)%Level_Pressure(k) <= Atmosphere(1)%Level_Pressure(k-1) ) then
               write(message(1),'(a,2i5.0,a)') ' Skipping the pixel at loc ', i, j,  &
                  ' where the pressure is not monotonic'
               call da_warning("da_get_innov_vector_crtm.inc",433,message(1:1))
               iv%instid(inst)%tb_inv(:,n) = missing_r
               cycle pixel_loop
            end if
         end do

         ! convert cloud content unit from kg/kg to kg/m^2        
         if (crtm_cloud) then
            do k=kts,kte
               do icld=1,Atmosphere(1)%n_Clouds
                  Atmosphere(1)%Cloud(icld)%Water_Content(k)= Atmosphere(1)%Cloud(icld)%Water_Content(k)* &
                     (Atmosphere(1)%Level_Pressure(k)- Atmosphere(1)%Level_Pressure(k-1))*100.0/gravity 
               end do
            end do
         end if

         if ( model_isflg == 0 ) then   ! over sea using SST
            call da_interp_2d_partial (grid%xb % tgrn, iv%instid(inst)%info, 1, n, n, model_ts(n:n))
         else
            call da_interp_2d_partial (grid%xb % tsk,  iv%instid(inst)%info, 1, n, n, model_ts(n:n))
         end if

         call da_interp_2d_partial (grid%xb % terr, iv%instid(inst)%info, 1, n, n, model_elv(n:n))

         ! variables for emissivity calculations
         !---------------------------------------- 
         call da_interp_2d_partial (grid%xb % smois,  iv%instid(inst)%info, 1, n, n, model_smois(n:n) )
         call da_interp_2d_partial (grid%xb % tslb,   iv%instid(inst)%info, 1, n, n, model_tslb(n:n) )
         call da_interp_2d_partial (grid%xb % snowh,  iv%instid(inst)%info, 1, n, n, model_snowh(n:n) )
         call da_interp_2d_partial (grid%xb % vegfra, iv%instid(inst)%info, 1, n, n, model_vegfra(n:n) )

         ! model_snowh(n) = model_snowh(n)*100.0   ! convert from m to mm
         model_vegfra(n) = 0.01*model_vegfra(n)  ! convert range to 0~1

         ! ADD for computing cloud liquid water path (mm) from guess
         clwp = 0.0
         do k = kts,kte ! from top to bottom
            dpf(k) = 100.0*(Atmosphere(1)%level_pressure(k) - Atmosphere(1)%level_pressure(k-1))
            clw  (k) = model_qcw(k)*dpf(k)/gravity ! kg/m2 or mm
            if (Atmosphere(1)%pressure(k)<100.0) clw(k) = 0.0
            clwp  = clwp + clw(k)
         end do

         ! 1 GeometryInfo Structure
         GeometryInfo(1)%Sensor_Zenith_Angle=iv%instid(inst)%satzen(n)
         GeometryInfo(1)%Source_Zenith_Angle=iv%instid(inst)%solzen(n)
         GeometryInfo(1)%iFOV=iv%instid(inst)%scanpos(n)
         ! GeometryInfo(1)%Satellite_Height=830.0
         ! GeometryInfo(1)%Sensor_Scan_Angle=
         ! GeometryInfo(1)%Sensor_Zenith_Angle=
         ! GeometryInfo(1)%Sensor_Scan_Angle=
         ! GeometryInfo(1)%Source_Zenith_Angle=

         ! 1 Surface parameter data

         if ( use_mspps_ts(inst) ) then
            ! only for AMSU-A over land
            if ( trim(crtm_sensor_name(rtminit_sensor(inst))) == 'amsua' &
                 .and. model_isflg == 2 ) then
               call da_mspps_ts(ob%instid(inst)%tb(1:nchanl,n), nchanl,  &
                                iv%instid(inst)%satzen(n), ts_mspps)
               ! ts_mspps is initilaized as negative values in the
               ! da_mspps_ts subroutine.  Apply only valid values here.
               if ( ts_mspps > 0.0 ) then
                  model_ts(n) = ts_mspps
               end if
            end if
         end if

         if (Surface(1)%Land_Coverage > 0.0) then
            Surface(1)%Land_Type       = nint(model_ivgtyp)
            Surface(1)%Vegetation_Type = max(1,wrf_to_crtm_mw(nint(model_ivgtyp)))
            Surface(1)%Soil_Type       = max(1,wrf_to_crtm_soil(nint(model_isltyp)))
            ! glacial land ice soil type and glacial vegetation type are not valid land types
            if ( surface(1)%Soil_Type == 9 .or. surface(1)%Vegetation_Type == 13 ) then
               surface(1)%ice_coverage = min(surface(1)%ice_coverage + surface(1)%land_coverage, 1.0)
               surface(1)%land_coverage = 0
            end if
            ! reset water-dominated land type
            if ( Surface(1)%Land_Type == grid%iswater ) Surface(1)%Land_Type = 1
            Surface(1)%Land_Temperature=model_ts(n)      ! K
            Surface(1)%Soil_Moisture_Content= model_smois(n) !0.05    ! volumetric water content (g/cm**3)
            ! Surface(1)%Canopy_Water_Content=0.05      ! gravimetric water content
            Surface(1)%Vegetation_Fraction=model_vegfra(n)
            Surface(1)%Soil_Temperature=model_tslb(n)
         end if
         if (Surface(1)%Water_Coverage > 0.0) then
            ! Surface%Water_Type=SEA_WATER          ! (Currently NOT used)
            Surface(1)%Water_Temperature=model_ts(n)     ! K
            Surface(1)%Wind_Speed=sqrt(model_u10(n)**2+model_v10(n)**2)  ! m/sec
            ! surface(1)%Wind_Direction=0.0            ! NOT used
            Surface(1)%Salinity=33.0                   ! ppmv
         end if

         if (Surface(1)%Snow_Coverage > 0.0) then
            Surface(1)%Snow_Type=2 !NEW_SNOW
            Surface(1)%Snow_Temperature=model_ts(n)      ! K
            Surface(1)%Snow_Depth=model_snowh(n)         ! mm
            ! Surface(1)%Snow_Density=0.2               ! g/cm**3
            ! Surface(1)%Snow_Grain_Size=2.0            ! mm
         end if
         if (Surface(1)%Ice_Coverage > 0.0) then
            ! Surface(1)%Ice_Type=FRESH_ICE             ! NO Table offered, single example is FRESH_ICE
            Surface(1)%Ice_Temperature=model_ts(n)       ! K
            Surface(1)%Ice_Thickness=10.0              ! mm
            ! Surface(1)%Ice_Density=0.9                ! g/cm**3
            ! Surface(1)%Ice_Roughness=0.0               ! NO Table offered, single example is ZERO

         end if
         if (nchanl > 0) then
            Surface(1)%SensorData%n_channels        = nchanl
            Surface(1)%SensorData%Sensor_Channel    = ChannelInfo(inst)%Sensor_Channel
            Surface(1)%SensorData%Sensor_Id         = ChannelInfo(inst)%Sensor_Id
            Surface(1)%SensorData%WMO_Satellite_Id  = ChannelInfo(inst)%WMO_Satellite_Id
            Surface(1)%SensorData%WMO_Sensor_Id     = ChannelInfo(inst)%WMO_Sensor_Id
            Surface(1)%SensorData%Tb(1:nchanl)      = ob%instid(inst)%tb(1:nchanl,n)
         end if

         ! user emissivity options
         Options(1)%use_emissivity = .false.
         if ( use_mspps_emis(inst) ) then
            ! Only for AMSU-A over land
            if ( trim(crtm_sensor_name(rtminit_sensor(inst))) == 'amsua' &
                 .and. model_isflg == 2 ) then
               call da_mspps_emis(ob%instid(inst)%tb(1:nchanl,n), nchanl, em_mspps)
               Options(1)%use_emissivity = .true.
               Options(1)%emissivity(1:nchanl) = em_mspps(1:nchanl)
            end if
         end if

         if (use_crtm_kmatrix) then

            ! 1 surface/atmosphere K initialization
            do l = 1, ChannelInfo(inst)%n_Channels
               ! -- Copy the adjoint atmosphere structure
               Atmosphere_K(l,1) = Atmosphere(1)
 
               ! -- Copy the adjoint surface structure
               Surface_K(l,1) = Surface(1)
            end do

            ! -- Zero the Adjoint outputs
            ! Important: adjoint variables must be initialized
            call CRTM_Atmosphere_Zero( Atmosphere_K )
            call CRTM_Surface_Zero( Surface_K )

            ! Assign tb = R^-1 Re :
            RTSolution_K(:,1)%brightness_temperature = 1.
            RTSolution_K(:,1)%radiance = 0.
	    
            ! [1.3] Call RTM K-Matrix model
            call da_crtm_k(1, nchanl, 1, Atmosphere,   &
                               Surface,      &
                               RTSolution_K,&
                               GeometryInfo, &
                               ChannelInfo(inst),  &
                               Atmosphere_K,&
                               Surface_K,   &
                               RTSolution,  &
                               Options)
         else
	 
            ! [1.3] Call RTM forward model
            call da_crtm_direct(1, nchanl, 1, Atmosphere,   &
               Surface,      &
               GeometryInfo, &
               ChannelInfo(inst:inst),  &
               RTSolution,              &
               Options)
	    
	 end if

         ! Compute Overcast Radiances for AIRS Cloud Detection (MMR)
         !----------------------------------------------------------------
	
   !---------------------------------------------------------------------------
   ! Precondition satellite control variable (Cloud Cover(s)):
   !---------------------------------------------------------------------------
         if (.false.) then  
!         if (use_satcv(2)) then  
            nclouds = iv%instid(inst)%cv_index(n)%nclouds
	    allocate (hessian(nclouds,nclouds))
	    allocate (eignvec(nclouds,nclouds))
	    allocate (eignval(nclouds))
            hessian(:,:)= 0.0      
            do ilev=1, nclouds
               do jlev=ilev, nclouds
	          do k = 1, nchanl
                     if (ALL(iv%instid(inst)%ichan(k) /=  Bands(:,1))) cycle   ! Only Channels in Band 1 	       
                        rad_clr      = iv%instid(inst)%rad_xb(k,n)             
                        rad_ovc_ilev = iv%instid(inst)%rad_ovc(k,kte-nclouds+ilev,n)
                        rad_ovc_jlev = iv%instid(inst)%rad_ovc(k,kte-nclouds+jlev,n)
                        hessian(ilev,jlev) = hessian(ilev,jlev) + (rad_ovc_ilev-rad_clr) * (rad_ovc_jlev-rad_clr)
      	          end do			     	 
                 hessian(jlev,ilev) = hessian(ilev,jlev)   
               end do
            end do  
	      
            call da_eof_decomposition(nclouds, hessian, eignvec, eignval)
	      
!	    if (ANY(eignval <= 0)) write(unit=stdout,fmt='(3A,I4,A,100F18.5)')      &
!	       'SATCV: non-positive Hessian for ', trim(iv%instid(inst)%rttovid_string), ' ,pixel ',n,'--> Eigenvalues =',eignval 
 	      
	    do ilev = 1, nclouds
	       do jlev = ilev, nclouds
                  iv%instid(inst)%cv_index(n)%vtox(ilev,jlev) = &
		     sum( eignvec(ilev,:) * sqrt(1.0/eignval(:)) * eignvec(jlev,:), mask = eignval >0 )
		       
		  iv%instid(inst)%cv_index(n)%vtox(jlev,ilev) = iv%instid(inst)%cv_index(n)%vtox(ilev,jlev) 
               end do
	    end do
	      	 
            deallocate(hessian,eignvec,eignval)
         end if   

         !----------------------------------------------------------------
         ! [2.0] calculate components of innovation vector:
         !----------------------------------------------------------------
         do k = 1, nchanl
            iv%instid(inst)%tb_xb(k,n)  = RTSolution(k,1)%Brightness_Temperature
            iv%instid(inst)%rad_xb(k,n) = RTSolution(k,1)%Radiance
            iv%instid(inst)%emiss(k,n)  = RTSolution(k,1)%surface_emissivity
      
            if (use_pseudo_rad .or. use_simulated_rad) then ! input is innovation
              ob%instid(inst)%tb(k,n) = missing_r
	      if ( iv%instid(inst)%tb_inv(k,n) > missing_r ) &
                ob%instid(inst)%tb(k,n) = RTSolution(k,1)%Brightness_Temperature + iv%instid(inst)%tb_inv(k,n)
            else
               if ( iv%instid(inst)%tb_inv(k,n) > missing_r ) & 
                  iv%instid(inst)%tb_inv(k,n) = ob%instid(inst)%tb(k,n) - iv%instid(inst)%tb_xb(k,n)	    
	    end if

            if (use_crtm_kmatrix) then
               ! surface Jacobian
               iv%instid(inst)%ts_jacobian(k,n) = Surface_k(k,1)%water_temperature
               iv%instid(inst)%windspeed_jacobian(k,n) = Surface_k(k,1)%wind_speed
               iv%instid(inst)%emiss_jacobian(k,n) = RTSolution_k(k,1)%surface_emissivity
	    end if   
         end do

         !----------------------------------------------------------------
         ! [3.0] store base state (and Jacobian) to innovation structure
         !----------------------------------------------------------------
         ! full level pressures
         iv%instid(inst)%pf(0,n)  = Atmosphere(1)%level_pressure(0)
         if (use_crtm_kmatrix) then
            ! PS Jacobian
            do l=1,nchanl
               iv%instid(inst)%ps_jacobian(l,n)  = Atmosphere_k(l,1)%level_pressure(Atmosphere(1)%n_layers)
            end do
	 end if
         do k=1,Atmosphere(1)%n_layers
            iv%instid(inst)%pm(k,n)  = Atmosphere(1)%pressure(k)
            iv%instid(inst)%pf(k,n)  = Atmosphere(1)%level_pressure(k)
            iv%instid(inst)%tm(k,n)  = Atmosphere(1)%temperature(k)
            iv%instid(inst)%qm(k,n)  = Atmosphere(1)%absorber(k,1)

            if (use_crtm_kmatrix) then
               ! T, Q Jacobian
               do l=1,nchanl
                  iv%instid(inst)%t_jacobian(l,k,n) = Atmosphere_k(l,1)%temperature(k)
                  iv%instid(inst)%q_jacobian(l,k,n) = Atmosphere_k(l,1)%absorber(k,1)
               end do
	    end if
         end do
	       
         if (crtm_cloud) then
            do k=1,Atmosphere(1)%n_layers
               iv%instid(inst)%qcw(k,n) = Atmosphere(1)%cloud(1)%water_content(k)
               iv%instid(inst)%qci(k,n) = Atmosphere(1)%cloud(2)%water_content(k)
               iv%instid(inst)%qrn(k,n) = Atmosphere(1)%cloud(3)%water_content(k)
               iv%instid(inst)%qsn(k,n) = Atmosphere(1)%cloud(4)%water_content(k)
               iv%instid(inst)%qgr(k,n) = Atmosphere(1)%cloud(5)%water_content(k)
               iv%instid(inst)%qhl(k,n) = Atmosphere(1)%cloud(6)%water_content(k)
               iv%instid(inst)%rcw(k,n) = Atmosphere(1)%cloud(1)%effective_radius(k)
               iv%instid(inst)%rci(k,n) = Atmosphere(1)%cloud(2)%effective_radius(k)
               iv%instid(inst)%rrn(k,n) = Atmosphere(1)%cloud(3)%effective_radius(k)
               iv%instid(inst)%rsn(k,n) = Atmosphere(1)%cloud(4)%effective_radius(k)
               iv%instid(inst)%rgr(k,n) = Atmosphere(1)%cloud(5)%effective_radius(k)
               iv%instid(inst)%rhl(k,n) = Atmosphere(1)%cloud(6)%effective_radius(k)
            
               if (use_crtm_kmatrix) then
	          ! Cloud Jacobian
                  do l=1,nchanl
                     iv%instid(inst)%water_jacobian(l,k,n)    = Atmosphere_k(l,1)%cloud(1)%water_content(k)
                     iv%instid(inst)%ice_jacobian(l,k,n)      = Atmosphere_k(l,1)%cloud(2)%water_content(k)
                     iv%instid(inst)%rain_jacobian(l,k,n)     = Atmosphere_k(l,1)%cloud(3)%water_content(k)
                     iv%instid(inst)%snow_jacobian(l,k,n)     = Atmosphere_k(l,1)%cloud(4)%water_content(k)
                     iv%instid(inst)%graupel_jacobian(l,k,n)  = Atmosphere_k(l,1)%cloud(5)%water_content(k)
                     iv%instid(inst)%hail_jacobian(l,k,n)     = Atmosphere_k(l,1)%cloud(6)%water_content(k)

                     iv%instid(inst)%water_r_jacobian(l,k,n)  = Atmosphere_k(l,1)%cloud(1)%effective_radius(k)
                     iv%instid(inst)%ice_r_jacobian(l,k,n)    = Atmosphere_k(l,1)%cloud(2)%effective_radius(k)
                     iv%instid(inst)%rain_r_jacobian(l,k,n)   = Atmosphere_k(l,1)%cloud(3)%effective_radius(k)
                     iv%instid(inst)%snow_r_jacobian(l,k,n)   = Atmosphere_k(l,1)%cloud(4)%effective_radius(k)
                     iv%instid(inst)%graupel_r_jacobian(l,k,n)= Atmosphere_k(l,1)%cloud(5)%effective_radius(k)
                     iv%instid(inst)%hail_r_jacobian(l,k,n)   = Atmosphere_k(l,1)%cloud(6)%effective_radius(k)
                  end do
	       end if	  
            end do
         end if

         if ( calc_weightfunc .and. use_crtm_kmatrix ) then
            do l = 1, nchanl
               total_od = 0.0
               do k = 1, Atmosphere(1)%n_layers
                  iv%instid(inst)%lod(l,k,n) = RTSolution(l,1)%layer_optical_depth(k)
                  iv%instid(inst)%lod_jacobian(l,k,n) = RTSolution_K(l,1)%layer_optical_depth(k)
                  total_od = total_od + RTSolution(l,1)%layer_optical_depth(k)
                  iv%instid(inst)%trans(l,k,n) = exp(-total_od)
                  iv%instid(inst)%trans_jacobian(l,k,n) = -iv%instid(inst)%lod_jacobian(l,k,n)*exp(iv%instid(inst)%lod(l,k,n))
               end do
               iv%instid(inst)%der_trans(l,1,n) = 0.0
               do k= 2, Atmosphere(1)%n_layers
                  iv%instid(inst)%der_trans(l,k,n) =  abs( (iv%instid(inst)%trans(l,k-1,n)-iv%instid(inst)%trans(l,k,n))/ &
                                                           (log(iv%instid(inst)%pm(k-1,n))-log(iv%instid(inst)%pm(k,n))) )
               end do
            end do
         end if
         !----------------------------------------------
         ! [4.0] store surface information to innovation structure
         !----------------------------------------------
         iv%instid(inst)%u10(n)       = model_u10(n)
         iv%instid(inst)%v10(n)       = model_v10(n)
         iv%instid(inst)%t2m(n)       = 0.01*missing_r !model_t2m
         iv%instid(inst)%mr2m(n)      = 0.01*missing_r !model_mr2m
         iv%instid(inst)%ps(n)        = model_psfc(n)
         iv%instid(inst)%ts(n)        = model_ts(n)
         iv%instid(inst)%smois(n)     = model_smois(n)
         iv%instid(inst)%tslb(n)      = model_tslb(n)
         iv%instid(inst)%snowh(n)     = model_snowh(n)
         iv%instid(inst)%isflg(n)     = model_isflg
         iv%instid(inst)%elevation(n) = model_elv(n)
         iv%instid(inst)%soiltyp(n)   = model_isltyp
         iv%instid(inst)%vegtyp(n)    = model_ivgtyp
         iv%instid(inst)%vegfra(n)    = model_vegfra(n)
         iv%instid(inst)%clwp(n)      = clwp
         iv%instid(inst)%water_coverage(n) = Surface(1)%water_coverage
         iv%instid(inst)%land_coverage(n)  = Surface(1)%land_coverage
         iv%instid(inst)%ice_coverage(n)   = Surface(1)%ice_coverage                              
         iv%instid(inst)%snow_coverage(n)  = Surface(1)%snow_coverage
      end do pixel_loop      !  end loop for pixels

      if ( use_mspps_emis(inst) ) deallocate(em_mspps)

      deallocate (model_u10)
      deallocate (model_v10)
      deallocate (model_psfc)
      deallocate (model_ts)
      deallocate (model_tslb)
      deallocate (model_snowh)
      deallocate (model_snow)
      deallocate (model_elv)
      deallocate (model_vegfra)
      deallocate (model_smois)

      call CRTM_RTSolution_Destroy(RTSolution)
      deallocate( RTSolution, STAT = Allocate_Status )
      if (Allocate_Status /= 0) &
         call da_error("da_get_innov_vector_crtm.inc",841, &
            (/"Error in deallocating RTSolution"/))

      if (use_crtm_kmatrix) then
         call CRTM_RTSolution_Destroy(RTSolution_K)
         deallocate( RTSolution_K, STAT = Allocate_Status )
         if (Allocate_Status /= 0) &
            call da_error("da_get_innov_vector_crtm.inc",848, &
               (/"Error in deallocating RTSolution_K"/))
      end if	 

      call CRTM_Options_Destroy(Options)
      IF ( CRTM_Options_Associated( Options(1) ) ) &
         call da_error("da_get_innov_vector_crtm.inc",854, &
            (/"Error in deallocating CRTM Options Structure"/))

      call CRTM_Surface_Destroy(Surface)
      IF ( CRTM_Surface_Associated( Surface(1) ) ) &
         call da_error("da_get_innov_vector_crtm.inc",859, &
            (/"Error in deallocating CRTM Surface Structure"/))

      if (use_crtm_kmatrix) then
         call CRTM_Surface_Destroy(Surface_K)
         IF ( ANY( CRTM_Surface_Associated( Surface_K(1,:)) ) ) &
            call da_error("da_get_innov_vector_crtm.inc",865, &
               (/"Error in deallocatting CRTM Surface_K Structure"/))
      endif

      if (use_crtm_kmatrix) then
         call CRTM_Atmosphere_Destroy( Atmosphere_K )
         IF ( ANY( CRTM_Atmosphere_Associated( Atmosphere_K(1,:)) ) ) &
            call da_error("da_get_innov_vector_crtm.inc",872, &
               (/"Error in deallocatting CRTM Atmosphere_K Structure"/))
      endif

      if (use_crtm_kmatrix) then
         deallocate( Atmosphere_K, Surface_K, STAT = Allocate_Status )
         if ( Allocate_Status /= 0 ) &
            call da_error("da_get_innov_vector_crtm.inc",879, &
               (/"Error in deallocatting CRTM Surface_K Structure"/))
      endif

   end do        ! end loop for sensor

   deallocate (wrf_to_crtm_mw)

   call CRTM_Atmosphere_Destroy (Atmosphere)
   IF ( CRTM_Atmosphere_Associated( Atmosphere(1) ) ) &
       call da_error("da_get_innov_vector_crtm.inc",889, &
         (/"Error in deallocating CRTM Atmosphere Structure"/))

   if (trace_use) call da_trace_exit("da_get_innov_vector_crtm")
 
end subroutine da_get_innov_vector_crtm

subroutine da_crtm_tl(  nsensor, nchan, nprof, Atmosphere,   &
                            Surface,      &
                            Atmosphere_TL,&
                            Surface_TL,   &
                            GeometryInfo, &
                            ChannelInfo,  &
                            RTSolution,   &
                            RTSolution_TL,&
                            Options)

   integer, intent(in)            :: nsensor, nchan, nprof
   type (CRTM_RTSolution_type ),  intent(inout)  :: RTSolution(nchan,nprof),RTSolution_TL(nchan,nprof)
   type (CRTM_ChannelInfo_type),  intent(in)  :: ChannelInfo(nsensor)
   type( CRTM_Atmosphere_type ),  intent(in)  :: Atmosphere(nprof), Atmosphere_TL(nprof)
   type( CRTM_Surface_type ),     intent(in)  :: Surface(nprof),Surface_TL(nprof)
   type( CRTM_Geometry_type ),    intent(inout)  :: GeometryInfo(nprof)
   type (CRTM_Options_type),      intent(in)     :: Options(nprof)

   integer :: Error_Status

   if (trace_use) call da_trace_entry("da_crtm_tl")

         Error_Status = CRTM_Tangent_Linear(Atmosphere,   &
                            Surface,      &
                            Atmosphere_TL,&
                            Surface_TL,   &
                            GeometryInfo, &
                            ChannelInfo,  &
                            RTSolution,   &
                            RTSolution_TL,&
                            Options) 

         if ( Error_Status /= 0 ) then
              call da_error("da_crtm_tl.inc",35, &
                 (/"Error in calling CRTM_Tangent_Linear"/))
         end if


   if (trace_use) call da_trace_exit("da_crtm_tl")

end subroutine da_crtm_tl
subroutine da_crtm_k(  nsensor, nchan, nprof, &
                            Atmosphere,   &
                            Surface,      &
                            RTSolution_K, &
                            GeometryInfo, &
                            ChannelInfo,  &
                            Atmosphere_K,   &
                            Surface_K,   &
                            RTSolution,   &
                            Options)

   integer, intent(in)            :: nsensor, nchan, nprof
   type (CRTM_RTSolution_type ),  intent(inout)  :: RTSolution(nchan,nprof)
   type (CRTM_RTSolution_type ),  intent(inout)  :: RTSolution_K(nchan,nprof)
   type (CRTM_ChannelInfo_type),  intent(in)     :: ChannelInfo(nsensor)
   type( CRTM_Atmosphere_type ),  intent(in)     :: Atmosphere(nprof)
   type( CRTM_Atmosphere_type ),  intent(inout)  :: Atmosphere_K(nchan,nprof)
   type( CRTM_Surface_type ),     intent(in)     :: Surface(nprof)
   type( CRTM_Surface_type ),     intent(inout)  :: Surface_K(nchan,nprof)
   type( CRTM_Geometry_type ),    intent(inout)  :: GeometryInfo(nprof)
   type (CRTM_Options_type),      intent(in)     :: Options(nprof)

   integer :: Error_Status

   if (trace_use) call da_trace_entry("da_crtm_k")

         Error_Status = CRTM_K_Matrix(Atmosphere,   &
                            Surface,      &
                            RTSolution_K,&
                            GeometryInfo, &
                            ChannelInfo,  &
                            Atmosphere_K,&
                            Surface_K,   &
                            RTSolution,   &
                            Options)
         if ( Error_Status /= 0 ) then
              call da_error("da_crtm_k.inc",38, &
                 (/"Error in calling CRTM_K_Matrix"/))
         end if

   if (trace_use) call da_trace_exit("da_crtm_k")

end subroutine da_crtm_k
subroutine da_crtm_direct(  nsensor, nchan, nprof, Atmosphere,   &
                            Surface,      &
                            GeometryInfo, &
                            ChannelInfo,  &
                            RTSolution,   &
                            Options)

   integer, intent(in)            :: nsensor, nchan, nprof
   type (CRTM_RTSolution_type ),  intent(inout)  :: RTSolution(nchan,nprof)
   type (CRTM_ChannelInfo_type),  intent(in)  :: ChannelInfo(nsensor)
   type( CRTM_Atmosphere_type ),  intent(in)  :: Atmosphere(nprof)
   type( CRTM_Surface_type ),     intent(in)  :: Surface(nprof)
   type( CRTM_Geometry_type ), intent(inout)  :: GeometryInfo(nprof)
   type (CRTM_Options_type), optional, intent(in) :: Options(nprof)

   integer :: Error_Status

   if (trace_use) call da_trace_entry("da_crtm_direct")

         Error_Status = CRTM_Forward (Atmosphere,   &
                            Surface,      &
                            GeometryInfo, &
                            ChannelInfo,  &
                            RTSolution,   &
                            Options)

         if ( Error_Status /= 0 ) then
              call da_error("da_crtm_direct.inc",29, &
                 (/"Error in calling CRTM_Forward"/))
         end if

  if (trace_use) call da_trace_exit("da_crtm_direct")

end subroutine da_crtm_direct
subroutine da_crtm_ad( nsensor, nchan, nprof, Atmosphere,   &
                            Surface,      &
                            RTSolution_AD, &
                            GeometryInfo, &
                            ChannelInfo,  &
                            Atmosphere_AD,   &
                            Surface_AD,   &
                            RTSolution,   &
                            Options)

   integer, intent(in)            :: nsensor, nchan, nprof
   type (CRTM_RTSolution_type ),  intent(inout)  :: RTSolution(nchan,nprof)
   type (CRTM_RTSolution_type ),  intent(inout)  :: RTSolution_AD(nchan,nprof)
   type (CRTM_ChannelInfo_type),  intent(in)  :: ChannelInfo(nsensor)
   type( CRTM_Atmosphere_type ),  intent(in)  :: Atmosphere(nprof)
   type( CRTM_Atmosphere_type ),  intent(inout)  :: Atmosphere_AD(nprof)
   type( CRTM_Surface_type ),     intent(in)  :: Surface(nprof)
   type( CRTM_Surface_type ),     intent(inout)  :: Surface_AD(nprof)
   type( CRTM_Geometry_type ),    intent(inout)  :: GeometryInfo(nprof)
   type (CRTM_Options_type),      intent(in)     :: Options(nprof)

   integer :: Error_Status

   if (trace_use) call da_trace_entry("da_crtm_ad")

         Error_Status = CRTM_Adjoint(Atmosphere,   &
                            Surface,      &
                            RTSolution_AD,&
                            GeometryInfo, &
                            ChannelInfo,  &
                            Atmosphere_AD,&
                            Surface_AD,   &
                            RTSolution,   &
                            Options)
         if ( Error_Status /= 0 ) then
              call da_error("da_crtm_ad.inc",37, &
                 (/"Error in calling CRTM_Adjoint"/))
         end if

   if (trace_use) call da_trace_exit("da_crtm_ad")

end subroutine da_crtm_ad
subroutine da_crtm_init(iv,ob, nsensor)
!------------------------------------------------------------------------------
!  PURPOSE: interface to the initialization subroutine of 1
!
!  METHOD:  read 1 coefs files
!
!  HISTORY: 10/15/2006  added crtm initialization    Tomislava Vukicevic, ATOC, University of Colorado
!           11/09/2006  Updated                      Zhiquan Liu
!           10/24/2007  limit to 1 init           Tom Auligne
!------------------------------------------------------------------------------

 implicit none 

 type (iv_type), intent (inout) :: iv
 type (y_type) , intent (inout) :: ob
 integer ,       intent (in)    :: nsensor

!
!  local arguments
!------------------- 
 integer   :: n, j, ichan

!
! 1 local ---------------------------------------------------
!
  integer :: Error_Status
!  character( 256 ) :: SpcCoeff_File
!  character( 256 ) :: TauCoeff_File
  character( 256 ) :: AerosolCoeff_File
  character( 256 ) :: CloudCoeff_File
  character( 256 ) :: File_Path
!  character( 80 ), pointer :: Sensor_Descriptor(:)
!
! end of 1 local

  call da_trace_entry("da_crtm_init")

!---------------------------------------------------------------------
! 1.0 get 1 sensor descriptor
!---------------------------------------------------------------------
  allocate(Sensor_Descriptor(nsensor))
  call da_crtm_sensor_descriptor(nsensor,Sensor_Descriptor)
  allocate(ChannelInfo(nsensor))

! 1 load coefficients
!-----------------------------------------------------------
! 1.1 call CRTM_Init to load coefficients and fill ChannelInfo structure
!-----------------------------------------------------------
  ! input: 
     AerosolCoeff_File = 'AerosolCoeff.bin'
     CloudCoeff_File   = 'CloudCoeff.bin'
     File_Path         = trim(crtm_coef_path)//'/'
  !----------------------------------------------------------------
  ! ChannelInfo structure contains on output: 
  !
  ! n_channels - integer, total number of channels
  ! Sensor_Index - integer
  ! Channel_Index - integer pointer, index of the channels loaded during initialization
  ! Sensor_Channel - integer pointer, the sensor channel #
  ! Sensor_ID - character pointer, character string containing satellite and sensor descr
  !                                        example: amsre_aqua (Appendix B in User guide)
  ! WMO_Satellite_ID - integer pointer
  ! WMO_Sensor_ID - integer pointer
  !----------------------------------------------------------------- 

     Error_Status = CRTM_Init(Sensor_Descriptor, &
                              ChannelInfo, &
                              AerosolCoeff_File = AerosolCoeff_File, &
                              CloudCoeff_File = CloudCoeff_File, &
                              IRwaterCoeff_File = crtm_irwater_coef, &
                              MWwaterCoeff_File = crtm_mwwater_coef, &
                              IRlandCoeff_File  = crtm_irland_coef,  &
                              VISlandCoeff_File = crtm_visland_coef, &
                              File_Path = File_Path) 

     if ( Error_Status /= 0 ) then 
       call da_error("da_crtm_init.inc",77, &
         (/"Error in initializing CRTM"/))
     END IF

     iv%instid(1:nsensor)%nlevels = kme-kms+1

     if (print_detail_rad) then
        do n = 1, nsensor
           write (message(1),*) 'in da_crtm_init: ChannelInfo content'
           write (message(2),*) 'Sensor_Index ',ChannelInfo(n)%Sensor_Index
           write (message(3),*) 'n_channels ',ChannelInfo(n)%n_channels
           write (message(4),*) 'Channel_Index ',ChannelInfo(n)%Channel_Index(:)
           write (message(5),*) 'Sensor_Channel ',ChannelInfo(n)%Sensor_Channel(:)
           write (message(6),*) 'Sensor_ID ',ChannelInfo(n)%Sensor_ID
           write (message(7),*) 'WMO_Satellite_ID ',ChannelInfo(n)%WMO_Satellite_ID
           write (message(8),*) 'WMO_Sensor_ID ',ChannelInfo(n)%WMO_Sensor_ID
           call da_message(message(1:8))
       end do
    end if

  call da_trace_exit("da_crtm_init")

end subroutine da_crtm_init
subroutine da_crtm_sensor_descriptor(nsensor,sensor_descriptor)

 integer,        intent(in)  :: nsensor
 character(len=80), intent(inout) :: sensor_descriptor(nsensor)

 integer :: i,platform_id,satellite_id,sensor_id
 character (len=80) :: crtm_sat, crtm_sensor

   if (trace_use) call da_trace_entry("da_crtm_sensor_descriptor")  

  do i=1,nsensor
     
     platform_id  = rtminit_platform(i)
     satellite_id = rtminit_satid(i)
     sensor_id    = rtminit_sensor(i)

     if (trim(crtm_platform_name(platform_id)) == 'eos') then
        if (satellite_id == 2) crtm_sat='aqua'
        if (satellite_id == 1) crtm_sat='terra'
     else if (trim(crtm_platform_name(platform_id)) == 'metop') then
        if (satellite_id == 1) crtm_sat='metop-b'
        if (satellite_id == 2) crtm_sat='metop-a'
        if (satellite_id == 3) crtm_sat='metop-c'
     else if (trim(crtm_platform_name(platform_id)) == 'tiros') then
        if (satellite_id == 0) crtm_sat='tirosn'
     else if (trim(crtm_platform_name(platform_id)) == 'fy3') then
        if (satellite_id == 1) crtm_sat='fy3a'
        if (satellite_id == 2) crtm_sat='fy3b'
     else if (trim(crtm_platform_name(platform_id)) == 'npp') then
        if (satellite_id == 0) crtm_sat='npp'
     else if (trim(crtm_platform_name(platform_id)) == 'msg') then
        if (satellite_id == 1) crtm_sat='m08'
        if (satellite_id == 2) crtm_sat='m09'
        if (satellite_id == 3) crtm_sat='m10'
     else if (trim(crtm_platform_name(platform_id)) == 'gcom-w') then
        if (satellite_id == 1) crtm_sat='gcom-w1'
     else
        write(crtm_sat, '(a,i2.2)')  &
             trim( crtm_platform_name(platform_id) ),satellite_id
     end if

     if ( trim(crtm_sensor_name(sensor_id)) == 'airs' ) then
        crtm_sensor='airs281'
     elseif ( trim(crtm_sensor_name(sensor_id)) == 'iasi' ) then
        crtm_sensor='iasi616'		
     elseif ( trim(crtm_sensor_name(sensor_id)) == 'hirs' ) then
        if (satellite_id <= 14) crtm_sensor='hirs2'
        if (satellite_id >= 15 .and. satellite_id <= 17) crtm_sensor='hirs3'
        if (satellite_id == 18) crtm_sensor='hirs4'
     elseif ( trim(crtm_sensor_name(sensor_id)) == 'avhrr' ) then
        if (satellite_id <= 14) crtm_sensor='avhrr2'
        if (satellite_id >= 15 .and. satellite_id <= 17) crtm_sensor='avhrr3'
        if (satellite_id == 18) crtm_sensor='avhrr4'
     else
        crtm_sensor=crtm_sensor_name(sensor_id)
     end if

     sensor_descriptor(i)=trim(crtm_sensor)//'_'//trim(crtm_sat)

  end do

   if (trace_use) call da_trace_exit("da_crtm_sensor_descriptor")  

end subroutine da_crtm_sensor_descriptor
subroutine da_det_crtm_climat (lat, month, crtm_climat)

! determine the 1 climatology model according to the 
! latitude and month
! 1 climatology model is categorized in CRTM_Atmosphere_Define.f90 as
! INTEGER, PARAMETER :: INVALID_MODEL          = 0
! INTEGER, PARAMETER :: TROPICAL               = 1
! INTEGER, PARAMETER :: MIDLATITUDE_SUMMER     = 2
! INTEGER, PARAMETER :: MIDLATITUDE_WINTER     = 3
! INTEGER, PARAMETER :: SUBARCTIC_SUMMER       = 4
! INTEGER, PARAMETER :: SUBARCTIC_WINTER       = 5
! INTEGER, PARAMETER :: US_STANDARD_ATMOSPHERE = 6

   implicit none

   real,    intent(in)  :: lat
   integer, intent(in)  :: month
   integer, intent(out) :: crtm_climat

   crtm_climat = 6  ! initialized to be us_standard_atmosphere

   ! 1: Tropical           (-23.4378 ~ 23.4378)
   if ( lat >= -23.4378 .and. lat <= 23.4378  ) crtm_climat = 1

   ! 2: Midlatitude summer (23.4378 ~ 66.561)
   ! 3: Midlatitude winter
   if ( lat > 23.4378 .and. lat <= 66.561  ) then ! North Mid-Lat
      if ( month >= 3 .and. month <= 8 ) then     
         crtm_climat = 2  ! Summer
      else 
         crtm_climat = 3  ! Winter
      end if
   end if

   if ( lat < -23.4378 .and. lat >= -66.561  ) then ! South Mid-Lat
      if ( month >= 3 .and. month <= 8 ) then
         crtm_climat = 3  ! Winter
      else
         crtm_climat = 2  ! Summer
      end if
   end if

   ! 4: Subarctic summer   ( > 66.561)
   ! 5: Subarctic winter
   if ( lat > 66.561  ) then ! Subarctic
      if ( month >= 3 .and. month <= 8 ) then
         crtm_climat = 4  ! Summer
      else
         crtm_climat = 5  ! Winter
      end if
   end if

   if ( lat < -66.561  ) then ! Subantarctic
      if ( month >= 3 .and. month <= 8 ) then
         crtm_climat = 5  ! Winter
      else
         crtm_climat = 4  ! Summer
      end if
   end if

end subroutine da_det_crtm_climat


end module da_crtm

