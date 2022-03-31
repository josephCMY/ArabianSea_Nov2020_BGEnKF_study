












module da_define_structures

   
   
   
   
   
   
   

   use module_domain, only: vp_type, x_type

   use da_control, only : anal_type_randomcv, stdout, max_fgat_time, &
      vert_corr, global, vert_evalue,print_detail_be, maxsensor, &
      max_ob_levels, trace_use, num_ob_indexes, kms, kme, &
      vert_corr_1, vert_corr_2, vert_evalue_global, cv_options, do_normalize, use_rf, &
      put_rand_seed, seed_array1, seed_array2, missing_r, &
      sound, synop, pilot, satem, geoamv, polaramv, airep, gpspw, gpsref, &
      metar, ships, ssmi_rv, ssmi_tb, ssmt1, ssmt2, qscat, profiler, buoy, bogus, &
      mtgirs, tamdar, tamdar_sfc, pseudo, radar, radiance, airsr, sonde_sfc, rain, &
      trace_use_dull,comm, num_pseudo

   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_tools_serial, only : da_array_print

   use da_reporting, only : da_error, da_warning, da_message, message
   use da_wavelet, only : nij,ws

   implicit none
   
   
   
   

   type xbx_type
      character (len=256) :: mminlu

      integer          :: fft_pad_i          
      integer          :: fft_pad_j          

      integer          :: pad_num            
      integer          :: pad_inc            
      integer, pointer :: pad_loc(:)         
      integer, pointer :: pad_pos(:)         

      integer          :: fft_ix             
      integer          :: fft_jy             

      integer, pointer :: fft_factors_x(:)   
      integer, pointer :: fft_factors_y(:)   

      real, pointer    :: trig_functs_x(:)   
      real, pointer    :: trig_functs_y(:)   

      real             :: psac_mean          
      real, pointer    :: latc_mean(:)       

      real, pointer    :: fft_coeffs(:,:)    

      real             :: fft_adjoint_factor 
      
      integer          :: inc                
      integer          :: ni
      integer          :: nj
      integer          :: nk
      integer          :: max_wavenumber
      integer          :: lenr
      integer          :: lensav
      integer          :: lenwrk
      integer          :: alp_size
      real, pointer       :: wsave(:)          
      real, pointer       :: lon(:)            
      real, pointer       :: sinlon(:)         
      real, pointer       :: coslon(:)         
      real, pointer       :: lat(:)            
      real, pointer       :: sinlat(:)         
      real, pointer       :: coslat(:)         
      real, pointer       :: int_wgts(:)       
      real, pointer       :: alp(:)            
   end type xbx_type

   
   
   

   

   type field_type
      real                   :: inv             
      integer                :: qc              
      real                   :: error           
      real                   :: sens            
      real                   :: imp             
   end type field_type

   type model_loc_type
      type (field_type)       :: slp            
      
      
      type (field_type)       :: pw             

      real                    :: x
      real                    :: y
      integer                 :: i
      integer                 :: j
      real                    :: dx
      real                    :: dxm
      real                    :: dy
      real                    :: dym
      logical                 :: proc_domain
      
      
      
      
      integer                 :: obs_global_index
   end type model_loc_type

   type each_level_type
      real                    :: height         
      integer                 :: height_qc      
      real                    :: zk             
      type (field_type)       :: u              
      type (field_type)       :: v              
      type (field_type)       :: p              
      type (field_type)       :: t              
      type (field_type)       :: q              
      type (field_type)       :: rh             
      type (field_type)       :: td             
      type (field_type)       :: Speed          
   end type each_level_type

   type radar_each_level_type
      real                   :: height         
      integer                :: height_qc      
      real                   :: zk             
      type (field_type)      :: rv
      type (field_type)      :: rf
   end type radar_each_level_type

   type info_type
      character (len = 40)   :: name          
      character (len = 12)   :: platform      
      character (len = 40)   :: id            
      character (len = 19)   :: date_char     
      integer                :: levels        
      real                   :: lat           
      real                   :: lon           
      real                   :: elv           
      real                   :: pstar         
      real                   :: dhr           
   end type info_type

   type infa_type
      integer                             :: max_lev
      integer                             :: nlocal
      integer                             :: ntotal
      integer                             :: thin_nlocal
      integer                             :: thin_ntotal
      integer                             :: plocal(0:max_fgat_time)
      integer                             :: ptotal(0:max_fgat_time)
      integer                             :: thin_plocal(0:max_fgat_time)
      integer                             :: thin_ptotal(0:max_fgat_time)
      integer                             :: n1
      integer                             :: n2
      character (len = 40) , allocatable  :: name(:)       
      character (len = 12), allocatable   :: platform(:)   
      character (len = 40), allocatable   :: id(:)         
      character (len = 19), allocatable   :: date_char(:)  
      integer, allocatable                :: levels(:)     
      real, allocatable                   :: lat(:,:)      
      real, allocatable                   :: lon(:,:)      
      real, allocatable                   :: elv(:)        
      real, allocatable                   :: pstar(:)      
      type (field_type), allocatable :: slp(:)         
      
      
      type (field_type), allocatable :: pw(:)          

      real, allocatable       :: x  (:,:)
      real, allocatable       :: y  (:,:)
      integer, allocatable    :: i  (:,:)
      integer, allocatable    :: j  (:,:)
      integer, allocatable    :: k  (:,:)
      real, allocatable       :: dx (:,:)
      real, allocatable       :: dxm(:,:)
      real, allocatable       :: dy (:,:)
      real, allocatable       :: dym(:,:)
      real, allocatable       :: dz (:,:)
      real, allocatable       :: dzm(:,:)
      real, allocatable       :: zk(:,:)
      logical, allocatable    :: proc_domain(:,:)
      logical, allocatable    :: thinned(:,:)
      
      
      
      
      integer, allocatable                 :: obs_global_index(:)
   end type infa_type

   type stn_loc_type
      real                    :: lon                  
      real                    :: lat                  
      real                    :: elv                  
      real                    :: x                    
      real                    :: y                    
      real                    :: zk                   
   end type stn_loc_type
 
   type radar_type
      type (stn_loc_type)     :: stn_loc

      real, pointer           :: model_p(:)
      real, pointer           :: model_t(:)
      real, pointer           :: model_rho(:)
      real, pointer           :: model_qrn(:)
      real, pointer           :: model_qcl(:)
      real, pointer           :: model_qci(:)
      real, pointer           :: model_qsn(:)
      real, pointer           :: model_qgr(:)
      real                    :: model_ps

      real                  , pointer :: height   (:) 
      integer               , pointer :: height_qc(:) 

      type (field_type)     , pointer :: rv       (:) 
      type (field_type)     , pointer :: rf       (:) 
      type (field_type)     , pointer :: rrn      (:) 
      type (field_type)     , pointer :: rcl      (:) 
      type (field_type)     , pointer :: rci      (:) 
      type (field_type)     , pointer :: rsn      (:) 
      type (field_type)     , pointer :: rgr      (:) 
      type (field_type)     , pointer :: rqv      (:) 
      real                  , pointer :: rrno     (:)
      real                  , pointer :: rclo     (:)
      real                  , pointer :: rcio     (:)
      real                  , pointer :: rsno     (:)
      real                  , pointer :: rgro     (:)
      real                  , pointer :: rqvo     (:)
   end type radar_type

   type multi_level_type
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (each_level_type)                  :: each(max_ob_levels)
   end type multi_level_type

   type multi_level_type_BUFR
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (each_level_type), pointer         :: each(:)
   end type multi_level_type_BUFR

   type radar_stn_type
      character (len = 5)    :: platform      
      character (len = 12)   :: name          
      character (len = 19)   :: date_char     
      integer                :: numobs        
      integer                :: levels        
      real                   :: lat           
      real                   :: lon           
      real                   :: elv           
   end type radar_stn_type

   type radar_multi_level_type
      type (radar_stn_type)                   :: stn
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (radar_each_level_type)            :: each(max_ob_levels)
   end type radar_multi_level_type

   type rain_stn_type
      character (len = 5)    :: platform      
      character (len = 12)   :: name          
      character (len = 19)   :: date_char     
      integer                :: numobs        
      integer                :: levels        
      real                   :: lat           
      real                   :: lon           
      real                   :: elv           
   end type rain_stn_type

   type rain_type
      real                    :: height    
      integer                 :: height_qc 
      type (stn_loc_type)     :: stn_loc
      type (field_type)       :: model_rainc    
      type (field_type)       :: model_rainnc   
      type (field_type)       :: rain             
   end type rain_type

   type rain_each_type
      real                   :: height         
      integer                :: height_qc      
      real                   :: zk             
      type (field_type)      :: rain
   end type rain_each_type

   type rain_single_level_type
      type (rain_stn_type)                    :: stn
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (rain_each_type)                   :: each(1)
   end type rain_single_level_type

   

   type airep_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
      type (field_type)     , pointer :: t        (:) 
      type (field_type)     , pointer :: q        (:) 
   end type airep_type

   type pilot_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
   end type pilot_type

   type bogus_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
      type (field_type)     , pointer :: t        (:) 
      type (field_type)     , pointer :: q        (:) 
      type (field_type)               :: slp          
   end type bogus_type

   type satem_type
      real                            :: ref_p        
      real                  , pointer :: p        (:) 

      type (field_type)     , pointer :: thickness(:)     
      type (field_type)     , pointer :: org_thickness(:) 
   end type satem_type

   type geoamv_type
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
   end type geoamv_type

   type polaramv_type
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
   end type polaramv_type

   type gpsref_type
      real             , pointer :: h  (:)      
      type (field_type), pointer :: ref(:)      
      type (field_type), pointer :: p  (:)      
      type (field_type), pointer :: t  (:)      
      type (field_type), pointer :: q  (:)      
   end type gpsref_type

   type synop_type
      real                    :: h              
      type (field_type)       :: u              
      type (field_type)       :: v              
      type (field_type)       :: t              
      type (field_type)       :: p              
      type (field_type)       :: q              
   end type synop_type

   type sound_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 

      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
      type (field_type)     , pointer :: t        (:) 
      type (field_type)     , pointer :: q        (:) 
   end type sound_type
     
   type mtgirs_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 

      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
      type (field_type)     , pointer :: t        (:) 
      type (field_type)     , pointer :: q        (:) 
   end type mtgirs_type

   type tamdar_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 

      type (field_type)     , pointer :: u        (:) 
      type (field_type)     , pointer :: v        (:) 
      type (field_type)     , pointer :: t        (:) 
      type (field_type)     , pointer :: q        (:) 
   end type tamdar_type

   type airsr_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: t        (:) 
      type (field_type)     , pointer :: q        (:) 
   end type airsr_type

   type gpspw_type
      type (field_type)       :: tpw  
  end type gpspw_type

   type ssmi_rv_type
      type (field_type)       :: Speed          
      type (field_type)       :: tpw            
   end type ssmi_rv_type

   type ssmi_tb_type

      type (field_type)       :: tb19v          
      type (field_type)       :: tb19h          
      type (field_type)       :: tb22v          
      type (field_type)       :: tb37v          
      type (field_type)       :: tb37h          
      type (field_type)       :: tb85v          
      type (field_type)       :: tb85h          
   end type ssmi_tb_type
   
   type ssmt1_type   
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: t        (:) 
   end type ssmt1_type

   type ssmt2_type
      real                  , pointer :: h        (:) 
      real                  , pointer :: p        (:) 
      type (field_type)     , pointer :: rh       (:) 
   end type ssmt2_type

   type pseudo_type
      type (field_type)       :: u              
      type (field_type)       :: v              
      type (field_type)       :: t              
      type (field_type)       :: p              
      type (field_type)       :: q              
   end type pseudo_type

   type qscat_type
      real                    :: h              
      type (field_type)       :: u              
      type (field_type)       :: v              
   end type qscat_type

   type varbc_info_type
      integer              :: platform_id, satellite_id, sensor_id
      integer              :: npredmax
      integer              :: gammapred
      integer              :: nchanl
      integer, pointer     :: nbgerr(:) 
      real,    pointer     :: pred(:,:)
      real,    pointer     :: pred_mean(:)
      real,    pointer     :: pred_std(:)
   end type varbc_info_type
   
   type varbc_type
      integer              :: nobs
      integer              :: npred 
      integer              :: ichanl
      integer, pointer     :: pred_use(:)
      integer, pointer     :: ipred(:)
      integer, pointer     :: index(:)
      real,    pointer     :: param(:)
      real,    pointer     :: bgerr(:) 
      real,    pointer     :: vtox(:,:)
   end type varbc_type
   
   type cv_index_type
      integer              :: ts
      integer              :: nclouds
      integer              :: ncv
      integer, pointer     :: cc(:)
      real, pointer        :: vtox(:,:)
   end type cv_index_type

   type instid_type
      
      integer              :: platform_id, satellite_id, sensor_id
      integer              :: rad_monitoring 
                                             
                                             
                                             
      character(len=20)    :: rttovid_string
      character(len=20)    :: rttovid_string_coef
      integer              :: num_rad, nchan, nlevels
      integer              :: num_rad_glo
      integer, pointer     :: ichan(:)
      real,    pointer     :: tb_inv(:,:)
      integer, pointer     :: tb_qc(:,:)
      real,    pointer     :: tb_error(:,:)
      real,    pointer     :: tb_xb(:,:) 
      real,    pointer     :: tb_sens(:,:)
      real,    pointer     :: tb_imp(:,:)
      real,    pointer     :: rad_xb(:,:)
      real,    pointer     :: rad_obs(:,:)
      real,    pointer     :: rad_ovc(:,:,:)
      integer, pointer     :: scanpos(:)
      integer, pointer     :: scanline(:)
      integer, pointer     :: cloud_flag(:,:)
      integer, pointer     :: rain_flag(:)
      real,    pointer     :: satzen(:) 
      real,    pointer     :: satazi(:) 
      real,    pointer     :: solzen(:) 
      real,    pointer     :: solazi(:) 
      real,    pointer     :: t(:,:)
      real,    pointer     :: q(:,:)
      real,    pointer     :: mr(:,:)
      real,    pointer     :: tm(:,:)
      real,    pointer     :: qm(:,:)
      real,    pointer     :: lod(:,:,:)       
      real,    pointer     :: trans(:,:,:)     
      real,    pointer     :: der_trans(:,:,:) 
      real,    pointer     :: kmin_t(:)	
      real,    pointer     :: kmax_p(:)	  
      real,    pointer     :: sensitivity_ratio(:,:,:)	  
      real,    pointer     :: p_chan_level(:,:)	  
      real,    pointer     :: qrn(:,:)
      real,    pointer     :: qcw(:,:)
      real,    pointer     :: qci(:,:)
      real,    pointer     :: qsn(:,:)
      real,    pointer     :: qgr(:,:)
      real,    pointer     :: qhl(:,:)
      real,    pointer     :: pm(:,:)
      real,    pointer     :: rcw(:,:) 
      real,    pointer     :: rci(:,:) 
      real,    pointer     :: rrn(:,:) 
      real,    pointer     :: rsn(:,:) 
      real,    pointer     :: rgr(:,:) 
      real,    pointer     :: rhl(:,:) 
      real,    pointer     :: pf(:,:)  
      real,    pointer     :: emiss(:,:)
      real,    pointer     :: u10(:)
      real,    pointer     :: v10(:)
      real,    pointer     :: t2m(:)
      real,    pointer     :: q2m(:)
      real,    pointer     :: mr2m(:)
      real,    pointer     :: psfc(:)
      real,    pointer     :: ps(:)
      real,    pointer     :: ts(:)
      real,    pointer     :: smois(:)
      real,    pointer     :: tslb(:)
      real,    pointer     :: snowh(:)
      integer, pointer     :: isflg(:)
      integer, pointer     :: ifgat(:)
      integer, pointer     :: landsea_mask(:)
      integer, pointer     :: surftype(:)     
      real,    pointer     :: snow_frac(:)    
      real,    pointer     :: elevation(:)
      real,    pointer     :: soiltyp(:)
      real,    pointer     :: vegtyp(:)
      real,    pointer     :: vegfra(:)
      real,    pointer     :: clwp(:) 
      real,    pointer     :: clw(:)  
      real,    pointer     :: ps_jacobian(:,:) 
      real,    pointer     :: ts_jacobian(:,:) 
      real,    pointer     :: windspeed_jacobian(:,:) 
      real,    pointer     :: emiss_jacobian(:,:)
      real,    pointer     :: gamma_jacobian(:,:)
      real,    pointer     :: t_jacobian(:,:,:)
      real,    pointer     :: q_jacobian(:,:,:)
      real,    pointer     :: lod_jacobian(:,:,:)
      real,    pointer     :: trans_jacobian(:,:,:)
      real,    pointer     :: water_jacobian(:,:,:) 
      real,    pointer     :: ice_jacobian(:,:,:)
      real,    pointer     :: rain_jacobian(:,:,:)
      real,    pointer     :: snow_jacobian(:,:,:)
      real,    pointer     :: graupel_jacobian(:,:,:)
      real,    pointer     :: hail_jacobian(:,:,:)
      real,    pointer     :: water_r_jacobian(:,:,:) 
      real,    pointer     :: ice_r_jacobian(:,:,:)
      real,    pointer     :: rain_r_jacobian(:,:,:)
      real,    pointer     :: snow_r_jacobian(:,:,:)
      real,    pointer     :: graupel_r_jacobian(:,:,:)
      real,    pointer     :: hail_r_jacobian(:,:,:)
      real,    pointer     :: water_coverage(:)
      real,    pointer     :: land_coverage(:)
      real,    pointer     :: ice_coverage(:)
      real,    pointer     :: snow_coverage(:)
      integer, pointer     :: crtm_climat(:) 

      type (varbc_info_type)        :: varbc_info
      type (varbc_type),pointer     :: varbc(:)
      type (cv_index_type), pointer :: cv_index(:)
      type (infa_type)              :: info
   end type instid_type

   type iv_type
      integer :: nstats(num_ob_indexes)

      integer :: time

      integer :: num_inst, total_rad_pixel, total_rad_channel

      real    :: synop_ef_u, synop_ef_v, synop_ef_t, synop_ef_p, synop_ef_q
      real    :: metar_ef_u, metar_ef_v, metar_ef_t, metar_ef_p, metar_ef_q
      real    :: ships_ef_u, ships_ef_v, ships_ef_t, ships_ef_p, ships_ef_q
      real    :: geoamv_ef_u, geoamv_ef_v
      real    :: polaramv_ef_u, polaramv_ef_v
      real    :: gpspw_ef_tpw
      real    :: sound_ef_u, sound_ef_v, sound_ef_t, sound_ef_q
      real    :: mtgirs_ef_u, mtgirs_ef_v, mtgirs_ef_t, mtgirs_ef_q
      real    :: tamdar_ef_u, tamdar_ef_v, tamdar_ef_t, tamdar_ef_q
      real    :: tamdar_sfc_ef_u, tamdar_sfc_ef_v, tamdar_sfc_ef_t, tamdar_sfc_ef_p, tamdar_sfc_ef_q
      real    :: airep_ef_u, airep_ef_v, airep_ef_t, airep_ef_q
      real    :: pilot_ef_u, pilot_ef_v
      real    :: ssmir_ef_speed, ssmir_ef_tpw
      real    :: satem_ef_thickness, ssmt1_ef_t, ssmt2_ef_rh
      real    :: gpsref_ef_ref, gpsref_ef_p, gpsref_ef_t, gpsref_ef_q
      real    :: qscat_ef_u, qscat_ef_v
      real    :: profiler_ef_u, profiler_ef_v
      real    :: buoy_ef_u, buoy_ef_v, buoy_ef_t, buoy_ef_p, buoy_ef_q
      real    :: radar_ef_rv, radar_ef_rf, radar_ef_rr
      real    :: bogus_ef_u, bogus_ef_v, bogus_ef_t, bogus_ef_p, bogus_ef_q, bogus_ef_slp
      real    :: airsr_ef_t,  airsr_ef_q
      real    :: rain_ef_r

      type (infa_type) :: info(num_ob_indexes)

      type (airsr_type)    , pointer :: airsr(:)
      type (sound_type)    , pointer :: sound(:)
      type (synop_type)    , pointer :: sonde_sfc(:)
      type (airep_type)    , pointer :: airep(:)
      type (pilot_type)    , pointer :: pilot(:)
      type (satem_type)    , pointer :: satem(:)
      type (geoamv_type)   , pointer :: geoamv(:)
      type (polaramv_type) , pointer :: polaramv(:)
      type (synop_type)    , pointer :: synop(:)
      type (synop_type)    , pointer :: metar(:)
      type (synop_type)    , pointer :: ships(:)
      type (gpspw_type)    , pointer :: gpspw(:)
      type (gpsref_type)   , pointer :: gpsref(:)
      type (ssmi_tb_type)  , pointer :: ssmi_tb(:)
      type (ssmi_rv_type)  , pointer :: ssmi_rv(:)
      type (ssmt1_type)    , pointer :: ssmt1(:)
      type (ssmt2_type)    , pointer :: ssmt2(:)
      type (pseudo_type)   , pointer :: pseudo(:)
      type (qscat_type)    , pointer :: qscat(:)
      type (synop_type)    , pointer :: buoy(:)
      type (pilot_type)    , pointer :: profiler(:)
      type (bogus_type)    , pointer :: bogus(:)
      type (radar_type)    , pointer :: radar(:)
      type (instid_type)   , pointer :: instid(:)
      type (mtgirs_type)   , pointer :: mtgirs(:)
      type (tamdar_type)   , pointer :: tamdar(:)
      type (synop_type)    , pointer :: tamdar_sfc(:)
      type (rain_type)     , pointer :: rain(:)

      real :: missing
      real :: ptop
   end type iv_type

   type number_type
      integer                    :: bad
      integer                    :: miss
      integer                    :: use
   end type number_type

   type bad_info_type
      type (number_type)         :: num
      integer                    :: nn(100000)
      integer                    :: kk(100000)
   end type bad_info_type

   type bad_data_type
      type (bad_info_type)       :: u
      type (bad_info_type)       :: v
      type (bad_info_type)       :: t
      type (bad_info_type)       :: p
      type (bad_info_type)       :: q
      type (bad_info_type)       :: tpw
      type (bad_info_type)       :: Speed
      type (bad_info_type)       :: gpsref
      type (bad_info_type)       :: thickness
      type (bad_info_type)       :: rh
      type (bad_info_type)       :: rv
      type (bad_info_type)       :: rf
      type (bad_info_type)       :: rrn
      type (bad_info_type)       :: rsn
      type (bad_info_type)       :: rgr
      type (bad_info_type)       :: rcl
      type (bad_info_type)       :: rci
      type (bad_info_type)       :: rqv
      type (bad_info_type)       :: slp
      type (bad_info_type)       :: rad
      type (bad_info_type)       :: rain
   end type bad_data_type

   type count_obs_number_type
      integer                                 :: num_used
      integer                                 :: num_outside_iyjx
      integer                                 :: num_max_err_chk
      integer                                 :: num_missing
   end type count_obs_number_type

   
   
   

   type residual_synop_type
      real :: u                                 
      real :: v                                 
      real :: t                                 
      real :: p                                 
      real :: q                                 
   end type residual_synop_type

   type residual_qscat_type
      real :: u                                 
      real :: v                                 
   end type residual_qscat_type

   type residual_geoamv_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
   end type residual_geoamv_type

   type residual_polaramv_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
   end type residual_polaramv_type

   type residual_gpspw_type
      real :: tpw                               
   end type residual_gpspw_type

   type residual_sound_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
      real, pointer :: t(:)                     
      real, pointer :: q(:)                     
   end type residual_sound_type
     
   type residual_mtgirs_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
      real, pointer :: t(:)                     
      real, pointer :: q(:)                     
   end type residual_mtgirs_type

   type residual_tamdar_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
      real, pointer :: t(:)                     
      real, pointer :: q(:)                     
   end type residual_tamdar_type

   type residual_airsr_type
      real, pointer :: t(:)                     
      real, pointer :: q(:)                     
   end type residual_airsr_type

   type residual_airep_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
      real, pointer :: t(:)                     
      real, pointer :: q(:)                     
   end type residual_airep_type

   type residual_pilot_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
   end type residual_pilot_type

   type residual_bogus_type
      real, pointer :: u(:)                     
      real, pointer :: v(:)                     
      real, pointer :: t(:)                     
      real, pointer :: q(:)                     
      real          :: slp                      
   end type residual_bogus_type

   type residual_satem_type
      real, pointer :: thickness(:)             
   end type residual_satem_type

   type residual_gpsref_type
      real, pointer :: ref(:)         
      real, pointer :: p  (:)         
      real, pointer :: t  (:)         
      real, pointer :: q  (:)         
   end type residual_gpsref_type

   type residual_ssmi_rv_type
      real                    :: tpw      
      real                    :: Speed    
   end type residual_ssmi_rv_type

   type residual_ssmi_tb_type
      real                    :: tb19v          
      real                    :: tb19h          
      real                    :: tb22v          
      real                    :: tb37v          
      real                    :: tb37h          
      real                    :: tb85v          
      real                    :: tb85h          
   end type residual_ssmi_tb_type
   
   type residual_ssmt1_type
      real, pointer :: t(:)                       
   end type residual_ssmt1_type
   
   type residual_ssmt2_type
      real, pointer :: rh(:)                      
   end type residual_ssmt2_type

   type residual_pseudo_type
      real :: u                                   
      real :: v                                   
      real :: t                                   
      real :: p                                   
      real :: q                                   
   end type residual_pseudo_type

   type residual_radar_type
      real, pointer :: rv(:)                    
      real, pointer :: rf(:)                    
      real, pointer :: rrn(:)                   
      real, pointer :: rcl(:)                   
      real, pointer :: rci(:)                   
      real, pointer :: rsn(:)                   
      real, pointer :: rgr(:)                   
      real, pointer :: rqv(:) 
   end type residual_radar_type

   type residual_instid_type
      integer                          :: num_rad
      integer                          :: nchan
      integer, pointer                 :: ichan (:)
      real, pointer                    :: tb(:,:)
   end type residual_instid_type

   type residual_rain_type
      real :: rain
   end type residual_rain_type 

   type y_type
      integer :: nlocal(num_ob_indexes)
      integer :: ntotal(num_ob_indexes)

      integer :: num_inst

      type (residual_synop_type),    pointer :: synop(:)
      type (residual_synop_type),    pointer :: metar(:) 
      type (residual_synop_type),    pointer :: ships(:) 
      type (residual_geoamv_type),   pointer :: geoamv(:)
      type (residual_polaramv_type), pointer :: polaramv(:)
      type (residual_gpspw_type),    pointer :: gpspw (:)
      type (residual_gpsref_type),   pointer :: gpsref(:)
      type (residual_sound_type),    pointer :: sound(:)
      type (residual_mtgirs_type),   pointer :: mtgirs(:)
      type (residual_tamdar_type),   pointer :: tamdar(:)
      type (residual_synop_type),    pointer :: tamdar_sfc(:)
      type (residual_airsr_type),    pointer :: airsr(:)
      type (residual_bogus_type),    pointer :: bogus(:)
      type (residual_synop_type),    pointer :: sonde_sfc(:) 
      type (residual_airep_type),    pointer :: airep(:)
      type (residual_pilot_type),    pointer :: pilot(:)
      type (residual_satem_type),    pointer :: satem(:)
      type (residual_ssmi_tb_type),  pointer :: ssmi_tb(:)
      type (residual_ssmi_rv_type),  pointer :: ssmi_rv(:)
      type (residual_ssmt1_type),    pointer :: ssmt1(:)
      type (residual_ssmt2_type),    pointer :: ssmt2(:)
      type (residual_pseudo_type),   pointer :: pseudo(:)
      type (residual_qscat_type),    pointer :: qscat(:)
      type (residual_synop_type),    pointer :: buoy(:) 
      type (residual_pilot_type),    pointer :: profiler(:) 
      type (residual_radar_type),    pointer :: radar(:)
      type (residual_instid_type),   pointer :: instid(:)
      type (residual_rain_type),     pointer :: rain(:)
   end type y_type

   
   
   

   

   type maxmin_type
      real                       :: value
      integer                    :: n, l
   end type maxmin_type

   
   
   
   
   type jo_type_rad
      integer, pointer :: num_ichan(:)
      real, pointer    :: jo_ichan(:)
   end type jo_type_rad

   type jo_type
      real                :: total
      real                :: synop_u, synop_v, synop_t, synop_p, synop_q
      real                :: metar_u, metar_v, metar_t, metar_p, metar_q
      real                :: ships_u, ships_v, ships_t, ships_p, ships_q
      real                :: geoamv_u, geoamv_v
      real                :: polaramv_u, polaramv_v
      real                :: gpspw_tpw, satem_thickness, gpsref_ref
      real                :: sound_u, sound_v, sound_t, sound_q
      real                :: sonde_sfc_u, sonde_sfc_v, sonde_sfc_t, &
                             sonde_sfc_p, sonde_sfc_q
      real                :: mtgirs_u, mtgirs_v, mtgirs_t, mtgirs_q
      real                :: tamdar_u, tamdar_v, tamdar_t, tamdar_q
      real                :: tamdar_sfc_u, tamdar_sfc_v, tamdar_sfc_t, &
                             tamdar_sfc_p, tamdar_sfc_q
      real                :: airep_u, airep_v, airep_t, airep_q
      real                :: pilot_u, pilot_v
      real                :: ssmir_speed, ssmir_tpw
      real                :: ssmi_tb19v, ssmi_tb19h, ssmi_tb22v, ssmi_tb37v, &
                             ssmi_tb37h, ssmi_tb85v, ssmi_tb85h
      real                :: ssmt1_t, ssmt2_rh
      real                :: pseudo_u, pseudo_v, pseudo_t, pseudo_p, pseudo_q
      real                :: qscat_u, qscat_v
      real                :: profiler_u, profiler_v
      real                :: buoy_u, buoy_v, buoy_t, buoy_p, buoy_q
      real                :: radar_rv, radar_rf, radar_rrn,radar_rsn,radar_rgr,radar_rcl,radar_rci,radar_rqv
      real                :: bogus_u, bogus_v, bogus_t, bogus_q, bogus_slp
      real                :: airsr_t, airsr_q
      real                :: rain_r
      type(jo_type_rad), pointer       :: rad(:)
   end type jo_type

   type j_type
      real             :: total
      real             :: jb
      real             :: jc
      real             :: je
      real             :: jp
      real             :: js
      real             :: jl
      real             :: jd
      type (jo_type)   :: jo
   end type j_type

   type cv_type
      integer :: size        
      integer :: size_jb     
      integer :: size_je     
      integer :: size_jp     
      integer :: size_js     
      integer :: size_jl     
      integer :: size1c      
      integer :: size2c      
      integer :: size3c      
      integer :: size4c      
      integer :: size5c      
      integer :: size_alphac 
      integer :: size1       
      integer :: size2       
      integer :: size3       
      integer :: size4       
      integer :: size5       
      integer :: size1l      
      integer :: size2l      
      integer :: size3l      
      integer :: size4l      
      integer :: size5l      
   end type cv_type

   type qhat_type
      integer          :: i
      real, allocatable:: values(:) 
   end type qhat_type

   type be_subtype
      integer           :: mz          
      integer           :: max_wave    
      character*5       :: name        
      real*8, pointer   :: rf_alpha(:) 
      real*8, pointer   :: val(:,:)    
      real*8, pointer   :: evec(:,:,:) 
      real*8, pointer   :: val_g(:)    
      real*8, pointer   :: evec_g(:,:) 
      real*8, pointer   :: power(:,:)  

      REAL, POINTER     ::sd(:,:,:)    
      REAL, POINTER     ::wsd(:,:,:)   
   end type be_subtype

   type be_type
      integer           :: ne
      integer           :: max_wave           
      integer           :: mix
      integer           :: mjy
      type (be_subtype) :: v1
      type (be_subtype) :: v2
      type (be_subtype) :: v3
      type (be_subtype) :: v4
      type (be_subtype) :: v5
      type (be_subtype) :: alpha
      real*8, pointer     :: pb_vert_reg(:,:,:)

      
      type (cv_type)    :: cv

      real, pointer :: reg_psi_chi  (:,:)
      real, pointer :: reg_psi_t (:,:,:)
      real, pointer :: reg_psi_ps   (:,:)
      real, pointer :: reg_psi_rh   (:,:,:)
      real, pointer :: reg_chi_u_t   (:,:,:)
      real, pointer :: reg_chi_u_ps     (:,:)
      real, pointer :: reg_chi_u_rh     (:,:,:)
      real, pointer :: reg_t_u_rh    (:,:,:)
      real, pointer :: reg_ps_u_rh      (:,:)


      INTEGER          :: ndeg,nta
      REAL             :: swidth
      REAL, POINTER    :: be(:)
      REAL, POINTER    :: rate(:)
      REAL, POINTER    :: table(:,:)
      REAL, POINTER    :: agvz(:,:,:,:)
      REAL, POINTER    :: bvz(:,:,:)
      REAL, POINTER    :: wgvz(:,:,:)
      REAL, POINTER    :: slix(:,:,:,:)
      REAL, POINTER    :: slipx(:,:)
      REAL, POINTER    :: sljy(:,:,:,:)
      REAL, POINTER    :: sljpy(:,:)
      REAL, POINTER    :: vz(:,:,:,:)
      REAL, POINTER    :: corz(:,:,:,:)
      REAL, POINTER    :: corp(:,:)


      REAL, POINTER    ::sd( :,:,:)
      REAL, POINTER    ::wsd(:,:,:)
   end type be_type

   

   type maxmin_field_type
      real                         :: value
      integer                      :: i, j
   end type maxmin_field_type

   
   
   
   
   
   

contains

subroutine da_allocate_background_errors  (jy, kz, l, e, be_eval_loc, &
                                           be_evec_loc, be_sub)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate components of wrfvar background errors.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)              :: jy                 ! i/y dimension.
   integer, intent(in)              :: kz                 ! k/z dimension.
   real*8, intent(in)               :: l(:)               ! Global eigenvalue.
   real*8, intent(in)               :: e(:,:)             ! Global eigenvectors.
   real*8, intent(in)               :: be_eval_loc(:,:)   ! Std dev/local evalue.
   real*8, intent(in)               :: be_evec_loc(:,:,:) ! Local eigenvectors.
   type (be_subtype), intent(inout) :: be_sub             ! Backgrd error struct.
    
   integer                          :: mz                 ! Vertical truncation.
   integer                          :: j, m, k            ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_allocate_background_errors")

   !--------------------------------------------------------------------------
   ! [1.0] Initialise:
   !--------------------------------------------------------------------------

   mz = be_sub % mz
   
   !--------------------------------------------------------------------------
   ! [2.0] Allocate components of be_sub structure:
   !--------------------------------------------------------------------------

   if (mz > 0) then
      allocate  (be_sub % val(1:jy,1:mz))
      
      if (vert_corr == vert_corr_2) then

         !--------------------------------------------------------------------
         ! [3.0] Allocate eigenvalues of vertical error covariance matrix:
         !--------------------------------------------------------------------

         if (vert_evalue == vert_evalue_global) then
            ! use global eigenvalues:
            do m = 1, mz
               be_sub % val(1:jy,m) = sqrt (l(m))
            end do
         else
            ! use eigenvalues varying with j-direction.
            do j = 1, jy
               do k = 1, mz
                  if (be_eval_loc(j,k) <=0) then
                     write (unit=message(1),fmt='(A,I5,A,I5,A,F10.2)') &
                        "At lat= ",j," For mode = ",k," local eigen value= ",be_eval_loc(j,k)
                     call da_error("da_allocate_background_errors.inc",54,message(1:1))
                  end if
               end do
            end do
            be_sub % val(1:jy,1:mz) = sqrt (be_eval_loc(1:jy,1:mz))            
         end if
 
         if (print_detail_be) then
            write(unit=message(1),fmt='(A,A)') 'j*k Eigenvalues for ', be_sub % name
            call da_array_print(2, be_sub % val(1:jy,1:mz), message(1))
         end if

         !----------------------------------------------------------------------- 
         ! [4.0] Allocate global eigenvectors of vertical error cov.:
         !-----------------------------------------------------------------------

         allocate  (be_sub % evec(1:jy,1:kz,1:mz))
         
         if (vert_evalue == vert_evalue_global) then
            ! use global eigenvectors:
            do j = 1, jy
               be_sub % evec(j,1:kz,1:mz) = e(1:kz,1:mz)
            end do
         else
            ! use eigenvectors varying with i-direction.
            be_sub % evec(1:jy,1:kz,1:mz) =  be_evec_loc(1:jy,1:kz,1:mz)
         end if
         
         if (print_detail_be) then      
            write(unit=message(1),fmt='(A,A)') 'k*k Eigenvectors for j = 1 ', be_sub % name
            call da_array_print (2, be_sub % evec(1,1:kz,1:mz), message(1))
         
            write(unit=message(1),fmt='(A,A)') 'k*k Eigenvectors for j = jy ', be_sub % name
            call da_array_print (2, be_sub % evec(jy,1:kz,1:mz), message(1))
         end if

         allocate (be_sub%val_g(1:mz))
         allocate (be_sub%evec_g(1:kz,1:mz))
  
         be_sub % val_g(1:mz) = l(1:mz)
         be_sub % evec_g(1:kz,1:mz) = e(1:kz,1:mz)
      else if (vert_corr == vert_corr_1) then
         if (print_detail_be) then
           write(unit=message(1),fmt='(A)') 'Change BE std dev to variance in NMC code'
           call da_message(message(1:1))
         end if
         if (vert_evalue == vert_evalue_global) then
            ! use global eigenvalues:
          do m = 1, mz
          be_sub % val(1:jy,m) = l(m)
          end do 
         else
          be_sub % val(1:jy,1:mz) = be_eval_loc(1:jy,1:mz)
         end if
      end if

      !-----------------------------------------------------------------------
      ! [2.2] Allocate recursive filter lengthscales and variance scaling factors:
      !-----------------------------------------------------------------------

      allocate (be_sub % rf_alpha(1:mz))

      be_sub % rf_alpha(1:mz) = 1.0    
   end if

   if (trace_use_dull) call da_trace_exit("da_allocate_background_errors")

end subroutine da_allocate_background_errors


subroutine da_allocate_observations (iv)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate components of observation structure.
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv     ! Observation structure.

   integer :: i

   if (trace_use) call da_trace_entry("da_allocate_observations")

   if (iv%info(sound)%nlocal     > 0) allocate(iv%sound    (1:iv%info(sound)%nlocal))
   if (iv%info(sonde_sfc)%nlocal > 0) allocate(iv%sonde_sfc(1:iv%info(sonde_sfc)%nlocal))
   if (iv%info(mtgirs)%nlocal    > 0) allocate(iv%mtgirs   (1:iv%info(mtgirs)%nlocal))
   if (iv%info(tamdar)%nlocal    > 0) allocate(iv%tamdar   (1:iv%info(tamdar)%nlocal))
   if (iv%info(tamdar_sfc)%nlocal > 0) allocate(iv%tamdar_sfc (1:iv%info(tamdar_sfc)%nlocal))
   if (iv%info(synop)%nlocal     > 0) allocate(iv%synop    (1:iv%info(synop)%nlocal))
   if (iv%info(airep)%nlocal     > 0) allocate(iv%airep    (1:iv%info(airep)%nlocal))
   if (iv%info(geoamv)%nlocal    > 0) allocate(iv%geoamv   (1:iv%info(geoamv)%nlocal))
   if (iv%info(polaramv)%nlocal  > 0) allocate(iv%polaramv (1:iv%info(polaramv)%nlocal))
   if (iv%info(satem)%nlocal     > 0) allocate(iv%satem    (1:iv%info(satem)%nlocal))
   if (iv%info(metar)%nlocal     > 0) allocate(iv%metar    (1:iv%info(metar)%nlocal))
   if (iv%info(ships)%nlocal     > 0) allocate(iv%ships    (1:iv%info(ships)%nlocal))
   if (iv%info(pilot)%nlocal     > 0) allocate(iv%pilot    (1:iv%info(pilot)%nlocal))
   if (iv%info(gpspw)%nlocal     > 0) allocate(iv%gpspw    (1:iv%info(gpspw)%nlocal))
   if (iv%info(gpsref)%nlocal    > 0) allocate(iv%gpsref   (1:iv%info(gpsref)%nlocal))
   if (iv%info(ssmi_tb)%nlocal   > 0) allocate(iv%ssmi_tb  (1:iv%info(ssmi_tb)%nlocal))
   if (iv%info(ssmi_rv)%nlocal   > 0) allocate(iv%ssmi_rv  (1:iv%info(ssmi_rv)%nlocal))
   if (iv%info(ssmt1)%nlocal     > 0) allocate(iv%ssmt1    (1:iv%info(ssmt1)%nlocal))
   if (iv%info(ssmt2)%nlocal     > 0) allocate(iv%ssmt2    (1:iv%info(ssmt2)%nlocal))
   if (iv%info(qscat)%nlocal     > 0) allocate(iv%qscat    (1:iv%info(qscat)%nlocal))
   if (iv%info(profiler)%nlocal  > 0) allocate(iv%profiler (1:iv%info(profiler)%nlocal))
   if (iv%info(buoy)%nlocal      > 0) allocate(iv%buoy     (1:iv%info(buoy)%nlocal))
   if (iv%info(radar)%nlocal     > 0) allocate(iv%radar    (1:iv%info(radar)%nlocal))
   if (iv%info(bogus)%nlocal     > 0) allocate(iv%bogus    (1:iv%info(bogus)%nlocal))
   if (iv%info(airsr)%nlocal     > 0) allocate(iv%airsr    (1:iv%info(airsr)%nlocal))
   if (iv%info(pseudo)%nlocal    > 0) allocate(iv%pseudo   (1:iv%info(pseudo)%nlocal))
   if (iv%info(rain)%nlocal      > 0) allocate(iv%rain     (1:iv%info(rain)%nlocal))
   do i=1,num_ob_indexes

      ! Radar structures now allocated in their own subroutine (da_setup_obs_structures_radar)
      if (i == radar) cycle

      if (iv%info(i)%nlocal > 0) then
         allocate (iv%info(i)%name(iv%info(i)%nlocal))     
         allocate (iv%info(i)%platform(iv%info(i)%nlocal)) 
         allocate (iv%info(i)%id(iv%info(i)%nlocal))       
         allocate (iv%info(i)%date_char(iv%info(i)%nlocal))
         allocate (iv%info(i)%levels(iv%info(i)%nlocal))   
         allocate (iv%info(i)%lat(iv%info(i)%max_lev,iv%info(i)%nlocal))    
         allocate (iv%info(i)%lon(iv%info(i)%max_lev,iv%info(i)%nlocal))    
         allocate (iv%info(i)%elv(iv%info(i)%nlocal))      
         allocate (iv%info(i)%pstar(iv%info(i)%nlocal))    

         allocate (iv%info(i)%slp(iv%info(i)%nlocal))   
         allocate (iv%info(i)%pw(iv%info(i)%nlocal))    

         allocate (iv%info(i)%x  (kms:kme,iv%info(i)%nlocal))   
         allocate (iv%info(i)%y  (kms:kme,iv%info(i)%nlocal))   
         allocate (iv%info(i)%i  (kms:kme,iv%info(i)%nlocal))   
         allocate (iv%info(i)%j  (kms:kme,iv%info(i)%nlocal))      
         allocate (iv%info(i)%dx (kms:kme,iv%info(i)%nlocal))  
         allocate (iv%info(i)%dxm(kms:kme,iv%info(i)%nlocal)) 
         allocate (iv%info(i)%dy (kms:kme,iv%info(i)%nlocal))  
         allocate (iv%info(i)%dym(kms:kme,iv%info(i)%nlocal)) 
         allocate (iv%info(i)%k  (iv%info(i)%max_lev,iv%info(i)%nlocal))
         allocate (iv%info(i)%dz (iv%info(i)%max_lev,iv%info(i)%nlocal))  
         allocate (iv%info(i)%dzm(iv%info(i)%max_lev,iv%info(i)%nlocal)) 
         allocate (iv%info(i)%zk (iv%info(i)%max_lev,iv%info(i)%nlocal)) 
         allocate (iv%info(i)%proc_domain(iv%info(i)%max_lev,iv%info(i)%nlocal)) 
         allocate (iv%info(i)%thinned(iv%info(i)%max_lev,iv%info(i)%nlocal)) 
         allocate (iv%info(i)%obs_global_index(iv%info(i)%nlocal)) 

         iv%info(i)%proc_domain(:,:)  = .false.
         iv%info(i)%thinned(:,:)      = .false.
         iv%info(i)%zk(:,:)           = missing_r
      end if
   end do

   if (trace_use) call da_trace_exit("da_allocate_observations")

end subroutine da_allocate_observations


subroutine da_allocate_observations_rain (iv)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate components of rain observation structure.
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv     ! Observation structure.

   integer :: i

   if (trace_use) call da_trace_entry("da_allocate_observations_rain")

   if (iv%info(rain)%nlocal    > 0) allocate(iv%rain   (1:iv%info(rain)%nlocal))
   if (iv%info(rain)%nlocal > 0) then
      allocate (iv%info(rain)%name(iv%info(rain)%nlocal))     
      allocate (iv%info(rain)%platform(iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%id(iv%info(rain)%nlocal))       
      allocate (iv%info(rain)%date_char(iv%info(rain)%nlocal))
      allocate (iv%info(rain)%levels(iv%info(rain)%nlocal))   
      allocate (iv%info(rain)%lat(iv%info(rain)%max_lev,iv%info(rain)%nlocal))    
      allocate (iv%info(rain)%lon(iv%info(rain)%max_lev,iv%info(rain)%nlocal))    
      allocate (iv%info(rain)%elv(iv%info(rain)%nlocal))      
      allocate (iv%info(rain)%pstar(iv%info(rain)%nlocal))    

      allocate (iv%info(rain)%slp(iv%info(rain)%nlocal))   
      allocate (iv%info(rain)%pw(iv%info(rain)%nlocal))    

      allocate (iv%info(rain)%x  (kms:kme,iv%info(rain)%nlocal))   
      allocate (iv%info(rain)%y  (kms:kme,iv%info(rain)%nlocal))   
      allocate (iv%info(rain)%i  (kms:kme,iv%info(rain)%nlocal))   
      allocate (iv%info(rain)%j  (kms:kme,iv%info(rain)%nlocal))      
      allocate (iv%info(rain)%dx (kms:kme,iv%info(rain)%nlocal))  
      allocate (iv%info(rain)%dxm(kms:kme,iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%dy (kms:kme,iv%info(rain)%nlocal))  
      allocate (iv%info(rain)%dym(kms:kme,iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%k  (iv%info(rain)%max_lev,iv%info(rain)%nlocal))
      allocate (iv%info(rain)%dz (iv%info(rain)%max_lev,iv%info(rain)%nlocal))  
      allocate (iv%info(rain)%dzm(iv%info(rain)%max_lev,iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%zk (iv%info(rain)%max_lev,iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%proc_domain(iv%info(rain)%max_lev,iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%thinned(iv%info(rain)%max_lev,iv%info(rain)%nlocal)) 
      allocate (iv%info(rain)%obs_global_index(iv%info(rain)%nlocal)) 

      iv%info(rain)%proc_domain(:,:)  = .false.
      iv%info(rain)%thinned(:,:)      = .false.
      iv%info(rain)%zk(:,:)           = missing_r
   end if

   if (trace_use) call da_trace_exit("da_allocate_observations_rain")

end subroutine da_allocate_observations_rain


subroutine da_allocate_y (iv, y)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate arrays used in y and residual obs structures.
   !---------------------------------------------------------------------------

   implicit none
   
   type (iv_type), intent(in)    :: iv      ! Ob type input.
   type (y_type),  intent(inout) :: y       ! Residual type structure.

   integer :: n, i    ! Loop counter.
   integer :: nlevels ! Number of levels.

   !---------------------------------------------------------------------------
   !  [1.0] Copy number of observations:
   !---------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_allocate_y")

   y % nlocal(:) = iv%info(:)%nlocal
   y % ntotal(:) = iv%info(:)%ntotal

   y % num_inst     = iv % num_inst

  !---------------------------------------------------------------------------
  ! [2.0] Allocate:
  !---------------------------------------------------------------------------

   if (y % nlocal(synop) > 0) then
      allocate (y % synop(1:y % nlocal(synop)))
      y % synop(1:y % nlocal(synop)) % u = 0.0
      y % synop(1:y % nlocal(synop)) % v = 0.0
      y % synop(1:y % nlocal(synop)) % t = 0.0
      y % synop(1:y % nlocal(synop)) % p = 0.0
      y % synop(1:y % nlocal(synop)) % q = 0.0
   end if

   if (y % nlocal(ships) > 0) then
      allocate (y % ships(1:y % nlocal(ships)))
      y % ships(1:y % nlocal(ships)) % u = 0.0
      y % ships(1:y % nlocal(ships)) % v = 0.0
      y % ships(1:y % nlocal(ships)) % t = 0.0
      y % ships(1:y % nlocal(ships)) % p = 0.0
      y % ships(1:y % nlocal(ships)) % q = 0.0
   end if

   if (y % nlocal(metar) > 0) then
      allocate (y % metar(1:y % nlocal(metar)))
      y % metar(1:y % nlocal(metar)) % u = 0.0
      y % metar(1:y % nlocal(metar)) % v = 0.0
      y % metar(1:y % nlocal(metar)) % t = 0.0
      y % metar(1:y % nlocal(metar)) % p = 0.0
      y % metar(1:y % nlocal(metar)) % q = 0.0
   end if

   if (y % nlocal(geoamv) > 0) then
      allocate (y % geoamv(1:y % nlocal(geoamv)))
      do n = 1, y % nlocal(geoamv)
         nlevels = iv%info(geoamv)%levels(n)
         allocate (y % geoamv(n)%u(1:nlevels))
         allocate (y % geoamv(n)%v(1:nlevels))
         y % geoamv(n) % u(1:nlevels) = 0.0
         y % geoamv(n) % v(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(polaramv) > 0) then
      allocate (y % polaramv(1:y % nlocal(polaramv)))
      do n = 1, y % nlocal(polaramv)
         nlevels = iv%info(polaramv)%levels(n)
         allocate (y % polaramv(n)%u(1:nlevels))
         allocate (y % polaramv(n)%v(1:nlevels))
         y % polaramv(n) % u(1:nlevels) = 0.0
         y % polaramv(n) % v(1:nlevels) = 0.0
      end do
   end if
 
   if (y % nlocal(gpspw) > 0) then
      allocate (y % gpspw(1:y % nlocal(gpspw)))
      y % gpspw(1:y % nlocal(gpspw)) % tpw = 0.0
   end if

   if (y % nlocal(gpsref) > 0) then
      allocate (y % gpsref(1:y % nlocal(gpsref)))
      do n = 1, y % nlocal(gpsref)
         nlevels = iv%info(gpsref)%levels(n)
         allocate (y % gpsref(n)%ref(1:nlevels))
         allocate (y % gpsref(n)%  p(1:nlevels))
         allocate (y % gpsref(n)%  t(1:nlevels))
         allocate (y % gpsref(n)%  q(1:nlevels))

         y % gpsref(n) % ref(1:nlevels) = 0.0
         y % gpsref(n) %   p(1:nlevels) = 0.0
         y % gpsref(n) %   t(1:nlevels) = 0.0
         y % gpsref(n) %   q(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(sound) > 0) then
      allocate (y % sound(1:y % nlocal(sound)))
      do n = 1, y % nlocal(sound)
         nlevels = max(1,iv%info(sound)%levels(n))
         allocate (y % sound(n)%u(1:nlevels))
         allocate (y % sound(n)%v(1:nlevels))
         allocate (y % sound(n)%t(1:nlevels))
         allocate (y % sound(n)%q(1:nlevels))
         y % sound(n) % u(1:nlevels) = 0.0
         y % sound(n) % v(1:nlevels) = 0.0
         y % sound(n) % t(1:nlevels) = 0.0
         y % sound(n) % q(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(sonde_sfc) > 0) then
      allocate (y % sonde_sfc(1:y % nlocal(sonde_sfc)))

      y % sonde_sfc(1:y % nlocal(sonde_sfc)) % u = 0.0
      y % sonde_sfc(1:y % nlocal(sonde_sfc)) % v = 0.0
      y % sonde_sfc(1:y % nlocal(sonde_sfc)) % t = 0.0
      y % sonde_sfc(1:y % nlocal(sonde_sfc)) % p = 0.0
      y % sonde_sfc(1:y % nlocal(sonde_sfc)) % q = 0.0
   end if
     
   if (y % nlocal(mtgirs) > 0) then
      allocate (y % mtgirs(1:y % nlocal(mtgirs)))
      do n = 1, y % nlocal(mtgirs)
         nlevels = max(1,iv%info(mtgirs)%levels(n))
         allocate (y % mtgirs(n)%u(1:nlevels))
         allocate (y % mtgirs(n)%v(1:nlevels))
         allocate (y % mtgirs(n)%t(1:nlevels))
         allocate (y % mtgirs(n)%q(1:nlevels))
         y % mtgirs(n) % u(1:nlevels) = 0.0
         y % mtgirs(n) % v(1:nlevels) = 0.0
         y % mtgirs(n) % t(1:nlevels) = 0.0
         y % mtgirs(n) % q(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(tamdar) > 0) then
      allocate (y % tamdar(1:y % nlocal(tamdar)))
      do n = 1, y % nlocal(tamdar)
         nlevels = max(1,iv%info(tamdar)%levels(n))
         allocate (y % tamdar(n)%u(1:nlevels))
         allocate (y % tamdar(n)%v(1:nlevels))
         allocate (y % tamdar(n)%t(1:nlevels))
         allocate (y % tamdar(n)%q(1:nlevels))
         y % tamdar(n) % u(1:nlevels) = 0.0
         y % tamdar(n) % v(1:nlevels) = 0.0
         y % tamdar(n) % t(1:nlevels) = 0.0
         y % tamdar(n) % q(1:nlevels) = 0.0
      end do
   end if
   if (y % nlocal(tamdar_sfc) > 0) then
      allocate (y % tamdar_sfc(1:y % nlocal(tamdar_sfc)))

      y % tamdar_sfc(1:y % nlocal(tamdar_sfc)) % u = 0.0
      y % tamdar_sfc(1:y % nlocal(tamdar_sfc)) % v = 0.0
      y % tamdar_sfc(1:y % nlocal(tamdar_sfc)) % t = 0.0
      y % tamdar_sfc(1:y % nlocal(tamdar_sfc)) % p = 0.0
      y % tamdar_sfc(1:y % nlocal(tamdar_sfc)) % q = 0.
   end if

   if (y % nlocal(pilot) > 0) then
      allocate (y % pilot(1:y % nlocal(pilot)))
      do n = 1, y % nlocal(pilot)
         nlevels = iv%info(pilot)%levels(n)
         allocate (y % pilot(n)%u(1:nlevels))
         allocate (y % pilot(n)%v(1:nlevels))
         y % pilot(n) % u(1:nlevels) = 0.0
         y % pilot(n) % v(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(radar) > 0) then
      allocate (y % radar(1:y % nlocal(radar)))
      do n = 1, y % nlocal(radar)
         nlevels = iv%info(radar)%levels(n)
         allocate (y % radar(n)%rv(1:nlevels))
         allocate (y % radar(n)%rf(1:nlevels))
         allocate (y % radar(n)%rrn(1:nlevels))
         allocate (y % radar(n)%rsn(1:nlevels))
         allocate (y % radar(n)%rgr(1:nlevels))
         allocate (y % radar(n)%rcl(1:nlevels))
         allocate (y % radar(n)%rci(1:nlevels))
         allocate (y % radar(n)%rqv(1:nlevels))

         y % radar(n) % rv(1:nlevels)  = 0.0
         y % radar(n) % rf(1:nlevels)  = 0.0
         y % radar(n) % rrn(1:nlevels) = 0.0
         y % radar(n) % rsn(1:nlevels) = 0.0
         y % radar(n) % rgr(1:nlevels) = 0.0
         y % radar(n) % rcl(1:nlevels) = 0.0
         y % radar(n) % rci(1:nlevels) = 0.0
         y % radar(n) % rqv(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(airep) > 0) then
      allocate (y % airep(1:y % nlocal(airep)))
      do n = 1, y % nlocal(airep)
         nlevels = iv%info(airep)%levels(n)
         allocate (y % airep(n)%u(1:nlevels))
         allocate (y % airep(n)%v(1:nlevels))
         allocate (y % airep(n)%t(1:nlevels))
         allocate (y % airep(n)%q(1:nlevels))
         y % airep(n) % u(1:nlevels) = 0.0
         y % airep(n) % v(1:nlevels) = 0.0
         y % airep(n) % t(1:nlevels) = 0.0
         y % airep(n) % q(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(bogus) > 0) then
      allocate (y % bogus(1:y % nlocal(bogus)))
      do n = 1, y % nlocal(bogus)
         nlevels = iv%info(bogus)%levels(n)
         allocate (y % bogus(n)%u(1:nlevels))
         allocate (y % bogus(n)%v(1:nlevels))
         allocate (y % bogus(n)%t(1:nlevels))
         allocate (y % bogus(n)%q(1:nlevels))
         y % bogus(n) % u(1:nlevels) = 0.0
         y % bogus(n) % v(1:nlevels) = 0.0
         y % bogus(n) % t(1:nlevels) = 0.0
         y % bogus(n) % q(1:nlevels) = 0.0
      end do

      y % bogus(1:y % nlocal(bogus)) % slp = 0.0
   end if

   if (y % nlocal(satem) > 0) then
      allocate (y % satem(1:y % nlocal(satem)))
      do n = 1, y % nlocal(satem)
         nlevels = iv%info(satem)%levels(n)
         allocate (y % satem(n) % thickness(1:nlevels))
         y % satem(n) % thickness(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(ssmi_tb) > 0) then
      allocate (y % ssmi_tb(1:y % nlocal(ssmi_tb)))
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb19v = 0.0
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb19h = 0.0
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb22v = 0.0
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb37v = 0.0
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb37h = 0.0
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb85v = 0.0
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb85h = 0.0
   end if

   if (y % nlocal(ssmi_rv) > 0) then
        allocate (y % ssmi_rv(1:y % nlocal(ssmi_rv)))
        y % ssmi_rv(1:y % nlocal(ssmi_rv)) % tpw = 0.0
        y % ssmi_rv(1:y % nlocal(ssmi_rv)) % Speed = 0.0
   end if
   
   if (y % nlocal(ssmt1) > 0) then
      allocate (y % ssmt1(1:y % nlocal(ssmt1)))
      do n = 1, y % nlocal(ssmt1)
         nlevels = iv%info(ssmt1)%levels(n)
         allocate (y % ssmt1(n) % t(1:nlevels))
         y % ssmt1(n) % t(1:nlevels) = 0.0
      end do
   end if
   
   if (y % nlocal(ssmt2) > 0) then
      allocate (y % ssmt2(1:y % nlocal(ssmt2)))
      do n = 1, y % nlocal(ssmt2)
         nlevels=iv%info(ssmt2)%levels(n)
         allocate (y % ssmt2(n) % rh(1:nlevels))
         y % ssmt2(n) % rh(1:nlevels) = 0.0
      end do
   end if
   
   if (y % nlocal(pseudo) > 0) then
        allocate (y % pseudo(1:y % nlocal(pseudo)))
        y % pseudo(1:y % nlocal(pseudo)) % u = 0.0
        y % pseudo(1:y % nlocal(pseudo)) % v = 0.0
        y % pseudo(1:y % nlocal(pseudo)) % t = 0.0
        y % pseudo(1:y % nlocal(pseudo)) % p = 0.0
        y % pseudo(1:y % nlocal(pseudo)) % q = 0.0
   end if

   if (y % nlocal(qscat) > 0) then
      allocate (y % qscat(1:y % nlocal(qscat)))
      y % qscat(1:y % nlocal(qscat)) % u = 0.0
      y % qscat(1:y % nlocal(qscat)) % v = 0.0
   end if
      
   if (y % nlocal(profiler) > 0) then
      allocate (y % profiler(1:y % nlocal(profiler)))
      do n = 1, y % nlocal(profiler)
         nlevels = iv%info(profiler)%levels(n)
         allocate (y % profiler(n)%u(1:nlevels))
         allocate (y % profiler(n)%v(1:nlevels))
         y % profiler(n) % u(1:nlevels) = 0.0
         y % profiler(n) % v(1:nlevels) = 0.0
      end do
   end if

   if (y % nlocal(buoy) > 0) then
      allocate (y % buoy(1:y % nlocal(buoy)))
      y % buoy(1:y % nlocal(buoy)) % u = 0.0
      y % buoy(1:y % nlocal(buoy)) % v = 0.0
      y % buoy(1:y % nlocal(buoy)) % t = 0.0
      y % buoy(1:y % nlocal(buoy)) % p = 0.0
      y % buoy(1:y % nlocal(buoy)) % q = 0.0
   end if

   if (y % nlocal(rain) > 0) then
      allocate (y % rain(1:y % nlocal(rain)))
      y % rain(1:y % nlocal(rain)) % rain = 0.0
   end if

   if (y % num_inst > 0) then
      allocate (y % instid(1:y % num_inst))
      do i = 1,  y % num_inst
         y % instid(i) % num_rad = iv % instid(i) % num_rad
         y % instid(i) % nchan   = iv % instid(i) % nchan
         ! allocate (y % instid(i) % ichan(1:y % instid(i) % nchan))
         ! do n = 1, y % instid(i) % nchan
         !     y % instid(i) % ichan(n) = n
         ! end do
         if (y % instid(i) % num_rad < 1)  then
            nullify (y % instid(i) % tb)
            cycle
         end if
         allocate (y % instid(i) % tb(1:y % instid(i) % nchan, y % instid(i) % num_rad))
         y % instid(i) % tb(:,:) = 0.0
      end do
   end if

   if (y % nlocal(airsr) > 0) then
      allocate (y % airsr(1:y % nlocal(airsr)))
      do n = 1, y % nlocal(airsr)
         nlevels = iv%info(airsr)%levels(n)
         allocate (y % airsr(n)%t(1:nlevels))
         allocate (y % airsr(n)%q(1:nlevels))
         y % airsr(n) % t(1:nlevels) = 0.0
         y % airsr(n) % q(1:nlevels) = 0.0
      end do
   end if

   if (trace_use) call da_trace_exit("da_allocate_y")

end subroutine da_allocate_y


subroutine da_allocate_y_radar (iv, y)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate arrays used in y and residual obs structures.
   !---------------------------------------------------------------------------

   implicit none
   
   type (iv_type), intent(in)    :: iv      ! Ob type input.
   type (y_type),  intent(inout) :: y       ! Residual type structure.

   integer                       :: n       ! Loop counter.
   integer                       :: nlevels ! Number of levels.

   !---------------------------------------------------------------------------
   !  [1.0] Copy number of observations:
   !---------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_allocate_y_radar")

   y % nlocal(radar) = iv%info(radar)%nlocal
   y % ntotal(radar) = iv%info(radar)%ntotal

  !---------------------------------------------------------------------------
  ! [2.0] Allocate:
  !---------------------------------------------------------------------------

   if (y % nlocal(radar) > 0) then
      allocate (y % radar(1:y % nlocal(radar)))
      do n = 1, y % nlocal(radar)
         nlevels = iv%info(radar)%levels(n)
         allocate (y % radar(n)%rv(1:nlevels))
         allocate (y % radar(n)%rf(1:nlevels))
         allocate (y % radar(n)%rrn(1:nlevels))
         allocate (y % radar(n)%rsn(1:nlevels))
         allocate (y % radar(n)%rgr(1:nlevels))
         allocate (y % radar(n)%rcl(1:nlevels))
         allocate (y % radar(n)%rci(1:nlevels))
         allocate (y % radar(n)%rqv(1:nlevels))

         y % radar(n) % rv(1:nlevels)  = 0.0
         y % radar(n) % rf(1:nlevels)  = 0.0
         y % radar(n) % rrn(1:nlevels) = 0.0
         y % radar(n) % rsn(1:nlevels) = 0.0
         y % radar(n) % rgr(1:nlevels) = 0.0
         y % radar(n) % rcl(1:nlevels) = 0.0
         y % radar(n) % rci(1:nlevels) = 0.0
         y % radar(n) % rqv(1:nlevels) = 0.0
      end do
   end if

   if (trace_use) call da_trace_exit("da_allocate_y_radar")

end subroutine da_allocate_y_radar


subroutine da_allocate_y_rain (iv, y)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate arrays used in y and residual obs structures.
   !---------------------------------------------------------------------------

   implicit none
   
   type (iv_type), intent(in)    :: iv      ! Ob type input.
   type (y_type),  intent(inout) :: y       ! Residual type structure.

   !---------------------------------------------------------------------------
   !  [1.0] Copy number of observations:
   !---------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_allocate_y_rain")

   y % nlocal(rain) = iv%info(rain)%nlocal
   y % ntotal(rain) = iv%info(rain)%ntotal

  !---------------------------------------------------------------------------
  ! [2.0] Allocate:
  !---------------------------------------------------------------------------

   if (y % nlocal(rain) > 0) then
      allocate (y % rain(1:y % nlocal(rain)))
      y % rain(1:y % nlocal(rain)) % rain = 0.0
   end if

   if (trace_use) call da_trace_exit("da_allocate_y_rain")

end subroutine da_allocate_y_rain


subroutine da_deallocate_background_errors (be)

   !---------------------------------------------------------------------------
   ! Purpose: Deallocate components of wrfvar background errors.
   !
   !  Update: Multivariate BE option (cv_options=6)
   !          Syed RH Rizvi (MMM/NESL/NCAR)   Date: 02/01/2010
   !
   !  Note: Please acknowledge author/institute in work that uses this code.
   !---------------------------------------------------------------------------

   implicit none

   type (be_type), intent(inout)        :: be     ! Background error structure.
   
   if (trace_use) call da_trace_entry("da_deallocate_background_errors")

   if (cv_options /= 3) then

      ! Deallocate gridpoint errors:

      if (be % v1 % mz > 0) deallocate (be % v1 % val)
      if (be % v2 % mz > 0) deallocate (be % v2 % val)
      if (be % v3 % mz > 0) deallocate (be % v3 % val)
      if (be % v4 % mz > 0) deallocate (be % v4 % val)
      if (be % v5 % mz > 0 .and. .not. global) deallocate (be % v5 % val) 
      if (be % v1 % mz > 0) deallocate (be % v1 % rf_alpha)
      if (be % v2 % mz > 0) deallocate (be % v2 % rf_alpha)
      if (be % v3 % mz > 0) deallocate (be % v3 % rf_alpha)
      if (be % v4 % mz > 0) deallocate (be % v4 % rf_alpha)
      if (be % v5 % mz > 0 .and. .not. global) deallocate (be % v5 % rf_alpha)
      if (global) then
         if (be % v1 % mz > 0) deallocate (be % v1 % power)
         if (be % v2 % mz > 0) deallocate (be % v2 % power)
         if (be % v3 % mz > 0) deallocate (be % v3 % power)
         if (be % v4 % mz > 0) deallocate (be % v4 % power)
         if (be % v5 % mz > 0) deallocate (be % v5 % power) 
      end if

      ! Deallocate eigenvectors of vertical error covariance:

      if (vert_corr == vert_corr_2) then
         if (be % v1 % mz > 0) deallocate (be % v1 % evec)
         if (be % v2 % mz > 0) deallocate (be % v2 % evec)
         if (be % v3 % mz > 0) deallocate (be % v3 % evec)
         if (be % v4 % mz > 0) deallocate (be % v4 % evec)
         if (be % v5 % mz > 0 .and. .not. global) deallocate (be % v5 % evec)
         if (be % v1 % mz > 0) deallocate (be % v1 % evec_g)
         if (be % v2 % mz > 0) deallocate (be % v2 % evec_g)
         if (be % v3 % mz > 0) deallocate (be % v3 % evec_g)
         if (be % v4 % mz > 0) deallocate (be % v4 % evec_g)
         if (be % v5 % mz > 0 .and. .not. global) deallocate (be % v5 % evec_g)
         if (be % v1 % mz > 0) deallocate (be % v1 % val_g)
         if (be % v2 % mz > 0) deallocate (be % v2 % val_g)
         if (be % v3 % mz > 0) deallocate (be % v3 % val_g)
         if (be % v4 % mz > 0) deallocate (be % v4 % val_g)
         if (be % v5 % mz > 0 .and. .not. global) deallocate (be % v5 % val_g)
      end if

      if ( cv_options /= 7 ) then
         deallocate (be % reg_psi_chi)
         deallocate (be % reg_psi_t)
         deallocate (be % reg_psi_ps)
         if ( cv_options == 6 ) then
            deallocate (be % reg_psi_rh)
            deallocate (be % reg_chi_u_t)
            deallocate (be % reg_chi_u_ps)
            deallocate (be % reg_chi_u_rh)
            deallocate (be % reg_t_u_rh)
            deallocate (be % reg_ps_u_rh)
         end if
      end if

      ! Deallocate control variable errors (in future uncomment use these to allow 
      ! eg NMC error correlations).

      ! deallocate (be % cv % val)

   else ! for cv_options = 3
    
      deallocate (be % corz)
      deallocate (be % corp)
      deallocate (be % vz)
      deallocate (be % agvz)
      deallocate (be % bvz)
      deallocate (be % wgvz)
      deallocate (be % be)
      deallocate (be % rate)
      deallocate (be % table)
      deallocate (be % slix)
      deallocate (be % slipx)
      deallocate (be % sljy)
      deallocate (be % sljpy)
   
   end if

   ! Deallocate wavelet parameters:
   if( .not. use_rf )deallocate(be%wsd,ws)
   if( do_normalize )deallocate(be%sd)
   if( do_normalize .or. .not. use_rf )deallocate(nij)

   if (trace_use) call da_trace_exit("da_deallocate_background_errors")

end subroutine da_deallocate_background_errors


subroutine da_deallocate_observations (iv)

   !---------------------------------------------------------------------------
   ! Purpose: Deallocate components of observation structure.
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout)        :: iv     ! Observation structure.
   integer   :: n

   if (trace_use) call da_trace_entry("da_deallocate_observations")

   !---------------------------------------------------------------------------
   ! [1.0] Deallocate:
   !---------------------------------------------------------------------------

   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%levels(n) > 0) then
            deallocate (iv%sound (n) % h)
            deallocate (iv%sound (n) % p)
            deallocate (iv%sound (n) % u)
            deallocate (iv%sound (n) % v)
            deallocate (iv%sound (n) % t)
            deallocate (iv%sound (n) % q)
         end if
      end do

      deallocate (iv%sound)
   end if

   if (iv%info(sonde_sfc)%nlocal > 0) deallocate (iv%sonde_sfc)
      
   if (iv%info(mtgirs)%nlocal > 0) then
      do n = 1, iv%info(mtgirs)%nlocal
         if (iv%info(mtgirs)%levels(n) > 0) then
            deallocate (iv%mtgirs (n) % h)
            deallocate (iv%mtgirs (n) % p)
            deallocate (iv%mtgirs (n) % u)
            deallocate (iv%mtgirs (n) % v)
            deallocate (iv%mtgirs (n) % t)
            deallocate (iv%mtgirs (n) % q)
         end if
      end do

      deallocate (iv%mtgirs)

   end if

   if (iv%info(tamdar)%nlocal > 0) then
      do n = 1, iv%info(tamdar)%nlocal
         if (iv%info(tamdar)%levels(n) > 0) then
            deallocate (iv%tamdar (n) % h)
            deallocate (iv%tamdar (n) % p)
            deallocate (iv%tamdar (n) % u)
            deallocate (iv%tamdar (n) % v)
            deallocate (iv%tamdar (n) % t)
            deallocate (iv%tamdar (n) % q)
         end if
      end do

      deallocate (iv%tamdar)
   end if
   if (iv%info(tamdar_sfc)%nlocal > 0) deallocate (iv%tamdar_sfc)
   if (iv%info(synop)%nlocal > 0) deallocate (iv%synop)

   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         deallocate (iv%airep (n) % h)
         deallocate (iv%airep (n) % p)
         deallocate (iv%airep (n) % u)
         deallocate (iv%airep (n) % v)
         deallocate (iv%airep (n) % t)
         deallocate (iv%airep (n) % q)
      end do

      deallocate (iv%airep)
   end if

   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         deallocate (iv%satem(n) % p)
         deallocate (iv%satem(n) % thickness)
         deallocate (iv%satem(n) % org_thickness)
      end do
      deallocate (iv%satem)
   end if

   if (iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
         deallocate (iv%geoamv(n) % p)
         deallocate (iv%geoamv(n) % u)
         deallocate (iv%geoamv(n) % v)
      end do
      deallocate (iv%geoamv)
   end if


   if (iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
         deallocate (iv%polaramv(n) % p)
         deallocate (iv%polaramv(n) % u)
         deallocate (iv%polaramv(n) % v)
      end do
      deallocate (iv%polaramv)
   end if

   if (iv%info(metar)%nlocal > 0) deallocate (iv%metar)
   if (iv%info(ships)%nlocal > 0) deallocate (iv%ships)

   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
         deallocate (iv%pilot (n) % h)
         deallocate (iv%pilot (n) % p)
         deallocate (iv%pilot (n) % u)
         deallocate (iv%pilot (n) % v)
      end do

      deallocate (iv%pilot)
   end if

   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         deallocate (iv%bogus (n) % h)
         deallocate (iv%bogus (n) % p)
         deallocate (iv%bogus (n) % u)
         deallocate (iv%bogus (n) % v)
         deallocate (iv%bogus (n) % t)
         deallocate (iv%bogus (n) % q)
      end do

      deallocate (iv%bogus)
   end if

   if (iv%info(radar)%nlocal > 0) then
      do n = 1, iv%info(radar)%nlocal
         deallocate (iv%radar (n) % model_p)
         deallocate (iv%radar (n) % model_rho)
         deallocate (iv%radar (n) % model_qrn)
         deallocate (iv%radar (n) % height  )
         deallocate (iv%radar (n) % height_qc)
         deallocate (iv%radar (n) % rv      )
         deallocate (iv%radar (n) % rf      )
      end do

      deallocate (iv%radar)
   end if

   if (iv%info(rain)%nlocal > 0) deallocate (iv%rain)
   if (iv%info(gpspw)%nlocal > 0) deallocate (iv%gpspw)

   if (iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
         deallocate (iv%gpsref(n) %  h)
         deallocate (iv%gpsref(n) % ref)
         deallocate (iv%gpsref(n) %   p)
         deallocate (iv%gpsref(n) %   t)
         deallocate (iv%gpsref(n) %   q)
      end do
      deallocate (iv%gpsref)
   end if

   if (iv%info(ssmi_tb)%nlocal > 0) deallocate (iv%ssmi_tb)
   if (iv%info(ssmi_rv)%nlocal > 0) deallocate (iv%ssmi_rv)

   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
         deallocate (iv%ssmt1(n) % h)
         deallocate (iv%ssmt1(n) % p)
         deallocate (iv%ssmt1(n) % t)
      end do
   
      deallocate (iv%ssmt1)
   end if
   
   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         deallocate (iv%ssmt2(n) % h)
         deallocate (iv%ssmt2(n) % p)
         deallocate (iv%ssmt2(n) % rh)
      end do
   
      deallocate (iv%ssmt2)
   end if

   if (iv%info(qscat)%nlocal > 0) deallocate (iv%qscat)

   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         deallocate (iv%profiler(n)%h)
         deallocate (iv%profiler(n)%p)
         deallocate (iv%profiler(n)%u)
         deallocate (iv%profiler(n)%v)
      end do

      deallocate(iv%profiler)
   end if

   if (iv%info(buoy)%nlocal     > 0) deallocate(iv%buoy)

   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         deallocate (iv%airsr (n) % h)
         deallocate (iv%airsr (n) % p)
         deallocate (iv%airsr (n) % t)
         deallocate (iv%airsr (n) % q)
      end do

      deallocate (iv%airsr)
   end if

   do n = 1, num_ob_indexes
      if (n .ne. radiance .and. iv%info(n)%nlocal > 0) then
         deallocate (iv%info(n)%name)
         deallocate (iv%info(n)%platform)
         deallocate (iv%info(n)%id)
         deallocate (iv%info(n)%date_char)
         deallocate (iv%info(n)%levels)
         deallocate (iv%info(n)%lat)
         deallocate (iv%info(n)%lon)
         deallocate (iv%info(n)%elv)
         deallocate (iv%info(n)%pstar)
         deallocate (iv%info(n)%slp)
         deallocate (iv%info(n)%pw)
         deallocate (iv%info(n)%x)
         deallocate (iv%info(n)%y)
         deallocate (iv%info(n)%i)
         deallocate (iv%info(n)%j)
         deallocate (iv%info(n)%dx)
         deallocate (iv%info(n)%dxm)
         deallocate (iv%info(n)%dy)
         deallocate (iv%info(n)%dym)
         deallocate (iv%info(n)%k)
         deallocate (iv%info(n)%dz)
         deallocate (iv%info(n)%dzm)
         deallocate (iv%info(n)%zk)
         deallocate (iv%info(n)%proc_domain)
         deallocate (iv%info(n)%thinned)
         deallocate (iv%info(n)%obs_global_index)
      end if
   end do

   if (trace_use) call da_trace_exit("da_deallocate_observations")

end subroutine da_deallocate_observations


subroutine da_deallocate_y(y)

   !---------------------------------------------------------------------------
   ! Purpose: Deallocate arrays used in y and residual obs structures.
   !
   ! Method:  Deallocate component in turn.
   !---------------------------------------------------------------------------

   implicit none
   
   type (y_type), intent(inout)          :: y      ! residual type structure.
   integer                               :: n,i  ! Loop counter.


   if (trace_use) call da_trace_entry("da_deallocate_y")

   !---------------------------------------------------------------------------
   ! [1.0] Deallocate:
   !---------------------------------------------------------------------------

   if (y % nlocal(synop) > 0) deallocate (y % synop)

   if (y % nlocal(ships) > 0) deallocate (y % ships)

   if (y % nlocal(metar) > 0) deallocate (y % metar)

   if (y % nlocal(rain) > 0)  deallocate (y % rain)

   if (y % nlocal(sound) > 0) then
      do n = 1, y % nlocal(sound)
         deallocate (y % sound(n)%u)
         deallocate (y % sound(n)%v)
         deallocate (y % sound(n)%t)
         deallocate (y % sound(n)%q)
      end do
      deallocate (y % sound)
   end if

   if (y % nlocal(sonde_sfc) > 0) deallocate (y % sonde_sfc)
      
   if (y % nlocal(mtgirs) > 0) then
      do n = 1, y % nlocal(mtgirs)
         deallocate (y % mtgirs(n)%u)
         deallocate (y % mtgirs(n)%v)
         deallocate (y % mtgirs(n)%t)
         deallocate (y % mtgirs(n)%q)
      end do

      deallocate (y % mtgirs)

   end if

   if (y % nlocal(tamdar) > 0) then
      do n = 1, y % nlocal(tamdar)
         deallocate (y % tamdar(n)%u)
         deallocate (y % tamdar(n)%v)
         deallocate (y % tamdar(n)%t)
         deallocate (y % tamdar(n)%q)
      end do

      deallocate (y % tamdar)
   end if

   if (y % nlocal(tamdar_sfc) > 0) deallocate (y % tamdar_sfc)

   if (y % nlocal(pilot) > 0) then
      do n = 1, y % nlocal(pilot)
         deallocate (y % pilot(n)%u)
         deallocate (y % pilot(n)%v)
      end do
      deallocate (y % pilot)
   end if

   if (y % nlocal(bogus) > 0) then
      do n = 1, y % nlocal(bogus)
         deallocate (y % bogus(n)%u)
         deallocate (y % bogus(n)%v)
         deallocate (y % bogus(n)%t)
         deallocate (y % bogus(n)%q)
      end do
      deallocate (y % bogus)
   end if

    if (y % nlocal(radar) > 0) then
       do n = 1, y % nlocal(radar)
          deallocate (y % radar(n)%rv)
          deallocate (y % radar(n)%rf)
       end do
       deallocate (y % radar)
    end if


   if (y % nlocal(airep) > 0) then
      do n = 1, y % nlocal(airep)
         deallocate (y % airep(n)%u)
         deallocate (y % airep(n)%v)
         deallocate (y % airep(n)%t)
         deallocate (y % airep(n)%q)
      end do
      deallocate (y % airep)
   end if

   if (y % nlocal(geoamv) > 0) then
      do n=1, y % nlocal(geoamv)
         deallocate (y % geoamv(n) % u)
         deallocate (y % geoamv(n) % v)
      end do
      deallocate (y % geoamv)
   end if

   if (y % nlocal(polaramv) > 0) then
      do n=1, y % nlocal(polaramv)
         deallocate (y % polaramv(n) % u)
         deallocate (y % polaramv(n) % v)
      end do
      deallocate (y % polaramv)
   end if

   if (y % nlocal(gpspw) > 0) deallocate (y % gpspw)

   if (y % nlocal(gpsref) > 0) then
      do n = 1, y % nlocal(gpsref)
         deallocate (y % gpsref(n)%ref)
         deallocate (y % gpsref(n)%  p)
         deallocate (y % gpsref(n)%  t)
         deallocate (y % gpsref(n)%  q)
      end do
      deallocate (y % gpsref)
   end if

   if (y % nlocal(satem) > 0) then
      do n = 1, y % nlocal(satem)
         deallocate (y % satem(n) % thickness)
      end do
      deallocate (y % satem)
   end if

   if (y % nlocal(ssmi_tb) > 0) deallocate (y % ssmi_tb)
   if (y % nlocal(ssmi_rv) > 0) deallocate (y % ssmi_rv)
   if (y % nlocal(pseudo)  > 0) deallocate (y % pseudo)

   if (y % nlocal(ssmt1) > 0) then
      do n = 1, y % nlocal(ssmt1)
         deallocate (y % ssmt1(n) % t)
      end do
      deallocate (y % ssmt1)
   end if

   if (y % nlocal(ssmt2) > 0) then
      do n = 1, y % nlocal(ssmt2)
         deallocate (y % ssmt2(n) % rh)
      end do
      deallocate (y % ssmt2)
   end if

   if (y % nlocal(qscat) > 0) deallocate (y % qscat)

   if (y % nlocal(profiler) > 0) then
      do n = 1, y % nlocal(profiler)
         deallocate (y % profiler(n)%u)
         deallocate (y % profiler(n)%v)
      end do
      deallocate (y % profiler)
   end if

   if (y % nlocal(buoy)  > 0) deallocate (y % buoy)

   !  radiance:
 
   if (y % num_inst > 0) then
      do i = 1,  y % num_inst
        if (y % instid(i) % num_rad < 1) cycle
        ! deallocate (y % instid(i) % ichan)
        deallocate ( y % instid(i) % tb )
      end do
      deallocate (y % instid)
   end if
   if (y % nlocal(airsr) > 0) then
      do n = 1, y % nlocal(airsr)
         deallocate (y % airsr(n)%t)
         deallocate (y % airsr(n)%q)
      end do
      deallocate (y % airsr)
   end if

   if (trace_use) call da_trace_exit("da_deallocate_y")

end subroutine da_deallocate_y


subroutine da_zero_x ( x )

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (x_type), intent(inout)        :: x      ! Analysis incrs structure.

   if (trace_use_dull) call da_trace_entry("da_zero_x")

   x % u(:,:,:) = 0.0
   x % v(:,:,:) = 0.0
   x % w(:,:,:) = 0.0
   x % t(:,:,:) = 0.0
   x % q(:,:,:) = 0.0
   x % p(:,:,:) = 0.0
   x % geoh(:,:,:) = 0.0
   x % rh(:,:,:) = 0.0
   x % wh(:,:,:) = 0.0
   x % rho(:,:,:) = 0.0
   x % ref(:,:,:) = 0.0

   x % qcw(:,:,:) = 0.0
   x % qrn(:,:,:) = 0.0
   x % qt (:,:,:) = 0.0
   x % qci(:,:,:) = 0.0
   x % qsn(:,:,:) = 0.0
   x % qgr(:,:,:) = 0.0

   x % tgrn(:,:) = 0.0
   x % psfc(:,:) = 0.0
   x % mu(:,:) = 0.0
   x % u10(:,:) = 0.0
   x % v10(:,:) = 0.0
   x % t2(:,:) = 0.0
   x % q2(:,:) = 0.0

   x % ztd(:,:) = 0.0
   x % tpw(:,:) = 0.0
   x % speed(:,:) = 0.0
   x % tb19v(:,:) = 0.0
   x % tb19h(:,:) = 0.0
   x % tb22v(:,:) = 0.0
   x % tb37v(:,:) = 0.0
   x % tb37h(:,:) = 0.0
   x % tb85v(:,:) = 0.0
   x % tb85h(:,:) = 0.0

   if (trace_use_dull) call da_trace_exit("da_zero_x")

end subroutine da_zero_x


subroutine da_zero_y( iv, y, value )

   !---------------------------------------------------------------------------
   ! Purpose: Initialises the Y-array
   !---------------------------------------------------------------------------

   implicit none
   
   type (iv_type), intent(in)            :: iv      ! Ob type input.
   type (y_type),  intent(inout)         :: y       ! Residual type structure.
   real, optional, intent(inout)         :: value

   integer                               :: n       ! Loop counter.
   integer                               :: nlevels ! Number of levels.

   if (trace_use_dull) call da_trace_entry("da_zero_y")
 
   if (.not.(present(value))) value = 0.0
   !---------------------------------------------------------------------------
   ! [1.0] Copy number of observations:
   !---------------------------------------------------------------------------

   y % nlocal(:) = iv % info(:) % nlocal

   !---------------------------------------------------------------------------
   ! [2.0] Allocate:
   !---------------------------------------------------------------------------

   ! Initialize synops:

   if ( y % nlocal(synop) > 0 ) then
      y % synop(1:y % nlocal(synop)) % u = value
      y % synop(1:y % nlocal(synop)) % v = value
      y % synop(1:y % nlocal(synop)) % t = value
      y % synop(1:y % nlocal(synop)) % p = value
      y % synop(1:y % nlocal(synop)) % q = value
   end if

   ! Initialize ships:

   if ( y % nlocal(ships) > 0 ) then
      y % ships(1:y % nlocal(ships)) % u = value
      y % ships(1:y % nlocal(ships)) % v = value
      y % ships(1:y % nlocal(ships)) % t = value
      y % ships(1:y % nlocal(ships)) % p = value
      y % ships(1:y % nlocal(ships)) % q = value
   end if

   ! Initialize metars:

   if ( y % nlocal(metar) > 0 ) then
      y % metar(1:y % nlocal(metar)) % u = value
      y % metar(1:y % nlocal(metar)) % v = value
      y % metar(1:y % nlocal(metar)) % t = value
      y % metar(1:y % nlocal(metar)) % p = value
      y % metar(1:y % nlocal(metar)) % q = value
   end if

   ! Initialize Geo. AMV's:

   if ( y % nlocal(geoamv) > 0 ) then
      do n = 1, y % nlocal(geoamv)
       nlevels = iv%info(geoamv)%levels(n)
       y % geoamv(n) % u(1:nlevels) = value
       y % geoamv(n) % v(1:nlevels) = value
      end do
   end if

   ! Initialize Polat AMVs:

   if ( y % nlocal(polaramv) > 0 ) then
      do n = 1, y % nlocal(polaramv)
       nlevels = iv%info(polaramv)%levels(n)
       y % polaramv(n) % u(1:nlevels) = value
       y % polaramv(n) % v(1:nlevels) = value
      end do
   end if

   ! Initialize GPS TPW:

   if ( y % nlocal(gpspw) > 0 ) then
      y % gpspw(1:y % nlocal(gpspw)) % tpw = value
   end if

   ! Initialize GPS REFRACTIVITY:

   if ( y % nlocal(gpsref) > 0 ) then
      do n = 1, y % nlocal(gpsref)
         nlevels = iv % info(gpsref) % levels(n)
         y % gpsref(n) % ref(1:nlevels) = value
         y % gpsref(n) %   p(1:nlevels) = value
         y % gpsref(n) %   t(1:nlevels) = value
         y % gpsref(n) %   q(1:nlevels) = value
      end do
   end if

   ! Initialize sondes:

   if ( y % nlocal(sound) > 0 ) then
      do n = 1, y % nlocal(sound)
         nlevels = iv%info(sound)%levels(n)

         y % sound(n) % u(1:nlevels) = value
         y % sound(n) % v(1:nlevels) = value
         y % sound(n) % t(1:nlevels) = value
         y % sound(n) % q(1:nlevels) = value
      end do
   end if

   ! Initialize sonde_sfc
   if ( y % nlocal(sonde_sfc) > 0 ) then
      do n = 1, y % nlocal(sonde_sfc)
         y % sonde_sfc(n) % u = value
         y % sonde_sfc(n) % v = value
         y % sonde_sfc(n) % t = value
         y % sonde_sfc(n) % p = value
         y % sonde_sfc(n) % q = value
      end do
   end if
      
   if ( y % nlocal(mtgirs) > 0 ) then
      do n = 1, y % nlocal(mtgirs)
         nlevels = iv%info(mtgirs)%levels(n)

         y % mtgirs(n) % u(1:nlevels) = value
         y % mtgirs(n) % v(1:nlevels) = value
         y % mtgirs(n) % t(1:nlevels) = value
         y % mtgirs(n) % q(1:nlevels) = value

      end do
   end if

   if ( y % nlocal(tamdar) > 0 ) then
      do n = 1, y % nlocal(tamdar)
         nlevels = iv%info(tamdar)%levels(n)

         y % tamdar(n) % u(1:nlevels) = 0.0
         y % tamdar(n) % v(1:nlevels) = 0.0
         y % tamdar(n) % t(1:nlevels) = 0.0
         y % tamdar(n) % q(1:nlevels) = 0.0

      end do
   end if

! Initialize tamdar_sfc
   if ( y % nlocal(tamdar_sfc) > 0 ) then

         y % tamdar_sfc(n) % u = 0.0
         y % tamdar_sfc(n) % v = 0.0
         y % tamdar_sfc(n) % t = 0.0
         y % tamdar_sfc(n) % p = 0.0
         y % tamdar_sfc(n) % q = 0.0
   end if

   if ( y % nlocal(bogus) > 0 ) then
      do n = 1, y % nlocal(bogus)
         nlevels = iv % info(bogus) % levels(n)

         y % bogus(n) % u(1:nlevels) = value
         y % bogus(n) % v(1:nlevels) = value
         y % bogus(n) % t(1:nlevels) = value
         y % bogus(n) % q(1:nlevels) = value
         y % bogus(n) % slp          = value
      end do
   end if

   ! Initialize pilots:

   if ( y % nlocal(pilot) > 0 ) then
      do n = 1, y % nlocal(pilot)
         nlevels = iv % info(pilot) % levels(n)

         y % pilot(n) % u(1:nlevels) = value
         y % pilot(n) % v(1:nlevels) = value
      end do
   end if

   ! Initialize AIREPs:

   if ( y % nlocal(airep) > 0 ) then
      do n = 1, y % nlocal(airep)
         nlevels = iv%info(airep)%levels(n)

         y % airep(n) % u(1:nlevels) = value
         y % airep(n) % v(1:nlevels) = value
         y % airep(n) % t(1:nlevels) = value
         y % airep(n) % q(1:nlevels) = value
      end do
   end if

   ! Initialize satem:

   if ( y % nlocal(satem) > 0 ) then
      do n = 1, y % nlocal(satem)
         nlevels = iv % info(satem) % levels(n)

         y % satem(n) % thickness(1:nlevels) = value
      end do
   end if

   if ( y % nlocal(ssmi_tb) > 0 ) then
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb19v = value
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb19h = value
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb22v = value
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb37v = value
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb37h = value
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb85v = value
      y % ssmi_tb(1:y % nlocal(ssmi_tb)) % tb85h = value
   end if

   if ( y % nlocal(ssmi_rv) > 0 ) then
        y % ssmi_rv(1:y % nlocal(ssmi_rv)) % tpw = value
        y % ssmi_rv(1:y % nlocal(ssmi_rv)) % Speed = value
   end if
   
   if ( y % nlocal(ssmt1) > 0 ) then
      do n = 1, y % nlocal(ssmt1)
         nlevels = iv % info(ssmt1) % levels(n)
         y % ssmt1(n) % t(1:nlevels) = value
      end do
   end if
   
   if ( y % nlocal(ssmt2) > 0 ) then
      do n = 1, y % nlocal(ssmt2)
         nlevels = iv % info(ssmt2) % levels(n)
         y % ssmt2(n) % rh(1:nlevels) = value
      end do
   end if
   
   if ( num_pseudo > 0 ) then
        y % pseudo(1:num_pseudo) % u = value
        y % pseudo(1:num_pseudo) % v = value
        y % pseudo(1:num_pseudo) % t = value
        y % pseudo(1:num_pseudo) % p = value
        y % pseudo(1:num_pseudo) % q = value
   end if

   !  Initialize Quikscat:

   if ( y % nlocal(qscat) > 0 ) then
      y % qscat(1:y % nlocal(qscat)) % u = value
      y % qscat(1:y % nlocal(qscat)) % v = value
   end if
      
   ! Initialize profilers:

   if ( y % nlocal(profiler) > 0 ) then
      do n = 1, y % nlocal(profiler)
         nlevels = iv % info(profiler) % levels(n)

         y % profiler(n) % u(1:nlevels) = value
         y % profiler(n) % v(1:nlevels) = value
      end do
   end if

   ! Initialize buoy:

   if ( y % nlocal(buoy) > 0 ) then
      y % buoy(1:y % nlocal(buoy)) % u = value
      y % buoy(1:y % nlocal(buoy)) % v = value
      y % buoy(1:y % nlocal(buoy)) % t = value
      y % buoy(1:y % nlocal(buoy)) % p = value
      y % buoy(1:y % nlocal(buoy)) % q = value
   end if

   ! Initialize radar:
   if ( y % nlocal(radar) > 0 ) then
      do n = 1, y % nlocal(radar)
         nlevels = iv % info(radar) % levels(n)

         y % radar(n) % rv(1:nlevels) = value
         y % radar(n) % rf(1:nlevels) = value
      end do
   end if

   ! Initialize rain:
   if ( y % nlocal(rain) > 0 ) then
          y % rain(1:y % nlocal(rain)) % rain = value
   end if

   ! Initialize AIRS retrievals:

   if ( y % nlocal(airsr) > 0 ) then
      do n = 1, y % nlocal(airsr)
         nlevels = iv % info(airsr) % levels(n)

         y % airsr(n) % t(1:nlevels) = value
         y % airsr(n) % q(1:nlevels) = value
      end do
   end if


   if (trace_use_dull) call da_trace_exit("da_zero_y")

end subroutine da_zero_y        


subroutine da_zero_vp_type( vp )

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (vp_type), intent(inout) :: vp

   if (trace_use_dull) call da_trace_entry("da_zero_vp_type")
 
   !  Standard fields:
   if (associated(vp % v1)) vp % v1(:,:,:) = 0.0
   if (associated(vp % v2)) vp % v2(:,:,:) = 0.0
   if (associated(vp % v3)) vp % v3(:,:,:) = 0.0
   if (associated(vp % v4)) vp % v4(:,:,:) = 0.0
   if (associated(vp % v5)) vp % v5(:,:,:) = 0.0
   ! Flow-dependent control variables:
   if (associated(vp % alpha) ) vp % alpha(:,:,:,:) = 0.0

   if (trace_use_dull) call da_trace_exit("da_zero_vp_type")

end subroutine da_zero_vp_type


subroutine da_initialize_cv(cv_size, cv)

   !---------------------------------------------------------------------------
   ! Purpose: Initialize components of control variable.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)   :: cv_size
   real,    intent(out)  :: cv(1:cv_size)    ! Control variable structure.

   integer                              :: i
   real                                 :: z, mean_cv, rms_cv, std_dev_cv

   if (trace_use) call da_trace_entry("da_initialize_cv")

   !---------------------------------------------------------------------------
   ! [1.0] Initialize cv:
   !---------------------------------------------------------------------------

   if (anal_type_randomcv) then
   
      ! [2.1] Initialize random number generator and scalars:

      call da_random_seed
      
      if( use_rf )then

         ! [2.2] Calculate random numbers with Gaussian distribution:

         do i = 1, cv_size
            call da_gauss_noise(z)
            cv(i) = z
         end do
      else
         write(unit=message(1),fmt='(a)')'Need to inject CV into wavelet space'
         call da_error("da_initialize_cv.inc",37,message(1:1))
      endif

      mean_cv = sum(cv) / real(cv_size)
      rms_cv = sqrt(sum(cv*cv) / real(cv_size))
      std_dev_cv = sqrt(rms_cv * rms_cv - mean_cv * mean_cv)

      write(unit=message(1),fmt='(a)')' Gaussian (Normal) noise statistics:'
      write(unit=message(2),fmt='(a,f15.5)')' Mean = ',mean_cv
      write(unit=message(3),fmt='(a,f15.5)')' RMS = ', rms_cv
      write(unit=message(4),fmt='(a,f15.5)')' STD DEV = ', std_dev_cv
      call da_message(message(1:4))
   else
      cv = 0.0
   end if

   if (trace_use) call da_trace_exit("da_initialize_cv")

end subroutine da_initialize_cv


subroutine da_random_seed

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   INCLUDE 'mpif.h'

   integer              :: seed_size
   integer, allocatable :: seed_array(:)

   integer              :: myproc,ierr,i

   if (trace_use) call da_trace_entry("da_random_seed")

   !----------------------------------------------------------------------------
   !  Check that right seed_size is being used:
   !----------------------------------------------------------------------------

   myproc=0
   call wrf_get_dm_communicator (comm)
   call mpi_comm_rank (comm, myproc, ierr)

   call random_seed(size=seed_size)              ! Get size of seed array.
   allocate(seed_array(1:seed_size))
   seed_array(1:seed_size) = 1

   if (put_rand_seed) then            ! Manually set random seed.

      if ( (seed_array1 == 0) .or. (seed_array2 == 0) ) then
         write(unit=message(1),fmt='(a)') ' Error: can not use "0" as a random seed!'
         write(unit=message(2),fmt='(a,i16)') ' seed_array1 = ',seed_array1
         write(unit=message(3),fmt='(a,i16)') ' seed_array2 = ',seed_array2
         call da_error("da_random_seed.inc",40,message(1:3))
      end if

      if (seed_size == 1) then
         write(unit=message(1),fmt='(a)') &
            ' Warning: this compiler only supports a single random seed; only using seed_array1!'
         call da_warning("da_random_seed.inc",46,message(1:1))
         seed_array(1) = seed_array1
         write(unit=message(1),fmt='(a,i16)')' Setting seed_array(1) = ', seed_array(1)
      else if (seed_size > 2) then
         write(unit=message(1),fmt='(a,i2,a)') &
            ' Note: this compiler expects an array of ',seed_size,' integers to the "random_seed" function; '
         write(unit=message(2),fmt='(a)') &
            ' filling the rest of the array with copies of seed_array1 and seed_array2'
         call da_warning("da_random_seed.inc",54,message(1:2))
         do i = 1,seed_size
            if ( mod (i,2) == 1 ) then
               seed_array(i) = seed_array1
            else
               seed_array(i) = seed_array2 * seed_array1 + myproc*10000000
            end if
            write(unit=message(1),fmt='(a,i0,a,i16)')' Setting seed_array(',i,') = ', seed_array(i)
            call da_message(message(1:1))
         end do
      else if (seed_size == 2) then
         seed_array(1) = seed_array1
         seed_array(2) = seed_array2 * seed_array1 + myproc*10000000
         write(unit=message(1),fmt='(a,i16)')' Setting seed_array(1) = ', seed_array(1)
         write(unit=message(2),fmt='(a,i16)')' Setting seed_array(2) = ', seed_array(2)
         call da_message(message(1:2))
      else
         write(unit=message(1),fmt='(a)') ' Error: failure in random number generator'
         write(unit=message(1),fmt='(a)') ' Your compiler does not follow the Fortran 95 standard!'
         call da_error("da_random_seed.inc",73,message(1:2))
      end if
      call random_seed(put=seed_array(1:seed_size)) ! Set random seed.
     
   else                                 ! Random seed set "randomly"
      call random_seed
      call random_seed(get=seed_array(1:seed_size))
      write(unit=message(1),fmt='(a,10i16)') 'Random number seed array = ', seed_array
      call da_message(message(1:1))
   end if
   
   deallocate(seed_array)

   if (trace_use) call da_trace_exit("da_random_seed")

end subroutine da_random_seed


subroutine da_gauss_noise( z)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------
      
   implicit none
   
   real,intent(out)     :: z
   real                 :: x, y, r, coeff

   if (trace_use) call da_trace_entry("da_gauss_noise")

   ! [2.1] Get two uniform variate random numbers in range 0 to 1:

   do
      call random_number( x)
      call random_number( y)

      ! [2.2] Transform to range -1 to 1 and calculate sum of squares:

      x = 2.0 * x - 1.0
      y = 2.0 * y - 1.0
      r = x * x + y * y
      
      if (r > 0.0 .and. r < 1.0) exit        
   end do

   ! [2.3] use Box-Muller transformation to get normal deviates:

   coeff = sqrt( -2.0 * log(r) / r)         
   z = coeff * x

   if (trace_use) call da_trace_exit("da_gauss_noise")
      
end subroutine da_gauss_noise



end module da_define_structures

