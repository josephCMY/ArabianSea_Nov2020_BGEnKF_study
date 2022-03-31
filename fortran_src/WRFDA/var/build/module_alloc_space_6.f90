

























MODULE module_alloc_space_6
CONTAINS










   SUBROUTINE alloc_space_field_core_6 ( grid,   id, setinitval_in ,  tl_in , inter_domain_in , okay_to_alloc_in, num_bytes_allocated , &
                                  sd31, ed31, sd32, ed32, sd33, ed33, &
                                  sm31 , em31 , sm32 , em32 , sm33 , em33 , &
                                  sp31 , ep31 , sp32 , ep32 , sp33 , ep33 , &
                                  sp31x, ep31x, sp32x, ep32x, sp33x, ep33x, &
                                  sp31y, ep31y, sp32y, ep32y, sp33y, ep33y, &
                                  sm31x, em31x, sm32x, em32x, sm33x, em33x, &
                                  sm31y, em31y, sm32y, em32y, sm33y, em33y )

      USE module_domain_type
      USE module_configure, ONLY : model_config_rec, grid_config_rec_type, in_use_for_config, model_to_grid_config_rec

      USE module_scalar_tables 

      IMPLICIT NONE

      

      TYPE(domain)               , POINTER          :: grid
      INTEGER , INTENT(IN)            :: id
      INTEGER , INTENT(IN)            :: setinitval_in   
      INTEGER , INTENT(IN)            :: sd31, ed31, sd32, ed32, sd33, ed33
      INTEGER , INTENT(IN)            :: sm31, em31, sm32, em32, sm33, em33
      INTEGER , INTENT(IN)            :: sp31, ep31, sp32, ep32, sp33, ep33
      INTEGER , INTENT(IN)            :: sp31x, ep31x, sp32x, ep32x, sp33x, ep33x
      INTEGER , INTENT(IN)            :: sp31y, ep31y, sp32y, ep32y, sp33y, ep33y
      INTEGER , INTENT(IN)            :: sm31x, em31x, sm32x, em32x, sm33x, em33x
      INTEGER , INTENT(IN)            :: sm31y, em31y, sm32y, em32y, sm33y, em33y

      
      
      
      
      INTEGER , INTENT(IN)            :: tl_in
 
      
      
      LOGICAL , INTENT(IN)            :: inter_domain_in, okay_to_alloc_in

      INTEGER(KIND=8) , INTENT(INOUT)         :: num_bytes_allocated


      
      INTEGER idum1, idum2, spec_bdy_width
      REAL    initial_data_value
      CHARACTER (LEN=256) message
      INTEGER tl
      LOGICAL inter_domain, okay_to_alloc
      INTEGER setinitval
      INTEGER sr_x, sr_y

      
      INTEGER ierr

      INTEGER                              :: loop

   

      TYPE ( grid_config_rec_type ) :: config_flags

      INTEGER                         :: k_start , k_end, its, ite, jts, jte
      INTEGER                         :: ids , ide , jds , jde , kds , kde , &
                                         ims , ime , jms , jme , kms , kme , &
                                         ips , ipe , jps , jpe , kps , kpe

      INTEGER                         :: sids , side , sjds , sjde , skds , skde , &
                                         sims , sime , sjms , sjme , skms , skme , &
                                         sips , sipe , sjps , sjpe , skps , skpe


      INTEGER ::              imsx, imex, jmsx, jmex, kmsx, kmex,    &
                              ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                              imsy, imey, jmsy, jmey, kmsy, kmey,    &
                              ipsy, ipey, jpsy, jpey, kpsy, kpey

      data_ordering : SELECT CASE ( model_data_order )
         CASE  ( DATA_ORDER_XYZ )
             ids = sd31 ; ide = ed31 ; jds = sd32 ; jde = ed32 ; kds = sd33 ; kde = ed33 ;
             ims = sm31 ; ime = em31 ; jms = sm32 ; jme = em32 ; kms = sm33 ; kme = em33 ;
             ips = sp31 ; ipe = ep31 ; jps = sp32 ; jpe = ep32 ; kps = sp33 ; kpe = ep33 ;
             imsx = sm31x ; imex = em31x ; jmsx = sm32x ; jmex = em32x ; kmsx = sm33x ; kmex = em33x ;
             ipsx = sp31x ; ipex = ep31x ; jpsx = sp32x ; jpex = ep32x ; kpsx = sp33x ; kpex = ep33x ;
             imsy = sm31y ; imey = em31y ; jmsy = sm32y ; jmey = em32y ; kmsy = sm33y ; kmey = em33y ;
             ipsy = sp31y ; ipey = ep31y ; jpsy = sp32y ; jpey = ep32y ; kpsy = sp33y ; kpey = ep33y ;
         CASE  ( DATA_ORDER_YXZ )
             ids = sd32  ; ide = ed32  ; jds = sd31  ; jde = ed31  ; kds = sd33  ; kde = ed33  ;
             ims = sm32  ; ime = em32  ; jms = sm31  ; jme = em31  ; kms = sm33  ; kme = em33  ;
             ips = sp32  ; ipe = ep32  ; jps = sp31  ; jpe = ep31  ; kps = sp33  ; kpe = ep33  ;
             imsx = sm32x  ; imex = em32x  ; jmsx = sm31x  ; jmex = em31x  ; kmsx = sm33x  ; kmex = em33x  ;
             ipsx = sp32x  ; ipex = ep32x  ; jpsx = sp31x  ; jpex = ep31x  ; kpsx = sp33x  ; kpex = ep33x  ;
             imsy = sm32y  ; imey = em32y  ; jmsy = sm31y  ; jmey = em31y  ; kmsy = sm33y  ; kmey = em33y  ;
             ipsy = sp32y  ; ipey = ep32y  ; jpsy = sp31y  ; jpey = ep31y  ; kpsy = sp33y  ; kpey = ep33y  ;
         CASE  ( DATA_ORDER_ZXY )
             ids = sd32  ; ide = ed32  ; jds = sd33  ; jde = ed33  ; kds = sd31  ; kde = ed31  ;
             ims = sm32  ; ime = em32  ; jms = sm33  ; jme = em33  ; kms = sm31  ; kme = em31  ;
             ips = sp32  ; ipe = ep32  ; jps = sp33  ; jpe = ep33  ; kps = sp31  ; kpe = ep31  ;
             imsx = sm32x  ; imex = em32x  ; jmsx = sm33x  ; jmex = em33x  ; kmsx = sm31x  ; kmex = em31x  ;
             ipsx = sp32x  ; ipex = ep32x  ; jpsx = sp33x  ; jpex = ep33x  ; kpsx = sp31x  ; kpex = ep31x  ;
             imsy = sm32y  ; imey = em32y  ; jmsy = sm33y  ; jmey = em33y  ; kmsy = sm31y  ; kmey = em31y  ;
             ipsy = sp32y  ; ipey = ep32y  ; jpsy = sp33y  ; jpey = ep33y  ; kpsy = sp31y  ; kpey = ep31y  ;
         CASE  ( DATA_ORDER_ZYX )
             ids = sd33  ; ide = ed33  ; jds = sd32  ; jde = ed32  ; kds = sd31  ; kde = ed31  ;
             ims = sm33  ; ime = em33  ; jms = sm32  ; jme = em32  ; kms = sm31  ; kme = em31  ;
             ips = sp33  ; ipe = ep33  ; jps = sp32  ; jpe = ep32  ; kps = sp31  ; kpe = ep31  ;
             imsx = sm33x  ; imex = em33x  ; jmsx = sm32x  ; jmex = em32x  ; kmsx = sm31x  ; kmex = em31x  ;
             ipsx = sp33x  ; ipex = ep33x  ; jpsx = sp32x  ; jpex = ep32x  ; kpsx = sp31x  ; kpex = ep31x  ;
             imsy = sm33y  ; imey = em33y  ; jmsy = sm32y  ; jmey = em32y  ; kmsy = sm31y  ; kmey = em31y  ;
             ipsy = sp33y  ; ipey = ep33y  ; jpsy = sp32y  ; jpey = ep32y  ; kpsy = sp31y  ; kpey = ep31y  ;
         CASE  ( DATA_ORDER_XZY )
             ids = sd31  ; ide = ed31  ; jds = sd33  ; jde = ed33  ; kds = sd32  ; kde = ed32  ;
             ims = sm31  ; ime = em31  ; jms = sm33  ; jme = em33  ; kms = sm32  ; kme = em32  ;
             ips = sp31  ; ipe = ep31  ; jps = sp33  ; jpe = ep33  ; kps = sp32  ; kpe = ep32  ;
             imsx = sm31x  ; imex = em31x  ; jmsx = sm33x  ; jmex = em33x  ; kmsx = sm32x  ; kmex = em32x  ;
             ipsx = sp31x  ; ipex = ep31x  ; jpsx = sp33x  ; jpex = ep33x  ; kpsx = sp32x  ; kpex = ep32x  ;
             imsy = sm31y  ; imey = em31y  ; jmsy = sm33y  ; jmey = em33y  ; kmsy = sm32y  ; kmey = em32y  ;
             ipsy = sp31y  ; ipey = ep31y  ; jpsy = sp33y  ; jpey = ep33y  ; kpsy = sp32y  ; kpey = ep32y  ;
         CASE  ( DATA_ORDER_YZX )
             ids = sd33  ; ide = ed33  ; jds = sd31  ; jde = ed31  ; kds = sd32  ; kde = ed32  ;
             ims = sm33  ; ime = em33  ; jms = sm31  ; jme = em31  ; kms = sm32  ; kme = em32  ;
             ips = sp33  ; ipe = ep33  ; jps = sp31  ; jpe = ep31  ; kps = sp32  ; kpe = ep32  ;
             imsx = sm33x  ; imex = em33x  ; jmsx = sm31x  ; jmex = em31x  ; kmsx = sm32x  ; kmex = em32x  ;
             ipsx = sp33x  ; ipex = ep33x  ; jpsx = sp31x  ; jpex = ep31x  ; kpsx = sp32x  ; kpex = ep32x  ;
             imsy = sm33y  ; imey = em33y  ; jmsy = sm31y  ; jmey = em31y  ; kmsy = sm32y  ; kmey = em32y  ;
             ipsy = sp33y  ; ipey = ep33y  ; jpsy = sp31y  ; jpey = ep31y  ; kpsy = sp32y  ; kpey = ep32y  ;
      END SELECT data_ordering

      CALL model_to_grid_config_rec ( id , model_config_rec , config_flags )

      CALL nl_get_sr_x( id , sr_x )
      CALL nl_get_sr_y( id , sr_y )

      tl = tl_in
      inter_domain = inter_domain_in
      okay_to_alloc = okay_to_alloc_in

      initial_data_value = 0.

      setinitval = setinitval_in

      CALL nl_get_spec_bdy_width( 1, spec_bdy_width )







IF ( setinitval .EQ. 3 ) grid%var4d=.FALSE.
IF ( setinitval .EQ. 3 ) grid%var4d_bin=0
IF ( setinitval .EQ. 3 ) grid%var4d_bin_rain=0
IF ( setinitval .EQ. 3 ) grid%var4d_lbc=.FALSE.
IF ( setinitval .EQ. 3 ) grid%multi_inc=0
IF ( setinitval .EQ. 3 ) grid%print_detail_radar=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_rain=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_xa=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_xb=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_f_obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_map=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_grad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_regression=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_spectral=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_testing=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_parallel=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_be=.FALSE.
IF ( setinitval .EQ. 3 ) grid%print_detail_outerloop=.FALSE.
IF ( setinitval .EQ. 3 ) grid%check_max_iv_print=.FALSE.
IF ( setinitval .EQ. 3 ) grid%check_buddy_print=.FALSE.
IF ( setinitval .EQ. 3 ) grid%analysis_accu=0
IF ( setinitval .EQ. 3 ) grid%calc_w_increment=.FALSE.
IF ( setinitval .EQ. 3 ) grid%dt_cloud_model=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_mod_filtered_obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_buoy=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_synop=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_ships=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_metar=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_sound=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_pilot=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_airep=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_qscat=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_tamdar=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_geoamv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_mtgirs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_polaramv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_sd_profiler=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wind_stats_sd=.FALSE.
IF ( setinitval .EQ. 3 ) grid%qc_rej_both=.FALSE.
IF ( setinitval .EQ. 3 ) grid%fg_format=0
IF ( setinitval .EQ. 3 ) grid%ob_format=0
IF ( setinitval .EQ. 3 ) grid%ob_format_gpsro=0
IF ( setinitval .EQ. 3 ) grid%num_fgat_time=0
IF ( setinitval .EQ. 3 ) grid%thin_conv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%thin_conv_ascii=.FALSE.
IF ( setinitval .EQ. 3 ) grid%thin_mesh_conv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%thin_rainobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_synopobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_shipsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_metarobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_soundobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mtgirsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_tamdarobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_pilotobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_airepobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_geoamvobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_polaramvobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_bogusobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_buoyobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_profilerobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_satemobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_gpsztdobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_gpspwobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_gpsrefobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%top_km_gpsro=initial_data_value
IF ( setinitval .EQ. 3 ) grid%bot_km_gpsro=initial_data_value
IF ( setinitval .EQ. 3 ) grid%use_ssmiretrievalobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_ssmitbobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_ssmt1obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_ssmt2obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_qscatobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_radarobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_radar_rv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_radar_rf=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_radar_rqv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_radar_rhv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_3dvar_phy=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_rainobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_hirs2obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_hirs3obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_hirs4obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mhsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_msuobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_amsuaobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_amsubobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_airsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_airsretobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_eos_amsuaobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_hsbobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_ssmisobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_iasiobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_seviriobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_amsr2obs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_kma1dvar=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_filtered_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_obs_errfac=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_atmsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mwtsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mwhsobs=.FALSE.
IF ( setinitval .EQ. 3 ) grid%check_max_iv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%max_error_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_spd=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_dir=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_omb_spd=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_omb_dir=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_pw=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_ref=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_rh=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_p=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_tb=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_thickness=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_rv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_rf=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_rain=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_buv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_bt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_bq=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_error_slp=initial_data_value
IF ( setinitval .EQ. 3 ) grid%check_buddy=.FALSE.
IF ( setinitval .EQ. 3 ) grid%put_rand_seed=.FALSE.
IF ( setinitval .EQ. 3 ) grid%omb_set_rand=.FALSE.
IF ( setinitval .EQ. 3 ) grid%omb_add_noise=.FALSE.
IF ( setinitval .EQ. 3 ) grid%position_lev_dependant=.FALSE.
IF ( setinitval .EQ. 3 ) grid%obs_qc_pointer=0
IF ( setinitval .EQ. 3 ) grid%qmarker_retain=0
IF ( setinitval .EQ. 3 ) grid%max_sound_input=0
IF ( setinitval .EQ. 3 ) grid%max_mtgirs_input=0
IF ( setinitval .EQ. 3 ) grid%max_tamdar_input=0
IF ( setinitval .EQ. 3 ) grid%max_synop_input=0
IF ( setinitval .EQ. 3 ) grid%max_geoamv_input=0
IF ( setinitval .EQ. 3 ) grid%max_polaramv_input=0
IF ( setinitval .EQ. 3 ) grid%max_airep_input=0
IF ( setinitval .EQ. 3 ) grid%max_satem_input=0
IF ( setinitval .EQ. 3 ) grid%max_pilot_input=0
IF ( setinitval .EQ. 3 ) grid%max_radar_input=0
IF ( setinitval .EQ. 3 ) grid%max_rain_input=0
IF ( setinitval .EQ. 3 ) grid%max_metar_input=0
IF ( setinitval .EQ. 3 ) grid%max_gpspw_input=0
IF ( setinitval .EQ. 3 ) grid%max_ships_input=0
IF ( setinitval .EQ. 3 ) grid%max_profiler_input=0
IF ( setinitval .EQ. 3 ) grid%max_bogus_input=0
IF ( setinitval .EQ. 3 ) grid%max_buoy_input=0
IF ( setinitval .EQ. 3 ) grid%max_ssmi_rv_input=0
IF ( setinitval .EQ. 3 ) grid%max_ssmi_tb_input=0
IF ( setinitval .EQ. 3 ) grid%max_ssmt1_input=0
IF ( setinitval .EQ. 3 ) grid%max_ssmt2_input=0
IF ( setinitval .EQ. 3 ) grid%max_qscat_input=0
IF ( setinitval .EQ. 3 ) grid%max_gpsref_input=0
IF ( setinitval .EQ. 3 ) grid%max_airsr_input=0
IF ( setinitval .EQ. 3 ) grid%max_tovs_input=0
IF ( setinitval .EQ. 3 ) grid%max_ssmis_input=0
IF ( setinitval .EQ. 3 ) grid%report_start=0
IF ( setinitval .EQ. 3 ) grid%report_end=0
IF ( setinitval .EQ. 3 ) grid%tovs_start=0
IF ( setinitval .EQ. 3 ) grid%tovs_end=0
IF ( setinitval .EQ. 3 ) grid%gpsref_thinning=.FALSE.
IF ( setinitval .EQ. 3 ) grid%outer_loop_restart=.FALSE.
IF ( setinitval .EQ. 3 ) grid%max_ext_its=0
IF ( setinitval .EQ. 3 ) grid%ntmax=0
IF ( setinitval .EQ. 3 ) grid%nsave=0
IF ( setinitval .EQ. 3 ) grid%write_interval=0
IF ( setinitval .EQ. 3 ) grid%eps=initial_data_value
IF ( setinitval .EQ. 3 ) grid%precondition_cg=.FALSE.
IF ( setinitval .EQ. 3 ) grid%precondition_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%use_lanczos=.FALSE.
IF ( setinitval .EQ. 3 ) grid%read_lanczos=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_lanczos=.FALSE.
IF ( setinitval .EQ. 3 ) grid%orthonorm_gradient=.FALSE.
IF ( setinitval .EQ. 3 ) grid%cv_options=0
IF ( setinitval .EQ. 3 ) grid%cloud_cv_options=0
IF ( setinitval .EQ. 3 ) grid%as1=initial_data_value
IF ( setinitval .EQ. 3 ) grid%as2=initial_data_value
IF ( setinitval .EQ. 3 ) grid%as3=initial_data_value
IF ( setinitval .EQ. 3 ) grid%as4=initial_data_value
IF ( setinitval .EQ. 3 ) grid%as5=initial_data_value
IF ( setinitval .EQ. 3 ) grid%do_normalize=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_rf=.FALSE.
IF ( setinitval .EQ. 3 ) grid%rf_passes=0
IF ( setinitval .EQ. 3 ) grid%var_scaling1=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling2=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling3=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling4=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling5=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling6=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling7=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling8=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling9=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling10=initial_data_value
IF ( setinitval .EQ. 3 ) grid%var_scaling11=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling1=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling2=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling3=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling4=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling5=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling6=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling7=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling8=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling9=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling10=initial_data_value
IF ( setinitval .EQ. 3 ) grid%len_scaling11=initial_data_value
IF ( setinitval .EQ. 3 ) grid%je_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%power_truncation=initial_data_value
IF ( setinitval .EQ. 3 ) grid%def_sub_domain=.FALSE.
IF ( setinitval .EQ. 3 ) grid%x_start_sub_domain=initial_data_value
IF ( setinitval .EQ. 3 ) grid%y_start_sub_domain=initial_data_value
IF ( setinitval .EQ. 3 ) grid%x_end_sub_domain=initial_data_value
IF ( setinitval .EQ. 3 ) grid%y_end_sub_domain=initial_data_value
IF ( setinitval .EQ. 3 ) grid%stdout=0
IF ( setinitval .EQ. 3 ) grid%stderr=0
IF ( setinitval .EQ. 3 ) grid%trace_unit=0
IF ( setinitval .EQ. 3 ) grid%trace_pe=0
IF ( setinitval .EQ. 3 ) grid%trace_repeat_head=0
IF ( setinitval .EQ. 3 ) grid%trace_repeat_body=0
IF ( setinitval .EQ. 3 ) grid%trace_max_depth=0
IF ( setinitval .EQ. 3 ) grid%trace_use=.FALSE.
IF ( setinitval .EQ. 3 ) grid%trace_use_frequent=.FALSE.
IF ( setinitval .EQ. 3 ) grid%trace_use_dull=.FALSE.
IF ( setinitval .EQ. 3 ) grid%trace_memory=.FALSE.
IF ( setinitval .EQ. 3 ) grid%trace_all_pes=.FALSE.
IF ( setinitval .EQ. 3 ) grid%trace_csv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_html=.FALSE.
IF ( setinitval .EQ. 3 ) grid%warnings_are_fatal=.FALSE.
IF ( setinitval .EQ. 3 ) grid%test_transforms=.FALSE.
IF ( setinitval .EQ. 3 ) grid%test_gradient=.FALSE.
IF ( setinitval .EQ. 3 ) grid%test_statistics=.FALSE.
IF ( setinitval .EQ. 3 ) grid%interpolate_stats=.FALSE.
IF ( setinitval .EQ. 3 ) grid%be_eta=initial_data_value
IF ( setinitval .EQ. 3 ) grid%test_dm_exact=.FALSE.
IF ( setinitval .EQ. 3 ) grid%cv_options_hum=0
IF ( setinitval .EQ. 3 ) grid%check_rh=0
IF ( setinitval .EQ. 3 ) grid%set_omb_rand_fac=initial_data_value
IF ( setinitval .EQ. 3 ) grid%seed_array1=0
IF ( setinitval .EQ. 3 ) grid%seed_array2=0
IF ( setinitval .EQ. 3 ) grid%sfc_assi_options=0
IF ( setinitval .EQ. 3 ) grid%psfc_from_slp=.FALSE.
IF ( setinitval .EQ. 3 ) grid%calculate_cg_cost_fn=.FALSE.
IF ( setinitval .EQ. 3 ) grid%lat_stats_option=.FALSE.
IF ( setinitval .EQ. 3 ) grid%interp_option=0
IF ( setinitval .EQ. 3 ) grid%balance_type=0
IF ( setinitval .EQ. 3 ) grid%use_wpec=.FALSE.
IF ( setinitval .EQ. 3 ) grid%wpec_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%vert_corr=0
IF ( setinitval .EQ. 3 ) grid%vertical_ip=0
IF ( setinitval .EQ. 3 ) grid%vert_evalue=0
IF ( setinitval .EQ. 3 ) grid%max_vert_var1=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var2=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var3=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var4=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var5=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var6=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var7=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var8=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var9=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var10=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var11=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_vert_var_alpha=initial_data_value
IF ( setinitval .EQ. 3 ) grid%psi_chi_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%psi_t_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%psi_ps_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%psi_rh_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%chi_u_t_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%chi_u_ps_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%chi_u_rh_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%t_u_rh_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%ps_u_rh_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%rttov_emis_atlas_ir=0
IF ( setinitval .EQ. 3 ) grid%rttov_emis_atlas_mw=0
IF ( setinitval .EQ. 3 ) grid%rtminit_print=0
IF ( setinitval .EQ. 3 ) grid%rtminit_nsensor=0
IF ( setinitval .EQ. 3 ) grid%rtminit_platform=0
IF ( setinitval .EQ. 3 ) grid%rtminit_satid=0
IF ( setinitval .EQ. 3 ) grid%rtminit_sensor=0
IF ( setinitval .EQ. 3 ) grid%rad_monitoring=0
IF ( setinitval .EQ. 3 ) grid%thinning_mesh=initial_data_value
IF ( setinitval .EQ. 3 ) grid%thinning=.FALSE.
IF ( setinitval .EQ. 3 ) grid%read_biascoef=.FALSE.
IF ( setinitval .EQ. 3 ) grid%biascorr=.FALSE.
IF ( setinitval .EQ. 3 ) grid%biasprep=.FALSE.
IF ( setinitval .EQ. 3 ) grid%rttov_scatt=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_profile=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_jacobian=.FALSE.
IF ( setinitval .EQ. 3 ) grid%qc_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_iv_rad_ascii=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_oa_rad_ascii=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_filtered_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_error_factor_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_landem=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_antcorr=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mspps_emis=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mspps_ts=.FALSE.
IF ( setinitval .EQ. 3 ) grid%mw_emis_sea=0
IF ( setinitval .EQ. 3 ) grid%tovs_min_transfer=0
IF ( setinitval .EQ. 3 ) grid%tovs_batch=.FALSE.
IF ( setinitval .EQ. 3 ) grid%rtm_option=0
IF ( setinitval .EQ. 3 ) grid%use_crtm_kmatrix=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_rttov_kmatrix=.FALSE.
IF ( setinitval .EQ. 3 ) grid%crtm_cloud=.FALSE.
IF ( setinitval .EQ. 3 ) grid%only_sea_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_pseudo_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_platid=0
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_satid=0
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_senid=0
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_ichan=0
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_lat=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_lon=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_inv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_rad_err=initial_data_value
IF ( setinitval .EQ. 3 ) grid%use_simulated_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%simulated_rad_io=.FALSE.
IF ( setinitval .EQ. 3 ) grid%simulated_rad_ngrid=0
IF ( setinitval .EQ. 3 ) grid%use_varbc=.FALSE.
IF ( setinitval .EQ. 3 ) grid%freeze_varbc=.FALSE.
IF ( setinitval .EQ. 3 ) grid%varbc_factor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%varbc_nbgerr=0
IF ( setinitval .EQ. 3 ) grid%varbc_nobsmin=0
IF ( setinitval .EQ. 3 ) grid%use_clddet_mmr=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_clddet_ecmwf=.FALSE.
IF ( setinitval .EQ. 3 ) grid%airs_warmest_fov=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_satcv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_blacklist_rad=.FALSE.
IF ( setinitval .EQ. 3 ) grid%calc_weightfunc=.FALSE.
IF ( setinitval .EQ. 3 ) grid%num_pseudo=0
IF ( setinitval .EQ. 3 ) grid%pseudo_x=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_y=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_z=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_val=initial_data_value
IF ( setinitval .EQ. 3 ) grid%pseudo_err=initial_data_value
IF ( setinitval .EQ. 3 ) grid%alphacv_method=0
IF ( setinitval .EQ. 3 ) grid%ensdim_alpha=0
IF ( setinitval .EQ. 3 ) grid%alpha_truncation=0
IF ( setinitval .EQ. 3 ) grid%alpha_corr_type=0
IF ( setinitval .EQ. 3 ) grid%alpha_corr_scale=initial_data_value
IF ( setinitval .EQ. 3 ) grid%alpha_std_dev=initial_data_value
IF ( setinitval .EQ. 3 ) grid%alpha_vertloc=.FALSE.


   END SUBROUTINE alloc_space_field_core_6

END MODULE module_alloc_space_6

