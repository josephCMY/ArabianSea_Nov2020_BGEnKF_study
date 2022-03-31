

























MODULE module_alloc_space_1
CONTAINS










   SUBROUTINE alloc_space_field_core_1 ( grid,   id, setinitval_in ,  tl_in , inter_domain_in , okay_to_alloc_in, num_bytes_allocated , &
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







IF ( setinitval .EQ. 3 ) grid%num_moves=0
IF ( setinitval .EQ. 3 ) grid%ts_buf_size=0
IF ( setinitval .EQ. 3 ) grid%max_ts_locs=0
IF ( setinitval .EQ. 3 ) grid%vortex_interval=0
IF ( setinitval .EQ. 3 ) grid%max_vortex_speed=0
IF ( setinitval .EQ. 3 ) grid%corral_dist=0
IF ( setinitval .EQ. 3 ) grid%track_level=0
IF ( setinitval .EQ. 3 ) grid%time_to_move=initial_data_value
IF ( setinitval .EQ. 3 ) grid%move_id=0
IF ( setinitval .EQ. 3 ) grid%move_interval=0
IF ( setinitval .EQ. 3 ) grid%move_cd_x=0
IF ( setinitval .EQ. 3 ) grid%move_cd_y=0
IF ( setinitval .EQ. 3 ) grid%swap_x=.FALSE.
IF ( setinitval .EQ. 3 ) grid%swap_y=.FALSE.
IF ( setinitval .EQ. 3 ) grid%cycle_x=.FALSE.
IF ( setinitval .EQ. 3 ) grid%cycle_y=.FALSE.
IF ( setinitval .EQ. 3 ) grid%reorder_mesh=.FALSE.
IF ( setinitval .EQ. 3 ) grid%perturb_input=.FALSE.
IF ( setinitval .EQ. 3 ) grid%eta_levels=initial_data_value
IF ( setinitval .EQ. 3 ) grid%max_dz=initial_data_value
IF ( setinitval .EQ. 3 ) grid%ocean_levels=0
IF ( setinitval .EQ. 3 ) grid%ocean_z=initial_data_value
IF ( setinitval .EQ. 3 ) grid%ocean_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%ocean_s=initial_data_value
IF ( setinitval .EQ. 3 ) grid%num_traj=0
IF ( setinitval .EQ. 3 ) grid%max_ts_level=0
IF ( setinitval .EQ. 3 ) grid%track_loc_in=0
IF ( setinitval .EQ. 3 ) grid%num_ext_model_couple_dom=0
IF ( setinitval .EQ. 3 ) grid%insert_bogus_storm=.FALSE.
IF ( setinitval .EQ. 3 ) grid%remove_storm=.FALSE.
IF ( setinitval .EQ. 3 ) grid%num_storm=0
IF ( setinitval .EQ. 3 ) grid%latc_loc=initial_data_value
IF ( setinitval .EQ. 3 ) grid%lonc_loc=initial_data_value
IF ( setinitval .EQ. 3 ) grid%vmax_meters_per_second=initial_data_value
IF ( setinitval .EQ. 3 ) grid%rmax=initial_data_value
IF ( setinitval .EQ. 3 ) grid%vmax_ratio=initial_data_value
IF ( setinitval .EQ. 3 ) grid%rankine_lid=initial_data_value
IF ( setinitval .EQ. 3 ) grid%force_read_thompson=.FALSE.
IF ( setinitval .EQ. 3 ) grid%write_thompson_tables=.FALSE.
IF ( setinitval .EQ. 3 ) grid%mp_physics=0
IF ( setinitval .EQ. 3 ) grid%nssl_cccn=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_alphah=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_alphahl=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_cnoh=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_cnohl=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_cnor=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_cnos=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_rho_qh=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_rho_qhl=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nssl_rho_qs=initial_data_value
IF ( setinitval .EQ. 3 ) grid%nudge_lightning=0
IF ( setinitval .EQ. 3 ) grid%nudge_light_times=0
IF ( setinitval .EQ. 3 ) grid%nudge_light_timee=0
IF ( setinitval .EQ. 3 ) grid%nudge_light_int=0
IF ( setinitval .EQ. 3 ) grid%gsfcgce_hail=0
IF ( setinitval .EQ. 3 ) grid%gsfcgce_2ice=0
IF ( setinitval .EQ. 3 ) grid%progn=0
IF ( setinitval .EQ. 3 ) grid%accum_mode=initial_data_value
IF ( setinitval .EQ. 3 ) grid%aitken_mode=initial_data_value
IF ( setinitval .EQ. 3 ) grid%coarse_mode=initial_data_value
IF ( setinitval .EQ. 3 ) grid%do_radar_ref=0
IF ( setinitval .EQ. 3 ) grid%compute_radar_ref=0
IF ( setinitval .EQ. 3 ) grid%ra_lw_physics=0
IF ( setinitval .EQ. 3 ) grid%ra_sw_physics=0
IF ( setinitval .EQ. 3 ) grid%radt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%naer=initial_data_value
IF ( setinitval .EQ. 3 ) grid%sf_sfclay_physics=0
IF ( setinitval .EQ. 3 ) grid%sf_surface_physics=0
IF ( setinitval .EQ. 3 ) grid%bl_pbl_physics=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_tkebudget=0
IF ( setinitval .EQ. 3 ) grid%ysu_topdown_pblmix=0
IF ( setinitval .EQ. 3 ) grid%shinhong_tke_diag=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_tkeadvect=.FALSE.
IF ( setinitval .EQ. 3 ) grid%bl_mynn_cloudpdf=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_mixlength=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_edmf=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_edmf_mom=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_edmf_tke=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_edmf_part=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_cloudmix=0
IF ( setinitval .EQ. 3 ) grid%bl_mynn_mixqt=0
IF ( setinitval .EQ. 3 ) grid%icloud_bl=0
IF ( setinitval .EQ. 3 ) grid%mfshconv=0
IF ( setinitval .EQ. 3 ) grid%sf_urban_physics=0
IF ( setinitval .EQ. 3 ) grid%bldt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%cu_physics=0
IF ( setinitval .EQ. 3 ) grid%shcu_physics=0
IF ( setinitval .EQ. 3 ) grid%cu_diag=0
IF ( setinitval .EQ. 3 ) grid%kf_edrates=0
IF ( setinitval .EQ. 3 ) grid%kfeta_trigger=0
IF ( setinitval .EQ. 3 ) grid%nsas_dx_factor=0
IF ( setinitval .EQ. 3 ) grid%cudt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%gsmdt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%isfflx=0
IF ( setinitval .EQ. 3 ) grid%ifsnow=0
IF ( setinitval .EQ. 3 ) grid%icloud=0
IF ( setinitval .EQ. 3 ) grid%ideal_xland=0
IF ( setinitval .EQ. 3 ) grid%swrad_scat=initial_data_value
IF ( setinitval .EQ. 3 ) grid%surface_input_source=0
IF ( setinitval .EQ. 3 ) grid%num_soil_layers=0
IF ( setinitval .EQ. 3 ) grid%maxpatch=0
IF ( setinitval .EQ. 3 ) grid%num_snow_layers=0
IF ( setinitval .EQ. 3 ) grid%num_snso_layers=0
IF ( setinitval .EQ. 3 ) grid%num_urban_layers=0
IF ( setinitval .EQ. 3 ) grid%num_urban_hi=0
IF ( setinitval .EQ. 3 ) grid%num_months=0
IF ( setinitval .EQ. 3 ) grid%sf_surface_mosaic=0
IF ( setinitval .EQ. 3 ) grid%mosaic_cat=0
IF ( setinitval .EQ. 3 ) grid%mosaic_cat_soil=0
IF ( setinitval .EQ. 3 ) grid%mosaic_lu=0
IF ( setinitval .EQ. 3 ) grid%mosaic_soil=0
IF ( setinitval .EQ. 3 ) grid%maxiens=0
IF ( setinitval .EQ. 3 ) grid%maxens=0
IF ( setinitval .EQ. 3 ) grid%maxens2=0
IF ( setinitval .EQ. 3 ) grid%maxens3=0
IF ( setinitval .EQ. 3 ) grid%ensdim=0
IF ( setinitval .EQ. 3 ) grid%cugd_avedx=0
IF ( setinitval .EQ. 3 ) grid%clos_choice=0
IF ( setinitval .EQ. 3 ) grid%imomentum=0
IF ( setinitval .EQ. 3 ) grid%ishallow=0
IF ( setinitval .EQ. 3 ) grid%convtrans_avglen_m=initial_data_value
IF ( setinitval .EQ. 3 ) grid%num_land_cat=0
IF ( setinitval .EQ. 3 ) grid%num_soil_cat=0
IF ( setinitval .EQ. 3 ) grid%mp_zero_out=0
IF ( setinitval .EQ. 3 ) grid%mp_zero_out_thresh=initial_data_value
IF ( setinitval .EQ. 3 ) grid%seaice_threshold=initial_data_value
IF ( setinitval .EQ. 3 ) grid%sst_update=0
IF ( setinitval .EQ. 3 ) grid%sst_skin=0
IF ( setinitval .EQ. 3 ) grid%tmn_update=0
IF ( setinitval .EQ. 3 ) grid%usemonalb=.FALSE.
IF ( setinitval .EQ. 3 ) grid%rdmaxalb=.FALSE.
IF ( setinitval .EQ. 3 ) grid%rdlai2d=.FALSE.
IF ( setinitval .EQ. 3 ) grid%ua_phys=.FALSE.
IF ( setinitval .EQ. 3 ) grid%opt_thcnd=0
IF ( setinitval .EQ. 3 ) grid%co2tf=0
IF ( setinitval .EQ. 3 ) grid%ra_call_offset=0
IF ( setinitval .EQ. 3 ) grid%cam_abs_freq_s=initial_data_value
IF ( setinitval .EQ. 3 ) grid%levsiz=0
IF ( setinitval .EQ. 3 ) grid%paerlev=0
IF ( setinitval .EQ. 3 ) grid%cam_abs_dim1=0
IF ( setinitval .EQ. 3 ) grid%cam_abs_dim2=0
IF ( setinitval .EQ. 3 ) grid%lagday=0
IF ( setinitval .EQ. 3 ) grid%no_src_types=0
IF ( setinitval .EQ. 3 ) grid%alevsiz=0
IF ( setinitval .EQ. 3 ) grid%o3input=0
IF ( setinitval .EQ. 3 ) grid%aer_opt=0
IF ( setinitval .EQ. 3 ) grid%swint_opt=0
IF ( setinitval .EQ. 3 ) grid%aer_type=0
IF ( setinitval .EQ. 3 ) grid%aer_aod550_opt=0
IF ( setinitval .EQ. 3 ) grid%aer_angexp_opt=0
IF ( setinitval .EQ. 3 ) grid%aer_ssa_opt=0
IF ( setinitval .EQ. 3 ) grid%aer_asy_opt=0
IF ( setinitval .EQ. 3 ) grid%aer_aod550_val=initial_data_value
IF ( setinitval .EQ. 3 ) grid%aer_angexp_val=initial_data_value
IF ( setinitval .EQ. 3 ) grid%aer_ssa_val=initial_data_value
IF ( setinitval .EQ. 3 ) grid%aer_asy_val=initial_data_value
IF ( setinitval .EQ. 3 ) grid%cu_rad_feedback=.FALSE.
IF ( setinitval .EQ. 3 ) grid%shallowcu_forced_ra=.FALSE.
IF ( setinitval .EQ. 3 ) grid%numbins=0
IF ( setinitval .EQ. 3 ) grid%thbinsize=initial_data_value
IF ( setinitval .EQ. 3 ) grid%rbinsize=initial_data_value
IF ( setinitval .EQ. 3 ) grid%mindeepfreq=initial_data_value
IF ( setinitval .EQ. 3 ) grid%minshallowfreq=initial_data_value
IF ( setinitval .EQ. 3 ) grid%shcu_aerosols_opt=0
IF ( setinitval .EQ. 3 ) grid%icloud_cu=0
IF ( setinitval .EQ. 3 ) grid%pxlsm_smois_init=0
IF ( setinitval .EQ. 3 ) grid%omlcall=0
IF ( setinitval .EQ. 3 ) grid%sf_ocean_physics=0
IF ( setinitval .EQ. 3 ) grid%traj_opt=0
IF ( setinitval .EQ. 3 ) grid%tracercall=0
IF ( setinitval .EQ. 3 ) grid%omdt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%oml_hml0=initial_data_value
IF ( setinitval .EQ. 3 ) grid%oml_gamma=initial_data_value
IF ( setinitval .EQ. 3 ) grid%oml_relaxation_time=initial_data_value
IF ( setinitval .EQ. 3 ) grid%isftcflx=0
IF ( setinitval .EQ. 3 ) grid%iz0tlnd=0
IF ( setinitval .EQ. 3 ) grid%shadlen=initial_data_value
IF ( setinitval .EQ. 3 ) grid%slope_rad=0
IF ( setinitval .EQ. 3 ) grid%topo_shading=0
IF ( setinitval .EQ. 3 ) grid%topo_wind=0
IF ( setinitval .EQ. 3 ) grid%no_mp_heating=0
IF ( setinitval .EQ. 3 ) grid%fractional_seaice=0
IF ( setinitval .EQ. 3 ) grid%seaice_snowdepth_opt=0
IF ( setinitval .EQ. 3 ) grid%seaice_snowdepth_max=initial_data_value
IF ( setinitval .EQ. 3 ) grid%seaice_snowdepth_min=initial_data_value
IF ( setinitval .EQ. 3 ) grid%seaice_albedo_opt=0
IF ( setinitval .EQ. 3 ) grid%seaice_albedo_default=initial_data_value
IF ( setinitval .EQ. 3 ) grid%seaice_thickness_opt=0
IF ( setinitval .EQ. 3 ) grid%seaice_thickness_default=initial_data_value
IF ( setinitval .EQ. 3 ) grid%tice2tsk_if2cold=.FALSE.
IF ( setinitval .EQ. 3 ) grid%bucket_mm=initial_data_value
IF ( setinitval .EQ. 3 ) grid%bucket_j=initial_data_value
IF ( setinitval .EQ. 3 ) grid%mp_tend_lim=initial_data_value
IF ( setinitval .EQ. 3 ) grid%prec_acc_dt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%prec_acc_opt=0
IF ( setinitval .EQ. 3 ) grid%bucketr_opt=0
IF ( setinitval .EQ. 3 ) grid%process_time_series=0
IF ( setinitval .EQ. 3 ) grid%grav_settling=0
IF ( setinitval .EQ. 3 ) grid%sas_pgcon=initial_data_value
IF ( setinitval .EQ. 3 ) grid%scalar_pblmix=0
IF ( setinitval .EQ. 3 ) grid%tracer_pblmix=0
IF ( setinitval .EQ. 3 ) grid%use_aero_icbc=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_rap_aero_icbc=.FALSE.
IF ( setinitval .EQ. 3 ) grid%use_mp_re=0
IF ( setinitval .EQ. 3 ) grid%ccn_conc=initial_data_value
IF ( setinitval .EQ. 3 ) grid%hail_opt=0
IF ( setinitval .EQ. 3 ) grid%dveg=0
IF ( setinitval .EQ. 3 ) grid%opt_crs=0
IF ( setinitval .EQ. 3 ) grid%opt_btr=0
IF ( setinitval .EQ. 3 ) grid%opt_run=0
IF ( setinitval .EQ. 3 ) grid%opt_sfc=0
IF ( setinitval .EQ. 3 ) grid%opt_frz=0
IF ( setinitval .EQ. 3 ) grid%opt_inf=0
IF ( setinitval .EQ. 3 ) grid%opt_rad=0
IF ( setinitval .EQ. 3 ) grid%opt_alb=0
IF ( setinitval .EQ. 3 ) grid%opt_snf=0
IF ( setinitval .EQ. 3 ) grid%opt_tbot=0
IF ( setinitval .EQ. 3 ) grid%opt_stc=0
IF ( setinitval .EQ. 3 ) grid%opt_gla=0
IF ( setinitval .EQ. 3 ) grid%opt_rsf=0
IF ( setinitval .EQ. 3 ) grid%wtddt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%wrf_hydro=0
IF ( setinitval .EQ. 3 ) grid%fgdt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%fgdtzero=0
IF ( setinitval .EQ. 3 ) grid%grid_fdda=0
IF ( setinitval .EQ. 3 ) grid%grid_sfdda=0
IF ( setinitval .EQ. 3 ) grid%if_no_pbl_nudging_uv=0
IF ( setinitval .EQ. 3 ) grid%if_no_pbl_nudging_t=0
IF ( setinitval .EQ. 3 ) grid%if_no_pbl_nudging_ph=0
IF ( setinitval .EQ. 3 ) grid%if_no_pbl_nudging_q=0
IF ( setinitval .EQ. 3 ) grid%if_zfac_uv=0
IF ( setinitval .EQ. 3 ) grid%k_zfac_uv=0
IF ( setinitval .EQ. 3 ) grid%if_zfac_t=0
IF ( setinitval .EQ. 3 ) grid%k_zfac_t=0
IF ( setinitval .EQ. 3 ) grid%if_zfac_ph=0
IF ( setinitval .EQ. 3 ) grid%k_zfac_ph=0
IF ( setinitval .EQ. 3 ) grid%if_zfac_q=0
IF ( setinitval .EQ. 3 ) grid%k_zfac_q=0
IF ( setinitval .EQ. 3 ) grid%dk_zfac_uv=0
IF ( setinitval .EQ. 3 ) grid%dk_zfac_t=0
IF ( setinitval .EQ. 3 ) grid%dk_zfac_ph=0
IF ( setinitval .EQ. 3 ) grid%guv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%guv_sfc=initial_data_value
IF ( setinitval .EQ. 3 ) grid%gt=initial_data_value
IF ( setinitval .EQ. 3 ) grid%gt_sfc=initial_data_value
IF ( setinitval .EQ. 3 ) grid%gq=initial_data_value
IF ( setinitval .EQ. 3 ) grid%gq_sfc=initial_data_value
IF ( setinitval .EQ. 3 ) grid%gph=initial_data_value
IF ( setinitval .EQ. 3 ) grid%dtramp_min=initial_data_value
IF ( setinitval .EQ. 3 ) grid%if_ramping=0
IF ( setinitval .EQ. 3 ) grid%rinblw=initial_data_value
IF ( setinitval .EQ. 3 ) grid%xwavenum=0
IF ( setinitval .EQ. 3 ) grid%ywavenum=0
IF ( setinitval .EQ. 3 ) grid%pxlsm_soil_nudge=0
IF ( setinitval .EQ. 3 ) grid%fasdas=0
IF ( setinitval .EQ. 3 ) grid%obs_nudge_opt=0
IF ( setinitval .EQ. 3 ) grid%max_obs=0
IF ( setinitval .EQ. 3 ) grid%fdda_start=initial_data_value
IF ( setinitval .EQ. 3 ) grid%fdda_end=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudge_wind=0
IF ( setinitval .EQ. 3 ) grid%obs_coef_wind=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudge_temp=0
IF ( setinitval .EQ. 3 ) grid%obs_coef_temp=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudge_mois=0
IF ( setinitval .EQ. 3 ) grid%obs_coef_mois=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudge_pstr=0
IF ( setinitval .EQ. 3 ) grid%obs_coef_pstr=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_no_pbl_nudge_uv=0
IF ( setinitval .EQ. 3 ) grid%obs_no_pbl_nudge_t=0
IF ( setinitval .EQ. 3 ) grid%obs_no_pbl_nudge_q=0
IF ( setinitval .EQ. 3 ) grid%obs_sfc_scheme_horiz=0
IF ( setinitval .EQ. 3 ) grid%obs_sfc_scheme_vert=0
IF ( setinitval .EQ. 3 ) grid%obs_max_sndng_gap=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr1_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr1_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr2_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr2_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr4_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr4_uv=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr1_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr1_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr2_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr2_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr4_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr4_t=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr1_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr1_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr2_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr2_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullr4_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampr4_q=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezfullmin=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezrampmin=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_nudgezmax=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_sfcfact=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_sfcfacr=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_dpsmx=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_rinxy=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_rinsig=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_twindo=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_npfi=0
IF ( setinitval .EQ. 3 ) grid%obs_ionf=0
IF ( setinitval .EQ. 3 ) grid%obs_idynin=0
IF ( setinitval .EQ. 3 ) grid%obs_dtramp=initial_data_value
IF ( setinitval .EQ. 3 ) grid%obs_prt_max=0
IF ( setinitval .EQ. 3 ) grid%obs_prt_freq=0
IF ( setinitval .EQ. 3 ) grid%obs_ipf_in4dob=.FALSE.
IF ( setinitval .EQ. 3 ) grid%obs_ipf_errob=.FALSE.
IF ( setinitval .EQ. 3 ) grid%obs_ipf_nudob=.FALSE.
IF ( setinitval .EQ. 3 ) grid%obs_ipf_init=.FALSE.
IF ( setinitval .EQ. 3 ) grid%obs_scl_neg_qv_innov=0
IF ( setinitval .EQ. 3 ) grid%scm_force=0
IF ( setinitval .EQ. 3 ) grid%scm_force_dx=initial_data_value
IF ( setinitval .EQ. 3 ) grid%num_force_layers=0
IF ( setinitval .EQ. 3 ) grid%scm_lu_index=0
IF ( setinitval .EQ. 3 ) grid%scm_isltyp=0
IF ( setinitval .EQ. 3 ) grid%scm_vegfra=initial_data_value
IF ( setinitval .EQ. 3 ) grid%scm_canwat=initial_data_value
IF ( setinitval .EQ. 3 ) grid%scm_lat=initial_data_value
IF ( setinitval .EQ. 3 ) grid%scm_lon=initial_data_value
IF ( setinitval .EQ. 3 ) grid%scm_th_t_tend=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_qv_t_tend=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_th_adv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_wind_adv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_qv_adv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_ql_adv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_vert_adv=.FALSE.
IF ( setinitval .EQ. 3 ) grid%num_force_soil_layers=0
IF ( setinitval .EQ. 3 ) grid%scm_soilt_force=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_soilq_force=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_force_th_largescale=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_force_qv_largescale=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_force_ql_largescale=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_force_wind_largescale=.FALSE.
IF ( setinitval .EQ. 3 ) grid%scm_force_skintemp=0
IF ( setinitval .EQ. 3 ) grid%scm_force_flux=0
IF ( setinitval .EQ. 3 ) grid%dyn_opt=0
IF ( setinitval .EQ. 3 ) grid%rk_ord=0
IF ( setinitval .EQ. 3 ) grid%w_damping=0
IF ( setinitval .EQ. 3 ) grid%diff_opt=0
IF ( setinitval .EQ. 3 ) grid%diff_opt_dfi=0
IF ( setinitval .EQ. 3 ) grid%km_opt=0
IF ( setinitval .EQ. 3 ) grid%km_opt_dfi=0


   END SUBROUTINE alloc_space_field_core_1

END MODULE module_alloc_space_1

