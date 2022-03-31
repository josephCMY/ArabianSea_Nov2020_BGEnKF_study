












module da_wrfvar_top

   
   
   

   use module_configure, only : grid_config_rec_type,model_config_rec, &
      model_to_grid_config_rec, get_config_as_buffer,set_config_as_buffer, &
      initial_config
   use module_domain, only : domain,alloc_and_configure_domain, head_grid, &
      program_name, domain_clock_get, domain_clock_set, x_type, dealloc_space_domain, &
      domain_destroy
   use module_driver_constants, only : max_comms
   use module_symbols_util, only : wrfu_finalize, wrfu_initialize, &
      wrfu_cal_gregorian
   use module_io_domain, only : close_dataset

   use module_radiance, only : satinfo

   use module_state_description, only : num_moist, num_a_moist, num_g_moist, &
      num_scalar, num_a_scalar, num_g_scalar, &
      num_chem, PARAM_FIRST_SCALAR, num_tracer 
   use module_tiles, only : set_tiles


   use module_dm, only : local_communicator, local_communicator_x, &
      local_communicator_y, ntasks_x, ntasks_y, data_order_xyz, mytask, &
      ntasks, data_order_xy,wrf_dm_initialize
   use module_comm_dm, only : halo_radar_xa_w_sub,halo_ssmi_xa_sub, &
      halo_sfc_xa_sub, halo_xa_sub, halo_psichi_uv_adj_sub, halo_bal_eqn_adj_sub, &
      halo_psichi_uv_sub, halo_init_sub, halo_psichi_uv_adj_sub, halo_2d_work_sub,&
      halo_wpec_sub, halo_wpec_adj_sub, halo_xa_all_sub, halo_xb_all_sub

   
   use da_control
   use da_define_structures, only : y_type, j_type, iv_type, be_type, &
      xbx_type,da_deallocate_background_errors,da_initialize_cv, &
      da_zero_vp_type,da_allocate_y,da_deallocate_observations, &
      da_deallocate_y, da_zero_x
   use da_minimisation, only : da_get_innov_vector,da_minimise_cg, &
      da_minimise_lz, da_write_diagnostics, da_calculate_residual, &
      da_calculate_grady, da_sensitivity, da_lanczos_io, da_calculate_j, &
      da_kmat_mul
   use da_obs, only : da_transform_xtoy_adj 
   use da_obs_io, only : da_write_filtered_obs, da_write_obs, da_final_write_obs , &
                         da_write_obs_etkf, da_write_modified_filtered_obs
   use da_par_util, only : da_system,da_copy_tile_dims,da_copy_dims
   use da_physics, only : da_uvprho_to_w_lin
   use da_radiance, only : da_deallocate_radiance
   use da_radiance1, only : num_tovs_before, tovs_recv_pe,tovs_copy_count, &
      tovs_send_pe,tovs_send_count,tovs_recv_start, num_tovs_after, &
      tovs_send_start, da_oi_stats_rad, da_write_oa_rad_ascii, da_setup_satcv
   use da_varbc, only : da_varbc_init,da_varbc_update
   use da_reporting, only : message, da_warning, da_error, da_message
   use da_setup_structures, only : da_setup_obs_structures, &
      da_setup_background_errors,da_setup_flow_predictors, &
      da_setup_cv, da_scale_background_errors, da_scale_background_errors_cv3
   use da_test, only : da_check, da_check_gradient
   use da_tools_serial, only : da_get_unit, da_free_unit
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace, da_trace_report
   use da_transfer_model, only : da_transfer_xatoanalysis,da_setup_firstguess, &
       da_transfer_wrftltoxa_adj
   use da_vtox_transforms, only : da_transform_vtox, da_transform_xtoxa, &
      da_transform_xtoxa_adj
   use da_wrfvar_io, only : da_med_initialdata_input, da_update_firstguess
   use da_tools, only : da_set_randomcv, da_get_julian_time

   use da_tools, only : map_info,map_info_ens,proj_merc, proj_ps,proj_lc,proj_latlon, &
      da_llxy_default,da_llxy_wrf,da_xyll,da_diff_seconds,da_map_set, &
      da_set_boundary_xb,da_togrid

   use module_radiance, only : crtm_destroy
   use da_crtm, only : channelinfo, sensor_descriptor

   use da_airep, only : da_oi_stats_airep
   use da_airsr , only : da_oi_stats_airsr
   use da_bogus, only : da_oi_stats_bogus
   use da_buoy , only : da_oi_stats_buoy
   use da_geoamv, only : da_oi_stats_geoamv
   use da_gpspw, only : da_oi_stats_gpspw
   use da_gpsref, only : da_oi_stats_gpsref
   use da_metar, only : da_oi_stats_metar
   use da_pilot, only : da_oi_stats_pilot
   use da_polaramv, only : da_oi_stats_polaramv
   use da_profiler, only : da_oi_stats_profiler
   use da_qscat, only : da_oi_stats_qscat
   use da_mtgirs, only : da_oi_stats_mtgirs
   use da_radar, only : da_oi_stats_radar, da_write_oa_radar_ascii
   use da_satem, only : da_oi_stats_satem
   use da_ships, only : da_oi_stats_ships
   use da_sound, only : da_oi_stats_sound, da_oi_stats_sonde_sfc
   use da_ssmi, only : da_oi_stats_ssmt1, da_oi_stats_ssmt2, da_oi_stats_ssmi_tb, da_oi_stats_ssmi_rv
   use da_synop, only : da_oi_stats_synop  
   use da_rain, only : da_oi_stats_rain

   use da_wrf_interfaces

   use da_netcdf_interface, only : da_get_var_2d_real_cdf

   implicit none

   integer :: loop, levels_to_process

   type (domain) , pointer :: keep_grid, grid_ptr, null_domain
   type (domain) , pointer :: another_grid, parent_grid, ensemble_grid, input_grid
   type (grid_config_rec_type), save :: config_flags
   integer                 :: number_at_same_level
   integer                 :: time_step_begin_restart

   integer :: domain_id , fid , oid , idum1 , idum2

   integer                 :: nbytes
   integer, parameter      :: configbuflen = 4* 65536
   integer                 :: configbuf( configbuflen )

   character (len=80)      :: rstname


contains

subroutine da_wrfvar_init1(no_init1)

   !-----------------------------------------------------------------------
   ! Purpose: WRFVAR initialization routine, part 1
   !-----------------------------------------------------------------------

   implicit none

   logical, optional, intent(in) :: no_init1

   !<DESCRIPTION>
   ! Program_name, a global variable defined in frame/module_domain.F, is
   ! set, then a routine <a href=init_modules.html>init_modules</a> is
   ! called. This calls all the init programs that are provided by the
   ! modules that are linked into WRF.  These include initialization of
   ! external I/O packages.   Also, some key initializations for
   ! distributed-memory parallelism occur here if 1 is specified
   ! in the compile: setting up I/O quilt processes to act as I/O servers
   ! and dividing up MPI communicators among those as well as initializing
   ! external communication packages.
   !
   !</DESCRIPTION>

   ! Set "program_name", which will be printed to output and netcdf metadata
   program_name = "WRFDA "//release_version

   ! Initialize WRF modules:  
   ! Phase 1 returns after mpi_init() (if it is called)
   if (.NOT. present(no_init1)) then
      call init_modules (1)
      ! Initialize utilities (time manager, etc.)
      call wrfu_initialize (defaultCalKind=WRFU_CAL_GREGORIAN)
   end if
   ! Phase 2 resumes after mpi_init() (if it is called)
   call init_modules (2)

   !<DESCRIPTION>
   ! The wrf namelist.input file is read and stored in the use associated
   ! structure model_config_rec, defined in frame/module_configure.F, by the
   ! call to <a href=initial_config.html>initial_config</a>.  On distributed
   ! memory parallel runs this is done only on one processor, and then
   ! broadcast as a buffer.  For distributed-memory, the broadcast of the
   ! configuration information is accomplished by first putting the
   ! configuration information into a buffer (<a
   ! href=get_config_as_buffer.html>get_config_as_buffer</a>), broadcasting
   ! the buffer, then setting the configuration information (<a
   ! href=set_config_as_buffer.html>set_config_as_buffer</a>).
   !
   !</DESCRIPTION>

   ! Don't use stdout here, too early
   write(unit=6,fmt='(/A)') "*** VARIATIONAL ANALYSIS ***"
   write(unit=6,fmt='(A,A/)') "    ",program_name

   call wrf_get_dm_communicator (comm)
   call mpi_comm_size (comm, num_procs, ierr)
   call mpi_comm_rank (comm, myproc, ierr)

   if (myproc==0) then
      rootproc=.true.
   else
      rootproc=.false.
   end if

   if (rootproc) then
      call initial_config
   end if
   call get_config_as_buffer (configbuf, configbuflen, nbytes)
   call wrf_dm_bcast_bytes (configbuf, nbytes)
   call set_config_as_buffer (configbuf, configbuflen)
   call wrf_dm_initialize

   ! Copy namelist variables to da_control


!STARTOFREGISTRYGENERATEDINCLUDE 'inc/config_assigns.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
! Contains config assign statements for module_domain.F.
  run_days                   = model_config_rec% run_days 
  run_hours                  = model_config_rec% run_hours 
  run_minutes                = model_config_rec% run_minutes 
  run_seconds                = model_config_rec% run_seconds 
  start_year                 = model_config_rec% start_year 
  start_month                = model_config_rec% start_month 
  start_day                  = model_config_rec% start_day 
  start_hour                 = model_config_rec% start_hour 
  start_minute               = model_config_rec% start_minute 
  start_second               = model_config_rec% start_second 
  end_year                   = model_config_rec% end_year 
  end_month                  = model_config_rec% end_month 
  end_day                    = model_config_rec% end_day 
  end_hour                   = model_config_rec% end_hour 
  end_minute                 = model_config_rec% end_minute 
  end_second                 = model_config_rec% end_second 
  interval_seconds           = model_config_rec% interval_seconds 
  input_from_file            = model_config_rec% input_from_file 
  fine_input_stream          = model_config_rec% fine_input_stream 
  input_from_hires           = model_config_rec% input_from_hires 
  rsmas_data_path            = model_config_rec% rsmas_data_path 
  all_ic_times               = model_config_rec% all_ic_times 
  julyr                      = model_config_rec% julyr 
  julday                     = model_config_rec% julday 
  gmt                        = model_config_rec% gmt 
  input_inname               = model_config_rec% input_inname 
  input_outname              = model_config_rec% input_outname 
  bdy_inname                 = model_config_rec% bdy_inname 
  bdy_outname                = model_config_rec% bdy_outname 
  rst_inname                 = model_config_rec% rst_inname 
  rst_outname                = model_config_rec% rst_outname 
  write_input                = model_config_rec% write_input 
  write_restart_at_0h        = model_config_rec% write_restart_at_0h 
  write_hist_at_0h_rst       = model_config_rec% write_hist_at_0h_rst 
  adjust_output_times        = model_config_rec% adjust_output_times 
  adjust_input_times         = model_config_rec% adjust_input_times 
  diag_print                 = model_config_rec% diag_print 
  nocolons                   = model_config_rec% nocolons 
  cycling                    = model_config_rec% cycling 
  output_diagnostics         = model_config_rec% output_diagnostics 
  nwp_diagnostics            = model_config_rec% nwp_diagnostics 
  output_ready_flag          = model_config_rec% output_ready_flag 
  usepio                     = model_config_rec% usepio 
  pioprocs                   = model_config_rec% pioprocs 
  piostart                   = model_config_rec% piostart 
  piostride                  = model_config_rec% piostride 
  pioshift                   = model_config_rec% pioshift 
  dfi_opt                    = model_config_rec% dfi_opt 
  dfi_savehydmeteors         = model_config_rec% dfi_savehydmeteors 
  dfi_nfilter                = model_config_rec% dfi_nfilter 
  dfi_write_filtered_input   = model_config_rec% dfi_write_filtered_input 
  dfi_write_dfi_history      = model_config_rec% dfi_write_dfi_history 
  dfi_cutoff_seconds         = model_config_rec% dfi_cutoff_seconds 
  dfi_time_dim               = model_config_rec% dfi_time_dim 
  dfi_fwdstop_year           = model_config_rec% dfi_fwdstop_year 
  dfi_fwdstop_month          = model_config_rec% dfi_fwdstop_month 
  dfi_fwdstop_day            = model_config_rec% dfi_fwdstop_day 
  dfi_fwdstop_hour           = model_config_rec% dfi_fwdstop_hour 
  dfi_fwdstop_minute         = model_config_rec% dfi_fwdstop_minute 
  dfi_fwdstop_second         = model_config_rec% dfi_fwdstop_second 
  dfi_bckstop_year           = model_config_rec% dfi_bckstop_year 
  dfi_bckstop_month          = model_config_rec% dfi_bckstop_month 
  dfi_bckstop_day            = model_config_rec% dfi_bckstop_day 
  dfi_bckstop_hour           = model_config_rec% dfi_bckstop_hour 
  dfi_bckstop_minute         = model_config_rec% dfi_bckstop_minute 
  dfi_bckstop_second         = model_config_rec% dfi_bckstop_second 
  time_step                  = model_config_rec% time_step 
  time_step_fract_num        = model_config_rec% time_step_fract_num 
  time_step_fract_den        = model_config_rec% time_step_fract_den 
  time_step_dfi              = model_config_rec% time_step_dfi 
  min_time_step              = model_config_rec% min_time_step 
  min_time_step_den          = model_config_rec% min_time_step_den 
  max_time_step              = model_config_rec% max_time_step 
  max_time_step_den          = model_config_rec% max_time_step_den 
  target_cfl                 = model_config_rec% target_cfl 
  target_hcfl                = model_config_rec% target_hcfl 
  max_step_increase_pct      = model_config_rec% max_step_increase_pct 
  starting_time_step         = model_config_rec% starting_time_step 
  starting_time_step_den     = model_config_rec% starting_time_step_den 
  step_to_output_time        = model_config_rec% step_to_output_time 
  adaptation_domain          = model_config_rec% adaptation_domain 
  use_adaptive_time_step     = model_config_rec% use_adaptive_time_step 
  use_adaptive_time_step_dfi = model_config_rec% use_adaptive_time_step_dfi 
  max_dom                    = model_config_rec% max_dom 
  lats_to_mic                = model_config_rec% lats_to_mic 
  s_we                       = model_config_rec% s_we 
  e_we                       = model_config_rec% e_we 
  s_sn                       = model_config_rec% s_sn 
  e_sn                       = model_config_rec% e_sn 
  s_vert                     = model_config_rec% s_vert 
  e_vert                     = model_config_rec% e_vert 
  num_metgrid_levels         = model_config_rec% num_metgrid_levels 
  num_metgrid_soil_levels    = model_config_rec% num_metgrid_soil_levels 
  p_top_requested            = model_config_rec% p_top_requested 
  interp_theta               = model_config_rec% interp_theta 
  interp_type                = model_config_rec% interp_type 
  rebalance                  = model_config_rec% rebalance 
  vert_refine_method         = model_config_rec% vert_refine_method 
  vert_refine_fact           = model_config_rec% vert_refine_fact 
  extrap_type                = model_config_rec% extrap_type 
  t_extrap_type              = model_config_rec% t_extrap_type 
  hypsometric_opt            = model_config_rec% hypsometric_opt 
  lowest_lev_from_sfc        = model_config_rec% lowest_lev_from_sfc 
  use_levels_below_ground    = model_config_rec% use_levels_below_ground 
  use_tavg_for_tsk           = model_config_rec% use_tavg_for_tsk 
  use_surface                = model_config_rec% use_surface 
  lagrange_order             = model_config_rec% lagrange_order 
  force_sfc_in_vinterp       = model_config_rec% force_sfc_in_vinterp 
  zap_close_levels           = model_config_rec% zap_close_levels 
  maxw_horiz_pres_diff       = model_config_rec% maxw_horiz_pres_diff 
  trop_horiz_pres_diff       = model_config_rec% trop_horiz_pres_diff 
  maxw_above_this_level      = model_config_rec% maxw_above_this_level 
  use_maxw_level             = model_config_rec% use_maxw_level 
  use_trop_level             = model_config_rec% use_trop_level 
  sfcp_to_sfcp               = model_config_rec% sfcp_to_sfcp 
  adjust_heights             = model_config_rec% adjust_heights 
  smooth_cg_topo             = model_config_rec% smooth_cg_topo 
  nest_interp_coord          = model_config_rec% nest_interp_coord 
  interp_method_type         = model_config_rec% interp_method_type 
  aggregate_lu               = model_config_rec% aggregate_lu 
  rh2qv_wrt_liquid           = model_config_rec% rh2qv_wrt_liquid 
  rh2qv_method               = model_config_rec% rh2qv_method 
  qv_max_p_safe              = model_config_rec% qv_max_p_safe 
  qv_max_flag                = model_config_rec% qv_max_flag 
  qv_max_value               = model_config_rec% qv_max_value 
  qv_min_p_safe              = model_config_rec% qv_min_p_safe 
  qv_min_flag                = model_config_rec% qv_min_flag 
  qv_min_value               = model_config_rec% qv_min_value 
  ideal_init_method          = model_config_rec% ideal_init_method 
  dx                         = model_config_rec% dx 
  dy                         = model_config_rec% dy 
  grid_id                    = model_config_rec% grid_id 
  grid_allowed               = model_config_rec% grid_allowed 
  parent_id                  = model_config_rec% parent_id 
  i_parent_start             = model_config_rec% i_parent_start 
  j_parent_start             = model_config_rec% j_parent_start 
  parent_grid_ratio          = model_config_rec% parent_grid_ratio 
  parent_time_step_ratio     = model_config_rec% parent_time_step_ratio 
  feedback                   = model_config_rec% feedback 
  smooth_option              = model_config_rec% smooth_option 
  blend_width                = model_config_rec% blend_width 
  ztop                       = model_config_rec% ztop 
  moad_grid_ratio            = model_config_rec% moad_grid_ratio 
  moad_time_step_ratio       = model_config_rec% moad_time_step_ratio 
  shw                        = model_config_rec% shw 
  tile_sz_x                  = model_config_rec% tile_sz_x 
  tile_sz_y                  = model_config_rec% tile_sz_y 
  numtiles                   = model_config_rec% numtiles 
  numtiles_inc               = model_config_rec% numtiles_inc 
  numtiles_x                 = model_config_rec% numtiles_x 
  numtiles_y                 = model_config_rec% numtiles_y 
  tile_strategy              = model_config_rec% tile_strategy 
  nproc_x                    = model_config_rec% nproc_x 
  nproc_y                    = model_config_rec% nproc_y 
  irand                      = model_config_rec% irand 
  dt                         = model_config_rec% dt 
  fft_used                   = model_config_rec% fft_used 
  cu_used                    = model_config_rec% cu_used 
  shcu_used                  = model_config_rec% shcu_used 
  cam_used                   = model_config_rec% cam_used 
  alloc_qndropsource         = model_config_rec% alloc_qndropsource 
  num_moves                  = model_config_rec% num_moves 
  ts_buf_size                = model_config_rec% ts_buf_size 
  max_ts_locs                = model_config_rec% max_ts_locs 
  vortex_interval            = model_config_rec% vortex_interval 
  max_vortex_speed           = model_config_rec% max_vortex_speed 
  corral_dist                = model_config_rec% corral_dist 
  track_level                = model_config_rec% track_level 
  time_to_move               = model_config_rec% time_to_move 
  move_id                    = model_config_rec% move_id 
  move_interval              = model_config_rec% move_interval 
  move_cd_x                  = model_config_rec% move_cd_x 
  move_cd_y                  = model_config_rec% move_cd_y 
  swap_x                     = model_config_rec% swap_x 
  swap_y                     = model_config_rec% swap_y 
  cycle_x                    = model_config_rec% cycle_x 
  cycle_y                    = model_config_rec% cycle_y 
  reorder_mesh               = model_config_rec% reorder_mesh 
  perturb_input              = model_config_rec% perturb_input 
  eta_levels                 = model_config_rec% eta_levels 
  max_dz                     = model_config_rec% max_dz 
  ocean_levels               = model_config_rec% ocean_levels 
  ocean_z                    = model_config_rec% ocean_z 
  ocean_t                    = model_config_rec% ocean_t 
  ocean_s                    = model_config_rec% ocean_s 
  num_traj                   = model_config_rec% num_traj 
  max_ts_level               = model_config_rec% max_ts_level 
  track_loc_in               = model_config_rec% track_loc_in 
  num_ext_model_couple_dom   = model_config_rec% num_ext_model_couple_dom 
  insert_bogus_storm         = model_config_rec% insert_bogus_storm 
  remove_storm               = model_config_rec% remove_storm 
  num_storm                  = model_config_rec% num_storm 
  latc_loc                   = model_config_rec% latc_loc 
  lonc_loc                   = model_config_rec% lonc_loc 
  vmax_meters_per_second     = model_config_rec% vmax_meters_per_second 
  rmax                       = model_config_rec% rmax 
  vmax_ratio                 = model_config_rec% vmax_ratio 
  rankine_lid                = model_config_rec% rankine_lid 
  force_read_thompson        = model_config_rec% force_read_thompson 
  write_thompson_tables      = model_config_rec% write_thompson_tables 
  mp_physics                 = model_config_rec% mp_physics 
  nssl_cccn                  = model_config_rec% nssl_cccn 
  nssl_alphah                = model_config_rec% nssl_alphah 
  nssl_alphahl               = model_config_rec% nssl_alphahl 
  nssl_cnoh                  = model_config_rec% nssl_cnoh 
  nssl_cnohl                 = model_config_rec% nssl_cnohl 
  nssl_cnor                  = model_config_rec% nssl_cnor 
  nssl_cnos                  = model_config_rec% nssl_cnos 
  nssl_rho_qh                = model_config_rec% nssl_rho_qh 
  nssl_rho_qhl               = model_config_rec% nssl_rho_qhl 
  nssl_rho_qs                = model_config_rec% nssl_rho_qs 
  nudge_lightning            = model_config_rec% nudge_lightning 
  nudge_light_times          = model_config_rec% nudge_light_times 
  nudge_light_timee          = model_config_rec% nudge_light_timee 
  nudge_light_int            = model_config_rec% nudge_light_int 
  path_to_files              = model_config_rec% path_to_files 
  gsfcgce_hail               = model_config_rec% gsfcgce_hail 
  gsfcgce_2ice               = model_config_rec% gsfcgce_2ice 
  progn                      = model_config_rec% progn 
  accum_mode                 = model_config_rec% accum_mode 
  aitken_mode                = model_config_rec% aitken_mode 
  coarse_mode                = model_config_rec% coarse_mode 
  do_radar_ref               = model_config_rec% do_radar_ref 
  compute_radar_ref          = model_config_rec% compute_radar_ref 
  ra_lw_physics              = model_config_rec% ra_lw_physics 
  ra_sw_physics              = model_config_rec% ra_sw_physics 
  radt                       = model_config_rec% radt 
  naer                       = model_config_rec% naer 
  sf_sfclay_physics          = model_config_rec% sf_sfclay_physics 
  sf_surface_physics         = model_config_rec% sf_surface_physics 
  bl_pbl_physics             = model_config_rec% bl_pbl_physics 
  bl_mynn_tkebudget          = model_config_rec% bl_mynn_tkebudget 
  ysu_topdown_pblmix         = model_config_rec% ysu_topdown_pblmix 
  shinhong_tke_diag          = model_config_rec% shinhong_tke_diag 
  bl_mynn_tkeadvect          = model_config_rec% bl_mynn_tkeadvect 
  bl_mynn_cloudpdf           = model_config_rec% bl_mynn_cloudpdf 
  bl_mynn_mixlength          = model_config_rec% bl_mynn_mixlength 
  bl_mynn_edmf               = model_config_rec% bl_mynn_edmf 
  bl_mynn_edmf_mom           = model_config_rec% bl_mynn_edmf_mom 
  bl_mynn_edmf_tke           = model_config_rec% bl_mynn_edmf_tke 
  bl_mynn_edmf_part          = model_config_rec% bl_mynn_edmf_part 
  bl_mynn_cloudmix           = model_config_rec% bl_mynn_cloudmix 
  bl_mynn_mixqt              = model_config_rec% bl_mynn_mixqt 
  icloud_bl                  = model_config_rec% icloud_bl 
  mfshconv                   = model_config_rec% mfshconv 
  sf_urban_physics           = model_config_rec% sf_urban_physics 
  bldt                       = model_config_rec% bldt 
  cu_physics                 = model_config_rec% cu_physics 
  shcu_physics               = model_config_rec% shcu_physics 
  cu_diag                    = model_config_rec% cu_diag 
  kf_edrates                 = model_config_rec% kf_edrates 
  kfeta_trigger              = model_config_rec% kfeta_trigger 
  nsas_dx_factor             = model_config_rec% nsas_dx_factor 
  cudt                       = model_config_rec% cudt 
  gsmdt                      = model_config_rec% gsmdt 
  isfflx                     = model_config_rec% isfflx 
  ifsnow                     = model_config_rec% ifsnow 
  icloud                     = model_config_rec% icloud 
  ideal_xland                = model_config_rec% ideal_xland 
  swrad_scat                 = model_config_rec% swrad_scat 
  surface_input_source       = model_config_rec% surface_input_source 
  num_soil_layers            = model_config_rec% num_soil_layers 
  maxpatch                   = model_config_rec% maxpatch 
  num_snow_layers            = model_config_rec% num_snow_layers 
  num_snso_layers            = model_config_rec% num_snso_layers 
  num_urban_layers           = model_config_rec% num_urban_layers 
  num_urban_hi               = model_config_rec% num_urban_hi 
  num_months                 = model_config_rec% num_months 
  sf_surface_mosaic          = model_config_rec% sf_surface_mosaic 
  mosaic_cat                 = model_config_rec% mosaic_cat 
  mosaic_cat_soil            = model_config_rec% mosaic_cat_soil 
  mosaic_lu                  = model_config_rec% mosaic_lu 
  mosaic_soil                = model_config_rec% mosaic_soil 
  maxiens                    = model_config_rec% maxiens 
  maxens                     = model_config_rec% maxens 
  maxens2                    = model_config_rec% maxens2 
  maxens3                    = model_config_rec% maxens3 
  ensdim                     = model_config_rec% ensdim 
  cugd_avedx                 = model_config_rec% cugd_avedx 
  clos_choice                = model_config_rec% clos_choice 
  imomentum                  = model_config_rec% imomentum 
  ishallow                   = model_config_rec% ishallow 
  convtrans_avglen_m         = model_config_rec% convtrans_avglen_m 
  num_land_cat               = model_config_rec% num_land_cat 
  num_soil_cat               = model_config_rec% num_soil_cat 
  mp_zero_out                = model_config_rec% mp_zero_out 
  mp_zero_out_thresh         = model_config_rec% mp_zero_out_thresh 
  seaice_threshold           = model_config_rec% seaice_threshold 
  sst_update                 = model_config_rec% sst_update 
  sst_skin                   = model_config_rec% sst_skin 
  tmn_update                 = model_config_rec% tmn_update 
  usemonalb                  = model_config_rec% usemonalb 
  rdmaxalb                   = model_config_rec% rdmaxalb 
  rdlai2d                    = model_config_rec% rdlai2d 
  ua_phys                    = model_config_rec% ua_phys 
  opt_thcnd                  = model_config_rec% opt_thcnd 
  co2tf                      = model_config_rec% co2tf 
  ra_call_offset             = model_config_rec% ra_call_offset 
  cam_abs_freq_s             = model_config_rec% cam_abs_freq_s 
  levsiz                     = model_config_rec% levsiz 
  paerlev                    = model_config_rec% paerlev 
  cam_abs_dim1               = model_config_rec% cam_abs_dim1 
  cam_abs_dim2               = model_config_rec% cam_abs_dim2 
  lagday                     = model_config_rec% lagday 
  no_src_types               = model_config_rec% no_src_types 
  alevsiz                    = model_config_rec% alevsiz 
  o3input                    = model_config_rec% o3input 
  aer_opt                    = model_config_rec% aer_opt 
  swint_opt                  = model_config_rec% swint_opt 
  aer_type                   = model_config_rec% aer_type 
  aer_aod550_opt             = model_config_rec% aer_aod550_opt 
  aer_angexp_opt             = model_config_rec% aer_angexp_opt 
  aer_ssa_opt                = model_config_rec% aer_ssa_opt 
  aer_asy_opt                = model_config_rec% aer_asy_opt 
  aer_aod550_val             = model_config_rec% aer_aod550_val 
  aer_angexp_val             = model_config_rec% aer_angexp_val 
  aer_ssa_val                = model_config_rec% aer_ssa_val 
  aer_asy_val                = model_config_rec% aer_asy_val 
  cu_rad_feedback            = model_config_rec% cu_rad_feedback 
  shallowcu_forced_ra        = model_config_rec% shallowcu_forced_ra 
  numbins                    = model_config_rec% numbins 
  thbinsize                  = model_config_rec% thbinsize 
  rbinsize                   = model_config_rec% rbinsize 
  mindeepfreq                = model_config_rec% mindeepfreq 
  minshallowfreq             = model_config_rec% minshallowfreq 
  shcu_aerosols_opt          = model_config_rec% shcu_aerosols_opt 
  icloud_cu                  = model_config_rec% icloud_cu 
  pxlsm_smois_init           = model_config_rec% pxlsm_smois_init 
  omlcall                    = model_config_rec% omlcall 
  sf_ocean_physics           = model_config_rec% sf_ocean_physics 
  traj_opt                   = model_config_rec% traj_opt 
  tracercall                 = model_config_rec% tracercall 
  omdt                       = model_config_rec% omdt 
  oml_hml0                   = model_config_rec% oml_hml0 
  oml_gamma                  = model_config_rec% oml_gamma 
  oml_relaxation_time        = model_config_rec% oml_relaxation_time 
  isftcflx                   = model_config_rec% isftcflx 
  iz0tlnd                    = model_config_rec% iz0tlnd 
  shadlen                    = model_config_rec% shadlen 
  slope_rad                  = model_config_rec% slope_rad 
  topo_shading               = model_config_rec% topo_shading 
  topo_wind                  = model_config_rec% topo_wind 
  no_mp_heating              = model_config_rec% no_mp_heating 
  fractional_seaice          = model_config_rec% fractional_seaice 
  seaice_snowdepth_opt       = model_config_rec% seaice_snowdepth_opt 
  seaice_snowdepth_max       = model_config_rec% seaice_snowdepth_max 
  seaice_snowdepth_min       = model_config_rec% seaice_snowdepth_min 
  seaice_albedo_opt          = model_config_rec% seaice_albedo_opt 
  seaice_albedo_default      = model_config_rec% seaice_albedo_default 
  seaice_thickness_opt       = model_config_rec% seaice_thickness_opt 
  seaice_thickness_default   = model_config_rec% seaice_thickness_default 
  tice2tsk_if2cold           = model_config_rec% tice2tsk_if2cold 
  bucket_mm                  = model_config_rec% bucket_mm 
  bucket_j                   = model_config_rec% bucket_j 
  mp_tend_lim                = model_config_rec% mp_tend_lim 
  prec_acc_dt                = model_config_rec% prec_acc_dt 
  prec_acc_opt               = model_config_rec% prec_acc_opt 
  bucketr_opt                = model_config_rec% bucketr_opt 
  process_time_series        = model_config_rec% process_time_series 
  grav_settling              = model_config_rec% grav_settling 
  sas_pgcon                  = model_config_rec% sas_pgcon 
  scalar_pblmix              = model_config_rec% scalar_pblmix 
  tracer_pblmix              = model_config_rec% tracer_pblmix 
  use_aero_icbc              = model_config_rec% use_aero_icbc 
  use_rap_aero_icbc          = model_config_rec% use_rap_aero_icbc 
  use_mp_re                  = model_config_rec% use_mp_re 
  ccn_conc                   = model_config_rec% ccn_conc 
  hail_opt                   = model_config_rec% hail_opt 
  dveg                       = model_config_rec% dveg 
  opt_crs                    = model_config_rec% opt_crs 
  opt_btr                    = model_config_rec% opt_btr 
  opt_run                    = model_config_rec% opt_run 
  opt_sfc                    = model_config_rec% opt_sfc 
  opt_frz                    = model_config_rec% opt_frz 
  opt_inf                    = model_config_rec% opt_inf 
  opt_rad                    = model_config_rec% opt_rad 
  opt_alb                    = model_config_rec% opt_alb 
  opt_snf                    = model_config_rec% opt_snf 
  opt_tbot                   = model_config_rec% opt_tbot 
  opt_stc                    = model_config_rec% opt_stc 
  opt_gla                    = model_config_rec% opt_gla 
  opt_rsf                    = model_config_rec% opt_rsf 
  wtddt                      = model_config_rec% wtddt 
  wrf_hydro                  = model_config_rec% wrf_hydro 
  fgdt                       = model_config_rec% fgdt 
  fgdtzero                   = model_config_rec% fgdtzero 
  grid_fdda                  = model_config_rec% grid_fdda 
  grid_sfdda                 = model_config_rec% grid_sfdda 
  if_no_pbl_nudging_uv       = model_config_rec% if_no_pbl_nudging_uv 
  if_no_pbl_nudging_t        = model_config_rec% if_no_pbl_nudging_t 
  if_no_pbl_nudging_ph       = model_config_rec% if_no_pbl_nudging_ph 
  if_no_pbl_nudging_q        = model_config_rec% if_no_pbl_nudging_q 
  if_zfac_uv                 = model_config_rec% if_zfac_uv 
  k_zfac_uv                  = model_config_rec% k_zfac_uv 
  if_zfac_t                  = model_config_rec% if_zfac_t 
  k_zfac_t                   = model_config_rec% k_zfac_t 
  if_zfac_ph                 = model_config_rec% if_zfac_ph 
  k_zfac_ph                  = model_config_rec% k_zfac_ph 
  if_zfac_q                  = model_config_rec% if_zfac_q 
  k_zfac_q                   = model_config_rec% k_zfac_q 
  dk_zfac_uv                 = model_config_rec% dk_zfac_uv 
  dk_zfac_t                  = model_config_rec% dk_zfac_t 
  dk_zfac_ph                 = model_config_rec% dk_zfac_ph 
  guv                        = model_config_rec% guv 
  guv_sfc                    = model_config_rec% guv_sfc 
  gt                         = model_config_rec% gt 
  gt_sfc                     = model_config_rec% gt_sfc 
  gq                         = model_config_rec% gq 
  gq_sfc                     = model_config_rec% gq_sfc 
  gph                        = model_config_rec% gph 
  dtramp_min                 = model_config_rec% dtramp_min 
  if_ramping                 = model_config_rec% if_ramping 
  rinblw                     = model_config_rec% rinblw 
  xwavenum                   = model_config_rec% xwavenum 
  ywavenum                   = model_config_rec% ywavenum 
  pxlsm_soil_nudge           = model_config_rec% pxlsm_soil_nudge 
  fasdas                     = model_config_rec% fasdas 
  obs_nudge_opt              = model_config_rec% obs_nudge_opt 
  max_obs                    = model_config_rec% max_obs 
  fdda_start                 = model_config_rec% fdda_start 
  fdda_end                   = model_config_rec% fdda_end 
  obs_nudge_wind             = model_config_rec% obs_nudge_wind 
  obs_coef_wind              = model_config_rec% obs_coef_wind 
  obs_nudge_temp             = model_config_rec% obs_nudge_temp 
  obs_coef_temp              = model_config_rec% obs_coef_temp 
  obs_nudge_mois             = model_config_rec% obs_nudge_mois 
  obs_coef_mois              = model_config_rec% obs_coef_mois 
  obs_nudge_pstr             = model_config_rec% obs_nudge_pstr 
  obs_coef_pstr              = model_config_rec% obs_coef_pstr 
  obs_no_pbl_nudge_uv        = model_config_rec% obs_no_pbl_nudge_uv 
  obs_no_pbl_nudge_t         = model_config_rec% obs_no_pbl_nudge_t 
  obs_no_pbl_nudge_q         = model_config_rec% obs_no_pbl_nudge_q 
  obs_sfc_scheme_horiz       = model_config_rec% obs_sfc_scheme_horiz 
  obs_sfc_scheme_vert        = model_config_rec% obs_sfc_scheme_vert 
  obs_max_sndng_gap          = model_config_rec% obs_max_sndng_gap 
  obs_nudgezfullr1_uv        = model_config_rec% obs_nudgezfullr1_uv 
  obs_nudgezrampr1_uv        = model_config_rec% obs_nudgezrampr1_uv 
  obs_nudgezfullr2_uv        = model_config_rec% obs_nudgezfullr2_uv 
  obs_nudgezrampr2_uv        = model_config_rec% obs_nudgezrampr2_uv 
  obs_nudgezfullr4_uv        = model_config_rec% obs_nudgezfullr4_uv 
  obs_nudgezrampr4_uv        = model_config_rec% obs_nudgezrampr4_uv 
  obs_nudgezfullr1_t         = model_config_rec% obs_nudgezfullr1_t 
  obs_nudgezrampr1_t         = model_config_rec% obs_nudgezrampr1_t 
  obs_nudgezfullr2_t         = model_config_rec% obs_nudgezfullr2_t 
  obs_nudgezrampr2_t         = model_config_rec% obs_nudgezrampr2_t 
  obs_nudgezfullr4_t         = model_config_rec% obs_nudgezfullr4_t 
  obs_nudgezrampr4_t         = model_config_rec% obs_nudgezrampr4_t 
  obs_nudgezfullr1_q         = model_config_rec% obs_nudgezfullr1_q 
  obs_nudgezrampr1_q         = model_config_rec% obs_nudgezrampr1_q 
  obs_nudgezfullr2_q         = model_config_rec% obs_nudgezfullr2_q 
  obs_nudgezrampr2_q         = model_config_rec% obs_nudgezrampr2_q 
  obs_nudgezfullr4_q         = model_config_rec% obs_nudgezfullr4_q 
  obs_nudgezrampr4_q         = model_config_rec% obs_nudgezrampr4_q 
  obs_nudgezfullmin          = model_config_rec% obs_nudgezfullmin 
  obs_nudgezrampmin          = model_config_rec% obs_nudgezrampmin 
  obs_nudgezmax              = model_config_rec% obs_nudgezmax 
  obs_sfcfact                = model_config_rec% obs_sfcfact 
  obs_sfcfacr                = model_config_rec% obs_sfcfacr 
  obs_dpsmx                  = model_config_rec% obs_dpsmx 
  obs_rinxy                  = model_config_rec% obs_rinxy 
  obs_rinsig                 = model_config_rec% obs_rinsig 
  obs_twindo                 = model_config_rec% obs_twindo 
  obs_npfi                   = model_config_rec% obs_npfi 
  obs_ionf                   = model_config_rec% obs_ionf 
  obs_idynin                 = model_config_rec% obs_idynin 
  obs_dtramp                 = model_config_rec% obs_dtramp 
  obs_prt_max                = model_config_rec% obs_prt_max 
  obs_prt_freq               = model_config_rec% obs_prt_freq 
  obs_ipf_in4dob             = model_config_rec% obs_ipf_in4dob 
  obs_ipf_errob              = model_config_rec% obs_ipf_errob 
  obs_ipf_nudob              = model_config_rec% obs_ipf_nudob 
  obs_ipf_init               = model_config_rec% obs_ipf_init 
  obs_scl_neg_qv_innov       = model_config_rec% obs_scl_neg_qv_innov 
  scm_force                  = model_config_rec% scm_force 
  scm_force_dx               = model_config_rec% scm_force_dx 
  num_force_layers           = model_config_rec% num_force_layers 
  scm_lu_index               = model_config_rec% scm_lu_index 
  scm_isltyp                 = model_config_rec% scm_isltyp 
  scm_vegfra                 = model_config_rec% scm_vegfra 
  scm_canwat                 = model_config_rec% scm_canwat 
  scm_lat                    = model_config_rec% scm_lat 
  scm_lon                    = model_config_rec% scm_lon 
  scm_th_t_tend              = model_config_rec% scm_th_t_tend 
  scm_qv_t_tend              = model_config_rec% scm_qv_t_tend 
  scm_th_adv                 = model_config_rec% scm_th_adv 
  scm_wind_adv               = model_config_rec% scm_wind_adv 
  scm_qv_adv                 = model_config_rec% scm_qv_adv 
  scm_ql_adv                 = model_config_rec% scm_ql_adv 
  scm_vert_adv               = model_config_rec% scm_vert_adv 
  num_force_soil_layers      = model_config_rec% num_force_soil_layers 
  scm_soilt_force            = model_config_rec% scm_soilt_force 
  scm_soilq_force            = model_config_rec% scm_soilq_force 
  scm_force_th_largescale    = model_config_rec% scm_force_th_largescale 
  scm_force_qv_largescale    = model_config_rec% scm_force_qv_largescale 
  scm_force_ql_largescale    = model_config_rec% scm_force_ql_largescale 
  scm_force_wind_largescale  = model_config_rec% scm_force_wind_largescale 
  scm_force_skintemp         = model_config_rec% scm_force_skintemp 
  scm_force_flux             = model_config_rec% scm_force_flux 
  dyn_opt                    = model_config_rec% dyn_opt 
  rk_ord                     = model_config_rec% rk_ord 
  w_damping                  = model_config_rec% w_damping 
  diff_opt                   = model_config_rec% diff_opt 
  diff_opt_dfi               = model_config_rec% diff_opt_dfi 
  km_opt                     = model_config_rec% km_opt 
  km_opt_dfi                 = model_config_rec% km_opt_dfi 
  damp_opt                   = model_config_rec% damp_opt 
  rad_nudge                  = model_config_rec% rad_nudge 
  gwd_opt                    = model_config_rec% gwd_opt 
  zdamp                      = model_config_rec% zdamp 
  dampcoef                   = model_config_rec% dampcoef 
  khdif                      = model_config_rec% khdif 
  kvdif                      = model_config_rec% kvdif 
  diff_6th_factor            = model_config_rec% diff_6th_factor 
  diff_6th_opt               = model_config_rec% diff_6th_opt 
  use_theta_m                = model_config_rec% use_theta_m 
  use_q_diabatic             = model_config_rec% use_q_diabatic 
  c_s                        = model_config_rec% c_s 
  c_k                        = model_config_rec% c_k 
  smdiv                      = model_config_rec% smdiv 
  emdiv                      = model_config_rec% emdiv 
  epssm                      = model_config_rec% epssm 
  non_hydrostatic            = model_config_rec% non_hydrostatic 
  use_input_w                = model_config_rec% use_input_w 
  time_step_sound            = model_config_rec% time_step_sound 
  h_mom_adv_order            = model_config_rec% h_mom_adv_order 
  v_mom_adv_order            = model_config_rec% v_mom_adv_order 
  h_sca_adv_order            = model_config_rec% h_sca_adv_order 
  v_sca_adv_order            = model_config_rec% v_sca_adv_order 
  momentum_adv_opt           = model_config_rec% momentum_adv_opt 
  moist_adv_opt              = model_config_rec% moist_adv_opt 
  moist_adv_dfi_opt          = model_config_rec% moist_adv_dfi_opt 
  chem_adv_opt               = model_config_rec% chem_adv_opt 
  tracer_adv_opt             = model_config_rec% tracer_adv_opt 
  scalar_adv_opt             = model_config_rec% scalar_adv_opt 
  tke_adv_opt                = model_config_rec% tke_adv_opt 
  top_radiation              = model_config_rec% top_radiation 
  mix_isotropic              = model_config_rec% mix_isotropic 
  mix_upper_bound            = model_config_rec% mix_upper_bound 
  top_lid                    = model_config_rec% top_lid 
  tke_upper_bound            = model_config_rec% tke_upper_bound 
  tke_drag_coefficient       = model_config_rec% tke_drag_coefficient 
  tke_heat_flux              = model_config_rec% tke_heat_flux 
  pert_coriolis              = model_config_rec% pert_coriolis 
  coriolis2d                 = model_config_rec% coriolis2d 
  mix_full_fields            = model_config_rec% mix_full_fields 
  base_pres                  = model_config_rec% base_pres 
  base_temp                  = model_config_rec% base_temp 
  base_lapse                 = model_config_rec% base_lapse 
  iso_temp                   = model_config_rec% iso_temp 
  base_pres_strat            = model_config_rec% base_pres_strat 
  base_lapse_strat           = model_config_rec% base_lapse_strat 
  use_baseparam_fr_nml       = model_config_rec% use_baseparam_fr_nml 
  fft_filter_lat             = model_config_rec% fft_filter_lat 
  coupled_filtering          = model_config_rec% coupled_filtering 
  pos_def                    = model_config_rec% pos_def 
  swap_pole_with_next_j      = model_config_rec% swap_pole_with_next_j 
  actual_distance_average    = model_config_rec% actual_distance_average 
  rotated_pole               = model_config_rec% rotated_pole 
  do_coriolis                = model_config_rec% do_coriolis 
  do_curvature               = model_config_rec% do_curvature 
  do_gradp                   = model_config_rec% do_gradp 
  tracer_opt                 = model_config_rec% tracer_opt 
  tenddiag                   = model_config_rec% tenddiag 
  spec_bdy_width             = model_config_rec% spec_bdy_width 
  spec_zone                  = model_config_rec% spec_zone 
  relax_zone                 = model_config_rec% relax_zone 
  specified                  = model_config_rec% specified 
  constant_bc                = model_config_rec% constant_bc 
  periodic_x                 = model_config_rec% periodic_x 
  symmetric_xs               = model_config_rec% symmetric_xs 
  symmetric_xe               = model_config_rec% symmetric_xe 
  open_xs                    = model_config_rec% open_xs 
  open_xe                    = model_config_rec% open_xe 
  periodic_y                 = model_config_rec% periodic_y 
  symmetric_ys               = model_config_rec% symmetric_ys 
  symmetric_ye               = model_config_rec% symmetric_ye 
  open_ys                    = model_config_rec% open_ys 
  open_ye                    = model_config_rec% open_ye 
  polar                      = model_config_rec% polar 
  nested                     = model_config_rec% nested 
  spec_exp                   = model_config_rec% spec_exp 
  spec_bdy_final_mu          = model_config_rec% spec_bdy_final_mu 
  real_data_init_type        = model_config_rec% real_data_init_type 
  have_bcs_moist             = model_config_rec% have_bcs_moist 
  have_bcs_scalar            = model_config_rec% have_bcs_scalar 
  background_proc_id         = model_config_rec% background_proc_id 
  forecast_proc_id           = model_config_rec% forecast_proc_id 
  production_status          = model_config_rec% production_status 
  compression                = model_config_rec% compression 
  nobs_ndg_vars              = model_config_rec% nobs_ndg_vars 
  nobs_err_flds              = model_config_rec% nobs_err_flds 
  cen_lat                    = model_config_rec% cen_lat 
  cen_lon                    = model_config_rec% cen_lon 
  truelat1                   = model_config_rec% truelat1 
  truelat2                   = model_config_rec% truelat2 
  moad_cen_lat               = model_config_rec% moad_cen_lat 
  stand_lon                  = model_config_rec% stand_lon 
  pole_lat                   = model_config_rec% pole_lat 
  pole_lon                   = model_config_rec% pole_lon 
  flag_metgrid               = model_config_rec% flag_metgrid 
  flag_snow                  = model_config_rec% flag_snow 
  flag_psfc                  = model_config_rec% flag_psfc 
  flag_sm000010              = model_config_rec% flag_sm000010 
  flag_sm010040              = model_config_rec% flag_sm010040 
  flag_sm040100              = model_config_rec% flag_sm040100 
  flag_sm100200              = model_config_rec% flag_sm100200 
  flag_st000010              = model_config_rec% flag_st000010 
  flag_st010040              = model_config_rec% flag_st010040 
  flag_st040100              = model_config_rec% flag_st040100 
  flag_st100200              = model_config_rec% flag_st100200 
  flag_soil_layers           = model_config_rec% flag_soil_layers 
  flag_slp                   = model_config_rec% flag_slp 
  flag_soilhgt               = model_config_rec% flag_soilhgt 
  flag_mf_xy                 = model_config_rec% flag_mf_xy 
  flag_um_soil               = model_config_rec% flag_um_soil 
  bdyfrq                     = model_config_rec% bdyfrq 
  mminlu                     = model_config_rec% mminlu 
  iswater                    = model_config_rec% iswater 
  islake                     = model_config_rec% islake 
  isice                      = model_config_rec% isice 
  isurban                    = model_config_rec% isurban 
  isoilwater                 = model_config_rec% isoilwater 
  map_proj                   = model_config_rec% map_proj 
  use_wps_input              = model_config_rec% use_wps_input 
  dfi_stage                  = model_config_rec% dfi_stage 
  mp_physics_dfi             = model_config_rec% mp_physics_dfi 
  bl_pbl_physics_dfi         = model_config_rec% bl_pbl_physics_dfi 
  windfarm_opt               = model_config_rec% windfarm_opt 
  windfarm_ij                = model_config_rec% windfarm_ij 
  lightning_option           = model_config_rec% lightning_option 
  lightning_dt               = model_config_rec% lightning_dt 
  lightning_start_seconds    = model_config_rec% lightning_start_seconds 
  flashrate_factor           = model_config_rec% flashrate_factor 
  iccg_method                = model_config_rec% iccg_method 
  iccg_prescribed_num        = model_config_rec% iccg_prescribed_num 
  iccg_prescribed_den        = model_config_rec% iccg_prescribed_den 
  cellcount_method           = model_config_rec% cellcount_method 
  cldtop_adjustment          = model_config_rec% cldtop_adjustment 
  sf_lake_physics            = model_config_rec% sf_lake_physics 
  auxinput1_inname           = model_config_rec% auxinput1_inname 
  io_form_auxinput1          = model_config_rec% io_form_auxinput1 
  override_restart_timers    = model_config_rec% override_restart_timers 
  auxhist1_inname            = model_config_rec% auxhist1_inname 
  auxhist1_outname           = model_config_rec% auxhist1_outname 
  auxhist1_interval_y        = model_config_rec% auxhist1_interval_y 
  auxhist1_interval_d        = model_config_rec% auxhist1_interval_d 
  auxhist1_interval_h        = model_config_rec% auxhist1_interval_h 
  auxhist1_interval_m        = model_config_rec% auxhist1_interval_m 
  auxhist1_interval_s        = model_config_rec% auxhist1_interval_s 
  auxhist1_interval          = model_config_rec% auxhist1_interval 
  auxhist1_begin_y           = model_config_rec% auxhist1_begin_y 
  auxhist1_begin_d           = model_config_rec% auxhist1_begin_d 
  auxhist1_begin_h           = model_config_rec% auxhist1_begin_h 
  auxhist1_begin_m           = model_config_rec% auxhist1_begin_m 
  auxhist1_begin_s           = model_config_rec% auxhist1_begin_s 
  auxhist1_begin             = model_config_rec% auxhist1_begin 
  auxhist1_end_y             = model_config_rec% auxhist1_end_y 
  auxhist1_end_d             = model_config_rec% auxhist1_end_d 
  auxhist1_end_h             = model_config_rec% auxhist1_end_h 
  auxhist1_end_m             = model_config_rec% auxhist1_end_m 
  auxhist1_end_s             = model_config_rec% auxhist1_end_s 
  auxhist1_end               = model_config_rec% auxhist1_end 
  io_form_auxhist1           = model_config_rec% io_form_auxhist1 
  frames_per_auxhist1        = model_config_rec% frames_per_auxhist1 
  auxhist2_inname            = model_config_rec% auxhist2_inname 
  auxhist2_outname           = model_config_rec% auxhist2_outname 
  auxhist2_interval_y        = model_config_rec% auxhist2_interval_y 
  auxhist2_interval_d        = model_config_rec% auxhist2_interval_d 
  auxhist2_interval_h        = model_config_rec% auxhist2_interval_h 
  auxhist2_interval_m        = model_config_rec% auxhist2_interval_m 
  auxhist2_interval_s        = model_config_rec% auxhist2_interval_s 
  auxhist2_interval          = model_config_rec% auxhist2_interval 
  auxhist2_begin_y           = model_config_rec% auxhist2_begin_y 
  auxhist2_begin_d           = model_config_rec% auxhist2_begin_d 
  auxhist2_begin_h           = model_config_rec% auxhist2_begin_h 
  auxhist2_begin_m           = model_config_rec% auxhist2_begin_m 
  auxhist2_begin_s           = model_config_rec% auxhist2_begin_s 
  auxhist2_begin             = model_config_rec% auxhist2_begin 
  auxhist2_end_y             = model_config_rec% auxhist2_end_y 
  auxhist2_end_d             = model_config_rec% auxhist2_end_d 
  auxhist2_end_h             = model_config_rec% auxhist2_end_h 
  auxhist2_end_m             = model_config_rec% auxhist2_end_m 
  auxhist2_end_s             = model_config_rec% auxhist2_end_s 
  auxhist2_end               = model_config_rec% auxhist2_end 
  io_form_auxhist2           = model_config_rec% io_form_auxhist2 
  frames_per_auxhist2        = model_config_rec% frames_per_auxhist2 
  auxhist3_inname            = model_config_rec% auxhist3_inname 
  auxhist3_outname           = model_config_rec% auxhist3_outname 
  auxhist3_interval_y        = model_config_rec% auxhist3_interval_y 
  auxhist3_interval_d        = model_config_rec% auxhist3_interval_d 
  auxhist3_interval_h        = model_config_rec% auxhist3_interval_h 
  auxhist3_interval_m        = model_config_rec% auxhist3_interval_m 
  auxhist3_interval_s        = model_config_rec% auxhist3_interval_s 
  auxhist3_interval          = model_config_rec% auxhist3_interval 
  auxhist3_begin_y           = model_config_rec% auxhist3_begin_y 
  auxhist3_begin_d           = model_config_rec% auxhist3_begin_d 
  auxhist3_begin_h           = model_config_rec% auxhist3_begin_h 
  auxhist3_begin_m           = model_config_rec% auxhist3_begin_m 
  auxhist3_begin_s           = model_config_rec% auxhist3_begin_s 
  auxhist3_begin             = model_config_rec% auxhist3_begin 
  auxhist3_end_y             = model_config_rec% auxhist3_end_y 
  auxhist3_end_d             = model_config_rec% auxhist3_end_d 
  auxhist3_end_h             = model_config_rec% auxhist3_end_h 
  auxhist3_end_m             = model_config_rec% auxhist3_end_m 
  auxhist3_end_s             = model_config_rec% auxhist3_end_s 
  auxhist3_end               = model_config_rec% auxhist3_end 
  io_form_auxhist3           = model_config_rec% io_form_auxhist3 
  frames_per_auxhist3        = model_config_rec% frames_per_auxhist3 
  auxhist4_inname            = model_config_rec% auxhist4_inname 
  auxhist4_outname           = model_config_rec% auxhist4_outname 
  auxhist4_interval_y        = model_config_rec% auxhist4_interval_y 
  auxhist4_interval_d        = model_config_rec% auxhist4_interval_d 
  auxhist4_interval_h        = model_config_rec% auxhist4_interval_h 
  auxhist4_interval_m        = model_config_rec% auxhist4_interval_m 
  auxhist4_interval_s        = model_config_rec% auxhist4_interval_s 
  auxhist4_interval          = model_config_rec% auxhist4_interval 
  auxhist4_begin_y           = model_config_rec% auxhist4_begin_y 
  auxhist4_begin_d           = model_config_rec% auxhist4_begin_d 
  auxhist4_begin_h           = model_config_rec% auxhist4_begin_h 
  auxhist4_begin_m           = model_config_rec% auxhist4_begin_m 
  auxhist4_begin_s           = model_config_rec% auxhist4_begin_s 
  auxhist4_begin             = model_config_rec% auxhist4_begin 
  auxhist4_end_y             = model_config_rec% auxhist4_end_y 
  auxhist4_end_d             = model_config_rec% auxhist4_end_d 
  auxhist4_end_h             = model_config_rec% auxhist4_end_h 
  auxhist4_end_m             = model_config_rec% auxhist4_end_m 
  auxhist4_end_s             = model_config_rec% auxhist4_end_s 
  auxhist4_end               = model_config_rec% auxhist4_end 
  io_form_auxhist4           = model_config_rec% io_form_auxhist4 
  frames_per_auxhist4        = model_config_rec% frames_per_auxhist4 
  auxhist5_inname            = model_config_rec% auxhist5_inname 
  auxhist5_outname           = model_config_rec% auxhist5_outname 
  auxhist5_interval_y        = model_config_rec% auxhist5_interval_y 
  auxhist5_interval_d        = model_config_rec% auxhist5_interval_d 
  auxhist5_interval_h        = model_config_rec% auxhist5_interval_h 
  auxhist5_interval_m        = model_config_rec% auxhist5_interval_m 
  auxhist5_interval_s        = model_config_rec% auxhist5_interval_s 
  auxhist5_interval          = model_config_rec% auxhist5_interval 
  auxhist5_begin_y           = model_config_rec% auxhist5_begin_y 
  auxhist5_begin_d           = model_config_rec% auxhist5_begin_d 
  auxhist5_begin_h           = model_config_rec% auxhist5_begin_h 
  auxhist5_begin_m           = model_config_rec% auxhist5_begin_m 
  auxhist5_begin_s           = model_config_rec% auxhist5_begin_s 
  auxhist5_begin             = model_config_rec% auxhist5_begin 
  auxhist5_end_y             = model_config_rec% auxhist5_end_y 
  auxhist5_end_d             = model_config_rec% auxhist5_end_d 
  auxhist5_end_h             = model_config_rec% auxhist5_end_h 
  auxhist5_end_m             = model_config_rec% auxhist5_end_m 
  auxhist5_end_s             = model_config_rec% auxhist5_end_s 
  auxhist5_end               = model_config_rec% auxhist5_end 
  io_form_auxhist5           = model_config_rec% io_form_auxhist5 
  frames_per_auxhist5        = model_config_rec% frames_per_auxhist5 
  auxhist6_inname            = model_config_rec% auxhist6_inname 
  auxhist6_outname           = model_config_rec% auxhist6_outname 
  auxhist6_interval_y        = model_config_rec% auxhist6_interval_y 
  auxhist6_interval_d        = model_config_rec% auxhist6_interval_d 
  auxhist6_interval_h        = model_config_rec% auxhist6_interval_h 
  auxhist6_interval_m        = model_config_rec% auxhist6_interval_m 
  auxhist6_interval_s        = model_config_rec% auxhist6_interval_s 
  auxhist6_interval          = model_config_rec% auxhist6_interval 
  auxhist6_begin_y           = model_config_rec% auxhist6_begin_y 
  auxhist6_begin_d           = model_config_rec% auxhist6_begin_d 
  auxhist6_begin_h           = model_config_rec% auxhist6_begin_h 
  auxhist6_begin_m           = model_config_rec% auxhist6_begin_m 
  auxhist6_begin_s           = model_config_rec% auxhist6_begin_s 
  auxhist6_begin             = model_config_rec% auxhist6_begin 
  auxhist6_end_y             = model_config_rec% auxhist6_end_y 
  auxhist6_end_d             = model_config_rec% auxhist6_end_d 
  auxhist6_end_h             = model_config_rec% auxhist6_end_h 
  auxhist6_end_m             = model_config_rec% auxhist6_end_m 
  auxhist6_end_s             = model_config_rec% auxhist6_end_s 
  auxhist6_end               = model_config_rec% auxhist6_end 
  io_form_auxhist6           = model_config_rec% io_form_auxhist6 
  frames_per_auxhist6        = model_config_rec% frames_per_auxhist6 
  auxhist7_inname            = model_config_rec% auxhist7_inname 
  auxhist7_outname           = model_config_rec% auxhist7_outname 
  auxhist7_interval_y        = model_config_rec% auxhist7_interval_y 
  auxhist7_interval_d        = model_config_rec% auxhist7_interval_d 
  auxhist7_interval_h        = model_config_rec% auxhist7_interval_h 
  auxhist7_interval_m        = model_config_rec% auxhist7_interval_m 
  auxhist7_interval_s        = model_config_rec% auxhist7_interval_s 
  auxhist7_interval          = model_config_rec% auxhist7_interval 
  auxhist7_begin_y           = model_config_rec% auxhist7_begin_y 
  auxhist7_begin_d           = model_config_rec% auxhist7_begin_d 
  auxhist7_begin_h           = model_config_rec% auxhist7_begin_h 
  auxhist7_begin_m           = model_config_rec% auxhist7_begin_m 
  auxhist7_begin_s           = model_config_rec% auxhist7_begin_s 
  auxhist7_begin             = model_config_rec% auxhist7_begin 
  auxhist7_end_y             = model_config_rec% auxhist7_end_y 
  auxhist7_end_d             = model_config_rec% auxhist7_end_d 
  auxhist7_end_h             = model_config_rec% auxhist7_end_h 
  auxhist7_end_m             = model_config_rec% auxhist7_end_m 
  auxhist7_end_s             = model_config_rec% auxhist7_end_s 
  auxhist7_end               = model_config_rec% auxhist7_end 
  io_form_auxhist7           = model_config_rec% io_form_auxhist7 
  frames_per_auxhist7        = model_config_rec% frames_per_auxhist7 
  auxhist8_inname            = model_config_rec% auxhist8_inname 
  auxhist8_outname           = model_config_rec% auxhist8_outname 
  auxhist8_interval_y        = model_config_rec% auxhist8_interval_y 
  auxhist8_interval_d        = model_config_rec% auxhist8_interval_d 
  auxhist8_interval_h        = model_config_rec% auxhist8_interval_h 
  auxhist8_interval_m        = model_config_rec% auxhist8_interval_m 
  auxhist8_interval_s        = model_config_rec% auxhist8_interval_s 
  auxhist8_interval          = model_config_rec% auxhist8_interval 
  auxhist8_begin_y           = model_config_rec% auxhist8_begin_y 
  auxhist8_begin_d           = model_config_rec% auxhist8_begin_d 
  auxhist8_begin_h           = model_config_rec% auxhist8_begin_h 
  auxhist8_begin_m           = model_config_rec% auxhist8_begin_m 
  auxhist8_begin_s           = model_config_rec% auxhist8_begin_s 
  auxhist8_begin             = model_config_rec% auxhist8_begin 
  auxhist8_end_y             = model_config_rec% auxhist8_end_y 
  auxhist8_end_d             = model_config_rec% auxhist8_end_d 
  auxhist8_end_h             = model_config_rec% auxhist8_end_h 
  auxhist8_end_m             = model_config_rec% auxhist8_end_m 
  auxhist8_end_s             = model_config_rec% auxhist8_end_s 
  auxhist8_end               = model_config_rec% auxhist8_end 
  io_form_auxhist8           = model_config_rec% io_form_auxhist8 
  frames_per_auxhist8        = model_config_rec% frames_per_auxhist8 
  auxhist9_inname            = model_config_rec% auxhist9_inname 
  auxhist9_outname           = model_config_rec% auxhist9_outname 
  auxhist9_interval_y        = model_config_rec% auxhist9_interval_y 
  auxhist9_interval_d        = model_config_rec% auxhist9_interval_d 
  auxhist9_interval_h        = model_config_rec% auxhist9_interval_h 
  auxhist9_interval_m        = model_config_rec% auxhist9_interval_m 
  auxhist9_interval_s        = model_config_rec% auxhist9_interval_s 
  auxhist9_interval          = model_config_rec% auxhist9_interval 
  auxhist9_begin_y           = model_config_rec% auxhist9_begin_y 
  auxhist9_begin_d           = model_config_rec% auxhist9_begin_d 
  auxhist9_begin_h           = model_config_rec% auxhist9_begin_h 
  auxhist9_begin_m           = model_config_rec% auxhist9_begin_m 
  auxhist9_begin_s           = model_config_rec% auxhist9_begin_s 
  auxhist9_begin             = model_config_rec% auxhist9_begin 
  auxhist9_end_y             = model_config_rec% auxhist9_end_y 
  auxhist9_end_d             = model_config_rec% auxhist9_end_d 
  auxhist9_end_h             = model_config_rec% auxhist9_end_h 
  auxhist9_end_m             = model_config_rec% auxhist9_end_m 
  auxhist9_end_s             = model_config_rec% auxhist9_end_s 
  auxhist9_end               = model_config_rec% auxhist9_end 
  io_form_auxhist9           = model_config_rec% io_form_auxhist9 
  frames_per_auxhist9        = model_config_rec% frames_per_auxhist9 
  auxhist10_inname           = model_config_rec% auxhist10_inname 
  auxhist10_outname          = model_config_rec% auxhist10_outname 
  auxhist10_interval_y       = model_config_rec% auxhist10_interval_y 
  auxhist10_interval_d       = model_config_rec% auxhist10_interval_d 
  auxhist10_interval_h       = model_config_rec% auxhist10_interval_h 
  auxhist10_interval_m       = model_config_rec% auxhist10_interval_m 
  auxhist10_interval_s       = model_config_rec% auxhist10_interval_s 
  auxhist10_interval         = model_config_rec% auxhist10_interval 
  auxhist10_begin_y          = model_config_rec% auxhist10_begin_y 
  auxhist10_begin_d          = model_config_rec% auxhist10_begin_d 
  auxhist10_begin_h          = model_config_rec% auxhist10_begin_h 
  auxhist10_begin_m          = model_config_rec% auxhist10_begin_m 
  auxhist10_begin_s          = model_config_rec% auxhist10_begin_s 
  auxhist10_begin            = model_config_rec% auxhist10_begin 
  auxhist10_end_y            = model_config_rec% auxhist10_end_y 
  auxhist10_end_d            = model_config_rec% auxhist10_end_d 
  auxhist10_end_h            = model_config_rec% auxhist10_end_h 
  auxhist10_end_m            = model_config_rec% auxhist10_end_m 
  auxhist10_end_s            = model_config_rec% auxhist10_end_s 
  auxhist10_end              = model_config_rec% auxhist10_end 
  io_form_auxhist10          = model_config_rec% io_form_auxhist10 
  frames_per_auxhist10       = model_config_rec% frames_per_auxhist10 
  auxhist11_inname           = model_config_rec% auxhist11_inname 
  auxhist11_outname          = model_config_rec% auxhist11_outname 
  auxhist11_interval_y       = model_config_rec% auxhist11_interval_y 
  auxhist11_interval_d       = model_config_rec% auxhist11_interval_d 
  auxhist11_interval_h       = model_config_rec% auxhist11_interval_h 
  auxhist11_interval_m       = model_config_rec% auxhist11_interval_m 
  auxhist11_interval_s       = model_config_rec% auxhist11_interval_s 
  auxhist11_interval         = model_config_rec% auxhist11_interval 
  auxhist11_begin_y          = model_config_rec% auxhist11_begin_y 
  auxhist11_begin_d          = model_config_rec% auxhist11_begin_d 
  auxhist11_begin_h          = model_config_rec% auxhist11_begin_h 
  auxhist11_begin_m          = model_config_rec% auxhist11_begin_m 
  auxhist11_begin_s          = model_config_rec% auxhist11_begin_s 
  auxhist11_begin            = model_config_rec% auxhist11_begin 
  auxhist11_end_y            = model_config_rec% auxhist11_end_y 
  auxhist11_end_d            = model_config_rec% auxhist11_end_d 
  auxhist11_end_h            = model_config_rec% auxhist11_end_h 
  auxhist11_end_m            = model_config_rec% auxhist11_end_m 
  auxhist11_end_s            = model_config_rec% auxhist11_end_s 
  auxhist11_end              = model_config_rec% auxhist11_end 
  io_form_auxhist11          = model_config_rec% io_form_auxhist11 
  frames_per_auxhist11       = model_config_rec% frames_per_auxhist11 
  auxhist12_inname           = model_config_rec% auxhist12_inname 
  auxhist12_outname          = model_config_rec% auxhist12_outname 
  auxhist12_interval_y       = model_config_rec% auxhist12_interval_y 
  auxhist12_interval_d       = model_config_rec% auxhist12_interval_d 
  auxhist12_interval_h       = model_config_rec% auxhist12_interval_h 
  auxhist12_interval_m       = model_config_rec% auxhist12_interval_m 
  auxhist12_interval_s       = model_config_rec% auxhist12_interval_s 
  auxhist12_interval         = model_config_rec% auxhist12_interval 
  auxhist12_begin_y          = model_config_rec% auxhist12_begin_y 
  auxhist12_begin_d          = model_config_rec% auxhist12_begin_d 
  auxhist12_begin_h          = model_config_rec% auxhist12_begin_h 
  auxhist12_begin_m          = model_config_rec% auxhist12_begin_m 
  auxhist12_begin_s          = model_config_rec% auxhist12_begin_s 
  auxhist12_begin            = model_config_rec% auxhist12_begin 
  auxhist12_end_y            = model_config_rec% auxhist12_end_y 
  auxhist12_end_d            = model_config_rec% auxhist12_end_d 
  auxhist12_end_h            = model_config_rec% auxhist12_end_h 
  auxhist12_end_m            = model_config_rec% auxhist12_end_m 
  auxhist12_end_s            = model_config_rec% auxhist12_end_s 
  auxhist12_end              = model_config_rec% auxhist12_end 
  io_form_auxhist12          = model_config_rec% io_form_auxhist12 
  frames_per_auxhist12       = model_config_rec% frames_per_auxhist12 
  auxhist13_inname           = model_config_rec% auxhist13_inname 
  auxhist13_outname          = model_config_rec% auxhist13_outname 
  auxhist13_interval_y       = model_config_rec% auxhist13_interval_y 
  auxhist13_interval_d       = model_config_rec% auxhist13_interval_d 
  auxhist13_interval_h       = model_config_rec% auxhist13_interval_h 
  auxhist13_interval_m       = model_config_rec% auxhist13_interval_m 
  auxhist13_interval_s       = model_config_rec% auxhist13_interval_s 
  auxhist13_interval         = model_config_rec% auxhist13_interval 
  auxhist13_begin_y          = model_config_rec% auxhist13_begin_y 
  auxhist13_begin_d          = model_config_rec% auxhist13_begin_d 
  auxhist13_begin_h          = model_config_rec% auxhist13_begin_h 
  auxhist13_begin_m          = model_config_rec% auxhist13_begin_m 
  auxhist13_begin_s          = model_config_rec% auxhist13_begin_s 
  auxhist13_begin            = model_config_rec% auxhist13_begin 
  auxhist13_end_y            = model_config_rec% auxhist13_end_y 
  auxhist13_end_d            = model_config_rec% auxhist13_end_d 
  auxhist13_end_h            = model_config_rec% auxhist13_end_h 
  auxhist13_end_m            = model_config_rec% auxhist13_end_m 
  auxhist13_end_s            = model_config_rec% auxhist13_end_s 
  auxhist13_end              = model_config_rec% auxhist13_end 
  io_form_auxhist13          = model_config_rec% io_form_auxhist13 
  frames_per_auxhist13       = model_config_rec% frames_per_auxhist13 
  auxhist14_inname           = model_config_rec% auxhist14_inname 
  auxhist14_outname          = model_config_rec% auxhist14_outname 
  auxhist14_interval_y       = model_config_rec% auxhist14_interval_y 
  auxhist14_interval_d       = model_config_rec% auxhist14_interval_d 
  auxhist14_interval_h       = model_config_rec% auxhist14_interval_h 
  auxhist14_interval_m       = model_config_rec% auxhist14_interval_m 
  auxhist14_interval_s       = model_config_rec% auxhist14_interval_s 
  auxhist14_interval         = model_config_rec% auxhist14_interval 
  auxhist14_begin_y          = model_config_rec% auxhist14_begin_y 
  auxhist14_begin_d          = model_config_rec% auxhist14_begin_d 
  auxhist14_begin_h          = model_config_rec% auxhist14_begin_h 
  auxhist14_begin_m          = model_config_rec% auxhist14_begin_m 
  auxhist14_begin_s          = model_config_rec% auxhist14_begin_s 
  auxhist14_begin            = model_config_rec% auxhist14_begin 
  auxhist14_end_y            = model_config_rec% auxhist14_end_y 
  auxhist14_end_d            = model_config_rec% auxhist14_end_d 
  auxhist14_end_h            = model_config_rec% auxhist14_end_h 
  auxhist14_end_m            = model_config_rec% auxhist14_end_m 
  auxhist14_end_s            = model_config_rec% auxhist14_end_s 
  auxhist14_end              = model_config_rec% auxhist14_end 
  io_form_auxhist14          = model_config_rec% io_form_auxhist14 
  frames_per_auxhist14       = model_config_rec% frames_per_auxhist14 
  auxhist15_inname           = model_config_rec% auxhist15_inname 
  auxhist15_outname          = model_config_rec% auxhist15_outname 
  auxhist15_interval_y       = model_config_rec% auxhist15_interval_y 
  auxhist15_interval_d       = model_config_rec% auxhist15_interval_d 
  auxhist15_interval_h       = model_config_rec% auxhist15_interval_h 
  auxhist15_interval_m       = model_config_rec% auxhist15_interval_m 
  auxhist15_interval_s       = model_config_rec% auxhist15_interval_s 
  auxhist15_interval         = model_config_rec% auxhist15_interval 
  auxhist15_begin_y          = model_config_rec% auxhist15_begin_y 
  auxhist15_begin_d          = model_config_rec% auxhist15_begin_d 
  auxhist15_begin_h          = model_config_rec% auxhist15_begin_h 
  auxhist15_begin_m          = model_config_rec% auxhist15_begin_m 
  auxhist15_begin_s          = model_config_rec% auxhist15_begin_s 
  auxhist15_begin            = model_config_rec% auxhist15_begin 
  auxhist15_end_y            = model_config_rec% auxhist15_end_y 
  auxhist15_end_d            = model_config_rec% auxhist15_end_d 
  auxhist15_end_h            = model_config_rec% auxhist15_end_h 
  auxhist15_end_m            = model_config_rec% auxhist15_end_m 
  auxhist15_end_s            = model_config_rec% auxhist15_end_s 
  auxhist15_end              = model_config_rec% auxhist15_end 
  io_form_auxhist15          = model_config_rec% io_form_auxhist15 
  frames_per_auxhist15       = model_config_rec% frames_per_auxhist15 
  auxhist16_inname           = model_config_rec% auxhist16_inname 
  auxhist16_outname          = model_config_rec% auxhist16_outname 
  auxhist16_interval_y       = model_config_rec% auxhist16_interval_y 
  auxhist16_interval_d       = model_config_rec% auxhist16_interval_d 
  auxhist16_interval_h       = model_config_rec% auxhist16_interval_h 
  auxhist16_interval_m       = model_config_rec% auxhist16_interval_m 
  auxhist16_interval_s       = model_config_rec% auxhist16_interval_s 
  auxhist16_interval         = model_config_rec% auxhist16_interval 
  auxhist16_begin_y          = model_config_rec% auxhist16_begin_y 
  auxhist16_begin_d          = model_config_rec% auxhist16_begin_d 
  auxhist16_begin_h          = model_config_rec% auxhist16_begin_h 
  auxhist16_begin_m          = model_config_rec% auxhist16_begin_m 
  auxhist16_begin_s          = model_config_rec% auxhist16_begin_s 
  auxhist16_begin            = model_config_rec% auxhist16_begin 
  auxhist16_end_y            = model_config_rec% auxhist16_end_y 
  auxhist16_end_d            = model_config_rec% auxhist16_end_d 
  auxhist16_end_h            = model_config_rec% auxhist16_end_h 
  auxhist16_end_m            = model_config_rec% auxhist16_end_m 
  auxhist16_end_s            = model_config_rec% auxhist16_end_s 
  auxhist16_end              = model_config_rec% auxhist16_end 
  io_form_auxhist16          = model_config_rec% io_form_auxhist16 
  frames_per_auxhist16       = model_config_rec% frames_per_auxhist16 
  auxhist17_inname           = model_config_rec% auxhist17_inname 
  auxhist17_outname          = model_config_rec% auxhist17_outname 
  auxhist17_interval_y       = model_config_rec% auxhist17_interval_y 
  auxhist17_interval_d       = model_config_rec% auxhist17_interval_d 
  auxhist17_interval_h       = model_config_rec% auxhist17_interval_h 
  auxhist17_interval_m       = model_config_rec% auxhist17_interval_m 
  auxhist17_interval_s       = model_config_rec% auxhist17_interval_s 
  auxhist17_interval         = model_config_rec% auxhist17_interval 
  auxhist17_begin_y          = model_config_rec% auxhist17_begin_y 
  auxhist17_begin_d          = model_config_rec% auxhist17_begin_d 
  auxhist17_begin_h          = model_config_rec% auxhist17_begin_h 
  auxhist17_begin_m          = model_config_rec% auxhist17_begin_m 
  auxhist17_begin_s          = model_config_rec% auxhist17_begin_s 
  auxhist17_begin            = model_config_rec% auxhist17_begin 
  auxhist17_end_y            = model_config_rec% auxhist17_end_y 
  auxhist17_end_d            = model_config_rec% auxhist17_end_d 
  auxhist17_end_h            = model_config_rec% auxhist17_end_h 
  auxhist17_end_m            = model_config_rec% auxhist17_end_m 
  auxhist17_end_s            = model_config_rec% auxhist17_end_s 
  auxhist17_end              = model_config_rec% auxhist17_end 
  io_form_auxhist17          = model_config_rec% io_form_auxhist17 
  frames_per_auxhist17       = model_config_rec% frames_per_auxhist17 
  auxhist18_inname           = model_config_rec% auxhist18_inname 
  auxhist18_outname          = model_config_rec% auxhist18_outname 
  auxhist18_interval_y       = model_config_rec% auxhist18_interval_y 
  auxhist18_interval_d       = model_config_rec% auxhist18_interval_d 
  auxhist18_interval_h       = model_config_rec% auxhist18_interval_h 
  auxhist18_interval_m       = model_config_rec% auxhist18_interval_m 
  auxhist18_interval_s       = model_config_rec% auxhist18_interval_s 
  auxhist18_interval         = model_config_rec% auxhist18_interval 
  auxhist18_begin_y          = model_config_rec% auxhist18_begin_y 
  auxhist18_begin_d          = model_config_rec% auxhist18_begin_d 
  auxhist18_begin_h          = model_config_rec% auxhist18_begin_h 
  auxhist18_begin_m          = model_config_rec% auxhist18_begin_m 
  auxhist18_begin_s          = model_config_rec% auxhist18_begin_s 
  auxhist18_begin            = model_config_rec% auxhist18_begin 
  auxhist18_end_y            = model_config_rec% auxhist18_end_y 
  auxhist18_end_d            = model_config_rec% auxhist18_end_d 
  auxhist18_end_h            = model_config_rec% auxhist18_end_h 
  auxhist18_end_m            = model_config_rec% auxhist18_end_m 
  auxhist18_end_s            = model_config_rec% auxhist18_end_s 
  auxhist18_end              = model_config_rec% auxhist18_end 
  io_form_auxhist18          = model_config_rec% io_form_auxhist18 
  frames_per_auxhist18       = model_config_rec% frames_per_auxhist18 
  auxhist19_inname           = model_config_rec% auxhist19_inname 
  auxhist19_outname          = model_config_rec% auxhist19_outname 
  auxhist19_interval_y       = model_config_rec% auxhist19_interval_y 
  auxhist19_interval_d       = model_config_rec% auxhist19_interval_d 
  auxhist19_interval_h       = model_config_rec% auxhist19_interval_h 
  auxhist19_interval_m       = model_config_rec% auxhist19_interval_m 
  auxhist19_interval_s       = model_config_rec% auxhist19_interval_s 
  auxhist19_interval         = model_config_rec% auxhist19_interval 
  auxhist19_begin_y          = model_config_rec% auxhist19_begin_y 
  auxhist19_begin_d          = model_config_rec% auxhist19_begin_d 
  auxhist19_begin_h          = model_config_rec% auxhist19_begin_h 
  auxhist19_begin_m          = model_config_rec% auxhist19_begin_m 
  auxhist19_begin_s          = model_config_rec% auxhist19_begin_s 
  auxhist19_begin            = model_config_rec% auxhist19_begin 
  auxhist19_end_y            = model_config_rec% auxhist19_end_y 
  auxhist19_end_d            = model_config_rec% auxhist19_end_d 
  auxhist19_end_h            = model_config_rec% auxhist19_end_h 
  auxhist19_end_m            = model_config_rec% auxhist19_end_m 
  auxhist19_end_s            = model_config_rec% auxhist19_end_s 
  auxhist19_end              = model_config_rec% auxhist19_end 
  io_form_auxhist19          = model_config_rec% io_form_auxhist19 
  frames_per_auxhist19       = model_config_rec% frames_per_auxhist19 
  auxhist20_inname           = model_config_rec% auxhist20_inname 
  auxhist20_outname          = model_config_rec% auxhist20_outname 
  auxhist20_interval_y       = model_config_rec% auxhist20_interval_y 
  auxhist20_interval_d       = model_config_rec% auxhist20_interval_d 
  auxhist20_interval_h       = model_config_rec% auxhist20_interval_h 
  auxhist20_interval_m       = model_config_rec% auxhist20_interval_m 
  auxhist20_interval_s       = model_config_rec% auxhist20_interval_s 
  auxhist20_interval         = model_config_rec% auxhist20_interval 
  auxhist20_begin_y          = model_config_rec% auxhist20_begin_y 
  auxhist20_begin_d          = model_config_rec% auxhist20_begin_d 
  auxhist20_begin_h          = model_config_rec% auxhist20_begin_h 
  auxhist20_begin_m          = model_config_rec% auxhist20_begin_m 
  auxhist20_begin_s          = model_config_rec% auxhist20_begin_s 
  auxhist20_begin            = model_config_rec% auxhist20_begin 
  auxhist20_end_y            = model_config_rec% auxhist20_end_y 
  auxhist20_end_d            = model_config_rec% auxhist20_end_d 
  auxhist20_end_h            = model_config_rec% auxhist20_end_h 
  auxhist20_end_m            = model_config_rec% auxhist20_end_m 
  auxhist20_end_s            = model_config_rec% auxhist20_end_s 
  auxhist20_end              = model_config_rec% auxhist20_end 
  io_form_auxhist20          = model_config_rec% io_form_auxhist20 
  frames_per_auxhist20       = model_config_rec% frames_per_auxhist20 
  auxhist21_inname           = model_config_rec% auxhist21_inname 
  auxhist21_outname          = model_config_rec% auxhist21_outname 
  auxhist21_interval_y       = model_config_rec% auxhist21_interval_y 
  auxhist21_interval_d       = model_config_rec% auxhist21_interval_d 
  auxhist21_interval_h       = model_config_rec% auxhist21_interval_h 
  auxhist21_interval_m       = model_config_rec% auxhist21_interval_m 
  auxhist21_interval_s       = model_config_rec% auxhist21_interval_s 
  auxhist21_interval         = model_config_rec% auxhist21_interval 
  auxhist21_begin_y          = model_config_rec% auxhist21_begin_y 
  auxhist21_begin_d          = model_config_rec% auxhist21_begin_d 
  auxhist21_begin_h          = model_config_rec% auxhist21_begin_h 
  auxhist21_begin_m          = model_config_rec% auxhist21_begin_m 
  auxhist21_begin_s          = model_config_rec% auxhist21_begin_s 
  auxhist21_begin            = model_config_rec% auxhist21_begin 
  auxhist21_end_y            = model_config_rec% auxhist21_end_y 
  auxhist21_end_d            = model_config_rec% auxhist21_end_d 
  auxhist21_end_h            = model_config_rec% auxhist21_end_h 
  auxhist21_end_m            = model_config_rec% auxhist21_end_m 
  auxhist21_end_s            = model_config_rec% auxhist21_end_s 
  auxhist21_end              = model_config_rec% auxhist21_end 
  io_form_auxhist21          = model_config_rec% io_form_auxhist21 
  frames_per_auxhist21       = model_config_rec% frames_per_auxhist21 
  auxhist22_inname           = model_config_rec% auxhist22_inname 
  auxhist22_outname          = model_config_rec% auxhist22_outname 
  auxhist22_interval_y       = model_config_rec% auxhist22_interval_y 
  auxhist22_interval_d       = model_config_rec% auxhist22_interval_d 
  auxhist22_interval_h       = model_config_rec% auxhist22_interval_h 
  auxhist22_interval_m       = model_config_rec% auxhist22_interval_m 
  auxhist22_interval_s       = model_config_rec% auxhist22_interval_s 
  auxhist22_interval         = model_config_rec% auxhist22_interval 
  auxhist22_begin_y          = model_config_rec% auxhist22_begin_y 
  auxhist22_begin_d          = model_config_rec% auxhist22_begin_d 
  auxhist22_begin_h          = model_config_rec% auxhist22_begin_h 
  auxhist22_begin_m          = model_config_rec% auxhist22_begin_m 
  auxhist22_begin_s          = model_config_rec% auxhist22_begin_s 
  auxhist22_begin            = model_config_rec% auxhist22_begin 
  auxhist22_end_y            = model_config_rec% auxhist22_end_y 
  auxhist22_end_d            = model_config_rec% auxhist22_end_d 
  auxhist22_end_h            = model_config_rec% auxhist22_end_h 
  auxhist22_end_m            = model_config_rec% auxhist22_end_m 
  auxhist22_end_s            = model_config_rec% auxhist22_end_s 
  auxhist22_end              = model_config_rec% auxhist22_end 
  io_form_auxhist22          = model_config_rec% io_form_auxhist22 
  frames_per_auxhist22       = model_config_rec% frames_per_auxhist22 
  auxhist23_inname           = model_config_rec% auxhist23_inname 
  auxhist23_outname          = model_config_rec% auxhist23_outname 
  auxhist23_interval_y       = model_config_rec% auxhist23_interval_y 
  auxhist23_interval_d       = model_config_rec% auxhist23_interval_d 
  auxhist23_interval_h       = model_config_rec% auxhist23_interval_h 
  auxhist23_interval_m       = model_config_rec% auxhist23_interval_m 
  auxhist23_interval_s       = model_config_rec% auxhist23_interval_s 
  auxhist23_interval         = model_config_rec% auxhist23_interval 
  auxhist23_begin_y          = model_config_rec% auxhist23_begin_y 
  auxhist23_begin_d          = model_config_rec% auxhist23_begin_d 
  auxhist23_begin_h          = model_config_rec% auxhist23_begin_h 
  auxhist23_begin_m          = model_config_rec% auxhist23_begin_m 
  auxhist23_begin_s          = model_config_rec% auxhist23_begin_s 
  auxhist23_begin            = model_config_rec% auxhist23_begin 
  auxhist23_end_y            = model_config_rec% auxhist23_end_y 
  auxhist23_end_d            = model_config_rec% auxhist23_end_d 
  auxhist23_end_h            = model_config_rec% auxhist23_end_h 
  auxhist23_end_m            = model_config_rec% auxhist23_end_m 
  auxhist23_end_s            = model_config_rec% auxhist23_end_s 
  auxhist23_end              = model_config_rec% auxhist23_end 
  io_form_auxhist23          = model_config_rec% io_form_auxhist23 
  frames_per_auxhist23       = model_config_rec% frames_per_auxhist23 
  auxhist24_inname           = model_config_rec% auxhist24_inname 
  auxhist24_outname          = model_config_rec% auxhist24_outname 
  auxhist24_interval_y       = model_config_rec% auxhist24_interval_y 
  auxhist24_interval_d       = model_config_rec% auxhist24_interval_d 
  auxhist24_interval_h       = model_config_rec% auxhist24_interval_h 
  auxhist24_interval_m       = model_config_rec% auxhist24_interval_m 
  auxhist24_interval_s       = model_config_rec% auxhist24_interval_s 
  auxhist24_interval         = model_config_rec% auxhist24_interval 
  auxhist24_begin_y          = model_config_rec% auxhist24_begin_y 
  auxhist24_begin_d          = model_config_rec% auxhist24_begin_d 
  auxhist24_begin_h          = model_config_rec% auxhist24_begin_h 
  auxhist24_begin_m          = model_config_rec% auxhist24_begin_m 
  auxhist24_begin_s          = model_config_rec% auxhist24_begin_s 
  auxhist24_begin            = model_config_rec% auxhist24_begin 
  auxhist24_end_y            = model_config_rec% auxhist24_end_y 
  auxhist24_end_d            = model_config_rec% auxhist24_end_d 
  auxhist24_end_h            = model_config_rec% auxhist24_end_h 
  auxhist24_end_m            = model_config_rec% auxhist24_end_m 
  auxhist24_end_s            = model_config_rec% auxhist24_end_s 
  auxhist24_end              = model_config_rec% auxhist24_end 
  io_form_auxhist24          = model_config_rec% io_form_auxhist24 
  frames_per_auxhist24       = model_config_rec% frames_per_auxhist24 
  auxinput1_outname          = model_config_rec% auxinput1_outname 
  auxinput1_interval_y       = model_config_rec% auxinput1_interval_y 
  auxinput1_interval_d       = model_config_rec% auxinput1_interval_d 
  auxinput1_interval_h       = model_config_rec% auxinput1_interval_h 
  auxinput1_interval_m       = model_config_rec% auxinput1_interval_m 
  auxinput1_interval_s       = model_config_rec% auxinput1_interval_s 
  auxinput1_interval         = model_config_rec% auxinput1_interval 
  auxinput1_begin_y          = model_config_rec% auxinput1_begin_y 
  auxinput1_begin_d          = model_config_rec% auxinput1_begin_d 
  auxinput1_begin_h          = model_config_rec% auxinput1_begin_h 
  auxinput1_begin_m          = model_config_rec% auxinput1_begin_m 
  auxinput1_begin_s          = model_config_rec% auxinput1_begin_s 
  auxinput1_begin            = model_config_rec% auxinput1_begin 
  auxinput1_end_y            = model_config_rec% auxinput1_end_y 
  auxinput1_end_d            = model_config_rec% auxinput1_end_d 
  auxinput1_end_h            = model_config_rec% auxinput1_end_h 
  auxinput1_end_m            = model_config_rec% auxinput1_end_m 
  auxinput1_end_s            = model_config_rec% auxinput1_end_s 
  auxinput1_end              = model_config_rec% auxinput1_end 
  frames_per_auxinput1       = model_config_rec% frames_per_auxinput1 
  auxinput2_inname           = model_config_rec% auxinput2_inname 
  auxinput2_outname          = model_config_rec% auxinput2_outname 
  auxinput2_interval_y       = model_config_rec% auxinput2_interval_y 
  auxinput2_interval_d       = model_config_rec% auxinput2_interval_d 
  auxinput2_interval_h       = model_config_rec% auxinput2_interval_h 
  auxinput2_interval_m       = model_config_rec% auxinput2_interval_m 
  auxinput2_interval_s       = model_config_rec% auxinput2_interval_s 
  auxinput2_interval         = model_config_rec% auxinput2_interval 
  auxinput2_begin_y          = model_config_rec% auxinput2_begin_y 
  auxinput2_begin_d          = model_config_rec% auxinput2_begin_d 
  auxinput2_begin_h          = model_config_rec% auxinput2_begin_h 
  auxinput2_begin_m          = model_config_rec% auxinput2_begin_m 
  auxinput2_begin_s          = model_config_rec% auxinput2_begin_s 
  auxinput2_begin            = model_config_rec% auxinput2_begin 
  auxinput2_end_y            = model_config_rec% auxinput2_end_y 
  auxinput2_end_d            = model_config_rec% auxinput2_end_d 
  auxinput2_end_h            = model_config_rec% auxinput2_end_h 
  auxinput2_end_m            = model_config_rec% auxinput2_end_m 
  auxinput2_end_s            = model_config_rec% auxinput2_end_s 
  auxinput2_end              = model_config_rec% auxinput2_end 
  io_form_auxinput2          = model_config_rec% io_form_auxinput2 
  frames_per_auxinput2       = model_config_rec% frames_per_auxinput2 
  auxinput3_inname           = model_config_rec% auxinput3_inname 
  auxinput3_outname          = model_config_rec% auxinput3_outname 
  auxinput3_interval_y       = model_config_rec% auxinput3_interval_y 
  auxinput3_interval_d       = model_config_rec% auxinput3_interval_d 
  auxinput3_interval_h       = model_config_rec% auxinput3_interval_h 
  auxinput3_interval_m       = model_config_rec% auxinput3_interval_m 
  auxinput3_interval_s       = model_config_rec% auxinput3_interval_s 
  auxinput3_interval         = model_config_rec% auxinput3_interval 
  auxinput3_begin_y          = model_config_rec% auxinput3_begin_y 
  auxinput3_begin_d          = model_config_rec% auxinput3_begin_d 
  auxinput3_begin_h          = model_config_rec% auxinput3_begin_h 
  auxinput3_begin_m          = model_config_rec% auxinput3_begin_m 
  auxinput3_begin_s          = model_config_rec% auxinput3_begin_s 
  auxinput3_begin            = model_config_rec% auxinput3_begin 
  auxinput3_end_y            = model_config_rec% auxinput3_end_y 
  auxinput3_end_d            = model_config_rec% auxinput3_end_d 
  auxinput3_end_h            = model_config_rec% auxinput3_end_h 
  auxinput3_end_m            = model_config_rec% auxinput3_end_m 
  auxinput3_end_s            = model_config_rec% auxinput3_end_s 
  auxinput3_end              = model_config_rec% auxinput3_end 
  io_form_auxinput3          = model_config_rec% io_form_auxinput3 
  frames_per_auxinput3       = model_config_rec% frames_per_auxinput3 
  auxinput4_inname           = model_config_rec% auxinput4_inname 
  auxinput4_outname          = model_config_rec% auxinput4_outname 
  auxinput4_interval_y       = model_config_rec% auxinput4_interval_y 
  auxinput4_interval_d       = model_config_rec% auxinput4_interval_d 
  auxinput4_interval_h       = model_config_rec% auxinput4_interval_h 
  auxinput4_interval_m       = model_config_rec% auxinput4_interval_m 
  auxinput4_interval_s       = model_config_rec% auxinput4_interval_s 
  auxinput4_interval         = model_config_rec% auxinput4_interval 
  auxinput4_begin_y          = model_config_rec% auxinput4_begin_y 
  auxinput4_begin_d          = model_config_rec% auxinput4_begin_d 
  auxinput4_begin_h          = model_config_rec% auxinput4_begin_h 
  auxinput4_begin_m          = model_config_rec% auxinput4_begin_m 
  auxinput4_begin_s          = model_config_rec% auxinput4_begin_s 
  auxinput4_begin            = model_config_rec% auxinput4_begin 
  auxinput4_end_y            = model_config_rec% auxinput4_end_y 
  auxinput4_end_d            = model_config_rec% auxinput4_end_d 
  auxinput4_end_h            = model_config_rec% auxinput4_end_h 
  auxinput4_end_m            = model_config_rec% auxinput4_end_m 
  auxinput4_end_s            = model_config_rec% auxinput4_end_s 
  auxinput4_end              = model_config_rec% auxinput4_end 
  io_form_auxinput4          = model_config_rec% io_form_auxinput4 
  frames_per_auxinput4       = model_config_rec% frames_per_auxinput4 
  auxinput5_inname           = model_config_rec% auxinput5_inname 
  auxinput5_outname          = model_config_rec% auxinput5_outname 
  auxinput5_interval_y       = model_config_rec% auxinput5_interval_y 
  auxinput5_interval_d       = model_config_rec% auxinput5_interval_d 
  auxinput5_interval_h       = model_config_rec% auxinput5_interval_h 
  auxinput5_interval_m       = model_config_rec% auxinput5_interval_m 
  auxinput5_interval_s       = model_config_rec% auxinput5_interval_s 
  auxinput5_interval         = model_config_rec% auxinput5_interval 
  auxinput5_begin_y          = model_config_rec% auxinput5_begin_y 
  auxinput5_begin_d          = model_config_rec% auxinput5_begin_d 
  auxinput5_begin_h          = model_config_rec% auxinput5_begin_h 
  auxinput5_begin_m          = model_config_rec% auxinput5_begin_m 
  auxinput5_begin_s          = model_config_rec% auxinput5_begin_s 
  auxinput5_begin            = model_config_rec% auxinput5_begin 
  auxinput5_end_y            = model_config_rec% auxinput5_end_y 
  auxinput5_end_d            = model_config_rec% auxinput5_end_d 
  auxinput5_end_h            = model_config_rec% auxinput5_end_h 
  auxinput5_end_m            = model_config_rec% auxinput5_end_m 
  auxinput5_end_s            = model_config_rec% auxinput5_end_s 
  auxinput5_end              = model_config_rec% auxinput5_end 
  io_form_auxinput5          = model_config_rec% io_form_auxinput5 
  frames_per_auxinput5       = model_config_rec% frames_per_auxinput5 
  auxinput6_inname           = model_config_rec% auxinput6_inname 
  auxinput6_outname          = model_config_rec% auxinput6_outname 
  auxinput6_interval_y       = model_config_rec% auxinput6_interval_y 
  auxinput6_interval_d       = model_config_rec% auxinput6_interval_d 
  auxinput6_interval_h       = model_config_rec% auxinput6_interval_h 
  auxinput6_interval_m       = model_config_rec% auxinput6_interval_m 
  auxinput6_interval_s       = model_config_rec% auxinput6_interval_s 
  auxinput6_interval         = model_config_rec% auxinput6_interval 
  auxinput6_begin_y          = model_config_rec% auxinput6_begin_y 
  auxinput6_begin_d          = model_config_rec% auxinput6_begin_d 
  auxinput6_begin_h          = model_config_rec% auxinput6_begin_h 
  auxinput6_begin_m          = model_config_rec% auxinput6_begin_m 
  auxinput6_begin_s          = model_config_rec% auxinput6_begin_s 
  auxinput6_begin            = model_config_rec% auxinput6_begin 
  auxinput6_end_y            = model_config_rec% auxinput6_end_y 
  auxinput6_end_d            = model_config_rec% auxinput6_end_d 
  auxinput6_end_h            = model_config_rec% auxinput6_end_h 
  auxinput6_end_m            = model_config_rec% auxinput6_end_m 
  auxinput6_end_s            = model_config_rec% auxinput6_end_s 
  auxinput6_end              = model_config_rec% auxinput6_end 
  io_form_auxinput6          = model_config_rec% io_form_auxinput6 
  frames_per_auxinput6       = model_config_rec% frames_per_auxinput6 
  auxinput7_inname           = model_config_rec% auxinput7_inname 
  auxinput7_outname          = model_config_rec% auxinput7_outname 
  auxinput7_interval_y       = model_config_rec% auxinput7_interval_y 
  auxinput7_interval_d       = model_config_rec% auxinput7_interval_d 
  auxinput7_interval_h       = model_config_rec% auxinput7_interval_h 
  auxinput7_interval_m       = model_config_rec% auxinput7_interval_m 
  auxinput7_interval_s       = model_config_rec% auxinput7_interval_s 
  auxinput7_interval         = model_config_rec% auxinput7_interval 
  auxinput7_begin_y          = model_config_rec% auxinput7_begin_y 
  auxinput7_begin_d          = model_config_rec% auxinput7_begin_d 
  auxinput7_begin_h          = model_config_rec% auxinput7_begin_h 
  auxinput7_begin_m          = model_config_rec% auxinput7_begin_m 
  auxinput7_begin_s          = model_config_rec% auxinput7_begin_s 
  auxinput7_begin            = model_config_rec% auxinput7_begin 
  auxinput7_end_y            = model_config_rec% auxinput7_end_y 
  auxinput7_end_d            = model_config_rec% auxinput7_end_d 
  auxinput7_end_h            = model_config_rec% auxinput7_end_h 
  auxinput7_end_m            = model_config_rec% auxinput7_end_m 
  auxinput7_end_s            = model_config_rec% auxinput7_end_s 
  auxinput7_end              = model_config_rec% auxinput7_end 
  io_form_auxinput7          = model_config_rec% io_form_auxinput7 
  frames_per_auxinput7       = model_config_rec% frames_per_auxinput7 
  auxinput8_inname           = model_config_rec% auxinput8_inname 
  auxinput8_outname          = model_config_rec% auxinput8_outname 
  auxinput8_interval_y       = model_config_rec% auxinput8_interval_y 
  auxinput8_interval_d       = model_config_rec% auxinput8_interval_d 
  auxinput8_interval_h       = model_config_rec% auxinput8_interval_h 
  auxinput8_interval_m       = model_config_rec% auxinput8_interval_m 
  auxinput8_interval_s       = model_config_rec% auxinput8_interval_s 
  auxinput8_interval         = model_config_rec% auxinput8_interval 
  auxinput8_begin_y          = model_config_rec% auxinput8_begin_y 
  auxinput8_begin_d          = model_config_rec% auxinput8_begin_d 
  auxinput8_begin_h          = model_config_rec% auxinput8_begin_h 
  auxinput8_begin_m          = model_config_rec% auxinput8_begin_m 
  auxinput8_begin_s          = model_config_rec% auxinput8_begin_s 
  auxinput8_begin            = model_config_rec% auxinput8_begin 
  auxinput8_end_y            = model_config_rec% auxinput8_end_y 
  auxinput8_end_d            = model_config_rec% auxinput8_end_d 
  auxinput8_end_h            = model_config_rec% auxinput8_end_h 
  auxinput8_end_m            = model_config_rec% auxinput8_end_m 
  auxinput8_end_s            = model_config_rec% auxinput8_end_s 
  auxinput8_end              = model_config_rec% auxinput8_end 
  io_form_auxinput8          = model_config_rec% io_form_auxinput8 
  frames_per_auxinput8       = model_config_rec% frames_per_auxinput8 
  auxinput9_inname           = model_config_rec% auxinput9_inname 
  auxinput9_outname          = model_config_rec% auxinput9_outname 
  auxinput9_interval_y       = model_config_rec% auxinput9_interval_y 
  auxinput9_interval_d       = model_config_rec% auxinput9_interval_d 
  auxinput9_interval_h       = model_config_rec% auxinput9_interval_h 
  auxinput9_interval_m       = model_config_rec% auxinput9_interval_m 
  auxinput9_interval_s       = model_config_rec% auxinput9_interval_s 
  auxinput9_interval         = model_config_rec% auxinput9_interval 
  auxinput9_begin_y          = model_config_rec% auxinput9_begin_y 
  auxinput9_begin_d          = model_config_rec% auxinput9_begin_d 
  auxinput9_begin_h          = model_config_rec% auxinput9_begin_h 
  auxinput9_begin_m          = model_config_rec% auxinput9_begin_m 
  auxinput9_begin_s          = model_config_rec% auxinput9_begin_s 
  auxinput9_begin            = model_config_rec% auxinput9_begin 
  auxinput9_end_y            = model_config_rec% auxinput9_end_y 
  auxinput9_end_d            = model_config_rec% auxinput9_end_d 
  auxinput9_end_h            = model_config_rec% auxinput9_end_h 
  auxinput9_end_m            = model_config_rec% auxinput9_end_m 
  auxinput9_end_s            = model_config_rec% auxinput9_end_s 
  auxinput9_end              = model_config_rec% auxinput9_end 
  io_form_auxinput9          = model_config_rec% io_form_auxinput9 
  frames_per_auxinput9       = model_config_rec% frames_per_auxinput9 
  auxinput10_inname          = model_config_rec% auxinput10_inname 
  auxinput10_outname         = model_config_rec% auxinput10_outname 
  auxinput10_interval_y      = model_config_rec% auxinput10_interval_y 
  auxinput10_interval_d      = model_config_rec% auxinput10_interval_d 
  auxinput10_interval_h      = model_config_rec% auxinput10_interval_h 
  auxinput10_interval_m      = model_config_rec% auxinput10_interval_m 
  auxinput10_interval_s      = model_config_rec% auxinput10_interval_s 
  auxinput10_interval        = model_config_rec% auxinput10_interval 
  auxinput10_begin_y         = model_config_rec% auxinput10_begin_y 
  auxinput10_begin_d         = model_config_rec% auxinput10_begin_d 
  auxinput10_begin_h         = model_config_rec% auxinput10_begin_h 
  auxinput10_begin_m         = model_config_rec% auxinput10_begin_m 
  auxinput10_begin_s         = model_config_rec% auxinput10_begin_s 
  auxinput10_begin           = model_config_rec% auxinput10_begin 
  auxinput10_end_y           = model_config_rec% auxinput10_end_y 
  auxinput10_end_d           = model_config_rec% auxinput10_end_d 
  auxinput10_end_h           = model_config_rec% auxinput10_end_h 
  auxinput10_end_m           = model_config_rec% auxinput10_end_m 
  auxinput10_end_s           = model_config_rec% auxinput10_end_s 
  auxinput10_end             = model_config_rec% auxinput10_end 
  io_form_auxinput10         = model_config_rec% io_form_auxinput10 
  frames_per_auxinput10      = model_config_rec% frames_per_auxinput10 
  auxinput11_inname          = model_config_rec% auxinput11_inname 
  auxinput11_outname         = model_config_rec% auxinput11_outname 
  auxinput11_interval_y      = model_config_rec% auxinput11_interval_y 
  auxinput11_interval_d      = model_config_rec% auxinput11_interval_d 
  auxinput11_interval_h      = model_config_rec% auxinput11_interval_h 
  auxinput11_interval_m      = model_config_rec% auxinput11_interval_m 
  auxinput11_interval_s      = model_config_rec% auxinput11_interval_s 
  auxinput11_interval        = model_config_rec% auxinput11_interval 
  auxinput11_begin_y         = model_config_rec% auxinput11_begin_y 
  auxinput11_begin_d         = model_config_rec% auxinput11_begin_d 
  auxinput11_begin_h         = model_config_rec% auxinput11_begin_h 
  auxinput11_begin_m         = model_config_rec% auxinput11_begin_m 
  auxinput11_begin_s         = model_config_rec% auxinput11_begin_s 
  auxinput11_begin           = model_config_rec% auxinput11_begin 
  auxinput11_end_y           = model_config_rec% auxinput11_end_y 
  auxinput11_end_d           = model_config_rec% auxinput11_end_d 
  auxinput11_end_h           = model_config_rec% auxinput11_end_h 
  auxinput11_end_m           = model_config_rec% auxinput11_end_m 
  auxinput11_end_s           = model_config_rec% auxinput11_end_s 
  auxinput11_end             = model_config_rec% auxinput11_end 
  io_form_auxinput11         = model_config_rec% io_form_auxinput11 
  frames_per_auxinput11      = model_config_rec% frames_per_auxinput11 
  auxinput12_inname          = model_config_rec% auxinput12_inname 
  auxinput12_outname         = model_config_rec% auxinput12_outname 
  auxinput12_interval_y      = model_config_rec% auxinput12_interval_y 
  auxinput12_interval_d      = model_config_rec% auxinput12_interval_d 
  auxinput12_interval_h      = model_config_rec% auxinput12_interval_h 
  auxinput12_interval_m      = model_config_rec% auxinput12_interval_m 
  auxinput12_interval_s      = model_config_rec% auxinput12_interval_s 
  auxinput12_interval        = model_config_rec% auxinput12_interval 
  auxinput12_begin_y         = model_config_rec% auxinput12_begin_y 
  auxinput12_begin_d         = model_config_rec% auxinput12_begin_d 
  auxinput12_begin_h         = model_config_rec% auxinput12_begin_h 
  auxinput12_begin_m         = model_config_rec% auxinput12_begin_m 
  auxinput12_begin_s         = model_config_rec% auxinput12_begin_s 
  auxinput12_begin           = model_config_rec% auxinput12_begin 
  auxinput12_end_y           = model_config_rec% auxinput12_end_y 
  auxinput12_end_d           = model_config_rec% auxinput12_end_d 
  auxinput12_end_h           = model_config_rec% auxinput12_end_h 
  auxinput12_end_m           = model_config_rec% auxinput12_end_m 
  auxinput12_end_s           = model_config_rec% auxinput12_end_s 
  auxinput12_end             = model_config_rec% auxinput12_end 
  io_form_auxinput12         = model_config_rec% io_form_auxinput12 
  frames_per_auxinput12      = model_config_rec% frames_per_auxinput12 
  auxinput13_inname          = model_config_rec% auxinput13_inname 
  auxinput13_outname         = model_config_rec% auxinput13_outname 
  auxinput13_interval_y      = model_config_rec% auxinput13_interval_y 
  auxinput13_interval_d      = model_config_rec% auxinput13_interval_d 
  auxinput13_interval_h      = model_config_rec% auxinput13_interval_h 
  auxinput13_interval_m      = model_config_rec% auxinput13_interval_m 
  auxinput13_interval_s      = model_config_rec% auxinput13_interval_s 
  auxinput13_interval        = model_config_rec% auxinput13_interval 
  auxinput13_begin_y         = model_config_rec% auxinput13_begin_y 
  auxinput13_begin_d         = model_config_rec% auxinput13_begin_d 
  auxinput13_begin_h         = model_config_rec% auxinput13_begin_h 
  auxinput13_begin_m         = model_config_rec% auxinput13_begin_m 
  auxinput13_begin_s         = model_config_rec% auxinput13_begin_s 
  auxinput13_begin           = model_config_rec% auxinput13_begin 
  auxinput13_end_y           = model_config_rec% auxinput13_end_y 
  auxinput13_end_d           = model_config_rec% auxinput13_end_d 
  auxinput13_end_h           = model_config_rec% auxinput13_end_h 
  auxinput13_end_m           = model_config_rec% auxinput13_end_m 
  auxinput13_end_s           = model_config_rec% auxinput13_end_s 
  auxinput13_end             = model_config_rec% auxinput13_end 
  io_form_auxinput13         = model_config_rec% io_form_auxinput13 
  frames_per_auxinput13      = model_config_rec% frames_per_auxinput13 
  auxinput14_inname          = model_config_rec% auxinput14_inname 
  auxinput14_outname         = model_config_rec% auxinput14_outname 
  auxinput14_interval_y      = model_config_rec% auxinput14_interval_y 
  auxinput14_interval_d      = model_config_rec% auxinput14_interval_d 
  auxinput14_interval_h      = model_config_rec% auxinput14_interval_h 
  auxinput14_interval_m      = model_config_rec% auxinput14_interval_m 
  auxinput14_interval_s      = model_config_rec% auxinput14_interval_s 
  auxinput14_interval        = model_config_rec% auxinput14_interval 
  auxinput14_begin_y         = model_config_rec% auxinput14_begin_y 
  auxinput14_begin_d         = model_config_rec% auxinput14_begin_d 
  auxinput14_begin_h         = model_config_rec% auxinput14_begin_h 
  auxinput14_begin_m         = model_config_rec% auxinput14_begin_m 
  auxinput14_begin_s         = model_config_rec% auxinput14_begin_s 
  auxinput14_begin           = model_config_rec% auxinput14_begin 
  auxinput14_end_y           = model_config_rec% auxinput14_end_y 
  auxinput14_end_d           = model_config_rec% auxinput14_end_d 
  auxinput14_end_h           = model_config_rec% auxinput14_end_h 
  auxinput14_end_m           = model_config_rec% auxinput14_end_m 
  auxinput14_end_s           = model_config_rec% auxinput14_end_s 
  auxinput14_end             = model_config_rec% auxinput14_end 
  io_form_auxinput14         = model_config_rec% io_form_auxinput14 
  frames_per_auxinput14      = model_config_rec% frames_per_auxinput14 
  auxinput15_inname          = model_config_rec% auxinput15_inname 
  auxinput15_outname         = model_config_rec% auxinput15_outname 
  auxinput15_interval_y      = model_config_rec% auxinput15_interval_y 
  auxinput15_interval_d      = model_config_rec% auxinput15_interval_d 
  auxinput15_interval_h      = model_config_rec% auxinput15_interval_h 
  auxinput15_interval_m      = model_config_rec% auxinput15_interval_m 
  auxinput15_interval_s      = model_config_rec% auxinput15_interval_s 
  auxinput15_interval        = model_config_rec% auxinput15_interval 
  auxinput15_begin_y         = model_config_rec% auxinput15_begin_y 
  auxinput15_begin_d         = model_config_rec% auxinput15_begin_d 
  auxinput15_begin_h         = model_config_rec% auxinput15_begin_h 
  auxinput15_begin_m         = model_config_rec% auxinput15_begin_m 
  auxinput15_begin_s         = model_config_rec% auxinput15_begin_s 
  auxinput15_begin           = model_config_rec% auxinput15_begin 
  auxinput15_end_y           = model_config_rec% auxinput15_end_y 
  auxinput15_end_d           = model_config_rec% auxinput15_end_d 
  auxinput15_end_h           = model_config_rec% auxinput15_end_h 
  auxinput15_end_m           = model_config_rec% auxinput15_end_m 
  auxinput15_end_s           = model_config_rec% auxinput15_end_s 
  auxinput15_end             = model_config_rec% auxinput15_end 
  io_form_auxinput15         = model_config_rec% io_form_auxinput15 
  frames_per_auxinput15      = model_config_rec% frames_per_auxinput15 
  auxinput16_inname          = model_config_rec% auxinput16_inname 
  auxinput16_outname         = model_config_rec% auxinput16_outname 
  auxinput16_interval_y      = model_config_rec% auxinput16_interval_y 
  auxinput16_interval_d      = model_config_rec% auxinput16_interval_d 
  auxinput16_interval_h      = model_config_rec% auxinput16_interval_h 
  auxinput16_interval_m      = model_config_rec% auxinput16_interval_m 
  auxinput16_interval_s      = model_config_rec% auxinput16_interval_s 
  auxinput16_interval        = model_config_rec% auxinput16_interval 
  auxinput16_begin_y         = model_config_rec% auxinput16_begin_y 
  auxinput16_begin_d         = model_config_rec% auxinput16_begin_d 
  auxinput16_begin_h         = model_config_rec% auxinput16_begin_h 
  auxinput16_begin_m         = model_config_rec% auxinput16_begin_m 
  auxinput16_begin_s         = model_config_rec% auxinput16_begin_s 
  auxinput16_begin           = model_config_rec% auxinput16_begin 
  auxinput16_end_y           = model_config_rec% auxinput16_end_y 
  auxinput16_end_d           = model_config_rec% auxinput16_end_d 
  auxinput16_end_h           = model_config_rec% auxinput16_end_h 
  auxinput16_end_m           = model_config_rec% auxinput16_end_m 
  auxinput16_end_s           = model_config_rec% auxinput16_end_s 
  auxinput16_end             = model_config_rec% auxinput16_end 
  io_form_auxinput16         = model_config_rec% io_form_auxinput16 
  frames_per_auxinput16      = model_config_rec% frames_per_auxinput16 
  auxinput17_inname          = model_config_rec% auxinput17_inname 
  auxinput17_outname         = model_config_rec% auxinput17_outname 
  auxinput17_interval_y      = model_config_rec% auxinput17_interval_y 
  auxinput17_interval_d      = model_config_rec% auxinput17_interval_d 
  auxinput17_interval_h      = model_config_rec% auxinput17_interval_h 
  auxinput17_interval_m      = model_config_rec% auxinput17_interval_m 
  auxinput17_interval_s      = model_config_rec% auxinput17_interval_s 
  auxinput17_interval        = model_config_rec% auxinput17_interval 
  auxinput17_begin_y         = model_config_rec% auxinput17_begin_y 
  auxinput17_begin_d         = model_config_rec% auxinput17_begin_d 
  auxinput17_begin_h         = model_config_rec% auxinput17_begin_h 
  auxinput17_begin_m         = model_config_rec% auxinput17_begin_m 
  auxinput17_begin_s         = model_config_rec% auxinput17_begin_s 
  auxinput17_begin           = model_config_rec% auxinput17_begin 
  auxinput17_end_y           = model_config_rec% auxinput17_end_y 
  auxinput17_end_d           = model_config_rec% auxinput17_end_d 
  auxinput17_end_h           = model_config_rec% auxinput17_end_h 
  auxinput17_end_m           = model_config_rec% auxinput17_end_m 
  auxinput17_end_s           = model_config_rec% auxinput17_end_s 
  auxinput17_end             = model_config_rec% auxinput17_end 
  io_form_auxinput17         = model_config_rec% io_form_auxinput17 
  frames_per_auxinput17      = model_config_rec% frames_per_auxinput17 
  auxinput18_inname          = model_config_rec% auxinput18_inname 
  auxinput18_outname         = model_config_rec% auxinput18_outname 
  auxinput18_interval_y      = model_config_rec% auxinput18_interval_y 
  auxinput18_interval_d      = model_config_rec% auxinput18_interval_d 
  auxinput18_interval_h      = model_config_rec% auxinput18_interval_h 
  auxinput18_interval_m      = model_config_rec% auxinput18_interval_m 
  auxinput18_interval_s      = model_config_rec% auxinput18_interval_s 
  auxinput18_interval        = model_config_rec% auxinput18_interval 
  auxinput18_begin_y         = model_config_rec% auxinput18_begin_y 
  auxinput18_begin_d         = model_config_rec% auxinput18_begin_d 
  auxinput18_begin_h         = model_config_rec% auxinput18_begin_h 
  auxinput18_begin_m         = model_config_rec% auxinput18_begin_m 
  auxinput18_begin_s         = model_config_rec% auxinput18_begin_s 
  auxinput18_begin           = model_config_rec% auxinput18_begin 
  auxinput18_end_y           = model_config_rec% auxinput18_end_y 
  auxinput18_end_d           = model_config_rec% auxinput18_end_d 
  auxinput18_end_h           = model_config_rec% auxinput18_end_h 
  auxinput18_end_m           = model_config_rec% auxinput18_end_m 
  auxinput18_end_s           = model_config_rec% auxinput18_end_s 
  auxinput18_end             = model_config_rec% auxinput18_end 
  io_form_auxinput18         = model_config_rec% io_form_auxinput18 
  frames_per_auxinput18      = model_config_rec% frames_per_auxinput18 
  auxinput19_inname          = model_config_rec% auxinput19_inname 
  auxinput19_outname         = model_config_rec% auxinput19_outname 
  auxinput19_interval_y      = model_config_rec% auxinput19_interval_y 
  auxinput19_interval_d      = model_config_rec% auxinput19_interval_d 
  auxinput19_interval_h      = model_config_rec% auxinput19_interval_h 
  auxinput19_interval_m      = model_config_rec% auxinput19_interval_m 
  auxinput19_interval_s      = model_config_rec% auxinput19_interval_s 
  auxinput19_interval        = model_config_rec% auxinput19_interval 
  auxinput19_begin_y         = model_config_rec% auxinput19_begin_y 
  auxinput19_begin_d         = model_config_rec% auxinput19_begin_d 
  auxinput19_begin_h         = model_config_rec% auxinput19_begin_h 
  auxinput19_begin_m         = model_config_rec% auxinput19_begin_m 
  auxinput19_begin_s         = model_config_rec% auxinput19_begin_s 
  auxinput19_begin           = model_config_rec% auxinput19_begin 
  auxinput19_end_y           = model_config_rec% auxinput19_end_y 
  auxinput19_end_d           = model_config_rec% auxinput19_end_d 
  auxinput19_end_h           = model_config_rec% auxinput19_end_h 
  auxinput19_end_m           = model_config_rec% auxinput19_end_m 
  auxinput19_end_s           = model_config_rec% auxinput19_end_s 
  auxinput19_end             = model_config_rec% auxinput19_end 
  io_form_auxinput19         = model_config_rec% io_form_auxinput19 
  frames_per_auxinput19      = model_config_rec% frames_per_auxinput19 
  auxinput20_inname          = model_config_rec% auxinput20_inname 
  auxinput20_outname         = model_config_rec% auxinput20_outname 
  auxinput20_interval_y      = model_config_rec% auxinput20_interval_y 
  auxinput20_interval_d      = model_config_rec% auxinput20_interval_d 
  auxinput20_interval_h      = model_config_rec% auxinput20_interval_h 
  auxinput20_interval_m      = model_config_rec% auxinput20_interval_m 
  auxinput20_interval_s      = model_config_rec% auxinput20_interval_s 
  auxinput20_interval        = model_config_rec% auxinput20_interval 
  auxinput20_begin_y         = model_config_rec% auxinput20_begin_y 
  auxinput20_begin_d         = model_config_rec% auxinput20_begin_d 
  auxinput20_begin_h         = model_config_rec% auxinput20_begin_h 
  auxinput20_begin_m         = model_config_rec% auxinput20_begin_m 
  auxinput20_begin_s         = model_config_rec% auxinput20_begin_s 
  auxinput20_begin           = model_config_rec% auxinput20_begin 
  auxinput20_end_y           = model_config_rec% auxinput20_end_y 
  auxinput20_end_d           = model_config_rec% auxinput20_end_d 
  auxinput20_end_h           = model_config_rec% auxinput20_end_h 
  auxinput20_end_m           = model_config_rec% auxinput20_end_m 
  auxinput20_end_s           = model_config_rec% auxinput20_end_s 
  auxinput20_end             = model_config_rec% auxinput20_end 
  io_form_auxinput20         = model_config_rec% io_form_auxinput20 
  frames_per_auxinput20      = model_config_rec% frames_per_auxinput20 
  auxinput21_inname          = model_config_rec% auxinput21_inname 
  auxinput21_outname         = model_config_rec% auxinput21_outname 
  auxinput21_interval_y      = model_config_rec% auxinput21_interval_y 
  auxinput21_interval_d      = model_config_rec% auxinput21_interval_d 
  auxinput21_interval_h      = model_config_rec% auxinput21_interval_h 
  auxinput21_interval_m      = model_config_rec% auxinput21_interval_m 
  auxinput21_interval_s      = model_config_rec% auxinput21_interval_s 
  auxinput21_interval        = model_config_rec% auxinput21_interval 
  auxinput21_begin_y         = model_config_rec% auxinput21_begin_y 
  auxinput21_begin_d         = model_config_rec% auxinput21_begin_d 
  auxinput21_begin_h         = model_config_rec% auxinput21_begin_h 
  auxinput21_begin_m         = model_config_rec% auxinput21_begin_m 
  auxinput21_begin_s         = model_config_rec% auxinput21_begin_s 
  auxinput21_begin           = model_config_rec% auxinput21_begin 
  auxinput21_end_y           = model_config_rec% auxinput21_end_y 
  auxinput21_end_d           = model_config_rec% auxinput21_end_d 
  auxinput21_end_h           = model_config_rec% auxinput21_end_h 
  auxinput21_end_m           = model_config_rec% auxinput21_end_m 
  auxinput21_end_s           = model_config_rec% auxinput21_end_s 
  auxinput21_end             = model_config_rec% auxinput21_end 
  io_form_auxinput21         = model_config_rec% io_form_auxinput21 
  frames_per_auxinput21      = model_config_rec% frames_per_auxinput21 
  auxinput22_inname          = model_config_rec% auxinput22_inname 
  auxinput22_outname         = model_config_rec% auxinput22_outname 
  auxinput22_interval_y      = model_config_rec% auxinput22_interval_y 
  auxinput22_interval_d      = model_config_rec% auxinput22_interval_d 
  auxinput22_interval_h      = model_config_rec% auxinput22_interval_h 
  auxinput22_interval_m      = model_config_rec% auxinput22_interval_m 
  auxinput22_interval_s      = model_config_rec% auxinput22_interval_s 
  auxinput22_interval        = model_config_rec% auxinput22_interval 
  auxinput22_begin_y         = model_config_rec% auxinput22_begin_y 
  auxinput22_begin_d         = model_config_rec% auxinput22_begin_d 
  auxinput22_begin_h         = model_config_rec% auxinput22_begin_h 
  auxinput22_begin_m         = model_config_rec% auxinput22_begin_m 
  auxinput22_begin_s         = model_config_rec% auxinput22_begin_s 
  auxinput22_begin           = model_config_rec% auxinput22_begin 
  auxinput22_end_y           = model_config_rec% auxinput22_end_y 
  auxinput22_end_d           = model_config_rec% auxinput22_end_d 
  auxinput22_end_h           = model_config_rec% auxinput22_end_h 
  auxinput22_end_m           = model_config_rec% auxinput22_end_m 
  auxinput22_end_s           = model_config_rec% auxinput22_end_s 
  auxinput22_end             = model_config_rec% auxinput22_end 
  io_form_auxinput22         = model_config_rec% io_form_auxinput22 
  frames_per_auxinput22      = model_config_rec% frames_per_auxinput22 
  auxinput23_inname          = model_config_rec% auxinput23_inname 
  auxinput23_outname         = model_config_rec% auxinput23_outname 
  auxinput23_interval_y      = model_config_rec% auxinput23_interval_y 
  auxinput23_interval_d      = model_config_rec% auxinput23_interval_d 
  auxinput23_interval_h      = model_config_rec% auxinput23_interval_h 
  auxinput23_interval_m      = model_config_rec% auxinput23_interval_m 
  auxinput23_interval_s      = model_config_rec% auxinput23_interval_s 
  auxinput23_interval        = model_config_rec% auxinput23_interval 
  auxinput23_begin_y         = model_config_rec% auxinput23_begin_y 
  auxinput23_begin_d         = model_config_rec% auxinput23_begin_d 
  auxinput23_begin_h         = model_config_rec% auxinput23_begin_h 
  auxinput23_begin_m         = model_config_rec% auxinput23_begin_m 
  auxinput23_begin_s         = model_config_rec% auxinput23_begin_s 
  auxinput23_begin           = model_config_rec% auxinput23_begin 
  auxinput23_end_y           = model_config_rec% auxinput23_end_y 
  auxinput23_end_d           = model_config_rec% auxinput23_end_d 
  auxinput23_end_h           = model_config_rec% auxinput23_end_h 
  auxinput23_end_m           = model_config_rec% auxinput23_end_m 
  auxinput23_end_s           = model_config_rec% auxinput23_end_s 
  auxinput23_end             = model_config_rec% auxinput23_end 
  io_form_auxinput23         = model_config_rec% io_form_auxinput23 
  frames_per_auxinput23      = model_config_rec% frames_per_auxinput23 
  auxinput24_inname          = model_config_rec% auxinput24_inname 
  auxinput24_outname         = model_config_rec% auxinput24_outname 
  auxinput24_interval_y      = model_config_rec% auxinput24_interval_y 
  auxinput24_interval_d      = model_config_rec% auxinput24_interval_d 
  auxinput24_interval_h      = model_config_rec% auxinput24_interval_h 
  auxinput24_interval_m      = model_config_rec% auxinput24_interval_m 
  auxinput24_interval_s      = model_config_rec% auxinput24_interval_s 
  auxinput24_interval        = model_config_rec% auxinput24_interval 
  auxinput24_begin_y         = model_config_rec% auxinput24_begin_y 
  auxinput24_begin_d         = model_config_rec% auxinput24_begin_d 
  auxinput24_begin_h         = model_config_rec% auxinput24_begin_h 
  auxinput24_begin_m         = model_config_rec% auxinput24_begin_m 
  auxinput24_begin_s         = model_config_rec% auxinput24_begin_s 
  auxinput24_begin           = model_config_rec% auxinput24_begin 
  auxinput24_end_y           = model_config_rec% auxinput24_end_y 
  auxinput24_end_d           = model_config_rec% auxinput24_end_d 
  auxinput24_end_h           = model_config_rec% auxinput24_end_h 
  auxinput24_end_m           = model_config_rec% auxinput24_end_m 
  auxinput24_end_s           = model_config_rec% auxinput24_end_s 
  auxinput24_end             = model_config_rec% auxinput24_end 
  io_form_auxinput24         = model_config_rec% io_form_auxinput24 
  frames_per_auxinput24      = model_config_rec% frames_per_auxinput24 
  history_interval           = model_config_rec% history_interval 
  frames_per_outfile         = model_config_rec% frames_per_outfile 
  restart                    = model_config_rec% restart 
  restart_interval           = model_config_rec% restart_interval 
  io_form_input              = model_config_rec% io_form_input 
  io_form_history            = model_config_rec% io_form_history 
  io_form_restart            = model_config_rec% io_form_restart 
  io_form_boundary           = model_config_rec% io_form_boundary 
  debug_level                = model_config_rec% debug_level 
  self_test_domain           = model_config_rec% self_test_domain 
  history_outname            = model_config_rec% history_outname 
  history_inname             = model_config_rec% history_inname 
  use_netcdf_classic         = model_config_rec% use_netcdf_classic 
  history_interval_d         = model_config_rec% history_interval_d 
  history_interval_h         = model_config_rec% history_interval_h 
  history_interval_m         = model_config_rec% history_interval_m 
  history_interval_s         = model_config_rec% history_interval_s 
  inputout_interval_d        = model_config_rec% inputout_interval_d 
  inputout_interval_h        = model_config_rec% inputout_interval_h 
  inputout_interval_m        = model_config_rec% inputout_interval_m 
  inputout_interval_s        = model_config_rec% inputout_interval_s 
  inputout_interval          = model_config_rec% inputout_interval 
  restart_interval_d         = model_config_rec% restart_interval_d 
  restart_interval_h         = model_config_rec% restart_interval_h 
  restart_interval_m         = model_config_rec% restart_interval_m 
  restart_interval_s         = model_config_rec% restart_interval_s 
  history_begin_y            = model_config_rec% history_begin_y 
  history_begin_d            = model_config_rec% history_begin_d 
  history_begin_h            = model_config_rec% history_begin_h 
  history_begin_m            = model_config_rec% history_begin_m 
  history_begin_s            = model_config_rec% history_begin_s 
  history_begin              = model_config_rec% history_begin 
  inputout_begin_y           = model_config_rec% inputout_begin_y 
  inputout_begin_d           = model_config_rec% inputout_begin_d 
  inputout_begin_h           = model_config_rec% inputout_begin_h 
  inputout_begin_m           = model_config_rec% inputout_begin_m 
  inputout_begin_s           = model_config_rec% inputout_begin_s 
  restart_begin_y            = model_config_rec% restart_begin_y 
  restart_begin_d            = model_config_rec% restart_begin_d 
  restart_begin_h            = model_config_rec% restart_begin_h 
  restart_begin_m            = model_config_rec% restart_begin_m 
  restart_begin_s            = model_config_rec% restart_begin_s 
  restart_begin              = model_config_rec% restart_begin 
  history_end_y              = model_config_rec% history_end_y 
  history_end_d              = model_config_rec% history_end_d 
  history_end_h              = model_config_rec% history_end_h 
  history_end_m              = model_config_rec% history_end_m 
  history_end_s              = model_config_rec% history_end_s 
  history_end                = model_config_rec% history_end 
  inputout_end_y             = model_config_rec% inputout_end_y 
  inputout_end_d             = model_config_rec% inputout_end_d 
  inputout_end_h             = model_config_rec% inputout_end_h 
  inputout_end_m             = model_config_rec% inputout_end_m 
  inputout_end_s             = model_config_rec% inputout_end_s 
  simulation_start_year      = model_config_rec% simulation_start_year 
  simulation_start_month     = model_config_rec% simulation_start_month 
  simulation_start_day       = model_config_rec% simulation_start_day 
  simulation_start_hour      = model_config_rec% simulation_start_hour 
  simulation_start_minute    = model_config_rec% simulation_start_minute 
  simulation_start_second    = model_config_rec% simulation_start_second 
  reset_simulation_start     = model_config_rec% reset_simulation_start 
  sr_x                       = model_config_rec% sr_x 
  sr_y                       = model_config_rec% sr_y 
  sgfdda_inname              = model_config_rec% sgfdda_inname 
  gfdda_inname               = model_config_rec% gfdda_inname 
  sgfdda_interval_d          = model_config_rec% sgfdda_interval_d 
  sgfdda_interval_h          = model_config_rec% sgfdda_interval_h 
  sgfdda_interval_m          = model_config_rec% sgfdda_interval_m 
  sgfdda_interval_s          = model_config_rec% sgfdda_interval_s 
  sgfdda_interval_y          = model_config_rec% sgfdda_interval_y 
  sgfdda_interval            = model_config_rec% sgfdda_interval 
  gfdda_interval_d           = model_config_rec% gfdda_interval_d 
  gfdda_interval_h           = model_config_rec% gfdda_interval_h 
  gfdda_interval_m           = model_config_rec% gfdda_interval_m 
  gfdda_interval_s           = model_config_rec% gfdda_interval_s 
  gfdda_interval_y           = model_config_rec% gfdda_interval_y 
  gfdda_interval             = model_config_rec% gfdda_interval 
  sgfdda_begin_y             = model_config_rec% sgfdda_begin_y 
  sgfdda_begin_d             = model_config_rec% sgfdda_begin_d 
  sgfdda_begin_h             = model_config_rec% sgfdda_begin_h 
  sgfdda_begin_m             = model_config_rec% sgfdda_begin_m 
  sgfdda_begin_s             = model_config_rec% sgfdda_begin_s 
  gfdda_begin_y              = model_config_rec% gfdda_begin_y 
  gfdda_begin_d              = model_config_rec% gfdda_begin_d 
  gfdda_begin_h              = model_config_rec% gfdda_begin_h 
  gfdda_begin_m              = model_config_rec% gfdda_begin_m 
  gfdda_begin_s              = model_config_rec% gfdda_begin_s 
  sgfdda_end_y               = model_config_rec% sgfdda_end_y 
  sgfdda_end_d               = model_config_rec% sgfdda_end_d 
  sgfdda_end_h               = model_config_rec% sgfdda_end_h 
  sgfdda_end_m               = model_config_rec% sgfdda_end_m 
  sgfdda_end_s               = model_config_rec% sgfdda_end_s 
  gfdda_end_y                = model_config_rec% gfdda_end_y 
  gfdda_end_d                = model_config_rec% gfdda_end_d 
  gfdda_end_h                = model_config_rec% gfdda_end_h 
  gfdda_end_m                = model_config_rec% gfdda_end_m 
  gfdda_end_s                = model_config_rec% gfdda_end_s 
  io_form_sgfdda             = model_config_rec% io_form_sgfdda 
  io_form_gfdda              = model_config_rec% io_form_gfdda 
  iofields_filename          = model_config_rec% iofields_filename 
  ignore_iofields_warning    = model_config_rec% ignore_iofields_warning 
  ncd_nofill                 = model_config_rec% ncd_nofill 
  update_sfcdiags            = model_config_rec% update_sfcdiags 
  use_wrf_sfcinfo            = model_config_rec% use_wrf_sfcinfo 
  use_background_errors      = model_config_rec% use_background_errors 
  write_increments           = model_config_rec% write_increments 
  var4d                      = model_config_rec% var4d 
  var4d_bin                  = model_config_rec% var4d_bin 
  var4d_bin_rain             = model_config_rec% var4d_bin_rain 
  var4d_lbc                  = model_config_rec% var4d_lbc 
  multi_inc                  = model_config_rec% multi_inc 
  print_detail_radar         = model_config_rec% print_detail_radar 
  print_detail_rain          = model_config_rec% print_detail_rain 
  print_detail_rad           = model_config_rec% print_detail_rad 
  print_detail_xa            = model_config_rec% print_detail_xa 
  print_detail_xb            = model_config_rec% print_detail_xb 
  print_detail_obs           = model_config_rec% print_detail_obs 
  print_detail_f_obs         = model_config_rec% print_detail_f_obs 
  print_detail_map           = model_config_rec% print_detail_map 
  print_detail_grad          = model_config_rec% print_detail_grad 
  print_detail_regression    = model_config_rec% print_detail_regression 
  print_detail_spectral      = model_config_rec% print_detail_spectral 
  print_detail_testing       = model_config_rec% print_detail_testing 
  print_detail_parallel      = model_config_rec% print_detail_parallel 
  print_detail_be            = model_config_rec% print_detail_be 
  print_detail_outerloop     = model_config_rec% print_detail_outerloop 
  check_max_iv_print         = model_config_rec% check_max_iv_print 
  check_buddy_print          = model_config_rec% check_buddy_print 
  analysis_accu              = model_config_rec% analysis_accu 
  calc_w_increment           = model_config_rec% calc_w_increment 
  dt_cloud_model             = model_config_rec% dt_cloud_model 
  write_mod_filtered_obs     = model_config_rec% write_mod_filtered_obs 
  wind_sd                    = model_config_rec% wind_sd 
  wind_sd_buoy               = model_config_rec% wind_sd_buoy 
  wind_sd_synop              = model_config_rec% wind_sd_synop 
  wind_sd_ships              = model_config_rec% wind_sd_ships 
  wind_sd_metar              = model_config_rec% wind_sd_metar 
  wind_sd_sound              = model_config_rec% wind_sd_sound 
  wind_sd_pilot              = model_config_rec% wind_sd_pilot 
  wind_sd_airep              = model_config_rec% wind_sd_airep 
  wind_sd_qscat              = model_config_rec% wind_sd_qscat 
  wind_sd_tamdar             = model_config_rec% wind_sd_tamdar 
  wind_sd_geoamv             = model_config_rec% wind_sd_geoamv 
  wind_sd_mtgirs             = model_config_rec% wind_sd_mtgirs 
  wind_sd_polaramv           = model_config_rec% wind_sd_polaramv 
  wind_sd_profiler           = model_config_rec% wind_sd_profiler 
  wind_stats_sd              = model_config_rec% wind_stats_sd 
  qc_rej_both                = model_config_rec% qc_rej_both 
  fg_format                  = model_config_rec% fg_format 
  ob_format                  = model_config_rec% ob_format 
  ob_format_gpsro            = model_config_rec% ob_format_gpsro 
  num_fgat_time              = model_config_rec% num_fgat_time 
  thin_conv                  = model_config_rec% thin_conv 
  thin_conv_ascii            = model_config_rec% thin_conv_ascii 
  thin_mesh_conv             = model_config_rec% thin_mesh_conv 
  thin_rainobs               = model_config_rec% thin_rainobs 
  use_synopobs               = model_config_rec% use_synopobs 
  use_shipsobs               = model_config_rec% use_shipsobs 
  use_metarobs               = model_config_rec% use_metarobs 
  use_soundobs               = model_config_rec% use_soundobs 
  use_mtgirsobs              = model_config_rec% use_mtgirsobs 
  use_tamdarobs              = model_config_rec% use_tamdarobs 
  use_pilotobs               = model_config_rec% use_pilotobs 
  use_airepobs               = model_config_rec% use_airepobs 
  use_geoamvobs              = model_config_rec% use_geoamvobs 
  use_polaramvobs            = model_config_rec% use_polaramvobs 
  use_bogusobs               = model_config_rec% use_bogusobs 
  use_buoyobs                = model_config_rec% use_buoyobs 
  use_profilerobs            = model_config_rec% use_profilerobs 
  use_satemobs               = model_config_rec% use_satemobs 
  use_gpsztdobs              = model_config_rec% use_gpsztdobs 
  use_gpspwobs               = model_config_rec% use_gpspwobs 
  use_gpsrefobs              = model_config_rec% use_gpsrefobs 
  top_km_gpsro               = model_config_rec% top_km_gpsro 
  bot_km_gpsro               = model_config_rec% bot_km_gpsro 
  use_ssmiretrievalobs       = model_config_rec% use_ssmiretrievalobs 
  use_ssmitbobs              = model_config_rec% use_ssmitbobs 
  use_ssmt1obs               = model_config_rec% use_ssmt1obs 
  use_ssmt2obs               = model_config_rec% use_ssmt2obs 
  use_qscatobs               = model_config_rec% use_qscatobs 
  use_radarobs               = model_config_rec% use_radarobs 
  use_radar_rv               = model_config_rec% use_radar_rv 
  use_radar_rf               = model_config_rec% use_radar_rf 
  use_radar_rqv              = model_config_rec% use_radar_rqv 
  use_radar_rhv              = model_config_rec% use_radar_rhv 
  use_3dvar_phy              = model_config_rec% use_3dvar_phy 
  use_rainobs                = model_config_rec% use_rainobs 
  use_hirs2obs               = model_config_rec% use_hirs2obs 
  use_hirs3obs               = model_config_rec% use_hirs3obs 
  use_hirs4obs               = model_config_rec% use_hirs4obs 
  use_mhsobs                 = model_config_rec% use_mhsobs 
  use_msuobs                 = model_config_rec% use_msuobs 
  use_amsuaobs               = model_config_rec% use_amsuaobs 
  use_amsubobs               = model_config_rec% use_amsubobs 
  use_airsobs                = model_config_rec% use_airsobs 
  use_airsretobs             = model_config_rec% use_airsretobs 
  use_eos_amsuaobs           = model_config_rec% use_eos_amsuaobs 
  use_hsbobs                 = model_config_rec% use_hsbobs 
  use_ssmisobs               = model_config_rec% use_ssmisobs 
  use_iasiobs                = model_config_rec% use_iasiobs 
  use_seviriobs              = model_config_rec% use_seviriobs 
  use_amsr2obs               = model_config_rec% use_amsr2obs 
  use_kma1dvar               = model_config_rec% use_kma1dvar 
  use_filtered_rad           = model_config_rec% use_filtered_rad 
  use_obs_errfac             = model_config_rec% use_obs_errfac 
  use_atmsobs                = model_config_rec% use_atmsobs 
  use_mwtsobs                = model_config_rec% use_mwtsobs 
  use_mwhsobs                = model_config_rec% use_mwhsobs 
  check_max_iv               = model_config_rec% check_max_iv 
  max_error_t                = model_config_rec% max_error_t 
  max_error_uv               = model_config_rec% max_error_uv 
  max_error_spd              = model_config_rec% max_error_spd 
  max_error_dir              = model_config_rec% max_error_dir 
  max_omb_spd                = model_config_rec% max_omb_spd 
  max_omb_dir                = model_config_rec% max_omb_dir 
  max_error_pw               = model_config_rec% max_error_pw 
  max_error_ref              = model_config_rec% max_error_ref 
  max_error_rh               = model_config_rec% max_error_rh 
  max_error_q                = model_config_rec% max_error_q 
  max_error_p                = model_config_rec% max_error_p 
  max_error_tb               = model_config_rec% max_error_tb 
  max_error_thickness        = model_config_rec% max_error_thickness 
  max_error_rv               = model_config_rec% max_error_rv 
  max_error_rf               = model_config_rec% max_error_rf 
  max_error_rain             = model_config_rec% max_error_rain 
  max_error_buv              = model_config_rec% max_error_buv 
  max_error_bt               = model_config_rec% max_error_bt 
  max_error_bq               = model_config_rec% max_error_bq 
  max_error_slp              = model_config_rec% max_error_slp 
  check_buddy                = model_config_rec% check_buddy 
  put_rand_seed              = model_config_rec% put_rand_seed 
  omb_set_rand               = model_config_rec% omb_set_rand 
  omb_add_noise              = model_config_rec% omb_add_noise 
  position_lev_dependant     = model_config_rec% position_lev_dependant 
  obs_qc_pointer             = model_config_rec% obs_qc_pointer 
  qmarker_retain             = model_config_rec% qmarker_retain 
  max_sound_input            = model_config_rec% max_sound_input 
  max_mtgirs_input           = model_config_rec% max_mtgirs_input 
  max_tamdar_input           = model_config_rec% max_tamdar_input 
  max_synop_input            = model_config_rec% max_synop_input 
  max_geoamv_input           = model_config_rec% max_geoamv_input 
  max_polaramv_input         = model_config_rec% max_polaramv_input 
  max_airep_input            = model_config_rec% max_airep_input 
  max_satem_input            = model_config_rec% max_satem_input 
  max_pilot_input            = model_config_rec% max_pilot_input 
  max_radar_input            = model_config_rec% max_radar_input 
  max_rain_input             = model_config_rec% max_rain_input 
  max_metar_input            = model_config_rec% max_metar_input 
  max_gpspw_input            = model_config_rec% max_gpspw_input 
  max_ships_input            = model_config_rec% max_ships_input 
  max_profiler_input         = model_config_rec% max_profiler_input 
  max_bogus_input            = model_config_rec% max_bogus_input 
  max_buoy_input             = model_config_rec% max_buoy_input 
  max_ssmi_rv_input          = model_config_rec% max_ssmi_rv_input 
  max_ssmi_tb_input          = model_config_rec% max_ssmi_tb_input 
  max_ssmt1_input            = model_config_rec% max_ssmt1_input 
  max_ssmt2_input            = model_config_rec% max_ssmt2_input 
  max_qscat_input            = model_config_rec% max_qscat_input 
  max_gpsref_input           = model_config_rec% max_gpsref_input 
  max_airsr_input            = model_config_rec% max_airsr_input 
  max_tovs_input             = model_config_rec% max_tovs_input 
  max_ssmis_input            = model_config_rec% max_ssmis_input 
  report_start               = model_config_rec% report_start 
  report_end                 = model_config_rec% report_end 
  tovs_start                 = model_config_rec% tovs_start 
  tovs_end                   = model_config_rec% tovs_end 
  gpsref_thinning            = model_config_rec% gpsref_thinning 
  outer_loop_restart         = model_config_rec% outer_loop_restart 
  max_ext_its                = model_config_rec% max_ext_its 
  ntmax                      = model_config_rec% ntmax 
  nsave                      = model_config_rec% nsave 
  write_interval             = model_config_rec% write_interval 
  eps                        = model_config_rec% eps 
  precondition_cg            = model_config_rec% precondition_cg 
  precondition_factor        = model_config_rec% precondition_factor 
  use_lanczos                = model_config_rec% use_lanczos 
  read_lanczos               = model_config_rec% read_lanczos 
  write_lanczos              = model_config_rec% write_lanczos 
  orthonorm_gradient         = model_config_rec% orthonorm_gradient 
  cv_options                 = model_config_rec% cv_options 
  cloud_cv_options           = model_config_rec% cloud_cv_options 
  as1                        = model_config_rec% as1 
  as2                        = model_config_rec% as2 
  as3                        = model_config_rec% as3 
  as4                        = model_config_rec% as4 
  as5                        = model_config_rec% as5 
  do_normalize               = model_config_rec% do_normalize 
  use_rf                     = model_config_rec% use_rf 
  rf_passes                  = model_config_rec% rf_passes 
  var_scaling1               = model_config_rec% var_scaling1 
  var_scaling2               = model_config_rec% var_scaling2 
  var_scaling3               = model_config_rec% var_scaling3 
  var_scaling4               = model_config_rec% var_scaling4 
  var_scaling5               = model_config_rec% var_scaling5 
  var_scaling6               = model_config_rec% var_scaling6 
  var_scaling7               = model_config_rec% var_scaling7 
  var_scaling8               = model_config_rec% var_scaling8 
  var_scaling9               = model_config_rec% var_scaling9 
  var_scaling10              = model_config_rec% var_scaling10 
  var_scaling11              = model_config_rec% var_scaling11 
  len_scaling1               = model_config_rec% len_scaling1 
  len_scaling2               = model_config_rec% len_scaling2 
  len_scaling3               = model_config_rec% len_scaling3 
  len_scaling4               = model_config_rec% len_scaling4 
  len_scaling5               = model_config_rec% len_scaling5 
  len_scaling6               = model_config_rec% len_scaling6 
  len_scaling7               = model_config_rec% len_scaling7 
  len_scaling8               = model_config_rec% len_scaling8 
  len_scaling9               = model_config_rec% len_scaling9 
  len_scaling10              = model_config_rec% len_scaling10 
  len_scaling11              = model_config_rec% len_scaling11 
  je_factor                  = model_config_rec% je_factor 
  power_truncation           = model_config_rec% power_truncation 
  def_sub_domain             = model_config_rec% def_sub_domain 
  x_start_sub_domain         = model_config_rec% x_start_sub_domain 
  y_start_sub_domain         = model_config_rec% y_start_sub_domain 
  x_end_sub_domain           = model_config_rec% x_end_sub_domain 
  y_end_sub_domain           = model_config_rec% y_end_sub_domain 
  stdout                     = model_config_rec% stdout 
  stderr                     = model_config_rec% stderr 
  trace_unit                 = model_config_rec% trace_unit 
  trace_pe                   = model_config_rec% trace_pe 
  trace_repeat_head          = model_config_rec% trace_repeat_head 
  trace_repeat_body          = model_config_rec% trace_repeat_body 
  trace_max_depth            = model_config_rec% trace_max_depth 
  trace_use                  = model_config_rec% trace_use 
  trace_use_frequent         = model_config_rec% trace_use_frequent 
  trace_use_dull             = model_config_rec% trace_use_dull 
  trace_memory               = model_config_rec% trace_memory 
  trace_all_pes              = model_config_rec% trace_all_pes 
  trace_csv                  = model_config_rec% trace_csv 
  use_html                   = model_config_rec% use_html 
  warnings_are_fatal         = model_config_rec% warnings_are_fatal 
  test_transforms            = model_config_rec% test_transforms 
  test_gradient              = model_config_rec% test_gradient 
  test_statistics            = model_config_rec% test_statistics 
  interpolate_stats          = model_config_rec% interpolate_stats 
  be_eta                     = model_config_rec% be_eta 
  test_dm_exact              = model_config_rec% test_dm_exact 
  cv_options_hum             = model_config_rec% cv_options_hum 
  check_rh                   = model_config_rec% check_rh 
  set_omb_rand_fac           = model_config_rec% set_omb_rand_fac 
  seed_array1                = model_config_rec% seed_array1 
  seed_array2                = model_config_rec% seed_array2 
  sfc_assi_options           = model_config_rec% sfc_assi_options 
  psfc_from_slp              = model_config_rec% psfc_from_slp 
  calculate_cg_cost_fn       = model_config_rec% calculate_cg_cost_fn 
  lat_stats_option           = model_config_rec% lat_stats_option 
  interp_option              = model_config_rec% interp_option 
  balance_type               = model_config_rec% balance_type 
  use_wpec                   = model_config_rec% use_wpec 
  wpec_factor                = model_config_rec% wpec_factor 
  vert_corr                  = model_config_rec% vert_corr 
  vertical_ip                = model_config_rec% vertical_ip 
  vert_evalue                = model_config_rec% vert_evalue 
  max_vert_var1              = model_config_rec% max_vert_var1 
  max_vert_var2              = model_config_rec% max_vert_var2 
  max_vert_var3              = model_config_rec% max_vert_var3 
  max_vert_var4              = model_config_rec% max_vert_var4 
  max_vert_var5              = model_config_rec% max_vert_var5 
  max_vert_var6              = model_config_rec% max_vert_var6 
  max_vert_var7              = model_config_rec% max_vert_var7 
  max_vert_var8              = model_config_rec% max_vert_var8 
  max_vert_var9              = model_config_rec% max_vert_var9 
  max_vert_var10             = model_config_rec% max_vert_var10 
  max_vert_var11             = model_config_rec% max_vert_var11 
  max_vert_var_alpha         = model_config_rec% max_vert_var_alpha 
  psi_chi_factor             = model_config_rec% psi_chi_factor 
  psi_t_factor               = model_config_rec% psi_t_factor 
  psi_ps_factor              = model_config_rec% psi_ps_factor 
  psi_rh_factor              = model_config_rec% psi_rh_factor 
  chi_u_t_factor             = model_config_rec% chi_u_t_factor 
  chi_u_ps_factor            = model_config_rec% chi_u_ps_factor 
  chi_u_rh_factor            = model_config_rec% chi_u_rh_factor 
  t_u_rh_factor              = model_config_rec% t_u_rh_factor 
  ps_u_rh_factor             = model_config_rec% ps_u_rh_factor 
  rttov_emis_atlas_ir        = model_config_rec% rttov_emis_atlas_ir 
  rttov_emis_atlas_mw        = model_config_rec% rttov_emis_atlas_mw 
  rtminit_print              = model_config_rec% rtminit_print 
  rtminit_nsensor            = model_config_rec% rtminit_nsensor 
  rtminit_platform           = model_config_rec% rtminit_platform 
  rtminit_satid              = model_config_rec% rtminit_satid 
  rtminit_sensor             = model_config_rec% rtminit_sensor 
  rad_monitoring             = model_config_rec% rad_monitoring 
  thinning_mesh              = model_config_rec% thinning_mesh 
  thinning                   = model_config_rec% thinning 
  read_biascoef              = model_config_rec% read_biascoef 
  biascorr                   = model_config_rec% biascorr 
  biasprep                   = model_config_rec% biasprep 
  rttov_scatt                = model_config_rec% rttov_scatt 
  write_profile              = model_config_rec% write_profile 
  write_jacobian             = model_config_rec% write_jacobian 
  qc_rad                     = model_config_rec% qc_rad 
  write_iv_rad_ascii         = model_config_rec% write_iv_rad_ascii 
  write_oa_rad_ascii         = model_config_rec% write_oa_rad_ascii 
  write_filtered_rad         = model_config_rec% write_filtered_rad 
  use_error_factor_rad       = model_config_rec% use_error_factor_rad 
  use_landem                 = model_config_rec% use_landem 
  use_antcorr                = model_config_rec% use_antcorr 
  use_mspps_emis             = model_config_rec% use_mspps_emis 
  use_mspps_ts               = model_config_rec% use_mspps_ts 
  mw_emis_sea                = model_config_rec% mw_emis_sea 
  tovs_min_transfer          = model_config_rec% tovs_min_transfer 
  tovs_batch                 = model_config_rec% tovs_batch 
  rtm_option                 = model_config_rec% rtm_option 
  use_crtm_kmatrix           = model_config_rec% use_crtm_kmatrix 
  use_rttov_kmatrix          = model_config_rec% use_rttov_kmatrix 
  crtm_cloud                 = model_config_rec% crtm_cloud 
  only_sea_rad               = model_config_rec% only_sea_rad 
  use_pseudo_rad             = model_config_rec% use_pseudo_rad 
  pseudo_rad_platid          = model_config_rec% pseudo_rad_platid 
  pseudo_rad_satid           = model_config_rec% pseudo_rad_satid 
  pseudo_rad_senid           = model_config_rec% pseudo_rad_senid 
  pseudo_rad_ichan           = model_config_rec% pseudo_rad_ichan 
  pseudo_rad_lat             = model_config_rec% pseudo_rad_lat 
  pseudo_rad_lon             = model_config_rec% pseudo_rad_lon 
  pseudo_rad_inv             = model_config_rec% pseudo_rad_inv 
  pseudo_rad_err             = model_config_rec% pseudo_rad_err 
  use_simulated_rad          = model_config_rec% use_simulated_rad 
  simulated_rad_io           = model_config_rec% simulated_rad_io 
  simulated_rad_ngrid        = model_config_rec% simulated_rad_ngrid 
  use_varbc                  = model_config_rec% use_varbc 
  freeze_varbc               = model_config_rec% freeze_varbc 
  varbc_factor               = model_config_rec% varbc_factor 
  varbc_nbgerr               = model_config_rec% varbc_nbgerr 
  varbc_nobsmin              = model_config_rec% varbc_nobsmin 
  use_clddet_mmr             = model_config_rec% use_clddet_mmr 
  use_clddet_ecmwf           = model_config_rec% use_clddet_ecmwf 
  airs_warmest_fov           = model_config_rec% airs_warmest_fov 
  use_satcv                  = model_config_rec% use_satcv 
  use_blacklist_rad          = model_config_rec% use_blacklist_rad 
  calc_weightfunc            = model_config_rec% calc_weightfunc 
  crtm_coef_path             = model_config_rec% crtm_coef_path 
  crtm_irwater_coef          = model_config_rec% crtm_irwater_coef 
  crtm_mwwater_coef          = model_config_rec% crtm_mwwater_coef 
  crtm_irland_coef           = model_config_rec% crtm_irland_coef 
  crtm_visland_coef          = model_config_rec% crtm_visland_coef 
  num_pseudo                 = model_config_rec% num_pseudo 
  pseudo_x                   = model_config_rec% pseudo_x 
  pseudo_y                   = model_config_rec% pseudo_y 
  pseudo_z                   = model_config_rec% pseudo_z 
  pseudo_val                 = model_config_rec% pseudo_val 
  pseudo_err                 = model_config_rec% pseudo_err 
  alphacv_method             = model_config_rec% alphacv_method 
  ensdim_alpha               = model_config_rec% ensdim_alpha 
  alpha_truncation           = model_config_rec% alpha_truncation 
  alpha_corr_type            = model_config_rec% alpha_corr_type 
  alpha_corr_scale           = model_config_rec% alpha_corr_scale 
  alpha_std_dev              = model_config_rec% alpha_std_dev 
  alpha_vertloc              = model_config_rec% alpha_vertloc 
  alpha_hydrometeors         = model_config_rec% alpha_hydrometeors 
  hybrid_dual_res            = model_config_rec% hybrid_dual_res 
  dual_res_upscale_opt       = model_config_rec% dual_res_upscale_opt 
  analysis_type              = model_config_rec% analysis_type 
  sensitivity_option         = model_config_rec% sensitivity_option 
  adj_sens                   = model_config_rec% adj_sens 
  analysis_date              = model_config_rec% analysis_date 
  pseudo_var                 = model_config_rec% pseudo_var 
  documentation_url          = model_config_rec% documentation_url 
  time_window_min            = model_config_rec% time_window_min 
  time_window_max            = model_config_rec% time_window_max 
  jcdfi_use                  = model_config_rec% jcdfi_use 
  jcdfi_diag                 = model_config_rec% jcdfi_diag 
  jcdfi_penalty              = model_config_rec% jcdfi_penalty 
  enable_identity            = model_config_rec% enable_identity 
  trajectory_io              = model_config_rec% trajectory_io 
  var4d_detail_out           = model_config_rec% var4d_detail_out 
  var4d_run                  = model_config_rec% var4d_run 
  mp_physics_ad              = model_config_rec% mp_physics_ad 
  mp_physics_4dvar           = model_config_rec% mp_physics_4dvar 
  chem_opt                   = model_config_rec% chem_opt 
!ENDOFREGISTRYGENERATEDINCLUDE

   if (.NOT. rootproc) check_max_iv_print=.false.

end subroutine da_wrfvar_init1


subroutine da_wrfvar_init2

   !-------------------------------------------------------------------------
   ! Purpose: WRFVAR initialization routine, part 2
   !-------------------------------------------------------------------------

   implicit none

   integer :: i
   character(len=80) :: filename
   logical           :: isfile
   logical           :: ex


   if (trace_use) call da_trace_entry("da_wrfvar_init2")

! Override the start time with the "analysis_date":
      read(analysis_date, fmt='(i4,5(1x,i2))') &
           start_year(1), start_month(1), start_day(1), start_hour(1), &
           start_minute(1), start_second(1)
      model_config_rec% start_year   = start_year
      model_config_rec% start_month  = start_month
      model_config_rec% start_day    = start_day
      model_config_rec% start_hour   = start_hour
      model_config_rec% start_minute = start_minute
      model_config_rec% start_second = start_second

   if (analysis_type(1:6) == "VERIFY" .or. analysis_type(1:6) == "verify") then
      anal_type_verify=.true.
   else
      anal_type_verify=.false.
   end if

   if (analysis_type(1:8) == "RANDOMCV" .or. analysis_type(1:8) == "randomcv") then
      anal_type_randomcv=.true.
   else
      anal_type_randomcv=.false.
   end if

   if (analysis_type(1:6) == "QC-OBS" .or. analysis_type(1:6) == "qc-obs") then
      anal_type_qcobs=.true.
   else
      anal_type_qcobs=.false.
   end if

   if (use_gpspwObs .and. use_gpsztdObs ) then
      call da_error("da_wrfvar_init2.inc",47, (/'can not assimilate gpspw and gpsztd simultaneously'/))
   end if

   if (fg_format==fg_format_kma_global .or. fg_format==fg_format_wrf_arw_global) then
      global = .true.
      nproc_x = 1
   else
      global = .false.
   end if

   anal_type_hybrid_dual_res = .false.
   if ( hybrid_dual_res ) then
      if ( ensdim_alpha >= 1 ) then
         if ( max_dom /= 2 ) then
            call da_error("da_wrfvar_init2.inc",61, (/'max_dom has to be 2 for hybrid_dual_res application'/))
         end if
         anal_type_hybrid_dual_res = .true.
      else
         write(unit=message(1),fmt='(A,2(I4))') "ensdim_alpha has to be non-zero for hybrid_dual_res application, resetting hybrid_dual_res=.false."
         call da_warning("da_wrfvar_init2.inc",66,message(1:1))
         anal_type_hybrid_dual_res = .false.
      end if
   endif

   if ( anal_type_hybrid_dual_res ) then
      call nl_set_shw( 1 , 0 )
      call nl_set_shw( 2 , 0 )
      !write(unit=message(1),fmt='(A,2(I4))') "Resetting shw for dual-res hybrid to shw = ",shw(1),shw(2)
      write(unit=message(1),fmt='(A,2(I4))') "Running WRFDA in dual-resolution hybrid mode"
      call da_message(message(1:1))
   endif


   !<DESCRIPTION>
   ! Among the configuration variables read from the namelist is
   ! debug_level. This is retrieved using nl_get_debug_level (Registry
   ! generated and defined in frame/module_configure.F).  The value is then
   ! used to set the debug-print information level for use by <a
   ! href=wrf_debug.html>wrf_debug</a> throughout the code. Debug_level
   ! of zero (the default) causes no information to be printed when the
   ! model runs. The higher the number (up to 1000) the more information is
   ! printed.
   ! 
   !</DESCRIPTION>

   call nl_get_debug_level (1, debug_level)
   call set_wrf_debug_level (debug_level)

   nullify(null_domain)


!  if (max_dom > 1 ) then 
   if (max_dom > 1 .and. ( .not. anal_type_hybrid_dual_res) ) then
      call da_error("da_wrfvar_init2.inc",100, (/'WRFDA does not handle nests (max_domain > 1)'/))
   end if

   !<DESCRIPTION>
   ! The top-most domain in the simulation is then allocated and configured
   ! by calling <a href=alloc_and_configure_domain.html>alloc_and_configure_domain</a>.
   ! Here, in the case of this root domain, the routine is passed the
   ! globally accessible pointer to type(domain), head_grid, defined in
   ! frame/module_domain.F.  The parent is null and the child index is given
   ! as negative, signifying none.  Afterwards, because the call to
   ! alloc_and_configure_domain may modify the model configuration data
   ! stored in model_config_rec, the configuration information is again
   ! repacked into a buffer, broadcast, and unpacked on each task (for
   ! 1 compiles). The call to <a
   ! href=setup_timekeeping.html>setup_timekeeping</a> for head_grid relies
   ! on this configuration information, and it must occur after the second
   ! broadcast of the configuration information.
   ! 
   !</DESCRIPTION>

   call da_trace("da_wrfvar_init2",message="calling alloc_and_configure_domain")

   call alloc_and_configure_domain (domain_id=1, grid=head_grid, parent=null_domain, kid=-1)  

   call da_trace("da_wrfvar_init2",message="calling model_to_grid_config_rec")
   call model_to_grid_config_rec (head_grid%id, model_config_rec, config_flags)  

   call da_trace("da_wrfvar_init2",message="calling set_scalar_indices_from_config")
   call set_scalar_indices_from_config (head_grid%id , idum1, idum2) 

   call da_trace("da_wrfvar_init2",message="calling init_wrfio")
   call init_wrfio

   call get_config_as_buffer (configbuf, configbuflen, nbytes)
   call wrf_dm_bcast_bytes (configbuf, nbytes)
   call set_config_as_buffer (configbuf, configbuflen)

   call setup_timekeeping (head_grid) 

   if ( anal_type_hybrid_dual_res ) then
      ! input_file_ens is 'fg_ens', set in da_control.f90
      inquire(file=trim(input_file_ens), exist=isfile)
      if ( .not. isfile ) then
         write(unit=message(1),fmt='(a,a,a)') 'File ',trim(input_file_ens),' (low-resolution ensemble file) is missing.'
         call da_error("da_wrfvar_init2.inc",146,message(1:1))
      endif
      call da_med_initialdata_input (head_grid, config_flags, trim(input_file_ens))
      parent_grid => head_grid
      call alloc_and_configure_domain (domain_id=2, grid=another_grid, parent=parent_grid, kid=1)
      call model_to_grid_config_rec (another_grid%id, model_config_rec, config_flags)
      call set_scalar_indices_from_config (another_grid%id , idum1, idum2)
      call init_wrfio
	 call get_config_as_buffer (configbuf, configbuflen, nbytes)
	 call wrf_dm_bcast_bytes (configbuf, nbytes)
	 call set_config_as_buffer (configbuf, configbuflen)
      call setup_timekeeping (another_grid)

      input_grid => another_grid
      ensemble_grid => head_grid
   else
      input_grid => head_grid 
      ensemble_grid => head_grid
   endif

   !<DESCRIPTION>
   ! The head grid is initialized with read-in data through the call to <a
   ! href=med_initialdata_input.html>med_initialdata_input</a>, which is
   ! passed the pointer head_grid and a locally declared configuration data
   ! structure, config_flags, that is set by a call to <a
   ! href=model_to_grid_config_rec.html>model_to_grid_config_rec</a>.  It is
   ! also necessary that the indices into the 4d tracer arrays such as
   ! moisture be set with a call to <a
   ! href=set_scalar_indices_from_config.html>set_scalar_indices_from_config</a>
   ! prior to the call to initialize the domain.  Both of these calls are
   ! told which domain they are setting up for by passing in the integer id
   ! of the head domain as <tt>head_grid%id</tt>, which is 1 for the
   ! top-most domain.
   ! 
   ! In the case that write_restart_at_0h is set to true in the namelist,
   ! the model simply generates a restart file using the just read-in data
   ! and then shuts down. This is used for ensemble breeding, and is not
   ! typically enabled.
   ! 
   !</DESCRIPTION>

   ! call med_initialdata_input(head_grid , config_flags,'fg')

   if ((config_flags%real_data_init_type == 1) .or. &
       (config_flags%real_data_init_type == 3)) then
!     call da_med_initialdata_input (head_grid, config_flags, 'fg')
      call da_med_initialdata_input (input_grid, config_flags, 'fg')
      if ( var4d ) then
         call med_latbound_in ( head_grid, config_flags )
         call close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )
      end if
   end if

   ! FIX?
   ! call da_warning("da_wrfvar_init2.inc",232,(/"Fix me"/))
   ! head_grid%start_subtime = head_grid%start_time
   ! head_grid%stop_subtime = head_grid%stop_time

   if (rootproc) then
      call da_get_unit (cost_unit)
      call da_get_unit (grad_unit)
      call da_get_unit (jo_unit)
      call da_get_unit (check_max_iv_unit)
      call da_get_unit (check_buddy_unit)
      open(unit=cost_unit,file="cost_fn",status="replace")
      open(unit=grad_unit,file="grad_fn",status="replace")
      if (.not. print_detail_outerloop) then
         call da_get_unit (stats_unit)
         open(unit=stats_unit,file="statistics",status="replace")
      end if
      open(unit=jo_unit,file="jo",status="replace")
      open(unit=check_max_iv_unit,file="check_max_iv",status="replace")
      open(unit=check_buddy_unit ,file="buddy_check" ,status="replace")
   end if

   if (trace_use) call da_trace_exit("da_wrfvar_init2")

end subroutine da_wrfvar_init2


subroutine da_wrfvar_run()

   !-------------------------------------------------------------------------
   ! Purpose: run wrfvar
   !-------------------------------------------------------------------------

   implicit none

   if (trace_use) call da_trace_entry("da_wrfvar_run")

   call da_wrfvar_interface (input_grid, config_flags)

   if (trace_use) call da_trace_exit("da_wrfvar_run")

end subroutine da_wrfvar_run


subroutine da_wrfvar_interface (grid, config_flags) 

   !------------------------------------------------------------------------
   ! Purpose: TBD
   !------------------------------------------------------------------------

   implicit none

   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   ! integer :: idum1, idum2

   call da_trace_entry("da_wrfvar_interface")

   ! call mediation_setup_step (grid , config_flags , 1 , 1 , 1)

   call set_scalar_indices_from_config (grid%id , idum1 , idum2)

   call model_to_grid_config_rec (grid%id , model_config_rec , config_flags)

   grid%itimestep = 1

   call da_solve (grid , config_flags)

   call da_trace_exit("da_wrfvar_interface")

end subroutine da_wrfvar_interface


subroutine da_wrfvar_finalize

   !-------------------------------------------------------------------------
   ! Purpose: Tidy up at the end
   !-------------------------------------------------------------------------

   implicit none

   integer               :: i
   type(domain), pointer :: grid

   if (trace_use) call da_trace_entry ("da_wrfvar_finalize")


!  grid => head_grid
   grid => input_grid

   if ( .not. anal_type_hybrid_dual_res ) then
      call dealloc_space_domain(grid%id)
   else
      !=== deallocate fine grid
      ! subroutine dealloc_space_domain checks for head_grid
      ! when anal_type_hybrid_dual_res, input_grid is another_grid
      ! dealloc_space_domain does not work in this case
      call domain_destroy(grid)
      !=== deallocate coarse grid
      grid => ensemble_grid
      call dealloc_space_domain(grid%id)
   end if


   if (allocated(num_tovs_before)) deallocate (num_tovs_before)
   if (allocated(num_tovs_after))  deallocate (num_tovs_after)
   if (allocated(tovs_copy_count)) deallocate (tovs_copy_count)
   if (allocated(tovs_send_pe))    deallocate (tovs_send_pe)
   if (allocated(tovs_recv_pe))    deallocate (tovs_recv_pe)
   if (allocated(tovs_send_start)) deallocate (tovs_send_start)
   if (allocated(tovs_send_count)) deallocate (tovs_send_count)
   if (allocated(tovs_recv_start)) deallocate (tovs_recv_start)

   if (rootproc) then
      close (cost_unit)
      close (grad_unit)
      if (.not. print_detail_outerloop) then
         close (stats_unit)
         call da_free_unit (stats_unit)
      end if
      close (jo_unit)
      close (check_max_iv_unit)
      close (check_buddy_unit)
      call da_free_unit (cost_unit)
      call da_free_unit (grad_unit)
      call da_free_unit (jo_unit)
      call da_free_unit (check_max_iv_unit)
      call da_free_unit (check_buddy_unit )
   end if

   if (use_rad .and. rtm_option == rtm_option_crtm) then
      ierr = CRTM_Destroy(ChannelInfo)
      deallocate(Sensor_Descriptor)
   end if

   do i=unit_start,unit_end
      if (unit_used(i)) then
         write(unit=stderr,FMT=*) "unit",i,"still used"
      end if
   end do

   if (trace_use) call da_trace_exit ("da_wrfvar_finalize")

end subroutine da_wrfvar_finalize


  subroutine da_solve ( grid , config_flags)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! 
   ! Edited 09/06/2012: Allow for variable ntmax for each outer loop (Mike Kavulich)
   !-----------------------------------------------------------------------

   implicit none



   type (domain),               intent(inout) :: grid
   type (grid_config_rec_type), intent(inout) :: config_flags

   type (xbx_type)              :: xbx         ! For header & non-grid arrays.
   type (be_type)               :: be          ! Background error structure.
   real, allocatable            :: cvt(:)      ! Control variable structure.
   real, allocatable            :: xhat(:)     ! Control variable structure.
   real, allocatable            :: qhat(:,:)   ! Control variable structure.
   real*8, allocatable          :: eignvec(:,:)
   real*8, allocatable          :: eignval(:)   
!   real, allocatable            :: full_eignvec(:)   
   type (y_type)                :: ob          ! Observation structure.
   type (iv_type)               :: iv          ! Obs. increment structure.
   type (y_type)                :: re          ! Residual (o-a) structure.
   type (y_type)                :: y           ! y = H(x_inc) structure.
   integer                      :: it          ! External loop counter.
   integer                      :: neign       
   type (j_type)                :: j           ! Cost function.
   type (y_type)                :: jo_grad_y   ! Grad_y(jo)

   integer                      :: cv_size, i, ichan, k, n
   real                         :: j_grad_norm_target ! Target j norm.

   character (len=3)            :: ci
   character (len=2)            :: outerloop
   character (len=256)          :: timestr
   integer   :: min_yyyy,min_mm,min_dd,min_hh,min_mn,min_ss
   integer   :: max_yyyy,max_mm,max_dd,max_hh,max_mn,max_ss
   character :: s
   real*8    :: time_min, time_max 
   integer   :: jl_start, jl_end
   character(len=256) :: timestr1
   type(x_type) :: shuffle

   real, allocatable :: grid_box_area(:,:), mapfac(:,:)
   character (len=10)    :: variable_name

   integer :: xx, yy, ii, jj,  xy, kk, jjj
   real    :: lat, lon, x_pos, y_pos, dx, dy, dxm, dym
   logical :: outside

   integer   :: cvt_unit, iost
   character(len=8) :: cvtfile
   logical :: ex
 
   if (trace_use) call da_trace_entry("da_solve")

   call mpi_barrier(comm,ierr)

   if ( config_flags%use_baseparam_fr_nml ) then
      call nl_get_base_pres  ( 1 , base_pres )
      call nl_get_base_temp  ( 1 , base_temp )
      call nl_get_base_lapse ( 1 , base_lapse )
      call nl_get_iso_temp   ( 1 , iso_temp )
      if ( iso_temp .NE. grid%tiso ) THEN
         write(unit=message(1),fmt='(A)') &
           'Namelist iso_temp does not equal iso_temp from fg. Reset namelist value and rerun.'
         write(unit=message(2),fmt='(A,F10.5)')'Namelist iso_temp   = ',iso_temp
         write(unit=message(3),fmt='(A,F10.5)')'Background iso_temp = ',grid%tiso
         call da_error("da_solve.inc",74,message(1:3))
      end if
      call nl_get_base_pres_strat  ( 1, base_pres_strat )
      call nl_get_base_lapse_strat ( 1, base_lapse_strat )

      grid%p00   = base_pres
      grid%t00   = base_temp
      grid%tlp   = base_lapse
      grid%tiso  = iso_temp
      grid%p_strat   = base_pres_strat
      grid%tlp_strat = base_lapse_strat
   else
      base_pres  = grid%p00
      base_temp  = grid%t00
      base_lapse = grid%tlp
      iso_temp   = grid%tiso
      base_pres_strat  = grid%p_strat
      base_lapse_strat = grid%tlp_strat
      if ( base_temp < 100.0 .or. base_pres < 10000.0 ) then
         write(unit=message(1),fmt='(A)') &
         'did not find base state parameters in fg. Add use_baseparam_fr_nml = .t. in &dynamics and rerun '
         call da_error("da_solve.inc",95,message(1:1))
      end if
   end if

   ! Calculate the num_fgat_time based on time_window_min, time_window_max
    if ( var4d ) then
       if (time_step == 0) then
          write(unit=message(1),fmt='(A)') &
          'For 4DVAR, in the &domains namelist, "time_step" must be set to a non-zero value'
          call da_error("da_solve.inc",104,message(1:1))
       endif
       read(unit=time_window_min,fmt='(i4,5(a1,i2))') min_yyyy,s,min_mm,s,min_dd,s,min_hh,s,min_mn,s,min_ss
       read(unit=time_window_max,fmt='(i4,5(a1,i2))') max_yyyy,s,max_mm,s,max_dd,s,max_hh,s,max_mn,s,max_ss
       call da_get_julian_time(min_yyyy,min_mm,min_dd,min_hh,min_mn,time_min)
       call da_get_julian_time(max_yyyy,max_mm,max_dd,max_hh,max_mn,time_max)
       if ( var4d_bin < time_step ) call nl_set_var4d_bin (1, time_step)
       time_max = (time_max - time_min) * 60   ! unit is : seconds
       num_fgat_time = NINT(time_max/var4d_bin)
       if ( NINT(time_max/var4d_bin)*var4d_bin .ne. NINT(time_max) ) then
          write(unit=message(1),fmt='(A)') &
          '4DVAR assimilation window must be evenly divisible by var4d_bin!'
          write(unit=message(2),fmt='(A,I7)') &
          'var4d_bin       = ',var4d_bin
          write(unit=message(3),fmt='(A,A)') &
          'time_window_max = ',time_window_max
          write(unit=message(4),fmt='(A,A)') &
          'time_window_min = ',time_window_min
          write(unit=message(5),fmt='(A,F10.4)') &
          'time_window_max - time_window_min = ',time_max
          write(unit=message(6),fmt='(A)')'Change var4d_bin, time_window_max, or time_window_min in namelist and rerun'
          call da_error("da_solve.inc",125,message(1:6))
       endif
       if ( var4d_bin/time_step*time_step .ne. var4d_bin ) then
          write(unit=message(1),fmt='(A)') &
          'var4d_bin must be evenly divisible by time_step!'
          write(unit=message(2),fmt='(A,I7)') &
          'var4d_bin = ',var4d_bin
          write(unit=message(3),fmt='(A,I7)') &
          'time_step = ',time_step
          write(unit=message(4),fmt='(A)')'Change var4d_bin or time_step in namelist and rerun'
          call da_error("da_solve.inc",135,message(1:4))
       endif

       num_fgat_time = num_fgat_time + 1
       write(unit=message(1),fmt='(a,i10)') 'num_fgat_time is: ', num_fgat_time
       call da_message(message(1:1))
       if ( use_rainobs ) then
          allocate (fgat_rain_flags(1:num_fgat_time))
          fgat_rain_flags = .false.
          if ( INT(var4d_bin_rain/var4d_bin)*var4d_bin .ne. INT(var4d_bin_rain) ) then
             write(unit=message(1),fmt='(A,A,2I7)') &
             'Please change var4d_bin_rain in namelist and rerun==>', 'var4d_bin_rain, var4d_bin:',var4d_bin_rain,var4d_bin
             call da_error("da_solve.inc",147,message(1:1))
          endif
          do n = 1, num_fgat_time, INT(var4d_bin_rain/var4d_bin)
             fgat_rain_flags(n) = .true.
          end do
       end if
    endif

   !---------------------------------------------------------------------------
   ! [1.0] Initial checks
   !---------------------------------------------------------------------------

   if (cv_options_hum /= cv_options_hum_specific_humidity .and. &
       cv_options_hum /= cv_options_hum_relative_humidity) then
      write(unit=message(1),fmt='(A,I3)') &
         'Invalid cv_options_hum = ', cv_options_hum
      call da_error("da_solve.inc",163,message(1:1))
   end if

   if (vert_corr == vert_corr_2) then
      if (vertical_ip < vertical_ip_0 .or. vertical_ip > vertical_ip_delta_p) then
         write (unit=message(1),fmt='(A,I3)') &
           'Invalid vertical_ip = ', vertical_ip
         call da_error("da_solve.inc",170,message(1:1))
      end if
   end if

   if( use_rf )then
      if (0.5 * real(rf_passes) /= real(rf_passes / 2)) then
         write(unit=stdout,fmt='(A,I4,A)') &
            'rf_passes = ', rf_passes, ' .Should be even.'
         rf_passes = int(real(rf_passes / 2))
         write(unit=stdout,fmt='(A,I4)') 'Resetting rf_passes = ', rf_passes
      end if
   else
      write(stdout,'("da_solve: using wavelet transform")')
   endif

   if ( anal_type_hybrid_dual_res .and. alphacv_method .ne. alphacv_method_xa ) then
      write (unit=message(1),fmt='(A)') &
        'Dual-res hybrid only with alphacv_method = 2'
      call da_error("da_solve.inc",188,message(1:1))
   endif

   if (anal_type_randomcv) then
      ntmax = 0
      write(unit=stdout,fmt='(a)') &
         ' Resetting ntmax = 0 for analysis_type = randomcv' 
   end if


   !---------------------------------------------------------------------------
   ! [2.0] Initialise wrfvar parameters:
   !---------------------------------------------------------------------------

   if ( anal_type_hybrid_dual_res ) then

      !---------------------------------
      ! Get full ensemble grid dimensions
      !---------------------------------
      call da_solve_init(ensemble_grid &
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/actual_new_args.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
,grid%moist,grid%moist_bxs,grid%moist_bxe,grid%moist_bys,grid%moist_bye,grid%moist_btxs,grid%moist_btxe,grid%moist_btys, &
grid%moist_btye,grid%scalar,grid%scalar_bxs,grid%scalar_bxe,grid%scalar_bys,grid%scalar_bye,grid%scalar_btxs,grid%scalar_btxe, &
grid%scalar_btys,grid%scalar_btye,grid%a_moist,grid%g_moist,grid%a_scalar,grid%g_scalar,grid%chem,grid%tracer,grid%tracer_bxs, &
grid%tracer_bxe,grid%tracer_bys,grid%tracer_bye,grid%tracer_btxs,grid%tracer_btxe,grid%tracer_btys,grid%tracer_btye &
!ENDOFREGISTRYGENERATEDINCLUDE
)
      ide_ens = ide ! these are unstaggered dimensions of the full ensemble domain
      jde_ens = jde
      kde_ens = kde

      !---------------------------------------
      ! Get "intermediate" grid sizes and tiles
      !---------------------------------------
      call da_solve_init(grid%intermediate_grid &
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/actual_new_args.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
,grid%moist,grid%moist_bxs,grid%moist_bxe,grid%moist_bys,grid%moist_bye,grid%moist_btxs,grid%moist_btxe,grid%moist_btys, &
grid%moist_btye,grid%scalar,grid%scalar_bxs,grid%scalar_bxe,grid%scalar_bys,grid%scalar_bye,grid%scalar_btxs,grid%scalar_btxe, &
grid%scalar_btys,grid%scalar_btye,grid%a_moist,grid%g_moist,grid%a_scalar,grid%g_scalar,grid%chem,grid%tracer,grid%tracer_bxs, &
grid%tracer_bxe,grid%tracer_bys,grid%tracer_bye,grid%tracer_btxs,grid%tracer_btxe,grid%tracer_btys,grid%tracer_btye &
!ENDOFREGISTRYGENERATEDINCLUDE
)

      ! these are unstaggered dimensions of the "intermediate" ensemble domain
      !  The intermediate grid is the coarse (ensemble) domain that is co-located with the
      !  hi-resolution (analysis) grid

      ids_int = ids ; jds_int = jds ; kds_int = kds
      ide_int = ide ; jde_int = jde ; kde_int = kde
      
      its_int = its ; ite_int = ite
      jts_int = jts ; jte_int = jte
      kts_int = kts ; kte_int = kte

      ims_int = ims ; ime_int = ime
      jms_int = jms ; jme_int = jme
      kms_int = kms ; kme_int = kme

      ips_int = ips ; ipe_int = ipe
      jps_int = jps ; jpe_int = jpe
      kps_int = kps ; kpe_int = kpe


      grid%imask_xstag = 1   ; grid%imask_ystag = 1
      grid%imask_nostag = 1  ; grid%imask_xystag = 1

      !---------------------------------------------------------------------------
      ! De-allocate parts of grid and replace with grid%intermediate_grid dimensions
      !---------------------------------------------------------------------------
      call reallocate_analysis_grid(grid)

      !----------------------------------------------------------
      ! Allocate and initialize some of grid%intermediate_grid
      !----------------------------------------------------------
      call allocate_intermediate_grid(grid%intermediate_grid)


      !---------------------------------------
      ! Get map projection information for the ensemble
      !---------------------------------------

       call da_setup_firstguess(xbx, ensemble_grid, config_flags, .true. )
       map_info_ens = map_info ! map_info is read in from da_tools.f90, call it something else

   endif

   call da_solve_init(grid &
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/actual_new_args.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
,grid%moist,grid%moist_bxs,grid%moist_bxe,grid%moist_bys,grid%moist_bye,grid%moist_btxs,grid%moist_btxe,grid%moist_btys, &
grid%moist_btye,grid%scalar,grid%scalar_bxs,grid%scalar_bxe,grid%scalar_bys,grid%scalar_bye,grid%scalar_btxs,grid%scalar_btxe, &
grid%scalar_btys,grid%scalar_btye,grid%a_moist,grid%g_moist,grid%a_scalar,grid%g_scalar,grid%chem,grid%tracer,grid%tracer_bxs, &
grid%tracer_bxe,grid%tracer_bys,grid%tracer_bye,grid%tracer_btxs,grid%tracer_btxe,grid%tracer_btys,grid%tracer_btye &
!ENDOFREGISTRYGENERATEDINCLUDE
)

   if ( .not. anal_type_hybrid_dual_res ) then
      ide_ens = ide ; jde_ens = jde ; kde_ens = kde

      ids_int = ids ; ide_int = ide 
      jds_int = jds ; jde_int = jde 
      kds_int = kds ; kde_int = kde

      its_int = its ; ite_int = ite 
      jts_int = jts ; jte_int = jte
      kts_int = kts ; kte_int = kte

      ims_int = ims ; ime_int = ime
      jms_int = jms ; jme_int = jme
      kms_int = kms ; kme_int = kme

      ips_int = ips ; ipe_int = ipe
      jps_int = jps ; jpe_int = jpe
      kps_int = kps ; kpe_int = kpe
   endif

   !---------------------------------------------------------------------------
   ! [3.0] Set up first guess field (grid%xb):
   !---------------------------------------------------------------------------

   call da_setup_firstguess(xbx, grid, config_flags, .false.)

   if ( anal_type_hybrid_dual_res ) then

      ! 
      ! Get ensemble grid mapfactor on entire coarse grid
      ! 
      variable_name = 'MAPFAC_M'

      allocate( grid_box_area(1:ide_ens,1:jde_ens), mapfac(1:ide_ens,1:jde_ens) )  
      call da_get_var_2d_real_cdf( input_file_ens, variable_name, mapfac, ide_ens, jde_ens, 1, .false. )
      grid_box_area(:,:) = ( (ensemble_grid%dx)/mapfac(:,:) )**2 
      grid%intermediate_grid%xb%grid_box_area(its_int:ite_int,jts_int:jte_int) = grid_box_area(its_int:ite_int,jts_int:jte_int)
      deallocate(mapfac,grid_box_area)

      !
      ! Make a list of "observations" from the the fine grid lat/lon
      !

      xy = ( ime - ims + 1 ) * ( jme - jms + 1 )

      allocate(ob_locs(1:xy)) ! From da_control

      kk = 0 ; jjj = 0
      do xx = ims, ime !ids, ide 
         do yy = jms, jme !jds, jde 

          outside = .false.

          lat = grid%xb%lat(xx,yy)
          lon = grid%xb%lon(xx,yy)

          x_pos = -1.0 ; y_pos = -1.0
          call da_llxy_wrf(map_info_ens,lat,lon,x_pos,y_pos)
          call da_togrid(x_pos,its_int-2, ite_int+2, ii, dx, dxm )
          call da_togrid(y_pos,jts_int-2, jte_int+2, jj, dy, dym )

          if ((int(x_pos) < ids_int) .or. (int(x_pos) >= ide_int) .or. &
             (int(y_pos) < jds_int) .or. (int(y_pos) >= jde_int)) then
             outside = .true.
          endif

	  if ((ii < ids_int) .or. (ii >= ide_int) .or. &
	     (jj < jds_int) .or. (jj >= jde_int)) then
	     outside     = .true.
	  endif

         if ((ii < its_int-1) .or. (ii > ite_int) .or. &
            (jj < jts_int-1) .or. (jj > jte_int)) then
            outside = .true.
         endif

         if ( .not. outside ) then
              kk = kk + 1
              ob_locs(kk)%x = x_pos
              ob_locs(kk)%y = y_pos
              ob_locs(kk)%i = ii
              ob_locs(kk)%j = jj
              ob_locs(kk)%dx = dx
              ob_locs(kk)%dy = dy
              ob_locs(kk)%dxm = dxm
              ob_locs(kk)%dym = dym
              ob_locs(kk)%xx = xx
              ob_locs(kk)%yy = yy
          else
             jjj = jjj + 1
          endif

        enddo
     enddo

      total_here = kk ! from da_control

   endif 

   !---------------------------------------------------------------------------
   ! [4.0] Set up observations (ob):
   !---------------------------------------------------------------------------
   call da_setup_obs_structures (grid, ob, iv, j)
   if (use_rad) then
      allocate (j % jo % rad(1:iv%num_inst))
      do i=1,iv%num_inst
         allocate (j % jo % rad(i) % jo_ichan(iv%instid(i)%nchan))
         allocate (j % jo % rad(i) % num_ichan(iv%instid(i)%nchan))
      end do
   end if

   !---------------------------------------------------------------------------
   ! [4.1] Observer (ANAL_TYPE="VERIFY")
   !---------------------------------------------------------------------------

   if (anal_type_verify) then
      check_max_iv = .false.
      ntmax=0
      it = 1
      num_qcstat_conv=0

      if (use_rad .and. (use_varbc.or.freeze_varbc)) call da_varbc_init(iv, be)

      call da_get_innov_vector (it, num_qcstat_conv, ob, iv, grid , config_flags)
      call da_allocate_y (iv, re)
 
      ! write out O-B statistics
      call da_write_diagnostics(it, grid,num_qcstat_conv, ob, iv, re, y, j)

      ! write out Gradient of Jo for adjoint sensitivity
      if (adj_sens) then
         cv_size = 1
         allocate (xhat(cv_size))
         call da_allocate_y (iv, y)
         call da_allocate_y (iv, jo_grad_y)

         call da_calculate_residual(iv, y, re)
         call da_calculate_grady(iv, re, jo_grad_y)
         call da_zero_x(grid%xa)

         call da_transform_xtoy_adj(cv_size, xhat, grid, iv, jo_grad_y, grid%xa)
         call da_transform_xtoxa_adj(grid)
         call da_transfer_wrftltoxa_adj(grid, config_flags, 'fcst', timestr)

         call da_deallocate_y (y)
         call da_deallocate_y (jo_grad_y)
      end if

      call da_deallocate_y(re)
      call da_deallocate_observations (iv)
      if (trace_use) call da_trace_exit ("da_solve")
      return
   end if
   
   !---------------------------------------------------------------------------
   ! [5.0] Set up control variable:
   !---------------------------------------------------------------------------
   be % cv % size_jb = 0
   be % cv % size_je = 0
   be % cv % size_jp = 0
   be % cv % size_js = 0
   be % cv % size_jl = 0
   
   !---------------------------------------------------------------------------
   ! [5.1] Set up background errors (be):
   !---------------------------------------------------------------------------
   if (use_background_errors .and. multi_inc /= 1) then
      call da_setup_background_errors (grid, be)
   else
      be % ne = ensdim_alpha
      be % v1 % mz = 0
      be % v2 % mz = 0
      be % v3 % mz = 0
      be % v4 % mz = 0
      be % v5 % mz = 0
   end if

    ! overwrite variables defined in da_setup_cv.inc set in the call to da_setup_background_errors
    if ( anal_type_hybrid_dual_res ) then
       be % cv % size_alphac = (ite_int - its_int + 1) * (jte_int - jts_int + 1) * be % alpha % mz * be % ne
       be % cv % size_je = be % cv % size_alphac
       cv_size_domain_je = (ide_int - ids_int + 1) * (jde_int - jds_int + 1) * be % alpha % mz * be % ne
    endif
   
   !---------------------------------------------------------------------------
   ! [5.2] Set up observation bias correction (VarBC):
   !---------------------------------------------------------------------------
   if (use_rad .and. (use_varbc.or.freeze_varbc)) call da_varbc_init(iv, be)

   !---------------------------------------------------------------------------
   ! [5.3] Set up satellite control variable:
   !---------------------------------------------------------------------------
   if (ANY(use_satcv)) call da_setup_satcv(iv, be)
   
   !---------------------------------------------------------------------------
   ! [5.4] Total control variable:
   !---------------------------------------------------------------------------   
   be % cv % size = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_js + be%cv%size_jl
   cv_size = be % cv % size

   !---------------------------------------------------------------------------
   ! [6.0] Set up ensemble perturbation input:
   !---------------------------------------------------------------------------

      grid % ep % ne = be % ne
      if (use_background_errors .and. be % ne > 0) then
!        call da_setup_flow_predictors ( ide, jde, kde, be % ne, grid%ep, &
!                                        its, ite, jts, jte, kts, kte )
         call da_setup_flow_predictors ( ide_ens, jde_ens, kde_ens, be % ne, grid%ep, &
                                         its_int, ite_int, jts_int, jte_int, kts_int, kte_int )
      end if

   !---------------------------------------------------------------------------
   ! [7.0] Setup control variable (cv):
   !---------------------------------------------------------------------------

!  Dynamically allocate the variables which don't rely on ntmax
      allocate (cvt(1:cv_size))
      allocate (xhat(1:cv_size))
!      if (use_lanczos) then
!        allocate (full_eignvec(cv_size))
!      end if

      if ( outer_loop_restart ) then
         !call da_get_unit(cvt_unit)
         cvt_unit=600
         if ( max_ext_its > 1 ) then
           max_ext_its=1
           write(unit=message(1),fmt='(a)') "Re-set max_ext_its = 1 for outer_loop_restart"
           call da_message(message(1:1))
         end if
         write(unit=cvtfile,fmt='(a,i4.4)') 'cvt_',myproc
         inquire(file=trim(cvtfile), exist=ex)
         if ( ex ) then
           open(unit=cvt_unit,file=trim(cvtfile),iostat=iost,form='UNFORMATTED',status='OLD')
           if (iost /= 0) then
             write(unit=message(1),fmt='(A,I5,A)') &
                "Error ",iost," opening cvt file "//trim(cvtfile)
             call da_error("da_solve.inc",527,message(1:1))
           end if
           write(unit=message(1),fmt='(a)') 'Reading cvt from : '//trim(cvtfile) 
           call da_message(message(1:1))
           read(cvt_unit) cvt
           close(cvt_unit)
         else
           write(unit=message(1),fmt='(a)') "cvt file '"//trim(cvtfile)//"' does not exists, initiallizing cvt."
           call da_message(message(1:1))
           call da_initialize_cv (cv_size, cvt)
         end if
      else
         call da_initialize_cv (cv_size, cvt)
      end if

      call da_zero_vp_type (grid%vv)
      call da_zero_vp_type (grid%vp)
  
      if ( var4d ) then
         call da_zero_vp_type (grid%vv6)
         call da_zero_vp_type (grid%vp6)
      end if

   !---------------------------------------------------------------------------
   ! [8] Outerloop
   !---------------------------------------------------------------------------

   j_grad_norm_target = 1.0
   do it = 1, max_ext_its

!  Dynamically allocate the variables which depend on ntmax
      if (use_lanczos) then
         allocate (qhat(1:cv_size, 0:ntmax(it)))
         allocate (eignvec(ntmax(it), ntmax(it)))
         allocate (eignval(ntmax(it)))
      end if

! Re-scale the variances and the scale-length for outer-loop > 1:
      if (it > 1 .and. (cv_options == 5 .or. cv_options == 7)) then
         print '(/10X,"===> Re-set BE SCALINGS for outer-loop=",i2)', it
         call  da_scale_background_errors ( be, it )
      else if (it > 1 .and. cv_options == 3) then
         print '(/10X,"===> Re-set CV3 BE SCALINGS for outer-loop=",i2)', it
         call  da_scale_background_errors_cv3 ( grid, be, it )
      endif

      call da_initialize_cv (cv_size, xhat)

      ! [8.1] Calculate nonlinear model trajectory 

!     if (var4d .and. multi_inc /= 2 ) then
      if (var4d) then
         write(unit=message(1),fmt='(A)')'Please re-compile the code with 4dvar option'
         call da_error("da_solve.inc",676,message(1:1))
      end if

      ! [8.2] Calculate innovation vector (O-B):

      num_qcstat_conv=0
      call da_get_innov_vector (it, num_qcstat_conv, ob, iv, grid , config_flags)
      if ( multi_inc == 1 ) then 
         if (trace_use) call da_trace_exit ("da_solve")
         return
      end if

      if (test_transforms .or. test_gradient) then

         if (test_gradient) then
            call da_allocate_y (iv, re)
            call da_allocate_y (iv, y)
            call da_check_gradient (grid, config_flags, cv_size, xhat, cvt, 1.0e-10, 8, &
                 xbx, be, iv, y, re, j)
            call da_deallocate_y (re)
            call da_deallocate_y (y)
         endif

         if (test_transforms) then
            call da_check (grid, config_flags, cv_size, xbx, be, grid%ep, iv, &
                        grid%vv, grid%vp, y)
         endif

         if (trace_use) call da_trace_exit("da_solve")
         return

      end if

      ! [8.4] Minimize cost function:

      call da_allocate_y (iv, re)
      call da_allocate_y (iv, y)

      if (use_lanczos) then
         if (read_lanczos) then
            call da_lanczos_io('r',cv_size,ntmax(it),neign,eignvec,eignval,qhat)

            call da_kmat_mul(grid,config_flags,it,cv_size,xbx, &
                             be,iv,xhat,qhat,cvt,re,y,j,eignvec,eignval,neign)
          ! Output Cost Function
            call da_calculate_j(it, 1, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                             be%cv%size_jl, xbx, be, iv, xhat, cvt, re, y, j, xhat, grid, config_flags )

         else
            call da_minimise_lz(grid, config_flags, it, cv_size, xbx,& 
                             be, iv, j_grad_norm_target, xhat, qhat, cvt, re, y, j, eignvec, eignval, neign )


         end if

         if (write_lanczos) call da_lanczos_io('w',cv_size,ntmax(it),neign,eignvec,eignval,qhat)

         if (adj_sens) call da_sensitivity(grid,config_flags,it,cv_size,xbx, &
                                           be,iv,xhat,qhat,cvt,y,eignvec,eignval,neign )

	 
      else

         call da_minimise_cg( grid, config_flags, it, be % cv % size, & 
              xbx, be, iv, j_grad_norm_target, xhat, cvt, re, y, j)
      end if
 
      ! Update outer-loop control variable
      cvt = cvt + xhat
      if ( outer_loop_restart ) then
        open(unit=cvt_unit,status='unknown',file=trim(cvtfile),iostat=iost,form='UNFORMATTED')
        if (iost /= 0) then
           write(unit=message(1),fmt='(A,I5,A)') &
             "Error ",iost," opening cvt file "//trim(cvtfile)
          call da_error("da_solve.inc",751,message(1:1))
        end if
        write(unit=message(1),fmt='(a)') 'Writing cvt to : '//trim(cvtfile) 
        call da_message(message(1:1))
        write(cvt_unit) cvt
        close(cvt_unit)
        !call da_free_unit(cvt_unit)
      end if
      !------------------------------------------------------------------------

      ! reset cv to random noise
      if (anal_type_randomcv) then
         call da_set_randomcv (cv_size, xhat)
      end if

      ! [8.5] Update latest analysis solution:
  
      if (.not. var4d) then
         call da_transform_vtox (grid,cv_size,xbx,be,grid%ep,xhat,grid%vv,grid%vp)
      else
         call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,xhat(1:be%cv%size_jb),grid%vv,grid%vp)
      endif
      call da_transform_xtoxa (grid)

      ! [8.6] Only when use_radarobs = .false. and calc_w_increment =.true.,
      !       the w_increment need to be diagnosed:

      if (calc_w_increment .and. .not. use_radarobs .and. .not. var4d) then
         call da_uvprho_to_w_lin (grid)

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_RADAR_XA_W.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_RADAR_XA_W_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
      end if

      ! [8.7] Write out diagnostics

      call da_write_diagnostics (it, grid, num_qcstat_conv, ob, iv, re, y, j)

      ! Write "clean" QCed observations if requested:
      if (anal_type_qcobs) then
         ! if (it == 1) then
          if (write_mod_filtered_obs) then
            call da_write_modified_filtered_obs (grid, ob, iv, &
               coarse_ix, coarse_jy, start_x, start_y)
          else
            call da_write_filtered_obs (it, grid, ob, iv, &
               coarse_ix, coarse_jy, start_x, start_y)
          end if     
         ! end if     
      end if

      ! [8.7.1] Write Ascii radar OMB and OMA file

      if (use_radarobs) then
         call da_write_oa_radar_ascii (ob,iv,re,it)
      end if

      ! [8.3] Interpolate x_g to low resolution grid

      ! [8.8] Write Ascii radiance OMB and OMA file

      if (use_rad .and. write_oa_rad_ascii) then
         call da_write_oa_rad_ascii (it,ob,iv,re)
      end if

      ! [8.9] Update VarBC parameters and write output file
      if ( use_rad .and. (use_varbc.or.freeze_varbc) ) &
                call da_varbc_update(it, cv_size, xhat, iv)


      !------------------------------------------------------------------------
      ! [8.10] Output WRFVAR analysis and analysis increments:
      !------------------------------------------------------------------------

      call da_transfer_xatoanalysis (it, xbx, grid, config_flags)

     if ( it < max_ext_its .and. print_detail_outerloop ) then
        write(outerloop,'(i2.2)') it
        call da_update_firstguess(grid,'wrfvar_output_'//outerloop)
     end if

     call da_deallocate_y (re)
     call da_deallocate_y (y)


   ! Deallocate arrays which depend on ntmax
     if (use_lanczos) then
        deallocate (qhat)
        deallocate (eignvec)
        deallocate (eignval)
     end if


   end do

   ! output wrfvar analysis

   if ((config_flags%real_data_init_type == 1) .or. &
       (config_flags%real_data_init_type == 3)) then
      call da_update_firstguess(input_grid)
      call med_shutdown_io (input_grid, config_flags)
   end if

   !---------------------------------------------------------------------------
   ! [9.0] Tidy up:
   !---------------------------------------------------------------------------

   deallocate (cvt)
   deallocate (xhat)
!   if (use_lanczos) then
!      deallocate (full_eignvec)
!   end if

   ! clean up radiance related arrays
   if (use_rad) then
      call da_deallocate_radiance (ob, iv, j)
      deallocate (time_slots)
   end if

   if (var4d .and. use_rainobs) deallocate(fgat_rain_flags)
   call da_deallocate_observations (iv)
   call da_deallocate_y (ob)
   if (use_background_errors) call da_deallocate_background_errors (be)

   if (xbx%pad_num > 0) then
      deallocate (xbx%pad_loc)
      deallocate (xbx%pad_pos)
   end if

   deallocate (xbx % fft_factors_x)
   deallocate (xbx % fft_factors_y)
   deallocate (xbx % fft_coeffs)
   deallocate (xbx % trig_functs_x)
   deallocate (xbx % trig_functs_y)

   if (global) then
      deallocate (xbx%coslat)
      deallocate (xbx%sinlat)
      deallocate (xbx%coslon)
      deallocate (xbx%sinlon)
      deallocate (xbx%int_wgts)
      deallocate (xbx%alp)
      deallocate (xbx%wsave)
      if (jts == jds) then
         deallocate (cos_xls)
         deallocate (sin_xls)
      end if
                                                                                
      if (jte == jde) then
         deallocate (cos_xle)
         deallocate (sin_xle)
      end if
   end if

   if ( anal_type_hybrid_dual_res ) deallocate(ob_locs)



   call mpi_barrier (comm,ierr)

   if (trace_use) call da_trace_exit ("da_solve")


contains

subroutine da_solve_init(grid &
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/dummy_new_args.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,scalar,scalar_bxs,scalar_bxe, &
scalar_bys,scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,a_moist,g_moist,a_scalar,g_scalar,chem,tracer,tracer_bxs, &
tracer_bxe,tracer_bys,tracer_bye,tracer_btxs,tracer_btxe,tracer_btys,tracer_btye &
!ENDOFREGISTRYGENERATEDINCLUDE
)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)      :: grid

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/dummy_new_decl.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_a_moist)           :: a_moist
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_g_moist)           :: g_moist
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_a_scalar)           :: a_scalar
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_g_scalar)           :: g_scalar
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_tracer)           :: tracer
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btye
!ENDOFREGISTRYGENERATEDINCLUDE

   integer :: ii

   integer :: sm31,sm32,sm33,sm31x,sm32x,sm33x,sm31y,sm32y,sm33y

   ! if (dwordsize != rwordsize)
   ! else
   !    define add_msg_xpose_real add_msg_xpose_doubleprecision
   ! end if

   if (trace_use) call da_trace_entry("da_solve_init")

   ! De-reference dimension information stored in the grid data structure.

   call da_copy_dims(grid)

   ! Compute these starting and stopping locations for each tile and number 
   ! of tiles.

   call set_tiles (grid , ids , ide , jds , jde , ips , ipe , jps , jpe)

   call da_copy_tile_dims(grid)

   sm31             = grid%sm31
   sm32             = grid%sm32
   sm33             = grid%sm33
   sm31x            = grid%sm31x
   sm32x            = grid%sm32x
   sm33x            = grid%sm33x
   sm31y            = grid%sm31y
   sm32y            = grid%sm32y
   sm33y            = grid%sm33y

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/data_calls.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use) call da_trace("da_solve_init", &
      Message="Setup halo region communication")

   ! Define halo region communication.
   !-----------------------------------------------------------------------
   !  Stencils for patch communications
   !                           * * * * *
   !         *        * * *    * * * * *
   !       * + *      * + *    * * + * *
   !         *        * * *    * * * * *
   !                           * * * * *
   !ij vp%v1            x
   !ij xb%cori          x
   !ij xb%rho           x
   !ij xa%u             x
   !ij xa%v             x
   !--------------------------------------------------------------
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_INIT.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_INIT_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_PSICHI_UV.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_PSICHI_UV_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_BAL_EQN_ADJ.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_BAL_EQN_ADJ_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_PSICHI_UV_ADJ.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_PSICHI_UV_ADJ_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_SFC_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_SFC_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_SSMI_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_SSMI_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_2D_WORK.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_2D_WORK_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_RADAR_XA_W.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_RADAR_XA_W_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   if (trace_use) call da_trace("da_solve_init", &
      Message="Copy domain and transpose descriptors")

   ! Copy domain and transpose descriptors.

   grid%xp%domdesc = grid%domdesc
   do ii = 1, max_comms
     grid%xp%comms(ii) = grid%comms(ii)
   end do


   ! Fill background scalars:

   grid%xb%ids = grid%xp%ids 
   grid%xb%ide = grid%xp%ide
   grid%xb%jds = grid%xp%jds 
   grid%xb%jde = grid%xp%jde
   grid%xb%kds = grid%xp%kds 
   grid%xb%kde = grid%xp%kde 

   grid%xb%ims = grid%xp%ims 
   grid%xb%ime = grid%xp%ime
   grid%xb%jms = grid%xp%jms 
   grid%xb%jme = grid%xp%jme
   grid%xb%kms = grid%xp%kms 
   grid%xb%kme = grid%xp%kme 

   grid%xb%its = grid%xp%its 
   grid%xb%ite = grid%xp%ite
   grid%xb%jts = grid%xp%jts 
   grid%xb%jte = grid%xp%jte 
   grid%xb%kts = grid%xp%kts
   grid%xb%kte = grid%xp%kte 

   if (trace_use) call da_trace_exit("da_solve_init")

end subroutine da_solve_init


subroutine reallocate_analysis_grid(grid)

   implicit none

   type(domain), intent(inout)      :: grid
   integer :: sm31, em31, sm32, em32, sm33, em33

   if (trace_use) call da_trace_entry("da_solve_dual_res_init")

   !
   ! First deallocate the arrays associated with alpha and ensemble perturbations
   !

   if (trace_use) call da_trace("da_solve_dual_res_init", &
      Message="Deallocating arrays")


   IF ( ASSOCIATED( grid%ep%v1 ) ) THEN
     DEALLOCATE(grid%ep%v1,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%ep%v1. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%ep%v2 ) ) THEN
     DEALLOCATE(grid%ep%v2,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%ep%v2. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%ep%v3 ) ) THEN
     DEALLOCATE(grid%ep%v3,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%ep%v3. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%ep%v4 ) ) THEN
     DEALLOCATE(grid%ep%v4,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%ep%v4. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%ep%v5 ) ) THEN
     DEALLOCATE(grid%ep%v5,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%ep%v5. '
    endif
   ENDIF
  IF ( ASSOCIATED( grid%vp%alpha ) ) THEN
     DEALLOCATE(grid%vp%alpha,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%vp%alpha. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%vv%alpha ) ) THEN
     DEALLOCATE(grid%vv%alpha,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%vv%alpha. '
    endif
   ENDIF


   !
   ! Now, reallocate the arrays with the intermediate grid dimensions
   !

   if (trace_use) call da_trace("da_solve_dual_res_init", &
      Message="Reallocating arrays")

   sm31 = grid%intermediate_grid%sm31
   em31 = grid%intermediate_grid%em31
   sm32 = grid%intermediate_grid%sm32
   em32 = grid%intermediate_grid%em32
   sm33 = grid%intermediate_grid%sm33
   em33 = grid%intermediate_grid%em33


   !
   ! allocate grid%vp%alpha
   !
   ALLOCATE(grid%vp%alpha(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%vp%alpha(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%vp%alpha=0.

   !
   ! allocate grid%vv%alpha
   !
   ALLOCATE(grid%vv%alpha(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%vv%alpha(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%vv%alpha=0.

   !
   ! allocate grid%ep%v1
   !
   ALLOCATE(grid%ep%v1(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%ep%v1(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%ep%v1=0.

   !
   ! allocate grid%ep%v2
   !
   ALLOCATE(grid%ep%v2(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%ep%v2(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%ep%v2=0.


   !
   ! allocate grid%ep%v3
   !
   ALLOCATE(grid%ep%v3(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%ep%v3(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%ep%v3=0.

   !
   ! allocate grid%ep%v4
   !
   ALLOCATE(grid%ep%v4(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%ep%v4(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%ep%v4=0.

   !
   ! allocate grid%ep%v5
   !
   ALLOCATE(grid%ep%v5(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha),STAT=ierr)
   if (ierr.ne.0) then
     print *,' Failed to allocate grid%ep%v5(sm31:em31,sm32:em32,sm33:em33,1:config_flags%ensdim_alpha). '
   endif
   grid%ep%v5=0.
 
   if (trace_use) call da_trace_exit("da_solve_dual_res_init")

end subroutine reallocate_analysis_grid

!!!!!!!!!!!!!!

subroutine allocate_intermediate_grid(grid)

   type(domain), intent(inout)      :: grid 

   integer :: sm31, em31, sm32, em32, sm33, em33
   integer :: sm31x, em31x, sm32x, em32x, sm33x, em33x
   integer :: sm31y, em31y, sm32y, em32y, sm33y, em33y

   !
   ! First deallocate the arrays
   !

   IF ( ASSOCIATED( grid%xp%vxy ) ) THEN
     DEALLOCATE(grid%xp%vxy,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%vxy. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xp%v1z ) ) THEN
     DEALLOCATE(grid%xp%v1z,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%v1z. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xp%v1x ) ) THEN
     DEALLOCATE(grid%xp%v1x,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%v1x. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xp%v1y ) ) THEN
     DEALLOCATE(grid%xp%v1y,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%v1y. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xp%v2z ) ) THEN
     DEALLOCATE(grid%xp%v2z,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%v2z. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xp%v2x ) ) THEN
     DEALLOCATE(grid%xp%v2x,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%v2x. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xp%v2y ) ) THEN
     DEALLOCATE(grid%xp%v2y,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xp%v2y. '
    endif
   ENDIF
   IF ( ASSOCIATED( grid%xb%grid_box_area ) ) THEN
     DEALLOCATE(grid%xb%grid_box_area,STAT=ierr)
    if (ierr.ne.0) then
       print *, ' Failed to deallocate grid%xb%grid_box_area. '
    endif
   ENDIF

   sm31 = grid%sm31
   em31 = grid%em31
   sm32 = grid%sm32
   em32 = grid%em32
   sm33 = grid%sm33
   em33 = grid%em33

   sm31x = grid%sm31x
   em31x = grid%em31x
   sm32x = grid%sm32x
   em32x = grid%em32x
   sm33x = grid%sm33x
   em33x = grid%em33x

   sm31y = grid%sm31y
   em31y = grid%em31y
   sm32y = grid%sm32y
   em32y = grid%em32y
   sm33y = grid%sm33y
   em33y = grid%em33y

!
! allocate grid%xp%vxy
!
  ALLOCATE(grid%xp%vxy(sm31:em31,sm32:em32),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%vxy'
  endif
  grid%xp%vxy=0.

!
! allocate grid%xp%v1z
!
  ALLOCATE(grid%xp%v1z(sm31:em31,sm32:em32,sm33:em33),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%v1z'
  endif
  grid%xp%v1z=0.

!
! allocate grid%xp%v1x
!
  ALLOCATE(grid%xp%v1x(sm31x:em31x,sm32x:em32x,sm33x:em33x),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%v1x'
  endif
  grid%xp%v1x=0.

!
! allocate grid%xp%v1y
!
 ALLOCATE(grid%xp%v1y(sm31y:em31y,sm32y:em32y,sm33y:em33y),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%v1y'
  endif
  grid%xp%v1y=0.

!
! allocate grid%xp%v2z
!
  ALLOCATE(grid%xp%v2z(sm31:em31,sm32:em32,sm33:em33),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%v2z'
  endif
  grid%xp%v2z=0.


!
! allocate grid%xp%v2x
!
  ALLOCATE(grid%xp%v2x(sm31x:em31x,sm32x:em32x,sm33x:em33x),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%v2x'
  endif
  grid%xp%v2x=0.


!
! allocate grid%xp%v2y
!
  ALLOCATE(grid%xp%v2y(sm31y:em31y,sm32y:em32y,sm33y:em33y),STAT=ierr)
  if (ierr.ne.0) then
      print *,' Failed to allocate grid%xp%v2y'
  endif
  grid%xp%v2y=0.

!
! allocate grid%xb%grid_box_area
!
   ALLOCATE(grid%xb%grid_box_area(sm31:em31,sm32:em32),STAT=ierr)
    if (ierr.ne.0) then
      print *,' Failed to allocate grid%xb%grid_box_area(sm31:em31,sm32:em32)'
    endif
    grid%xb%grid_box_area=0.
  
end subroutine allocate_intermediate_grid



end subroutine da_solve


end module da_wrfvar_top
