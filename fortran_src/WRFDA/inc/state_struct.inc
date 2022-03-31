!STARTOFREGISTRYGENERATEDINCLUDE 'inc/state_struct.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
real                                     :: cfn
real                                     :: cfn1
integer                                  :: step_number
logical                                  :: this_is_an_ideal_run
logical                                  :: stepping_to_time
integer                                  :: last_step_updated
real                                     :: rdx
real                                     :: rdy
real                                     :: dts
real                                     :: dtseps
real                                     :: resm
real                                     :: zetatop
real                                     :: cf1
real                                     :: cf2
real                                     :: cf3
integer                                  :: itimestep
real                                     :: xtime
real                                     :: julian
integer                                  :: lbc_fid
logical                                  :: tiled
logical                                  :: patched
logical                                  :: press_adj
real                                     :: dtbc
real                                     :: u_frame
real                                     :: v_frame
real                                     :: p_top
real                                     :: t00
real                                     :: p00
real                                     :: tlp
real                                     :: tiso
real                                     :: tlp_strat
real                                     :: p_strat
real                                     :: max_msftx
real                                     :: max_msfty
integer                                  :: stepra
logical                                  :: warm_rain
logical                                  :: adv_moist_cond
integer                                  :: save_topo_from_real
integer                                  :: run_days
integer                                  :: run_hours
integer                                  :: run_minutes
integer                                  :: run_seconds
integer                                  :: start_year
integer                                  :: start_month
integer                                  :: start_day
integer                                  :: start_hour
integer                                  :: start_minute
integer                                  :: start_second
integer                                  :: end_year
integer                                  :: end_month
integer                                  :: end_day
integer                                  :: end_hour
integer                                  :: end_minute
integer                                  :: end_second
integer                                  :: interval_seconds
logical                                  :: input_from_file
integer                                  :: fine_input_stream
logical                                  :: input_from_hires
character*256                               :: rsmas_data_path
logical                                  :: all_ic_times
integer                                  :: julyr
integer                                  :: julday
real                                     :: gmt
character*256                               :: input_inname
character*256                               :: input_outname
character*256                               :: bdy_inname
character*256                               :: bdy_outname
character*256                               :: rst_inname
character*256                               :: rst_outname
logical                                  :: write_input
logical                                  :: write_restart_at_0h
logical                                  :: write_hist_at_0h_rst
logical                                  :: adjust_output_times
logical                                  :: adjust_input_times
integer                                  :: diag_print
logical                                  :: nocolons
logical                                  :: cycling
integer                                  :: output_diagnostics
integer                                  :: nwp_diagnostics
logical                                  :: output_ready_flag
logical                                  :: usepio
integer                                  :: pioprocs
integer                                  :: piostart
integer                                  :: piostride
integer                                  :: pioshift
integer                                  :: dfi_opt
integer                                  :: dfi_savehydmeteors
integer                                  :: dfi_nfilter
logical                                  :: dfi_write_filtered_input
logical                                  :: dfi_write_dfi_history
integer                                  :: dfi_cutoff_seconds
integer                                  :: dfi_time_dim
integer                                  :: dfi_fwdstop_year
integer                                  :: dfi_fwdstop_month
integer                                  :: dfi_fwdstop_day
integer                                  :: dfi_fwdstop_hour
integer                                  :: dfi_fwdstop_minute
integer                                  :: dfi_fwdstop_second
integer                                  :: dfi_bckstop_year
integer                                  :: dfi_bckstop_month
integer                                  :: dfi_bckstop_day
integer                                  :: dfi_bckstop_hour
integer                                  :: dfi_bckstop_minute
integer                                  :: dfi_bckstop_second
integer                                  :: time_step
integer                                  :: time_step_fract_num
integer                                  :: time_step_fract_den
integer                                  :: time_step_dfi
integer                                  :: min_time_step
integer                                  :: min_time_step_den
integer                                  :: max_time_step
integer                                  :: max_time_step_den
real                                     :: target_cfl
real                                     :: target_hcfl
integer                                  :: max_step_increase_pct
integer                                  :: starting_time_step
integer                                  :: starting_time_step_den
logical                                  :: step_to_output_time
integer                                  :: adaptation_domain
logical                                  :: use_adaptive_time_step
logical                                  :: use_adaptive_time_step_dfi
integer                                  :: max_dom
integer                                  :: lats_to_mic
integer                                  :: s_we
integer                                  :: e_we
integer                                  :: s_sn
integer                                  :: e_sn
integer                                  :: s_vert
integer                                  :: e_vert
integer                                  :: num_metgrid_levels
integer                                  :: num_metgrid_soil_levels
real                                     :: p_top_requested
logical                                  :: interp_theta
integer                                  :: interp_type
integer                                  :: rebalance
integer                                  :: vert_refine_method
integer                                  :: vert_refine_fact
integer                                  :: extrap_type
integer                                  :: t_extrap_type
integer                                  :: hypsometric_opt
logical                                  :: lowest_lev_from_sfc
logical                                  :: use_levels_below_ground
logical                                  :: use_tavg_for_tsk
logical                                  :: use_surface
integer                                  :: lagrange_order
integer                                  :: force_sfc_in_vinterp
real                                     :: zap_close_levels
real                                     :: maxw_horiz_pres_diff
real                                     :: trop_horiz_pres_diff
real                                     :: maxw_above_this_level
integer                                  :: use_maxw_level
integer                                  :: use_trop_level
logical                                  :: sfcp_to_sfcp
logical                                  :: adjust_heights
logical                                  :: smooth_cg_topo
integer                                  :: nest_interp_coord
integer                                  :: interp_method_type
logical                                  :: aggregate_lu
logical                                  :: rh2qv_wrt_liquid
integer                                  :: rh2qv_method
real                                     :: qv_max_p_safe
real                                     :: qv_max_flag
real                                     :: qv_max_value
real                                     :: qv_min_p_safe
real                                     :: qv_min_flag
real                                     :: qv_min_value
integer                                  :: ideal_init_method
real                                     :: dx
real                                     :: dy
integer                                  :: grid_id
logical                                  :: grid_allowed
integer                                  :: parent_id
integer                                  :: i_parent_start
integer                                  :: j_parent_start
integer                                  :: parent_grid_ratio
integer                                  :: parent_time_step_ratio
integer                                  :: feedback
integer                                  :: smooth_option
integer                                  :: blend_width
real                                     :: ztop
integer                                  :: moad_grid_ratio
integer                                  :: moad_time_step_ratio
integer                                  :: shw
integer                                  :: tile_sz_x
integer                                  :: tile_sz_y
integer                                  :: numtiles
integer                                  :: numtiles_inc
integer                                  :: numtiles_x
integer                                  :: numtiles_y
integer                                  :: tile_strategy
integer                                  :: nproc_x
integer                                  :: nproc_y
integer                                  :: irand
real                                     :: dt
integer                                  :: fft_used
integer                                  :: cu_used
integer                                  :: shcu_used
integer                                  :: cam_used
integer                                  :: alloc_qndropsource
integer                                  :: num_moves
integer                                  :: ts_buf_size
integer                                  :: max_ts_locs
integer                                  :: vortex_interval
integer                                  :: max_vortex_speed
integer                                  :: corral_dist
integer                                  :: track_level
real                                     :: time_to_move
integer                                  :: move_id
integer                                  :: move_interval
integer                                  :: move_cd_x
integer                                  :: move_cd_y
logical                                  :: swap_x
logical                                  :: swap_y
logical                                  :: cycle_x
logical                                  :: cycle_y
logical                                  :: reorder_mesh
logical                                  :: perturb_input
real                                     :: eta_levels
real                                     :: max_dz
integer                                  :: ocean_levels
real                                     :: ocean_z
real                                     :: ocean_t
real                                     :: ocean_s
integer                                  :: num_traj
integer                                  :: max_ts_level
integer                                  :: track_loc_in
integer                                  :: num_ext_model_couple_dom
logical                                  :: insert_bogus_storm
logical                                  :: remove_storm
integer                                  :: num_storm
real                                     :: latc_loc
real                                     :: lonc_loc
real                                     :: vmax_meters_per_second
real                                     :: rmax
real                                     :: vmax_ratio
real                                     :: rankine_lid
logical                                  :: force_read_thompson
logical                                  :: write_thompson_tables
integer                                  :: mp_physics
real                                     :: nssl_cccn
real                                     :: nssl_alphah
real                                     :: nssl_alphahl
real                                     :: nssl_cnoh
real                                     :: nssl_cnohl
real                                     :: nssl_cnor
real                                     :: nssl_cnos
real                                     :: nssl_rho_qh
real                                     :: nssl_rho_qhl
real                                     :: nssl_rho_qs
integer                                  :: nudge_lightning
integer                                  :: nudge_light_times
integer                                  :: nudge_light_timee
integer                                  :: nudge_light_int
character*256                               :: path_to_files
integer                                  :: gsfcgce_hail
integer                                  :: gsfcgce_2ice
integer                                  :: progn
real                                     :: accum_mode
real                                     :: aitken_mode
real                                     :: coarse_mode
integer                                  :: do_radar_ref
integer                                  :: compute_radar_ref
integer                                  :: ra_lw_physics
integer                                  :: ra_sw_physics
real                                     :: radt
real                                     :: naer
integer                                  :: sf_sfclay_physics
integer                                  :: sf_surface_physics
integer                                  :: bl_pbl_physics
integer                                  :: bl_mynn_tkebudget
integer                                  :: ysu_topdown_pblmix
integer                                  :: shinhong_tke_diag
logical                                  :: bl_mynn_tkeadvect
integer                                  :: bl_mynn_cloudpdf
integer                                  :: bl_mynn_mixlength
integer                                  :: bl_mynn_edmf
integer                                  :: bl_mynn_edmf_mom
integer                                  :: bl_mynn_edmf_tke
integer                                  :: bl_mynn_edmf_part
integer                                  :: bl_mynn_cloudmix
integer                                  :: bl_mynn_mixqt
integer                                  :: icloud_bl
integer                                  :: mfshconv
integer                                  :: sf_urban_physics
real                                     :: bldt
integer                                  :: cu_physics
integer                                  :: shcu_physics
integer                                  :: cu_diag
integer                                  :: kf_edrates
integer                                  :: kfeta_trigger
integer                                  :: nsas_dx_factor
real                                     :: cudt
real                                     :: gsmdt
integer                                  :: isfflx
integer                                  :: ifsnow
integer                                  :: icloud
integer                                  :: ideal_xland
real                                     :: swrad_scat
integer                                  :: surface_input_source
integer                                  :: num_soil_layers
integer                                  :: maxpatch
integer                                  :: num_snow_layers
integer                                  :: num_snso_layers
integer                                  :: num_urban_layers
integer                                  :: num_urban_hi
integer                                  :: num_months
integer                                  :: sf_surface_mosaic
integer                                  :: mosaic_cat
integer                                  :: mosaic_cat_soil
integer                                  :: mosaic_lu
integer                                  :: mosaic_soil
integer                                  :: maxiens
integer                                  :: maxens
integer                                  :: maxens2
integer                                  :: maxens3
integer                                  :: ensdim
integer                                  :: cugd_avedx
integer                                  :: clos_choice
integer                                  :: imomentum
integer                                  :: ishallow
real                                     :: convtrans_avglen_m
integer                                  :: num_land_cat
integer                                  :: num_soil_cat
integer                                  :: mp_zero_out
real                                     :: mp_zero_out_thresh
real                                     :: seaice_threshold
integer                                  :: sst_update
integer                                  :: sst_skin
integer                                  :: tmn_update
logical                                  :: usemonalb
logical                                  :: rdmaxalb
logical                                  :: rdlai2d
logical                                  :: ua_phys
integer                                  :: opt_thcnd
integer                                  :: co2tf
integer                                  :: ra_call_offset
real                                     :: cam_abs_freq_s
integer                                  :: levsiz
integer                                  :: paerlev
integer                                  :: cam_abs_dim1
integer                                  :: cam_abs_dim2
integer                                  :: lagday
integer                                  :: no_src_types
integer                                  :: alevsiz
integer                                  :: o3input
integer                                  :: aer_opt
integer                                  :: swint_opt
integer                                  :: aer_type
integer                                  :: aer_aod550_opt
integer                                  :: aer_angexp_opt
integer                                  :: aer_ssa_opt
integer                                  :: aer_asy_opt
real                                     :: aer_aod550_val
real                                     :: aer_angexp_val
real                                     :: aer_ssa_val
real                                     :: aer_asy_val
logical                                  :: cu_rad_feedback
logical                                  :: shallowcu_forced_ra
integer                                  :: numbins
real                                     :: thbinsize
real                                     :: rbinsize
real                                     :: mindeepfreq
real                                     :: minshallowfreq
integer                                  :: shcu_aerosols_opt
integer                                  :: icloud_cu
integer                                  :: pxlsm_smois_init
integer                                  :: omlcall
integer                                  :: sf_ocean_physics
integer                                  :: traj_opt
integer                                  :: tracercall
real                                     :: omdt
real                                     :: oml_hml0
real                                     :: oml_gamma
real                                     :: oml_relaxation_time
integer                                  :: isftcflx
integer                                  :: iz0tlnd
real                                     :: shadlen
integer                                  :: slope_rad
integer                                  :: topo_shading
integer                                  :: topo_wind
integer                                  :: no_mp_heating
integer                                  :: fractional_seaice
integer                                  :: seaice_snowdepth_opt
real                                     :: seaice_snowdepth_max
real                                     :: seaice_snowdepth_min
integer                                  :: seaice_albedo_opt
real                                     :: seaice_albedo_default
integer                                  :: seaice_thickness_opt
real                                     :: seaice_thickness_default
logical                                  :: tice2tsk_if2cold
real                                     :: bucket_mm
real                                     :: bucket_j
real                                     :: mp_tend_lim
real                                     :: prec_acc_dt
integer                                  :: prec_acc_opt
integer                                  :: bucketr_opt
integer                                  :: process_time_series
integer                                  :: grav_settling
real                                     :: sas_pgcon
integer                                  :: scalar_pblmix
integer                                  :: tracer_pblmix
logical                                  :: use_aero_icbc
logical                                  :: use_rap_aero_icbc
integer                                  :: use_mp_re
real                                     :: ccn_conc
integer                                  :: hail_opt
integer                                  :: dveg
integer                                  :: opt_crs
integer                                  :: opt_btr
integer                                  :: opt_run
integer                                  :: opt_sfc
integer                                  :: opt_frz
integer                                  :: opt_inf
integer                                  :: opt_rad
integer                                  :: opt_alb
integer                                  :: opt_snf
integer                                  :: opt_tbot
integer                                  :: opt_stc
integer                                  :: opt_gla
integer                                  :: opt_rsf
real                                     :: wtddt
integer                                  :: wrf_hydro
real                                     :: fgdt
integer                                  :: fgdtzero
integer                                  :: grid_fdda
integer                                  :: grid_sfdda
integer                                  :: if_no_pbl_nudging_uv
integer                                  :: if_no_pbl_nudging_t
integer                                  :: if_no_pbl_nudging_ph
integer                                  :: if_no_pbl_nudging_q
integer                                  :: if_zfac_uv
integer                                  :: k_zfac_uv
integer                                  :: if_zfac_t
integer                                  :: k_zfac_t
integer                                  :: if_zfac_ph
integer                                  :: k_zfac_ph
integer                                  :: if_zfac_q
integer                                  :: k_zfac_q
integer                                  :: dk_zfac_uv
integer                                  :: dk_zfac_t
integer                                  :: dk_zfac_ph
real                                     :: guv
real                                     :: guv_sfc
real                                     :: gt
real                                     :: gt_sfc
real                                     :: gq
real                                     :: gq_sfc
real                                     :: gph
real                                     :: dtramp_min
integer                                  :: if_ramping
real                                     :: rinblw
integer                                  :: xwavenum
integer                                  :: ywavenum
integer                                  :: pxlsm_soil_nudge
integer                                  :: fasdas
integer                                  :: obs_nudge_opt
integer                                  :: max_obs
real                                     :: fdda_start
real                                     :: fdda_end
integer                                  :: obs_nudge_wind
real                                     :: obs_coef_wind
integer                                  :: obs_nudge_temp
real                                     :: obs_coef_temp
integer                                  :: obs_nudge_mois
real                                     :: obs_coef_mois
integer                                  :: obs_nudge_pstr
real                                     :: obs_coef_pstr
integer                                  :: obs_no_pbl_nudge_uv
integer                                  :: obs_no_pbl_nudge_t
integer                                  :: obs_no_pbl_nudge_q
integer                                  :: obs_sfc_scheme_horiz
integer                                  :: obs_sfc_scheme_vert
real                                     :: obs_max_sndng_gap
real                                     :: obs_nudgezfullr1_uv
real                                     :: obs_nudgezrampr1_uv
real                                     :: obs_nudgezfullr2_uv
real                                     :: obs_nudgezrampr2_uv
real                                     :: obs_nudgezfullr4_uv
real                                     :: obs_nudgezrampr4_uv
real                                     :: obs_nudgezfullr1_t
real                                     :: obs_nudgezrampr1_t
real                                     :: obs_nudgezfullr2_t
real                                     :: obs_nudgezrampr2_t
real                                     :: obs_nudgezfullr4_t
real                                     :: obs_nudgezrampr4_t
real                                     :: obs_nudgezfullr1_q
real                                     :: obs_nudgezrampr1_q
real                                     :: obs_nudgezfullr2_q
real                                     :: obs_nudgezrampr2_q
real                                     :: obs_nudgezfullr4_q
real                                     :: obs_nudgezrampr4_q
real                                     :: obs_nudgezfullmin
real                                     :: obs_nudgezrampmin
real                                     :: obs_nudgezmax
real                                     :: obs_sfcfact
real                                     :: obs_sfcfacr
real                                     :: obs_dpsmx
real                                     :: obs_rinxy
real                                     :: obs_rinsig
real                                     :: obs_twindo
integer                                  :: obs_npfi
integer                                  :: obs_ionf
integer                                  :: obs_idynin
real                                     :: obs_dtramp
integer                                  :: obs_prt_max
integer                                  :: obs_prt_freq
logical                                  :: obs_ipf_in4dob
logical                                  :: obs_ipf_errob
logical                                  :: obs_ipf_nudob
logical                                  :: obs_ipf_init
integer                                  :: obs_scl_neg_qv_innov
integer                                  :: scm_force
real                                     :: scm_force_dx
integer                                  :: num_force_layers
integer                                  :: scm_lu_index
integer                                  :: scm_isltyp
real                                     :: scm_vegfra
real                                     :: scm_canwat
real                                     :: scm_lat
real                                     :: scm_lon
logical                                  :: scm_th_t_tend
logical                                  :: scm_qv_t_tend
logical                                  :: scm_th_adv
logical                                  :: scm_wind_adv
logical                                  :: scm_qv_adv
logical                                  :: scm_ql_adv
logical                                  :: scm_vert_adv
integer                                  :: num_force_soil_layers
logical                                  :: scm_soilt_force
logical                                  :: scm_soilq_force
logical                                  :: scm_force_th_largescale
logical                                  :: scm_force_qv_largescale
logical                                  :: scm_force_ql_largescale
logical                                  :: scm_force_wind_largescale
integer                                  :: scm_force_skintemp
integer                                  :: scm_force_flux
integer                                  :: dyn_opt
integer                                  :: rk_ord
integer                                  :: w_damping
integer                                  :: diff_opt
integer                                  :: diff_opt_dfi
integer                                  :: km_opt
integer                                  :: km_opt_dfi
integer                                  :: damp_opt
integer                                  :: rad_nudge
integer                                  :: gwd_opt
real                                     :: zdamp
real                                     :: dampcoef
real                                     :: khdif
real                                     :: kvdif
real                                     :: diff_6th_factor
integer                                  :: diff_6th_opt
integer                                  :: use_theta_m
integer                                  :: use_q_diabatic
real                                     :: c_s
real                                     :: c_k
real                                     :: smdiv
real                                     :: emdiv
real                                     :: epssm
logical                                  :: non_hydrostatic
logical                                  :: use_input_w
integer                                  :: time_step_sound
integer                                  :: h_mom_adv_order
integer                                  :: v_mom_adv_order
integer                                  :: h_sca_adv_order
integer                                  :: v_sca_adv_order
integer                                  :: momentum_adv_opt
integer                                  :: moist_adv_opt
integer                                  :: moist_adv_dfi_opt
integer                                  :: chem_adv_opt
integer                                  :: tracer_adv_opt
integer                                  :: scalar_adv_opt
integer                                  :: tke_adv_opt
logical                                  :: top_radiation
integer                                  :: mix_isotropic
real                                     :: mix_upper_bound
logical                                  :: top_lid
real                                     :: tke_upper_bound
real                                     :: tke_drag_coefficient
real                                     :: tke_heat_flux
logical                                  :: pert_coriolis
logical                                  :: coriolis2d
logical                                  :: mix_full_fields
real                                     :: base_pres
real                                     :: base_temp
real                                     :: base_lapse
real                                     :: iso_temp
real                                     :: base_pres_strat
real                                     :: base_lapse_strat
logical                                  :: use_baseparam_fr_nml
real                                     :: fft_filter_lat
logical                                  :: coupled_filtering
logical                                  :: pos_def
logical                                  :: swap_pole_with_next_j
logical                                  :: actual_distance_average
logical                                  :: rotated_pole
logical                                  :: do_coriolis
logical                                  :: do_curvature
logical                                  :: do_gradp
integer                                  :: tracer_opt
integer                                  :: tenddiag
integer                                  :: spec_bdy_width
integer                                  :: spec_zone
integer                                  :: relax_zone
logical                                  :: specified
logical                                  :: constant_bc
logical                                  :: periodic_x
logical                                  :: symmetric_xs
logical                                  :: symmetric_xe
logical                                  :: open_xs
logical                                  :: open_xe
logical                                  :: periodic_y
logical                                  :: symmetric_ys
logical                                  :: symmetric_ye
logical                                  :: open_ys
logical                                  :: open_ye
logical                                  :: polar
logical                                  :: nested
real                                     :: spec_exp
integer                                  :: spec_bdy_final_mu
integer                                  :: real_data_init_type
logical                                  :: have_bcs_moist
logical                                  :: have_bcs_scalar
integer                                  :: background_proc_id
integer                                  :: forecast_proc_id
integer                                  :: production_status
integer                                  :: compression
integer                                  :: nobs_ndg_vars
integer                                  :: nobs_err_flds
real                                     :: cen_lat
real                                     :: cen_lon
real                                     :: truelat1
real                                     :: truelat2
real                                     :: moad_cen_lat
real                                     :: stand_lon
real                                     :: pole_lat
real                                     :: pole_lon
integer                                  :: flag_metgrid
integer                                  :: flag_snow
integer                                  :: flag_psfc
integer                                  :: flag_sm000010
integer                                  :: flag_sm010040
integer                                  :: flag_sm040100
integer                                  :: flag_sm100200
integer                                  :: flag_st000010
integer                                  :: flag_st010040
integer                                  :: flag_st040100
integer                                  :: flag_st100200
integer                                  :: flag_soil_layers
integer                                  :: flag_slp
integer                                  :: flag_soilhgt
integer                                  :: flag_mf_xy
integer                                  :: flag_um_soil
real                                     :: bdyfrq
character*256                               :: mminlu
integer                                  :: iswater
integer                                  :: islake
integer                                  :: isice
integer                                  :: isurban
integer                                  :: isoilwater
integer                                  :: map_proj
integer                                  :: use_wps_input
integer                                  :: dfi_stage
integer                                  :: mp_physics_dfi
integer                                  :: bl_pbl_physics_dfi
integer                                  :: windfarm_opt
integer                                  :: windfarm_ij
integer                                  :: lightning_option
real                                     :: lightning_dt
real                                     :: lightning_start_seconds
real                                     :: flashrate_factor
integer                                  :: iccg_method
real                                     :: iccg_prescribed_num
real                                     :: iccg_prescribed_den
integer                                  :: cellcount_method
real                                     :: cldtop_adjustment
integer                                  :: sf_lake_physics
character*256                               :: auxinput1_inname
integer                                  :: io_form_auxinput1
logical                                  :: override_restart_timers
integer                                  :: auxhist1_oid
character*256                               :: auxhist1_inname
character*256                               :: auxhist1_outname
integer                                  :: auxhist1_interval_y
integer                                  :: auxhist1_interval_d
integer                                  :: auxhist1_interval_h
integer                                  :: auxhist1_interval_m
integer                                  :: auxhist1_interval_s
integer                                  :: auxhist1_interval
integer                                  :: auxhist1_begin_y
integer                                  :: auxhist1_begin_d
integer                                  :: auxhist1_begin_h
integer                                  :: auxhist1_begin_m
integer                                  :: auxhist1_begin_s
integer                                  :: auxhist1_begin
integer                                  :: auxhist1_end_y
integer                                  :: auxhist1_end_d
integer                                  :: auxhist1_end_h
integer                                  :: auxhist1_end_m
integer                                  :: auxhist1_end_s
integer                                  :: auxhist1_end
integer                                  :: io_form_auxhist1
integer                                  :: frames_per_auxhist1
integer                                  :: auxhist2_oid
character*256                               :: auxhist2_inname
character*256                               :: auxhist2_outname
integer                                  :: auxhist2_interval_y
integer                                  :: auxhist2_interval_d
integer                                  :: auxhist2_interval_h
integer                                  :: auxhist2_interval_m
integer                                  :: auxhist2_interval_s
integer                                  :: auxhist2_interval
integer                                  :: auxhist2_begin_y
integer                                  :: auxhist2_begin_d
integer                                  :: auxhist2_begin_h
integer                                  :: auxhist2_begin_m
integer                                  :: auxhist2_begin_s
integer                                  :: auxhist2_begin
integer                                  :: auxhist2_end_y
integer                                  :: auxhist2_end_d
integer                                  :: auxhist2_end_h
integer                                  :: auxhist2_end_m
integer                                  :: auxhist2_end_s
integer                                  :: auxhist2_end
integer                                  :: io_form_auxhist2
integer                                  :: frames_per_auxhist2
integer                                  :: auxhist3_oid
character*256                               :: auxhist3_inname
character*256                               :: auxhist3_outname
integer                                  :: auxhist3_interval_y
integer                                  :: auxhist3_interval_d
integer                                  :: auxhist3_interval_h
integer                                  :: auxhist3_interval_m
integer                                  :: auxhist3_interval_s
integer                                  :: auxhist3_interval
integer                                  :: auxhist3_begin_y
integer                                  :: auxhist3_begin_d
integer                                  :: auxhist3_begin_h
integer                                  :: auxhist3_begin_m
integer                                  :: auxhist3_begin_s
integer                                  :: auxhist3_begin
integer                                  :: auxhist3_end_y
integer                                  :: auxhist3_end_d
integer                                  :: auxhist3_end_h
integer                                  :: auxhist3_end_m
integer                                  :: auxhist3_end_s
integer                                  :: auxhist3_end
integer                                  :: io_form_auxhist3
integer                                  :: frames_per_auxhist3
integer                                  :: auxhist4_oid
character*256                               :: auxhist4_inname
character*256                               :: auxhist4_outname
integer                                  :: auxhist4_interval_y
integer                                  :: auxhist4_interval_d
integer                                  :: auxhist4_interval_h
integer                                  :: auxhist4_interval_m
integer                                  :: auxhist4_interval_s
integer                                  :: auxhist4_interval
integer                                  :: auxhist4_begin_y
integer                                  :: auxhist4_begin_d
integer                                  :: auxhist4_begin_h
integer                                  :: auxhist4_begin_m
integer                                  :: auxhist4_begin_s
integer                                  :: auxhist4_begin
integer                                  :: auxhist4_end_y
integer                                  :: auxhist4_end_d
integer                                  :: auxhist4_end_h
integer                                  :: auxhist4_end_m
integer                                  :: auxhist4_end_s
integer                                  :: auxhist4_end
integer                                  :: io_form_auxhist4
integer                                  :: frames_per_auxhist4
integer                                  :: auxhist5_oid
character*256                               :: auxhist5_inname
character*256                               :: auxhist5_outname
integer                                  :: auxhist5_interval_y
integer                                  :: auxhist5_interval_d
integer                                  :: auxhist5_interval_h
integer                                  :: auxhist5_interval_m
integer                                  :: auxhist5_interval_s
integer                                  :: auxhist5_interval
integer                                  :: auxhist5_begin_y
integer                                  :: auxhist5_begin_d
integer                                  :: auxhist5_begin_h
integer                                  :: auxhist5_begin_m
integer                                  :: auxhist5_begin_s
integer                                  :: auxhist5_begin
integer                                  :: auxhist5_end_y
integer                                  :: auxhist5_end_d
integer                                  :: auxhist5_end_h
integer                                  :: auxhist5_end_m
integer                                  :: auxhist5_end_s
integer                                  :: auxhist5_end
integer                                  :: io_form_auxhist5
integer                                  :: frames_per_auxhist5
integer                                  :: auxhist6_oid
character*256                               :: auxhist6_inname
character*256                               :: auxhist6_outname
integer                                  :: auxhist6_interval_y
integer                                  :: auxhist6_interval_d
integer                                  :: auxhist6_interval_h
integer                                  :: auxhist6_interval_m
integer                                  :: auxhist6_interval_s
integer                                  :: auxhist6_interval
integer                                  :: auxhist6_begin_y
integer                                  :: auxhist6_begin_d
integer                                  :: auxhist6_begin_h
integer                                  :: auxhist6_begin_m
integer                                  :: auxhist6_begin_s
integer                                  :: auxhist6_begin
integer                                  :: auxhist6_end_y
integer                                  :: auxhist6_end_d
integer                                  :: auxhist6_end_h
integer                                  :: auxhist6_end_m
integer                                  :: auxhist6_end_s
integer                                  :: auxhist6_end
integer                                  :: io_form_auxhist6
integer                                  :: frames_per_auxhist6
integer                                  :: auxhist7_oid
character*256                               :: auxhist7_inname
character*256                               :: auxhist7_outname
integer                                  :: auxhist7_interval_y
integer                                  :: auxhist7_interval_d
integer                                  :: auxhist7_interval_h
integer                                  :: auxhist7_interval_m
integer                                  :: auxhist7_interval_s
integer                                  :: auxhist7_interval
integer                                  :: auxhist7_begin_y
integer                                  :: auxhist7_begin_d
integer                                  :: auxhist7_begin_h
integer                                  :: auxhist7_begin_m
integer                                  :: auxhist7_begin_s
integer                                  :: auxhist7_begin
integer                                  :: auxhist7_end_y
integer                                  :: auxhist7_end_d
integer                                  :: auxhist7_end_h
integer                                  :: auxhist7_end_m
integer                                  :: auxhist7_end_s
integer                                  :: auxhist7_end
integer                                  :: io_form_auxhist7
integer                                  :: frames_per_auxhist7
integer                                  :: auxhist8_oid
character*256                               :: auxhist8_inname
character*256                               :: auxhist8_outname
integer                                  :: auxhist8_interval_y
integer                                  :: auxhist8_interval_d
integer                                  :: auxhist8_interval_h
integer                                  :: auxhist8_interval_m
integer                                  :: auxhist8_interval_s
integer                                  :: auxhist8_interval
integer                                  :: auxhist8_begin_y
integer                                  :: auxhist8_begin_d
integer                                  :: auxhist8_begin_h
integer                                  :: auxhist8_begin_m
integer                                  :: auxhist8_begin_s
integer                                  :: auxhist8_begin
integer                                  :: auxhist8_end_y
integer                                  :: auxhist8_end_d
integer                                  :: auxhist8_end_h
integer                                  :: auxhist8_end_m
integer                                  :: auxhist8_end_s
integer                                  :: auxhist8_end
integer                                  :: io_form_auxhist8
integer                                  :: frames_per_auxhist8
integer                                  :: auxhist9_oid
character*256                               :: auxhist9_inname
character*256                               :: auxhist9_outname
integer                                  :: auxhist9_interval_y
integer                                  :: auxhist9_interval_d
integer                                  :: auxhist9_interval_h
integer                                  :: auxhist9_interval_m
integer                                  :: auxhist9_interval_s
integer                                  :: auxhist9_interval
integer                                  :: auxhist9_begin_y
integer                                  :: auxhist9_begin_d
integer                                  :: auxhist9_begin_h
integer                                  :: auxhist9_begin_m
integer                                  :: auxhist9_begin_s
integer                                  :: auxhist9_begin
integer                                  :: auxhist9_end_y
integer                                  :: auxhist9_end_d
integer                                  :: auxhist9_end_h
integer                                  :: auxhist9_end_m
integer                                  :: auxhist9_end_s
integer                                  :: auxhist9_end
integer                                  :: io_form_auxhist9
integer                                  :: frames_per_auxhist9
integer                                  :: auxhist10_oid
character*256                               :: auxhist10_inname
character*256                               :: auxhist10_outname
integer                                  :: auxhist10_interval_y
integer                                  :: auxhist10_interval_d
integer                                  :: auxhist10_interval_h
integer                                  :: auxhist10_interval_m
integer                                  :: auxhist10_interval_s
integer                                  :: auxhist10_interval
integer                                  :: auxhist10_begin_y
integer                                  :: auxhist10_begin_d
integer                                  :: auxhist10_begin_h
integer                                  :: auxhist10_begin_m
integer                                  :: auxhist10_begin_s
integer                                  :: auxhist10_begin
integer                                  :: auxhist10_end_y
integer                                  :: auxhist10_end_d
integer                                  :: auxhist10_end_h
integer                                  :: auxhist10_end_m
integer                                  :: auxhist10_end_s
integer                                  :: auxhist10_end
integer                                  :: io_form_auxhist10
integer                                  :: frames_per_auxhist10
integer                                  :: auxhist11_oid
character*256                               :: auxhist11_inname
character*256                               :: auxhist11_outname
integer                                  :: auxhist11_interval_y
integer                                  :: auxhist11_interval_d
integer                                  :: auxhist11_interval_h
integer                                  :: auxhist11_interval_m
integer                                  :: auxhist11_interval_s
integer                                  :: auxhist11_interval
integer                                  :: auxhist11_begin_y
integer                                  :: auxhist11_begin_d
integer                                  :: auxhist11_begin_h
integer                                  :: auxhist11_begin_m
integer                                  :: auxhist11_begin_s
integer                                  :: auxhist11_begin
integer                                  :: auxhist11_end_y
integer                                  :: auxhist11_end_d
integer                                  :: auxhist11_end_h
integer                                  :: auxhist11_end_m
integer                                  :: auxhist11_end_s
integer                                  :: auxhist11_end
integer                                  :: io_form_auxhist11
integer                                  :: frames_per_auxhist11
integer                                  :: auxhist12_oid
character*256                               :: auxhist12_inname
character*256                               :: auxhist12_outname
integer                                  :: auxhist12_interval_y
integer                                  :: auxhist12_interval_d
integer                                  :: auxhist12_interval_h
integer                                  :: auxhist12_interval_m
integer                                  :: auxhist12_interval_s
integer                                  :: auxhist12_interval
integer                                  :: auxhist12_begin_y
integer                                  :: auxhist12_begin_d
integer                                  :: auxhist12_begin_h
integer                                  :: auxhist12_begin_m
integer                                  :: auxhist12_begin_s
integer                                  :: auxhist12_begin
integer                                  :: auxhist12_end_y
integer                                  :: auxhist12_end_d
integer                                  :: auxhist12_end_h
integer                                  :: auxhist12_end_m
integer                                  :: auxhist12_end_s
integer                                  :: auxhist12_end
integer                                  :: io_form_auxhist12
integer                                  :: frames_per_auxhist12
integer                                  :: auxhist13_oid
character*256                               :: auxhist13_inname
character*256                               :: auxhist13_outname
integer                                  :: auxhist13_interval_y
integer                                  :: auxhist13_interval_d
integer                                  :: auxhist13_interval_h
integer                                  :: auxhist13_interval_m
integer                                  :: auxhist13_interval_s
integer                                  :: auxhist13_interval
integer                                  :: auxhist13_begin_y
integer                                  :: auxhist13_begin_d
integer                                  :: auxhist13_begin_h
integer                                  :: auxhist13_begin_m
integer                                  :: auxhist13_begin_s
integer                                  :: auxhist13_begin
integer                                  :: auxhist13_end_y
integer                                  :: auxhist13_end_d
integer                                  :: auxhist13_end_h
integer                                  :: auxhist13_end_m
integer                                  :: auxhist13_end_s
integer                                  :: auxhist13_end
integer                                  :: io_form_auxhist13
integer                                  :: frames_per_auxhist13
integer                                  :: auxhist14_oid
character*256                               :: auxhist14_inname
character*256                               :: auxhist14_outname
integer                                  :: auxhist14_interval_y
integer                                  :: auxhist14_interval_d
integer                                  :: auxhist14_interval_h
integer                                  :: auxhist14_interval_m
integer                                  :: auxhist14_interval_s
integer                                  :: auxhist14_interval
integer                                  :: auxhist14_begin_y
integer                                  :: auxhist14_begin_d
integer                                  :: auxhist14_begin_h
integer                                  :: auxhist14_begin_m
integer                                  :: auxhist14_begin_s
integer                                  :: auxhist14_begin
integer                                  :: auxhist14_end_y
integer                                  :: auxhist14_end_d
integer                                  :: auxhist14_end_h
integer                                  :: auxhist14_end_m
integer                                  :: auxhist14_end_s
integer                                  :: auxhist14_end
integer                                  :: io_form_auxhist14
integer                                  :: frames_per_auxhist14
integer                                  :: auxhist15_oid
character*256                               :: auxhist15_inname
character*256                               :: auxhist15_outname
integer                                  :: auxhist15_interval_y
integer                                  :: auxhist15_interval_d
integer                                  :: auxhist15_interval_h
integer                                  :: auxhist15_interval_m
integer                                  :: auxhist15_interval_s
integer                                  :: auxhist15_interval
integer                                  :: auxhist15_begin_y
integer                                  :: auxhist15_begin_d
integer                                  :: auxhist15_begin_h
integer                                  :: auxhist15_begin_m
integer                                  :: auxhist15_begin_s
integer                                  :: auxhist15_begin
integer                                  :: auxhist15_end_y
integer                                  :: auxhist15_end_d
integer                                  :: auxhist15_end_h
integer                                  :: auxhist15_end_m
integer                                  :: auxhist15_end_s
integer                                  :: auxhist15_end
integer                                  :: io_form_auxhist15
integer                                  :: frames_per_auxhist15
integer                                  :: auxhist16_oid
character*256                               :: auxhist16_inname
character*256                               :: auxhist16_outname
integer                                  :: auxhist16_interval_y
integer                                  :: auxhist16_interval_d
integer                                  :: auxhist16_interval_h
integer                                  :: auxhist16_interval_m
integer                                  :: auxhist16_interval_s
integer                                  :: auxhist16_interval
integer                                  :: auxhist16_begin_y
integer                                  :: auxhist16_begin_d
integer                                  :: auxhist16_begin_h
integer                                  :: auxhist16_begin_m
integer                                  :: auxhist16_begin_s
integer                                  :: auxhist16_begin
integer                                  :: auxhist16_end_y
integer                                  :: auxhist16_end_d
integer                                  :: auxhist16_end_h
integer                                  :: auxhist16_end_m
integer                                  :: auxhist16_end_s
integer                                  :: auxhist16_end
integer                                  :: io_form_auxhist16
integer                                  :: frames_per_auxhist16
integer                                  :: auxhist17_oid
character*256                               :: auxhist17_inname
character*256                               :: auxhist17_outname
integer                                  :: auxhist17_interval_y
integer                                  :: auxhist17_interval_d
integer                                  :: auxhist17_interval_h
integer                                  :: auxhist17_interval_m
integer                                  :: auxhist17_interval_s
integer                                  :: auxhist17_interval
integer                                  :: auxhist17_begin_y
integer                                  :: auxhist17_begin_d
integer                                  :: auxhist17_begin_h
integer                                  :: auxhist17_begin_m
integer                                  :: auxhist17_begin_s
integer                                  :: auxhist17_begin
integer                                  :: auxhist17_end_y
integer                                  :: auxhist17_end_d
integer                                  :: auxhist17_end_h
integer                                  :: auxhist17_end_m
integer                                  :: auxhist17_end_s
integer                                  :: auxhist17_end
integer                                  :: io_form_auxhist17
integer                                  :: frames_per_auxhist17
integer                                  :: auxhist18_oid
character*256                               :: auxhist18_inname
character*256                               :: auxhist18_outname
integer                                  :: auxhist18_interval_y
integer                                  :: auxhist18_interval_d
integer                                  :: auxhist18_interval_h
integer                                  :: auxhist18_interval_m
integer                                  :: auxhist18_interval_s
integer                                  :: auxhist18_interval
integer                                  :: auxhist18_begin_y
integer                                  :: auxhist18_begin_d
integer                                  :: auxhist18_begin_h
integer                                  :: auxhist18_begin_m
integer                                  :: auxhist18_begin_s
integer                                  :: auxhist18_begin
integer                                  :: auxhist18_end_y
integer                                  :: auxhist18_end_d
integer                                  :: auxhist18_end_h
integer                                  :: auxhist18_end_m
integer                                  :: auxhist18_end_s
integer                                  :: auxhist18_end
integer                                  :: io_form_auxhist18
integer                                  :: frames_per_auxhist18
integer                                  :: auxhist19_oid
character*256                               :: auxhist19_inname
character*256                               :: auxhist19_outname
integer                                  :: auxhist19_interval_y
integer                                  :: auxhist19_interval_d
integer                                  :: auxhist19_interval_h
integer                                  :: auxhist19_interval_m
integer                                  :: auxhist19_interval_s
integer                                  :: auxhist19_interval
integer                                  :: auxhist19_begin_y
integer                                  :: auxhist19_begin_d
integer                                  :: auxhist19_begin_h
integer                                  :: auxhist19_begin_m
integer                                  :: auxhist19_begin_s
integer                                  :: auxhist19_begin
integer                                  :: auxhist19_end_y
integer                                  :: auxhist19_end_d
integer                                  :: auxhist19_end_h
integer                                  :: auxhist19_end_m
integer                                  :: auxhist19_end_s
integer                                  :: auxhist19_end
integer                                  :: io_form_auxhist19
integer                                  :: frames_per_auxhist19
integer                                  :: auxhist20_oid
character*256                               :: auxhist20_inname
character*256                               :: auxhist20_outname
integer                                  :: auxhist20_interval_y
integer                                  :: auxhist20_interval_d
integer                                  :: auxhist20_interval_h
integer                                  :: auxhist20_interval_m
integer                                  :: auxhist20_interval_s
integer                                  :: auxhist20_interval
integer                                  :: auxhist20_begin_y
integer                                  :: auxhist20_begin_d
integer                                  :: auxhist20_begin_h
integer                                  :: auxhist20_begin_m
integer                                  :: auxhist20_begin_s
integer                                  :: auxhist20_begin
integer                                  :: auxhist20_end_y
integer                                  :: auxhist20_end_d
integer                                  :: auxhist20_end_h
integer                                  :: auxhist20_end_m
integer                                  :: auxhist20_end_s
integer                                  :: auxhist20_end
integer                                  :: io_form_auxhist20
integer                                  :: frames_per_auxhist20
integer                                  :: auxhist21_oid
character*256                               :: auxhist21_inname
character*256                               :: auxhist21_outname
integer                                  :: auxhist21_interval_y
integer                                  :: auxhist21_interval_d
integer                                  :: auxhist21_interval_h
integer                                  :: auxhist21_interval_m
integer                                  :: auxhist21_interval_s
integer                                  :: auxhist21_interval
integer                                  :: auxhist21_begin_y
integer                                  :: auxhist21_begin_d
integer                                  :: auxhist21_begin_h
integer                                  :: auxhist21_begin_m
integer                                  :: auxhist21_begin_s
integer                                  :: auxhist21_begin
integer                                  :: auxhist21_end_y
integer                                  :: auxhist21_end_d
integer                                  :: auxhist21_end_h
integer                                  :: auxhist21_end_m
integer                                  :: auxhist21_end_s
integer                                  :: auxhist21_end
integer                                  :: io_form_auxhist21
integer                                  :: frames_per_auxhist21
integer                                  :: auxhist22_oid
character*256                               :: auxhist22_inname
character*256                               :: auxhist22_outname
integer                                  :: auxhist22_interval_y
integer                                  :: auxhist22_interval_d
integer                                  :: auxhist22_interval_h
integer                                  :: auxhist22_interval_m
integer                                  :: auxhist22_interval_s
integer                                  :: auxhist22_interval
integer                                  :: auxhist22_begin_y
integer                                  :: auxhist22_begin_d
integer                                  :: auxhist22_begin_h
integer                                  :: auxhist22_begin_m
integer                                  :: auxhist22_begin_s
integer                                  :: auxhist22_begin
integer                                  :: auxhist22_end_y
integer                                  :: auxhist22_end_d
integer                                  :: auxhist22_end_h
integer                                  :: auxhist22_end_m
integer                                  :: auxhist22_end_s
integer                                  :: auxhist22_end
integer                                  :: io_form_auxhist22
integer                                  :: frames_per_auxhist22
integer                                  :: auxhist23_oid
character*256                               :: auxhist23_inname
character*256                               :: auxhist23_outname
integer                                  :: auxhist23_interval_y
integer                                  :: auxhist23_interval_d
integer                                  :: auxhist23_interval_h
integer                                  :: auxhist23_interval_m
integer                                  :: auxhist23_interval_s
integer                                  :: auxhist23_interval
integer                                  :: auxhist23_begin_y
integer                                  :: auxhist23_begin_d
integer                                  :: auxhist23_begin_h
integer                                  :: auxhist23_begin_m
integer                                  :: auxhist23_begin_s
integer                                  :: auxhist23_begin
integer                                  :: auxhist23_end_y
integer                                  :: auxhist23_end_d
integer                                  :: auxhist23_end_h
integer                                  :: auxhist23_end_m
integer                                  :: auxhist23_end_s
integer                                  :: auxhist23_end
integer                                  :: io_form_auxhist23
integer                                  :: frames_per_auxhist23
integer                                  :: auxhist24_oid
character*256                               :: auxhist24_inname
character*256                               :: auxhist24_outname
integer                                  :: auxhist24_interval_y
integer                                  :: auxhist24_interval_d
integer                                  :: auxhist24_interval_h
integer                                  :: auxhist24_interval_m
integer                                  :: auxhist24_interval_s
integer                                  :: auxhist24_interval
integer                                  :: auxhist24_begin_y
integer                                  :: auxhist24_begin_d
integer                                  :: auxhist24_begin_h
integer                                  :: auxhist24_begin_m
integer                                  :: auxhist24_begin_s
integer                                  :: auxhist24_begin
integer                                  :: auxhist24_end_y
integer                                  :: auxhist24_end_d
integer                                  :: auxhist24_end_h
integer                                  :: auxhist24_end_m
integer                                  :: auxhist24_end_s
integer                                  :: auxhist24_end
integer                                  :: io_form_auxhist24
integer                                  :: frames_per_auxhist24
integer                                  :: auxinput1_oid
character*256                               :: auxinput1_outname
integer                                  :: auxinput1_interval_y
integer                                  :: auxinput1_interval_d
integer                                  :: auxinput1_interval_h
integer                                  :: auxinput1_interval_m
integer                                  :: auxinput1_interval_s
integer                                  :: auxinput1_interval
integer                                  :: auxinput1_begin_y
integer                                  :: auxinput1_begin_d
integer                                  :: auxinput1_begin_h
integer                                  :: auxinput1_begin_m
integer                                  :: auxinput1_begin_s
integer                                  :: auxinput1_begin
integer                                  :: auxinput1_end_y
integer                                  :: auxinput1_end_d
integer                                  :: auxinput1_end_h
integer                                  :: auxinput1_end_m
integer                                  :: auxinput1_end_s
integer                                  :: auxinput1_end
integer                                  :: frames_per_auxinput1
integer                                  :: auxinput2_oid
character*256                               :: auxinput2_inname
character*256                               :: auxinput2_outname
integer                                  :: auxinput2_interval_y
integer                                  :: auxinput2_interval_d
integer                                  :: auxinput2_interval_h
integer                                  :: auxinput2_interval_m
integer                                  :: auxinput2_interval_s
integer                                  :: auxinput2_interval
integer                                  :: auxinput2_begin_y
integer                                  :: auxinput2_begin_d
integer                                  :: auxinput2_begin_h
integer                                  :: auxinput2_begin_m
integer                                  :: auxinput2_begin_s
integer                                  :: auxinput2_begin
integer                                  :: auxinput2_end_y
integer                                  :: auxinput2_end_d
integer                                  :: auxinput2_end_h
integer                                  :: auxinput2_end_m
integer                                  :: auxinput2_end_s
integer                                  :: auxinput2_end
integer                                  :: io_form_auxinput2
integer                                  :: frames_per_auxinput2
integer                                  :: auxinput3_oid
character*256                               :: auxinput3_inname
character*256                               :: auxinput3_outname
integer                                  :: auxinput3_interval_y
integer                                  :: auxinput3_interval_d
integer                                  :: auxinput3_interval_h
integer                                  :: auxinput3_interval_m
integer                                  :: auxinput3_interval_s
integer                                  :: auxinput3_interval
integer                                  :: auxinput3_begin_y
integer                                  :: auxinput3_begin_d
integer                                  :: auxinput3_begin_h
integer                                  :: auxinput3_begin_m
integer                                  :: auxinput3_begin_s
integer                                  :: auxinput3_begin
integer                                  :: auxinput3_end_y
integer                                  :: auxinput3_end_d
integer                                  :: auxinput3_end_h
integer                                  :: auxinput3_end_m
integer                                  :: auxinput3_end_s
integer                                  :: auxinput3_end
integer                                  :: io_form_auxinput3
integer                                  :: frames_per_auxinput3
integer                                  :: auxinput4_oid
character*256                               :: auxinput4_inname
character*256                               :: auxinput4_outname
integer                                  :: auxinput4_interval_y
integer                                  :: auxinput4_interval_d
integer                                  :: auxinput4_interval_h
integer                                  :: auxinput4_interval_m
integer                                  :: auxinput4_interval_s
integer                                  :: auxinput4_interval
integer                                  :: auxinput4_begin_y
integer                                  :: auxinput4_begin_d
integer                                  :: auxinput4_begin_h
integer                                  :: auxinput4_begin_m
integer                                  :: auxinput4_begin_s
integer                                  :: auxinput4_begin
integer                                  :: auxinput4_end_y
integer                                  :: auxinput4_end_d
integer                                  :: auxinput4_end_h
integer                                  :: auxinput4_end_m
integer                                  :: auxinput4_end_s
integer                                  :: auxinput4_end
integer                                  :: io_form_auxinput4
integer                                  :: frames_per_auxinput4
integer                                  :: auxinput5_oid
character*256                               :: auxinput5_inname
character*256                               :: auxinput5_outname
integer                                  :: auxinput5_interval_y
integer                                  :: auxinput5_interval_d
integer                                  :: auxinput5_interval_h
integer                                  :: auxinput5_interval_m
integer                                  :: auxinput5_interval_s
integer                                  :: auxinput5_interval
integer                                  :: auxinput5_begin_y
integer                                  :: auxinput5_begin_d
integer                                  :: auxinput5_begin_h
integer                                  :: auxinput5_begin_m
integer                                  :: auxinput5_begin_s
integer                                  :: auxinput5_begin
integer                                  :: auxinput5_end_y
integer                                  :: auxinput5_end_d
integer                                  :: auxinput5_end_h
integer                                  :: auxinput5_end_m
integer                                  :: auxinput5_end_s
integer                                  :: auxinput5_end
integer                                  :: io_form_auxinput5
integer                                  :: frames_per_auxinput5
integer                                  :: auxinput6_oid
character*256                               :: auxinput6_inname
character*256                               :: auxinput6_outname
integer                                  :: auxinput6_interval_y
integer                                  :: auxinput6_interval_d
integer                                  :: auxinput6_interval_h
integer                                  :: auxinput6_interval_m
integer                                  :: auxinput6_interval_s
integer                                  :: auxinput6_interval
integer                                  :: auxinput6_begin_y
integer                                  :: auxinput6_begin_d
integer                                  :: auxinput6_begin_h
integer                                  :: auxinput6_begin_m
integer                                  :: auxinput6_begin_s
integer                                  :: auxinput6_begin
integer                                  :: auxinput6_end_y
integer                                  :: auxinput6_end_d
integer                                  :: auxinput6_end_h
integer                                  :: auxinput6_end_m
integer                                  :: auxinput6_end_s
integer                                  :: auxinput6_end
integer                                  :: io_form_auxinput6
integer                                  :: frames_per_auxinput6
integer                                  :: auxinput7_oid
character*256                               :: auxinput7_inname
character*256                               :: auxinput7_outname
integer                                  :: auxinput7_interval_y
integer                                  :: auxinput7_interval_d
integer                                  :: auxinput7_interval_h
integer                                  :: auxinput7_interval_m
integer                                  :: auxinput7_interval_s
integer                                  :: auxinput7_interval
integer                                  :: auxinput7_begin_y
integer                                  :: auxinput7_begin_d
integer                                  :: auxinput7_begin_h
integer                                  :: auxinput7_begin_m
integer                                  :: auxinput7_begin_s
integer                                  :: auxinput7_begin
integer                                  :: auxinput7_end_y
integer                                  :: auxinput7_end_d
integer                                  :: auxinput7_end_h
integer                                  :: auxinput7_end_m
integer                                  :: auxinput7_end_s
integer                                  :: auxinput7_end
integer                                  :: io_form_auxinput7
integer                                  :: frames_per_auxinput7
integer                                  :: auxinput8_oid
character*256                               :: auxinput8_inname
character*256                               :: auxinput8_outname
integer                                  :: auxinput8_interval_y
integer                                  :: auxinput8_interval_d
integer                                  :: auxinput8_interval_h
integer                                  :: auxinput8_interval_m
integer                                  :: auxinput8_interval_s
integer                                  :: auxinput8_interval
integer                                  :: auxinput8_begin_y
integer                                  :: auxinput8_begin_d
integer                                  :: auxinput8_begin_h
integer                                  :: auxinput8_begin_m
integer                                  :: auxinput8_begin_s
integer                                  :: auxinput8_begin
integer                                  :: auxinput8_end_y
integer                                  :: auxinput8_end_d
integer                                  :: auxinput8_end_h
integer                                  :: auxinput8_end_m
integer                                  :: auxinput8_end_s
integer                                  :: auxinput8_end
integer                                  :: io_form_auxinput8
integer                                  :: frames_per_auxinput8
integer                                  :: auxinput9_oid
character*256                               :: auxinput9_inname
character*256                               :: auxinput9_outname
integer                                  :: auxinput9_interval_y
integer                                  :: auxinput9_interval_d
integer                                  :: auxinput9_interval_h
integer                                  :: auxinput9_interval_m
integer                                  :: auxinput9_interval_s
integer                                  :: auxinput9_interval
integer                                  :: auxinput9_begin_y
integer                                  :: auxinput9_begin_d
integer                                  :: auxinput9_begin_h
integer                                  :: auxinput9_begin_m
integer                                  :: auxinput9_begin_s
integer                                  :: auxinput9_begin
integer                                  :: auxinput9_end_y
integer                                  :: auxinput9_end_d
integer                                  :: auxinput9_end_h
integer                                  :: auxinput9_end_m
integer                                  :: auxinput9_end_s
integer                                  :: auxinput9_end
integer                                  :: io_form_auxinput9
integer                                  :: frames_per_auxinput9
integer                                  :: auxinput10_oid
character*256                               :: auxinput10_inname
character*256                               :: auxinput10_outname
integer                                  :: auxinput10_interval_y
integer                                  :: auxinput10_interval_d
integer                                  :: auxinput10_interval_h
integer                                  :: auxinput10_interval_m
integer                                  :: auxinput10_interval_s
integer                                  :: auxinput10_interval
integer                                  :: auxinput10_begin_y
integer                                  :: auxinput10_begin_d
integer                                  :: auxinput10_begin_h
integer                                  :: auxinput10_begin_m
integer                                  :: auxinput10_begin_s
integer                                  :: auxinput10_begin
integer                                  :: auxinput10_end_y
integer                                  :: auxinput10_end_d
integer                                  :: auxinput10_end_h
integer                                  :: auxinput10_end_m
integer                                  :: auxinput10_end_s
integer                                  :: auxinput10_end
integer                                  :: io_form_auxinput10
integer                                  :: frames_per_auxinput10
integer                                  :: auxinput11_oid
character*256                               :: auxinput11_inname
character*256                               :: auxinput11_outname
integer                                  :: auxinput11_interval_y
integer                                  :: auxinput11_interval_d
integer                                  :: auxinput11_interval_h
integer                                  :: auxinput11_interval_m
integer                                  :: auxinput11_interval_s
integer                                  :: auxinput11_interval
integer                                  :: auxinput11_begin_y
integer                                  :: auxinput11_begin_d
integer                                  :: auxinput11_begin_h
integer                                  :: auxinput11_begin_m
integer                                  :: auxinput11_begin_s
integer                                  :: auxinput11_begin
integer                                  :: auxinput11_end_y
integer                                  :: auxinput11_end_d
integer                                  :: auxinput11_end_h
integer                                  :: auxinput11_end_m
integer                                  :: auxinput11_end_s
integer                                  :: auxinput11_end
integer                                  :: io_form_auxinput11
integer                                  :: frames_per_auxinput11
integer                                  :: auxinput12_oid
character*256                               :: auxinput12_inname
character*256                               :: auxinput12_outname
integer                                  :: auxinput12_interval_y
integer                                  :: auxinput12_interval_d
integer                                  :: auxinput12_interval_h
integer                                  :: auxinput12_interval_m
integer                                  :: auxinput12_interval_s
integer                                  :: auxinput12_interval
integer                                  :: auxinput12_begin_y
integer                                  :: auxinput12_begin_d
integer                                  :: auxinput12_begin_h
integer                                  :: auxinput12_begin_m
integer                                  :: auxinput12_begin_s
integer                                  :: auxinput12_begin
integer                                  :: auxinput12_end_y
integer                                  :: auxinput12_end_d
integer                                  :: auxinput12_end_h
integer                                  :: auxinput12_end_m
integer                                  :: auxinput12_end_s
integer                                  :: auxinput12_end
integer                                  :: io_form_auxinput12
integer                                  :: frames_per_auxinput12
integer                                  :: auxinput13_oid
character*256                               :: auxinput13_inname
character*256                               :: auxinput13_outname
integer                                  :: auxinput13_interval_y
integer                                  :: auxinput13_interval_d
integer                                  :: auxinput13_interval_h
integer                                  :: auxinput13_interval_m
integer                                  :: auxinput13_interval_s
integer                                  :: auxinput13_interval
integer                                  :: auxinput13_begin_y
integer                                  :: auxinput13_begin_d
integer                                  :: auxinput13_begin_h
integer                                  :: auxinput13_begin_m
integer                                  :: auxinput13_begin_s
integer                                  :: auxinput13_begin
integer                                  :: auxinput13_end_y
integer                                  :: auxinput13_end_d
integer                                  :: auxinput13_end_h
integer                                  :: auxinput13_end_m
integer                                  :: auxinput13_end_s
integer                                  :: auxinput13_end
integer                                  :: io_form_auxinput13
integer                                  :: frames_per_auxinput13
integer                                  :: auxinput14_oid
character*256                               :: auxinput14_inname
character*256                               :: auxinput14_outname
integer                                  :: auxinput14_interval_y
integer                                  :: auxinput14_interval_d
integer                                  :: auxinput14_interval_h
integer                                  :: auxinput14_interval_m
integer                                  :: auxinput14_interval_s
integer                                  :: auxinput14_interval
integer                                  :: auxinput14_begin_y
integer                                  :: auxinput14_begin_d
integer                                  :: auxinput14_begin_h
integer                                  :: auxinput14_begin_m
integer                                  :: auxinput14_begin_s
integer                                  :: auxinput14_begin
integer                                  :: auxinput14_end_y
integer                                  :: auxinput14_end_d
integer                                  :: auxinput14_end_h
integer                                  :: auxinput14_end_m
integer                                  :: auxinput14_end_s
integer                                  :: auxinput14_end
integer                                  :: io_form_auxinput14
integer                                  :: frames_per_auxinput14
integer                                  :: auxinput15_oid
character*256                               :: auxinput15_inname
character*256                               :: auxinput15_outname
integer                                  :: auxinput15_interval_y
integer                                  :: auxinput15_interval_d
integer                                  :: auxinput15_interval_h
integer                                  :: auxinput15_interval_m
integer                                  :: auxinput15_interval_s
integer                                  :: auxinput15_interval
integer                                  :: auxinput15_begin_y
integer                                  :: auxinput15_begin_d
integer                                  :: auxinput15_begin_h
integer                                  :: auxinput15_begin_m
integer                                  :: auxinput15_begin_s
integer                                  :: auxinput15_begin
integer                                  :: auxinput15_end_y
integer                                  :: auxinput15_end_d
integer                                  :: auxinput15_end_h
integer                                  :: auxinput15_end_m
integer                                  :: auxinput15_end_s
integer                                  :: auxinput15_end
integer                                  :: io_form_auxinput15
integer                                  :: frames_per_auxinput15
integer                                  :: auxinput16_oid
character*256                               :: auxinput16_inname
character*256                               :: auxinput16_outname
integer                                  :: auxinput16_interval_y
integer                                  :: auxinput16_interval_d
integer                                  :: auxinput16_interval_h
integer                                  :: auxinput16_interval_m
integer                                  :: auxinput16_interval_s
integer                                  :: auxinput16_interval
integer                                  :: auxinput16_begin_y
integer                                  :: auxinput16_begin_d
integer                                  :: auxinput16_begin_h
integer                                  :: auxinput16_begin_m
integer                                  :: auxinput16_begin_s
integer                                  :: auxinput16_begin
integer                                  :: auxinput16_end_y
integer                                  :: auxinput16_end_d
integer                                  :: auxinput16_end_h
integer                                  :: auxinput16_end_m
integer                                  :: auxinput16_end_s
integer                                  :: auxinput16_end
integer                                  :: io_form_auxinput16
integer                                  :: frames_per_auxinput16
integer                                  :: auxinput17_oid
character*256                               :: auxinput17_inname
character*256                               :: auxinput17_outname
integer                                  :: auxinput17_interval_y
integer                                  :: auxinput17_interval_d
integer                                  :: auxinput17_interval_h
integer                                  :: auxinput17_interval_m
integer                                  :: auxinput17_interval_s
integer                                  :: auxinput17_interval
integer                                  :: auxinput17_begin_y
integer                                  :: auxinput17_begin_d
integer                                  :: auxinput17_begin_h
integer                                  :: auxinput17_begin_m
integer                                  :: auxinput17_begin_s
integer                                  :: auxinput17_begin
integer                                  :: auxinput17_end_y
integer                                  :: auxinput17_end_d
integer                                  :: auxinput17_end_h
integer                                  :: auxinput17_end_m
integer                                  :: auxinput17_end_s
integer                                  :: auxinput17_end
integer                                  :: io_form_auxinput17
integer                                  :: frames_per_auxinput17
integer                                  :: auxinput18_oid
character*256                               :: auxinput18_inname
character*256                               :: auxinput18_outname
integer                                  :: auxinput18_interval_y
integer                                  :: auxinput18_interval_d
integer                                  :: auxinput18_interval_h
integer                                  :: auxinput18_interval_m
integer                                  :: auxinput18_interval_s
integer                                  :: auxinput18_interval
integer                                  :: auxinput18_begin_y
integer                                  :: auxinput18_begin_d
integer                                  :: auxinput18_begin_h
integer                                  :: auxinput18_begin_m
integer                                  :: auxinput18_begin_s
integer                                  :: auxinput18_begin
integer                                  :: auxinput18_end_y
integer                                  :: auxinput18_end_d
integer                                  :: auxinput18_end_h
integer                                  :: auxinput18_end_m
integer                                  :: auxinput18_end_s
integer                                  :: auxinput18_end
integer                                  :: io_form_auxinput18
integer                                  :: frames_per_auxinput18
integer                                  :: auxinput19_oid
character*256                               :: auxinput19_inname
character*256                               :: auxinput19_outname
integer                                  :: auxinput19_interval_y
integer                                  :: auxinput19_interval_d
integer                                  :: auxinput19_interval_h
integer                                  :: auxinput19_interval_m
integer                                  :: auxinput19_interval_s
integer                                  :: auxinput19_interval
integer                                  :: auxinput19_begin_y
integer                                  :: auxinput19_begin_d
integer                                  :: auxinput19_begin_h
integer                                  :: auxinput19_begin_m
integer                                  :: auxinput19_begin_s
integer                                  :: auxinput19_begin
integer                                  :: auxinput19_end_y
integer                                  :: auxinput19_end_d
integer                                  :: auxinput19_end_h
integer                                  :: auxinput19_end_m
integer                                  :: auxinput19_end_s
integer                                  :: auxinput19_end
integer                                  :: io_form_auxinput19
integer                                  :: frames_per_auxinput19
integer                                  :: auxinput20_oid
character*256                               :: auxinput20_inname
character*256                               :: auxinput20_outname
integer                                  :: auxinput20_interval_y
integer                                  :: auxinput20_interval_d
integer                                  :: auxinput20_interval_h
integer                                  :: auxinput20_interval_m
integer                                  :: auxinput20_interval_s
integer                                  :: auxinput20_interval
integer                                  :: auxinput20_begin_y
integer                                  :: auxinput20_begin_d
integer                                  :: auxinput20_begin_h
integer                                  :: auxinput20_begin_m
integer                                  :: auxinput20_begin_s
integer                                  :: auxinput20_begin
integer                                  :: auxinput20_end_y
integer                                  :: auxinput20_end_d
integer                                  :: auxinput20_end_h
integer                                  :: auxinput20_end_m
integer                                  :: auxinput20_end_s
integer                                  :: auxinput20_end
integer                                  :: io_form_auxinput20
integer                                  :: frames_per_auxinput20
integer                                  :: auxinput21_oid
character*256                               :: auxinput21_inname
character*256                               :: auxinput21_outname
integer                                  :: auxinput21_interval_y
integer                                  :: auxinput21_interval_d
integer                                  :: auxinput21_interval_h
integer                                  :: auxinput21_interval_m
integer                                  :: auxinput21_interval_s
integer                                  :: auxinput21_interval
integer                                  :: auxinput21_begin_y
integer                                  :: auxinput21_begin_d
integer                                  :: auxinput21_begin_h
integer                                  :: auxinput21_begin_m
integer                                  :: auxinput21_begin_s
integer                                  :: auxinput21_begin
integer                                  :: auxinput21_end_y
integer                                  :: auxinput21_end_d
integer                                  :: auxinput21_end_h
integer                                  :: auxinput21_end_m
integer                                  :: auxinput21_end_s
integer                                  :: auxinput21_end
integer                                  :: io_form_auxinput21
integer                                  :: frames_per_auxinput21
integer                                  :: auxinput22_oid
character*256                               :: auxinput22_inname
character*256                               :: auxinput22_outname
integer                                  :: auxinput22_interval_y
integer                                  :: auxinput22_interval_d
integer                                  :: auxinput22_interval_h
integer                                  :: auxinput22_interval_m
integer                                  :: auxinput22_interval_s
integer                                  :: auxinput22_interval
integer                                  :: auxinput22_begin_y
integer                                  :: auxinput22_begin_d
integer                                  :: auxinput22_begin_h
integer                                  :: auxinput22_begin_m
integer                                  :: auxinput22_begin_s
integer                                  :: auxinput22_begin
integer                                  :: auxinput22_end_y
integer                                  :: auxinput22_end_d
integer                                  :: auxinput22_end_h
integer                                  :: auxinput22_end_m
integer                                  :: auxinput22_end_s
integer                                  :: auxinput22_end
integer                                  :: io_form_auxinput22
integer                                  :: frames_per_auxinput22
integer                                  :: auxinput23_oid
character*256                               :: auxinput23_inname
character*256                               :: auxinput23_outname
integer                                  :: auxinput23_interval_y
integer                                  :: auxinput23_interval_d
integer                                  :: auxinput23_interval_h
integer                                  :: auxinput23_interval_m
integer                                  :: auxinput23_interval_s
integer                                  :: auxinput23_interval
integer                                  :: auxinput23_begin_y
integer                                  :: auxinput23_begin_d
integer                                  :: auxinput23_begin_h
integer                                  :: auxinput23_begin_m
integer                                  :: auxinput23_begin_s
integer                                  :: auxinput23_begin
integer                                  :: auxinput23_end_y
integer                                  :: auxinput23_end_d
integer                                  :: auxinput23_end_h
integer                                  :: auxinput23_end_m
integer                                  :: auxinput23_end_s
integer                                  :: auxinput23_end
integer                                  :: io_form_auxinput23
integer                                  :: frames_per_auxinput23
integer                                  :: auxinput24_oid
character*256                               :: auxinput24_inname
character*256                               :: auxinput24_outname
integer                                  :: auxinput24_interval_y
integer                                  :: auxinput24_interval_d
integer                                  :: auxinput24_interval_h
integer                                  :: auxinput24_interval_m
integer                                  :: auxinput24_interval_s
integer                                  :: auxinput24_interval
integer                                  :: auxinput24_begin_y
integer                                  :: auxinput24_begin_d
integer                                  :: auxinput24_begin_h
integer                                  :: auxinput24_begin_m
integer                                  :: auxinput24_begin_s
integer                                  :: auxinput24_begin
integer                                  :: auxinput24_end_y
integer                                  :: auxinput24_end_d
integer                                  :: auxinput24_end_h
integer                                  :: auxinput24_end_m
integer                                  :: auxinput24_end_s
integer                                  :: auxinput24_end
integer                                  :: io_form_auxinput24
integer                                  :: frames_per_auxinput24
integer                                  :: oid
integer                                  :: history_interval
integer                                  :: frames_per_outfile
logical                                  :: restart
integer                                  :: restart_interval
integer                                  :: io_form_input
integer                                  :: io_form_history
integer                                  :: io_form_restart
integer                                  :: io_form_boundary
integer                                  :: debug_level
logical                                  :: self_test_domain
character*256                               :: history_outname
character*256                               :: history_inname
logical                                  :: use_netcdf_classic
integer                                  :: history_interval_d
integer                                  :: history_interval_h
integer                                  :: history_interval_m
integer                                  :: history_interval_s
integer                                  :: inputout_interval_d
integer                                  :: inputout_interval_h
integer                                  :: inputout_interval_m
integer                                  :: inputout_interval_s
integer                                  :: inputout_interval
integer                                  :: restart_interval_d
integer                                  :: restart_interval_h
integer                                  :: restart_interval_m
integer                                  :: restart_interval_s
integer                                  :: history_begin_y
integer                                  :: history_begin_d
integer                                  :: history_begin_h
integer                                  :: history_begin_m
integer                                  :: history_begin_s
integer                                  :: history_begin
integer                                  :: inputout_begin_y
integer                                  :: inputout_begin_d
integer                                  :: inputout_begin_h
integer                                  :: inputout_begin_m
integer                                  :: inputout_begin_s
integer                                  :: restart_begin_y
integer                                  :: restart_begin_d
integer                                  :: restart_begin_h
integer                                  :: restart_begin_m
integer                                  :: restart_begin_s
integer                                  :: restart_begin
integer                                  :: history_end_y
integer                                  :: history_end_d
integer                                  :: history_end_h
integer                                  :: history_end_m
integer                                  :: history_end_s
integer                                  :: history_end
integer                                  :: inputout_end_y
integer                                  :: inputout_end_d
integer                                  :: inputout_end_h
integer                                  :: inputout_end_m
integer                                  :: inputout_end_s
integer                                  :: simulation_start_year
integer                                  :: simulation_start_month
integer                                  :: simulation_start_day
integer                                  :: simulation_start_hour
integer                                  :: simulation_start_minute
integer                                  :: simulation_start_second
logical                                  :: reset_simulation_start
integer                                  :: sr_x
integer                                  :: sr_y
character*256                               :: sgfdda_inname
character*256                               :: gfdda_inname
integer                                  :: sgfdda_interval_d
integer                                  :: sgfdda_interval_h
integer                                  :: sgfdda_interval_m
integer                                  :: sgfdda_interval_s
integer                                  :: sgfdda_interval_y
integer                                  :: sgfdda_interval
integer                                  :: gfdda_interval_d
integer                                  :: gfdda_interval_h
integer                                  :: gfdda_interval_m
integer                                  :: gfdda_interval_s
integer                                  :: gfdda_interval_y
integer                                  :: gfdda_interval
integer                                  :: sgfdda_begin_y
integer                                  :: sgfdda_begin_d
integer                                  :: sgfdda_begin_h
integer                                  :: sgfdda_begin_m
integer                                  :: sgfdda_begin_s
integer                                  :: gfdda_begin_y
integer                                  :: gfdda_begin_d
integer                                  :: gfdda_begin_h
integer                                  :: gfdda_begin_m
integer                                  :: gfdda_begin_s
integer                                  :: sgfdda_end_y
integer                                  :: sgfdda_end_d
integer                                  :: sgfdda_end_h
integer                                  :: sgfdda_end_m
integer                                  :: sgfdda_end_s
integer                                  :: gfdda_end_y
integer                                  :: gfdda_end_d
integer                                  :: gfdda_end_h
integer                                  :: gfdda_end_m
integer                                  :: gfdda_end_s
integer                                  :: io_form_sgfdda
integer                                  :: io_form_gfdda
character*256                               :: iofields_filename
logical                                  :: ignore_iofields_warning
logical                                  :: ncd_nofill
logical                                  :: update_sfcdiags
logical                                  :: use_wrf_sfcinfo
logical                                  :: use_background_errors
logical                                  :: write_increments
logical                                  :: var4d
integer                                  :: var4d_bin
integer                                  :: var4d_bin_rain
logical                                  :: var4d_lbc
integer                                  :: multi_inc
logical                                  :: print_detail_radar
logical                                  :: print_detail_rain
logical                                  :: print_detail_rad
logical                                  :: print_detail_xa
logical                                  :: print_detail_xb
logical                                  :: print_detail_obs
logical                                  :: print_detail_f_obs
logical                                  :: print_detail_map
logical                                  :: print_detail_grad
logical                                  :: print_detail_regression
logical                                  :: print_detail_spectral
logical                                  :: print_detail_testing
logical                                  :: print_detail_parallel
logical                                  :: print_detail_be
logical                                  :: print_detail_outerloop
logical                                  :: check_max_iv_print
logical                                  :: check_buddy_print
integer                                  :: analysis_accu
logical                                  :: calc_w_increment
logical                                  :: dt_cloud_model
logical                                  :: write_mod_filtered_obs
logical                                  :: wind_sd
logical                                  :: wind_sd_buoy
logical                                  :: wind_sd_synop
logical                                  :: wind_sd_ships
logical                                  :: wind_sd_metar
logical                                  :: wind_sd_sound
logical                                  :: wind_sd_pilot
logical                                  :: wind_sd_airep
logical                                  :: wind_sd_qscat
logical                                  :: wind_sd_tamdar
logical                                  :: wind_sd_geoamv
logical                                  :: wind_sd_mtgirs
logical                                  :: wind_sd_polaramv
logical                                  :: wind_sd_profiler
logical                                  :: wind_stats_sd
logical                                  :: qc_rej_both
integer                                  :: fg_format
integer                                  :: ob_format
integer                                  :: ob_format_gpsro
integer                                  :: num_fgat_time
logical                                  :: thin_conv
logical                                  :: thin_conv_ascii
real                                     :: thin_mesh_conv
logical                                  :: thin_rainobs
logical                                  :: use_synopobs
logical                                  :: use_shipsobs
logical                                  :: use_metarobs
logical                                  :: use_soundobs
logical                                  :: use_mtgirsobs
logical                                  :: use_tamdarobs
logical                                  :: use_pilotobs
logical                                  :: use_airepobs
logical                                  :: use_geoamvobs
logical                                  :: use_polaramvobs
logical                                  :: use_bogusobs
logical                                  :: use_buoyobs
logical                                  :: use_profilerobs
logical                                  :: use_satemobs
logical                                  :: use_gpsztdobs
logical                                  :: use_gpspwobs
logical                                  :: use_gpsrefobs
real                                     :: top_km_gpsro
real                                     :: bot_km_gpsro
logical                                  :: use_ssmiretrievalobs
logical                                  :: use_ssmitbobs
logical                                  :: use_ssmt1obs
logical                                  :: use_ssmt2obs
logical                                  :: use_qscatobs
logical                                  :: use_radarobs
logical                                  :: use_radar_rv
logical                                  :: use_radar_rf
logical                                  :: use_radar_rqv
logical                                  :: use_radar_rhv
logical                                  :: use_3dvar_phy
logical                                  :: use_rainobs
logical                                  :: use_hirs2obs
logical                                  :: use_hirs3obs
logical                                  :: use_hirs4obs
logical                                  :: use_mhsobs
logical                                  :: use_msuobs
logical                                  :: use_amsuaobs
logical                                  :: use_amsubobs
logical                                  :: use_airsobs
logical                                  :: use_airsretobs
logical                                  :: use_eos_amsuaobs
logical                                  :: use_hsbobs
logical                                  :: use_ssmisobs
logical                                  :: use_iasiobs
logical                                  :: use_seviriobs
logical                                  :: use_amsr2obs
logical                                  :: use_kma1dvar
logical                                  :: use_filtered_rad
logical                                  :: use_obs_errfac
logical                                  :: use_atmsobs
logical                                  :: use_mwtsobs
logical                                  :: use_mwhsobs
logical                                  :: check_max_iv
real                                     :: max_error_t
real                                     :: max_error_uv
real                                     :: max_error_spd
real                                     :: max_error_dir
real                                     :: max_omb_spd
real                                     :: max_omb_dir
real                                     :: max_error_pw
real                                     :: max_error_ref
real                                     :: max_error_rh
real                                     :: max_error_q
real                                     :: max_error_p
real                                     :: max_error_tb
real                                     :: max_error_thickness
real                                     :: max_error_rv
real                                     :: max_error_rf
real                                     :: max_error_rain
real                                     :: max_error_buv
real                                     :: max_error_bt
real                                     :: max_error_bq
real                                     :: max_error_slp
logical                                  :: check_buddy
logical                                  :: put_rand_seed
logical                                  :: omb_set_rand
logical                                  :: omb_add_noise
logical                                  :: position_lev_dependant
integer                                  :: obs_qc_pointer
integer                                  :: qmarker_retain
integer                                  :: max_sound_input
integer                                  :: max_mtgirs_input
integer                                  :: max_tamdar_input
integer                                  :: max_synop_input
integer                                  :: max_geoamv_input
integer                                  :: max_polaramv_input
integer                                  :: max_airep_input
integer                                  :: max_satem_input
integer                                  :: max_pilot_input
integer                                  :: max_radar_input
integer                                  :: max_rain_input
integer                                  :: max_metar_input
integer                                  :: max_gpspw_input
integer                                  :: max_ships_input
integer                                  :: max_profiler_input
integer                                  :: max_bogus_input
integer                                  :: max_buoy_input
integer                                  :: max_ssmi_rv_input
integer                                  :: max_ssmi_tb_input
integer                                  :: max_ssmt1_input
integer                                  :: max_ssmt2_input
integer                                  :: max_qscat_input
integer                                  :: max_gpsref_input
integer                                  :: max_airsr_input
integer                                  :: max_tovs_input
integer                                  :: max_ssmis_input
integer                                  :: report_start
integer                                  :: report_end
integer                                  :: tovs_start
integer                                  :: tovs_end
logical                                  :: gpsref_thinning
logical                                  :: outer_loop_restart
integer                                  :: max_ext_its
integer                                  :: ntmax
integer                                  :: nsave
integer                                  :: write_interval
real                                     :: eps
logical                                  :: precondition_cg
real                                     :: precondition_factor
logical                                  :: use_lanczos
logical                                  :: read_lanczos
logical                                  :: write_lanczos
logical                                  :: orthonorm_gradient
integer                                  :: cv_options
integer                                  :: cloud_cv_options
real                                     :: as1
real                                     :: as2
real                                     :: as3
real                                     :: as4
real                                     :: as5
logical                                  :: do_normalize
logical                                  :: use_rf
integer                                  :: rf_passes
real                                     :: var_scaling1
real                                     :: var_scaling2
real                                     :: var_scaling3
real                                     :: var_scaling4
real                                     :: var_scaling5
real                                     :: var_scaling6
real                                     :: var_scaling7
real                                     :: var_scaling8
real                                     :: var_scaling9
real                                     :: var_scaling10
real                                     :: var_scaling11
real                                     :: len_scaling1
real                                     :: len_scaling2
real                                     :: len_scaling3
real                                     :: len_scaling4
real                                     :: len_scaling5
real                                     :: len_scaling6
real                                     :: len_scaling7
real                                     :: len_scaling8
real                                     :: len_scaling9
real                                     :: len_scaling10
real                                     :: len_scaling11
real                                     :: je_factor
real                                     :: power_truncation
logical                                  :: def_sub_domain
real                                     :: x_start_sub_domain
real                                     :: y_start_sub_domain
real                                     :: x_end_sub_domain
real                                     :: y_end_sub_domain
integer                                  :: stdout
integer                                  :: stderr
integer                                  :: trace_unit
integer                                  :: trace_pe
integer                                  :: trace_repeat_head
integer                                  :: trace_repeat_body
integer                                  :: trace_max_depth
logical                                  :: trace_use
logical                                  :: trace_use_frequent
logical                                  :: trace_use_dull
logical                                  :: trace_memory
logical                                  :: trace_all_pes
logical                                  :: trace_csv
logical                                  :: use_html
logical                                  :: warnings_are_fatal
logical                                  :: test_transforms
logical                                  :: test_gradient
logical                                  :: test_statistics
logical                                  :: interpolate_stats
real                                     :: be_eta
logical                                  :: test_dm_exact
integer                                  :: cv_options_hum
integer                                  :: check_rh
real                                     :: set_omb_rand_fac
integer                                  :: seed_array1
integer                                  :: seed_array2
integer                                  :: sfc_assi_options
logical                                  :: psfc_from_slp
logical                                  :: calculate_cg_cost_fn
logical                                  :: lat_stats_option
integer                                  :: interp_option
integer                                  :: balance_type
logical                                  :: use_wpec
real                                     :: wpec_factor
integer                                  :: vert_corr
integer                                  :: vertical_ip
integer                                  :: vert_evalue
real                                     :: max_vert_var1
real                                     :: max_vert_var2
real                                     :: max_vert_var3
real                                     :: max_vert_var4
real                                     :: max_vert_var5
real                                     :: max_vert_var6
real                                     :: max_vert_var7
real                                     :: max_vert_var8
real                                     :: max_vert_var9
real                                     :: max_vert_var10
real                                     :: max_vert_var11
real                                     :: max_vert_var_alpha
real                                     :: psi_chi_factor
real                                     :: psi_t_factor
real                                     :: psi_ps_factor
real                                     :: psi_rh_factor
real                                     :: chi_u_t_factor
real                                     :: chi_u_ps_factor
real                                     :: chi_u_rh_factor
real                                     :: t_u_rh_factor
real                                     :: ps_u_rh_factor
integer                                  :: rttov_emis_atlas_ir
integer                                  :: rttov_emis_atlas_mw
integer                                  :: rtminit_print
integer                                  :: rtminit_nsensor
integer                                  :: rtminit_platform
integer                                  :: rtminit_satid
integer                                  :: rtminit_sensor
integer                                  :: rad_monitoring
real                                     :: thinning_mesh
logical                                  :: thinning
logical                                  :: read_biascoef
logical                                  :: biascorr
logical                                  :: biasprep
logical                                  :: rttov_scatt
logical                                  :: write_profile
logical                                  :: write_jacobian
logical                                  :: qc_rad
logical                                  :: write_iv_rad_ascii
logical                                  :: write_oa_rad_ascii
logical                                  :: write_filtered_rad
logical                                  :: use_error_factor_rad
logical                                  :: use_landem
logical                                  :: use_antcorr
logical                                  :: use_mspps_emis
logical                                  :: use_mspps_ts
integer                                  :: mw_emis_sea
integer                                  :: tovs_min_transfer
logical                                  :: tovs_batch
integer                                  :: rtm_option
logical                                  :: use_crtm_kmatrix
logical                                  :: use_rttov_kmatrix
logical                                  :: crtm_cloud
logical                                  :: only_sea_rad
logical                                  :: use_pseudo_rad
integer                                  :: pseudo_rad_platid
integer                                  :: pseudo_rad_satid
integer                                  :: pseudo_rad_senid
integer                                  :: pseudo_rad_ichan
real                                     :: pseudo_rad_lat
real                                     :: pseudo_rad_lon
real                                     :: pseudo_rad_inv
real                                     :: pseudo_rad_err
logical                                  :: use_simulated_rad
logical                                  :: simulated_rad_io
integer                                  :: simulated_rad_ngrid
logical                                  :: use_varbc
logical                                  :: freeze_varbc
real                                     :: varbc_factor
integer                                  :: varbc_nbgerr
integer                                  :: varbc_nobsmin
logical                                  :: use_clddet_mmr
logical                                  :: use_clddet_ecmwf
logical                                  :: airs_warmest_fov
logical                                  :: use_satcv
logical                                  :: use_blacklist_rad
logical                                  :: calc_weightfunc
character*256                               :: crtm_coef_path
character*256                               :: crtm_irwater_coef
character*256                               :: crtm_mwwater_coef
character*256                               :: crtm_irland_coef
character*256                               :: crtm_visland_coef
integer                                  :: num_pseudo
real                                     :: pseudo_x
real                                     :: pseudo_y
real                                     :: pseudo_z
real                                     :: pseudo_val
real                                     :: pseudo_err
integer                                  :: alphacv_method
integer                                  :: ensdim_alpha
integer                                  :: alpha_truncation
integer                                  :: alpha_corr_type
real                                     :: alpha_corr_scale
real                                     :: alpha_std_dev
logical                                  :: alpha_vertloc
logical                                  :: alpha_hydrometeors
logical                                  :: hybrid_dual_res
integer                                  :: dual_res_upscale_opt
character*256                               :: analysis_type
integer                                  :: sensitivity_option
logical                                  :: adj_sens
character*256                               :: analysis_date
character*256                               :: pseudo_var
character*256                               :: documentation_url
character*256                               :: time_window_min
character*256                               :: time_window_max
logical                                  :: jcdfi_use
integer                                  :: jcdfi_diag
real                                     :: jcdfi_penalty
logical                                  :: enable_identity
logical                                  :: trajectory_io
logical                                  :: var4d_detail_out
logical                                  :: var4d_run
integer                                  :: mp_physics_ad
integer                                  :: mp_physics_4dvar
TYPE(vp_type)                               :: vv
TYPE(vp_type)                               :: vv6
TYPE(vp_type)                               :: vp
TYPE(vp_type)                               :: vp6
TYPE(ep_type)                               :: ep
TYPE(xb_type)                               :: xb
TYPE(x_type)                               :: xa
TYPE(x_type)                               :: x6a
TYPE(xpose_type)                               :: xp
integer                                  :: chem_opt
real      ,DIMENSION(:,:)     ,POINTER   :: lu_index
real      ,DIMENSION(:,:)     ,POINTER   :: lu_mask
real      ,DIMENSION(:)       ,POINTER   :: znu
real      ,DIMENSION(:)       ,POINTER   :: znw
real      ,DIMENSION(:)       ,POINTER   :: zs
real      ,DIMENSION(:)       ,POINTER   :: dzs
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_bxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_bxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_bys
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_bye
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_btxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_btxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_btys
real      ,DIMENSION(:,:,:)   ,POINTER   :: u_btye
real      ,DIMENSION(:,:,:)   ,POINTER   :: ru
real      ,DIMENSION(:,:,:)   ,POINTER   :: ru_m
real      ,DIMENSION(:,:,:)   ,POINTER   :: ru_tend
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_bxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_bxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_bys
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_bye
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_btxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_btxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_btys
real      ,DIMENSION(:,:,:)   ,POINTER   :: v_btye
real      ,DIMENSION(:,:,:)   ,POINTER   :: rv
real      ,DIMENSION(:,:,:)   ,POINTER   :: rv_m
real      ,DIMENSION(:,:,:)   ,POINTER   :: rv_tend
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_bxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_bxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_bys
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_bye
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_btxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_btxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_btys
real      ,DIMENSION(:,:,:)   ,POINTER   :: w_btye
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_bxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_bxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_bys
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_bye
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_btxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_btxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_btys
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph_btye
real      ,DIMENSION(:,:,:)   ,POINTER   :: phb
real      ,DIMENSION(:,:,:)   ,POINTER   :: ph0
real      ,DIMENSION(:,:,:)   ,POINTER   :: php
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_bxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_bxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_bys
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_bye
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_btxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_btxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_btys
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_btye
real      ,DIMENSION(:,:,:)   ,POINTER   :: t_init
real      ,DIMENSION(:,:)     ,POINTER   :: mu_1
real      ,DIMENSION(:,:)     ,POINTER   :: mu_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_bxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_bxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_bys
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_bye
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_btxs
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_btxe
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_btys
real      ,DIMENSION(:,:,:)   ,POINTER   :: mu_btye
real      ,DIMENSION(:,:)     ,POINTER   :: mub
real      ,DIMENSION(:,:)     ,POINTER   :: mu0
real      ,DIMENSION(:,:)     ,POINTER   :: mudf
real      ,DIMENSION(:,:)     ,POINTER   :: muu
real      ,DIMENSION(:,:)     ,POINTER   :: muv
real      ,DIMENSION(:,:)     ,POINTER   :: mut
real      ,DIMENSION(:,:)     ,POINTER   :: muts
real      ,DIMENSION(:,:,:)   ,POINTER   :: p
real      ,DIMENSION(:,:,:)   ,POINTER   :: al
real      ,DIMENSION(:,:,:)   ,POINTER   :: alt
real      ,DIMENSION(:,:,:)   ,POINTER   :: alb
real      ,DIMENSION(:,:,:)   ,POINTER   :: zx
real      ,DIMENSION(:,:,:)   ,POINTER   :: zy
real      ,DIMENSION(:,:,:)   ,POINTER   :: rdz
real      ,DIMENSION(:,:,:)   ,POINTER   :: rdzw
real      ,DIMENSION(:,:,:)   ,POINTER   :: pb
real      ,DIMENSION(:)       ,POINTER   :: fnm
real      ,DIMENSION(:)       ,POINTER   :: fnp
real      ,DIMENSION(:)       ,POINTER   :: rdnw
real      ,DIMENSION(:)       ,POINTER   :: rdn
real      ,DIMENSION(:)       ,POINTER   :: dnw
real      ,DIMENSION(:)       ,POINTER   :: dn
real      ,DIMENSION(:)       ,POINTER   :: t_base
real      ,DIMENSION(:,:,:)   ,POINTER   :: z
real      ,DIMENSION(:,:,:)   ,POINTER   :: z_at_w
real      ,DIMENSION(:,:,:)   ,POINTER   :: p_hyd
real      ,DIMENSION(:,:,:)   ,POINTER   :: p_hyd_w
real      ,DIMENSION(:,:)     ,POINTER   :: q2
real      ,DIMENSION(:,:)     ,POINTER   :: t2
real      ,DIMENSION(:,:)     ,POINTER   :: th2
real      ,DIMENSION(:,:)     ,POINTER   :: psfc
real      ,DIMENSION(:,:)     ,POINTER   :: qsfc
real      ,DIMENSION(:,:)     ,POINTER   :: u10
real      ,DIMENSION(:,:)     ,POINTER   :: v10
integer   ,DIMENSION(:,:)     ,POINTER   :: imask_nostag
integer   ,DIMENSION(:,:)     ,POINTER   :: imask_xstag
integer   ,DIMENSION(:,:)     ,POINTER   :: imask_ystag
integer   ,DIMENSION(:,:)     ,POINTER   :: imask_xystag
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_bxs
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_bxe
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_bys
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_bye
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_btxs
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_btxe
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_btys
real      ,DIMENSION(:,:,:,:) ,POINTER   :: moist_btye
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_bxs
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_bxe
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_bys
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_bye
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_btxs
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_btxe
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_btys
real      ,DIMENSION(:,:,:,:) ,POINTER   :: scalar_btye
real      ,DIMENSION(:)       ,POINTER   :: fcx
real      ,DIMENSION(:)       ,POINTER   :: gcx
real      ,DIMENSION(:,:,:)   ,POINTER   :: smois
real      ,DIMENSION(:,:,:)   ,POINTER   :: sh2o
real      ,DIMENSION(:,:,:)   ,POINTER   :: smcrel
real      ,DIMENSION(:,:)     ,POINTER   :: xice
integer   ,DIMENSION(:,:)     ,POINTER   :: ivgtyp
integer   ,DIMENSION(:,:)     ,POINTER   :: isltyp
real      ,DIMENSION(:,:)     ,POINTER   :: vegfra
real      ,DIMENSION(:,:)     ,POINTER   :: snow
real      ,DIMENSION(:,:)     ,POINTER   :: snowh
real      ,DIMENSION(:,:)     ,POINTER   :: htop
real      ,DIMENSION(:,:)     ,POINTER   :: hbot
real      ,DIMENSION(:,:,:)   ,POINTER   :: tslb
real      ,DIMENSION(:,:)     ,POINTER   :: msft
real      ,DIMENSION(:,:)     ,POINTER   :: msfu
real      ,DIMENSION(:,:)     ,POINTER   :: msfv
real      ,DIMENSION(:,:)     ,POINTER   :: msftx
real      ,DIMENSION(:,:)     ,POINTER   :: msfty
real      ,DIMENSION(:,:)     ,POINTER   :: msfux
real      ,DIMENSION(:,:)     ,POINTER   :: msfuy
real      ,DIMENSION(:,:)     ,POINTER   :: msfvx
real      ,DIMENSION(:,:)     ,POINTER   :: msfvx_inv
real      ,DIMENSION(:,:)     ,POINTER   :: msfvy
real      ,DIMENSION(:,:)     ,POINTER   :: f
real      ,DIMENSION(:,:)     ,POINTER   :: e
real      ,DIMENSION(:,:)     ,POINTER   :: sina
real      ,DIMENSION(:,:)     ,POINTER   :: cosa
real      ,DIMENSION(:,:)     ,POINTER   :: ht
real      ,DIMENSION(:,:)     ,POINTER   :: ht_fine
real      ,DIMENSION(:,:)     ,POINTER   :: ht_int
real      ,DIMENSION(:,:)     ,POINTER   :: ht_input
real      ,DIMENSION(:,:)     ,POINTER   :: tsk
real      ,DIMENSION(:)       ,POINTER   :: u_base
real      ,DIMENSION(:)       ,POINTER   :: v_base
real      ,DIMENSION(:)       ,POINTER   :: qv_base
real      ,DIMENSION(:)       ,POINTER   :: z_base
real      ,DIMENSION(:,:)     ,POINTER   :: rainc
real      ,DIMENSION(:,:)     ,POINTER   :: rainnc
real      ,DIMENSION(:,:)     ,POINTER   :: raincv
real      ,DIMENSION(:,:)     ,POINTER   :: rainncv
real      ,DIMENSION(:,:)     ,POINTER   :: xlat
real      ,DIMENSION(:,:)     ,POINTER   :: xlong
real      ,DIMENSION(:,:)     ,POINTER   :: xlat_u
real      ,DIMENSION(:,:)     ,POINTER   :: xlong_u
real      ,DIMENSION(:,:)     ,POINTER   :: xlat_v
real      ,DIMENSION(:,:)     ,POINTER   :: xlong_v
real      ,DIMENSION(:,:)     ,POINTER   :: albedo
real      ,DIMENSION(:,:)     ,POINTER   :: clat
real      ,DIMENSION(:,:)     ,POINTER   :: albbck
real      ,DIMENSION(:,:)     ,POINTER   :: embck
real      ,DIMENSION(:,:)     ,POINTER   :: emiss
real      ,DIMENSION(:,:)     ,POINTER   :: snotime
real      ,DIMENSION(:,:)     ,POINTER   :: noahres
real      ,DIMENSION(:,:)     ,POINTER   :: cldefi
real      ,DIMENSION(:,:)     ,POINTER   :: tmn
real      ,DIMENSION(:,:)     ,POINTER   :: xland
real      ,DIMENSION(:,:)     ,POINTER   :: snowc
real      ,DIMENSION(:,:)     ,POINTER   :: znt
real      ,DIMENSION(:,:)     ,POINTER   :: ust
real      ,DIMENSION(:,:)     ,POINTER   :: mol
real      ,DIMENSION(:,:)     ,POINTER   :: pblh
real      ,DIMENSION(:,:)     ,POINTER   :: hfx
real      ,DIMENSION(:,:)     ,POINTER   :: qfx
real      ,DIMENSION(:,:)     ,POINTER   :: regime
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_u_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_u_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_u_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_u_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_v_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_v_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_v_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_v_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_w_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_w_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_w_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_w_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_ph_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_ph_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_ph_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_ph_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_t_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_t_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_t_1
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_t_2
real      ,DIMENSION(:,:)     ,POINTER   :: a_mu_1
real      ,DIMENSION(:,:)     ,POINTER   :: a_mu_2
real      ,DIMENSION(:,:)     ,POINTER   :: g_mu_1
real      ,DIMENSION(:,:)     ,POINTER   :: g_mu_2
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_p
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_p
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_z
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_z
real      ,DIMENSION(:)       ,POINTER   :: kma_a
real      ,DIMENSION(:)       ,POINTER   :: kma_b
real      ,DIMENSION(:,:,:,:) ,POINTER   :: a_moist
real      ,DIMENSION(:,:,:,:) ,POINTER   :: g_moist
real      ,DIMENSION(:,:,:)   ,POINTER   :: g_h_diabatic
real      ,DIMENSION(:,:,:)   ,POINTER   :: a_h_diabatic
real      ,DIMENSION(:,:)     ,POINTER   :: g_rainc
real      ,DIMENSION(:,:)     ,POINTER   :: a_rainc
real      ,DIMENSION(:,:)     ,POINTER   :: g_rainnc
real      ,DIMENSION(:,:)     ,POINTER   :: a_rainnc
real      ,DIMENSION(:,:)     ,POINTER   :: g_raincv
real      ,DIMENSION(:,:)     ,POINTER   :: a_raincv
real      ,DIMENSION(:,:)     ,POINTER   :: g_rainncv
real      ,DIMENSION(:,:)     ,POINTER   :: a_rainncv
real      ,DIMENSION(:)       ,POINTER   :: dummy
real      ,DIMENSION(:,:,:,:) ,POINTER   :: a_scalar
real      ,DIMENSION(:,:,:,:) ,POINTER   :: g_scalar
real      ,DIMENSION(:,:)     ,POINTER   :: landmask
real      ,DIMENSION(:,:)     ,POINTER   :: sst
real      ,DIMENSION(:,:,:,:) ,POINTER   :: chem
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_bxs
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_bxe
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_bys
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_bye
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_btxs
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_btxe
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_btys
real      ,DIMENSION(:,:,:,:) ,POINTER   :: tracer_btye
!ENDOFREGISTRYGENERATEDINCLUDE
