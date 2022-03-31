












module da_setup_structures

   
   
   

   use da_wavelet, only: lf,namw,nb,nij,ws
   use module_domain, only : xb_type, ep_type, domain

   use da_define_structures, only : xbx_type,be_subtype, be_type, y_type, j_type, &
      iv_type,da_allocate_background_errors,da_allocate_observations, &
      multi_level_type,each_level_type, da_allocate_observations_rain
   use da_wrf_interfaces, only : wrf_debug
   use da_control, only : trace_use,vert_evalue,stdout,rootproc, &
      analysis_date,coarse_ix,coarse_ds,map_projection,coarse_jy, c2,dsm,phic, &
      pole, cone_factor, start_x,base_pres,ptop,psi1,start_y, base_lapse,base_temp,truelat2_3dv, &
      truelat1_3dv,xlonc,t0,num_fft_factors,pi,print_detail_spectral, global, print_detail_obs, &
      use_radar_rf, use_radar_rhv, use_radar_rqv, use_3dvar_phy, &
      num_ob_indexes,kts, kte, time_window_max, time_window_min, &
      max_fgat_time, num_fgat_time, dt_cloud_model, &
      use_ssmiretrievalobs,use_radarobs,use_ssmitbobs,use_qscatobs, num_procs, use_rainobs, &
      num_pseudo, missing, ob_format, ob_format_bufr,ob_format_ascii, ob_format_madis, ob_format_gpsro, &
      use_airepobs, use_tamdarobs, test_dm_exact, use_amsuaobs, use_amsubobs, &
      use_airsobs, use_bogusobs, sfc_assi_options, use_eos_amsuaobs, &
      use_filtered_rad, use_gpsrefobs, use_hirs2obs, &
      use_hsbobs,use_hirs3obs, use_gpspwobs, use_gpsztdobs, use_metarobs, use_msuobs, &
      use_kma1dvar,use_pilotobs, use_polaramvobs, use_rad, crtm_cloud, use_soundobs,use_mtgirsobs, &
      use_ssmt1obs,use_ssmt2obs, use_shipsobs, use_satemobs, use_synopobs, &
      use_radar_rv,use_profilerobs, use_obsgts, use_geoamvobs, use_buoyobs, &
      jb_factor, je_factor, alphacv_method,its,ite,jts,jte,cv_size_domain_jb, cv_size_domain_jl, &
      cv_size_domain_je, cv_size_domain,ensdim_alpha, alpha_vertloc, alpha_hydrometeors, &
      lat_stats_option,alpha_std_dev,sigma_alpha,alpha_corr_scale, &
      len_scaling1, len_scaling2, len_scaling3, len_scaling4, len_scaling5,&
      len_scaling6, len_scaling7, len_scaling8, len_scaling9, &
      len_scaling10, len_scaling11, &
      max_vert_var1, max_vert_var2, max_vert_var3, max_vert_var4, max_vert_var5, &
      max_vert_var6, max_vert_var7, max_vert_var8, max_vert_var9, max_vert_var10,&
      max_vert_var11, max_vert_var_alpha, &
      print_detail_be, test_statistics, do_normalize, use_rf, &
      var_scaling1, var_scaling2, var_scaling3, var_scaling4, &
      var_scaling5, var_scaling6, var_scaling7, var_scaling8, &
      var_scaling9, var_scaling10, var_scaling11,&
      vert_corr,max_vert_var5,power_truncation,alpha_truncation, &
      print_detail_regression,gas_constant, use_airsretobs, &
      filename_len, use_ssmisobs, gravity, t_triple, use_hirs4obs, use_mhsobs, &
      use_mwtsobs, use_mwhsobs, use_atmsobs,    &
      vert_corr_2, alphacv_method_xa, vert_evalue_global, &
      vert_evalue_local, obs_names, thin_conv, thin_conv_ascii, &
      sound, sonde_sfc, mtgirs, tamdar, tamdar_sfc, synop, profiler, gpsref, gpspw, polaramv, geoamv, ships, metar, &
      satem, radar, ssmi_rv, ssmi_tb, ssmt1, ssmt2, airsr, pilot, airep, rain, &
      bogus, buoy, qscat, radiance, pseudo, trace_use_dull, kts,kte, &
      use_simulated_rad, use_pseudo_rad, pseudo_rad_platid, pseudo_rad_satid, &
      pseudo_rad_senid, rtminit_nsensor, rtminit_platform, rtminit_satid, &
      rtminit_sensor, thinning, qc_rad, var4d, & 
      num_pseudo,pseudo_x, pseudo_y, pseudo_z, pseudo_var,pseudo_val, pseudo_err,&
      fg_format, fg_format_wrf_arw_regional,fg_format_wrf_nmm_regional, &
      fg_format_wrf_arw_global, fg_format_kma_global, deg_to_rad, rad_to_deg, &
      sonde_sfc, missing_data, missing_r, qc_good, thin_mesh_conv, time_slots, &
      cv_options, cloud_cv_options, cv_size, as1, as2, as3, as4, as5, print_detail_be, &
      ids,ide,jds,jde,kds,kde, ims,ime,jms,jme,kms,kme, &
      its,ite,jts,jte,kts,kte, ips,ipe,jps,jpe,kps,kpe, root, comm, ierr, &
      fmt_info, fmt_srfc, fmt_each, unit_end, max_ext_its, &  
      psi_chi_factor, psi_t_factor, psi_ps_factor, psi_rh_factor, &
      chi_u_t_factor, chi_u_ps_factor,chi_u_rh_factor, t_u_rh_factor, ps_u_rh_factor, &
      interpolate_stats, be_eta, thin_rainobs, fgat_rain_flags, use_iasiobs, &
      use_seviriobs, jds_int, jde_int, anal_type_hybrid_dual_res, use_amsr2obs, nrange

   use da_obs, only : da_fill_obs_structures, da_store_obs_grid_info, da_store_obs_grid_info_rad, &
                      da_fill_obs_structures_rain, da_fill_obs_structures_radar, da_set_obs_missing,da_set_3d_obs_missing
   use da_obs_io, only : da_read_obs_bufr,da_read_obs_radar, &
      da_scan_obs_radar,da_scan_obs_ascii,da_read_obs_ascii, &
      da_read_obs_bufrgpsro, da_scan_obs_rain, da_read_obs_rain
   use da_par_util1, only : da_proc_sum_real, da_proc_sum_int, da_proc_sum_ints
   use da_par_util, only : da_patch_to_global
   use da_lapack, only : dsyev
   use da_radiance, only : da_setup_radiance_structures
   use da_reporting, only : da_error, message, da_warning, da_message
   use da_recursive_filter, only : da_calculate_rf_factors
   use da_spectral, only : da_initialize_h,da_calc_power_spectrum
   use da_ssmi, only : da_read_obs_ssmi,da_scan_obs_ssmi
   use da_tools_serial, only : da_get_unit, da_free_unit, da_array_print, da_find_fft_factors, &
      da_find_fft_trig_funcs
   use da_tools, only: da_get_time_slots, da_1d_eigendecomposition
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_vtox_transforms, only : da_check_eof_decomposition
   use da_rfz_cv3, only : da_rfz0
   use da_rf_cv3, only : RFDPAR1, RFDPAR2, RFDPARV
   use module_radiance, only : init_constants_derived
   use gsi_thinning, only : r999,r360,rlat_min,rlat_max,rlon_min,rlon_max, &
                            dlat_grid,dlon_grid,thinning_grid_conv,thinning_grid, &
                            make3grids, makegrids, destroygrids, destroygrids_conv, cleangrids_conv

   use da_par_util, only : true_mpi_real

   implicit none

   include 'mpif.h'

contains

subroutine da_get_vertical_truncation( max_vert_var, eigenval, be_sub)

   !---------------------------------------------------------------------------
   !  Purpose: Calculate vertical mode truncation from explained variance.
   !
   !  Method:  Calculate cumulative variance and compare with limit.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)                 :: max_vert_var    ! Vertical variance limit.
   real*8, intent(in)               :: eigenval(:)     ! Global eigenvaluess.
   type(be_subtype), intent(inout) :: be_sub          ! Back. error sub type.
    
   integer                          :: kz              ! # vertical levels.
   integer                          :: k               ! Loop counter.
   real*8                           :: tot_variance    ! Total variance.
   real*8                           :: cum_variance    ! Cumulative variance.
   character(LEN = 6)              :: name            ! Variable name.

   if (trace_use_dull) call da_trace_entry("da_get_vertical_truncation")


   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   kz = size( eigenval(:))
   name = trim(be_sub % name)

   !---------------------------------------------------------------------------
   ! [2.0] Calculate vertical truncation:
   !---------------------------------------------------------------------------

   if (max_vert_var >= 100.0) then
   
      ! [2.1] No truncation: 
      be_sub % mz = kz

      ! Disregard zero/-ve eigenvalues(which should be very small and only 
      ! appear if statistics have been interpolated between domains):

      do k = 1, kz
         if (eigenval(k) <= 0.0) then
            be_sub % mz = k - 1
            exit
         end if
      end do      
   else
   
      ! [2.2] Calculate cumulative variance and truncate:

      tot_variance = sum( eigenval(1:kz))
      cum_variance = 0.0
      
      do k = 1, kz
         cum_variance = cum_variance + eigenval(k)
         
         if (eigenval(k) <= 0.0) then
            be_sub % mz = k - 1
            exit
         end if
         
         if (cum_variance/tot_variance >= 0.01 * max_vert_var) then
            be_sub % mz = k
            exit
         end if  
      end do
      
      if (max_vert_var == 0.0) be_sub % mz = 0 

   end if

   write(unit=stdout,fmt='(A,A6,A3,I3,A1,f7.2,A2)') &
      'Vertical truncation for ', name, &
      ' = ', be_sub % mz, '(', &
      max_vert_var, '%)'
   write (unit=stdout,fmt='(A)') " "

   if (trace_use_dull) call da_trace_exit("da_get_vertical_truncation")
                                       
end subroutine da_get_vertical_truncation


subroutine da_interpolate_regcoeff (iy, iys, kz, kzs, meanl_stats, meanl_xb, meanp_stats, meanp_xb, &
   pb_vert_reg_stats, pb_vert_reg)

   !---------------------------------------------------------------------------
   ! Purpose: Interpolate statistical regression coefficient to new domain.
   !
   ! Method:  i,k,k Interpolation.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: iy                       ! Number of rows in  xb.
   integer, intent(in)  :: iys                      ! Number of rows in stats.
   integer, intent(in)  :: kz                       ! Number of levels in xb.
   integer, intent(in)  :: kzs                      ! Number of levels in stats.
   real,    intent(in)  :: meanl_stats(:)           ! Mean latitude on stats rows.
   real,    intent(in)  :: meanl_xb(:)              ! Mean latitude on xb rows.
   real,    intent(in)  :: meanp_stats(:)           ! Mean pressure on stats levs.
   real,    intent(in)  :: meanp_xb(:)              ! Mean pressure on xb levs.
   real*8,  intent(in)  :: pb_vert_reg_stats(:,:,:) ! Coefficient on stats grid.
   real*8,  intent(out) :: pb_vert_reg(:,:,:)       ! Coefficient on xb grid.
     
   integer :: i, is, i_south           ! Loop counters.
   integer :: k1, k2, k, ks            ! Loop counters.
   integer :: k1s, k2s
   real    :: lat_wgt

   integer :: k_above(1:kz)
   real    :: pb_vert_reg_temp(1:iys,1:kz,1:kz)
   real    :: weight(1:kz)

   if (trace_use) call da_trace_entry("da_interpolate_regcoeff")

   pb_vert_reg = 0.0

   !------------------------------------------------------------------------
   ! [1.0] Find xb levels/rows bounded by stats domain:
   !------------------------------------------------------------------------

   do k = 1, kz
      if (meanp_xb(k) <= meanp_stats(1)) then
         weight(k) = 1.0e-6
         k_above(k) = 1
      else if (meanp_xb(k) >= meanp_stats(kzs)) then
         weight(k) = 1.0-1.0e-6
         k_above(k) = kzs-1
      else
         do ks = 1, kzs-1
            if (meanp_xb(k) >= meanp_stats(ks) .AND. meanp_xb(k) <= meanp_stats(ks+1)) then
               weight(k) = (meanp_xb(k) - meanp_stats(ks)) / (meanp_stats(ks+1) - meanp_stats(ks))
               k_above(k) = ks
               exit
            end if
         end do
      end if
   end do

   !------------------------------------------------------------------------   
   ! [3.0] Interpolate regression coefficient from stats to xb levels:
   !------------------------------------------------------------------------

   pb_vert_reg_temp(1:iys,1:kz,1:kz) = 0.0

   do is = 1, iys
      do k1 = 1, kz
         k1s = k_above(k1)
         do k2 = 1, kz
            k2s = k_above(k2)

            pb_vert_reg_temp(is,k1,k2) = (1.0-weight(k1)) * (1.0-weight(k2)) * &
                                         pb_vert_reg_stats(is,k1s,k2s) + &
                                         weight(k1) * (1.0-weight(k2)) * &
                                         pb_vert_reg_stats(is,k1s+1,k2s) + &
                                         weight(k2) * (1.0-weight(k1)) * &
                                         pb_vert_reg_stats(is,k1s,k2s+1) + &
                                         weight(k2) * weight(k1) * &
                                         pb_vert_reg_stats(is,k1s+1,k2s+1)
         end do
      end do     
   end do
         
   !------------------------------------------------------------------------
   ! [4.0] Interpolate to from statistics latitudes to xb latitudes:
   !------------------------------------------------------------------------

   i_south = 2

   do i = 1, iy
   
      ! Find position of xb latitude in statistics rows:

      if (meanl_xb(i) <= meanl_stats(2)) then
         i_south = 2
         lat_wgt = 0.5
      else if (meanl_xb(i) >= meanl_stats(iys-1)) then
         i_south = iys-2
         lat_wgt = 0.5
      else
         do is = 1, iys-1
            if (meanl_xb(i) >= meanl_stats(is) .AND. meanl_xb(i) <= meanl_stats(is+1)) then

               lat_wgt = (meanl_xb(i) - meanl_stats(is)) / (meanl_stats(is+1) - meanl_stats(is))
               i_south = is
               exit
            end if
         end do
      end if
   
      do k1 = 1, kz
         do k2 = 1, kz
            pb_vert_reg(i,k1,k2) = lat_wgt * pb_vert_reg_temp(i_south+1,k1,k2) + &
               (1.0 - lat_wgt) * pb_vert_reg_temp(i_south,k1,k2)
         end do
      end do     
   end do

   if (print_detail_regression) then
      call da_array_print(1, pb_vert_reg_stats(1,:,:), 'pb_vert_reg_stats(1,:,:)')
      call da_array_print(1, pb_vert_reg(1,:,:),       'pb_vert_reg(1,:,:)')
      call da_array_print(1, pb_vert_reg_stats(:,1,:), 'pb_vert_reg_stats(:,1,:)')
      call da_array_print(1, pb_vert_reg(:,1,:),       'pb_vert_reg(:,1,:)')
   end if

   if (trace_use) call da_trace_exit("da_interpolate_regcoeff")
   
end subroutine da_interpolate_regcoeff


subroutine da_rescale_background_errors (var_scaling, len_scaling, &
                                          ds, s, be_sub)

   !---------------------------------------------------------------------------
   ! Purpose: Rescale wrfvar background errors.
   !
   ! Method:  Empirical scaling and inclusion of recursive filter rescaling.
   !---------------------------------------------------------------------------

   implicit none

   real, intent(in)                 :: var_scaling       ! Variance factor.
   real, intent(in)                 :: len_scaling       ! Lengthscale factor.
   real, intent(in)                 :: ds                ! Resolution (m)
   real*8, intent(inout)            :: s(:)              ! RF lengthscale.
   type (be_subtype), intent(inout) :: be_sub            ! Backgrd error struct.
    
   integer                          :: mz                ! Vertical truncation.
   integer                          :: m
   real*8, allocatable              :: rf_scale_factor(:)! RF rescaling.

   if (trace_use_dull) call da_trace_entry("da_rescale_background_errors")

   write(unit=stdout,fmt='(3x,"Scaling: var, len, ds:",3e15.6 )') &
                                     var_scaling, len_scaling, ds

   !--------------------------------------------------------------------------
   ! [1.0] Initialise:
   !--------------------------------------------------------------------------

   mz = be_sub % mz

   !--------------------------------------------------------------------------
   ! [2.0] Perform various rescalings:
   !--------------------------------------------------------------------------

   if (mz > 0) then

      ! [2.1] Empirical rescaling of lengthscales:
      s(1:mz) = len_scaling * s(1:mz)
   
      if (print_detail_be) then
         write(unit=stdout,fmt='(a,a)')trim(be_sub % name), ' Error Lengthscales (m):'
         do m = 1, mz
            write(unit=stdout,fmt='(a,i4,1pe13.5)')be_sub % name, m, s(m)
         end do
      end if
      
      ! [2.2] Make lengthscale nondimensional:
      s(1:mz) = s(1:mz) / ds

      ! [2.3] Empirical rescaling of variances:
      be_sub % val(:,:) = var_scaling * be_sub % val(:,:)

      ! Calculate recursive filter rescaling:

      allocate(rf_scale_factor(1:mz))

      call da_calculate_rf_factors(s(:), be_sub % rf_alpha(:), &
                                    rf_scale_factor(:))

      do m = 1, mz
         be_sub % val(:,m) = rf_scale_factor(m) * be_sub % val(:,m)
      end do
                                       
      deallocate (rf_scale_factor)   

   end if

   if (trace_use_dull) call da_trace_exit("da_rescale_background_errors")

end subroutine da_rescale_background_errors


subroutine da_scale_background_errors ( be, it )
!
   TYPE (be_type), INTENT(INOUT) :: be     ! Back. errors structure
   INTEGER,        INTENT(IN)    :: it     ! outer-loop index
!
   real, allocatable, dimension(:,:) :: v1_val , v2_val , v3_val , &
                                        v4_val , v5_val
   real*8, allocatable, dimension(:) :: rf_len1, rf_len2, rf_len3, &
                                        rf_len4, rf_len5
!
   integer                     :: be_rf_unit, be_print_unit
   integer  :: i, ix, jy, kz, v1_mz, v2_mz, v3_mz, v4_mz, v5_mz
   real     :: ds

   if ( jb_factor <= 0.0 ) return
!
! Rewind the unit:
    be_rf_unit    = unit_end + 1
    be_print_unit = unit_end + 2
    rewind (be_rf_unit)
!
! Read the dimensions and allocate the arrays:
    read(be_rf_unit) kz, jy, ix, v1_mz, v2_mz, v3_mz, v4_mz, v5_mz, ds
!
    allocate ( v1_val (1:jy,1:v1_mz) )
    allocate ( rf_len1(1:kz) )
    allocate ( v2_val (1:jy,1:v2_mz) )
    allocate ( rf_len2(1:kz) )
    allocate ( v3_val (1:jy,1:v3_mz) )
    allocate ( rf_len3(1:kz) )
    allocate ( v4_val (1:jy,1:v4_mz) )
    allocate ( rf_len4(1:kz) )
    allocate ( v5_val (1:jy,1:v5_mz) )
    allocate ( rf_len5(1:1) )
!
! Read the variances and scale-lengths and restore them to be:
    read(be_rf_unit) v1_val , v2_val , v3_val , v4_val , v5_val , &
                     rf_len1, rf_len2, rf_len3, rf_len4, rf_len5
!
    be % v1 % val = v1_val
    be % v2 % val = v2_val
    be % v3 % val = v3_val
    be % v4 % val = v4_val
    be % v5 % val = v5_val
!
! Rescale the scale-lengths and variances:
   CALL da_rescale_background_errors( var_scaling1(it), len_scaling1(it), &
                                      ds, rf_len1, be % v1 )
!  .........................................................    
   CALL da_rescale_background_errors( var_scaling2(it), len_scaling2(it), &
                                      ds, rf_len2, be % v2 )
! ..........................................................
   CALL da_rescale_background_errors( var_scaling3(it), len_scaling3(it), &
                                      ds, rf_len3, be % v3 )
! ...............................................................
   CALL da_rescale_background_errors( var_scaling4(it), len_scaling4(it), &
                                      ds, rf_len4, be % v4 )
! ..............................................................
   CALL da_rescale_background_errors( var_scaling5(it), len_scaling5(it), &
                                      ds, rf_len5, be % v5 )
!
! Print the variances and RF (Recursive Filter) factors rf_alpha:
    write(unit=stdout,fmt='(/5x,"Complete the Rescale BES in outer-loop:" i2)') it
    
    if ( print_detail_be ) then
       write(be_print_unit,'(/"============================================================")')
       write(be_print_unit,'("For outer loop ",i2)') it
       write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                                      it, kz, jy, ix, ds
       write(be_print_unit,'("Namelist options specified for this iteration:")')
       write(be_print_unit,'("var_scaling1(it) = ",e12.5,2x,"len_scaling1(it) = "e12.5)')var_scaling1(it),len_scaling1(it)
       write(be_print_unit,'("var_scaling2(it) = ",e12.5,2x,"len_scaling2(it) = "e12.5)')var_scaling2(it),len_scaling2(it)
       write(be_print_unit,'("var_scaling3(it) = ",e12.5,2x,"len_scaling3(it) = "e12.5)')var_scaling3(it),len_scaling3(it)
       write(be_print_unit,'("var_scaling4(it) = ",e12.5,2x,"len_scaling4(it) = "e12.5)')var_scaling4(it),len_scaling4(it)
       write(be_print_unit,'("var_scaling5(it) = ",e12.5,2x,"len_scaling5(it) = "e12.5)')var_scaling5(it),len_scaling5(it)
       write(be_print_unit,'("Background error statistics for this iteration:")')
       write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)
       write(be_print_unit,'("be%v1%rf_alpha:"/(10e12.5))') be % v1 % rf_alpha(:)
       write(be_print_unit,'("be%v2%rf_alpha:"/(10e12.5))') be % v2 % rf_alpha(:)
       write(be_print_unit,'("be%v3%rf_alpha:"/(10e12.5))') be % v3 % rf_alpha(:)
       write(be_print_unit,'("be%v4%rf_alpha:"/(10e12.5))') be % v4 % rf_alpha(:)
       write(be_print_unit,'("be%v5%rf_alpha:"/(10e12.5))') be % v5 % rf_alpha(:)
       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       do i = 1,kz 
          if (i == 1) then
             write(be_print_unit,'(i3,2x,5e15.5)') i, rf_len1(i), rf_len2(i), rf_len3(i), rf_len4(i), rf_len5(i)
          else
             write(be_print_unit,'(i3,2x,4e15.5)') i, rf_len1(i), rf_len2(i), rf_len3(i), rf_len4(i)
          endif
       enddo
    endif
!
! Deallocate the arrays:
    deallocate ( v1_val )
    deallocate ( rf_len1 )
    deallocate ( v2_val )
    deallocate ( rf_len2 )
    deallocate ( v3_val )
    deallocate ( rf_len3 )
    deallocate ( v4_val )
    deallocate ( rf_len4 )
    deallocate ( v5_val )
    deallocate ( rf_len5 )
!
end subroutine da_scale_background_errors

subroutine da_scale_background_errors_cv3 ( grid, be, it )
!
   type (domain), intent(in)     :: grid
   TYPE (be_type), INTENT(INOUT) :: be     ! Back. errors structure
   INTEGER,        INTENT(IN)    :: it     ! outer-loop index
!
   REAL, DIMENSION( ids: ide,  jds: jde,  kds: kde, 1:4) :: hwll
   REAL, DIMENSION( ids: ide,  jds: jde)                 :: hwllp, global_fac
   REAL, DIMENSION( kts: kte,  kts: kte)                 :: vv
   integer                     :: n, i, j, k, ic, jc, ii, ij, ijk, &
                                  iis, iie, jjs, jje, kks, kke
   real, dimension(1)          :: xsum
   real, dimension(ids: ide,  jds: jde)                   :: global_2d
   real, dimension(ids: ide,  jds: jde,  kds: kde)        :: global_3d
   character(len=6)            :: vname
   character(len=9)            :: chr
   character(len=8)            :: i_char

   INTEGER :: nta, ndeg, ku, kz
   real    :: samp,s2u,tin,as(5),as0(5),slim
   character(len=256) :: mesg

   integer                     :: ier, be_rf_unit, be_print_unit
!
! 1, get global_fac from grid:

   ij = ( ide- ids+1)*( jde- jds+1)
   ijk = ij * ( kde - kds + 1)
!  Collect xb component of fac into global buffer.
   call da_patch_to_global( grid, grid%xb%map_factor, global_fac )
   call wrf_dm_bcast_real( global_fac, ij )

! 2, Read in the variance and horizontal scales
!
! 2.1 assign the units for reading and printing:

! Rewind the unit:
    be_rf_unit    = unit_end + 3
    be_print_unit = unit_end + 5
    if (rootproc) rewind (be_rf_unit)

    if (rootproc .and. print_detail_be) &
    write(be_print_unit,'(//10x,a,i2,a)') '====== SCALE CV3 BE for IT =', it, ' ======' 
!
! 2.2 Read the dimensions
!
    if (rootproc) &
       read(be_rf_unit) iis, iie, jjs, jje, kks, kke, ku, samp
    call wrf_dm_bcast_integer( iis , 1 )
    call wrf_dm_bcast_integer( iie , 1 )
    call wrf_dm_bcast_integer( jjs , 1 )
    call wrf_dm_bcast_integer( jje , 1 )
    call wrf_dm_bcast_integer( kks , 1 )
    call wrf_dm_bcast_integer( kke , 1 )
    call wrf_dm_bcast_integer( ku , 1 )
    call wrf_dm_bcast_real( samp , 1 )
    if (rootproc .and. print_detail_be) &
    write (be_print_unit,'(/a,7i5,f15.5)') &
              "Read in: ids, ide, jds, jde, kds, kde, ku, samp:", &
                        iis, iie, jjs, jje, kks, kke, ku, samp
!
! 2.3 Read in the variances before normalization and write to "be_print_unit":

!    write (*,'(3x,a)')  'READ IN VARIANCE BEFORE NORMALIZATION:'
! Print out the sum(x*X) for verification:
!
   if (rootproc .and. print_detail_be) &
   write (be_print_unit,'(3x,a)')  'VARIANCE:'
   do ii = 1, 4
     if (rootproc) &
        read (be_rf_unit) chr, vname, i, global_3d
     call wrf_dm_bcast_string( chr, 9 )
     call wrf_dm_bcast_string( vname, 6 )
     call wrf_dm_bcast_integer( i, 1 )
     call wrf_dm_bcast_real( global_3d, ijk )
     be%corz(its:ite,jts:jte,kts:kte,i) = global_3d(its:ite,jts:jte,kts:kte)
     xsum(1) = sum (be%corz(its:ite,jts:jte,kts:kte,i)*be%corz(its:ite,jts:jte,kts:kte,i))
     call da_proc_sum_real(xsum)
     if (rootproc .and. print_detail_be) &
     write (be_print_unit,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, vname, xsum(1)
!     write (*,'(3x,a9,2x,a6,2x,i3)') chr, vname, i
   enddo

!  Psfc Variance before the normalization.
   if (rootproc) &
      read (be_rf_unit) chr, vname, ii, global_2d
   call wrf_dm_bcast_string( chr, 9 )
   call wrf_dm_bcast_string( vname, 6 )
   call wrf_dm_bcast_integer( ii, 1 )
   call wrf_dm_bcast_real( global_2d, ij )
   be%corp(its:ite,jts:jte) = global_2d(its:ite,jts:jte)
   xsum(1) = sum (be%corp(its:ite,jts:jte)*be%corp(its:ite,jts:jte))
   call da_proc_sum_real(xsum)
   if (rootproc .and. print_detail_be) &
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'PSFC_u', xsum(1)
!   write (*,'(3x,a9,2x,a6,2x,i3)') chr, vname, ii
!
! 2.4 Read in the vertical scales to "be_print_unit":
!
!   write (*,'(3x,a)')  'READ IN VERTICAL SCALE BEFORE TUNING:'
!
   if (rootproc .and. print_detail_be) &
   write (be_print_unit,'(3x,a)')  'VERTICAL SCALE:'
   do ii = 1, 4
     if (rootproc) &
         read (be_rf_unit) chr, vname, n, global_3d
     call wrf_dm_bcast_string( chr, 9 )
     call wrf_dm_bcast_string( vname, 6 )
     call wrf_dm_bcast_integer( n, 1 )
     call wrf_dm_bcast_real( global_3d, ijk )
     do i = its, ite
     do j = jts, jte
     do k = kts, kte
        be%vz(k,i,j,n) = global_3d(i,j,k)
     enddo
     enddo
     enddo
     xsum(1) = sum (be%vz(kts:kte,its:ite,jts:jte,n)*be%vz(kts:kte,its:ite,jts:jte,n))
     call da_proc_sum_real(xsum)
     if (rootproc .and. print_detail_be) &
     write (be_print_unit,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, vname, xsum(1)
!     write (*,'(3x,a9,2x,a6,2x,i3)') chr, vname, n
   enddo
!
! 2.5 Read in the Horizontal scales, and write out vertical scales to "be_print_unit":

   if (rootproc .and. print_detail_be) &
   write (be_print_unit,'(3x,a)')  'READ IN THE HORIZONTAL SCALES:'
   if (rootproc) &
      read (be_rf_unit) hwll
   call wrf_dm_bcast_real( hwll, 4*ijk )
   if (rootproc .and. print_detail_be) then
   write (be_print_unit,'(3x,a)')  'HORIZONTAL SCALES:'
   xsum(1) = sum (hwll*hwll)
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PSI, CHI_u, T_u, RH_u SCALES:', xsum(1)
   endif

   if (rootproc) &
      read (be_rf_unit) hwllp
   call wrf_dm_bcast_real( hwllp, ij )
   if (rootproc .and. print_detail_be) then
   xsum(1) = sum (hwllp*hwllp)
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PS_u SCALES:', xsum(1)
!   write (*,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PS_u SCALES:', xsum(1)
   endif

! 2.6 Write out the regression coefficients to "be_print_unit":

   if (print_detail_be) then

   if (rootproc) &
   write (be_print_unit,'(3x,a)')  'REGRESSION COEFF. T, CHI, PSFC:'
   do ii = kts, kte
     xsum(1) = sum (be%agvz(:,:,:,ii)*be%agvz(:,:,:,ii))
     call da_proc_sum_real(xsum)
     if (rootproc) &
     write (be_print_unit,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, 'TMP-PSI', xsum(1)
   enddo

!  Rg. Coeff. already stored in "be" data structure, no need to save.
   xsum(1) = sum (be%bvz*be%bvz)
   call da_proc_sum_real(xsum)
   if (rootproc) &
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'CHI-PSI', xsum(1)

!  Rg. Coeff. already stored in "be" data structure, no need to save.
   xsum(1) = sum (be%wgvz*be%wgvz)
   call da_proc_sum_real(xsum)
   if (rootproc) &
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'PSF-PSI', xsum(1)
!
   endif

! 3.0 Re-scale the scale-lengths
!
! 3.1 horizontal scales

   as0(1)=as1(2)
   as0(2)=as2(2)
   as0(3)=as3(2)
   as0(4)=as4(2)

   as(1)=as1(2+(it-1)*3)
   as(2)=as2(2+(it-1)*3)
   as(3)=as3(2+(it-1)*3)
   as(4)=as4(2+(it-1)*3)

   do n=1, 4
      hwll(:,:,:,n) = hwll(:,:,:,n) * as(n) / as0(n)
   enddo

   hwllp = hwllp * as5(2+(it-1)*3) / as5(2)

! 3.2 vertical scales

   as0(1)=as1(3)
   as0(2)=as2(3)
   as0(3)=as3(3)
   as0(4)=as4(3)

   as(1)=as1(3+(it-1)*3)
   as(2)=as2(3+(it-1)*3)
   as(3)=as3(3+(it-1)*3)
   as(4)=as4(3+(it-1)*3)

   do n=1, 4
      be%vz(kts:kte,its:ite,jts:jte,n) = be%vz(kts:kte,its:ite,jts:jte,n) * as(n) / as0(n)
   enddo

! 4, Normalize the variances
 
   kz =  kde - kds + 1
!
! 4.1 Get the square-root of the variance tuning factor:

   n = 1+(it-1)*3
   if ( as1(n) <= 1.0E-7 .or. as2(n) <= 1.0E-7 .or. as3(n) <= 1.0E-7 .or. as4(n) <= 1.0E-7 .or. as5(n) <= 1.0E-7 ) then
     Write (i_char, '(i8)') it
     call da_error("da_scale_background_errors_cv3.inc",227, &
                    (/"Missing or invalid AS1 ~ AS5 value for IT = "//i_char/))
   endif

   as(1)=sqrt(as1(1+(it-1)*3))
   as(2)=sqrt(as2(1+(it-1)*3))
   as(3)=sqrt(as3(1+(it-1)*3))
   as(4)=sqrt(as4(1+(it-1)*3))
   as(5)=sqrt(as5(1+(it-1)*3))

! 4.2 Scale the horizontal scale in unit of grid-point:

   s2u= 1./grid%xb%ds
   hwll=hwll*s2u
   hwllp=hwllp*s2u
 
! 4.3 Normalize the covariance for psi, chi, t, and rh:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, n, j, i, vv, k)
   do ij = 1, grid%num_tiles
   do n=1,4

      do j=grid%j_start(ij), grid%j_end(ij)
         do i=its,ite
     ! Initialize the matrix vv:
            vv=0.
            do k=kts,kte
               vv(k,k)=1.
            enddo

     ! Apply pseudo-Gaussian "right-left" smoother in vertical with 
     ! the vertical scale length be%vz: 
            call da_rfz0(vv,kz,kz,be%ndeg,&
                         be%vz(kts:kte,i,j,n),be%be,be%table,be%nta,be%swidth)

     ! Normalize the covariance for psi, chi, t, and rh:
            do k=kts,kte
               be % corz(i,j,k,n)=be % corz(i,j,k,n)*as(n) &
                                  *samp/hwll(i,j,k,n)/vv(k,k)/global_fac(i,j)
            enddo

         enddo
      enddo
   enddo
   enddo
   !$OMP END PARALLEL DO

! 4.4 Normalize the covariance for Psfc:

!    xsum(1) = sum (be%corp(its:ite,jts:jte)*be%corp(its:ite,jts:jte))
!    call da_proc_sum_real(xsum)
!    if (rootproc) &
!    write (*,'("BEFORE: ",a,2x,"sum^2=",e20.12)') 'PSFC_u', xsum(1)

!   if (rootproc) then
!   xsum(1) = sum (hwllp*hwllp)
!   write (*,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PS_u SCALES:', xsum(1)
!   endif
!   write (*,'("as(5)=",f15.5,2x,"samp=",f15.5)') as(5), samp

    be % corp(its:ite,jts:jte)=be % corp(its:ite,jts:jte)*as(5) &
         *samp/hwllp(its:ite,jts:jte)/global_fac(its:ite,jts:jte)

!    xsum(1) = sum (be%corp(its:ite,jts:jte)*be%corp(its:ite,jts:jte))
!    call da_proc_sum_real(xsum)
!    if (rootproc) &
!    write (*,'(" AFTER: ",a,2x,"sum^2=",e20.12)') 'PSFC_u', xsum(1)

    write(mesg,*) 'Normalized covariance for Psfc: sum(be%corp*be%corp)=', &
                  sum(be%corp*be%corp)
    call wrf_debug ( 1 , mesg )
!
! 5, Assign the inverse of scale length fields for recursive filter:

   ic = 1 + ( ide- ids)/2
   jc = 1 + ( jde- jds)/2

!   if (rootproc) &
!   write (*,'("i=",i4,2x,"jc=",i4)') ic, jc
!
! 5.1 allocate the arrays for be scales: y-direction and x-direction:

    ALLOCATE ( be % sljpy ( grid%xp%ipsy: grid%xp%ipey, grid%xp%jpsy: grid%xp%jpey) )
    ALLOCATE ( be % sljy ( grid%xp%ipsy: grid%xp%ipey, grid%xp%jpsy: grid%xp%jpey, grid%xp%kpsy: grid%xp%kpey,1:4) )

    ALLOCATE ( be % slipx ( grid%xp%ipsx: grid%xp%ipex, grid%xp%jpsx: grid%xp%jpex) )
    ALLOCATE ( be % slix ( grid%xp%ipsx: grid%xp%ipex, grid%xp%jpsx: grid%xp%jpex, grid%xp%kpsx: grid%xp%kpex,1:4) )

! 5.2 Y-direction:

! 5.2.1 3-D fields: psi, chi, t, and rh:

    do n=1,4
       do k= grid%xp%kpsy, grid%xp%kpey
          do j= grid%xp%jpsy, grid%xp%jpey
             do i= grid%xp%ipsy, grid%xp%ipey
                be%sljy(i,j,k,n)=1./global_fac(i,j)/hwll(i,j,k,n)
             enddo
          enddo
       enddo
    enddo
!    write (*,'("a, be%sljy:sum^2=",e20.12)') sum(be%sljy*be%sljy)

    ! Above level ku, the sljy fields are set to a constant 
    ! for psi and chi, i.e. homogenous:
    do n=1,2
       do k=max(ku, grid%xp%kpsy), grid%xp%kpey
          slim=1./global_fac(ic,jc)/hwll(ic,jc,k,n)
          do j= grid%xp%jpsy, grid%xp%jpey
             do i= grid%xp%ipsy, grid%xp%ipey
                be%sljy(i,j,k,n)=slim
             enddo
          enddo
       enddo
    enddo
!    write (*,'("b, be%sljy:sum^2=",e20.12)') sum(be%sljy*be%sljy)

! 5.2.2 2-D field: Psfc:
    do j= grid%xp%jpsy, grid%xp%jpey
       do i= grid%xp%ipsy, grid%xp%ipey
          be%sljpy(i,j)=1./global_fac(i,j)/hwllp(i,j)
       enddo
    enddo
!    write (*,'("   be%sljpy:sum^2=",e20.12)') sum(be%sljpy*be%sljpy)

! 5.3 X-direction:

! 5.3.1 3-D fields: psi, chi, t, and rh:

    do n=1,4
      do k= grid%xp%kpsx, grid%xp%kpex
         do j= grid%xp%jpsx, grid%xp%jpex
            do i= grid%xp%ipsx, grid%xp%ipex
               be%slix(i,j,k,n)=1./global_fac(i,j)/hwll(i,j,k,n)
            enddo
         enddo
      enddo
    enddo
!    write (*,'("a, be%slix:sum^2=",e20.12)') sum(be%slix*be%slix)

   ! Above level ku, the sljy fields are set to a constant 
   ! for psi and chi, i.e. homogenous:
    do n=1,2
      do k=max(ku, grid%xp%kpsx), grid%xp%kpex
         slim=1./global_fac(ic,jc)/hwll(ic,jc,k,n)
         do j= grid%xp%jpsx, grid%xp%jpex
            do i= grid%xp%ipsx, grid%xp%ipex
               be%slix(i,j,k,n)=slim
            enddo
         enddo
      enddo
    enddo
!    write (*,'("b, be%slix:sum^2=",e20.12)') sum(be%slix*be%slix)

! 5.3.2 2-D field: Psfc:
    do j= grid%xp%jpsx, grid%xp%jpex
      do i= grid%xp%ipsx, grid%xp%ipex
         be%slipx(i,j)=1./global_fac(i,j)/hwllp(i,j)
      enddo
    enddo
!    write (*,'("   be%slipx:sum^2=",e20.12)') sum(be%slipx*be%slipx)

end subroutine da_scale_background_errors_cv3
subroutine da_setup_background_errors(grid, be)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors.
   !          Wrapper subroutine.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(in)   :: grid
!  type (xb_type), intent(in)  :: xb       ! First guess structure.
   type (be_type), intent(out) :: be       ! Back. errors structure.

   if (trace_use) call da_trace_entry("da_setup_background_errors")

!  Hybrid parameters:
   be % ne = ensdim_alpha                          ! Size of ensemble.

   if (be % ne > 0) then     ! Calculation to preserve total variance.
      if ( je_factor > alpha_std_dev**2 ) then
         jb_factor   = je_factor / ( je_factor - alpha_std_dev**2 )
      else
         jb_factor   = -999.
         write(6,*) 'Full ensemble mode: deactivating Jb control variable'
            max_vert_var1 = 0
            max_vert_var2 = 0
            max_vert_var3 = 0
            max_vert_var4 = 0
            max_vert_var5 = 0
      end if
      sigma_alpha = alpha_std_dev
      write(6,'(a,4f15.5)')' jb_factor, je_factor, alpha_std_dev, sigma_alpha = ', &
                    jb_factor, je_factor, alpha_std_dev, sigma_alpha
   else
      jb_factor = 1.0
   end if

   be % v1 % mz = 0
   be % v2 % mz = 0
   be % v3 % mz = 0
   be % v4 % mz = 0
   be % v5 % mz = 0
   if (global) then
      call da_setup_be_global(be)
   else if(fg_format == fg_format_wrf_arw_regional) then    
      if ( (cv_options == 5) .or. (cv_options == 6) .or. (cv_options == 7) ) then
         call da_setup_be_regional (grid%xb, be, grid)
      else if(cv_options == 3 ) then
         call da_setup_be_ncep_gfs (grid, be)
      else 
         write(unit=message(1),fmt='(A,I4)') &
             'Invalid CV option chosen:  cv_options = ', cv_options
         call da_error("da_setup_background_errors.inc",69,message(1:1))
      endif
   else if(fg_format == fg_format_wrf_nmm_regional ) then
!rizvi TBD            call da_setup_be_regional (grid%xb, be)
            call da_setup_be_nmm_regional (grid%xb, be)
   end if

   call da_setup_cv (be)

   if (trace_use) call da_trace_exit("da_setup_background_errors")

end subroutine da_setup_background_errors


subroutine da_setup_be_global (be)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors
   ! 
   ! Updates: 
   !
   !       Implementation of multi-variate BE  
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !---------------------------------------------------------------------------

   implicit none

   type (be_type), intent(out) :: be                 ! Back. errors structure

   integer                     :: nrec
   integer                     :: max_wave_in        ! Dimension of input power
   integer                     :: i, j, k, n ! Loop counters
   ! real, allocatable   :: height(:,:,:)      ! Height field.
   integer, allocatable:: bin(:,:,:)      ! Bin assigned to each 3D point
   integer, allocatable:: bin2d(:,:)      ! Bin assigned to each 2D point
   integer             :: bin_type        ! Type of bin to average over.
   integer             :: num_bins        ! Number of bins (3D fields).
   integer             :: num_bins2d      ! Number of bins (3D fields).
   logical             :: dummy

   real                :: binwidth_lat       ! Used if bin_type = 2 (degrees) 
   real                :: binwidth_hgt       ! Used if bin_type = 2 (m)
   real                :: hgt_min, hgt_max   ! Used if bin_type = 2 (m)
   real                :: lat_min, lat_max   ! Used if bin_type = 2 (degrees)

   character*10         :: variable
   integer              :: ni, nj, nk, b, be_unit
   integer, allocatable :: max_wave(:)
   real*8, allocatable  :: evec_g(:,:)
   real*8, allocatable  :: eval_g(:)  
   real*8, allocatable  :: evec_loc(:,:,:)
   real*8, allocatable  :: eval_loc(:,:)  
   real, allocatable   :: regcoeff_psi_chi(:)        ! psi/chi    regression cooefficient.
   real, allocatable   :: regcoeff_psi_t(:,:,:)      ! psi/t   regression cooefficient.
   real, allocatable   :: regcoeff_psi_ps(:,:)       ! psi/ps     regression cooefficient.
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     ! psi/rh     regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    ! chi_u/t regression coefficient
   real, allocatable   :: regcoeff_chi_u_ps(:,:)     ! chi_u/ps   regression coefficient
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   ! chi_u/rh   regression coefficient
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     ! t_u/rh  regression coefficient
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      ! ps_u/rh    regression coefficient

   real*8, allocatable :: power(:)               ! Temporary power spectrum.
   real, allocatable    :: power2d(:,:)           ! Temporary power spectrum.

   if (trace_use) call da_trace_entry("da_setup_be_global")

   call da_message((/"[3.0] Set up background errors (be) for global WRF-Var"/))

   be % max_wave = ide / 2 - 1

   !---------------------------------------------------------------------
   ! [1] Read in be metadata:
   !---------------------------------------------------------------------

   call da_get_unit(be_unit)
   open (unit=be_unit,file="be.dat",status="old",form="unformatted")

   read (be_unit, end= 99, err = 100) ni, nj, nk 
   read (be_unit, err = 100) bin_type
   read (be_unit, err = 100) lat_min, lat_max, binwidth_lat
   read (be_unit, err = 100) hgt_min, hgt_max, binwidth_hgt
   read (be_unit, err = 100) num_bins, num_bins2d

   allocate (bin(1:ni,1:nj,1:nk))
   allocate (bin2d(1:ni,1:nj))
   read (be_unit, err = 100) bin(1:ni,1:nj,1:nk)
   read (be_unit, err = 100) bin2d(1:ni,1:nj)

   if (ni /= ide .or. nj /= jde .or. nk /= kde) then
      call da_error("da_setup_be_global.inc",77, &
         (/"Cannot generate BE at this resolution"/))
   end if

   !---------------------------------------------------------------------
   ! [2] Read in be regression coefficients:
   !---------------------------------------------------------------------


   allocate  (regcoeff_psi_chi(1:num_bins))
   allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))

   regcoeff_psi_chi = 0.
   regcoeff_psi_t   = 0.
   regcoeff_psi_ps  = 0.
   regcoeff_psi_rh  = 0.
   regcoeff_chi_u_t = 0.
   regcoeff_chi_u_ps= 0.
   regcoeff_chi_u_rh= 0.
   regcoeff_t_u_rh  = 0.
   regcoeff_ps_u_rh = 0.

   read (be_unit,end= 99,  err = 100) regcoeff_psi_chi
   read (be_unit,end= 99,  err = 100) regcoeff_psi_ps 
   read (be_unit,end= 99,  err = 100) regcoeff_psi_t

   ! Fill regression coeff. array
   allocate (be%reg_psi_chi(1:jde,1:nk))
   allocate (be%reg_psi_t  (1:jde,1:nk,1:nk))
   allocate (be%reg_psi_ps (1:jde,1:nk))
   allocate (be%reg_psi_rh (1:jde,1:nk,1:nk))
   allocate (be%reg_chi_u_t(1:jde,1:nk,1:nk))
   allocate (be%reg_chi_u_ps(1:jde,1:nk))
   allocate (be%reg_chi_u_rh(1:jde,1:nk,1:nk))
   allocate (be%reg_t_u_rh (1:jde,1:nk,1:nk))
   allocate (be%reg_ps_u_rh(1:jde,1:nk))


   be%reg_psi_chi = 0.
   be%reg_psi_t   = 0.
   be%reg_psi_ps  = 0.
   be%reg_psi_rh  = 0.
   be%reg_chi_u_t = 0.
   be%reg_chi_u_ps= 0.
   be%reg_chi_u_rh= 0.
   be%reg_t_u_rh  = 0.


   do k=1,nk
      do j =1, jde
         b = bin(1,j,k)
         be%reg_psi_chi(j,k) = psi_chi_factor * regcoeff_psi_chi(b)
      end do
   end do

   do j=1,jde
      b = bin2d(1,j)
      do k=1,nk
         be%reg_psi_ps(j,k)   = psi_ps_factor   * regcoeff_psi_ps(k,b)
         be%reg_ps_u_rh(j,k)  = ps_u_rh_factor  * regcoeff_ps_u_rh(k,b)
         be%reg_chi_u_ps(j,k) = chi_u_ps_factor * regcoeff_chi_u_ps(k,b)
      end do
   end do

   do j=1,jde
      b = bin2d(1,j)
      do i=1,nk
         do k=1,nk
            be%reg_psi_t(j,i,k)   = psi_t_factor      *  regcoeff_psi_t(i,k,b)
            be%reg_psi_rh(j,i,k)  = psi_rh_factor     *  regcoeff_psi_rh(i,k,b)
            be%reg_chi_u_t(j,i,k) = chi_u_t_factor    *  regcoeff_chi_u_t(i,k,b)
            be%reg_chi_u_rh(j,i,k)= chi_u_rh_factor   *  regcoeff_chi_u_rh(i,k,b)
            be%reg_t_u_rh(j,i,k)  = t_u_rh_factor     *  regcoeff_t_u_rh(i,k,b)
         end do
      end do
   end do

   deallocate (regcoeff_psi_chi)
   deallocate (regcoeff_psi_t)
   deallocate (regcoeff_psi_ps)
   deallocate (regcoeff_psi_rh)
   deallocate (regcoeff_chi_u_t)
   deallocate (regcoeff_chi_u_ps)
   deallocate (regcoeff_chi_u_rh)
   deallocate (regcoeff_t_u_rh)
   deallocate (regcoeff_ps_u_rh)

   !---------------------------------------------------------------------
   ! [3] Read in be vertical eigenmodes:
   !---------------------------------------------------------------------

   do nrec = 1, 4
      read (be_unit,end= 99,  err = 100) variable   
      read (be_unit,end= 99,  err = 100) nk, num_bins2d 

      allocate (evec_g(1:nk,1:nk))
      allocate (eval_g(1:nk))
      allocate (evec_loc(1:nk,1:nk,num_bins2d))
      allocate (eval_loc(1:nk,num_bins2d))

      read (be_unit,end= 99, err = 100) evec_g     
      read (be_unit,end= 99, err = 100) eval_g     
      read (be_unit,end= 99, err = 100) evec_loc     
      read (be_unit,end= 99, err = 100) eval_loc    

      if (nrec == 1) then
         be % v1 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var1, var_scaling1(1), be%v1)

      else if (nrec == 2) then
         be % v2 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var2, var_scaling2(1), be%v2)

      else if (nrec == 3) then
         be % v3 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var3, var_scaling3(1), be%v3)

      else if (nrec == 4) then
         be % v4 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var4, var_scaling4(1), be%v4)
      end if 

      deallocate (evec_g)     
      deallocate (eval_g)     
      deallocate (evec_loc)     
      deallocate (eval_loc)     

   end do ! loop nrec

   deallocate (bin)
   deallocate (bin2d)

   !---------------------------------------------------------------------
   ! [4] Read in be power spectra:
   !---------------------------------------------------------------------

   allocate(max_wave(1:nk))

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format 
      if (k == 1) then
          allocate (power(0:max_wave_in))                      ! Temporary.
          allocate (power2d(0:max_wave_in,1:nk))               ! Temporary.
      end if
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra(be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v1 % max_wave = maxval(max_wave(1:nk))
   write (unit=stdout,fmt='(/3x,3a,i6)') &
      'Horizontal truncation for ', be % v1 % name, ' = ', be % v1 % max_wave
   allocate (be % v1 % power(0:be % v1 % max_wave,1:nk))
   be % v1 % power(0:be % v1 % max_wave,1:nk) = power2d(0:be % v1 % max_wave,1:nk)
   be % v1 % power(0,1:nk) = len_scaling1(1) * be % v1 % power(0,1:nk) 

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format    
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v2 % max_wave = maxval(max_wave(1:nk))
   write (unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v2 % name, ' = ', be % v2 % max_wave
   allocate (be % v2 % power(0:be % v2 % max_wave,1:nk))
   be % v2 % power(0:be % v2 % max_wave,1:nk) = power2d(0:be % v2 % max_wave,1:nk)
   be % v2 % power(0,1:nk) = len_scaling2(1) * be % v2 % power(0,1:nk) 

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format    
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v3 % max_wave = maxval(max_wave(1:nk))
   write(unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v3 % name, ' = ', be % v3 % max_wave
   allocate (be % v3 % power(0:be % v3 % max_wave,1:nk))
   be % v3 % power(0:be % v3 % max_wave,1:nk) = power2d(0:be % v3 % max_wave,1:nk)
   be % v3 % power(0,1:nk) = len_scaling3(1) * be % v3 % power(0,1:nk) 

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v4 % max_wave = maxval(max_wave(1:nk))
   write (unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v4 % name, ' = ', be % v4 % max_wave
   allocate (be % v4 % power(0:be % v4 % max_wave,1:nk))
   be % v4 % power(0:be % v4 % max_wave,1:nk) = power2d(0:be % v4 % max_wave,1:nk)
   be % v4 % power(0,1:nk) = len_scaling4(1) * be % v4 % power(0,1:nk) 

   ! ps_u:
   read (be_unit) variable
   be % v5 % name = variable
   be % v5 % mz = 1
   if (max_vert_var5 <=  0.0) be % v5 % mz = 0                         
   read (be_unit) max_wave_in, nrec
   read (be_unit) dummy ! use to preserve file format
   read (be_unit) power(0:max_wave_in) 

   ! Truncate power spectra:
   call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                             power, be % v5 % max_wave)

   write (unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v5 % name, ' = ', be % v5 % max_wave
   allocate (be % v5 % power(0:be % v5 % max_wave,1))
   be % v5 % power(0:be % v5 % max_wave,1) = power(0:be % v5 % max_wave)
   be % v5 % power(0,1) = len_scaling5(1) * be%v5%power(0,1) 

   deallocate(power)

   !--------------------------------------------------------------
   ! [5] Perform checks on eigenvectors:
   !--------------------------------------------------------------

   if (test_statistics) then
      call da_check_eof_decomposition(be%v1%val_g(:), be%v1%evec_g(:,:),&
                                     be%v1%name)
      call da_check_eof_decomposition(be%v2%val_g(:), be%v2%evec_g(:,:),&
                                     be%v2%name)
      call da_check_eof_decomposition(be%v3%val_g(:), be%v3%evec_g(:,:),&
                                     be%v3%name)
      call da_check_eof_decomposition(be%v4%val_g(:), be%v4%evec_g(:,:),&
                                     be%v4%name)
   end if
 
   !--------------------------------------------------------------
   ! [6] Set up alpha (flow-dependent) control variable "errors": 
   !--------------------------------------------------------------

   if (be % ne > 0) then
      be % alpha % mz = be % ne
      be % alpha % name = 'alpha'
      be % alpha % max_wave = alpha_truncation
      if (print_detail_be) then
         write(unit=stdout,fmt='(3x,3a,i6)') &
            'Horizontal truncation for ', be % alpha % name, ' = ', &
            be % alpha % max_wave
      end if
      allocate (power(0:be % alpha % max_wave))
      call da_calc_power_spectrum(be % alpha % max_wave, power)

      allocate (be % alpha % power(0:be % alpha % max_wave, be % ne))
      do n = 1, be % ne
         be % alpha % power(0:be % alpha % max_wave,n) = power(0:be % alpha % max_wave)
      end do
   end if

   deallocate (max_wave)

   write (unit=stdout,fmt='(A)') " "

   close (be_unit)
   call da_free_unit(be_unit)

   if (trace_use) call da_trace_exit("da_setup_be_global")

   return

99 write (unit=message(1),fmt='(a, i5)')' Unexpected end on BE-unit = ',be_unit
   call da_error("da_setup_be_global.inc",376,message(1:1))

100 write (unit=message(1),fmt='(a, i5)')' Read error on BE-unit = ',be_unit
   call da_error("da_setup_be_global.inc",379,message(1:1))

   ! FIX? daft having these here

   deallocate(power)
   deallocate(power2d)

end subroutine da_setup_be_global


subroutine da_setup_be_ncep_gfs( grid, be )
!------------------------------------------------------------------------------
!  PURPOSE: Define and allocate components of background errors for cv_option 3.
!
!  METHOD:  Allocate components in turn.
!
!  HISTORY: 08/02/2002 - Creation of F90 version.           Wan-Shu Wu
!
!  PARENT_MODULE: DA_Setup_Structures
!
!  Modified by Yong-Run Guo,  03/07/2009, for WRFVar 3.1
!
!------------------------------------------------------------------------------

   IMPLICIT NONE

   type (domain), intent(in)     :: grid
   TYPE (be_type), INTENT(INOUT) :: be                    ! Back. errors structure.

   INTEGER                     :: ij,ijk                ! Scalar.
   INTEGER                     :: i, j, k, ic, jc, ii   ! Loop counters.
   INTEGER                     :: ier, be_unit          ! error index


! added for AVN
   integer                     :: nlath
   integer                     :: nsig
   integer                     :: m,n,m1,n1,n4            ! loop counter
   integer                     :: msig,mlath,nmdszh,kcap  ! dummy variables
   REAL, ALLOCATABLE           :: corz_kz(:,:)
   REAL, ALLOCATABLE           :: cord_kz(:,:)
   REAL, ALLOCATABLE           :: corh_kz(:,:)
   REAL, ALLOCATABLE           :: corq_kz(:,:)
   REAL, ALLOCATABLE           :: corz_avn(:,:)
   REAL, ALLOCATABLE           :: cord_avn(:,:)
   REAL, ALLOCATABLE           :: corh_avn(:,:)
   REAL, ALLOCATABLE           :: corq_avn(:,:)
   REAL, ALLOCATABLE           :: corp_avn(:)
   REAL, ALLOCATABLE           :: clat_avn(:),sigma_avn(:)
   REAL, ALLOCATABLE           :: hwll_avn(:,:,:),hwllp_avn(:),hwll_kz(:,:,:)
   REAL, ALLOCATABLE           :: vztdq_avn(:,:,:),vztdq_kz(:,:,:)

   REAL, ALLOCATABLE           :: agv_avn(:,:,:),agv_kz(:,:,:)
   REAL, ALLOCATABLE           :: bv_avn(:,:),wgv_avn(:,:),bv_kz(:,:),wgv_kz(:,:)
   REAL, ALLOCATABLE           :: dsh(:),turn(:,:)

   REAL, DIMENSION( kts: kte,  kts: kte) :: vv
   REAL, DIMENSION( ids: ide,  jds: jde,  kds: kde, 1:4) :: hwll
   REAL, DIMENSION( ids: ide,  jds: jde)                 :: hwllp, &
                                                            coef1, coef2, &
                                                            global_lat, global_fac
   INTEGER, DIMENSION( ids: ide,  jds: jde)              :: mlat

   INTEGER :: nta,ndeg,ku,kz
   real    :: samp,s2u,tin,as(5),slim
   character(len=256) :: mesg

   integer                     :: be_rf_unit, be_print_unit, it
   real                        :: xsum
   real, dimension(1)          :: xxsum
   real, dimension(ids: ide,  jds: jde)                   :: global_2d
   real, dimension(ids: ide,  jds: jde,  kds: kde)        :: global_3d
   real, dimension(ims: ime,  jms: jme,  kms: kme)        :: corz_3d, vz_3d
   real, dimension(ims: ime,  jms: jme)                   :: corp_2d
   character(len=6)            :: vname

   write (6,'(A)') ' ----------------------------------------------------------'
   write (6,'(A,I3)') ' [3.0] Set up background errors (be) for cv_option:', cv_options
   write (6,'(A)') ' ----------------------------------------------------------'
   write (6,*)

   if (cv_options /= 3) then
      write(unit=message(1),fmt='(A)') 'Something has gone horribly wrong here'
      write(unit=message(2),fmt='(A,I4)') &
          'This subroutine is for cv_options = 3, yet cv_options = ', cv_options
      call da_error("da_setup_be_ncep_gfs.inc",76,message(1:2))
   endif

!!!!!!!!! cv_options=3
   be % v1 % name = 'psi  '           ! Streamfunction
   be % v2 % name = 'chi_u'           ! Uncorrelated velocity potential.
   be % v3 % name = 't_u'             ! Unbalanced temperature.
   be % v4 % name = 'q/qsg'
   be % v5 % name = 'psfc'            ! surface pressure
   write(6,'(3x,A)')' DA_Setup_Background_Errors: 3DVAR dry control variables are:'
   write(6,'(4x,7A)')TRIM(be % v1 % name), ', ', TRIM(be % v2 % name), ', ', &
                  TRIM(be % v3 % name), ' and ', TRIM(be % v5 % name)

   write(6,'(3x,A,A)')' DA_Setup_Background_Errors: 3DVAR humidity control variable is ',&
                     TRIM(be % v4 % name)

   write(6,*)

   be % mix =  ide -  ids + 1
   be % mjy =  jde -  jds + 1

   ij = ( ite- its+1)*( jte- jts+1)
   ijk = ij * ( kte- kts+1)

   be % v1 % mz =  kde- kds+1
   be % v2 % mz =  kde- kds+1
   be % v3 % mz =  kde- kds+1
   be % v4 % mz =  kde- kds+1
   be % v5 % mz = 1           

   be % cv % size1  = ijk
   be % cv % size2  = ijk
   be % cv % size3  = ijk
   be % cv % size4  = ijk
   be % cv % size5  = ij

   be % cv % size = be % cv % size1 + be % cv % size2 + be % cv % size3 + &
                    be % cv % size4 + be % cv % size5

   cv_size = be % cv % size


   call da_get_unit(be_unit)
   open (unit=be_unit,file="be.dat",status="old",form="unformatted")

   rewind(be_unit)

   be_rf_unit    = unit_end + 3

   read(be_unit, iostat= ier) nsig,nlath

   if (ier /= 0) then
      write(unit=message(1), fmt='(A,I4,A,I4)') &
           'cv_options:', cv_options,' Reading error in unit=',be_unit
      write(unit=message(2), fmt='(A)') 'There is likely a problem with your be.dat file'
      call da_error("da_setup_be_ncep_gfs.inc",131,message(1:2))
   endif

   write(unit=message(1),fmt='(A,I10)') 'Number of vertical level for stats = ', nsig
   write(unit=message(2),fmt='(A,I10)') 'Number of latitude           nlath = ', nlath
   call da_message(message(1:2))

   kz =  kde- kds+1

   if(nsig.ne.kz)then
      write(unit=message(1),fmt='(A,I6)') 'Number of vertical level for WRFVar=', kz
      call da_message(message(1:1))
   end if

! 1, Allocate the arrays 

!   1.1 for BK STATS at WRF eta-levels:

  ! Variances for psi(corz), chi_u(cord), t_u(corh), and p-rh(corq)  
   ALLOCATE ( corz_kz(1:2*nlath+1,1:kz),cord_kz(1:2*nlath+1,1:kz) )
   ALLOCATE ( corh_kz(1:2*nlath+1,1:kz),corq_kz(1:2*nlath+1,1:kz) )
  ! Scale lengths: horizontal (hwll), vertical (vztdq)
   ALLOCATE ( hwll_kz(0:nlath*2+1,1:kz,1:4)                      )
   ALLOCATE ( vztdq_kz(1:kz,0:nlath*2+1,1:4)                     )
  ! Regression coefficients: t(agv), chi(bv), psfc(wgv)
   ALLOCATE ( agv_kz(0:nlath*2+1,1:kz,1:kz)                      )
   ALLOCATE ( bv_kz(0:nlath*2+1,1:kz),wgv_kz(0:nlath*2+1,1:kz)   )

!   1.2 for BK STATS inputed from NCEP GFS sigma levels:

  ! Variances for psi(corz), chi_u(cord), t_u(corh), and p-rh(corq)
   ALLOCATE ( corz_avn(1:2*nlath+1,1:nsig),cord_avn(1:2*nlath+1,1:nsig) )
   ALLOCATE ( corh_avn(1:2*nlath+1,1:nsig),corq_avn(1:2*nlath+1,1:nsig) )
  ! Variance for Psfc:
   ALLOCATE ( corp_avn(1:2*nlath+1),clat_avn(1:2*nlath),sigma_avn(1:nsig) )
  ! Scale lengths: horizontal (hwll), vertical (vztdq)
   ALLOCATE ( hwll_avn(0:nlath*2+1,1:nsig,1:4),hwllp_avn(0:nlath*2+1) )
   ALLOCATE ( vztdq_avn(1:nsig,0:nlath*2+1,1:4)                     )
  ! Regression coefficients: t(agv), chi(bv), psfc(wgv)
   ALLOCATE ( agv_avn(0:nlath*2+1,1:nsig,1:nsig)                    )
   ALLOCATE ( bv_avn(0:nlath*2+1,1:nsig),wgv_avn(0:nlath*2+1,1:nsig) )

!   1.3 for BK STATS at the WRF model grids:

   ALLOCATE ( be % corz(its:ite,jts:jte,kts:kte,1:4) )
   ALLOCATE ( be % corp(its:ite,jts:jte) )

   ALLOCATE ( be % vz(kts:kte,its:ite,jts:jte,1:4) )

   ALLOCATE ( be % agvz(its:ite,jts:jte,kts:kte,kts:kte) )
   ALLOCATE ( be % bvz(its:ite,jts:jte,kts:kte) )
   ALLOCATE ( be % wgvz(its:ite,jts:jte,kts:kte) )
!

! 2, load the WFR model latitude and map factor:

   ij = ( ide- ids+1)*( jde- jds+1)
!  Collect xb component of lat into global buffer.
   call da_patch_to_global( grid, grid%xb%lat, global_lat )
   call wrf_dm_bcast_real( global_lat, ij )
!  Collect xb component of fac into global buffer.
   call da_patch_to_global( grid, grid%xb%map_factor, global_fac )
   call wrf_dm_bcast_real( global_fac, ij )

! 3, Read in the NCEP GFS BK STATS:

!   3.1 Latitude and sigma values:

   read(be_unit, iostat=ier) clat_avn,(sigma_avn(k),k=1,nsig)
   if (ier /= 0) then
      write(unit=message(1), fmt='(A,I4,A,I4)') &
           'cv_options:', cv_options,' Reading error in unit=',be_unit
      write(unit=message(2), fmt='(A)') 'There is likely a problem with your be.dat file'
      call da_error("da_setup_be_ncep_gfs.inc",213,message(1:2))
   endif

!   3.2 Variances:

   m=2*nlath+1
   read (be_unit) &
                   ((corz_avn(i,k),i=1,m),k=1,nsig),   &
                   ((cord_avn(i,k),i=1,m),k=1,nsig),   &
                   ((corh_avn(i,k),i=1,m),k=1,nsig),   &
                   ((corq_avn(i,k),i=1,m),k=1,nsig),corp_avn

!   3.3 Scale lengths

   ! horizontal
   read(be_unit) (((hwll_avn(i,k,m),i=0,nlath*2+1),k=1,nsig),m=1,4),   &
                     hwllp_avn
   ! vertical:
   read(be_unit) (((vztdq_avn(k,i,m),k=1,nsig),i=0,nlath*2+1),m=1,4)

!   3.4 Regression coefficients:
  
   read(be_unit) (((agv_avn(i,k,m),i=0,nlath*2+1),k=1,nsig),m=1,nsig), &
                         ((bv_avn(i,k),i=0,nlath*2+1),k=1,nsig), &
                         ((wgv_avn(i,k),i=0,nlath*2+1),k=1,nsig)

! 4, incorporate the tuning factors:

!   4.1 Horizontal scales:
 
   as(1)=as1(2)
   as(2)=as2(2)
   as(3)=as3(2)
   as(4)=as4(2)
   do m=1,4
      do k=1,nsig
         do i=0,nlath*2+1
            hwll_avn(i,k,m)=hwll_avn(i,k,m)*as(m)
         enddo
      enddo
   enddo
   do i=0,nlath*2+1
      hwllp_avn(i)=hwllp_avn(i)*as5(2)
   enddo

!   4.2 Vertical scales:
   as(1)=as1(3)
   as(2)=as2(3)
   as(3)=as3(3)
   as(4)=as4(3)
   do m=1,4
      do i=0,nlath*2+1
         do k=1,nsig
            vztdq_avn(k,i,m)=vztdq_avn(k,i,m)*as(m)
         enddo
       enddo
   enddo

! 5, determine the level ku, which locates above and nearest sigma=0.15:

   ku= kde+1
   k_loop: do k= kds+1, kde
      if (grid%xb%sigmah(k-1)>0.15 .and. grid%xb%sigmah(k)<=0.15) then
         ku=k
         exit k_loop
      endif
   end do k_loop
   write(mesg,'(a,i5,a)') "==> level ku =", ku ," locates above and nearest sigma=0.15"
   call wrf_debug ( 1 , mesg )

! 6, Vertical interpolation of BK STATS: 
!    to convert the NCEP GFS sigma levels to WRF model eta levels

   call da_chgvres(nlath,nsig,kz,grid%xb%sigmah,sigma_avn,&
                   corz_avn,cord_avn,corh_avn,corq_avn,hwll_avn,vztdq_avn,agv_avn,bv_avn,wgv_avn,&
                   corz_kz, cord_kz, corh_kz, corq_kz, hwll_kz, vztdq_kz, agv_kz, bv_kz, wgv_kz)

! 7, Horizontal interpolation

!   7.1 Calculate the interpolation coefficients:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( j, i, m, m1 )
   do j= jds, jde
      do i= ids,  ide

         if (global_lat(i,j).ge.clat_avn(2*nlath)) then
            ! 7.1.1 Model lat >= max AVN lat: 
            mlat(i,j)=nlath*2-1
            coef1(i,j)=0.
            coef2(i,j)=1.
         else
            ! 7.1.2 Model lat < max AVN lat: 
            do m=1,2*nlath-1
               m1=m+1
               if ((global_lat(i,j).ge.clat_avn(m)).and.  &
                  (global_lat(i,j).lt.clat_avn(m1))) then
                  mlat(i,j)=m
                  exit
               end if
            end do

            coef2(i,j)=(global_lat(i,j)-clat_avn(m))/(clat_avn(m1)-clat_avn(m))
            coef1(i,j)=1.-coef2(i,j)
         endif

      end do
   end do
   !$OMP END PARALLEL DO

!   7.2 interpolation of the covariance

   ! Psfc:
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, j, i, m, m1 )
   do ij = 1, grid%num_tiles

   do j=grid%j_start(ij), grid%j_end(ij)
      do i=its,ite
         m=mlat(i,j)
         m1=m+1
         be%corp(i,j)=corp_avn(m)*coef1(i,j)+corp_avn(m1)*coef2(i,j)
      enddo
   enddo

   enddo
   !$OMP END PARALLEL DO

   ! psi, chi, t, and rh:
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, j, i, k, m, m1 )
   do ij = 1, grid%num_tiles

   do k=kts,kte
      do j=grid%j_start(ij), grid%j_end(ij)
         do i=its,ite
            m=mlat(i,j)
            m1=m+1
            be%corz(i,j,k,1)=corz_kz(m,k)*coef1(i,j)+corz_kz(m1,k)*coef2(i,j)
            be%corz(i,j,k,2)=cord_kz(m,k)*coef1(i,j)+cord_kz(m1,k)*coef2(i,j)
            be%corz(i,j,k,3)=corh_kz(m,k)*coef1(i,j)+corh_kz(m1,k)*coef2(i,j)
            be%corz(i,j,k,4)=corq_kz(m,k)*coef1(i,j)+corq_kz(m1,k)*coef2(i,j)
         end do
      end do
   end do

   end do
   !$OMP END PARALLEL DO
  
!   7.3 interpolation of the horizontal scale lengths

   ic = 1 + ( ide- ids)/2
   jc = 1 + ( jde- jds)/2

   ! Psfc:
   do j= jds,  jde
      do i= ids,  ide
         m=mlat(i,j)
         m1=m+1
         hwllp(i,j)=hwllp_avn(m)*coef1(i,j)+hwllp_avn(m1)*coef2(i,j)
      end do
   end do
   write(mesg,'(a,2i5,a,f20.10)') 'Horizontal scale length (m) for Psfc at (',ic/2,jc/2,') hwllp=',hwllp(ic/2,jc/2)
   call wrf_debug ( 1 , mesg )

   ! psi, chi, t, and rh:
   do n4=1,4
      do k= kds,  kde
         do j= jds,  jde
            do i= ids,  ide
               m=mlat(i,j)
               m1=m+1
               hwll(i,j,k,n4)=hwll_kz(m,k,n4)*coef1(i,j)+hwll_kz(m1,k,n4)*coef2(i,j)
            end do
         end do
      end do
   end do

!   7.4 interpolation of the vertical scale lengths

   do n4=1,4
      do j=jts,jte
         do i=its,ite
            m=mlat(i,j)
            m1=m+1
            do k=kts,kte
               be%vz(k,i,j,n4)=vztdq_kz(k,m,n4)*coef1(i,j)+vztdq_kz(k,m1,n4)*coef2(i,j)
            end do
         end do
      end do
   end do
   write(mesg,'(a,3i5,a,4f20.16)') 'Vertical scale length for Psi, Chi, t, and rh at (k,i,j)=', &
                                  10,its+2,jts+5,'  vz=',(be%vz(10,its+2,jts+5,n4),n4=1,4)
   call wrf_debug ( 1 , mesg )

!   7.5 interpolation of the regression coefficients

   ! Temperature:
   do k=kts,kte
      do n=kts,kte
         do j=jts,jte
            do i=its,ite
               m=mlat(i,j)
               m1=m+1
               be%agvz(i,j,n,k)=agv_kz(m,n,k)*coef1(i,j)+agv_kz(m1,n,k)*coef2(i,j)
            end do
         end do
      end do
   end do

   ! potentail velocity:
   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            m=mlat(i,j)
            m1=m+1
            be%bvz(i,j,k)=bv_kz(m,k)*coef1(i,j)+bv_kz(m1,k)*coef2(i,j)
         end do
      end do
   end do

   ! Surface pressure:
   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            m=mlat(i,j)
            m1=m+1
            be%wgvz(i,j,k)=wgv_kz(m,k)*coef1(i,j)+wgv_kz(m1,k)*coef2(i,j)
         end do
      end do
   end do

!   7.6 Deallocate the arrays:

   ! For NCEP GFS BK STATS:
   DEALLOCATE ( corz_avn,cord_avn )
   DEALLOCATE ( corh_avn,corq_avn )
   DEALLOCATE ( corp_avn,clat_avn,sigma_avn )
   DEALLOCATE ( hwll_avn,hwllp_avn )
   DEALLOCATE ( vztdq_avn )
   DEALLOCATE ( agv_avn )
   DEALLOCATE ( bv_avn,wgv_avn )

   ! For WRF model levels:
   DEALLOCATE ( corz_kz,cord_kz )
   DEALLOCATE ( corh_kz,corq_kz )
   DEALLOCATE ( hwll_kz )
   DEALLOCATE ( vztdq_kz )
   DEALLOCATE ( agv_kz )
   DEALLOCATE ( bv_kz,wgv_kz )

!
! 8, Create the parameters for recursive filter

!  call da_prerf(grid%xb,be)
!
!   8.1 Set the constants and generate be%be, samp, be%table:

   ! Constant-1: ndeg:
   ndeg=4

   be%ndeg=ndeg
 
   ALLOCATE ( turn (1:ndeg,1:ndeg) )
   ALLOCATE ( be % be (1:ndeg) )
   ALLOCATE ( be % rate (1:ndeg) )

   CALL RFDPAR1(be%BE,be%RATE,ndeg)
   CALL RFDPAR2(be%BE,be%RATE,TURn,SAMP,ndeg)

    ! Constant-2: nta:
   nta=5600

   be%nta=nta
   allocate (dsh(1:nta)        )
   ALLOCATE ( be % table (1:nta,1:ndeg) )

    ! Constant-3: be%swidth:
   be%swidth=10.

   tin=be%swidth/dble(nta)
   do i=1,nta
      dsh(i)=dble(i-1)*tin
   enddo

   call  RFDPARV(DSH,be%RATE,be%table,nta,ndeg )

!   8.2 Deallocate the working arrays:

   deallocate (dsh )
   deallocate (turn )
   be_print_unit = unit_end + 4
   
!  8.2.1 Save the CV3 BE at model grids before tuned:

if ( max_ext_its > 1 ) then
   write (*,'(/a,i3)') 'mpp code ==> Write CV3 BE to fort.',be_rf_unit

   if (rootproc) then
   write (be_rf_unit) ids, ide, jds, jde, kds, kde, ku, samp
   write (*,'(a,7i5,f15.5)') "ids, ide, jds, jde, kds, kde, ku, samp:", &
                        ids, ide, jds, jde, kds, kde, ku, samp
   endif

   ij = ( ide - ids + 1)*( jde - jds + 1)
   ijk = ij * ( kde - kds + 1)
!  Collect the variance and write out:
   do ii = 1, 4
! tile --> patch: 
     corz_3d(its:ite,jts:jte,kts:kte) = be%corz(its:ite,jts:jte,kts:kte,ii)
! patch -->global:
     call da_patch_to_global(grid, corz_3d, global_3d)
! broadcast:
     call wrf_dm_bcast_real(global_3d, ijk)

     if (ii == 1) vname = 'PSI   '
     if (ii == 2) vname = 'CHI_u '
     if (ii == 3) vname = 'TMP_u '
     if (ii == 4) vname = 'PSD_RH'

     xxsum(1) = sum (be%corz(:,:,:,ii)*be%corz(:,:,:,ii))
     call da_proc_sum_real(xxsum)
     if (rootproc) then
       write (be_rf_unit) 'VARIANCE:', vname, ii, global_3d
!       write (*,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, vname, xxsum(1)
     endif
   enddo
! tile --> patch:
   corp_2d(its:ite,jts:jte) = be%corp(its:ite,jts:jte)
! patch --> global:
   call da_patch_to_global(grid, corp_2d, global_2d)
! broadcast:
   call wrf_dm_bcast_real(global_2d, ij)
   xxsum(1) = sum (be%corp*be%corp)
   call da_proc_sum_real(xxsum)
   if (rootproc) then
     write (be_rf_unit) 'VARIANCE:', 'PSFC_u',  1,  global_2d
!     write (*,'(9x,a,2x,"sum^2=",e20.12)') 'PSFC_u', xxsum(1)
   endif
!
   do ii = 1, 4
! tile --> patch: 
     do i = its, ite
     do j = jts, jte
     do k = kts,kte
        vz_3d(i,j,k) = be%vz(k,i,j,ii)
     enddo
     enddo
     enddo
! patch -->global:
     call da_patch_to_global(grid, vz_3d, global_3d)
! broadcast:
     call wrf_dm_bcast_real(global_3d, ijk)

     if (ii == 1) vname = 'PSI   '
     if (ii == 2) vname = 'CHI_u '
     if (ii == 3) vname = 'TMP_u '
     if (ii == 4) vname = 'PSD_RH'

     xxsum(1) = sum (be%vz(:,:,:,ii)*be%vz(:,:,:,ii))
     call da_proc_sum_real(xxsum)
     if (rootproc) then
       write (be_rf_unit) 'VZ-SCALE:', vname, ii, global_3d
!       write (*,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, vname, xxsum(1)
     endif
   enddo
!
!  Horizotal scales:
    if (rootproc) then
      write (be_rf_unit) hwll
      xsum = sum (hwll*hwll)
!      write (*,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PSI, CHI_u, T_u, RH_u SCALES:', xsum
      write (be_rf_unit) hwllp
      xsum = sum (hwllp*hwllp)
!      write (*,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PS_u SCALES:', xsum
    endif


!  8.2.2 Write out the CV3 BE information to "be_print_unit":

   if (PRINT_DETAIL_BE) then

   write (be_print_unit,'(a,7i5,f15.5)') "ids, ide, jds, jde, kds, kde, ku, samp:", &
                                          ids, ide, jds, jde, kds, kde, ku, samp

   write (be_print_unit,'(3x,a)')  'VARIANCE:'
   do ii = 1, 4
     if (ii == 1) vname = 'PSI   '
     if (ii == 2) vname = 'CHI_u '
     if (ii == 3) vname = 'TMP_u '
     if (ii == 4) vname = 'PSD_RH'
     xxsum(1) = sum (be%corz(its:ite,jts:jte,kts:kte,ii)*be%corz(its:ite,jts:jte,kts:kte,ii))
     call da_proc_sum_real(xxsum)
     if (rootproc) &
     write (be_print_unit,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, vname, xxsum(1)
   enddo

!  Pscf Variance before the normalization.
   xxsum(1) = sum (be%corp(its:ite,jts:jte)*be%corp(its:ite,jts:jte))
   call da_proc_sum_real(xxsum)
   if (rootproc) &
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'PSFC_u', xxsum(1)
!
!  Collect the vertical scales and write out:
   write (be_print_unit,'(3x,a)')  'VERTICAL SCALES:'
   do ii = 1, 4
     if (ii == 1) vname = 'PSI   '
     if (ii == 2) vname = 'CHI_u '
     if (ii == 3) vname = 'TMP_u '
     if (ii == 4) vname = 'PSD_RH'
     xxsum(1) = sum (be%vz(:,:,:,ii)*be%vz(:,:,:,ii))
     call da_proc_sum_real(xxsum)
     if (rootproc) &
     write (be_print_unit,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, vname, xxsum(1)
   enddo

!  Horizotal scales:
   if (rootproc) &
   write (be_print_unit,'(3x,a)')  'READ IN THE HORIZONTAL SCALES:'
   if (rootproc) then
   write (be_print_unit,'(3x,a)')  'HORIZONTAL SCALES:'
   xxsum(1) = sum (hwll*hwll)
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PSI, CHI_u, T_u, RH_u SCALES:', xxsum(1)
   endif

   if (rootproc) then
   xxsum(1) = sum (hwllp*hwllp)
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'HORIZONTAL PS_u SCALES:', xxsum(1)
   endif

!  Collect the regression coefficients:
   write (be_print_unit,'(3x,a)')  'REGRESSION COEFF. T, CHI, PSFC:'
   do ii = kts, kte
     xxsum(1) = sum (be%agvz(:,:,:,ii)*be%agvz(:,:,:,ii))
     call da_proc_sum_real(xxsum)
     if (rootproc) &
     write (be_print_unit,'(5x,i3,1x,a,2x,"sum^2=",e20.12)') ii, 'TMP-PSI', xxsum(1)
   enddo

!  Rg. Coeff. has already stored in "be" data structure, not need to be save.
   xxsum(1) = sum (be%bvz*be%bvz)
   call da_proc_sum_real(xxsum)
   if (rootproc) &
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'CHI-PSI', xxsum(1)

!  Rg. Coeff. has already stored in "be" data structure, not need to be save.
   xxsum(1) = sum (be%wgvz*be%wgvz)
   call da_proc_sum_real(xxsum)
   if (rootproc) &
   write (be_print_unit,'(9x,a,2x,"sum^2=",e20.12)') 'PSF-PSI', xxsum(1)

   endif
end if   ! max_ext_its > 1

! 9, Incorporate the tuning factors for covariance

! sli in scale  unit (map_factor come with ds )
!           variance* amp for 3d/2d RF

   as(1)=sqrt(as1(1))
   as(2)=sqrt(as2(1))
   as(3)=sqrt(as3(1))
   as(4)=sqrt(as4(1))
   as(5)=sqrt(as5(1))

!   9.1 Scale the horizontal scale in unit of grid-point:

   s2u= 1./grid%xb%ds
   hwll=hwll*s2u
   hwllp=hwllp*s2u
 
!   9.2 Re-scale the covariance for psi, chi, t, and rh:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, n, j, i, vv, k)
   do ij = 1, grid%num_tiles
   do n=1,4
      do j=grid%j_start(ij), grid%j_end(ij)
         do i=its,ite

            vv=0.
            do k=kts,kte
               vv(k,k)=1.
            enddo

     ! Recursive filter routie applied in vertical with 
     ! the vertical scale length be%vz: 
            call da_rfz0(vv,kz,kz,be%ndeg,&
                         be%vz(kts:kte,i,j,n),be%be,be%table,be%nta,be%swidth)

     ! Re-scale the covariance for psi, chi, t, and rh:
            do k=kts,kte
               be % corz(i,j,k,n)=be % corz(i,j,k,n)*as(n) &
                                  *samp/hwll(i,j,k,n)/vv(k,k)/global_fac(i,j)
            enddo

         enddo
      enddo
   enddo
   enddo
   !$OMP END PARALLEL DO

!   9.3 Re-scale the covariance for Psfc:

    be % corp(its:ite,jts:jte)=be % corp(its:ite,jts:jte)*as(5) &
         *samp/hwllp(its:ite,jts:jte)/global_fac(its:ite,jts:jte)

    write(mesg,*) 'Re-scaled covariance for Psfc: sum(be%corp*be%corp)=', &
                  sum(be%corp*be%corp)
    call wrf_debug ( 1 , mesg )

!
! 10, Assign the inverse of scale length fields for recursive filter:
!
!   10.1 allocate the arrays for be scales: y-direction and x-direction:

    ALLOCATE ( be % sljpy ( grid%xp%ipsy: grid%xp%ipey, grid%xp%jpsy: grid%xp%jpey) )
    ALLOCATE ( be % sljy ( grid%xp%ipsy: grid%xp%ipey, grid%xp%jpsy: grid%xp%jpey, grid%xp%kpsy: grid%xp%kpey,1:4) )

    ALLOCATE ( be % slipx ( grid%xp%ipsx: grid%xp%ipex, grid%xp%jpsx: grid%xp%jpex) )
    ALLOCATE ( be % slix ( grid%xp%ipsx: grid%xp%ipex, grid%xp%jpsx: grid%xp%jpex, grid%xp%kpsx: grid%xp%kpex,1:4) )

!   10.2 Y-direction:

    ! 3-D fields: psi, chi, t, and rh:
    do n=1,4
       do k= grid%xp%kpsy, grid%xp%kpey
          do j= grid%xp%jpsy, grid%xp%jpey
             do i= grid%xp%ipsy, grid%xp%ipey
                be%sljy(i,j,k,n)=1./global_fac(i,j)/hwll(i,j,k,n)
             enddo
          enddo
       enddo
    enddo

    ! Above level ku,the sljy fields are set to a constant 
    ! for psi and chi, i.e. homogenous:
    do n=1,2
       do k=max(ku, grid%xp%kpsy), grid%xp%kpey
          slim=1./global_fac(ic,jc)/hwll(ic,jc,k,n)
          do j= grid%xp%jpsy, grid%xp%jpey
             do i= grid%xp%ipsy, grid%xp%ipey
                be%sljy(i,j,k,n)=slim
             enddo
          enddo
       enddo
    enddo

    ! 2-D field: Psfc:
    do j= grid%xp%jpsy, grid%xp%jpey
       do i= grid%xp%ipsy, grid%xp%ipey
          be%sljpy(i,j)=1./global_fac(i,j)/hwllp(i,j)
       enddo
    enddo

!   10.3 X-direction:

   ! 3-D fields: psi, chi, t, and rh:
   do n=1,4
      do k= grid%xp%kpsx, grid%xp%kpex
         do j= grid%xp%jpsx, grid%xp%jpex
            do i= grid%xp%ipsx, grid%xp%ipex
               be%slix(i,j,k,n)=1./global_fac(i,j)/hwll(i,j,k,n)
            enddo
         enddo
      enddo
   enddo

   ! Above level ku,the sljy fields are set to a constant 
   ! for psi and chi, i.e. homogenous:
   do n=1,2
      do k=max(ku, grid%xp%kpsx), grid%xp%kpex
         slim=1./global_fac(ic,jc)/hwll(ic,jc,k,n)
         do j= grid%xp%jpsx, grid%xp%jpex
            do i= grid%xp%ipsx, grid%xp%ipex
               be%slix(i,j,k,n)=slim
            enddo
         enddo
      enddo
   enddo

   ! 2-D field: Psfc:
   do j= grid%xp%jpsx, grid%xp%jpex
      do i= grid%xp%ipsx, grid%xp%ipex
         be%slipx(i,j)=1./global_fac(i,j)/hwllp(i,j)
      enddo
   enddo

   close(be_unit)
   call da_free_unit(be_unit)

end subroutine da_setup_be_ncep_gfs

subroutine da_setup_be_regional(xb, be, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors
   !
   ! Updates:
   !
   !       Implementation of multi-variate BE
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb                    ! First guess structure.
   type (be_type), intent(inout) :: be                    ! Back. errors structure.
   type (domain),  intent(in)    :: grid

   logical                do_normalize1      ! Test do_normalize.
   integer             :: i, j, k, m         ! Loop counters.
   integer, allocatable:: bin(:,:,:)         ! Bin assigned to each 3D point
   integer, allocatable:: bin2d(:,:)         ! Bin assigned to each 2D point
   integer             :: bin_type           ! Type of bin to average over.
   integer             :: num_bins           ! Number of bins (3D fields).
   integer             :: num_bins2d         ! Number of bins (3D fields).
   real*8  :: junk                           ! For reading truncated sd, wsd
   real*8  :: lat_min, lat_max, binwidth_lat ! Used if bin_type = 2 (degrees)..
   real*8  :: hgt_min, hgt_max, binwidth_hgt ! Used if bin_type = 2 (m). .

   real*8, allocatable           :: be1_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be2_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be3_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be4_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be5_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be1_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be2_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be3_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be4_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be5_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: alpha_val(:)          ! Global Eigenvalues.

   real*8, allocatable           :: be1_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be2_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be3_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be4_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be5_evec_loc(:,:,:)   ! Local Eigenvectors.

   real*8, allocatable           :: be1_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be2_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be3_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be4_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be5_evec_glo(:,:)     ! Global Eigenvectors.

   real*8, allocatable           :: alpha_evec(:,:)       ! Global Eigenvectors.

   real*8, allocatable           :: be1_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be2_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be3_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be4_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be5_rf_lengthscale(:)
   real*8, allocatable           :: alpha_rf_lengthscale(:)
   real*8, allocatable           :: alpha_rf_scale_factor(:)

   real*8, allocatable           :: evec_loc(:,:,:)        ! Latitudinally varying eigenvectors.
   real*8, allocatable           :: eval_loc(:,:)          ! Latitudinally varying eigenvalues.

   character*5                 :: namv(5)=(/"psi  ","chi_u","t_u  ","rh   ","ps_u "/)
   character*10                :: variable
   character*80                :: var80c
   character(len=15)           :: write_fmt
   integer                     :: ni, nj, nk, nk_2d, b
   integer                     :: ix, jy, kz, kzs(5), mz, mzs(5)
   real, allocatable   :: regcoeff_psi_chi(:)        ! psi/chi    regression cooefficient.
   real, allocatable   :: regcoeff_psi_t(:,:,:)      ! psi/t      regression cooefficient.
   real, allocatable   :: regcoeff_psi_ps(:,:)       ! psi/ps     regression cooefficient.
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     ! psi/rh     regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    ! chi_u/t    regression coefficient
   real, allocatable   :: regcoeff_chi_u_ps(:,:)     ! chi_u/ps   regression coefficient
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   ! chi_u/rh   regression coefficient
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     ! t_u/rh     regression coefficient
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      ! ps_u/rh    regression coefficient
   real                :: qrain_th_low, qrain_th_high

   integer                     :: be_unit, ier, be_rf_unit, be_print_unit, it, idummy

!-----------for interpolating CV5--------------------------------------------------------------
   REAL, ALLOCATABLE           :: reg_psi_ps0(:,:), reg_psi_chi0(:,:), reg_psi_t0(:,:,:), &
                                  reg_psi_ps   (:), reg_psi_chi   (:), reg_psi_t   (:,:), &
                                  reg_psi_ps_be(:), reg_psi_chi_be(:), reg_psi_t_be(:,:)
   REAL, DIMENSION(:,:), ALLOCATABLE:: covm1_be, covm2_be, covm3_be, covm4_be,&
                                       covm1   , covm2   , covm3   , covm4
   REAL, DIMENSION(:),   ALLOCATABLE:: rfls1_be, rfls2_be, rfls3_be, rfls4_be
   REAL, DIMENSION(:),   ALLOCATABLE:: eta_be
! For write out the Interpolated CV5 BE (bin_type=5) file:
   integer              :: bin_type_out, num_bins_out, num_bins2d_out
   integer, allocatable :: bin_out(:,:,:), bin2d_out(:,:)
   integer, parameter  :: be_out_unit   = 135

!===========for interpolating CV5==============================================================

   if (trace_use) call da_trace_entry("da_setup_be_regional")
   if( (cv_options == 5) .or. (cv_options == 6) .or. (cv_options == 7) ) then
      write (unit=message(1),fmt='(A,I3)') &
         'Set up background errors for regional application for cv_options = ',cv_options
      write (unit=message(2),fmt='(a)') ' '
      call da_message(message(1:2))
      if ( (cv_options == 7 .or. cv_options == 6) .and. interpolate_stats ) then
         write (unit=message(1),fmt='(A,I3)') 'interpolate_stats is not implemented for cv_options=',cv_options
         call da_error("da_setup_be_regional.inc",151,message(1:1))
      end if
   else 
      write (unit=message(1),fmt='(a)') 'cv_options should be 5, 6, or 7'
      call da_error("da_setup_be_regional.inc",155,message(1:1))
   end if

   ix = xb % mix
   jy = xb % mjy
   kz = xb % mkz

   be_rf_unit    = unit_end + 1
   be_print_unit = unit_end + 2

   call da_get_unit(be_unit)
   open(unit=be_unit,file="be.dat", status="old",form="unformatted")

   rewind (be_unit)
   read (be_unit, iostat=ier) ni, nj, nk
   if (ier /= 0) then
      write (unit=message(1),fmt='(a,i3)') 'Error in reading be.dat, unit= ',be_unit
      call da_error("da_setup_be_regional.inc",172,message(1:1))
   end if
   read (be_unit) bin_type

!-----------for interpolating CV5--------------------------------------------------------------
   if ( .not. interpolate_stats ) then 
      if (ni /= ix .or. nj /= jy .or. nk /= kz) then
         write(unit=message(1),fmt='(a)')  &
            'Dimensions of the assimilation domain do not match those in the BE file.'
         write(unit=message(2),fmt='(3x,a,3i4)') "in be.dat: ni, nj, nk = ",ni, nj, nk
         write(unit=message(3),fmt='(3x,a,3i4)') "in fg:     ix, jy, kz = ",ix, jy, kz
         write(message(4),'(a)') 'be.dat need to be generated using data consistent with the assimilation domain.'
         call da_error("da_setup_be_regional.inc",184,message(1:4))
      end if
   else  ! when interpolate_stats = true
      ! Must use the domain averaged Reg. Coeff.:
      lat_stats_option = .FALSE.
      ! Must use the global Eigenvector/eigenvalue
      vert_evalue = 1
      write (unit=message(1),fmt='(A,I3,A)') 'interpolate_stats is .true., BE for cv_options = ',cv_options,' will be interpolated.'
      call da_message(message(1:1))
      if (rootproc) write(be_out_unit) ix, jy, kz
      write(6,'(a/("k=",i3,2x,f10.4))') 'xb%sigmah=',(k,xb%sigmah(k),k=1,kz)
      allocate (eta_be(1:nk))
      read(be_unit,iostat=ier) eta_be
      !  If not available from BE file, get it from namelist
      if (ier /= 0) then
         write(message(1),'("ier=",i5,2x,a)') ier, " NO Eta values existed in BE file."
         write(message(2),'(a)') "Will use BE Eta values from the namelist be_eta."
         call da_message(message(1:2))
         rewind(be_unit)
         read(be_unit) idummy, idummy, idummy
         read(be_unit) idummy
         eta_be = be_eta(1:nk)
      endif
      write(6,'(a/("k=",i3,2x,f10.4))') 'eta_be=',(k,eta_be(k),k=1,nk)
      !  If eta values not defined correctly, stop
      if (eta_be(1) == eta_be(nk)) then
         write(message(1),'(a)') 'Wrong eta values found.'
         write(message(2),'(a)') 'The eta levels used in generating be.dat should be specified in namelist &wrfvar10 be_eta.'
         call da_error("da_setup_be_regional.inc",212,message(1:2))
      endif
   end if

!===========for interpolating CV5==============================================================

   allocate (bin(1:ni,1:nj,1:nk))
   allocate (bin2d(1:ni,1:nj))

   read (be_unit)lat_min, lat_max, binwidth_lat
   read (be_unit)hgt_min, hgt_max, binwidth_hgt
   read (be_unit)num_bins, num_bins2d
   read (be_unit)bin(1:ni,1:nj,1:nk)
   read (be_unit)bin2d(1:ni,1:nj)
!
!-----------for interpolating CV5--------------------------------------------------------------
! Write out the interpolated CV5 BE:
    if (interpolate_stats) then
 
! No matter which bin_type BE inputed, after Interpolated, the BE always
! became bin_type = 1 BE, i.e. the global BE.
     bin_type_out   = 5
     num_bins_out   = kz
     num_bins2d_out = 1
     if (rootproc) then 
     write(be_out_unit)bin_type_out
     write(be_out_unit)lat_min, lat_max, binwidth_lat
     write(be_out_unit)hgt_min, hgt_max, binwidth_hgt
     write(be_out_unit)num_bins_out, num_bins2d_out
     endif
      allocate (bin_out(1:ix,1:jy,1:kz))
      allocate (bin2d_out(1:ix,1:jy))
        do k = 1, kz
           bin_out (:,:,k) = k
        enddo
        bin2d_out = 1
      if (rootproc) then
      write(be_out_unit)bin_out
      write(be_out_unit)bin2d_out
      endif
      deallocate (bin_out)
      deallocate (bin2d_out)
    endif
!===========for interpolating CV5==============================================================
!
   if ( cv_options /= 7 ) then   ! No regression coefficients for cv_options == 7

      ! 1.1 Read in regression coefficients
      allocate  (regcoeff_psi_chi(1:num_bins))
      allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
      if ( cv_options == 6 ) then
         allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
         allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
         allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
         allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
         allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
         allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))
      end if

      if ( cv_options == 5 ) then
         read (be_unit) regcoeff_psi_chi
         read (be_unit) regcoeff_psi_ps 
         read (be_unit) regcoeff_psi_t
      else
         do i = 1 , 9
            read (be_unit) var80c  
            select case( trim(adjustl(var80c)) )
            case ('regcoeff_psi_chi')
               read (be_unit) regcoeff_psi_chi
            case ('regcoeff_psi_t')
               read (be_unit) regcoeff_psi_t
            case ('regcoeff_psi_ps')
               read (be_unit) regcoeff_psi_ps
            case ('regcoeff_psi_rh')
               read (be_unit) regcoeff_psi_rh
            case ('regcoeff_chi_u_t')
               read (be_unit) regcoeff_chi_u_t
            case ('regcoeff_chi_u_ps')
               read (be_unit) regcoeff_chi_u_ps
            case ('regcoeff_chi_u_rh')
               read (be_unit) regcoeff_chi_u_rh
            case ('regcoeff_t_u_rh')
               read (be_unit) regcoeff_t_u_rh
            case ('regcoeff_ps_u_rh')
               read (be_unit) regcoeff_ps_u_rh
            case default;
               message(1)=' Read problem in regression coefficients in BE file '
               write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
               call da_error("da_setup_be_regional.inc",323,message(1:2))
            end select
         end do
      end if
      ! 1.2 Fill regression coeff. array for model BE, vertical dimension is kz:

      allocate (be%reg_psi_chi  (1:jy,1:kz))
      allocate (be%reg_psi_t    (1:jy,1:kz,1:kz))
      allocate (be%reg_psi_ps   (1:jy,1:kz))
      if ( cv_options == 6 ) then
         allocate (be%reg_psi_rh   (1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_t  (1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_ps (1:jy,1:kz))
         allocate (be%reg_chi_u_rh (1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_rh   (1:jy,1:kz,1:kz))
         allocate (be%reg_ps_u_rh  (1:jy,1:kz))
      end if

      !-----------for interpolating CV5--------------------------------------------------------------
      if ( interpolate_stats ) then
         allocate (reg_psi_chi0(1:nj,1:nk))
         allocate (reg_psi_ps0 (1:nj,1:nk))
         allocate (reg_psi_t0  (1:nj,1:nk,1:nk))
      end if

      be%reg_psi_chi    = 0.
      be%reg_psi_t      = 0.
      be%reg_psi_ps     = 0.
      if ( cv_options == 6 ) then
         be%reg_psi_rh     = 0.
         be%reg_chi_u_t    = 0.
         be%reg_chi_u_ps   = 0.
         be%reg_chi_u_rh   = 0.
         be%reg_t_u_rh     = 0.
      end if

      if ( interpolate_stats ) then
         do k=1,nk
            do j=1,nj
               reg_psi_chi0(j,k) = regcoeff_psi_chi(bin(1,j,k))
               b = bin2d(1,j)
               reg_psi_ps0(j,k)  = regcoeff_psi_ps(k,b)
            end do
         end do
         do k=1,nk
            do i=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  reg_psi_t0(j,i,k) = regcoeff_psi_t(i,k,b)
               end do
            end do
         end do
      end if   ! if interpolate_stats

      if ( .not. interpolate_stats ) then
         do k=1,nk
            do j=1,nj
               be%reg_psi_chi(j,k) = psi_chi_factor * regcoeff_psi_chi(bin(1,j,k))
               b = bin2d(1,j)
               be%reg_psi_ps(j,k)  = psi_ps_factor  * regcoeff_psi_ps(k,b)
            end do
         end do
         do k=1,nk
            do i=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  be%reg_psi_t(j,i,k) = psi_t_factor * regcoeff_psi_t(i,k,b)
               end do
            end do
         end do

         if ( cv_options == 6 ) then
            do k=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  be%reg_ps_u_rh(j,k)  = ps_u_rh_factor  * regcoeff_ps_u_rh(k,b)
                  be%reg_chi_u_ps(j,k) = chi_u_ps_factor * regcoeff_chi_u_ps(k,b)
               end do
            end do
            do k=1,nk
               do i=1,nk
                  do j=1,nj
                     b = bin2d(1,j)
                     be%reg_psi_rh(j,i,k)   = psi_rh_factor   * regcoeff_psi_rh(i,k,b)
                     be%reg_chi_u_t(j,i,k)  = chi_u_t_factor  * regcoeff_chi_u_t(i,k,b)
                     be%reg_chi_u_rh(j,i,k) = chi_u_rh_factor * regcoeff_chi_u_rh(i,k,b)
                     be%reg_t_u_rh(j,i,k)   = t_u_rh_factor   * regcoeff_t_u_rh(i,k,b)
                  end do
               end do
            end do
         end if   ! if cv_options 6
      end if   ! if not interpolate_stats

      deallocate (regcoeff_psi_chi)
      deallocate (regcoeff_psi_t)
      deallocate (regcoeff_psi_ps)
      if ( cv_options == 6 ) then
         deallocate (regcoeff_psi_rh)
         deallocate (regcoeff_chi_u_t)
         deallocate (regcoeff_chi_u_ps)
         deallocate (regcoeff_chi_u_rh)
         deallocate (regcoeff_t_u_rh)
         deallocate (regcoeff_ps_u_rh)
      end if
      ! 1.3 Domain_averaged regression coefficients

      if (.not.lat_stats_option) then
         write (unit=message(1), fmt='(a)') ' '
         write (unit=message(2), fmt='(3x, a)') &
            'Using the averaged regression coefficients for unbalanced part'

         !-----------for interpolating CV5--------------------------------------------------------------
         if ( interpolate_stats ) then
            allocate (reg_psi_ps_be (1:nk))
            allocate (reg_psi_chi_be(1:nk))
            allocate (reg_psi_t_be  (1:nk,1:nk))
            do k=1,nk
               reg_psi_ps_be(k)  = sum(reg_psi_ps0 (:,k))/dble(nj)
               reg_psi_chi_be(k) = sum(reg_psi_chi0(:,k))/dble(nj)
            end do
            do m=1,nk
               do k=1,nk
                  reg_psi_t_be(k,m)=sum(reg_psi_t0(:,k,m))/dble(nj)
               end do
            end do
            deallocate (reg_psi_chi0)
            deallocate (reg_psi_ps0 )
            deallocate (reg_psi_t0  )
         end if   ! if interpolate_stats
         !===========for interpolating CV5==============================================================

         if ( (.not. interpolate_stats) .and. (bin_type /= 5) ) then
            do k=1,kz
                be%reg_psi_ps  (:,k) = sum(be%reg_psi_ps (:,k))/dble(jy)
                be%reg_psi_chi (:,k) = sum(be%reg_psi_chi(:,k))/dble(jy)
            end do
            do m=1,kz
               do k=1,kz
                   be%reg_psi_t (:,k,m) = sum(be%reg_psi_t(:,k,m))/dble(jy)
                end do
            end do
            if ( cv_options == 6 ) then
               do k=1,kz
                  be%reg_ps_u_rh (:,k) = sum(be%reg_ps_u_rh (:,k))/dble(jy)
                  be%reg_chi_u_ps(:,k) = sum(be%reg_chi_u_ps(:,k))/dble(jy)
               end do
               do m=1,kz
                  do k=1,kz
                     be%reg_psi_rh  (:,k,m) = sum(be%reg_psi_rh  (:,k,m))/dble(jy)
                     be%reg_chi_u_t (:,k,m) = sum(be%reg_chi_u_t (:,k,m))/dble(jy)
                     be%reg_chi_u_rh(:,k,m) = sum(be%reg_chi_u_rh(:,k,m))/dble(jy)
                     be%reg_t_u_rh  (:,k,m) = sum(be%reg_t_u_rh  (:,k,m))/dble(jy)
                  end do
               end do
            end if   ! if cv_options 6
         end if   ! if not interpolate_stats

      else
         write (unit=message(1), fmt='(a)') ' '
         write (unit=message(2), fmt='(3x, a)') &
            'Using the geographically-dependent regression coefficients for unbalanced part'
      end if

      call da_message(message(1:2))

   end if ! cv_options /= 7

   ! 2.0 Load the eigenvector and eigenvalue

   allocate (be1_eval_loc (1:nj,1:nk))
   allocate (be2_eval_loc (1:nj,1:nk))
   allocate (be3_eval_loc (1:nj,1:nk))
   allocate (be4_eval_loc (1:nj,1:nk))
   allocate (be5_eval_loc (1:nj,1:1))

   if (vert_corr == vert_corr_2) then

      allocate (be1_eval_glo(1:nk))
      allocate (be2_eval_glo(1:nk))
      allocate (be3_eval_glo(1:nk))
      allocate (be4_eval_glo(1:nk))
      allocate (be5_eval_glo(1:1))
      allocate (be1_evec_loc(1:nj,1:nk,1:nk))
      allocate (be2_evec_loc(1:nj,1:nk,1:nk))
      allocate (be3_evec_loc(1:nj,1:nk,1:nk))
      allocate (be4_evec_loc(1:nj,1:nk,1:nk))
      allocate (be5_evec_loc(1:nj,1: 1,1: 1))
      allocate (be1_evec_glo(1:nk,1:nk))
      allocate (be2_evec_glo(1:nk,1:nk))
      allocate (be3_evec_glo(1:nk,1:nk))
      allocate (be4_evec_glo(1:nk,1:nk))
      allocate (be5_evec_glo(1:1,1:1))
   end if
   ! 2.2 Read in the eigenvector and eigenvalue 

   do i = 1 , 4
      read (be_unit) variable
      select case( trim(adjustl(variable)) )
      case ('psi', 'u')
         be % v1 % name = trim(adjustl(variable))
         read (be_unit) nk, num_bins2d
         read (be_unit)  be1_evec_glo
         read (be_unit)  be1_eval_glo
         if( i == 1) then
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
         end if
         read (be_unit)  evec_loc
         read (be_unit)  eval_loc
         do j=1,nj
            b = bin2d(1,j)
            be1_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
            be1_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
         end do

      case ('chi_u', 'v')
         be % v2 % name = trim(adjustl(variable))
         read (be_unit) nk, num_bins2d
         read (be_unit)  be2_evec_glo
         read (be_unit)  be2_eval_glo
         read (be_unit)  evec_loc
         read (be_unit)  eval_loc
         do j=1,nj
            b = bin2d(1,j)
            be2_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
            be2_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
         end do

      case ('t_u', 't')
         be % v3 % name = trim(adjustl(variable))
         read (be_unit) nk, num_bins2d
         read (be_unit)  be3_evec_glo
         read (be_unit)  be3_eval_glo
         read (be_unit)  evec_loc
         read (be_unit)  eval_loc
         do j=1,nj
            b = bin2d(1,j)
            be3_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
            be3_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
         end do

      case ('rh_u' , 'rh' )
         be % v4 % name = trim(adjustl(variable))
         read (be_unit) nk, num_bins2d
         read (be_unit)  be4_evec_glo
         read (be_unit)  be4_eval_glo
         read (be_unit)  evec_loc
         read (be_unit)  eval_loc
         do j=1,nj
            b = bin2d(1,j)
            be4_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
            be4_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
         end do


      case default;
         message(1)=' Read problem in eigen vectors/values in BE file '
         write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
         write (unit=message(3),fmt='(A)') ' Make sure you are using the correct be.dat file for your cv_options setting!'
         call da_error("da_setup_be_regional.inc",742,message(1:3))
      end select
   end do

   deallocate (evec_loc)
   deallocate (eval_loc)
   ! 2.2.5 Control variable ps_u
   read (be_unit) variable
   read (be_unit) nk_2d, num_bins2d
   allocate (evec_loc(1:nk_2d,1:nk_2d,1:num_bins2d))
   allocate (eval_loc(1:nk_2d,        1:num_bins2d))
   read (be_unit)  be5_evec_glo
   read (be_unit)  be5_eval_glo
   read (be_unit)  evec_loc
   read (be_unit)  eval_loc
   if( (trim(adjustl(variable)) /= 'ps_u') .and. (trim(adjustl(variable)) /= 'ps') ) then
      message(1)=' Variable mismatch while transfering eigen values from BE file '
      if (cv_options == 7) then
         write (unit=message(2),fmt='(A,A)') ' Expected ps but got ',trim(adjustl(variable))
      else
         write (unit=message(2),fmt='(A,A)') ' Expected ps_u but got ',trim(adjustl(variable))
      endif
      call da_error("da_setup_be_regional.inc",767,message(1:2))
   end if
   be % v5 % name = trim(adjustl(variable))
   do j=1,nj
      b = bin2d(1,j)
      be5_evec_loc(j,1:1,1:1) = evec_loc(1:1,1:1,b)
      be5_eval_loc(j,1:1 ) = eval_loc(1:1,b)
   end do

   deallocate (evec_loc)
   deallocate (eval_loc)

!
   if(use_radarobs .and. use_radar_rf .or. use_rad .and. crtm_cloud) then  
      if ( cloud_cv_options == 1 ) be % v4 % name = 'qt   '
   end if

   write (unit=message(1),fmt='(3x,A,7A)') 'WRFDA dry control variables are: ', &
      trim(be % v1 % name), ', ', trim(be % v2 % name), ', ', &
      trim(be % v3 % name), ' and ', trim(be % v5 % name)
   write (unit=message(2),fmt='(3x,A,A)') &
      'Humidity control variable is ', trim(be % v4 % name)

   call da_message(message(1:2))
!
  if (use_rf) then

! 3.0 Load the scale lengths
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! 3.1 allocate the array for scale lengths
   allocate (rfls1_be(1:nk))
   allocate (rfls2_be(1:nk))
   allocate (rfls3_be(1:nk))
   allocate (rfls4_be(1:nk))
   allocate (be1_rf_lengthscale(1:kz))
   allocate (be2_rf_lengthscale(1:kz))
   allocate (be3_rf_lengthscale(1:kz))
   allocate (be4_rf_lengthscale(1:kz))
   allocate (be5_rf_lengthscale(1:1)) ! YRG, 01/09/2012

! ===================================================================
! This is used for do_normalize = .T., I am not sure if it is working 
! properly (???), more testing needed. (YRG, 01/12/2012)
!
      b=0			! pointer to variable*mode
      allocate(nij(0:0,0:1,0:0))

      kzs=(/kz,kz,kz,kz,1/)
      be%v1%mz=1; be%v2%mz=1;be%v3%mz=1;be%v4%mz=1;be%v5%mz=1; 
      mzs=(/be%v1%mz,be%v2%mz,be%v3%mz,be%v4%mz,be%v5%mz/)
! ====================================================================
   ! 3.2 read in the scale lengths
   do i = 1 , 5              ! variables loop:
      read (be_unit, IOSTAT=k) variable
      if (print_detail_be)&
           write(unit=stdout,fmt='(2a)')' Reading lengthscale for: ',trim(adjustl(variable))
      if( k/=0 ) then
         write(unit=message(1),fmt='(A,I0,A,I0,A,A)') 'Error ',k,' reading namv(',i,')=',namv(i)
         call da_warning("da_setup_be_regional.inc",872,message(1:1))
      end if
      select case( trim(adjustl(variable)) )
         case ('psi', 'u')
            read(be_unit) rfls1_be
            if (.not.interpolate_stats) be1_rf_lengthscale = rfls1_be
         case ('chi_u', 'v')
            read(be_unit) rfls2_be
            if (.not.interpolate_stats) be2_rf_lengthscale = rfls2_be
         case ('t_u', 't')
            read(be_unit) rfls3_be
            if (.not.interpolate_stats) be3_rf_lengthscale = rfls3_be
         case ('rh_u' , 'rh')
            read(be_unit) rfls4_be
            if (.not.interpolate_stats) be4_rf_lengthscale = rfls4_be

         case ('ps_u', 'ps')
            read(be_unit) be5_rf_lengthscale(1:1)
         case default;
            message(1)='Read problem in lengthscales in be.dat'
            write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
            call da_error("da_setup_be_regional.inc",947,message(1:2))
      endselect

         if( do_normalize )then
            read(be_unit)do_normalize1
            if( do_normalize.neqv.do_normalize1 ) &
               call da_error("da_setup_be_regional.inc",953,(/namv(i)//": do_normalize.neqv.do_normalize1"/))
            read(be_unit)nij
            if( i==1 )allocate(be%sd(nij(0,1,0),nij(0,0,0),sum(mzs)))
            do k=1,mzs(i)
               read(be_unit)be%sd(:,:,b+k)
            enddo
            do k=mzs(i)+1,kzs(i)! read and discard truncated modes:
               read(be_unit)(junk,j=1,nij(0,1,0)*nij(0,0,0))
            enddo
            write(*,'(A,": |sd[",A,"]|=",es9.3)')"da_setup_be_regional.inc",namv(i),sqrt(sum(be%sd(:,:,b+1:b+mzs(i))**2))
            b=b+mzs(i)		! point to next variable.
         endif
      enddo
      if (print_detail_be) then
! Lengthscale variables (except 5!) are arrays of length "kz", so need to account for this in fmt statement
         write(write_fmt,'(A,I0,A)')'(A,',kz,'F10.5)'
         write(unit=message(1),fmt='(A)') 'BE Length scales:'
         write(unit=message(2),fmt=write_fmt) 'be1= ',be1_rf_lengthscale
         write(unit=message(3),fmt=write_fmt) 'be2= ',be2_rf_lengthscale
         write(unit=message(4),fmt=write_fmt) 'be3= ',be3_rf_lengthscale
         write(unit=message(5),fmt=write_fmt) 'be4= ',be4_rf_lengthscale
         write(unit=message(6),fmt='(A,F10.5)') 'be5= ',be5_rf_lengthscale
         call da_message(message(1:6))
      end if

  endif
!================== for interpolating CV5 =======================================
   ! 3.3 Interpolate the CV5 BE to the model resolution

   if (interpolate_stats) then
 
   ! 3.3.1 Compute the Vertical covariance matrix arrays for input CV5 BE:
 
   ! 3.3.1.1 Allocate the covariance arrays for input CV5 BE
 
      allocate (covm1_be(1:nk,1:nk))
      allocate (covm2_be(1:nk,1:nk))
      allocate (covm3_be(1:nk,1:nk))
      allocate (covm4_be(1:nk,1:nk))
   ! 3.3.1.2 Compute the vertical covariance matrix from Input CV5 BE
   !       eigenvector/eigenvalue:
      call da_Eigen_to_covmatrix(nk, be1_evec_glo, be1_eval_glo, covm1_be)
      call da_Eigen_to_covmatrix(nk, be2_evec_glo, be2_eval_glo, covm2_be)
      call da_Eigen_to_covmatrix(nk, be3_evec_glo, be3_eval_glo, covm3_be)
      call da_Eigen_to_covmatrix(nk, be4_evec_glo, be4_eval_glo, covm4_be)
   ! Deallocate the input CV5 BE Global eigenvector/eigenvalue:
      DEALLOCATE ( be1_eval_glo )
      DEALLOCATE ( be2_eval_glo )
      DEALLOCATE ( be3_eval_glo )
      DEALLOCATE ( be4_eval_glo )
 
      DEALLOCATE ( be1_evec_glo )
      DEALLOCATE ( be2_evec_glo )
      DEALLOCATE ( be3_evec_glo )
      DEALLOCATE ( be4_evec_glo )

   ! Allocate the model CV5 BE Global eigenvector/eigenvalue:
      ALLOCATE ( be1_eval_glo(1:kz) )
      ALLOCATE ( be2_eval_glo(1:kz) )
      ALLOCATE ( be3_eval_glo(1:kz) )
      ALLOCATE ( be4_eval_glo(1:kz) )
 
      ALLOCATE ( be1_evec_glo(1:kz,1:kz) )
      ALLOCATE ( be2_evec_glo(1:kz,1:kz) )
      ALLOCATE ( be3_evec_glo(1:kz,1:kz) )
      ALLOCATE ( be4_evec_glo(1:kz,1:kz) )
 
   ! 3.3.2 Convert the Vertical resolution from input CV5 BE to Model:

   ! 3.3.2.1 Allocate the Covariance matrix arrays for model
 
      allocate (covm1(1:kz,1:kz))
      allocate (covm2(1:kz,1:kz))
      allocate (covm3(1:kz,1:kz))
      allocate (covm4(1:kz,1:kz))
 
   ! 3.3.2.2 Allocate the regression coeff. arrays for model BE
 
      allocate (reg_psi_ps (1:kz))
      allocate (reg_psi_chi(1:kz))
      allocate (reg_psi_t  (1:kz,1:kz))
   ! 3.3.2 Vertical resolution conversion from input BE(nk) to model(kz) BE:
      call da_chg_be_vres(kz, nk, xb%sigmah, eta_be,&
                        reg_psi_chi_be,reg_psi_ps_be,reg_psi_t_be, &
                        reg_psi_chi   ,reg_psi_ps   ,reg_psi_t   , &
                        covm1_be, covm2_be, covm3_be, covm4_be, &
                        covm1   , covm2   , covm3   , covm4   , &
                        rfls1_be, rfls2_be, rfls3_be, rfls4_be, &
                        be1_rf_lengthscale, be2_rf_lengthscale, &
                        be3_rf_lengthscale, be4_rf_lengthscale)
       
   ! 3.3.3 Deallocate the arrays for input CV5 BE:
      deallocate (reg_psi_ps_be )
      deallocate (reg_psi_chi_be)
      deallocate (reg_psi_t_be  )
 
      deallocate ( rfls1_be )
      deallocate ( rfls2_be )
      deallocate ( rfls3_be )
      deallocate ( rfls4_be )
   ! 3.3.4 Assign the Regression coefficients for model BE:
 
   ! 3.3.4.1 Ps and Chi regression coefficients
      do k=1,kz
        do j=1,jy
          be%reg_psi_ps(j,k)  = reg_psi_ps(k)
          be%reg_psi_chi (j,k)= reg_psi_chi(k)
        enddo
      enddo
!          write(6,*)'be%reg_psi_ps(10,10)=',be%reg_psi_ps(10,10)
!          write(6,*)'be%reg_psi_chi(10,10)=',be%reg_psi_chi(10,10)
   ! 3.3.4.2 Temp. regression coefficients
      do m=1,kz
      do k=1,kz
        do j=1,jy
          be%reg_psi_t(j,k,m) = reg_psi_t(k,m)
        enddo
      enddo
      enddo
!          write(6,*)'be%reg_psi_t(10,10,10)=',be%reg_psi_t(10,10,10)
   ! 3.3.4.3 RWite out the Interpolated BE: Regression coefficients
        if (rootproc) then 
        write(be_out_unit) reg_psi_chi
        write(be_out_unit) reg_psi_ps
        write(be_out_unit) reg_psi_t
        endif
   ! 3.3.4.4 Deallocate the domain-averaged Reg. Coeff. arrays:
 
      deallocate (reg_psi_ps )
      deallocate (reg_psi_chi)
      deallocate (reg_psi_t  )
   ! 3.3.5 Re-compute the Eigenvector/eigenvalue from the model covariance matrix
      call da_gen_eigen(kz, covm1, be1_evec_glo, be1_eval_glo)
      call da_gen_eigen(kz, covm2, be2_evec_glo, be2_eval_glo)
      call da_gen_eigen(kz, covm3, be3_evec_glo, be3_eval_glo)
      call da_gen_eigen(kz, covm4, be4_evec_glo, be4_eval_glo)
   ! 3.3.6 Write out the interpolated CV5 BE
 
   ! 3.3.6.1 Eigenvector/eigenvalue:
      variable = be % v1 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be1_evec_glo
         write(be_out_unit) be1_eval_glo
      ! - Because num_bins2d = 1, be1_evec_loc = be1_evec_glo:
         write(be_out_unit) be1_evec_glo
         write(be_out_unit) be1_eval_glo
      endif
      variable = be % v2 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be2_evec_glo
         write(be_out_unit) be2_eval_glo
      ! - Because num_bins2d = 1, be2_evec_loc = be2_evec_glo:
         write(be_out_unit) be2_evec_glo
         write(be_out_unit) be2_eval_glo
      endif
      variable = be % v3 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be3_evec_glo
         write(be_out_unit) be3_eval_glo
      ! - Because num_bins2d = 1, be3_evec_loc = be3_evec_glo:
         write(be_out_unit) be3_evec_glo
         write(be_out_unit) be3_eval_glo
      endif
      variable = be % v4 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be4_evec_glo
         write(be_out_unit) be4_eval_glo
      ! - Because num_bins2d = 1, be4_evec_loc = be4_evec_glo:
         write(be_out_unit) be4_evec_glo
         write(be_out_unit) be4_eval_glo
      endif
      variable = be % v5 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) nk_2d, num_bins2d_out
         write(be_out_unit) be5_evec_glo
         write(be_out_unit) be5_eval_glo
      ! - Because num_bins2d = 1, be5_evec_loc = be5_evec_glo:
         write(be_out_unit) be5_evec_glo
         write(be_out_unit) be5_eval_glo
      endif

   ! 3.3.6.2 Scale-length
      variable = be % v1 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be1_rf_lengthscale
      endif
      variable = be % v2 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be2_rf_lengthscale
      endif
      variable = be % v3 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be3_rf_lengthscale
      endif
      variable = be % v4 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be4_rf_lengthscale
      endif
      variable = be % v5 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be5_rf_lengthscale
      endif 
   ! - Close the output unit:
 
      close(be_out_unit)
   ! 3.3.7 Deallocate the arrays for reading the input CV5 BE:
 
      DEALLOCATE ( be1_eval_loc )
      DEALLOCATE ( be2_eval_loc )
      DEALLOCATE ( be3_eval_loc )
      DEALLOCATE ( be4_eval_loc )
      DEALLOCATE ( be5_eval_loc )
      DEALLOCATE ( be1_evec_loc )
      DEALLOCATE ( be2_evec_loc )
      DEALLOCATE ( be3_evec_loc )
      DEALLOCATE ( be4_evec_loc )
      DEALLOCATE ( be5_evec_loc )
 
   ! Allocate the arrays for Model local BE:
      ALLOCATE ( be1_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be2_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be3_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be4_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be5_eval_loc (1:jy,1:1) )
      ALLOCATE ( be1_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be2_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be3_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be4_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be5_evec_loc(1:jy,1:1,1:1) )
 
   endif ! if interpolate_stats

   ! 4.0 Check and get the truncated number of the vertical modes
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (vert_corr == vert_corr_2) then

      ! 4.1 Perform checks on eigenvectors:

      if (test_statistics) then
         call da_check_eof_decomposition(be1_eval_glo(:), be1_evec_glo(:,:), be % v1 % name)
         call da_check_eof_decomposition(be2_eval_glo(:), be2_evec_glo(:,:), be % v2 % name)
         call da_check_eof_decomposition(be3_eval_glo(:), be3_evec_glo(:,:), be % v3 % name)
         call da_check_eof_decomposition(be4_eval_glo(:), be4_evec_glo(:,:), be % v4 % name)
      end if

      ! 4.2 Truncate in vertical:

      call da_get_vertical_truncation(max_vert_var1, be1_eval_glo(:), be % v1)
      call da_get_vertical_truncation(max_vert_var2, be2_eval_glo(:), be % v2)
      call da_get_vertical_truncation(max_vert_var3, be3_eval_glo(:), be % v3)
      call da_get_vertical_truncation(max_vert_var4, be4_eval_glo(:), be % v4)
      if (max_vert_var5 == 0.0) then
         be % v5 % mz = 0
      else
         be % v5 % mz = 1
      end if

      write (unit=stdout,fmt=*) ' '

   else

      ! 4.3 no truncated

      be % v1 % mz = xb % mkz
      be % v2 % mz = xb % mkz
      be % v3 % mz = xb % mkz
      be % v4 % mz = xb % mkz
      be % v5 % mz = xb % mkz
   end if

   kzs=(/kz,kz,kz,kz,1/)
   mzs=(/be%v1%mz,be%v2%mz,be%v3%mz,be%v4%mz,be%v5%mz/)

   ! 5.0 Initialise control variable space components of header:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! 5.1 Compute the size of the control variables

   be % mix = ix
   be % mjy = jy

   ! 5.2 Transfer errors to error structure:

   call da_allocate_background_errors(jy, kz, be1_eval_glo, be1_evec_glo, be1_eval_loc, &
                                       be1_evec_loc, be % v1)
   call da_allocate_background_errors(jy, kz, be2_eval_glo, be2_evec_glo, be2_eval_loc, &
                                       be2_evec_loc, be % v2)
   call da_allocate_background_errors(jy, kz, be3_eval_glo, be3_evec_glo, be3_eval_loc, &
                                       be3_evec_loc, be % v3)
   call da_allocate_background_errors(jy, kz, be4_eval_glo, be4_evec_glo, be4_eval_loc, &
                                       be4_evec_loc, be % v4)

   ! 5.2.1 transfer the ps_u variance to be % v5:

   call da_allocate_background_errors(jy,  1, be5_eval_glo, be5_evec_glo, be5_eval_loc, &
                                       be5_evec_loc, be % v5)
   if (print_detail_be) then
      write (unit=stdout,fmt='(3x,a,i10)') "b) Finished eigenvector processing!"
   end if

   if( use_rf )then		! use recursive filters:

   ! 5.3 Convert the scale lengths in the real distance (meter)

      be1_rf_lengthscale(1:kz) = be1_rf_lengthscale(1:kz) * xb%ds
      be2_rf_lengthscale(1:kz) = be2_rf_lengthscale(1:kz) * xb%ds
      be3_rf_lengthscale(1:kz) = be3_rf_lengthscale(1:kz) * xb%ds
      be4_rf_lengthscale(1:kz) = be4_rf_lengthscale(1:kz) * xb%ds
      be5_rf_lengthscale(1:1)  = be5_rf_lengthscale(1:1)  * xb%ds
   else				! use wavelets:
      read(be_unit)k		! trimmed string length
      if( k>len(variable) )then
         write(message(1),'(i0,">",i0," for use_rf=.false.")')k,len(variable)
         call da_error("da_setup_be_regional.inc",1365,message(1:1))
      endif
      read(be_unit)variable(1:k)
      if( variable(1:k) /= 'wavelet' )then
         write(message(1),'(A,"/=''wavelet''")')variable(1:k)
         call da_error("da_setup_be_regional.inc",1370,message(1:1))
      endif
      b = 0			! pointer to variable*mode
      do m=1,5			! variables loop:
         read(be_unit)k,nk
         if( k>len(variable) )then
            write(message(1),'(i0,">",i0," for ",A)')k,len(variable),namv(m)
            call da_error("da_setup_be_regional.inc",1377,message(1:1))
         elseif( nk /= kzs(m) )then
            write(message(1),'(i0,"=nk/=kzs(",i0,")=",i0," for ",A)')nk,m,kzs(m),namv(m)
            call da_error("da_setup_be_regional.inc",1380,message(1:1))
         endif
         read(be_unit)variable(1:k)
         if( variable(1:k) /= trim(namv(m)) )then
            write(message(1),'(A,"/=",A)')variable(1:k),trim(namv(m))
            call da_error("da_setup_be_regional.inc",1385,message(1:1))
         endif
         if( m==1 )then
!           Possibly reassign namelist do_normalize value:
            read(be_unit)do_normalize,namw,lf,nb
            allocate(nij(0:nb,0:1,0:2))
            read(be_unit)nij	! wavelet indexes
            write(*,'(A,": ")',advance="no")"da_setup_be_regional.inc"
            do i=0,nb		! wavelet-band loop:
               write(*,'(i2,"{",2(i3,","),i3,";",2(i3,","),i3,"}")',advance="no")&
                  i,transpose(nij(i,:,:))
            enddo
            allocate(be%wsd(nij(0,1,2),nij(0,0,2),sum(mzs)), &
!____________________max of {,i}dwtai wavelet scratch workspace sizes:
                     ws(max(maxval(nij(0,:,2)), &
                            2*floor(.5*(real(maxval(nij(0,:,0)))+lf))+lf-2)))
            if( do_normalize )allocate(be%sd(nij(0,1,0),nij(0,0,0),sum(mzs)))
         endif			! if( m==1 )
         do k=1,mzs(m)
!___________mode-k field & wavelet-coefficient std. devs.:
            read(be_unit)be%wsd(:,:,b+k)
            if( do_normalize )read(be_unit)be%sd(:,:,b+k)
         enddo			! mode-k loop.
         do k=mzs(m)+1,kzs(m)	! read and discard truncated modes:
            read(be_unit)(junk,j=1,nij(0,1,2)*nij(0,0,2))
            if( do_normalize )read(be_unit)(junk,j=1,nij(0,1,0)*nij(0,0,0))
         enddo			! mode-k loop.
         write(*,'(A,": |wsd[",A,"]|=",es9.3)')"da_setup_be_regional.inc",namv(m),sqrt(sum(be%wsd(:,:,b+1:b+mzs(m))**2))
         if( do_normalize )write(*,'(A,": | sd[",A,"]|=",es9.3)')"da_setup_be_regional.inc",namv(m),sqrt(sum(be%sd(:,:,b+1:b+mzs(m))**2))
         b = b+mzs(m)		! point to next variable.
      enddo			! variables loop
   endif			! if( use_rf )

   ! 6.0 Perform checks on eigenvectors with be data structure:
   if (jb_factor > 0.0 .and. test_statistics) then
      call da_check_eof_decomposition(be%v1%val_g(:), be%v1%evec_g(:,:),&
                                     be%v1%name)
      call da_check_eof_decomposition(be%v2%val_g(:), be%v2%evec_g(:,:),&
                                     be%v2%name)
      call da_check_eof_decomposition(be%v3%val_g(:), be%v3%evec_g(:,:),&
                                     be%v3%name)
      call da_check_eof_decomposition(be%v4%val_g(:), be%v4%evec_g(:,:),&
                                     be%v4%name)
   end if

   ! 6.1 Close the be unit

   close(be_unit)
   call da_free_unit(be_unit)

   if( use_rf )then
   ! 6.2 Keep the original be % v1, be % v2,...., and lengthscale in the first loop
   !     for the rescaling in the later loops:

      it = 1

      if (max_ext_its > 1 .and. jb_factor > 0.0) then

          write(unit=message(1),fmt='(A,I4)') '>>> Save the variances and scale-lengths in outer-loop', it
          call da_message(message(1:1))
          write(be_rf_unit)  kz, jy, ix, be % v1 % mz, be % v2 % mz, be% v3 % mz, &
                                     be % v4 % mz, be % v5 % mz, xb % ds
          write(be_rf_unit) be % v1 % val, be % v2 % val, be% v3 % val, &
                                     be % v4 % val, be % v5 % val, &
             be1_rf_lengthscale, be2_rf_lengthscale, be3_rf_lengthscale, &
             be4_rf_lengthscale, be5_rf_lengthscale
     
          if (print_detail_be ) then
             write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                               it, kz, jy, ix, xb % ds
             write(be_print_unit,'("Original val and rf, and mz:",5i5)') &
                      be % v1 % mz, be % v2 % mz, be% v3 % mz, be % v4 % mz, be % v5 % mz
             write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)
             write(be_print_unit,'(/"scale-length: kz=",i3)') kz
             do i = 1,kz 
               if (i == 1) then
                 write(be_print_unit,'(i3,2x,5e15.5)') i,be1_rf_lengthscale(i), &
                   be2_rf_lengthscale(i), be3_rf_lengthscale(i), be4_rf_lengthscale(i), &
                   be5_rf_lengthscale(i)
               else
                 write(be_print_unit,'(i3,2x,5e15.5)') i,be1_rf_lengthscale(i), &
                   be2_rf_lengthscale(i), be3_rf_lengthscale(i), be4_rf_lengthscale(i)
               endif
             enddo
    
          endif

      endif

      ! 7.0 Apply empirical and recursive filter rescaling factor:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call da_rescale_background_errors(var_scaling1(1), len_scaling1(1), &
                                         xb % ds, be1_rf_lengthscale, be % v1)
      call da_rescale_background_errors(var_scaling2(1), len_scaling2(1), &
                                         xb % ds, be2_rf_lengthscale, be % v2)
      call da_rescale_background_errors(var_scaling3(1), len_scaling3(1), &
                                         xb % ds, be3_rf_lengthscale, be % v3)
      call da_rescale_background_errors(var_scaling4(1), len_scaling4(1), &
                                         xb % ds, be4_rf_lengthscale, be % v4)
      call da_rescale_background_errors(var_scaling5(1), len_scaling5(1), &
                                         xb % ds, be5_rf_lengthscale, be % v5)
   endif

   if (print_detail_be .and. jb_factor > 0.0) then
       write(be_print_unit,'(/"============================================================")')
       write(be_print_unit,'("For outer loop ",i2)') it
       write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                                      it, kz, jy, ix, xb % ds
       write(be_print_unit,'("Namelist options specified for this iteration:")')
       write(be_print_unit,'("var_scaling1(it) = ",e12.5,2x,"len_scaling1(it) = "e12.5)')var_scaling1(it),len_scaling1(it)
       write(be_print_unit,'("var_scaling2(it) = ",e12.5,2x,"len_scaling2(it) = "e12.5)')var_scaling2(it),len_scaling2(it)
       write(be_print_unit,'("var_scaling3(it) = ",e12.5,2x,"len_scaling3(it) = "e12.5)')var_scaling3(it),len_scaling3(it)
       write(be_print_unit,'("var_scaling4(it) = ",e12.5,2x,"len_scaling4(it) = "e12.5)')var_scaling4(it),len_scaling4(it)
       write(be_print_unit,'("var_scaling5(it) = ",e12.5,2x,"len_scaling5(it) = "e12.5)')var_scaling5(it),len_scaling5(it)
       write(be_print_unit,'("Background error statistics for this iteration:")')
       write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)
       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       write(be_print_unit,'("be%v1%rf_alpha:"/(10e12.5))') be % v1 % rf_alpha(:)
       write(be_print_unit,'("be%v2%rf_alpha:"/(10e12.5))') be % v2 % rf_alpha(:)
       write(be_print_unit,'("be%v3%rf_alpha:"/(10e12.5))') be % v3 % rf_alpha(:)
       write(be_print_unit,'("be%v4%rf_alpha:"/(10e12.5))') be % v4 % rf_alpha(:)
       write(be_print_unit,'("be%v5%rf_alpha:"/(10e12.5))') be % v5 % rf_alpha(:)
       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       do i = 1,kz
          if (i == 1) then
             write(be_print_unit,'(i3,2x,5e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                  be3_rf_lengthscale(i), be4_rf_lengthscale(i), be5_rf_lengthscale(i)
          else
             write(be_print_unit,'(i3,2x,4e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                  be3_rf_lengthscale(i), be4_rf_lengthscale(i)
          endif
       enddo

   endif

   ! 8.0 deallocate input model state:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   deallocate (be1_eval_loc)
   deallocate (be2_eval_loc)
   deallocate (be3_eval_loc)
   deallocate (be4_eval_loc)
   deallocate (be5_eval_loc)
   if( use_rf )then
      deallocate (be1_rf_lengthscale)
      deallocate (be2_rf_lengthscale)
      deallocate (be3_rf_lengthscale)
      deallocate (be4_rf_lengthscale)
      deallocate (be5_rf_lengthscale)
   endif
   if (vert_corr == vert_corr_2) then
      deallocate (be1_eval_glo)
      deallocate (be2_eval_glo)
      deallocate (be3_eval_glo)
      deallocate (be4_eval_glo)
      deallocate (be5_eval_glo)
      deallocate (be1_evec_loc)
      deallocate (be2_evec_loc)
      deallocate (be3_evec_loc)
      deallocate (be4_evec_loc)
      deallocate (be5_evec_loc)
      deallocate (be1_evec_glo)
      deallocate (be2_evec_glo)
      deallocate (be3_evec_glo)
      deallocate (be4_evec_glo)
      deallocate (be5_evec_glo)
   end if

   deallocate (bin)
   deallocate (bin2d)

   if (be % ne > 0) then
      be % alpha % name = 'alpha'
      allocate (alpha_val(1:nk)) ! Not using jy dimension yet.
      allocate (alpha_evec(1:nk,1:nk)) ! Not using jy dimension yet.

      if ( alpha_vertloc ) then ! Use vertical localization:
         call da_get_unit(be_unit)
         open(unit=be_unit,file='be.vertloc.dat', status='old', form='unformatted')

         read (be_unit)nk
         read (be_unit)alpha_val(1:nk)
         read (be_unit)alpha_evec(1:nk,1:nk)
         close(be_unit)
         call da_free_unit(be_unit)
         be % alpha % mz = nk

         call da_get_vertical_truncation(max_vert_var_alpha, alpha_val, be % alpha)
      else
         be % alpha % mz = 1 ! No vertical localization.
         alpha_val(1) = 1.0
         alpha_val(2:kz) = 0.0
         alpha_evec(:,:) = 1.0
      end if
      mz = be % alpha % mz

!     Alpha eigenvalues and eigenvectors:
      allocate (be % alpha % val(1:jy,1:mz)) ! Not using jy dimension but here for consistency.
      allocate (be % alpha % evec(1:nj,1:nk,1:mz))

      if ( anal_type_hybrid_dual_res ) then
         deallocate(be % alpha % val, be % alpha % evec)
         allocate (be % alpha % val(jds_int:jde_int,1:mz))
         allocate (be % alpha % evec(jds_int:jde_int,1:nk,1:mz))
      endif


      do m = 1, mz
         be % alpha % val(:,m) = sigma_alpha * alpha_val(m)
         do k = 1, nk
            be % alpha % evec(:,k,m) = alpha_evec(k,m)
         end do
      end do

!     Alpha RF lengthscales and variance scaling factors:
      allocate (alpha_rf_lengthscale(1:mz))
      allocate (be % alpha % rf_alpha(1:mz))
      allocate (alpha_rf_scale_factor(1:mz))


      if ( anal_type_hybrid_dual_res ) then
         alpha_rf_lengthscale(1:mz) = 1000.0 * alpha_corr_scale / (grid%intermediate_grid % dx )
      else
         alpha_rf_lengthscale(1:mz) = 1000.0 * alpha_corr_scale / xb % ds ! Convert km to grid spacings.
      endif


      call da_calculate_rf_factors( alpha_rf_lengthscale(:), be % alpha % rf_alpha(:), &
                                    alpha_rf_scale_factor(:) )
      do m = 1, mz
         be % alpha % val(:,m) = alpha_rf_scale_factor(m) * be % alpha % val(:,m)
      end do

      if( .not. use_rf ) then
         allocate(be%alpha%wsd(nij(0,1,2),nij(0,0,2),mz))
         call random_number(be%alpha%wsd)			! need to parallelize
         if( do_normalize )then
            allocate(be%alpha%sd(nij(0,1,0),nij(0,0,0),mz))
            call random_number(be%alpha%sd)			! need to parallelize
         endif
      endif

      deallocate(alpha_val)
      deallocate(alpha_evec)
      deallocate(alpha_rf_lengthscale)
      deallocate(alpha_rf_scale_factor)

   else
      be % alpha % mz = 0
   end if

   if (trace_use) call da_trace_exit("da_setup_be_regional")

end subroutine da_setup_be_regional


subroutine da_setup_be_nmm_regional(xb, be)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors
   ! 
   ! Updates: 
   !
   !       Implementation of multi-variate BE  
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb                    ! First guess structure.
   type (be_type), intent(inout) :: be                    ! Back. errors structure.

   integer                     :: i, j, k, m       ! Loop counters.
   integer, allocatable:: bin(:,:,:)         ! Bin assigned to each 3D point
   integer, allocatable:: bin2d(:,:)         ! Bin assigned to each 2D point
   integer             :: bin_type           ! Type of bin to average over.
   integer             :: num_bins           ! Number of bins (3D fields).
   integer             :: num_bins2d         ! Number of bins (3D fields).
   real    :: lat_min, lat_max, binwidth_lat ! Used if bin_type = 2 (degrees)..
   real    :: hgt_min, hgt_max, binwidth_hgt ! Used if bin_type = 2 (m). .

   real*8, allocatable         :: be1_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable         :: be2_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable         :: be3_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable         :: be4_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable         :: be5_eval_loc(:,:)     ! Temp arrays.

   real*8, allocatable         :: be1_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable         :: be2_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable         :: be3_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable         :: be4_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable         :: be5_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable         :: alpha_val(:)          ! Global Eigenvalues.

   real*8, allocatable         :: be1_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable         :: be2_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable         :: be3_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable         :: be4_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable         :: be5_evec_loc(:,:,:)   ! Local Eigenvectors.

   real*8, allocatable         :: be1_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable         :: be2_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable         :: be3_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable         :: be4_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable         :: be5_evec_glo(:,:)     ! Global Eigenvectors.
   real, allocatable           :: alpha_evec(:,:)       ! Global Eigenvectors.

   real*8, allocatable         :: be1_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable         :: be2_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable         :: be3_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable         :: be4_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable         :: be5_rf_lengthscale(:)
   real*8, allocatable         :: alpha_rf_lengthscale(:)
   real*8, allocatable         :: alpha_rf_scale_factor(:)

   real, allocatable           :: evec_loc(:,:,:)        ! Latitudinally varying eigenvectors.
   real, allocatable           :: eval_loc(:,:)          ! Latitudinally varying eigenvalues.

   character*10                :: variable
   character*80                :: var80c    
   integer                     :: ni, nj, nk, nk_2d, b         
   integer                     :: ix, jy, kz, mz
   real, allocatable   :: regcoeff_psi_chi(:)        ! psi/chi    regression cooefficient.
   real, allocatable   :: regcoeff_psi_t(:,:,:)      ! psi/t   regression cooefficient.
   real, allocatable   :: regcoeff_psi_ps(:,:)       ! psi/ps     regression cooefficient.
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     ! psi/rh     regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    ! chi_u/t regression coefficient
   real, allocatable   :: regcoeff_chi_u_ps(:,:)     ! chi_u/ps   regression coefficient
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   ! chi_u/rh   regression coefficient
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     ! t_u/rh  regression coefficient
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      ! ps_u/rh    regression coefficient

   integer                     :: be_unit, ier, be_rf_unit, be_print_unit, it

   if (trace_use) call da_trace_entry("da_setup_be_nmm_regional")
   
   if( (cv_options == 5) .or. (cv_options == 6) .or. (cv_options == 7) ) then
      write (unit=message(1),fmt='(A,i3)') &
          'Set up background errors for regional application for cv_option= ',cv_options
      call da_message(message(1:1))
   else 
      write(unit=message(1),fmt='(A,I4)') &
          'This subroutine is for cv_options = 5, 6, or 7; you have selected cv_options = ', cv_options
      call da_error("da_setup_be_nmm_regional.inc",88,message(1:1))
   end if

   ix = xb % mix
   jy = xb % mjy
   kz = xb % mkz

   be_rf_unit    = unit_end + 1
   be_print_unit = unit_end + 2

   call da_get_unit(be_unit)
   open(unit=be_unit,file="be.dat", status="old",form="unformatted")

   rewind (be_unit)
   read (be_unit, iostat=ier) ni, nj, nk
   if (ier /= 0) then
      write(unit=message(1),fmt='(A,I4,A,I4)') 'cv_options:', cv_options,' Reading error in unit=',be_unit
      call da_error("da_setup_be_nmm_regional.inc",105,message(1:1))
   else
      if ( nk /= kz ) then
         call da_error("da_setup_be_nmm_regional.inc",108,  &
            (/"Vertical levels in fg and be.dat do not match."/))
      end if
   endif

   allocate (bin(1:ni,1:nj,1:nk))
   allocate (bin2d(1:ni,1:nj))

   read (be_unit)bin_type
   read (be_unit)lat_min, lat_max, binwidth_lat
   read (be_unit)hgt_min, hgt_max, binwidth_hgt
   read (be_unit)num_bins, num_bins2d
   read (be_unit)bin(1:ni,1:nj,1:nk)
   read (be_unit)bin2d(1:ni,1:nj)

   ! 1.1 Read in regression coefficients
   allocate  (regcoeff_psi_chi(1:num_bins))
   allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))

   regcoeff_psi_chi    = 0.
   regcoeff_psi_t      = 0.
   regcoeff_psi_ps     = 0.
   regcoeff_psi_rh     = 0.
   regcoeff_chi_u_t    = 0.
   regcoeff_chi_u_ps   = 0.
   regcoeff_chi_u_rh   = 0.
   regcoeff_t_u_rh     = 0.
   regcoeff_ps_u_rh    = 0.

   if ( (cv_options == 5) .or. (cv_options == 7) ) then
   read (be_unit) regcoeff_psi_chi
   read (be_unit) regcoeff_psi_ps
   read (be_unit) regcoeff_psi_t
   else
   do i = 1 , 9
   read (be_unit) var80c
   select case( trim(adjustl(var80c)) )

   case ('regcoeff_psi_chi')
   read (be_unit) regcoeff_psi_chi
   case ('regcoeff_psi_t')
   read (be_unit) regcoeff_psi_t
   case ('regcoeff_psi_ps')
   read (be_unit) regcoeff_psi_ps
   case ('regcoeff_psi_rh')
   read (be_unit) regcoeff_psi_rh
   case ('regcoeff_chi_u_t')
   read (be_unit) regcoeff_chi_u_t
   case ('regcoeff_chi_u_ps')
   read (be_unit) regcoeff_chi_u_ps
   case ('regcoeff_chi_u_rh')
   read (be_unit) regcoeff_chi_u_rh
   case ('regcoeff_t_u_rh')
   read (be_unit) regcoeff_t_u_rh
   case ('regcoeff_ps_u_rh')
   read (be_unit) regcoeff_ps_u_rh
   case default;
      message(1)=' Read problem in regression coefficients in BE file '
      write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
      call da_error("da_setup_be_nmm_regional.inc",174,message(1:2))
   end select
   end do
   end if

   ! 1.2 Fill regression coeff. array
   allocate (be%reg_psi_chi (1:jy,1:nk))
   allocate (be%reg_psi_t   (1:jy,1:nk,1:nk))
   allocate (be%reg_psi_ps  (1:jy,1:nk))
   allocate (be%reg_psi_rh  (1:jy,1:nk,1:nk))
   allocate (be%reg_chi_u_t (1:jy,1:nk,1:nk))
   allocate (be%reg_chi_u_ps(1:jy,1:nk))
   allocate (be%reg_chi_u_rh(1:jy,1:nk,1:nk))
   allocate (be%reg_t_u_rh  (1:jy,1:nk,1:nk))
   allocate (be%reg_ps_u_rh (1:jy,1:nk))

   be%reg_psi_chi = 0.
   be%reg_psi_t   = 0.
   be%reg_psi_ps  = 0.
   be%reg_psi_rh  = 0.
   be%reg_chi_u_t = 0.
   be%reg_chi_u_ps= 0.
   be%reg_chi_u_rh= 0.
   be%reg_t_u_rh  = 0.

   do k=1,nk
      do j =1, jy
         b = bin(1,1,k)
         be%reg_psi_chi(j,k) = psi_chi_factor * regcoeff_psi_chi(b)
      end do
   end do

   do j=1,jy
      b = bin2d(1,1)
      do k=1,nk
         be%reg_psi_ps(j,k)   = psi_ps_factor   * regcoeff_psi_ps(k,b)
         be%reg_ps_u_rh(j,k)  = ps_u_rh_factor  * regcoeff_ps_u_rh(k,b)
         be%reg_chi_u_ps(j,k) = chi_u_ps_factor * regcoeff_chi_u_ps(k,b)
      end do
   end do

   do j=1,jy
      b = bin2d(1,1)
      do i=1,nk
         do k=1,nk
            be%reg_psi_t(j,i,k)   = psi_t_factor   *  regcoeff_psi_t(i,k,b)
            be%reg_psi_rh(j,i,k)  = psi_rh_factor  *  regcoeff_psi_rh(i,k,b)
            be%reg_chi_u_t(j,i,k) = chi_u_t_factor *  regcoeff_chi_u_t(i,k,b)
            be%reg_chi_u_rh(j,i,k)= chi_u_rh_factor*  regcoeff_chi_u_rh(i,k,b)
            be%reg_t_u_rh(j,i,k)  = t_u_rh_factor  *  regcoeff_t_u_rh(i,k,b)
         end do
      end do
   end do

   deallocate (regcoeff_psi_chi)
   deallocate (regcoeff_psi_t)
   deallocate (regcoeff_psi_ps)
   deallocate (regcoeff_psi_rh)
   deallocate (regcoeff_chi_u_t)
   deallocate (regcoeff_chi_u_ps)
   deallocate (regcoeff_chi_u_rh)
   deallocate (regcoeff_t_u_rh)
   deallocate (regcoeff_ps_u_rh)

   ! 1.3 Domain_averaged regression coefficients

  if (.not.lat_stats_option) then
      write (unit=message(4), fmt='(3x, a)') &
         'Using the averaged regression coefficients for unbalanced part'

      do k=1,nk
         be%reg_psi_ps  (1:num_bins2d,k) = sum(be%reg_psi_ps  (1:num_bins2d,k))/dble(num_bins2d)
         be%reg_ps_u_rh (1:num_bins2d,k) = sum(be%reg_ps_u_rh (1:num_bins2d,k))/dble(num_bins2d)
         be%reg_chi_u_ps(1:num_bins2d,k) = sum(be%reg_chi_u_ps(1:num_bins2d,k))/dble(num_bins2d)
      end do

      do m=1,nk
         do k=1,nk

           be%reg_psi_t  (1:num_bins2d,k,m)= sum(be%reg_psi_t(1:num_bins2d,k,m))/dble(num_bins2d)
           be%reg_psi_rh (1:num_bins2d,k,m)= sum(be%reg_psi_rh  (1:num_bins2d,k,m))/dble(num_bins2d)
           be%reg_chi_u_t(1:num_bins2d,k,m)= sum(be%reg_chi_u_t(1:num_bins2d,k,m))/dble(num_bins2d)
           be%reg_chi_u_rh(1:num_bins2d,k,m)= sum(be%reg_chi_u_rh  (1:num_bins2d,k,m))/dble(num_bins2d)
           be%reg_t_u_rh (1:num_bins2d,k,m)= sum(be%reg_t_u_rh (1:num_bins2d,k,m))/dble(num_bins2d)

         end do
      end do
   else
      write (unit=message(4), fmt='(3x, a)') &
         'Using the latitude-dependent regression coefficients for unbalanced part'
   end if

   call da_message(message(1:4))

   ! 2.0 Load the eigenvector and eigenvalue

   allocate (be1_eval_loc (1:jy,1:kz))
   allocate (be2_eval_loc (1:jy,1:kz))
   allocate (be3_eval_loc (1:jy,1:kz))
   allocate (be4_eval_loc (1:jy,1:kz))
   allocate (be5_eval_loc (1:jy,1:1))

   if (vert_corr == vert_corr_2) then

      allocate (be1_eval_glo(1:kz))
      allocate (be2_eval_glo(1:kz))
      allocate (be3_eval_glo(1:kz))
      allocate (be4_eval_glo(1:kz))
      allocate (be5_eval_glo(1:1))

      allocate (be1_evec_loc(1:jy,1:kz,1:kz))
      allocate (be2_evec_loc(1:jy,1:kz,1:kz))
      allocate (be3_evec_loc(1:jy,1:kz,1:kz))
      allocate (be4_evec_loc(1:jy,1:kz,1:kz))
      allocate (be5_evec_loc(1:jy,1: 1,1: 1))

      allocate (be1_evec_glo(1:kz,1:kz))
      allocate (be2_evec_glo(1:kz,1:kz))
      allocate (be3_evec_glo(1:kz,1:kz))
      allocate (be4_evec_glo(1:kz,1:kz))
      allocate (be5_evec_glo(1:1,1:1))
   end if

   ! 2.2 Read in the eigenvector and eigenvalue 

   do i = 1 , 4
    read (be_unit) variable
   select case( trim(adjustl(variable)) )
   case ('psi')
   read (be_unit) nk, num_bins2d
   read (be_unit)  be1_evec_glo
   read (be_unit)  be1_eval_glo
   if( i == 1) then
   allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
   allocate (eval_loc(1:nk,     1:num_bins2d))
   end if
   read (be_unit)  evec_loc
   read (be_unit)  eval_loc
   do j=1,jy
      b = bin2d(1,1)
      be1_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
      be1_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
   end do
   be % v1 % name = trim(adjustl(variable))

   case ('chi_u')
   read (be_unit) nk, num_bins2d
   read (be_unit)  be2_evec_glo
   read (be_unit)  be2_eval_glo
   read (be_unit)  evec_loc
   read (be_unit)  eval_loc
   do j=1,jy
      b = bin2d(1,1)
      be2_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
      be2_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
   end do
   be % v2 % name = trim(adjustl(variable))

   case ('t_u')
   read (be_unit) nk, num_bins2d
   read (be_unit)  be3_evec_glo
   read (be_unit)  be3_eval_glo
   read (be_unit)  evec_loc
   read (be_unit)  eval_loc
   do j=1,jy
      b = bin2d(1,1)
      be3_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
      be3_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
   end do
   be % v3 % name = trim(adjustl(variable))

   case ('rh_u' , 'rh' )
   read (be_unit) nk, num_bins2d
   read (be_unit)  be4_evec_glo
   read (be_unit)  be4_eval_glo
   read (be_unit)  evec_loc
   read (be_unit)  eval_loc
   do j=1,jy
      b = bin2d(1,1)
      be4_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
      be4_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
   end do
   be % v4 % name = trim(adjustl(variable))
   case default;
      message(1)=' Read problem in eigen vaectors/values in BE file '
      write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
      call da_error("da_setup_be_nmm_regional.inc",360,message(1:2))
   end select
   end do

   deallocate (evec_loc)
   deallocate (eval_loc)

   ! 2.2.5 Control variable ps_u
   read (be_unit) variable
   read (be_unit) nk_2d, num_bins2d
   allocate (evec_loc(1:nk_2d,1:nk_2d,1:num_bins2d))
   allocate (eval_loc(1:nk_2d,        1:num_bins2d))
   read (be_unit)  be5_evec_glo
   read (be_unit)  be5_eval_glo
   read (be_unit)  evec_loc
   read (be_unit)  eval_loc
   if( trim(adjustl(variable)) /= 'ps_u' ) then
      message(1)=' Variable mismatch while transfering eigen values from BE file '
      write (unit=message(2),fmt='(A,A)') ' Expected ps_u but got ',trim(adjustl(variable))
      call da_error("da_setup_be_nmm_regional.inc",379,message(1:2))
   end if
   be % v5 % name = trim(adjustl(variable))
   do j=1,jy
      b = bin2d(1,1)
      be5_evec_loc(j,1:1,1:1) = evec_loc(1:1,1:1,b)
      be5_eval_loc(j,1:1 ) = eval_loc(1:1,b)
   end do

   deallocate (evec_loc)
   deallocate (eval_loc)

!
   if(use_radarobs .and. use_radar_rf .or. use_rad .and. crtm_cloud) then  
      if ( cloud_cv_options == 1 ) be % v4 % name = 'qt   '
   end if

   write (unit=message(2),fmt='(3x,A,7A)') 'WRF-Var dry control variables are:', &
      trim(be % v1 % name), ', ', trim(be % v2 % name), ', ', &
      trim(be % v3 % name), ' and ', trim(be % v5 % name)

   write (unit=message(3),fmt='(3x,A,A)') &
      'Humidity control variable is ', trim(be % v4 % name)
!
   ! 3.0 Check and get the truncated number of the vertical modes
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (vert_corr == vert_corr_2) then

      ! 3.1 Perform checks on eigenvectors:

      if (test_statistics) then
         call da_check_eof_decomposition(be1_eval_glo(:), be1_evec_glo(:,:), be % v1 % name)
         call da_check_eof_decomposition(be2_eval_glo(:), be2_evec_glo(:,:), be % v2 % name)
         call da_check_eof_decomposition(be3_eval_glo(:), be3_evec_glo(:,:), be % v3 % name)
         call da_check_eof_decomposition(be4_eval_glo(:), be4_evec_glo(:,:), be % v4 % name)
      end if

      ! 3.2 Truncate in vertical:

      call da_get_vertical_truncation(max_vert_var1, be1_eval_glo(:), be % v1)
      call da_get_vertical_truncation(max_vert_var2, be2_eval_glo(:), be % v2)
      call da_get_vertical_truncation(max_vert_var3, be3_eval_glo(:), be % v3)
      call da_get_vertical_truncation(max_vert_var4, be4_eval_glo(:), be % v4)

      if (max_vert_var5 == 0.0) then
         be % v5 % mz = 0
      else
         be % v5 % mz = 1
      end if

      write (unit=stdout,fmt=*) ' '

   else

      ! 3.3 no truncated

      be % v1 % mz = xb % mkz
      be % v2 % mz = xb % mkz
      be % v3 % mz = xb % mkz
      be % v4 % mz = xb % mkz
      be % v5 % mz = xb % mkz

   end if

   ! 4.0 Initialise control variable space components of header:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! 4.1 Compute the size of the control variables

   be % mix = ix
   be % mjy = jy

   ! 4.2 Transfer errors to error structure:

   call da_allocate_background_errors(jy, kz, be1_eval_glo, be1_evec_glo, be1_eval_loc, &
                                       be1_evec_loc, be % v1)
   call da_allocate_background_errors(jy, kz, be2_eval_glo, be2_evec_glo, be2_eval_loc, &
                                       be2_evec_loc, be % v2)
   call da_allocate_background_errors(jy, kz, be3_eval_glo, be3_evec_glo, be3_eval_loc, &
                                       be3_evec_loc, be % v3)
   call da_allocate_background_errors(jy, kz, be4_eval_glo, be4_evec_glo, be4_eval_loc, &
                                       be4_evec_loc, be % v4)

   ! 4.2.1 transfer the ps_u variance to be % v5:

   call da_allocate_background_errors(jy,  1, be5_eval_glo, be5_evec_glo, be5_eval_loc, &
                                       be5_evec_loc, be % v5)

   if (print_detail_be) then
      write (unit=stdout,fmt='(3x,a,i10)') "b) Finished eigenvector processing!"
   end if

   ! 5.0 Load the scale lengths
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! 5.1 allocate the array for scale lengths

   allocate (be1_rf_lengthscale(1:nk))
   allocate (be2_rf_lengthscale(1:nk))
   allocate (be3_rf_lengthscale(1:nk))
   allocate (be4_rf_lengthscale(1:nk))
   allocate (be5_rf_lengthscale(1:nk))

   ! 5.2 read in the scale lengths
   do i = 1 , 5
    read (be_unit) variable
   select case( trim(adjustl(variable)) )

   case ('psi')
     read(be_unit) be1_rf_lengthscale
   case ('chi_u')
     read(be_unit) be2_rf_lengthscale
   case ('t_u')
     read(be_unit) be3_rf_lengthscale
   case ('rh_u' , 'rh')
     read(be_unit) be4_rf_lengthscale
   case ('ps_u')
     read(be_unit) be5_rf_lengthscale(1:1)
   case default;
      message(1)=' Read problem in lengthscales in BE file '
      write (unit=message(2),fmt='(A,A)') ' Trying to read lengthscales for variable: ',trim(adjustl(variable))
      call da_error("da_setup_be_nmm_regional.inc",501,message(1:2))
   end select
   end do


   ! 5.3 Convert the scale lengths in the real distance (meter)

   be1_rf_lengthscale(1:nk) = be1_rf_lengthscale(1:nk) * xb%ds
   be2_rf_lengthscale(1:nk) = be2_rf_lengthscale(1:nk) * xb%ds
   be3_rf_lengthscale(1:nk) = be3_rf_lengthscale(1:nk) * xb%ds
   be4_rf_lengthscale(1:nk) = be4_rf_lengthscale(1:nk) * xb%ds
   be5_rf_lengthscale(1:1)  = be5_rf_lengthscale(1:1)  * xb%ds

   ! 6.0 Perform checks on eigenvectors with be data structure:
   if (test_statistics) then
      call da_check_eof_decomposition(be%v1%val_g(:), be%v1%evec_g(:,:),&
                                     be%v1%name)
      call da_check_eof_decomposition(be%v2%val_g(:), be%v2%evec_g(:,:),&
                                     be%v2%name)
      call da_check_eof_decomposition(be%v3%val_g(:), be%v3%evec_g(:,:),&
                                     be%v3%name)
      call da_check_eof_decomposition(be%v4%val_g(:), be%v4%evec_g(:,:),&
                                     be%v4%name)
   end if

   ! 6.2 Close the be unit

   close(be_unit)
   call da_free_unit(be_unit)

   ! 6.3 Keep the original be % v1, be % v2,...., and lengthscale in the first loop
   !     for the rescaling in the later loops:

   it = 1

   if (max_ext_its > 1) then

       write(unit=message(1),fmt='(A,I4)') '>>> Save the variances and scale-lengths in outer-loop', it
       call da_message(message(1:1))
       write(be_rf_unit)  kz, jy, ix, be % v1 % mz, be % v2 % mz, be% v3 % mz, &
                                     be % v4 % mz, be % v5 % mz, xb % ds
       write(be_rf_unit) be % v1 % val, be % v2 % val, be% v3 % val, &
                                     be % v4 % val, be % v5 % val, &
             be1_rf_lengthscale, be2_rf_lengthscale, be3_rf_lengthscale, &
             be4_rf_lengthscale, be5_rf_lengthscale
     
       if (print_detail_be ) then
       write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                               it, kz, jy, ix, xb % ds
       write(be_print_unit,'("Original val and rf, and mz:",5i5)') &
                      be % v1 % mz, be % v2 % mz, be% v3 % mz, be % v4 % mz, be % v5 % mz
       write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)
       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       do i = 1,kz 
         if (i == 1) then
           write(be_print_unit,'(i3,2x,5e15.5)') i,be1_rf_lengthscale(i), &
                be2_rf_lengthscale(i), be3_rf_lengthscale(i), be4_rf_lengthscale(i), &
                be5_rf_lengthscale(i)
         else
           write(be_print_unit,'(i3,2x,5e15.5)') i,be1_rf_lengthscale(i), &
                be2_rf_lengthscale(i), be3_rf_lengthscale(i), be4_rf_lengthscale(i)
       endif
       enddo
    
       endif

   endif

   ! 7.0 Apply empirical and recursive filter rescaling factor:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call da_rescale_background_errors(var_scaling1(1), len_scaling1(1), &
                                      xb % ds, be1_rf_lengthscale, be % v1)
   call da_rescale_background_errors(var_scaling2(1), len_scaling2(1), &
                                      xb % ds, be2_rf_lengthscale, be % v2)
   call da_rescale_background_errors(var_scaling3(1), len_scaling3(1), &
                                      xb % ds, be3_rf_lengthscale, be % v3)
   call da_rescale_background_errors(var_scaling4(1), len_scaling4(1), &
                                      xb % ds, be4_rf_lengthscale, be % v4)
   call da_rescale_background_errors(var_scaling5(1), len_scaling5(1), &
                                      xb % ds, be5_rf_lengthscale, be % v5)
   if (print_detail_be ) then

       write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                                      it, kz, jy, ix, xb % ds
       write(be_print_unit,'("Loop it=",i2," val and rf, and mz:",5i5)') &
                  it, be % v1 % mz, be % v2 % mz, be% v3 % mz, be % v4 % mz, be % v5 % mz
       write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)
       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       do i = 1,kz 
         if (i == 1) then
           write(be_print_unit,'(i3,2x,5e15.5)') i, be % v1 % rf_alpha(i), &
                be % v2 % rf_alpha(i), be % v3 % rf_alpha(i), be % v4 % rf_alpha(i), &
                be % v5 % rf_alpha(i)
         else
           write(be_print_unit,'(i3,2x,5e15.5)') i, be % v1 % rf_alpha(i), &
                be % v2 % rf_alpha(i), be % v3 % rf_alpha(i), be % v4 % rf_alpha(i)
       endif
       enddo

   endif

   ! 8.0 deallocate input model state:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   deallocate (be1_eval_loc)
   deallocate (be2_eval_loc)
   deallocate (be3_eval_loc)
   deallocate (be4_eval_loc)
   deallocate (be5_eval_loc)

   deallocate (be1_rf_lengthscale)
   deallocate (be2_rf_lengthscale)
   deallocate (be3_rf_lengthscale)
   deallocate (be4_rf_lengthscale)
   deallocate (be5_rf_lengthscale)

   if (vert_corr == vert_corr_2) then
      deallocate (be1_eval_glo)
      deallocate (be2_eval_glo)
      deallocate (be3_eval_glo)
      deallocate (be4_eval_glo)
      deallocate (be5_eval_glo)

      deallocate (be1_evec_loc)
      deallocate (be2_evec_loc)
      deallocate (be3_evec_loc)
      deallocate (be4_evec_loc)
      deallocate (be5_evec_loc)

      deallocate (be1_evec_glo)
      deallocate (be2_evec_glo)
      deallocate (be3_evec_glo)
      deallocate (be4_evec_glo)
      deallocate (be5_evec_glo)

   end if

   deallocate (bin)
   deallocate (bin2d)

   if (be % ne > 0) then
      be % alpha % name = 'alpha'
      allocate (alpha_val(1:kz)) ! Not using jy dimension yet.
      allocate (alpha_evec(1:kz,1:kz)) ! Not using jy dimension yet.

      if ( alpha_vertloc ) then ! Use vertical localization:
         call da_get_unit(be_unit)
         open(unit=be_unit,file='be.vertloc.dat', status='old', form='unformatted')

         read (be_unit)nk
         read (be_unit)alpha_val(1:nk)
         read (be_unit)alpha_evec(1:nk,1:nk)
         close(be_unit)

         call da_get_vertical_truncation(max_vert_var_alpha, alpha_val, be % alpha)
      else
         be % alpha % mz = 1 ! No vertical localization.
         alpha_val(1) = 1.0
         alpha_val(2:kz) = 0.0
         alpha_evec(:,:) = 1.0
      end if
      mz = be % alpha % mz

!     Alpha eigenvalues and eigenvectors:
      allocate (be % alpha % val(1:jy,1:mz)) ! Not using jy dimension but here for consistency.
      allocate (be % alpha % evec(1:jy,1:kz,1:mz))
      do m = 1, mz
         be % alpha % val(:,m) = sigma_alpha * alpha_val(m)
         do k = 1, nk
            be % alpha % evec(:,k,m) = alpha_evec(k,m)
         end do
      end do

!     Alpha RF lengthscales and variance scaling factors:
      allocate (alpha_rf_lengthscale(1:mz))
      allocate (be % alpha % rf_alpha(1:mz))
      allocate (alpha_rf_scale_factor(1:mz))

      alpha_rf_lengthscale(1:mz) = 1000.0 * alpha_corr_scale / xb % ds ! Convert km to grid spacings.

      call da_calculate_rf_factors( alpha_rf_lengthscale(:), be % alpha % rf_alpha(:), &
                                    alpha_rf_scale_factor(:) )
      do m = 1, mz
         be % alpha % val(:,m) = alpha_rf_scale_factor(m) * be % alpha % val(:,m)
      end do

   else
      be % alpha % mz = 0
   end if

   if (trace_use) call da_trace_exit("da_setup_be_nmm_regional")

end subroutine da_setup_be_nmm_regional


subroutine da_setup_cv(be)

   !---------------------------------------------------------------------------
   ! Purpose: Calculate the size of the 1-dimensional control variable array.        
   !---------------------------------------------------------------------------

   implicit none

   type (be_type), intent(inout) :: be       ! background error.

   integer :: iy, jx   ! Local horizontal domain dimensions.
   integer :: ij       ! Product of horizontal dims.

   if (trace_use) call da_trace_entry("da_setup_cv")

   !--------------------------------------------------------------
   ! [1] Define standard control variable size:
   !--------------------------------------------------------------

   if (global) then
      be % cv % size1c = (be % v1 % max_wave+1) * (be % v1 % max_wave+2)/2
      be % cv % size2c = (be % v2 % max_wave+1) * (be % v2 % max_wave+2)/2
      be % cv % size3c = (be % v3 % max_wave+1) * (be % v3 % max_wave+2)/2
      be % cv % size4c = (be % v4 % max_wave+1) * (be % v4 % max_wave+2)/2
      be % cv % size5c = (be % v5 % max_wave+1) * (be % v5 % max_wave+2)/2

      be % cv % size1  = 2 * be % cv % size1c * be % v1 % mz
      be % cv % size2  = 2 * be % cv % size2c * be % v2 % mz
      be % cv % size3  = 2 * be % cv % size3c * be % v3 % mz
      be % cv % size4  = 2 * be % cv % size4c * be % v4 % mz
      be % cv % size5  = 2 * be % cv % size5c * be % v5 % mz
   else
      iy = ite - its + 1
      jx = jte - jts + 1
      if( use_rf )then		! Use recursive filters:
         ij = iy * jx		! cv space size same as state space size
      else			! cv space size slightly larger for dwtai2():
         ij = size(be%wsd,1)*size(be%wsd,2)
!        call da_setup_cv_iyjx2ij(iy,jx,ij,nij,nb,lf,"tile")
         write(unit=message(1),fmt='("da_setup_cv: using ",i0," 2D ",A1,i0," wavelets in ",i0," bands.")') ij,namw,lf,nb
         call da_message(message(1:1))
      endif			! if( use_rf )
      be % cv % size1 = ij * be % v1 % mz
      be % cv % size2 = ij * be % v2 % mz
      be % cv % size3 = ij * be % v3 % mz
      be % cv % size4 = ij * be % v4 % mz
      be % cv % size5 = ij * be % v5 % mz
   end if
   be % cv % size_jb = be % cv % size1 + be % cv % size2 + be % cv % size3 + &
      be % cv % size4 + be % cv % size5
   !--------------------------------------------------------------
   ! [1.1] Define 4D-Var lateral boundary condition control variable size:
   !--------------------------------------------------------------

   be % cv % size_jl = 0

   if ( .not. global .and. var4d ) then
      be % cv % size1l = be % cv % size1
      be % cv % size2l = be % cv % size2
      be % cv % size3l = be % cv % size3
      be % cv % size4l = be % cv % size4
      be % cv % size5l = be % cv % size5

      be % cv % size_jl = be % cv % size_jb 
   endif

   !--------------------------------------------------------------
   ! [2] Define flow-dependent control variable size:
   !--------------------------------------------------------------

   if ( be % ne > 0) then
      if (global) then
         be % cv % size_alphac = (be % alpha % max_wave + 1) * &
                                 (be % alpha % max_wave + 2)  / 2
         be % cv % size_je  = 2 * be % cv % size_alphac * be % alpha % mz
      else
         be % cv % size_alphac = ij * be % alpha % mz * be % ne
         be % cv % size_je = be % cv % size_alphac
      end if
   end if
   
   !--------------------------------------------------------------
   ! [3] Define domain-wide cv sizes for bit-repro option:
   !--------------------------------------------------------------

   cv_size_domain_jb = 0
   cv_size_domain_je = 0
   cv_size_domain_jl = 0

   if (.not. global) then
      iy = ide - ids + 1
      jx = jde - jds + 1
      if( use_rf )then
         ij = iy * jx
      else
         ij = size(be%wsd,1)*size(be%wsd,2)
      endif
      cv_size_domain_jb = ij * (be % v1 % mz + be % v2 % mz + be % v3 % mz + &
                                 be % v4 % mz + be % v5 % mz)
      cv_size_domain_je = ij * be % alpha % mz * be % ne
      if ( var4d ) then
         cv_size_domain_jl = cv_size_domain_jb 
      end if
   end if

   if (trace_use) call da_trace_exit("da_setup_cv")

end subroutine da_setup_cv


subroutine da_chgvres(nlath,nsig,kz,sigmah,sigma_avn,&
   corz_avn,cord_avn,corh_avn,corq_avn,hwll_avn,vztdq_avn,agv_avn,bv_avn,wgv_avn,&
   corz_kz, cord_kz, corh_kz, corq_kz, hwll_kz, vztdq_kz, agv_kz, bv_kz, wgv_kz)

   !---------------------------------------------------------------------------
   ! Purpose: Change vertical resolution of background stats for cv_options=3
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: nlath,nsig,kz
   real, intent(in)    ::  sigmah(kz),sigma_avn(1:nsig)

   real, intent(out) :: corz_kz(1:2*nlath+1,1:kz),cord_kz(1:2*nlath+1,1:kz) 
   real, intent(out) :: corh_kz(1:2*nlath+1,1:kz),corq_kz(1:2*nlath+1,1:kz)
   real, intent(out) :: hwll_kz(0:nlath*2+1,1:kz,1:4)
   real, intent(out) :: vztdq_kz(1:kz,0:nlath*2+1,1:4)
   real, intent(out) :: agv_kz(0:nlath*2+1,1:kz,1:kz)
   real, intent(out) :: bv_kz(0:nlath*2+1,1:kz),wgv_kz(0:nlath*2+1,1:kz)

   real, intent(in) :: corz_avn(1:2*nlath+1,1:nsig),cord_avn(1:2*nlath+1,1:nsig)
   real, intent(in) :: corh_avn(1:2*nlath+1,1:nsig),corq_avn(1:2*nlath+1,1:nsig)
   real, intent(in) :: hwll_avn(0:nlath*2+1,1:nsig,1:4)
   real, intent(in) :: vztdq_avn(1:nsig,0:nlath*2+1,1:4)
   real, intent(in) :: agv_avn(0:nlath*2+1,1:nsig,1:nsig)
   real, intent(in) :: bv_avn(0:nlath*2+1,1:nsig),wgv_avn(0:nlath*2+1,1:nsig)

   integer :: i,j,k,m,l,l1,m1,n
   real    :: rsigo(nsig),rsig(kz)
   real    :: coef1(kz),coef2(kz)
   integer :: lsig(kz)

   if (trace_use) call da_trace_entry("da_chgvres")

   if (kz==nsig) then
      corz_kz=corz_avn
      cord_kz=cord_avn
      corh_kz=corh_avn
      corq_kz=corq_avn
      hwll_kz=hwll_avn
      vztdq_kz=vztdq_avn
      agv_kz=agv_avn
      bv_kz=bv_avn
      wgv_kz=wgv_avn
      if (trace_use) call da_trace_exit("da_chgvres")
      return
   end if

   do k=1,kz
      rsig(k)=log(sigmah(k))
   end do
   do k=1,nsig
      rsigo(k)=log(sigma_avn(k))
   end do

   do k=1,kz
      if (rsig(k).ge.rsigo(1)) then
        m=1
        m1=2
        lsig(k)=1
        coef1(k)=1.0
           coef2(k)=0.0
      else if (rsig(k).lt.rsigo(nsig)) then
         m=nsig-1
         m1=nsig
         lsig(k)=nsig-1
         coef1(k)=0.0
         coef2(k)=1.0
      else
         do m=1,nsig
            m1=m+1
            if ((rsig(k).le.rsigo(m))   .and.  &
                (rsig(k).gt.rsigo(m1))    )then
               lsig(k)=m
               go to 2345
            end if
         end do
2345     continue
         coef1(k)=(rsigo(m1)-rsig(k))/(rsigo(m1)-rsigo(m))
         coef2(k)=1.0-coef1(k)
         if (lsig(k)==nsig) then
            lsig(k)=nsig-1
            coef2(k)=1.0
            coef1(k)=0.0
         end if
      end if
   end do

   ! agv wgv bv
   do k=1,kz
      m=lsig(k)
      m1=m+1
      do i=1,nlath*2
         wgv_kz(i,k)=wgv_avn(i,m)*coef1(k)+wgv_avn(i,m1)*coef2(k)
         bv_kz(i,k)=bv_avn(i,m)*coef1(k)+bv_avn(i,m1)*coef2(k)
      end do

      do j=1,kz
         l=lsig(j)
         l1=l+1
         do i=1,nlath*2
            agv_kz(i,j,k)=(agv_avn(i,l,m)*coef1(j)+agv_avn(i,l1,m)*coef2(j))*coef1(k) &
                    +(agv_avn(i,l,m1)*coef1(j)+agv_avn(i,l1,m1)*coef2(j))*coef2(k)
         end do
      end do
   end do

   agv_kz(0,:,:)=agv_kz(1,:,:)
   wgv_kz(0,:)=wgv_kz(1,:)
   bv_kz(0,:)=bv_kz(1,:)
   agv_kz(nlath*2+1,:,:)=agv_kz(nlath*2,:,:)
   wgv_kz(nlath*2+1,:)=wgv_kz(nlath*2,:)
   bv_kz(nlath*2+1,:)=bv_kz(nlath*2,:)

   do k=1,kz
      m=lsig(k)
      m1=m+1

      ! corz,cord,corh,corq
      do i=1,nlath*2
         corz_kz(i,k)=corz_avn(i,m)*coef1(k)+corz_avn(i,m1)*coef2(k)
         cord_kz(i,k)=cord_avn(i,m)*coef1(k)+cord_avn(i,m1)*coef2(k)
         corh_kz(i,k)=corh_avn(i,m)*coef1(k)+corh_avn(i,m1)*coef2(k)
         corq_kz(i,k)=corq_avn(i,m)*coef1(k)+corq_avn(i,m1)*coef2(k)
      end do

      do n=1,4 
         do i=1,nlath*2
            ! hwll
            hwll_kz(i,k,n)=hwll_avn(i,m,n)*coef1(k)+hwll_avn(i,m1,n)*coef2(k)
            ! vztdq
            vztdq_kz(k,i,n)=vztdq_avn(m,i,n)*coef1(k)+vztdq_avn(m1,i,n)*coef2(k)
          end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_chgvres")

end subroutine da_chgvres


subroutine da_setup_flow_predictors( ix, jy, kz, ne, ep, its, ite, jts, jte, kts, kte )

   !------------------------------------------------------------------------------
   ! Purpose: Setup structures for flow-dependent information and read it in.
   !------------------------------------------------------------------------------

   implicit none

   integer, intent(in)         :: ix, jy, kz            ! EP grid dimensions.
   integer, intent(in)         :: its, jts, kts         ! Tile start.
   integer, intent(in)         :: ite, jte, kte         ! Tile end.
   integer, intent(in)         :: ne                    ! Ensemble size.
   type (ep_type), intent(inout):: ep                   ! Flow-dependent info.

   character*10                :: ce                    ! Member index -> character.
   character(len=filename_len) :: filename              ! Input filename.
   character*10                :: var(1:5)              ! Variable name.
   integer                     :: ni, nj, nk            ! Grid dimensions.
   integer                     :: e                     ! Loop counter
   logical                     :: ldum1, ldum2,nkdum    ! Dummy.
   real                        :: temp3d(1:ix,1:jy,1:kz)! Temporary, real*4 array.
   real                        :: temp2d(1:ix,1:jy)     ! Temporary, real*4 array.

   real                        :: ens_scaling_inv       ! Ensemble scaling of perturbations.

   integer                     :: ep_unit

   if (trace_use) call da_trace_entry("da_setup_flow_predictors")

   call da_get_unit(ep_unit)

   call da_message(&
      (/"Setup structures for flow-dependent information and read in"/))

   ep % ne = ne

   ens_scaling_inv = 1.0
   if (ne > 1) ens_scaling_inv = 1.0 / sqrt(real(ne-1))

   ! Decide which space we are introducing flow-dependence:
   if (alphacv_method == alphacv_method_xa) then    ! xa space
      var(1) = 'u'
      var(2) = 'v'
      var(3) = 't'
      var(4) = 'q'
      var(5) = 'ps'
   else                               ! vp space (default)
      var(1) = 'psi'
      var(2) = 'chi_u'
      var(3) = 't_u'
      var(4) = 'rh'
      var(5) = 'ps_u'
   end if

   !---------------------------------------------------------------------------
   ! Input ensemble perturbations
   !---------------------------------------------------------------------------

   do e = 1, ne

      write(unit=ce,fmt='(i3.3)')e

      ! v1:
      filename = 'ep/'//trim(var(1))//'.e'//trim(ce)
      open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
      read(unit=ep_unit) ni, nj, nk

      if (ni /= ix .or. nj /= jy .or. nk /= kz) then
         write(unit=message(1),fmt='(a)') &
            'Inconsistent grid dimensions'
         write(unit=message(2),fmt='(a,3i6)') &
            ' Grid dims for analysis grid: ', ix, jy
         write(unit=message(3),fmt='(a,3i6)') &
            ' Grid dims for flow predictors: ', ni, nj
         call da_warning("da_setup_flow_predictors.inc",75,message(1:3))
      end if

      read(unit=ep_unit) temp3d(1:ix,1:jy,1:kz)
      close(unit=ep_unit)
      ep % v1(its:ite,jts:jte,kts:kte,e) = ens_scaling_inv * &
                                           temp3d(its:ite,jts:jte,kts:kte)

      ! v2:
      filename = 'ep/'//trim(var(2))//'.e'//trim(ce)
      open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
      read(unit=ep_unit) ni, nj, nk
      read(unit=ep_unit) temp3d(1:ix,1:jy,1:kz)
      ep % v2(its:ite,jts:jte,kts:kte,e) = ens_scaling_inv * &
                                           temp3d(its:ite,jts:jte,kts:kte)
      close(unit=ep_unit)

      ! v3:
      filename = 'ep/'//trim(var(3))//'.e'//trim(ce)
      open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
      read(unit=ep_unit) ni, nj, nk
      read(unit=ep_unit) temp3d(1:ix,1:jy,1:kz)
      ep % v3(its:ite,jts:jte,kts:kte,e) = ens_scaling_inv * &
                                           temp3d(its:ite,jts:jte,kts:kte)
      close(unit=ep_unit)

      ! v4:
      filename = 'ep/'//trim(var(4))//'.e'//trim(ce)
      open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
      read(unit=ep_unit) ni, nj, nk
      read(unit=ep_unit) temp3d(1:ix,1:jy,1:kz)
      ep % v4(its:ite,jts:jte,kts:kte,e) = ens_scaling_inv * &
                                           temp3d(its:ite,jts:jte,kts:kte)
      close(unit=ep_unit)

      ! v5:
      filename = 'ep/'//trim(var(5))//'.e'//trim(ce)
      open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
      read(unit=ep_unit) ni, nj, nkdum
      read(unit=ep_unit) temp2d(1:ix,1:jy)
      ep % v5(its:ite,jts:jte,1,e) = ens_scaling_inv * temp2d(its:ite,jts:jte)
      close(unit=ep_unit)

   end do 

!  Optional include hydrometeors:

   if( use_radarobs .and. use_radar_rf .or. use_rad .and. crtm_cloud & ! qt control variable
       .and. alphacv_method == alphacv_method_xa .and. alpha_hydrometeors ) then    ! xa space

      do e = 1, ne
         write(unit=ce,fmt='(i3.3)')e

         filename = 'ep/'//'qcloud'//'.e'//trim(ce)
         open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
         read(unit=ep_unit) ni, nj, nk
         read(unit=ep_unit) temp3d(1:ix,1:jy,1:kz)
         ep % v4(its:ite,jts:jte,kts:kte,e) = ep % v4(its:ite,jts:jte,kts:kte,e) + &
                                              ens_scaling_inv * temp3d(its:ite,jts:jte,kts:kte)
         close(unit=ep_unit)

         filename = 'ep/'//'qrain'//'.e'//trim(ce)
         open(unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
         read(unit=ep_unit) ni, nj, nk
         read(unit=ep_unit) temp3d(1:ix,1:jy,1:kz)
         ep % v4(its:ite,jts:jte,kts:kte,e) = ep % v4(its:ite,jts:jte,kts:kte,e) + &
                                              ens_scaling_inv * temp3d(its:ite,jts:jte,kts:kte)
         close(unit=ep_unit)
      end do
   end if

   call da_free_unit(ep_unit)

   if (trace_use) call da_trace_exit("da_setup_flow_predictors")

end subroutine da_setup_flow_predictors


subroutine da_setup_obs_structures( grid, ob, iv, j_cost)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate and read in components of observation structure.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain),     intent(inout) :: grid   ! Model data
   type (y_type),     intent(out)   :: ob     ! Observation structure.
   type (iv_type),    intent(out)   :: iv     ! O-B structure.
   type (j_type),     intent(inout) :: j_cost ! Cost function.

   integer :: i,j
   character (len=14)   :: ob_name

   if (trace_use) call da_trace_entry("da_setup_obs_structures")

! Default values for fmt_info, fmt_srfc, and fmt_each(YRG, 02/25/2009):

   fmt_info = '(a12,1x,a19,1x,a40,1x,i6,3(f12.3,11x),6x,a40)'
   fmt_srfc = '(7(:,f12.3,i4,f7.2))'
   fmt_each = '(3(f12.3,i4,f7.2),11x,3(f12.3,i4,f7.2),11x,1(f12.3,i4,f7.2))'

   call da_message((/'Set up observations (ob)'/))

   ! Adjust obs switches

   if (use_synopobs .OR. use_shipsobs .OR. use_metarobs .OR. use_pilotobs .OR. &
      use_profilerobs .OR. use_buoyobs .OR. use_soundobs .OR. use_mtgirsobs .OR. use_bogusobs .OR. &
      use_satemobs .OR. use_geoamvobs .OR. use_polaramvobs .OR. use_airepobs .OR. use_tamdarobs .OR. &
      use_gpspwobs .OR. use_gpsztdobs .OR. use_gpsrefobs .OR. use_ssmiretrievalobs .OR. &
      use_ssmitbobs .OR. use_ssmt1obs .OR. use_ssmt2obs .OR. use_qscatobs .OR. &
      use_airsretobs) then
 
      use_obsgts = .true.
   else
      use_obsgts = .false.
   end if

   if (use_hirs2obs .OR. use_hirs3obs .OR. use_msuobs .OR. use_amsuaobs .OR. &
      use_amsubobs .OR. use_airsobs .OR. use_eos_amsuaobs .OR. &
      use_hsbobs .OR. use_kma1dvar .OR. use_filtered_rad .OR. &
      use_ssmisobs .OR. use_hirs4obs .OR. use_mhsobs .OR. use_pseudo_rad .OR. &
      use_mwtsobs  .OR. use_mwhsobs .OR. use_atmsobs .OR. use_simulated_rad .OR. &
      use_iasiobs .OR. use_seviriobs .OR. use_amsr2obs) then
      use_rad = .true.
   else
      use_rad = .false.
   end if

   ! test_dm_exact can be set to .true. to force 1 runs 
   ! to produce results that are bitwise-identical regardless of the number of 
   ! MPI tasks used.  This is useful for validating that multi-processor runs 
   ! are not a source of bugs.  Runtime will be much longer.  This option is 
   ! automatically overridden to .false. for serial or 1-MPI-task runs.  

   if (test_dm_exact) then
      if (num_procs == 1) then
         test_dm_exact = .false.
         write(unit=stdout,fmt='(A)') &
            ' test_dm_exact overridden to .false. for serial or 1-MPI-task run'
      end if
      ! only implmenented for Sound and Synop, so switch other types off
      use_shipsobs         = .false.
      use_metarobs         = .false.
      use_bogusobs         = .false.
      use_pilotobs         = .false.
      use_airepobs         = .false.
      use_geoamvobs        = .false.
      use_polaramvobs      = .false.
      use_buoyobs          = .false.
      use_profilerobs      = .false.
      use_satemobs         = .false.
      use_gpspwobs         = .false.
      use_gpsztdobs        = .false.
      use_gpsrefobs        = .false.
      use_ssmiretrievalobs = .false.
      use_ssmitbobs        = .false.
      use_ssmt1obs         = .false.
      use_ssmt2obs         = .false.
      use_qscatobs         = .false.
      use_hirs2obs         = .false.
      use_hirs3obs         = .false.
      use_hirs4obs         = .false.
      use_mhsobs           = .false.
      use_msuobs           = .false.
      use_amsuaobs         = .false.
      use_amsubobs         = .false.
      use_mwtsobs          = .false.
      use_mwhsobs          = .false.
      use_atmsobs          = .false.
      use_airsobs          = .false.
      use_eos_amsuaobs     = .false.
      use_hsbobs           = .false.
      use_obsgts           = .false.
      use_rad              = .false.
      use_airsretobs       = .false.
      use_rainobs          = .false.
      use_iasiobs          = .false.	  
      use_radarobs         = .false.
      use_radar_rv         = .false.
      use_radar_rf         = .false.
   end if
    
   if (use_pseudo_rad .or. num_pseudo > 0) then
      call da_message((/"Single OBS Test:: Turn off all the OBS switches ***"/))
      use_synopobs         = .false.
      use_shipsobs         = .false.
      use_metarobs         = .false.
      use_soundobs         = .false.
      use_mtgirsobs        = .false.
      use_tamdarobs        = .false.
      use_bogusobs         = .false.
      use_pilotobs         = .false.
      use_airepobs         = .false.
      use_geoamvobs        = .false.
      use_polaramvobs      = .false.
      use_buoyobs          = .false.
      use_profilerobs      = .false.
      use_satemobs         = .false.
      use_gpspwobs         = .false.
      use_gpsztdobs        = .false.
      use_gpsrefobs        = .false.
      use_ssmiretrievalobs = .false.
      use_ssmitbobs        = .false.
      use_ssmt1obs         = .false.
      use_ssmt2obs         = .false.
      use_qscatobs         = .false.
      use_hirs2obs         = .false.
      use_hirs3obs         = .false.
      use_hirs4obs         = .false.
      use_mhsobs           = .false.
      use_msuobs           = .false.
      use_amsuaobs         = .false.
      use_amsubobs         = .false.
      use_mwtsobs          = .false.
      use_mwhsobs          = .false.
      use_atmsobs          = .false.
      use_airsobs          = .false.
      use_eos_amsuaobs     = .false.
      use_hsbobs           = .false.
      use_obsgts           = .true.
      use_rad              = .false.
      use_airsretobs       = .false.
      use_rainobs          = .false.
      use_iasiobs          = .false.	  
      use_radarobs         = .false.
      use_radar_rv         = .false.
      use_radar_rf         = .false.
      ob_format = ob_format_ascii
   end if

   if ( use_pseudo_rad ) then
      use_rad = .true.
      thinning         = .false.
      qc_rad           = .false.
      rtminit_nsensor  = 1
      rtminit_platform = pseudo_rad_platid
      rtminit_satid    = pseudo_rad_satid
      rtminit_sensor   = pseudo_rad_senid
   end if

   if (sfc_assi_options < 1 .OR. sfc_assi_options > 2) then
      write(unit=message(1),fmt='(A,I3)') &
         'Invalid sfc_assi_option = ', sfc_assi_options
      call da_error("da_setup_obs_structures.inc",167,message(1:1))
   end if

   if ( use_gpsrefobs ) then
      if ( ob_format == ob_format_bufr ) then
          ! when main conv input is from bufr, gpsro has to be from bufr too
          ! overwrite the namelist setting
          ob_format_gpsro = ob_format_bufr
      else 
         if ( ob_format == ob_format_ascii .and. &
              ob_format_gpsro == ob_format_bufr ) then
            call da_message((/'ob_format=2 while ob_format_gpsro=1, will be reading GPSRO data from gpsro.bufr'/))
         end if
      end if
   end if

   !---------------------------------------------------------------------------      
   ! [1.0] Setup and read in fields from first guess:
   !----------------------------------------------------------------------------     

   iv%missing = missing
   ! iv%ptop    = grid%xb%ptop

   iv%total_rad_pixel   = 0
   iv%total_rad_channel = 0

   iv%info(:)%nlocal = 0
   iv%info(:)%ntotal = 0
   iv%info(:)%thin_nlocal = 0
   iv%info(:)%thin_ntotal = 0
   do i=1,num_ob_indexes
      iv%info(i)%plocal(:) = 0
      iv%info(i)%ptotal(:) = 0
      iv%info(i)%thin_plocal(:) = 0
      iv%info(i)%thin_ptotal(:) = 0
   end do
   iv%num_inst  = 0 

   ob%nlocal(:) = 0
   ob%ntotal(:) = 0
   ob%num_inst  = 0

   iv%info(:)%max_lev = 1

   ! get FGAT time slots

   allocate ( time_slots(0:num_fgat_time) )
   call da_get_time_slots(num_fgat_time,time_window_min,analysis_date, &
                          time_window_max, time_slots)

   if (use_obsgts) then
      ! Conventional obs can be in 1 or ascii format
      if (ob_format == ob_format_bufr) then
         call da_message((/'Using PREPBUFR format observation input'/))
         call da_setup_obs_structures_bufr (grid, ob, iv)
      else if (ob_format == ob_format_ascii) then
         call da_message((/'Using ASCII format observation input'/))
         call da_setup_obs_structures_ascii (ob, iv, grid)
      else
         write(unit=message(1),fmt='(A,I3)') 'Invalid ob_format = ', ob_format
         call da_error("da_setup_obs_structures.inc",227,message(1:1))
      end if
   end if


   if (use_radarobs) then
      ! Radar obs are read from separate file(s)
         call da_message((/'Using ASCII format radar observation input'/))
         call da_setup_obs_structures_radar (grid, ob, iv)
   end if

   if (use_rainobs .and. var4d) then
      call da_message((/'Using ASCII format precipitation observation input'/))
      call da_setup_obs_structures_rain (grid, ob, iv) 
   end if

   if (use_rad) then
      ! Radiance files can only be in 1
      call da_message((/'Using NCEP BUFR radiance 1b input'/))
      if ( use_amsr2obs ) then
         call da_message((/'Using AMSR2 radiance input in HDF5 format'/))
      end if
      call da_setup_radiance_structures(grid, ob, iv)
   end if

   ! Summarize observations 

   write(unit=stdout, fmt='(a)')  'Observation summary'
   do i=1,num_fgat_time
      write(unit=stdout, fmt='(3x,a,i2)') 'ob time ', i
      do j=1,num_ob_indexes
         if(j.eq.27)cycle 
         if (use_gpsztdobs .and.obs_names(j) == 'gpspw         ' ) then
            ob_name = 'gpsztd        '
         else
            ob_name = obs_names(j)
         endif
         if(iv%info(j)%ptotal(i) - iv%info(j)%ptotal(i-1) > 0) then          
         if ( thin_conv_ascii .and. ob_format == ob_format_ascii ) then
            if ( j == radiance ) cycle !temporary fix to not print incorrect counts for radiance
            write(unit=stdout, fmt='(2x,2(a,i10),a,i8,a)') &
               ob_name, iv%info(j)%ptotal(i) - iv%info(j)%ptotal(i-1), ' global, ', &
               iv%info(j)%thin_ptotal(i) - iv%info(j)%thin_ptotal(i-1), ' global (post-thinning), ', &
               iv%info(j)%thin_plocal(i) - iv%info(j)%thin_plocal(i-1), ' local (post-thinning)'
         else
            write(unit=stdout, fmt='(6x,a,i10,a,i8,a)') &
               ob_name, iv%info(j)%ptotal(i) - iv%info(j)%ptotal(i-1), ' global,', &
               iv%info(j)%plocal(i) - iv%info(j)%plocal(i-1), ' local'
         end if ! thin_conv_ascii
         end if
      end do
   end do
   write(unit=stdout, fmt='(a)') ' '
  
   ! Get horizontal interpolation weights.

   call da_setup_obs_interp_wts (iv) 

   if (trace_use) call da_trace_exit("da_setup_obs_structures")    

end subroutine da_setup_obs_structures


subroutine da_setup_obs_structures_ascii( ob, iv, grid )

   !-------------------------------------------------------------------------
   ! Purpose: Define, allocate and read of observation structure.
   ! Updates:
   !          Syed RH Rizvi NCAR/NESL/MMM/DAS Date:  02/21/2013 
   !          Updated with thinning option
   !-------------------------------------------------------------------------

   implicit none

   type (y_type),  intent(out)   :: ob  ! Observation structure.
   type (iv_type), intent(inout) :: iv  ! O-B structure.
   type (domain),  intent(inout) :: grid  ! First guess structure

   character(len=filename_len)  :: filename
   integer                      :: n, i, j, k
   logical                      :: outside, thin_3d
   logical                      :: uvq_direct=.false.
   integer  :: istart,iend,jstart,jend
   real     :: rlonlat(4)

   if (trace_use) call da_trace_entry("da_setup_obs_structures_ascii")
   !-------------------------------
   ! 0.0 Make thinning grids
   !------------------------------
   call init_constants_derived
   thin_3d=.false.
   if ( thin_conv_ascii ) then
      rlat_min =  r999
      rlat_max = -r999
      rlon_min =  r999
      rlon_max = -r999

      istart=MINVAL( grid%i_start(1:grid%num_tiles) )
      iend  =MAXVAL( grid%i_end  (1:grid%num_tiles) )
      jstart=MINVAL( grid%j_start(1:grid%num_tiles) )
      jend  =MAXVAL( grid%j_end  (1:grid%num_tiles) )

      do i = istart, iend
         do j = jstart, jend
            rlat_min=min(rlat_min, grid%xb%lat(i,j))
            rlat_max=max(rlat_max, grid%xb%lat(i,j))
            if( grid%xb%lon(i,j) < 0.0) then
              rlon_min=min(rlon_min, (r360+grid%xb%lon(i,j)))
              rlon_max=max(rlon_max, (r360+grid%xb%lon(i,j)))
            else
              rlon_min=min(rlon_min, grid%xb%lon(i,j))
              rlon_max=max(rlon_max, grid%xb%lon(i,j))
            endif
         enddo
      enddo

      call mpi_reduce(rlat_min, rlonlat(1), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlon_min, rlonlat(2), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlat_max, rlonlat(3), 1, true_mpi_real, mpi_max, root, comm, ierr)
      call mpi_reduce(rlon_max, rlonlat(4), 1, true_mpi_real, mpi_max, root, comm, ierr)

      CALL mpi_bcast( rlonlat, 4 , true_mpi_real , root , comm, ierr )

      rlat_min = rlonlat(1)
      rlon_min = rlonlat(2)
      rlat_max = rlonlat(3)
      rlon_max = rlonlat(4)

      dlat_grid = rlat_max - rlat_min
      dlon_grid = rlon_max - rlon_min

      allocate(thinning_grid_conv(num_ob_indexes))
      do n = 1, num_ob_indexes
         if (n == radar) cycle
         if( n == airep .or. n == tamdar ) then
            thin_3d=.true.
            call make3grids (n,thin_mesh_conv(n), thin_3d)
         else
            call make3grids (n,thin_mesh_conv(n))
         end if 
      end do
   end if


   !-------------------------------


   ! Find out if pseudo ob is local
   !-------------------------------
   if (num_pseudo > 0) then
      iv%time = 1
      outside = .false.
      i = int(pseudo_x)
      j = int(pseudo_y)
      if (fg_format == fg_format_kma_global) then
         if ((j < jts-1) .or. (j > jte)) outside = .true.
      else
         if ((i < ids)   .or. (i >= ide) .or. (j < jds)   .or. (j >= jde)) outside = .true.
         if ((i < its-1) .or. (i >  ite) .or. (j < jts-1) .or. (j >  jte)) outside = .true.
      end if
   else
      iv%time = 1
      outside = .true.
   end if

   if (num_pseudo > 0) then
   if (outside) then
      iv%info(pseudo)%nlocal          = 0
      ob%nlocal(pseudo)               = 0
      iv%info(pseudo)%ntotal          = 0
      iv%info(pseudo)%plocal(iv%time) = 0
   else
      iv%info(pseudo)%nlocal          = num_pseudo
      ob%nlocal(pseudo)               = num_pseudo
      iv%info(pseudo)%ntotal          = num_pseudo
      iv%info(pseudo)%max_lev         = 1            
   end if
   end if

   !--------------------------------------------------------------------------
   ! [1.0] Scan GTS observation header and get idea of number of obs:
   !--------------------------------------------------------------------------
  
   if (num_fgat_time > 1) then
!      filename = ' '

      do n=1, num_fgat_time

         iv%time = n
         filename = ' '

        write(filename(1:10), fmt='(a, i2.2, a)') 'ob', n,'.ascii'

         ! scan main body of gts observation file
         call da_scan_obs_ascii (iv, filename,grid)

         if (use_ssmiretrievalobs .or. use_ssmitbobs) then
            ! scan SSMI observation file
            write(filename(1:9), fmt='(a, i2.2, a)') 'ob', n,'.ssmi'
            call da_scan_obs_ssmi (iv, filename)
         end if

         iv%info(:)%plocal(n) = iv%info(:)%nlocal
         iv%info(:)%ptotal(n) = iv%info(:)%ntotal
      end do
   else
      iv%time = 1
      call da_scan_obs_ascii(iv, 'ob.ascii', grid)
      !-----------------------------------------------------------------------
      ! read header of ssmi observation file
      !-----------------------------------------------------------------------
      if (use_ssmiretrievalobs .or. use_ssmitbobs) then
         call da_scan_obs_ssmi(iv, 'ob.ssmi')
      end if

      do i=1,num_ob_indexes
         if (i == radar) cycle
         iv%info(i)%plocal(iv%time) = iv%info(i)%nlocal
         iv%info(i)%ptotal(iv%time) = iv%info(i)%ntotal
      end do
   end if

   !--------------------------------------------------------------------------
   ! Allocate the ob based on input number of obs:
   !--------------------------------------------------------------------------
   call da_allocate_observations (iv)

   if (num_fgat_time > 1) then

      do n=1, num_fgat_time
         iv%time = n
         filename = ' '  

         write(filename(1:10), fmt='(a, i2.2, a)') 'ob', n,'.ascii'

         ! Read gts observation file
         call da_read_obs_ascii (iv, filename, uvq_direct, grid)

         if (use_ssmiretrievalobs .or. use_ssmitbobs) then
            ! read ssmi observation file
            write(filename(1:9), fmt='(a, i2.2, a)') 'ob', n,'.ssmi'
            call da_read_obs_ssmi (iv, filename)
         end if

         do i=1,num_ob_indexes
            if (i == radar) cycle
            iv%info(i)%thin_ptotal(n) = iv%info(i)%thin_ntotal
            iv%info(i)%thin_plocal(n) = iv%info(i)%thin_nlocal
         end do
      end do
   else
      iv%time = 1

      call da_read_obs_ascii(iv, 'ob.ascii', uvq_direct, grid)

      if (use_ssmiretrievalobs .or. use_ssmitbobs) then
         ! read ssmi observation file
         call da_read_obs_ssmi (iv, 'ob.ssmi')
      end if

      do i=1,num_ob_indexes
         if (i == radar) cycle
         iv%info(i)%thin_ptotal(iv%time) = iv%info(i)%thin_ntotal
         iv%info(i)%thin_plocal(iv%time) = iv%info(i)%thin_nlocal
      end do
   end if

   if ( use_gpsrefobs .and. (ob_format_gpsro == ob_format_bufr) ) then
      call da_read_obs_bufrgpsro(iv)
   end if

   !--------------------------------------------------------------------------
   ! [2.5] Set all thinned obs missing  
   !--------------------------------------------------------------------------
    if ( thin_conv_ascii ) then
       do i = 1, num_ob_indexes
          if (i == radar) cycle
          if ( iv%info(i)%ntotal > 0 ) then
             if ( iv%info(i)%nlocal > 0 ) then
                if ( ANY(iv%info(i)%thinned(:,:)) ) then
                   if( i == airep .or. i==tamdar ) then
                      call da_set_3d_obs_missing(iv,i)  ! assign missing values level-by-level if thinned=true data
                   else
                      call da_set_obs_missing(iv,i)  ! assign missing values to those thinned=true data
                   end if
                end if
             end if
          end if
       end do
    end if ! thin_conv_ascii
  
   !--------------------------------------------------------------------------
   ! [3.0] Calculate innovation vector (O-B) and create (smaller) ob structure:
   !--------------------------------------------------------------------------

   if (uvq_direct) then
      call da_fill_obs_structures(iv, ob, uvq_direct)
   else
      call da_fill_obs_structures(iv, ob)
   endif



   iv%time = 1

   if ( thin_conv_ascii ) then
      do n = 1, num_ob_indexes
         if (n == radar) cycle
         if( n == airep .or. n==tamdar ) then
            thin_3d=.true.
            call destroygrids_conv (n, thin_3d)
         else
            call destroygrids_conv (n)
         end if
      end do
      deallocate(thinning_grid_conv)
   end if

   if (trace_use) call da_trace_exit("da_setup_obs_structures_ascii")
end subroutine da_setup_obs_structures_ascii


subroutine da_setup_obs_structures_bufr(grid, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: Define, allocate and read observation structure.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain) , intent(in)    :: grid        ! model data
   type (y_type),  intent(out)   :: ob          ! Observation structure.
   type (iv_type), intent(inout) :: iv          ! O-B structure.


   character(len=filename_len) :: filename
   integer                     :: n,i,j
   
   ! thinning variables
   integer  :: istart,iend,jstart,jend
   real     :: rlonlat(4)

   if (trace_use) call da_trace_entry("da_setup_obs_structures_bufr")

   !-------------------------------
   ! 0.0 Make thinning grids
   !------------------------------
   call init_constants_derived

   if ( thin_conv ) then
      rlat_min =  r999
      rlat_max = -r999
      rlon_min =  r999
      rlon_max = -r999

      istart=MINVAL( grid%i_start(1:grid%num_tiles) )
      iend  =MAXVAL( grid%i_end  (1:grid%num_tiles) )
      jstart=MINVAL( grid%j_start(1:grid%num_tiles) )
      jend  =MAXVAL( grid%j_end  (1:grid%num_tiles) ) 

      do i = istart, iend
         do j = jstart, jend
            rlat_min=min(rlat_min, grid%xb%lat(i,j))
            rlat_max=max(rlat_max, grid%xb%lat(i,j))
            if( grid%xb%lon(i,j) < 0.0) then
              rlon_min=min(rlon_min, (r360+grid%xb%lon(i,j)))
              rlon_max=max(rlon_max, (r360+grid%xb%lon(i,j)))
            else
              rlon_min=min(rlon_min, grid%xb%lon(i,j))
              rlon_max=max(rlon_max, grid%xb%lon(i,j))
            endif
         enddo
      enddo

      call mpi_reduce(rlat_min, rlonlat(1), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlon_min, rlonlat(2), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlat_max, rlonlat(3), 1, true_mpi_real, mpi_max, root, comm, ierr)
      call mpi_reduce(rlon_max, rlonlat(4), 1, true_mpi_real, mpi_max, root, comm, ierr)

      CALL mpi_bcast( rlonlat, 4 , true_mpi_real , root , comm, ierr )

      rlat_min = rlonlat(1)
      rlon_min = rlonlat(2)
      rlat_max = rlonlat(3)
      rlon_max = rlonlat(4)

      dlat_grid = rlat_max - rlat_min
      dlon_grid = rlon_max - rlon_min

      allocate(thinning_grid_conv(num_ob_indexes))
      do n = 1, num_ob_indexes
         call make3grids (n,thin_mesh_conv(n))
      end do
   end if

   !--------------------------------------------------------------------------
   ! [1.0] Read data
   !--------------------------------------------------------------------------
    
      iv%time = 1
      call da_read_obs_bufr(iv)
! 
!for gps

   if ( use_gpsrefobs ) then
      call da_read_obs_bufrgpsro(iv)
   end if

   if ( thin_conv ) then
      do n = 1, num_ob_indexes
         call destroygrids_conv (n)
      end do
      deallocate(thinning_grid_conv)
   end if

   !--------------------------------------------------------------------------
   ! [3.0] Calculate innovation vector (O-B) and create (smaller) ob structure:
   !--------------------------------------------------------------------------

   call da_fill_obs_structures(iv, ob)

   iv%time = 1


   if (trace_use) call da_trace_exit("da_setup_obs_structures_bufr")

end subroutine da_setup_obs_structures_bufr


subroutine da_setup_obs_structures_madis( ob, iv)

!------------------------------------------------------------------------------
! PURPOSE: Define, allocate and read of observation structure.
!
! METHOD:  Define, allocate and read of observation structure.
!
! HISTORY: 05-21-03  M. Barth         The setup_obs routine to ingest data from
!                                     the FSL MADIS netCDF files.
!          10-20-03  M. Barth         - Changed calls to get SATOB data for the
!                                       changes made in the MADIS API when the
!                                       true SATWND dataset was added.  (This code
!                                       originally accessed a dummy SATW dataset.)
!                                     - Changed max_map from 150 to 2000 to 
!                                       accommodate full time window options.
!          04-01-04  M. Barth         Added MADIS COOP surface dataset.
!          08-27-04  M. Barth         - Added sonde_sfc stuff.
!                                     - Initialize new count_obs entries that we
!                                       don't use (radar, profiler, buoy).  Note
!                                       that, as with the 1 interface, we
!                                       continue to lump profiler in with "sound"
!                                       and buoy with "ships" -- as the downstream
!                                       processing is identical regardless of
!                                       the designation, and we offer MADIS-level
!                                       control of these datasets, this is simpler
!                                       than restructuring this code.
!          02-06-06 M. Barth          Changes for WRFVAR version 2.1.
!                                     - Renamed all "fsl" references to "madis".
!                                       (Exception:  "FSL" is still used as an
!                                        argument to MINIT.)
!                                     - Changed xb to xp.  
!                                     - Removed xbx arg, so analysis time is
!                                       now determined by the NAMELIST ANALYSIS_DATE
!                                       variable.
!                                     - Changed satob to geoamv.
!                                     - Changed max_sfc from 50K to 200K.
!                                     - Changed max_raoblev from 256 to 405.
!                                     - Miscellaneous other minor changes to
!                                       conform to current paradigm.
!                                     - Added protection against invalid lat/lons,
!                                       as well as protection against latitude
!                                       at the poles -- some projections can't
!                                       handle this, and the WRFVAR routines
!                                       that do the coordinate conversions aren't
!                                       smart enough to trap this case.  [Are you
!                                       surprised?]
!          08-07-07 M. Barth          - Fixed problem introduced in 2005 version
!                                       of WRFVAR by filling in the ob_num 
!                                       and current_ob_time data structures in iv.
!          02-06-06 M. Barth          Changes for WRFVAR version 3:
!                                     - Removed xp argument.
!                                     - Changed DA_ll_to_xy call to da_llxy.
!                                     - Added trace usage.
!                                     - Changed ob_type structure to iv_type.
!                                     - Changes in iv structure:
!                                       - current_ob_time -> time
!                                       - New info structure
!                                     - Fill in global vs. local obs counts
!                                       into iv%info.  This is *presumably* the
!                                       replacement for iv%ob_num.
!                                     - Fill in iv%info with stuff that used
!                                       to go into iv%<dataset>%info.
!                                     - Since DA_Constants and DA_Read_Namelist no
!                                       longer exist, the declaration of, and
!                                       input of, the MADIS namelist variables
!                                       have been moved to this routine.
!                                     - Changed print_detail to print_detail_obs.
!                                     - Added duplication logic (see explanation
!                                       in PROCESS_DATASET).
!          11-28-08 M. Barth          Since I finally got the darn thing to
!                                     compile, I debugged the massive changes
!                                     in version 3.
!
! PARENT_MODULE: da_setup_structures
!------------------------------------------------------------------------------

   use da_tools, only : da_llxy
   IMPLICIT NONE
   
! ARGUMENTS:

   TYPE ( y_type), INTENT(OUT)    :: ob          ! Observation structure.
   TYPE (iv_type), INTENT(INOUT)  :: iv          ! O-B structure.

   call da_error("da_setup_obs_structures_madis.inc",2106,(/"Needs to be compiled with a MADIS library"/))

   if (trace_use) CALL da_trace_exit("da_setup_obs_structures_madis")

END SUBROUTINE da_setup_obs_structures_madis

subroutine da_setup_obs_structures_rain(grid, ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: Define, allocate and read observation structure.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain) , intent(in)    :: grid        ! model data
   type (y_type),  intent(out)   :: ob          ! Observation structure.
   type (iv_type), intent(inout) :: iv          ! O-B structure.

   character(len=filename_len)  :: filename
   integer                     :: n,i,j,k
   
   ! thinning variables
   integer  :: istart,iend,jstart,jend
   real     :: rlonlat(4)

   if (trace_use) call da_trace_entry("da_setup_obs_structures_rain")

   !-------------------------------
   ! 0.0 Make thinning grids
   !------------------------------
   call init_constants_derived

   if ( thin_rainobs ) then
      rlat_min =  r999
      rlat_max = -r999
      rlon_min =  r999
      rlon_max = -r999

      istart=MINVAL( grid%i_start(1:grid%num_tiles) )
      iend  =MAXVAL( grid%i_end  (1:grid%num_tiles) )
      jstart=MINVAL( grid%j_start(1:grid%num_tiles) )
      jend  =MAXVAL( grid%j_end  (1:grid%num_tiles) ) 

      do i = istart, iend
         do j = jstart, jend
            rlat_min=min(rlat_min, grid%xb%lat(i,j))
            rlat_max=max(rlat_max, grid%xb%lat(i,j))
            if( grid%xb%lon(i,j) < 0.0) then
              rlon_min=min(rlon_min, (r360+grid%xb%lon(i,j)))
              rlon_max=max(rlon_max, (r360+grid%xb%lon(i,j)))
            else
              rlon_min=min(rlon_min, grid%xb%lon(i,j))
              rlon_max=max(rlon_max, grid%xb%lon(i,j))
            endif
         enddo
      enddo

      call mpi_reduce(rlat_min, rlonlat(1), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlon_min, rlonlat(2), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlat_max, rlonlat(3), 1, true_mpi_real, mpi_max, root, comm, ierr)
      call mpi_reduce(rlon_max, rlonlat(4), 1, true_mpi_real, mpi_max, root, comm, ierr)

      CALL mpi_bcast( rlonlat, 4 , true_mpi_real , root , comm, ierr )

      rlat_min = rlonlat(1)
      rlon_min = rlonlat(2)
      rlat_max = rlonlat(3)
      rlon_max = rlonlat(4)

      dlat_grid = rlat_max - rlat_min
      dlon_grid = rlon_max - rlon_min

      allocate(thinning_grid(num_ob_indexes,num_fgat_time))
      do n=1, num_fgat_time
         call makegrids (rain,thin_mesh_conv(rain), n)
      end do
   end if

   !--------------------------------------------------------------------------
   ! [1.0] Scan data
   !--------------------------------------------------------------------------

   do n=2, num_fgat_time
      if ( .not. fgat_rain_flags(n) ) cycle
      iv%time = n
      filename = ' '

      ! scan rainfall observation file
      write(filename(1:9), fmt='(a, i2.2, a)') 'ob', n,'.rain'
      call da_scan_obs_rain(iv, filename, iv%time)

      iv%info(rain)%plocal(n) = iv%info(rain)%nlocal
      iv%info(rain)%ptotal(n) = iv%info(rain)%ntotal
   end do

   !--------------------------------------------------------------------------
   ! Allocate the ob based on input number of obs:
   !--------------------------------------------------------------------------

   call da_allocate_observations_rain (iv)

   !--------------------------------------------------------------------------
   ! [2.0] Read data
   !--------------------------------------------------------------------------
    
   do n=2, num_fgat_time
      if ( .not. fgat_rain_flags(n) ) cycle
      iv%time = n
      filename = ' '

      ! read rainfall observation file
      write(filename(1:9), fmt='(a, i2.2, a)') 'ob', n,'.rain'
     call da_read_obs_rain(iv, filename, iv%time)
   end do

   if ( thin_rainobs ) then
     do n = 1, num_fgat_time
        call destroygrids (rain,n)
     end do
     deallocate(thinning_grid)
   end if

   !--------------------------------------------------------------------------
   ! [3.0] Calculate innovation vector (O-B) and create (smaller) ob structure:
   !--------------------------------------------------------------------------

   call da_fill_obs_structures_rain(iv, ob)

   iv%time = 1

   if (trace_use) call da_trace_exit("da_setup_obs_structures_rain")

end subroutine da_setup_obs_structures_rain


subroutine da_setup_obs_structures_radar( grid, ob, iv )

   !-------------------------------------------------------------------------
   ! Purpose: Define, allocate and read radar observation structure.
   !-------------------------------------------------------------------------

   implicit none

   type (y_type),  intent(out)   :: ob  ! Observation structure.
   type (iv_type), intent(inout) :: iv  ! O-B structure.
   type (domain),  intent(inout) :: grid  ! First guess structure

   character(len=filename_len)  :: filename
   integer                      :: n, i, j, k
   integer  :: istart,iend,jstart,jend
   real     :: rlonlat(4)

   if (trace_use) call da_trace_entry("da_setup_obs_structures_radar")

   call init_constants_derived

   !--------------------------------------------------------------------------
   ! [1.0] Scan RADAR observation header and get number of obs:
   !--------------------------------------------------------------------------
   if (num_fgat_time > 1) then
      do n=1, num_fgat_time

         iv%time = n
         filename = ' '

         ! scan radar observation file
         write(filename(1:10), fmt='(a, i2.2, a)') 'ob', n,'.radar'
         call da_scan_obs_radar(iv, filename, grid)

         iv%info(radar)%plocal(n) = iv%info(radar)%nlocal
         iv%info(radar)%ptotal(n) = iv%info(radar)%ntotal
      end do
   else
      iv%time = 1
      ! scan main body of radar observation file
      call da_scan_obs_radar(iv, 'ob.radar', grid)
      iv%info(radar)%plocal(iv%time) = iv%info(radar)%nlocal
      iv%info(radar)%ptotal(iv%time) = iv%info(radar)%ntotal
   end if

   !--------------------------------------------------------------------------
   ! Allocate based on input number of obs:
   !--------------------------------------------------------------------------
   ! This logic was originally found in da_allocate_observations; moved here
   if (iv%info(radar)%nlocal > 0) allocate(iv%radar (1:iv%info(radar)%nlocal))
   if (iv%info(radar)%nlocal > 0) then
      allocate (iv%info(radar)%name(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%platform(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%id(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%date_char(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%levels(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%lat(iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%lon(iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%elv(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%pstar(iv%info(radar)%nlocal))

      allocate (iv%info(radar)%slp(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%pw(iv%info(radar)%nlocal))

      allocate (iv%info(radar)%x  (kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%y  (kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%i  (kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%j  (kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%dx (kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%dxm(kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%dy (kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%dym(kms:kme,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%k  (iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%dz (iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%dzm(iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%zk (iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%proc_domain(iv%info(radar)%max_lev,iv%info(radar)%nlocal))
      allocate (iv%info(radar)%obs_global_index(iv%info(radar)%nlocal))
      allocate (iv%info(radar)%thinned(iv%info(radar)%max_lev,iv%info(radar)%nlocal))

      iv%info(radar)%proc_domain(:,:)  = .false.
      iv%info(radar)%thinned(:,:)      = .false.
      iv%info(radar)%zk(:,:)           = missing_r
   end if

   if (num_fgat_time > 1) then

      do n=1, num_fgat_time
         iv%time = n
         filename = ' '  

         ! read radar observation file
         write(filename(1:10), fmt='(a, i2.2, a)') 'ob', n,'.radar'
         call da_read_obs_radar(iv, filename, grid)

      end do
   else
      iv%time = 1

      ! read radar observation file
      call da_read_obs_radar(iv, 'ob.radar', grid)
   end if


   ! Calculate DT for RF DA

   if (use_radar_rf) then
      if (.not. DT_cloud_model) then
         do j = jts,jte
            do i = its, ite
               do k = kts, kte
                   grid%xb%delt(i,j,k) = 0.0
               end do
            end do
         end do

         do n = 1, iv%info(radar)%nlocal
            do i=int(iv%info(radar)%i(1,n)), int(iv%info(radar)%i(1,n))+1
               do j=int(iv%info(radar)%j(1,n)), int(iv%info(radar)%j(1,n))+1
                  do k=kts, kte
                     grid%xb%delt(i,j,k) = 1800.0
                     grid%xb%qrn(i,j,k) = amax1(5.0E-5, grid%xb%qrn(i,j,k))
                     grid%xb%qcw(i,j,k) = amax1(5.0E-12, grid%xb%qcw(i,j,k))
                  end do
               end do
            end do
         end do
      end if
   end if

  
   !--------------------------------------------------------------------------
   ! [3.0] Calculate innovation vector (O-B) and create (smaller) ob structure:
   !--------------------------------------------------------------------------

   call da_fill_obs_structures_radar(iv, ob)

   iv%time = 1

   if (trace_use) call da_trace_exit("da_setup_obs_structures_radar")
end subroutine da_setup_obs_structures_radar

subroutine da_setup_obs_interp_wts (iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout)  :: iv         ! Innovation vector (O-B).

   integer                        :: i         ! Loop counter.

   if (trace_use) call da_trace_entry("da_setup_obs_interp_wts")

   do i=1,num_ob_indexes
      if (i /= radiance .and. iv%info(i)%nlocal > 0) then ! i=22 is radiance , should be excluded here
         if ( ob_format == ob_format_ascii ) then
            call da_store_obs_grid_info (iv%info(i))
         else if ( ob_format == ob_format_bufr ) then
            call da_store_obs_grid_info (iv%info(i))
         end if
      end if
   end do

   do i = 1, iv % num_inst
      if (iv % instid(i) % num_rad < 1) cycle
      call da_store_obs_grid_info_rad (iv%instid(i)%info)
   end do

   if (trace_use) call da_trace_exit("da_setup_obs_interp_wts")

end subroutine da_setup_obs_interp_wts


subroutine da_setup_runconstants(grid, xbx)

   !---------------------------------------------------------------------------
   !  Purpose: Define constants used.
   !---------------------------------------------------------------------------

   implicit none

   type (domain),   intent(inout) :: grid                
   type (xbx_type), intent(inout) :: xbx     ! Header & non-gridded vars.

   integer                  :: n             ! Loop counter.

   integer                  :: fft_size_i    ! Efficient FFT 1st dimension.
   integer                  :: fft_size_j    ! Efficient FFT 2nd dimension.
   logical                  :: found_magic   ! true when efficient FFT found
   logical                  :: need_pad      ! True when need pad_i > 0

   integer                  :: i, j
   integer                  :: nx            ! nx + 1 = ix + pad_i
   integer                  :: ny            ! ny + 1 = jy + pad_j
   real                     :: const         ! Multiplicative constants.
   real                     :: coeff_nx      ! Multiplicative coefficients.
   real                     :: coeff_ny      ! Multiplicative coefficients.
   real                     :: cos_coeff_nx  ! Multiplicative coefficients.
   real                     :: cos_coeff_ny  ! Multiplicative coefficients.

   if (trace_use) call da_trace_entry("da_setup_runconstants")

   !---------------------------------------------------------------------------
   ! [1.0] Calculate padded grid-size required for efficient FFTs and factors:
   !---------------------------------------------------------------------------

   ! [1.1] In x direction:

   nx = ide - ids

   allocate(xbx % fft_factors_x(num_fft_factors))

   do n = nx, nx+nrange
      call da_find_fft_factors(n, found_magic, xbx % fft_factors_x)

      if (found_magic) then
         if (mod(n, 2) == 0) then
            fft_size_i = n
            xbx % fft_pad_i = n - nx
            exit
         end if
      end if
   end do

   if (.NOT. found_magic) then
     call da_error("da_setup_runconstants.inc",53, &
       (/"No FFT factor found: invalid e_we value"/))
   end if

   allocate(xbx % trig_functs_x(1:3*fft_size_i))  

   call da_find_fft_trig_funcs(fft_size_i, xbx % trig_functs_x(:))

   ! [1.2] In y direction:

   ny = jde - jds

   allocate(xbx % fft_factors_y(num_fft_factors))

   do n = ny, ny+nrange
      call da_find_fft_factors(n, found_magic, xbx % fft_factors_y)

      if (found_magic) then
         if (mod(n, 2) == 0) then
            fft_size_j = n
            xbx % fft_pad_j = n - ny
            exit
         end if
      end if
   end do
      
   if (.NOT. found_magic) then
     call da_error("da_setup_runconstants.inc",80, &
       (/"No FFT factor found: invalid e_sn value"/))
   end if

   allocate(xbx % trig_functs_y(1:3*fft_size_j))  

   call da_find_fft_trig_funcs(fft_size_j, xbx % trig_functs_y(:))

   !-----Multiplicative coefficent for solution of spectral Poission eqn:

   !mslee.Bgrid
   ! const = -0.5 * grid%xb%ds * grid%xb%ds
   const = -2.0 * grid%xb%ds * grid%xb%ds

   nx = ide - ids + xbx%fft_pad_i
   ny = jde - jds + xbx%fft_pad_j
   ! YRG: A-grid:
   coeff_nx = 2.0 * pi / real(nx)
   coeff_ny = 2.0 * pi / real(ny)

   ! YRG: B-grid::
   ! coeff_nx = pi / real(nx)
   ! coeff_ny = pi / real(ny)

   xbx%fft_ix = nx + 1
   xbx%fft_jy = ny + 1

   allocate(xbx%fft_coeffs(1:xbx%fft_ix,1:xbx%fft_jy))

   do j = 2, ny
      cos_coeff_ny = COS(coeff_ny * real(j-1))
      do i = 2, nx
         cos_coeff_nx = COS(coeff_nx * real(i-1))
         !mslee.Bgrid
         ! xbx%fft_coeffs(i,j) = const / (1.0 - cos_coeff_nx * cos_coeff_ny)
         if (cos_coeff_nx.eq.1.and.cos_coeff_ny.eq.1) then
            xbx%fft_coeffs(i,j) = 0.0
         else
            xbx%fft_coeffs(i,j) = const / (2.0 - cos_coeff_nx - cos_coeff_ny)
         end if
      end do
   end do

   ! Set to zero all coefficients which are multiplied by sin(0.0) in FST:

   !mslee      xbx%fft_coeffs(1,1:xbx%fft_jy)  = 0.0
   !mslee      xbx%fft_coeffs(xbx%fft_ix,1:xbx%fft_jy) = 0.0
   !mslee      xbx%fft_coeffs(1:xbx%fft_ix,1)  = 0.0
   !mslee      xbx%fft_coeffs(1:xbx%fft_ix,xbx%fft_jy) = 0.0
   !mslee. we need to check the following

   xbx%fft_adjoint_factor = 4.0 / real(nx * ny)

   !---------------------------------------------------------------------------

   ! Calculate i increment for distributing pad region across processors.

   need_pad = .true.
   if (xbx%fft_pad_i <= 0) then
      need_pad = .false.
   else if (xbx%fft_pad_i > (ide-ids+1)) then
      write(unit=message(1), fmt='(a)') &
       'FFT xbx%fft_pad_i is too large!'
      write(unit=message(2), fmt='(2x,a,i4)') &
         '(ide-ids+1) = ', (ide-ids+1), &
         'xbx%fft_pad_i = ', xbx%fft_pad_i
      call da_error ("da_setup_runconstants.inc",146,message(1:2))
   else
      xbx%pad_inc = (ide-ids+1)/xbx%fft_pad_i
   end if

   ! Calculate number of pad vectors in x to be done on this processor.
   ! Need to save xbx%pad_num, xbx%pad_inc, and xbx%pad_loc

   xbx%pad_num = 0
   if (need_pad) then
      do n=1, xbx%fft_pad_i
         i = (n-1)*xbx%pad_inc + 1
         if ((i >= grid%xp%itsy) .and. (i <= grid%xp%itey)) then
            xbx%pad_num = xbx%pad_num + 1
         end if
      end do

      if (xbx%pad_num > 0) then
         allocate(xbx%pad_loc(1:xbx%pad_num))
         allocate(xbx%pad_pos(1:xbx%pad_num))
      end if

      xbx%pad_num = 0
      do n=1, xbx%fft_pad_i
         i = (n-1)*xbx%pad_inc + 1
         if ((i >= grid%xp%itsy) .and. (i <= grid%xp%itey)) then
            xbx%pad_num = xbx%pad_num + 1
            xbx%pad_loc(xbx%pad_num) = i
            xbx%pad_pos(xbx%pad_num) = grid%xp%ide + n
         end if
      end do
   end if
   
   !---------------------------------------------------------------------------

   if (global) then
      ! Set up Spectral transform constants:
      xbx%inc    = 1
      xbx%ni     = ide-ids+1
      xbx%nj     = jde-jds+1
      xbx%lenr   = xbx%inc * (xbx%ni - 1) + 1
      xbx%lensav = xbx%ni + int(log(real(xbx%ni))) + 4
      xbx%lenwrk = xbx%ni
      xbx%max_wavenumber = xbx%ni/2 -1
      xbx%alp_size = (xbx%nj+ 1)*(xbx%max_wavenumber+1)*(xbx%max_wavenumber+2)/4

      if (print_detail_spectral) then
         write (unit=stdout,fmt='(a)') "Spectral transform constants"
         write (unit=stdout, fmt='(a, i8)') &
            '   inc            =', xbx%inc   , &
            '   ni             =', xbx%ni    , &
            '   nj             =', xbx%nj    , &
            '   lenr           =', xbx%lenr  , &
            '   lensav         =', xbx%lensav, &
            '   lenwrk         =', xbx%lenwrk, &
            '   max_wavenumber =', xbx%max_wavenumber, &
            '   alp_size       =', xbx%alp_size
         write (unit=stdout,fmt='(a)') " "
      end if

      allocate(xbx%coslat(1:xbx%nj))
      allocate(xbx%sinlat(1:xbx%nj))
      allocate(xbx%coslon(1:xbx%ni))
      allocate(xbx%sinlon(1:xbx%ni))
      allocate(xbx%int_wgts(1:xbx%nj))      ! Interpolation weights
      allocate(xbx%alp(1:xbx%alp_size))
      allocate(xbx%wsave(1:xbx%lensav))

      call da_initialize_h(xbx%ni, xbx%nj, xbx%max_wavenumber, &
                           xbx%lensav, xbx%alp_size,           &
                           xbx%wsave, grid%xb%lon, xbx%sinlon,      &
                           xbx%coslon, grid%xb%lat, xbx%sinlat,     &
                           xbx%coslat, xbx%int_wgts, xbx%alp)


   end if

   if (trace_use) call da_trace_exit("da_setup_runconstants")
   
end subroutine da_setup_runconstants


subroutine da_cloud_model (TB, PB, QB, QCWB, QRNB, ZB, ZFB, DT, kts, kte)

   !-----------------------------------------------------------------
   ! Purpose: Calculate DT (=dz/w) using cumulus parameterization 
   !          of a one-dimensional cloud model.
   !-----------------------------------------------------------------

   ! Calculate DT

   implicit none

   integer, intent(in)                     :: kts, kte
   real, intent(in),  dimension(kts:kte)   :: TB, PB, QB, QCWB, QRNB, ZB
   real, intent(in),  dimension(kts:kte+1) :: ZFB
   real, intent(out), dimension(kts:kte)   :: DT

   integer                    :: k
   real                       :: P0, Z0, T0, Q0
   real                       :: PLCL, ZLCL, TLCL, QLCL
   integer                    :: KCB, KCT
   real                       :: PCT, ZCT
   real, dimension(kts:kte)   :: ZC, TC, QC, PP, QT
   real, dimension(kts:kte)   :: TCV, TBV, B
   real                       :: ALPHA, RC, MU, XX, YY
   real, dimension(kts:kte+1) :: W0, W

   if (trace_use) call da_trace_entry("da_cloud_model")

   ALPHA=0.5
   RC=100.0
   MU=0.183/RC

   do k = kts, kte+1
      W0(k)=0.0  
      W(k)=0.0  
   end do

   do k = kts, kte
      PP(k)=PB(k)/100.0
      DT(k)=0.0
   end do

   P0 = PP(kts)
   Z0 = ZB(kts)
   T0 = MAX(TB(kts),303.0)

   call da_qfrmrh (P0, T0, 95.0, Q0)

   call da_lcl (P0, Z0, T0, Q0, PLCL, ZLCL, TLCL, QLCL)

   call da_qfrmrh (PLCL, TLCL, 95.0, QLCL)

   call da_cumulus (ZLCL, TLCL, QLCL, PLCL, PP, TB,            &
                  ZC, TC, QC, KCB, KCT, PCT, ZCT, kts, kte)

   do k = KCB, KCT
      TCV(k) = TC(k) * (1.0 + 0.608 * QC(k))
      TBV(k) = TB(k) * (1.0 + 0.608 * QB(k))
   
      B(k) = (TCV(k)-TBV(k)) / TBV(k)

      QT(k) = QC(k) + QCWB(k) + QRNB(k)
   end do

   W0(KCB) = 0.0
   do k = KCB+1, KCT+1
      XX = 1.0+2.0*MU*(ZFB(k)-ZFB(k-1))
      YY = 2.0*gravity*(B(k-1)/(1.0+ALPHA) - QT(k-1)) * (ZFB(k)-ZFB(k-1))
      W0(k) =  (W0(k-1)+YY) / XX
   end do
     
   do k = KCB, KCT+1
      if (W0(k) >= 0.0) then
         W(k) = sqrt(W0(k))
      end if
   end do


   do k = KCT, KCB+1, -1
      if (W(k) >= 0.01) then
         DT(k) = (ZB(k)-ZB(k-1))/W(k)
      else
         DT(k) = 0.0
      end if
   end do

   if (trace_use) call da_trace_exit("da_cloud_model")

end subroutine da_cloud_model


subroutine da_lcl(p0, z0, t0, q0, plcl, zlcl, tlcl, qlcl)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)   :: p0, z0, t0, q0
   real, intent(out)  :: plcl, zlcl, tlcl, qlcl

   real   :: dp, qs, eps

   if (trace_use) call da_trace_entry("da_lcl")

   dp=5.0
   plcl=300.0

   do
      tlcl=t0*((plcl/p0)**0.286)

      call da_qfrmrh(plcl, tlcl, 100.0, qs)

      eps=qs-q0

      if (eps >= 0.0) then
         zlcl=(1004.0/gravity)*(t0-tlcl)+z0
         qlcl=qs
         if (trace_use) call da_trace_exit("da_lcl")
         return
      else
         plcl=plcl+dp

         if (plcl >= p0) then
            zlcl=z0
            qlcl=q0
            plcl=p0
            if (trace_use) call da_trace_exit("da_lcl")
            return
        end if
      end if
   end do

   if (trace_use) call da_trace_exit("da_lcl")

end subroutine da_lcl


subroutine da_cumulus (zcb, tcb, qcb, pcb, pk, te, z, t, q, lcb, lct, pct, zct, kts, kte)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: kts, kte
   real,    intent(inout) :: zcb, tcb, qcb, pcb
   real,    intent(in)    :: pk(kts:kte)
   real,    intent(in)    :: te(kts:kte)
   real,    intent(out)   :: z(kts:kte)
   real,    intent(out)   :: t(kts:kte)
   real,    intent(out)   :: q(kts:kte)
   integer, intent(out)   :: lcb, lct
   real,    intent(out)   :: pct, zct

   integer   :: k, ia, l, ncb
   real      :: cp, r, hl, em, et, p
   real      :: tll, qll, pll, zll, tbar, pbar, qbar
   real      :: dp, dz, ddt, dt 

   if (trace_use) call da_trace_entry("da_cumulus")    

   cp=1004.0
   r=2000.0/7.0
   hl=2.49e06
   dt=0.1
   ia=1000

   do k = kts, kte
      z(k) = 0.0
      t(k) = 0.0
      q(k) = 0.0
   end do

   em=gravity*zcb+cp*tcb+hl*qcb

   ncb=kts

   if (pk(kte) > pcb) then
      ncb=kte
   end if

   do l=kte-1,kts,-1
      if (pk(l) > pcb) then
         ncb=l+1
         exit
      end if
   end do

   do l=ncb,kte
      p=pk(l)
      do k=1,ia
         if (l == ncb) then
            tll=tcb
            qll=qcb
            pll=pcb
            zll=zcb
         else
            tll=t(l-1)
            qll=q(l-1)
            pll=pk(l-1)
            zll=z(l-1)
         end if

         t(l)=tll-(k*dt)

         call da_qfrmrh(p, t(l), 100.0, q(l))

         tbar=0.5*(t(l)+tll)
         qbar=0.5*(q(l)+qll)
         pbar=0.5*(p+pll)
         dp=pll-p
         dz=(r*tbar*(1.0+0.61*qbar)*dp)/(gravity*pbar)
         z(l)=zll+dz
         et=gravity*z(l)+cp*t(l)+hl*q(l)
         if ((et-em) <= 0.0) exit
      end do
   end do

   lct=ncb

   do k=kte,ncb+1,-1
      ddt=t(k)-te(k)

      if (ddt >= 0.0) then
         lct=k
         exit
      end if
   end do

   lcb=lct

   do k=ncb,kte
      ddt=t(k)-te(k)
      if (ddt >= 0.0) then
         lcb=k
         exit
      end if
   end do

   pct=pk(lct)
   zct=z(lct)
   pcb=pk(lcb)
   zcb=z(lcb)

   if (trace_use) call da_trace_exit("da_cumulus")    

end subroutine da_cumulus


subroutine da_qfrmrh ( p, t, rh, q )

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)  :: p, t, rh
   real, intent(out) :: q

   real   :: a, b, e, qs

   if (trace_use) call da_trace_entry("da_qfrmrh")

   a=17.26
   b=35.86
   if (t <= 263.0) a=21.87
   if (t <= 263.0) b= 7.66
   e  = 6.11*exp(a*(t-t_triple)/(t-b))
   qs = 0.622*e/(p-0.378*e)
   q  = qs*rh/100.0

   if (trace_use) call da_trace_exit("da_qfrmrh")

end subroutine da_qfrmrh


subroutine da_write_increments (grid, q_cgrid, mu_cgrid, ph_cgrid)

   !----------------------------------------------------------------------
   ! Purpose: Write analysis increments
   !----------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)                         :: grid
   real,intent(in) :: q_cgrid(ims:ime,jms:jme,kms:kme)
   real,intent(in) :: ph_cgrid(ims:ime,jms:jme,kms:kme)
   real,intent(in) :: mu_cgrid(ims:ime,jms:jme)

   ! Arrays for write out increments:
   integer                                          :: ix, jy, kz
   real, dimension(1:grid%xb%mix,1:grid%xb%mjy)               ::     gbuf_2d
   real, dimension(1:grid%xb%mix+1,1:grid%xb%mjy+1)           ::     gbuf_2dd
   real, dimension(1:grid%xb%mix,1:grid%xb%mjy,1:grid%xb%mkz)      ::     gbuf

   real, dimension(1:grid%xb%mix,1:grid%xb%mjy,1:grid%xb%mkz+1)    ::    wgbuf
   real, dimension(:,:,:), allocatable :: u_global, v_global, w_global, &
                                          p_global, t_global, q_global, &
                                         ph_global
   real, dimension(:,:)  , allocatable :: mu_global, psfc_global, &
                       psac_global, tgrn_global, terr_global, snow_global,&
                        lat_global,  lon_global, lanu_global,             &
                 map_factor_global, cori_global, landmask_global

   integer :: anl_inc_unit

   if (trace_use) call da_trace_entry("da_write_increments")


   ! Dimension of the domain:
   ix = grid%xb%mix
   jy = grid%xb%mjy
   kz = grid%xb%mkz

 
   ! 3-d and 2-d increments:

   allocate (   p_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (   t_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (   q_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (   u_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (   v_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (   w_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (  ph_global (1:ix+1,1:jy+1,1:kz+1))
   allocate (psfc_global (1:ix+1,1:jy+1))
   allocate (  mu_global (1:ix+1,1:jy+1))
   call da_patch_to_global(grid, grid%xa % p, gbuf) 
   if (rootproc) then 
      p_global(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 
   call da_patch_to_global(grid, grid%xa % t, gbuf) 
   if (rootproc) then 
      t_global(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 
   call da_patch_to_global(grid, q_cgrid, gbuf) 
   if (rootproc) then 
      q_global(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 
   call da_patch_to_global(grid, grid%xa % u, gbuf) 
   if (rootproc) then 
      u_global(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 
   call da_patch_to_global(grid, grid%xa % v, gbuf) 
   if (rootproc) then 
      v_global(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if

   ! One more level for w and ph:
   grid%xp%kde=grid%xp%kde+1
   kde=kde+1
   call da_patch_to_global(grid, grid%xa % w, wgbuf) 
   if (rootproc) then 
      w_global(1:ix,1:jy,1:kz+1) = wgbuf(1:ix,1:jy,1:kz+1) 
   end if 
   call da_patch_to_global(grid, ph_cgrid, wgbuf) 
   if (rootproc) then 
      ph_global(1:ix,1:jy,1:kz+1) = wgbuf(1:ix,1:jy,1:kz+1) 
   end if 
   kde=kde-1
   grid%xp%kde=grid%xp%kde-1
 
   call da_patch_to_global(grid, grid%xa % psfc, gbuf_2d) 
   if (rootproc) then 
      psfc_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if 
   call da_patch_to_global(grid, mu_cgrid, gbuf_2d) 
   if (rootproc) then 
      mu_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if 

   ! 2d constant fields:

   allocate (      psac_global (1:ix+1,1:jy+1))
   allocate (      tgrn_global (1:ix+1,1:jy+1))
   allocate (      terr_global (1:ix+1,1:jy+1))
   allocate (      snow_global (1:ix+1,1:jy+1))
   allocate (       lat_global (1:ix+1,1:jy+1))
   allocate (       lon_global (1:ix+1,1:jy+1))
   allocate (      lanu_global (1:ix+1,1:jy+1))
   allocate (map_factor_global (1:ix+1,1:jy+1))
   allocate (      cori_global (1:ix+1,1:jy+1))
   allocate (  landmask_global (1:ix+1,1:jy+1))

   call da_patch_to_global(grid, grid%xb%psac, gbuf_2d) 
   if (rootproc) then 
      psac_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%tgrn, gbuf_2d) 
   if (rootproc) then 
      tgrn_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%terr, gbuf_2d) 
   if (rootproc) then 
      terr_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%snow, gbuf_2d) 
   if (rootproc) then 
      snow_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%lat , gbuf_2d) 
   if (rootproc) then 
      lat_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%lon , gbuf_2d) 
   if (rootproc) then 
      lon_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%lanu, gbuf_2d) 
   if (rootproc) then 
      lanu_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if
   call da_patch_to_global(grid, grid%xb%map_factor, gbuf_2d) 
   if (rootproc) then 
      map_factor_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if

   ! temporary increase to dimensions for cori
   ide=ide+1
   jde=jde+1
   grid%xp%ide=grid%xp%ide+1
   grid%xp%jde=grid%xp%jde+1
   call da_patch_to_global(grid, grid%xb%cori, gbuf_2dd) 
   if (rootproc) then
      cori_global(1:ix+1,1:jy+1) = gbuf_2dd(1:ix+1,1:jy+1) 
   end if
   ide=ide-1
   jde=jde-1
   grid%xp%ide=grid%xp%ide-1
   grid%xp%jde=grid%xp%jde-1

   call da_patch_to_global(grid, grid%xb%landmask, gbuf_2d)
   if (rootproc) then 
      landmask_global(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if


   if (rootproc) then
      call da_get_unit(anl_inc_unit)
      open(unit=anl_inc_unit, file='analysis_increments', form='unformatted')

      write (unit=anl_inc_unit) ANALYSIS_DATE

      write (unit=anl_inc_unit) 1, ix, 1, jy, 1, kz 

      ! Map projection information:
      write (unit=anl_inc_unit) map_projection, coarse_ix, coarse_jy
      write (unit=anl_inc_unit) &
         coarse_ds, start_x, start_y, &
         phic, xlonc, cone_factor, truelat1_3dv, truelat2_3dv, pole, dsm,   &
         psi1, c2, ptop, base_pres, t0, base_lapse, base_temp

      ! 1d constant fields:

      write (unit=anl_inc_unit) grid%xb%sigmah, grid%xb%sigmaf


      ! 3d- and 2d-increments:
      write (unit=anl_inc_unit) u_global, v_global, w_global, p_global, &
         t_global, q_global, ph_global, mu_global, psfc_global

      ! 2d-constant fields:
      write (unit=anl_inc_unit) psac_global, tgrn_global, terr_global, &
         snow_global, lat_global, lon_global, lanu_global, map_factor_global, &
         cori_global, landmask_global
      close(anl_inc_unit)
      call da_free_unit(anl_inc_unit)

   end if

   if (trace_use) call da_trace_exit("da_write_increments")

end subroutine da_write_increments


subroutine da_write_increments_for_wrf_nmm_regional (grid)                                 

   !----------------------------------------------------------------------
   ! Purpose: Write analysis increments
   !----------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)                                    :: grid
   integer                                                         :: ix, jy, kz
   real, dimension(1:grid%xb%mix,1:grid%xb%mjy)                    :: gbuf_2d
   real, dimension(1:grid%xb%mix,1:grid%xb%mjy,1:grid%xb%mkz)      :: gbuf

   integer :: anl_inc_unit

   if (trace_use) call da_trace_entry("da_write_increments_for_wrf_nmm_regional")


   ! Dimension of the domain:
   ix = grid%xb%mix
   jy = grid%xb%mjy
   kz = grid%xb%mkz

   if (rootproc) then
      call da_get_unit(anl_inc_unit)
      open(unit=anl_inc_unit, file='analysis_increments_for_wrf-nmm', form='unformatted')

      write (unit=anl_inc_unit) analysis_date 

      write (unit=anl_inc_unit) ix,jy, kz 
   end if


   call da_patch_to_global(grid, grid%xa % u, gbuf) 
   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa % v, gbuf) 
   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa % t, gbuf) 
   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa % q, gbuf) 
   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf(1:ix,1:jy,1:kz) 
   end if 


   call da_patch_to_global(grid, grid%xa % psfc, gbuf_2d) 

   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf_2d(1:ix,1:jy) 
   end if 

   call da_patch_to_global(grid, grid%xb % lat, gbuf_2d) 

   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf_2d(1:ix,1:jy) 
   end if 

   call da_patch_to_global(grid, grid%xb % lon, gbuf_2d) 

   if (rootproc) then 
      write (unit=anl_inc_unit) gbuf_2d(1:ix,1:jy) 
   end if 

      close(anl_inc_unit)
      call da_free_unit(anl_inc_unit)

   if (trace_use) call da_trace_exit("da_write_increments_for_wrf_nmm_regional")

end subroutine da_write_increments_for_wrf_nmm_regional
subroutine da_write_kma_increments(grid)

   !---------------------------------------------------------------------------
   ! Purpose: Gathers KMA analysis increments and writes 
   !           on "anl_inc_unit" unit 
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(in) :: grid

   ! Arrays for write out increments:
   integer                                     :: ix, jy, kz

   real, dimension(1:grid%xb%mix,1:grid%xb%mjy)          :: gbuf_2d
   real, dimension(1:grid%xb%mix,1:grid%xb%mjy,1:grid%xb%mkz) :: gbuf
   real, dimension(:,:)  , allocatable :: psfc_g
   real, dimension(:,:,:), allocatable :: u_g, v_g, t_g, q_g, p_g

   integer                                     :: i, j, k,anl_inc_unit

   if (trace_use) call da_trace_entry("da_write_kma_increments")

   ! Dimension of the domain:
   ix = grid%xb%mix
   jy = grid%xb%mjy
   kz = grid%xb%mkz

 
   ! 3-d and 2-d increments:

   allocate (psfc_g (1:ix,1:jy))
   allocate (   u_g (1:ix,1:jy,1:kz))
   allocate (   v_g (1:ix,1:jy,1:kz))
   allocate (   t_g (1:ix,1:jy,1:kz))
   allocate (   q_g (1:ix,1:jy,1:kz))
   allocate (   p_g (1:ix,1:jy,1:kz))

   call da_patch_to_global(grid, grid%xa%psfc, gbuf_2d) 
   if (rootproc) then 
      psfc_g(1:ix,1:jy) = gbuf_2d(1:ix,1:jy) 
   end if 

   call da_patch_to_global(grid, grid%xa%u, gbuf) 
   if (rootproc) then 
      u_g(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa%v, gbuf) 
   if (rootproc) then 
      v_g(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa%t, gbuf) 
   if (rootproc) then 
      t_g(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa%q, gbuf) 
   if (rootproc) then 
      q_g(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 

   call da_patch_to_global(grid, grid%xa%p, gbuf) 
   if (rootproc) then 
      p_g(1:ix,1:jy,1:kz) = gbuf(1:ix,1:jy,1:kz) 
   end if 

   if (rootproc) then
      ! 3d- and 2d-increments:

      call da_get_unit(anl_inc_unit)
      open(unit=anl_inc_unit,file="analysis_increments_kma",status="replace", &
         form="unformatted")
      write(anl_inc_unit) ((psfc_g(i,j),i=ids,ide),j=jds,jde)
      write(anl_inc_unit) (((u_g(i,j,k),i=ids,ide),j=ids,jde),k=kds,kde)
      write(anl_inc_unit) (((v_g(i,j,k),i=ids,ide),j=ids,jde),k=kds,kde)
      write(anl_inc_unit) (((t_g(i,j,k),i=ids,ide),j=ids,jde),k=kds,kde)
      write(anl_inc_unit) (((q_g(i,j,k),i=ids,ide),j=ids,jde),k=kds,kde)
      write(anl_inc_unit) (((p_g(i,j,k),i=ids,ide),j=ids,jde),k=kds,kde)
      close(anl_inc_unit)
      call da_free_unit(anl_inc_unit)
   end if

   if (trace_use) call da_trace_exit("da_write_kma_increments")

end subroutine da_write_kma_increments 


subroutine da_get_bins_info(nj, nk, bin2d, evec_g, eval_g,&
                  evec_loc, eval_loc, max_vert_var, var_scaling, be_sub)

   !---------------------------------------------------------------------------
   !  Purpose: Extracts Eigen vectors/values info from bins
   !           and builds up background error structure 
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: nj, nk       ! W-E, S-N and Vert. Dims 
   integer, intent(in)    :: bin2d(:,:)       ! Bin assigned to each 2D
   real*8,  intent(in)    :: evec_g(:,:)      ! Global Eig. vectors 
   real*8,  intent(in)    :: eval_g(:)        ! Global Eig. values  
   real*8,  intent(in)    :: evec_loc(:,:,:)  ! Local Eig. vectors   
   real*8,  intent(in)    :: eval_loc(:,:)    ! Local Eig. values  
   real,    intent(in)    :: max_vert_var     ! Vertical variance
   real,    intent(in)    :: var_scaling      ! Variance re-scaling factor
   type(be_subtype), intent(inout) :: be_sub ! Background error structure
  
   real*8, allocatable    :: e_vec_loc(:,:,:)
   real*8, allocatable    :: e_val_loc(:,:)
  
   integer               :: i, j, k, b

   !---------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_get_bins_info")

   allocate(e_vec_loc(1:nj,1:nk,1:nk))
   allocate(e_val_loc(1:nj,1:nk))
   if (vert_evalue == vert_evalue_global) then       ! use global eigen vectors
      do k = 1, nk
         e_val_loc(1:nj,k)     = eval_g(k)
         do i =1,nk
            e_vec_loc(1:nj,k,i) = evec_g(k,i)  
         end do  
      end do  
  
   else if (vert_evalue == vert_evalue_local) then  ! use local  eigen vectors
      do j = 1, nj
         b = bin2d(1,j)         
         do k=1,nk
           e_val_loc(j,k)   = eval_loc(k,b)
           do i = 1,nk
              e_vec_loc(j,k,i) = evec_loc(k,i,b)
            end do
         end do
      end do
   else
      write(unit=message(1),fmt='(A,I5)') &
         "Invalid value of vert_evalue=",vert_evalue
      call da_error("da_get_bins_info.inc",53,message(1:1))
   end if

   call da_get_vertical_truncation(max_vert_var, eval_g, be_sub)
   call da_allocate_background_errors(nj, nk, eval_g, evec_g, &
      e_val_loc, e_vec_loc, be_sub)
   if (be_sub%mz > 0) be_sub%val(:,:) = var_scaling * be_sub%val(:,:)
   deallocate( e_vec_loc)
   deallocate( e_val_loc)

   if (trace_use) call da_trace_exit("da_get_bins_info")
  
end subroutine da_get_bins_info


subroutine da_truncate_spectra(max_wave, nw, power_trunc, power, max_wave_trunc)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: max_wave         ! Smallest wave for domain.
   integer, intent(in)  :: nw               ! Dimension of input power spectrum.
   real,    intent(in)  :: power_trunc      ! Power truncation (fraction).
   real*8,  intent(in)  :: power(0:nw)      ! Power spectrum.
   integer, intent(out) :: max_wave_trunc   ! Smallest wave after truncation.

   integer :: l                ! Loop counter.
   real    :: truncated_power  ! Truncated power.
   real    :: cumul_power      ! Cumulative power.

   if (trace_use) call da_trace_entry("da_truncate_spectra")

   truncated_power = power_trunc * sum(power(0:nw))

   cumul_power = 0.0
   max_wave_trunc = max_wave
   do l = 0, nw - 1 
      cumul_power = cumul_power + power(l)
      if (cumul_power > truncated_power) then
         max_wave_trunc = l - 1
         exit
      end if
   end do

   if (max_wave_trunc > max_wave) then
      write(unit=message(1),fmt='(a)') &
         'da_truncate_spectra: Power requested needs higher resolution.'     
      write(unit=message(2),fmt='(a,i8)') &
         'Maximum grid wavenumber =  ', max_wave
      write(unit=message(3),fmt='(a,i8)') &
         'Truncating to wavenumber = ', max_wave_trunc
      call da_warning("da_truncate_spectra.inc",40,message(1:3))
      max_wave_trunc = max_wave
   end if

   if (trace_use) call da_trace_exit("da_truncate_spectra")

end subroutine da_truncate_spectra


SUBROUTINE da_chg_be_Vres(kz, nk, eta_h, eta_be,&
                        reg_chi_be,reg_ps_be,reg_t_be, &
                        reg_chi   ,reg_ps   ,reg_t   , &
                        covm1_be, covm2_be, covm3_be, covm4_be, &
                        covm1   , covm2   , covm3   , covm4   , &
                        rfls1_be, rfls2_be, rfls3_be, rfls4_be, &
                        rfls1, rfls2, rfls3, rfls4)

!------------------------------------------------------------------------------
!  PURPOSE: Change vertical resolution of background stats for cv_options=3
!------------------------------------------------------------------------------
   integer, intent(in) :: kz, nk  ! model V-resolution and BE V-resolution
   real   , intent(in) :: eta_h(kz), eta_be(nk) ! Eta level definition
!
!  Regression coefficient's arrays:
   real, dimension(1:nk),      intent(in) :: reg_chi_be, reg_ps_be 
   real, dimension(1:nk,1:nk), intent(in) :: reg_t_be
   real, dimension(1:kz),      intent(out):: reg_chi   , reg_ps 
   real, dimension(1:kz,1:kz), intent(out):: reg_t

! Vertical Convarince matrix arrays:
   real, dimension(1:nk,1:nk), intent(in) :: covm1_be, covm2_be, &
                                             covm3_be, covm4_be
   real, dimension(1:kz,1:kz), intent(out):: covm1   , covm2   , &
                                             covm3   , covm4 

! Recursive filter scale-length arrays:
   real, dimension(nk),        intent(in) :: rfls1_be, rfls2_be, &
                                             rfls3_be, rfls4_be
   real, dimension(kz),        intent(out):: rfls1   , rfls2   , &
                                             rfls3   , rfls4  


   integer             :: i,j,k,m,l,l1,m1,n
   integer             :: lsig(kz)
   real                :: rsig(kz), rsigo(nk)
   real                :: coef1(kz),coef2(kz)
   logical             :: NO_INTERP
! ---------------------------------------------------------------------

   NO_INTERP = .FALSE.

! Check if the # of levels and the values of eta are same: 
!
   NO_INTERP = (kz == nk)
   if (NO_INTERP) then
      do k = 1, nk
        if (abs(eta_h(k)-eta_be(k)) > 1.0e-6) then
          NO_INTERP = .FALSE.
          exit
        endif
       enddo
   endif

   if(NO_INTERP )then
! Regression coefficients:
     reg_chi = reg_chi_be
     reg_ps  = reg_ps_be
     reg_t   = reg_t_be

! Vertical covarince matrix:
     covm1 = covm1_be
     covm2 = covm2_be
     covm3 = covm3_be
     covm4 = covm4_be

! Recursive filter scale-length:
     rfls1 = rfls1_be
     rfls2 = rfls2_be
     rfls3 = rfls3_be
     rfls4 = rfls4_be

     return
   endif

!   if (.not.NO_INTERP) then
!     write(6,'(/10X,a/2X,"Model eta levels:")') &
!             "Vertical resolution conversion needed for CV5 BE:::"
!     write(6,'(2X,I3,2X,f10.5)') ( k,eta_h(k),k=1,kz)
!   endif

! Convert Eta to log(eta):
    do k=1,kz
      rsig(k)=log(eta_h(k))
    enddo
    do k=1,nk
      rsigo(k)=log(eta_be(k))
    enddo

! Find the coef1 and coef2 for the vertical interpolation:
!  
  do k=1,kz

! Model levels below the lowest BE level:
  if(rsig(k).ge.rsigo(1))then
     m=1
     m1=2
     lsig(k)=1
     coef1(k)=1.
     coef2(k)=0.

! Model levels above the highest BE level:
  else if(rsig(k).lt.rsigo(nk))then
     m=nk-1
     m1=nk
     lsig(k)=nk-1
     coef1(k)=0.
     coef2(k)=1

! Model levels located within the BE levels:
  else
     do m=1,nk
       m1=m+1
       if((rsig(k).le.rsigo(m))   .and.  &
          (rsig(k).gt.rsigo(m1))     )then
         lsig(k)=m
        go to 2345
       end if
     end do

2345    continue
    coef1(k)=(rsigo(m1)-rsig(k))/(rsigo(m1)-rsigo(m))
    coef2(k)=1.-coef1(k)
     if(lsig(k)==nk)then
     lsig(k)=nk-1
     coef2(k)=1.
     coef1(k)=0.
     endif
 endif

   end do

  do k=1,kz
    m=lsig(k)
    m1=m+1
! Interpolation for Regression coefficients:
      reg_chi(k)=reg_chi_be(m)*coef1(k)+reg_chi_be(m1)*coef2(k)
      reg_ps(k) =reg_ps_be (m)*coef1(k)+reg_ps_be(m1) *coef2(k)

! Recursive filter scale-lengths:
      rfls1(k) =rfls1_be (m)*coef1(k)+rfls1_be(m1) *coef2(k)
      rfls2(k) =rfls2_be (m)*coef1(k)+rfls2_be(m1) *coef2(k)
      rfls3(k) =rfls3_be (m)*coef1(k)+rfls3_be(m1) *coef2(k)
      rfls4(k) =rfls4_be (m)*coef1(k)+rfls4_be(m1) *coef2(k)

    do j=1,kz
      l=lsig(j)
      l1=l+1
! Interpolation for Regression coefficients:
        reg_t(j,k)=(reg_t_be(l,m)*coef1(j)+reg_t_be(l1,m)*coef2(j))*coef1(k) &
                  +(reg_t_be(l,m1)*coef1(j)+reg_t_be(l1,m1)*coef2(j))*coef2(k)
! Vertical covariance matrix:
        covm1(j,k)=(covm1_be(l,m)*coef1(j)+covm1_be(l1,m)*coef2(j))*coef1(k) &
                  +(covm1_be(l,m1)*coef1(j)+covm1_be(l1,m1)*coef2(j))*coef2(k)
        covm2(j,k)=(covm2_be(l,m)*coef1(j)+covm2_be(l1,m)*coef2(j))*coef1(k) &
                  +(covm2_be(l,m1)*coef1(j)+covm2_be(l1,m1)*coef2(j))*coef2(k)
        covm3(j,k)=(covm3_be(l,m)*coef1(j)+covm3_be(l1,m)*coef2(j))*coef1(k) &
                  +(covm3_be(l,m1)*coef1(j)+covm3_be(l1,m1)*coef2(j))*coef2(k)
        covm4(j,k)=(covm4_be(l,m)*coef1(j)+covm4_be(l1,m)*coef2(j))*coef1(k) &
                  +(covm4_be(l,m1)*coef1(j)+covm4_be(l1,m1)*coef2(j))*coef2(k)
    enddo
  enddo

!--------------------------------------------------------------------
end subroutine da_chg_be_Vres
Subroutine da_gen_eigen(kz, covm, eve, val)

  implicit none

  integer, intent(in) :: kz
  real   , dimension(1:kz,1:kz), intent(inout) :: covm
  real   , dimension(1:kz,1:kz), intent(out)   :: eve
  real   , dimension(1:kz)     , intent(out)   :: val

  integer  :: k1, k2

!
      call da_1d_eigendecomposition( covm, eve, val )

      do k1 = 1, kz
        if ( val(k1) < 0.0 ) val(k1) = 0.0
      end do

      do k1 = 1, kz
        do k2 = k1, kz
          covm(k1,k2) = SUM( eve(k1,:) * val(:) * eve(k2,:) )
        end do
      end do
   
      do k1 = 2, kz
        do k2 = 1, k1-1
          covm(k1,k2) = covm(k2,k1)
        end do
      end do

      call da_1d_eigendecomposition( covm, eve, val )

      do k1 = 1, kz
        if ( val(k1) < 0.0 ) val(k1) = 0.0
      end do
!
end Subroutine da_gen_eigen
Subroutine da_Eigen_to_covmatrix(kzs, be_evec_glo, be_eval_glo, b_native)
!---------------------------------------------------------------------------
!   Recalculate native global vertical background error cov matrix:
!---------------------------------------------------------------------------
   implicit none
   
   integer, intent(in) :: kzs
   real   , dimension(1:kzs,1:kzs), intent(in) :: be_evec_glo
   real   , dimension(1:kzs)                   :: be_eval_glo
   real   , dimension(1:kzs,1:kzs), intent(out):: b_native

   integer    :: k1, k2

     do k1 = 1, kzs
        do k2 = k1, kzs
          b_native(k1,k2) = SUM( be_evec_glo(k1,:) * be_eval_glo(:) * &
                                be_evec_glo(k2,:) )
        end do
     end do
   
     do k1 = 2, kzs
        do k2 = 1, k1-1
          b_native(k1,k2) = b_native(k2,k1)
        end do
     end do

end  Subroutine da_Eigen_to_covmatrix

end module da_setup_structures
