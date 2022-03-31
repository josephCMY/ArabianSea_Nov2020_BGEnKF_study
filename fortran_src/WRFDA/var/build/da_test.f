












module da_test

   
   
   

   use module_configure, only : grid_config_rec_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer

   use module_dm, only : local_communicator, &
      ntasks_x, ntasks_y, data_order_xyz, mytask, &
      ntasks, data_order_xy
   use module_comm_dm, only : halo_psichi_uv_adj_sub, halo_xa_sub, &
      halo_sfc_xa_sub, halo_ssmi_xa_sub, halo_radar_xa_w_sub


   use da_control, only : num_procs, var4d_bin, var4d_lbc
   use module_domain, only : vp_type, xb_type, x_type, ep_type, &
      domain, domain_clock_get, domain_clock_set, domain_clockprint, domain_clockadvance
   use module_state_description, only : dyn_em,dyn_em_tl,dyn_em_ad,p_a_qv

   use da_control, only : trace_use,ierr, trace_use_dull, comm,global,stdout,rootproc, &
      sfc_assi_options,typical_qrn_rms,typical_qci_rms,typical_qsn_rms,typical_qgr_rms,jcdfi_use, jcdfi_diag, &
      typical_u_rms,typical_v_rms,typical_w_rms,typical_t_rms, typical_p_rms, typical_rain_rms, &
      typical_q_rms,typical_qcw_rms,print_detail_testing,typical_rh_rms, &
      fg_format, fg_format_wrf_arw_global, fg_format_wrf_arw_regional,fg_format_wrf_nmm_regional, &
      typical_rf_rms,typical_rv_rms, typical_thickness_rms, typical_tb19v_rms,typical_tb37h_rms, &
      typical_tb85h_rms,typical_tb37v_rms,typical_tb85v_rms,typical_tb22v_rms, &
      typical_tb19h_rms,typical_speed_rms,typical_tpw_rms,typical_ref_rms, &
      cv_options_hum,inv_typ_vp5_sumsq,inv_typ_vp1_sumsq, trajectory_io, &
      inv_typ_vp3_sumsq,inv_typ_vp2_sumsq,inv_typ_vpalpha_sumsq, &
      inv_typ_vp4_sumsq,typical_rho_rms,balance_geo,balance_cyc,balance_type, &
      balance_geocyc, var4d, num_fgat_time,cv_options_hum_specific_humidity, &
      cv_options_hum_relative_humidity, ids, ide, jds, jde, kds, kde, &
      sound, sonde_sfc, mtgirs, synop, profiler, gpsref, gpspw, polaramv, geoamv, ships, metar, &
      satem, radar, ssmi_rv, ssmi_tb, ssmt1, ssmt2, airsr, pilot, airep, tamdar, tamdar_sfc, rain, &
      bogus, buoy, qscat, pseudo, radiance, use_radarobs, use_ssmiretrievalobs,use_rainobs, &
      use_gpsrefobs, use_ssmt1obs, use_ssmitbobs, use_ssmt2obs, use_gpspwobs, &
      use_gpsztdobs, use_radar_rf, use_radar_rhv, use_rad, crtm_cloud, cloud_cv_options, &
      ids,ide,jds,jde,kds,kde, ims,ime,jms,jme,kms,kme, fgat_rain_flags, &
      its,ite,jts,jte,kts,kte, ips,ipe,jps,jpe,kps,kpe, cv_options, cv_size, &
      cloud_cv_options, cp, gas_constant, test_dm_exact, cv_size_domain, &
      its_int, ite_int, jts_int, jte_int, kts_int, kte_int, &
      ims_int, ime_int, jms_int, jme_int, kms_int, kme_int

   use da_define_structures, only : da_zero_x,da_zero_vp_type,da_allocate_y, &
      da_deallocate_y,be_type, xbx_type, iv_type, y_type, j_type, da_initialize_cv
   use da_dynamics, only : da_uv_to_divergence,da_uv_to_vorticity, &
      da_psichi_to_uv, da_psichi_to_uv_adj
   use da_ffts, only : da_solve_poissoneqn_fct
   use da_minimisation, only : da_transform_vtoy_adj,da_transform_vtoy, da_swap_xtraj, &
       da_read_basicstates, da_calculate_j
   use da_obs, only : da_transform_xtoy,da_transform_xtoy_adj
   use da_par_util, only : da_patch_to_global, da_system, da_cv_to_global
   use da_par_util1, only : true_mpi_real
   use da_physics, only : da_transform_xtopsfc,da_transform_xtopsfc_adj, &
      da_pt_to_rho_lin,da_transform_xtotpw,da_transform_xtogpsref_lin, &
      da_transform_xtowtq, da_transform_xtowtq_adj,da_pt_to_rho_adj, &
      da_transform_xtotpw_adj, da_transform_xtoztd_lin, da_transform_xtoztd_adj, &
      da_moist_phys_lin, da_moist_phys_adj, da_uvprho_to_w_lin, da_uvprho_to_w_adj
   use da_reporting, only : da_error, message, da_message
   use da_spectral, only : da_test_spectral
   use da_ssmi, only : da_transform_xtoseasfcwind_lin, &
      da_transform_xtoseasfcwind_adj
   use da_statistics, only : da_correlation_coeff1d,da_correlation_coeff2d
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tracing, only : da_trace_entry,da_trace_exit
   use da_transfer_model, only : da_transfer_wrftltoxa,da_transfer_xatowrftl, &
      da_transfer_xatowrftl_adj,da_transfer_wrftltoxa_adj,da_transfer_wrftoxb
   
   
   use da_wrf_interfaces, only : wrf_debug, wrf_shutdown
   use da_wrfvar_io, only : da_med_initialdata_output,da_med_initialdata_input
   use da_vtox_transforms, only : da_transform_xtotb_lin, &
      da_transform_xtotb_adj, da_vertical_transform, da_transform_vptox, &
      da_transform_xtogpsref_adj,da_transform_vptox_adj,da_transform_vtox, &
      da_transform_vtox_adj,da_transform_vtovv,da_transform_vtovv_global, &
      da_transform_vtovv_global_adj, da_transform_vtovv_adj, da_transform_xtoxa, &
      da_transform_xtoxa_adj, da_apply_be, da_apply_be_adj, da_transform_bal, &
      da_transform_bal_adj

   implicit none

   private :: da_dot_cv, da_dot

   include 'mpif.h'

contains

subroutine da_check_balance(phi, phi_u)

   !---------------------------------------------------------------------------
   ! Purpose: Compare balanced mass (phi_b - function of wind) and actual phi.
   !
   ! Method:  Calculate correlation between balanced and actual phi.
   !---------------------------------------------------------------------------

   implicit none
      
   real, intent(in)             :: phi(:,:,:)      ! Total phi.
   real, intent(in)             :: phi_u(:,:,:)    ! Unbalanced phi.

   integer                      :: iy              ! Size of 1st dimension.
   integer                      :: jx              ! Size of 2nd dimension.
   integer                      :: kz              ! Size of 3rd dimension.
   integer                      :: i, k            ! Loop counters
   real                         :: corr_coef       ! Correlation coefficient.
   real                         :: accurac         ! Accuracy.
   real, allocatable            :: phi_b1(:)       ! 1D balanced phi.
   real, allocatable            :: phi_b2(:,:)     ! 2D balanced phi.
   real, allocatable            :: corr_coeff(:,:) ! Correlation coefficient.
   real, allocatable            :: accuracy(:,:)   ! Accuracy.

   if (trace_use) call da_trace_entry("da_check_balance")
          
   if (balance_type == balance_geo) then
      write(unit=stdout, fmt='(a)') ' da_check_balance: Balance is geostrophic.'
   else if (balance_type == balance_cyc) then
      write(unit=stdout, fmt='(a)') &
         ' da_check_balance: Balance is cyclostrophic.'
   else if (balance_type == balance_geocyc) then
      write(unit=stdout, fmt='(a)') &
         ' da_check_balance: Balance is geo/cyclostrophic.'
   end if
      
   write(unit=stdout, fmt='(a)') ' da_check_balance: Correlation/accuracy: '
      
   !-------------------------------------------------------------------------
   ! [1.0]: Initialise:
   !-------------------------------------------------------------------------  

   iy = size(phi_u, DIM=1)
   jx = size(phi_u, DIM=2)
   kz = size(phi_u, DIM=3)
      
   allocate(phi_b1(1:jx))
   allocate(phi_b2(1:iy,1:jx))

   allocate(corr_coeff(1:kz,1:iy))
   corr_coeff(1:kz,1:iy) = 0.0

   allocate(accuracy(1:kz,1:iy))
   accuracy(1:kz,1:iy) = 0.0
      
   !-------------------------------------------------------------------------
   ! [2.0]: Calculate correlations/accuracy:
   !-------------------------------------------------------------------------  

   do k = 1, kz
      do i = 1, iy

         phi_b1(2:jx-1) = phi(i,2:jx-1,k) - phi_u(i,2:jx-1,k)
            
         call da_correlation_coeff1d(phi_b1(2:jx-1), phi(i,2:jx-1,k), &
                                      corr_coeff(k,i), accuracy(k,i))
     
         ! write(58,*) corr_coeff(k,i), accuracy(k,i)
      end do
         
      phi_b2(2:iy-1,2:jx-1) = phi(2:iy-1,2:jx-1,k) - phi_u(2:iy-1,2:jx-1,k)
      call da_correlation_coeff2d(phi_b2(2:iy-1,2:jx-1), &
                                   phi(2:iy-1,2:jx-1,k), &
                                   corr_coef, accurac)

      write(unit=stdout, fmt='(i6,1pe9.2,1pe9.2)') &
            k, corr_coef, accurac

   end do

   !-------------------------------------------------------------------------
   ! [3.0]: Tidy up:
   !-------------------------------------------------------------------------  

   deallocate(phi_b1)
   deallocate(phi_b2)
   deallocate(corr_coeff)
   deallocate(accuracy)

   if (trace_use) call da_trace_exit("da_check_balance")

end subroutine da_check_balance


subroutine da_check_cvtovv_adjoint(grid, cv_size, xbx, be, cv, vv)

   !---------------------------------------------------------------------------
   ! Purpose: Test vtovv routine and adjoint for compatibility.
   !
   ! Method:  Standard adjoint test: < vv, vv > = < cv_adj, cv >.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   integer, intent(in)               :: cv_size ! Size of cv array.  
   type (xbx_type),intent(in)        :: xbx   ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be    ! background error structure.
   real, intent(in)                  :: cv(1:cv_size) ! control variable.
   type (vp_type), intent(inout)     :: vv    ! CV(i,j,m).

   real                              :: adj_par_lhs ! < Vv, Vv >
   real                              :: adj_par_rhs ! < cv_adj, cv >
   real                              :: adj_sum_lhs ! < Vv, Vv >
   real                              :: adj_sum_rhs ! < cv_adj, cv >
   real                              :: cv2(1:cv_size)! control variable.
   integer                           :: cv_size_tmp

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_check_cvtovv_adjoint")

   if (cv_options == 3 ) then
      write(unit=stdout, fmt='(/a,i2,a/)') 'cv_options =',cv_options, &
                     '   no da_check_cvtovv_adjoint check...'
      goto 1234
   end if

   write(unit=stdout, fmt='(/a/)') 'da_check_cvtovv_adjoint: Test Results:'
      
   !-------------------------------------------------------------------------
   ! [2.0] Perform Vp = U_v Vv transform:
   !-------------------------------------------------------------------------

   if (global) then
      call da_transform_vtovv_global(cv_size, xbx, be, cv, vv)
   else
      call da_transform_vtovv(grid, cv_size, be, cv, vv)
   end if

   !----------------------------------------------------------------------
   ! [3.0] Calculate LHS of adjoint test equation:
   !----------------------------------------------------------------------
   adj_par_lhs = sum(vv % v1(its:ite,jts:jte,1:be%v1%mz)**2) &
               + sum(vv % v2(its:ite,jts:jte,1:be%v2%mz)**2) &
               + sum(vv % v3(its:ite,jts:jte,1:be%v3%mz)**2) &
               + sum(vv % v4(its:ite,jts:jte,1:be%v4%mz)**2) &
               + sum(vv % v5(its:ite,jts:jte,1:be%v5%mz)**2) 

   if (be % ne > 0) then
!     adj_par_lhs = adj_par_lhs + sum(vv % alpha(its:ite,jts:jte,1:be%alpha%mz,1:be%ne)**2)
      adj_par_lhs = adj_par_lhs + sum(vv % alpha(its_int:ite_int,jts_int:jte_int,1:be%alpha%mz,1:be%ne)**2)
   end if

   !----------------------------------------------------------------------
   ! [4.0] Calculate RHS of adjoint test equation:
   !----------------------------------------------------------------------

   if (global) then
      call da_transform_vtovv_global_adj(cv_size, xbx, be, cv2, vv)
   else
      call da_transform_vtovv_adj(grid, cv_size, be, cv2, vv)
   end if

   cv_size_tmp = cv_size - be%cv%size_jp - be%cv%size_js - be%cv%size_jl
   adj_par_rhs = sum(cv(1:cv_size_tmp) * cv2(1:cv_size_tmp))

   !----------------------------------------------------------------------
   ! [5.0] Print output:
   !----------------------------------------------------------------------

   if (.not. global ) then
    if( num_procs == 1) then
      write(unit=stdout, fmt='(a,e22.14)') &
         'Single Domain: < Vv, Vv >     = ', adj_par_lhs, &
         'Single Domain: < cv_adj, cv > = ', adj_par_rhs
    else
      write(unit=stdout, fmt='(/a/,a/)')&
        'It is Multi Processor Run: ',&
        'For Single Domain: da_check_cvtovv_adjoint Test: Not Performed'
    endif
   end if

   adj_sum_lhs = wrf_dm_sum_real(adj_par_lhs)

   if (global) then
      adj_sum_rhs = adj_par_rhs
   else
      adj_sum_rhs = wrf_dm_sum_real(adj_par_rhs)
   end if  

   if (rootproc) then
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
           'Whole  Domain: < Vv, Vv >     = ', adj_sum_lhs, &
           'Whole  Domain: < cv_adj, cv > = ', adj_sum_rhs
   end if
      
   write(unit=stdout, fmt='(/a/)') &
      'da_check_cvtovv_adjoint: Test Finished.'

1234 continue

   if (trace_use) call da_trace_exit("da_check_cvtovv_adjoint")

end subroutine da_check_cvtovv_adjoint


subroutine da_check_vtox_adjoint(grid, cv_size, xbx, be, ep, cv1, vv, vp)

   !---------------------------------------------------------------------------
   ! Purpose: Test V to X routine and adjoint for compatibility.
   !
   ! Method:  Standard adjoint test: < x, x > = < v_adj, v >.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   integer, intent(in)             :: cv_size
   type (xbx_type),intent(in)      :: xbx    ! For header & non-grid arrays.
   type (be_type), intent(in)      :: be     ! background error structure.
   type (ep_type), intent(in)      :: ep     ! Ensemble perturbation structure.
   real, intent(in)                :: cv1(1:cv_size) ! control variable (local).
   type (vp_type), intent(inout)   :: vv     ! Grdipt/EOF CV.
   type (vp_type), intent(inout)   :: vp     ! Grdipt/level CV.

   real                   :: cv2(1:cv_size)    ! control variable (local).
   real                   :: adj_par_lhs ! < x, x >
   real                   :: adj_par_rhs ! < x, x >
   real                   :: adj_sum_lhs ! < x, x >
   real                   :: adj_sum_rhs ! < v_adj, v >

   if (trace_use) call da_trace_entry("da_check_vtox_adjoint")

   write(unit=stdout, fmt='(/a/)') &
      'da_check_vtox_adjoint: Adjoint Test Results:'

   !----------------------------------------------------------------------
   ! [1.0] Initialise:
   !----------------------------------------------------------------------

   cv2(:) = 0.0
      
   !----------------------------------------------------------------------
   ! [2.0] Perform x = U v transform:
   !----------------------------------------------------------------------

   call da_zero_x (grid%xa)

   call da_transform_vtox(grid,cv_size, xbx, be, ep, cv1, vv, vp)
   call da_transform_xtoxa(grid)

   !----------------------------------------------------------------------
   ! [3.0] Calculate LHS of adjoint test equation: 
   !----------------------------------------------------------------------

   adj_par_lhs = sum(grid%xa % u(its:ite, jts:jte, kts:kte)**2) / typical_u_rms**2 &
               + sum(grid%xa % v(its:ite, jts:jte, kts:kte)**2) / typical_v_rms**2 &     
               + sum(grid%xa % p(its:ite, jts:jte, kts:kte)**2) / typical_p_rms**2 &     
               + sum(grid%xa % t(its:ite, jts:jte, kts:kte)**2) / typical_t_rms**2 &     
               + sum(grid%xa % q(its:ite, jts:jte, kts:kte)**2) / typical_q_rms**2 &     
               + sum(grid%xa % rho(its:ite,jts:jte,kts:kte)**2)/ typical_rho_rms**2 & 
               + sum(grid%xa % psfc(its:ite, jts:jte)**2) / typical_p_rms**2             

   if (cv_options_hum == cv_options_hum_relative_humidity) then
      adj_par_lhs = adj_par_lhs &
              + sum(grid%xa % rh(its:ite, jts:jte, kts:kte)**2) / typical_rh_rms**2
   end if
!
   if (use_radar_rf .or. use_radar_rhv .or.  crtm_cloud ) then
      adj_par_lhs = adj_par_lhs &
             + sum(grid%xa % qcw(its:ite, jts:jte, kts:kte)**2)/typical_qcw_rms**2 &
             + sum(grid%xa % qrn(its:ite, jts:jte, kts:kte)**2)/typical_qrn_rms**2
      if ( cloud_cv_options /= 1 ) then
         adj_par_lhs = adj_par_lhs &
             + sum(grid%xa % qci(its:ite, jts:jte, kts:kte)**2)/typical_qci_rms**2 &
             + sum(grid%xa % qsn(its:ite, jts:jte, kts:kte)**2)/typical_qsn_rms**2 &
             + sum(grid%xa % qgr(its:ite, jts:jte, kts:kte)**2)/typical_qgr_rms**2
      end if
   end if

   if (use_radarobs) then
      adj_par_lhs = adj_par_lhs &
         + sum(grid%xa % wh (its:ite, jts:jte, kts:kte)**2)/typical_w_rms**2
   else
      adj_par_lhs = adj_par_lhs &
         + sum(grid%xa % w  (its:ite, jts:jte, kts:kte)**2)/typical_w_rms**2
   end if

   !-------------------------------------------------------------------------
   ! [4.0] Rescale input to adjoint routine:
   !-------------------------------------------------------------------------

   grid%xa % u(:,:,:) = grid%xa % u(:,:,:) / typical_u_rms**2
   grid%xa % v(:,:,:) = grid%xa % v(:,:,:) / typical_v_rms**2
   grid%xa % p(:,:,:) = grid%xa % p(:,:,:) / typical_p_rms**2
   grid%xa % t(:,:,:) = grid%xa % t(:,:,:) / typical_t_rms**2
   grid%xa % q(:,:,:) = grid%xa % q(:,:,:) / typical_q_rms**2
   grid%xa % rho(:,:,:) = grid%xa % rho(:,:,:) / typical_rho_rms**2

   grid%xa % psfc(:,:) = grid%xa % psfc(:,:) / typical_p_rms**2

   if (cv_options_hum == cv_options_hum_relative_humidity) then
      grid%xa % rh(:,:,:) = grid%xa % rh(:,:,:) / typical_rh_rms**2
   end if

!
   if (use_radar_rf .or. use_radar_rhv .or. crtm_cloud ) then
      grid%xa % qcw(:,:,:) = grid%xa % qcw(:,:,:) / typical_qcw_rms**2
      grid%xa % qrn(:,:,:) = grid%xa % qrn(:,:,:) / typical_qrn_rms**2
      if ( cloud_cv_options /= 1 ) then
         grid%xa % qci(:,:,:) = grid%xa % qci(:,:,:) / typical_qci_rms**2
         grid%xa % qsn(:,:,:) = grid%xa % qsn(:,:,:) / typical_qsn_rms**2
         grid%xa % qgr(:,:,:) = grid%xa % qgr(:,:,:) / typical_qgr_rms**2
      end if
   end if

   if (use_radarobs) then
      grid%xa %wh(:,:,:) = grid%xa %wh(:,:,:) / typical_w_rms**2
      grid%xa % w(:,:,:) = 0.0
   else
      grid%xa %w (:,:,:) = grid%xa %w (:,:,:) / typical_w_rms**2
   end if

   !-------------------------------------------------------------------------
   ! [5.0] Perform adjoint operation:
   !-------------------------------------------------------------------------

   call da_transform_xtoxa_adj(grid)
   call da_transform_vtox_adj(grid, cv_size, xbx, be, ep, vp, vv, cv2)

   !-------------------------------------------------------------------------
   ! [6.0] Calculate RHS of adjoint test equation:
   !-------------------------------------------------------------------------

   adj_par_rhs = sum(cv1(1:cv_size) * cv2(1:cv_size))

   !-------------------------------------------------------------------------
   ! [7.0] Print output:
   !-------------------------------------------------------------------------

   if (.not. global ) then
    if( num_procs == 1) then
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
         'Single Domain: < x, x >     = ', adj_par_lhs, &
         'Single Domain: < v_adj, v > = ', adj_par_rhs
    else
      write(unit=stdout, fmt='(/a/,a/)')&
        'It is Multi Processor Run: ',&
        'For Single Domain: da_check_vtox_adjoint Test: Not Performed'
    endif
   end if

   adj_sum_lhs = wrf_dm_sum_real(adj_par_lhs)

   if (global) then
      adj_sum_rhs = adj_par_rhs
   else
      adj_sum_rhs = wrf_dm_sum_real(adj_par_rhs)
   end if

   if (rootproc) then
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
         'Whole  Domain: < x, x >     = ', adj_sum_lhs, &
         'Whole  Domain: < v_adj, v > = ', adj_sum_rhs
   end if

   write(unit=stdout, fmt='(/a/)') 'da_check_vtox_adjoint: Finished'

   if (trace_use) call da_trace_exit("da_check_vtox_adjoint")

end subroutine da_check_vtox_adjoint


subroutine da_check_vptox_adjoint(grid, ne, be, ep, vp, cv_size)

   !---------------------------------------------------------------------------
   ! Purpose: Test Vp to X routine and adjoint for compatibility.
   !
   ! Method:  Standard adjoint test: < x, x > = < v_adj, v >.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)     :: grid

   integer, intent(in)              :: ne   ! Ensemble size.
   type (be_type), intent(in)       :: be   ! background errors.
   type (ep_type), intent(in)       :: ep   ! Ensemble perturbation type.
   type (vp_type),intent(inout)     :: vp   ! grdipt/level cv (local)
   integer, intent(in)              :: cv_size

! control variable
   real                             :: cv(1:cv_size), cv_2(1:cv_size) 

   real                             :: adj_par_lhs ! < x, x >
   real                             :: adj_par_rhs ! < v_adj, v >
   real                             :: adj_sum_lhs ! < x, x >
   real                             :: adj_sum_rhs ! < v_adj, v >
   real                             :: vp2_v1(ims:ime,jms:jme,kms:kme)
   real                             :: vp2_v2(ims:ime,jms:jme,kms:kme)
   real                             :: vp2_v3(ims:ime,jms:jme,kms:kme)
   real                             :: vp2_v4(ims:ime,jms:jme,kms:kme)
   real                             :: vp2_v5(ims:ime,jms:jme,kms:kme)
!  real                             :: vp2_alpha(ims:ime,jms:jme,kms:kme,1:ne)
   real                             :: vp2_alpha(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int,1:ne)

   if (trace_use) call da_trace_entry("da_check_vptox_adjoint")

   !-------------------------------------------------------------------------
   ! [1.0] Initialise:
   !-------------------------------------------------------------------------


   call da_zero_x(grid%xa)

   vp2_v1(:,:,:) = vp % v1(:,:,:)
   vp2_v2(:,:,:) = vp % v2(:,:,:)

   call da_psichi_to_uv(vp % v1, vp % v2, grid%xb % coefx, &
                        grid%xb % coefy , grid%xa % u, grid%xa % v)

   adj_par_lhs = sum(grid%xa % u(its:ite,jts:jte,:)**2) / typical_u_rms**2
   adj_par_lhs = sum(grid%xa % v(its:ite,jts:jte,:)**2) / typical_v_rms**2 + &
      adj_par_lhs

   grid%xa % u(:,:,:) = grid%xa % u(:,:,:) / typical_u_rms**2
   grid%xa % v(:,:,:) = grid%xa % v(:,:,:) / typical_v_rms**2

   vp%v1(:,:,:)=0.0
   vp%v2(:,:,:)=0.0

   call da_psichi_to_uv_adj(grid%xa % u, grid%xa % v, grid%xb % coefx,   &
                             grid%xb % coefy, vp % v1, vp % v2)

   adj_par_rhs = sum(vp % v1(its:ite,jts:jte,:) * vp2_v1(its:ite,jts:jte,:))
   adj_par_rhs = sum(vp % v2(its:ite,jts:jte,:) * vp2_v2(its:ite,jts:jte,:)) + &
      adj_par_rhs
   
      write(unit=stdout, fmt='(/a/)') &
          'da_check_da_psichi_to_uv_adjoint: Test Results:'
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
          'Single Domain: < u_v,     u_v         > = ', adj_par_lhs, &
          'Single Domain: < psi_chi, psi_chi_adj > = ', adj_par_rhs

   adj_sum_lhs = wrf_dm_sum_real(adj_par_lhs)
   adj_sum_rhs = wrf_dm_sum_real(adj_par_rhs)

   if (rootproc) then

      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
          'Whole  Domain: < u_v,     u_v         > = ', adj_sum_lhs, &
          'Whole  Domain: < psi_chi, psi_chi_adj > = ', adj_sum_rhs
   end if
      write(unit=stdout, fmt='(/a/)') &
          'da_check_da_psichi_to_uv_adjoint: Test Finished:'

   vp%v1(:,:,:) = vp2_v1(:,:,:)
   vp%v2(:,:,:) = vp2_v2(:,:,:)

   call da_zero_x(grid%xa)

   vp2_v1(:,:,:) = vp % v1(:,:,:)
   vp2_v2(:,:,:) = vp % v2(:,:,:)
   vp2_v3(:,:,:) = vp % v3(:,:,:)
   vp2_v4(:,:,:) = vp % v4(:,:,:)
   vp2_v5(:,:,:) = vp % v5(:,:,:)

   if (be % ne > 0) vp2_alpha(:,:,:,:) = vp % alpha(:,:,:,:)

   !-------------------------------------------------------------------------
   ! [2.0] Perform x = U vp transform:
   !-------------------------------------------------------------------------

   if ( cv_options == 3 ) then
      
      call random_number(cv(:))
      cv(:) = cv(:) - 0.5
      cv_2 = cv

      call da_apply_be( be, cv, vp, grid)
      call da_transform_bal( vp, be, grid)
 
   else

      call da_transform_vptox(grid, vp, be, ep)

   end if

   !-------------------------------------------------------------------------
   ! [3.0] Calculate LHS of adjoint test equation:
   !-------------------------------------------------------------------------

   !  grid%xa % u(:,:,:) = 0.0
   !  grid%xa % v(:,:,:) = 0.0
   !  grid%xa % t(:,:,:) = 0.0
   !  grid%xa % q(:,:,:) = 0.0
   !  grid%xa%psfc(:,:) = 0.0

   !  grid%xa % p(:,:,:) = 0.0
   !  grid%xa % rho(:,:,:) = 0.0
   !  grid%xa % w(:,:,:) = 0.0
   !  grid%xa % wh(:,:,:) = 0.0
   !  grid%xa % rh(:,:,:) = 0.0
   !  grid%xa % qt(:,:,:) = 0.0
   !  grid%xa % qcw(:,:,:) = 0.0
   !  grid%xa % qrn(:,:,:) = 0.0

   adj_par_lhs = sum(grid%xa%u(its:ite,jts:jte,:)**2)/typical_u_rms**2
   adj_par_lhs = sum(grid%xa%v(its:ite,jts:jte,:)**2)/typical_v_rms**2 + adj_par_lhs
   adj_par_lhs = sum(grid%xa%t(its:ite,jts:jte,:)**2)/typical_t_rms**2 + adj_par_lhs
   if ( (use_radar_rf .or. crtm_cloud) .and. (cloud_cv_options == 1) ) then
      adj_par_lhs = sum(grid%xa%qt(its:ite,jts:jte,:)**2)/typical_q_rms**2 + adj_par_lhs
   else
      adj_par_lhs = sum(grid%xa%q(its:ite,jts:jte,:)**2)/typical_q_rms**2 + adj_par_lhs
   end if
   adj_par_lhs = sum(grid%xa%psfc(its:ite,jts:jte)**2)/typical_p_rms**2 + adj_par_lhs

   adj_par_lhs = sum(grid%xa%p(its:ite,jts:jte,:)**2)/typical_p_rms**2 + adj_par_lhs
   adj_par_lhs = sum(grid%xa%rho(its:ite,jts:jte,:)**2)/typical_rho_rms**2 + &
      adj_par_lhs

   if (use_radarobs) then
      adj_par_lhs = adj_par_lhs &
                + sum(grid%xa % wh (its:ite, jts:jte, kts:kte)**2)/typical_w_rms**2
   else
      adj_par_lhs = adj_par_lhs &
                + sum(grid%xa % w  (its:ite, jts:jte, kts:kte)**2)/typical_w_rms**2
   end if

   if (cv_options_hum == cv_options_hum_relative_humidity) then
      adj_par_lhs = sum(grid%xa % rh(its:ite,jts:jte,:)**2) / &
         typical_rh_rms**2 + adj_par_lhs
   end if


   if (use_radar_rf .or. crtm_cloud) then
      if ( cloud_cv_options == 1 ) then
         adj_par_lhs = sum(grid%xa % qcw(its:ite,jts:jte,kts:kte)**2) / &
            typical_qcw_rms**2 + adj_par_lhs
         adj_par_lhs = sum(grid%xa % qrn(its:ite,jts:jte,kts:kte)**2) / &
            typical_qrn_rms**2 + adj_par_lhs
      end if
   end if

   !-------------------------------------------------------------------------
   ! [4.0] Rescale input to adjoint routine:
   !-------------------------------------------------------------------------
      
   grid%xa % u(:,:,:) = grid%xa % u(:,:,:) / typical_u_rms**2
   grid%xa % v(:,:,:) = grid%xa % v(:,:,:) / typical_v_rms**2
   grid%xa % t(:,:,:) = grid%xa % t(:,:,:) / typical_t_rms**2
   if ( (use_radar_rf .or. crtm_cloud) .and. (cloud_cv_options == 1) ) then
      grid%xa % qt(:,:,:) = grid%xa % qt(:,:,:) / typical_q_rms**2
   else
      grid%xa % q(:,:,:) = grid%xa % q(:,:,:) / typical_q_rms**2
   end if
   grid%xa%psfc(:,:) = grid%xa%psfc(:,:) / typical_p_rms**2
   
   grid%xa % p(:,:,:) = grid%xa % p(:,:,:) / typical_p_rms**2
   grid%xa % rho(:,:,:) = grid%xa % rho(:,:,:) / typical_rho_rms**2

   if (use_radarobs) then
      grid%xa %wh(:,:,:) = grid%xa %wh(:,:,:) / typical_w_rms**2
      grid%xa % w(:,:,:) = 0.0
   else
      grid%xa %w (:,:,:) = grid%xa %w (:,:,:) / typical_w_rms**2
   end if

   if (cv_options_hum == cv_options_hum_relative_humidity) then
      grid%xa % rh(:,:,:) = grid%xa % rh(:,:,:) / typical_rh_rms**2
   end if


   if (use_radar_rf .or. crtm_cloud) then
      if ( cloud_cv_options == 1 ) then
         grid%xa % qcw(:,:,:) = grid%xa % qcw(:,:,:) / typical_qcw_rms**2
         grid%xa % qrn(:,:,:) = grid%xa % qrn(:,:,:) / typical_qrn_rms**2
      end if
   end if
   
   !-------------------------------------------------------------------------
   ! [5.0] Perform adjoint operation:
   !-------------------------------------------------------------------------

   call da_zero_vp_type (vp)

   if (cv_options == 3 ) then

      cv = 0.0

      call da_transform_bal_adj( vp, be, grid)
      call da_apply_be_adj( be, cv, vp, grid)

   else

       call da_transform_vptox_adj(grid, vp, be, ep)

   end if

   !-------------------------------------------------------------------------
   ! [6.0] Calculate RHS of adjoint test equation:
   !-------------------------------------------------------------------------

   adj_par_rhs = sum(vp % v1(its:ite,jts:jte,:) * vp2_v1(its:ite,jts:jte,:))
   adj_par_rhs = sum(vp % v2(its:ite,jts:jte,:) * vp2_v2(its:ite,jts:jte,:)) + &
      adj_par_rhs
   adj_par_rhs = sum(vp % v3(its:ite,jts:jte,:) * vp2_v3(its:ite,jts:jte,:)) + &
      adj_par_rhs
   adj_par_rhs = sum(vp % v4(its:ite,jts:jte,:) * vp2_v4(its:ite,jts:jte,:)) + &
      adj_par_rhs
   adj_par_rhs = sum(vp % v5(its:ite,jts:jte,:) * vp2_v5(its:ite,jts:jte,:)) + &
      adj_par_rhs
   if (be % ne > 0) then
      adj_par_rhs = sum(vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,:) * &
         vp2_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,:)) + adj_par_rhs
   end if

   if ( cv_options == 3 ) adj_par_rhs = sum (cv_2*cv)

   !-------------------------------------------------------------------------
   ! [7.0] Print output:
   !-------------------------------------------------------------------------

   adj_sum_lhs = wrf_dm_sum_real(adj_par_lhs)
   adj_sum_rhs = wrf_dm_sum_real(adj_par_rhs)

   write(unit=stdout, fmt='(/a/)') 'da_check_vptox_adjoint: Test Results:'
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
         'Single Domain: < x, x >       = ', adj_par_lhs, &
         'Single Domain: < vp_adj, vp > = ', adj_par_rhs

   if (rootproc) then
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
         'Whole  Domain: < x, x >       = ', adj_sum_lhs, &
         'Whole  Domain: < vp_adj, vp > = ', adj_sum_rhs
   end if

   vp % v1(:,:,:) = vp2_v1(:,:,:)
   vp % v2(:,:,:) = vp2_v2(:,:,:)
   vp % v3(:,:,:) = vp2_v3(:,:,:)
   vp % v4(:,:,:) = vp2_v4(:,:,:)
   vp % v5(:,:,:) = vp2_v5(:,:,:)
   if (be % ne > 0) vp % alpha(:,:,:,:) = vp2_alpha(:,:,:,:)

   write(unit=stdout, fmt='(/a/)') 'da_check_vptox_adjoint: Test Finished:'

   if (trace_use) call da_trace_exit("da_check_vptox_adjoint")
      
end subroutine da_check_vptox_adjoint


subroutine da_check_vp_errors(vp1, vp2, ne, &
                               its,ite, jts,jte, kts,kte)

   !---------------------------------------------------------------------------
   ! Purpose: Test invertibility of transform to/from Vp or Vv
   !
   ! Method:  Perform statistics on differences in initial and final Vv or Vp
   !---------------------------------------------------------------------------

   implicit none

   type (vp_type), intent(in)     :: vp1         ! Test input
   type (vp_type), intent(in)     :: vp2         ! Test output.
   integer, intent(in)            :: ne          ! Ensemble size.
   integer, intent(in)            :: its,ite, jts,jte, kts,kte ! tile   dims.

   real                           :: inv_size    ! 1/size of array.
   real                           :: rms_fild    ! RMS of field.
   real                           :: rms_diff    ! RMS of differnce.

   real, dimension(its:ite, jts:jte, kts:kte) :: diff ! Difference
!  real                           :: diff_alpha(its:ite,jts:jte,kts:kte,1:ne)
   real                           :: diff_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne)
   real                           :: inv_size_ens    ! 1/size of ensemble array.

   if (trace_use) call da_trace_entry("da_check_vp_errors")

   inv_size = 1.0 / real((ite-its+1) * (jte-jts+1) * (kte-kts+1))
   inv_size_ens = 1.0 / real((ite_int-its_int+1) * (jte_int-jts_int+1) * (kte_int-kts_int+1)*ne)

   !-------------------------------------------------------------------------
   ! [1.0]: Check v1 differences:
   !-------------------------------------------------------------------------

   diff(its:ite,jts:jte, kts:kte) = vp2 % v1(its:ite,jts:jte,kts:kte) - &
                                    vp1 % v1(its:ite,jts:jte,kts:kte)

   rms_fild = sqrt(sum(vp1 % v1(its:ite, jts:jte,kts:kte) &
                       * vp1 % v1(its:ite, jts:jte,kts:kte)) * inv_size)
   rms_diff = sqrt(sum(diff(its:ite, jts:jte,kts:kte) &
                       * diff(its:ite, jts:jte,kts:kte)) * inv_size)
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' v1 is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') &
           ' v1 RMS error/RMS field = ', rms_diff/rms_fild
   end if      

   !-------------------------------------------------------------------------
   ! [2.0]: Check v2 differences:
   !-------------------------------------------------------------------------

   diff(its:ite,jts:jte, kts:kte) = vp2 % v2(its:ite,jts:jte,kts:kte) - &
                                    vp1 % v2(its:ite,jts:jte,kts:kte)

   rms_fild = sqrt(sum(vp1 % v2(its:ite, jts:jte,kts:kte) &
                       * vp1 % v2(its:ite, jts:jte,kts:kte)) * inv_size)
   rms_diff = sqrt(sum(diff(its:ite, jts:jte,kts:kte) &
                       * diff(its:ite, jts:jte,kts:kte)) * inv_size)
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' v2 is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') &
           ' v2 RMS error/RMS field = ', rms_diff/rms_fild
   end if      

   !-------------------------------------------------------------------------
   ! [3.0]: Check v3 differences:
   !-------------------------------------------------------------------------

   diff(its:ite,jts:jte, kts:kte) = vp2 % v3(its:ite,jts:jte,kts:kte) - &
                                    vp1 % v3(its:ite,jts:jte,kts:kte)

   rms_fild = sqrt(sum(vp1 % v3(its:ite, jts:jte,kts:kte) &
                       * vp1 % v3(its:ite, jts:jte,kts:kte)) * inv_size)
   rms_diff = sqrt(sum(diff(its:ite, jts:jte,kts:kte) &
                       * diff(its:ite, jts:jte,kts:kte)) * inv_size)
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' v3 is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') &
         ' v3 RMS error/RMS field = ', rms_diff/rms_fild
   end if      

   !-------------------------------------------------------------------------
   ! [4.0]: Check v4 differences:
   !-------------------------------------------------------------------------

   diff(its:ite,jts:jte, kts:kte) = vp2 % v4(its:ite,jts:jte,kts:kte) - &
                                    vp1 % v4(its:ite,jts:jte,kts:kte)

   rms_fild = sqrt(sum(vp1 % v4(its:ite, jts:jte,kts:kte) &
                       * vp1 % v4(its:ite, jts:jte,kts:kte)) * inv_size)
   rms_diff = sqrt(sum(diff(its:ite, jts:jte,kts:kte) &
                       * diff(its:ite, jts:jte,kts:kte)) * inv_size)
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' v4 is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') &
         ' v4 RMS error/RMS field = ', rms_diff/rms_fild
   end if
      
   !-------------------------------------------------------------------------
   ! [5.0]: Check v5 differences:
   !-------------------------------------------------------------------------

   inv_size = 1.0 / real((ite-its+1) * (jte-jts+1))

   diff(its:ite, jts:jte,kts:kte) = vp2 % v5(its:ite, jts:jte,kts:kte) - &
      vp1 % v5(its:ite, jts:jte,kts:kte)

   rms_fild = sqrt(sum(vp1 % v5(its:ite, jts:jte,kts:kte) * &
      vp1 % v5(its:ite, jts:jte,kts:kte)) * inv_size)
   rms_diff = sqrt(sum(diff(its:ite, jts:jte,1) * &
      diff(its:ite, jts:jte,1)) * inv_size)
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' v5 is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') &
         ' v5 RMS error/RMS field = ', rms_diff/rms_fild
   end if    
      
   !-------------------------------------------------------------------------
   ! [6.0]: Check alpha differences:
   !-------------------------------------------------------------------------
 
   ! changed all dimensions below to reflect ensemble tiles
   if (ne > 0) then
!     inv_size = 1.0 / real((ite-its+1) * (jte-jts+1) * ne)
      diff_alpha(its:ite,jts:jte,kts:kte,1:ne) = vp2 % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne) - &
                                         vp1 % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne)
      rms_fild = sqrt(sum(vp1 % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne) &
                          * vp1 % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne)) * inv_size_ens)
      rms_diff = sqrt(sum(diff_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne) &
                          * diff_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:ne)) * inv_size_ens)
     
      if (rms_fild /= 0.0) then
         write(unit=stdout, fmt='(a,1pe10.4)') ' alpha RMS error/RMS field = ',&
            rms_diff/rms_fild
      end if
   end if

   if (trace_use) call da_trace_exit("da_check_vp_errors")

end subroutine da_check_vp_errors


subroutine da_check_vvtovp_adjoint(grid, ne, xb, be, vv, vp)

   !---------------------------------------------------------------------------
   ! Purpose: Test Vv to Vp routine and adjoint for compatibility.
   !
   ! Method:  Standard adjoint test: < Vp, Vp > = < Vv_adj, Vv >.
   !---------------------------------------------------------------------------

   implicit none

   type (domain), intent(in)         :: grid
   integer, intent(in)               :: ne    ! Ensemble size.
   type (xb_type), intent(in)        :: xb    ! first guess (local).
   type (be_type), intent(in)        :: be    ! background error structure.
   type (vp_type), intent(inout)     :: vv    ! CV(i,j,m).
   type (vp_type), intent(inout)     :: vp    ! CV(i,j,k)

   real                              :: adj_par_lhs ! < x, x >
   real                              :: adj_par_rhs ! < v_adj, v >
   real                              :: adj_sum_lhs ! < x, x >
   real                              :: adj_sum_rhs ! < v_adj, v >

   real                              :: vv2_v1(ims:ime,jms:jme,kms:kme)
   real                              :: vv2_v2(ims:ime,jms:jme,kms:kme)
   real                              :: vv2_v3(ims:ime,jms:jme,kms:kme)
   real                              :: vv2_v4(ims:ime,jms:jme,kms:kme)
   real                              :: vv2_v5(ims:ime,jms:jme,kms:kme)
!  real                              :: vv2_alpha(ims:ime,jms:jme,kts:kte,1:ne)
   real                              :: vv2_alpha(ims_int:ime_int,jms_int:jme_int,kts_int:kte_int,1:ne)

   if (trace_use) call da_trace_entry("da_check_vvtovp_adjoint")

   if (cv_options == 3 ) then
      write(unit=stdout, fmt='(/a,i2,a/)') 'cv_options =',cv_options, &
                     '   no da_check_vvtovp_adjoint check...'
      goto 1235
   end if

   !----------------------------------------------------------------------
   ! [1.0] Initialise:
   !----------------------------------------------------------------------

   write(unit=stdout, fmt='(/a/)') 'da_check_vvtovp_adjoint: Test Results:'
      
   !----------------------------------------------------------------------
   ! [2.0] Perform Vp = U_v Vv transform:
   !----------------------------------------------------------------------

   call da_vertical_transform(grid, 'u', be, &
                               xb % vertical_inner_product, &
                               vv, vp)

   !----------------------------------------------------------------------
   ! [3.0] Calculate LHS of adjoint test equation:
   !----------------------------------------------------------------------

   adj_par_lhs = sum(vp % v1(its:ite,jts:jte,kts:kte)**2) * inv_typ_vp1_sumsq &
               + sum(vp % v2(its:ite,jts:jte,kts:kte)**2) * inv_typ_vp2_sumsq &
               + sum(vp % v3(its:ite,jts:jte,kts:kte)**2) * inv_typ_vp3_sumsq &
               + sum(vp % v4(its:ite,jts:jte,kts:kte)**2) * inv_typ_vp4_sumsq &
               + sum(vp % v5(its:ite,jts:jte,kts:kte)**2) * inv_typ_vp5_sumsq

   if (be % ne > 0) then
      adj_par_lhs = adj_par_lhs + &
!        sum(vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne)**2) * inv_typ_vpalpha_sumsq
         sum(vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)**2) * inv_typ_vpalpha_sumsq
   end if

   !----------------------------------------------------------------------
   ! [4.0] Rescale input to adjoint routine:
   !----------------------------------------------------------------------

   vp % v1(its:ite,jts:jte,kts:kte) = vp % v1(its:ite,jts:jte,kts:kte) * &
      inv_typ_vp1_sumsq
   vp % v2(its:ite,jts:jte,kts:kte) = vp % v2(its:ite,jts:jte,kts:kte) * &
      inv_typ_vp2_sumsq
   vp % v3(its:ite,jts:jte,kts:kte) = vp % v3(its:ite,jts:jte,kts:kte) * &
      inv_typ_vp3_sumsq
   vp % v4(its:ite,jts:jte,kts:kte) = vp % v4(its:ite,jts:jte,kts:kte) * &
      inv_typ_vp4_sumsq
   vp % v5(its:ite,jts:jte,kts:kte) = vp % v5(its:ite,jts:jte,kts:kte) * &
      inv_typ_vp5_sumsq

   if (be % ne > 0) then
!     vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne) = &
!        vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne) * inv_typ_vpalpha_sumsq
      vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) = &
         vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) * inv_typ_vpalpha_sumsq
   end if

   !----------------------------------------------------------------------
   ! [5.0] Perform adjoint operation:
   !----------------------------------------------------------------------

   vv2_v1(its:ite,jts:jte,1:be%v1%mz) = vv % v1(its:ite,jts:jte,1:be%v1%mz)
   vv2_v2(its:ite,jts:jte,1:be%v2%mz) = vv % v2(its:ite,jts:jte,1:be%v2%mz)
   vv2_v3(its:ite,jts:jte,1:be%v3%mz) = vv % v3(its:ite,jts:jte,1:be%v3%mz)
   vv2_v4(its:ite,jts:jte,1:be%v4%mz) = vv % v4(its:ite,jts:jte,1:be%v4%mz)
   vv2_v5(its:ite,jts:jte,1:be%v5%mz) = vv % v5(its:ite,jts:jte,1:be%v5%mz)

   if (be % ne > 0) then
!     vv2_alpha(its:ite,jts:jte,kts:kte,1:be%ne) = vv % alpha(its:ite,jts:jte,kts:kte,1:be%ne)
      vv2_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) = &
           vv % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)
   end if      

   call da_vertical_transform(grid, 'u_adj', be, &
                               xb % vertical_inner_product, &
                               vv, vp)

   !----------------------------------------------------------------------
   ! [6.0] Calculate RHS of adjoint test equation:
   !----------------------------------------------------------------------

   adj_par_rhs = 0.0
   if (be % v1 % mz > 0) &
      adj_par_rhs = sum(vv % v1(its:ite,jts:jte,1:be%v1%mz) * &
         vv2_v1(its:ite,jts:jte,1:be%v1%mz)) + adj_par_rhs
   if (be % v2 % mz > 0) &
      adj_par_rhs = sum(vv % v2(its:ite,jts:jte,1:be%v2%mz) * &
         vv2_v2(its:ite,jts:jte,1:be%v2%mz)) + adj_par_rhs
   if (be % v3 % mz > 0) &
      adj_par_rhs = sum(vv % v3(its:ite,jts:jte,1:be%v3%mz) * &
         vv2_v3(its:ite,jts:jte,1:be%v3%mz)) + adj_par_rhs
   if (be % v4 % mz > 0) &
      adj_par_rhs = sum(vv % v4(its:ite,jts:jte,1:be%v4%mz) * &
         vv2_v4(its:ite,jts:jte,1:be%v4%mz)) + adj_par_rhs
   if (be % v5 % mz == 1) &
      adj_par_rhs = sum(vv % v5(its:ite,jts:jte,1:be%v5%mz) * &
         vv2_v5(its:ite,jts:jte,1:be%v5%mz)) + adj_par_rhs
!  if (be % ne > 0) &
!     adj_par_rhs = sum(vv % alpha(its:ite,jts:jte,kts:kte,1:be%ne) * &
!        vv2_alpha(its:ite,jts:jte,kts:kte,1:be%ne)) + adj_par_rhs
   if (be % ne > 0) &
      adj_par_rhs = sum(vv % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) * &
         vv2_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)) + adj_par_rhs

   !----------------------------------------------------------------------
   ! [7.0] Print output:
   !----------------------------------------------------------------------

   write(unit=stdout, fmt='(a,1pe22.14)') &
        'Single domain < vp,     vp > = ', adj_par_lhs, &
        'Single domain < Vv_adj, Vv > = ', adj_par_rhs

   adj_sum_lhs = wrf_dm_sum_real(adj_par_lhs)
   adj_sum_rhs = wrf_dm_sum_real(adj_par_rhs)


   if (rootproc) then
      write(unit=stdout, fmt='(/)')
      write(unit=stdout, fmt='(a,1pe22.14)') &
         'Whole  Domain: < Vp, Vp >     = ', adj_sum_lhs, &
         'Whole  Domain: < Vv_adj, Vv > = ', adj_sum_rhs
   end if
      
   vv % v1(its:ite,jts:jte,1:be%v1%mz) = vv2_v1(its:ite,jts:jte,1:be%v1%mz)
   vv % v2(its:ite,jts:jte,1:be%v2%mz) = vv2_v2(its:ite,jts:jte,1:be%v2%mz)
   vv % v3(its:ite,jts:jte,1:be%v3%mz) = vv2_v3(its:ite,jts:jte,1:be%v3%mz)
   vv % v4(its:ite,jts:jte,1:be%v4%mz) = vv2_v4(its:ite,jts:jte,1:be%v4%mz)
   vv % v5(its:ite,jts:jte,1:be%v5%mz) = vv2_v5(its:ite,jts:jte,1:be%v5%mz)

   if (be % ne > 0) then
!     vv % alpha(its:ite,jts:jte,kts:kte,1:be%ne) = vv2_alpha(its:ite,jts:jte,kts:kte,1:be%ne)
      vv % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) = &
                 vv2_alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)
   end if

   write(unit=stdout, fmt='(/a/)') 'da_check_vvtovp_adjoint: Test Finished.'

1235 continue

   if (trace_use) call da_trace_exit("da_check_vvtovp_adjoint")

end subroutine da_check_vvtovp_adjoint


subroutine da_check_xtovptox_errors(xa, xa2_u, xa2_v, xa2_w, xa2_t, &
                                     xa2_p, xa2_q, xa2_rho, &
                                     xa2_qt, xa2_qcw, xa2_qrn)

   !---------------------------------------------------------------------------
   ! Purpose: Test invertibility of v = U^{-1} x followed by x = Uv.
   !
   !  Method:  Perform statistics on differences in initial and final x.
   !---------------------------------------------------------------------------

   implicit none
      
   type (x_type), intent(in)      :: xa          ! Test input

   real, dimension(ims:ime, jms:jme, kms:kme), &
                 intent(in)      :: xa2_u, xa2_v, xa2_t, &
                                    xa2_p, xa2_q, xa2_rho, &
                                    xa2_qt, xa2_qcw, xa2_qrn
   real, dimension(ims:ime, jms:jme, kms:kme), &
                 intent(in)      :: xa2_w    !xiao


   real                           :: rms_fild    ! RMS of field.
   real                           :: rms_diff    ! RMS of differnce.

   real, dimension(ims:ime, jms:jme, kms:kme) :: diff ! Difference

   if (trace_use) call da_trace_entry("da_check_xtovpx_errors")

   !----------------------------------------------------------------------
   ! [1.0]: Check u differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte) = xa2_u(its:ite, jts:jte, kts:kte) &
                                   - xa% u(its:ite, jts:jte, kts:kte)
   
   rms_fild = sqrt(sum(xa % u(its:ite, jts:jte, kts:kte) &
                       * xa % u(its:ite, jts:jte, kts:kte)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte) &
                       * diff(its:ite, jts:jte, kts:kte)))
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' u is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' u RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' u RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' u RMS error/RMS field = ', &
         rms_diff/rms_fild
   end if        
     
   !----------------------------------------------------------------------
   ! [2.0]: Check v differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte) = xa2_v(its:ite, jts:jte, kts:kte) &
                                   - xa% v(its:ite, jts:jte, kts:kte)
   
   rms_fild = sqrt(sum(xa % v(its:ite, jts:jte, kts:kte) &
                       * xa % v(its:ite, jts:jte, kts:kte)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte) &
                       * diff(its:ite, jts:jte, kts:kte)))
     
   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' v is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' v RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' v RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' v RMS error/RMS field = ', &
      rms_diff/rms_fild
   end if    
      
   !----------------------------------------------------------------------
   ! [3.0]: Check t differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte) = xa2_t(its:ite, jts:jte, kts:kte) &
                                   - xa% t(its:ite, jts:jte, kts:kte)

   rms_fild = sqrt(sum(xa % t(its:ite, jts:jte, kts:kte) &
                       * xa % t(its:ite, jts:jte, kts:kte)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte) &
                       * diff(its:ite, jts:jte, kts:kte)))

   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' t is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' t RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' t RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' t RMS error/RMS field = ', &
         rms_diff/rms_fild
   end if         
        
   !----------------------------------------------------------------------
   ! [4.0]: Check p differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte) = xa2_p(its:ite, jts:jte, kts:kte) &
                                   - xa% p(its:ite, jts:jte, kts:kte)

   rms_fild = sqrt(sum(xa % p(its:ite, jts:jte, kts:kte) &
                       * xa % p(its:ite, jts:jte, kts:kte)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte) &
                       * diff(its:ite, jts:jte, kts:kte)))

   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' p is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' p RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' p RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' p RMS error/RMS field = ', &
         rms_diff/rms_fild
   end if           

   !----------------------------------------------------------------------
   ! [5.0]: Check q differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte) = xa2_q(its:ite, jts:jte, kts:kte) &
                                   - xa% q(its:ite, jts:jte, kts:kte)

   rms_fild = sqrt(sum(xa % q(its:ite, jts:jte, kts:kte) &
                       * xa % q(its:ite, jts:jte, kts:kte)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte) &
                       * diff(its:ite, jts:jte, kts:kte)))

   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' q is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' q RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' q RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' q RMS error/RMS field = ', &
         rms_diff/rms_fild
   end if        

   !----------------------------------------------------------------------
   ! [6.0]: Check rho differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte) = xa2_rho(its:ite, jts:jte, kts:kte) &
                                   - xa% rho(its:ite, jts:jte, kts:kte)

   rms_fild = sqrt(sum(xa % rho(its:ite, jts:jte, kts:kte) &
                       * xa % rho(its:ite, jts:jte, kts:kte)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte) &
                       * diff(its:ite, jts:jte, kts:kte)))

   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' rho is zero ' 
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' rho RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' rho RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' rho RMS error/RMS field = ', &
         rms_diff/rms_fild
   end if        

   !----------------------------------------------------------------------
   ! [7.0]: Check w differences:
   !----------------------------------------------------------------------

   diff(its:ite, jts:jte, kts:kte+1) = xa2_w(its:ite, jts:jte, kts:kte+1) &
                                     - xa% w(its:ite, jts:jte, kts:kte+1)

   rms_fild = sqrt(sum(xa % w(its:ite, jts:jte, kts:kte+1) &
                       * xa % w(its:ite, jts:jte, kts:kte+1)))
   rms_diff = sqrt(sum(diff(its:ite, jts:jte, kts:kte+1) &
                       * diff(its:ite, jts:jte, kts:kte+1)))

   if (rms_fild == 0.0) then
      write(unit=stdout, fmt='(a)') ' w is zero '
   else
      write(unit=stdout, fmt='(a,1pe10.4)') ' w RMS error = ', rms_diff
      write(unit=stdout, fmt='(a,1pe10.4)') ' w RMS field = ', rms_fild
      write(unit=stdout, fmt='(a,1pe10.4)') ' w RMS error/RMS field = ', &
         rms_diff/rms_fild
   end if

   if (trace_use) call da_trace_exit("da_check_xtovpx_errors")
         
end subroutine da_check_xtovptox_errors


subroutine da_check_xtoy_adjoint(cv_size, cv, xbx, be, grid, config_flags, iv, y)
   
   !--------------------------------------------------------------------------
   ! Purpose: Test observation operator transform and adjoint for compatibility.
   !
   ! Method:  Standard adjoint test: < y, y > = < x, x_adj >.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------------
   
   implicit none
   
   integer, intent(in)                       :: cv_size ! Size of cv array.
   type (be_type),             intent(in)    :: be    ! background error structure.
   real, intent(inout)                       :: cv(1:cv_size)   ! control variables.
   type (xbx_type),            intent(inout) :: xbx   ! Header & non-gridded vars.
   type (domain),              intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   type (iv_type),             intent(inout) :: iv    ! ob. increment vector.
   type (y_type),              intent(inout) :: y     ! y = h (grid%xa)

   real                           :: adj_ttl_lhs   ! < y, y >
   real                           :: adj_ttl_rhs   ! < x, x_adj >

   real                           :: partial_lhs   ! < y, y >
   real                           :: partial_rhs   ! < x, x_adj >

   real                           :: pertile_lhs   ! < y, y >
   real                           :: pertile_rhs   ! < x, x_adj >
 
   real, dimension(ims:ime, jms:jme, kms:kme) :: xa2_u, xa2_v, xa2_t, &
                                                 xa2_p, xa2_q, xa2_rh
   real, dimension(ims:ime, jms:jme, kms:kme) :: xa2_w
   real, dimension(ims:ime, jms:jme)          :: xa2_psfc
   real, dimension(ims:ime, jms:jme, kms:kme) :: xa2_qcw, xa2_qci, xa2_qrn, xa2_qsn, xa2_qgr
   real, dimension(ims:ime, jms:jme, kms:kme) :: x6a2_u, x6a2_v, x6a2_t, &
                                                 x6a2_p, x6a2_q, x6a2_rh
   real, dimension(ims:ime, jms:jme, kms:kme) :: x6a2_w
   real, dimension(ims:ime, jms:jme)          :: x6a2_psfc
   real, dimension(ims:ime, jms:jme, kms:kme) :: x6a2_qcw, x6a2_qci, x6a2_qrn, x6a2_qsn, x6a2_qgr
   real, dimension(:,:,:), allocatable :: a_hr_rainc, a_hr_rainnc
   integer :: nobwin, i, j, k, fgat_rain
   character(len=4) :: filnam
   character(len=256) :: timestr
   integer :: time_step_seconds
   type(x_type) :: shuffle
   real             :: subarea, whole_area

   if (trace_use) call da_trace_entry("da_check_xtoy_adjoint")

   write (unit=stdout, fmt='(/a/)') 'da_check_xtoy_adjoint: Test Results:'

   !----------------------------------------------------------------------
   ! [1.0] Initialise:
   !----------------------------------------------------------------------

   partial_lhs = 0.0
   pertile_lhs = 0.0


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

   xa2_u(ims:ime, jms:jme, kms:kme) = grid%xa%u(ims:ime, jms:jme, kms:kme)
   xa2_v(ims:ime, jms:jme, kms:kme) = grid%xa%v(ims:ime, jms:jme, kms:kme)
   xa2_t(ims:ime, jms:jme, kms:kme) = grid%xa%t(ims:ime, jms:jme, kms:kme)
   xa2_p(ims:ime, jms:jme, kms:kme) = grid%xa%p(ims:ime, jms:jme, kms:kme)
   xa2_q(ims:ime, jms:jme, kms:kme) = grid%xa%q(ims:ime, jms:jme, kms:kme)
   xa2_w(ims:ime, jms:jme, kms:kme) = grid%xa%w(ims:ime, jms:jme, kms:kme)
   xa2_rh(ims:ime, jms:jme, kms:kme)= grid%xa%rh(ims:ime, jms:jme, kms:kme)
   xa2_psfc(ims:ime, jms:jme)       = grid%xa%psfc(ims:ime, jms:jme)

   xa2_qcw(ims:ime, jms:jme, kms:kme) = grid%xa%qcw(ims:ime, jms:jme, kms:kme)
   xa2_qci(ims:ime, jms:jme, kms:kme) = grid%xa%qci(ims:ime, jms:jme, kms:kme)
   xa2_qrn(ims:ime, jms:jme, kms:kme) = grid%xa%qrn(ims:ime, jms:jme, kms:kme)
   xa2_qsn(ims:ime, jms:jme, kms:kme) = grid%xa%qsn(ims:ime, jms:jme, kms:kme)
   xa2_qgr(ims:ime, jms:jme, kms:kme) = grid%xa%qgr(ims:ime, jms:jme, kms:kme)

   x6a2_u = 0.0
   x6a2_v = 0.0
   x6a2_t = 0.0
   x6a2_p = 0.0
   x6a2_q = 0.0
   x6a2_w = 0.0
   x6a2_rh = 0.0
   x6a2_psfc = 0.0

   x6a2_qcw = 0.0
   x6a2_qci = 0.0
   x6a2_qrn = 0.0
   x6a2_qsn = 0.0
   x6a2_qgr = 0.0

   if (var4d) then
      write(unit=message(1),fmt='(A)')'Please recompile the code with 4dvar option' 
      call da_error("da_check_xtoy_adjoint.inc",205,message(1:1))
   end if

   if ( num_fgat_time > 1 ) then
      call domain_clock_get (grid, stop_timestr=timestr)
      call domain_clock_set( grid, current_timestr=timestr )
      call domain_clock_set (grid, time_step_seconds=-1*var4d_bin)
      call domain_clockprint(150, grid, 'get CurrTime from clock,')
   endif

   fgat_rain = num_fgat_time
   do nobwin= num_fgat_time, 1, -1

      iv%time = nobwin
      iv%info(:)%n1 = iv%info(:)%plocal(iv%time-1) + 1
      iv%info(:)%n2 = iv%info(:)%plocal(iv%time)

      if (var4d) then
      end if

      call da_pt_to_rho_lin(grid)
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

      if (sfc_assi_options == 2) then
         call da_transform_xtowtq (grid)
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
      end if

      if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. &
          use_gpsztdobs .or. use_gpsrefobs .or.                 &
          use_ssmitbobs .or. use_ssmiretrievalobs) then

         ! Now do something for PW
         call da_transform_xtotpw(grid)

         ! GPS Refractivity:
         if (use_gpsrefobs .or. use_gpsztdobs) then
            call da_transform_xtogpsref_lin(grid)
            if (use_gpsztdobs) call da_transform_xtoztd_lin(grid)
         end if

         if (use_ssmt1obs .or. use_ssmt2obs .or. &
             use_ssmitbobs .or. use_ssmiretrievalobs) then
            if (global) then
              call da_error("da_check_xtoy_adjoint.inc",270, &
                (/"grid%xb%speed is not available, see da_transfer_kmatoxb.inc"/))
            end if
            call da_transform_xtoseasfcwind_lin(grid)
         end if

         if (use_ssmitbobs) call da_transform_xtotb_lin (grid)

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
      end if

   ! Compute w increments using Richardson's eqn.

   if ( Use_RadarObs ) then
      if ( .not. var4d ) call da_uvprho_to_w_lin(grid)

      do k=kts,kte
         do j=jts,jte
            do i=its,ite
               grid%xa%wh(i,j,k)=0.5*(grid%xa%w(i,j,k)+grid%xa%w(i,j,k+1))
            end do
         end do
      end do

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

   if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) ) then
    if ( cloud_cv_options == 1 )then
      ! Partition of hydrometeor increments via warm rain process
      call da_moist_phys_lin(grid)
    end if
   end if

      !----------------------------------------------------------------------
      ! [2.0] Perform y = Hx transform:
      !----------------------------------------------------------------------
      call da_transform_xtoy (cv_size, cv, grid, iv, y)


      !----------------------------------------------------------------------
      ! [3.0] Calculate LHS of adjoint test equation and
      !       Rescale input to adjoint routine :
      !----------------------------------------------------------------------

      if (iv%info(sound)%nlocal > 0) call da_check_xtoy_adjoint_sound(iv, y, partial_lhs, pertile_lhs)
      if (iv%info(sonde_sfc)%nlocal > 0) call da_check_xtoy_adjoint_sonde_sfc (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(mtgirs)%nlocal   > 0) call da_check_xtoy_adjoint_mtgirs   (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(tamdar)%nlocal   > 0) call da_check_xtoy_adjoint_tamdar   (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(tamdar_sfc)%nlocal   > 0) call da_check_xtoy_adjoint_tamdar_sfc(iv, y, partial_lhs, pertile_lhs)
      if (iv%info(synop)%nlocal    > 0) call da_check_xtoy_adjoint_synop    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(geoamv)%nlocal   > 0) call da_check_xtoy_adjoint_geoamv   (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(polaramv)%nlocal > 0) call da_check_xtoy_adjoint_polaramv (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(airep)%nlocal    > 0) call da_check_xtoy_adjoint_airep    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(pilot)%nlocal    > 0) call da_check_xtoy_adjoint_pilot    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(radar)%nlocal    > 0) call da_check_xtoy_adjoint_radar    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(satem)%nlocal    > 0) call da_check_xtoy_adjoint_satem    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(metar)%nlocal    > 0) call da_check_xtoy_adjoint_metar    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(ships)%nlocal    > 0) call da_check_xtoy_adjoint_ships    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(gpspw)%nlocal    > 0) call da_check_xtoy_adjoint_gpspw    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(gpsref)%nlocal   > 0) call da_check_xtoy_adjoint_gpsref   (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(ssmi_tb)%nlocal  > 0) call da_check_xtoy_adjoint_ssmi_tb  (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(ssmi_rv)%nlocal  > 0) call da_check_xtoy_adjoint_ssmi_rv  (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(ssmt2)%nlocal    > 0) call da_check_xtoy_adjoint_ssmt1    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(ssmt2)%nlocal    > 0) call da_check_xtoy_adjoint_ssmt2    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(qscat)%nlocal    > 0) call da_check_xtoy_adjoint_qscat    (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(profiler)%nlocal > 0) call da_check_xtoy_adjoint_profiler (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(buoy)%nlocal     > 0) call da_check_xtoy_adjoint_buoy     (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(bogus)%nlocal    > 0) call da_check_xtoy_adjoint_bogus    (iv, y, partial_lhs, pertile_lhs)
      if (iv%num_inst              > 0) call da_check_xtoy_adjoint_rad      (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(rain)%nlocal     > 0) call da_check_xtoy_adjoint_rain     (iv, y, partial_lhs, pertile_lhs)
      if (iv%info(pseudo)%nlocal   > 0) call da_check_xtoy_adjoint_pseudo   (iv, y, partial_lhs, pertile_lhs)

      !----------------------------------------------------------------------
      ! [5.0] Perform adjoint operation:
      !----------------------------------------------------------------------
      call da_zero_x (grid%xa)

      if (use_rainobs .and. num_fgat_time > 1) then
         a_hr_rainc(:,:,nobwin) = 0.0
         a_hr_rainnc(:,:,nobwin) = 0.0
      endif


      call da_transform_xtoy_adj (cv_size, cv, grid, iv, y, grid%xa)


   ! Compute w increments using Richardson's eqn.
   if ( Use_RadarObs)  then
      do k=kts,kte
         do j=jts,jte
            do i=its,ite
               grid%xa%w(i,j,k)=grid%xa%w(i,j,k)+0.5*grid%xa%wh(i,j,k)
               grid%xa%w(i,j,k+1)=grid%xa%w(i,j,k+1)+0.5*grid%xa%wh(i,j,k)
               grid%xa%wh(i,j,k)=0.0
            end do
         end do
      end do

      if ( .not. var4d ) call da_uvprho_to_w_adj(grid)
   end if

   if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) ) then
     if ( cloud_cv_options == 1) then
      ! Partition of hydrometeor increments via warm rain process
      call da_moist_phys_adj(grid)
     end if
   end if

      if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. &
          use_gpsztdobs .or. use_gpsrefobs .or.                 &
          use_ssmitbobs .or. use_ssmiretrievalobs) then

         if (use_ssmitbobs) call da_transform_xtotb_adj (grid)

         ! for PW
         call da_transform_xtotpw_adj (grid)

         ! GPS Refractivity:
         if (use_gpsrefobs .or. use_gpsztdobs) then
            if (use_gpsztdobs) call da_transform_xtoztd_adj(grid)
            call da_transform_xtogpsref_adj (grid)
         end if

         if (use_ssmt1obs .or. use_ssmt2obs .or. &
             use_ssmitbobs .or. use_ssmiretrievalobs) then
            if (global) then
               call da_error("da_check_xtoy_adjoint.inc",416, &
                  (/"grid%xb%speed is not available, see da_transfer_kmatoxb.inc"/))
            end if
            call da_transform_xtoseasfcwind_adj (grid)
         end if
      end if

      ! Now do something for surface variables
      if (sfc_assi_options == 2) then
         call da_transform_xtowtq_adj (grid)

      end if

      call da_pt_to_rho_adj (grid)

      if (var4d) then
      end if

      if ( nobwin > 1 ) call domain_clockadvance (grid)
      call domain_clockprint(150, grid, 'DEBUG Adjoint Check:  get CurrTime from clock,')

   end do

   if ( num_fgat_time > 1 ) then
      call nl_get_time_step ( grid%id, time_step_seconds)
      call domain_clock_set (grid, time_step_seconds=time_step_seconds)
      call domain_clockprint(150, grid, 'get CurrTime from clock,')
   endif

   if (var4d) then
   end if


   pertile_rhs = sum (grid%xa%u(ims:ime, jms:jme, kms:kme) * xa2_u(ims:ime, jms:jme, kms:kme)) &
      + sum (grid%xa%v(ims:ime, jms:jme, kms:kme) * xa2_v(ims:ime, jms:jme, kms:kme))          &
      + sum (grid%xa%w(ims:ime, jms:jme, kms:kme) * xa2_w(ims:ime, jms:jme, kms:kme))          &
      + sum (grid%xa%t(ims:ime, jms:jme, kms:kme) * xa2_t(ims:ime, jms:jme, kms:kme))          &
      + sum (grid%xa%p(ims:ime, jms:jme, kms:kme) * xa2_p(ims:ime, jms:jme, kms:kme))          &
      + sum (grid%xa%q(ims:ime, jms:jme, kms:kme) * xa2_q(ims:ime, jms:jme, kms:kme))          &
      + sum (grid%xa%rh(ims:ime, jms:jme, kms:kme)* xa2_rh(ims:ime, jms:jme, kms:kme))         &
      + sum (grid%xa%psfc(ims:ime, jms:jme) * xa2_psfc(ims:ime, jms:jme))                      &
      + sum (grid%x6a%u(ims:ime, jms:jme, kms:kme) * x6a2_u(ims:ime, jms:jme, kms:kme))        &
      + sum (grid%x6a%v(ims:ime, jms:jme, kms:kme) * x6a2_v(ims:ime, jms:jme, kms:kme))        &
      + sum (grid%x6a%w(ims:ime, jms:jme, kms:kme) * x6a2_w(ims:ime, jms:jme, kms:kme))        &
      + sum (grid%x6a%t(ims:ime, jms:jme, kms:kme) * x6a2_t(ims:ime, jms:jme, kms:kme))        &
      + sum (grid%x6a%p(ims:ime, jms:jme, kms:kme) * x6a2_p(ims:ime, jms:jme, kms:kme))        &
      + sum (grid%x6a%q(ims:ime, jms:jme, kms:kme) * x6a2_q(ims:ime, jms:jme, kms:kme))        &
      + sum (grid%x6a%rh(ims:ime, jms:jme, kms:kme)* x6a2_rh(ims:ime, jms:jme, kms:kme))       &
      + sum (grid%x6a%psfc(ims:ime, jms:jme) * x6a2_psfc(ims:ime, jms:jme))
   pertile_rhs = pertile_rhs &
      + sum (grid%xa%qcw(ims:ime, jms:jme, kms:kme) * xa2_qcw(ims:ime, jms:jme, kms:kme))      &
      + sum (grid%xa%qci(ims:ime, jms:jme, kms:kme) * xa2_qci(ims:ime, jms:jme, kms:kme))      &
      + sum (grid%xa%qrn(ims:ime, jms:jme, kms:kme) * xa2_qrn(ims:ime, jms:jme, kms:kme))      &
      + sum (grid%xa%qsn(ims:ime, jms:jme, kms:kme) * xa2_qsn(ims:ime, jms:jme, kms:kme))      &
      + sum (grid%xa%qgr(ims:ime, jms:jme, kms:kme) * xa2_qgr(ims:ime, jms:jme, kms:kme))      &
      + sum (grid%x6a%qcw(ims:ime, jms:jme, kms:kme) * x6a2_qcw(ims:ime, jms:jme, kms:kme))    &
      + sum (grid%x6a%qci(ims:ime, jms:jme, kms:kme) * x6a2_qci(ims:ime, jms:jme, kms:kme))    &
      + sum (grid%x6a%qrn(ims:ime, jms:jme, kms:kme) * x6a2_qrn(ims:ime, jms:jme, kms:kme))    &
      + sum (grid%x6a%qsn(ims:ime, jms:jme, kms:kme) * x6a2_qsn(ims:ime, jms:jme, kms:kme))    &
      + sum (grid%x6a%qgr(ims:ime, jms:jme, kms:kme) * x6a2_qgr(ims:ime, jms:jme, kms:kme))


   !----------------------------------------------------------------------
   ! [6.0] Calculate RHS of adjoint test equation:
   !----------------------------------------------------------------------
   
   partial_rhs = sum (grid%xa%u(its:ite, jts:jte, kts:kte) * xa2_u(its:ite, jts:jte, kts:kte)) &
      + sum (grid%xa%v(its:ite, jts:jte, kts:kte) * xa2_v(its:ite, jts:jte, kts:kte))          &
      + sum (grid%xa%w(its:ite, jts:jte, kts:kte+1) * xa2_w(its:ite, jts:jte, kts:kte+1))      &
      + sum (grid%xa%t(its:ite, jts:jte, kts:kte) * xa2_t(its:ite, jts:jte, kts:kte))          &
      + sum (grid%xa%p(its:ite, jts:jte, kts:kte) * xa2_p(its:ite, jts:jte, kts:kte))          &
      + sum (grid%xa%q(its:ite, jts:jte, kts:kte) * xa2_q(its:ite, jts:jte, kts:kte))          &
      + sum (grid%xa%rh(its:ite, jts:jte, kts:kte)* xa2_rh(its:ite, jts:jte, kts:kte))         &
      + sum (grid%xa%psfc(its:ite, jts:jte) * xa2_psfc(its:ite, jts:jte))                      &
      + sum (grid%x6a%u(its:ite, jts:jte, kts:kte) * x6a2_u(its:ite, jts:jte, kts:kte))        &
      + sum (grid%x6a%v(its:ite, jts:jte, kts:kte) * x6a2_v(its:ite, jts:jte, kts:kte))        &
      + sum (grid%x6a%w(its:ite, jts:jte, kts:kte+1) * x6a2_w(its:ite, jts:jte, kts:kte+1))    &
      + sum (grid%x6a%t(its:ite, jts:jte, kts:kte) * x6a2_t(its:ite, jts:jte, kts:kte))        &
      + sum (grid%x6a%p(its:ite, jts:jte, kts:kte) * x6a2_p(its:ite, jts:jte, kts:kte))        &
      + sum (grid%x6a%q(its:ite, jts:jte, kts:kte) * x6a2_q(its:ite, jts:jte, kts:kte))        &
      + sum (grid%x6a%rh(its:ite, jts:jte, kts:kte)* x6a2_rh(its:ite, jts:jte, kts:kte))       &
      + sum (grid%x6a%psfc(its:ite, jts:jte) * x6a2_psfc(its:ite, jts:jte)) 

   partial_rhs = partial_rhs &
      + sum (grid%xa%qcw(its:ite, jts:jte, kts:kte) * xa2_qcw(its:ite, jts:jte, kts:kte))   &
      + sum (grid%xa%qci(its:ite, jts:jte, kts:kte) * xa2_qci(its:ite, jts:jte, kts:kte))   & 
      + sum (grid%xa%qrn(its:ite, jts:jte, kts:kte) * xa2_qrn(its:ite, jts:jte, kts:kte))   & 
      + sum (grid%xa%qsn(its:ite, jts:jte, kts:kte) * xa2_qsn(its:ite, jts:jte, kts:kte))   & 
      + sum (grid%xa%qgr(its:ite, jts:jte, kts:kte) * xa2_qgr(its:ite, jts:jte, kts:kte))   & 
      + sum (grid%x6a%qcw(its:ite, jts:jte, kts:kte) * x6a2_qcw(its:ite, jts:jte, kts:kte)) &
      + sum (grid%x6a%qci(its:ite, jts:jte, kts:kte) * x6a2_qci(its:ite, jts:jte, kts:kte)) &
      + sum (grid%x6a%qrn(its:ite, jts:jte, kts:kte) * x6a2_qrn(its:ite, jts:jte, kts:kte)) &
      + sum (grid%x6a%qsn(its:ite, jts:jte, kts:kte) * x6a2_qsn(its:ite, jts:jte, kts:kte)) &
      + sum (grid%x6a%qgr(its:ite, jts:jte, kts:kte) * x6a2_qgr(its:ite, jts:jte, kts:kte))


   !----------------------------------------------------------------------
   !  [7.0] Print output:
   !----------------------------------------------------------------------
   write (unit=stdout, fmt='(A,1pe22.14)') ' Single Domain < y, y     > = ', pertile_lhs
   write (unit=stdout, fmt='(A,1pe22.14)') ' Single Domain < x, x_adj > = ', pertile_rhs

   adj_ttl_lhs = wrf_dm_sum_real (partial_lhs)
   adj_ttl_rhs = wrf_dm_sum_real (partial_rhs)
   
   if (rootproc) then
      write(unit=stdout, fmt='(/)')
      write (unit=stdout, fmt='(A,1pe22.14)') ' Whole Domain < y, y     > = ', adj_ttl_lhs
      write (unit=stdout, fmt='(A,1pe22.14)') ' Whole Domain < x, x_adj > = ', adj_ttl_rhs
   end if

   write (unit=stdout, fmt='(/a/)') 'da_check_xtoy_adjoint: Test Finished:'
   if (trace_use) call da_trace_exit("da_check_xtoy_adjoint")
   
end subroutine da_check_xtoy_adjoint


subroutine da_check_xtoy_adjoint_airep(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n, k          ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_airep")

   do n=iv%info(airep)%n1, iv%info(airep)%n2
      do k=1, iv%info(airep)%levels(n)
         if (iv%info(airep)%proc_domain(k,n)) then
            adjtest_lhs = adjtest_lhs + &
               (y%airep(n)%u(k) / typical_u_rms)**2 + &
               (y%airep(n)%v(k) / typical_v_rms)**2 + &
               (y%airep(n)%t(k) / typical_t_rms)**2 + & 
               (y%airep(n)%q(k) / typical_q_rms)**2 
         end if
      end do

      do k=1, iv%info(airep)%levels(n)
         pertile_lhs = pertile_lhs + &
            (y%airep(n)%u(k) / typical_u_rms)**2 + &
            (y%airep(n)%v(k) / typical_v_rms)**2 + &
            (y%airep(n)%t(k) / typical_t_rms)**2 + &
            (y%airep(n)%q(k) / typical_q_rms)**2

        y%airep(n)%u(k) = y%airep(n)%u(k) / typical_u_rms ** 2
        y%airep(n)%v(k) = y%airep(n)%v(k) / typical_v_rms ** 2
        y%airep(n)%t(k) = y%airep(n)%t(k) / typical_t_rms ** 2
        y%airep(n)%q(k) = y%airep(n)%q(k) / typical_q_rms ** 2
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_airep")

end subroutine da_check_xtoy_adjoint_airep


subroutine da_check_xtoy_adjoint_gpspw(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_gpspw")

   do n=iv%info(gpspw)%n1, iv%info(gpspw)%n2
      if (iv%info(gpspw)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs + (y%gpspw(n) %tpw/typical_tpw_rms) ** 2
      end if

      pertile_lhs = pertile_lhs + (y%gpspw(n) %tpw/typical_tpw_rms) ** 2

      y%gpspw (n)%tpw = y%gpspw (n)%tpw/typical_tpw_rms ** 2
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_gpspw")

end subroutine da_check_xtoy_adjoint_gpspw


subroutine da_check_xtoy_adjoint_gpsref(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_gpsref")

   do n=iv%info(gpsref)%n1, iv%info(gpsref)%n2
      if (iv%info(gpsref)%proc_domain(1,n)) then
         do k=1, iv%info(gpsref)%levels(n)
            adjtest_lhs = adjtest_lhs + (y%gpsref(n)%ref(k) / typical_ref_rms)**2 
         end do
      end if

      do k=1, iv%info(gpsref)%levels(n)
         pertile_lhs = pertile_lhs + (y%gpsref(n)%ref(k) / typical_ref_rms)**2
         y%gpsref(n)%ref(k) = y%gpsref(n)%ref(k) / typical_ref_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_gpsref")

end subroutine da_check_xtoy_adjoint_gpsref


subroutine da_check_xtoy_adjoint_metar(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type(y_type) , intent(inout)  :: y             ! y = h (xa)
   real,          intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_metar")

   do n=iv%info(metar)%n1, iv%info(metar)%n2
      if (iv%info(metar)%zk(1,n) < 1.0 .and. sfc_assi_options /= 2) then
         y%metar(n)%u = 0.0
         y%metar(n)%v = 0.0
         y%metar(n)%t = 0.0
         y%metar(n)%p = 0.0
         y%metar(n)%q = 0.0

         cycle
      end if

      if ( sfc_assi_options == 2 ) then
          if (iv%metar(n)%u%qc < 0) y%metar(n)%u = 0.0
          if (iv%metar(n)%v%qc < 0) y%metar(n)%v = 0.0
          if (iv%metar(n)%t%qc < 0) y%metar(n)%t = 0.0
          if (iv%metar(n)%p%qc < 0) y%metar(n)%p = 0.0
          if (iv%metar(n)%q%qc < 0) y%metar(n)%q = 0.0
      end if

      y%metar(n)%u = y%metar(n)%u/typical_u_rms
      y%metar(n)%v = y%metar(n)%v/typical_v_rms
      y%metar(n)%t = y%metar(n)%t/typical_t_rms
      y%metar(n)%p = y%metar(n)%p/typical_p_rms
      y%metar(n)%q = y%metar(n)%q/typical_q_rms

      if (iv%info(metar)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs  &
                     + (y%metar(n)%u)**2 &
                     + (y%metar(n)%v)**2 &
                     + (y%metar(n)%t)**2 &
                     + (y%metar(n)%p)**2 &
                     + (y%metar(n)%q)**2
      end if

      pertile_lhs = pertile_lhs &
                  + (y%metar(n)%u)**2 &
                  + (y%metar(n)%v)**2 &
                  + (y%metar(n)%t)**2 &
                  + (y%metar(n)%p)**2 &
                  + (y%metar(n)%q)**2

      y%metar(n)%u = y%metar(n)%u/typical_u_rms
      y%metar(n)%v = y%metar(n)%v/typical_v_rms
      y%metar(n)%t = y%metar(n)%t/typical_t_rms
      y%metar(n)%p = y%metar(n)%p/typical_p_rms
      y%metar(n)%q = y%metar(n)%q/typical_q_rms
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_metar")

end subroutine da_check_xtoy_adjoint_metar


subroutine da_check_xtoy_adjoint_pilot(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs  

   integer  :: n, k          ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_pilot")

   do n=iv%info(pilot)%n1, iv%info(pilot)%n2
      if (iv%info(pilot)%proc_domain(1,n)) then
         do k=1, iv%info(pilot)%levels(n)
            adjtest_lhs = adjtest_lhs + &
                (y%pilot(n)%u(k)/typical_u_rms)**2 + (y%pilot(n)%v(k)/typical_v_rms)**2
         end do
      end if

      do k=1, iv%info(pilot)%levels(n)
         pertile_lhs = pertile_lhs + &
            (y%pilot(n)%u(k)/typical_u_rms)**2 + (y%pilot(n)%v(k)/typical_v_rms)**2

         y%pilot(n)%u(k)= y%pilot(n)%u(k) / typical_u_rms ** 2
         y%pilot(n)%v(k)= y%pilot(n)%v(k) / typical_v_rms ** 2
      end do
   end do

   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_pilot")

end subroutine da_check_xtoy_adjoint_pilot


subroutine da_check_xtoy_adjoint_ssmi_rv(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! loop counter.
   real    :: var

   if (trace_use) call da_trace_entry("da_check_xtoy_adjoint_ssmi_rv")

   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n=iv%info(ssmi_rv)%n1, iv%info(ssmi_rv)%n2
         y%ssmi_rv(n)%speed = y%ssmi_rv(n)%speed/typical_speed_rms
         y%ssmi_rv(n)%tpw   = y%ssmi_rv(n)%tpw/typical_tpw_rms

         var = (y%ssmi_rv(n)%speed) ** 2 + (y%ssmi_rv(n)%tpw) ** 2

         pertile_lhs = pertile_lhs + var

         if (iv%info(ssmi_rv)%proc_domain(1,n)) then
            adjtest_lhs = adjtest_lhs + var
         end if

         y%ssmi_rv(n)%speed = y%ssmi_rv(n)%speed/typical_speed_rms
         y%ssmi_rv(n)%tpw   = y%ssmi_rv(n)%tpw/typical_tpw_rms
      end do
   end if

   if (trace_use) call da_trace_exit("da_check_xtoy_adjoint_ssmi_rv")

end subroutine da_check_xtoy_adjoint_ssmi_rv


subroutine da_check_xtoy_adjoint_ssmi_tb (iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! loop counter.
   real    :: var

   if (trace_use) call da_trace_entry("da_check_xtoy_adjoint_ssmi_tb")

   do n=iv%info(ssmi_tb)%n1, iv%info(ssmi_tb)%n2
      y%ssmi_tb(n)%tb19v = y%ssmi_tb(n)%tb19v/typical_tb19v_rms
      y%ssmi_tb(n)%tb19h = y%ssmi_tb(n)%tb19h/typical_tb19h_rms
      y%ssmi_tb(n)%tb22v = y%ssmi_tb(n)%tb22v/typical_tb22v_rms
      y%ssmi_tb(n)%tb37v = y%ssmi_tb(n)%tb37v/typical_tb37v_rms
      y%ssmi_tb(n)%tb37h = y%ssmi_tb(n)%tb37h/typical_tb37h_rms
      y%ssmi_tb(n)%tb85v = y%ssmi_tb(n)%tb85v/typical_tb85v_rms
      y%ssmi_tb(n)%tb85h = y%ssmi_tb(n)%tb85h/typical_tb85h_rms

       var = (y%ssmi_tb(n)%tb19v) ** 2 &
           + (y%ssmi_tb(n)%tb19h) ** 2 &
           + (y%ssmi_tb(n)%tb22v) ** 2 &
           + (y%ssmi_tb(n)%tb37v) ** 2 &
           + (y%ssmi_tb(n)%tb37h) ** 2 &
           + (y%ssmi_tb(n)%tb85v) ** 2 &
           + (y%ssmi_tb(n)%tb85h) ** 2 

      pertile_lhs = pertile_lhs + var

      if (iv%info(ssmi_tb)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs + var
      end if
      y%ssmi_tb(n)%tb19v = y%ssmi_tb(n)%tb19v/typical_tb19v_rms
      y%ssmi_tb(n)%tb19h = y%ssmi_tb(n)%tb19h/typical_tb19h_rms
      y%ssmi_tb(n)%tb22v = y%ssmi_tb(n)%tb22v/typical_tb22v_rms
      y%ssmi_tb(n)%tb37v = y%ssmi_tb(n)%tb37v/typical_tb37v_rms
      y%ssmi_tb(n)%tb37h = y%ssmi_tb(n)%tb37h/typical_tb37h_rms
      y%ssmi_tb(n)%tb85v = y%ssmi_tb(n)%tb85v/typical_tb85v_rms
      y%ssmi_tb(n)%tb85h = y%ssmi_tb(n)%tb85h/typical_tb85h_rms
   end do

   if (trace_use) call da_trace_exit("da_check_xtoy_adjoint_ssmi_tb")

end subroutine da_check_xtoy_adjoint_ssmi_tb


subroutine da_check_xtoy_adjoint_satem(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_satem")

   do n=iv%info(satem)%n1, iv%info(satem)%n2
      if (iv%info(satem)%proc_domain(1,n)) then
         do k=1, iv%info(satem)%levels(n)
            adjtest_lhs = adjtest_lhs + (y%satem(n)%thickness(k)/typical_thickness_rms)**2
         end do
      end if

      do k=1, iv%info(satem)%levels(n)
         pertile_lhs = pertile_lhs + (y%satem(n)%thickness(k)/typical_thickness_rms)**2
         y%satem(n)%thickness(k) = y%satem(n)%thickness(k) / typical_thickness_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_satem")

end subroutine da_check_xtoy_adjoint_satem


subroutine da_check_xtoy_adjoint_geoamv (iv, y, adjtest_lhs, pertile_lhs)

   !-------------------------------------------------------------------------
   ! Purpose:  Adjoint Test for Geo. AMVs Obs
   !-------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type(y_type),  intent(inout)  :: y             ! y = h (xa)
   real,          intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_geoamv")

   do n=iv%info(geoamv)%n1, iv%info(geoamv)%n2
      if (iv%info(geoamv)%proc_domain(1,n)) then
         do k=1, iv%info(geoamv)%levels(n)
            adjtest_lhs = adjtest_lhs + &
                        (y%geoamv(n)%u(k)/typical_u_rms)**2 + &
                        (y%geoamv(n)%v(k)/typical_v_rms)**2
         end do
      end if

      do k=1, iv%info(geoamv)%levels(n)
         pertile_lhs = pertile_lhs + &
                     (y%geoamv(n)%u(k)/typical_u_rms)**2 + &
                     (y%geoamv(n)%v(k)/typical_v_rms)**2

         y%geoamv(n)%u(k)= y%geoamv(n)%u(k) / typical_u_rms ** 2
         y%geoamv(n)%v(k)= y%geoamv(n)%v(k) / typical_v_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_geoamv")

end subroutine da_check_xtoy_adjoint_geoamv


subroutine da_check_xtoy_adjoint_polaramv (iv, y, adjtest_lhs, pertile_lhs)

   !-------------------------------------------------------------------------
   ! Purpose:  Adjoint Test for Polar AMVs Obs
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n, k          ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_polaramv")

   do n=iv%info(polaramv)%n1, iv%info(polaramv)%n2
      do k=1, iv%info(polaramv)%levels(n)
         if (iv%info(polaramv)%proc_domain(k,n)) then
            adjtest_lhs = adjtest_lhs + (y%polaramv(n)%u(k)/typical_u_rms)**2 + (y%polaramv(n)%v(k)/typical_v_rms)**2
         end if
      end do

      do k=1, iv%info(polaramv)%levels(n)
         pertile_lhs = pertile_lhs + &
            (y%polaramv(n)%u(k)/typical_u_rms)**2 + (y%polaramv(n)%v(k)/typical_v_rms)**2

         y%polaramv(n)%u(k)= y%polaramv(n)%u(k) / typical_u_rms ** 2
         y%polaramv(n)%v(k)= y%polaramv(n)%v(k) / typical_v_rms ** 2
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_polaramv")

end subroutine da_check_xtoy_adjoint_polaramv


subroutine da_check_xtoy_adjoint_ships(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_ships")

   do n=iv%info(ships)%n1, iv%info(ships)%n2
      if (iv%info(ships)%zk(1,n) < 1.0 .and. sfc_assi_options == 1) then
         y%ships(n)%u = 0.0
         y%ships(n)%v = 0.0
         y%ships(n)%t = 0.0
         y%ships(n)%p = 0.0
         y%ships(n)%q = 0.0

         cycle
      end if
 
      if ( sfc_assi_options == 2 ) then
          if (iv%ships(n)%u%qc < 0) y%ships(n)%u = 0.0
          if (iv%ships(n)%v%qc < 0) y%ships(n)%v = 0.0
          if (iv%ships(n)%t%qc < 0) y%ships(n)%t = 0.0
          if (iv%ships(n)%p%qc < 0) y%ships(n)%p = 0.0
          if (iv%ships(n)%q%qc < 0) y%ships(n)%q = 0.0
      end if

      y%ships(n)%u = y%ships(n)%u/typical_u_rms
      y%ships(n)%v = y%ships(n)%v/typical_v_rms
      y%ships(n)%t = y%ships(n)%t/typical_t_rms
      y%ships(n)%p = y%ships(n)%p/typical_p_rms
      y%ships(n)%q = y%ships(n)%q/typical_q_rms

      if (iv%info(ships)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs  &
            + (y%ships(n)%u)**2 &
            + (y%ships(n)%v)**2 &
            + (y%ships(n)%t)**2 &
            + (y%ships(n)%p)**2 &
            + (y%ships(n)%q)**2
      end if

      pertile_lhs = pertile_lhs &
         + (y%ships(n)%u)**2 &
         + (y%ships(n)%v)**2 &
         + (y%ships(n)%t)**2 &
         + (y%ships(n)%p)**2 &
         + (y%ships(n)%q)**2

      y%ships(n)%u = y%ships(n)%u/typical_u_rms
      y%ships(n)%v = y%ships(n)%v/typical_v_rms
      y%ships(n)%t = y%ships(n)%t/typical_t_rms
      y%ships(n)%p = y%ships(n)%p/typical_p_rms
      y%ships(n)%q = y%ships(n)%q/typical_q_rms
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_ships")

end subroutine da_check_xtoy_adjoint_ships


subroutine da_check_xtoy_adjoint_radar(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs

   integer :: n, k          ! Loop counter.

   if (trace_use) call da_trace_entry("da_check_xtoy_adjoint_radar")

   do n=iv%info(radar)%n1, iv%info(radar)%n2
      if (iv%info(radar)%proc_domain(1,n)) then
         do k=1, iv%info(radar)%levels(n)
            adjtest_lhs = adjtest_lhs + &
               (y%radar(n)%rv(k)/typical_rv_rms)**2   + (y%radar(n)%rf(k)/typical_rf_rms)**2   + &
               (y%radar(n)%rrn(k)/typical_qrn_rms)**2 + (y%radar(n)%rsn(k)/typical_qsn_rms)**2 + &
               (y%radar(n)%rgr(k)/typical_qgr_rms)**2 + (y%radar(n)%rqv(k)/typical_q_rms)**2 
         end do
      end if

      do k=1, iv%info(radar)%levels(n)
         pertile_lhs = pertile_lhs + &
            (y%radar(n)%rv(k)/typical_rv_rms)**2   + (y%radar(n)%rf(k)/typical_rf_rms)**2 + &
            (y%radar(n)%rrn(k)/typical_qrn_rms)**2 + (y%radar(n)%rsn(k)/typical_qsn_rms)**2 + &
            (y%radar(n)%rgr(k)/typical_qgr_rms)**2 + (y%radar(n)%rqv(k)/typical_q_rms)**2

         y%radar(n)%rv(k)= y%radar(n)%rv(k) / typical_rv_rms ** 2
         y%radar(n)%rf(k)= y%radar(n)%rf(k) / typical_rf_rms ** 2
         y%radar(n)%rrn(k)= y%radar(n)%rrn(k) / typical_qrn_rms ** 2
         y%radar(n)%rsn(k)= y%radar(n)%rsn(k) / typical_qsn_rms ** 2
         y%radar(n)%rgr(k)= y%radar(n)%rgr(k) / typical_qgr_rms ** 2
         y%radar(n)%rqv(k)= y%radar(n)%rqv(k) / typical_q_rms ** 2
      end do
   end do

end subroutine da_check_xtoy_adjoint_radar


subroutine da_check_xtoy_adjoint_rain(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! loop counter.
   real    :: var

   if (trace_use) call da_trace_entry("da_check_xtoy_adjoint_rain")

   if (iv%info(rain)%nlocal > 0) then
      do n=iv%info(rain)%n1, iv%info(rain)%n2
         y%rain(n)%rain = y%rain(n)%rain/typical_rain_rms

         var = (y%rain(n)%rain) ** 2

         pertile_lhs = pertile_lhs + var

         if (iv%info(rain)%proc_domain(1,n)) then
            adjtest_lhs = adjtest_lhs + var
         end if

         y%rain(n)%rain = y%rain(n)%rain/typical_rain_rms
      end do
   end if

   if (trace_use) call da_trace_exit("da_check_xtoy_adjoint_rain")

end subroutine da_check_xtoy_adjoint_rain

subroutine da_check_xtoy_adjoint_bogus(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer  :: n, k          ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_bogus")

   do n=iv%info(bogus)%n1, iv%info(bogus)%n2
      if (iv%info(bogus)%proc_domain(1,n)) then
         do k=1, iv%info(bogus)%levels(n)
            adjtest_lhs = adjtest_lhs + &
                          (y%bogus(n)%u(k)/typical_u_rms)**2 + &
                          (y%bogus(n)%v(k)/typical_v_rms)**2 + &
                          (y%bogus(n)%t(k)/typical_t_rms)**2 + &
                          (y%bogus(n)%q(k)/typical_q_rms)**2 
         end do
         adjtest_lhs = adjtest_lhs + (y%bogus(n)%slp/typical_p_rms)**2
      end if

      do k=1, iv%info(bogus)%levels(n)
         pertile_lhs = pertile_lhs + &
                       (y%bogus(n)%u(k)/typical_u_rms)**2 + &
                       (y%bogus(n)%v(k)/typical_v_rms)**2 + &
                       (y%bogus(n)%t(k)/typical_t_rms)**2 + &
                       (y%bogus(n)%q(k)/typical_q_rms)**2

         y%bogus(n)%u(k) = y%bogus(n)%u(k) / typical_u_rms ** 2
         y%bogus(n)%v(k) = y%bogus(n)%v(k) / typical_v_rms ** 2
         y%bogus(n)%t(k) = y%bogus(n)%t(k) / typical_t_rms ** 2
         y%bogus(n)%q(k) = y%bogus(n)%q(k) / typical_q_rms ** 2
      end do
      pertile_lhs = pertile_lhs + (y%bogus(n)%slp/typical_p_rms)**2
      y%bogus(n)%slp = y%bogus(n)%slp / typical_p_rms ** 2
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_bogus")

end subroutine da_check_xtoy_adjoint_bogus


subroutine da_check_xtoy_adjoint_sound (iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv            ! obs. inc. vector (o-b).
   type(y_type) , intent(inout) :: y             ! y = h (xa)
   real,          intent(inout) :: adjtest_lhs, pertile_lhs   

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_sound")

   do n=iv%info(sound)%n1, iv%info(sound)%n2
      do k=1, iv%info(sound)%levels(n)
         if (iv%info(sound)%proc_domain(k,n)) then
            adjtest_lhs = adjtest_lhs + &
                          (y%sound(n)%u(k)/typical_u_rms)**2 + &
                          (y%sound(n)%v(k)/typical_v_rms)**2 + &
                          (y%sound(n)%t(k)/typical_t_rms)**2 + &
                          (y%sound(n)%q(k)/typical_q_rms)**2
         end if
      end do

      do k=1, iv%info(sound)%levels(n)
         pertile_lhs = pertile_lhs + &
                       (y%sound(n)%u(k)/typical_u_rms)**2 + &
                       (y%sound(n)%v(k)/typical_v_rms)**2 + &
                       (y%sound(n)%t(k)/typical_t_rms)**2 + &
                       (y%sound(n)%q(k)/typical_q_rms)**2

         y%sound(n)%u(k) = y%sound(n)%u(k) / typical_u_rms ** 2
         y%sound(n)%v(k) = y%sound(n)%v(k) / typical_v_rms ** 2
         y%sound(n)%t(k) = y%sound(n)%t(k) / typical_t_rms ** 2
         y%sound(n)%q(k) = y%sound(n)%q(k) / typical_q_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_sound")

end subroutine da_check_xtoy_adjoint_sound


subroutine da_check_xtoy_adjoint_sonde_sfc(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_sonde_sfc")

   do n=iv%info(sound)%n1, iv%info(sound)%n2
      if (iv%info(sound)%zk(1,n) < 1.0 .and. sfc_assi_options /= 2) then
         y%sonde_sfc(n)%u = 0.0
         y%sonde_sfc(n)%v = 0.0
         y%sonde_sfc(n)%t = 0.0
         y%sonde_sfc(n)%p = 0.0
         y%sonde_sfc(n)%q = 0.0

         cycle
      end if

      if ( sfc_assi_options == 2 ) then
          if (iv%sonde_sfc(n)%u%qc < 0) y%sonde_sfc(n)%u = 0.0
          if (iv%sonde_sfc(n)%v%qc < 0) y%sonde_sfc(n)%v = 0.0
          if (iv%sonde_sfc(n)%t%qc < 0) y%sonde_sfc(n)%t = 0.0
          if (iv%sonde_sfc(n)%p%qc < 0) y%sonde_sfc(n)%p = 0.0
          if (iv%sonde_sfc(n)%q%qc < 0) y%sonde_sfc(n)%q = 0.0
      end if

      y%sonde_sfc(n)%u = y%sonde_sfc(n)%u/typical_u_rms
      y%sonde_sfc(n)%v = y%sonde_sfc(n)%v/typical_v_rms
      y%sonde_sfc(n)%t = y%sonde_sfc(n)%t/typical_t_rms
      y%sonde_sfc(n)%p = y%sonde_sfc(n)%p/typical_p_rms
      y%sonde_sfc(n)%q = y%sonde_sfc(n)%q/typical_q_rms

      if (iv%info(sound)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs  &
                     + (y%sonde_sfc(n)%u)**2 &
                     + (y%sonde_sfc(n)%v)**2 &
                     + (y%sonde_sfc(n)%t)**2 &
                     + (y%sonde_sfc(n)%p)**2 &
                     + (y%sonde_sfc(n)%q)**2
      end if

      pertile_lhs = pertile_lhs &
                  + (y%sonde_sfc(n)%u)**2 &
                  + (y%sonde_sfc(n)%v)**2 &
                  + (y%sonde_sfc(n)%t)**2 &
                  + (y%sonde_sfc(n)%p)**2 &
                  + (y%sonde_sfc(n)%q)**2

      y%sonde_sfc(n)%u = y%sonde_sfc(n)%u/typical_u_rms
      y%sonde_sfc(n)%v = y%sonde_sfc(n)%v/typical_v_rms
      y%sonde_sfc(n)%t = y%sonde_sfc(n)%t/typical_t_rms
      y%sonde_sfc(n)%p = y%sonde_sfc(n)%p/typical_p_rms
      y%sonde_sfc(n)%q = y%sonde_sfc(n)%q/typical_q_rms
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_sonde_sfc")

end subroutine da_check_xtoy_adjoint_sonde_sfc


subroutine da_check_xtoy_adjoint_mtgirs (iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv            ! obs. inc. vector (o-b).
   type(y_type) , intent(inout) :: y             ! y = h (xa)
   real,          intent(inout) :: adjtest_lhs, pertile_lhs   

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_mtgirs")

   do n=iv%info(mtgirs)%n1, iv%info(mtgirs)%n2
      do k=1, iv%info(mtgirs)%levels(n)
         if (iv%info(mtgirs)%proc_domain(k,n)) then
            adjtest_lhs = adjtest_lhs + &
                          (y%mtgirs(n)%u(k)/typical_u_rms)**2 + &
                          (y%mtgirs(n)%v(k)/typical_v_rms)**2 + &
                          (y%mtgirs(n)%t(k)/typical_t_rms)**2 + &
                          (y%mtgirs(n)%q(k)/typical_q_rms)**2
         end if
      end do

      do k=1, iv%info(mtgirs)%levels(n)
         pertile_lhs = pertile_lhs + &
                       (y%mtgirs(n)%u(k)/typical_u_rms)**2 + &
                       (y%mtgirs(n)%v(k)/typical_v_rms)**2 + &
                       (y%mtgirs(n)%t(k)/typical_t_rms)**2 + &
                       (y%mtgirs(n)%q(k)/typical_q_rms)**2

         y%mtgirs(n)%u(k) = y%mtgirs(n)%u(k) / typical_u_rms ** 2
         y%mtgirs(n)%v(k) = y%mtgirs(n)%v(k) / typical_v_rms ** 2
         y%mtgirs(n)%t(k) = y%mtgirs(n)%t(k) / typical_t_rms ** 2
         y%mtgirs(n)%q(k) = y%mtgirs(n)%q(k) / typical_q_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_mtgirs")

end subroutine da_check_xtoy_adjoint_mtgirs


subroutine da_check_xtoy_adjoint_tamdar (iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv            ! obs. inc. vector (o-b).
   type(y_type) , intent(inout) :: y             ! y = h (xa)
   real,          intent(inout) :: adjtest_lhs, pertile_lhs   

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_tamdar")

   do n=iv%info(tamdar)%n1, iv%info(tamdar)%n2
      do k=1, iv%info(tamdar)%levels(n)
         if (iv%info(tamdar)%proc_domain(k,n)) then
            adjtest_lhs = adjtest_lhs + &
                          (y%tamdar(n)%u(k)/typical_u_rms)**2 + &
                          (y%tamdar(n)%v(k)/typical_v_rms)**2 + &
                          (y%tamdar(n)%t(k)/typical_t_rms)**2 + &
                          (y%tamdar(n)%q(k)/typical_q_rms)**2
         end if
      end do

      do k=1, iv%info(tamdar)%levels(n)
         pertile_lhs = pertile_lhs + &
                       (y%tamdar(n)%u(k)/typical_u_rms)**2 + &
                       (y%tamdar(n)%v(k)/typical_v_rms)**2 + &
                       (y%tamdar(n)%t(k)/typical_t_rms)**2 + &
                       (y%tamdar(n)%q(k)/typical_q_rms)**2

         y%tamdar(n)%u(k) = y%tamdar(n)%u(k) / typical_u_rms ** 2
         y%tamdar(n)%v(k) = y%tamdar(n)%v(k) / typical_v_rms ** 2
         y%tamdar(n)%t(k) = y%tamdar(n)%t(k) / typical_t_rms ** 2
         y%tamdar(n)%q(k) = y%tamdar(n)%q(k) / typical_q_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_tamdar")

end subroutine da_check_xtoy_adjoint_tamdar


subroutine da_check_xtoy_adjoint_tamdar_sfc(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer  :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_tamdar_sfc")

   do n=iv%info(tamdar_sfc)%n1, iv%info(tamdar_sfc)%n2
      if (iv%info(tamdar_sfc)%zk(1,n) < 1.0 .and. sfc_assi_options /= 2) then
         y%tamdar_sfc(n)%u = 0.0
         y%tamdar_sfc(n)%v = 0.0
         y%tamdar_sfc(n)%t = 0.0
         y%tamdar_sfc(n)%p = 0.0
         y%tamdar_sfc(n)%q = 0.0
         cycle
      end if

      if ( sfc_assi_options == 2 ) then
          if (iv%tamdar_sfc(n)%u%qc < 0) y%tamdar_sfc(n)%u = 0.0
          if (iv%tamdar_sfc(n)%v%qc < 0) y%tamdar_sfc(n)%v = 0.0
          if (iv%tamdar_sfc(n)%t%qc < 0) y%tamdar_sfc(n)%t = 0.0
          if (iv%tamdar_sfc(n)%p%qc < 0) y%tamdar_sfc(n)%p = 0.0
          if (iv%tamdar_sfc(n)%q%qc < 0) y%tamdar_sfc(n)%q = 0.0
      end if

      y%tamdar_sfc(n)%u = y%tamdar_sfc(n)%u/typical_u_rms
      y%tamdar_sfc(n)%v = y%tamdar_sfc(n)%v/typical_v_rms
      y%tamdar_sfc(n)%t = y%tamdar_sfc(n)%t/typical_t_rms
      y%tamdar_sfc(n)%p = y%tamdar_sfc(n)%p/typical_p_rms
      y%tamdar_sfc(n)%q = y%tamdar_sfc(n)%q/typical_q_rms

      if (iv%info(tamdar_sfc)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs + (y%tamdar_sfc(n)%u)**2 + (y%tamdar_sfc(n)%v)**2 + (y%tamdar_sfc(n)%t)**2 &
            + (y%tamdar_sfc(n)%p)**2 + (y%tamdar_sfc(n)%q)**2
      end if

      pertile_lhs = pertile_lhs + (y%tamdar_sfc(n)%u)**2 + (y%tamdar_sfc(n)%v)**2 + (y%tamdar_sfc(n)%t)**2 &
         + (y%tamdar_sfc(n)%p)**2 + (y%tamdar_sfc(n)%q)**2

      y%tamdar_sfc(n)%u = y%tamdar_sfc(n)%u/typical_u_rms
      y%tamdar_sfc(n)%v = y%tamdar_sfc(n)%v/typical_v_rms
      y%tamdar_sfc(n)%t = y%tamdar_sfc(n)%t/typical_t_rms
      y%tamdar_sfc(n)%p = y%tamdar_sfc(n)%p/typical_p_rms
      y%tamdar_sfc(n)%q = y%tamdar_sfc(n)%q/typical_q_rms
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_tamdar_sfc")

end subroutine da_check_xtoy_adjoint_tamdar_sfc


subroutine da_check_xtoy_adjoint_synop(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer  :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_synop")

   do n=iv%info(synop)%n1, iv%info(synop)%n2
      if (iv%info(synop)%zk(1,n) < 1.0 .and. sfc_assi_options /= 2) then
         y%synop(n)%u = 0.0
         y%synop(n)%v = 0.0
         y%synop(n)%t = 0.0
         y%synop(n)%p = 0.0
         y%synop(n)%q = 0.0
         cycle
      end if

      if ( sfc_assi_options == 2 ) then
          if (iv%synop(n)%u%qc < 0) y%synop(n)%u = 0.0
          if (iv%synop(n)%v%qc < 0) y%synop(n)%v = 0.0
          if (iv%synop(n)%t%qc < 0) y%synop(n)%t = 0.0
          if (iv%synop(n)%p%qc < 0) y%synop(n)%p = 0.0
          if (iv%synop(n)%q%qc < 0) y%synop(n)%q = 0.0
      end if

      y%synop(n)%u = y%synop(n)%u/typical_u_rms
      y%synop(n)%v = y%synop(n)%v/typical_v_rms
      y%synop(n)%t = y%synop(n)%t/typical_t_rms
      y%synop(n)%p = y%synop(n)%p/typical_p_rms
      y%synop(n)%q = y%synop(n)%q/typical_q_rms

      if (iv%info(synop)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs + (y%synop(n)%u)**2 + (y%synop(n)%v)**2 + (y%synop(n)%t)**2 &
            + (y%synop(n)%p)**2 + (y%synop(n)%q)**2
      end if

      pertile_lhs = pertile_lhs + (y%synop(n)%u)**2 + (y%synop(n)%v)**2 + (y%synop(n)%t)**2 &
         + (y%synop(n)%p)**2 + (y%synop(n)%q)**2

      y%synop(n)%u = y%synop(n)%u/typical_u_rms
      y%synop(n)%v = y%synop(n)%v/typical_v_rms
      y%synop(n)%t = y%synop(n)%t/typical_t_rms
      y%synop(n)%p = y%synop(n)%p/typical_p_rms
      y%synop(n)%q = y%synop(n)%q/typical_q_rms
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_synop")

end subroutine da_check_xtoy_adjoint_synop


subroutine da_check_xtoy_adjoint_rad(iv, y, adjtest_lhs, pertile_lhs)

   !------------------------------------------------------------------------------
   ! Purpose: Calculate innovation vector for radiance data.
   !------------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type(y_type),  intent(inout)  :: y             ! y = h (xa)
   real,          intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: inst, n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_rad")

   do inst = 1, iv%num_inst                 ! loop for sensor
      if (iv%instid(inst)%num_rad < 1) cycle

      do n= 1, iv%instid(inst)%num_rad       ! loop for pixel
         ! if (iv%instid(inst)%rad(n)%loc%proc_domain_with_halo) then
         if (iv%instid(inst)%info%proc_domain(1,n)) then
            do k = 1, iv%instid(inst)%nchan
               ! if ( iv%instid(inst)%tb_qc(k,n) >= obs_qc_pointer ) &
               adjtest_lhs = adjtest_lhs + &
                 ( y%instid(inst)%tb(k,n)/iv%instid(inst)%tb_error(k,n) )**2
            end do
         end if

         do k=1, iv%instid(inst)%nchan
            ! if ( iv%instid(inst)%tb_qc(k,n) >= obs_qc_pointer ) &
            pertile_lhs = pertile_lhs + &
             ( y%instid(inst)%tb(k,n)/iv%instid(inst)%tb_error(k,n) )**2

            y%instid(inst)%tb(k,n) =  &
              y%instid(inst)%tb(k,n) / (iv%instid(inst)%tb_error(k,n))**2
         end do
         ! end if
      end do
   end do 

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_rad")

end subroutine da_check_xtoy_adjoint_rad


subroutine da_transform_xtovp(grid, xb, xbx, xa, vp, be)

   !---------------------------------------------------------------------------
   ! Purpose: Transforms analysis to control variables (Vp-type)
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !
   ! Updates:
   !
   !       Implementation of multi-variate BE
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !---------------------------------------------------------------------------

   implicit none

   type(domain),            intent(inout) :: grid
   type(xb_type),           intent(in)    :: xb         ! First guess structure.
   type(xbx_type),          intent(in)    :: xbx        ! Header/non-gridded vars.
   type(x_type),            intent(inout) :: xa         ! Analysis increments.
   type(vp_type),           intent(out)   :: vp         ! CV on grid structure.
   type(be_type), optional, intent(in)    :: be         ! Background errors.

   real    :: vor(ims:ime,jms:jme,kms:kme) ! Vorticity.
   real    :: div(ims:ime,jms:jme,kms:kme) ! Divergence.

   real    :: one_over_m2(ims:ime,jms:jme) !   Multiplicative coeff.

   integer :: i, j, k , k1                 ! Loop counters.

   if (trace_use) call da_trace_entry("da_transform_xtovp")

   if ( (cv_options == 3) .or. (cv_options == 7) ) then
      write(unit=message(1),fmt='(A,I3)') 'Cannot perform transform_xtovp for cv_options == ',cv_options
      call da_error("da_transform_xtovp.inc",34,message(1:1))
   endif

   !----------------------------------------------------------------
   ! [1.0] Perform transform v = U^{-1} x
   !----------------------------------------------------------------      

   call da_zero_vp_type (vp)

   ! [2.2] Transform u, v to streamfunction via vorticity:


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



   call da_uv_to_vorticity(xb, xa % u, xa % v, vor)

   ! Convert vorticity to Del**2 psi:
   if (.not. global) then               
      if (fg_format == fg_format_wrf_arw_regional) then
         one_over_m2(its:ite,jts:jte) = 1.0 / (xb % map_factor(its:ite,jts:jte) * &
                                        xb % map_factor(its:ite,jts:jte))
         do k = kts, kte
            vor(its:ite,jts:jte,k) = &
              one_over_m2(its:ite,jts:jte)*vor(its:ite,jts:jte,k)
         end do
      else if (fg_format == fg_format_wrf_nmm_regional) then
         write(unit=message(1),fmt='(A,I5)') &
         "Needs to be developed for fg_format_nmm_regional = ",fg_format
         call da_error("da_transform_xtovp.inc",92,message(1:1))
      else
         write(unit=message(1),fmt='(A,I5,A,L10)') &
            ' Wrong choice of fg_format= ',fg_format,' with global = ',global
         call da_error("da_transform_xtovp.inc",96,message(1:1))
      end if
   end if

   ! Calculate psi:
   write (unit=stdout,fmt=*) ' calling Solve_PoissonEquation for Psi'
   call da_solve_poissoneqn_fct(grid,xbx, vor, vp%v1)

   ! [2.3] Transform u, v to velocity potential via divergence:

   call da_message((/'calling UV_To_Divergence'/))
   call da_uv_to_divergence(xb, xa % u, xa % v, div)

   ! Convert divergence to Del**2 chi:
   if (.not. global)  then              
      if (fg_format == fg_format_wrf_arw_regional) then
         do k = kts, kte
            div(its:ite,jts:jte,k) = &
               one_over_m2(its:ite,jts:jte) * div(its:ite,jts:jte,k)
         end do
      else if (fg_format == fg_format_wrf_nmm_regional) then
         write(unit=message(1),fmt='(A,I5)') &
         "Needs to be developed for fg_format_nmm_regional = ",fg_format
         call da_error("da_transform_xtovp.inc",119,message(1:1))
      else
         write(unit=message(1),fmt='(A,I5,A,L10)') &
            ' Wrong choice of fg_format= ',fg_format,' with global = ',global
         call da_error("da_transform_xtovp.inc",123,message(1:1))
      end if
   end if

   ! Calculate chi:

   call da_message((/' calling Solve_PoissonEquation for Chi'/))
   call da_solve_poissoneqn_fct(grid,xbx, div, vp%v2)

   ! [2.4] Transform chi to chi_u:
   call da_message((/' calculating chi_u'/))
   do k=kts,kte
      do j=jts,jte
         vp%v2(its:ite,j,k) = vp%v2(its:ite,j,k) - &
            be%reg_psi_chi(j,k)*vp%v1(its:ite,j,k)
      end do
   end do

   ! [2.5] Compute t_u:
   call da_message((/' calculating t_u'/))
   do k1=kts,kte
    do k=kts,kte
      do j=jts,jte
            vp%v3(its:ite,j,k) = vp%v3(its:ite,j,k) - be%reg_psi_t(j,k,k1)*vp%v1(its:ite,j,k1)
      end do
    end do
   end do
   if ( cv_options == 6 ) then
      do k1=kts,kte
         do k=kts,kte
            do j=jts,jte
               vp%v3(its:ite,j,k) = vp%v3(its:ite,j,k) + be%reg_chi_u_t(j,k,k1)*vp%v2(its:ite,j,k1)
            end do
         end do
      end do
   end if

   ! [2.6] Choice of moisture control variable:
   call da_message((/' calculating psudo rh'/))
   vp % v4(its:ite,jts:jte,kts:kte) = xa % q  (its:ite,jts:jte,kts:kte) /   &
      xb % qs (its:ite,jts:jte,kts:kte)

   if ( cv_options == 6 ) then
      do k1 = kts, kte
         do k = kts, kte
            do j = jts, jte
               vp%v4(its:ite,j,k) = vp%v4(its:ite,j,k)    - &
                                    be%reg_psi_rh(j,k,k1)*vp%v1(its:ite,j,k1)    + &
                                    be%reg_chi_u_rh(j,k,k1)*vp%v2(its:ite,j,k1)  + &
                                    be%reg_t_u_rh(j,k,k1)*vp%v3(its:ite,j,k1) + &
                                    be%reg_ps_u_rh(j,k1)*vp%v5(its:ite,j,1)
            end do
         end do
      end do
   end if

   ! [2.7] compute psfc_u:
   call da_message((/' calculating psfc_u '/))
   do j=jts,jte
      do i=its,ite
         vp % v5(i,j,1) = xa%psfc(i,j) - be%reg_psi_ps(j,k)*vp%v1(i,j,k)
      end do
   end do
   if ( cv_options == 6 ) then
      do j=jts,jte
         do i=its,ite
            vp % v5(i,j,1) = xa%psfc(i,j) + be%reg_chi_u_ps(j,k)*vp%v2(i,j,k)
         end do
      end do
   end if

   if (trace_use) call da_trace_exit("da_transform_xtovp")

end subroutine da_transform_xtovp


subroutine da_check(grid, config_flags, cv_size, xbx, be, ep, iv, vv, vp, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   integer, intent(in)              :: cv_size ! Size of cv array.
   type (xbx_type),   intent(inout) :: xbx   ! Header & non-gridded vars.
   type (be_type),    intent(in)    :: be    ! background error structure.
   type (ep_type),    intent(in)    :: ep    ! Ensemble perturbation structure.
   type (iv_type),    intent(inout) :: iv    ! ob. increment vector.
   type (vp_type),    intent(inout) :: vv    ! Grdipt/EOF CV.
   type (vp_type),    intent(inout) :: vp    ! Grdipt/level CV.
   type (y_type),     intent(inout) :: y             ! y = h (xa)

   integer :: sizec
   real    :: cvtest(1:cv_size)    ! background error structure.
   real    :: field(its:ite,jts:jte) ! Field for spectral transform test.

   call da_trace_entry("da_check")

   !----------------------------------------------------------------------------
   ! [1] Set up test data:
   !----------------------------------------------------------------------------

   ! Initialize cv values with random data:
   call random_number(cvtest(:))
   cvtest(:) = cvtest(:) - 0.5

   ! vv arrays initialized already.
   ! vp arrays initialized already.

   !----------------------------------------------------------------------------
   ! [2] Perform vtox adjoint tests:
   !----------------------------------------------------------------------------

   call da_message((/"Performing vtox adjoint tests"/))

   ! v_to_vv adjoint test:

   call da_check_cvtovv_adjoint(grid, cv_size, xbx, be, cvtest, vv)

   !-------------------------------------------------------------------------
   ! vv_to_vp adjoint test:
   !-------------------------------------------------------------------------

   call da_check_vvtovp_adjoint(grid, be % ne, grid%xb, be, vv, vp)

   !-------------------------------------------------------------------------
   ! vptox adjoint test:
   !-------------------------------------------------------------------------

   call da_check_vptox_adjoint(grid, be % ne, be, ep, vp, cv_size)

   !-------------------------------------------------------------------------
   ! vtox adjoint test: <x,x> = <v_adj,v>
   !-------------------------------------------------------------------------

   call da_check_vtox_adjoint(grid, cv_size, xbx, be, ep, cvtest, vv, vp)

   !----------------------------------------------------------------------------
   ! [2] Perform xtoy adjoint tests:
   !----------------------------------------------------------------------------

   call da_message((/"Performing xtoy adjoint tests"/))

   call da_allocate_y(iv, y)
   call da_zero_x(grid%xa)

   call da_setup_testfield(grid)

   ! WHY?
   ! Make cv_array random.

   ! call random_number(cvtest(1:cv_size))
   ! cvtest(1:cv_size) = cvtest(1:cv_size) - 0.5

   ! call da_transform_vtox(grid, cv_size, xbx, be, ep, cvtest, vv, vp)

   call da_check_xtoy_adjoint(cv_size, cvtest, xbx, be, grid, config_flags, iv, y)

   call da_deallocate_y(y)

   !----------------------------------------------------------------------------
   ! [4] Perform spectral test:
   !----------------------------------------------------------------------------

   if (global) then

      call da_message((/"Performing spectral tests"/))

      call random_number(field(:,:))
      field(:,:) = field(:,:) - 0.5

      sizec = (be % max_wave+1) * (be % max_wave+2)/2
      call da_test_spectral(be % max_wave, sizec, xbx, field)

   end if

   call da_trace_exit("da_check")

end subroutine da_check


real function da_dot(n,x,y)

   !-----------------------------------------------------------------------
   ! Purpose: forms the dot product of two vectors.
   ! uses unrolled loops for increments equal to one.
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in) :: n
   real,    intent(in) :: x(n)
   real,    intent(in) :: y(n)

   real    :: dtemp1
   integer :: i,m,mp1

   da_dot = 0.0
   if (n <= 0) return

   if (trace_use) call da_trace_entry("da_dot")    

   dtemp1 = 0.0

   ! code for both increments equal to 1

   if (n > 0) then
      m = mod(n,5)
      if (m /= 0) then
         do i = 1,m
            dtemp1 = dtemp1 + x(i)*y(i)
         end do
      end if
      if (n >= 5) then
         mp1 = m + 1
         do i = mp1,n,5
            dtemp1 = dtemp1 + x(i   )*y(i   ) + x(i + 1)*y(i + 1) + &
                              x(i + 2)*y(i + 2) + x(i + 3)*y(i + 3) + &
                              x(i + 4)*y(i + 4)
         end do
      end if
   end if

   da_dot = dtemp1

   if (trace_use) call da_trace_exit("da_dot")    


end function da_dot


real function da_dot_cv(cv_size, x, y, grid, mzs, jp_start, jp_end)

   !-----------------------------------------------------------------------
   ! Purpose: Forms the dot product of two vectors that are organized in the 
   ! format of a "cv_type".  
   !
   ! Capable of producing bitwise-exact results for distributed-memory 
   ! parallel runs for testing.  This feature is very slow and consumes 
   ! lots of memory. 
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in) :: cv_size           ! Size of array (tile).
   real,             intent(in) :: x(cv_size)        ! 1st vector.
   real,             intent(in) :: y(cv_size)        ! 1st vector.
   type(domain),     intent(in) :: grid              ! decomposed dimensions
   integer,          intent(in) :: mzs(7)            ! mz for each variable
                                                     ! (to identify 2D arrays)
   integer, optional,intent(in) :: jp_start, jp_end

   logical       :: lvarbc
   real, pointer :: xx(:), yy(:)                     ! Temporary vectors.
   real, pointer :: xg(:), yg(:)                     ! Serial data arrays.
   real          :: dtemp1(1), dtmp                  ! Temporary.

   if (trace_use) call da_trace_entry("da_dot_cv")

   allocate(xx(1:cv_size))
   allocate(yy(1:cv_size))
   xx = x
   yy = y

 ! VarBC parameters are global (no summation across processors)
 !-------------------------------------------------------------
   lvarbc = present(jp_start) .and. present(jp_end)
   if (lvarbc) lvarbc = lvarbc .and. (jp_end >= jp_start)   
   if (lvarbc .and. .not. rootproc) then
      xx(jp_start:jp_end) = 0.
      yy(jp_start:jp_end) = 0.
   end if
      
   ! Bitwise-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not
   ! testing.  This will always be .false. for a serial run or 
   ! one-processor 1 run.

   if (test_dm_exact) then

      ! Collect local cv arrays x and y to globally-sized serially-ordered 
      ! arrays xg and yg.  Note that xg and yg will only exist on the 
      ! monitor task.  

      if (rootproc) then
!         cv_size_domain = cv_size_domain_jb+cv_size_domain_je+cv_size_domain_jp+cv_size_domain_js   
         cv_size_domain = wrf_dm_sum_integer(cv_size)
         allocate(xg(1:cv_size_domain))
         allocate(yg(1:cv_size_domain))
      end if

      call da_cv_to_global(cv_size, cv_size_domain, xx, grid, mzs, xg)
      call da_cv_to_global(cv_size, cv_size_domain, yy, grid, mzs, yg)

      if (rootproc) then
         dtemp1(1) = da_dot(cv_size_domain, xg, yg)
         deallocate(xg, yg)
      end if

      ! Broadcast result from monitor to other tasks.  
      call wrf_dm_bcast_real(dtemp1, 1)

   else

      dtemp1(1) = da_dot(cv_size, xx, yy)
      
      !if (.not. global) then
         dtmp = dtemp1(1)
         ! summation across processors:
         dtemp1(1) = wrf_dm_sum_real(dtmp)
      !end if

   end if
  
   deallocate(xx,yy)

   da_dot_cv = dtemp1(1)

   if (trace_use) call da_trace_exit("da_dot_cv")

end function da_dot_cv
subroutine da_check_xtoy_adjoint_pseudo(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_pseudo")

   do n=iv%info(pseudo)%n1, iv%info(pseudo)%n2
      if (iv%info(pseudo)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs + &
            (y%pseudo(n)%u/typical_u_rms)**2 + &
            (y%pseudo(n)%v/typical_v_rms)**2 + &
            (y%pseudo(n)%t/typical_t_rms)**2 + &
            (y%pseudo(n)%p/typical_p_rms)**2 + &
            (y%pseudo(n)%q/typical_q_rms)**2
      end if

      pertile_lhs = pertile_lhs + &
         (y%pseudo(n)%u/typical_u_rms)**2 + &
         (y%pseudo(n)%v/typical_v_rms)**2 + &
         (y%pseudo(n)%t/typical_t_rms)**2 + &
         (y%pseudo(n)%p/typical_p_rms)**2 + &
         (y%pseudo(n)%q/typical_q_rms)**2

      y%pseudo(n)%u = y%pseudo(n)%u/typical_u_rms**2
      y%pseudo(n)%v = y%pseudo(n)%v/typical_v_rms**2
      y%pseudo(n)%t = y%pseudo(n)%t/typical_t_rms**2
      y%pseudo(n)%p = y%pseudo(n)%p/typical_p_rms**2
      y%pseudo(n)%q = y%pseudo(n)%q/typical_q_rms**2
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_pseudo")

end subroutine da_check_xtoy_adjoint_pseudo


subroutine da_check_xtoy_adjoint_qscat(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_qscat")

   do n=iv%info(qscat)%n1, iv%info(qscat)%n2
      if (iv%info(qscat)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs + (y%qscat(n)%u/typical_u_rms)**2 + (y%qscat(n)%v/typical_v_rms)**2
      end if

      pertile_lhs = pertile_lhs +(y%qscat(n)%u/typical_u_rms)**2 + (y%qscat(n)%v/typical_v_rms)**2

      y%qscat (n)%u = y%qscat (n)%u/typical_u_rms ** 2
      y%qscat (n)%v = y%qscat (n)%v/typical_v_rms ** 2
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_qscat")

end subroutine da_check_xtoy_adjoint_qscat


subroutine da_check_xtoy_adjoint_ssmt1(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer  :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_ssmt1")

   do n=iv%info(ssmt1)%n1, iv%info(ssmt1)%n2
      if (iv%info(ssmt1)%proc_domain(1,n)) then
         do k=1, iv%info(ssmt1)%levels(n)
            adjtest_lhs = adjtest_lhs + (y%ssmt1(n)%t(k)/typical_t_rms)**2
         end do
      end if

      do k=1, iv%info(ssmt1)%levels(n)
         pertile_lhs = pertile_lhs + (y%ssmt1(n)%t(k)/typical_t_rms)**2

         y%ssmt1(n)%t(k) = y%ssmt1(n)%t(k) / typical_t_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_ssmt1")

end subroutine da_check_xtoy_adjoint_ssmt1


subroutine da_check_xtoy_adjoint_ssmt2(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n, k          ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_ssmt2")

   do n=iv%info(ssmt2)%n1, iv%info(ssmt2)%n2
      if (iv%info(ssmt2)%proc_domain(1,n)) then
         do k=1, iv%info(ssmt2)%levels(n)
            adjtest_lhs = adjtest_lhs + (y%ssmt2(n)%rh(k)/typical_rh_rms)**2
         end do
      end if

      do k=1, iv%info(ssmt2)%levels(n)
         pertile_lhs = pertile_lhs + (y%ssmt2(n)%rh(k)/typical_rh_rms)**2

         y%ssmt2(n)%rh(k) = y%ssmt2(n)%rh(k) / typical_rh_rms ** 2
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_ssmt2")

end subroutine da_check_xtoy_adjoint_ssmt2


subroutine da_check_xtoy_adjoint_profiler(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs  

   integer :: n, k          ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_profiler")

   do n=iv%info(profiler)%n1, iv%info(profiler)%n2
      if (iv%info(profiler)%proc_domain(1,n)) then
         do k=1, iv%info(profiler)%levels(n)
            adjtest_lhs = adjtest_lhs &
              + (y%profiler(n)%u(k)/typical_u_rms)**2 + (y%profiler(n)%v(k)/typical_v_rms)**2

         end do
      end if

      do k=1, iv%info(profiler)%levels(n)
         pertile_lhs = pertile_lhs &
            + (y%profiler(n)%u(k)/typical_u_rms)**2 + (y%profiler(n)%v(k)/typical_v_rms)**2

         y%profiler(n)%u(k)= y%profiler(n)%u(k) / typical_u_rms ** 2
         y%profiler(n)%v(k)= y%profiler(n)%v(k) / typical_v_rms ** 2
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_profiler")

end subroutine da_check_xtoy_adjoint_profiler


subroutine da_check_xtoy_adjoint_buoy(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.
   
   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_buoy")

   do n=iv%info(buoy)%n1, iv%info(buoy)%n2
      if (iv%info(buoy)%zk(1,n) < 1.0 .and. sfc_assi_options /= 2) then
         y%buoy(n)%u = 0.0
         y%buoy(n)%v = 0.0
         y%buoy(n)%t = 0.0
         y%buoy(n)%p = 0.0
         y%buoy(n)%q = 0.0

         cycle
      end if

      if ( sfc_assi_options == 2 ) then
          if (iv%buoy(n)%u%qc < 0) y%buoy(n)%u = 0.0
          if (iv%buoy(n)%v%qc < 0) y%buoy(n)%v = 0.0
          if (iv%buoy(n)%t%qc < 0) y%buoy(n)%t = 0.0
          if (iv%buoy(n)%p%qc < 0) y%buoy(n)%p = 0.0
          if (iv%buoy(n)%q%qc < 0) y%buoy(n)%q = 0.0
      end if

      y%buoy(n)%u = y%buoy(n)%u/typical_u_rms
      y%buoy(n)%v = y%buoy(n)%v/typical_v_rms
      y%buoy(n)%t = y%buoy(n)%t/typical_t_rms
      y%buoy(n)%p = y%buoy(n)%p/typical_p_rms
      y%buoy(n)%q = y%buoy(n)%q/typical_q_rms

      if (iv%info(buoy)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs  &
                     + (y%buoy(n)%u)**2 &
                     + (y%buoy(n)%v)**2 &
                     + (y%buoy(n)%t)**2 &
                     + (y%buoy(n)%p)**2 &
                     + (y%buoy(n)%q)**2
      end if

      pertile_lhs = pertile_lhs &
                  + (y%buoy(n)%u)**2 &
                  + (y%buoy(n)%v)**2 &
                  + (y%buoy(n)%t)**2 &
                  + (y%buoy(n)%p)**2 &
                  + (y%buoy(n)%q)**2

      y%buoy(n)%u = y%buoy(n)%u/typical_u_rms
      y%buoy(n)%v = y%buoy(n)%v/typical_v_rms
      y%buoy(n)%t = y%buoy(n)%t/typical_t_rms
      y%buoy(n)%p = y%buoy(n)%p/typical_p_rms
      y%buoy(n)%q = y%buoy(n)%q/typical_q_rms
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_buoy")

end subroutine da_check_xtoy_adjoint_buoy


subroutine da_setup_testfield(grid)

   !----------------------------------------------------------------------------
   ! Purpose: produce test increment field based on grid%xb field.
   !
   ! Method:  pass through x=uv transfom to ensure satisfies boundary conditions
   !----------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)   :: grid

   integer :: i, j

   if (trace_use) call da_trace_entry("da_setup_testfield")

   !-------------------------------------------------------------------------
   ! [1.0]: initialise:
   !-------------------------------------------------------------------------

   write(unit=stdout, fmt='(/a/)') &
      'Starting da_setup_testfield ...'

   !-------------------------------------------------------------------------
   ! [2.0]: set up test increment field structure:
   !-------------------------------------------------------------------------

   ! [2.1] test wind, temperature, pressure, humidity and cloud parameters

   call da_set_tst_trnsf_fld(grid, grid%xa%u, grid%xb%u, typical_u_rms)
   call da_set_tst_trnsf_fld(grid, grid%xa%v, grid%xb%v, typical_v_rms)
   call da_set_tst_trnsf_fld(grid, grid%xa%w, grid%xb%w, typical_w_rms)
   call da_set_tst_trnsf_fld(grid, grid%xa%t, grid%xb%t, typical_t_rms)
   call da_set_tst_trnsf_fld(grid, grid%xa%p, grid%xb%p, typical_p_rms)
   call da_set_tst_trnsf_fld(grid, grid%xa%q, grid%xb%q, typical_q_rms)
   if ( ( use_rad .and. crtm_cloud ) .or. use_radar_rf .or. use_radar_rhv ) then
      call da_set_tst_trnsf_fld(grid, grid%xa%qcw, grid%xb%qcw, typical_qcw_rms)
      call da_set_tst_trnsf_fld(grid, grid%xa%qrn, grid%xb%qrn, typical_qrn_rms)
      call da_set_tst_trnsf_fld(grid, grid%xa%qci, grid%xb%qci, typical_qci_rms)
      call da_set_tst_trnsf_fld(grid, grid%xa%qsn, grid%xb%qsn, typical_qsn_rms)
      call da_set_tst_trnsf_fld(grid, grid%xa%qgr, grid%xb%qgr, typical_qgr_rms)
   end if

   ! [2.5] get test density increment from linearised ideal gas law:

   call da_pt_to_rho_lin(grid)

   grid%xa%psfc(grid%xp%its:grid%xp%ite, grid%xp%jts:grid%xp%jte) = &
   grid%xa%p   (grid%xp%its:grid%xp%ite, grid%xp%jts:grid%xp%jte, grid%xp%kts)

   if (print_detail_testing) then
      write(unit=stdout, fmt='(2a,4x,a,i8)') &
         'file:', "da_setup_testfield.inc", 'line:', 53

      write(unit=stdout, fmt=*) 'grid%xp%its, grid%xp%ite, grid%xp%jts, grid%xp%jte) =', &
         grid%xp%its, grid%xp%ite, grid%xp%jts, grid%xp%jte
   
      do j=grid%xp%jts, grid%xp%jte
         do i=grid%xp%its, grid%xp%ite
            if (i == j) then
               write(unit=stdout, fmt='(2(a,i4),a,f14.6)') &
                  'grid%xa%psfc(', i, ',', j, ') =', grid%xa%psfc(i, j)
            end if
         end do
      end do
   end if


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



   if (sfc_assi_options == 2) then
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
   end if


   if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. &
        use_ssmitbobs .or. use_ssmiretrievalobs) then

      ! Now do something for PW
      call da_transform_xtotpw(grid)

      ! GPS Refractivity:
      if (use_gpsrefobs) &
         call da_transform_xtogpsref_lin(grid)

      if (use_ssmt1obs .or. use_ssmt2obs .or. use_ssmitbobs .or. use_ssmiretrievalobs) then
         call da_transform_xtoseasfcwind_lin(grid)
      end if

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
   end if

   write(unit=stdout, fmt='(/a/)') 'End of da_setup_testfield.'

   if (trace_use) call da_trace_exit("da_setup_testfield")

end subroutine da_setup_testfield


subroutine da_check_sfc_assi(grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------
   
   implicit none

   type (domain), intent(inout)     :: grid
  
   type (iv_type),    intent(inout) :: iv    ! ob. increment vector.
   type (y_type),     intent(inout) :: y     ! y = h (grid%xa)

   real                           :: adj_ttl_lhs   ! < y, y >
   real                           :: adj_ttl_rhs   ! < x, x_adj >

   real                           :: partial_lhs   ! < y, y >
   real                           :: partial_rhs   ! < x, x_adj >

   real                           :: pertile_lhs   ! < y, y >
   real                           :: pertile_rhs   ! < x, x_adj >

   integer                        :: n

   real, dimension(ims:ime, jms:jme, kms:kme) :: xa2_u, xa2_v, xa2_t, &
                                                 xa2_p, xa2_q
 
   real, dimension(ims:ime, jms:jme)          :: xa2_u10, xa2_v10, xa2_t2, &
                                                 xa2_q2, xa2_tgrn, xa2_psfc


   if (trace_use) call da_trace_entry("da_check_sfc_assi")
  
   call da_message((/"check_sfc_assi: Adjoint Test Results:"/))
    
   xa2_u(ims:ime, jms:jme, kms:kme) = grid%xa%u(ims:ime, jms:jme, kms:kme)
   xa2_v(ims:ime, jms:jme, kms:kme) = grid%xa%v(ims:ime, jms:jme, kms:kme)
   xa2_t(ims:ime, jms:jme, kms:kme) = grid%xa%t(ims:ime, jms:jme, kms:kme)
   xa2_p(ims:ime, jms:jme, kms:kme) = grid%xa%p(ims:ime, jms:jme, kms:kme)
   xa2_q(ims:ime, jms:jme, kms:kme) = grid%xa%q(ims:ime, jms:jme, kms:kme)

   xa2_psfc(ims:ime, jms:jme) = grid%xa%psfc(ims:ime, jms:jme)
   xa2_tgrn(ims:ime, jms:jme) = grid%xa%tgrn(ims:ime, jms:jme)
   xa2_u10 (ims:ime, jms:jme) = grid%xa%u10 (ims:ime, jms:jme)
   xa2_v10 (ims:ime, jms:jme) = grid%xa%v10 (ims:ime, jms:jme)
   xa2_t2  (ims:ime, jms:jme) = grid%xa%t2  (ims:ime, jms:jme)
   xa2_q2  (ims:ime, jms:jme) = grid%xa%q2  (ims:ime, jms:jme)

   ! WHY?
   ! call check_psfc(grid, iv, y)

   call da_transform_xtowtq (grid)

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


   partial_lhs = 0.0
   pertile_lhs = 0.0

   do n=1, iv%info(synop)%nlocal
      call da_transform_xtopsfc(grid, iv, synop, iv%synop(:), y%synop(:))


      pertile_lhs = pertile_lhs &
                  + y%synop(n)%u * y%synop(n)%u &
                  + y%synop(n)%v * y%synop(n)%v &
                  + y%synop(n)%t * y%synop(n)%t &
                  + y%synop(n)%p * y%synop(n)%p &
                  + y%synop(n)%q * y%synop(n)%q

      if (iv%info(synop)%proc_domain(1,n)) then
         partial_lhs = partial_lhs &
                     + y%synop(n)%u * y%synop(n)%u &
                     + y%synop(n)%v * y%synop(n)%v &
                     + y%synop(n)%t * y%synop(n)%t &
                     + y%synop(n)%p * y%synop(n)%p &
                     + y%synop(n)%q * y%synop(n)%q
      end if
   end do

   !----------------------------------------------------------------------
   ! [5.0] Perform adjoint operation:
   !----------------------------------------------------------------------

   call da_zero_x(grid%xa)

   do n=1, iv%info(synop)%nlocal
      call da_transform_xtopsfc_adj(grid,iv, synop,iv%synop(:),y%synop(:),grid%xa)
   end do

   call da_transform_xtowtq_adj (grid)
   
   pertile_rhs = sum(grid%xa%u(ims:ime, jms:jme, kms:kme) * &
      xa2_u(ims:ime, jms:jme, kms:kme)) + &
                 sum(grid%xa%v(ims:ime, jms:jme, kms:kme) * &
      xa2_v(ims:ime, jms:jme, kms:kme)) + &
                 sum(grid%xa%t(ims:ime, jms:jme, kms:kme) * &
      xa2_t(ims:ime, jms:jme, kms:kme)) + &
                 sum(grid%xa%p(ims:ime, jms:jme, kms:kme) * &
      xa2_p(ims:ime, jms:jme, kms:kme)) + &
                 sum(grid%xa%q(ims:ime, jms:jme, kms:kme) * &
      xa2_q(ims:ime, jms:jme, kms:kme)) + &
                 sum(grid%xa%psfc(ims:ime, jms:jme) * xa2_psfc(ims:ime, jms:jme))

   !-------------------------------------------------------------------------
   ! [6.0] Calculate RHS of adjivnt test equation:
   !-------------------------------------------------------------------------
   
   partial_rhs = &
      sum(grid%xa%u(its:ite, jts:jte, kts:kte) * xa2_u(its:ite,jts:jte,kts:kte)) + &
      sum(grid%xa%v(its:ite, jts:jte, kts:kte) * xa2_v(its:ite,jts:jte,kts:kte)) + &
      sum(grid%xa%t(its:ite, jts:jte, kts:kte) * xa2_t(its:ite,jts:jte,kts:kte)) + &
      sum(grid%xa%p(its:ite, jts:jte, kts:kte) * xa2_p(its:ite,jts:jte,kts:kte)) + &
      sum(grid%xa%q(its:ite, jts:jte, kts:kte) * xa2_q(its:ite,jts:jte,kts:kte)) + &
      sum(grid%xa%psfc(its:ite, jts:jte) * xa2_psfc(its:ite, jts:jte))
   
   !-------------------------------------------------------------------------
   ! [7.0] Print output:
   !-------------------------------------------------------------------------
   
   write(unit=stdout, fmt='(A,1pe22.14)') &
        ' Tile < y, y     > = ', pertile_lhs, &
        ' Tile < x, x_adj > = ', pertile_rhs

   adj_ttl_lhs = wrf_dm_sum_real(partial_lhs)
   adj_ttl_rhs = wrf_dm_sum_real(partial_rhs)
   write (unit=stdout,fmt='(A,2F10.2)') &
      'TEST_COVERAGE_check_sfc_assi_A:  adj_ttl_lhs,adj_ttl_rhs = ', &
      adj_ttl_lhs,adj_ttl_rhs
   if (rootproc) then
      write(unit=stdout, fmt='(A,1pe22.14)') &
         ' Whole Domain < y, y     > = ', adj_ttl_lhs
      write(unit=stdout, fmt='(A,1pe22.14)') &
         ' Whole Domain < x, x_adj > = ', adj_ttl_rhs
   end if

   ! recover grid%xa
   grid%xa%u(ims:ime, jms:jme, kms:kme) = xa2_u(ims:ime, jms:jme, kms:kme)
   grid%xa%v(ims:ime, jms:jme, kms:kme) = xa2_v(ims:ime, jms:jme, kms:kme)
   grid%xa%t(ims:ime, jms:jme, kms:kme) = xa2_t(ims:ime, jms:jme, kms:kme)
   grid%xa%p(ims:ime, jms:jme, kms:kme) = xa2_p(ims:ime, jms:jme, kms:kme)
   grid%xa%q(ims:ime, jms:jme, kms:kme) = xa2_q(ims:ime, jms:jme, kms:kme)

   grid%xa%psfc(ims:ime, jms:jme) = xa2_psfc(ims:ime, jms:jme)
   grid%xa%tgrn(ims:ime, jms:jme) = xa2_tgrn(ims:ime, jms:jme)
   grid%xa%u10 (ims:ime, jms:jme) = xa2_u10 (ims:ime, jms:jme)
   grid%xa%v10 (ims:ime, jms:jme) = xa2_v10 (ims:ime, jms:jme)
   grid%xa%t2  (ims:ime, jms:jme) = xa2_t2  (ims:ime, jms:jme)
   grid%xa%q2  (ims:ime, jms:jme) = xa2_q2  (ims:ime, jms:jme)

   call wrf_shutdown

   if (trace_use) call da_trace_exit("da_check_sfc_assi")
   
end subroutine da_check_sfc_assi



subroutine da_check_psfc(grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------
     
   implicit none
   
   type (domain),  intent(inout)     :: grid
   type (iv_type), intent(inout)     :: iv  ! ob. increment vector.
   type (y_type),  intent(inout)     :: y   ! residual

   real                           :: adj_ttl_lhs   ! < y, y >
   real                           :: adj_ttl_rhs   ! < x, x_adj >

   real                           :: partial_lhs   ! < y, y >
   real                           :: partial_rhs   ! < x, x_adj >

   real                           :: pertile_lhs   ! < y, y >
   real                           :: pertile_rhs   ! < x, x_adj >

   integer                        :: n

   real, dimension(ims:ime, jms:jme) :: xa2_u10, xa2_v10, xa2_t2, &
                                        xa2_q2, xa2_psfc

   if (trace_use) call da_trace_entry("da_check_psfc")
   
   write(unit=stdout, fmt='(/3a,i6/a)') &
        'File: ', "da_check_psfc.inc", ', line:', 30, &
        'Adjoint Test Results:'

   ! save input

   xa2_psfc(ims:ime, jms:jme) = grid%xa%p   (ims:ime, jms:jme, kts)
   xa2_u10 (ims:ime, jms:jme) = grid%xa%u10 (ims:ime, jms:jme)
   xa2_v10 (ims:ime, jms:jme) = grid%xa%v10 (ims:ime, jms:jme)
   xa2_t2  (ims:ime, jms:jme) = grid%xa%t2  (ims:ime, jms:jme)
   xa2_q2  (ims:ime, jms:jme) = grid%xa%q2  (ims:ime, jms:jme)

   !----------------------------------------------------------------------

   partial_lhs = 0.0
   pertile_lhs = 0.0

   do n=1, iv%info(synop)%nlocal
      call da_transform_xtopsfc(grid, iv, synop, iv%synop(:), y%synop(:))
      pertile_lhs = pertile_lhs &
                  + y%synop(n)%u * y%synop(n)%u &
                  + y%synop(n)%v * y%synop(n)%v &
                  + y%synop(n)%t * y%synop(n)%t &
                  + y%synop(n)%p * y%synop(n)%p &
                  + y%synop(n)%q * y%synop(n)%q

      if (iv%info(synop)%proc_domain(1,n)) then
         partial_lhs = partial_lhs & 
                     + y%synop(n)%u * y%synop(n)%u &
                     + y%synop(n)%v * y%synop(n)%v &
                     + y%synop(n)%t * y%synop(n)%t &
                     + y%synop(n)%p * y%synop(n)%p &
                     + y%synop(n)%q * y%synop(n)%q
      end if
   end do

   !-------------------------------------------------------------------------
   ! [5.0] Perform adjoint operation:
   !-------------------------------------------------------------------------

   grid%xa%psfc(ims:ime, jms:jme) = 0.0
   grid%xa%tgrn(ims:ime, jms:jme) = 0.0
   grid%xa%u10 (ims:ime, jms:jme) = 0.0
   grid%xa%v10 (ims:ime, jms:jme) = 0.0
   grid%xa%t2  (ims:ime, jms:jme) = 0.0
   grid%xa%q2  (ims:ime, jms:jme) = 0.0
   
   do n=1, iv%info(synop)%nlocal
      call da_transform_xtopsfc_adj(grid,iv,synop,iv%synop(:),y%synop(:),grid%xa)
   end do

   pertile_rhs = sum(grid%xa%u10 (ims:ime, jms:jme) * xa2_u10 (ims:ime, jms:jme)) &
               + sum(grid%xa%v10 (ims:ime, jms:jme) * xa2_v10 (ims:ime, jms:jme)) &
               + sum(grid%xa%t2  (ims:ime, jms:jme) * xa2_t2  (ims:ime, jms:jme)) &
               + sum(grid%xa%q2  (ims:ime, jms:jme) * xa2_q2  (ims:ime, jms:jme)) &
               + sum(grid%xa%psfc(ims:ime, jms:jme) * xa2_psfc(ims:ime, jms:jme))

   partial_rhs = sum(grid%xa%u10 (its:ite, jts:jte) * xa2_u10 (its:ite, jts:jte)) &
               + sum(grid%xa%v10 (its:ite, jts:jte) * xa2_v10 (its:ite, jts:jte)) &
               + sum(grid%xa%t2  (its:ite, jts:jte) * xa2_t2  (its:ite, jts:jte)) &
               + sum(grid%xa%q2  (its:ite, jts:jte) * xa2_q2  (its:ite, jts:jte)) &
               + sum(grid%xa%psfc(its:ite, jts:jte) * xa2_psfc(its:ite, jts:jte))
   
   !----------------------------------------------------------------------
   ! [6.0] Calculate RHS of adjoint test equation:
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   ! [7.0] Print output:
   !----------------------------------------------------------------------
   
   write(unit=stdout, fmt='(A,1pe22.14)') &
      ' Tile < y, y     > = ', pertile_lhs, &
      ' Tile < x, x_adj > = ', pertile_rhs

   adj_ttl_lhs = wrf_dm_sum_real(partial_lhs)
   adj_ttl_rhs = wrf_dm_sum_real(partial_rhs)
   write (unit=stdout,fmt='(A,2F10.2)') &
      'TEST_COVERAGE_check_sfc_assi_B:  adj_ttl_lhs,adj_ttl_rhs = ', &
      adj_ttl_lhs,adj_ttl_rhs
   if (rootproc) then
      write(unit=stdout, fmt='(A,1pe22.14)') ' Whole Domain < y, y     > = ', &
         adj_ttl_lhs
      write(unit=stdout, fmt='(A,1pe22.14)') ' Whole Domain < x, x_adj > = ', &
         adj_ttl_rhs
   end if

   ! recover
   grid%xa%psfc(ims:ime, jms:jme) = xa2_psfc(ims:ime, jms:jme)
   grid%xa%u10 (ims:ime, jms:jme) = xa2_u10 (ims:ime, jms:jme)
   grid%xa%v10 (ims:ime, jms:jme) = xa2_v10 (ims:ime, jms:jme)
   grid%xa%t2  (ims:ime, jms:jme) = xa2_t2  (ims:ime, jms:jme)
   grid%xa%q2  (ims:ime, jms:jme) = xa2_q2  (ims:ime, jms:jme)

   if (trace_use) call da_trace_exit("da_check_psfc")

end subroutine da_check_psfc


subroutine da_set_tst_trnsf_fld(grid, va, vb, typical_rms)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain), intent(in)  :: grid   ! model state
   real,         intent(out) :: va(ims:ime, jms:jme, kms:kme)
   real,         intent(in)  :: vb(ims:ime, jms:jme, kms:kme)
   real,         intent(in)  :: typical_rms

   real    :: gbuf(ids:ide, jds:jde, kds:kde) 
   integer :: size3d
   integer :: k             ! loop counter.
   real    :: scale         ! scale factor.
   real    :: xy_dim_inv    ! 1 / real(iy * jx)
   real    :: field_mean    ! mean field.
   real    :: field_rms     ! rms field.

   if (trace_use) call da_trace_entry("da_set_tst_trnsf_fld")

   size3d = (ide-ids+1)*(jde-jds+1)*(kde-kds+1)

   call da_patch_to_global(grid, vb, gbuf)
   call wrf_dm_bcast_real(gbuf, size3d)

   xy_dim_inv  = 1.0 / real((ide-ids+1)*(jde-jds+1))

   do k = kts,kte
      field_mean = sum(gbuf(ids:ide, jds:jde, k)) * xy_dim_inv

      gbuf(ids:ide, jds:jde, k) = gbuf(ids:ide, jds:jde, k) - field_mean

      field_rms = sqrt(1.0e-6+sum(gbuf(ids:ide, jds:jde, k)**2) * xy_dim_inv)

      scale = typical_rms / field_rms

      va(its:ite, jts:jte, k) = scale * gbuf(its:ite, jts:jte, k)
   end do

   if (trace_use) call da_trace_exit("da_set_tst_trnsf_fld")

end subroutine da_set_tst_trnsf_fld


subroutine da_check_vtoy_adjoint(cv_size,grid, config_flags, vp, vv, xbx, be, ep, iv, y)

   !---------------------------------------------------------------------------
   !  Purpose: Perform V To Y Adjoint transform test                             
   !
   !  Method:  Perform adjoint test on complete transform: <y,y> = <v_adj,v>.!
   !---------------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size
   type(grid_config_rec_type), intent(inout) :: config_flags
   type(domain),               intent(inout) :: grid 
   type(vp_type),              intent(inout) :: vv    ! Grdipt/EOF CV.
   type(vp_type),              intent(inout) :: vp    ! Grdipt/level CV.
   type(xbx_type),             intent(inout) :: xbx   ! Header & non-gridded vars.
   type(be_type),              intent(in)    :: be    ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! ensemble perturbation structure.
   type(iv_type),              intent(inout) :: iv    ! ob. increment vector.
   type(y_type),               intent(inout) :: y     ! y = h (xa)

   real    :: cv(1:cv_size)          ! Test control variable.
   real    :: cv_2(1:cv_size)
   real    :: adj_sum_lhs               ! < y, y >
   real    :: adj_rhs,adj_sum_rhs       ! < v, v_adj >
   real    :: partial_lhs   ! < y, y >
   real    :: pertile_lhs   ! < y, y >  

   if (trace_use_dull) call da_trace_entry("da_check_vtoy_adjoint")

   write(unit=stdout, fmt='(/a/a)') &
        '       da_check_vtoy_adjoint:',&
        '---------------------------------------'

   call random_number(cv(:))
   cv(:) = cv(:) - 0.5

   cv_2(1:cv_size) = cv(1:cv_size)

   call da_zero_x(grid%xa)
   call da_zero_vp_type(vp)
   call da_zero_vp_type(vv)

   call da_transform_vtoy(cv_size, be, ep, cv, iv, vp, vv, vp, vv, xbx, y, grid, config_flags)

   !-------------------------------------------------------------------------
   ! [3.0] Calculate LHS of adjoint test equation and
   !        Rescale input to adjoint routine :
   !-------------------------------------------------------------------------

   call da_get_y_lhs_value(iv, y, partial_lhs, pertile_lhs, adj_sum_lhs)

   cv = 0.0

   ! WHY?
   ! call da_zero_vp_type(vp)
   ! call da_zero_vp_type(vv)
   ! call da_zero_x(grid%xa)      

   call da_transform_vtoy_adj(cv_size, be, ep, cv, iv, vp, vv, vp, vv, xbx, y, grid, &
      config_flags, .true.)

   adj_rhs = sum(cv(1:cv_size) * cv_2(1:cv_size))

   !-------------------------------------------------------------------------
   !      Print output:
   !-------------------------------------------------------------------------

   if (global) then
      adj_sum_rhs = adj_rhs
   else
      call mpi_allreduce(adj_rhs, adj_sum_rhs, 1, true_mpi_real, mpi_sum, &
                       comm, ierr)
   end if

   if (rootproc) then
      write(unit=stdout, fmt='(A,1pe22.14)') &
      'Whole Domain  < y, y     > = ', adj_sum_lhs
      write(unit=stdout, fmt='(A,1pe22.14)') &
         'Whole Domain  < v, v_adj > = ', adj_sum_rhs
   end if

   if (trace_use_dull) call da_trace_exit("da_check_vtoy_adjoint")

end subroutine da_check_vtoy_adjoint


subroutine da_get_y_lhs_value (iv, y, partial_lhs, pertile_lhs, adj_ttl_lhs) 

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------
    
   implicit none
   
   type(iv_type), intent(in)    :: iv    ! ob. increment vector.
   type(y_type),  intent(inout) :: y     ! y = h(xa)
   real,          intent(out)   :: partial_lhs, pertile_lhs, adj_ttl_lhs

   if (trace_use) call da_trace_entry("da_get_y_lhs_value")

   partial_lhs = 0.0
   pertile_lhs = 0.0

   if (iv%info(sound)%nlocal > 0)  call da_check_xtoy_adjoint_sound( iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(sonde_sfc)%nlocal > 0) call da_check_xtoy_adjoint_sonde_sfc( iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(mtgirs)%nlocal         > 0) call da_check_xtoy_adjoint_mtgirs   (iv, y, partial_lhs, pertile_lhs)
   if (iv%info(tamdar)%nlocal         > 0) call da_check_xtoy_adjoint_tamdar   (iv, y, partial_lhs, pertile_lhs)
   if (iv%info(tamdar_sfc)%nlocal     > 0) call da_check_xtoy_adjoint_tamdar_sfc (iv, y, partial_lhs, pertile_lhs)
   if (iv%info(synop)%nlocal          > 0) call da_check_xtoy_adjoint_synop    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(geoamv)%nlocal         > 0) call da_check_xtoy_adjoint_geoamv   (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(polaramv)%nlocal       > 0) call da_check_xtoy_adjoint_polaramv (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(airep)%nlocal          > 0) call da_check_xtoy_adjoint_airep    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(pilot)%nlocal          > 0) call da_check_xtoy_adjoint_pilot    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(radar)%nlocal          > 0) call da_check_xtoy_adjoint_radar    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(satem)%nlocal          > 0) call da_check_xtoy_adjoint_satem    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(metar)%nlocal          > 0) call da_check_xtoy_adjoint_metar    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(ships)%nlocal          > 0) call da_check_xtoy_adjoint_ships    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(gpspw)%nlocal          > 0) call da_check_xtoy_adjoint_gpspw    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(gpsref)%nlocal         > 0) call da_check_xtoy_adjoint_gpsref   (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(ssmi_tb)%nlocal        > 0) call da_check_xtoy_adjoint_ssmi_tb  (iv, y, partial_lhs, pertile_lhs)
   if (iv%info(ssmi_rv)%nlocal        > 0) call da_check_xtoy_adjoint_ssmi_rv  (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(ssmt2)%nlocal          > 0) call da_check_xtoy_adjoint_ssmt1    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(ssmt2)%nlocal          > 0) call da_check_xtoy_adjoint_ssmt2    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(qscat)%nlocal          > 0) call da_check_xtoy_adjoint_qscat    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(profiler)%nlocal       > 0) call da_check_xtoy_adjoint_profiler (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(buoy)%nlocal           > 0) call da_check_xtoy_adjoint_buoy     (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(bogus)%nlocal          > 0) call da_check_xtoy_adjoint_bogus    (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(radiance)%nlocal       > 0) call da_check_xtoy_adjoint_rad      (iv, y, partial_lhs, pertile_lhs) 
   if (iv%info(rain)%nlocal           > 0) call da_check_xtoy_adjoint_rain     (iv, y, partial_lhs, pertile_lhs)
   ! FIX? consider using dm_sum_real
   call mpi_allreduce( partial_lhs, adj_ttl_lhs, 1, true_mpi_real, mpi_sum, comm, ierr) 

   if (trace_use) call da_trace_exit("da_get_y_lhs_value")
   
end subroutine da_get_y_lhs_value


subroutine da_check_gradient(grid, config_flags, cv_size, xhat, cv, pdx, itertest, xbx, be, iv, y, re, j_cost)

   !-----------------------------------------------------------------------
   ! Purpose: Gradient test
   ! Adopted from grtest.f90 of GSI by Xin Zhang , September, 2011
   !-----------------------------------------------------------------------

   implicit none

   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   real,    intent(in)    ::   pdx
   integer, intent(in   ) :: itertest
   integer, intent(in)           :: cv_size ! Size of cv array.
   real,    intent(inout)        :: xhat(1:cv_size)  ! control variable (local).
   real,    intent(inout)           :: cv(1:cv_size)    ! control variable (local).
   type (xbx_type),   intent(inout) :: xbx   ! Header & non-gridded vars.
   type (be_type),    intent(in) :: be    ! background error structure.
   type (iv_type),    intent(inout) :: iv    ! ob. increment vector.
   type (y_type),     intent(inout) :: y             ! y = h (xa)
   type (y_type), intent(inout)  :: re    ! residual (o-a) structure.
   type (j_type), intent(inout)  :: j_cost                 ! cost function

   real    :: xdir(1:cv_size) , grad(1:cv_size), yhat(1:cv_size)
   real    :: zfy,zf0,zdf0,za,zfa,zdfa
   real    :: zabuf(itertest),zfabuf(itertest),ztf2buf(itertest)
   real    :: ZT1,ZB,ZFB,ztco,ZTC1,ZT2,ZTC1A,ZTC2,ZTF2
   real    :: ZAL,ZFAL,ZBL,ZFBL,ZTF2L
   real    :: ZTC00,ZTC02,ZTC10,ZTC12
   real    :: ZERMIN,ZT1TST,ZREF
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                           :: mz(7)
   integer                           :: sz(5)
   integer :: ibest,idig
   integer :: ii, jj

   call da_trace_entry("da_check_gradient")

   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be % ne /)
   sz = (/ be%cv%size1, be%cv%size2, be%cv%size3, be%cv%size4, be%cv%size5 /)
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

   call da_message((/' gradient test starting'/))

   if (pdx<=EPSILON(pdx)) then
      if (rootproc) write(6,*)'grtest, pdx=',pdx
      write(6,*)'grtest: pdx too small',pdx
      call da_trace_exit("da_check_gradient")
      return
   endif

   !----------------------------------------------------------------------------
   ! [1] Initial point
   !----------------------------------------------------------------------------

   call da_initialize_cv (cv_size, cv)
   call da_initialize_cv (cv_size, xhat)

   ! Initialize cv values with random data:
   call random_number(xhat(:))
   xhat(:) = xhat(:) - 0.5
   if (rootproc) write(6,*)'grtest: use random_number(xhat)'

   yhat=xhat

   call da_calculate_j(1, 1, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, yhat, cv, re, y, j_cost, grad, grid, config_flags)

   zfy = j_cost%total

   !----------------------------------------------------------------------------
   ! [1.1] Define perturbation direction ZH
   !----------------------------------------------------------------------------

   call da_message((/' The test direction is the opposite of the gradient '/))

   xdir  = -1.0  * grad
   !----------------------------------------------------------------------------
   ! [1.2] Set function f value and derivative at origin
   !----------------------------------------------------------------------------

   zf0=zfy
   zdf0=da_dot_cv(cv_size,grad,xdir,grid,mz,jp_start,jp_end)
   if (rootproc) write(6,*)'grtest: F(0)=',zf0,' DF(0)=',zdf0

   IF (ZDF0>0.0) write(6,*) 'GRTEST Warning, DF should be negative'
   IF (ABS(ZDF0) < SQRT(EPSILON(ZDF0))) THEN
      if (rootproc) write(6,*) 'GRTEST WARNING, DERIVATIVE IS TOO SMALL'
   ENDIF

   !----------------------------------------------------------------------------
   ! [2] Loop on test point
   !----------------------------------------------------------------------------

   ztf2buf(1)=0.0

   DO jj=1,itertest

      za=pdx*(10.0**(jj-1))

      if (rootproc) write(6,*)'grtest iter=',jj,' alpha=',za
   
      !----------------------------------------------------------------------------
      ! [2.1] Compute f and df at new point y=x+a.h
      !----------------------------------------------------------------------------

      do ii=1,cv_size
         yhat(ii) = xhat(ii) + za * xdir(ii)
      end do

      call da_calculate_j(1, 1, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                          be%cv%size_jl, xbx, be, iv, yhat, cv, re, y, j_cost, grad, grid, config_flags)
      zfy = j_cost%total

      zfa=zfy
      zdfa=da_dot_cv(cv_size,grad,xdir,grid,mz,jp_start,jp_end)

      if (rootproc) write(6,*)'grtest: alpha=',za,' F(a)=',zfa,' DF(a)=',zdfa
   
      zabuf(jj)=za
      zfabuf(jj)=zfa

      !----------------------------------------------------------------------------
      ! [2.2] Quantity TC0=f(a)/f(0)-1 
      !----------------------------------------------------------------------------

!         if f is continuous then TC0->1 at origin,
!         at least linearly with a.

      IF (ABS(zf0)<=TINY(zf0)) THEN
!           do not compute T1 in this unlikely case
         if (rootproc) write(6,*) 'grtest: Warning: zf0 is suspiciously small.'
         if (rootproc) write(6,*) 'grtest: F(a)-F(0)=',zfa-zf0
      ELSE
         ztco=zfa/zf0-1.0
         if (rootproc) write(6,*)'grtest: continuity TC0=',ztco
      ENDIF

      !----------------------------------------------------------------------------
      !                     f(a)-f(0)
      ! [2.3] Quantity T1=-----------
      !                      a.df(0)
      !----------------------------------------------------------------------------

!         if df is the gradient then T1->1 at origin,
!         linearly with a. T1 is undefined if df(0)=0.

      IF (ABS(za*zdf0)<=SQRT(TINY(zf0))) THEN
         if (rootproc) write(6,*)'grtest: Warning: could not compute ',&
          & 'gradient test T1, a.df(0)=',za*zdf0
      ELSE
         zt1=(zfa-zf0)/(za*zdf0)
         if (rootproc) write(6,*)'grtest: gradient T1=',zt1
      ENDIF

      !----------------------------------------------------------------------------
      ! [2.4] Quantity TC1=( f(a)-f(0)-a.df(0) )/a
      !----------------------------------------------------------------------------

!         if df is the gradient and df is continuous,
!         then TC1->0 linearly with a.
      ZTC1=(ZFA-ZF0-ZA*ZDF0)/ZA
      if (rootproc) write(6,*)'grtest: grad continuity TC1=',ZTC1

      !----------------------------------------------------------------------------
      ! [2.5] Quantity T2=( f(a)-f(0)-a.df(0) )*2/a**2
      !----------------------------------------------------------------------------

!         if d2f exists then T2 -> d2f(0) linearly with a.
      ZT2=(ZFA-ZF0-ZA*ZDF0)*2.0/(ZA**2)
      if (rootproc) write(6,*)'grtest: second derivative T2=',ZT2

      !----------------------------------------------------------------------------
      ! [2.6] Quantity TC1A=df(a)-df(0)
      !----------------------------------------------------------------------------

!         if df is the gradient in a and df is continuous,
!         then TC1A->0 linearly with a.
      ZTC1A=ZDFA-ZDF0
      if (rootproc) write(6,*)'grtest: a-grad continuity TC1A=',ZTC1A

      !----------------------------------------------------------------------------
      ! [2.7] Quantity TC2=( 2(f(0)-f(a))+ a(df(0)+df(a))/a**2
      !----------------------------------------------------------------------------

!         if f is exactly quadratic, then TC2=0, always: numerically
!         it has to -> 0 when a is BIG. Otherwise TC2->0 linearly for
!         small a is trivially implied by TC1A and T2.
      ZTC2=(2.0*(ZF0-ZFA)+ZA*(ZDF0+ZDFA))/(ZA**2)
      if (rootproc) write(6,*)'grtest: quadraticity TC2=',ZTC2

      !----------------------------------------------------------------------------
      !                     2   f(0)-f(b)   f(a)-f(b)
      ! [2.8] Quantity TF2=---( --------- + --------- )
      !                     a       b          a-b
      !----------------------------------------------------------------------------

!         if 0, a and b are distinct and f is quadratic then
!         TF2=d2f, always. The estimate is most stable when a,b are big.
!         This test works only after two loops, but it is immune against
!         gradient bugs. 

      IF (jj>=2) THEN
         ZB =ZABUF (jj-1)
         ZFB=ZFABUF(jj-1)
         ZTF2=2.0/ZA*((ZF0-ZFB)/ZB+(ZFA-ZFB)/(ZA-ZB))
         if (rootproc) write(6,*)'grtest: convexity ZTF2=',ZTF2
         ztf2buf(jj)=ztf2
      ENDIF

! End loop
   ENDDO

   !----------------------------------------------------------------------------
   ! [3] Comment on the results
   !----------------------------------------------------------------------------

!       TC0(0)/TC0(2)<.011 -> df looks continuous
!       item with (T1<1 and 1-T1 is min) = best grad test item
!       reldif(TF2(last),TF2(last-1)) = precision on quadraticity
   
   !----------------------------------------------------------------------------
   !       3.1 Fundamental checks
   !----------------------------------------------------------------------------

   if (rootproc) then
      write(6,*) 'GRTEST: TENTATIVE CONCLUSIONS :'

      ZTC00=ABS(zfabuf(1)-zf0)
      ZTC02=ABS(zfabuf(3)-zf0)
      IF( ZTC00/zabuf(1)  <=  1.5*(ZTC02/zabuf(3)) )THEN
         write(6,*) 'GRTEST: function f looks continous.'
      ELSE
         write(6,*) 'GRTEST: WARNING f does not look continuous',&
          & ' (perhaps truncation problem)'
      ENDIF
   
   !----------------------------------------------------------------------------
   !       3.2 Gradient quality
   !----------------------------------------------------------------------------

      IF (ABS(zdf0)<=SQRT(TINY(zf0))) THEN
         write(6,*) 'GRTEST: The gradient is 0, which is unusual !'
         ZTC10=ABS(zfabuf(1)-zf0)
         ZTC12=ABS(zfabuf(3)-zf0)
         IF( ZTC10/zabuf(1)**2  <=  1.1*ZTC12/zabuf(3)**2)THEN
            write(6,*)'GRTEST: The gradient looks good anyway.'
         ENDIF
      ELSE
!        Find best gradient test index
         ZERMIN=HUGE(0.0)
         ibest=-1
         DO jj=1,itertest
            ZT1TST=(zfabuf(jj)-zf0)/(zabuf(jj)*zdf0)
            ZT1TST=ABS(ZT1TST-1.0)
            IF (ZT1TST<ZERMIN) THEN
               ibest=jj
               ZERMIN=ZT1TST
            ENDIF
         ENDDO
         IF(ibest == -1)THEN
            write(6,*)'GRTEST: gradient test problem : bad ',&
             & 'gradient, non-convex cost, or truncation errors ?'
         ELSE
            idig=INT(-LOG(ZERMIN+TINY(ZERMIN))/LOG(10.0))
            write(6,*)'GRTEST: the best gradient test found has ',&
             & idig,' satisfactory digits.'
            IF(idig <= 1)THEN
               write(6,*)'GRTEST: SAYS: THE GRADIENT IS VERY BAD.'
            ELSEIF(idig <= 3)THEN
               write(6,*)'GRTEST: SAYS: THE GRADIENT IS SUSPICIOUS.'
            ELSEIF(idig <= 5)THEN
               write(6,*)'GRTEST: SAYS: THE GRADIENT IS ACCEPTABLE.'
            ELSE
               write(6,*)'GRTEST: SAYS: THE GRADIENT IS EXCELLENT.'
            ENDIF

            IF (ibest<=itertest-2) THEN
               ZTC10=ABS(zfabuf(ibest         )-zf0-zabuf(ibest         )*zdf0)/zabuf(ibest         )
               ZTC12=ABS(zfabuf(ibest+2)-zf0-zabuf(ibest+2)*zdf0)/zabuf(ibest+2)
               IF(ZTC10/zabuf(ibest) <=  1.1*ZTC12/zabuf(ibest+2) )THEN
                  write(6,*)'GRTEST: Gradient convergence looks good.'
               ELSE
                  write(6,*)'GRTEST: Gradient convergence is suspicious.'
               ENDIF
            ELSE
               write(6,*)'GRTEST: could not check grad convergence.'
            ENDIF
         ENDIF
      ENDIF

   !----------------------------------------------------------------------------
   !         3.3 Quadraticity
   !      finite difference quadraticity test (gradient-free)
   !----------------------------------------------------------------------------

      ZTF2=ztf2buf(itertest)
      ZTF2L=ztf2buf(itertest-1)
      write(6,*) 'GRTEST: finite diff. d2f estimate no1:',ZTF2
      write(6,*) 'GRTEST: finite diff. d2f estimate no2:',ZTF2L
      ZREF=(ABS(ZTF2L)+ABS(ZTF2))/2.0
      IF (ZREF<=TINY(ZREF)) THEN
         write(6,*) 'GRTEST: they are too small to decide whether ',&
          & 'they agree or not.'
      ELSE
         idig=INT(-LOG(ABS(ZTF2L-ZTF2)/ZREF+TINY(ZTF2))/LOG(10.0))
         write(6,*) 'GRTEST: the fin.dif. estimates of d2f ',&
          & 'have ',idig,' satisfactory digits.'
         IF(idig <= 1)THEN
            write(6,*) 'GRTEST: THE FD-QUADRATICITY IS BAD.'
         ELSEIF(idig <= 3)THEN
            write(6,*) 'GRTEST:: THE FD-QUADRATICITY IS SUSPICIOUS.'
         ELSEIF(idig <= 5)THEN
            write(6,*) 'GRTEST: THE FD-QUADRATICITY IS ACCEPTABLE.'
         ELSE
            write(6,*) 'GRTEST: THE FD-QUADRATICITY IS EXCELLENT.'
         ENDIF
      ENDIF

      write(6,*) 'grtest: Goodbye.'
   endif

   call da_trace_exit("da_check_gradient")

end subroutine da_check_gradient



end module da_test
