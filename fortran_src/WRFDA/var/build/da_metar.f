












module da_metar

   use module_domain, only : domain

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, &
      missing, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv, fails_error_max, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, metar,max_error_thickness, &
      sfc_assi_options, trace_use_dull, anal_type_verify, kts,kte, &
      kms,kme,sfc_assi_options_1,sfc_assi_options_2, max_ext_its, &
      qcstat_conv_unit,ob_vars, &
      convert_fd2uv, convert_uv2fd, max_error_spd, max_error_dir, &
      max_omb_spd, max_omb_dir, pi, qc_rej_both, & 
      wind_sd_metar, wind_stats_sd 
   use da_grid_definitions, only : da_ffdduv, da_ffdduv_model, da_ffdduv_diagnose
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      maxmin_type
   use da_interpolation, only : da_to_zk,da_interp_lin_2d, &
      da_interp_lin_3d,da_interp_lin_3d_adj,da_interp_lin_2d_adj
      
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_physics, only : da_sfc_pre,da_transform_xtopsfc_adj, &
      da_transform_xtopsfc, da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual,da_obs_sfc_correction, da_convert_zk
   use da_tracing, only : da_trace_entry, da_trace_exit

  

  type residual_metar1_type
    real          :: u                        
    real          :: v                        
    real          :: t                        
    real          :: p                        
    real          :: q                        
  end type residual_metar1_type

  type maxmin_metar_stats_type
     type (maxmin_type)         :: u, v, t, p, q
  end type maxmin_metar_stats_type

  type stats_metar_type
     type (maxmin_metar_stats_type)  :: maximum, minimum
     type (residual_metar1_type)     :: average, rms_err
  end type stats_metar_type

contains

subroutine da_ao_stats_metar (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O
   type(y_type),   intent (in)   :: ob            ! Observation structure.

   type (stats_metar_type) :: stats
   integer                 :: nu, nv, nt, np, nq
   integer                 :: n
   real                    :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_metar")


   nu = 0
   nv = 0
   nt = 0
   np = 0
   nq = 0

   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%maximum%t = maxmin_type (missing_r, 0, 0)
   stats%maximum%p = maxmin_type (missing_r, 0, 0)
   stats%maximum%q = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%p = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_metar1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(metar)%nlocal
       if (iv%info(metar)%proc_domain(1,n)) then

          u_inc = re%metar(n)%u
          v_inc = re%metar(n)%v
          u_obs = ob%metar(n)%u
          v_obs = ob%metar(n)%v

          if (.not. wind_sd_metar .and. wind_stats_sd) &
             call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                     iv%metar(n)%u%qc, iv%metar(n)%v%qc, convert_uv2fd)
          if (wind_sd_metar .and. .not. wind_stats_sd) &
             call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                     iv%metar(n)%u%qc, iv%metar(n)%v%qc, convert_fd2uv)

         call da_stats_calculate (n, 0, iv%metar(n)%u%qc, & 
            u_inc, nu, &
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%metar(n)%v%qc, & 
            v_inc, nv, &
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate (n, 0, iv%metar(n)%t%qc, & 
            re%metar(n)%t, nt, &
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate (n, 0, iv%metar(n)%p%qc, & 
            re%metar(n)%p, np, &
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate (n, 0, iv%metar(n)%q%qc, & 
            re%metar(n)%q, nq, &
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)

      end if    ! end if (iv%info(metar)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (np)
   call da_proc_sum_int (nq)
   iv%nstats(metar) = nu + nv + nt + np + nq

   call da_proc_stats_combine(stats%average%u, stats%rms_err%u, &
      stats%minimum%u%value, stats%maximum%u%value, &
      stats%minimum%u%n, stats%maximum%u%n, &
      stats%minimum%u%l, stats%maximum%u%l)
   call da_proc_stats_combine(stats%average%v, stats%rms_err%v, &
      stats%minimum%v%value, stats%maximum%v%value, &
      stats%minimum%v%n, stats%maximum%v%n, &
      stats%minimum%v%l, stats%maximum%v%l)
   call da_proc_stats_combine(stats%average%t, stats%rms_err%t, &
      stats%minimum%t%value, stats%maximum%t%value, &
      stats%minimum%t%n, stats%maximum%t%n, &
      stats%minimum%t%l, stats%maximum%t%l)
   call da_proc_stats_combine(stats%average%p, stats%rms_err%p, &
      stats%minimum%p%value, stats%maximum%p%value, &
      stats%minimum%p%n, stats%maximum%p%n, &
      stats%minimum%p%l, stats%maximum%p%l)
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
      stats%minimum%q%value, stats%maximum%q%value, &
      stats%minimum%q%n, stats%maximum%q%n, &
      stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. np /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for metar'
         call da_print_stats_metar(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_metar")

end subroutine da_ao_stats_metar


subroutine da_jo_and_grady_metar(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type(jo_type), intent(inout) :: jo          ! Obs cost function.

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_metar")


   jo % metar_u = 0.0
   jo % metar_v = 0.0
   jo % metar_t = 0.0
   jo % metar_p = 0.0
   jo % metar_q = 0.0

   if (iv%info(metar)%nlocal > 0) then
      do n=1, iv%info(metar)%nlocal
         jo_grad_y%metar(n)%u = -re%metar(n)%u / &
                              (iv%metar(n)%u%error * iv%metar(n)%u%error)
         jo_grad_y%metar(n)%v = -re%metar(n)%v / &
                              (iv%metar(n)%v%error * iv%metar(n)%v%error)
         jo_grad_y%metar(n)%t = -re%metar(n)%t / &
                              (iv%metar(n)%t%error * iv%metar(n)%t%error)
         jo_grad_y%metar(n)%p = -re%metar(n)%p / &
                              (iv%metar(n)%p%error * iv%metar(n)%p%error)
         jo_grad_y%metar(n)%q = -re%metar(n)%q / &
                              (iv%metar(n)%q%error * iv%metar(n)%q%error)

         if (iv%info(metar)%proc_domain(1,n)) then
            jo % metar_u = jo % metar_u - re%metar(n)%u * jo_grad_y%metar(n)%u
            jo % metar_v = jo % metar_v - re%metar(n)%v * jo_grad_y%metar(n)%v
            jo % metar_t = jo % metar_t - re%metar(n)%t * jo_grad_y%metar(n)%t
            jo % metar_p = jo % metar_p - re%metar(n)%p * jo_grad_y%metar(n)%p
            jo % metar_q = jo % metar_q - re%metar(n)%q * jo_grad_y%metar(n)%q
         end if
      end do

      jo % metar_u = 0.5 * jo % metar_u
      jo % metar_v = 0.5 * jo % metar_v
      jo % metar_t = 0.5 * jo % metar_t
      jo % metar_p = 0.5 * jo % metar_p
      jo % metar_q = 0.5 * jo % metar_q
   end if

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_metar")

end subroutine da_jo_and_grady_metar


subroutine da_residual_metar(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !----------------------------------------------------------------------
   ! Purpose: Calculate residuals for metars
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual vector (O-A).

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type)              :: n_obs_bad
   integer                           :: n

   if (trace_use_dull) call da_trace_entry("da_residual_metar")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % p % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(metar)%nlocal
      np_available = np_available + 5
      re%metar(n)%u = da_residual(n, 0, y%metar(n)%u, iv%metar(n)%u, n_obs_bad % u)
      re%metar(n)%v = da_residual(n, 0, y%metar(n)%v, iv%metar(n)%v, n_obs_bad % v)
      re%metar(n)%t = da_residual(n, 0, y%metar(n)%t, iv%metar(n)%t, n_obs_bad % t)
      re%metar(n)%p = da_residual(n, 0, y%metar(n)%p, iv%metar(n)%p, n_obs_bad % p)
      re%metar(n)%q = da_residual(n, 0, y%metar(n)%q, iv%metar(n)%q, n_obs_bad % q)
   end do

   np_missing = np_missing + n_obs_bad % u % num % miss + &
      n_obs_bad % v % num % miss + n_obs_bad % t % num % miss + &
      n_obs_bad % p % num % miss + n_obs_bad % q % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad + &
      n_obs_bad % v % num % bad + n_obs_bad % t % num % bad + &
      n_obs_bad % p % num % bad + n_obs_bad % q % num % bad
   np_obs_used = np_obs_used + n_obs_bad % u % num % use + &
      n_obs_bad % v % num % use + n_obs_bad % t % num % use + &
      n_obs_bad % p % num % use + n_obs_bad % q % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_metar")

end subroutine da_residual_metar


subroutine da_oi_stats_metar (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_metar_type)  :: stats
   integer                  :: nu, nv, nt, np, nq
   integer                  :: n
   real                     :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_metar")

   nu = 0
   nv = 0
   nt = 0
   np = 0
   nq = 0

   stats%maximum%u = maxmin_type(missing_r, 0, 0)
   stats%maximum%v = maxmin_type(missing_r, 0, 0)
   stats%maximum%t = maxmin_type(missing_r, 0, 0)
   stats%maximum%p = maxmin_type(missing_r, 0, 0)
   stats%maximum%q = maxmin_type(missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%p = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_metar1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(metar)%nlocal
      if (iv%info(metar)%proc_domain(1,n)) then

         u_inv = iv%metar(n)%u%inv
         v_inv = iv%metar(n)%v%inv
         u_obs = ob%metar(n)%u
         v_obs = ob%metar(n)%v

        if (.not. wind_sd_metar .and. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%metar(n)%u%qc, iv%metar(n)%v%qc, convert_uv2fd)
        if (wind_sd_metar .and. .not. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%metar(n)%u%qc, iv%metar(n)%v%qc, convert_fd2uv)

         call da_stats_calculate(iv%info(metar)%obs_global_index(n), &
            0, iv%metar(n)%u%qc, &
            u_inv, nu, &
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate(iv%info(metar)%obs_global_index(n), &
            0, iv%metar(n)%v%qc, &
            v_inv, nv, &
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate(iv%info(metar)%obs_global_index(n), &
             0, iv%metar(n)%t%qc, &
             iv%metar(n)%t%inv, nt, &
             stats%minimum%t, stats%maximum%t, &
             stats%average%t, stats%rms_err%t)
         call da_stats_calculate(iv%info(metar)%obs_global_index(n), &
            0, iv%metar(n)%p%qc, &
            iv%metar(n)%p%inv, np, &
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate(iv%info(metar)%obs_global_index(n), &
            0, iv%metar(n)%q%qc, &
            iv%metar(n)%q%inv, nq, &
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)
      end if
   end do

   ! Do inter-processor communication to gather statistics.

   call da_proc_sum_int(nu)
   call da_proc_sum_int(nv)
   call da_proc_sum_int(nt)
   call da_proc_sum_int(np)
   call da_proc_sum_int(nq)

   call da_proc_stats_combine(stats%average%u, stats%rms_err%u, &
      stats%minimum%u%value, stats%maximum%u%value, &
      stats%minimum%u%n, stats%maximum%u%n, &
      stats%minimum%u%l, stats%maximum%u%l)
   call da_proc_stats_combine(stats%average%v, stats%rms_err%v, &
      stats%minimum%v%value, stats%maximum%v%value, &
      stats%minimum%v%n, stats%maximum%v%n, &
      stats%minimum%v%l, stats%maximum%v%l)
   call da_proc_stats_combine(stats%average%t, stats%rms_err%t, &
      stats%minimum%t%value, stats%maximum%t%value, &
      stats%minimum%t%n, stats%maximum%t%n, &
      stats%minimum%t%l, stats%maximum%t%l)
   call da_proc_stats_combine(stats%average%p, stats%rms_err%p, &
      stats%minimum%p%value, stats%maximum%p%value, &
      stats%minimum%p%n, stats%maximum%p%n, &
      stats%minimum%p%l, stats%maximum%p%l)
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
       stats%minimum%q%value, stats%maximum%q%value, &
       stats%minimum%q%n, stats%maximum%q%n, &
       stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. np /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for metar'
         call da_print_stats_metar(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_metar")

end subroutine da_oi_stats_metar


subroutine da_print_stats_metar(stats_unit, nu, nv, nt, np, nq, metar)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,           intent(in)    :: stats_unit
   integer,           intent(inout) :: nu, nv, nt, np, nq
   type (stats_metar_type), intent(in)    :: metar

   if (trace_use_dull) call da_trace_entry("da_print_stats_metar")

   write(unit=stats_unit, fmt='(6a/)') &
      '   var             ', &
      'u (m/s)     n    k    ', &
      'v (m/s)     n    k    ', &
      't (K)       n    k    ', &
      'p (Pa)      n    k    ', &
      'q (kg/kg)   n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
      '  Number: ', nu, nv, nt, np, nq

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1
   if (nt < 1) nt = 1
   if (np < 1) np = 1
   if (nq < 1) nq = 1

   write(unit=stats_unit, fmt='((a,4(f12.4,2i5),e12.4,2i5))') &
      ' Minimum(n,k): ', metar%minimum%u, metar%minimum%v, metar%minimum%t, &
                       metar%minimum%p, metar%minimum%q, &
      ' Maximum(n,k): ', metar%maximum%u, metar%maximum%v, metar%maximum%t, &
                       metar%maximum%p, metar%maximum%q
   write(unit=stats_unit, fmt='((a,4(f12.4,10x),e12.4,10x))') &
      ' Average     : ', metar%average%u/real(nu), metar%average%v/real(nv), &
                       metar%average%t/real(nt), metar%average%p/real(np), &
                       metar%average%q/real(nq), &
      '    RMSE     : ', sqrt(metar%rms_err%u/real(nu)), &
                       sqrt(metar%rms_err%v/real(nv)), &
                       sqrt(metar%rms_err%t/real(nt)), &
                       sqrt(metar%rms_err%p/real(np)), &
                       sqrt(metar%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_metar")

end subroutine da_print_stats_metar


subroutine da_transform_xtoy_metar (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer                      :: n   ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_p(:)

   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_metar")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_v(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_t(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_q(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_p(iv%info(metar)%n1:iv%info(metar)%n2))

      allocate (ub(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (vb(1,iv%info(metar)%n1:iv%info(metar)%n2))


      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xa%u, iv%info(metar),model_u)
      call da_interp_lin_3d (grid%xa%v, iv%info(metar),model_v)
      call da_interp_lin_3d (grid%xa%t, iv%info(metar),model_t)
      call da_interp_lin_3d (grid%xa%q, iv%info(metar),model_q)

      call da_interp_lin_2d (grid%xa%psfc, iv%info(metar), 1, model_p)

      call da_interp_lin_3d (grid%xb%u, iv%info(metar), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(metar), vb)

      do n=iv%info(metar)%n1,iv%info(metar)%n2

         if (wind_sd_metar) then
            call da_uv_to_sd_lin(y%metar(n)%u,y%metar(n)%v,model_u(1,n),model_v(1,n),ub(1,n),vb(1,n))
         else
            y%metar(n)%u = model_u(1,n)
            y%metar(n)%v = model_v(1,n)
         end if

         y%metar(n)%t = model_t(1,n)
         y%metar(n)%q = model_q(1,n)
         y%metar(n)%p = model_p(n)

      end do

      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_p)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      ! [2.0] Surface assmiilation approach 2
      call da_transform_xtopsfc(grid,iv,metar,iv%metar(:),y%metar(:))
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_metar")

end subroutine da_transform_xtoy_metar


subroutine da_transform_xtoy_metar_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid          ! first guess state.
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_psfc(:)

   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_metar_adj")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_v(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_t(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_q(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_psfc(iv%info(metar)%n1:iv%info(metar)%n2))

      allocate (ub(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (vb(1,iv%info(metar)%n1:iv%info(metar)%n2))

      call da_interp_lin_3d (grid%xb%u, iv%info(metar), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(metar), vb)


      ! [1.2] Interpolate horizontally:
      do n=iv%info(metar)%n1, iv%info(metar)%n2

         if (wind_sd_metar) then
             call da_uv_to_sd_adj(jo_grad_y%metar(n)%u, &
                                  jo_grad_y%metar(n)%v, model_u(1,n), model_v(1,n), ub(1,n), vb(1,n))
         else
            model_u(1,n)  = jo_grad_y%metar(n)%u
            model_v(1,n)  = jo_grad_y%metar(n)%v
         end if

         model_t(1,n)  = jo_grad_y%metar(n)%t
         model_q(1,n)  = jo_grad_y%metar(n)%q
         model_psfc(n) = jo_grad_y%metar(n)%p

      end do
      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(metar), model_u)
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(metar), model_v)
      call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(metar), model_t)
      call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(metar), model_q)
      call da_interp_lin_2d_adj (jo_grad_x%psfc, iv%info(metar), 1, model_psfc)

      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc_adj(grid,iv,metar,iv%metar(:), jo_grad_y%metar(:),jo_grad_x)
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_metar_adj")

end subroutine da_transform_xtoy_metar_adj


subroutine da_check_max_iv_metar(iv,ob, it, num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it    ! outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)
   type(y_type),  intent(in)    :: ob      ! Observation structure.


   logical :: failed,failed1,failed2
   integer :: n

  if (trace_use_dull) call da_trace_entry("da_check_max_iv_metar")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(metar)%n1,iv%info(metar)%n2
      if(.not. qc_rej_both)then
         if(wind_sd_metar)then
            failed=.false.
            if( iv%metar(n)%u%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%u, max_error_spd,failed)
                if( iv%info(metar)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,metar,1,1) = num_qcstat_conv(1,metar,1,1) + 1
                    if(failed) then
                       num_qcstat_conv(2,metar,1,1) = num_qcstat_conv(2,metar,1,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                       'metar',ob_vars(1),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
                    end if
                end if
            end if

            failed=.false.
            if( iv%metar(n)%v%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%v, max_error_dir,failed)
                if( iv%info(metar)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,metar,2,1) = num_qcstat_conv(1,metar,2,1) + 1
                    if(failed)then
                       num_qcstat_conv(2,metar,2,1) = num_qcstat_conv(2,metar,2,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                            'metar',ob_vars(2),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
                    end if
                end if
            end if
         else
            failed=.false.
            if( iv%metar(n)%u%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%u, max_error_uv,failed)
                if( iv%info(metar)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,metar,1,1) = num_qcstat_conv(1,metar,1,1) + 1
                    if(failed) then
                       num_qcstat_conv(2,metar,1,1) = num_qcstat_conv(2,metar,1,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                             'metar',ob_vars(1),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
                    end if
                end if
             end if
             failed=.false.
             if( iv%metar(n)%v%qc >= obs_qc_pointer ) then
                 call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%v, max_error_uv,failed)
                 if( iv%info(metar)%proc_domain(1,n) ) then
                     num_qcstat_conv(1,metar,2,1) = num_qcstat_conv(1,metar,2,1) + 1
                     if(failed)then
                        num_qcstat_conv(2,metar,2,1) = num_qcstat_conv(2,metar,2,1) + 1
                        write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                              'metar',ob_vars(2),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
                     end if
                 end if
             end if
          end if

          if(wind_sd_metar)then
             if(iv%metar(n)%u%qc == fails_error_max .or.  abs(iv%metar(n)%u%inv) >= max_omb_spd) then
                iv%metar(n)%u%qc = fails_error_max
                iv%metar(n)%u%inv = 0.0
             endif
             if(iv%metar(n)%v%qc == fails_error_max .or.  abs(iv%metar(n)%v%inv) >= max_omb_dir) then
                iv%metar(n)%v%qc = fails_error_max
                iv%metar(n)%v%inv = 0.0
             endif
          endif
       else
          failed1=.false.
          failed2=.false.
          if( iv%metar(n)%v%qc >= obs_qc_pointer .or. iv%metar(n)%u%qc >= obs_qc_pointer )  then
              if(wind_sd_metar)then
                 call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%u, max_error_spd,failed1)
                 call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%v, max_error_dir,failed2)
              else
                 call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%u, max_error_uv,failed1)
                 call da_max_error_qc (it,iv%info(metar), n, iv%metar(n)%v, max_error_uv,failed2)
              endif
          endif

          if( iv%info(metar)%proc_domain(1,n) ) then
              num_qcstat_conv(1,metar,1,1) = num_qcstat_conv(1,metar,1,1) + 1
              num_qcstat_conv(1,metar,2,1) = num_qcstat_conv(1,metar,2,1) + 1

              if(failed1 .or. failed2) then
                 num_qcstat_conv(2,metar,1,1) = num_qcstat_conv(2,metar,1,1) + 1
                 write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                       'metar',ob_vars(1),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
                 num_qcstat_conv(2,metar,2,1) = num_qcstat_conv(2,metar,2,1) + 1
                 write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                       'metar',ob_vars(2),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
              endif
           endif

           if(wind_sd_metar)then
              if(iv%metar(n)%u%qc == fails_error_max .or. iv%metar(n)%v%qc == fails_error_max .or. &
                 abs(iv%metar(n)%v%inv) >= max_omb_dir .or. abs(iv%metar(n)%u%inv) >= max_omb_spd )then
                 iv%metar(n)%u%qc = fails_error_max
                 iv%metar(n)%v%qc = fails_error_max
                 iv%metar(n)%u%inv = 0.0
                 iv%metar(n)%v%inv = 0.0
              endif
           else
              if(iv%metar(n)%u%qc == fails_error_max .or. iv%metar(n)%v%qc == fails_error_max ) then
                 iv%metar(n)%u%qc = fails_error_max
                 iv%metar(n)%v%qc = fails_error_max
                 iv%metar(n)%u%inv = 0.0
                 iv%metar(n)%v%inv = 0.0
              endif
           endif
        endif

      failed=.false.
      if( iv%metar(n)%t%qc >= obs_qc_pointer )  then
      call da_max_error_qc (it, iv%info(metar), n, iv%metar(n)%t, max_error_t , failed)
      if( iv%info(metar)%proc_domain(1,n) ) then
      num_qcstat_conv(1,metar,3,1)= num_qcstat_conv(1,metar,3,1) + 1
      if(failed) then
      num_qcstat_conv(2,metar,3,1)= num_qcstat_conv(2,metar,3,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
      'metar',ob_vars(3),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%metar(n)%p%qc >= obs_qc_pointer )  then 
      call da_max_error_qc (it, iv%info(metar), n, iv%metar(n)%p, max_error_p , failed)         
      if( iv%info(metar)%proc_domain(1,n) ) then
      num_qcstat_conv(1,metar,5,1)= num_qcstat_conv(1,metar,5,1) + 1
      if(failed) then
      num_qcstat_conv(2,metar,5,1)= num_qcstat_conv(2,metar,5,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
      'metar',ob_vars(5),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%metar(n)%q%qc >= obs_qc_pointer ) then
       if( iv%metar(n)%t%qc == fails_error_max .or. iv%metar(n)%p%qc == fails_error_max) then
       failed=.true.
       iv%metar(n)%q%qc  = fails_error_max
       iv%metar(n)%q%inv = 0.0
       else
       call da_max_error_qc (it, iv%info(metar), n, iv%metar(n)%q, max_error_q , failed)
       endif
      if( iv%info(metar)%proc_domain(1,n) ) then
      num_qcstat_conv(1,metar,4,1)= num_qcstat_conv(1,metar,4,1) + 1
      if(failed) then
      num_qcstat_conv(2,metar,4,1)= num_qcstat_conv(2,metar,4,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
      'metar',ob_vars(4),iv%info(metar)%lat(1,n),iv%info(metar)%lon(1,n),0.01*ob%metar(n)%p
      end if
      end if
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_metar")

end subroutine da_check_max_iv_metar
subroutine da_get_innov_vector_metar(it,num_qcstat_conv, grid, ob, iv)

   !-------------------------------------------------------------------------
   ! Purpose: TBD   
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-------------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid    ! first guess state.
   type(y_type),     intent(inout) :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)
   real  :: speed, direction
   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real, allocatable :: model_u(:,:)  ! Model value u at oblocation.
   real, allocatable :: model_v(:,:)  ! Model value v at oblocation.
   real, allocatable :: model_t(:,:)  ! Model value t at oblocation.
   real, allocatable :: model_p(:,:)  ! Model value p at oblocation.
   real, allocatable :: model_q(:,:)  ! Model value q at oblocation.
   real, allocatable :: model_hsm(:,:)

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.

   real    :: hd, psfcm

   real    :: ho, to, qo

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_metar")

   allocate (model_u  (1,iv%info(metar)%n1:iv%info(metar)%n2))
   allocate (model_v  (1,iv%info(metar)%n1:iv%info(metar)%n2))
   allocate (model_t  (1,iv%info(metar)%n1:iv%info(metar)%n2))
   allocate (model_q  (1,iv%info(metar)%n1:iv%info(metar)%n2))
   allocate (model_p  (1,iv%info(metar)%n1:iv%info(metar)%n2))
   allocate (model_hsm(1,iv%info(metar)%n1:iv%info(metar)%n2))

   if ( it > 1 ) then
      do n=iv%info(metar)%n1,iv%info(metar)%n2
         if (iv%metar(n)%u%qc == fails_error_max) iv%metar(n)%u%qc = 0
         if (iv%metar(n)%v%qc == fails_error_max) iv%metar(n)%v%qc = 0
         if (iv%metar(n)%t%qc == fails_error_max) iv%metar(n)%t%qc = 0
         if (iv%metar(n)%p%qc == fails_error_max) iv%metar(n)%p%qc = 0
         if (iv%metar(n)%q%qc == fails_error_max) iv%metar(n)%q%qc = 0
      end do
   end if

   if (sfc_assi_options == sfc_assi_options_1) then

      do n=iv%info(metar)%n1,iv%info(metar)%n2
         ! [1.1] Get horizontal interpolation weights:

         i   = iv%info(metar)%i(1,n)
         j   = iv%info(metar)%j(1,n)
         dx  = iv%info(metar)%dx(1,n)
         dy  = iv%info(metar)%dy(1,n)
         dxm = iv%info(metar)%dxm(1,n)
         dym = iv%info(metar)%dym(1,n)

         ! Surface correction

         iv%metar(n)%p%inv = ob%metar(n)%p
         iv%metar(n)%t%inv = ob%metar(n)%t
         iv%metar(n)%q%inv = ob%metar(n)%q
         iv%metar(n)%u%inv = ob%metar(n)%u
         iv%metar(n)%v%inv = ob%metar(n)%v

         if (iv % metar(n) % h > missing_r) then
            do k=kts,kte
               v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                  + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
            end do

            hd = v_h(kts) - iv % metar(n) % h

            if (abs(hd) <= Max_StHeight_Diff .or. anal_type_verify) then
               if (iv % metar(n) % h < v_h(kts)) then
                  iv%info(metar)%zk(:,n) = 1.0+1.0e-6
                  call da_obs_sfc_correction(iv%info(metar), iv%metar(n), n, grid%xb)
               else
                  call da_to_zk(iv % metar(n) % h, v_h, v_interp_h, iv%info(metar)%zk(1,n))
               end if
            end if
         else if (ob % metar(n) % p > 1.0) then
            do k=kts,kte
               v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                  + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
            end do

            call da_to_zk(ob % metar(n) % p, v_p, v_interp_p, iv%info(metar)%zk(1,n))

            if (iv%info(metar)%zk(1,n) < 0.0 .and.  .not.anal_type_verify) then
               iv % metar(n) % p % inv = missing_r
               iv % metar(n) % p % qc  = missing_data
               iv%info(metar)%zk(:,n) = 1.0+1.0e-6
            end if
         end if
      end do

      call da_convert_zk(iv%info(metar))

      !------------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !------------------------------------------------------------------------

      if (.not.anal_type_verify ) then
         do n=iv%info(metar)%n1,iv%info(metar)%n2
            if (iv%info(metar)%zk(1,n) < 0.0) then
               iv % metar(n) % u % qc = missing_data
               iv % metar(n) % v % qc = missing_data
               iv % metar(n) % t % qc = missing_data
               iv % metar(n) % q % qc = missing_data
               iv % metar(n) % p % qc = missing_data
            end if
         end do
      end if

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xb%u, iv%info(metar),model_u)
      call da_interp_lin_3d (grid%xb%v, iv%info(metar),model_v)
      call da_interp_lin_3d (grid%xb%t, iv%info(metar),model_t)
      call da_interp_lin_3d (grid%xb%q, iv%info(metar),model_q)
      call da_interp_lin_3d (grid%xb%p, iv%info(metar),model_p)

   else if (sfc_assi_options == sfc_assi_options_2) then
      ! 1.2.1 Surface assimilation approach 2 
      !(10-m u, v, 2-m t, q, and sfc_p)

      call da_interp_lin_2d (grid%xb%u10,  iv%info(metar), 1, model_u)
      call da_interp_lin_2d (grid%xb%v10,  iv%info(metar), 1, model_v)
      call da_interp_lin_2d (grid%xb%t2,   iv%info(metar), 1, model_t)
      call da_interp_lin_2d (grid%xb%q2,   iv%info(metar), 1, model_q)
      call da_interp_lin_2d (grid%xb%psfc, iv%info(metar), 1, model_p)
      call da_interp_lin_2d (grid%xb%terr, iv%info(metar), 1, model_hsm)

      do n=iv%info(metar)%n1,iv%info(metar)%n2

         iv%metar(n)%p%inv = ob%metar(n)%p
         iv%metar(n)%t%inv = ob%metar(n)%t
         iv%metar(n)%q%inv = ob%metar(n)%q
         iv%metar(n)%u%inv = ob%metar(n)%u
         iv%metar(n)%v%inv = ob%metar(n)%v

         if (iv%metar(n)%p%qc >= 0) then
            ! model surface p, t, q, h at observed site:

            ho = iv%metar(n)%h
            to = -888888.0
            qo = -888888.0

            if (iv%metar(n)%t%qc >= 0 .and. iv%metar(n)%q%qc >= 0) then
               to = ob%metar(n)%t
               qo = ob%metar(n)%q
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to, qo)
            else if (iv%metar(n)%t%qc >= 0 .and. iv%metar(n)%q%qc < 0) then
               to = ob%metar(n)%t
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to)
            else
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho)
            end if

            ! Pressure at the observed height:
            model_p(1,n) = psfcm
         end if
      end do
   end if

   do n=iv%info(metar)%n1,iv%info(metar)%n2

      !-----------------------------------------------------------------------
      ! [3.0] Fast interpolation:
      !-----------------------------------------------------------------------
      if (wind_sd_metar) then
         call da_ffdduv_model (speed,direction,model_u(1,n), model_v(1,n), convert_uv2fd)

         if (ob%metar(n)%u > missing_r .AND. iv%metar(n)%u%qc >= obs_qc_pointer) then
             iv%metar(n)%u%inv = iv%metar(n)%u%inv - speed
         else
             iv % metar(n) % u % inv = 0.0
         end if

         if (ob%metar(n)%v > missing_r .AND. iv%metar(n)%v%qc >= obs_qc_pointer) then
             iv%metar(n)%v%inv = iv%metar(n)%v%inv - direction
             if (iv%metar(n)%v%inv > 180.0 ) iv%metar(n)%v%inv = iv%metar(n)%v%inv - 360.0
             if (iv%metar(n)%v%inv < -180.0 ) iv%metar(n)%v%inv = iv%metar(n)%v%inv + 360.0
         else
             iv % metar(n) % v % inv = 0.0
         end if
      else
         if (ob % metar(n) % u > missing_r .AND. iv % metar(n) % u % qc >= obs_qc_pointer) then
             iv % metar(n) % u % inv = iv%metar(n)%u%inv - model_u(1,n)
         else
             iv % metar(n) % u % inv = 0.0
         end if

         if (ob % metar(n) % v > missing_r .AND. iv % metar(n) % v % qc >= obs_qc_pointer) then
             iv % metar(n) % v % inv = iv%metar(n)%v%inv - model_v(1,n)
         else
             iv % metar(n) % v % inv = 0.0
         end if
      end if

      !if (ob % metar(n) % p > 0.0 .AND. iv % metar(n) % p % qc >= obs_qc_pointer) then
      if ( iv % metar(n) % p % qc >= obs_qc_pointer ) then
         iv % metar(n) % p % inv = iv%metar(n)%p%inv - model_p(1,n)
      else
         iv % metar(n) % p % inv = 0.0
      end if

      if (ob % metar(n) % t > 0.0 .AND. &
         iv % metar(n) % t % qc >= obs_qc_pointer) then
         iv % metar(n) % t % inv = iv%metar(n)%t%inv - model_t(1,n)
      else
         iv % metar(n) % t % inv = 0.0
      end if

      if (ob % metar(n) % q > 0.0 .AND. &
         iv % metar(n) % q % qc >= obs_qc_pointer) then
         iv % metar(n) % q % inv = iv%metar(n)%q%inv - model_q(1,n)
      else
         iv % metar(n) % q % inv = 0.0
      end if
   end do

   !-----------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_metar(iv,ob, it, num_qcstat_conv)

   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_q)
   deallocate (model_p)
   deallocate (model_hsm)

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_metar")

end subroutine da_get_innov_vector_metar


subroutine da_calculate_grady_metar(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(inout) :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_metar")

   do n=1, iv%info(metar)%nlocal
      if (iv%metar(n)%u%qc < obs_qc_pointer) re%metar(n)%u = 0.0
      if (iv%metar(n)%v%qc < obs_qc_pointer) re%metar(n)%v = 0.0
      if (iv%metar(n)%t%qc < obs_qc_pointer) re%metar(n)%t = 0.0
      if (iv%metar(n)%p%qc < obs_qc_pointer) re%metar(n)%p = 0.0
      if (iv%metar(n)%q%qc < obs_qc_pointer) re%metar(n)%q = 0.0

      jo_grad_y%metar(n)%u = -re%metar(n)%u / (iv%metar(n)%u%error * iv%metar(n)%u%error)
      jo_grad_y%metar(n)%v = -re%metar(n)%v / (iv%metar(n)%v%error * iv%metar(n)%v%error)
      jo_grad_y%metar(n)%t = -re%metar(n)%t / (iv%metar(n)%t%error * iv%metar(n)%t%error)
      jo_grad_y%metar(n)%p = -re%metar(n)%p / (iv%metar(n)%p%error * iv%metar(n)%p%error)
      jo_grad_y%metar(n)%q = -re%metar(n)%q / (iv%metar(n)%q%error * iv%metar(n)%q%error)
   end do
  
if (trace_use_dull) call da_trace_exit("da_calculate_grady_metar")
     
end subroutine da_calculate_grady_metar




end module da_metar

