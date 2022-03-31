












module da_buoy 

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, max_error_uv, max_error_t, rootproc, &
      buoy, max_error_p,max_error_q, trace_use_dull,fails_error_max, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv,sfc_assi_options, anal_type_verify, &
      kms,kme,kts,kte,sfc_assi_options_1,sfc_assi_options_2, max_ext_its,&
      qcstat_conv_unit,ob_vars, &
      convert_fd2uv, convert_uv2fd, max_error_spd, max_error_dir, &
      max_omb_spd, max_omb_dir, pi, qc_rej_both, &
      wind_sd_buoy, wind_stats_sd
   use da_grid_definitions, only : da_ffdduv, da_ffdduv_model, da_ffdduv_diagnose
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_to_zk, &
      da_interp_lin_3d,da_interp_lin_3d_adj, &
      da_interp_lin_2d, da_interp_lin_2d_adj, da_interp_lin_2d_partial
   use da_par_util1, only : da_proc_sum_int
   use da_par_util, only : da_proc_stats_combine
   use da_physics, only : da_sfc_pre, da_transform_xtopsfc, da_transform_xtopsfc_adj, &
                          da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_obs_sfc_correction, da_convert_zk
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_buoy1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: p                        
      real          :: q                        
   end type residual_buoy1_type

   type maxmin_buoy_stats_type
      type (maxmin_type)         :: u, v, t, p, q
   end type maxmin_buoy_stats_type

   type stats_buoy_type
      type (maxmin_buoy_stats_type)  :: maximum, minimum
      type (residual_buoy1_type)     :: average, rms_err
   end type stats_buoy_type

contains

subroutine da_ao_stats_buoy (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O
   type(y_type),   intent (in)   :: ob            ! Observation structure.

   type (stats_buoy_type) :: stats
   integer                :: nu, nv, nt, np, nq
   integer                :: n
   real                   :: u_inc, v_inc, u_obs, v_obs
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_buoy")    

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

   stats%average = residual_buoy1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(buoy)%nlocal
      if (iv%info(buoy)%proc_domain(1,n)) then

          u_inc = re%buoy(n)%u
          v_inc = re%buoy(n)%v
          u_obs = ob%buoy(n)%u
          v_obs = ob%buoy(n)%v

          if (.not. wind_sd_buoy .and. wind_stats_sd) &
             call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                     iv%buoy(n)%u%qc, iv%buoy(n)%v%qc, convert_uv2fd)
          if (wind_sd_buoy .and. .not. wind_stats_sd) &
             call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                     iv%buoy(n)%u%qc, iv%buoy(n)%v%qc, convert_fd2uv)

         call da_stats_calculate (n, 0, iv%buoy(n)%u%qc, & 
            u_inc, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%buoy(n)%v%qc, & 
            v_inc, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate (n, 0, iv%buoy(n)%t%qc, & 
            re%buoy(n)%t, nt, & 
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate (n, 0, iv%buoy(n)%p%qc, & 
            re%buoy(n)%p, np, & 
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate (n, 0, iv%buoy(n)%q%qc, & 
            re%buoy(n)%q, nq, & 
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)

      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (np)
   call da_proc_sum_int (nq)
   iv%nstats(buoy) = nu + nv + nt + np + nq

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for buoy'
         call da_print_stats_buoy(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_buoy")    

 end subroutine da_ao_stats_buoy


subroutine da_jo_and_grady_buoy(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(in)    :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_buoy")

   jo % buoy_u = 0.0
   jo % buoy_v = 0.0
   jo % buoy_t = 0.0
   jo % buoy_p = 0.0
   jo % buoy_q = 0.0

   do n=1, iv%info(buoy)%nlocal
      jo_grad_y%buoy(n)%u = -re%buoy(n)%u / (iv%buoy(n)%u%error * iv%buoy(n)%u%error)
      jo_grad_y%buoy(n)%v = -re%buoy(n)%v / (iv%buoy(n)%v%error * iv%buoy(n)%v%error)
      jo_grad_y%buoy(n)%t = -re%buoy(n)%t / (iv%buoy(n)%t%error * iv%buoy(n)%t%error)
      jo_grad_y%buoy(n)%p = -re%buoy(n)%p / (iv%buoy(n)%p%error * iv%buoy(n)%p%error)
      jo_grad_y%buoy(n)%q = -re%buoy(n)%q / (iv%buoy(n)%q%error * iv%buoy(n)%q%error)

      if (iv%info(buoy)%proc_domain(1,n)) then
         jo % buoy_u = jo % buoy_u - re%buoy(n)%u * jo_grad_y%buoy(n)%u
         jo % buoy_v = jo % buoy_v - re%buoy(n)%v * jo_grad_y%buoy(n)%v
         jo % buoy_t = jo % buoy_t - re%buoy(n)%t * jo_grad_y%buoy(n)%t
         jo % buoy_p = jo % buoy_p - re%buoy(n)%p * jo_grad_y%buoy(n)%p
         jo % buoy_q = jo % buoy_q - re%buoy(n)%q * jo_grad_y%buoy(n)%q
      end if
   end do

   jo % buoy_u = 0.5 * jo % buoy_u
   jo % buoy_v = 0.5 * jo % buoy_v
   jo % buoy_t = 0.5 * jo % buoy_t
   jo % buoy_p = 0.5 * jo % buoy_p
   jo % buoy_q = 0.5 * jo % buoy_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_buoy")
     
end subroutine da_jo_and_grady_buoy


subroutine da_residual_buoy(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
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

   if (trace_use_dull) call da_trace_entry("da_residual_buoy")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % p % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(buoy)%nlocal
      np_available = np_available + 5
      re%buoy(n)%u = da_residual(n, 0, y%buoy(n)%u, iv%buoy(n)%u, n_obs_bad % u)
      re%buoy(n)%v = da_residual(n, 0, y%buoy(n)%v, iv%buoy(n)%v, n_obs_bad % v)
      re%buoy(n)%t = da_residual(n, 0, y%buoy(n)%t, iv%buoy(n)%t, n_obs_bad % t)
      re%buoy(n)%p = da_residual(n, 0, y%buoy(n)%p, iv%buoy(n)%p, n_obs_bad % p)
      re%buoy(n)%q = da_residual(n, 0, y%buoy(n)%q, iv%buoy(n)%q, n_obs_bad % q)
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

   if (trace_use_dull) call da_trace_exit("da_residual_buoy")

end subroutine da_residual_buoy


subroutine da_oi_stats_buoy (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_buoy_type) :: stats      
   integer                :: nu, nv, nt, np, nq
   integer                :: n
   real                   :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_buoy")

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

   stats%average = residual_buoy1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(buoy)%nlocal
      if (iv%info(buoy)%proc_domain(1,n)) then

          u_inv = iv%buoy(n)%u%inv
          v_inv = iv%buoy(n)%v%inv
          u_obs = ob%buoy(n)%u
          v_obs = ob%buoy(n)%v

          if (.not. wind_sd_buoy .and. wind_stats_sd) &
             call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                     iv%buoy(n)%u%qc, iv%buoy(n)%v%qc, convert_uv2fd)
          if (wind_sd_buoy .and. .not. wind_stats_sd) &
             call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                     iv%buoy(n)%u%qc, iv%buoy(n)%v%qc, convert_fd2uv)

         call da_stats_calculate(iv%info(buoy)%obs_global_index(n), &
            0, iv%buoy(n)%u%qc, &
            u_inv, nu, &
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate(iv%info(buoy)%obs_global_index(n), &
            0, iv%buoy(n)%v%qc, &
            v_inv, nv, &
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate(iv%info(buoy)%obs_global_index(n), &
             0, iv%buoy(n)%t%qc, &
             iv%buoy(n)%t%inv, nt, &
             stats%minimum%t, stats%maximum%t, &
             stats%average%t, stats%rms_err%t)
         call da_stats_calculate(iv%info(buoy)%obs_global_index(n), &
            0, iv%buoy(n)%p%qc, &
            iv%buoy(n)%p%inv, np, &
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate(iv%info(buoy)%obs_global_index(n), &
             0, iv%buoy(n)%q%qc, &
             iv%buoy(n)%q%inv, nq, &
             stats%minimum%q, stats%maximum%q, &
             stats%average%q, stats%rms_err%q)
      end if    ! end if (iv%info(buoy)%proc_domain(1,n))
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for buoy'
         call da_print_stats_buoy(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_buoy")

end subroutine da_oi_stats_buoy


subroutine da_print_stats_buoy(stats_unit, nu, nv, nt, np, nq, buoy)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                intent(in)    :: stats_unit
   integer,                intent(inout) :: nu, nv, nt, np, nq
   type (stats_buoy_type), intent(in)    :: buoy

   if (trace_use_dull) call da_trace_entry("da_print_stats_buoy")

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
      ' Minimum(n,k): ', buoy%minimum%u, buoy%minimum%v, buoy%minimum%t, &
                         buoy%minimum%p, buoy%minimum%q, &
      ' Maximum(n,k): ', buoy%maximum%u, buoy%maximum%v, buoy%maximum%t, &
                         buoy%maximum%p, buoy%maximum%q
   write(unit=stats_unit, fmt='((a,4(f12.4,10x),e12.4,10x))') &
      ' Average     : ', buoy%average%u/real(nu), buoy%average%v/real(nv), &
                         buoy%average%t/real(nt), buoy%average%p/real(np), &
                         buoy%average%q/real(nq), &
      '    RMSE     : ', sqrt(buoy%rms_err%u/real(nu)), &
                         sqrt(buoy%rms_err%v/real(nv)), &
                         sqrt(buoy%rms_err%t/real(nt)), &
                         sqrt(buoy%rms_err%p/real(np)), &
                         sqrt(buoy%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_buoy")

end subroutine da_print_stats_buoy


subroutine da_transform_xtoy_buoy (grid, iv, y)

   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !--------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n        ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_psfc(:)

   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_buoy")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_v(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_t(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_q(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_psfc(iv%info(buoy)%n1:iv%info(buoy)%n2))

      allocate (ub(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (vb(1,iv%info(buoy)%n1:iv%info(buoy)%n2))

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xa%u, iv%info(buoy), model_u)
      call da_interp_lin_3d (grid%xa%v, iv%info(buoy), model_v)
      call da_interp_lin_3d (grid%xa%t, iv%info(buoy), model_t)
      call da_interp_lin_3d (grid%xa%q, iv%info(buoy), model_q)

      call da_interp_lin_2d (grid%xa%psfc, iv%info(buoy), 1, model_psfc)

      call da_interp_lin_3d (grid%xb%u, iv%info(buoy), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(buoy), vb)

      do n=iv%info(buoy)%n1,iv%info(buoy)%n2
         if(wind_sd_buoy)then
            call da_uv_to_sd_lin(y%buoy(n)%u,y%buoy(n)%v,model_u(1,n),model_v(1,n),ub(1,n),vb(1,n))
         else
            y%buoy(n)%u = model_u(1,n)
            y%buoy(n)%v = model_v(1,n)
         end if
         y%buoy(n)%t = model_t(1,n)
         y%buoy(n)%q = model_q(1,n)
         y%buoy(n)%p = model_psfc(n)
      end do
      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      ! [2.0] Surface assmiilation approach 2
      call da_transform_xtopsfc(grid,iv,buoy,iv%buoy(:),y%buoy(:))
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_buoy")

end subroutine da_transform_xtoy_buoy


subroutine da_transform_xtoy_buoy_adj(grid, iv, jo_grad_y, jo_grad_x)

   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !--------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_buoy_adj")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_v(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_t(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_q(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (model_psfc(iv%info(buoy)%n1:iv%info(buoy)%n2))

      allocate (ub(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
      allocate (vb(1,iv%info(buoy)%n1:iv%info(buoy)%n2))

      call da_interp_lin_3d (grid%xb%u, iv%info(buoy), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(buoy), vb)

      ! [1.2] Interpolate horizontally:
      do n=iv%info(buoy)%n1,iv%info(buoy)%n2
         if(wind_sd_buoy)then
            call da_uv_to_sd_adj(jo_grad_y%buoy(n)%u, &
                                 jo_grad_y%buoy(n)%v, model_u(1,n), model_v(1,n), ub(1,n), vb(1,n))
         else
            model_u(1,n)  = jo_grad_y%buoy(n)%u
            model_v(1,n)  = jo_grad_y%buoy(n)%v
         end if

         model_t(1,n)  = jo_grad_y%buoy(n)%t
         model_q(1,n)  = jo_grad_y%buoy(n)%q
         model_psfc(n) = jo_grad_y%buoy(n)%p
      end do
      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(buoy), model_u)
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(buoy), model_v)
      call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(buoy), model_t)
      call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(buoy), model_q)

      call da_interp_lin_2d_adj (jo_grad_x%psfc, iv%info(buoy), 1, model_psfc)
      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc_adj(grid,iv, buoy, iv%buoy(:), jo_grad_y%buoy(:), jo_grad_x)
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_buoy_adj")

end subroutine da_transform_xtoy_buoy_adj


subroutine da_check_max_iv_buoy(iv,ob, it, num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration 
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)
   type(y_type),  intent(in)    :: ob      ! Observation structure.

   logical :: failed,failed1,failed2
   integer :: n
   
   if (trace_use_dull) call da_trace_entry("da_check_max_iv_buoy")       


   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

      do n=iv%info(buoy)%n1,iv%info(buoy)%n2
         if(.not. qc_rej_both)then
            if(wind_sd_buoy)then
               failed=.false.
               if( iv%buoy(n)%u%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%u, max_error_spd,failed)
                   if( iv%info(buoy)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,buoy,1,1) = num_qcstat_conv(1,buoy,1,1) + 1
                       if(failed) then
                          num_qcstat_conv(2,buoy,1,1) = num_qcstat_conv(2,buoy,1,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'buoy',ob_vars(1),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
                       end if
                   end if
                end if

                failed=.false.
                if( iv%buoy(n)%v%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%v, max_error_dir,failed)
                    if( iv%info(buoy)%proc_domain(1,n) ) then
                        num_qcstat_conv(1,buoy,2,1) = num_qcstat_conv(1,buoy,2,1) + 1
                        if(failed)then
                           num_qcstat_conv(2,buoy,2,1) = num_qcstat_conv(2,buoy,2,1) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'buoy',ob_vars(2),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
                        end if
                    end if
                end if

             else

                failed=.false.
                if( iv%buoy(n)%u%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%u, max_error_uv,failed)
                    if( iv%info(buoy)%proc_domain(1,n) ) then
                        num_qcstat_conv(1,buoy,1,1) = num_qcstat_conv(1,buoy,1,1) + 1
                        if(failed) then
                           num_qcstat_conv(2,buoy,1,1) = num_qcstat_conv(2,buoy,1,1) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'buoy',ob_vars(1),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
                        end if
                    end if
                 end if

                 failed=.false.
                 if( iv%buoy(n)%v%qc >= obs_qc_pointer ) then
                     call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%v, max_error_uv,failed)
                     if( iv%info(buoy)%proc_domain(1,n) ) then
                         num_qcstat_conv(1,buoy,2,1) = num_qcstat_conv(1,buoy,2,1) + 1
                         if(failed)then
                            num_qcstat_conv(2,buoy,2,1) = num_qcstat_conv(2,buoy,2,1) + 1
                            write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                            'buoy',ob_vars(2),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
                         end if
                     end if
                 end if
              end if
              if(wind_sd_buoy)then
                 if(iv%buoy(n)%u%qc == fails_error_max .or.  abs(iv%buoy(n)%u%inv) >= max_omb_spd) then
                    iv%buoy(n)%u%qc = fails_error_max
                    iv%buoy(n)%u%inv = 0.0
                 endif
                 if(iv%buoy(n)%v%qc == fails_error_max .or.  abs(iv%buoy(n)%v%inv) >= max_omb_dir) then
                    iv%buoy(n)%v%qc = fails_error_max
                    iv%buoy(n)%v%inv = 0.0
                 endif
              endif

           else
              failed1=.false.
              failed2=.false.

              if( iv%buoy(n)%v%qc >= obs_qc_pointer .or. iv%buoy(n)%u%qc >= obs_qc_pointer )  then
                  if(wind_sd_buoy)then
                     call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%u, max_error_spd,failed1)
                     call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%v, max_error_dir,failed2)
                  else
                     call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%u, max_error_uv,failed1)
                     call da_max_error_qc (it,iv%info(buoy), n, iv%buoy(n)%v, max_error_uv,failed2)
                  endif
              endif
                        
              if( iv%info(buoy)%proc_domain(1,n) ) then
                  num_qcstat_conv(1,buoy,1,1) = num_qcstat_conv(1,buoy,1,1) + 1
                  num_qcstat_conv(1,buoy,2,1) = num_qcstat_conv(1,buoy,2,1) + 1

                  if(failed1 .or. failed2) then
                     num_qcstat_conv(2,buoy,1,1) = num_qcstat_conv(2,buoy,1,1) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'buoy',ob_vars(1),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
                     num_qcstat_conv(2,buoy,2,1) = num_qcstat_conv(2,buoy,2,1) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'buoy',ob_vars(2),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
                  endif
               endif

               if(wind_sd_buoy)then
                  if(iv%buoy(n)%u%qc == fails_error_max .or. iv%buoy(n)%v%qc == fails_error_max .or. &
                     abs(iv%buoy(n)%v%inv) >= max_omb_dir .or. abs(iv%buoy(n)%u%inv) >= max_omb_spd )then
                     iv%buoy(n)%u%qc = fails_error_max
                     iv%buoy(n)%v%qc = fails_error_max
                     iv%buoy(n)%u%inv = 0.0
                     iv%buoy(n)%v%inv = 0.0
                  endif
               else
                  if(iv%buoy(n)%u%qc == fails_error_max .or. iv%buoy(n)%v%qc == fails_error_max ) then
                     iv%buoy(n)%u%qc = fails_error_max
                     iv%buoy(n)%v%qc = fails_error_max
                     iv%buoy(n)%u%inv = 0.0
                     iv%buoy(n)%v%inv = 0.0
                  endif
               endif
            endif

      failed=.false.
      if( iv%buoy(n)%t%qc >= obs_qc_pointer )  then 
      call da_max_error_qc (it, iv%info(buoy), n, iv%buoy(n)%t, max_error_t , failed)
      if( iv%info(buoy)%proc_domain(1,n) ) then
      num_qcstat_conv(1,buoy,3,1)= num_qcstat_conv(1,buoy,3,1) + 1
      if(failed) then
      num_qcstat_conv(2,buoy,3,1)= num_qcstat_conv(2,buoy,3,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'buoy',ob_vars(3),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%buoy(n)%p%qc >= obs_qc_pointer )  then 
      call da_max_error_qc (it, iv%info(buoy), n, iv%buoy(n)%p, max_error_p , failed)         
      if( iv%info(buoy)%proc_domain(1,n) ) then
      num_qcstat_conv(1,buoy,5,1)= num_qcstat_conv(1,buoy,5,1) + 1
      if(failed) then
      num_qcstat_conv(2,buoy,5,1)= num_qcstat_conv(2,buoy,5,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'buoy',ob_vars(5),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
      end if
      end if
      end if


      failed=.false.
      if( iv%buoy(n)%q%qc >= obs_qc_pointer ) then
       if( iv%buoy(n)%t%qc == fails_error_max .or. iv%buoy(n)%p%qc == fails_error_max) then
       failed=.true.
       iv%buoy(n)%q%qc  = fails_error_max
       iv%buoy(n)%q%inv = 0.0
       else
       call da_max_error_qc (it, iv%info(buoy), n, iv%buoy(n)%q, max_error_q , failed)
       endif
      if( iv%info(buoy)%proc_domain(1,n) ) then
      num_qcstat_conv(1,buoy,4,1)= num_qcstat_conv(1,buoy,4,1) + 1
      if(failed) then
      num_qcstat_conv(2,buoy,4,1)= num_qcstat_conv(2,buoy,4,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'buoy',ob_vars(4),iv%info(buoy)%lat(1,n),iv%info(buoy)%lon(1,n),0.01*ob%buoy(n)%p
      end if
      end if
      end if

   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_max_iv_buoy")       

end subroutine da_check_max_iv_buoy
subroutine da_get_innov_vector_buoy( it,num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD    
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid      ! first guess state.
   type(y_type),     intent(inout) :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real    :: speed, direction

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

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_buoy")

   allocate (model_u(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
   allocate (model_v(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
   allocate (model_t(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
   allocate (model_p(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
   allocate (model_q(1,iv%info(buoy)%n1:iv%info(buoy)%n2))
   allocate (model_hsm(1,iv%info(buoy)%n1:iv%info(buoy)%n2))

   if ( it > 1 ) then
      do n=iv%info(buoy)%n1,iv%info(buoy)%n2
         if (iv%buoy(n)%u%qc == fails_error_max) iv%buoy(n)%u%qc = 0
         if (iv%buoy(n)%v%qc == fails_error_max) iv%buoy(n)%v%qc = 0
         if (iv%buoy(n)%t%qc == fails_error_max) iv%buoy(n)%t%qc = 0
         if (iv%buoy(n)%p%qc == fails_error_max) iv%buoy(n)%p%qc = 0
         if (iv%buoy(n)%q%qc == fails_error_max) iv%buoy(n)%q%qc = 0
      end do
   end if

   if (sfc_assi_options == sfc_assi_options_1) then
      do n=iv%info(buoy)%n1,iv%info(buoy)%n2
         ! [1.1] Get horizontal interpolation weights:

         i   = iv%info(buoy)%i(1,n)
         j   = iv%info(buoy)%j(1,n)
         dx  = iv%info(buoy)%dx(1,n)
         dy  = iv%info(buoy)%dy(1,n)
         dxm = iv%info(buoy)%dxm(1,n)
         dym = iv%info(buoy)%dym(1,n)

         ! Surface correction

         iv%buoy(n)%p%inv = ob%buoy(n)%p
         iv%buoy(n)%t%inv = ob%buoy(n)%t
         iv%buoy(n)%q%inv = ob%buoy(n)%q
         iv%buoy(n)%u%inv = ob%buoy(n)%u
         iv%buoy(n)%v%inv = ob%buoy(n)%v

         if (iv % buoy(n) % h > missing_r) then
            do k=kts,kte
               v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                      + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
            end do

            hd = v_h(kts) - iv % buoy(n) % h
            if (abs(hd) <= Max_StHeight_Diff .or. anal_type_verify) then
               if (iv % buoy(n) % h < v_h(kts)) then
                  iv%info(buoy)%zk(:,n) = 1.0+1.0e-6
                  call da_obs_sfc_correction(iv%info(buoy), iv%buoy(n), n, grid%xb)

               else
                  call da_to_zk(iv % buoy(n) % h, v_h, v_interp_h, iv%info(buoy)%zk(1,n))
               end if
            end if
         else if (ob % buoy(n) % p > 1.0) then
            do k=kts,kte
               v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                       + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
            end do

            call da_to_zk(ob % buoy(n) % p, v_p, v_interp_p, iv%info(buoy)%zk(1,n))

            if (iv%info(buoy)%zk(1,n) < 0.0 .and.  .not.anal_type_verify) then
               iv % buoy(n) % p % inv = missing_r
               iv % buoy(n) % p % qc  = missing_data
               iv%info(buoy)%zk(:,n) = 1.0+1.0e-6
            end if
         end if
      end do

      call da_convert_zk (iv%info(buoy))

      if (.not.anal_type_verify) then
         do n=iv%info(buoy)%n1,iv%info(buoy)%n2
            if (iv%info(buoy)%zk(1,n) < 0.0) then
               iv % buoy(n) % u % qc = missing_data
               iv % buoy(n) % v % qc = missing_data
               iv % buoy(n) % t % qc = missing_data
               iv % buoy(n) % q % qc = missing_data
               iv % buoy(n) % p % qc = missing_data
            end if
         end do
      end if

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xb%u, iv%info(buoy),model_u)
      call da_interp_lin_3d (grid%xb%v, iv%info(buoy),model_v)
      call da_interp_lin_3d (grid%xb%t, iv%info(buoy),model_t)
      call da_interp_lin_3d (grid%xb%q, iv%info(buoy),model_q)
      call da_interp_lin_3d (grid%xb%p, iv%info(buoy),model_p)
   else if (sfc_assi_options == sfc_assi_options_2) then

      ! Surface data assimilation approach 2
      ! -----------------------------------

      ! 1.2.1 Surface assmiilation approach 2(10-m u, v, 2-m t, q, 
      ! and sfc_p)

      call da_interp_lin_2d (grid%xb%u10,  iv%info(buoy), 1,model_u)
      call da_interp_lin_2d (grid%xb%v10,  iv%info(buoy), 1,model_v)
      call da_interp_lin_2d (grid%xb%t2,   iv%info(buoy), 1,model_t)
      call da_interp_lin_2d (grid%xb%q2,   iv%info(buoy), 1,model_q)
      call da_interp_lin_2d (grid%xb%psfc, iv%info(buoy), 1,model_p)

      do n=iv%info(buoy)%n1,iv%info(buoy)%n2

         iv%buoy(n)%p%inv = ob%buoy(n)%p
         iv%buoy(n)%t%inv = ob%buoy(n)%t
         iv%buoy(n)%q%inv = ob%buoy(n)%q
         iv%buoy(n)%u%inv = ob%buoy(n)%u
         iv%buoy(n)%v%inv = ob%buoy(n)%v

         if (iv%buoy(n)%p%qc >= 0) then
            ! model surface p, t, q, h at observed site:

            call da_interp_lin_2d_partial (grid%xb%terr, iv%info(buoy), 1, n, n, model_hsm(:,n))

            ho = iv%buoy(n)%h
            to = -888888.0
            qo = -888888.0

            if (iv%buoy(n)%t%qc >= 0 .and. iv%buoy(n)%q%qc >= 0) then
               to = ob%buoy(n)%t
               qo = ob%buoy(n)%q
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to, qo)
            else if (iv%buoy(n)%t%qc >= 0 .and. iv%buoy(n)%q%qc < 0) then
               to = ob%buoy(n)%t
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to)
            else
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho)
            end if

            ! Pressure at the observed height:
            model_p(1,n) = psfcm
         end if
      end do
   end if

   do n=iv%info(buoy)%n1,iv%info(buoy)%n2

      !-----------------------------------------------------------------------
      ! [3.0] Fast interpolation: 
      !-----------------------------------------------------------------------
   if(wind_sd_buoy)then
      call da_ffdduv_model (speed,direction,model_u(1,n), model_v(1,n), convert_uv2fd)

      if (ob%buoy(n)%u > missing_r .AND. iv%buoy(n)%u%qc >= obs_qc_pointer) then
          iv%buoy(n)%u%inv = iv%buoy(n)%u%inv - speed
      else
         iv % buoy(n) % u % inv = 0.0
      end if

      if (ob%buoy(n)%v > missing_r .AND. iv%buoy(n)%v%qc >= obs_qc_pointer) then
          iv%buoy(n)%v%inv = iv%buoy(n)%v%inv - direction
          if (iv%buoy(n)%v%inv > 180.0 ) iv%buoy(n)%v%inv = iv%buoy(n)%v%inv - 360.0
          if (iv%buoy(n)%v%inv < -180.0 ) iv%buoy(n)%v%inv = iv%buoy(n)%v%inv + 360.0
      else
         iv % buoy(n) % v % inv = 0.0
      end if
    else
      if (ob % buoy(n) % u > missing_r .AND. iv % buoy(n) % u % qc >= obs_qc_pointer) then
         iv % buoy(n) % u % inv = iv%buoy(n)%u%inv - model_u(1,n)
      else
         iv % buoy(n) % u % inv = 0.0
      end if

      if (ob % buoy(n) % v > missing_r .AND. iv % buoy(n) % v % qc >= obs_qc_pointer) then
         iv % buoy(n) % v % inv = iv%buoy(n)%v%inv - model_v(1,n)
      else
         iv % buoy(n) % v % inv = 0.0
      end if
    end if

      !if (ob % buoy(n) % p > 0.0 .AND. iv % buoy(n) % p % qc >= obs_qc_pointer) then
      if ( iv % buoy(n) % p % qc >= obs_qc_pointer ) then
         iv % buoy(n) % p % inv = iv%buoy(n)%p%inv - model_p(1,n)
      else
         iv % buoy(n) % p % inv = 0.0
      end if

      if (ob % buoy(n) % t > 0.0 .AND. iv % buoy(n) % t % qc >= obs_qc_pointer) then
         iv % buoy(n) % t % inv = iv%buoy(n)%t%inv - model_t(1,n)
      else
         iv % buoy(n) % t % inv = 0.0
      end if

      if (ob % buoy(n) % q > 0.0 .AND. iv % buoy(n) % q % qc >= obs_qc_pointer) then
         iv % buoy(n) % q % inv = iv%buoy(n)%q%inv - model_q(1,n)
      else
         iv % buoy(n) % q % inv = 0.0
      end if
   end do

   ! -----------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_buoy(iv,ob, it, num_qcstat_conv)

   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_p)
   deallocate (model_q)
   deallocate (model_hsm)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_buoy")

end subroutine da_get_innov_vector_buoy


subroutine da_calculate_grady_buoy(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_buoy")       

   do n=1, iv%info(buoy)%nlocal
      if (iv%buoy(n)%u%qc < obs_qc_pointer) re%buoy(n)%u = 0.0
      if (iv%buoy(n)%v%qc < obs_qc_pointer) re%buoy(n)%v = 0.0
      if (iv%buoy(n)%t%qc < obs_qc_pointer) re%buoy(n)%t = 0.0
      if (iv%buoy(n)%p%qc < obs_qc_pointer) re%buoy(n)%p = 0.0
      if (iv%buoy(n)%q%qc < obs_qc_pointer) re%buoy(n)%q = 0.0

      jo_grad_y%buoy(n)%u = -re%buoy(n)%u / (iv%buoy(n)%u%error * iv%buoy(n)%u%error)
      jo_grad_y%buoy(n)%v = -re%buoy(n)%v / (iv%buoy(n)%v%error * iv%buoy(n)%v%error)
      jo_grad_y%buoy(n)%t = -re%buoy(n)%t / (iv%buoy(n)%t%error * iv%buoy(n)%t%error)
      jo_grad_y%buoy(n)%p = -re%buoy(n)%p / (iv%buoy(n)%p%error * iv%buoy(n)%p%error)
      jo_grad_y%buoy(n)%q = -re%buoy(n)%q / (iv%buoy(n)%q%error * iv%buoy(n)%q%error)
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_buoy")  
     
end subroutine da_calculate_grady_buoy




end module da_buoy 

