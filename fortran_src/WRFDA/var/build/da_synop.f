












module da_synop

   use module_domain, only : domain

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, missing_data, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, sfc_assi_options, no_buddies, fails_error_max, &
      fails_buddy_check, check_buddy, check_buddy_print, check_buddy_unit, &
      buddy_weight , max_buddy_uv, max_buddy_t, max_buddy_p, max_buddy_rh, &
      max_stheight_diff,test_dm_exact, anal_type_verify, &
      kts,kte,kms,kme,sfc_assi_options_1,sfc_assi_options_2 , &
      trace_use_dull, synop, max_ext_its,qcstat_conv_unit,ob_vars, &
      convert_fd2uv, convert_uv2fd, max_error_spd, max_error_dir, &
      max_omb_spd, max_omb_dir, pi, qc_rej_both, &
      wind_sd_synop, wind_stats_sd 
   use da_grid_definitions, only : da_ffdduv, da_ffdduv_model, da_ffdduv_diagnose 
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_to_zk, &
      da_interp_lin_3d,da_interp_lin_3d_adj, &
      da_interp_lin_2d, da_interp_lin_2d_adj
   use da_par_util1, only : da_proc_sum_int
   use da_par_util, only : da_proc_stats_combine, &
      da_deallocate_global_synop, da_to_global_synop
   use da_physics, only : da_sfc_pre, da_transform_xtopsfc, &
      da_transform_xtopsfc_adj, da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_obs_sfc_correction, &
                        da_buddy_qc, da_convert_zk
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_synop1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: p                        
      real          :: q                        
   end type residual_synop1_type

   type maxmin_synop_stats_type
      type (maxmin_type)         :: u, v, t, p, q
   end type maxmin_synop_stats_type

   type stats_synop_type
      type (maxmin_synop_stats_type)  :: maximum, minimum
      type (residual_synop1_type)     :: average, rms_err
   end type stats_synop_type

contains

subroutine da_ao_stats_synop (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type  (y_type), intent (in)    :: re            ! A - O
   type(y_type),   intent (in)    :: ob            ! Observation structure.

   type (stats_synop_type) :: stats
   integer                 :: nu, nv, nt, np, nq
   integer                 :: n
   real                    :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_synop")

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

   stats%average = residual_synop1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(synop)%nlocal
      if (iv%info(synop)%proc_domain(1,n)) then

         u_inc = re%synop(n)%u
         v_inc = re%synop(n)%v
         u_obs = ob%synop(n)%u
         v_obs = ob%synop(n)%v

         if (.not. wind_sd_synop .and. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%synop(n)%u%qc, iv%synop(n)%v%qc, convert_uv2fd)
         if (wind_sd_synop .and. .not. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%synop(n)%u%qc, iv%synop(n)%v%qc, convert_fd2uv)

         call da_stats_calculate (n, 0, iv%synop(n)%u%qc,  & 
            u_inc, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%synop(n)%v%qc,  & 
            v_inc, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate (n, 0, iv%synop(n)%t%qc,  & 
            re%synop(n)%t, nt, & 
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate (n, 0, iv%synop(n)%p%qc,  & 
            re%synop(n)%p, np, & 
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate (n, 0, iv%synop(n)%q%qc,  & 
            re%synop(n)%q, nq, & 
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)
      end if    ! end if (iv%info(synop)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (np)
   call da_proc_sum_int (nq)
   iv%nstats(synop) = nu + nv + nt + np + nq

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for synop'
         call da_print_stats_synop(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_synop")

end subroutine da_ao_stats_synop


subroutine da_jo_and_grady_synop(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type),  intent(in)   :: re          ! Residual vector.
   type (y_type),  intent(inout):: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout):: jo          ! Obs cost function.

   integer        :: n

   ! the following "global" objects are used only when testing
   type (iv_type) :: iv_glob         ! Global Innovation vector (O-B).
   type (y_type)  :: re_glob         ! Global Residual vector (O-A).
   type (y_type)  :: jo_grad_y_glob  ! Global Grad_y(Jo)

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_synop")

   jo % synop_u = 0.0
   jo % synop_v = 0.0
   jo % synop_t = 0.0
   jo % synop_p = 0.0
   jo % synop_q = 0.0

   if (test_dm_exact) then
      if (iv%info(synop)%ntotal == 0) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_synop")
         return
      end if
   else
      if (iv%info(synop)%nlocal < 1) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_synop")
         return
      end if
   end if

   do n=1, iv%info(synop)%nlocal
      jo_grad_y%synop(n)%u = -re%synop(n)%u / (iv%synop(n)%u%error * iv%synop(n)%u%error)
      jo_grad_y%synop(n)%v = -re%synop(n)%v / (iv%synop(n)%v%error * iv%synop(n)%v%error)
      jo_grad_y%synop(n)%t = -re%synop(n)%t / (iv%synop(n)%t%error * iv%synop(n)%t%error)
      jo_grad_y%synop(n)%p = -re%synop(n)%p / (iv%synop(n)%p%error * iv%synop(n)%p%error)
      jo_grad_y%synop(n)%q = -re%synop(n)%q / (iv%synop(n)%q%error * iv%synop(n)%q%error)
   end do

   ! Bitwise-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not
   ! testing.  This will always be .false. for a serial run.
   if (test_dm_exact) then
      ! collect all obs in serial order and allocate global objects
      call da_to_global_synop(iv, re, jo_grad_y, iv_glob, re_glob, jo_grad_y_glob)
      ! perform remaining computations
      call da_jo_synop_uvtq(iv_glob, re_glob, jo_grad_y_glob, jo)
      ! free global objects
      call da_deallocate_global_synop(iv_glob, re_glob, jo_grad_y_glob)
   else
      ! perform remaining computations
      call da_jo_synop_uvtq(iv, re, jo_grad_y, jo)
   end if

   jo % synop_u = 0.5 * jo % synop_u
   jo % synop_v = 0.5 * jo % synop_v
   jo % synop_t = 0.5 * jo % synop_t
   jo % synop_p = 0.5 * jo % synop_p
   jo % synop_q = 0.5 * jo % synop_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_synop")

end subroutine da_jo_and_grady_synop


subroutine da_jo_synop_uvtq(iv, re, jo_grad_y, jo)

   !-----------------------------------------------------------------------
   ! Purpose: Ensures that exactly the same code is used for both
   ! serial and parallel computations when testing_dm_bitwise_exact is set.
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv         ! Innovation vector.
   type (y_type),  intent(in)    :: re         ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y  ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo         ! Obs cost function.

   integer :: n

   if (trace_use_dull) call da_trace_entry("da_jo_synop_uvtq")

   do n=1, iv%info(synop)%nlocal
      if (iv%info(synop)%proc_domain(1,n)) then
         jo % synop_u = jo % synop_u - re%synop(n)%u * jo_grad_y%synop(n)%u
         jo % synop_v = jo % synop_v - re%synop(n)%v * jo_grad_y%synop(n)%v
         jo % synop_t = jo % synop_t - re%synop(n)%t * jo_grad_y%synop(n)%t
         jo % synop_p = jo % synop_p - re%synop(n)%p * jo_grad_y%synop(n)%p
         jo % synop_q = jo % synop_q - re%synop(n)%q * jo_grad_y%synop(n)%q
     end if
  end do

   if (trace_use_dull) call da_trace_exit("da_jo_synop_uvtq")

end subroutine da_jo_synop_uvtq


subroutine da_residual_synop(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for synop obs
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

   if (trace_use_dull) call da_trace_entry("da_residual_synop")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % p % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(synop)%nlocal
      np_available = np_available + 5

      re%synop(n)%u = da_residual(n, 0, y%synop(n)%u, iv%synop(n)%u, n_obs_bad % u)
      re%synop(n)%v = da_residual(n, 0, y%synop(n)%v, iv%synop(n)%v, n_obs_bad % v)
      re%synop(n)%t = da_residual(n, 0, y%synop(n)%t, iv%synop(n)%t, n_obs_bad % t)
      re%synop(n)%p = da_residual(n, 0, y%synop(n)%p, iv%synop(n)%p, n_obs_bad % p)
      re%synop(n)%q = da_residual(n, 0, y%synop(n)%q, iv%synop(n)%q, n_obs_bad % q)
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

   if (trace_use_dull) call da_trace_exit("da_residual_synop")

end subroutine da_residual_synop


subroutine da_oi_stats_synop (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_synop_type) :: stats
   integer                 :: nu, nv, nt, np, nq
   integer                 :: n
   real                    :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_synop")

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

   stats%average = residual_synop1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(synop)%nlocal
      if (iv%info(synop)%proc_domain(1,n)) then

        u_inv = iv%synop(n)%u%inv
        v_inv = iv%synop(n)%v%inv
        u_obs = ob%synop(n)%u
        v_obs = ob%synop(n)%v

        if (.not. wind_sd_synop .and. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%synop(n)%u%qc, iv%synop(n)%v%qc, convert_uv2fd)
        if (wind_sd_synop .and. .not. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%synop(n)%u%qc, iv%synop(n)%v%qc, convert_fd2uv)

        call da_stats_calculate(iv%info(synop)%obs_global_index(n), &
           0, iv%synop(n)%u%qc, &
           u_inv, nu, &
           stats%minimum%u, stats%maximum%u, &
           stats%average%u, stats%rms_err%u)
        call da_stats_calculate(iv%info(synop)%obs_global_index(n), &
           0, iv%synop(n)%v%qc, &
           v_inv, nv, &
           stats%minimum%v, stats%maximum%v, &
           stats%average%v, stats%rms_err%v)
        call da_stats_calculate(iv%info(synop)%obs_global_index(n), &
           0, iv%synop(n)%t%qc, &
           iv%synop(n)%t%inv, nt, &
           stats%minimum%t, stats%maximum%t, &
           stats%average%t, stats%rms_err%t)
        call da_stats_calculate(iv%info(synop)%obs_global_index(n), &
           0, iv%synop(n)%p%qc, &
           iv%synop(n)%p%inv, np, &
           stats%minimum%p, stats%maximum%p, &
           stats%average%p, stats%rms_err%p)
        call da_stats_calculate(iv%info(synop)%obs_global_index(n), &
            0, iv%synop(n)%q%qc, &
            iv%synop(n)%q%inv, nq, &
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for synop'
         call da_print_stats_synop(stats_unit, nu, nv, nt, np, nq, stats)

         ! min_val = -50.0
         ! max_val = 50.0
         ! bin_width = 1.0
         ! call da_data_distribution('synop u ', iv%total(synop), min_val, max_val, &
         !                            bin_width, iv%synop(:)%u%inv)
         ! call da_data_distribution('synop v ', iv%total(synop), min_val, max_val, &
         !                            bin_width, iv%synop(:)%v%inv)
         ! call da_data_distribution('synop t ', iv%total(synop), min_val, max_val, &
         !                            bin_width, iv%synop(:)%t%inv)
         ! call da_data_distribution('synop p ', iv%total(synop), -1000.0, 1000.0, &
         !                            50.0, iv%synop(:)%p%inv)
         ! call da_data_distribution('synop q ', iv%total(synop), -0.03, 0.03, &
         !                               0.001, iv%synop(:)%q%inv)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_synop")

end subroutine da_oi_stats_synop


subroutine da_print_stats_synop(stats_unit, nu, nv, nt, np, nq, stats)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nu, nv, nt, np, nq
   type (stats_synop_type), intent(in)    :: stats

   if (trace_use_dull) call da_trace_entry("da_print_stats_synop")

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
     ' Minimum(n,k): ', stats%minimum%u, stats%minimum%v, stats%minimum%t, &
                        stats%minimum%p, stats%minimum%q, &
     ' Maximum(n,k): ', stats%maximum%u, stats%maximum%v, stats%maximum%t, &
                        stats%maximum%p, stats%maximum%q
   write(unit=stats_unit, fmt='((a,4(f12.4,10x),e12.4,10x))') &
      ' Average     : ', stats%average%u/real(nu), stats%average%v/real(nv), &
                         stats%average%t/real(nt), stats%average%p/real(np), &
                         stats%average%q/real(nq), &
      '    RMSE     : ', sqrt(stats%rms_err%u/real(nu)), &
                         sqrt(stats%rms_err%v/real(nv)), &
                         sqrt(stats%rms_err%t/real(nt)), &
                         sqrt(stats%rms_err%p/real(np)), &
                         sqrt(stats%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_synop")

end subroutine da_print_stats_synop


subroutine da_transform_xtoy_synop (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n       ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_psfc(:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)
   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_synop")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_v(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_t(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_q(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_psfc(iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (ub(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (vb(1,iv%info(synop)%n1:iv%info(synop)%n2))

      call da_interp_lin_3d (grid%xb%u, iv%info(synop), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(synop), vb)

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xa%u, iv%info(synop), model_u)
      call da_interp_lin_3d (grid%xa%v, iv%info(synop), model_v)
      call da_interp_lin_3d (grid%xa%t, iv%info(synop), model_t)
      call da_interp_lin_3d (grid%xa%q, iv%info(synop), model_q)

      call da_interp_lin_2d(grid%xa%psfc, iv%info(synop), 1, model_psfc)

      do n=iv%info(synop)%n1,iv%info(synop)%n2
         if(wind_sd_synop)then
            call da_uv_to_sd_lin(y%synop(n)%u,y%synop(n)%v,model_u(1,n),model_v(1,n),ub(1,n),vb(1,n))
         else
            y%synop(n)%u = model_u(1,n)
            y%synop(n)%v = model_v(1,n)
         end if
         y%synop(n)%t = model_t(1,n)
         y%synop(n)%q = model_q(1,n)
         y%synop(n)%p = model_psfc(n)
      end do
      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc(grid, iv, synop, iv%synop(:), y%synop(:))
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_synop")

end subroutine da_transform_xtoy_synop


subroutine da_transform_xtoy_synop_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

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
   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_synop_adj")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_v(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_t(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_q(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (model_psfc(iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (ub(1,iv%info(synop)%n1:iv%info(synop)%n2))
      allocate (vb(1,iv%info(synop)%n1:iv%info(synop)%n2))

      call da_interp_lin_3d (grid%xb%u, iv%info(synop), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(synop), vb)

      ! [1.2] Interpolate horizontally:
      do n=iv%info(synop)%n1,iv%info(synop)%n2
         if(wind_sd_synop)then
            call da_uv_to_sd_adj(jo_grad_y%synop(n)%u, &
                                 jo_grad_y%synop(n)%v, model_u(1,n), model_v(1,n), ub(1,n), vb(1,n))
         else
            model_u(1,n)  = jo_grad_y%synop(n)%u
            model_v(1,n)  = jo_grad_y%synop(n)%v
         end if
         model_t(1,n)  = jo_grad_y%synop(n)%t
         model_q(1,n)  = jo_grad_y%synop(n)%q
         model_psfc(n) = jo_grad_y%synop(n)%p
      end do
      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(synop),model_u)
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(synop),model_v)
      call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(synop),model_t)
      call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(synop),model_q)

      call da_interp_lin_2d_adj (jo_grad_x%psfc, iv%info(synop), 1, model_psfc)
      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)
   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc_adj(grid,iv, synop, iv%synop(:), jo_grad_y%synop(:),jo_grad_x)
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_synop_adj")

end subroutine da_transform_xtoy_synop_adj


subroutine da_get_innov_vector_synop( it,num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   integer,       intent(in)    :: it      ! External iteration.
   type(domain),  intent(in)    :: grid    ! model state
   type(y_type),  intent(inout) :: ob      ! Observation structure.
   type(iv_type), intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.

   real    :: v_h(kms:kme)      ! Model value h at ob hor. loc
   real    :: v_p(kms:kme)      ! Model value p at ob hor. loc

   real    :: hd, psfcm

   real    :: ho, to, qo
   real  :: speed, direction

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_p(:,:)
   real, allocatable :: model_hsm(:,:)
    
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_synop")

   allocate (model_u  (1,iv%info(synop)%n1:iv%info(synop)%n2))
   allocate (model_v  (1,iv%info(synop)%n1:iv%info(synop)%n2))
   allocate (model_t  (1,iv%info(synop)%n1:iv%info(synop)%n2))
   allocate (model_q  (1,iv%info(synop)%n1:iv%info(synop)%n2))
   allocate (model_p  (1,iv%info(synop)%n1:iv%info(synop)%n2))
   allocate (model_hsm(1,iv%info(synop)%n1:iv%info(synop)%n2))

   if ( it > 1 ) then
      do n=iv%info(synop)%n1,iv%info(synop)%n2
         if (iv%synop(n)%u%qc == fails_error_max) iv%synop(n)%u%qc = 0
         if (iv%synop(n)%v%qc == fails_error_max) iv%synop(n)%v%qc = 0
         if (iv%synop(n)%t%qc == fails_error_max) iv%synop(n)%t%qc = 0
         if (iv%synop(n)%p%qc == fails_error_max) iv%synop(n)%p%qc = 0
         if (iv%synop(n)%q%qc == fails_error_max) iv%synop(n)%q%qc = 0
      end do
   end if

   if (sfc_assi_options == sfc_assi_options_1) then

      do n=iv%info(synop)%n1,iv%info(synop)%n2
         ! [1.1] Get horizontal interpolation weights:
         i   = iv%info(synop)%i(1,n)
         j   = iv%info(synop)%j(1,n)
         dx  = iv%info(synop)%dx(1,n)
         dy  = iv%info(synop)%dy(1,n)
         dxm = iv%info(synop)%dxm(1,n)
         dym = iv%info(synop)%dym(1,n)

         ! Surface correction

         iv%synop(n)%p%inv = ob%synop(n)%p
         iv%synop(n)%t%inv = ob%synop(n)%t
         iv%synop(n)%q%inv = ob%synop(n)%q
         iv%synop(n)%u%inv = ob%synop(n)%u
         iv%synop(n)%v%inv = ob%synop(n)%v

         if (iv % synop(n) % h > missing_r) then
            do k=kts,kte
               v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                       + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
            end do

            hd = v_h(kts) - iv % synop(n) % h
            if (abs(hd) <= Max_StHeight_Diff .or.  anal_type_verify ) then
               if (iv % synop(n) % h < v_h(kts)) then
                  iv%info(synop)%zk(:,n) = 1.0+1.0e-6
                  call da_obs_sfc_correction(iv%info(synop), iv%synop(n), n, grid%xb)
               else
                  call da_to_zk(iv % synop(n) % h, v_h, v_interp_h, iv%info(synop)%zk(1,n))
               end if
            end if
         else if (ob % synop(n) % p > 1.0) then
            do k=kts,kte
              v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                       + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
            end do

            call da_to_zk(ob % synop(n) % p, v_p, v_interp_p, iv%info(synop)%zk(1,n))
            if (iv%info(synop)%zk(1,n) < 0.0 .and. .not.anal_type_verify ) then
               iv % synop(n) % p % inv = missing_r
               iv % synop(n) % p % qc  = missing_data
               iv%info(synop)%zk(:,n) = 1.0+1.0e-6
            end if
         end if
      end do

      call da_convert_zk (iv%info(synop))

      if (.not.anal_type_verify ) then
         do n=iv%info(synop)%n1,iv%info(synop)%n2
            if (iv%info(synop)%zk(1,n) < 0.0) then
               iv % synop(n) % u % qc = missing_data
               iv % synop(n) % v % qc = missing_data
               iv % synop(n) % t % qc = missing_data
               iv % synop(n) % q % qc = missing_data
               iv % synop(n) % p % qc = missing_data
            end if
         end do
      end if
            
      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xb%u, iv%info(synop), model_u)
      call da_interp_lin_3d (grid%xb%v, iv%info(synop), model_v)
      call da_interp_lin_3d (grid%xb%t, iv%info(synop), model_t)
      call da_interp_lin_3d (grid%xb%q, iv%info(synop), model_q)
      call da_interp_lin_3d (grid%xb%p, iv%info(synop), model_p)
   else if (sfc_assi_options == sfc_assi_options_2) then
      ! Surface data assimilation approach 2
      !------------------------------------

      ! 1.2.1 Surface assmiilation approach 2(10-m u, v, 2-m t, q, and sfc_p)

      call da_interp_lin_2d (grid%xb%u10,  iv%info(synop), 1, model_u(1,:))
      call da_interp_lin_2d (grid%xb%v10,  iv%info(synop), 1, model_v(1,:))
      call da_interp_lin_2d (grid%xb%t2,   iv%info(synop), 1, model_t(1,:))
      call da_interp_lin_2d (grid%xb%q2,   iv%info(synop), 1, model_q(1,:))
      call da_interp_lin_2d (grid%xb%psfc, iv%info(synop), 1, model_p(1,:))
      call da_interp_lin_2d (grid%xb%terr, iv%info(synop), 1, model_hsm(1,:))

      do n=iv%info(synop)%n1,iv%info(synop)%n2

         iv%synop(n)%p%inv = ob%synop(n)%p
         iv%synop(n)%t%inv = ob%synop(n)%t
         iv%synop(n)%q%inv = ob%synop(n)%q
         iv%synop(n)%u%inv = ob%synop(n)%u
         iv%synop(n)%v%inv = ob%synop(n)%v

         if (iv%synop(n)%p%qc >= 0) then
            ho = iv%synop(n)%h
            to = -888888.0
            qo = -888888.0

            if (iv%synop(n)%t%qc >= 0 .and. iv%synop(n)%q%qc >= 0) then
               to = ob%synop(n)%t
               qo = ob%synop(n)%q
               call da_sfc_pre (psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to, qo)
            else if (iv%synop(n)%t%qc >= 0 .and. iv%synop(n)%q%qc < 0) then
               to = ob%synop(n)%t
               call da_sfc_pre (psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to)
            else
               call da_sfc_pre (psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho)
            end if

            ! Pressure at the observed height:
            model_p(1,n) = psfcm
         end if
      end do
   end if

   do n=iv%info(synop)%n1,iv%info(synop)%n2

      !--------------------------------------------------------------------
      !     [3.0] Fast interpolation:
      !--------------------------------------------------------------------
      if(wind_sd_synop)then
         call da_ffdduv_model (speed,direction,model_u(1,n), model_v(1,n), convert_uv2fd)

         if (ob%synop(n)%u > missing_r .AND. iv%synop(n)%u%qc >= obs_qc_pointer) then
             iv%synop(n)%u%inv = iv%synop(n)%u%inv - speed
         else
             iv % synop(n) % u % inv = 0.0
         end if

         if (ob%synop(n)%v > missing_r .AND. iv%synop(n)%v%qc >= obs_qc_pointer) then
             iv%synop(n)%v%inv = iv%synop(n)%v%inv - direction
             if (iv%synop(n)%v%inv > 180.0 ) iv%synop(n)%v%inv = iv%synop(n)%v%inv - 360.0
             if (iv%synop(n)%v%inv < -180.0 ) iv%synop(n)%v%inv = iv%synop(n)%v%inv + 360.0
         else
             iv % synop(n) % v % inv = 0.0
         end if
      else
         if (ob % synop(n) % u > missing_r .AND. iv % synop(n) % u % qc >= obs_qc_pointer) then
             iv % synop(n) % u % inv = iv%synop(n)%u%inv - model_u(1,n)
         else
             iv % synop(n) % u % inv = 0.0
         end if

         if (ob % synop(n) % v > missing_r .AND. iv % synop(n) % v % qc >= obs_qc_pointer) then
             iv % synop(n) % v % inv = iv%synop(n)%v%inv - model_v(1,n)
         else
             iv % synop(n) % v % inv = 0.0
         end if
      end if

      !if (ob % synop(n) % p > 0.0 .AND. iv % synop(n) % p % qc >= obs_qc_pointer) then
      if ( iv % synop(n) % p % qc >= obs_qc_pointer ) then
         iv % synop(n) % p % inv = iv%synop(n)%p%inv - model_p(1,n)
      else
         iv % synop(n) % p % inv = 0.0
      end if

      if (ob % synop(n) % t > 0.0 .AND. iv % synop(n) % t % qc >= obs_qc_pointer) then
         iv % synop(n) % t % inv = iv%synop(n)%t%inv - model_t(1,n)
      else
         iv % synop(n) % t % inv = 0.0
      end if

      if (ob % synop(n) % q > 0.0 .AND. iv % synop(n) % q % qc >= obs_qc_pointer) then
         iv % synop(n) % q % inv = iv%synop(n)%q%inv - model_q(1,n)
      else
         iv % synop(n) % q % inv = 0.0
      end if
   end do

   !--------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !--------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_synop(iv,ob, it,num_qcstat_conv)

   if (check_buddy) call da_check_buddy_synop(iv, ob, grid%dx, it)
 !
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_q)
   deallocate (model_p)
   deallocate (model_hsm)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_synop")

end subroutine da_get_innov_vector_synop


subroutine da_check_max_iv_synop(iv,ob, it, num_qcstat_conv)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)
   type(y_type),  intent(in)    :: ob      ! Observation structure.


   logical :: failed,failed1,failed2
   integer :: n

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_synop")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(synop)%n1,iv%info(synop)%n2
         if(.not. qc_rej_both)then
            if(wind_sd_synop)then
               failed=.false.
               if( iv%synop(n)%u%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%u, max_error_spd,failed)
                   if( iv%info(synop)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,synop,1,1) = num_qcstat_conv(1,synop,1,1) + 1
                       if(failed) then
                          num_qcstat_conv(2,synop,1,1) = num_qcstat_conv(2,synop,1,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'synop',ob_vars(1),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
                       end if
                   end if
               end if

               failed=.false.
               if( iv%synop(n)%v%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%v, max_error_dir,failed)
                   if( iv%info(synop)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,synop,2,1) = num_qcstat_conv(1,synop,2,1) + 1
                       if(failed)then
                          num_qcstat_conv(2,synop,2,1) = num_qcstat_conv(2,synop,2,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'synop',ob_vars(2),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
                       end if
                   end if
               end if
            else
               failed=.false.
               if( iv%synop(n)%u%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%u, max_error_uv,failed)
                   if( iv%info(synop)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,synop,1,1) = num_qcstat_conv(1,synop,1,1) + 1
                       if(failed) then
                          num_qcstat_conv(2,synop,1,1) = num_qcstat_conv(2,synop,1,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'synop',ob_vars(1),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
                       end if
                   end if
               end if

               failed=.false.
               if( iv%synop(n)%v%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%v, max_error_uv,failed)
                   if( iv%info(synop)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,synop,2,1) = num_qcstat_conv(1,synop,2,1) + 1
                       if(failed)then
                          num_qcstat_conv(2,synop,2,1) = num_qcstat_conv(2,synop,2,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'synop',ob_vars(2),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
                       end if
                   end if
               end if
            end if

            if(wind_sd_synop)then
               if(iv%synop(n)%u%qc == fails_error_max .or.  abs(iv%synop(n)%u%inv) >= max_omb_spd) then
                  iv%synop(n)%u%qc = fails_error_max
                  iv%synop(n)%u%inv = 0.0
               endif
               if(iv%synop(n)%v%qc == fails_error_max .or.  abs(iv%synop(n)%v%inv) >= max_omb_dir) then
                  iv%synop(n)%v%qc = fails_error_max
                  iv%synop(n)%v%inv = 0.0
               endif
            endif
         else
            failed1=.false.
            failed2=.false.

            if( iv%synop(n)%v%qc >= obs_qc_pointer .or. iv%synop(n)%u%qc >= obs_qc_pointer )  then
                if(wind_sd_synop)then
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%u, max_error_spd,failed1)
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%v, max_error_dir,failed2)
                else
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%u, max_error_uv,failed1)
                   call da_max_error_qc (it,iv%info(synop), n, iv%synop(n)%v, max_error_uv,failed2)
                endif
            endif

            if( iv%info(synop)%proc_domain(1,n) ) then
                num_qcstat_conv(1,synop,1,1) = num_qcstat_conv(1,synop,1,1) + 1
                num_qcstat_conv(1,synop,2,1) = num_qcstat_conv(1,synop,2,1) + 1

                if(failed1 .or. failed2) then
                   num_qcstat_conv(2,synop,1,1) = num_qcstat_conv(2,synop,1,1) + 1
                   write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'synop',ob_vars(1),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
                   num_qcstat_conv(2,synop,2,1) = num_qcstat_conv(2,synop,2,1) + 1
                   write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'synop',ob_vars(2),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
                endif
             endif

             if(wind_sd_synop)then
                if(iv%synop(n)%u%qc == fails_error_max .or. iv%synop(n)%v%qc == fails_error_max .or. &
                   abs(iv%synop(n)%v%inv) >= max_omb_dir .or. abs(iv%synop(n)%u%inv) >= max_omb_spd )then
                   iv%synop(n)%u%qc = fails_error_max
                   iv%synop(n)%v%qc = fails_error_max
                   iv%synop(n)%u%inv = 0.0
                   iv%synop(n)%v%inv = 0.0
                endif
             else
                if(iv%synop(n)%u%qc == fails_error_max .or. iv%synop(n)%v%qc == fails_error_max ) then
                   iv%synop(n)%u%qc = fails_error_max
                   iv%synop(n)%v%qc = fails_error_max
                   iv%synop(n)%u%inv = 0.0
                   iv%synop(n)%v%inv = 0.0
                endif
             endif
          endif


      failed=.false.
      if( iv%synop(n)%t%qc >= obs_qc_pointer )  then
      call da_max_error_qc (it, iv%info(synop), n, iv%synop(n)%t, max_error_t , failed)
      if( iv%info(synop)%proc_domain(1,n) ) then
      num_qcstat_conv(1,synop,3,1)= num_qcstat_conv(1,synop,3,1) + 1
      if(failed) then
      num_qcstat_conv(2,synop,3,1)= num_qcstat_conv(2,synop,3,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'synop',ob_vars(3),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%synop(n)%p%qc >= obs_qc_pointer )  then
      call da_max_error_qc (it, iv%info(synop), n, iv%synop(n)%p, max_error_p , failed)         
      if( iv%info(synop)%proc_domain(1,n) ) then
      num_qcstat_conv(1,synop,5,1)= num_qcstat_conv(1,synop,5,1) + 1
      if(failed) then
      num_qcstat_conv(2,synop,5,1)= num_qcstat_conv(2,synop,5,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'synop',ob_vars(5),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%synop(n)%q%qc >= obs_qc_pointer ) then
       if( iv%synop(n)%t%qc == fails_error_max .or. iv%synop(n)%p%qc == fails_error_max) then
       failed=.true.
       iv%synop(n)%q%qc  = fails_error_max
       iv%synop(n)%q%inv = 0.0
       else
       call da_max_error_qc (it, iv%info(synop), n, iv%synop(n)%q, max_error_q , failed)
       endif
      if( iv%info(synop)%proc_domain(1,n) ) then
      num_qcstat_conv(1,synop,4,1)= num_qcstat_conv(1,synop,4,1) + 1
      if(failed) then
      num_qcstat_conv(2,synop,4,1)= num_qcstat_conv(2,synop,4,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'synop',ob_vars(4),iv%info(synop)%lat(1,n),iv%info(synop)%lon(1,n),0.01*ob%synop(n)%p
      end if
      end if
      end if 
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_synop")

end subroutine da_check_max_iv_synop


subroutine da_calculate_grady_synop(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_synop")

   do n=1, iv%info(synop)%nlocal
      if (iv%synop(n)%u%qc < obs_qc_pointer) re%synop(n)%u = 0.0
      if (iv%synop(n)%v%qc < obs_qc_pointer) re%synop(n)%v = 0.0
      if (iv%synop(n)%t%qc < obs_qc_pointer) re%synop(n)%t = 0.0
      if (iv%synop(n)%p%qc < obs_qc_pointer) re%synop(n)%p = 0.0
      if (iv%synop(n)%q%qc < obs_qc_pointer) re%synop(n)%q = 0.0

      if (iv%synop(n)%u%qc < obs_qc_pointer) re%synop(n)%u = 0.0
      if (iv%synop(n)%v%qc < obs_qc_pointer) re%synop(n)%v = 0.0
      if (iv%synop(n)%t%qc < obs_qc_pointer) re%synop(n)%t = 0.0
      if (iv%synop(n)%p%qc < obs_qc_pointer) re%synop(n)%p = 0.0
      if (iv%synop(n)%q%qc < obs_qc_pointer) re%synop(n)%q = 0.0

      jo_grad_y%synop(n)%u = -re%synop(n)%u / (iv%synop(n)%u%error * iv%synop(n)%u%error)
      jo_grad_y%synop(n)%v = -re%synop(n)%v / (iv%synop(n)%v%error * iv%synop(n)%v%error)
      jo_grad_y%synop(n)%t = -re%synop(n)%t / (iv%synop(n)%t%error * iv%synop(n)%t%error)
      jo_grad_y%synop(n)%p = -re%synop(n)%p / (iv%synop(n)%p%error * iv%synop(n)%p%error)
      jo_grad_y%synop(n)%q = -re%synop(n)%q / (iv%synop(n)%q%error * iv%synop(n)%q%error)
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_synop")
     
end subroutine da_calculate_grady_synop


subroutine da_check_buddy_synop(iv, ob, dx, it)

   !-----------------------------------------------------------------------
   ! Purpose: Buddy check for SYNOP observations.
   !
   ! For SYNOP, there may not need the binning procedure before going
   ! into the da_buddy_qc. So  bottom_pressure = 30000.0 nad num_bins_p = 1.
   ! If you want to do binning, minor modifications needed.
   !
   !                       Yong-Run Guo, 10/10/2008
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   type(y_type),  intent(in)    :: ob      ! Observation structure
   integer,       intent(in)    :: it      ! Outer iteration
   real,          intent(in)    :: dx

   integer :: k, n, bin, i, j, m_max, kk, nn, numobs
   real    :: dx_km, Press_mb

! All data in one bin:
   integer, parameter               :: num_bins_p = 1
   real, parameter                  :: bottom_pressure = 30000.0
! 
! Total of 13 bins used:
!   integer, parameter               :: num_bins_p = 13
!   real, parameter                  :: bottom_pressure = 100000.0

   real, parameter                  :: bin_width_p = 10000.0
   real   , dimension(0:num_bins_p) :: bin_start_p, pressure, bin_width
   integer, dimension(0:num_bins_p) :: num
!
   integer, allocatable, dimension(:,:) :: n_iv

   integer,          allocatable, dimension(:) :: qc_flag_small
   real,             allocatable, dimension(:) :: xob, yob, obs
   character(len=5), allocatable, dimension(:) :: station_id
!-----------------------------------------------------------------------------
   
!   if (trace_use_dull) call da_trace_entry("da_check_buddy_synop")

   !--------------------------------------------------------------------------- 
   ! [1.0] Open diagnostic file:
   !---------------------------------------------------------------------------

   if (rootproc .and. check_buddy_print) then
      write (check_buddy_unit,'(/A)')  &
         '================================================================'
      write (unit = check_buddy_unit, fmt = '(A,i4,A,i4/)') &
            'SYNOP BUDDY TEST QC:  no_buddies_qc=',no_buddies,&
            '  fails_buddy_check_qc=',fails_buddy_check
   end if

   !---------------------------------------------------------------------------
   ! [2.0] Bin the data vertically based on the obs p::
   !---------------------------------------------------------------------------

!   print*,'==> Synop Buddy check: num_bins_p = ',num_bins_p
   dx_km = dx / 1000.0
!  
   bin_start_p(0) = bottom_pressure
   pressure   (0) = bin_start_p(0)
   bin_width      (0) = 0.0     
   do n = 1, num_bins_p
      bin_start_p(n) = bin_start_p(n-1) - bin_width(n-1)
      if (bin_start_p(n) > 30000.0) then
         bin_width(n) = bin_width_p
      else
         bin_width(n) = bin_width_p / 2.0
      endif
      pressure(n) = bin_start_p(n) - bin_width(n)/2.0
   enddo
   bin_start_p(0) = bottom_pressure + 10.0
!
! Only 1 bin=0 used, if you want to do the normal binning, comment out 
! the line below:
   pressure   (0) = 100000.0

!   print '(I3,2x,"start_p=",f10.1," mid-pressure=",f10.1," width=",f10.1)', &
!        (n, bin_start_p(n), pressure(n), bin_width(n), n=0, num_bins_p)
!
! 2.1 Get the maximum dimension for all the bins:
!
   num = 0
   do n = iv%info(synop)%n1,iv%info(synop)%n2
         if (ob%synop(n)%p > missing_r) then
           do i = 0, num_bins_p - 1
              if (iv%synop(n)%p%qc >=0 .and.             &
                  (ob%synop(n)%p <= bin_start_p(i) .and. &
                   ob%synop(n)%p >  bin_start_p(i+1)) ) then
                 bin = i
                 exit
              endif
           enddo
!           bin = int( (bottom_pressure - ob%synop(n)%p)/bin_width(n) ) + 1
           if (ob%synop(n)%p > bottom_pressure) bin = 0
           if (ob%synop(n)%p <=  bin_start_p(num_bins_p)) bin = num_bins_p
           num(bin) = num(bin) + 1
         endif
   enddo
   m_max = maxval(num)
!   print *,(i,num(i),i=0,num_bins_p)
!   print *,"m_max=", m_max
!
! 2.2 Save the location indices (n,k) for each of bins:
!
!   print '("Synop n1=",i5,"  n2=",i5)',iv%info(synop)%n1,iv%info(synop)%n2
   allocate ( n_iv( 0: num_bins_p,1:m_max+10 ) )

   num = 0
   do n = iv%info(synop)%n1,iv%info(synop)%n2
         if (ob%synop(n)%p > missing_r) then
           do i = 0, num_bins_p - 1
              if (iv%synop(n)%p%qc >=0 .and.             &
                  (ob%synop(n)%p <= bin_start_p(i) .and. &
                   ob%synop(n)%p >  bin_start_p(i+1)) ) then
                 bin = i
                 exit
              endif
           enddo
!           bin = int( (bottom_pressure - ob%synop(n)%p)/bin_width(n) ) + 1
           if (ob%synop(n)%p > bottom_pressure) bin = 0
           if (ob%synop(n)%p <=  bin_start_p(num_bins_p)) bin = num_bins_p

           num(bin) = num(bin) + 1
           n_iv(bin,num(bin)) = n

         endif
   end do
!
! 2.3 Print out the binned results:
!
!   do i = 0, num_bins_p
!      print '("bin:",I2,"  start_p=",f8.1," num=",i5)', &
!                      i, bin_start_p(i), num(i)
!      do j = 1, num(i)
!         n = n_iv(i,j)
!         print '("j, n:",2i5,2x,"p=",f10.1)', &
!                  j, n, ob%synop(n)%p
!      enddo
!   enddo
   !---------------------------------------------------------------------------
   ! [3.0] Buddy check for each of the pressure-bins::
   !---------------------------------------------------------------------------

   do i = 0, num_bins_p

     if (num(i) <= 1) cycle

! 3.1 Get the Station locations:
   
     ! Pressure level:
     Press_mb = pressure(i) / 100.0
     numobs = num(i)

     allocate(xob(1:numobs))
     allocate(yob(1:numobs))
     allocate(obs(1:numobs))
     allocate(qc_flag_small(1:numobs))
     allocate(station_id   (1:numobs))

     if (rootproc .and. check_buddy_print) then
         write (check_buddy_unit,'(5X,A)')  &
         '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write (unit = check_buddy_unit, fmt = '(5X,A,I3,2X,A,I6)') &
             'BIN =', i, 'NUMOBS =', numobs
     end if
!     print *,'SYNOP: BIN=', i, '  numobs=',num(i)

     ! Station locations

     do n = 1, numobs
        nn = n_iv(i,n)
!
        station_id(n)        = iv%info(synop)%id(nn)
        xob(n)               = iv%info(synop)%x(1,nn)
        yob(n)               = iv%info(synop)%y(1,nn)
     enddo
 
! 3.2 U-component buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'UU      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_uv=', max_buddy_uv 

     obs = 0.0
     qc_flag_small(n) = missing
     do n = 1, numobs
        nn = n_iv(i,n)
        obs(n)               = iv%synop(nn)%u%inv
        qc_flag_small(n)     = iv%synop(nn)%u%qc
     enddo

     call da_buddy_qc (numobs, m_max, station_id, xob, yob, obs, qc_flag_small,&
                       'UU      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_uv , check_buddy_unit, check_buddy_print )

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        iv%synop(nn)%u%qc = qc_flag_small(n)
     enddo

! 3.2 V-component buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'VV      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_uv=', max_buddy_uv
 

     obs = 0.0
     qc_flag_small(n) = missing
     do n = 1, numobs
        nn = n_iv(i,n)
        obs(n)               = iv%synop(nn)%v%inv
        qc_flag_small(n)     = iv%synop(nn)%v%qc
     enddo

     call da_buddy_qc (numobs, m_max, station_id, xob, yob, obs, qc_flag_small,&
                       'VV      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_uv , check_buddy_unit, check_buddy_print )

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        iv%synop(nn)%v%qc = qc_flag_small(n)
     enddo

! 3.3 Temperature buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'TT      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_t=', max_buddy_t 

     obs = 0.0
     qc_flag_small(n) = missing
     do n = 1, numobs
        nn = n_iv(i,n)
        obs(n)               = iv%synop(nn)%t%inv
        qc_flag_small(n)     = iv%synop(nn)%t%qc
     enddo

     call da_buddy_qc (numobs, m_max, station_id, xob, yob, obs, qc_flag_small,&
                       'TT      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_t , check_buddy_unit, check_buddy_print )

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        iv%synop(nn)%t%qc = qc_flag_small(n)
     enddo

! 3.3 Specific humidity buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'QQ      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_rh=', max_buddy_rh 

     obs = 0.0
     qc_flag_small(n) = missing
     do n = 1, numobs
        nn = n_iv(i,n)
        obs(n)               = iv%synop(nn)%q%inv
        qc_flag_small(n)     = iv%synop(nn)%q%qc
     enddo

     call da_buddy_qc (numobs, m_max, station_id, xob, yob, obs, qc_flag_small,&
                       'QQ      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_rh , check_buddy_unit, check_buddy_print )

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        iv%synop(nn)%q%qc = qc_flag_small(n)
     enddo

! 3.4 Pressure buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'PMSL    ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_p=', max_buddy_p 

     obs = 0.0
     qc_flag_small(n) = missing
     do n = 1, numobs
        nn = n_iv(i,n)
        obs(n)               = iv%synop(nn)%p%inv
        qc_flag_small(n)     = iv%synop(nn)%p%qc
     enddo

     call da_buddy_qc (numobs, m_max, station_id, xob, yob, obs, qc_flag_small,&
                       'PMSL    ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_p , check_buddy_unit, check_buddy_print )

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        iv%synop(nn)%p%qc = qc_flag_small(n)
     enddo

! 3.5 Deallocate arrays

1234 continue

     deallocate(xob)
     deallocate(yob)
     deallocate(obs)
     deallocate(qc_flag_small)
     deallocate(station_id   )

   enddo

   deallocate ( n_iv )

   if (trace_use_dull) call da_trace_exit("da_check_buddy_synop")

end subroutine da_check_buddy_synop

end module da_synop

