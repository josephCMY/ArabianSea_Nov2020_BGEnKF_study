












module da_mtgirs

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing_data, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, fails_error_max, &
      max_stheight_diff,test_dm_exact, anal_type_verify, &
      kms,kme,kts,kte, &
      trace_use_dull, mtgirs, position_lev_dependant, max_ext_its, qcstat_conv_unit, ob_vars, &
      convert_fd2uv, convert_uv2fd, max_error_spd, max_error_dir, max_omb_spd, max_omb_dir, pi, qc_rej_both,&
      wind_sd_mtgirs, wind_stats_sd 
   use da_grid_definitions, only : da_ffdduv,da_ffdduv_model, da_ffdduv_diagnose
   use da_physics, only : da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use module_domain, only : domain
   use da_interpolation, only : da_to_zk, da_interp_lin_3d, &
      da_interp_lin_3d_adj, da_interp_lin_2d, da_interp_lin_2d_adj, da_interp_lin_2d_partial
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk, da_get_print_lvl
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_mtgirs1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: q                        
   end type residual_mtgirs1_type

   type maxmin_mtgirs_stats_type
      type (maxmin_type)         :: u, v, t, q
   end type maxmin_mtgirs_stats_type

   type stats_mtgirs_type
      type (maxmin_mtgirs_stats_type)  :: maximum, minimum
      type (residual_mtgirs1_type)     :: average, rms_err
   end type stats_mtgirs_type

contains

subroutine da_ao_stats_mtgirs (stats_unit, iv, re, ob)

   !-------------------------------------------------------------
   ! Purpose: TBD   
   !-------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type (y_type),  intent (in)    :: re            ! A - O
   type(y_type),   intent (in)    :: ob            ! Observation structure.

   type (stats_mtgirs_type) :: stats
   integer                 :: nu, nv, nt, nq
   integer                 :: n, k
   real                    :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_mtgirs")

   nu = 0
   nv = 0
   nt = 0
   nq = 0

   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%maximum%t = maxmin_type (missing_r, 0, 0)
   stats%maximum%q = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_mtgirs1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(mtgirs)%nlocal
      do k=1, iv%info(mtgirs)%levels(n)
         if (iv%info(mtgirs)%proc_domain(1,n)) then

            u_inc = re%mtgirs(n)%u(k)
            v_inc = re%mtgirs(n)%v(k)
            u_obs = ob%mtgirs(n)%u(k)
            v_obs = ob%mtgirs(n)%v(k)

            if (.not. wind_sd_mtgirs .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%mtgirs(n)%u(k)%qc, iv%mtgirs(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_mtgirs .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%mtgirs(n)%u(k)%qc, iv%mtgirs(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate (n, k, iv%mtgirs(n)%u(k)%qc, & 
               u_inc, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%mtgirs(n)%v(k)%qc, & 
               v_inc, nu, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
            call da_stats_calculate (n, k, iv%mtgirs(n)%t(k)%qc, & 
               re%mtgirs(n)%t(k), nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate (n, k, iv%mtgirs(n)%q(k)%qc, & 
               re%mtgirs(n)%q(k), nq, &
               stats%minimum%q, stats%maximum%q, &
               stats%average%q, stats%rms_err%q)
         end if
      end do
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (nq)
   iv%nstats(mtgirs) = nu + nv + nt + nq

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
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
      stats%minimum%q%value, stats%maximum%q%value, &
      stats%minimum%q%n, stats%maximum%q%n, &
      stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for mtgirs'
         call da_print_stats_mtgirs(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_mtgirs")

end subroutine da_ao_stats_mtgirs


subroutine da_jo_and_grady_mtgirs(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv         ! Innovation vector.
   type (y_type),  intent(in)    :: re         ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y  ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo         ! Obs cost function.

   integer                      :: n, k
   ! the following "global" objects are used only when testing
   type (iv_type) :: iv_glob         ! Global Innovation vector (O-B).
   type (y_type)  :: re_glob         ! Global Residual vector (O-A).
   type (y_type)  :: jo_grad_y_glob  ! Global Grad_y(Jo)
   
   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_mtgirs")

   jo % mtgirs_u = 0.0
   jo % mtgirs_v = 0.0
   jo % mtgirs_t = 0.0
   jo % mtgirs_q = 0.0

   if (test_dm_exact) then
      if (iv%info(mtgirs)%ntotal == 0) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_mtgirs")
         return
      end if
   else
      if (iv%info(mtgirs)%nlocal < 1) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_mtgirs")
         return
      end if
   end if

   do n=1, iv%info(mtgirs)%nlocal
       do k=1, iv%info(mtgirs)%levels(n)
          jo_grad_y%mtgirs(n)%u(k) = -re%mtgirs(n)%u(k) / (iv%mtgirs(n)%u(k)%error * iv%mtgirs(n)%u(k)%error)
          jo_grad_y%mtgirs(n)%v(k) = -re%mtgirs(n)%v(k) / (iv%mtgirs(n)%v(k)%error * iv%mtgirs(n)%v(k)%error)
          jo_grad_y%mtgirs(n)%t(k) = -re%mtgirs(n)%t(k) / (iv%mtgirs(n)%t(k)%error * iv%mtgirs(n)%t(k)%error)
          jo_grad_y%mtgirs(n)%q(k) = -re%mtgirs(n)%q(k) / (iv%mtgirs(n)%q(k)%error * iv%mtgirs(n)%q(k)%error)
      end do
   end do

   ! Bitwise-exact reduction preserves operation order of serial code for 
   ! testing, at the cost of much-increased run-time.  Turn it off when not 
   ! testing.  This will always be .false. for a serial or 1-MPI-process run.  
   if (test_dm_exact) then
      ! collect all obs in serial order and allocate global objects
!      call da_to_global_mtgirs(iv, re, jo_grad_y, iv_glob, re_glob, jo_grad_y_glob)
      ! perform remaining computations
      call da_jo_mtgirs_uvtq(iv_glob, re_glob, jo_grad_y_glob, jo)
      ! free global objects
!      call da_deallocate_global_mtgirs(iv_glob, re_glob, jo_grad_y_glob)
   else
      ! perform remaining computations
      call da_jo_mtgirs_uvtq(iv, re, jo_grad_y, jo)
   end if

   jo % mtgirs_u = 0.5 * jo % mtgirs_u
   jo % mtgirs_v = 0.5 * jo % mtgirs_v
   jo % mtgirs_t = 0.5 * jo % mtgirs_t
   jo % mtgirs_q = 0.5 * jo % mtgirs_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_mtgirs")

end subroutine da_jo_and_grady_mtgirs


subroutine da_jo_mtgirs_uvtq (iv, re, jo_grad_y, jo)

   !-----------------------------------------------------------------------
   ! Purpose: Ensures that exactly the same code is used for both 
   ! serial and parallel computations when testing_dm_bitwise_exact is set.
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv         ! Innovation vector.
   type (y_type),  intent(in)    :: re         ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y  ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo         ! Obs cost function.
 
   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_mtgirs_uvtq")

   do n=1, iv%info(mtgirs)%nlocal
      do k=1, iv%info(mtgirs)%levels(n)
         if (iv%info(mtgirs)%proc_domain(1,n)) then
            jo % mtgirs_u = jo % mtgirs_u - re%mtgirs(n)%u(k) * jo_grad_y%mtgirs(n)%u(k)
            jo % mtgirs_v = jo % mtgirs_v - re%mtgirs(n)%v(k) * jo_grad_y%mtgirs(n)%v(k)
            jo % mtgirs_t = jo % mtgirs_t - re%mtgirs(n)%t(k) * jo_grad_y%mtgirs(n)%t(k)
            jo % mtgirs_q = jo % mtgirs_q - re%mtgirs(n)%q(k) * jo_grad_y%mtgirs(n)%q(k)
         end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_jo_mtgirs_uvtq")

end subroutine da_jo_mtgirs_uvtq


subroutine da_residual_mtgirs(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for mtgirs obs
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual structure.

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type)              :: n_obs_bad
   integer                           :: n, k

   if (trace_use_dull) call da_trace_entry("da_residual_mtgirs")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(mtgirs)%nlocal
      do k=1, iv%info(mtgirs)%levels(n)
         np_available = np_available + 4

         re%mtgirs(n)%u(k) = da_residual(n, k, y%mtgirs(n)%u(k), iv%mtgirs(n)%u(k), n_obs_bad % u)
         re%mtgirs(n)%v(k) = da_residual(n, k, y%mtgirs(n)%v(k), iv%mtgirs(n)%v(k), n_obs_bad % v)
         re%mtgirs(n)%t(k) = da_residual(n, k, y%mtgirs(n)%t(k), iv%mtgirs(n)%t(k), n_obs_bad % t)
         re%mtgirs(n)%q(k) = da_residual(n, k, y%mtgirs(n)%q(k), iv%mtgirs(n)%q(k), n_obs_bad % q)
      end do
   end do

   np_missing = np_missing + n_obs_bad % u % num % miss + &
      n_obs_bad % v % num % miss + n_obs_bad % t % num % miss + &
      n_obs_bad % q % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad + &
      n_obs_bad % v % num % bad + n_obs_bad % t % num % bad + &
      n_obs_bad % q % num % bad
   np_obs_used = np_obs_used + n_obs_bad % u % num % use + &
      n_obs_bad % v % num % use + n_obs_bad % t % num % use + &
      n_obs_bad % q % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_mtgirs")

end subroutine da_residual_mtgirs


subroutine da_oi_stats_mtgirs (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! OI
   type(y_type),   intent (in)      :: ob            ! Observation structure.

   type (stats_mtgirs_type)         :: stats
   integer                          :: nu, nv, nt, nq
   integer                          :: n, k
   real                             :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_mtgirs")

   nu = 0
   nv = 0
   nt = 0
   nq = 0

   stats%maximum%u = maxmin_type(missing_r, 0, 0)
   stats%maximum%v = maxmin_type(missing_r, 0, 0)
   stats%maximum%t = maxmin_type(missing_r, 0, 0)
   stats%maximum%q = maxmin_type(missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_mtgirs1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(mtgirs)%nlocal
      do k=1, iv%info(mtgirs)%levels(n)
         if (iv%info(mtgirs)%proc_domain(1,n)) then

            u_inv = iv%mtgirs(n)%u(k)%inv
            v_inv = iv%mtgirs(n)%v(k)%inv
            u_obs = ob%mtgirs(n)%u(k)
            v_obs = ob%mtgirs(n)%v(k)

            if (.not. wind_sd_mtgirs .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%mtgirs(n)%u(k)%qc, iv%mtgirs(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_mtgirs .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%mtgirs(n)%u(k)%qc, iv%mtgirs(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate(iv%info(mtgirs)%obs_global_index(n), &
               k, iv%mtgirs(n)%u(k)%qc, u_inv, nu, &
               stats%minimum%u, stats%maximum%u, stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(mtgirs)%obs_global_index(n), &
               k, iv%mtgirs(n)%v(k)%qc, v_inv, nv, &
               stats%minimum%v, stats%maximum%v, stats%average%v, stats%rms_err%v)
            call da_stats_calculate(iv%info(mtgirs)%obs_global_index(n), &
               k, iv%mtgirs(n)%t(k)%qc, iv%mtgirs(n)%t(k)%inv, nt, &
               stats%minimum%t, stats%maximum%t, stats%average%t, stats%rms_err%t)
            call da_stats_calculate(iv%info(mtgirs)%obs_global_index(n), &
               k, iv%mtgirs(n)%q(k)%qc, iv%mtgirs(n)%q(k)%inv, nq, &
               stats%minimum%q, stats%maximum%q, stats%average%q, stats%rms_err%q)
         end if
      end do
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nu)
   call da_proc_sum_int(nv)
   call da_proc_sum_int(nt)
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
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
      stats%minimum%q%value, stats%maximum%q%value, &
      stats%minimum%q%n, stats%maximum%q%n, &
      stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for mtgirs'
         call da_print_stats_mtgirs(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_mtgirs")

end subroutine da_oi_stats_mtgirs


subroutine da_print_stats_mtgirs(stats_unit, nu, nv, nt, nq, mtgirs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nu, nv, nt, nq
   type (stats_mtgirs_type), intent(in)    :: mtgirs

   if (trace_use_dull) call da_trace_entry("da_print_stats_mtgirs")

   write(unit=stats_unit, fmt='(5a/)') &
      '   var             ', &
      'u (m/s)     n    k    ', &
      'v (m/s)     n    k    ', &
      't (K)       n    k    ', &
      'q (kg/kg)   n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
      '  Number: ', nu, nv, nt, nq

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1
   if (nt < 1) nt = 1
   if (nq < 1) nq = 1

   write(unit=stats_unit, fmt='((a,3(f12.4,2i5),e12.4,2i5))') &
      ' Minimum(n,k): ', mtgirs%minimum%u, mtgirs%minimum%v, &
                         mtgirs%minimum%t, mtgirs%minimum%q, &
      ' Maximum(n,k): ', mtgirs%maximum%u, mtgirs%maximum%v, &
                         mtgirs%maximum%t, mtgirs%maximum%q
   write(unit=stats_unit, fmt='((a,3(f12.4,10x),e12.4,10x))') &
      ' Average     : ', mtgirs%average%u/real(nu), &
                         mtgirs%average%v/real(nv), &
                         mtgirs%average%t/real(nt), &
                         mtgirs%average%q/real(nq), &
      '    RMSE     : ', sqrt(mtgirs%rms_err%u/real(nu)), &
                         sqrt(mtgirs%rms_err%v/real(nv)), &
                         sqrt(mtgirs%rms_err%t/real(nt)), &
                         sqrt(mtgirs%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_mtgirs")

end subroutine da_print_stats_mtgirs


subroutine da_transform_xtoy_mtgirs (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),     intent(in)    :: grid
   type (iv_type),    intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: y        ! y = h (grid%xa) (linear)

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   integer :: n,k

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_mtgirs")

   allocate (u(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (v(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (t(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (q(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
  
   allocate (ub(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (vb(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))

   call da_interp_lin_3d (grid%xa%u, iv%info(mtgirs), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(mtgirs), v)
   call da_interp_lin_3d (grid%xa%t, iv%info(mtgirs), t)
   call da_interp_lin_3d (grid%xa%q, iv%info(mtgirs), q)

   call da_interp_lin_3d (grid%xb%u, iv%info(mtgirs), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(mtgirs), vb)

   do n=iv%info(mtgirs)%n1,iv%info(mtgirs)%n2
      do k = 1, iv%info(mtgirs)%levels(n)
         if(wind_sd_mtgirs) then
            call da_uv_to_sd_lin(y%mtgirs(n)%u(k),y%mtgirs(n)%v(k),u(k,n),v(k,n),ub(k,n),vb(k,n))
         else
            y%mtgirs(n)%u(k) = u(k,n)
            y%mtgirs(n)%v(k) = v(k,n)
         end if
      end do
      y%mtgirs(n)%t(:) = t(1:size(y%mtgirs(n)%t),n)
      y%mtgirs(n)%q(:) = q(1:size(y%mtgirs(n)%q),n)
   end do

   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_mtgirs")

end subroutine da_transform_xtoy_mtgirs


subroutine da_transform_xtoy_mtgirs_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none
   type(domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n,k

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_mtgirs_adj")

   allocate (u(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (v(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (t(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (q(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))

   allocate (ub(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (vb(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   call da_interp_lin_3d (grid%xb%u, iv%info(mtgirs), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(mtgirs), vb)

   do n=iv%info(mtgirs)%n1,iv%info(mtgirs)%n2
      do k = 1, iv%info(mtgirs)%levels(n)
         if(wind_sd_mtgirs) then
            call da_uv_to_sd_adj(jo_grad_y%mtgirs(n)%u(k), &
                                 jo_grad_y%mtgirs(n)%v(k), u(k,n), v(k,n), ub(k,n), vb(k,n))
         else
            u(k,n) = jo_grad_y%mtgirs(n)%u(k)
            v(k,n) = jo_grad_y%mtgirs(n)%v(k)
         end if
      end do
      t(1:size(jo_grad_y%mtgirs(n)%t),n) = jo_grad_y%mtgirs(n)%t(:)
      q(1:size(jo_grad_y%mtgirs(n)%q),n) = jo_grad_y%mtgirs(n)%q(:)
   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(mtgirs), u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(mtgirs), v)
   call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(mtgirs), t)
   call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(mtgirs), q)
   call da_interp_lin_3d (grid%xb%u, iv%info(mtgirs), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(mtgirs), vb)

   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)
   deallocate (ub)
   deallocate (vb)
   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_mtgirs_adj")

end subroutine da_transform_xtoy_mtgirs_adj


subroutine da_check_max_iv_mtgirs(iv, it, num_qcstat_conv)

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

   integer :: k,n, ipr
   logical :: failed,failed1,failed2

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_mtgirs")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(mtgirs)%n1,iv%info(mtgirs)%n2
    do k = 1, iv%info(mtgirs)%levels(n)
         call da_get_print_lvl(iv%mtgirs(n)%p(k),ipr)
         if(.not. qc_rej_both)then
             if(wind_sd_mtgirs)then
               failed=.false.
               if( iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%u(k), max_error_spd,failed)
                   if( iv%info(mtgirs)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,mtgirs,1,ipr) = num_qcstat_conv(1,mtgirs,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,mtgirs,1,ipr) = num_qcstat_conv(2,mtgirs,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'mtgirs',ob_vars(1),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
                       end if
                   end if
                end if

                failed=.false.
                if( iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%v(k), max_error_dir,failed)
                    if( iv%info(mtgirs)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,mtgirs,2,ipr) = num_qcstat_conv(1,mtgirs,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,mtgirs,2,ipr) = num_qcstat_conv(2,mtgirs,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'mtgirs',ob_vars(2),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
                        end if
                    end if
                end if
             else
                failed=.false.
                if( iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%u(k), max_error_uv,failed)
                    if( iv%info(mtgirs)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,mtgirs,1,ipr) = num_qcstat_conv(1,mtgirs,1,ipr) + 1
                        if(failed) then
                           num_qcstat_conv(2,mtgirs,1,ipr) = num_qcstat_conv(2,mtgirs,1,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'mtgirs',ob_vars(1),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
                        end if
                    end if
                end if

                failed=.false.
                if( iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%v(k), max_error_uv,failed)
                    if( iv%info(mtgirs)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,mtgirs,2,ipr) = num_qcstat_conv(1,mtgirs,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,mtgirs,2,ipr) = num_qcstat_conv(2,mtgirs,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'mtgirs',ob_vars(2),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
                        end if
                    end if
                 end if
             end if

             if(wind_sd_mtgirs)then
                if(iv%mtgirs(n)%u(k)%qc == fails_error_max .or. abs(iv%mtgirs(n)%u(k)%inv) >= max_omb_spd) then
                   iv%mtgirs(n)%u(k)%qc = fails_error_max
                   iv%mtgirs(n)%u(k)%inv = 0.0
                endif
                if(iv%mtgirs(n)%v(k)%qc == fails_error_max .or. abs(iv%mtgirs(n)%v(k)%inv) >= max_omb_dir) then
                   iv%mtgirs(n)%v(k)%qc = fails_error_max
                   iv%mtgirs(n)%v(k)%inv = 0.0
                endif
             endif
          else
             failed1=.false.
             failed2=.false.

             if( iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer .or. iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer )  then
                 if(wind_sd_mtgirs)then
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%u(k), max_error_spd,failed1)
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%v(k), max_error_dir,failed2)
                 else
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%u(k), max_error_uv,failed1)
                    call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%v(k), max_error_uv,failed2)
                 endif
             endif

             if( iv%info(mtgirs)%proc_domain(k,n) ) then
                 num_qcstat_conv(1,mtgirs,1,ipr) = num_qcstat_conv(1,mtgirs,1,ipr) + 1
                 num_qcstat_conv(1,mtgirs,2,ipr) = num_qcstat_conv(1,mtgirs,2,ipr) + 1

                 if(failed1 .or. failed2) then
                    num_qcstat_conv(2,mtgirs,1,ipr) = num_qcstat_conv(2,mtgirs,1,ipr) + 1
                    write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'mtgirs',ob_vars(1),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
                    num_qcstat_conv(2,mtgirs,2,ipr) = num_qcstat_conv(2,mtgirs,2,ipr) + 1
                    write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'mtgirs',ob_vars(2),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
                 end if
             end if

	     if(wind_sd_mtgirs)then
                if(iv%mtgirs(n)%u(k)%qc == fails_error_max .or. iv%mtgirs(n)%v(k)%qc == fails_error_max .or. &
                   abs(iv%mtgirs(n)%v(k)%inv) >= max_omb_dir .or. abs(iv%mtgirs(n)%u(k)%inv) >= max_omb_spd )then
                   iv%mtgirs(n)%u(k)%qc = fails_error_max
                   iv%mtgirs(n)%v(k)%qc = fails_error_max
                   iv%mtgirs(n)%u(k)%inv = 0.0
                   iv%mtgirs(n)%v(k)%inv = 0.0
                endif
             else
                if(iv%mtgirs(n)%u(k)%qc == fails_error_max .or. iv%mtgirs(n)%v(k)%qc == fails_error_max ) then
                   iv%mtgirs(n)%u(k)%qc = fails_error_max
                   iv%mtgirs(n)%v(k)%qc = fails_error_max
                   iv%mtgirs(n)%u(k)%inv = 0.0
                   iv%mtgirs(n)%v(k)%inv = 0.0
                endif
             endif
         endif


         failed=.false.
         if( iv%mtgirs(n)%t(k)%qc >= obs_qc_pointer )  then
         call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%t(k), max_error_t ,failed)
         if( iv%info(mtgirs)%proc_domain(k,n) ) then
             num_qcstat_conv(1,mtgirs,3,ipr) = num_qcstat_conv(1,mtgirs,3,ipr) + 1
         if(failed) then
          num_qcstat_conv(2,mtgirs,3,ipr) = num_qcstat_conv(2,mtgirs,3,ipr) + 1
           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'mtgirs',ob_vars(3),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
         end if
         end if
         end if

         failed=.false.
         if( iv%mtgirs(n)%q(k)%qc >= obs_qc_pointer )  then 
          if( iv%mtgirs(n)%t(k)%qc == fails_error_max ) then
          failed=.true.
          iv%mtgirs(n)%q(k)%qc  = fails_error_max 
          iv%mtgirs(n)%q(k)%inv = 0.0
          else
          call da_max_error_qc (it,iv%info(mtgirs), n, iv%mtgirs(n)%q(k), max_error_q ,failed)
          endif
         if( iv%info(mtgirs)%proc_domain(k,n) ) then
            num_qcstat_conv(1,mtgirs,4,ipr) = num_qcstat_conv(1,mtgirs,4,ipr) + 1
         if(failed) then
         num_qcstat_conv(2,mtgirs,4,ipr) = num_qcstat_conv(2,mtgirs,4,ipr) + 1
           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'mtgirs',ob_vars(4),iv%info(mtgirs)%lat(k,n),iv%info(mtgirs)%lon(k,n),0.01*iv%mtgirs(n)%p(k)
         end if
         end if
         end if

    end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_mtgirs")

end subroutine da_check_max_iv_mtgirs
subroutine da_get_innov_vector_mtgirs (it, num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD  
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it       ! External iteration.
   type(domain),     intent(in)    :: grid     ! first guess state.
   type(y_type),     intent(inout) :: ob       ! Observation structure.
   type(iv_type),    intent(inout) :: iv       ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n, k        ! Loop counter.
   integer :: i  (kms:kme)
   integer :: j  (kms:kme)
   real    :: dx (kms:kme)
   real    :: dxm(kms:kme)  
   real    :: dy (kms:kme)
   real    :: dym(kms:kme)  
   real    :: speed, direction

   real, allocatable :: model_u(:,:)  ! Model value u at ob location.
   real, allocatable :: model_v(:,:)  ! Model value v at ob location.
   real, allocatable :: model_t(:,:)  ! Model value t at ob location.
   real, allocatable :: model_q(:,:)  ! Model value q at ob location.

   real    :: v_h(kts:kte)      ! Model value h at ob hor. location.
   real    :: v_p(kts:kte)      ! Model value p at ob hor. location.


   if (trace_use_dull) call da_trace_entry ("da_get_innov_vector_mtgirs")

   allocate (model_u(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (model_v(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (model_t(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))
   allocate (model_q(iv%info(mtgirs)%max_lev,iv%info(mtgirs)%n1:iv%info(mtgirs)%n2))

   model_u(:,:) = 0.0
   model_v(:,:) = 0.0
   model_t(:,:) = 0.0
   model_q(:,:) = 0.0

   if ( it > 1 ) then
      do n=iv%info(mtgirs)%n1, iv%info(mtgirs)%n2
         do k=1, iv%info(mtgirs)%levels(n)
            if (iv%mtgirs(n)%u(k)%qc == fails_error_max) iv%mtgirs(n)%u(k)%qc = 0
            if (iv%mtgirs(n)%v(k)%qc == fails_error_max) iv%mtgirs(n)%v(k)%qc = 0
            if (iv%mtgirs(n)%t(k)%qc == fails_error_max) iv%mtgirs(n)%t(k)%qc = 0
            if (iv%mtgirs(n)%q(k)%qc == fails_error_max) iv%mtgirs(n)%q(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(mtgirs)%n1, iv%info(mtgirs)%n2
      if (iv%info(mtgirs)%levels(n) < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      if (position_lev_dependant) then
         i(:)   = iv%info(mtgirs)%i(:,n)
         j(:)   = iv%info(mtgirs)%j(:,n)
         dx(:)  = iv%info(mtgirs)%dx(:,n)
         dy(:)  = iv%info(mtgirs)%dy(:,n)
         dxm(:) = iv%info(mtgirs)%dxm(:,n)
         dym(:) = iv%info(mtgirs)%dym(:,n)
         do k=kts,kte
            v_h(k) = dym(k)*(dxm(k)*grid%xb%h(i(k),j(k),k) + dx(k)*grid%xb%h(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%h(i(k),j(k)+1,k) + dx(k)*grid%xb%h(i(k)+1,j(k)+1,k))
            v_p(k) = dym(k)*(dxm(k)*grid%xb%p(i(k),j(k),k) + dx(k)*grid%xb%p(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%p(i(k),j(k)+1,k) + dx(k)*grid%xb%p(i(k)+1,j(k)+1,k))
         end do
      else
         i(1)   = iv%info(mtgirs)%i(1,n)
         j(1)   = iv%info(mtgirs)%j(1,n)
         dx(1)  = iv%info(mtgirs)%dx(1,n)
         dy(1)  = iv%info(mtgirs)%dy(1,n)
         dxm(1) = iv%info(mtgirs)%dxm(1,n)
         dym(1) = iv%info(mtgirs)%dym(1,n)

         v_h(kts:kte) = dym(1) * (dxm(1)*grid%xb%h(i(1),j(1),kts:kte)   + dx(1)*grid%xb%h(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%h(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%h(i(1)+1,j(1)+1,kts:kte))
         v_p(kts:kte) = dym(1) * (dxm(1)*grid%xb%p(i(1),j(1),kts:kte)   + dx(1)*grid%xb%p(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%p(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%p(i(1)+1,j(1)+1,kts:kte))
      end if

      do k=1, iv%info(mtgirs)%levels(n)
         if (iv%mtgirs(n)%p(k) > 1.0) then
            call da_to_zk (iv%mtgirs(n)%p(k), v_p, v_interp_p, iv%info(mtgirs)%zk(k,n))
         else if (iv%mtgirs(n)%h(k) > 0.0) then
            call da_to_zk (iv%mtgirs(n)%h(k), v_h, v_interp_h, iv%info(mtgirs)%zk(k,n))
         end if
      end do

   end do

   call da_convert_zk (iv%info(mtgirs))

   if (.not. anal_type_verify) then
      do n=iv%info(mtgirs)%n1,iv%info(mtgirs)%n2
         do k=1, iv%info(mtgirs)%levels(n)
            if (iv%info(mtgirs)%zk(k,n) < 0.0) then
               iv%mtgirs(n)%u(k)%qc = missing_data
               iv%mtgirs(n)%v(k)%qc = missing_data
               iv%mtgirs(n)%t(k)%qc = missing_data
               iv%mtgirs(n)%q(k)%qc = missing_data
            end if
         end do
      end do
   end if

   ! [1.2] Interpolate horizontally to ob:

   call da_interp_lin_3d (grid%xb%u, iv%info(mtgirs), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(mtgirs), model_v)
   call da_interp_lin_3d (grid%xb%t, iv%info(mtgirs), model_t)
   call da_interp_lin_3d (grid%xb%q, iv%info(mtgirs), model_q)

   do n=iv%info(mtgirs)%n1, iv%info(mtgirs)%n2
      !----------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !----------------------------------------------------------------------

      do k = 1, iv%info(mtgirs)%levels(n)
         iv%mtgirs(n)%u(k)%inv = 0.0
         iv%mtgirs(n)%v(k)%inv = 0.0
         iv%mtgirs(n)%t(k)%inv = 0.0
         iv%mtgirs(n)%q(k)%inv = 0.0

         !-------------------------------------------------------------------
         ! [3.0] Interpolation:
         !-------------------------------------------------------------------
         if (wind_sd_mtgirs) then
             call da_ffdduv_model (speed,direction,model_u(k,n), model_v(k,n), convert_uv2fd)

              if (ob%mtgirs(n)%u(k) > missing_r .AND. iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer) then
                  iv%mtgirs(n)%u(k)%inv = ob%mtgirs(n)%u(k) - speed
              end if

              if (ob%mtgirs(n)%v(k) > missing_r .AND. iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer) then
                  iv%mtgirs(n)%v(k)%inv = ob%mtgirs(n)%v(k) - direction
                  if (iv%mtgirs(n)%v(k)%inv > 180.0 ) iv%mtgirs(n)%v(k)%inv = iv%mtgirs(n)%v(k)%inv - 360.0
                  if (iv%mtgirs(n)%v(k)%inv < -180.0 ) iv%mtgirs(n)%v(k)%inv = iv%mtgirs (n)%v(k)%inv + 360.0
              end if
         else
              if (ob%mtgirs(n)%u(k) > missing_r .AND. iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer) then
                  iv%mtgirs(n)%u(k)%inv = ob%mtgirs(n)%u(k) - model_u(k,n)
              end if

              if (ob%mtgirs(n)%v(k) > missing_r .AND. iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer) then
                  iv%mtgirs(n)%v(k)%inv = ob%mtgirs(n)%v(k) - model_v(k,n)
              end if
          end if


         if (ob%mtgirs(n)%t(k) > missing_r .AND. iv%mtgirs(n)%t(k)%qc >= obs_qc_pointer) then
            iv%mtgirs(n)%t(k)%inv = ob%mtgirs(n)%t(k) - model_t(k,n)
         end if

         if (ob%mtgirs(n)%q(k) > missing_r .AND. iv%mtgirs(n)%q(k)%qc >= obs_qc_pointer) then
            iv%mtgirs(n)%q(k)%inv = ob%mtgirs(n)%q(k) - model_q(k,n)
         end if
      end do
   end do

   !----------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !----------------------------------------------------------------------

    if ( check_max_iv ) &
       call da_check_max_iv_mtgirs (iv, it, num_qcstat_conv)

   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_q)
   
   if (trace_use_dull) call da_trace_exit ("da_get_innov_vector_mtgirs")

end subroutine da_get_innov_vector_mtgirs


subroutine da_calculate_grady_mtgirs(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_mtgirs")

   do n=1, iv%info(mtgirs)%nlocal
      do k=1, iv%info(mtgirs)%levels(n)
         if (iv%mtgirs(n)%u(k)%qc < obs_qc_pointer) re%mtgirs(n)%u(k) = 0.0
         if (iv%mtgirs(n)%v(k)%qc < obs_qc_pointer) re%mtgirs(n)%v(k) = 0.0
         if (iv%mtgirs(n)%t(k)%qc < obs_qc_pointer) re%mtgirs(n)%t(k) = 0.0
         if (iv%mtgirs(n)%q(k)%qc < obs_qc_pointer) re%mtgirs(n)%q(k) = 0.0

         jo_grad_y%mtgirs(n)%u(k) = -re%mtgirs(n)%u(k) / (iv%mtgirs(n)%u(k)%error * iv%mtgirs(n)%u(k)%error)
         jo_grad_y%mtgirs(n)%v(k) = -re%mtgirs(n)%v(k) / (iv%mtgirs(n)%v(k)%error * iv%mtgirs(n)%v(k)%error)
         jo_grad_y%mtgirs(n)%t(k) = -re%mtgirs(n)%t(k) / (iv%mtgirs(n)%t(k)%error * iv%mtgirs(n)%t(k)%error)
         jo_grad_y%mtgirs(n)%q(k) = -re%mtgirs(n)%q(k) / (iv%mtgirs(n)%q(k)%error * iv%mtgirs(n)%q(k)%error)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_mtgirs")

end subroutine da_calculate_grady_mtgirs



end module da_mtgirs

