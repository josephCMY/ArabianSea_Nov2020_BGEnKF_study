












module da_tamdar

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing_data, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, sfc_assi_options, &



      max_stheight_diff,test_dm_exact, anal_type_verify, &
      kms,kme,kts,kte,sfc_assi_options_1,sfc_assi_options_2, num_procs, comm, &
      trace_use_dull, tamdar, tamdar_sfc, position_lev_dependant, max_ext_its, &
      qcstat_conv_unit,ob_vars, fails_error_max, &
      convert_fd2uv,convert_uv2fd,max_error_spd,max_error_dir,max_omb_spd,max_omb_dir,pi,qc_rej_both, &
      wind_sd_tamdar, wind_stats_sd
   use da_grid_definitions, only : da_ffdduv,da_ffdduv_model, da_ffdduv_diagnose
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use module_domain, only : domain
   use da_interpolation, only : da_to_zk, da_interp_lin_3d, &
      da_interp_lin_3d_adj, da_interp_lin_2d, da_interp_lin_2d_adj, da_interp_lin_2d_partial
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_obs_sfc_correction,da_convert_zk, &
                        da_buddy_qc, da_get_print_lvl

   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int

   use da_physics, only : da_sfc_pre, da_transform_xtopsfc, &
      da_transform_xtopsfc_adj, da_uv_to_sd_lin, da_uv_to_sd_adj

   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_tamdar1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: q                        
   end type residual_tamdar1_type

   type maxmin_tamdar_stats_type
      type (maxmin_type)         :: u, v, t, q
   end type maxmin_tamdar_stats_type

   type stats_tamdar_type
      type (maxmin_tamdar_stats_type)  :: maximum, minimum
      type (residual_tamdar1_type)     :: average, rms_err
   end type stats_tamdar_type

   

   type residual_tamdar_sfc1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: p                        
      real          :: q                        
   end type residual_tamdar_sfc1_type

   type maxmin_tamdar_sfc_stats_type
      type (maxmin_type)         :: u, v, t, p, q
   end type maxmin_tamdar_sfc_stats_type

   type stats_tamdar_sfc_type
      type (maxmin_tamdar_sfc_stats_type)  :: maximum, minimum
      type (residual_tamdar_sfc1_type)     :: average, rms_err
   end type stats_tamdar_sfc_type

contains

subroutine da_ao_stats_tamdar (stats_unit, iv, re, ob)

   !-------------------------------------------------------------
   ! Purpose: TBD   
   !-------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type (y_type),  intent (in)    :: re            ! A - O
   type(y_type),   intent (in)    :: ob            ! Observation structure.

   type (stats_tamdar_type) :: stats
   integer                 :: nu, nv, nt, nq
   integer                 :: n, k
   real                    :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_tamdar")

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

   stats%average = residual_tamdar1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(tamdar)%nlocal
      do k=1, iv%info(tamdar)%levels(n)
         if (iv%info(tamdar)%proc_domain(1,n)) then

            u_inc = re%tamdar(n)%u(k)
            v_inc = re%tamdar(n)%v(k)
            u_obs = ob%tamdar(n)%u(k)
            v_obs = ob%tamdar(n)%v(k)

            if (.not. wind_sd_tamdar .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%tamdar(n)%u(k)%qc, iv%tamdar(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_tamdar .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%tamdar(n)%u(k)%qc, iv%tamdar(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate (n, k, iv%tamdar(n)%u(k)%qc, & 
               u_inc, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%tamdar(n)%v(k)%qc, & 
               v_inc, nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
            call da_stats_calculate (n, k, iv%tamdar(n)%t(k)%qc, & 
               re%tamdar(n)%t(k), nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate (n, k, iv%tamdar(n)%q(k)%qc, & 
               re%tamdar(n)%q(k), nq, &
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
   iv%nstats(tamdar) = nu + nv + nt + nq

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for tamdar'
         call da_print_stats_tamdar(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_tamdar")

end subroutine da_ao_stats_tamdar


subroutine da_jo_and_grady_tamdar(iv, re, jo, jo_grad_y)

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
   
   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_tamdar")

   jo % tamdar_u = 0.0
   jo % tamdar_v = 0.0
   jo % tamdar_t = 0.0
   jo % tamdar_q = 0.0

   if (test_dm_exact) then
      if (iv%info(tamdar)%ntotal == 0) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_tamdar")
         return
      end if
   else
      if (iv%info(tamdar)%nlocal < 1) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_tamdar")
         return
      end if
   end if

   do n=1, iv%info(tamdar)%nlocal
       do k=1, iv%info(tamdar)%levels(n)
          jo_grad_y%tamdar(n)%u(k) = -re%tamdar(n)%u(k) / (iv%tamdar(n)%u(k)%error * iv%tamdar(n)%u(k)%error)
          jo_grad_y%tamdar(n)%v(k) = -re%tamdar(n)%v(k) / (iv%tamdar(n)%v(k)%error * iv%tamdar(n)%v(k)%error)
          jo_grad_y%tamdar(n)%t(k) = -re%tamdar(n)%t(k) / (iv%tamdar(n)%t(k)%error * iv%tamdar(n)%t(k)%error)
          jo_grad_y%tamdar(n)%q(k) = -re%tamdar(n)%q(k) / (iv%tamdar(n)%q(k)%error * iv%tamdar(n)%q(k)%error)
      end do
   end do

   ! Bitwise-exact reduction preserves operation order of serial code for 
   ! testing, at the cost of much-increased run-time.  Turn it off when not 
   ! testing.  This will always be .false. for a serial or 1-MPI-process run.  
   if (test_dm_exact) then
      ! perform remaining computations
      call da_jo_tamdar_uvtq(iv_glob, re_glob, jo_grad_y_glob, jo)
   else
      ! perform remaining computations
      call da_jo_tamdar_uvtq(iv, re, jo_grad_y, jo)
   end if

   jo % tamdar_u = 0.5 * jo % tamdar_u
   jo % tamdar_v = 0.5 * jo % tamdar_v
   jo % tamdar_t = 0.5 * jo % tamdar_t
   jo % tamdar_q = 0.5 * jo % tamdar_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_tamdar")

end subroutine da_jo_and_grady_tamdar


subroutine da_jo_tamdar_uvtq (iv, re, jo_grad_y, jo)

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

   if (trace_use_dull) call da_trace_entry("da_jo_tamdar_uvtq")

   do n=1, iv%info(tamdar)%nlocal
      do k=1, iv%info(tamdar)%levels(n)
         if (iv%info(tamdar)%proc_domain(1,n)) then
            jo % tamdar_u = jo % tamdar_u - re%tamdar(n)%u(k) * jo_grad_y%tamdar(n)%u(k)
            jo % tamdar_v = jo % tamdar_v - re%tamdar(n)%v(k) * jo_grad_y%tamdar(n)%v(k)
            jo % tamdar_t = jo % tamdar_t - re%tamdar(n)%t(k) * jo_grad_y%tamdar(n)%t(k)
            jo % tamdar_q = jo % tamdar_q - re%tamdar(n)%q(k) * jo_grad_y%tamdar(n)%q(k)
         end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_jo_tamdar_uvtq")

end subroutine da_jo_tamdar_uvtq


subroutine da_residual_tamdar(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for tamdar obs
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

   if (trace_use_dull) call da_trace_entry("da_residual_tamdar")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(tamdar)%nlocal
      do k=1, iv%info(tamdar)%levels(n)
         np_available = np_available + 4

         re%tamdar(n)%u(k) = da_residual(n, k, y%tamdar(n)%u(k), iv%tamdar(n)%u(k), n_obs_bad % u)
         re%tamdar(n)%v(k) = da_residual(n, k, y%tamdar(n)%v(k), iv%tamdar(n)%v(k), n_obs_bad % v)
         re%tamdar(n)%t(k) = da_residual(n, k, y%tamdar(n)%t(k), iv%tamdar(n)%t(k), n_obs_bad % t)
         re%tamdar(n)%q(k) = da_residual(n, k, y%tamdar(n)%q(k), iv%tamdar(n)%q(k), n_obs_bad % q)
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

   if (trace_use_dull) call da_trace_exit("da_residual_tamdar")

end subroutine da_residual_tamdar


subroutine da_oi_stats_tamdar (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! OI
   type(y_type),   intent (in)      :: ob            ! Observation structure.

   type (stats_tamdar_type)         :: stats
   integer                          :: nu, nv, nt, nq
   integer                          :: n, k
   real                             :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_tamdar")

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
   stats%average = residual_tamdar1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(tamdar)%nlocal
      do k=1, iv%info(tamdar)%levels(n)
!         if (iv%info(tamdar)%proc_domain(k,n)) then
         if (iv%info(tamdar)%proc_domain(1,n)) then

            u_inv = iv%tamdar(n)%u(k)%inv
            v_inv = iv%tamdar(n)%v(k)%inv
            u_obs = ob%tamdar(n)%u(k)
            v_obs = ob%tamdar(n)%v(k)

            if (.not. wind_sd_tamdar .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%tamdar(n)%u(k)%qc,iv%tamdar(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_tamdar .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%tamdar(n)%u(k)%qc,iv%tamdar(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate(iv%info(tamdar)%obs_global_index(n), &
               k, iv%tamdar(n)%u(k)%qc, u_inv, nu, &
               stats%minimum%u, stats%maximum%u, stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(tamdar)%obs_global_index(n), &
               k, iv%tamdar(n)%v(k)%qc, v_inv, nv, &
               stats%minimum%v, stats%maximum%v, stats%average%v, stats%rms_err%v)
            call da_stats_calculate(iv%info(tamdar)%obs_global_index(n), &
               k, iv%tamdar(n)%t(k)%qc, iv%tamdar(n)%t(k)%inv, nt, &
               stats%minimum%t, stats%maximum%t, stats%average%t, stats%rms_err%t)
            call da_stats_calculate(iv%info(tamdar)%obs_global_index(n), &
               k, iv%tamdar(n)%q(k)%qc, iv%tamdar(n)%q(k)%inv, nq, &
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for tamdar'
         call da_print_stats_tamdar(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_tamdar")

end subroutine da_oi_stats_tamdar


subroutine da_print_stats_tamdar(stats_unit, nu, nv, nt, nq, tamdar)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nu, nv, nt, nq
   type (stats_tamdar_type), intent(in)    :: tamdar

   if (trace_use_dull) call da_trace_entry("da_print_stats_tamdar")

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
      ' Minimum(n,k): ', tamdar%minimum%u, tamdar%minimum%v, &
                         tamdar%minimum%t, tamdar%minimum%q, &
      ' Maximum(n,k): ', tamdar%maximum%u, tamdar%maximum%v, &
                         tamdar%maximum%t, tamdar%maximum%q
   write(unit=stats_unit, fmt='((a,3(f12.4,10x),e12.4,10x))') &
      ' Average     : ', tamdar%average%u/real(nu), &
                         tamdar%average%v/real(nv), &
                         tamdar%average%t/real(nt), &
                         tamdar%average%q/real(nq), &
      '    RMSE     : ', sqrt(tamdar%rms_err%u/real(nu)), &
                         sqrt(tamdar%rms_err%v/real(nv)), &
                         sqrt(tamdar%rms_err%t/real(nt)), &
                         sqrt(tamdar%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_tamdar")

end subroutine da_print_stats_tamdar


subroutine da_transform_xtoy_tamdar (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
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

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_tamdar")

   allocate (u(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (v(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (t(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (q(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))

   allocate (ub(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (vb(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))

   call da_interp_lin_3d (grid%xa%u, iv%info(tamdar), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(tamdar), v)
   call da_interp_lin_3d (grid%xa%t, iv%info(tamdar), t)
   call da_interp_lin_3d (grid%xa%q, iv%info(tamdar), q)

   call da_interp_lin_3d (grid%xb%u, iv%info(tamdar), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(tamdar), vb)

   do n=iv%info(tamdar)%n1,iv%info(tamdar)%n2
      do k = 1, iv%info(tamdar)%levels(n)
         if(wind_sd_tamdar) then
            call da_uv_to_sd_lin(y%tamdar(n)%u(k),y%tamdar(n)%v(k),u(k,n),v(k,n),ub(k,n),vb(k,n))
         else
            y%tamdar(n)%u(k) = u(k,n)
            y%tamdar(n)%v(k) = v(k,n)
         end if
      end do
      y%tamdar(n)%t(:) = t(1:size(y%tamdar(n)%t),n)
      y%tamdar(n)%q(:) = q(1:size(y%tamdar(n)%q),n)
   end do

   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_tamdar")

end subroutine da_transform_xtoy_tamdar


subroutine da_transform_xtoy_tamdar_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n, k

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)

   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)


   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_tamdar_adj")

   allocate (u(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (v(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (t(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (q(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))

   allocate (ub(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (vb(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))

   call da_interp_lin_3d (grid%xb%u, iv%info(tamdar), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(tamdar), vb)

   do n=iv%info(tamdar)%n1,iv%info(tamdar)%n2
      do k = 1, iv%info(tamdar)%levels(n)
         if(wind_sd_tamdar) then
            call da_uv_to_sd_adj(jo_grad_y%tamdar(n)%u(k), &
                                 jo_grad_y%tamdar(n)%v(k), u(k,n), v(k,n), ub(k,n), vb(k,n))
         else
            u(k,n) = jo_grad_y%tamdar(n)%u(k)
            v(k,n) = jo_grad_y%tamdar(n)%v(k)
         end if
      end do
      t(1:size(jo_grad_y%tamdar(n)%t),n) = jo_grad_y%tamdar(n)%t(:)
      q(1:size(jo_grad_y%tamdar(n)%q),n) = jo_grad_y%tamdar(n)%q(:)
   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(tamdar), u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(tamdar), v)
   call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(tamdar), t)
   call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(tamdar), q)

   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_tamdar_adj")

end subroutine da_transform_xtoy_tamdar_adj


subroutine da_check_max_iv_tamdar(iv, it,num_qcstat_conv)

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

   integer :: k,n,ipr
   logical :: failed,failed1,failed2

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_tamdar")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(tamdar)%n1,iv%info(tamdar)%n2
      do k = 1, iv%info(tamdar)%levels(n)
         call da_get_print_lvl(iv%tamdar(n)%p(k),ipr) 

         if(.not. qc_rej_both)then
            if(wind_sd_tamdar)then
               failed=.false.
               if( iv%tamdar(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%u(k), max_error_spd,failed)
                   if( iv%info(tamdar)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,tamdar,1,ipr) = num_qcstat_conv(1,tamdar,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,tamdar,1,ipr) = num_qcstat_conv(2,tamdar,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'tamdar',ob_vars(1),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
                       end if
                   end if
                end if

                failed=.false.
                if( iv%tamdar(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%v(k), max_error_dir,failed)
                    if( iv%info(tamdar)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,tamdar,2,ipr) = num_qcstat_conv(1,tamdar,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,tamdar,2,ipr) = num_qcstat_conv(2,tamdar,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'tamdar',ob_vars(2),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
                        end if
                    end if
                end if
             else
                failed=.false.
                if( iv%tamdar(n)%u(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%u(k), max_error_uv,failed)
                    if( iv%info(tamdar)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,tamdar,1,ipr) = num_qcstat_conv(1,tamdar,1,ipr) + 1
                        if(failed) then
                           num_qcstat_conv(2,tamdar,1,ipr) = num_qcstat_conv(2,tamdar,1,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'tamdar',ob_vars(1),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
                         end if
                    end if
                 end if

                 failed=.false.
                 if( iv%tamdar(n)%v(k)%qc >= obs_qc_pointer ) then
                     call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%v(k), max_error_uv,failed)
                     if( iv%info(tamdar)%proc_domain(k,n) ) then
                         num_qcstat_conv(1,tamdar,2,ipr) = num_qcstat_conv(1,tamdar,2,ipr) + 1
                         if(failed)then
                            num_qcstat_conv(2,tamdar,2,ipr) = num_qcstat_conv(2,tamdar,2,ipr) + 1
                            write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                            'tamdar',ob_vars(2),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
                         end if
                     end if
                 end if
              end if

              if(wind_sd_tamdar)then
                 if(iv%tamdar(n)%u(k)%qc == fails_error_max .or. abs(iv%tamdar(n)%u(k)%inv) >= max_omb_spd) then
                    iv%tamdar(n)%u(k)%qc = fails_error_max
                    iv%tamdar(n)%u(k)%inv = 0.0
                 endif
                 if(iv%tamdar(n)%v(k)%qc == fails_error_max .or. abs(iv%tamdar(n)%v(k)%inv) >= max_omb_dir) then
                    iv%tamdar(n)%v(k)%qc = fails_error_max
                    iv%tamdar(n)%v(k)%inv = 0.0
                 endif
              endif
           else
              failed1=.false.
              failed2=.false.

              if( iv%tamdar(n)%v(k)%qc >= obs_qc_pointer .or. iv%tamdar(n)%u(k)%qc >= obs_qc_pointer )  then
                  if(wind_sd_tamdar)then
                     call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%u(k), max_error_spd,failed1)
                     call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%v(k), max_error_dir,failed2)
                  else
                     call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%u(k), max_error_uv,failed1)
                     call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%v(k), max_error_uv,failed2)
                  endif
              endif

              if( iv%info(tamdar)%proc_domain(k,n) ) then
                  num_qcstat_conv(1,tamdar,1,ipr) = num_qcstat_conv(1,tamdar,1,ipr) + 1
                  num_qcstat_conv(1,tamdar,2,ipr) = num_qcstat_conv(1,tamdar,2,ipr) + 1

                  if(failed1 .or. failed2) then
                     num_qcstat_conv(2,tamdar,1,ipr) = num_qcstat_conv(2,tamdar,1,ipr) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'tamdar',ob_vars(1),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
                     num_qcstat_conv(2,tamdar,2,ipr) = num_qcstat_conv(2,tamdar,2,ipr) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'tamdar',ob_vars(2),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
                  endif
               endif

               if(wind_sd_tamdar)then
                  if(iv%tamdar(n)%u(k)%qc == fails_error_max .or. iv%tamdar(n)%v(k)%qc == fails_error_max .or. &
                     abs(iv%tamdar(n)%v(k)%inv) >= max_omb_dir .or. abs(iv%tamdar(n)%u(k)%inv) >= max_omb_spd )then
                     iv%tamdar(n)%u(k)%qc = fails_error_max
                     iv%tamdar(n)%v(k)%qc = fails_error_max
                     iv%tamdar(n)%u(k)%inv = 0.0
                     iv%tamdar(n)%v(k)%inv = 0.0
                  endif
               else
                  if(iv%tamdar(n)%u(k)%qc == fails_error_max .or. iv%tamdar(n)%v(k)%qc == fails_error_max ) then
                     iv%tamdar(n)%u(k)%qc = fails_error_max
                     iv%tamdar(n)%v(k)%qc = fails_error_max
                     iv%tamdar(n)%u(k)%inv = 0.0
                     iv%tamdar(n)%v(k)%inv = 0.0
                  endif
               endif
            endif

         failed=.false.
         if( iv%tamdar(n)%t(k)%qc >= obs_qc_pointer )  then
         call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%t(k), max_error_t ,failed)
         if( iv%info(tamdar)%proc_domain(k,n) ) then
             num_qcstat_conv(1,tamdar,3,ipr) = num_qcstat_conv(1,tamdar,3,ipr) + 1
         if(failed) then
          num_qcstat_conv(2,tamdar,3,ipr) = num_qcstat_conv(2,tamdar,3,ipr) + 1
           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'tamdar',ob_vars(3),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
         end if
         end if
         end if

         failed=.false.
         if( iv%tamdar(n)%q(k)%qc >= obs_qc_pointer ) then
          if( iv%tamdar(n)%t(k)%qc == fails_error_max ) then
          failed=.true.
          iv%tamdar(n)%q(k)%qc  = fails_error_max
          iv%tamdar(n)%q(k)%inv = 0.0
          else
          call da_max_error_qc (it,iv%info(tamdar), n, iv%tamdar(n)%q(k), max_error_q ,failed)
          endif
         if( iv%info(tamdar)%proc_domain(k,n) ) then
                    num_qcstat_conv(1,tamdar,4,ipr) = num_qcstat_conv(1,tamdar,4,ipr) + 1
         if(failed) then
         num_qcstat_conv(2,tamdar,4,ipr) = num_qcstat_conv(2,tamdar,4,ipr) + 1
           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'tamdar',ob_vars(4),iv%info(tamdar)%lat(k,n),iv%info(tamdar)%lon(k,n),0.01*iv%tamdar(n)%p(k)
         end if
         end if
         end if

      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_tamdar")

end subroutine da_check_max_iv_tamdar
subroutine da_get_innov_vector_tamdar (it, num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD  
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

   real, allocatable :: model_u(:,:)  ! Model value u at ob location.
   real, allocatable :: model_v(:,:)  ! Model value v at ob location.
   real, allocatable :: model_t(:,:)  ! Model value t at ob location.
   real, allocatable :: model_q(:,:)  ! Model value q at ob location.

   real    :: speed, direction

   real    :: v_h(kts:kte)      ! Model value h at ob hor. location.
   real    :: v_p(kts:kte)      ! Model value p at ob hor. location.

   if (trace_use_dull) call da_trace_entry ("da_get_innov_vector_tamdar")

   allocate (model_u(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (model_v(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (model_t(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
   allocate (model_q(iv%info(tamdar)%max_lev,iv%info(tamdar)%n1:iv%info(tamdar)%n2))

   model_u(:,:) = 0.0
   model_v(:,:) = 0.0
   model_t(:,:) = 0.0
   model_q(:,:) = 0.0

   if ( it > 1 ) then
      do n=iv%info(tamdar)%n1,iv%info(tamdar)%n2
         do k=1, iv%info(tamdar)%levels(n)
            if (iv%tamdar(n)%u(k)%qc == fails_error_max) iv%tamdar(n)%u(k)%qc = 0
            if (iv%tamdar(n)%v(k)%qc == fails_error_max) iv%tamdar(n)%v(k)%qc = 0
            if (iv%tamdar(n)%t(k)%qc == fails_error_max) iv%tamdar(n)%t(k)%qc = 0
            if (iv%tamdar(n)%q(k)%qc == fails_error_max) iv%tamdar(n)%q(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(tamdar)%n1, iv%info(tamdar)%n2
      if (iv%info(tamdar)%levels(n) < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      if (position_lev_dependant) then
         i(:)   = iv%info(tamdar)%i(:,n)
         j(:)   = iv%info(tamdar)%j(:,n)
         dx(:)  = iv%info(tamdar)%dx(:,n)
         dy(:)  = iv%info(tamdar)%dy(:,n)
         dxm(:) = iv%info(tamdar)%dxm(:,n)
         dym(:) = iv%info(tamdar)%dym(:,n)
         do k=kts,kte
            v_h(k) = dym(k)*(dxm(k)*grid%xb%h(i(k),j(k),k) + dx(k)*grid%xb%h(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%h(i(k),j(k)+1,k) + dx(k)*grid%xb%h(i(k)+1,j(k)+1,k))
            v_p(k) = dym(k)*(dxm(k)*grid%xb%p(i(k),j(k),k) + dx(k)*grid%xb%p(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%p(i(k),j(k)+1,k) + dx(k)*grid%xb%p(i(k)+1,j(k)+1,k))
         end do
      else
         i(1)   = iv%info(tamdar)%i(1,n)
         j(1)   = iv%info(tamdar)%j(1,n)
         dx(1)  = iv%info(tamdar)%dx(1,n)
         dy(1)  = iv%info(tamdar)%dy(1,n)
         dxm(1) = iv%info(tamdar)%dxm(1,n)
         dym(1) = iv%info(tamdar)%dym(1,n)

         v_h(kts:kte) = dym(1) * (dxm(1)*grid%xb%h(i(1),j(1),kts:kte)   + dx(1)*grid%xb%h(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%h(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%h(i(1)+1,j(1)+1,kts:kte))
         v_p(kts:kte) = dym(1) * (dxm(1)*grid%xb%p(i(1),j(1),kts:kte)   + dx(1)*grid%xb%p(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%p(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%p(i(1)+1,j(1)+1,kts:kte))
      end if

      do k=1, iv%info(tamdar)%levels(n)
         if (iv%tamdar(n)%p(k) > 1.0) then
            call da_to_zk (iv%tamdar(n)%p(k), v_p, v_interp_p, iv%info(tamdar)%zk(k,n))
         else if (iv%tamdar(n)%h(k) > 0.0) then
            call da_to_zk (iv%tamdar(n)%h(k), v_h, v_interp_h, iv%info(tamdar)%zk(k,n))
         end if
      end do

   end do

   call da_convert_zk (iv%info(tamdar))

   if (.not. anal_type_verify) then
      do n=iv%info(tamdar)%n1,iv%info(tamdar)%n2
         do k=1, iv%info(tamdar)%levels(n)
            if (iv%info(tamdar)%zk(k,n) < 0.0) then
               iv%tamdar(n)%u(k)%qc = missing_data
               iv%tamdar(n)%v(k)%qc = missing_data
               iv%tamdar(n)%t(k)%qc = missing_data
               iv%tamdar(n)%q(k)%qc = missing_data
            end if
         end do
      end do
   end if

   ! [1.2] Interpolate horizontally to ob:
   call da_interp_lin_3d (grid%xb%u, iv%info(tamdar), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(tamdar), model_v)
   call da_interp_lin_3d (grid%xb%t, iv%info(tamdar), model_t)
   call da_interp_lin_3d (grid%xb%q, iv%info(tamdar), model_q)

   do n=iv%info(tamdar)%n1, iv%info(tamdar)%n2
      !----------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !----------------------------------------------------------------------

      do k = 1, iv%info(tamdar)%levels(n)
         iv%tamdar(n)%u(k)%inv = 0.0
         iv%tamdar(n)%v(k)%inv = 0.0
         iv%tamdar(n)%t(k)%inv = 0.0
         iv%tamdar(n)%q(k)%inv = 0.0

         !-------------------------------------------------------------------
         ! [3.0] Interpolation:
         !-----------------------------------------------------------------
          if (wind_sd_tamdar) then
              call da_ffdduv_model (speed,direction,model_u(k,n), model_v(k,n), convert_uv2fd)

              if (ob%tamdar(n)%u(k) > missing_r .AND. iv%tamdar(n)%u(k)%qc >= obs_qc_pointer) then
                  iv%tamdar(n)%u(k)%inv = ob%tamdar(n)%u(k) - speed
              end if

              if (ob%tamdar(n)%v(k) > missing_r .AND. iv%tamdar(n)%v(k)%qc >= obs_qc_pointer) then
                  iv%tamdar(n)%v(k)%inv = ob%tamdar(n)%v(k) - direction
                  if (iv%tamdar(n)%v(k)%inv > 180.0 ) iv%tamdar(n)%v(k)%inv = iv%tamdar(n)%v(k)%inv - 360.0
                  if (iv%tamdar(n)%v(k)%inv < -180.0 ) iv%tamdar(n)%v(k)%inv = iv%tamdar(n)%v(k)%inv + 360.0
              end if
          else
             if (ob%tamdar(n)%u(k) > missing_r .AND. iv%tamdar(n)%u(k)%qc >= obs_qc_pointer) then
                 iv%tamdar(n)%u(k)%inv = ob%tamdar(n)%u(k) - model_u(k,n)
             end if

             if (ob%tamdar(n)%v(k) > missing_r .AND. iv%tamdar(n)%v(k)%qc >= obs_qc_pointer) then
                 iv%tamdar(n)%v(k)%inv = ob%tamdar(n)%v(k) - model_v(k,n)
             end if
          end if


             if (ob%tamdar(n)%t(k) > missing_r .AND. iv%tamdar(n)%t(k)%qc >= obs_qc_pointer) then
                 iv%tamdar(n)%t(k)%inv = ob%tamdar(n)%t(k) - model_t(k,n)
             end if

             if (ob%tamdar(n)%q(k) > missing_r .AND. iv%tamdar(n)%q(k)%qc >= obs_qc_pointer) then
                 iv%tamdar(n)%q(k)%inv = ob%tamdar(n)%q(k) - model_q(k,n)
             end if
      end do
   end do

   !----------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !----------------------------------------------------------------------

   if (check_max_iv) call da_check_max_iv_tamdar (iv, it, num_qcstat_conv)

   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_q)
   if (trace_use_dull) call da_trace_exit ("da_get_innov_vector_tamdar")

end subroutine da_get_innov_vector_tamdar


subroutine da_calculate_grady_tamdar(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_tamdar")

   do n=1, iv%info(tamdar)%nlocal
      do k=1, iv%info(tamdar)%levels(n)
         if (iv%tamdar(n)%u(k)%qc < obs_qc_pointer) re%tamdar(n)%u(k) = 0.0
         if (iv%tamdar(n)%v(k)%qc < obs_qc_pointer) re%tamdar(n)%v(k) = 0.0
         if (iv%tamdar(n)%t(k)%qc < obs_qc_pointer) re%tamdar(n)%t(k) = 0.0
         if (iv%tamdar(n)%q(k)%qc < obs_qc_pointer) re%tamdar(n)%q(k) = 0.0

         jo_grad_y%tamdar(n)%u(k) = -re%tamdar(n)%u(k) / (iv%tamdar(n)%u(k)%error * iv%tamdar(n)%u(k)%error)
         jo_grad_y%tamdar(n)%v(k) = -re%tamdar(n)%v(k) / (iv%tamdar(n)%v(k)%error * iv%tamdar(n)%v(k)%error)
         jo_grad_y%tamdar(n)%t(k) = -re%tamdar(n)%t(k) / (iv%tamdar(n)%t(k)%error * iv%tamdar(n)%t(k)%error)
         jo_grad_y%tamdar(n)%q(k) = -re%tamdar(n)%q(k) / (iv%tamdar(n)%q(k)%error * iv%tamdar(n)%q(k)%error)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_tamdar")

end subroutine da_calculate_grady_tamdar




subroutine da_ao_stats_tamdar_sfc (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type  (y_type), intent (in)    :: re            ! A - O
   type(y_type),   intent (in)    :: ob            ! Observation structure.

   type (stats_tamdar_sfc_type)     :: stats
   integer                          :: nu, nv, nt, np, nq
   integer                          :: n
   real                             :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_tamdar_sfc")

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
   stats%average = residual_tamdar_sfc1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(tamdar_sfc)%nlocal
      if (iv%info(tamdar_sfc)%proc_domain(1,n)) then

         u_inc = re%tamdar_sfc(n)%u
         v_inc = re%tamdar_sfc(n)%v
         u_obs = ob%tamdar_sfc(n)%u
         v_obs = ob%tamdar_sfc(n)%v

         if (.not. wind_sd_tamdar .and. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%tamdar_sfc(n)%u%qc, iv%tamdar_sfc(n)%v%qc, convert_uv2fd)
         if (wind_sd_tamdar .and. .not. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%tamdar_sfc(n)%u%qc, iv%tamdar_sfc(n)%v%qc, convert_fd2uv)

         call da_stats_calculate (n, 0, iv%tamdar_sfc(n)%u%qc,  & 
            u_inc, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%tamdar_sfc(n)%v%qc,  & 
            v_inc, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate (n, 0, iv%tamdar_sfc(n)%t%qc,  & 
            re%tamdar_sfc(n)%t, nt, & 
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)  
         call da_stats_calculate (n, 0, iv%tamdar_sfc(n)%p%qc,  & 
            re%tamdar_sfc(n)%p, np, & 
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate (n, 0, iv%tamdar_sfc(n)%q%qc,  & 
            re%tamdar_sfc(n)%q, nq, & 
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
   iv%nstats(tamdar_sfc) =  nu + nv + nt + np + nq

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for tamdar_sfc'
         call da_print_stats_tamdar_sfc(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_tamdar_sfc")

end subroutine da_ao_stats_tamdar_sfc


subroutine da_jo_and_grady_tamdar_sfc( iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv          ! Innovation vector.
   type(y_type),  intent(in)    :: re          ! Residual vector.
   type(y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type(jo_type), intent(inout) :: jo          ! Obs cost function.

   integer                      :: n
   ! the following "global" objects are used only when testing
   type(iv_type) :: iv_glob         ! Global Innovation vector(O-B).
   type(y_type)  :: re_glob         ! Global Residual vector(O-A).
   type(y_type)  :: jo_grad_y_glob  ! Global Grad_y(Jo)

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_tamdar_sfc")

   jo % tamdar_sfc_u = 0.0
   jo % tamdar_sfc_v = 0.0
   jo % tamdar_sfc_t = 0.0
   jo % tamdar_sfc_p = 0.0
   jo % tamdar_sfc_q = 0.0

   if (test_dm_exact) then
      if (iv%info(tamdar_sfc)%ntotal == 0) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_tamdar_sfc")
         return
      end if
   else
      if (iv%info(tamdar_sfc)%nlocal < 1) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_tamdar_sfc")
         return
      end if
   end if

   do n=1, iv%info(tamdar_sfc)%nlocal
      jo_grad_y%tamdar_sfc(n)%u = -re%tamdar_sfc(n)%u / (iv%tamdar_sfc(n)%u%error * iv%tamdar_sfc(n)%u%error)
      jo_grad_y%tamdar_sfc(n)%v = -re%tamdar_sfc(n)%v / (iv%tamdar_sfc(n)%v%error * iv%tamdar_sfc(n)%v%error)
      jo_grad_y%tamdar_sfc(n)%t = -re%tamdar_sfc(n)%t / (iv%tamdar_sfc(n)%t%error * iv%tamdar_sfc(n)%t%error)
      jo_grad_y%tamdar_sfc(n)%p = -re%tamdar_sfc(n)%p / (iv%tamdar_sfc(n)%p%error * iv%tamdar_sfc(n)%p%error)
      jo_grad_y%tamdar_sfc(n)%q = -re%tamdar_sfc(n)%q / (iv%tamdar_sfc(n)%q%error * iv%tamdar_sfc(n)%q%error)
   end do

      call da_jo_tamdar_sfc_uvtq( iv, re, jo_grad_y, jo)

   jo % tamdar_sfc_u = 0.5 * jo % tamdar_sfc_u
   jo % tamdar_sfc_v = 0.5 * jo % tamdar_sfc_v
   jo % tamdar_sfc_t = 0.5 * jo % tamdar_sfc_t
   jo % tamdar_sfc_p = 0.5 * jo % tamdar_sfc_p
   jo % tamdar_sfc_q = 0.5 * jo % tamdar_sfc_q
   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_tamdar_sfc")

end subroutine da_jo_and_grady_tamdar_sfc


subroutine da_jo_tamdar_sfc_uvtq(iv, re, jo_grad_y, jo)

   !-----------------------------------------------------------------------
   ! Purpose: Ensures that exactly the same code is used for both
   ! serial and parallel computations when testing_dm_bitwise_exact is set.
   !-----------------------------------------------------------------------
 
   implicit none

   type (iv_type), intent(in  ) :: iv         ! Innovation vector.
   type (y_type),  intent(in  ) :: re         ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y  ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo         ! Obs cost function.
  
   integer :: n

   if (trace_use_dull) call da_trace_entry("da_jo_tamdar_sfc_uvtq")

   do n=1, iv%info(tamdar_sfc)%nlocal
      if (iv%info(tamdar_sfc)%proc_domain(1,n)) then
        jo % tamdar_sfc_u = jo % tamdar_sfc_u - re%tamdar_sfc(n)%u * jo_grad_y%tamdar_sfc(n)%u
        jo % tamdar_sfc_v = jo % tamdar_sfc_v - re%tamdar_sfc(n)%v * jo_grad_y%tamdar_sfc(n)%v
        jo % tamdar_sfc_t = jo % tamdar_sfc_t - re%tamdar_sfc(n)%t * jo_grad_y%tamdar_sfc(n)%t
        jo % tamdar_sfc_p = jo % tamdar_sfc_p - re%tamdar_sfc(n)%p * jo_grad_y%tamdar_sfc(n)%p
        jo % tamdar_sfc_q = jo % tamdar_sfc_q - re%tamdar_sfc(n)%q * jo_grad_y%tamdar_sfc(n)%q
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_jo_tamdar_sfc_uvtq")

end subroutine da_jo_tamdar_sfc_uvtq


subroutine da_residual_tamdar_sfc(iv, y, re,np_missing, np_bad_data,np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for tamdar surface obs
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual vector (O-A).

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type)  :: n_obs_bad
   integer               :: n

   if (trace_use_dull) call da_trace_entry("da_residual_tamdar_sfc")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % p % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(tamdar_sfc)%nlocal
      np_available = np_available + 5

      re%tamdar_sfc(n)%u = da_residual(n, 0, y%tamdar_sfc(n)%u, iv%tamdar_sfc(n)%u, n_obs_bad % u)
      re%tamdar_sfc(n)%v = da_residual(n, 0, y%tamdar_sfc(n)%v, iv%tamdar_sfc(n)%v, n_obs_bad % v)
      re%tamdar_sfc(n)%t = da_residual(n, 0, y%tamdar_sfc(n)%t, iv%tamdar_sfc(n)%t, n_obs_bad % t)
      re%tamdar_sfc(n)%p = da_residual(n, 0, y%tamdar_sfc(n)%p, iv%tamdar_sfc(n)%p, n_obs_bad % p)
      re%tamdar_sfc(n)%q = da_residual(n, 0, y%tamdar_sfc(n)%q, iv%tamdar_sfc(n)%q, n_obs_bad % q)
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

   if (trace_use_dull) call da_trace_exit("da_residual_tamdar_sfc")

end subroutine da_residual_tamdar_sfc


subroutine da_oi_stats_tamdar_sfc (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_tamdar_sfc_type):: stats
   integer                     :: nu, nv, nt, np, nq
   integer                     :: n
   real                        :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_tamdar_sfc")

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

   stats%average = residual_tamdar_sfc1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(tamdar_sfc)%nlocal
      if (iv%info(tamdar_sfc)%proc_domain(1,n)) then

        u_inv = iv%tamdar_sfc(n)%u%inv
        v_inv = iv%tamdar_sfc(n)%v%inv
        u_obs = ob%tamdar_sfc(n)%u
        v_obs = ob%tamdar_sfc(n)%v

        if (.not. wind_sd_tamdar .and. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%tamdar_sfc(n)%u%qc, iv%tamdar_sfc(n)%v%qc, convert_uv2fd)
        if (wind_sd_tamdar.and. .not. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%tamdar_sfc(n)%u%qc, iv%tamdar_sfc(n)%v%qc, convert_fd2uv)

         call da_stats_calculate(iv%info(tamdar_sfc)%obs_global_index(n), &
            0, iv%tamdar_sfc(n)%u%qc, &
            u_inv, nu, &
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate(iv%info(tamdar_sfc)%obs_global_index(n), &
            0, iv%tamdar_sfc(n)%v%qc, &
            v_inv, nv, &
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate(iv%info(tamdar_sfc)%obs_global_index(n), &
            0, iv%tamdar_sfc(n)%t%qc, &
            iv%tamdar_sfc(n)%t%inv, nt, &
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate(iv%info(tamdar_sfc)%obs_global_index(n), &
            0, iv%tamdar_sfc(n)%p%qc, &
            iv%tamdar_sfc(n)%p%inv, np, &
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)

         call da_stats_calculate(iv%info(tamdar_sfc)%obs_global_index(n), &
            0, iv%tamdar_sfc(n)%q%qc, &
            iv%tamdar_sfc(n)%q%inv, nq, &
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)
      end if    ! end if (iv%info(tamdar_sfc)%proc_domain(1,n))
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for tamdar_sfc'
         call da_print_stats_tamdar_sfc(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_tamdar_sfc")

 end subroutine da_oi_stats_tamdar_sfc


subroutine da_print_stats_tamdar_sfc(stats_unit, nu, nv, nt, np, nq, tamdar_sfc)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                     intent(in)    :: stats_unit
   integer,                     intent(inout) :: nu, nv, nt, np, nq
   type (stats_tamdar_sfc_type), intent(in)    :: tamdar_sfc

   if (trace_use_dull) call da_trace_entry("da_print_stats_tamdar_sfc")

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
      ' Minimum(n,k): ', tamdar_sfc%minimum%u, tamdar_sfc%minimum%v, &
      tamdar_sfc%minimum%t, tamdar_sfc%minimum%p, tamdar_sfc%minimum%q, &
      ' Maximum(n,k): ', tamdar_sfc%maximum%u, tamdar_sfc%maximum%v, &
      tamdar_sfc%maximum%t, &
                        tamdar_sfc%maximum%p, tamdar_sfc%maximum%q
   write(unit=stats_unit, fmt='((a,4(f12.4,10x),e12.4,10x))') &
      ' Average     : ', tamdar_sfc%average%u/real(nu), &
      tamdar_sfc%average%v/real(nv), &
      tamdar_sfc%average%t/real(nt), tamdar_sfc%average%p/real(np), &
      tamdar_sfc%average%q/real(nq), &
      '    RMSE     : ', sqrt(tamdar_sfc%rms_err%u/real(nu)), &
      sqrt(tamdar_sfc%rms_err%v/real(nv)), &
      sqrt(tamdar_sfc%rms_err%t/real(nt)), &
      sqrt(tamdar_sfc%rms_err%p/real(np)), &
      sqrt(tamdar_sfc%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_tamdar_sfc")

end subroutine da_print_stats_tamdar_sfc


subroutine da_transform_xtoy_tamdar_sfc (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),     intent(inout) :: grid
   type (iv_type),    intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n        ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_psfc(:)

   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_tamdar_sfc")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_v(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_t(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_q(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_psfc(iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))

      allocate (ub(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (vb(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xa%u, iv%info(tamdar_sfc), model_u)
      call da_interp_lin_3d (grid%xa%v, iv%info(tamdar_sfc), model_v)
      call da_interp_lin_3d (grid%xa%t, iv%info(tamdar_sfc), model_t)
      call da_interp_lin_3d (grid%xa%q, iv%info(tamdar_sfc), model_q)

      call da_interp_lin_2d (grid%xa%psfc, iv%info(tamdar_sfc), 1, model_psfc)

      call da_interp_lin_3d (grid%xb%u, iv%info(tamdar_sfc), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(tamdar_sfc), vb)

      do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
         if(wind_sd_tamdar)then
            call da_uv_to_sd_lin(y%tamdar_sfc(n)%u,y%tamdar_sfc(n)%v,model_u(1,n),model_v(1,n),ub(1,n),vb(1,n))
         else
            y%tamdar_sfc(n)%u = model_u(1,n)
            y%tamdar_sfc(n)%v = model_v(1,n)
         end if
         y%tamdar_sfc(n)%t = model_t(1,n)
         y%tamdar_sfc(n)%q = model_q(1,n)
         y%tamdar_sfc(n)%p = model_psfc(n)
      end do

      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc(grid, iv, tamdar_sfc, iv%tamdar_sfc(:), y%tamdar_sfc(:))
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_tamdar_sfc")

end subroutine da_transform_xtoy_tamdar_sfc


subroutine da_transform_xtoy_tamdar_sfc_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_tamdar_sfc_adj")

   if (sfc_assi_options == sfc_assi_options_1) then

      allocate (model_u(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_v(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_t(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_q(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (model_psfc(iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))

      allocate (ub(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
      allocate (vb(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))

      call da_interp_lin_3d (grid%xb%u, iv%info(tamdar_sfc), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(tamdar_sfc), vb)

      ! [1.2] Interpolate horizontally:
      do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
         if(wind_sd_tamdar)then
            call da_uv_to_sd_adj(jo_grad_y%tamdar_sfc(n)%u, &
                                 jo_grad_y%tamdar_sfc(n)%v, model_u(1,n), model_v(1,n), ub(1,n), vb(1,n))
         else
            model_u(1,n)  = jo_grad_y%tamdar_sfc(n)%u
            model_v(1,n)  = jo_grad_y%tamdar_sfc(n)%v
         end if
         model_t(1,n)  = jo_grad_y%tamdar_sfc(n)%t
         model_q(1,n)  = jo_grad_y%tamdar_sfc(n)%q
         model_psfc(n) = jo_grad_y%tamdar_sfc(n)%p
      end do

      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(tamdar_sfc), model_u)
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(tamdar_sfc), model_v)
      call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(tamdar_sfc), model_t)
      call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(tamdar_sfc), model_q)

      call da_interp_lin_2d_adj (jo_grad_x%psfc, iv%info(tamdar_sfc), 1, model_psfc)


      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc_adj(grid,iv,tamdar_sfc,iv%tamdar_sfc(:), jo_grad_y%tamdar_sfc(:),jo_grad_x)
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_tamdar_sfc_adj")

end subroutine da_transform_xtoy_tamdar_sfc_adj


subroutine da_get_innov_vector_tamdar_sfc( it, num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD    
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
   real, allocatable :: model_u(:,:)  ! Model value u at oblocation.
   real, allocatable :: model_v(:,:)  ! Model value v at oblocation.
   real, allocatable :: model_t(:,:)  ! Model value t at oblocation.
   real, allocatable :: model_p(:,:)  ! Model value p at oblocation.
   real, allocatable :: model_q(:,:)  ! Model value q at oblocation.
   real, allocatable :: model_hsm(:,:)

   real  :: speed, direction

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.

   real    :: hd, psfcm
   real    :: ho, to, qo
   
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_tamdar_sfc")


   allocate (model_u(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
   allocate (model_v(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
   allocate (model_t(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
   allocate (model_p(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
   allocate (model_q(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))
   allocate (model_hsm(1,iv%info(tamdar_sfc)%n1:iv%info(tamdar_sfc)%n2))

   if ( it > 1 ) then
      do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
         if (iv%tamdar_sfc(n)%u%qc == fails_error_max) iv%tamdar_sfc(n)%u%qc = 0
         if (iv%tamdar_sfc(n)%v%qc == fails_error_max) iv%tamdar_sfc(n)%v%qc = 0
         if (iv%tamdar_sfc(n)%t%qc == fails_error_max) iv%tamdar_sfc(n)%t%qc = 0
         if (iv%tamdar_sfc(n)%p%qc == fails_error_max) iv%tamdar_sfc(n)%p%qc = 0
         if (iv%tamdar_sfc(n)%q%qc == fails_error_max) iv%tamdar_sfc(n)%q%qc = 0
      end do
   end if

   if (sfc_assi_options == sfc_assi_options_1) then
      do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
         ! [1.1] Get horizontal interpolation weights:

         i   = iv%info(tamdar_sfc)%i(1,n)
         j   = iv%info(tamdar_sfc)%j(1,n)
         dx  = iv%info(tamdar_sfc)%dx(1,n)
         dy  = iv%info(tamdar_sfc)%dy(1,n)
         dxm = iv%info(tamdar_sfc)%dxm(1,n)
         dym = iv%info(tamdar_sfc)%dym(1,n)

         ! Surface correction

         iv%tamdar_sfc(n)%p%inv = ob%tamdar_sfc(n)%p
         iv%tamdar_sfc(n)%t%inv = ob%tamdar_sfc(n)%t
         iv%tamdar_sfc(n)%q%inv = ob%tamdar_sfc(n)%q
         iv%tamdar_sfc(n)%u%inv = ob%tamdar_sfc(n)%u
         iv%tamdar_sfc(n)%v%inv = ob%tamdar_sfc(n)%v

         if (iv % tamdar_sfc(n) % h > missing_r) then
            do k=kts,kte
               v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                  + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
            end do

            hd = v_h(kts) - iv % tamdar_sfc(n) % h

            if (abs(hd) <= Max_StHeight_Diff .or. anal_type_verify ) then
               if (iv % tamdar_sfc(n) % h < v_h(kts)) then
                  iv%info(tamdar_sfc)%zk(:,n) = 1.0+1.0e-6
                  call da_obs_sfc_correction(iv%info(tamdar_sfc), iv%tamdar_sfc(n), n, grid%xb)
               else
                  call da_to_zk(iv % tamdar_sfc(n) % h, v_h, v_interp_h, iv%info(tamdar_sfc)%zk(1,n))
               end if
            end if
         else if (ob % tamdar_sfc(n) % p > 1.0) then
            do k=kts,kte
               v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                      + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
            end do

            call da_to_zk(ob % tamdar_sfc(n) % p, v_p, v_interp_p, iv%info(tamdar_sfc)%zk(1,n))

            if (iv%info(tamdar_sfc)%zk(1,n) < 0.0 .and. .not. anal_type_verify) then
               iv % tamdar_sfc(n) % p % inv = missing_r
               iv % tamdar_sfc(n) % p % qc  = missing_data
               iv%info(tamdar_sfc)%zk(:,n) = 1.0+1.0e-6
            end if
         end if
      end do

      call da_convert_zk (iv%info(tamdar_sfc))

      if (.not.anal_type_verify ) then
         do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
            if (iv%info(tamdar_sfc)%zk(1,n) < 0.0) then
               iv % tamdar_sfc(n) % u % qc = missing_data
               iv % tamdar_sfc(n) % v % qc = missing_data
               iv % tamdar_sfc(n) % t % qc = missing_data
               iv % tamdar_sfc(n) % q % qc = missing_data
               iv % tamdar_sfc(n) % p % qc = missing_data
            end if
         end do
      end if

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xb%u, iv%info(tamdar_sfc), model_u)
      call da_interp_lin_3d (grid%xb%v, iv%info(tamdar_sfc), model_v)
      call da_interp_lin_3d (grid%xb%t, iv%info(tamdar_sfc), model_t)
      call da_interp_lin_3d (grid%xb%q, iv%info(tamdar_sfc), model_q)
      call da_interp_lin_3d (grid%xb%p, iv%info(tamdar_sfc), model_p)

   else if (sfc_assi_options == sfc_assi_options_2) then
      ! 1.2.1 Surface assimilation approach 2(10-m u, v, 2-m t, q, and sfc_p)
      call da_interp_lin_3d (grid%xb%u, iv%info(tamdar_sfc), model_u)
      call da_interp_lin_3d (grid%xb%v, iv%info(tamdar_sfc), model_v)
      call da_interp_lin_2d (grid%xb%t2,   iv%info(tamdar_sfc), 1,model_t)
      call da_interp_lin_2d (grid%xb%q2,   iv%info(tamdar_sfc), 1,model_q)
      call da_interp_lin_2d (grid%xb%psfc, iv%info(tamdar_sfc), 1,model_p)

      do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2

         iv%tamdar_sfc(n)%p%inv = ob%tamdar_sfc(n)%p
         iv%tamdar_sfc(n)%t%inv = ob%tamdar_sfc(n)%t
         iv%tamdar_sfc(n)%q%inv = ob%tamdar_sfc(n)%q
         iv%tamdar_sfc(n)%u%inv = ob%tamdar_sfc(n)%u
         iv%tamdar_sfc(n)%v%inv = ob%tamdar_sfc(n)%v
         if (iv%tamdar_sfc(n)%p%qc >= 0) then
            ! model surface p, t, q, h at observed site:

            call da_interp_lin_2d_partial (grid%xb%terr, iv%info(tamdar_sfc), 1, n, n, model_hsm(:,n))

            ho = iv%tamdar_sfc(n)%h
            to = -888888.0
            qo = -888888.0

            if (iv%tamdar_sfc(n)%t%qc >= 0 .and. iv%tamdar_sfc(n)%q%qc >= 0) then
               to = ob%tamdar_sfc(n)%t
               qo = ob%tamdar_sfc(n)%q
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to, qo)
            else if (iv%tamdar_sfc(n)%t%qc >= 0 .and. iv%tamdar_sfc(n)%q%qc < 0) then
               to = ob%tamdar_sfc(n)%t
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho,to)
            else
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho)
            end if

            ! Pressure at the observed height:
            model_p(1,n) = psfcm
         end if
      end do
   end if

   do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
      !-----------------------------------------------------------------------
      ! [3.0] Fast interpolation:
      !-----------------------------------------------------------------------

       if(wind_sd_tamdar)then
         call da_ffdduv_model (speed,direction,model_u(1,n), model_v(1,n), convert_uv2fd)

         if (ob%tamdar_sfc(n)%u > missing_r .AND. iv%tamdar_sfc(n)%u%qc >= obs_qc_pointer) then
             iv%tamdar_sfc(n)%u%inv = iv%tamdar_sfc(n)%u%inv - speed
         else
             iv % tamdar_sfc(n) % u % inv = 0.0
         end if

         if (ob%tamdar_sfc(n)%v > missing_r .AND. iv%tamdar_sfc(n)%v%qc >= obs_qc_pointer) then
             iv%tamdar_sfc(n)%v%inv = iv%tamdar_sfc(n)%v%inv - direction
             if (iv%tamdar_sfc(n)%v%inv > 180.0 ) iv%tamdar_sfc(n)%v%inv = iv%tamdar_sfc(n)%v%inv - 360.0
             if (iv%tamdar_sfc(n)%v%inv < -180.0 ) iv%tamdar_sfc(n)%v%inv = iv%tamdar_sfc(n)%v%inv + 360.0
         else
             iv % tamdar_sfc(n) % v % inv = 0.0
         end if
       else
          if (ob % tamdar_sfc(n) % u > missing_r .AND. iv % tamdar_sfc(n) % u % qc >= obs_qc_pointer) then
              iv % tamdar_sfc(n) % u % inv = iv%tamdar_sfc(n)%u%inv - model_u(1,n)
          else
              iv % tamdar_sfc(n) % u % inv = 0.0
          end if

          if (ob % tamdar_sfc(n) % v > missing_r .AND. iv % tamdar_sfc(n) % v % qc >= obs_qc_pointer) then
              iv % tamdar_sfc(n) % v % inv = iv%tamdar_sfc(n)%v%inv - model_v(1,n)
          else
              iv % tamdar_sfc(n) % v % inv = 0.0
          end if
       end if


      if (ob % tamdar_sfc(n) % p > 0.0 .AND. iv % tamdar_sfc(n) % p % qc >= obs_qc_pointer) then
         iv % tamdar_sfc(n) % p % inv = iv%tamdar_sfc(n)%p%inv - model_p(1,n)
      else
         iv % tamdar_sfc(n) % p % inv = 0.0
      end if

      if (ob % tamdar_sfc(n) % t > 0.0 .AND. iv % tamdar_sfc(n) % t % qc >= obs_qc_pointer) then
         iv % tamdar_sfc(n) % t % inv = iv%tamdar_sfc(n)%t%inv - model_t(1,n)
      else
         iv % tamdar_sfc(n) % t % inv = 0.0
      end if

      if (ob % tamdar_sfc(n) % q > 0.0 .AND. iv % tamdar_sfc(n) % q % qc >= obs_qc_pointer) then
         iv % tamdar_sfc(n) % q % inv = iv%tamdar_sfc(n)%q%inv - model_q(1,n)
      else
         iv % tamdar_sfc(n) % q % inv = 0.0
      end if
   end do

   !-----------------------------------------------------------------------
   !     [5.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_tamdar_sfc(iv,ob, it, num_qcstat_conv) 

   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_p)
   deallocate (model_q)
   deallocate (model_hsm)
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_tamdar_sfc")

end subroutine da_get_innov_vector_tamdar_sfc


subroutine da_check_max_iv_tamdar_sfc(iv,ob, it, num_qcstat_conv)                                 

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

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_tamdar_sfc")


   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(tamdar_sfc)%n1,iv%info(tamdar_sfc)%n2
         if(.not. qc_rej_both)then
            if(wind_sd_tamdar)then
               failed=.false.
               if( iv%tamdar_sfc(n)%u%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%u, max_error_spd,failed)
                   if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,tamdar_sfc,1,1) = num_qcstat_conv(1,tamdar_sfc,1,1) + 1
                       if(failed) then
                          num_qcstat_conv(2,tamdar_sfc,1,1) = num_qcstat_conv(2,tamdar_sfc,1,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'tamdar_sfc',ob_vars(1),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
                       end if
                   end if
               end if

               failed=.false.
               if( iv%tamdar_sfc(n)%v%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%v, max_error_dir,failed)
                   if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,tamdar_sfc,2,1) = num_qcstat_conv(1,tamdar_sfc,2,1) + 1
                       if(failed)then
                          num_qcstat_conv(2,tamdar_sfc,2,1) = num_qcstat_conv(2,tamdar_sfc,2,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'tamdar_sfc',ob_vars(2),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
                       end if
                   end if
               end if
            else
               failed=.false.
               if( iv%tamdar_sfc(n)%u%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%u, max_error_uv,failed)
                   if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,tamdar_sfc,1,1) = num_qcstat_conv(1,tamdar_sfc,1,1) + 1
                       if(failed) then
                          num_qcstat_conv(2,tamdar_sfc,1,1) = num_qcstat_conv(2,tamdar_sfc,1,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'tamdar_sfc',ob_vars(1),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
                       end if
                   end if
               end if

               failed=.false.
               if( iv%tamdar_sfc(n)%v%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%v, max_error_uv,failed)
                   if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,tamdar_sfc,2,1) = num_qcstat_conv(1,tamdar_sfc,2,1) + 1
                       if(failed)then
                          num_qcstat_conv(2,tamdar_sfc,2,1) = num_qcstat_conv(2,tamdar_sfc,2,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'tamdar_sfc',ob_vars(2),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
                       end if
                   end if
               end if
            end if

            if(wind_sd_tamdar)then
               if(iv%tamdar_sfc(n)%u%qc == fails_error_max .or.  abs(iv%tamdar_sfc(n)%u%inv) >= max_omb_spd) then
                  iv%tamdar_sfc(n)%u%qc = fails_error_max
                  iv%tamdar_sfc(n)%u%inv = 0.0
               endif
               if(iv%tamdar_sfc(n)%v%qc == fails_error_max .or.  abs(iv%tamdar_sfc(n)%v%inv) >= max_omb_dir) then
                  iv%tamdar_sfc(n)%v%qc = fails_error_max
                  iv%tamdar_sfc(n)%v%inv = 0.0
               endif
            endif
         else
            failed1=.false.
            failed2=.false.

            if( iv%tamdar_sfc(n)%v%qc >= obs_qc_pointer .or. iv%tamdar_sfc(n)%u%qc >= obs_qc_pointer )  then
                if(wind_sd_tamdar)then
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%u, max_error_spd,failed1)
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%v, max_error_dir,failed2)
                else
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%u, max_error_uv,failed1)
                   call da_max_error_qc (it,iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%v, max_error_uv,failed2)
                endif
            endif

            if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
                num_qcstat_conv(1,tamdar_sfc,1,1) = num_qcstat_conv(1,tamdar_sfc,1,1) + 1
                num_qcstat_conv(1,tamdar_sfc,2,1) = num_qcstat_conv(1,tamdar_sfc,2,1) + 1

                if(failed1 .or. failed2) then
                   num_qcstat_conv(2,tamdar_sfc,1,1) = num_qcstat_conv(2,tamdar_sfc,1,1) + 1
                   write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'tamdar_sfc',ob_vars(1),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
                   num_qcstat_conv(2,tamdar_sfc,2,1) = num_qcstat_conv(2,tamdar_sfc,2,1) + 1
                   write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'tamdar_sfc',ob_vars(2),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
                endif
            endif

            if(wind_sd_tamdar)then
               if(iv%tamdar_sfc(n)%u%qc == fails_error_max .or. iv%tamdar_sfc(n)%v%qc == fails_error_max .or. &
                  abs(iv%tamdar_sfc(n)%v%inv) >= max_omb_dir .or. abs(iv%tamdar_sfc(n)%u%inv) >= max_omb_spd )then
                  iv%tamdar_sfc(n)%u%qc = fails_error_max
                  iv%tamdar_sfc(n)%v%qc = fails_error_max
                  iv%tamdar_sfc(n)%u%inv = 0.0
                  iv%tamdar_sfc(n)%v%inv = 0.0
               endif
            else
               if(iv%tamdar_sfc(n)%u%qc == fails_error_max .or. iv%tamdar_sfc(n)%v%qc == fails_error_max ) then
                  iv%tamdar_sfc(n)%u%qc = fails_error_max
                  iv%tamdar_sfc(n)%v%qc = fails_error_max
                  iv%tamdar_sfc(n)%u%inv = 0.0
                  iv%tamdar_sfc(n)%v%inv = 0.0
               endif
            endif
         endif

      failed=.false.
      if( iv%tamdar_sfc(n)%t%qc >= obs_qc_pointer )  then 
      call da_max_error_qc (it, iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%t, max_error_t , failed)
      if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
      num_qcstat_conv(1,tamdar_sfc,3,1)= num_qcstat_conv(1,tamdar_sfc,3,1) + 1
      if(failed) then
      num_qcstat_conv(2,tamdar_sfc,3,1)= num_qcstat_conv(2,tamdar_sfc,3,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'tamdar_sfc',ob_vars(3),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%tamdar_sfc(n)%p%qc >= obs_qc_pointer ) then
      call da_max_error_qc (it, iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%p, max_error_p , failed)         
      if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
      num_qcstat_conv(1,tamdar_sfc,5,1)= num_qcstat_conv(1,tamdar_sfc,5,1) + 1
      if(failed) then
      num_qcstat_conv(2,tamdar_sfc,5,1)= num_qcstat_conv(2,tamdar_sfc,5,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'tamdar_sfc',ob_vars(5),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%tamdar_sfc(n)%q%qc >= obs_qc_pointer ) then
       if( iv%tamdar_sfc(n)%t%qc == fails_error_max .or. iv%tamdar_sfc(n)%p%qc == fails_error_max) then
       failed=.true.
       iv%tamdar_sfc(n)%q%qc  = fails_error_max
       iv%tamdar_sfc(n)%q%inv = 0.0
       else
       call da_max_error_qc (it, iv%info(tamdar_sfc), n, iv%tamdar_sfc(n)%q, max_error_q , failed)
       endif
      if( iv%info(tamdar_sfc)%proc_domain(1,n) ) then
      num_qcstat_conv(1,tamdar_sfc,4,1)= num_qcstat_conv(1,tamdar_sfc,4,1) + 1
      if(failed) then
      num_qcstat_conv(2,tamdar_sfc,4,1)= num_qcstat_conv(2,tamdar_sfc,4,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'tamdar_sfc',ob_vars(4),iv%info(tamdar_sfc)%lat(1,n),iv%info(tamdar_sfc)%lon(1,n),0.01*ob%tamdar_sfc(n)%p
      end if
      end if
      end if

   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_tamdar_sfc")

end subroutine da_check_max_iv_tamdar_sfc
subroutine da_calculate_grady_tamdar_sfc(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector              
   !-------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_tamdar_sfc")

   do n=1, iv%info(tamdar_sfc)%nlocal
      if (iv%tamdar_sfc(n)%u%qc < obs_qc_pointer) re%tamdar_sfc(n)%u = 0.0
      if (iv%tamdar_sfc(n)%v%qc < obs_qc_pointer) re%tamdar_sfc(n)%v = 0.0
      if (iv%tamdar_sfc(n)%t%qc < obs_qc_pointer) re%tamdar_sfc(n)%t = 0.0
      if (iv%tamdar_sfc(n)%p%qc < obs_qc_pointer) re%tamdar_sfc(n)%p = 0.0
      if (iv%tamdar_sfc(n)%q%qc < obs_qc_pointer) re%tamdar_sfc(n)%q = 0.0

      if (iv%tamdar_sfc(n)%u%qc < obs_qc_pointer) re%tamdar_sfc(n)%u = 0.0
      if (iv%tamdar_sfc(n)%v%qc < obs_qc_pointer) re%tamdar_sfc(n)%v = 0.0
      if (iv%tamdar_sfc(n)%t%qc < obs_qc_pointer) re%tamdar_sfc(n)%t = 0.0
      if (iv%tamdar_sfc(n)%p%qc < obs_qc_pointer) re%tamdar_sfc(n)%p = 0.0
      if (iv%tamdar_sfc(n)%q%qc < obs_qc_pointer) re%tamdar_sfc(n)%q = 0.0

      jo_grad_y%tamdar_sfc(n)%u = -re%tamdar_sfc(n)%u / &
          (iv%tamdar_sfc(n)%u%error * iv%tamdar_sfc(n)%u%error)
      jo_grad_y%tamdar_sfc(n)%v = -re%tamdar_sfc(n)%v / &
          (iv%tamdar_sfc(n)%v%error * iv%tamdar_sfc(n)%v%error)
      jo_grad_y%tamdar_sfc(n)%t = -re%tamdar_sfc(n)%t / &
          (iv%tamdar_sfc(n)%t%error * iv%tamdar_sfc(n)%t%error)
      jo_grad_y%tamdar_sfc(n)%p = -re%tamdar_sfc(n)%p / &
          (iv%tamdar_sfc(n)%p%error * iv%tamdar_sfc(n)%p%error)
      jo_grad_y%tamdar_sfc(n)%q = -re%tamdar_sfc(n)%q / &
          (iv%tamdar_sfc(n)%q%error * iv%tamdar_sfc(n)%q%error)
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_tamdar_sfc")
     
end subroutine da_calculate_grady_tamdar_sfc



end module da_tamdar

