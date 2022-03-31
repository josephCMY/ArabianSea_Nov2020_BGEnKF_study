












module da_sound

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, missing_data, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, sfc_assi_options, no_buddies, fails_error_max, &
      fails_buddy_check, check_buddy, check_buddy_print, check_buddy_unit, &
      buddy_weight , max_buddy_uv, max_buddy_t, max_buddy_p, max_buddy_rh, &
      max_stheight_diff,test_dm_exact, anal_type_verify, &
      kms,kme,kts,kte,sfc_assi_options_1,sfc_assi_options_2, num_procs, comm, &
      trace_use_dull, sound, sonde_sfc, position_lev_dependant, max_ext_its,qcstat_conv_unit,ob_vars, &
      convert_fd2uv,convert_uv2fd,max_error_spd,max_error_dir,max_omb_spd,max_omb_dir,pi,qc_rej_both, &
      wind_sd_sound, wind_stats_sd
   use da_grid_definitions, only : da_ffdduv,da_ffdduv_model, da_ffdduv_diagnose


   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use module_domain, only : domain
   use da_interpolation, only : da_to_zk, da_interp_lin_3d, &
      da_interp_lin_3d_adj, da_interp_lin_2d, da_interp_lin_2d_adj, da_interp_lin_2d_partial
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_obs_sfc_correction, da_convert_zk,&
                        da_buddy_qc, da_get_print_lvl
   use da_par_util, only : da_proc_stats_combine, &
      da_deallocate_global_sound, da_to_global_sound, da_to_global_sonde_sfc, &
      da_deallocate_global_sonde_sfc
   use da_par_util1, only : da_proc_sum_int 
   use da_physics, only : da_sfc_pre, da_transform_xtopsfc, &
      da_transform_xtopsfc_adj, da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_sound1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: q                        
   end type residual_sound1_type

   type maxmin_sound_stats_type
      type (maxmin_type)         :: u, v, t, q
   end type maxmin_sound_stats_type

   type stats_sound_type
      type (maxmin_sound_stats_type)  :: maximum, minimum
      type (residual_sound1_type)     :: average, rms_err
   end type stats_sound_type

   

   type residual_sonde_sfc1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: p                        
      real          :: q                        
   end type residual_sonde_sfc1_type

   type maxmin_sonde_sfc_stats_type
      type (maxmin_type)         :: u, v, t, p, q
   end type maxmin_sonde_sfc_stats_type

   type stats_sonde_sfc_type
      type (maxmin_sonde_sfc_stats_type)  :: maximum, minimum
      type (residual_sonde_sfc1_type)     :: average, rms_err
   end type stats_sonde_sfc_type

   include 'mpif.h'

contains

subroutine da_ao_stats_sound (stats_unit, iv, re, ob)

   !-------------------------------------------------------------
   ! Purpose: TBD   
   !-------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type (y_type),  intent (in)    :: re            ! A - O
   type(y_type),   intent (in)    :: ob            ! Observation structure.

   type (stats_sound_type) :: stats
   integer                 :: nu, nv, nt, nq
   integer                 :: n, k
   real                    :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_sound")

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

   stats%average = residual_sound1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(sound)%nlocal
      do k=1, iv%info(sound)%levels(n)

            u_inc = re%sound(n)%u(k)
            v_inc = re%sound(n)%v(k)
            u_obs = ob%sound(n)%u(k)
            v_obs = ob%sound(n)%v(k)

            if (.not. wind_sd_sound .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%sound(n)%u(k)%qc, iv%sound(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_sound .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%sound(n)%u(k)%qc, iv%sound(n)%v(k)%qc, convert_fd2uv)

         if (iv%info(sound)%proc_domain(1,n)) then
            call da_stats_calculate (n, k, iv%sound(n)%u(k)%qc, & 
               u_inc, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%sound(n)%v(k)%qc, & 
               v_inc, nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
            call da_stats_calculate (n, k, iv%sound(n)%t(k)%qc, & 
               re%sound(n)%t(k), nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate (n, k, iv%sound(n)%q(k)%qc, & 
               re%sound(n)%q(k), nq, &
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
   iv%nstats(sound) = nu + nv + nt + nq

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for sound'
         call da_print_stats_sound(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_sound")

end subroutine da_ao_stats_sound


subroutine da_jo_and_grady_sound(iv, re, jo, jo_grad_y)

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
   
   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_sound")

   jo % sound_u = 0.0
   jo % sound_v = 0.0
   jo % sound_t = 0.0
   jo % sound_q = 0.0

   if (test_dm_exact) then
      if (iv%info(sound)%ntotal == 0) then 
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_sound")
         return
      end if
   else
      if (iv%info(sound)%nlocal < 1) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_sound")
         return
      end if
   end if

   do n=1, iv%info(sound)%nlocal
       do k=1, iv%info(sound)%levels(n)
          jo_grad_y%sound(n)%u(k) = -re%sound(n)%u(k) / (iv%sound(n)%u(k)%error * iv%sound(n)%u(k)%error)
          jo_grad_y%sound(n)%v(k) = -re%sound(n)%v(k) / (iv%sound(n)%v(k)%error * iv%sound(n)%v(k)%error)
          jo_grad_y%sound(n)%t(k) = -re%sound(n)%t(k) / (iv%sound(n)%t(k)%error * iv%sound(n)%t(k)%error)
          jo_grad_y%sound(n)%q(k) = -re%sound(n)%q(k) / (iv%sound(n)%q(k)%error * iv%sound(n)%q(k)%error)
      end do
   end do

   ! Bitwise-exact reduction preserves operation order of serial code for 
   ! testing, at the cost of much-increased run-time.  Turn it off when not 
   ! testing.  This will always be .false. for a serial or 1-MPI-process run.  
   if (test_dm_exact) then
      ! collect all obs in serial order and allocate global objects
      call da_to_global_sound(iv, re, jo_grad_y, iv_glob, re_glob, jo_grad_y_glob)
      ! perform remaining computations
      call da_jo_sound_uvtq(iv_glob, re_glob, jo_grad_y_glob, jo)
      ! free global objects
      call da_deallocate_global_sound(iv_glob, re_glob, jo_grad_y_glob)
   else
      ! perform remaining computations
      call da_jo_sound_uvtq(iv, re, jo_grad_y, jo)
   end if

   jo % sound_u = 0.5 * jo % sound_u
   jo % sound_v = 0.5 * jo % sound_v
   jo % sound_t = 0.5 * jo % sound_t
   jo % sound_q = 0.5 * jo % sound_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_sound")

end subroutine da_jo_and_grady_sound


subroutine da_jo_sound_uvtq (iv, re, jo_grad_y, jo)

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

   if (trace_use_dull) call da_trace_entry("da_jo_sound_uvtq")

   do n=1, iv%info(sound)%nlocal
      do k=1, iv%info(sound)%levels(n)
         if (iv%info(sound)%proc_domain(1,n)) then
            jo % sound_u = jo % sound_u - re%sound(n)%u(k) * jo_grad_y%sound(n)%u(k)
            jo % sound_v = jo % sound_v - re%sound(n)%v(k) * jo_grad_y%sound(n)%v(k)
            jo % sound_t = jo % sound_t - re%sound(n)%t(k) * jo_grad_y%sound(n)%t(k)
            jo % sound_q = jo % sound_q - re%sound(n)%q(k) * jo_grad_y%sound(n)%q(k)
         end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_jo_sound_uvtq")

end subroutine da_jo_sound_uvtq


subroutine da_residual_sound(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for sound obs
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

   if (trace_use_dull) call da_trace_entry("da_residual_sound")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(sound)%nlocal
      do k=1, iv%info(sound)%levels(n)
         np_available = np_available + 4

         re%sound(n)%u(k) = da_residual(n, k, y%sound(n)%u(k), iv%sound(n)%u(k), n_obs_bad % u)
         re%sound(n)%v(k) = da_residual(n, k, y%sound(n)%v(k), iv%sound(n)%v(k), n_obs_bad % v)
         re%sound(n)%t(k) = da_residual(n, k, y%sound(n)%t(k), iv%sound(n)%t(k), n_obs_bad % t)
         re%sound(n)%q(k) = da_residual(n, k, y%sound(n)%q(k), iv%sound(n)%q(k), n_obs_bad % q)
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

   if (trace_use_dull) call da_trace_exit("da_residual_sound")

end subroutine da_residual_sound


subroutine da_oi_stats_sound (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! OI
   type(y_type),   intent (in)      :: ob            ! Observation structure.

   type (stats_sound_type)          :: stats
   integer                          :: nu, nv, nt, nq
   integer                          :: n, k
   real                             :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_sound")

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
   stats%average = residual_sound1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(sound)%nlocal
      do k=1, iv%info(sound)%levels(n)
!         if (iv%info(sound)%proc_domain(k,n)) then
         if (iv%info(sound)%proc_domain(1,n)) then

            u_inv = iv%sound(n)%u(k)%inv
            v_inv = iv%sound(n)%v(k)%inv
            u_obs = ob%sound(n)%u(k)
            v_obs = ob%sound(n)%v(k)

            if (.not. wind_sd_sound .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%sound(n)%u(k)%qc, iv%sound(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_sound .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%sound(n)%u(k)%qc, iv%sound(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate(iv%info(sound)%obs_global_index(n), &
               k, iv%sound(n)%u(k)%qc, u_inv, nu, &
               stats%minimum%u, stats%maximum%u, stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(sound)%obs_global_index(n), &
               k, iv%sound(n)%v(k)%qc, v_inv, nv, &
               stats%minimum%v, stats%maximum%v, stats%average%v, stats%rms_err%v)
            call da_stats_calculate(iv%info(sound)%obs_global_index(n), &
               k, iv%sound(n)%t(k)%qc, iv%sound(n)%t(k)%inv, nt, &
               stats%minimum%t, stats%maximum%t, stats%average%t, stats%rms_err%t)
            call da_stats_calculate(iv%info(sound)%obs_global_index(n), &
               k, iv%sound(n)%q(k)%qc, iv%sound(n)%q(k)%inv, nq, &
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for sound'
         call da_print_stats_sound(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_sound")

end subroutine da_oi_stats_sound


subroutine da_print_stats_sound(stats_unit, nu, nv, nt, nq, Sound)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nu, nv, nt, nq
   type (stats_sound_type), intent(in)    :: sound

   if (trace_use_dull) call da_trace_entry("da_print_stats_sound")

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
      ' Minimum(n,k): ', sound%minimum%u, sound%minimum%v, &
                         sound%minimum%t, sound%minimum%q, &
      ' Maximum(n,k): ', sound%maximum%u, sound%maximum%v, &
                         sound%maximum%t, sound%maximum%q
   write(unit=stats_unit, fmt='((a,3(f12.4,10x),e12.4,10x))') &
      ' Average     : ', sound%average%u/real(nu), &
                         sound%average%v/real(nv), &
                         sound%average%t/real(nt), &
                         sound%average%q/real(nq), &
      '    RMSE     : ', sqrt(sound%rms_err%u/real(nu)), &
                         sqrt(sound%rms_err%v/real(nv)), &
                         sqrt(sound%rms_err%t/real(nt)), &
                         sqrt(sound%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_sound")

end subroutine da_print_stats_sound


subroutine da_transform_xtoy_sound (grid, iv, y)

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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_sound")

   allocate (u(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (v(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (t(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (q(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))

   allocate (ub(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (vb(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
  
   call da_interp_lin_3d (grid%xa%u, iv%info(sound), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(sound), v)
   call da_interp_lin_3d (grid%xa%t, iv%info(sound), t)
   call da_interp_lin_3d (grid%xa%q, iv%info(sound), q)

   call da_interp_lin_3d (grid%xb%u, iv%info(sound), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(sound), vb)


   do n=iv%info(sound)%n1,iv%info(sound)%n2
      do k = 1, iv%info(sound)%levels(n)
         if(wind_sd_sound) then
            call da_uv_to_sd_lin(y%sound(n)%u(k),y%sound(n)%v(k),u(k,n),v(k,n),ub(k,n),vb(k,n))
         else
            y%sound(n)%u(k) = u(k,n)
            y%sound(n)%v(k) = v(k,n)
         end if
      end do
      y%sound(n)%t(:) = t(1:size(y%sound(n)%t),n)
      y%sound(n)%q(:) = q(1:size(y%sound(n)%q),n)
   end do

   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_sound")

end subroutine da_transform_xtoy_sound


subroutine da_transform_xtoy_sound_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_sound_adj")

   allocate (u(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (v(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (t(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (q(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))

   allocate (ub(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (vb(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))

   call da_interp_lin_3d (grid%xb%u, iv%info(sound), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(sound), vb)

   do n=iv%info(sound)%n1,iv%info(sound)%n2
       do k = 1, iv%info(sound)%levels(n)
         if(wind_sd_sound) then
            call da_uv_to_sd_adj(jo_grad_y%sound(n)%u(k), &
                                 jo_grad_y%sound(n)%v(k), u(k,n), v(k,n), ub(k,n), vb(k,n))
         else
            u(k,n) = jo_grad_y%sound(n)%u(k)
            v(k,n) = jo_grad_y%sound(n)%v(k)
         end if
      end do
      t(1:size(jo_grad_y%sound(n)%t),n) = jo_grad_y%sound(n)%t(:)
      q(1:size(jo_grad_y%sound(n)%q),n) = jo_grad_y%sound(n)%q(:)
   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(sound), u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(sound), v)
   call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(sound), t)
   call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(sound), q)
   
   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_sound_adj")

end subroutine da_transform_xtoy_sound_adj


subroutine da_check_max_iv_sound(iv, it,num_qcstat_conv)

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

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_sound")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(sound)%n1,iv%info(sound)%n2
      do k = 1, iv%info(sound)%levels(n)
         call da_get_print_lvl(iv%sound(n)%p(k),ipr) 

         if(.not. qc_rej_both)then
            if(wind_sd_sound)then
               failed=.false.
               if( iv%sound(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%u(k), max_error_spd,failed)
                   if( iv%info(sound)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,sound,1,ipr) = num_qcstat_conv(1,sound,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,sound,1,ipr) = num_qcstat_conv(2,sound,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'sound',ob_vars(1),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
                       end if
                   end if
                end if

                failed=.false.
                if( iv%sound(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%v(k), max_error_dir,failed)
                    if( iv%info(sound)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,sound,2,ipr) = num_qcstat_conv(1,sound,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,sound,2,ipr) = num_qcstat_conv(2,sound,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'sound',ob_vars(2),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
                        end if
                    end if
                end if
             else
                failed=.false.
                if( iv%sound(n)%u(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%u(k), max_error_uv,failed)
                    if( iv%info(sound)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,sound,1,ipr) = num_qcstat_conv(1,sound,1,ipr) + 1
                        if(failed) then
                           num_qcstat_conv(2,sound,1,ipr) = num_qcstat_conv(2,sound,1,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'sound',ob_vars(1),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
                        end if
                    end if
                 end if

                 failed=.false.
                 if( iv%sound(n)%v(k)%qc >= obs_qc_pointer ) then
                     call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%v(k), max_error_uv,failed)
                     if( iv%info(sound)%proc_domain(k,n) ) then
                         num_qcstat_conv(1,sound,2,ipr) = num_qcstat_conv(1,sound,2,ipr) + 1
                         if(failed)then
                            num_qcstat_conv(2,sound,2,ipr) = num_qcstat_conv(2,sound,2,ipr) + 1
                            write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                            'sound',ob_vars(2),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
                         end if
                     end if
                 end if
              end if

              if(wind_sd_sound)then
                 if(iv%sound(n)%u(k)%qc == fails_error_max .or. abs(iv%sound(n)%u(k)%inv) >= max_omb_spd) then
                    iv%sound(n)%u(k)%qc = fails_error_max
                    iv%sound(n)%u(k)%inv = 0.0
                 endif
                 if(iv%sound(n)%v(k)%qc == fails_error_max .or. abs(iv%sound(n)%v(k)%inv) >= max_omb_dir) then
                    iv%sound(n)%v(k)%qc = fails_error_max
                    iv%sound(n)%v(k)%inv = 0.0
                 endif
              endif
           else
              failed1=.false.
              failed2=.false.

              if( iv%sound(n)%v(k)%qc >= obs_qc_pointer .or. iv%sound(n)%u(k)%qc >= obs_qc_pointer )  then
                  if(wind_sd_sound)then
                     call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%u(k), max_error_spd,failed1)
                     call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%v(k), max_error_dir,failed2)
                  else
                     call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%u(k), max_error_uv,failed1)
                     call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%v(k), max_error_uv,failed2)
                  endif
              endif

              if( iv%info(sound)%proc_domain(k,n) ) then
                  num_qcstat_conv(1,sound,1,ipr) = num_qcstat_conv(1,sound,1,ipr) + 1
                  num_qcstat_conv(1,sound,2,ipr) = num_qcstat_conv(1,sound,2,ipr) + 1

                  if(failed1 .or. failed2) then
                     num_qcstat_conv(2,sound,1,ipr) = num_qcstat_conv(2,sound,1,ipr) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'sound',ob_vars(1),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
                    num_qcstat_conv(2,sound,2,ipr) = num_qcstat_conv(2,sound,2,ipr) + 1
                    write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'sound',ob_vars(2),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
                 endif
              endif

              if(wind_sd_sound)then
                 if(iv%sound(n)%u(k)%qc == fails_error_max .or. iv%sound(n)%v(k)%qc == fails_error_max .or. &
                    abs(iv%sound(n)%v(k)%inv) >= max_omb_dir .or. abs(iv%sound(n)%u(k)%inv) >= max_omb_spd )then
                    iv%sound(n)%u(k)%qc = fails_error_max
                    iv%sound(n)%v(k)%qc = fails_error_max
                    iv%sound(n)%u(k)%inv = 0.0
                    iv%sound(n)%v(k)%inv = 0.0
                 endif
              else
                 if(iv%sound(n)%u(k)%qc == fails_error_max .or. iv%sound(n)%v(k)%qc == fails_error_max ) then
                    iv%sound(n)%u(k)%qc = fails_error_max
                    iv%sound(n)%v(k)%qc = fails_error_max
                    iv%sound(n)%u(k)%inv = 0.0
                    iv%sound(n)%v(k)%inv = 0.0
                 endif
              endif
           endif


         failed=.false.
         if( iv%sound(n)%t(k)%qc >= obs_qc_pointer )  then
         call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%t(k), max_error_t ,failed)
         if( iv%info(sound)%proc_domain(k,n) ) then
                    num_qcstat_conv(1,sound,3,ipr) = num_qcstat_conv(1,sound,3,ipr) + 1
         if(failed) then
          num_qcstat_conv(2,sound,3,ipr) = num_qcstat_conv(2,sound,3,ipr) + 1
           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'sound',ob_vars(3),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
         end if
         end if
         end if

         failed=.false.
         if( iv%sound(n)%q(k)%qc >= obs_qc_pointer ) then 
          if( iv%sound(n)%t(k)%qc == fails_error_max ) then
          failed=.true.
          iv%sound(n)%q(k)%qc  = fails_error_max
          iv%sound(n)%q(k)%inv = 0.0
          else
          call da_max_error_qc (it,iv%info(sound), n, iv%sound(n)%q(k), max_error_q ,failed)
          endif
         if( iv%info(sound)%proc_domain(k,n) ) then
                    num_qcstat_conv(1,sound,4,ipr) = num_qcstat_conv(1,sound,4,ipr) + 1
         if(failed) then
         num_qcstat_conv(2,sound,4,ipr) = num_qcstat_conv(2,sound,4,ipr) + 1
           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'sound',ob_vars(4),iv%info(sound)%lat(k,n),iv%info(sound)%lon(k,n),0.01*iv%sound(n)%p(k)
         end if
         end if
         end if

      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_sound")

end subroutine da_check_max_iv_sound
subroutine da_get_innov_vector_sound (it,num_qcstat_conv, grid, ob, iv)

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

   real, allocatable :: model_u(:,:)  ! Model value u at ob location.
   real, allocatable :: model_v(:,:)  ! Model value v at ob location.
   real, allocatable :: model_t(:,:)  ! Model value t at ob location.
   real, allocatable :: model_q(:,:)  ! Model value q at ob location.

   real    :: speed, direction

   real    :: v_h(kts:kte)      ! Model value h at ob hor. location.
   real    :: v_p(kts:kte)      ! Model value p at ob hor. location.

   if (trace_use_dull) call da_trace_entry ("da_get_innov_vector_sound")

   allocate (model_u(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (model_v(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (model_t(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))
   allocate (model_q(iv%info(sound)%max_lev,iv%info(sound)%n1:iv%info(sound)%n2))

   model_u(:,:) = 0.0
   model_v(:,:) = 0.0
   model_t(:,:) = 0.0
   model_q(:,:) = 0.0

   if ( it > 1 ) then
     do n=iv%info(sound)%n1,iv%info(sound)%n2
        do k=1, iv%info(sound)%levels(n)
           if (iv%sound(n)%u(k)%qc == fails_error_max) iv%sound(n)%u(k)%qc = 0
           if (iv%sound(n)%v(k)%qc == fails_error_max) iv%sound(n)%v(k)%qc = 0
           if (iv%sound(n)%t(k)%qc == fails_error_max) iv%sound(n)%t(k)%qc = 0
           if (iv%sound(n)%q(k)%qc == fails_error_max) iv%sound(n)%q(k)%qc = 0
        end do
     end do
   end if

   do n=iv%info(sound)%n1, iv%info(sound)%n2
      if (iv%info(sound)%levels(n) < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      if (position_lev_dependant) then
         i(:)   = iv%info(sound)%i(:,n)
         j(:)   = iv%info(sound)%j(:,n)
         dx(:)  = iv%info(sound)%dx(:,n)
         dy(:)  = iv%info(sound)%dy(:,n)
         dxm(:) = iv%info(sound)%dxm(:,n)
         dym(:) = iv%info(sound)%dym(:,n)
         do k=kts,kte
            v_h(k) = dym(k)*(dxm(k)*grid%xb%h(i(k),j(k),k) + dx(k)*grid%xb%h(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%h(i(k),j(k)+1,k) + dx(k)*grid%xb%h(i(k)+1,j(k)+1,k))
            v_p(k) = dym(k)*(dxm(k)*grid%xb%p(i(k),j(k),k) + dx(k)*grid%xb%p(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%p(i(k),j(k)+1,k) + dx(k)*grid%xb%p(i(k)+1,j(k)+1,k))
         end do
      else
         i(1)   = iv%info(sound)%i(1,n)
         j(1)   = iv%info(sound)%j(1,n)
         dx(1)  = iv%info(sound)%dx(1,n)
         dy(1)  = iv%info(sound)%dy(1,n)
         dxm(1) = iv%info(sound)%dxm(1,n)
         dym(1) = iv%info(sound)%dym(1,n)

         v_h(kts:kte) = dym(1) * (dxm(1)*grid%xb%h(i(1),j(1),kts:kte)   + dx(1)*grid%xb%h(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%h(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%h(i(1)+1,j(1)+1,kts:kte))
         v_p(kts:kte) = dym(1) * (dxm(1)*grid%xb%p(i(1),j(1),kts:kte)   + dx(1)*grid%xb%p(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%p(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%p(i(1)+1,j(1)+1,kts:kte))
      end if

      do k=1, iv%info(sound)%levels(n)
         if (iv%sound(n)%p(k) > 1.0) then
            call da_to_zk (iv%sound(n)%p(k), v_p, v_interp_p, iv%info(sound)%zk(k,n))
         else if (iv%sound(n)%h(k) > 0.0) then
            call da_to_zk (iv%sound(n)%h(k), v_h, v_interp_h, iv%info(sound)%zk(k,n))
         end if

      end do

   end do

   call da_convert_zk (iv%info(sound))

   if (.not. anal_type_verify) then
      do n=iv%info(sound)%n1,iv%info(sound)%n2
         do k=1, iv%info(sound)%levels(n)
            if (iv%info(sound)%zk(k,n) < 0.0) then
               iv%sound(n)%u(k)%qc = missing_data
               iv%sound(n)%v(k)%qc = missing_data
               iv%sound(n)%t(k)%qc = missing_data
               iv%sound(n)%q(k)%qc = missing_data
            end if
         end do
      end do
   end if

   ! [1.2] Interpolate horizontally to ob:

   call da_interp_lin_3d (grid%xb%u, iv%info(sound), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(sound), model_v)
   call da_interp_lin_3d (grid%xb%t, iv%info(sound), model_t)
   call da_interp_lin_3d (grid%xb%q, iv%info(sound), model_q)

   do n=iv%info(sound)%n1, iv%info(sound)%n2
      !----------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !----------------------------------------------------------------------

      do k = 1, iv%info(sound)%levels(n)
         iv%sound(n)%u(k)%inv = 0.0
         iv%sound(n)%v(k)%inv = 0.0
         iv%sound(n)%t(k)%inv = 0.0
         iv%sound(n)%q(k)%inv = 0.0
         !-------------------------------------------------------------------
         ! [3.0] Interpolation:
         !-------------------------------------------------------------------

          if (wind_sd_sound) then
              call da_ffdduv_model (speed,direction,model_u(k,n), model_v(k,n), convert_uv2fd)

              if (ob%sound(n)%u(k) > missing_r .AND. iv%sound(n)%u(k)%qc >= obs_qc_pointer) then
                  iv%sound(n)%u(k)%inv = ob%sound(n)%u(k) - speed
              end if

              if (ob%sound(n)%v(k) > missing_r .AND. iv%sound(n)%v(k)%qc >= obs_qc_pointer) then
                  iv%sound(n)%v(k)%inv = ob%sound(n)%v(k) - direction
                  if (iv%sound(n)%v(k)%inv > 180.0 ) iv%sound(n)%v(k)%inv = iv%sound(n)%v(k)%inv - 360.0
                  if (iv%sound(n)%v(k)%inv < -180.0 ) iv%sound(n)%v(k)%inv = iv%sound(n)%v(k)%inv + 360.0
              end if
           else
             if (ob%sound(n)%u(k) > missing_r .AND. iv%sound(n)%u(k)%qc >= obs_qc_pointer) then
                 iv%sound(n)%u(k)%inv = ob%sound(n)%u(k) - model_u(k,n)
             end if

             if (ob%sound(n)%v(k) > missing_r .AND. iv%sound(n)%v(k)%qc >= obs_qc_pointer) then
                 iv%sound(n)%v(k)%inv = ob%sound(n)%v(k) - model_v(k,n)
             end if
          end if

         if (ob%sound(n)%t(k) > missing_r .AND. iv%sound(n)%t(k)%qc >= obs_qc_pointer) then
            iv%sound(n)%t(k)%inv = ob%sound(n)%t(k) - model_t(k,n)
         end if

         if (ob%sound(n)%q(k) > missing_r .AND. iv%sound(n)%q(k)%qc >= obs_qc_pointer) then
            iv%sound(n)%q(k)%inv = ob%sound(n)%q(k) - model_q(k,n)
         end if
      end do
   end do

   !----------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_sound (iv, it, num_qcstat_conv)

   if (check_buddy) call da_check_buddy_sound(iv, grid%dx, it)
!
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_q)
   
   if (trace_use_dull) call da_trace_exit ("da_get_innov_vector_sound")

end subroutine da_get_innov_vector_sound


subroutine da_calculate_grady_sound(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_sound")

   do n=1, iv%info(sound)%nlocal
      do k=1, iv%info(sound)%levels(n)
         if (iv%sound(n)%u(k)%qc < obs_qc_pointer) re%sound(n)%u(k) = 0.0
         if (iv%sound(n)%v(k)%qc < obs_qc_pointer) re%sound(n)%v(k) = 0.0
         if (iv%sound(n)%t(k)%qc < obs_qc_pointer) re%sound(n)%t(k) = 0.0
         if (iv%sound(n)%q(k)%qc < obs_qc_pointer) re%sound(n)%q(k) = 0.0

         jo_grad_y%sound(n)%u(k) = -re%sound(n)%u(k) / (iv%sound(n)%u(k)%error * iv%sound(n)%u(k)%error)
         jo_grad_y%sound(n)%v(k) = -re%sound(n)%v(k) / (iv%sound(n)%v(k)%error * iv%sound(n)%v(k)%error)
         jo_grad_y%sound(n)%t(k) = -re%sound(n)%t(k) / (iv%sound(n)%t(k)%error * iv%sound(n)%t(k)%error)
         jo_grad_y%sound(n)%q(k) = -re%sound(n)%q(k) / (iv%sound(n)%q(k)%error * iv%sound(n)%q(k)%error)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_sound")

end subroutine da_calculate_grady_sound


subroutine da_check_buddy_sound(iv, dx, it)

   !-----------------------------------------------------------------------
   ! Purpose: Buddy check for SOUND observations.
   !
   ! For SOUND, the binning procedure should be completed before going
   ! into the da_buddy_qc. 
   !
   !                       Yong-Run Guo, 10/10/2008
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration
   real,          intent(in)    :: dx

   integer :: k, n, bin, i, j, m_max, m_max_g, kk, nn, numobs, ierr
   real    :: dx_km, Press_mb

   integer, parameter               :: num_bins_p = 13
   real, parameter                  :: bin_width_p = 10000.0
   real, parameter                  :: bottom_pressure = 100000.0
   real   , dimension(0:num_bins_p) :: bin_start_p, pressure, bin_width
   integer, dimension(0:num_bins_p) :: num
!
   integer, allocatable, dimension(:,:) :: n_iv, k_iv

   integer,          allocatable, dimension(:,:,:) :: qc_flag_small
   real,             allocatable, dimension(:,:) :: xob, yob
   real,             allocatable, dimension(:,:,:) :: obs
   character(len=5), allocatable, dimension(:,:) :: station_id
   integer,          allocatable, dimension(:,:,:,:) :: qc_flag_small_g
   real,             allocatable, dimension(:,:,:) :: xob_g, yob_g
   real,             allocatable, dimension(:,:,:,:) :: obs_g
   integer,                       dimension(0:num_bins_p,0:num_procs-1) :: num_recv
!-----------------------------------------------------------------------------
   
   if (trace_use_dull) call da_trace_entry("da_check_buddy_sound")

   !--------------------------------------------------------------------------- 
   ! [1.0] Open diagnostic file:
   !---------------------------------------------------------------------------

   if (rootproc .and. check_buddy_print) then
      write (check_buddy_unit,'(/A)')  &
         '================================================================'
      write (unit = check_buddy_unit, fmt = '(A,i4,A,i4/)') &
            'SOUND BUDDY TEST QC:  no_buddies_qc=',no_buddies,&
            '  fails_buddy_check_qc=',fails_buddy_check
   end if

   !---------------------------------------------------------------------------
   ! [2.0] Bin the data vertically based on the obs p::
   !---------------------------------------------------------------------------

!   print*,'==> Sound Buddy check: num_bins_p = ',num_bins_p
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

!   print '(I3,2x,"start_p=",f10.1," mid-pressure=",f10.1," width=",f10.1)', &
!        (n, bin_start_p(n), pressure(n), bin_width(n), n=0, num_bins_p)
!
! 2.1 Get the maximum dimension for all the bins:
!
   num = 0
   do n = iv%info(sound)%n1,iv%info(sound)%n2
      do k =  1, iv%info(sound)%levels(n)
         if (iv%sound(n)%p(k) > missing_r) then
           do i = 0, num_bins_p - 1
              if (iv%sound(n)%p(k) <= bin_start_p(i) .and. &
                  iv%sound(n)%p(k) >  bin_start_p(i+1) ) then
                 bin = i
                 exit
              endif
           enddo
!           bin = int( (bottom_pressure - iv%sound(n)%p(k))/bin_width(n) ) + 1
           if (iv%sound(n)%p(k) > bottom_pressure) bin = 0
           if (iv%sound(n)%p(k) <=  bin_start_p(num_bins_p)) bin = num_bins_p
           num(bin) = num(bin) + 1
         endif
      enddo
   enddo
   m_max = maxval(num)
!   print *,(i,num(i),i=0,num_bins_p)
!   print *,"m_max=", m_max
!
! 2.2 Save the location indices (n,k) for each of bins:
!
!   print '("Sound n1=",i5,"  n2=",i5)',iv%info(sound)%n1,iv%info(sound)%n2
   allocate ( n_iv( 0: num_bins_p,1:m_max+10 ) )
   allocate ( k_iv( 0: num_bins_p,1:m_max+10 ) )

   num = 0
   do n = iv%info(sound)%n1,iv%info(sound)%n2
      do k =  1, iv%info(sound)%levels(n)
         if (iv%sound(n)%p(k) > missing_r) then
           do i = 0, num_bins_p - 1
              if (iv%sound(n)%p(k) <= bin_start_p(i) .and. &
                  iv%sound(n)%p(k) >  bin_start_p(i+1) ) then
                 bin = i
                 exit
              endif
           enddo
!           bin = int( (bottom_pressure - iv%sound(n)%p(k))/bin_width(n) ) + 1
           if (iv%sound(n)%p(k) > bottom_pressure) bin = 0
           if (iv%sound(n)%p(k) <=  bin_start_p(num_bins_p)) bin = num_bins_p

           num(bin) = num(bin) + 1
           n_iv(bin,num(bin)) = n
           k_iv(bin,num(bin)) = k

         endif
      enddo
   end do
!
! 2.3 Print out the binned results:
!
!   do i = 0, num_bins_p
!      print '("bin:",I2,"  start_p=",f8.1," num=",i5)', &
!                      i, bin_start_p(i), num(i)
!      do j = 1, num(i)
!         n = n_iv(i,j)
!         k = k_iv(i,j)
!         print '("j, n, k:",3i5,2x,"p=",f10.1)', &
!                  j, n, k, iv%sound(n)%p(k)
!      enddo
!   enddo
   !---------------------------------------------------------------------------
   ! [3.0] Buddy check for each of the pressure-bins::
   !---------------------------------------------------------------------------

   call mpi_allgather(num, num_bins_p+1, mpi_integer, num_recv, num_bins_p+1, mpi_integer, comm, ierr)
   call mpi_allreduce(m_max, m_max_g, 1, mpi_integer, mpi_max, comm, ierr)
   m_max = m_max_g
   allocate(xob_g(1:m_max,0:num_bins_p,0:num_procs-1))
   allocate(yob_g(1:m_max,0:num_bins_p,0:num_procs-1))
   allocate(obs_g(1:m_max,4,0:num_bins_p,0:num_procs-1))
   allocate(qc_flag_small_g(1:m_max,4,0:num_bins_p,0:num_procs-1))

   allocate(xob(1:m_max,0:num_bins_p))
   allocate(yob(1:m_max,0:num_bins_p))
   allocate(obs(1:m_max,4,0:num_bins_p))
   allocate(qc_flag_small(1:m_max,4,0:num_bins_p))
   allocate(station_id   (1:m_max,0:num_bins_p))

   obs = 0.0
   qc_flag_small = missing
   obs_g = 0.0
   qc_flag_small_g = missing

   do i = 0, num_bins_p
      
      numobs = num(i)
 
      do n = 1, numobs
         nn = n_iv(i,n)
         kk = k_iv(i,n)
 
         station_id(n,i)          = iv%info(sound)%id(nn)
         xob(n,i)                 = iv%info(sound)%x(1,nn)
         yob(n,i)                 = iv%info(sound)%y(1,nn)

         obs(n,1,i)               = iv%sound(nn)%u(kk)%inv
         qc_flag_small(n,1,i)     = iv%sound(nn)%u(kk)%qc
         obs(n,2,i)               = iv%sound(nn)%v(kk)%inv
         qc_flag_small(n,2,i)     = iv%sound(nn)%v(kk)%qc
         obs(n,3,i)               = iv%sound(nn)%t(kk)%inv
         qc_flag_small(n,3,i)     = iv%sound(nn)%t(kk)%qc
         obs(n,4,i)               = iv%sound(nn)%q(kk)%inv
         qc_flag_small(n,4,i)     = iv%sound(nn)%q(kk)%qc
      enddo

   enddo

   call mpi_allgather (xob, m_max*(num_bins_p+1), mpi_real8, xob_g, m_max*(num_bins_p+1), mpi_real8, comm, ierr)
   call mpi_allgather (yob, m_max*(num_bins_p+1), mpi_real8, yob_g, m_max*(num_bins_p+1), mpi_real8, comm, ierr)
   call mpi_allgather (obs, 4*m_max*(num_bins_p+1), mpi_real8, obs_g, 4*m_max*(num_bins_p+1), mpi_real8, comm, ierr)
   call mpi_allgather (qc_flag_small, 4*m_max*(num_bins_p+1), mpi_integer, qc_flag_small_g, 4*m_max*(num_bins_p+1), &
                        mpi_integer, comm, ierr)

   do i = 0, num_bins_p

     if ( num(i) <= 1 ) cycle


! 3.1 Get the Station locations:
   
     ! Pressure level:
     Press_mb = pressure(i) / 100.0
     numobs = num(i)

     if (rootproc .and. check_buddy_print) then
         write (check_buddy_unit,'(5X,A)')  &
         '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write (unit = check_buddy_unit, fmt = '(5X,A,I3,2X,A,I6)') &
             'BIN =', i, 'NUMOBS =', numobs
     end if


! 3.2 U-component buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'UU      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_uv=', max_buddy_uv 

     call da_buddy_qc (numobs, m_max, station_id(:,i), xob(:,i), yob(:,i), obs(:,1,i), qc_flag_small(:,1,i), &
                       'UU      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_uv , check_buddy_unit, check_buddy_print,  xob_g(:,i,:), &
                       yob_g(:,i,:), obs_g(:,1,i,:), qc_flag_small_g(:,1,i,:),num_recv(i,:))

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        kk = k_iv(i,n)
        iv%sound(nn)%u(kk)%qc = qc_flag_small(n,1,i)
     enddo

! 3.2 V-component buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'VV      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_uv=', max_buddy_uv
 
     call da_buddy_qc (numobs, m_max, station_id(:,i), xob(:,i), yob(:,i), obs(:,2,i), qc_flag_small(:,2,i), &
                       'UU      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_uv , check_buddy_unit, check_buddy_print,  xob_g(:,i,:), &
                       yob_g(:,i,:), obs_g(:,2,i,:), qc_flag_small_g(:,2,i,:),num_recv(i,:))

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        kk = k_iv(i,n)
        iv%sound(nn)%v(kk)%qc = qc_flag_small(n,2,i)
     enddo

! 3.3 Temperature buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'TT      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_t=', max_buddy_t 

     call da_buddy_qc (numobs, m_max, station_id(:,i), xob(:,i), yob(:,i), obs(:,3,i), qc_flag_small(:,3,i), &
                       'UU      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_uv , check_buddy_unit, check_buddy_print,  xob_g(:,i,:), &
                       yob_g(:,i,:), obs_g(:,3,i,:), qc_flag_small_g(:,3,i,:),num_recv(i,:))

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        kk = k_iv(i,n)
        iv%sound(nn)%t(kk)%qc = qc_flag_small(n,3,i)
     enddo

! 3.3 Specific humidity buddy check:

     if (rootproc .and. check_buddy_print) &
         write (check_buddy_unit,'(8X,A,A,f10.1,3(A,f6.1))')  &
                'QQ      ', ' Pressure(mb)=',Press_mb, ' ds(km)=',dx_km,&
                '  buddy_weight=', buddy_weight , &
                '  max_buddy_rh=', max_buddy_rh 

     call da_buddy_qc (numobs, m_max, station_id(:,i), xob(:,i), yob(:,i), obs(:,4,i), qc_flag_small(:,4,i), &
                       'UU      ', Press_mb, dx_km, buddy_weight , &
                       max_buddy_uv , check_buddy_unit, check_buddy_print,  xob_g(:,i,:), &
                       yob_g(:,i,:), obs_g(:,4,i,:), qc_flag_small_g(:,4,i,:),num_recv(i,:))

   !  Put the qc_flag back into the permanent space.
   
     do n = 1, numobs
        nn = n_iv(i,n)
        kk = k_iv(i,n)
        iv%sound(nn)%q(kk)%qc = qc_flag_small(n,4,i)
     enddo

! 3.4 Deallocate arrays

   enddo

   deallocate(xob_g)
   deallocate(yob_g)
   deallocate(obs_g)
   deallocate(qc_flag_small_g)
   deallocate(xob)
   deallocate(yob)
   deallocate(obs)
   deallocate(qc_flag_small)
   deallocate(station_id   )

   deallocate ( n_iv )
   deallocate ( k_iv )

   if (trace_use_dull) call da_trace_exit("da_check_buddy_sound")

end subroutine da_check_buddy_sound

subroutine da_ao_stats_sonde_sfc (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type  (y_type), intent (in)    :: re            ! A - O
   type(y_type),   intent (in)    :: ob            ! Observation structure.

   type (stats_sonde_sfc_type)      :: stats
   integer                          :: nu, nv, nt, np, nq
   integer                          :: n
   real                             :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_sonde_sfc")

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
   stats%average = residual_sonde_sfc1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(sound)%nlocal
      if (iv%info(sonde_sfc)%proc_domain(1,n)) then

         u_inc = re%sonde_sfc(n)%u
         v_inc = re%sonde_sfc(n)%v
         u_obs = ob%sonde_sfc(n)%u
         v_obs = ob%sonde_sfc(n)%v

         if (.not. wind_sd_sound .and. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%sonde_sfc(n)%u%qc, iv%sonde_sfc(n)%v%qc, convert_uv2fd)
         if (wind_sd_sound .and. .not. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%sonde_sfc(n)%u%qc, iv%sonde_sfc(n)%v%qc, convert_fd2uv)

         call da_stats_calculate (n, 0, iv%sonde_sfc(n)%u%qc,  & 
            u_inc, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%sonde_sfc(n)%v%qc,  & 
            v_inc, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate (n, 0, iv%sonde_sfc(n)%t%qc,  & 
            re%sonde_sfc(n)%t, nt, & 
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)  
         call da_stats_calculate (n, 0, iv%sonde_sfc(n)%p%qc,  & 
            re%sonde_sfc(n)%p, np, & 
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate (n, 0, iv%sonde_sfc(n)%q%qc,  & 
            re%sonde_sfc(n)%q, nq, & 
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
   iv%nstats(sonde_sfc) = nu + nv + nt + np + nq

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for sonde_sfc'
         call da_print_stats_sonde_sfc(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_sonde_sfc")

end subroutine da_ao_stats_sonde_sfc


subroutine da_jo_and_grady_sonde_sfc( iv, re, jo, jo_grad_y)

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

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_sonde_sfc")

   jo % sonde_sfc_u = 0.0
   jo % sonde_sfc_v = 0.0
   jo % sonde_sfc_t = 0.0
   jo % sonde_sfc_p = 0.0
   jo % sonde_sfc_q = 0.0

   if (test_dm_exact) then
      if (iv%info(sound)%ntotal == 0) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_sonde_sfc")
         return
      end if
   else
      if (iv%info(sound)%nlocal < 1) then
         if (trace_use_dull) call da_trace_exit("da_jo_and_grady_sonde_sfc")
         return
      end if
   end if

   do n=1, iv%info(sound)%nlocal
      jo_grad_y%sonde_sfc(n)%u = -re%sonde_sfc(n)%u / (iv%sonde_sfc(n)%u%error * iv%sonde_sfc(n)%u%error)
      jo_grad_y%sonde_sfc(n)%v = -re%sonde_sfc(n)%v / (iv%sonde_sfc(n)%v%error * iv%sonde_sfc(n)%v%error)
      jo_grad_y%sonde_sfc(n)%t = -re%sonde_sfc(n)%t / (iv%sonde_sfc(n)%t%error * iv%sonde_sfc(n)%t%error)
      jo_grad_y%sonde_sfc(n)%p = -re%sonde_sfc(n)%p / (iv%sonde_sfc(n)%p%error * iv%sonde_sfc(n)%p%error)
      jo_grad_y%sonde_sfc(n)%q = -re%sonde_sfc(n)%q / (iv%sonde_sfc(n)%q%error * iv%sonde_sfc(n)%q%error)
   end do

   ! Bitwise-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not
   ! testing.  This will always be .false. for a serial run.
   if (test_dm_exact) then
      ! collect all obs in serial order and allocate global objects
      call da_to_global_sonde_sfc( iv, re, jo_grad_y, iv_glob, re_glob, jo_grad_y_glob)
      ! perform remaining computations
      call da_jo_sonde_sfc_uvtq( iv_glob, re_glob, jo_grad_y_glob, jo)
      ! free global objects
      call da_deallocate_global_sonde_sfc( iv_glob, re_glob, jo_grad_y_glob)
   else
      ! perform remaining computations
      call da_jo_sonde_sfc_uvtq( iv, re, jo_grad_y, jo)
   end if

   jo % sonde_sfc_u = 0.5 * jo % sonde_sfc_u
   jo % sonde_sfc_v = 0.5 * jo % sonde_sfc_v
   jo % sonde_sfc_t = 0.5 * jo % sonde_sfc_t
   jo % sonde_sfc_p = 0.5 * jo % sonde_sfc_p
   jo % sonde_sfc_q = 0.5 * jo % sonde_sfc_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_sonde_sfc")

end subroutine da_jo_and_grady_sonde_sfc


subroutine da_jo_sonde_sfc_uvtq(iv, re, jo_grad_y, jo)

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

   if (trace_use_dull) call da_trace_entry("da_jo_sonde_sfc_uvtq")

   do n=1, iv%info(sound)%nlocal
      if (iv%info(sonde_sfc)%proc_domain(1,n)) then
        jo % sonde_sfc_u = jo % sonde_sfc_u - re%sonde_sfc(n)%u * jo_grad_y%sonde_sfc(n)%u
        jo % sonde_sfc_v = jo % sonde_sfc_v - re%sonde_sfc(n)%v * jo_grad_y%sonde_sfc(n)%v
        jo % sonde_sfc_t = jo % sonde_sfc_t - re%sonde_sfc(n)%t * jo_grad_y%sonde_sfc(n)%t
        jo % sonde_sfc_p = jo % sonde_sfc_p - re%sonde_sfc(n)%p * jo_grad_y%sonde_sfc(n)%p
        jo % sonde_sfc_q = jo % sonde_sfc_q - re%sonde_sfc(n)%q * jo_grad_y%sonde_sfc(n)%q
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_jo_sonde_sfc_uvtq")

end subroutine da_jo_sonde_sfc_uvtq


subroutine da_residual_sonde_sfc(iv, y, re,np_missing, np_bad_data,np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for sonde surface obs
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

   if (trace_use_dull) call da_trace_entry("da_residual_sonde_sfc")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % p % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(sound)%nlocal
      np_available = np_available + 5

      re%sonde_sfc(n)%u = da_residual(n, 0, y%sonde_sfc(n)%u, iv%sonde_sfc(n)%u, n_obs_bad % u)
      re%sonde_sfc(n)%v = da_residual(n, 0, y%sonde_sfc(n)%v, iv%sonde_sfc(n)%v, n_obs_bad % v)
      re%sonde_sfc(n)%t = da_residual(n, 0, y%sonde_sfc(n)%t, iv%sonde_sfc(n)%t, n_obs_bad % t)
      re%sonde_sfc(n)%p = da_residual(n, 0, y%sonde_sfc(n)%p, iv%sonde_sfc(n)%p, n_obs_bad % p)
      re%sonde_sfc(n)%q = da_residual(n, 0, y%sonde_sfc(n)%q, iv%sonde_sfc(n)%q, n_obs_bad % q)
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

   if (trace_use_dull) call da_trace_exit("da_residual_sonde_sfc")

end subroutine da_residual_sonde_sfc


subroutine da_oi_stats_sonde_sfc (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_sonde_sfc_type) :: stats
   integer                     :: nu, nv, nt, np, nq
   integer                     :: n
   real                        :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_sonde_sfc")

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

   stats%average = residual_sonde_sfc1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(sound)%nlocal
      if (iv%info(sonde_sfc)%proc_domain(1,n)) then

        u_inv = iv%sonde_sfc(n)%u%inv
        v_inv = iv%sonde_sfc(n)%v%inv
        u_obs = ob%sonde_sfc(n)%u
        v_obs = ob%sonde_sfc(n)%v

        if (.not. wind_sd_sound .and. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%sonde_sfc(n)%u%qc, iv%sonde_sfc(n)%v%qc, convert_uv2fd)
        if (wind_sd_sound .and. .not. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%sonde_sfc(n)%u%qc, iv%sonde_sfc(n)%v%qc, convert_fd2uv)

         call da_stats_calculate(iv%info(sonde_sfc)%obs_global_index(n), &
            0, iv%sonde_sfc(n)%u%qc, &
            u_inv, nu, &
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate(iv%info(sonde_sfc)%obs_global_index(n), &
            0, iv%sonde_sfc(n)%v%qc, &
            v_inv, nv, &
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate(iv%info(sonde_sfc)%obs_global_index(n), &
            0, iv%sonde_sfc(n)%t%qc, &
            iv%sonde_sfc(n)%t%inv, nt, &
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate(iv%info(sonde_sfc)%obs_global_index(n), &
            0, iv%sonde_sfc(n)%p%qc, &
            iv%sonde_sfc(n)%p%inv, np, &
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)

         call da_stats_calculate(iv%info(sonde_sfc)%obs_global_index(n), &
            0, iv%sonde_sfc(n)%q%qc, &
            iv%sonde_sfc(n)%q%inv, nq, &
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)
      end if    ! end if (iv%info(sonde_sfc)%proc_domain(1,n))
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for sonde_sfc'
         call da_print_stats_sonde_sfc(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_sonde_sfc")

 end subroutine da_oi_stats_sonde_sfc


subroutine da_print_stats_sonde_sfc(stats_unit, nu, nv, nt, np, nq, sonde_sfc)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                     intent(in)    :: stats_unit
   integer,                     intent(inout) :: nu, nv, nt, np, nq
   type (stats_sonde_sfc_type), intent(in)    :: sonde_sfc

   if (trace_use_dull) call da_trace_entry("da_print_stats_sonde_sfc")

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
      ' Minimum(n,k): ', sonde_sfc%minimum%u, sonde_sfc%minimum%v, &
      sonde_sfc%minimum%t, sonde_sfc%minimum%p, sonde_sfc%minimum%q, &
      ' Maximum(n,k): ', sonde_sfc%maximum%u, sonde_sfc%maximum%v, &
      sonde_sfc%maximum%t, &
                        sonde_sfc%maximum%p, sonde_sfc%maximum%q
   write(unit=stats_unit, fmt='((a,4(f12.4,10x),e12.4,10x))') &
      ' Average     : ', sonde_sfc%average%u/real(nu), &
      sonde_sfc%average%v/real(nv), &
      sonde_sfc%average%t/real(nt), sonde_sfc%average%p/real(np), &
      sonde_sfc%average%q/real(nq), &
      '    RMSE     : ', sqrt(sonde_sfc%rms_err%u/real(nu)), &
      sqrt(sonde_sfc%rms_err%v/real(nv)), &
      sqrt(sonde_sfc%rms_err%t/real(nt)), &
      sqrt(sonde_sfc%rms_err%p/real(np)), &
      sqrt(sonde_sfc%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_sonde_sfc")

end subroutine da_print_stats_sonde_sfc


subroutine da_transform_xtoy_sonde_sfc (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_sonde_sfc")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_v(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_t(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_q(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_psfc(iv%info(sound)%n1:iv%info(sound)%n2))

      allocate (ub(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (vb(1,iv%info(sound)%n1:iv%info(sound)%n2))

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xa%u, iv%info(sonde_sfc), model_u)
      call da_interp_lin_3d (grid%xa%v, iv%info(sonde_sfc), model_v)
      call da_interp_lin_3d (grid%xa%t, iv%info(sonde_sfc), model_t)
      call da_interp_lin_3d (grid%xa%q, iv%info(sonde_sfc), model_q)

      call da_interp_lin_2d (grid%xa%psfc, iv%info(sonde_sfc), 1, model_psfc)

      call da_interp_lin_3d (grid%xb%u, iv%info(sonde_sfc), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(sonde_sfc), vb)

      do n=iv%info(sound)%n1,iv%info(sound)%n2
         if(wind_sd_sound)then
            call da_uv_to_sd_lin(y%sonde_sfc(n)%u,y%sonde_sfc(n)%v,model_u(1,n),model_v(1,n),ub(1,n),vb(1,n))
         else
            y%sonde_sfc(n)%u = model_u(1,n)
            y%sonde_sfc(n)%v = model_v(1,n)
         end if
         y%sonde_sfc(n)%t = model_t(1,n)
         y%sonde_sfc(n)%q = model_q(1,n)
         y%sonde_sfc(n)%p = model_psfc(n)
      end do

      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc(grid, iv, sonde_sfc, iv%sonde_sfc(:), y%sonde_sfc(:))
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_sonde_sfc")

end subroutine da_transform_xtoy_sonde_sfc


subroutine da_transform_xtoy_sonde_sfc_adj(grid, iv, jo_grad_y, jo_grad_x)

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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_sonde_sfc_adj")

   if (sfc_assi_options == sfc_assi_options_1) then

      allocate (model_u(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_v(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_t(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_q(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (model_psfc(iv%info(sound)%n1:iv%info(sound)%n2))

      allocate (ub(1,iv%info(sound)%n1:iv%info(sound)%n2))
      allocate (vb(1,iv%info(sound)%n1:iv%info(sound)%n2))

      call da_interp_lin_3d (grid%xb%u, iv%info(sonde_sfc), ub)
      call da_interp_lin_3d (grid%xb%v, iv%info(sonde_sfc), vb)

      ! [1.2] Interpolate horizontally:
      do n=iv%info(sound)%n1,iv%info(sound)%n2
         if(wind_sd_sound)then
            call da_uv_to_sd_adj(jo_grad_y%sonde_sfc(n)%u, &
                                jo_grad_y%sonde_sfc(n)%v, model_u(1,n), model_v(1,n), ub(1,n), vb(1,n))
         else
            model_u(1,n)  = jo_grad_y%sonde_sfc(n)%u
            model_v(1,n)  = jo_grad_y%sonde_sfc(n)%v
         end if
         model_t(1,n)  = jo_grad_y%sonde_sfc(n)%t
         model_q(1,n)  = jo_grad_y%sonde_sfc(n)%q
         model_psfc(n) = jo_grad_y%sonde_sfc(n)%p
      end do

      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(sonde_sfc), model_u)
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(sonde_sfc), model_v)
      call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(sonde_sfc), model_t)
      call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(sonde_sfc), model_q)

      call da_interp_lin_2d_adj (jo_grad_x%psfc, iv%info(sonde_sfc), 1, model_psfc)
      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
      deallocate (ub)
      deallocate (vb)

   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc_adj(grid,iv,sonde_sfc,iv%sonde_sfc(:), jo_grad_y%sonde_sfc(:),jo_grad_x)
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_sonde_sfc_adj")

end subroutine da_transform_xtoy_sonde_sfc_adj


subroutine da_get_innov_vector_sonde_sfc( it, num_qcstat_conv,grid, ob, iv)

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
   
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_sonde_sfc")

   allocate (model_u(1,iv%info(sonde_sfc)%n1:iv%info(sonde_sfc)%n2))
   allocate (model_v(1,iv%info(sonde_sfc)%n1:iv%info(sonde_sfc)%n2))
   allocate (model_t(1,iv%info(sonde_sfc)%n1:iv%info(sonde_sfc)%n2))
   allocate (model_p(1,iv%info(sonde_sfc)%n1:iv%info(sonde_sfc)%n2))
   allocate (model_q(1,iv%info(sonde_sfc)%n1:iv%info(sonde_sfc)%n2))
   allocate (model_hsm(1,iv%info(sonde_sfc)%n1:iv%info(sonde_sfc)%n2))

   if ( it > 1 ) then
      do n=iv%info(sonde_sfc)%n1,iv%info(sonde_sfc)%n2
         if (iv%sonde_sfc(n)%u%qc == fails_error_max) iv%sonde_sfc(n)%u%qc = 0
         if (iv%sonde_sfc(n)%v%qc == fails_error_max) iv%sonde_sfc(n)%v%qc = 0
         if (iv%sonde_sfc(n)%t%qc == fails_error_max) iv%sonde_sfc(n)%t%qc = 0
         if (iv%sonde_sfc(n)%p%qc == fails_error_max) iv%sonde_sfc(n)%p%qc = 0
         if (iv%sonde_sfc(n)%q%qc == fails_error_max) iv%sonde_sfc(n)%q%qc = 0
      end do
   end if

   if (sfc_assi_options == sfc_assi_options_1) then
      do n=iv%info(sonde_sfc)%n1,iv%info(sonde_sfc)%n2
         ! [1.1] Get horizontal interpolation weights:

         i   = iv%info(sonde_sfc)%i(1,n)
         j   = iv%info(sonde_sfc)%j(1,n)
         dx  = iv%info(sonde_sfc)%dx(1,n)
         dy  = iv%info(sonde_sfc)%dy(1,n)
         dxm = iv%info(sonde_sfc)%dxm(1,n)
         dym = iv%info(sonde_sfc)%dym(1,n)

         ! Surface correction

         iv%sonde_sfc(n)%p%inv = ob%sonde_sfc(n)%p
         iv%sonde_sfc(n)%t%inv = ob%sonde_sfc(n)%t
         iv%sonde_sfc(n)%q%inv = ob%sonde_sfc(n)%q
         iv%sonde_sfc(n)%u%inv = ob%sonde_sfc(n)%u
         iv%sonde_sfc(n)%v%inv = ob%sonde_sfc(n)%v

         if (iv % sonde_sfc(n) % h > missing_r) then
            do k=kts,kte
               v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                  + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
            end do

            hd = v_h(kts) - iv % sonde_sfc(n) % h

            if (abs(hd) <= Max_StHeight_Diff .or. anal_type_verify ) then
               if (iv % sonde_sfc(n) % h < v_h(kts)) then
                  iv%info(sonde_sfc)%zk(:,n) = 1.0+1.0e-6
                  call da_obs_sfc_correction(iv%info(sonde_sfc), iv%sonde_sfc(n), n, grid%xb)

               else
                  call da_to_zk(iv % sonde_sfc(n) % h, v_h, v_interp_h, iv%info(sonde_sfc)%zk(1,n))
               end if
            end if
         else if (ob % sonde_sfc(n) % p > 1.0) then
            do k=kts,kte
               v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                      + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
            end do

            call da_to_zk(ob % sonde_sfc(n) % p, v_p, v_interp_p, iv%info(sonde_sfc)%zk(1,n))

            if (iv%info(sonde_sfc)%zk(1,n) < 0.0 .and. .not. anal_type_verify) then
               iv % sonde_sfc(n) % p % inv = missing_r
               iv % sonde_sfc(n) % p % qc  = missing_data
               iv%info(sonde_sfc)%zk(:,n) = 1.0+1.0e-6
            end if
         end if
      end do

      call da_convert_zk (iv%info(sonde_sfc))

      if (.not.anal_type_verify ) then
         do n=iv%info(sonde_sfc)%n1,iv%info(sonde_sfc)%n2
            if (iv%info(sonde_sfc)%zk(1,n) < 0.0) then
               iv % sonde_sfc(n) % u % qc = missing_data
               iv % sonde_sfc(n) % v % qc = missing_data
               iv % sonde_sfc(n) % t % qc = missing_data
               iv % sonde_sfc(n) % q % qc = missing_data
               iv % sonde_sfc(n) % p % qc = missing_data
            end if
         end do
      end if

      ! [1.2] Interpolate horizontally:
      call da_interp_lin_3d (grid%xb%u, iv%info(sonde_sfc), model_u)
      call da_interp_lin_3d (grid%xb%v, iv%info(sonde_sfc), model_v)
      call da_interp_lin_3d (grid%xb%t, iv%info(sonde_sfc), model_t)
      call da_interp_lin_3d (grid%xb%q, iv%info(sonde_sfc), model_q)
      call da_interp_lin_3d (grid%xb%p, iv%info(sonde_sfc), model_p)

   else if (sfc_assi_options == sfc_assi_options_2) then
      ! 1.2.1 Surface assimilation approach 2(10-m u, v, 2-m t, q, and sfc_p)

      call da_interp_lin_2d (grid%xb%u10,  iv%info(sonde_sfc), 1,model_u)
      call da_interp_lin_2d (grid%xb%v10,  iv%info(sonde_sfc), 1,model_v)
      call da_interp_lin_2d (grid%xb%t2,   iv%info(sonde_sfc), 1,model_t)
      call da_interp_lin_2d (grid%xb%q2,   iv%info(sonde_sfc), 1,model_q)
      call da_interp_lin_2d (grid%xb%psfc, iv%info(sonde_sfc), 1,model_p)

      do n=iv%info(sonde_sfc)%n1,iv%info(sonde_sfc)%n2

         iv%sonde_sfc(n)%p%inv = ob%sonde_sfc(n)%p
         iv%sonde_sfc(n)%t%inv = ob%sonde_sfc(n)%t
         iv%sonde_sfc(n)%q%inv = ob%sonde_sfc(n)%q
         iv%sonde_sfc(n)%u%inv = ob%sonde_sfc(n)%u
         iv%sonde_sfc(n)%v%inv = ob%sonde_sfc(n)%v

         if (iv%sonde_sfc(n)%p%qc >= 0) then
            ! model surface p, t, q, h at observed site:

            call da_interp_lin_2d_partial (grid%xb%terr, iv%info(sonde_sfc), 1, n, n, model_hsm(:,n))

            ho = iv%sonde_sfc(n)%h
            to = -888888.0
            qo = -888888.0

            if (iv%sonde_sfc(n)%t%qc >= 0 .and. iv%sonde_sfc(n)%q%qc >= 0) then
               to = ob%sonde_sfc(n)%t
               qo = ob%sonde_sfc(n)%q
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho, to, qo)
            else if (iv%sonde_sfc(n)%t%qc >= 0 .and. iv%sonde_sfc(n)%q%qc < 0) then
               to = ob%sonde_sfc(n)%t
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho,to)
            else
               call da_sfc_pre(psfcm, model_p(1,n), model_t(1,n), model_q(1,n), model_hsm(1,n), ho)
            end if

            ! Pressure at the observed height:
            model_p(1,n) = psfcm
         end if
      end do
   end if

   do n=iv%info(sonde_sfc)%n1,iv%info(sonde_sfc)%n2
      !-----------------------------------------------------------------------
      ! [3.0] Fast interpolation:
      !-----------------------------------------------------------------------
      if(wind_sd_sound)then
         call da_ffdduv_model (speed,direction,model_u(1,n), model_v(1,n), convert_uv2fd)

         if (ob%sonde_sfc(n)%u > missing_r .AND. iv%sonde_sfc(n)%u%qc >= obs_qc_pointer) then
             iv%sonde_sfc(n)%u%inv = iv%sonde_sfc(n)%u%inv - speed
         else
             iv % sonde_sfc(n) % u % inv = 0.0
         end if

         if (ob%sonde_sfc(n)%v > missing_r .AND. iv%sonde_sfc(n)%v%qc >= obs_qc_pointer) then
             iv%sonde_sfc(n)%v%inv = iv%sonde_sfc(n)%v%inv - direction
             if (iv%sonde_sfc(n)%v%inv > 180.0 ) iv%sonde_sfc(n)%v%inv = iv%sonde_sfc(n)%v%inv - 360.0
             if (iv%sonde_sfc(n)%v%inv < -180.0 ) iv%sonde_sfc(n)%v%inv = iv%sonde_sfc(n)%v%inv + 360.0
         else
             iv % sonde_sfc(n) % v % inv = 0.0
         end if
      else
         if (ob % sonde_sfc(n) % u > missing_r .AND. iv % sonde_sfc(n) % u % qc >= obs_qc_pointer) then
             iv % sonde_sfc(n) % u % inv = iv%sonde_sfc(n)%u%inv - model_u(1,n)
         else
             iv % sonde_sfc(n) % u % inv = 0.0
         end if

         if (ob % sonde_sfc(n) % v > missing_r .AND. iv % sonde_sfc(n) % v % qc >= obs_qc_pointer) then
             iv % sonde_sfc(n) % v % inv = iv%sonde_sfc(n)%v%inv - model_v(1,n)
         else
             iv % sonde_sfc(n) % v % inv = 0.0
         end if
      end if

      !if (ob % sonde_sfc(n) % p > 0.0 .AND. iv % sonde_sfc(n) % p % qc >= obs_qc_pointer) then
      if ( iv % sonde_sfc(n) % p % qc >= obs_qc_pointer ) then
         iv % sonde_sfc(n) % p % inv = iv%sonde_sfc(n)%p%inv - model_p(1,n)
      else
         iv % sonde_sfc(n) % p % inv = 0.0
      end if

      if (ob % sonde_sfc(n) % t > 0.0 .AND. iv % sonde_sfc(n) % t % qc >= obs_qc_pointer) then
         iv % sonde_sfc(n) % t % inv = iv%sonde_sfc(n)%t%inv - model_t(1,n)
      else
         iv % sonde_sfc(n) % t % inv = 0.0
      end if

      if (ob % sonde_sfc(n) % q > 0.0 .AND. iv % sonde_sfc(n) % q % qc >= obs_qc_pointer) then
         iv % sonde_sfc(n) % q % inv = iv%sonde_sfc(n)%q%inv - model_q(1,n)
      else
         iv % sonde_sfc(n) % q % inv = 0.0
      end if
   end do

   !-----------------------------------------------------------------------
   !     [5.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_sonde_sfc(iv,ob, it, num_qcstat_conv)
 
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_t)
   deallocate (model_p)
   deallocate (model_q)
   deallocate (model_hsm)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_sonde_sfc")

end subroutine da_get_innov_vector_sonde_sfc


subroutine da_check_max_iv_sonde_sfc(iv,ob, it, num_qcstat_conv)                                 

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

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_sonde_sfc")


   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(sonde_sfc)%n1,iv%info(sonde_sfc)%n2
         if(.not. qc_rej_both)then
            if(wind_sd_sound)then
               failed=.false.
               if( iv%sonde_sfc(n)%u%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%u, max_error_spd,failed)
                   if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
                       num_qcstat_conv(1,sonde_sfc,1,1) = num_qcstat_conv(1,sonde_sfc,1,1) + 1
                       if(failed) then
                          num_qcstat_conv(2,sonde_sfc,1,1) = num_qcstat_conv(2,sonde_sfc,1,1) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'sonde_sfc',ob_vars(1),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
                       end if
                   end if
                end if

                failed=.false.
                if( iv%sonde_sfc(n)%v%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%v, max_error_dir,failed)
                    if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
                        num_qcstat_conv(1,sonde_sfc,2,1) = num_qcstat_conv(1,sonde_sfc,2,1) + 1
                        if(failed)then
                           num_qcstat_conv(2,sonde_sfc,2,1) = num_qcstat_conv(2,sonde_sfc,2,1) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'sonde_sfc',ob_vars(2),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
                        end if
                    end if
                 end if
             else
                 failed=.false.
                 if( iv%sonde_sfc(n)%u%qc >= obs_qc_pointer ) then
                     call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%u, max_error_uv,failed)
                     if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
                         num_qcstat_conv(1,sonde_sfc,1,1) = num_qcstat_conv(1,sonde_sfc,1,1) + 1
                         if(failed) then
                            num_qcstat_conv(2,sonde_sfc,1,1) = num_qcstat_conv(2,sonde_sfc,1,1) + 1
                            write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                            'sonde_sfc',ob_vars(1),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
                         end if
                     end if
                  end if

                  failed=.false.
                  if( iv%sonde_sfc(n)%v%qc >= obs_qc_pointer ) then
                      call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%v, max_error_uv,failed)
                      if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
                          num_qcstat_conv(1,sonde_sfc,2,1) = num_qcstat_conv(1,sonde_sfc,2,1) + 1
                          if(failed)then
                             num_qcstat_conv(2,sonde_sfc,2,1) = num_qcstat_conv(2,sonde_sfc,2,1) + 1
                             write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                             'sonde_sfc',ob_vars(2),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
                          end if
                      end if
                  end if
              end if

              if(wind_sd_sound)then
                 if(iv%sonde_sfc(n)%u%qc == fails_error_max .or. abs(iv%sonde_sfc(n)%u%inv) >= max_omb_spd) then
                    iv%sonde_sfc(n)%u%qc = fails_error_max
                    iv%sonde_sfc(n)%u%inv = 0.0
                 endif
                 if(iv%sonde_sfc(n)%v%qc == fails_error_max .or. abs(iv%sonde_sfc(n)%v%inv) >= max_omb_dir) then
                    iv%sonde_sfc(n)%v%qc = fails_error_max
                    iv%sonde_sfc(n)%v%inv = 0.0
                 endif
              endif
           else
              failed1=.false.
              failed2=.false.

              if( iv%sonde_sfc(n)%v%qc >= obs_qc_pointer .or. iv%sonde_sfc(n)%u%qc >= obs_qc_pointer )  then
                  if(wind_sd_sound)then
                     call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%u, max_error_spd,failed1)
                     call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%v, max_error_dir,failed2)
                  else
                     call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%u, max_error_uv,failed1)
                     call da_max_error_qc (it,iv%info(sonde_sfc), n, iv%sonde_sfc(n)%v, max_error_uv,failed2)
                  endif
              endif

              if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
                  num_qcstat_conv(1,sonde_sfc,1,1) = num_qcstat_conv(1,sonde_sfc,1,1) + 1
                  num_qcstat_conv(1,sonde_sfc,2,1) = num_qcstat_conv(1,sonde_sfc,2,1) + 1
                  if(failed1 .or. failed2) then
                     num_qcstat_conv(2,sonde_sfc,1,1) = num_qcstat_conv(2,sonde_sfc,1,1) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'sonde_sfc',ob_vars(1),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
                     num_qcstat_conv(2,sonde_sfc,2,1) = num_qcstat_conv(2,sonde_sfc,2,1) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'sonde_sfc',ob_vars(2),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
                  endif
               endif

               if(wind_sd_sound)then
                  if(iv%sonde_sfc(n)%u%qc == fails_error_max .or. iv%sonde_sfc(n)%v%qc == fails_error_max .or. &
                     abs(iv%sonde_sfc(n)%v%inv) >= max_omb_dir .or. abs(iv%sonde_sfc(n)%u%inv) >= max_omb_spd )then
                     iv%sonde_sfc(n)%u%qc = fails_error_max
                     iv%sonde_sfc(n)%v%qc = fails_error_max
                     iv%sonde_sfc(n)%u%inv = 0.0
                     iv%sonde_sfc(n)%v%inv = 0.0
                  endif
               else
                  if(iv%sonde_sfc(n)%u%qc == fails_error_max .or. iv%sonde_sfc(n)%v%qc == fails_error_max ) then
                     iv%sonde_sfc(n)%u%qc = fails_error_max
                     iv%sonde_sfc(n)%v%qc = fails_error_max
                     iv%sonde_sfc(n)%u%inv = 0.0
                     iv%sonde_sfc(n)%v%inv = 0.0
                  endif
               endif
            endif

      failed=.false.
      if( iv%sonde_sfc(n)%t%qc >= obs_qc_pointer )  then
      call da_max_error_qc (it, iv%info(sonde_sfc), n, iv%sonde_sfc(n)%t, max_error_t , failed)
      if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
      num_qcstat_conv(1,sonde_sfc,3,1)= num_qcstat_conv(1,sonde_sfc,3,1) + 1
      if(failed) then
      num_qcstat_conv(2,sonde_sfc,3,1)= num_qcstat_conv(2,sonde_sfc,3,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'sonde_sfc',ob_vars(3),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%sonde_sfc(n)%p%qc >= obs_qc_pointer )  then
      call da_max_error_qc (it, iv%info(sonde_sfc), n, iv%sonde_sfc(n)%p, max_error_p , failed)         
      if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
      num_qcstat_conv(1,sonde_sfc,5,1)= num_qcstat_conv(1,sonde_sfc,5,1) + 1
      if(failed) then
      num_qcstat_conv(2,sonde_sfc,5,1)= num_qcstat_conv(2,sonde_sfc,5,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'sonde_sfc',ob_vars(5),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
      end if
      end if
      end if

      failed=.false.
      if( iv%sonde_sfc(n)%q%qc >= obs_qc_pointer ) then
       if( iv%sonde_sfc(n)%t%qc == fails_error_max .or. iv%sonde_sfc(n)%p%qc == fails_error_max) then
       failed=.true.
       iv%sonde_sfc(n)%q%qc  = fails_error_max
       iv%sonde_sfc(n)%q%inv = 0.0
       else
       call da_max_error_qc (it, iv%info(sonde_sfc), n, iv%sonde_sfc(n)%q, max_error_q , failed)
       endif
      if( iv%info(sonde_sfc)%proc_domain(1,n) ) then
      num_qcstat_conv(1,sonde_sfc,4,1)= num_qcstat_conv(1,sonde_sfc,4,1) + 1
      if(failed) then
      num_qcstat_conv(2,sonde_sfc,4,1)= num_qcstat_conv(2,sonde_sfc,4,1) + 1
      write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
           'sonde_sfc',ob_vars(4),iv%info(sonde_sfc)%lat(1,n),iv%info(sonde_sfc)%lon(1,n),0.01*ob%sonde_sfc(n)%p
      end if
      end if
      end if 
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_sonde_sfc")

end subroutine da_check_max_iv_sonde_sfc


subroutine da_calculate_grady_sonde_sfc(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector              
   !-------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_sonde_sfc")

   do n=1, iv%info(sound)%nlocal
      if (iv%sonde_sfc(n)%u%qc < obs_qc_pointer) re%sonde_sfc(n)%u = 0.0
      if (iv%sonde_sfc(n)%v%qc < obs_qc_pointer) re%sonde_sfc(n)%v = 0.0
      if (iv%sonde_sfc(n)%t%qc < obs_qc_pointer) re%sonde_sfc(n)%t = 0.0
      if (iv%sonde_sfc(n)%p%qc < obs_qc_pointer) re%sonde_sfc(n)%p = 0.0
      if (iv%sonde_sfc(n)%q%qc < obs_qc_pointer) re%sonde_sfc(n)%q = 0.0

      if (iv%sonde_sfc(n)%u%qc < obs_qc_pointer) re%sonde_sfc(n)%u = 0.0
      if (iv%sonde_sfc(n)%v%qc < obs_qc_pointer) re%sonde_sfc(n)%v = 0.0
      if (iv%sonde_sfc(n)%t%qc < obs_qc_pointer) re%sonde_sfc(n)%t = 0.0
      if (iv%sonde_sfc(n)%p%qc < obs_qc_pointer) re%sonde_sfc(n)%p = 0.0
      if (iv%sonde_sfc(n)%q%qc < obs_qc_pointer) re%sonde_sfc(n)%q = 0.0

      jo_grad_y%sonde_sfc(n)%u = -re%sonde_sfc(n)%u / &
          (iv%sonde_sfc(n)%u%error * iv%sonde_sfc(n)%u%error)
      jo_grad_y%sonde_sfc(n)%v = -re%sonde_sfc(n)%v / &
          (iv%sonde_sfc(n)%v%error * iv%sonde_sfc(n)%v%error)
      jo_grad_y%sonde_sfc(n)%t = -re%sonde_sfc(n)%t / &
          (iv%sonde_sfc(n)%t%error * iv%sonde_sfc(n)%t%error)
      jo_grad_y%sonde_sfc(n)%p = -re%sonde_sfc(n)%p / &
          (iv%sonde_sfc(n)%p%error * iv%sonde_sfc(n)%p%error)
      jo_grad_y%sonde_sfc(n)%q = -re%sonde_sfc(n)%q / &
          (iv%sonde_sfc(n)%q%error * iv%sonde_sfc(n)%q%error)
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_sonde_sfc")
     
end subroutine da_calculate_grady_sonde_sfc



end module da_sound

