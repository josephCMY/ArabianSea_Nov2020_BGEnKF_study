












module da_bogus

   use module_domain, only : domain

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, max_error_uv, max_error_t, rootproc, &
      bogus, max_error_p,max_error_q, trace_use_dull,fails_error_max, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, anal_type_verify, kms,kme,kts,kte, &
      ob_vars,qcstat_conv_unit
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_interp_lin_3d, da_to_zk, &
      da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_physics, only : da_tpq_to_slp_adj,da_tpq_to_slp_lin
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_bogus1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: q                        
      real          :: slp                      
   end type residual_bogus1_type

   type maxmin_bogus_stats_type
      type (maxmin_type)         :: u, v, t, q, slp 
   end type maxmin_bogus_stats_type

   type stats_bogus_type
      type (maxmin_bogus_stats_type)  :: maximum, minimum
      type (residual_bogus1_type)     :: average, rms_err
   end type stats_bogus_type

contains

subroutine da_ao_stats_bogus (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type (y_type),  intent (in)    :: re            ! A - O

   type (stats_bogus_type) :: stats
   integer                 :: nu, nv, nt, nq, nslp
   integer                 :: n, k
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_bogus")

   nu = 0
   nv = 0
   nt = 0
   nq = 0
   nslp = 0

   stats%maximum%u   = maxmin_type (missing_r, 0, 0)
   stats%maximum%v   = maxmin_type (missing_r, 0, 0)
   stats%maximum%t   = maxmin_type (missing_r, 0, 0)
   stats%maximum%q   = maxmin_type (missing_r, 0, 0)
   stats%maximum%slp = maxmin_type (missing_r, 0, 0)
   stats%minimum%u   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%slp = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_bogus1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(bogus)%nlocal
      if (iv%info(bogus)%proc_domain(1,n)) then
         call da_stats_calculate (n, 0, iv%bogus(n)%slp%qc,  &
            re%bogus(n)%slp, nslp,     &
            stats%minimum%slp, stats%maximum%slp,  &
            stats%average%slp, stats%rms_err%slp)

         do k=1, iv%info(bogus)%levels(n)
            call da_stats_calculate (n, k, iv%bogus(n)%u(k)%qc, & 
               re%bogus(n)%u(k), nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%bogus(n)%v(k)%qc, & 
               re%bogus(n)%v(k), nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
            call da_stats_calculate (n, k, iv%bogus(n)%t(k)%qc, & 
               re%bogus(n)%t(k), nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate (n, k, iv%bogus(n)%q(k)%qc, & 
               re%bogus(n)%q(k), nq, &
               stats%minimum%q, stats%maximum%q, &
               stats%average%q, stats%rms_err%q)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (nq)
   call da_proc_sum_int (nslp)
   iv%nstats(bogus) = nu + nv + nt + nq + nslp 

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
   call da_proc_stats_combine(stats%average%slp, stats%rms_err%slp, &
      stats%minimum%slp%value, stats%maximum%slp%value, &
      stats%minimum%slp%n, stats%maximum%slp%n, &
      stats%minimum%slp%l, stats%maximum%slp%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. nq /= 0 .or. nslp /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for bogus'
         call da_print_stats_bogus(stats_unit, nu, nv, nt, nq, nslp, stats)
      end if
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_bogus")

end subroutine da_ao_stats_bogus


subroutine da_jo_and_grady_bogus(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(in)    :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_bogus") 

   jo%bogus_u   = 0.0
   jo%bogus_v   = 0.0
   jo%bogus_t   = 0.0
   jo%bogus_q   = 0.0
   jo%bogus_slp = 0.0

   if (iv%info(bogus)%nlocal > 0) then
      do n=1, iv%info(bogus)%nlocal
         jo_grad_y%bogus(n)%slp = -re%bogus(n)%slp &
            / (iv%bogus(n)%slp%error * iv%bogus(n)%slp%error)
         do k=1, iv%info(bogus)%levels(n)
            jo_grad_y%bogus(n)%u(k) = -re%bogus(n)%u(k) &
               / (iv%bogus(n)%u(k)%error *iv%bogus(n)%u(k)%error)
            jo_grad_y%bogus(n)%v(k) = -re%bogus(n)%v(k) &
               /(iv%bogus(n)%v(k)%error * iv%bogus(n)%v(k)%error)
            jo_grad_y%bogus(n)%t(k) = -re%bogus(n)%t(k) &
               / (iv%bogus(n)%t(k)%error * iv%bogus(n)%t(k)%error)
            jo_grad_y%bogus(n)%q(k) = -re%bogus(n)%q(k) &
               / (iv%bogus(n)%q(k)%error * iv%bogus(n)%q(k)%error)
         end do

         if (iv%info(bogus)%proc_domain(1,n)) then
            jo%bogus_slp = jo%bogus_slp - re%bogus(n)%slp * jo_grad_y%bogus(n)%slp

            do k=1, iv%info(bogus)%levels(n)
               jo%bogus_u = jo%bogus_u - re%bogus(n)%u(k) * jo_grad_y%bogus(n)%u(k)
               jo%bogus_v = jo%bogus_v - re%bogus(n)%v(k) * jo_grad_y%bogus(n)%v(k)
               jo%bogus_t = jo%bogus_t - re%bogus(n)%t(k) * jo_grad_y%bogus(n)%t(k)
               jo%bogus_q = jo%bogus_q - re%bogus(n)%q(k) * jo_grad_y%bogus(n)%q(k)
            end do
         end if 
      end do

      jo % bogus_slp = 0.5 * jo % bogus_slp
      jo % bogus_u = 0.5 * jo % bogus_u
      jo % bogus_v = 0.5 * jo % bogus_v
      jo % bogus_t = 0.5 * jo % bogus_t
      jo % bogus_q = 0.5 * jo % bogus_q
   end if

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_bogus") 

end subroutine da_jo_and_grady_bogus


subroutine da_residual_bogus(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

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

   type (bad_data_type) :: n_obs_bad
   integer              :: n, k

   if (trace_use_dull) call da_trace_entry("da_residual_bogus")      

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)
   n_obs_bad % slp % num = number_type(0, 0, 0)

   do n=1, iv%info(bogus)%nlocal
      do k=1, iv%info(bogus)%levels(n)
         np_available = np_available + 4
         re%bogus(n)%u(k) = da_residual(n, k, y%bogus(n)%u(k), iv%bogus(n)%u(k), n_obs_bad % u)
         re%bogus(n)%v(k) = da_residual(n, k, y%bogus(n)%v(k), iv%bogus(n)%v(k), n_obs_bad % v)
         re%bogus(n)%t(k) = da_residual(n, k, y%bogus(n)%t(k), iv%bogus(n)%t(k), n_obs_bad % t)
         re%bogus(n)%q(k) = da_residual(n, k, y%bogus(n)%q(k), iv%bogus(n)%q(k), n_obs_bad % q)
      end do

      np_available = np_available + 1
      re%bogus(n)%slp = da_residual(n, 0, y%bogus(n)%slp, iv%bogus(n)%slp, n_obs_bad % slp) 
   end do

   np_missing = np_missing + n_obs_bad % u % num % miss + &
      n_obs_bad % v % num % miss + n_obs_bad % slp % num % miss + &
      n_obs_bad % t % num % miss + n_obs_bad % q % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad + &
      n_obs_bad % v % num % bad + n_obs_bad % slp % num % bad  + &
      n_obs_bad % t % num % bad + n_obs_bad % q % num % bad
   np_obs_used = np_obs_used + n_obs_bad % u % num % use + &
      n_obs_bad % v % num % use + n_obs_bad % slp % num % use  + &
      n_obs_bad % t % num % use + n_obs_bad % q % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_bogus")   

end subroutine da_residual_bogus


subroutine da_oi_stats_bogus (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_bogus_type) :: stats
   integer                 :: nu, nv, nt, nq, nslp
   integer                 :: n, k

   if (trace_use_dull) call da_trace_entry("da_oi_stats_bogus")

   nslp = 0
   nu = 0
   nv = 0
   nt = 0
   nq = 0

   stats%maximum%u   = maxmin_type(missing_r, 0, 0)
   stats%maximum%v   = maxmin_type(missing_r, 0, 0)
   stats%maximum%t   = maxmin_type(missing_r, 0, 0)
   stats%maximum%q   = maxmin_type(missing_r, 0, 0)
   stats%maximum%slp = maxmin_type(missing_r, 0, 0)
   stats%minimum%u   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%slp = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_bogus1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(bogus)%nlocal
      if (iv%info(bogus)%proc_domain(1,n)) then
         do k=1, iv%info(bogus)%levels(n)
            call da_stats_calculate(iv%info(bogus)%obs_global_index(n), &
               k, iv%bogus(n)%u(k)%qc, &
               iv%bogus(n)%u(k)%inv, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(bogus)%obs_global_index(n), &
               k, iv%bogus(n)%v(k)%qc, &
               iv%bogus(n)%v(k)%inv, nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
            call da_stats_calculate(iv%info(bogus)%obs_global_index(n), &
               k, iv%bogus(n)%t(k)%qc, &
               iv%bogus(n)%t(k)%inv, nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate(iv%info(bogus)%obs_global_index(n), &
               k, iv%bogus(n)%q(k)%qc, &
               iv%bogus(n)%q(k)%inv, nq, &
               stats%minimum%q, stats%maximum%q, &
               stats%average%q, stats%rms_err%q)
         end do
         call da_stats_calculate( iv%info(bogus)%obs_global_index(n), &
            0, iv%bogus(n)%slp%qc, &
            iv%bogus(n)%slp%inv, nslp,  &
            stats%minimum%slp, stats%maximum%slp, &
            stats%average%slp, stats%rms_err%slp)
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nu)
   call da_proc_sum_int(nv)
   call da_proc_sum_int(nt)
   call da_proc_sum_int(nq)
   call da_proc_sum_int(nslp)

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
   call da_proc_stats_combine(stats%average%slp, stats%rms_err%slp, &
      stats%minimum%slp%value, stats%maximum%slp%value, &
      stats%minimum%slp%n, stats%maximum%slp%n, &
      stats%minimum%slp%l, stats%maximum%slp%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. nq /= 0 .or. nslp /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for bogus'
         call da_print_stats_bogus(stats_unit, nu, nv, nt, nq, nslp, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_bogus")

end subroutine da_oi_stats_bogus


subroutine da_print_stats_bogus(stats_unit, nu, nv, nt, nq, nslp, Bogus)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nu, nv, nt, nq, nslp
   type (stats_bogus_type), intent(in)    :: bogus

   if (trace_use_dull) call da_trace_entry("da_print_stats_bogus") 

   write(unit=stats_unit, fmt='(6a/)') &
      '   var             ', &
      'u (m/s)     n    k    ', &
      'v (m/s)     n    k    ', &
      't (K)       n    k    ', &
      'q (kg/kg)   n    k    ', &
      'slp (pa)    n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
      '  Number: ', nu, nv, nt, nq, nslp

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1
   if (nt < 1) nt = 1
   if (nq < 1) nq = 1
   if (nslp < 1) nslp = 1

   write(unit=stats_unit, fmt='((a,3(f12.4,2i5),e12.4,2i5,f12.4,2i5))') &
      ' Minimum(n,k): ', Bogus%minimum%u, Bogus%minimum%v, Bogus%minimum%t, &
                         Bogus%minimum%q, Bogus%minimum%slp,                &
      ' Maximum(n,k): ', Bogus%maximum%u, Bogus%maximum%v, Bogus%maximum%t, &
                         Bogus%maximum%q, Bogus%maximum%slp
   write(unit=stats_unit, fmt='((a,3(f12.4,10x),e12.4,10x,f12.4,10x))') &
      ' Average     : ', Bogus%average%u/real(nu), &
                         Bogus%average%v/real(nv), &
                         Bogus%average%t/real(nt), &
                         Bogus%average%q/real(nq), &
                         Bogus%average%slp/real(nslp),  &
      '    RMSE     : ', sqrt(Bogus%rms_err%u/real(nu)), &
                         sqrt(Bogus%rms_err%v/real(nv)), &
                         sqrt(Bogus%rms_err%t/real(nt)), &
                         sqrt(Bogus%rms_err%q/real(nq)), &
                         sqrt(Bogus%rms_err%slp/real(nslp))

   if (trace_use_dull) call da_trace_exit("da_print_stats_bogus") 

end subroutine da_print_stats_bogus


subroutine da_transform_xtoy_bogus (grid, iv, y)

   !---------------------------------------------------------------------
   ! Purpose: the linearized bogus observation operators.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs ! Number of obs levels.
   real    :: dx, dxm 
   real    :: dy, dym
   real    :: model_t(kts:kte)
   real    :: model_t9(kts:kte)
   real    :: model_q(kts:kte)
   real    :: model_q9(kts:kte)
   real    :: model_p_c(kts:kte)
   real    :: model_p_c9(kts:kte)
   real    :: psfc_loc,terr_loc,psfc_loc9

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_bogus")      

   do n = iv%info(bogus)%n1,iv%info(bogus)%n2
      num_levs = iv%info(bogus)%levels(n)

      if (num_levs < 1) cycle

      ! [1.1] Get cross pt. horizontal interpolation weights:

      i   = iv%info(bogus)%i(1,n)
      dy  = iv%info(bogus)%dy(1,n)
      dym = iv%info(bogus)%dym(1,n)
      j   = iv%info(bogus)%j(1,n)
      dx  = iv%info(bogus)%dx(1,n)
      dxm = iv%info(bogus)%dxm(1,n)

      ! [1.2] Interpolate horizontally from cross points:

      do k = kts, kte
         model_t9(k) = dym*(dxm*grid%xa%t(i,j,k) + dx*grid%xa%t(i+1,j,k)) &
            + dy *(dxm*grid%xa%t(i,j+1,k) + dx*grid%xa%t(i+1,j+1,k))
         model_t(k) = dym*(dxm*grid%xb%t(i,j,k) + dx*grid%xb%t(i+1,j,k)) &
            + dy *(dxm*grid%xb%t(i,j+1,k) + dx*grid%xb%t(i+1,j+1,k))
         model_q9(k) = dym*(dxm*grid%xa%q(i,j,k) + dx*grid%xa%q(i+1,j,k)) &
            + dy *(dxm*grid%xa%q(i,j+1,k) + dx*grid%xa%q(i+1,j+1,k))
         model_q(k) = dym*(dxm*grid%xb%q(i,j,k) + dx*grid%xb%q(i+1,j,k)) &
            + dy *(dxm*grid%xb%q(i,j+1,k) + dx*grid%xb%q(i+1,j+1,k))
         model_p_c9(k) = dym*(dxm*grid%xa%p(i,j,k) + dx*grid%xa%p(i+1,j,k)) &
            + dy *(dxm*grid%xa%p(i,j+1,k) + dx*grid%xa%p(i+1,j+1,k))
         model_p_c(k) = dym*(dxm*grid%xb%p(i,j,k) + dx*grid%xb%p(i+1,j,k)) &
            + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
      end do

      terr_loc = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) &
         + dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))
      psfc_loc = dym*(dxm*grid%xb%psfc(i,j)   + dx*grid%xb%psfc(i+1,j)) &
         + dy *(dxm*grid%xb%psfc(i,j+1) + dx*grid%xb%psfc(i+1,j+1))
      psfc_loc9 = dym*(dxm*grid%xa%psfc(i,j)   + dx*grid%xa%psfc(i+1,j)) &
         + dy *(dxm*grid%xa%psfc(i,j+1) + dx*grid%xa%psfc(i+1,j+1))

      call da_tpq_to_slp_lin (model_t, model_q, model_p_c, terr_loc, psfc_loc,   &
         model_t9, model_q9, model_p_c9, psfc_loc9, y%bogus(n)%slp) 
   end do

   allocate (u(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (v(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (t(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (q(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   u(:,:)=0.0
   v(:,:)=0.0
   t(:,:)=0.0
   q(:,:)=0.0
  
   call da_interp_lin_3d (grid%xa%u, iv%info(bogus), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(bogus), v)
   call da_interp_lin_3d (grid%xa%t, iv%info(bogus), t)
   call da_interp_lin_3d (grid%xa%q, iv%info(bogus), q)

   do n=iv%info(bogus)%n1,iv%info(bogus)%n2
      y%bogus(n)%u(:) = u(1:size(y%bogus(n)%u),n)
      y%bogus(n)%v(:) = v(1:size(y%bogus(n)%v),n)
      y%bogus(n)%t(:) = t(1:size(y%bogus(n)%t),n)
      y%bogus(n)%q(:) = q(1:size(y%bogus(n)%q),n)
   end do

   deallocate (u)
   deallocate (v)
   deallocate (t)
   deallocate (q)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_bogus") 

end subroutine da_transform_xtoy_bogus


subroutine da_transform_xtoy_bogus_adj(grid, iv, jo_grad_y, jo_grad_x)

   !---------------------------------------------------------------------
   ! Purpose: the adjoint of bogus observation operators.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------

   implicit none

   type (domain),     intent(in)     :: grid
   type (iv_type),    intent(in)     :: iv          ! obs. inc vector (o-b).
   type (y_type) ,    intent(inout)  :: jo_grad_y   ! grad_y(jo)
   type (x_type) ,    intent(inout)  :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs
   real    :: dx, dxm
   real    :: dy, dym

   real    :: model_t(kms:kme)
   real    :: tt(kms:kme)
   real    :: model_q(kms:kme)
   real    :: qq(kms:kme)
   real    :: model_p_c(kms:kme)
   real    :: pp(kms:kme)
   real    :: psfc_loc,terr_loc,ppsfc   

   real, allocatable :: temp_u(:,:)
   real, allocatable :: temp_v(:,:)
   real, allocatable :: temp_t(:,:)
   real, allocatable :: temp_q(:,:)           

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_bogus_adj")

   allocate (temp_u(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (temp_v(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (temp_t(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (temp_q(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))

   do n=iv%info(bogus)%n1,iv%info(bogus)%n2
      temp_u(1:size(jo_grad_y%bogus(n)%u),n)  = jo_grad_y%bogus(n)%u
      temp_v(1:size(jo_grad_y%bogus(n)%u),n)  = jo_grad_y%bogus(n)%v
      temp_t(1:size(jo_grad_y%bogus(n)%u),n)  = jo_grad_y%bogus(n)%t
      temp_q(1:size(jo_grad_y%bogus(n)%u),n)  = jo_grad_y%bogus(n)%q
   end do

   ! [1.1] Adjoint feedback from Y to X for u and v:

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(bogus), temp_u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(bogus), temp_v)
   call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(bogus), temp_t)
   call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(bogus), temp_q)
   deallocate (temp_u)
   deallocate (temp_v)
   deallocate (temp_t)
   deallocate (temp_q)

   do n= iv%info(bogus)%n1, iv%info(bogus)%n2
      num_levs = iv%info(bogus)%levels(n)
      if (num_levs < 1) cycle

      ! [1.1] Get cross pt. horizontal interpolation weights:

      i   = iv%info(bogus)%i(1,n)
      dy  = iv%info(bogus)%dy(1,n)
      dym = iv%info(bogus)%dym(1,n)
      j   = iv%info(bogus)%j(1,n)
      dx  = iv%info(bogus)%dx(1,n)
      dxm = iv%info(bogus)%dxm(1,n)

      ! [1.2] Compute the feedback from SLP to t and q:

      ! 1.2.1 Background at the observation location:

      do k = kts, kte
         model_t(k) = dym*(dxm*grid%xb%t(i,j,k) + dx*grid%xb%t(i+1,j,k)) &
            + dy *(dxm*grid%xb%t(i,j+1,k) + dx*grid%xb%t(i+1,j+1,k))
         model_q(k) = dym*(dxm*grid%xb%q(i,j,k) + dx*grid%xb%q(i+1,j,k)) &
            + dy *(dxm*grid%xb%q(i,j+1,k) + dx*grid%xb%q(i+1,j+1,k))
         model_p_c(k) = dym*(dxm*grid%xb%p(i,j,k) + dx*grid%xb%p(i+1,j,k)) &
            + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
      end do

      terr_loc = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) &
         + dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))
      psfc_loc = dym*(dxm*grid%xb%psfc(i,j)   + dx*grid%xb%psfc(i+1,j)) &
         + dy *(dxm*grid%xb%psfc(i,j+1) + dx*grid%xb%psfc(i+1,j+1))

      ! 1.2.2 Compute the feedback from SLP to p, t, q, and psfc 
      !       at the observed location:

      call da_tpq_to_slp_adj(model_t, model_q ,model_p_c, terr_loc, psfc_loc, &
         tt, qq, pp, ppsfc, jo_grad_y%bogus(n)%slp)  

      ! 1.2.3 feedback from the observed location to grid space:

      ! 1.2.3.1 for Psfc

      jo_grad_x % psfc(i,j)     = jo_grad_x % psfc(i,j) + dym*dxm*ppsfc
      jo_grad_x % psfc(i+1,j)   = jo_grad_x % psfc(i+1,j)+ dym*dx*ppsfc
      jo_grad_x % psfc(i,j+1)   = jo_grad_x % psfc(i,j+1)+ dy*dxm*ppsfc
      jo_grad_x % psfc(i+1,j+1) = jo_grad_x % psfc(i+1,j+1)+dy*dx*ppsfc

      ! 1.2.3.2 for t, q, p

      do k = kts, kte
         jo_grad_x % t(i,j,k)     = jo_grad_x % t(i,j,k)+dym*dxm*tt(k)
         jo_grad_x % t(i+1,j,k)   = jo_grad_x % t(i+1,j,k)+ dym*dx*tt(k)
         jo_grad_x % t(i,j+1,k)   = jo_grad_x % t(i,j+1,k)+ dy*dxm*tt(k)
         jo_grad_x % t(i+1,j+1,k) = jo_grad_x % t(i+1,j+1,k)+dy*dx*tt(k)
         jo_grad_x % q(i,j,k)     = jo_grad_x % q(i,j,k)+dym*dxm*qq(k)
         jo_grad_x % q(i+1,j,k)   = jo_grad_x % q(i+1,j,k)+ dym*dx*qq(k)
         jo_grad_x % q(i,j+1,k)   = jo_grad_x % q(i,j+1,k)+ dy*dxm*qq(k)
         jo_grad_x % q(i+1,j+1,k) = jo_grad_x % q(i+1,j+1,k)+dy*dx*qq(k)
         jo_grad_x % p(i,j,k)     = jo_grad_x % p(i,j,k)+dym*dxm*pp(k)
         jo_grad_x % p(i+1,j,k)   = jo_grad_x % p(i+1,j,k)+ dym*dx*pp(k)
         jo_grad_x % p(i,j+1,k)   = jo_grad_x % p(i,j+1,k)+ dy*dxm*pp(k)
         jo_grad_x % p(i+1,j+1,k) = jo_grad_x % p(i+1,j+1,k)+dy*dx*pp(k)
      end do
   end do                  

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_bogus_adj")

end subroutine da_transform_xtoy_bogus_adj


subroutine da_check_max_iv_bogus(iv,ob, it, num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)
   type(y_type),  intent(in)    :: ob      ! Observation structure.



   integer :: k,n, ipr
   logical :: failed
   
   if (trace_use_dull) call da_trace_entry("da_check_max_iv_bogus")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(bogus)%n1,iv%info(bogus)%n2
      do k = 1, iv%info(bogus)%levels(n)
         call da_get_print_lvl(iv%bogus(n)%p(k),ipr)

        failed=.false.
        if( iv%bogus(n)%u(k)%qc >= obs_qc_pointer ) then
         call da_max_error_qc (it,iv%info(bogus), n, iv%bogus(n)%u(k), max_error_buv,failed)
        if( iv%info(bogus)%proc_domain(k,n) ) then
           num_qcstat_conv(1,bogus,1,ipr) = num_qcstat_conv(1,bogus,1,ipr) + 1
         if(failed)then
          num_qcstat_conv(2,bogus,1,ipr) = num_qcstat_conv(2,bogus,1,ipr) + 1
          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
          'bogus',ob_vars(1),iv%info(bogus)%lat(k,n),iv%info(bogus)%lon(k,n),0.01*iv%bogus(n)%p(k)
         end if
        end if
        end if

        failed=.false.
        if( iv%bogus(n)%v(k)%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it,iv%info(bogus), n, iv%bogus(n)%v(k), max_error_buv,failed)
        if( iv%info(bogus)%proc_domain(k,n) ) then
          num_qcstat_conv(1,bogus,2,ipr) = num_qcstat_conv(1,bogus,2,ipr) + 1
         if(failed)then
          num_qcstat_conv(2,bogus,2,ipr) = num_qcstat_conv(2,bogus,2,ipr) + 1
          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
          'bogus',ob_vars(2),iv%info(bogus)%lat(k,n),iv%info(bogus)%lon(k,n),0.01*iv%bogus(n)%p(k)
         end if
        end if
        end if

        failed=.false.
        if( iv%bogus(n)%t(k)%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it,iv%info(bogus), n, iv%bogus(n)%t(k), max_error_bt ,failed)
        if( iv%info(bogus)%proc_domain(k,n) ) then
           num_qcstat_conv(1,bogus,3,ipr) = num_qcstat_conv(1,bogus,3,ipr) + 1
         if(failed)then
          num_qcstat_conv(2,bogus,3,ipr) = num_qcstat_conv(2,bogus,3,ipr) + 1
          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
          'bogus',ob_vars(3),iv%info(bogus)%lat(k,n),iv%info(bogus)%lon(k,n),0.01*iv%bogus(n)%p(k)
         end if
        end if
        end if

        failed=.false.
        if( iv%bogus(n)%q(k)%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it,iv%info(bogus), n, iv%bogus(n)%q(k), max_error_bq ,failed)
        if( iv%info(bogus)%proc_domain(k,n) ) then
           num_qcstat_conv(1,bogus,4,ipr) = num_qcstat_conv(1,bogus,4,ipr) + 1
         if(failed)then
          num_qcstat_conv(2,bogus,4,ipr) = num_qcstat_conv(2,bogus,4,ipr) + 1
          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
          'bogus',ob_vars(4),iv%info(bogus)%lat(k,n),iv%info(bogus)%lon(k,n),0.01*iv%bogus(n)%p(k)
         end if
        end if
        end if

      end do 
      ! Sea Level Pressure

      if( iv%info(bogus)%proc_domain(1,n) ) then
        failed=.false.
        if( iv%bogus(n)%slp%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it,iv%info(bogus), n, iv%bogus(n)%slp, max_error_slp ,failed)
           num_qcstat_conv(1,bogus,5,1) = num_qcstat_conv(1,bogus,5,1) + 1
         if(failed) then
          num_qcstat_conv(2,bogus,5,1) = num_qcstat_conv(2,bogus,5,1) + 1
          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
          'bogus',ob_vars(5),iv%info(bogus)%lat(1,n),iv%info(bogus)%lon(1,n),ob%bogus(n)%slp
        endif
      endif
      endif

   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_max_iv_bogus")

end subroutine da_check_max_iv_bogus
subroutine da_get_innov_vector_bogus(it,num_qcstat_conv, grid, ob, iv)

   !------------------------------------------------------------------------------
   ! Purpose: calculate the innovations for the bogus data.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !------------------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it       ! External iteration.
   type(domain),     intent(in)    :: grid     ! first guess state.
   type(y_type),     intent(in)    :: ob       ! Observation structure.
   type(iv_type),    intent(inout) :: iv       ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs ! Number of obs levels.

   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.

   real, allocatable :: model_u(:,:)  ! Model value u at ob location.
   real, allocatable :: model_v(:,:)  ! Model value v at ob location.
   real, allocatable :: model_t(:,:)  ! Model value t at ob location.
   real, allocatable :: model_q(:,:)  ! Model value t at ob location.
   real    :: model_slp                  ! Model value slp at ob location.

   
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_bogus")

   allocate (model_u(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (model_v(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (model_t(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   allocate (model_q(iv%info(bogus)%max_lev,iv%info(bogus)%n1:iv%info(bogus)%n2))
   model_u(:,:) = 0.0
   model_v(:,:) = 0.0
   model_t(:,:) = 0.0
   model_q(:,:) = 0.0

   if ( it > 1 ) then
      do n=iv%info(bogus)%n1,iv%info(bogus)%n2
         do k=1, iv%info(bogus)%levels(n)
            if (iv%bogus(n)%u(k)%qc == fails_error_max) iv%bogus(n)%u(k)%qc = 0
            if (iv%bogus(n)%v(k)%qc == fails_error_max) iv%bogus(n)%v(k)%qc = 0
            if (iv%bogus(n)%t(k)%qc == fails_error_max) iv%bogus(n)%t(k)%qc = 0
            if (iv%bogus(n)%q(k)%qc == fails_error_max) iv%bogus(n)%q(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(bogus)%n1,iv%info(bogus)%n2
      num_levs = iv%info(bogus)%levels(n)

      if (num_levs < 1) cycle

      i   = iv%info(bogus)%i(1,n)
      j   = iv%info(bogus)%j(1,n)
      dx  = iv%info(bogus)%dx(1,n)
      dy  = iv%info(bogus)%dy(1,n)
      dxm = iv%info(bogus)%dxm(1,n)
      dym = iv%info(bogus)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j,k) + dx*grid%xb%h(i+1,j,k)) + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
         v_p(k) = dym*(dxm*grid%xb%p(i,j,k) + dx*grid%xb%p(i+1,j,k)) + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
      end do

      do k=1, iv%info(bogus)%levels(n)
         if (iv % bogus(n) % p(k) > 1.0) then
            call da_to_zk(iv % bogus(n) % p(k), v_p, v_interp_p, iv%info(bogus)%zk(k,n))
         else if (iv % bogus(n) % h(k) > 0.0) then
            call da_to_zk(iv % bogus(n) % h(k), v_h, v_interp_h, iv%info(bogus)%zk(k,n))
         end if

         if (iv%info(bogus)%zk(k,n) < 0.0 .and.  .not.anal_type_verify) then
            iv % bogus(n) % u(k) % qc = missing_data
            iv % bogus(n) % v(k) % qc = missing_data
            iv % bogus(n) % t(k) % qc = missing_data
            iv % bogus(n) % q(k) % qc = missing_data
         end if
      end do
   end do

   call da_convert_zk (iv%info(bogus))

   ! [1.4] Interpolate horizontally:

   call da_interp_lin_3d (grid%xb%u, iv%info(bogus), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(bogus), model_v)
   call da_interp_lin_3d (grid%xb%t, iv%info(bogus), model_t)
   call da_interp_lin_3d (grid%xb%q, iv%info(bogus), model_q)

   do n=iv%info(bogus)%n1,iv%info(bogus)%n2
      num_levs = iv%info(bogus)%levels(n)

      if (num_levs < 1) cycle

      i   = iv%info(bogus)%i(1,n)
      j   = iv%info(bogus)%j(1,n)
      dx  = iv%info(bogus)%dx(1,n)
      dy  = iv%info(bogus)%dy(1,n)
      dxm = iv%info(bogus)%dxm(1,n)
      dym = iv%info(bogus)%dym(1,n)

      model_slp = dym*(dxm*grid%xb%slp(i,j)   + dx*grid%xb%slp(i+1,j)) &
         + dy *(dxm*grid%xb%slp(i,j+1) + dx*grid%xb%slp(i+1,j+1))

      !------------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !------------------------------------------------------------------------

      iv % bogus(n) % slp % inv = 0.0

      if (ABS(ob % bogus(n) % slp - missing_r) > 1.0 .AND. &
           iv % bogus(n) % slp % qc >= obs_qc_pointer) then
        iv % bogus(n) % slp % inv = ob % bogus(n) % slp - model_slp
      end if

      do k = 1, iv%info(bogus)%levels(n)
         iv % bogus(n) % u(k) % inv = 0.0
         iv % bogus(n) % v(k) % inv = 0.0
         iv % bogus(n) % t(k) % inv = 0.0
         iv % bogus(n) % q(k) % inv = 0.0

         !------------------------------------------------------------------------
         ! [4.0] Fast interpolation:
         !------------------------------------------------------------------------

         if (ob % bogus(n) % u(k) > missing_r .AND. iv % bogus(n) % u(k) % qc >= obs_qc_pointer) then
           iv % bogus(n) % u(k) % inv = ob % bogus(n) % u(k) - model_u(k,n)
         end if

         if (ob % bogus(n) % v(k) > missing_r .AND. iv % bogus(n) % v(k) % qc >= obs_qc_pointer) then
           iv % bogus(n) % v(k) % inv = ob % bogus(n) % v(k) - model_v(k,n)
         end if

         if (ob % bogus(n) % t(k) > missing_r .AND. iv % bogus(n) % t(k) % qc >= obs_qc_pointer) then
            ! only for global Bogus(YRG 07/15/2005):
            if (iv%info(bogus)%platform(n)(8:12) /= 'TCBOG') then
               iv % bogus(n) % t(k) % inv = ob % bogus(n) % t(k) - model_t(k,n)
            else
               iv % bogus(n) % t(k) % inv = missing_r 
               iv % bogus(n) % t(k) % qc  = missing_data
            end if
         end if

         if (ob % bogus(n) % q(k) > missing_r .AND. iv % bogus(n) % q(k) % qc >= obs_qc_pointer) then
            ! only for global Bogus(YRG 07/15/2005):
            if (iv%info(bogus)%platform(n)(8:12) /= 'TCBOG') then
               iv % bogus(n) % q(k) % inv = ob % bogus(n) % q(k) - model_q(k,n)
            else
              iv % bogus(n) % q(k) % inv = missing_r 
              iv % bogus(n) % q(k) % qc  = missing_data
            end if
         end if
      end do
   end do

   !------------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_bogus(iv,ob, it, num_qcstat_conv)

    deallocate (model_u)
    deallocate (model_v)
    deallocate (model_t)
    deallocate (model_q)

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_bogus")

end subroutine da_get_innov_vector_bogus


subroutine da_calculate_grady_bogus(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n, k
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_bogus")

   do n=1, iv%info(bogus)%nlocal
      if (iv%bogus(n)%slp%qc < obs_qc_pointer) then
         re%bogus(n)%slp = 0.0
      end if

      jo_grad_y%bogus(n)%slp = -re%bogus(n)%slp / (iv%bogus(n)%slp%error * iv%bogus(n)%slp%error)

      do k=1, iv%info(bogus)%levels(n)
         if (iv%bogus(n)%u(k)%qc < obs_qc_pointer) re%bogus(n)%u(k) = 0.0
         if (iv%bogus(n)%v(k)%qc < obs_qc_pointer) re%bogus(n)%v(k) = 0.0
         if (iv%bogus(n)%t(k)%qc < obs_qc_pointer) re%bogus(n)%t(k) = 0.0
         if (iv%bogus(n)%q(k)%qc < obs_qc_pointer) re%bogus(n)%q(k) = 0.0

         jo_grad_y%bogus(n)%u(k) = -re%bogus(n)%u(k) / (iv%bogus(n)%u(k)%error * iv%bogus(n)%u(k)%error)
         jo_grad_y%bogus(n)%v(k) = -re%bogus(n)%v(k) / (iv%bogus(n)%v(k)%error * iv%bogus(n)%v(k)%error)
         jo_grad_y%bogus(n)%t(k) = -re%bogus(n)%t(k) / (iv%bogus(n)%t(k)%error * iv%bogus(n)%t(k)%error)
         jo_grad_y%bogus(n)%q(k) = -re%bogus(n)%q(k) / (iv%bogus(n)%q(k)%error * iv%bogus(n)%q(k)%error)
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_bogus")

end subroutine da_calculate_grady_bogus



end module da_bogus

