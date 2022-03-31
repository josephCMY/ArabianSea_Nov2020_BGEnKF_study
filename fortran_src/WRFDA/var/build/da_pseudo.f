












module da_pseudo

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, trace_use_dull, &
      missing, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, pseudo, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, pseudo_var
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use module_domain, only : domain
   use da_interpolation, only : da_interp_lin_3d,da_interp_lin_3d_adj
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_residual, da_convert_zk
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_pseudo1_type
      real          :: u                        
      real          :: v                        
      real          :: t                        
      real          :: p                        
      real          :: q                        
   end type residual_pseudo1_type

   type maxmin_pseudo_stats_type
      type (maxmin_type)         :: u, v, t, p, q
   end type maxmin_pseudo_stats_type

   type stats_pseudo_type
      type (maxmin_pseudo_stats_type)  :: maximum, minimum
      type (residual_pseudo1_type)     :: average, rms_err
   end type stats_pseudo_type

contains

subroutine da_jo_and_grady_pseudo(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout):: jo          ! Obs cost function.

   integer :: n

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_pseudo")

   jo % pseudo_u = 0.0
   jo % pseudo_v = 0.0
   jo % pseudo_t = 0.0
   jo % pseudo_p = 0.0
   jo % pseudo_q = 0.0

   do n=1, iv%info(pseudo)%nlocal
      jo_grad_y%pseudo(n)%u = -re%pseudo(n)%u / (iv%pseudo(n)%u%error * iv%pseudo(n)%u%error)
      jo_grad_y%pseudo(n)%v = -re%pseudo(n)%v / (iv%pseudo(n)%v%error * iv%pseudo(n)%v%error)
      jo_grad_y%pseudo(n)%t = -re%pseudo(n)%t / (iv%pseudo(n)%t%error * iv%pseudo(n)%t%error)
      jo_grad_y%pseudo(n)%p = -re%pseudo(n)%p / (iv%pseudo(n)%p%error * iv%pseudo(n)%p%error)
      jo_grad_y%pseudo(n)%q = -re%pseudo(n)%q / (iv%pseudo(n)%q%error * iv%pseudo(n)%q%error)

      if (iv%info(pseudo)%proc_domain(1,n)) then
         jo % pseudo_u = jo % pseudo_u - re%pseudo(n)%u * jo_grad_y%pseudo(n)%u
         jo % pseudo_v = jo % pseudo_v - re%pseudo(n)%v * jo_grad_y%pseudo(n)%v
         jo % pseudo_t = jo % pseudo_t - re%pseudo(n)%t * jo_grad_y%pseudo(n)%t
         jo % pseudo_p = jo % pseudo_p - re%pseudo(n)%p * jo_grad_y%pseudo(n)%p
         jo % pseudo_q = jo % pseudo_q - re%pseudo(n)%q * jo_grad_y%pseudo(n)%q
      end if
   end do

   jo % pseudo_u = 0.5 * jo % pseudo_u
   jo % pseudo_v = 0.5 * jo % pseudo_v
   jo % pseudo_t = 0.5 * jo % pseudo_t
   jo % pseudo_p = 0.5 * jo % pseudo_p
   jo % pseudo_q = 0.5 * jo % pseudo_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_pseudo")
  
end subroutine da_jo_and_grady_pseudo


subroutine da_residual_pseudo(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for pseudo obs
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
   integer              :: n

   if (trace_use_dull) call da_trace_entry("da_residual_pseudo")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)
   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % p % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(pseudo)%nlocal
      if (iv%info(pseudo)%proc_domain(1,n)) then
         np_available = np_available + 5
         re%pseudo(n)%u = da_residual(n, 0, y%pseudo(n)%u, iv%pseudo(n)%u, n_obs_bad % u)
         re%pseudo(n)%v = da_residual(n, 0, y%pseudo(n)%v, iv%pseudo(n)%v, n_obs_bad % v)
         re%pseudo(n)%t = da_residual(n, 0, y%pseudo(n)%t, iv%pseudo(n)%t, n_obs_bad % t)
         re%pseudo(n)%p = da_residual(n, 0, y%pseudo(n)%p, iv%pseudo(n)%p, n_obs_bad % p)
         re%pseudo(n)%q = da_residual(n, 0, y%pseudo(n)%q, iv%pseudo(n)%q, n_obs_bad % q)
      end if
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

   if (trace_use_dull) call da_trace_exit("da_residual_pseudo")

end subroutine da_residual_pseudo


subroutine da_get_innov_vector_pseudo(grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   type(domain),     intent(in)    :: grid        ! Background structure 
   type(y_type),     intent(inout) :: ob          ! Observation structure.
   type(iv_type),    intent(inout) :: iv          ! O-B structure.

   integer :: n        ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_p(:,:)
   real, allocatable :: model_t(:,:)
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_pseudo")


   allocate (model_u(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (model_v(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (model_q(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (model_p(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (model_t(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))

   call da_convert_zk (iv%info(pseudo))

   call da_interp_lin_3d (grid%xb%u, iv%info(pseudo), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(pseudo), model_v)
   call da_interp_lin_3d (grid%xb%t, iv%info(pseudo), model_t)
   call da_interp_lin_3d (grid%xb%p, iv%info(pseudo), model_p)
   call da_interp_lin_3d (grid%xb%q, iv%info(pseudo), model_q)
   do n=iv%info(pseudo)%n1,iv%info(pseudo)%n2
      !---------------------------------------------------------------
      ! [3.0] Calculate observation O = B +(O-B):
      !---------------------------------------------------------------

      select case(pseudo_var(1:1))
      case ('u', 'U')
         if (ob % pseudo(n) % u > missing_r) then
             iv % pseudo(n) % u % inv = ob%pseudo(n)%u - model_u(1,n)
         else
             ob % pseudo(n) % u = model_u(1,n) + iv % pseudo(n) % u % inv
	 endif   
      case ('v', 'V')
         if (ob % pseudo(n) % v > missing_r) then
             iv % pseudo(n) % v % inv = ob%pseudo(n)%v - model_v(1,n)
         else
             ob % pseudo(n) % v = model_v(1,n) + iv % pseudo(n) % v % inv
	 endif   
      case ('t', 'T')
         if (ob % pseudo(n) % t > missing_r) then
             iv % pseudo(n) % t % inv = ob%pseudo(n)%t - model_t(1,n)
         else
             ob % pseudo(n) % t = model_t(1,n) + iv % pseudo(n) % t % inv
	 endif   
      case ('p', 'P')
         if (ob % pseudo(n) % p > missing_r) then
             iv % pseudo(n) % p % inv = ob%pseudo(n)%p - model_p(1,n)
         else
             ob % pseudo(n) % p = model_p(1,n) + iv % pseudo(n) % p % inv
	 endif   
      case ('q', 'Q')
         if (ob % pseudo(n) % q > missing_r) then
             iv % pseudo(n) % q % inv = ob%pseudo(n)%q - model_q(1,n)
         else
             ob % pseudo(n) % q = model_q(1,n) + iv % pseudo(n) % q % inv
	 endif   
      end select 

   end do

   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_q)
   deallocate (model_p)
   deallocate (model_t)
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_pseudo")

end subroutine da_get_innov_vector_pseudo


subroutine da_ao_stats_pseudo  (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! O-B
   type  (y_type), intent(in)    :: re            ! O-A

   type (stats_pseudo_type) :: stats
   integer                  :: nu, nv, nt, np, nq
   integer                  :: n
   real                     :: o_minus_b, o_minus_a, sigma_o, sigma_b
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_pseudo")

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

   stats%average = residual_pseudo1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(pseudo)%nlocal
      if (iv%info(pseudo)%proc_domain(1,n)) then
         call da_stats_calculate (n, 0, iv%pseudo(n)%u%qc, & 
            re%pseudo(n)%u, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%pseudo(n)%v%qc, & 
            re%pseudo(n)%v, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate (n, 0, iv%pseudo(n)%t%qc, & 
            re%pseudo(n)%t, nt, & 
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate (n, 0, iv%pseudo(n)%p%qc, & 
            re%pseudo(n)%p, np, & 
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate (n, 0, iv%pseudo(n)%q%qc, & 
            re%pseudo(n)%q, nq, & 
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)

         if (nu > 0) then
            o_minus_b = iv%pseudo(n)%u%inv
            o_minus_a = re%pseudo(n)%u
            sigma_o   = iv%pseudo(n)%u%error
         else if (nv > 0) then
            o_minus_b = iv%pseudo(n)%v%inv
            o_minus_a = re%pseudo(n)%v
            sigma_o   = iv%pseudo(n)%v%error
         else if (nt > 0) then
            o_minus_b = iv%pseudo(n)%t%inv
            o_minus_a = re%pseudo(n)%t
            sigma_o   = iv%pseudo(n)%t%error
         else if (np > 0) then
            o_minus_b = iv%pseudo(n)%p%inv
            o_minus_a = re%pseudo(n)%p
            sigma_o   = iv%pseudo(n)%p%error
         else if (nq > 0) then
            o_minus_b = iv%pseudo(n)%q%inv
            o_minus_a = re%pseudo(n)%q
            sigma_o   = iv%pseudo(n)%q%error
         end if

         write(stats_unit,'(/a,a1,a,e15.5)')' pseudo ', pseudo_var, ' o-b: ', o_minus_b
         write(stats_unit,' (a,a1,a,e15.5)')' pseudo ', pseudo_var, ' o-a: ', o_minus_a
         write(stats_unit,' (a,a1,a,e15.5)')' pseudo ', pseudo_var, ' sigma_o: ', sigma_o

         ! calculate equivalent sigma_b using o-a=(o-b)*sigma_o/(sigma_o+sigma_b)
         sigma_b = sqrt ((o_minus_b - o_minus_a) / o_minus_a) * sigma_o
         write(stats_unit,'(a,a1,a,e15.5)')' pseudo ', pseudo_var, ' sigma_b: ', sigma_b
      end if
   end do    

   ! do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (np)
   call da_proc_sum_int (nq)
   iv%nstats(pseudo) = nu + nv + nt + np + nq
   
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
      if (nu /= 0 .or. nv /= 0 .OR. nt /= 0 .or. np /= 0 .OR. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' O-A Diagnostics for pseudo'
         call da_print_stats_pseudo (stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_pseudo")

end subroutine da_ao_stats_pseudo


subroutine da_oi_stats_pseudo (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! O-B

   type (stats_pseudo_type)         :: stats
   integer                          :: nu, nv, nt, np, nq
   integer                          :: n

   if (trace_use_dull) call da_trace_entry("da_oi_stats_pseudo")

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
   stats%average = residual_pseudo1_type(0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(pseudo)%nlocal
      if (iv%info(pseudo)%proc_domain(1,n)) then
         call da_stats_calculate(n, 0, iv%pseudo(n)%u%qc, & 
            iv%pseudo(n)%u%inv, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate(n, 0, iv%pseudo(n)%v%qc, & 
            iv%pseudo(n)%v%inv, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
         call da_stats_calculate(n, 0, iv%pseudo(n)%t%qc, & 
            iv%pseudo(n)%t%inv, nt, & 
            stats%minimum%t, stats%maximum%t, &
            stats%average%t, stats%rms_err%t)
         call da_stats_calculate(n, 0, iv%pseudo(n)%p%qc, & 
            iv%pseudo(n)%p%inv, np, & 
            stats%minimum%p, stats%maximum%p, &
            stats%average%p, stats%rms_err%p)
         call da_stats_calculate(n, 0, iv%pseudo(n)%q%qc, & 
            iv%pseudo(n)%q%inv, nq, & 
            stats%minimum%q, stats%maximum%q, &
            stats%average%q, stats%rms_err%q)
      end if    ! end if (iv%info(pseudo)%proc_domain(1,n))
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
      if (nu /= 0 .or. nv /= 0 .OR. nt /= 0 .or. np /= 0 .OR. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' O-B Diagnostics for pseudo'
         call da_print_stats_pseudo(stats_unit, nu, nv, nt, np, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_pseudo")

end subroutine da_oi_stats_pseudo


subroutine da_print_stats_pseudo(stats_unit, nu, nv, nt, np, nq, pseudo)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                  intent(in)    :: stats_unit
   integer,                  intent(inout) :: nu, nv, nt, np, nq
   type (stats_pseudo_type), intent(in)    :: pseudo

   if (trace_use_dull) call da_trace_entry("da_print_stats_pseudo")

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
        ' Minimum(n,k): ', pseudo%minimum%u, pseudo%minimum%v, pseudo%minimum%t, &
                      pseudo%minimum%p, pseudo%minimum%q, &
        ' Maximum(n,k): ', pseudo%maximum%u, pseudo%maximum%v, pseudo%maximum%t, &
                      pseudo%maximum%p, pseudo%maximum%q
   write(unit=stats_unit, fmt='((a,4(f12.4,10x),e12.4,10x))') &
        ' Average     : ', pseudo%average%u/real(nu), pseudo%average%v/real(nv), &
                      pseudo%average%t/real(nt), pseudo%average%p/real(np), &
                      pseudo%average%q/real(nq), &
        '    RMSE     : ', sqrt(pseudo%rms_err%u/real(nu)), &
                      sqrt(pseudo%rms_err%v/real(nv)), &
                      sqrt(pseudo%rms_err%t/real(nt)), &
                      sqrt(pseudo%rms_err%p/real(np)), &
                      sqrt(pseudo%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_pseudo")

end subroutine da_print_stats_pseudo


subroutine da_transform_xtoy_pseudo(grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer :: n        ! Loop counter.

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: q(:,:)
   real, allocatable :: p(:,:)
   real, allocatable :: t(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_pseudo")

   allocate (u(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (v(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (q(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (p(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (t(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))

   call da_interp_lin_3d(grid%xa%u, iv%info(pseudo), u)
   call da_interp_lin_3d(grid%xa%v, iv%info(pseudo), v)
   call da_interp_lin_3d(grid%xa%q, iv%info(pseudo), q)
   call da_interp_lin_3d(grid%xa%p, iv%info(pseudo), p)
   call da_interp_lin_3d(grid%xa%t, iv%info(pseudo), t)
   do n=iv%info(pseudo)%n1,iv%info(pseudo)%n2
      y%pseudo(n)%u = u(1,n)
      y%pseudo(n)%v = v(1,n)
      y%pseudo(n)%q = q(1,n)
      y%pseudo(n)%p = p(1,n)
      y%pseudo(n)%t = t(1,n)
   end do

   deallocate (u)
   deallocate (v)
   deallocate (q)
   deallocate (p)
   deallocate (t)
   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_pseudo")

end subroutine da_transform_xtoy_pseudo


subroutine da_transform_xtoy_pseudo_adj(iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: q(:,:)
   real, allocatable :: p(:,:)
   real, allocatable :: t(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_pseudo_adj")

   allocate (u(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (v(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (q(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (p(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))
   allocate (t(1,iv%info(pseudo)%n1:iv%info(pseudo)%n2))

   do n=iv%info(pseudo)%n1,iv%info(pseudo)%n2
      u(1,n) = jo_grad_y%pseudo(n)%u
      v(1,n) = jo_grad_y%pseudo(n)%v
      q(1,n) = jo_grad_y%pseudo(n)%q
      p(1,n) = jo_grad_y%pseudo(n)%p
      t(1,n) = jo_grad_y%pseudo(n)%t
   end do

   call da_interp_lin_3d_adj(jo_grad_x%u, iv%info(pseudo), u)
   call da_interp_lin_3d_adj(jo_grad_x%v, iv%info(pseudo), v)
   call da_interp_lin_3d_adj(jo_grad_x%q, iv%info(pseudo), q)
   call da_interp_lin_3d_adj(jo_grad_x%p, iv%info(pseudo), p)
   call da_interp_lin_3d_adj(jo_grad_x%t, iv%info(pseudo), t)


   deallocate (u)
   deallocate (v)
   deallocate (q)
   deallocate (p)
   deallocate (t)
   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_pseudo_adj")

end subroutine da_transform_xtoy_pseudo_adj


subroutine da_calculate_grady_pseudo(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_pseudo")

   do n=1, iv%info(pseudo)%nlocal
      if (iv%pseudo(n)%u%qc < obs_qc_pointer) re%pseudo(n)%u = 0.0
      if (iv%pseudo(n)%v%qc < obs_qc_pointer) re%pseudo(n)%v = 0.0
      if (iv%pseudo(n)%t%qc < obs_qc_pointer) re%pseudo(n)%t = 0.0
      if (iv%pseudo(n)%p%qc < obs_qc_pointer) re%pseudo(n)%p = 0.0
      if (iv%pseudo(n)%q%qc < obs_qc_pointer) re%pseudo(n)%q = 0.0

      jo_grad_y%pseudo(n)%u = -re%pseudo(n)%u / (iv%pseudo(n)%u%error * iv%pseudo(n)%u%error)
      jo_grad_y%pseudo(n)%v = -re%pseudo(n)%v / (iv%pseudo(n)%v%error * iv%pseudo(n)%v%error)
      jo_grad_y%pseudo(n)%t = -re%pseudo(n)%t / (iv%pseudo(n)%t%error * iv%pseudo(n)%t%error)
      jo_grad_y%pseudo(n)%p = -re%pseudo(n)%p / (iv%pseudo(n)%p%error * iv%pseudo(n)%p%error)
      jo_grad_y%pseudo(n)%q = -re%pseudo(n)%q / (iv%pseudo(n)%q%error * iv%pseudo(n)%q%error)
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_pseudo")
  
end subroutine da_calculate_grady_pseudo



end module da_pseudo

