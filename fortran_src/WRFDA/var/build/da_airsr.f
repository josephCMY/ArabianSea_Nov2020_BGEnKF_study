












module da_airsr

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, max_error_uv, max_error_t, rootproc, &
      airsr, max_error_p,max_error_q, trace_use_dull,fails_error_max, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, anal_type_verify, kms,kme,kts,kte, &
      ob_vars, qcstat_conv_unit, fails_error_max
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_interp_lin_3d, da_to_zk, &
      da_interp_lin_3d_adj
   use da_par_util1, only : da_proc_sum_int
   use da_par_util, only : da_proc_stats_combine
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk, da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_airsr1_type
      real          :: t                        
      real          :: q                        
   end type residual_airsr1_type

   type maxmin_airsr_stats_type
      type (maxmin_type)         :: t, q
   end type maxmin_airsr_stats_type

   type stats_airsr_type
      type (maxmin_airsr_stats_type)  :: maximum, minimum
      type (residual_airsr1_type)     :: average, rms_err
   end type stats_airsr_type

contains

subroutine da_ao_stats_airsr (stats_unit, iv, re)

   !-------------------------------------------------------------------------
   ! Purpose: Compute analysis increment at AIRS retrieval locations
   !-------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O

   type (stats_airsr_type) :: stats
   integer                 :: nt, nq
   integer                 :: n, k

   if (trace_use_dull) call da_trace_entry("da_ao_stats_airsr")

   nt = 0
   nq = 0

   stats%maximum%t = maxmin_type (missing_r, 0, 0)
   stats%maximum%q = maxmin_type (missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_airsr1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(airsr)%nlocal
      if (iv%info(airsr)%proc_domain(1,n)) then
         do k=1, iv%info(airsr)%levels(n)
            call da_stats_calculate (n, k, iv%airsr(n)%t(k)%qc, & 
               re%airsr(n)%t(k), nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate (n, k, iv%airsr(n)%q(k)%qc, & 
               re%airsr(n)%q(k), nq, &
               stats%minimum%q, stats%maximum%q, &
               stats%average%q, stats%rms_err%q)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nt)
   call da_proc_sum_int (nq)
   iv%nstats(airsr) = nt + nq

   call da_proc_stats_combine(stats%average%t, stats%rms_err%t, &
      stats%minimum%t%value, stats%maximum%t%value, &
      stats%minimum%t%n, stats%maximum%t%n, &
      stats%minimum%t%l, stats%maximum%t%l)
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
      stats%minimum%q%value, stats%maximum%q%value, &
      stats%minimum%q%n, stats%maximum%q%n, &
      stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nt /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for airs retrievals'
         call da_print_stats_airsr(stats_unit, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_airsr")

end subroutine da_ao_stats_airsr


subroutine da_jo_and_grady_airsr(iv, re, jo, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates cost function and its gradient at all AIRS 
   ! retrieval locations   
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv         ! Innovation vector.
   type (y_type),  intent(in)    :: re         ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y  ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo         ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_airsr")

   jo % airsr_t = 0.0
   jo % airsr_q = 0.0

   do n=1, iv%info(airsr)%nlocal
      do k=1, iv%info(airsr)%levels(n)
         jo_grad_y%airsr(n)%t(k) = -re%airsr(n)%t(k) / (iv%airsr(n)%t(k)%error * iv%airsr(n)%t(k)%error)
         jo_grad_y%airsr(n)%q(k) = -re%airsr(n)%q(k) / (iv%airsr(n)%q(k)%error * iv%airsr(n)%q(k)%error)
      end do

      if (iv%info(airsr)%proc_domain(1,n)) then
         do k=1, iv%info(airsr)%levels(n)
            jo % airsr_t = jo % airsr_t - re%airsr(n)%t(k) * jo_grad_y%airsr(n)%t(k)
            jo % airsr_q = jo % airsr_q - re%airsr(n)%q(k) * jo_grad_y%airsr(n)%q(k)
         end do
      end if
   end do

   jo % airsr_t = 0.5 * jo % airsr_t
   jo % airsr_q = 0.5 * jo % airsr_q

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_airsr")

 end subroutine da_jo_and_grady_airsr


subroutine da_jo_airsr_tq(iv, re, jo_grad_y, jo)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv         ! Innovation vector.
   type (y_type),  intent(in)    :: re         ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y  ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo         ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_airsr_tq")

   do n=1, iv%info(airsr)%nlocal
      if (iv%info(airsr)%proc_domain(1,n)) then
         do k=1, iv%info(airsr)%levels(n)
            jo % airsr_t = jo % airsr_t - re%airsr(n)%t(k) * jo_grad_y%airsr(n)%t(k)
            jo % airsr_q = jo % airsr_q - re%airsr(n)%q(k) * jo_grad_y%airsr(n)%q(k)
         end do
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_jo_airsr_tq")

end subroutine da_jo_airsr_tq


subroutine da_residual_airsr(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !----------------------------------------------------------------------
   ! Purpose: Calculates residual at AIRS retrieval locations                   
   !----------------------------------------------------------------------

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

   if (trace_use_dull) call da_trace_entry("da_residual_airsr")

   n_obs_bad % t % num = number_type(0, 0, 0)
   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(airsr)%nlocal
      do k=1, iv%info(airsr)%levels(n)
        np_available = np_available + 4
        re%airsr(n)%t(k) = da_residual(n, k, y%airsr(n)%t(k), iv%airsr(n)%t(k), n_obs_bad % t)
        re%airsr(n)%q(k) = da_residual(n, k, y%airsr(n)%q(k), iv%airsr(n)%q(k), n_obs_bad % q)
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % t % num % miss + n_obs_bad % q % num % miss
   np_bad_data = np_bad_data + n_obs_bad % t % num % bad  + n_obs_bad % q % num % bad
   np_obs_used = np_obs_used + n_obs_bad % t % num % use  + n_obs_bad % q % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_airsr")

end subroutine da_residual_airsr


subroutine da_oi_stats_airsr (stats_unit, iv)

  !----------------------------------------------------------------------
  ! Purpose: Prints out innov diagnostics for AIRS retrievals  
  !----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_airsr_type) :: stats
   integer                 :: nt, nq
   integer                 :: n, k

   if (trace_use_dull) call da_trace_entry("da_oi_stats_airsr")

   nt = 0
   nq = 0

   stats%maximum%t = maxmin_type(missing_r, 0, 0)
   stats%maximum%q = maxmin_type(missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_airsr1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(airsr)%nlocal
      if (iv%info(airsr)%proc_domain(1,n)) then
         do k=1, iv%info(airsr)%levels(n)
            call da_stats_calculate(iv%info(airsr)%obs_global_index(n), &
               k, iv%airsr(n)%t(k)%qc, &
               iv%airsr(n)%t(k)%inv, nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate(iv%info(airsr)%obs_global_index(n), &
               k, iv%airsr(n)%q(k)%qc, &
               iv%airsr(n)%q(k)%inv, nq, &
               stats%minimum%q, stats%maximum%q, &
               stats%average%q, stats%rms_err%q)
         end do
      end if    ! end if (iv%info(airsr)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nt)
   call da_proc_sum_int(nq)

   call da_proc_stats_combine(stats%average%t, stats%rms_err%t, &
      stats%minimum%t%value, stats%maximum%t%value, &
      stats%minimum%t%n, stats%maximum%t%n, &
      stats%minimum%t%l, stats%maximum%t%l)
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
      stats%minimum%q%value, stats%maximum%q%value, &
      stats%minimum%q%n, stats%maximum%q%n, &
      stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nt /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for airs retrievals'
         call da_print_stats_airsr(stats_unit, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_airsr")

end subroutine da_oi_stats_airsr


subroutine da_print_stats_airsr(stats_unit, nt, nq, airsr)

   !-----------------------------------------------------------------------
   ! Purpose: Prints out table for AIRS diagnostics
   !----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nt, nq
   type (stats_airsr_type), intent(in)    :: airsr

   if (trace_use_dull) call da_trace_entry("da_print_stats_airsr")
   
   write(unit=stats_unit, fmt='(5a/)') &
      '   var             ', &
      't (K)       n    k    ', 'q (kg/kg)   n    k'

   write(unit=stats_unit, fmt='(a,i16,i22)') &
      '  Number: ', nt, nq

   if (nt < 1) nt = 1
   if (nq < 1) nq = 1
   
   write(unit=stats_unit, fmt='((a,2(f12.4,2i5)))') &
      ' Minimum(n,k): ', airsr%minimum%t, airsr%minimum%q, &
      ' Maximum(n,k): ', airsr%maximum%t, airsr%maximum%q
   write(unit=stats_unit, fmt='((a,f12.4,10x,e12.4,10x))')  &
      ' Average     : ', airsr%average%t/real(nt), airsr%average%q/real(nq),      &
      '    RMSE     : ', sqrt(airsr%rms_err%t/real(nt)), sqrt(airsr%rms_err%q/real(nq))

   if (trace_use_dull) call da_trace_exit("da_print_stats_airsr")

end subroutine da_print_stats_airsr


subroutine da_transform_xtoy_airsr (grid, iv, y)

   !-------------------------------------------------------------------------
   ! Purpose: Does transforms from model space to AIRS locations
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-------------------------------------------------------------------------

   implicit none

   type (domain),     intent(in)    :: grid     ! gridded analysis increment.
   type (iv_type),    intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n  ! Loop counter.

   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_airsr")

   allocate (t(iv%info(airsr)%max_lev,iv%info(airsr)%n1:iv%info(airsr)%n2))
   allocate (q(iv%info(airsr)%max_lev,iv%info(airsr)%n1:iv%info(airsr)%n2))
  
   call da_interp_lin_3d (grid%xa%t, iv%info(airsr), t)
   call da_interp_lin_3d (grid%xa%q, iv%info(airsr), q)

   do n=iv%info(airsr)%n1,iv%info(airsr)%n2
      y%airsr(n)%t(:) = t(1:size(y%airsr(n)%t),n)
      y%airsr(n)%q(:) = q(1:size(y%airsr(n)%q),n)
   end do

   deallocate (t)
   deallocate (q)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_airsr")

end subroutine da_transform_xtoy_airsr


subroutine da_transform_xtoy_airsr_adj(iv, jo_grad_y, jo_grad_x)

   !----------------------------------------------------------------------
   ! Purpose: Does adjoint computation at AIRS retrieval locations
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n  ! Loop counter.

   real, allocatable :: t(:,:)
   real, allocatable :: q(:,:)           

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_airsr_adj")

   allocate (t(iv%info(airsr)%max_lev,iv%info(airsr)%n1:iv%info(airsr)%n2))
   allocate (q(iv%info(airsr)%max_lev,iv%info(airsr)%n1:iv%info(airsr)%n2))

   do n=iv%info(airsr)%n1,iv%info(airsr)%n2
      t(1:size(jo_grad_y%airsr(n)%t),n)  = jo_grad_y%airsr(n)%t
      q(1:size(jo_grad_y%airsr(n)%q),n)  = jo_grad_y%airsr(n)%q
   end do

   ! [1.1] Adjoint feedback from Y to X for u and v:

   call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(airsr), t)
   call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(airsr), q)

   deallocate (t)
   deallocate (q)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_airsr_adj")

end subroutine da_transform_xtoy_airsr_adj


subroutine da_check_max_iv_airsr(iv, it,num_qcstat_conv)     

   !-------------------------------------------------------------------------
   ! Purpose: Applies max error check on AIRS retrievals
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: k,n, ipr
   logical :: failed
   
   if (trace_use_dull) call da_trace_entry("da_check_max_iv_airsr")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(airsr)%n1,iv%info(airsr)%n2
      do k = 1, iv%info(airsr)%levels(n)
        call da_get_print_lvl(iv%airsr(n)%p(k),ipr)
        failed = .false.
        if( iv%airsr(n)%t(k)%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it,iv%info(airsr), n, iv%airsr(n)%t(k), max_error_t ,failed)
        if( iv%info(airsr)%proc_domain(k,n) ) then
          num_qcstat_conv(1,airsr,3,ipr) = num_qcstat_conv(1,airsr,3,ipr) + 1
         if(failed)then
          num_qcstat_conv(2,airsr,3,ipr) = num_qcstat_conv(2,airsr,3,ipr) + 1
        write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
        'airsr',ob_vars(3),iv%info(airsr)%lat(k,n),iv%info(airsr)%lon(k,n),0.01*iv%airsr(n)%p(k)
         endif
        endif
        endif

        failed = .false.
        if( iv%airsr(n)%q(k)%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it,iv%info(airsr), n, iv%airsr(n)%q(k), max_error_q ,failed)
         if( iv%info(airsr)%proc_domain(k,n) ) then
           num_qcstat_conv(1,airsr,4,ipr) = num_qcstat_conv(1,airsr,4,ipr) + 1
         if(failed)then
          num_qcstat_conv(2,airsr,4,ipr) = num_qcstat_conv(2,airsr,4,ipr) + 1
        write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
        'airsr',ob_vars(4),iv%info(airsr)%lat(k,n),iv%info(airsr)%lon(k,n),0.01*iv%airsr(n)%p(k)
         endif
        endif
        endif
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_max_iv_airsr")
end subroutine da_check_max_iv_airsr
subroutine da_get_innov_vector_airsr( it,num_qcstat_conv, grid, ob, iv)

   !----------------------------------------------------------------------
   ! Purpose:   a) Rcomputes innovation vecrot at
   !                   AIRS retrieval locations 
   !                b) Does quality control check on innovation vector
   !                   if required.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it       ! External iteration.
   type(domain),     intent(in)    :: grid       ! first guess state.
   type(y_type),     intent(inout) :: ob       ! Observation structure.
   type(iv_type),    intent(inout) :: iv       ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: n  ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs ! Number of obs levels.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.

   ! real    :: model_h(1:max_ob_levels)  ! Model value h at ob location.

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.

   real, allocatable :: model_q(:,:)  ! Model value q at ob location.
   real, allocatable :: model_t(:,:)  ! Model value t at ob location.

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_airsr")

   allocate (model_t(iv%info(airsr)%max_lev,iv%info(airsr)%n1:iv%info(airsr)%n2))
   allocate (model_q(iv%info(airsr)%max_lev,iv%info(airsr)%n1:iv%info(airsr)%n2))

   model_t(:,:) = 0.0
   model_q(:,:) = 0.0

   if ( it > 1 ) then
      do n = iv%info(airsr)%n1, iv%info(airsr)%n2
         do k = 1, iv%info(airsr)%levels(n)
            if (iv%airsr(n)%t(k)%qc == fails_error_max) iv%airsr(n)%t(k)%qc = 0
            if (iv%airsr(n)%q(k)%qc == fails_error_max) iv%airsr(n)%q(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(airsr)%n1, iv%info(airsr)%n2
      num_levs = iv%info(airsr)%levels(n)

      if (num_levs < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      i   = iv%info(airsr)%i(1,n)
      j   = iv%info(airsr)%j(1,n)
      dx  = iv%info(airsr)%dx(1,n)
      dy  = iv%info(airsr)%dy(1,n)
      dxm = iv%info(airsr)%dxm(1,n)
      dym = iv%info(airsr)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j,k) + dx*grid%xb%h(i+1,j,k)) + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
         v_p(k) = dym*(dxm*grid%xb%p(i,j,k) + dx*grid%xb%p(i+1,j,k)) + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
      end do

      do k=1, num_levs
         if (iv%airsr(n)%p(k) > 1.0) then
            call da_to_zk(iv%airsr(n)%p(k), v_p, v_interp_p, iv%info(airsr)%zk(k,n))
         else if (iv%airsr(n)%h(k) > 0.0) then
            call da_to_zk(iv%airsr(n)%h(k), v_h, v_interp_h, iv%info(airsr)%zk(k,n))
         end if

         if (iv%info(airsr)%zk(k,n) < 0.0 .and. .not. anal_type_verify) then
            iv%airsr(n)%t(k)%qc = missing_data
            iv%airsr(n)%q(k)%qc = missing_data
         end if
      end do
   end do

   call da_convert_zk (iv%info(airsr))

   ! [1.2] Interpolate horizontally to ob:
   call da_interp_lin_3d (grid%xb%t, iv%info(airsr), model_t)
   call da_interp_lin_3d (grid%xb%q, iv%info(airsr), model_q)

   do n=iv%info(airsr)%n1,iv%info(airsr)%n2
      !------------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !------------------------------------------------------------------------

      do k = 1, iv%info(airsr)%levels(n)
         iv%airsr(n)%t(k)%inv = 0.0
         iv%airsr(n)%q(k)%inv = 0.0

         !------------------------------------------------------------------------
         ! [3.0] Interpolation:
         !------------------------------------------------------------------------

         if (ob%airsr(n)%t(k) > missing_r .AND. iv%airsr(n)%t(k)%qc >= obs_qc_pointer) then
           iv%airsr(n)%t(k)%inv = ob%airsr(n)%t(k) - model_t(k,n)
         end if

         if (ob%airsr(n)%q(k) > missing_r .AND. iv%airsr(n)%q(k)%qc >= obs_qc_pointer) then
            iv%airsr(n)%q(k)%inv = ob%airsr(n)%q(k) - model_q(k,n)
         end if
      end do
   end do

   !------------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_airsr(iv, it, num_qcstat_conv)     
 
   deallocate (model_t)
   deallocate (model_q)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_airsr")

end subroutine da_get_innov_vector_airsr


subroutine da_calculate_grady_airsr(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on AIRS re-vector             
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer :: n, k
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_airsr")

   do n=1, iv%info(airsr)%nlocal
      do k=1, iv%info(airsr)%levels(n)
         if (iv%airsr(n)%t(k)%qc < obs_qc_pointer) re%airsr(n)%t(k) = 0.0
         if (iv%airsr(n)%q(k)%qc < obs_qc_pointer) re%airsr(n)%q(k) = 0.0

         jo_grad_y%airsr(n)%t(k) = -re%airsr(n)%t(k) / (iv%airsr(n)%t(k)%error * iv%airsr(n)%t(k)%error)
         jo_grad_y%airsr(n)%q(k) = -re%airsr(n)%q(k) / (iv%airsr(n)%q(k)%error * iv%airsr(n)%q(k)%error)
       end do
    end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_airsr")

end subroutine da_calculate_grady_airsr



end module da_airsr

