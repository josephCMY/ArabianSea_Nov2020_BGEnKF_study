












module da_polaramv

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, trace_use_dull, &
      missing, max_error_uv, max_error_t, rootproc, kms,kme,kts,kte, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv, fails_error_max,  &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, polaramv, anal_type_verify, &
      position_lev_dependant, qcstat_conv_unit,ob_vars, &
      convert_fd2uv,convert_uv2fd,max_error_spd,max_error_dir,max_omb_spd,max_omb_dir,pi,qc_rej_both, &
      wind_sd_polaramv, wind_stats_sd 
   use da_grid_definitions, only : da_ffdduv,da_ffdduv_model, da_ffdduv_diagnose
   use da_physics, only : da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_interp_lin_3d, da_to_zk, &
      da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk, da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_polaramv1_type
      real          :: u                        
      real          :: v                        
   end type residual_polaramv1_type

   type maxmin_polaramv_stats_type
      type (maxmin_type)         :: u, v, t, q
   end type maxmin_polaramv_stats_type

   type stats_polaramv_type
      type (maxmin_polaramv_stats_type)  :: maximum, minimum
      type (residual_polaramv1_type)     :: average, rms_err
   end type stats_polaramv_type


contains

subroutine da_ao_stats_polaramv (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type  (y_type), intent(in)    :: re            ! A - O
   type(y_type),   intent (in)   :: ob            ! Observation structure.

   type (stats_polaramv_type) :: stats
   integer                    :: nu, nv
   integer                    :: n, k
   real                       :: u_inc, v_inc, u_obs, v_obs
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_polaramv")

   nu = 0
   nv = 0
   
   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_polaramv1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(polaramv)%nlocal
      if (iv%info(polaramv)%proc_domain(1,n)) then
         do k=1, iv%info(polaramv)%levels(n)

            u_inc = re%polaramv(n)%u(k)
            v_inc = re%polaramv(n)%v(k)
            u_obs = ob%polaramv(n)%u(k)
            v_obs = ob%polaramv(n)%v(k)

            if (.not. wind_sd_polaramv .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%polaramv(n)%u(k)%qc, iv%polaramv(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_polaramv .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%polaramv(n)%u(k)%qc, iv%polaramv(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate (n, 0, iv%polaramv(n)%u(k)%qc, & 
               u_inc, nu, & 
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, 0, iv%polaramv(n)%v(k)%qc, & 
               v_inc, nv, & 
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   iv%nstats(polaramv) = nu + nv
   
   call da_proc_stats_combine(stats%average%u, stats%rms_err%u, &
      stats%minimum%u%value, stats%maximum%u%value, &
      stats%minimum%u%n, stats%maximum%u%n, &
      stats%minimum%u%l, stats%maximum%u%l)
   call da_proc_stats_combine(stats%average%v, stats%rms_err%v, &
      stats%minimum%v%value, stats%maximum%v%value, &
      stats%minimum%v%n, stats%maximum%v%n, &
      stats%minimum%v%l, stats%maximum%v%l)
      
   if (rootproc) then   
      if (nu /= 0 .or. nv /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for polaramv'
         call da_print_stats_polaramv (stats_unit, nu, nv, stats)
      end if
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_polaramv")

end subroutine da_ao_stats_polaramv  


subroutine da_jo_and_grady_polaramv( iv, re, jo, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose:  Calculates Cost function and Gradient for Polar AMVs Obs
   !----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv          ! Innovation vector.
   type(y_type),  intent(in)    :: re          ! Residual vector.
   type(y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type(jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_polaramv")

   jo % polaramv_u = 0.0
   jo % polaramv_v = 0.0

   do n=1, iv%info(polaramv)%nlocal
      do k=1, iv%info(polaramv)%levels(n)
         jo_grad_y%polaramv(n)%u(k) = -re%polaramv(n)%u(k) / &
            ( iv%polaramv(n)%u(k)%error * iv%polaramv(n)%u(k)%error)
         jo_grad_y%polaramv(n)%v(k) = -re%polaramv(n)%v(k) / &
            ( iv%polaramv(n)%v(k)%error * iv%polaramv(n)%v(k)%error)
      end do

      do k=1, iv%info(polaramv)%levels(n)
         if (iv%info(polaramv)%proc_domain(k,n)) then
            jo % polaramv_u = jo % polaramv_u - re%polaramv(n)%u(k) * jo_grad_y%polaramv(n)%u(k)
            jo % polaramv_v = jo % polaramv_v - re%polaramv(n)%v(k) * jo_grad_y%polaramv(n)%v(k)
         end if
      end do
   end do

   jo % polaramv_u = 0.5 * jo % polaramv_u
   jo % polaramv_v = 0.5 * jo % polaramv_v

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_polaramv")
     
end subroutine da_jo_and_grady_polaramv


subroutine da_residual_polaramv(iv, y, re,np_missing, np_bad_data,np_obs_used, np_available)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------------

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

   if (trace_use_dull) call da_trace_entry("da_residual_polaramv")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)

   do n=1, iv%info(polaramv)%nlocal
      do k=1, iv%info(polaramv)%levels(n)
         np_available = np_available + 2
         re%polaramv(n)%u(k) = da_residual(n, 0, y%polaramv(n)%u(k), iv%polaramv(n)%u(k), n_obs_bad % u)
         re%polaramv(n)%v (k)= da_residual(n, 0, y%polaramv(n)%v(k), iv%polaramv(n)%v(k), n_obs_bad % v)            
     end do
   end do

   np_missing  = np_missing  + n_obs_bad % u % num % miss + n_obs_bad % v % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad + n_obs_bad % v % num % bad
   np_obs_used = np_obs_used + n_obs_bad % u % num % use + n_obs_bad % v % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_polaramv")

end subroutine da_residual_polaramv


subroutine da_oi_stats_polaramv(stats_unit, iv, ob)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates (Obs - Analysis) statistics for Polar AMVs
   !-------------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! OI
   type(y_type),   intent (in)      :: ob            ! Observation structure.

   type (stats_polaramv_type)       :: stats
   integer                          :: nu, nv
   integer                          :: n, k
   real                             :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_polaramv")

   nu = 0
   nv = 0
   
   stats%maximum%u = maxmin_type(missing_r, 0, 0)
   stats%maximum%v = maxmin_type(missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_polaramv1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(polaramv)%nlocal
      if (iv%info(polaramv)%proc_domain(1,n)) then
         do k=1, iv%info(polaramv)%levels(n)

            u_inv = iv%polaramv(n)%u(k)%inv
            v_inv = iv%polaramv(n)%v(k)%inv
            u_obs = ob%polaramv(n)%u(k)
            v_obs = ob%polaramv(n)%v(k)

            if (.not. wind_sd_polaramv .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%polaramv(n)%u(k)%qc, iv%polaramv(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_polaramv .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%polaramv(n)%u(k)%qc, iv%polaramv(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate(iv%info(polaramv)%obs_global_index(n), &
               0, iv%polaramv(n)%u(k)%qc, &
               u_inv, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(polaramv)%obs_global_index(n), &
               0, iv%polaramv(n)%v(k)%qc, &
               v_inv, nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nu)
   call da_proc_sum_int(nv)
   call da_proc_stats_combine(stats%average%u, stats%rms_err%u, &
      stats%minimum%u%value, stats%maximum%u%value, &
      stats%minimum%u%n, stats%maximum%u%n, &
      stats%minimum%u%l, stats%maximum%u%l)
   call da_proc_stats_combine(stats%average%v, stats%rms_err%v, &
      stats%minimum%v%value, stats%maximum%v%value, &
      stats%minimum%v%n, stats%maximum%v%n, &
      stats%minimum%v%l, stats%maximum%v%l)
      
   if (rootproc) then
      if (nu /= 0 .or. nv /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for polaramv'
         call da_print_stats_polaramv(stats_unit, nu, nv, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_polaramv")


end subroutine da_oi_stats_polaramv


subroutine da_print_stats_polaramv(stats_unit, nu, nv, polaramv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: stats_unit
   integer,                    intent(inout) :: nu, nv
   type (stats_polaramv_type), intent(in)    :: polaramv

   if (trace_use_dull) call da_trace_entry("da_print_stats_polaramv")

   write(unit=stats_unit, fmt='(a/)') &
      '   var             u (m/s)     n    k    v (m/s)     n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
      '  Number: ', nu, nv

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1

   write(unit=stats_unit, fmt='((a,2(f12.4,2i5)))') &
      ' Minimum(n,k): ', polaramv%minimum%u, polaramv%minimum%v, &
      ' Maximum(n,k): ', polaramv%maximum%u, polaramv%maximum%v

   write(unit=stats_unit, fmt='((a,2(f12.4,10x)))') &
      ' Average     : ', polaramv%average%u/real(nu), &
      polaramv%average%v/real(nv), &
      '    RMSE     : ', sqrt(polaramv%rms_err%u/real(nu)), &
                      sqrt(polaramv%rms_err%v/real(nv))

   if (trace_use_dull) call da_trace_exit("da_print_stats_polaramv")

end subroutine da_print_stats_polaramv


subroutine da_transform_xtoy_polaramv (grid, iv, y)

   !----------------------------------------------------------------------
   ! Purpose: X to Y Transform operator for Polar AMV's               
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer           :: n, k ! Loop counter.
   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_polaramv")

   allocate (u(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (v(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (ub(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (vb(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))

   call da_interp_lin_3d (grid%xa%u, iv%info(polaramv), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(polaramv), v)
   call da_interp_lin_3d (grid%xb%u, iv%info(polaramv), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(polaramv), vb)

   do n=iv%info(polaramv)%n1,iv%info(polaramv)%n2
      do k = 1, iv%info(polaramv)%levels(n)

         if(wind_sd_polaramv) then
            call da_uv_to_sd_lin(y%polaramv(n)%u(k),y%polaramv(n)%v(k),u(k,n),v(k,n),ub(k,n),vb(k,n))
         else
            y%polaramv(n)%u(k) = u(k,n)
            y%polaramv(n)%v(k) = v(k,n)
         end if

      end do
   end do

   deallocate (u)
   deallocate (v)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_polaramv")

end subroutine da_transform_xtoy_polaramv


subroutine da_transform_xtoy_polaramv_adj (grid, iv, jo_grad_y, jo_grad_x)

   !-------------------------------------------------------------------------
   ! Purpose: X to Y Transpose operator for Polar AMVs 
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-------------------------------------------------------------------------

   implicit none
   type(domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer           :: n, k     ! Loop counter.
   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_polaramv_adj")

   allocate (u(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (v(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (ub(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (vb(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))

   call da_interp_lin_3d (grid%xb%u, iv%info(polaramv), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(polaramv), vb)

   do n=iv%info(polaramv)%n1,iv%info(polaramv)%n2
       do k = 1, iv%info(polaramv)%levels(n)

         if(wind_sd_polaramv) then
            call da_uv_to_sd_adj(jo_grad_y%polaramv(n)%u(k), &
                                 jo_grad_y%polaramv(n)%v(k), u(k,n), v(k,n), ub(k,n), vb(k,n))
         else
             u(k,n) = jo_grad_y%polaramv(n)%u(k)
             v(k,n) = jo_grad_y%polaramv(n)%v(k)
         end if

      end do
   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(polaramv), u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(polaramv), v)

   deallocate (u)
   deallocate (v)
   deallocate (ub)
   deallocate (vb)
   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_polaramv_adj")

end subroutine da_transform_xtoy_polaramv_adj


subroutine da_check_max_iv_polaramv(iv,it,num_qcstat_conv)   

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type),    intent(inout) :: iv
   integer,          intent(in)    :: it      ! Outer iteration
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer  :: k,n, ipr
   logical  :: failed,failed1,failed2
   
   if (trace_use_dull) call da_trace_entry("da_check_max_iv_polaramv")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(polaramv)%n1,iv%info(polaramv)%n2
      do k = 1, iv%info(polaramv)%levels(n)
        call da_get_print_lvl(iv%polaramv(n)%p(k),ipr)
         if(.not. qc_rej_both)then
             if(wind_sd_polaramv)then
               failed=.false.
               if( iv%polaramv(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%u(k), max_error_spd,failed)
                   if( iv%info(polaramv)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,polaramv,1,ipr) = num_qcstat_conv(1,polaramv,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,polaramv,1,ipr) = num_qcstat_conv(2,polaramv,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'polaramv',ob_vars(1),iv%info(polaramv)%lat(k,n),iv%info(polaramv)%lon(k,n),0.01*iv%polaramv(n)%p(k)
                       end if
                   end if
                end if

                failed=.false.
                if( iv%polaramv(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%v(k), max_error_dir,failed)
                    if( iv%info(polaramv)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,polaramv,2,ipr) = num_qcstat_conv(1,polaramv,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,polaramv,2,ipr) = num_qcstat_conv(2,polaramv,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'polaramv',ob_vars(2),iv%info(polaramv)%lat(k,n),iv%info(polaramv)%lon(k,n),0.01*iv%polaramv(n)%p(k)
                        end if
                    end if
                end if

             else

                failed=.false.
                if( iv%polaramv(n)%u(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%u(k), max_error_uv,failed)
                    if( iv%info(polaramv)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,polaramv,1,ipr) = num_qcstat_conv(1,polaramv,1,ipr) + 1
                        if(failed) then
                           num_qcstat_conv(2,polaramv,1,ipr) = num_qcstat_conv(2,polaramv,1,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'polaramv',ob_vars(1),iv%info(polaramv)%lat(k,n),iv%info(polaramv)%lon(k,n),0.01*iv%polaramv(n)%p(k)
                        end if
                    end if
                end if

                failed=.false.
                if( iv%polaramv(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%v(k), max_error_uv,failed)
                    if( iv%info(polaramv)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,polaramv,2,ipr) = num_qcstat_conv(1,polaramv,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,polaramv,2,ipr) = num_qcstat_conv(2,polaramv,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'polaramv',ob_vars(2),iv%info(polaramv)%lat(k,n),iv%info(polaramv)%lon(k,n),0.01*iv%polaramv(n)%p(k)
                        end if
                    end if
                 end if
             end if

             if(wind_sd_polaramv)then
                if(iv%polaramv(n)%u(k)%qc == fails_error_max .or. abs(iv%polaramv(n)%u(k)%inv) >= max_omb_spd) then
                   iv%polaramv(n)%u(k)%qc = fails_error_max
                   iv%polaramv(n)%u(k)%inv = 0.0
                endif
                if(iv%polaramv(n)%v(k)%qc == fails_error_max .or. abs(iv%polaramv(n)%v(k)%inv) >= max_omb_dir) then
                   iv%polaramv(n)%v(k)%qc = fails_error_max
                   iv%polaramv(n)%v(k)%inv = 0.0
                endif
             endif

          else


         failed1=.false.
         failed2=.false.

         if( iv%polaramv(n)%v(k)%qc >= obs_qc_pointer .or. iv%polaramv(n)%u(k)%qc >= obs_qc_pointer )  then
             if(wind_sd_polaramv)then
                call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%u(k), max_error_spd,failed1)
                call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%v(k), max_error_dir,failed2)
             else
                call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%u(k), max_error_uv,failed1)
                call da_max_error_qc (it,iv%info(polaramv), n, iv%polaramv(n)%v(k), max_error_uv,failed2)
             endif

             if( iv%info(polaramv)%proc_domain(k,n) ) then
                 num_qcstat_conv(1,polaramv,1,ipr) = num_qcstat_conv(1,polaramv,1,ipr) + 1
                 num_qcstat_conv(1,polaramv,2,ipr) = num_qcstat_conv(1,polaramv,2,ipr) + 1

                 if(failed1 .or. failed2) then
                    num_qcstat_conv(2,polaramv,1,ipr) = num_qcstat_conv(2,polaramv,1,ipr) + 1
                    write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'polaramv',ob_vars(1),iv%info(polaramv)%lat(k,n),iv%info(polaramv)%lon(k,n),0.01*iv%polaramv(n)%p(k)
                    num_qcstat_conv(2,polaramv,2,ipr) = num_qcstat_conv(2,polaramv,2,ipr) + 1
                    write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'polaramv',ob_vars(2),iv%info(polaramv)%lat(k,n),iv%info(polaramv)%lon(k,n),0.01*iv%polaramv(n)%p(k)
                 end if
             end if
          end if

	  if(wind_sd_polaramv)then
             if(iv%polaramv(n)%u(k)%qc == fails_error_max .or. iv%polaramv(n)%v(k)%qc == fails_error_max .or. &
                abs(iv%polaramv(n)%v(k)%inv) >= max_omb_dir .or. abs(iv%polaramv(n)%u(k)%inv) >= max_omb_spd )then
                iv%polaramv(n)%u(k)%qc = fails_error_max
                iv%polaramv(n)%v(k)%qc = fails_error_max
                iv%polaramv(n)%u(k)%inv = 0.0
                iv%polaramv(n)%v(k)%inv = 0.0
             endif
          else
             if(iv%polaramv(n)%u(k)%qc == fails_error_max .or. iv%polaramv(n)%v(k)%qc == fails_error_max ) then
                iv%polaramv(n)%u(k)%qc = fails_error_max
                iv%polaramv(n)%v(k)%qc = fails_error_max
                iv%polaramv(n)%u(k)%inv = 0.0
                iv%polaramv(n)%v(k)%inv = 0.0
              endif
          endif
       endif
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_max_iv_polaramv")

end subroutine da_check_max_iv_polaramv      
subroutine da_get_innov_vector_polaramv( it, num_qcstat_conv, grid, ob, iv)

   !----------------------------------------------------------------------
   ! Purpose: Calculates innovation vector, does QC for polaramv
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid      ! first guess state.
   type(y_type),     intent(in)    :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
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

   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.
   real    :: speed, direction
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_polaramv")

   allocate (model_u(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))
   allocate (model_v(iv%info(polaramv)%max_lev,iv%info(polaramv)%n1:iv%info(polaramv)%n2))

   model_u(:,:) = 0.0
   model_v(:,:) = 0.0

   if( it > 1 ) then
      do n=iv%info(polaramv)%n1, iv%info(polaramv)%n2
         do k=1, iv%info(polaramv)%levels(n)
            if (iv%polaramv(n)%u(k)%qc == fails_error_max) iv%polaramv(n)%u(k)%qc = 0
            if (iv%polaramv(n)%v(k)%qc == fails_error_max) iv%polaramv(n)%v(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(polaramv)%n1, iv%info(polaramv)%n2
      if (iv%info(polaramv)%levels(n) < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      if (position_lev_dependant) then
         i(:)   = iv%info(polaramv)%i(:,n)
         j(:)   = iv%info(polaramv)%j(:,n)
         dx(:)  = iv%info(polaramv)%dx(:,n)
         dy(:)  = iv%info(polaramv)%dy(:,n)
         dxm(:) = iv%info(polaramv)%dxm(:,n)
         dym(:) = iv%info(polaramv)%dym(:,n)
         do k=kts,kte
            v_p(k) = dym(k)*(dxm(k)*grid%xb%p(i(k),j(k),k) + dx(k)*grid%xb%p(i(k)+1,j(k),k)) &
               + dy(k) *(dxm(k)*grid%xb%p(i(k),j(k)+1,k) + dx(k)*grid%xb%p(i(k)+1,j(k)+1,k))
         end do
      else
         i(1)   = iv%info(polaramv)%i(1,n)
         j(1)   = iv%info(polaramv)%j(1,n)
         dx(1)  = iv%info(polaramv)%dx(1,n)
         dy(1)  = iv%info(polaramv)%dy(1,n)
         dxm(1) = iv%info(polaramv)%dxm(1,n)
         dym(1) = iv%info(polaramv)%dym(1,n)

         v_p(kts:kte) = dym(1) * (dxm(1)*grid%xb%p(i(1),j(1),kts:kte)   + dx(1)*grid%xb%p(i(1)+1,j(1),kts:kte)) &
                       + dy(1) * (dxm(1)*grid%xb%p(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%p(i(1)+1,j(1)+1,kts:kte))
      end if

      do k=1, iv%info(polaramv)%levels(n)
         if (iv%polaramv(n)%p(k) > 1.0) then
            call da_to_zk (iv%polaramv(n)%p(k), v_p,v_interp_p, iv%info(polaramv)%zk(k,n))
         end if
      end do
   end do

   call da_convert_zk(iv%info(polaramv))

   if (.not. anal_type_verify) then
      do n=iv%info(polaramv)%n1,iv%info(polaramv)%n2
         do k=1, iv%info(polaramv)%levels(n)
            if (iv%info(polaramv)%zk(k,n) < 0.0) then
               iv%polaramv(n)%u(k)%qc = missing_data
               iv%polaramv(n)%v(k)%qc = missing_data
            end if
         end do
      end do
   end if

   ! [1.2] Interpolate horizontally to ob:

   call da_interp_lin_3d (grid%xb%u, iv%info(polaramv), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(polaramv), model_v)

   do n=iv%info(polaramv)%n1,iv%info(polaramv)%n2
      do k = 1, iv%info(polaramv)%levels(n)

          if (wind_sd_polaramv) then
              call da_ffdduv_model (speed,direction,model_u(k,n), model_v(k,n), convert_uv2fd)

              if (ob%polaramv(n)%u(k) > missing_r .AND. iv%polaramv(n)%u(k)%qc >= obs_qc_pointer) then
                  iv%polaramv(n)%u(k)%inv = ob%polaramv(n)%u(k) - speed
              end if

              if (ob%polaramv(n)%v(k) > missing_r .AND. iv%polaramv(n)%v(k)%qc >= obs_qc_pointer) then
                  iv%polaramv(n)%v(k)%inv = ob%polaramv(n)%v(k) - direction
                  if (iv%polaramv(n)%v(k)%inv > 180.0 ) iv%polaramv(n)%v(k)%inv = iv%polaramv(n)%v(k)%inv - 360.0
                  if (iv%polaramv(n)%v(k)%inv < -180.0 ) iv%polaramv(n)%v(k)%inv = iv%polaramv(n)%v(k)%inv + 360.0
              end if
          else

             if (ob%polaramv(n)%u(k) > missing_r .AND. iv%polaramv(n)%u(k)%qc >= obs_qc_pointer) then
                 iv%polaramv(n)%u(k)%inv = ob%polaramv(n)%u(k) - model_u(k,n)
             end if

             if (ob%polaramv(n)%v(k) > missing_r .AND. iv%polaramv(n)%v(k)%qc >= obs_qc_pointer) then
                 iv%polaramv(n)%v(k)%inv = ob%polaramv(n)%v(k) - model_v(k,n)
             end if
          end if

      end do
   end do
   
   !------------------------------------------------------------------------
   ! Perform optional maximum error check:
   !------------------------------------------------------------------------

    if ( check_max_iv ) &
       call da_check_max_iv_polaramv(iv, it, num_qcstat_conv)     

   deallocate (model_u)
   deallocate (model_v)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_polaramv")

end subroutine da_get_innov_vector_polaramv


subroutine da_calculate_grady_polaramv(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates Gradient of Polar AMVs  Obs.          
   !-------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer :: n , k
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_polaramv")

   do n=1, iv%info(polaramv)%nlocal
      do k=1, iv%info(polaramv)%levels(n)
         if (iv%polaramv(n)%u(k)%qc < obs_qc_pointer) then
            re%polaramv(n)%u(k) = 0.0
         end if
         if (iv%polaramv(n)%v(k)%qc < obs_qc_pointer) then
            re%polaramv(n)%v(k) = 0.0
         end if

         jo_grad_y%polaramv(n)%u(k) = -re%polaramv(n)%u(k) / &
             (iv%polaramv(n)%u(k)%error * iv%polaramv(n)%u(k)%error)
         jo_grad_y%polaramv(n)%v(k) = -re%polaramv(n)%v(k) / &
            (iv%polaramv(n)%v(k)%error * iv%polaramv(n)%v(k)%error)
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_polaramv")

end subroutine da_calculate_grady_polaramv



end module da_polaramv     

