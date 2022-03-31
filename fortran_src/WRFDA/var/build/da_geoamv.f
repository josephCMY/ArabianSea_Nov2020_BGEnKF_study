












module da_geoamv

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, &
      missing, max_error_uv, max_error_t, rootproc, kms,kme,kts,kte, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv, trace_use_dull, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp,fails_error_max, &
      max_error_bt, max_error_buv, geoamv, anal_type_verify,qcstat_conv_unit,ob_vars, &
      convert_fd2uv, convert_uv2fd, max_error_spd, max_error_dir, max_omb_spd, max_omb_dir, pi, qc_rej_both,&
      wind_sd_geoamv, wind_stats_sd
   use da_grid_definitions, only : da_ffdduv, da_ffdduv_model,da_ffdduv_diagnose
   use da_physics, only : da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      maxmin_type
  use da_interpolation, only : da_interp_lin_3d, da_to_zk, &
      da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk, da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_geoamv1_type
      real          :: u                        
      real          :: v                        
   end type residual_geoamv1_type

   type maxmin_geoamv_stats_type
      type (maxmin_type)         :: u, v, t, q
   end type maxmin_geoamv_stats_type

   type stats_geoamv_type
      type (maxmin_geoamv_stats_type)  :: maximum, minimum
      type (residual_geoamv1_type)     :: average, rms_err
   end type stats_geoamv_type

contains

subroutine da_ao_stats_geoamv (stats_unit, iv, re, ob)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates (Obs - Analysis) statistics for Geo.  AMV's
   !
   !-------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O
   type(y_type),   intent (in)   :: ob            ! Observation structure.

   type (stats_geoamv_type) :: stats
   integer                  :: nu, nv
   integer                  :: n, k
   real                     :: u_inc, v_inc, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_ao_stats_geoamv")

   nu = 0
   nv = 0

   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_geoamv1_type (0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(geoamv)%nlocal
      if (iv%info(geoamv)%proc_domain(1,n)) then
         do k=1, iv%info(geoamv)%levels(n)

            u_inc = re%geoamv(n)%u(k)
            v_inc = re%geoamv(n)%v(k)
            u_obs = ob%geoamv(n)%u(k)
            v_obs = ob%geoamv(n)%v(k)

            if (.not. wind_sd_geoamv .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%geoamv(n)%u(k)%qc, iv%geoamv(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_geoamv .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%geoamv(n)%u(k)%qc, iv%geoamv(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate (n, 0, iv%geoamv(n)%u(k)%qc, & 
               u_inc, nu, & 
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, 0, iv%geoamv(n)%v(k)%qc, & 
               v_inc, nv, & 
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   iv%nstats(geoamv) = nu + nv

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for geoamv'
         call da_print_stats_geoamv (stats_unit, nu, nv, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_geoamv")

end subroutine da_ao_stats_geoamv  


subroutine da_jo_and_grady_geoamv( iv, re, jo, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose:  Calculates Cost function and Gradient for Geo. CVMs 
   !-------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv          ! Innovation vector.
   type(y_type),  intent(in)    :: re          ! Residual vector.
   type(y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type(jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_geoamv")

   jo % geoamv_u = 0.0
   jo % geoamv_v = 0.0

   if (iv%info(geoamv)%nlocal > 0) then
      do n=1, iv%info(geoamv)%nlocal
         do k=1, iv%info(geoamv)%levels(n)
            jo_grad_y%geoamv(n)%u(k) = -re%geoamv(n)%u(k) / ( iv%geoamv(n)%u(k)%error * iv%geoamv(n)%u(k)%error)
            jo_grad_y%geoamv(n)%v(k) = -re%geoamv(n)%v(k) / ( iv%geoamv(n)%v(k)%error * iv%geoamv(n)%v(k)%error)
         end do

         do k=1, iv%info(geoamv)%levels(n)
            if (iv%info(geoamv)%proc_domain(1,n)) then
               jo % geoamv_u = jo % geoamv_u - re%geoamv(n)%u(k) * jo_grad_y%geoamv(n)%u(k)
               jo % geoamv_v = jo % geoamv_v - re%geoamv(n)%v(k) * jo_grad_y%geoamv(n)%v(k)
            end if
         end do
      end do

      jo % geoamv_u = 0.5 * jo % geoamv_u
      jo % geoamv_v = 0.5 * jo % geoamv_v
   end if

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_geoamv")
     
end subroutine da_jo_and_grady_geoamv


subroutine da_residual_geoamv(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates residual vector for Geo. CMV's
   !               
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual vector (O-A).

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type)              :: n_obs_bad
   integer                           :: n, k

   if (trace_use_dull) call da_trace_entry("da_residual_geoamv")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)

   do n=1, iv%info(geoamv)%nlocal
      do k=1, iv%info(geoamv)%levels(n)
         np_available = np_available + 2
         re%geoamv(n)%u(k) = da_residual(n, 0, y%geoamv(n)%u(k), iv%geoamv(n)%u(k), n_obs_bad % u)
         re%geoamv(n)%v (k)= da_residual(n, 0, y%geoamv(n)%v(k), iv%geoamv(n)%v(k), n_obs_bad % v)
      end do
   end do

   np_missing = np_missing + n_obs_bad % u % num % miss + n_obs_bad % v % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad + n_obs_bad % v % num % bad
   np_obs_used = np_obs_used + n_obs_bad % u % num % use + n_obs_bad % v % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_geoamv")

end subroutine da_residual_geoamv


subroutine da_oi_stats_geoamv(stats_unit, iv, ob)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates (Obs - Analysis) statistics for Geo. AMV's
   !-------------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_geoamv_type) :: stats
   integer                  :: nu, nv
   integer                  :: n, k
   real                     :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_geoamv")

   nu = 0
   nv = 0

   stats%maximum%u = maxmin_type(missing_r, 0, 0)
   stats%maximum%v = maxmin_type(missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_geoamv1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(geoamv)%nlocal
      if (iv%info(geoamv)%proc_domain(1,n)) then
         do k=1, iv%info(geoamv)%levels(n)

            u_inv = iv%geoamv(n)%u(k)%inv
            v_inv = iv%geoamv(n)%v(k)%inv
            u_obs = ob%geoamv(n)%u(k)
            v_obs = ob%geoamv(n)%v(k)

            if (.not. wind_sd_geoamv .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%geoamv(n)%u(k)%qc, iv%geoamv(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_geoamv .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%geoamv(n)%u(k)%qc, iv%geoamv(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate(iv%info(geoamv)%obs_global_index(n), &
               0, iv%geoamv(n)%u(k)%qc, &
               u_inv, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(geoamv)%obs_global_index(n), &
               0, iv%geoamv(n)%v(k)%qc, &
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for geoamv'
         call da_print_stats_geoamv(stats_unit, nu, nv, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_geoamv")

end subroutine da_oi_stats_geoamv


subroutine da_print_stats_geoamv(stats_unit, nu, nv, geoamv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                  intent(in)    :: stats_unit
   integer,                  intent(inout) :: nu, nv
   type (stats_geoamv_type), intent(in)    :: geoamv

   if (trace_use_dull) call da_trace_entry("da_print_stats_geoamv")

   write(unit=stats_unit, fmt='(a/)') &
      '   var             u (m/s)     n    k    v (m/s)     n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
      '  Number: ', nu, nv

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1

   write(unit=stats_unit, fmt='((a,2(f12.4,2i5)))') &
      ' Minimum(n,k): ', geoamv%minimum%u, geoamv%minimum%v, &
      ' Maximum(n,k): ', geoamv%maximum%u, geoamv%maximum%v

   write(unit=stats_unit, fmt='((a,2(f12.4,10x)))') &
      ' Average     : ', geoamv%average%u/real(nu), geoamv%average%v/real(nv), &
      '    RMSE     : ', sqrt(geoamv%rms_err%u/real(nu)), &
                         sqrt(geoamv%rms_err%v/real(nv))

   if (trace_use_dull) call da_trace_exit("da_print_stats_geoamv")

end subroutine da_print_stats_geoamv


subroutine da_transform_xtoy_geoamv (grid, iv, y)

   !-------------------------------------------------------------------------
   ! Purpose: X to Y Transform operator for Geo. AMVs
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer           :: n,k
   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)

   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_geoamv")

   allocate (u(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))
   allocate (v(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))

   allocate (ub(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))
   allocate (vb(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))

   call da_interp_lin_3d (grid%xa%u, iv%info(geoamv), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(geoamv), v)

   call da_interp_lin_3d (grid%xb%u, iv%info(geoamv), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(geoamv), vb)

   do n=iv%info(geoamv)%n1,iv%info(geoamv)%n2
      do k = 1, iv%info(geoamv)%levels(n)
         if(wind_sd_geoamv) then
            call da_uv_to_sd_lin(y%geoamv(n)%u(k),y%geoamv(n)%v(k),u(k,n),v(k,n),ub(k,n),vb(k,n))
         else
            y%geoamv(n)%u(k) = u(k,n)
            y%geoamv(n)%v(k) = v(k,n)
         end if
       end do
   end do

   deallocate (u)
   deallocate (v)
   deallocate (vb)
   deallocate (ub)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_geoamv")

end subroutine da_transform_xtoy_geoamv


subroutine da_transform_xtoy_geoamv_adj (grid, iv, jo_grad_y, jo_grad_x)

   !-------------------------------------------------------------------------
   ! Purpose: X to Y Transpose operator for Geo. AMVs
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

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_geoamv_adj")
 
   allocate (u(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))
   allocate (v(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))

   allocate (ub(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))
   allocate (vb(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))

   call da_interp_lin_3d (grid%xb%u, iv%info(geoamv), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(geoamv), vb)

   do n=iv%info(geoamv)%n1,iv%info(geoamv)%n2
      do k = 1, iv%info(geoamv)%levels(n)
          if(wind_sd_geoamv) then
             call da_uv_to_sd_adj(jo_grad_y%geoamv(n)%u(k), &
                                  jo_grad_y%geoamv(n)%v(k), u(k,n), v(k,n), ub(k,n), vb(k,n))
          else
             u(k,n) = jo_grad_y%geoamv(n)%u(k)
             v(k,n) = jo_grad_y%geoamv(n)%v(k)
          end if
      end do
   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(geoamv), u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(geoamv), v)

   deallocate (u)
   deallocate (v)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_geoamv_adj")

end subroutine da_transform_xtoy_geoamv_adj


subroutine da_check_max_iv_geoamv(iv, it, num_qcstat_conv)    

   !-------------------------------------------------------------------------
   ! Purpose: Innovation vector check for Geo. AMVs               
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-------------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: k,n,ipr
   logical :: failed,failed1,failed2

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_geoamv")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(geoamv)%n1,iv%info(geoamv)%n2
      do k = 1, iv%info(geoamv)%levels(n)
         call da_get_print_lvl(iv%geoamv(n)%p(k),ipr)

         if(.not. qc_rej_both)then
            if(wind_sd_geoamv)then
               failed=.false.
               if( iv%geoamv(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%u(k), max_error_spd,failed)
                   if( iv%info(geoamv)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,geoamv,1,ipr) = num_qcstat_conv(1,geoamv,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,geoamv,1,ipr) = num_qcstat_conv(2,geoamv,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'geoamv',ob_vars(1),iv%info(geoamv)%lat(k,n),iv%info(geoamv)%lon(k,n),0.01*iv%geoamv(n)%p(k)
                       end if
                   end if
                end if

                failed=.false.
                if( iv%geoamv(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%v(k), max_error_dir,failed)
                    if( iv%info(geoamv)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,geoamv,2,ipr) = num_qcstat_conv(1,geoamv,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,geoamv,2,ipr) = num_qcstat_conv(2,geoamv,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'geoamv',ob_vars(2),iv%info(geoamv)%lat(k,n),iv%info(geoamv)%lon(k,n),0.01*iv%geoamv(n)%p(k)
                        end if
                    end if
                end if

            else

               failed=.false.
               if( iv%geoamv(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%u(k), max_error_uv,failed)
                   if( iv%info(geoamv)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,geoamv,1,ipr) = num_qcstat_conv(1,geoamv,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,geoamv,1,ipr) = num_qcstat_conv(2,geoamv,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'geoamv',ob_vars(1),iv%info(geoamv)%lat(k,n),iv%info(geoamv)%lon(k,n),0.01*iv%geoamv(n)%p(k)
                       end if
                   end if
               end if

               failed=.false.
               if( iv%geoamv(n)%v(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%v(k), max_error_uv,failed)
                   if( iv%info(geoamv)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,geoamv,2,ipr) = num_qcstat_conv(1,geoamv,2,ipr) + 1
                       if(failed)then
                          num_qcstat_conv(2,geoamv,2,ipr) = num_qcstat_conv(2,geoamv,2,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'geoamv',ob_vars(2),iv%info(geoamv)%lat(k,n),iv%info(geoamv)%lon(k,n),0.01*iv%geoamv(n)%p(k)
                       end if
                   end if
               end if

               if(wind_sd_geoamv)then
                  if(iv%geoamv(n)%u(k)%qc == fails_error_max .or.  abs(iv%geoamv(n)%u(k)%inv) >= max_omb_spd) then
                     iv%geoamv(n)%u(k)%qc = fails_error_max
                     iv%geoamv(n)%u(k)%inv = 0.0
                  endif
                  if(iv%geoamv(n)%v(k)%qc == fails_error_max .or.  abs(iv%geoamv(n)%v(k)%inv) >= max_omb_dir) then
                     iv%geoamv(n)%v(k)%qc = fails_error_max
                     iv%geoamv(n)%v(k)%inv = 0.0
                  endif
               endif
            endif
         else
            failed1=.false.
            failed2=.false.

            if( iv%geoamv(n)%v(k)%qc >= obs_qc_pointer .or. iv%geoamv(n)%u(k)%qc >= obs_qc_pointer )  then
                if(wind_sd_geoamv)then
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%u(k), max_error_spd,failed1)
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%v(k), max_error_dir,failed2)
                else
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%u(k), max_error_uv,failed1)
                   call da_max_error_qc (it,iv%info(geoamv), n, iv%geoamv(n)%v(k), max_error_uv,failed2)
                endif
            endif

            if( iv%info(geoamv)%proc_domain(k,n) ) then
                num_qcstat_conv(1,geoamv,1,ipr) = num_qcstat_conv(1,geoamv,1,ipr) + 1
                num_qcstat_conv(1,geoamv,2,ipr) = num_qcstat_conv(1,geoamv,2,ipr) + 1

                if(failed1 .or. failed2) then
                   num_qcstat_conv(2,geoamv,1,ipr) = num_qcstat_conv(2,geoamv,1,ipr) + 1
                   write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'geoamv',ob_vars(1),iv%info(geoamv)%lat(k,n),iv%info(geoamv)%lon(k,n),0.01*iv%geoamv(n)%p(k)
                   num_qcstat_conv(2,geoamv,2,ipr) = num_qcstat_conv(2,geoamv,2,ipr) + 1
                   write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                         'geoamv',ob_vars(2),iv%info(geoamv)%lat(k,n),iv%info(geoamv)%lon(k,n),0.01*iv%geoamv(n)%p(k)
                endif
            endif

            if(wind_sd_geoamv)then
               if(iv%geoamv(n)%u(k)%qc == fails_error_max .or. iv%geoamv(n)%v(k)%qc == fails_error_max .or. &
                  abs(iv%geoamv(n)%v(k)%inv) >= max_omb_dir .or. abs(iv%geoamv(n)%u(k)%inv) >= max_omb_spd )then
                  iv%geoamv(n)%u(k)%qc = fails_error_max
                  iv%geoamv(n)%v(k)%qc = fails_error_max
                  iv%geoamv(n)%u(k)%inv = 0.0
                  iv%geoamv(n)%v(k)%inv = 0.0
               endif
            else
               if(iv%geoamv(n)%u(k)%qc == fails_error_max .or. iv%geoamv(n)%v(k)%qc == fails_error_max ) then
                  iv%geoamv(n)%u(k)%qc = fails_error_max
                  iv%geoamv(n)%v(k)%qc = fails_error_max
                  iv%geoamv(n)%u(k)%inv = 0.0
                  iv%geoamv(n)%v(k)%inv = 0.0
               endif
            endif
         endif

      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_geoamv")

end subroutine da_check_max_iv_geoamv      

subroutine da_get_innov_vector_geoamv( it,num_qcstat_conv,  grid, ob, iv)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates innovation vector does QC for geoamv
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-------------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid    ! first guess state.
   type(y_type),     intent(in)    :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n              ! Loop counter.
   integer :: k              ! Index dimension.
   integer :: num_levs       ! Number of obs levels.
   real    :: speed, direction

   integer :: i  (kms:kme)
   integer :: j  (kms:kme)
   real    :: dx (kms:kme)
   real    :: dxm(kms:kme)  
   real    :: dy (kms:kme)
   real    :: dym(kms:kme) 
   real,allocatable :: model_u(:,:)
   real,allocatable :: model_v(:,:)

   real    :: v_p(kts:kte)      ! Model value p at ob hor. location.

   integer :: itu,ituf,itvv,itvvf

   if (trace_use_dull) call da_trace_entry ("da_get_innov_vector_geoamv")

   allocate (model_u(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))
   allocate (model_v(iv%info(geoamv)%max_lev,iv%info(geoamv)%n1:iv%info(geoamv)%n2))

   model_u(:,:) = 0.0
   model_v(:,:) = 0.0
   
   if ( it > 1 ) then
      do n = iv%info(geoamv)%n1, iv%info(geoamv)%n2
         do k = 1, iv%info(geoamv)%levels(n)
            if (iv%geoamv(n)%u(k)%qc == fails_error_max) iv%geoamv(n)%u(k)%qc = 0
            if (iv%geoamv(n)%v(k)%qc == fails_error_max) iv%geoamv(n)%v(k)%qc = 0
         end do
      end do
   end if

   do n = iv%info(geoamv)%n1, iv%info(geoamv)%n2
      ! [1.3] Get horizontal interpolation weights:

      num_levs = iv%info(geoamv)%levels(n)
      if (num_levs < 1) cycle

      ! slower
      ! i(:)   = iv%info(geoamv)%i(:,n)
      ! j(:)   = iv%info(geoamv)%j(:,n)
      ! dx(:)  = iv%info(geoamv)%dx(:,n)
      ! dy(:)  = iv%info(geoamv)%dy(:,n)
      ! dxm(:) = iv%info(geoamv)%dxm(:,n)
      ! dym(:) = iv%info(geoamv)%dym(:,n)

      ! faster
      i(1)   = iv%info(geoamv)%i(1,n)
      j(1)   = iv%info(geoamv)%j(1,n)
      dx(1)  = iv%info(geoamv)%dx(1,n)
      dy(1)  = iv%info(geoamv)%dy(1,n)
      dxm(1) = iv%info(geoamv)%dxm(1,n)
      dym(1) = iv%info(geoamv)%dym(1,n)

      ! if position varies with height, slower
      ! do k=kts,kte
      !    v_p(k) = dym(k)*(dxm(k)*grid%xb%p(i(k),j(k),k)+dx(k)*grid%xb%p(i(k)+1,j(k),k)) &
      !       + dy(k)*(dxm(k)*grid%xb%p(i(k),j(k)+1,k)+dx(k)*grid%xb%p(i(k)+1,j(k)+1,k))
      ! end do
 
      ! If position does not, faster
      v_p(kts:kte) = dym(1)*(dxm(1)*grid%xb%p(i(1),j(1),kts:kte) + dx(1)*grid%xb%p(i(1)+1,j(1),kts:kte)) &
         + dy(1)*(dxm(1)*grid%xb%p(i(1),j(1)+1,kts:kte) + dx(1)*grid%xb%p(i(1)+1,j(1)+1,kts:kte))

      do k=1, iv%info(geoamv)%levels(n)
         if (iv%geoamv(n)%p(k) > 1.0) then
            call da_to_zk (iv%geoamv(n)%p(k), v_p, v_interp_p, iv%info(geoamv)%zk(k,n))
         end if
      end do

   end do

   call da_convert_zk (iv%info(geoamv))

   if (.not. anal_type_verify) then
      do n = iv%info(geoamv)%n1, iv%info(geoamv)%n2
         do k=1, iv%info(geoamv)%levels(n)
            if (iv%info(geoamv)%zk(k,n) < 0.0) then
               iv%geoamv(n)%u(k)% qc = missing_data
               iv%geoamv(n)%v(k)% qc = missing_data
            end if
         end do
      end do
   end if

   call da_interp_lin_3d (grid%xb%u, iv%info(geoamv), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(geoamv), model_v)

   do n = iv%info(geoamv)%n1, iv%info(geoamv)%n2
      do k = 1, iv%info(geoamv)%levels(n)
         iv%geoamv(n)%u(k)%inv = 0.0
         iv%geoamv(n)%v(k)%inv = 0.0

         if (wind_sd_geoamv) then
             call da_ffdduv_model (speed,direction,model_u(k,n), model_v(k,n), convert_uv2fd)

             if (ob%geoamv(n)%u(k) > missing_r .AND. iv%geoamv(n)%u(k)%qc >= obs_qc_pointer) then
                 iv%geoamv(n)%u(k)%inv = ob%geoamv(n)%u(k) - speed
             end if

             if (ob%geoamv(n)%v(k) > missing_r .AND. iv%geoamv(n)%v(k)%qc >= obs_qc_pointer) then
                 iv%geoamv(n)%v(k)%inv = ob%geoamv(n)%v(k) - direction
                 if (iv%geoamv(n)%v(k)%inv > 180.0 ) iv%geoamv(n)%v(k)%inv = iv%geoamv(n)%v(k)%inv - 360.0
                 if (iv%geoamv(n)%v(k)%inv < -180.0 ) iv%geoamv(n)%v(k)%inv = iv%geoamv(n)%v(k)%inv + 360.0
             end if
         else
             if (ob%geoamv(n)%u(k) > missing_r .AND. iv%geoamv(n)%u(k)%qc >= obs_qc_pointer) then
                 iv%geoamv(n)%u(k)%inv = ob%geoamv(n)%u(k) - model_u(k,n)
             end if

             if (ob%geoamv(n)%v(k) > missing_r .AND. iv%geoamv(n)%v(k)%qc >= obs_qc_pointer) then
                 iv%geoamv(n)%v(k)%inv = ob%geoamv(n)%v(k) - model_v(k,n)
             end if
         end if

      end do
   end do

   !------------------------------------------------------------------------
   ! Perform optional maximum error check:
   !------------------------------------------------------------------------
   if ( check_max_iv ) &
      call da_check_max_iv_geoamv(iv,it,num_qcstat_conv)     

   deallocate (model_u)
   deallocate (model_v)

   if (trace_use_dull) call da_trace_exit ("da_get_innov_vector_geoamv")
   
end subroutine da_get_innov_vector_geoamv
subroutine da_calculate_grady_geoamv(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Calculates Gradient of Geo. CMVs Obs.
   !               
   !-------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n , k
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_geoamv")  

   do n=1, iv%info(geoamv)%nlocal
      do k=1, iv%info(geoamv)%levels(n)
         if (iv%geoamv(n)%u(k)%qc < obs_qc_pointer) re%geoamv(n)%u(k) = 0.0
         if (iv%geoamv(n)%v(k)%qc < obs_qc_pointer) re%geoamv(n)%v(k) = 0.0

         jo_grad_y%geoamv(n)%u(k) = -re%geoamv(n)%u(k) / (iv%geoamv(n)%u(k)%error * iv%geoamv(n)%u(k)%error)
         jo_grad_y%geoamv(n)%v(k) = -re%geoamv(n)%v(k) / (iv%geoamv(n)%v(k)%error * iv%geoamv(n)%v(k)%error)
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_geoamv")  

end subroutine da_calculate_grady_geoamv



end module da_geoamv     

