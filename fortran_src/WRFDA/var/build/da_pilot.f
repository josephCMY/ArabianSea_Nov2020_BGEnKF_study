












module da_pilot

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing_data, max_error_uv, max_error_t, rootproc, &
      pilot, max_error_p,max_error_q, trace_use_dull,fails_error_max, &
      max_stheight_diff, anal_type_verify, kms,kme,kts,kte,ob_vars,qcstat_conv_unit, &
      convert_fd2uv, convert_uv2fd, max_error_spd, max_error_dir, &
      max_omb_spd, max_omb_dir, pi, qc_rej_both, & 
      wind_sd_pilot, wind_stats_sd 
   use da_grid_definitions, only : da_ffdduv, da_ffdduv_model, da_ffdduv_diagnose
   use da_physics, only : da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_interp_lin_3d, da_to_zk, &
      da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_pilot1_type
      real          :: u                        
      real          :: v                        
   end type residual_pilot1_type

   type maxmin_pilot_stats_type
      type (maxmin_type)         :: u, v
   end type maxmin_pilot_stats_type

   type stats_pilot_type
      type (maxmin_pilot_stats_type)  :: maximum, minimum
      type (residual_pilot1_type)     :: average, rms_err
   end type stats_pilot_type

contains

subroutine da_ao_stats_pilot (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type  (y_type), intent(in)    :: re            ! A - O
   type(y_type),   intent (in)   :: ob            ! Observation structure.

   type (stats_pilot_type) :: stats
   integer                 :: nu, nv
   integer                 :: n, k
   real                    :: u_inc, v_inc, u_obs, v_obs
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_pilot")

   nu = 0
   nv = 0

   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_pilot1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(pilot)%nlocal
      if (iv%info(pilot)%proc_domain(1,n)) then
         do k=1, iv%info(pilot)%levels(n)

            u_inc = re%pilot(n)%u(k)
            v_inc = re%pilot(n)%v(k)
            u_obs = ob%pilot(n)%u(k)
            v_obs = ob%pilot(n)%v(k)

            if (.not. wind_sd_pilot .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%pilot(n)%u(k)%qc, iv%pilot(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_pilot .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                       iv%pilot(n)%u(k)%qc, iv%pilot(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate (n, k, iv%pilot(n)%u(k)%qc, & 
               u_inc, nu, & 
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%pilot(n)%v(k)%qc, & 
               v_inc, nv, & 
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
         end do
      end if    ! end if (iv%info(pilot)%proc_domain(1,n))
   end do

   ! do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   iv%nstats(pilot) = nu + nv

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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for pilot'
         call da_print_stats_pilot(stats_unit, nu, nv, stats)
      end if
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_pilot")

end subroutine da_ao_stats_pilot


subroutine da_jo_and_grady_pilot(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type(jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_pilot")

   jo % pilot_u = 0.0
   jo % pilot_v = 0.0

   do n=1, iv%info(pilot)%nlocal
      do k=1, iv%info(pilot)%levels(n)
         jo_grad_y%pilot(n)%u(k) = -re%pilot(n)%u(k) / (iv%pilot(n)%u(k)%error * iv%pilot(n)%u(k)%error)
         jo_grad_y%pilot(n)%v(k) = -re%pilot(n)%v(k) / (iv%pilot(n)%v(k)%error * iv%pilot(n)%v(k)%error)
      end do

     if (iv%info(pilot)%proc_domain(1,n)) then
        do k=1, iv%info(pilot)%levels(n)
           jo % pilot_u = jo % pilot_u - re%pilot(n)%u(k) * jo_grad_y%pilot(n)%u(k)
           jo % pilot_v = jo % pilot_v - re%pilot(n)%v(k) * jo_grad_y%pilot(n)%v(k)
        end do
      end if
   end do

   jo % pilot_u = 0.5 * jo % pilot_u
   jo % pilot_v = 0.5 * jo % pilot_v

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_pilot")

end subroutine da_jo_and_grady_pilot


subroutine da_residual_pilot(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for pilots
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

   if (trace_use_dull) call da_trace_entry("da_residual_pilot")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)

   do n=1, iv%info(pilot)%nlocal
      do k=1, iv%info(pilot)%levels(n)
         np_available = np_available + 2
         re%pilot(n)%u(k) = da_residual(n, k, y%pilot(n)%u(k), iv%pilot(n)%u(k), n_obs_bad % u)
         re%pilot(n)%v(k) = da_residual(n, k, y%pilot(n)%v(k), iv%pilot(n)%v(k), n_obs_bad % v)
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % u % num % miss + n_obs_bad % v % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad + n_obs_bad % v % num % bad 
   np_obs_used = np_obs_used + n_obs_bad % u % num % use + n_obs_bad % v % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_pilot")

end subroutine da_residual_pilot


subroutine da_oi_stats_pilot (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_pilot_type) :: stats
   integer                 :: nu, nv
   integer                 :: n, k
   real                    :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_pilot")

   nu = 0
   nv = 0

   stats%maximum%u = maxmin_type(missing_r, 0, 0)
   stats%maximum%v = maxmin_type(missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_pilot1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(pilot)%nlocal
      if (iv%info(pilot)%proc_domain(1,n)) then
         do k=1, iv%info(pilot)%levels(n)

            u_inv = iv%pilot(n)%u(k)%inv
            v_inv = iv%pilot(n)%v(k)%inv
            u_obs = ob%pilot(n)%u(k)
            v_obs = ob%pilot(n)%v(k)

            if (.not. wind_sd_pilot .and. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%pilot(n)%u(k)%qc, iv%pilot(n)%v(k)%qc, convert_uv2fd)
            if (wind_sd_pilot .and. .not. wind_stats_sd) &
               call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                       iv%pilot(n)%u(k)%qc, iv%pilot(n)%v(k)%qc, convert_fd2uv)

            call da_stats_calculate(iv%info(pilot)%obs_global_index(n), &
               k, iv%pilot(n)%u(k)%qc,  &
               u_inv, nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate(iv%info(pilot)%obs_global_index(n), &
               k, iv%pilot(n)%v(k)%qc, &
               v_inv, nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
         end do
      end if
   end do

   ! do inter-processor communication to gather statistics.
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for pilot'
         call da_print_stats_pilot(stats_unit, nu, nv, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_pilot")

 end subroutine da_oi_stats_pilot


subroutine da_print_stats_pilot(stats_unit, nu, nv, Pilot)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nu, nv
   type (stats_pilot_type), intent(in)    :: pilot

   if (trace_use_dull) call da_trace_entry("da_print_stats_pilot")

   write(unit=stats_unit, fmt='(a/)') &
      '   var             u (m/s)     n    k    v (m/s)     n    k  '

   write(unit=stats_unit, fmt='(a,i16,i22)') &
      '  Number: ', nu, nv

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1

   write(unit=stats_unit, fmt='((a,2(f12.4,2i5)))') &
       ' Minimum(n,k): ', pilot%minimum%u, pilot%minimum%v, &
       ' Maximum(n,k): ', pilot%maximum%u, pilot%maximum%v
   write(unit=stats_unit, fmt='((a,2(f12.4,10x)))')  &
      ' Average     : ', pilot%average%u/real(nu), &
                         pilot%average%v/real(nv),      &
      '    RMSE     : ', sqrt(pilot%rms_err%u/real(nu)), &
                         sqrt(pilot%rms_err%v/real(nv))

   if (trace_use_dull) call da_trace_exit("da_print_stats_pilot")

end subroutine da_print_stats_pilot


subroutine da_transform_xtoy_pilot (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n, k  ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_pilot")

   allocate (model_u(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (model_v(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (ub(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (vb(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))

   call da_interp_lin_3d (grid%xa%u, iv%info(pilot), model_u)
   call da_interp_lin_3d (grid%xa%v, iv%info(pilot), model_v)
   call da_interp_lin_3d (grid%xb%u, iv%info(pilot), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(pilot), vb)

   do n=iv%info(pilot)%n1,iv%info(pilot)%n2
      do k = 1, iv%info(pilot)%levels(n)

         if(wind_sd_pilot) then
            call da_uv_to_sd_lin(y%pilot(n)%u(k),y%pilot(n)%v(k),model_u(k,n),model_v(k,n),ub(k,n),vb(k,n))
         else
            y%pilot(n)%u(k) = model_u(k,n)
            y%pilot(n)%v(k) = model_v(k,n)
         end if

      end do

   end do

   deallocate (model_u)
   deallocate (model_v)
   deallocate (ub)
   deallocate (vb)
   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_pilot")

end subroutine da_transform_xtoy_pilot 


subroutine da_transform_xtoy_pilot_adj(grid, iv, jo_grad_y, jo_grad_x)

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

   integer :: n, k  ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_pilot_adj")

   allocate (model_u(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (model_v(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (ub(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (vb(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   call da_interp_lin_3d (grid%xb%u, iv%info(pilot), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(pilot), vb)

   do n=iv%info(pilot)%n1,iv%info(pilot)%n2
       do k = 1, iv%info(pilot)%levels(n)

         if(wind_sd_pilot) then
            call da_uv_to_sd_adj(jo_grad_y%pilot(n)%u(k), &
                                 jo_grad_y%pilot(n)%v(k), model_u(k,n), model_v(k,n), ub(k,n), vb(k,n))
         else
             model_u(k,n) = jo_grad_y%pilot(n)%u(k)
             model_v(k,n) = jo_grad_y%pilot(n)%v(k)
         end if

      end do

   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(pilot), model_u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(pilot), model_v)
   call da_interp_lin_3d (grid%xb%u, iv%info(pilot), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(pilot), vb)

   deallocate(model_u)
   deallocate(model_v)
   deallocate (ub)
   deallocate (vb)
   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_pilot_adj")

end subroutine da_transform_xtoy_pilot_adj


subroutine da_check_max_iv_pilot(iv, it, num_qcstat_conv)     

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
   
   if (trace_use_dull) call da_trace_entry("da_check_max_iv_pilot")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(pilot)%n1,iv%info(pilot)%n2
      do k = 1, iv%info(pilot)%levels(n)
         call da_get_print_lvl(iv%pilot(n)%p(k),ipr)

         if(.not. qc_rej_both)then
            if(wind_sd_pilot)then
               failed=.false.
               if( iv%pilot(n)%u(k)%qc >= obs_qc_pointer ) then
                   call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%u(k), max_error_spd,failed)
                   if( iv%info(pilot)%proc_domain(k,n) ) then
                       num_qcstat_conv(1,pilot,1,ipr) = num_qcstat_conv(1,pilot,1,ipr) + 1
                       if(failed) then
                          num_qcstat_conv(2,pilot,1,ipr) = num_qcstat_conv(2,pilot,1,ipr) + 1
                          write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                          'pilot',ob_vars(1),iv%info(pilot)%lat(k,n),iv%info(pilot)%lon(k,n),0.01*iv%pilot(n)%p(k)
                       end if
                   end if
                end if

                failed=.false.
                if( iv%pilot(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%v(k), max_error_dir,failed)
                    if( iv%info(pilot)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,pilot,2,ipr) = num_qcstat_conv(1,pilot,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,pilot,2,ipr) = num_qcstat_conv(2,pilot,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'pilot',ob_vars(2),iv%info(pilot)%lat(k,n),iv%info(pilot)%lon(k,n),0.01*iv%pilot(n)%p(k)
                        end if
                    end if
                end if
             else
                failed=.false.
                if( iv%pilot(n)%u(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%u(k), max_error_uv,failed)
                    if( iv%info(pilot)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,pilot,1,ipr) = num_qcstat_conv(1,pilot,1,ipr) + 1
                        if(failed) then
                           num_qcstat_conv(2,pilot,1,ipr) = num_qcstat_conv(2,pilot,1,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'pilot',ob_vars(1),iv%info(pilot)%lat(k,n),iv%info(pilot)%lon(k,n),0.01*iv%pilot(n)%p(k)
                        end if
                    end if
                end if

                failed=.false.
                if( iv%pilot(n)%v(k)%qc >= obs_qc_pointer ) then
                    call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%v(k), max_error_uv,failed)
                    if( iv%info(pilot)%proc_domain(k,n) ) then
                        num_qcstat_conv(1,pilot,2,ipr) = num_qcstat_conv(1,pilot,2,ipr) + 1
                        if(failed)then
                           num_qcstat_conv(2,pilot,2,ipr) = num_qcstat_conv(2,pilot,2,ipr) + 1
                           write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                           'pilot',ob_vars(2),iv%info(pilot)%lat(k,n),iv%info(pilot)%lon(k,n),0.01*iv%pilot(n)%p(k)
                        end if
                    end if
                end if
              end if

              if(wind_sd_pilot)then
                 if(iv%pilot(n)%u(k)%qc == fails_error_max .or.  abs(iv%pilot(n)%u(k)%inv) >= max_omb_spd) then
                 iv%pilot(n)%u(k)%qc = fails_error_max
                 iv%pilot(n)%u(k)%inv = 0.0
              endif
              if(iv%pilot(n)%v(k)%qc == fails_error_max .or.  abs(iv%pilot(n)%v(k)%inv) >= max_omb_dir) then
                 iv%pilot(n)%v(k)%qc = fails_error_max
                 iv%pilot(n)%v(k)%inv = 0.0
              endif
           endif

        else
           failed1=.false.
           failed2=.false.

           if( iv%pilot(n)%v(k)%qc >= obs_qc_pointer .or. iv%pilot(n)%u(k)%qc >= obs_qc_pointer )  then
               if(wind_sd_pilot)then
                  call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%u(k), max_error_spd,failed1)
                  call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%v(k), max_error_dir,failed2)
               else
                  call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%u(k), max_error_uv,failed1)
                  call da_max_error_qc (it,iv%info(pilot), n, iv%pilot(n)%v(k), max_error_uv,failed2)
               endif
           endif

           if( iv%info(pilot)%proc_domain(k,n) ) then
               num_qcstat_conv(1,pilot,1,ipr) = num_qcstat_conv(1,pilot,1,ipr) + 1
               num_qcstat_conv(1,pilot,2,ipr) = num_qcstat_conv(1,pilot,2,ipr) + 1

               if(failed1 .or. failed2) then
                  num_qcstat_conv(2,pilot,1,ipr) = num_qcstat_conv(2,pilot,1,ipr) + 1
                  write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                        'pilot',ob_vars(1),iv%info(pilot)%lat(k,n),iv%info(pilot)%lon(k,n),0.01*iv%pilot(n)%p(k)
                  num_qcstat_conv(2,pilot,2,ipr) = num_qcstat_conv(2,pilot,2,ipr) + 1
                  write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                        'pilot',ob_vars(2),iv%info(pilot)%lat(k,n),iv%info(pilot)%lon(k,n),0.01*iv%pilot(n)%p(k)
               endif
           endif

           if(wind_sd_pilot)then
              if(iv%pilot(n)%u(k)%qc == fails_error_max .or. iv%pilot(n)%v(k)%qc == fails_error_max .or. &
                 abs(iv%pilot(n)%v(k)%inv) >= max_omb_dir .or. abs(iv%pilot(n)%u(k)%inv) >= max_omb_spd )then
                 iv%pilot(n)%u(k)%qc = fails_error_max
                 iv%pilot(n)%v(k)%qc = fails_error_max
                 iv%pilot(n)%u(k)%inv = 0.0
                 iv%pilot(n)%v(k)%inv = 0.0
              endif
           else
              if(iv%pilot(n)%u(k)%qc == fails_error_max .or. iv%pilot(n)%v(k)%qc == fails_error_max ) then
                 iv%pilot(n)%u(k)%qc = fails_error_max
                 iv%pilot(n)%v(k)%qc = fails_error_max
                 iv%pilot(n)%u(k)%inv = 0.0
                 iv%pilot(n)%v(k)%inv = 0.0
              endif
           endif
       endif

      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_check_max_iv_pilot")

end subroutine da_check_max_iv_pilot
subroutine da_get_innov_vector_pilot( it,num_qcstat_conv, grid, ob, iv)

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

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.

   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.

   real, allocatable :: model_u(:,:)  ! Model value u at ob location.
   real, allocatable :: model_v(:,:)  ! Model value v at ob location.

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.
   real    :: speed, direction

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_pilot")

   allocate (model_u(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))
   allocate (model_v(iv%info(pilot)%max_lev,iv%info(pilot)%n1:iv%info(pilot)%n2))

   model_u(:,:) = 0.0
   model_v(:,:) = 0.0

   if ( it > 1 ) then
      do n=iv%info(pilot)%n1,iv%info(pilot)%n2
         do k=1, iv%info(pilot)%levels(n)
            if (iv%pilot(n)%u(k)%qc == fails_error_max) iv%pilot(n)%u(k)%qc = 0
            if (iv%pilot(n)%v(k)%qc == fails_error_max) iv%pilot(n)%v(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(pilot)%n1,iv%info(pilot)%n2
      ! [1.3] Get horizontal interpolation weights:

      i   = iv%info(pilot)%i(1,n)
      j   = iv%info(pilot)%j(1,n)
      dx  = iv%info(pilot)%dx(1,n)
      dy  = iv%info(pilot)%dy(1,n)
      dxm = iv%info(pilot)%dxm(1,n)
      dym = iv%info(pilot)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j,k) + dx*grid%xb%h(i+1,j,k)) + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
         v_p(k) = dym*(dxm*grid%xb%p(i,j,k) + dx*grid%xb%p(i+1,j,k)) + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
      end do

      do k=1, iv%info(pilot)%levels(n)
         if (iv % pilot(n) % p(k) > 1.0) then
            call da_to_zk(iv % pilot(n) % p(k), v_p, v_interp_p, iv%info(pilot)%zk(k,n))
         else if (iv%pilot(n) % h(k) > missing_r) then
            call da_to_zk(iv % pilot(n) % h(k), v_h, v_interp_h, iv%info(pilot)%zk(k,n))
         end if

         if (iv%info(pilot)%zk(k,n) < 0.0 .and.  .not.anal_type_verify) then
            iv % pilot(n) % u(k) % qc = missing_data
            iv % pilot(n) % v(k) % qc = missing_data
         end if
      end do
   end do

   call da_convert_zk (iv%info(pilot))

   ! [1.4] Interpolate horizontally:
   call da_interp_lin_3d (grid%xb%u, iv%info(pilot), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(pilot), model_v)

   do n=iv%info(pilot)%n1,iv%info(pilot)%n2
      !------------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !------------------------------------------------------------------------

      do k = 1, iv%info(pilot)%levels(n)
         iv % pilot(n) % u(k) % inv = 0.0
         iv % pilot(n) % v(k) % inv = 0.0

         !------------------------------------------------------------------------
         ! [4.0] Fast interpolation:
         !------------------------------------------------------------------------


          if (wind_sd_pilot) then
              call da_ffdduv_model (speed,direction,model_u(k,n), model_v(k,n), convert_uv2fd)

              if (ob%pilot(n)%u(k) > missing_r .AND. iv%pilot(n)%u(k)%qc >= obs_qc_pointer) then
                  iv%pilot(n)%u(k)%inv = ob%pilot(n)%u(k) - speed
              end if

              if (ob%pilot(n)%v(k) > missing_r .AND. iv%pilot(n)%v(k)%qc >= obs_qc_pointer) then
                  iv%pilot(n)%v(k)%inv = ob%pilot(n)%v(k) - direction
                  if (iv%pilot(n)%v(k)%inv > 180.0 ) iv%pilot(n)%v(k)%inv = iv%pilot(n)%v(k)%inv - 360.0
                  if (iv%pilot(n)%v(k)%inv < -180.0 ) iv%pilot(n)%v(k)%inv = iv%pilot(n)%v(k)%inv + 360.0
              end if
          else
             if (ob%pilot(n)%u(k) > missing_r .AND. iv%pilot(n)%u(k)%qc >= obs_qc_pointer) then
                 iv%pilot(n)%u(k)%inv = ob%pilot(n)%u(k) - model_u(k,n)
             end if

             if (ob%pilot(n)%v(k) > missing_r .AND. iv%pilot(n)%v(k)%qc >= obs_qc_pointer) then
                 iv%pilot(n)%v(k)%inv = ob%pilot(n)%v(k) - model_v(k,n)
             end if
          end if

      end do
   end do

   !------------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_pilot(iv, it, num_qcstat_conv)    

   deallocate (model_u)
   deallocate (model_v)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_pilot")

end subroutine da_get_innov_vector_pilot


subroutine da_calculate_grady_pilot(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer  :: n, k
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_pilot")

   do n=1, iv%info(pilot)%nlocal
      do k=1, iv%info(pilot)%levels(n)
         if (iv%pilot(n)%u(k)%qc < obs_qc_pointer) re%pilot(n)%u(k) = 0.0
         if (iv%pilot(n)%v(k)%qc < obs_qc_pointer) re%pilot(n)%v(k) = 0.0

         jo_grad_y%pilot(n)%u(k) = -re%pilot(n)%u(k) / (iv%pilot(n)%u(k)%error * iv%pilot(n)%u(k)%error)
         jo_grad_y%pilot(n)%v(k) = -re%pilot(n)%v(k) / (iv%pilot(n)%v(k)%error * iv%pilot(n)%v(k)%error)
      end do
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_pilot")

end subroutine da_calculate_grady_pilot



end module da_pilot

