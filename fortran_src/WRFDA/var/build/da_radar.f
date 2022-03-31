












module da_radar

   use module_domain, only : domain

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, trace_use, &
      missing, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv,  &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, radar,fails_error_max, &
      use_radar_rv, use_radar_rf, use_radar_rhv, use_radar_rqv, &
      below_model_surface,mkz,above_model_lid,&
      fg_format,fg_format_wrf_arw_regional,fg_format_wrf_nmm_regional,fg_format_wrf_arw_global,&
      fg_format_kma_global,max_error_rv,max_error_rf, &
      far_below_model_surface,kms,kme,kts,kte, trace_use_dull,filename_len,&
      myproc, analysis_date, num_procs , ierr, comm, es_beta, es_gamma, a_ew
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      infa_type, field_type
   use da_interpolation, only : da_to_zk, da_interp_lin_3d,da_interp_lin_3d_adj
   use da_par_util, only :da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_residual, map_info, da_llxy_wrf, da_llxy_default, da_convert_zk
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_reporting, only : da_error, da_warning, da_message, message
   use da_tools_serial, only : da_get_unit, da_free_unit

   

   type residual_radar1_type
      real                    :: rv
      real                    :: rf
      real                    :: rrn
      real                    :: rsn
      real                    :: rgr
      real                    :: rcl
      real                    :: rci
      real                    :: rqv
   end type residual_radar1_type

   type maxmin_radar_stats_type
      type (maxmin_type)         :: rv       
      type (maxmin_type)         :: rf       
      type (maxmin_type)         :: rrn
      type (maxmin_type)         :: rsn
      type (maxmin_type)         :: rgr
      type (maxmin_type)         :: rcl
      type (maxmin_type)         :: rci
      type (maxmin_type)         :: rqv
   end type maxmin_radar_stats_type

   type stats_radar_type
      type (maxmin_radar_stats_type)  :: maximum, minimum
      type (residual_radar1_type)     :: average, rms_err
   end type stats_radar_type

   real, parameter :: leh1=43.1
   real, parameter :: leh2=17.5

contains

subroutine da_ao_stats_radar (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type (y_type),  intent (in)    :: re            ! A - O

   type (stats_radar_type) :: stats
   integer                 :: nrv, nrf, nrrn, nrsn, nrgr,nrqv
   integer                 :: n, k

   if (trace_use) call da_trace_entry("da_ao_stats_radar")

   nrv = 0
   nrf = 0
   nrrn= 0
   nrsn= 0
   nrgr= 0
   nrqv= 0

   stats%maximum%rv = maxmin_type (missing_r, 0, 0)
   stats%maximum%rf = maxmin_type (missing_r, 0, 0)
   stats%minimum%rv = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rf = maxmin_type(-missing_r, 0, 0)

   stats%maximum%rrn = maxmin_type (missing_r, 0, 0)
   stats%maximum%rsn = maxmin_type (missing_r, 0, 0)
   stats%maximum%rgr = maxmin_type (missing_r, 0, 0)
   stats%maximum%rcl = maxmin_type (missing_r, 0, 0)
   stats%maximum%rci = maxmin_type (missing_r, 0, 0)
   stats%maximum%rqv = maxmin_type (missing_r, 0, 0)

   stats%minimum%rrn = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rsn = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rgr = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rcl = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rci = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rqv = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_radar1_type(0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(radar)%nlocal
      if (iv%info(radar)%proc_domain(1,n)) then
         do k=1, iv%info(radar)%levels(n)
            if (use_radar_rv) then
               call da_stats_calculate (n, k, iv%radar(n)%rv(k)%qc, & 
                  re%radar(n)%rv(k), nrv, & 
                  stats%minimum%rv, stats%maximum%rv, &
                  stats%average%rv, stats%rms_err%rv)
            end if

            if (use_radar_rf) then
               call da_stats_calculate (n, k, iv%radar(n)%rf(k)%qc, & 
                  re%radar(n)%rf(k), nrf, & 
                  stats%minimum%rf, stats%maximum%rf, &
                  stats%average%rf, stats%rms_err%rf)
            end if

            if (.not. use_radar_rf .and. use_radar_rhv) then
               call da_stats_calculate (n, k, iv%radar(n)%rrn(k)%qc, &
                  re%radar(n)%rrn(k), nrrn, &
                  stats%minimum%rrn, stats%maximum%rrn, &
                  stats%average%rrn, stats%rms_err%rrn)
              call da_stats_calculate (n, k, iv%radar(n)%rsn(k)%qc, &
                  re%radar(n)%rsn(k), nrsn, &
                  stats%minimum%rsn, stats%maximum%rsn, &
                  stats%average%rsn, stats%rms_err%rsn)
               call da_stats_calculate (n, k, iv%radar(n)%rgr(k)%qc, &
                  re%radar(n)%rgr(k), nrgr, &
                  stats%minimum%rgr, stats%maximum%rgr, &
                  stats%average%rgr, stats%rms_err%rgr)
            end if

            if (use_radar_rqv) then
               call da_stats_calculate (n, k, iv%radar(n)%rqv(k)%qc, &
                  re%radar(n)%rqv(k), nrqv, &
                  stats%minimum%rqv, stats%maximum%rqv, &
                  stats%average%rqv, stats%rms_err%rqv)
            end if
         end do
      end if    
   end do
   ! Do inter-processor communication to gather statistics.
   if (use_radar_rv) then
      call da_proc_sum_int (nrv)
      call da_proc_stats_combine(stats%average%rv, stats%rms_err%rv, &
         stats%minimum%rv%value, stats%maximum%rv%value, &
         stats%minimum%rv%n, stats%maximum%rv%n, &
         stats%minimum%rv%l, stats%maximum%rv%l)
   end if

   if (use_radar_rf) then
      call da_proc_sum_int (nrf)
      call da_proc_stats_combine(stats%average%rf, stats%rms_err%rf, &
          stats%minimum%rf%value, stats%maximum%rf%value, &
          stats%minimum%rf%n, stats%maximum%rf%n, &
          stats%minimum%rf%l, stats%maximum%rf%l)
   end if

   if (.not.use_radar_rf .and. use_radar_rhv) then
      call da_proc_sum_int (nrrn)
      call da_proc_stats_combine(stats%average%rrn, stats%rms_err%rrn, &
          stats%minimum%rrn%value, stats%maximum%rrn%value, &
          stats%minimum%rrn%n, stats%maximum%rrn%n, &
          stats%minimum%rrn%l, stats%maximum%rrn%l)
      call da_proc_sum_int (nrsn)
      call da_proc_stats_combine(stats%average%rsn, stats%rms_err%rsn, &
          stats%minimum%rsn%value, stats%maximum%rsn%value, &
          stats%minimum%rsn%n, stats%maximum%rsn%n, &
          stats%minimum%rsn%l, stats%maximum%rsn%l)
      call da_proc_sum_int (nrgr)
      call da_proc_stats_combine(stats%average%rgr, stats%rms_err%rgr, &
          stats%minimum%rgr%value, stats%maximum%rgr%value, &
          stats%minimum%rgr%n, stats%maximum%rgr%n, &
          stats%minimum%rgr%l, stats%maximum%rgr%l)
   end if

   if (use_radar_rqv) then
      call da_proc_sum_int (nrqv)
      call da_proc_stats_combine(stats%average%rqv, stats%rms_err%rqv, &
          stats%minimum%rqv%value, stats%maximum%rqv%value, &
          stats%minimum%rqv%n, stats%maximum%rqv%n, &
          stats%minimum%rqv%l, stats%maximum%rqv%l)
   end if
   iv%nstats(radar) = nrv + nrf + nrrn + nrsn + nrgr + nrqv

   if (rootproc) then
      if (nrv /= 0 .or. nrf /= 0 .or. nrrn /= 0 .or. nrsn /= 0 .or. nrgr /= 0 .or. nrqv /= 0) then 
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for radar'
         call da_print_stats_radar(stats_unit, nrv, nrf, nrrn, nrsn, nrgr, nrqv, stats)
      end if
   end if

   if (trace_use) call da_trace_exit("da_ao_stats_radar")

end subroutine da_ao_stats_radar


subroutine da_jo_and_grady_radar(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type),intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use) call da_trace_entry("da_jo_and_grady_radar")

   jo % radar_rv = 0.0
   jo % radar_rf = 0.0
   jo % radar_rrn = 0.0
   jo % radar_rsn = 0.0
   jo % radar_rgr = 0.0
   jo % radar_rci = 0.0
   jo % radar_rcl = 0.0
   jo % radar_rqv = 0.0 

   do n=1, iv%info(radar)%nlocal
      do k=1, iv%info(radar)%levels(n)
         if (use_radar_rv) then
            jo_grad_y%radar(n)%rv(k) = -re%radar(n)%rv(k) / (iv%radar(n)%rv(k)%error * iv%radar(n)%rv(k)%error) 
         end if

         if (use_radar_rf) then
            jo_grad_y%radar(n)%rf(k) = -re%radar(n)%rf(k) / (iv%radar(n)%rf(k)%error * iv%radar(n)%rf(k)%error) 
         end if

         if (use_radar_rhv) then
            jo_grad_y%radar(n)%rrn(k) = -re%radar(n)%rrn(k) / (iv%radar(n)%rrn(k)%error * iv%radar(n)%rrn(k)%error) 
            jo_grad_y%radar(n)%rsn(k) = -re%radar(n)%rsn(k) / (iv%radar(n)%rsn(k)%error * iv%radar(n)%rsn(k)%error) 
            jo_grad_y%radar(n)%rgr(k) = -re%radar(n)%rgr(k) / (iv%radar(n)%rgr(k)%error * iv%radar(n)%rgr(k)%error) 
         end if

         if (use_radar_rqv) then
            jo_grad_y%radar(n)%rqv(k) = -re%radar(n)%rqv(k) / (iv%radar(n)%rqv(k)%error * iv%radar(n)%rqv(k)%error) 
         end if
      end do

      if (iv%info(radar)%proc_domain(1,n)) then
         do k=1, iv%info(radar)%levels(n)
            if (use_radar_rv) then
               jo % radar_rv = jo % radar_rv - re%radar(n)%rv(k) * jo_grad_y%radar(n)%rv(k)
            end if

            if (use_radar_rf) then
               jo % radar_rf = jo % radar_rf - re%radar(n)%rf(k) * jo_grad_y%radar(n)%rf(k)
            end if
       
            if (use_radar_rhv) then
               jo % radar_rrn = jo % radar_rrn - re%radar(n)%rrn(k) * jo_grad_y%radar(n)%rrn(k)
               jo % radar_rsn = jo % radar_rsn - re%radar(n)%rsn(k) * jo_grad_y%radar(n)%rsn(k)
               jo % radar_rgr = jo % radar_rgr - re%radar(n)%rgr(k) * jo_grad_y%radar(n)%rgr(k)
            end if

            if (use_radar_rqv) then
               jo % radar_rqv = jo % radar_rqv - re%radar(n)%rqv(k) * jo_grad_y%radar(n)%rqv(k)
            end if
         end do
      end if
   end do
      
   jo % radar_rv = 0.5 * jo % radar_rv
   jo % radar_rf = 0.5 * jo % radar_rf
   jo % radar_rrn= 0.5 * jo % radar_rrn
   jo % radar_rsn= 0.5 * jo % radar_rsn
   jo % radar_rgr= 0.5 * jo % radar_rgr
   jo % radar_rqv= 0.5 * jo % radar_rqv

   if (trace_use) call da_trace_exit("da_jo_and_grady_radar")

end subroutine da_jo_and_grady_radar


subroutine da_residual_radar(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for radar obs
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

   if (trace_use) call da_trace_entry("da_residual_radar")

   n_obs_bad % rv % num = number_type(0, 0, 0)
   n_obs_bad % rf % num = number_type(0, 0, 0)
   n_obs_bad % rrn % num = number_type(0, 0, 0)
   n_obs_bad % rsn % num = number_type(0, 0, 0)
   n_obs_bad % rgr % num = number_type(0, 0, 0)
   n_obs_bad % rqv % num = number_type(0, 0, 0)

   do n=1, iv%info(radar)%nlocal
      do k=1, iv%info(radar)%levels(n)
         if (use_radar_rv) then
            np_available = np_available + 1
            re%radar(n)%rv(k) = da_residual(n, k, y%radar(n)%rv(k), iv%radar(n)%rv(k), n_obs_bad % rv)
         end if

         if (use_radar_rf) then
            np_available = np_available + 1
            if (.not. use_radar_rhv) then
               re%radar(n)%rf(k) = da_residual(n, k, y%radar(n)%rf(k), iv%radar(n)%rf(k), n_obs_bad % rf)
            end if
         end if

        if (.not.use_radar_rf .and. use_radar_rhv) then
           re%radar(n)%rrn(k) = da_residual(n, k, y%radar(n)%rrn(k), iv%radar(n)%rrn(k), n_obs_bad % rrn)
           re%radar(n)%rsn(k) = da_residual(n, k, y%radar(n)%rsn(k), iv%radar(n)%rsn(k), n_obs_bad % rsn)
           re%radar(n)%rgr(k) = da_residual(n, k, y%radar(n)%rgr(k), iv%radar(n)%rgr(k), n_obs_bad % rgr)
        end if

        if (use_radar_rqv) then
           re%radar(n)%rqv(k) = da_residual(n, k, y%radar(n)%rqv(k), iv%radar(n)%rqv(k), n_obs_bad % rqv)
        end if
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % rv % num % miss + n_obs_bad % rf % num % miss + &
                               n_obs_bad % rrn% num % miss + n_obs_bad % rsn% num % miss + &
                               n_obs_bad % rqv% num % miss + n_obs_bad % rqv% num % miss 
   np_bad_data = np_bad_data + n_obs_bad % rv % num % bad  + n_obs_bad % rf % num % bad  + &
                               n_obs_bad % rrn% num % bad  + n_obs_bad % rsn% num % bad  + &
                               n_obs_bad % rgr% num % bad  + n_obs_bad % rqv% num % bad
   np_obs_used = np_obs_used + n_obs_bad % rv % num % use  + n_obs_bad % rf % num % use  + &
                               n_obs_bad % rrn% num % use  + n_obs_bad % rsn% num % use  + &
                               n_obs_bad % rgr% num % use  + n_obs_bad % rqv% num % use 

   if (trace_use) call da_trace_exit("da_residual_radar")

end subroutine da_residual_radar


subroutine da_oi_stats_radar (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_radar_type) :: stats
   integer                 :: nrv, nrf, nrrn,nrsn,nrgr,nrqv
   integer                 :: n, k

   if (trace_use) call da_trace_entry("da_oi_stats_radar")

   nrv = 0
   nrf = 0
   nrrn= 0
   nrsn= 0
   nrgr= 0
   nrqv= 0
     
   stats%maximum%rv = maxmin_type(missing_r, 0, 0)
   stats%minimum%rv = maxmin_type(-missing_r, 0, 0)
   stats%maximum%rf = maxmin_type(missing_r, 0, 0)
   stats%minimum%rf = maxmin_type(-missing_r, 0, 0)

   stats%maximum%rrn = maxmin_type (missing_r, 0, 0)
   stats%maximum%rsn = maxmin_type (missing_r, 0, 0)
   stats%maximum%rgr = maxmin_type (missing_r, 0, 0)
   stats%maximum%rcl = maxmin_type (missing_r, 0, 0)
   stats%maximum%rci = maxmin_type (missing_r, 0, 0)
   stats%maximum%rqv = maxmin_type (missing_r, 0, 0)

   stats%minimum%rrn = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rsn = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rgr = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rcl = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rci = maxmin_type(-missing_r, 0, 0)
   stats%minimum%rqv = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_radar1_type(0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0)

   stats%rms_err = stats%average

   do n=1, iv%info(radar)%nlocal
      if (iv%info(radar)%proc_domain(1,n)) then
         do k=1, iv%info(radar)%levels(n)
            if (use_radar_rv) then
               call da_stats_calculate(iv%info(radar)%obs_global_index(n), &
                  k, iv%radar(n)%rv(k)%qc, &
                  iv%radar(n)%rv(k)%inv, nrv, &
                  stats%minimum%rv, stats%maximum%rv, &
                  stats%average%rv, stats%rms_err%rv)
            end if

            if (use_radar_rf) then
               call da_stats_calculate(iv%info(radar)%obs_global_index(n), &
                  k, iv%radar(n)%rf(k)%qc, &
                  iv%radar(n)%rf(k)%inv, nrf, &
                  stats%minimum%rf, stats%maximum%rf, &
                  stats%average%rf, stats%rms_err%rf)
            end if

            if (.not. use_radar_rf .and. use_radar_rhv) then
               call da_stats_calculate(iv%info(radar)%obs_global_index(n), &
                  k, iv%radar(n)%rrn(k)%qc, &
                  iv%radar(n)%rrn(k)%inv, nrrn, &
                  stats%minimum%rrn, stats%maximum%rrn, &
                  stats%average%rrn, stats%rms_err%rrn)
               call da_stats_calculate(iv%info(radar)%obs_global_index(n), &
                  k, iv%radar(n)%rsn(k)%qc, &
                  iv%radar(n)%rsn(k)%inv, nrsn, &
                  stats%minimum%rsn, stats%maximum%rsn, &
                  stats%average%rsn, stats%rms_err%rsn)
               call da_stats_calculate(iv%info(radar)%obs_global_index(n), &
                  k, iv%radar(n)%rgr(k)%qc, &
                  iv%radar(n)%rgr(k)%inv, nrgr, &
                  stats%minimum%rgr, stats%maximum%rgr, &
                  stats%average%rgr, stats%rms_err%rgr)
            end if

            if (use_radar_rqv) then
               call da_stats_calculate(iv%info(radar)%obs_global_index(n), &
                  k, iv%radar(n)%rqv(k)%qc, &
                  iv%radar(n)%rqv(k)%inv, nrqv, &
                  stats%minimum%rqv, stats%maximum%rqv, &
                  stats%average%rqv, stats%rms_err%rqv)
            end if

         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.

   if (use_radar_rv) then
      call da_proc_sum_int(nrv)
      call da_proc_stats_combine(stats%average%rv, stats%rms_err%rv, &
         stats%minimum%rv%value, stats%maximum%rv%value, &
         stats%minimum%rv%n, stats%maximum%rv%n, &
         stats%minimum%rv%l, stats%maximum%rv%l)
   end if
   
   if (use_radar_rf) then
      call da_proc_sum_int(nrf)
      call da_proc_stats_combine(stats%average%rf, stats%rms_err%rf, &
         stats%minimum%rf%value, stats%maximum%rf%value, &
         stats%minimum%rf%n, stats%maximum%rf%n, &
         stats%minimum%rf%l, stats%maximum%rf%l)
   end if

   if (.not.use_radar_rf .and. use_radar_rhv) then
      call da_proc_sum_int(nrrn)
      call da_proc_stats_combine(stats%average%rrn, stats%rms_err%rrn, &
         stats%minimum%rrn%value, stats%maximum%rrn%value, &
         stats%minimum%rrn%n, stats%maximum%rrn%n, &
         stats%minimum%rrn%l, stats%maximum%rrn%l)
      call da_proc_sum_int(nrsn)
      call da_proc_stats_combine(stats%average%rsn, stats%rms_err%rsn, &
         stats%minimum%rsn%value, stats%maximum%rsn%value, &
         stats%minimum%rsn%n, stats%maximum%rsn%n, &
         stats%minimum%rsn%l, stats%maximum%rsn%l)
      call da_proc_stats_combine(stats%average%rgr, stats%rms_err%rgr, &
         stats%minimum%rgr%value, stats%maximum%rgr%value, &
         stats%minimum%rgr%n, stats%maximum%rgr%n, &
         stats%minimum%rgr%l, stats%maximum%rgr%l)
   end if

   if (use_radar_rqv) then
      call da_proc_sum_int(nrqv)
      call da_proc_stats_combine(stats%average%rqv, stats%rms_err%rqv, &
         stats%minimum%rqv%value, stats%maximum%rqv%value, &
         stats%minimum%rqv%n, stats%maximum%rqv%n, &
         stats%minimum%rqv%l, stats%maximum%rqv%l)
   end if

   if (rootproc) then
      if (nrv /= 0 .or. nrf /= 0 .or. nrrn /= 0 .or. nrsn /= 0 .or. nrgr /= 0 .or. nrqv /= 0 ) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for radar'
         call da_print_stats_radar(stats_unit, nrv, nrf, nrrn, nrsn, nrgr, nrqv,stats)
      end if
   end if

   if (trace_use) call da_trace_exit("da_oi_stats_radar")

end subroutine da_oi_stats_radar


subroutine da_print_stats_radar(stats_unit, nrv, nrf, nrrn, nrsn, nrgr, nrqv, radar)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nrv, nrf, nrrn, nrsn, nrgr,nrqv
   type (stats_radar_type), intent(in)    :: radar

   if (trace_use) call da_trace_entry("da_print_stats_radar")
   
   write(unit=stats_unit, fmt='(a/)') &
      '   var             rv (m/s)        n    k         rf (dBZ)        n    k        rrn(kg/kg)       n    k        rsn(kg/kg)       n    k        rgr(kg/kg)       n    k        rqv(kg/kg)       n    k'

   write(unit=stats_unit, fmt='(a,i16,5(i31))') &
      '  Number: ', nrv, nrf, nrrn, nrsn, nrgr, nrqv

   if (nrv < 1) nrv = 1
   if (nrf < 1) nrf = 1
   if (nrrn < 1) nrrn = 1
   if (nrsn < 1) nrsn = 1
   if (nrgr < 1) nrgr = 1
   if (nrqv < 1) nrqv = 1
   
   write(unit=stats_unit, fmt='((a,f12.4,i9,i5, 5(f17.4,i9,i5)))')  &
      ' Minimum(n,k): ', radar%minimum%rv, radar%minimum%rf, radar%minimum%rrn, radar%minimum%rsn,radar%minimum%rgr, radar%minimum%rqv
   write(unit=stats_unit, fmt='((a,f12.4,i9,i5, 5(f17.4,i9,i5)))') &
      ' Maximum(n,k): ', radar%maximum%rv, radar%maximum%rf, radar%maximum%rrn, radar%maximum%rsn,radar%maximum%rgr, radar%maximum%rqv
   write(unit=stats_unit, fmt='((a,6(f12.4,19x)))')              &
      ' Average     : ', radar%average%rv/real(nrv), radar%average%rf/real(nrf),radar%average%rrn/real(nrrn), radar%average%rsn/real(nrsn),radar%average%rgr/real(nrgr), radar%average%rqv/real(nrqv) 
   write(unit=stats_unit, fmt='((a,6(f12.4,19x)))')              &
      '    RMSE     : ', sqrt(radar%rms_err%rv/real(nrv)), sqrt(radar%rms_err%rf/real(nrf)), &
                         sqrt(radar%rms_err%rrn/real(nrrn)),sqrt(radar%rms_err%rsn/real(nrsn)),sqrt(radar%rms_err%rgr/real(nrgr)),sqrt(radar%rms_err%rqv/real(nrqv))

   if (trace_use) call da_trace_exit("da_print_stats_radar")

end subroutine da_print_stats_radar


subroutine da_transform_xtoy_radar (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: calculate the Doppler radial velocity and 
   ! reflectivity at the observation location from the first guess.
   ! It is linearized. 
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------
 
   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n                   ! Loop counter.      
   integer :: k                         ! Index dimension.   

   real, allocatable :: model_p(:,:)
   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_w(:,:)
   real, allocatable :: model_rho(:,:)
   real, allocatable :: model_qrn(:,:)
   real, allocatable :: model_qrnb(:,:)
   real, allocatable :: model_qrnb2(:,:)
   real, allocatable :: model_ps(:)
   real, allocatable :: model_qsn(:,:)
   real, allocatable :: model_qgr(:,:)
   real, allocatable :: model_qsnb(:,:)
   real, allocatable :: model_qgrb(:,:)
   real, allocatable :: model_qv(:,:)
   real, allocatable :: model_qvb(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_tb(:,:)
   real, allocatable :: model_tc(:,:)

   real    :: xr,yr,zr

   real    :: alog_10

   if (trace_use) call da_trace_entry("da_transform_xtoy_radar")

   alog_10 = alog(10.0)

   allocate (model_p(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_u(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_v(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_w(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_rho(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrnb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrnb2(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_ps(iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qsn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgr(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qsnb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgrb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qv(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qvb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_t(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_tb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_tc(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))

   do n=iv%info(radar)%n1,iv%info(radar)%n2
      do k = 1, iv%info(radar)%levels(n)
         model_qrnb(k,n) = iv%radar(n)%model_qrn(k)
         model_p(k,n)    = iv%radar(n)%model_p(k)
      end do

      model_ps(n) = iv%radar(n)%model_ps
   end do

   ! [1.4] Interpolate horizontally from dot points:
   call da_interp_lin_3d (grid%xa%u,   iv%info(radar), model_u)
   call da_interp_lin_3d (grid%xa%v,   iv%info(radar), model_v)
   call da_interp_lin_3d (grid%xa%qrn, iv%info(radar), model_qrn)
   call da_interp_lin_3d (grid%xa%wh,  iv%info(radar), model_w)
   call da_interp_lin_3d (grid%xa%qsn, iv%info(radar), model_qsn)
   call da_interp_lin_3d (grid%xa%qgr, iv%info(radar), model_qgr)
   call da_interp_lin_3d (grid%xa%q,   iv%info(radar), model_qv)
   call da_interp_lin_3d (grid%xa%t,   iv%info(radar), model_t)
   !basic states
   call da_interp_lin_3d (grid%xb % rho, iv%info(radar), model_rho)
   call da_interp_lin_3d (grid%xb%t,   iv%info(radar), model_tb)
   call da_interp_lin_3d (grid%xb%q,   iv%info(radar), model_qvb)
   call da_interp_lin_3d (grid%xb%qrn, iv%info(radar), model_qrnb2)
   call da_interp_lin_3d (grid%xb%qsn, iv%info(radar), model_qsnb)
   call da_interp_lin_3d (grid%xb%qgr, iv%info(radar), model_qgrb)

   model_tc = model_tb - 273.15

   do n=iv%info(radar)%n1,iv%info(radar)%n2

      ! [1.7] Calculate rv and rf at OBS location

      xr = grid%xb%ds * (iv%info(radar)%x(1,n) - iv%radar(n)%stn_loc%x)
      yr = grid%xb%ds * (iv%info(radar)%y(1,n) - iv%radar(n)%stn_loc%y)

      do k = 1, iv%info(radar)%levels(n)
         if (iv % radar(n) % height_qc(k) /= below_model_surface .and.  &
              iv % radar(n) % height_qc(k) /= above_model_lid) then
            if (use_radar_rv) then
               if (iv % radar(n) % rv(k) % qc >= obs_qc_pointer) then
                  zr=iv%radar(n)%height(k) - iv%radar(n)%stn_loc%elv

                  call da_radial_velocity_lin(y%radar(n)%rv(k), &
                     model_p(k,n), &
                     model_u(k,n), model_v(k,n), model_w(k,n), model_qrn(k,n),    &
                     model_ps(n), xr, yr, zr, model_qrnb(k,n))
               end if
            end if

            if (use_radar_rf .and. .not. use_radar_rhv) then
               if (iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then
                  y%radar(n)%rf(k) = leh2 * model_qrn(k,n) /(model_qrnb(k,n)*alog_10) 
               end if
            end if

            if (.not.use_radar_rf .and. use_radar_rhv) then
               if (iv % radar(n) % rrn(k) % qc >= obs_qc_pointer) then
                  y%radar(n)%rrn(k) = model_qrn(k,n) 
               end if
               if (iv % radar(n) % rsn(k) % qc >= obs_qc_pointer) then
                  y%radar(n)%rsn(k) = model_qsn(k,n) 
               end if
               if (iv % radar(n) % rgr(k) % qc >= obs_qc_pointer) then
                  y%radar(n)%rgr(k) = model_qgr(k,n)
               end if
            end if

            if (use_radar_rqv) then
               !dqv=qs*drh+(c2*c3*/(T+c3)**2.0-c4/T)*qv**dT
               !c2=17.67 is es_beta
               !c3=243.5 is es_gamma
               !c4=0.622 is a_ew
               ! use qc from get_inv.
               if (iv % radar(n) % rqv(k) % qc >= obs_qc_pointer) then
                  y%radar(n)%rqv(k) = model_qv(k,n)
                  ! Wang JAMC
                  y%radar(n)%rqv(k) = y%radar(n)%rqv(k) + ( es_beta*es_gamma/(model_tb(k,n)+es_gamma)**2.0 )*model_qvb(k,n)*model_t(k,n)
               end if
            end if

         end if
      end do
   end do

   deallocate (model_p)
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_w)
   deallocate (model_qrn)
   deallocate (model_qrnb)
   deallocate (model_ps)

   if (trace_use) call da_trace_exit("da_transform_xtoy_radar")

end subroutine da_transform_xtoy_radar 


subroutine da_transform_xtoy_radar_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   !------------------------------------------------------------------------
   ! This subroutine is the adjoint of Doppler radar observation operators.
   !------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: k  ! Index dimension.

   integer :: n

   real, allocatable :: model_p(:,:)
   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_w(:,:)
   real, allocatable :: model_rho(:,:)
   real, allocatable :: model_qrn(:,:)
   real, allocatable :: model_qrnb(:,:)
   real, allocatable :: model_qrnb2(:,:)
   real, allocatable :: model_ps(:)

   real, allocatable :: model_qsn(:,:)
   real, allocatable :: model_qgr(:,:)
   real, allocatable :: model_qsnb(:,:)
   real, allocatable :: model_qgrb(:,:)
   real, allocatable :: model_qv(:,:)
   real, allocatable :: model_qvb(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_tb(:,:)
   real, allocatable :: model_tc(:,:)

   real    :: xr,yr,zr

   real    :: alog10,qrn1,qsn1,qgr1

   if (trace_use) call da_trace_entry("da_transform_xtoy_radar_adj")

   alog10= alog(10.0)

   allocate (model_p(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_u(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_v(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_w(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrnb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrnb2(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_ps(iv%info(radar)%n1:iv%info(radar)%n2))

   allocate (model_rho(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qsn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgr(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qsnb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgrb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qv(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qvb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_t(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_tb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_tc(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))

   !basic states
   call da_interp_lin_3d (grid%xb % rho, iv%info(radar), model_rho)
   call da_interp_lin_3d (grid%xb%t,   iv%info(radar), model_tb)
   call da_interp_lin_3d (grid%xb%q,   iv%info(radar), model_qvb)
   call da_interp_lin_3d (grid%xb%qrn, iv%info(radar), model_qrnb2)
   call da_interp_lin_3d (grid%xb%qsn, iv%info(radar), model_qsnb)
   call da_interp_lin_3d (grid%xb%qgr, iv%info(radar), model_qgrb)

   model_tc = model_tb - 273.15

   ! Needed
   model_u = 0.0
   model_v = 0.0
   model_w = 0.0
   model_qrn = 0.0
   model_qsn = 0.0
   model_qgr = 0.0
   model_qv  = 0.0
   model_t   = 0.0


   ! W_HALF is vertical velocity at half-sigma levels.

   model_ps(iv%info(radar)%n1:iv%info(radar)%n2) = iv%radar(iv%info(radar)%n1:iv%info(radar)%n2)%model_ps 

   do n=iv%info(radar)%n1,iv%info(radar)%n2

      ! [1.7] Calculate rv and rf at OBS location

      xr = grid%xb%ds * (iv%info(radar)%x(1,n) - iv%radar(n)%stn_loc%x)
      yr = grid%xb%ds * (iv%info(radar)%y(1,n) - iv%radar(n)%stn_loc%y)

      model_qrnb(1:iv%info(radar)%levels(n),n) = iv%radar(n)%model_qrn(1:iv%info(radar)%levels(n))
      model_p   (1:iv%info(radar)%levels(n),n) = iv%radar(n)%model_p(1:iv%info(radar)%levels(n))

      do k = 1,iv%info(radar)%levels(n)
         if (iv % radar(n) % height_qc(k) /= below_model_surface .and.  &
              iv % radar(n) % height_qc(k) /= above_model_lid) then

            if (use_radar_rf .and. .not.use_radar_rhv) then
               if (iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then
                  model_qrn(k,n) = model_qrn(k,n) + leh2/(model_qrnb(k,n)*alog10) * jo_grad_y%radar(n)%rf(k)
               end if
            end if

            if (.not.use_radar_rf .and. use_radar_rhv) then
               if (iv % radar(n) % rrn(k) % qc >= obs_qc_pointer) then
                  model_qrn(k,n) = model_qrn(k,n) + jo_grad_y%radar(n)%rrn(k)
               end if
               if (iv % radar(n) % rsn(k) % qc >= obs_qc_pointer) then
                  model_qsn(k,n) = model_qsn(k,n) + jo_grad_y%radar(n)%rsn(k)
               end if
               if (iv % radar(n) % rgr(k) % qc >= obs_qc_pointer) then
                  model_qgr(k,n) = model_qgr(k,n) + jo_grad_y%radar(n)%rgr(k)
               end if
            end if

            if (use_radar_rqv) then
               if (iv % radar(n) % rqv(k) % qc >= obs_qc_pointer) then
               !TL  y%radar(n)%rqv(k) = y%radar(n)%rqv(k) +( es_beta*es_gamma/(model_tb(k,n)+es_gamma)**2.0 )*model_qvb(k,n)*model_t(k,n)
                  model_qv(k,n) = model_qv(k,n) + jo_grad_y%radar(n)%rqv(k)
                  model_t(k,n)  = model_t(k,n)  + (es_beta*es_gamma/(model_tb(k,n)+es_gamma)**2.0)*model_qvb(k,n)*jo_grad_y%radar(n)%rqv(k)
               end if
            end if


            if (use_radar_rv) then
               if (iv % radar(n) % rv(k) % qc >= obs_qc_pointer) then
                  zr=iv%radar(n)%height(k) - iv%radar(n)%stn_loc%elv

                  call da_radial_velocity_adj(jo_grad_y%radar(n)%rv(k), &
                     model_p(k,n), model_u(k,n), model_v(k,n), model_w(k,n),  &
                     model_qrn(k,n), model_ps(n), xr, yr, zr, model_qrnb(k,n))

               end if
            end if
         end if
      end do
      jo_grad_y%radar(n)%rv(:) = 0.0
      jo_grad_y%radar(n)%rf(:) = 0.0
      jo_grad_y%radar(n)%rrn(:)= 0.0
      jo_grad_y%radar(n)%rsn(:)= 0.0
      jo_grad_y%radar(n)%rgr(:)= 0.0
      jo_grad_y%radar(n)%rqv(:)= 0.0
   end do ! n

   ! [1.6] Interpolate horizontally from crs points:

   call da_interp_lin_3d_adj (jo_grad_x % wh,  iv%info(radar), model_w)
   call da_interp_lin_3d_adj (jo_grad_x % qrn, iv%info(radar), model_qrn)
   call da_interp_lin_3d_adj (jo_grad_x % qsn, iv%info(radar), model_qsn)
   call da_interp_lin_3d_adj (jo_grad_x % qgr, iv%info(radar), model_qgr)
   call da_interp_lin_3d_adj (jo_grad_x % q,   iv%info(radar), model_qv)
   call da_interp_lin_3d_adj (jo_grad_x % t,   iv%info(radar), model_t)
   call da_interp_lin_3d_adj (jo_grad_x % v,   iv%info(radar), model_v)
   call da_interp_lin_3d_adj (jo_grad_x % u,   iv%info(radar), model_u)

   deallocate (model_p)
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_w)
   deallocate (model_rho)
   deallocate (model_qrn)
   deallocate (model_qrnb)
   deallocate (model_qrnb2)
   deallocate (model_ps)
   deallocate (model_qv)
   deallocate (model_qvb)
   deallocate (model_t)
   deallocate (model_tb)
   deallocate (model_qsn)
   deallocate (model_qsnb)
   deallocate (model_qgr)
   deallocate (model_qgrb)

   if (trace_use) call da_trace_exit("da_transform_xtoy_radar_adj")

end subroutine da_transform_xtoy_radar_adj


subroutine da_check_max_iv_radar(iv, it, irv, irf, irvf, irff)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it    
   integer,       intent(inout) :: irv, irf, irvf, irff    

   integer :: k,n
   logical :: failed 

   if (trace_use) call da_trace_entry("da_check_max_iv_radar")

   !---------------------------------------------------------------------------
   ! [1.0] Open diagnostic file:
   !---------------------------------------------------------------------------

   if (rootproc .and. check_max_iv_print) then
      write (check_max_iv_unit,'(a)')  &
         '----------------------------------------------------------------'
      write (unit = check_max_iv_unit, fmt = '(A,/)') 'MAX ERROR TEST QC:'

      write (unit = check_max_iv_unit, fmt = '(/,9(A,F3.0,/))')  &
         'Error max test ratio for radar_rv   = ',max_error_rv, &
         'Error max test ratio for radar_rf   = ',max_error_rf
   end if

   !------------------------------------------------------------------------
   ! [2.0] Perform maximum innovation vector check:
   !------------------------------------------------------------------------

   failed = .false.

   do n = iv%info(radar)%n1,iv%info(radar)%n2
      do k = 1, iv%info(radar)%levels(n)
         if (iv%radar(n)%height_qc(k) /= far_below_model_surface .and. &
              iv%radar(n)%height_qc(k) /= above_model_lid) then
            ! rv
            if (use_radar_rv) then
               call da_max_error_qc_radar(it, iv%info(radar), n, iv%radar(n)%rv(k), max_error_rv, irv, irvf, &
                  check_max_iv_unit, 'rv   ', failed, check_max_iv_print)
            end if

            ! rf
            if (use_radar_rf) then
               call da_max_error_qc_radar(it, iv%info(radar), n, iv%radar(n)%rf(k), max_error_rf, irf, irff, &
                  check_max_iv_unit, 'rf   ', failed, check_max_iv_print)
            end if
         end if
      end do
   end do

   if (trace_use) call da_trace_exit("da_check_max_iv_radar")

end subroutine da_check_max_iv_radar
subroutine da_get_innov_vector_radar (it, grid, ob, iv)

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

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.

   real, allocatable :: model_p(:,:)
   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_w(:,:)

   real, allocatable :: model_rv(:,:)
   real, allocatable :: model_rf(:,:)
   real, allocatable :: model_ps(:)

   real, allocatable :: model_qv(:,:)
   real, allocatable :: model_qs(:,:)
   real, allocatable :: model_tc(:,:)

   real, allocatable :: model_rho(:,:)
   real, allocatable :: model_qrn(:,:)
   real, allocatable :: model_qcl(:,:)
   real, allocatable :: model_qci(:,:)
   real, allocatable :: model_qsn(:,:)
   real, allocatable :: model_qgr(:,:)

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.

   real    :: xr,yr,zr
   integer :: irv, irvf
   integer :: irf, irff

   real    :: alog_10, czr,czs,czg, zrr,zds,zws,zg,rze
   alog_10 = alog(10.0)

   ! Ze=zv*(ro*v)**1.75
   ! Zdb=10*log10(Ze)
   zrr = 3.63*1.00e+9  ! rainwater
   zds = 9.80*1.00e+8  ! dry snow
   zws = 4.26*1.00e+11 ! wet snow
   zg  = 4.33*1.00e+10 ! grauple

   
   if (trace_use) call da_trace_entry("da_get_innov_vector_radar")

   irv = 0; irvf = 0; irf = 0; irff = 0


   ! No point in going through and allocating all these variables if we're just going to quit anyway

   if ( use_radar_rf .and. use_radar_rhv ) then
      write(unit=message(1),fmt='(A)') "Both 'use_radar_rf' and 'use_radar_rhv' are set to true"
      write(unit=message(2),fmt='(A)') "You must only choose one of these options"
      call da_error("da_get_innov_vector_radar.inc",68,message(1:2))
   end if


   allocate (model_p(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_u(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_v(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_w(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))

   allocate (model_rv(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_rf(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_ps(iv%info(radar)%n1:iv%info(radar)%n2))

   allocate (model_qv(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qs(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_tc(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))

   allocate (model_rho(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qcl(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qci(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qsn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgr(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))

   model_p(:,:)    = 0.
   model_u(:,:)    = 0.
   model_v(:,:)    = 0.
   model_w(:,:)    = 0.

   model_rv(:,:)   = 0.
   model_rf(:,:)   = 0.
   model_ps(:)     = 0.

   model_qv(:,:)   = 0.
   model_qs(:,:)   = 0.
   model_tc(:,:)   = 0.

   model_rho(:,:)  = 0.
   model_qrn(:,:)  = 0.
   model_qcl(:,:)  = 0.
   model_qci(:,:)  = 0.
   model_qsn(:,:)  = 0.
   model_qgr(:,:)  = 0.

   if ( it > 1 ) then
     do n=iv%info(radar)%n1,iv%info(radar)%n2
        do k=1,iv%info(radar)%levels(n)
           if (iv%radar(n)%rv(k)%qc == fails_error_max) iv%radar(n)%rv(k)%qc = 0
           if (iv%radar(n)%rf(k)%qc == fails_error_max) iv%radar(n)%rf(k)%qc = 0
        end do
     end do
   end if

   do n=iv%info(radar)%n1,iv%info(radar)%n2
      if (iv%info(radar)%levels(n) < 1) cycle

      ! [1.0] Get horizontal interpolation weights:

      i   = iv%info(radar)%i(1,n)
      j   = iv%info(radar)%j(1,n)
      dx  = iv%info(radar)%dx(1,n)
      dy  = iv%info(radar)%dy(1,n)
      dxm = iv%info(radar)%dxm(1,n)
      dym = iv%info(radar)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j,k)+dx*grid%xb%h(i+1,j,k)) + dy*(dxm*grid%xb%h(i,j+1,k)+dx*grid%xb%h(i+1,j+1,k))
      end do

      do k=1, iv%info(radar)%levels(n)
         call da_to_zk(iv%radar(n)%height(k), v_h, v_interp_h, iv%info(radar)%zk(k,n))

         if (iv%info(radar)%zk(k,n) < 1.0) then
            iv%radar(n)%height_qc(k) = below_model_surface
         else if (iv%info(radar)%zk(k,n) > mkz) then
            iv%radar(n)%height_qc(k) = above_model_lid
         end if
      end do
   end do

   call da_convert_zk (iv%info(radar))

   ! [2.0] Interpolate horizontally to ob points:

   call da_interp_lin_3d (grid%xb % p,   iv%info(radar), model_p)
   call da_interp_lin_3d (grid%xb % u,   iv%info(radar), model_u)
   call da_interp_lin_3d (grid%xb % v,   iv%info(radar), model_v)
   call da_interp_lin_3d (grid%xb % wh,  iv%info(radar), model_w)
   call da_interp_lin_3d (grid%xb % rho, iv%info(radar), model_rho)
   call da_interp_lin_3d (grid%xb % qrn, iv%info(radar), model_qrn)
   call da_interp_lin_3d (grid%xb % qcw, iv%info(radar), model_qcl)
   call da_interp_lin_3d (grid%xb % qci, iv%info(radar), model_qci)
   call da_interp_lin_3d (grid%xb % qsn, iv%info(radar), model_qsn)
IF ( ASSOCIATED( grid%xb%qgr ) ) THEN
   call da_interp_lin_3d (grid%xb % qgr, iv%info(radar), model_qgr)
END IF
   call da_interp_lin_3d (grid%xb % q,  iv%info(radar), model_qv)
   call da_interp_lin_3d (grid%xb % qs, iv%info(radar), model_qs)
   call da_interp_lin_3d (grid%xb % t,  iv%info(radar), model_tc)
   ! whl to TC
   model_tc = model_tc - 273.15

   ! Test 5.0E-5 as critical value. It can not be smaller.
   do n=iv%info(radar)%n1,iv%info(radar)%n2
      do k=1,iv%info(radar)%levels(n)
         ! for Xiao's default scheme
         if(use_radar_rf) model_qrn(k,n)=amax1(5.0E-5,model_qrn(k,n))
      end do

      i   = iv%info(radar)%i(1,n)
      j   = iv%info(radar)%j(1,n)
      dx  = iv%info(radar)%dx(1,n)
      dy  = iv%info(radar)%dy(1,n)
      dxm = iv%info(radar)%dxm(1,n)
      dym = iv%info(radar)%dym(1,n)


      model_ps(n) = dxm *(dym * grid%xb % psac(i,  j) + dy * grid%xb%psac(i+1,  j)) + &
                 dx  *(dym * grid%xb % psac(i,j+1) + dy * grid%xb%psac(i+1,j+1)) + &
                 grid%xb % ptop

      iv%radar(n)%model_p(1:iv%info(radar)%levels(n))   = model_p(1:iv%info(radar)%levels(n),n)
      iv%radar(n)%model_rho(1:iv%info(radar)%levels(n)) = model_rho(1:iv%info(radar)%levels(n),n)
      iv%radar(n)%model_qrn(1:iv%info(radar)%levels(n)) = model_qrn(1:iv%info(radar)%levels(n),n)
      iv%radar(n)%model_qcl(1:iv%info(radar)%levels(n)) = model_qcl(1:iv%info(radar)%levels(n),n)
      iv%radar(n)%model_qci(1:iv%info(radar)%levels(n)) = model_qci(1:iv%info(radar)%levels(n),n)
      iv%radar(n)%model_qsn(1:iv%info(radar)%levels(n)) = model_qsn(1:iv%info(radar)%levels(n),n)
      iv%radar(n)%model_qgr(1:iv%info(radar)%levels(n)) = model_qgr(1:iv%info(radar)%levels(n),n)

      iv%radar(n)%model_ps     = model_ps(n)

      ! [3.0] Calculate rv, rf at OBS location and initialise components of &
      ! innovation vector:

      if (fg_format == fg_format_wrf_arw_regional .or. &
          fg_format == fg_format_wrf_arw_global ) then
         call da_llxy_wrf(map_info, &
            iv%radar(n)%stn_loc%lat, iv%radar(n)%stn_loc%lon, &
            iv%radar(n)%stn_loc%x,   iv%radar(n)%stn_loc%y)
      else
         call da_llxy_default( iv%radar(n)%stn_loc%lat, iv%radar(n)%stn_loc%lon, &
            iv%radar(n)%stn_loc%x,   iv%radar(n)%stn_loc%y)
      end if

      xr = grid%xb%ds *(iv%info(radar)%x(1,n) - iv%radar(n)%stn_loc%x)
      yr = grid%xb%ds *(iv%info(radar)%y(1,n) - iv%radar(n)%stn_loc%y)

      do k=1, iv%info(radar)%levels(n)
         iv % radar(n) % rv(k) % inv = 0.0
         iv % radar(n) % rf(k) % inv = 0.0

         if (iv % radar(n) % height_qc(k) /= below_model_surface .and.  &
             iv % radar(n) % height_qc(k) /= above_model_lid) then

            if (use_radar_rv) then
               if (abs(iv % radar(n) % rv(k) % qc - missing_data) > 1) then
                  if (abs(ob % radar(n) % rv(k) - missing_r) > 1.0 .AND. &
                       iv % radar(n) % rv(k) % qc >= obs_qc_pointer) then
                  ! Rf can be used to control rv
                     zr=iv%radar(n)%height(k) - iv%radar(n)%stn_loc%elv

                     call da_radial_velocity(model_rv(k,n), model_p(k,n),  &
                        model_u(k,n), model_v(k,n), model_w(k,n),          &
                        model_qrn(k,n), model_ps(n), xr, yr, zr)

                     iv % radar(n) % rv(k) % inv = ob % radar(n) % rv(k) - model_rv(k,n)
                  end if
               end if
            end if

            if (use_radar_rf) then
               if (abs(iv % radar(n) % rf(k) % qc - missing_data) > 1) then
                  if (abs(ob % radar(n) % rf(k) - missing_r) > 1.0 .AND. &
                     iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then
                     model_rf(k,n) = leh1 + leh2 * alog10(model_rho(k,n) * model_qrn(k,n) * 1.0e+3)
                     iv % radar(n) % rf(k) % inv = ob % radar(n) % rf(k) - model_rf(k,n)
                  end if
               end if
            end if

            ! calculate retrieved hydrometeorological variables
            ! Jidong Gao JAS 2013
            if (use_radar_rhv) then
               if (abs(iv % radar(n) % rf(k) % qc - missing_data) > 1) then
                  if (abs(ob % radar(n) % rf(k) - missing_r) > 1.0 .AND. &
                     iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then

                     ! compute retrieved hydrometeoro varibles rhv 
                     if (it.eq.1) then
                        if (model_tc(k,n).ge.5.0) then
                           !JAMC 
                           !iv % radar(n) % rrno(k) = 10.0** ( (ob % radar(n) % rf(k)-leh1)/leh2 ) /model_rho(k,n)*0.001
                           !JAS
                           rze = 10.0**(ob % radar(n) % rf(k)/10.0)
                           iv % radar(n) % rrno(k) = exp ( log(rze/zrr)/1.75 )/model_rho(k,n)
                           ! rrn and rrno were assighed to missing values in read_obs_radar_ascii.inc
                           ! maximum value check, use the data under threshold 15g/kg
                           call da_radar_rf (model_qrn(k,n),model_qsn(k,n),model_qgr(k,n),model_tc(k,n),model_rho(k,n),rze)
                           if (iv % radar(n) % rrno(k) .ge. 0.015.or.ob % radar(n) % rf(k).ge.55) then
                              iv % radar(n) % rrn(k) % qc = -5  ! old -5 
                           else
                              iv % radar(n) % rrn(k) % qc = 0
                           end if
                        else if (model_tc(k,n).lt.5.0 .and. model_tc(k,n).ge.-5.0 ) then
                             czr=(model_tc(k,n)+5)/10.0
                             if (model_tc(k,n).le.0.0) then
                             czs = (1.0-czr)*zds/(zds+zg)
                             czg = (1.0-czr)*zg/(zds+zg)
                             else
                             czs = (1.0-czr)*zws/(zws+zg)
                             czg = (1.0-czr)*zg/(zws+zg)
                             end if
                             !iv % radar(n) % rrno(k) = 10.0** ( (czr*ob % radar(n) % rf(k)-leh1)/leh2 ) /model_rho(k,n)*0.001
                             rze = 10.0**(czr*ob % radar(n) % rf(k)/10.0)
                             iv % radar(n) % rrno(k) = exp ( log(rze/zrr)/1.75 )/model_rho(k,n)
                             rze = 10.0**(czs*ob % radar(n) % rf(k)/10.0)
                             iv % radar(n) % rsno(k) = exp ( log(rze/zds)/1.75 )/model_rho(k,n)
                             rze = 10.0**(czg*ob % radar(n) % rf(k)/10.0)
                             iv % radar(n) % rgro(k) = exp ( log(rze/zg )/1.75 )/model_rho(k,n)
                             iv % radar(n) % rrn(k) % qc = 0
                             iv % radar(n) % rsn(k) % qc = 0
                             iv % radar(n) % rgr(k) % qc = 0
                        else if (model_tc(k,n).le.-5.0) then
                             czs = zds/(zds+zg)
                             czg = 1.0 - czs
                             !   rra  = exp(log(rrze/zgr)/1.75)/ro1
                             rze = 10.0**(czs*ob % radar(n) % rf(k)/10.0) 
                             iv % radar(n) % rsno(k) = exp ( log(rze/zds)/1.75 )/model_rho(k,n)
                             rze = 10.0**(czg*ob % radar(n) % rf(k)/10.0)
                             iv % radar(n) % rgro(k) = exp ( log(rze/zg )/1.75 )/model_rho(k,n)
                             iv % radar(n) % rsn(k) % qc = 0
                             iv % radar(n) % rgr(k) % qc = 0
                        end if  ! temp
                     end if  ! it=1

                     if (iv % radar(n) % rrn(k) % qc >= obs_qc_pointer) then
                        ! x=b**y; y=logb(x)
                        iv % radar(n) % rrn(k) % inv = iv % radar(n) % rrno(k) - model_qrn(k,n)
                        iv % radar(n) % rrn(k) % error = iv % radar(n) % rf(k) % error * iv % radar(n) % rrno(k) * alog_10/leh2
                        ! rainwater error
                        iv % radar(n) % rrn(k) % error = amax1(0.0005,iv % radar(n) % rrn(k) % error)
                        iv % radar(n) % rrn(k) % error = amin1( 0.001,iv % radar(n) % rrn(k) % error)
                     end if
                     if (iv % radar(n) % rsn(k) % qc >= obs_qc_pointer) then
                        iv % radar(n) % rsn(k) % inv = iv % radar(n) % rsno(k) - model_qsn(k,n)
                        iv % radar(n) % rsn(k) % error = iv % radar(n) % rf(k) % error * iv % radar(n) % rsno(k) * alog_10/leh2
                        ! rainwater error
                        iv % radar(n) % rsn(k) % error = amax1(0.0005,iv % radar(n) % rrn(k) % error)
                        iv % radar(n) % rsn(k) % error = amin1( 0.001,iv % radar(n) % rrn(k) % error)
                     end if
                     if (iv % radar(n) % rgr(k) % qc >= obs_qc_pointer) then
                        iv % radar(n) % rgr(k) % inv = iv % radar(n) % rgro(k) - model_qgr(k,n)
                        iv % radar(n) % rgr(k) % error = iv % radar(n) % rf(k) % error * iv % radar(n) % rgro(k) * alog_10/leh2
                        ! rainwater error
                        iv % radar(n) % rgr(k) % error = amax1(0.0005,iv % radar(n) % rgr(k) % error)
                        iv % radar(n) % rgr(k) % error = amin1( 0.001,iv % radar(n) % rgr(k) % error)
                     end if
                  end if ! ob check
               end if ! iv qc and missing value check
            end if ! rhv

            ! retrieved water vapor
            if (use_radar_rqv) then
               !iv%%rqv and iv%%rqvo were assigned to missing values in read_obs_radar_ascii.inc
               !iter=1, rqv is missing; for second loop, dont change rqv value
               if (it .eq. 1) then
                  iv % radar(n) % rqvo(k) = 1.0*model_qs(k,n)
                  iv % radar(n) % rqv(k) % error = amax1(0.0005,0.1*iv % radar(n) % rqvo(k))
               end if

               !same condtions as use_radar_rf
               if (abs(iv % radar(n) % rf(k) % qc - missing_data) > 1) then
                  if (abs(ob % radar(n) % rf(k) - missing_r) > 1.0 .AND. &
                       iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then
                     zr=iv%radar(n)%height(k) - iv%radar(n)%stn_loc%elv
                     if ( ob % radar(n) % rf(k) .lt. 25.0  .or. zr.lt.1500)then 
                        iv % radar(n) % rqv(k) % qc = -5
                     else
                        iv % radar(n) % rqv(k) % qc = 0
                     end if
                     if (iv % radar(n) % rqv(k) % qc >= obs_qc_pointer) then                        
                        if (ob % radar(n) % rf(k) >= 25.0 .and. ob % radar(n) % rf(k) < 40.0) then
                           iv % radar(n) % rqvo(k)= 0.85*model_qs(k,n)
                        else if (ob % radar(n) % rf(k) >= 40.0 .and. ob % radar(n) % rf(k) < 50.0) then
                           iv % radar(n) % rqvo(k)= 0.95*model_qs(k,n)
                        end if
                        iv % radar(n) % rqv(k) % inv = iv % radar(n) % rqvo(k) - model_qv(k,n)
                        ! iv % radar(n) % rqv(k) % inv must be >= 0.0
                        iv % radar(n) % rqv(k) % inv = amax1(0.0,iv % radar(n) % rqv(k) % inv )
                        iv % radar(n) % rqv(k) % error = amax1(0.001,0.20*iv % radar(n) % rqvo(k))
                     end if
                  end if ! obs qc and mising
               end if  ! iv qc and missing
            end if  ! use rqv

         end if  ! not surface or model lid
      end do
   end do

   !------------------------------------------------------------------------
   ! [4.0] Perform optional maximum error check:  
   !------------------------------------------------------------------------

   if (check_max_iv)  then
      call da_check_max_iv_radar(iv, it, irv, irf, irvf, irff)
   end if

   if (rootproc .and. check_max_iv_print) then
      write(unit = check_max_iv_unit, fmt ='(/,A,i5,A)')   &
         'For outer iteration ', it, ', Total Rejections for radar follows:'

      if (use_radar_rv) then
          write( unit = check_max_iv_unit, fmt = '(/,2(A,I6))') &
            'Number of failed rv observations:     ',irvf, ' on ',irv
      end if

      if (use_radar_rf) then
         write( unit = check_max_iv_unit, fmt = '(/,2(A,I6))') &
            'Number of failed rf observations:     ',irff, ' on ',irf
      end if
   end if

   deallocate (model_p)
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_w)

   deallocate (model_rv)
   deallocate (model_rf)
   deallocate (model_ps)

   deallocate (model_qv)
   deallocate (model_qs)
   deallocate (model_tc)

   deallocate (model_qrn)
   deallocate (model_rho)
   deallocate (model_qcl) 
   deallocate (model_qci)
   deallocate (model_qsn)
   deallocate (model_qgr)
  
   if (trace_use) call da_trace_exit("da_get_innov_vector_radar")

end subroutine da_get_innov_vector_radar


subroutine da_radial_velocity(rv,p,u,v,w,qrn,ps,x,y,z)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)  :: x, y, z
   real, intent(in)  :: p, u, v, w, qrn, ps
   real, intent(out) :: rv

   real :: r, alpha, vt
   real    :: qrrc

   qrrc = 1.0e-3
   vt = 0.0

   if (trace_use) call da_trace_entry("da_radial_velocity")

   r=sqrt(x*x+y*y+z*z)
   alpha=(ps/p)**0.4

!   if (qrn <= 0.0) vt=0.0
!   if (qrn >  0.0) vt=5.4*alpha*qrn**0.125

   if (use_radar_rf .or. use_radar_rhv)then
      if (qrn <= qrrc)then
         vt=0.0
      else
         vt=5.4*alpha*qrn**0.125
      end if
   end if
   rv=u*x+v*y+(w-vt)*z
   rv=rv/r

   if (trace_use) call da_trace_exit("da_radial_velocity")

end subroutine da_radial_velocity


subroutine da_radial_velocity_lin(rv,p,u,v,w,qrn,ps,x,y,z,qrn9)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)  :: x, y, z
   real, intent(in)  :: p, u, v, w, qrn, ps
   real, intent(in)  :: qrn9
   real, intent(out) :: rv

   real    :: r, alpha, vt
   real    :: qrrc

   qrrc = 1.0e-3
   vt = 0.0

   if (trace_use) call da_trace_entry("da_radial_velocity_lin")

   r     = sqrt(x*x+y*y+z*z)
   alpha = (ps/p)**0.4


   if (use_radar_rf .or. use_radar_rhv)then
      if (qrn9 <= qrrc)then
         vt=0.0
      else
         vt=0.675*alpha*qrn9**(-0.875)*qrn
      end if
   end if

!   if (qrn9 <= 0.0) then
!      vt=0.0
!   end if

!   if (qrn9 >  0.0) then
!      vt=0.675*alpha*qrn9**(-0.875)*qrn
!   end if

   rv = u*x+v*y+(w-vt)*z
   rv = rv/r

   if (trace_use) call da_trace_exit("da_radial_velocity_lin")

end subroutine da_radial_velocity_lin


subroutine da_radial_velocity_adj(rv,p,u,v,w,qrn,ps,x,y,z,qrn9)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: x, y, z
   real, intent(in)    :: p
   real, intent(in)    :: qrn9
   real, intent(in)    :: ps
   real, intent(inout) :: rv
   real, intent(inout) :: u, v, w, qrn

   real :: r, alpha, vt
   real :: qrrc

   qrrc = 1.0e-3

   if (trace_use) call da_trace_entry("da_radial_velocity_adj")

   r     = sqrt(x*x+y*y+z*z)
   alpha = (ps/p)**0.4

   rv = rv/r
   u  = u + rv*x
   v  = v + rv*y
   w  = w + rv*z
   vt = -rv*z

  if (use_radar_rf .or. use_radar_rhv)then
!   if (qrn9 >  0.0) then
!      qrn = qrn + vt*0.675*alpha*qrn9**(-0.875)
!   end if
   if (qrn9 >  qrrc) then
      qrn = qrn + vt*0.675*alpha*qrn9**(-0.875)
   end if
  end if
   if (trace_use) call da_trace_exit("da_radial_velocity_adj")

end subroutine da_radial_velocity_adj


subroutine da_calculate_grady_radar(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer  :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_radar")

   do n=1, iv%info(radar)%nlocal
      do k=1, iv%info(radar)%levels(n)
         if (use_radar_rv) then
            if (iv%radar(n)%rv(k)%qc < obs_qc_pointer) then
               re%radar(n)%rv(k) = 0.0
            end if
            jo_grad_y%radar(n)%rv(k) = -re%radar(n)%rv(k) / (iv%radar(n)%rv(k)%error * iv%radar(n)%rv(k)%error) 
         end if

         if (use_radar_rf) then
            if (iv%radar(n)%rf(k)%qc < obs_qc_pointer) then
               re%radar(n)%rf(k) = 0.0
            end if
            jo_grad_y%radar(n)%rf(k) = -re%radar(n)%rf(k) / (iv%radar(n)%rf(k)%error * iv%radar(n)%rf(k)%error) 
         end if

         if (.not.use_radar_rf .and. use_radar_rhv) then
            if (iv%radar(n)%rrn(k)%qc < obs_qc_pointer) then
               re%radar(n)%rrn(k) = 0.0
            end if
            jo_grad_y%radar(n)%rrn(k) = -re%radar(n)%rrn(k) / (iv%radar(n)%rrn(k)%error * iv%radar(n)%rrn(k)%error)
            if (iv%radar(n)%rsn(k)%qc < obs_qc_pointer) then
               re%radar(n)%rsn(k) = 0.0
            end if
            jo_grad_y%radar(n)%rsn(k) = -re%radar(n)%rsn(k) / (iv%radar(n)%rsn(k)%error * iv%radar(n)%rsn(k)%error)
            if (iv%radar(n)%rgr(k)%qc < obs_qc_pointer) then
               re%radar(n)%rgr(k) = 0.0
            end if
            jo_grad_y%radar(n)%rgr(k) = -re%radar(n)%rgr(k) / (iv%radar(n)%rgr(k)%error * iv%radar(n)%rgr(k)%error)
         end if

         if (use_radar_rqv) then
            if (iv%radar(n)%rqv(k)%qc < obs_qc_pointer) then
               re%radar(n)%rqv(k) = 0.0
            end if
            jo_grad_y%radar(n)%rqv(k) = -re%radar(n)%rqv(k) / (iv%radar(n)%rqv(k)%error * iv%radar(n)%rqv(k)%error)
         end if

      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_radar")

end subroutine da_calculate_grady_radar


subroutine da_max_error_qc_radar (it, info, n,field, max_error, ix, ixf, iunit, var, failed, print_details)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer          ,   intent (in)   :: it
   type (infa_type) ,   intent(in)    :: info
   integer,             intent(in)    :: n
   type (field_type),   intent(inout) :: field
   real             ,   intent(in)    :: max_error
   integer          ,   intent(inout) :: ix, ixf
   integer          ,   intent(in)    :: iunit
   character (len=*),   intent(in)    :: var
   logical          ,   intent(out)   :: failed
   logical          ,   intent(in)    :: print_details

   real                               :: err, err_max
   integer                            :: qc_flag

   if (trace_use) call da_trace_entry("da_max_error_qc_radar")

   failed = .false.

   qc_flag = field % qc
   err_max = field % error * max_error
   err     = field % inv
   err     = ABS (err)

   ix     = ix + 1
   if (it > 1 .and. qc_flag == fails_error_max) field%qc = 0
   if (err > err_max) then
      if (field % qc > fails_error_max) field % qc = fails_error_max 

      ixf = ixf + 1
      failed = .true.

      if (print_details .and. failed) then
         if (err_max .LE. 0.0) then
            write (iunit , fmt = '(A,3(F12.1,1X),A,A,A,A,A,3f10.2)')   &
               "Err_max < 0 ==> ",err,err_max,max_error, " for ", var, &
               " OBS ID: ", info%platform(n),     &
               " LA/LON/ELV:", info%lat(1,n), info%lon(1,n), info%elv(n)
            ! call da_error("da_max_error_qc_radar.inc",46,(/"Erk"/))
         end if

         write (iunit , fmt = '(A,A,A,I5,A,I5,A,F4.1,A,A,A,2F12.1)') &
            "Err_max failed:ID=", info%platform(n),&
            "  Ix=", ix, "  Ixf=", ixf, " Err_max ratio =",err/err_max, &
            " for ", var, " inv, error:",field % inv, field % error
      end if
      field % inv = 0.0
   end if

   if (trace_use) call da_trace_exit("da_max_error_qc_radar")

end subroutine da_max_error_qc_radar


subroutine da_write_oa_radar_ascii ( ob, iv, re, it )

   !---------------------------------------------------------------------------
   ! Purpose: write out OMB and OMA vector structure for radar data.
   !---------------------------------------------------------------------------

   implicit none

   type (y_type),     intent(in)  :: ob       ! Observation structure.
   type (iv_type),    intent(in)  :: iv       ! O-B structure.
   type (y_type),     intent(in)  :: re       ! O-A structure.
   integer,           intent(in)  :: it       ! external iteration counter

   integer                        :: n , num_obs       ! Loop counter.
   integer                        :: i, k     ! Index dimension.
   integer                        :: nlevelss ! Number of obs levels.

   integer            :: ios, oma_radar_unit, omb_radar_unit, omb_radar_iter_unit 
   character(len=filename_len)  :: filename , file 
   integer            :: ndomain
   character(len=filename_len), allocatable     ::filename1(:) 
   INTEGER :: m, kk, l 
   REAL :: rv_obs, rv_inv, rv_error , rf_obs, rf_inv, rf_error, rf_inc, rv_inc 
   INTEGER :: rv_qc, rf_qc, levels, num, error 
   REAL :: lat,lon,press, zk 
   CHARACTER(len=5) :: stn_id 

    if (trace_use) call da_trace_entry("da_write_oa_radar_ascii")

    write(unit=message(1),fmt='(A)') 'Writing radar OMA ascii file'
    call da_message(message(1:1))

   write(unit=filename, fmt='(a,i2.2,a,i4.4)') 'radar_omb_oma_',it,'.', myproc
   call da_get_unit(oma_radar_unit)
   open(unit=oma_radar_unit,file=trim(filename),form='formatted',iostat=ios)
   if (ios /= 0) Then
       call da_error("da_write_oa_radar_ascii.inc",45, &
         (/"Cannot open oma radar file"//filename/))
   endif

   if (iv % info(radar)%nlocal  >0 ) then
      num_obs = 0
      do n = 1, iv% info(radar)%nlocal
         if (iv%info(radar)%proc_domain(1,n)) num_obs=num_obs+1    
      end do
      if (num_obs > 0) then
         write(oma_radar_unit,'(a20,i8)')'radar', num_obs
         num_obs = 0
        do n = 1, iv % info(radar)%nlocal
            if (iv%info(radar)%proc_domain(1,n)) then  
               num_obs = num_obs + 1
               write(oma_radar_unit,'(i8)') iv % info(radar) % levels(n)
               do k = 1, iv % info(radar) % levels(n)
                 write(oma_radar_unit,'(2i8,a5,2f9.2,f17.7,2(2f17.7,i8,2f17.7),f17.7)')&
                    num_obs , k, 'RADAR', &  
                    iv % info (radar)% lat(1,n), &       ! Latitude
                    iv % info (radar) % lon(1,n), &       ! Longitude
                    iv % radar(n) % height(k), &           ! Obs height in m
                    ob%radar(n)%rv(k),&
                    iv%radar(n)%rv(k)%inv,iv%radar(n)%rv(k)%qc,iv%radar(n)%rv(k)%error,&
                    re%radar(n)%rv(k), &! O, O-B, O-A rv
                    ob%radar(n)%rf(k), &
                    iv%radar(n)%rf(k)%inv,iv%radar(n)%rf(k)%qc,iv%radar(n)%rf(k)%error,&
                    re%radar(n)%rf(k),iv%info(radar)%zk(k,n)   ! O, O-B, O-A rf
               end do
            end if
         end do
      end if
   end if


   ! Wait to ensure all temporary files have been written
   call mpi_barrier(comm, ierr)


  close (oma_radar_unit)
  call da_free_unit(oma_radar_unit)

  
  IF (rootproc) THEN
  call da_get_unit(omb_radar_unit)
  allocate (filename1(0:num_procs-1)) 
      do k = 0,num_procs-1
         write(unit=filename1(k),fmt ='(a,i2.2,a,i4.4)')'radar_omb_oma_',it,'.',k
      end do
   call da_get_unit(omb_radar_iter_unit)
   write(unit=file,fmt ='(a,i2.2)')'radar_omb_oma_',it
   open(omb_radar_iter_unit,file=trim(file),form='formatted', status='replace', iostat=ios)
     if (ios /= 0) call da_error("da_write_oa_radar_ascii.inc",99, &
         (/"Cannot open file "//file/))
  ENDIF
 
  num_obs = 0
  IF (iv % info(radar)%nlocal  >0 ) then
      do n = 1, iv% info(radar)%nlocal
        if (iv%info(radar)%proc_domain(1,n)) num_obs=num_obs+1
      end do
  ENDIF
   call da_proc_sum_int(num_obs)    
    IF (num_obs > 0 .and. rootproc) then
      write(omb_radar_iter_unit,'(a20,i8)')'radar', num_obs  
      num_obs = 0
      num = 0
      do k = 0,num_procs-1

        open(omb_radar_unit,file=trim(filename1(k)),status='old',iostat=error)
          if (error /= 0) call da_error("da_write_oa_radar_ascii.inc",117, &
         (/"Cannot open file "//file/))

      read(omb_radar_unit, '(20x,i8)', iostat=error)num_obs
      IF(error /= 0)THEN
         write(unit=message(1),fmt='(A,A)') 'Nothing to read from ',filename1(k)
         call da_message(message(1:1))
         cycle
      ENDIF   
         if (num_obs > 0) then
           do  n = 1, num_obs    
               read(omb_radar_unit,'(i8)') levels
                  write(omb_radar_iter_unit,'(i8)')levels
                  num = num + 1 
               do m = 1, levels
                  read(omb_radar_unit,'(2i8,a5,2f9.2,f17.7,2(2f17.7,i8,2f17.7),f17.7)')&
                      kk,l, stn_id, &          ! Station
                      lat, lon, press, &       ! Lat/lon, dummy    
                      rv_obs, rv_inv, rv_qc, rv_error, rv_inc, & 
                      rf_obs, rf_inv, rf_qc, rf_error, rf_inc, zk
                     write(omb_radar_iter_unit,'(2i8,a5,2f9.2,f17.7,2(2f17.7,i8,2f17.7),f17.7)')&
                        num,m,stn_id, &          ! Station
                         lat, lon, press, &       ! Lat/lon, dummy    
                     rv_obs, rv_inv, rv_qc, rv_error, rv_inc, & 
                     rf_obs, rf_inv, rf_qc, rf_error, rf_inc, zk
               end do
            enddo
         endif                 
      end do
    ENDIF




   if (rootproc) then
      close(omb_radar_unit)
      close(omb_radar_iter_unit)
      call da_free_unit(omb_radar_unit)
      call da_free_unit(omb_radar_iter_unit)
      deallocate (filename1)
   end if


   if (trace_use) call da_trace_exit("da_write_oa_radar_ascii")

end subroutine da_write_oa_radar_ascii





subroutine da_radar_rf (ra,rs,rg,tc1,ro1,rze)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: ra, rs, rg, tc1, ro1
   real, intent(inout) :: rze

   real                :: zrr, zds, zws, zgr
   real                :: zerr, zews, zeds, zegr
   real                :: czr

   if (trace_use) call da_trace_entry("da_radar_rf")

   ! ro1 use wrfda value, do not need to scale by 0.001

   zrr = 3.63*1.00e+9  ! rainwater
   zws = 4.26*1.00e+11 ! wet snow
   zds = 9.80*1.00e+8  ! dry snow
   zgr = 4.33*1.00e+10 ! graupel

   zerr = zrr*(ro1*ra)**1.75
   zews = zws*(ro1*rs)**1.75
   zeds = zds*(ro1*rs)**1.75
   zegr = zgr*(ro1*rg)**1.75

   if (tc1.ge.5.0) then
      rze = zerr
   elseif (tc1.le.5.0 .and. tc1 .ge.-5.0) then
      czr = (tc1+5.0)/10.0
      if (tc1.le.0.0) then
         rze = czr*zerr + (1.0-czr)*(zeds+zegr)
      else
         rze = czr*zerr + (1.0-czr)*(zews+zegr)
      end if
   elseif (tc1.lt.-5.0) then
       rze = zeds + zegr
   end if

   if (rze.lt.1.0e-20) rze=1.0e-20

   if (trace_use) call da_trace_exit("da_radar_rf")

end subroutine da_radar_rf


end module da_radar

