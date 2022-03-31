












module da_gpspw

   use module_dm, only : wrf_dm_sum_real
   use module_domain, only : domain

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print,kts,kte, &
      missing, max_error_uv, max_error_t, rootproc, gpspw, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv,  &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, gpspw,max_error_thickness, &
      pseudo_var, num_pseudo, use_gpspwobs, use_gpsztdobs, max_error_pw,fails_error_max, &
      fails_error_max,pseudo_err,pseudo_x, pseudo_y, stdout, &
      pseudo_z,pseudo_val,max_error_ref, trace_use_dull, pseudo, its,ite,jts,jte,&
      ob_vars,qcstat_conv_unit
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      maxmin_type, da_allocate_observations
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_reporting, only : da_error, da_message, message
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_gpspw1_type
      real          :: tpw                      
   end type residual_gpspw1_type

   type maxmin_gpspw_stats_type
      type (maxmin_type)         :: tpw
   end type maxmin_gpspw_stats_type

   type stats_gpspw_type
      type (maxmin_gpspw_stats_type)  :: maximum, minimum
      type (residual_gpspw1_type)     :: average, rms_err
   end type stats_gpspw_type

contains

subroutine da_ao_stats_gpspw (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O

   type (stats_gpspw_type) :: stats
   integer                 :: ntpw
   integer                 :: n
   real                    :: o_minus_b, o_minus_a, sigma_o, sigma_b
   real                    :: o_minus_b_0, o_minus_a_0, sigma_o_0, sigma_b_0 

   if (trace_use_dull) call da_trace_entry("da_ao_stats_gpspw")

   ntpw = 0
   o_minus_b_0 = 0.0; o_minus_a_0 = 0.0; sigma_o_0 = 0.0; sigma_b_0 = 0.0

   stats%maximum%tpw = maxmin_type (missing_r, 0, 0)
   stats%minimum%tpw = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_gpspw1_type(0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(gpspw)%nlocal
      if (iv%info(gpspw)%proc_domain(1,n)) then
         call da_stats_calculate (n, 0, iv%gpspw(n)%tpw%qc, &
            re%gpspw(n)%tpw, ntpw,   &
            stats%minimum%tpw,  stats%maximum%tpw, &
            stats%average%tpw,  stats%rms_err%tpw)

         if ((pseudo_var(1:3) == 'tpw' .or.pseudo_var(1:3) == 'ztd' ) &
                        .and. num_pseudo > 0) then
            o_minus_b = iv%GPSpw(n)%tpw%inv
            o_minus_a = re%gpspw(n)%tpw
            sigma_o   = iv%gpspw(n)%tpw%error

            ! Calculate equivalent sigma_b using
            ! 
            ! O-A=(O-B)*sigma_o/(sigma_o+sigma_b)

            sigma_b = sqrt ((o_minus_b - o_minus_a) / o_minus_a) * sigma_o
            o_minus_b_0 = wrf_dm_sum_real (o_minus_b)
            o_minus_a_0 = wrf_dm_sum_real (o_minus_a)
            sigma_o_0 = wrf_dm_sum_real (sigma_o)
            sigma_b_0 = wrf_dm_sum_real (sigma_b)
!           write(unit=stdout,fmt='(A,F10.2)') &
!              'TEST_COVERAGE_da_ao_stats_gpspw:  o_minus_b_0 = ', o_minus_b_0
!           write(unit=stdout,fmt='(A,F10.2)') &
!              'TEST_COVERAGE_da_ao_stats_gpspw:  o_minus_a_0 = ', o_minus_a_0
!           write(unit=stdout,fmt='(A,F10.2)') &
!              'TEST_COVERAGE_da_ao_stats_gpspw:  sigma_o_0 = ', sigma_o_0
!           write(unit=stdout,fmt='(A,F10.2)') &
!              'TEST_COVERAGE_da_ao_stats_gpspw:  sigma_b_0 = ', sigma_b_0
            if (rootproc) then  
               write(stats_unit,'(/A,A3,A,f12.3)')  &  
                 ' Pseudo ', pseudo_var, ' O-B: ', o_minus_b_0  
               write(stats_unit,' (A,A3,A,f12.3)')  &  
                 ' Pseudo ', pseudo_var, ' O-A: ', o_minus_a_0  
               write(stats_unit,' (A,A3,A,f12.3)')  &  
                 ' Pseudo ', pseudo_var, ' sigma_o: ', sigma_o_0  
               write(stats_unit,'(A,A3,A,f12.3)')  &  
                 ' Pseudo ', pseudo_var, ' sigma_b: ', sigma_b_0 
            end if  
         end if
      end if    ! end if (iv%info(gpspw)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (ntpw)
   iv%nstats(gpspw) = ntpw

   call da_proc_stats_combine(stats%average%tpw, stats%rms_err%tpw, &
      stats%minimum%tpw%value, stats%maximum%tpw%value, &
      stats%minimum%tpw%n, stats%maximum%tpw%n, &
      stats%minimum%tpw%l, stats%maximum%tpw%l)

   if (rootproc) then
      if (ntpw /= 0) then
         if (use_gpspwObs) then
           write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for gpspw'
         else if (use_gpsztdObs) then
           write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for gpsztd'
         endif

         call da_print_stats_gpspw(stats_unit, ntpw, stats)
       end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_gpspw")

end subroutine da_ao_stats_gpspw


subroutine da_jo_and_grady_gpspw(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout):: jo          ! Obs cost function.

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_gpspw")

   jo % gpspw_tpw = 0.0

   if (iv%info(gpspw)%nlocal > 0) then
      do n=1, iv%info(gpspw)%nlocal
         jo_grad_y%gpspw(n)%tpw = -re%gpspw(n)%tpw / (iv%gpspw(n)%tpw%error * iv%gpspw(n)%tpw%error)

         if (iv%info(gpspw)%proc_domain(1,n)) then
            jo % gpspw_tpw = jo % gpspw_tpw - re%gpspw(n)%tpw * jo_grad_y%gpspw(n)%tpw
         end if
      end do

      jo % gpspw_tpw = 0.5 * jo % gpspw_tpw
   end if

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_gpspw")

end subroutine da_jo_and_grady_gpspw


subroutine da_residual_gpspw(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

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

   type (bad_data_type)              :: n_obs_bad
   integer                           :: n

   if (trace_use_dull) call da_trace_entry("da_residual_gpspw")

   n_obs_bad % q % num = number_type(0, 0, 0)

   do n=1, iv%info(gpspw)%nlocal
      np_available = np_available + 1
      re%gpspw(n)%tpw = da_residual(n, 0, y%gpspw(n)%tpw, iv%gpspw(n)%tpw, n_obs_bad % q)
   end do

   np_missing  = np_missing  + n_obs_bad % q % num % miss 
   np_bad_data = np_bad_data + n_obs_bad % q % num % bad
   np_obs_used = np_obs_used + n_obs_bad % q % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_gpspw")

end subroutine da_residual_gpspw



subroutine da_oi_stats_gpspw (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_gpspw_type) :: stats
   integer                 :: ntpw
   integer                 :: n

   if (trace_use_dull) call da_trace_entry("da_oi_stats_gpspw")

   ntpw = 0

   stats%maximum%tpw = maxmin_type(missing_r, 0, 0)
   stats%minimum%tpw = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_gpspw1_type(0.0)
   stats%rms_err = stats%average
   do n=1, iv%info(gpspw)%nlocal
      if (iv%info(gpspw)%proc_domain(1,n)) then
         call da_stats_calculate(iv%info(gpspw)%obs_global_index(n), &
            0, iv%gpspw(n)%tpw%qc, &
            iv%gpspw(n)%tpw%inv, ntpw, &
            stats%minimum%tpw  , stats%maximum%tpw, &
            stats%average%tpw  , stats%rms_err%tpw)

      end if    ! end if (iv%info(gpspw)%proc_domain(1,n))
   end do

   ! do inter-processor communication to gather statistics.

   call da_proc_sum_int(ntpw)

   call da_proc_stats_combine(stats%average%tpw, stats%rms_err%tpw, &
      stats%minimum%tpw%value, stats%maximum%tpw%value, &
      stats%minimum%tpw%n, stats%maximum%tpw%n, &
      stats%minimum%tpw%l, stats%maximum%tpw%l)

   if (rootproc) then
      if (ntpw /= 0) then
        if (use_gpspwObs) then
          write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for gpspw'
        else if (use_gpsztdObs) then
          write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for gpsztd'
        endif

         call da_print_stats_gpspw(stats_unit, ntpw, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_gpspw")

end subroutine da_oi_stats_gpspw


subroutine da_print_stats_gpspw(stats_unit, ntpw, Gpspw)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: ntpw
   type (stats_gpspw_type), intent(in)    :: gpspw

   if (trace_use_dull) call da_trace_entry("da_residual_gpspw")

   if (ntpw > 0) then
      if (use_gpspwObs) then
        write(unit=stats_unit, fmt='(a/)') '   var           tpw(cm)     n'
      else if (use_gpsztdObs) then
        write(unit=stats_unit, fmt='(a/)') '   var           ztd(cm)     n'
      endif
 
      write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntpw

      write(unit=stats_unit, fmt='(a, f12.4,i5)') &
        ' Minimum(n): ', gpspw%minimum%tpw%value, &
                         gpspw%minimum%tpw%n    , &
        ' Maximum(n): ', gpspw%maximum%tpw%value, &
                         gpspw%maximum%tpw%n
      write(unit=stats_unit, fmt='(a, f12.4,5x)') &
        ' Average   : ', gpspw%average%tpw/real(ntpw), &
        '    RMSE   : ', sqrt(gpspw%rms_err%tpw/real(ntpw))
   end if

   if (trace_use_dull) call da_trace_exit("da_residual_gpspw")

end subroutine da_print_stats_gpspw


subroutine da_transform_xtoy_gpspw (grid, iv, y)

   !----------------------------------------------------------------
   ! Purpose: Calculate the difference in the elevation between model surface and GPS site
   !----------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid 
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer :: n        ! Loop counter.
   integer :: i, j     ! Index dimension.
   real    :: dx, dxm  ! 
   real    :: dy, dym  !

   integer :: k          ! Index dimension.
   real    :: dpw, ddpw  ! adjustment pw [mm]
   real    :: obs_terr   ! real terrain height at GPS site [m]
   real    :: model_terr ! model terrain height at GPS site[m]
   real    :: model_q(kts:kte)    ! model q at GPS site [kg/kg]
   real    :: model_z(kts:kte)    ! model z at GPS site [m]
   real    :: model_rho(kts:kte)  ! model rho at GPS site [kg/m^3]
   real    :: model_dq(kts:kte)   ! model dq at GPS site [kg/kg]
   real    :: model_drho(kts:kte) ! model drho at GPS site [kg/m^3]

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpspw")

   do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2
      i   = iv%info(gpspw)%i(1,n)
      dy  = iv%info(gpspw)%dy(1,n)
      dym = iv%info(gpspw)%dym(1,n)
      j   = iv%info(gpspw)%j(1,n)
      dx  = iv%info(gpspw)%dx(1,n)
      dxm = iv%info(gpspw)%dxm(1,n)

      ! Mathematical transformation:

      ! kusaka 2004-04-08
      ! y % gpspw(n)% tpw = dym* (dxm * grid%xa%tpw(i,j) + &
      !    dx  * grid%xa%tpw(i+1,j)) + &
      !    dy * (dxm * grid%xa%tpw(i,j+1) + &
      !    dx  * grid%xa%tpw(i+1,j+1))
      dpw = 0.0
      obs_terr   = iv%info(gpspw)%elv(n)
      model_terr = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + &
                   dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))

      if (obs_terr <= model_terr) then 
         model_q(1) = dym*(dxm*grid%xb%q(i,j,1)   + dx*grid%xb%q(i+1,j,1)) + &
             dy *(dxm*grid%xb%q(i,j+1,1) + dx*grid%xb%q(i+1,j+1,1))
         model_rho(1) = dym*(dxm*grid%xb%rho(i,j,1)   + dx*grid%xb%rho(i+1,j,1)) + &
            dy *(dxm*grid%xb%rho(i,j+1,1) + dx*grid%xb%rho(i+1,j+1,1))

         model_dq(1) = dym*(dxm*grid%xa%q(i,j,1)   + dx*grid%xa%q(i+1,j,1)) + &
            dy *(dxm*grid%xa%q(i,j+1,1) + dx*grid%xa%q(i+1,j+1,1))
         model_drho(1) = dym*(dxm*grid%xa%rho(i,j,1)   + dx*grid%xa%rho(i+1,j,1)) + &
            dy *(dxm*grid%xa%rho(i,j+1,1) + dx*grid%xa%rho(i+1,j+1,1))

         dpw = (model_rho(1)*model_dq(1) + model_drho(1)*model_q(1)) &
            * (obs_terr - model_terr)
      else 
         model_z(1) = dym*(dxm*grid%xb%hf(i,j,1)   + dx*grid%xb%hf(i+1,j,1)) + &
            dy *(dxm*grid%xb%hf(i,j+1,1) + dx*grid%xb%hf(i+1,j+1,1))
         do k = kts, kte
            if (model_z(k) >= obs_terr) exit

            model_z(k+1) = dym*(dxm*grid%xb%hf(i,j,k+1)   + dx*grid%xb%hf(i+1,j,k+1)) + &
               dy *(dxm*grid%xb%hf(i,j+1,k+1) + dx*grid%xb%hf(i+1,j+1,k+1))

            model_q(k) = dym*(dxm*grid%xb%q(i,j,k)   + dx*grid%xb%q(i+1,j,k)) + &
               dy *(dxm*grid%xb%q(i,j+1,k) + dx*grid%xb%q(i+1,j+1,k))
            model_rho(k) = dym*(dxm*grid%xb%rho(i,j,k)   + dx*grid%xb%rho(i+1,j,k)) + &
               dy *(dxm*grid%xb%rho(i,j+1,k) + dx*grid%xb%rho(i+1,j+1,k))

            model_dq(k) = dym*(dxm*grid%xa%q(i,j,k)   + dx*grid%xa%q(i+1,j,k)) + &
               dy *(dxm*grid%xa%q(i,j+1,k) + dx*grid%xa%q(i+1,j+1,k))
            model_drho(k) = dym*(dxm*grid%xa%rho(i,j,k)   + dx*grid%xa%rho(i+1,j,k)) + &
               dy *(dxm*grid%xa%rho(i,j+1,k) + dx*grid%xa%rho(i+1,j+1,k))

            ! Here assumed that "model_z" is constant, i.e. grid%xa%hf=0.0. With MM5, 
            ! this is true, but with WRF, grid%xa%hf may not be zero (?). In the WRF 
            ! model space (x,y,znu), only the "znu" is constant, but all variables 
            ! including height could be changed at the "znu" levels. So here is only 
            ! an approximation for WRF. The following few lines of code is written
            ! by Y.-R. Guo 07/16/2004.

            if (model_z(k+1) <= obs_terr) then
               ddpw = (model_rho(k)*model_dq(k) + model_drho(k)*model_q(k)) * (model_z(k+1)-model_z(k))
            else 
               ddpw = (model_rho(k)*model_dq(k) + model_drho(k)*model_q(k)) * (obs_terr-model_z(k))
            end if

            dpw = dpw + ddpw
         end do 
      end if 
      y % gpspw(n)% tpw = dym* (dxm * grid%xa%tpw(i,j) + &
                                dx  * grid%xa%tpw(i+1,j)) + &
                         dy * (dxm * grid%xa%tpw(i,j+1) + &
                                dx  * grid%xa%tpw(i+1,j+1)) &
                         + 0.1*dpw
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpspw")

end subroutine da_transform_xtoy_gpspw


subroutine da_transform_xtoy_gpspw_adj(grid, iv, jo_grad_y, jo_grad_x)

   !----------------------------------------------------------------
   ! Purpose: Calculate the difference in the elevation between model surface and GPS site
   !----------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.
   integer :: i, j     ! Index dimension.
   real    :: dx, dxm  !
   real    :: dy, dym  !

   integer :: k        ! Index dimension.
   real    :: dpw,ddpw     ! adjustment pw [mm]
   real    :: obs_terr     ! real terrain height at GPS site [m]
   real    :: model_terr   ! model terrain height at GPS site[m]
   real    :: model_q(kts:kte)   ! model q at GPS site [kg/kg]
   real    :: model_z(kts:kte)   ! model z at GPS site [m]
   real    :: model_rho(kts:kte) ! model rho at GPS site [kg/m^3]
   real    :: model_dq(kts:kte)   ! model dq at GPS site [kg/kg]
   real    :: model_drho(kts:kte) ! model drho at GPS site [kg/m^3]

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpspw_adj")

   do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2
      i   = iv%info(gpspw)%i(1,n)
      dy  = iv%info(gpspw)%dy(1,n)
      dym = iv%info(gpspw)%dym(1,n)
      j   = iv%info(gpspw)%j(1,n)
      dx  = iv%info(gpspw)%dx(1,n)
      dxm = iv%info(gpspw)%dxm(1,n)

      !  Initialise the varibles
      dpw           = 0.0
      model_q(:)    = 0.0
      model_dq(:)   = 0.0
      model_rho(:)  = 0.0
      model_drho(:) = 0.0

      obs_terr   = iv%info(gpspw)%elv(n)
      model_terr = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + &
         dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))

      dpw = 0.1 *jo_grad_y%gpspw(n)%tpw

      jo_grad_x%tpw(i  ,j )  = jo_grad_x%tpw(i  ,j )  + dxm*dym*jo_grad_y%gpspw(n)%tpw
      jo_grad_x%tpw(i+1,j )  = jo_grad_x%tpw(i+1,j )  + dym*dx *jo_grad_y%gpspw(n)%tpw
      jo_grad_x%tpw(i  ,j+1) = jo_grad_x%tpw(i,j+1)   + dxm *dy*jo_grad_y%gpspw(n)%tpw
      jo_grad_x%tpw(i+1,j+1) = jo_grad_x%tpw(i+1,j+1) + dx *dy *jo_grad_y%gpspw(n)%tpw

      if (obs_terr <= model_terr) then 
         model_q(1)   = dym*(dxm*grid%xb%q(i,j,1)   + dx*grid%xb%q(i+1,j,1)) + &
            dy *(dxm*grid%xb%q(i,j+1,1) + dx*grid%xb%q(i+1,j+1,1))
         model_rho(1) = dym*(dxm*grid%xb%rho(i,j,1)   + dx*grid%xb%rho(i+1,j,1)) + &
            dy *(dxm*grid%xb%rho(i,j+1,1) + dx*grid%xb%rho(i+1,j+1,1))

         model_dq(1)   =  model_rho(1)*(obs_terr - model_terr)*dpw
         model_drho(1) =  model_q(1)  *(obs_terr - model_terr)*dpw

         jo_grad_x%q(i,j,1)       = jo_grad_x%q(i,j,1)       + dym*dxm*model_dq(1)
         jo_grad_x%q(i+1,j,1)     = jo_grad_x%q(i+1,j,1)     + dym*dx*model_dq(1)
         jo_grad_x%q(i,j+1,1)     = jo_grad_x%q(i,j+1,1)     + dy*dxm*model_dq(1)
         jo_grad_x%q(i+1,j+1,1)   = jo_grad_x%q(i+1,j+1,1)   + dy*dx*model_dq(1)

         jo_grad_x%rho(i,j,1)     = jo_grad_x%rho(i,j,1)     + dym*dxm*model_drho(1)
         jo_grad_x%rho(i+1,j,1)   = jo_grad_x%rho(i+1,j,1)   + dym*dx*model_drho(1)
         jo_grad_x%rho(i,j+1,1)   = jo_grad_x%rho(i,j+1,1)   + dy*dxm*model_drho(1)
         jo_grad_x%rho(i+1,j+1,1) = jo_grad_x%rho(i+1,j+1,1) + dy*dx*model_drho(1)
      else 
         model_z(1) = dym*(dxm*grid%xb%hf(i,j,1)   + dx*grid%xb%hf(i+1,j,1)) + &
            dy *(dxm*grid%xb%hf(i,j+1,1) + dx*grid%xb%hf(i+1,j+1,1))

         do k = kts, kte
           if (model_z(k) >= obs_terr) exit

           model_z(k+1) = dym*(dxm*grid%xb%hf(i,j,k+1)   + dx*grid%xb%hf(i+1,j,k+1)) + &
              dy *(dxm*grid%xb%hf(i,j+1,k+1) + dx*grid%xb%hf(i+1,j+1,k+1))
           model_q(k) = dym*(dxm*grid%xb%q(i,j,k)   + dx*grid%xb%q(i+1,j,k)) + &
              dy *(dxm*grid%xb%q(i,j+1,k) + dx*grid%xb%q(i+1,j+1,k))
           model_rho(k) = dym*(dxm*grid%xb%rho(i,j,k)   + dx*grid%xb%rho(i+1,j,k)) + &
              dy *(dxm*grid%xb%rho(i,j+1,k) + dx*grid%xb%rho(i+1,j+1,k))

           ddpw = dpw

           if (model_z(k+1) <= obs_terr) then
             model_dq  (k) = model_rho(k) *(model_z(k+1) - model_z(k))*ddpw 
             model_drho(k) = model_q(k)   *(model_z(k+1) - model_z(k))*ddpw 
           else
             model_dq  (k) = model_rho(k) *(obs_terr-model_z(k))*ddpw 
             model_drho(k) = model_q(k)   *(obs_terr-model_z(k))*ddpw 
           end if 

           ! No feedback to x%hf was considered here (Refer to comments in
           ! da_transform_xtoy_gpspw.inc).       Y.-R. Guo  07/15/2002
           !......................................................................... 

           jo_grad_x%q(i,j,k)       = jo_grad_x%q(i,j,k)       + dym*dxm*model_dq(k)
           jo_grad_x%q(i+1,j,k)     = jo_grad_x%q(i+1,j,k)     + dym*dx*model_dq(k)
           jo_grad_x%q(i,j+1,k)     = jo_grad_x%q(i,j+1,k)     + dy*dxm*model_dq(k)
           jo_grad_x%q(i+1,j+1,k)   = jo_grad_x%q(i+1,j+1,k)   + dy*dx*model_dq(k)

           jo_grad_x%rho(i,j,k)     = jo_grad_x%rho(i,j,k)     + dym*dxm*model_drho(k)
           jo_grad_x%rho(i+1,j,k)   = jo_grad_x%rho(i+1,j,k)   + dym*dx*model_drho(k)
           jo_grad_x%rho(i,j+1,k)   = jo_grad_x%rho(i,j+1,k)   + dy*dxm*model_drho(k)
           jo_grad_x%rho(i+1,j+1,k) = jo_grad_x%rho(i+1,j+1,k) + dy*dx*model_drho(k)
        end do
      end if 
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpspw_adj")

end subroutine da_transform_xtoy_gpspw_adj


subroutine da_transform_xtoy_gpsztd ( grid, iv, y )
!----------------------------------------------------------------
! Purpose: To calculate the observed ZTD with the height
!          correction in considering the differexbnce between
!          the model terrain height and receiver height.
!
!    The logic is similar to the Gpspw: 
!
!      when the receiver below the model surface, a correction,
!      "dzd" should be subtructed from the observed ZTD;
!      when the receiver higher than the model suface, a
!      correction "ztd" should be added to the observed ZTD.
!          
!----------------------------------------------------------------

   IMPLICIT NONE

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(INOUT) :: y        ! y = h (xa)

   INTEGER                      :: n        ! Loop counter.
   INTEGER                      :: i, j     ! Index dimension.
   REAL                         :: dx, dxm  ! 
   REAL                         :: dy, dym  !
!--   
   INTEGER                :: k          ! Index dimension.
   REAL                   :: dzd, ddzd  ! adjustment pw [mm]
   REAL                   :: obs_terr   ! real terrain height at GPS site [m]
   REAL                   :: model_terr ! model terrain height at GPS site[m]
   REAL,DIMENSION(kts:kte):: model_ztd   ! model ztd at GPS site [m]
   REAL,DIMENSION(kts:kte):: model_z     ! model height at GPS site [m]
!--
   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpsztd")

      y % gpspw(:)% tpw = 0.0

      do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2

         i   = iv%info(gpspw)%i(1,n)
         dy  = iv%info(gpspw)%dy(1,n)
         dym = iv%info(gpspw)%dym(1,n)
         j   = iv%info(gpspw)%j(1,n)
         dx  = iv%info(gpspw)%dx(1,n)
         dxm = iv%info(gpspw)%dxm(1,n)

! Mathematical transformation:

         dzd = 0.0
         obs_terr   = iv%info(gpspw)%elv(n)
         model_terr = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + &
                      dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))

         if ( obs_terr <= model_terr ) then 

            model_ztd(1) = dym*(dxm*grid%xa%ref(i,j,1)   + dx*grid%xa%ref(i+1,j,1)) + &
                           dy *(dxm*grid%xa%ref(i,j+1,1) + dx*grid%xa%ref(i+1,j+1,1))

            dzd =  model_ztd(1) * ( obs_terr - model_terr )

         else 

            model_z(1) = dym*(dxm*grid%xb%hf(i,j,1)   + dx*grid%xb%hf(i+1,j,1)) + &
                         dy *(dxm*grid%xb%hf(i,j+1,1) + dx*grid%xb%hf(i+1,j+1,1))

            do k = kts, kte

              if (model_z(k) >= obs_terr ) exit

              model_z(k+1) = dym*(dxm*grid%xb%hf(i,j,k+1)   + dx*grid%xb%hf(i+1,j,k+1)) + &
                             dy *(dxm*grid%xb%hf(i,j+1,k+1) + dx*grid%xb%hf(i+1,j+1,k+1))

              model_ztd(k) = dym*(dxm*grid%xa%ref(i,j,k)   + dx*grid%xa%ref(i+1,j,k)) + &
                             dy *(dxm*grid%xa%ref(i,j+1,k) + dx*grid%xa%ref(i+1,j+1,k))
!
! Here assumed that "model_z" is constant, i.e. grid%xa%hf=0.0. With MM5, 
! this is true, but with WRF, grid%xa%hf may not be zero (?). In the WRF 
! model space (x,y,znu), only the "znu" is constant, but all variables 
! including height could be changed at the "znu" levels. So here is only 
! an approximation for WRF. The following few lines of code is written
! by Y.-R. Guo 07/16/2004.
!
              if ( model_z(k+1) <= obs_terr ) then
                 ddzd = model_ztd(k) * ( model_z(k+1) - model_z(k) )
              else 
                 ddzd = model_ztd(k) * ( obs_terr -  model_z(k) )
              endif

              dzd = dzd + ddzd

            end do 
         end if 

         y % gpspw(n)% tpw = dym* ( dxm * grid%xa%ztd(i,j) + &
                                    dx  * grid%xa%ztd(i+1,j) ) + &
                             dy * ( dxm * grid%xa%ztd(i,j+1) + &
                                    dx  * grid%xa%ztd(i+1,j+1) ) &
                             + 1.e-4 * dzd

      end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpsztd")

end subroutine da_transform_xtoy_gpsztd

subroutine da_transform_xtoy_gpsztd_adj( grid, iv, jo_grad_y, jo_grad_x )
!----------------------------------------------------------------
!  HISTORY
!
!    Purpose:  Considering the difference in the elevation 
!              between model surface and GPS ZTD site
!
!                                   Y.-R. Guo 05/21/2008
!----------------------------------------------------------------

   IMPLICIT NONE

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   INTEGER                       :: n        ! Loop counter.
   INTEGER                       :: i, j     ! Index dimension.
   REAL                          :: dx, dxm  !
   REAL                          :: dy, dym  !

!-- 2004-04-08
   INTEGER                :: k        ! Index dimension.
   REAL                   :: dzd,ddzd     ! adjustment pw [mm]
   REAL                   :: obs_terr     ! real terrain height at GPS site [m]
   REAL                   :: model_terr   ! model terrain height at GPS site[m]
   REAL,DIMENSION(kts:kte):: model_z     ! model z at GPS site [m]
   REAL,DIMENSION(kts:kte):: model_dztd  ! model dq at GPS site [kg/kg]
!--
   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpspw_adj")

   do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2

         i   = iv%info(gpspw)%i(1,n)
         dy  = iv%info(gpspw)%dy(1,n)
         dym = iv%info(gpspw)%dym(1,n)
         j   = iv%info(gpspw)%j(1,n)
         dx  = iv%info(gpspw)%dx(1,n)
         dxm = iv%info(gpspw)%dxm(1,n)

!  Initialise the varibles
         dzd             = 0.
         model_dztd(:)   = 0.

         obs_terr   = iv%info(gpspw)%elv(n)
         model_terr = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + &
                      dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))
         
         dzd =       1.e-4 * jo_grad_y%gpspw(n)%tpw

         jo_grad_x%ztd(i  ,j  ) = jo_grad_x%ztd(i  ,j  ) &
                                + dxm*dym*jo_grad_y%gpspw(n)%tpw
         jo_grad_x%ztd(i+1,j  ) = jo_grad_x%ztd(i+1,j  ) &
                                + dym*dx *jo_grad_y%gpspw(n)%tpw
         jo_grad_x%ztd(i  ,j+1) = jo_grad_x%ztd(i  ,j+1) &
                                + dxm *dy*jo_grad_y%gpspw(n)%tpw
         jo_grad_x%ztd(i+1,j+1) = jo_grad_x%ztd(i+1,j+1) &
                                + dx *dy *jo_grad_y%gpspw(n)%tpw

         IF ( obs_terr <= model_terr ) THEN 

            model_dztd(1)  =  (obs_terr - model_terr) * dzd

            jo_grad_x%ref(i,j,1)     = jo_grad_x%ref(i,j,1)     + dym*dxm*model_dztd(1)
            jo_grad_x%ref(i+1,j,1)   = jo_grad_x%ref(i+1,j,1)   + dym*dx *model_dztd(1)
            jo_grad_x%ref(i,j+1,1)   = jo_grad_x%ref(i,j+1,1)   + dy *dxm*model_dztd(1)
            jo_grad_x%ref(i+1,j+1,1) = jo_grad_x%ref(i+1,j+1,1) + dy *dx *model_dztd(1)

         ELSE 

            model_z(1) = dym*(dxm*grid%xb%hf(i,j,1)   + dx*grid%xb%hf(i+1,j,1)) + &
                         dy *(dxm*grid%xb%hf(i,j+1,1) + dx*grid%xb%hf(i+1,j+1,1))

! ..............................................................................
! The following part of code is written by Y.-R. Guo             07/16/2004

            do k = kts, kte

               if ( model_z(k) >= obs_terr ) exit

               model_z(k+1) = dym*(dxm*grid%xb%hf(i,j,k+1)   + dx*grid%xb%hf(i+1,j,k+1)) + &
                              dy *(dxm*grid%xb%hf(i,j+1,k+1) + dx*grid%xb%hf(i+1,j+1,k+1))
               
               ddzd = dzd

               if ( model_z(k+1) <= obs_terr ) then
                 model_dztd(k) = (model_z(k+1) - model_z(k))*ddzd 
               else
                 model_dztd(k) = (obs_terr - model_z(k))*ddzd 
               end if 
!
! No feedback to x%hf was considered here (Refer to comments in
! DA_Transform_XToY_Gpspw.inc).       Y.-R. Guo  04/06/2005
! ..................................................................................... 
               jo_grad_x%ref(i,j,k)     = jo_grad_x%ref(i,j,k)     + dym*dxm*model_dztd(k)
               jo_grad_x%ref(i+1,j,k)   = jo_grad_x%ref(i+1,j,k)   + dym*dx *model_dztd(k)
               jo_grad_x%ref(i,j+1,k)   = jo_grad_x%ref(i,j+1,k)   + dy *dxm*model_dztd(k)
               jo_grad_x%ref(i+1,j+1,k) = jo_grad_x%ref(i+1,j+1,k) + dy *dx *model_dztd(k)

             end do
         end if
      end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpsztd_adj")

end subroutine da_transform_xtoy_gpsztd_adj
subroutine da_check_max_iv_gpspw(iv, it, num_qcstat_conv) 

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-------------------------------------------------------------------------

   implicit none


   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   logical :: failed 
   integer :: n

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_gpspw")

   ! TPW:

   do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2

      failed=.false.
      if( iv%gpspw(n)%tpw%qc >= obs_qc_pointer ) then 
      call da_max_error_qc(it, iv%info(gpspw), n, iv%gpspw(n)%tpw, max_error_pw,failed)
      if( iv%info(gpspw)%proc_domain(1,n) ) then
       num_qcstat_conv(1,gpspw,7,1) = num_qcstat_conv(1,gpspw,7,1) + 1
      if(failed) then
       num_qcstat_conv(2,gpspw,7,1) = num_qcstat_conv(2,gpspw,7,1) + 1
       write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
       'gpspw',ob_vars(7),iv%info(gpspw)%lat(1,n),iv%info(gpspw)%lon(1,n),'1013.25'                     
      end if
      end if
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_gpspw")

end subroutine da_check_max_iv_gpspw


subroutine da_get_innov_vector_gpspw (it,num_qcstat_conv, grid, ob, iv)

   !----------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !----------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid    ! first guess state.
   type(y_type),     intent(inout) :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: n        ! Loop counter.
   integer :: i, j     ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real    :: model_tpw! Model value u at oblocation.

   !--------------------------------------------------------------------------
   integer :: k            ! Index dimension
   real    :: dpw, ddpw    ! adjustment pw [mm]
   real    :: obs_terr     ! real terrain height at GPS site [m]
   real    :: model_terr   ! model terrain height at GPS site[m]
   real    :: model_q(kts:kte) ! model q at GPS site [kg/kg]
   real    :: model_z(kts:kte) ! model z at GPS site [m]
   real    :: model_rho(kts:kte) ! model rho at GPS site [kg/m^3]

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_gpspw")

   ! unit_gps = myproc + 140
   !---------------------------------------------------------------------------

   ! GPS TPW Pseudo OBS test:
   if ( (pseudo_var(1:3) == 'tpw' .and. num_pseudo > 0) .and. &
        it == 1 ) then
      ! Deallocate:
      if (iv%info(gpspw)%nlocal > 0) then
         write(unit=stdout, fmt='(a,i4)') 'iv%info(gpspw)%nlocal =', iv%info(gpspw)%nlocal
         deallocate(iv % gpspw)
         deallocate(ob % gpspw)
      end if

      use_gpspwobs = .true.

      iv%info(gpspw)%nlocal   = num_pseudo
      iv%info(gpspw)%plocal(1) = num_pseudo
      iv%info(gpspw)%ntotal   = num_pseudo
      iv%info(pseudo)%nlocal  = 0

      ! Allocate:
      call da_allocate_observations (iv)
      iv%info(gpspw)%n1 = 1
      iv%info(gpspw)%n2 = 1
      allocate (ob % gpspw (1:num_pseudo))
      ob % gpspw(1) % tpw   = 0.0
    
      write(unit=message(1),fmt='(a,i2)')'==> GPS TPW pseudo OBS test: num_pseudo=',num_pseudo
      call da_message(message(1:1))

      iv%info(gpspw)%x(:,1)   = pseudo_x
      iv%info(gpspw)%y(:,1)   = pseudo_y
      iv%info(gpspw)%i(:,1)   = int(pseudo_x)
      iv%info(gpspw)%j(:,1)   = int(pseudo_y)
      iv%info(gpspw)%dx(:,1)  = pseudo_x-real(iv%info(gpspw)%i(1,1))
      iv%info(gpspw)%dy(:,1)  = pseudo_y-real(iv%info(gpspw)%j(1,1))
      iv%info(gpspw)%dxm(:,1) = 1.0 - iv%info(gpspw)%dx(1,1)
      iv%info(gpspw)%dym(:,1) = 1.0 - iv%info(gpspw)%dy(1,1)

      iv % gpspw(1) % tpw % inv   = pseudo_val
      iv % gpspw(1) % tpw % qc    = 0
      iv % gpspw(1) % tpw % error = pseudo_err

      ! To consider the site elevation effect, set the model terrain height
      ! to elevation for pseudo OBS:

      i   = iv%info(gpspw)%i(1,1)
      j   = iv%info(gpspw)%j(1,1)
      dx  = iv%info(gpspw)%dx(1,1)
      dy  = iv%info(gpspw)%dy(1,1)
      dxm = iv%info(gpspw)%dxm(1,1)
      dym = iv%info(gpspw)%dym(1,1)

      iv%info(gpspw)%elv(1) = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + & 
                             dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))

      ! Set halo:
      if ((iv%info(gpspw)%i(1,1) < its-1) .or.(iv%info(gpspw)%i(1,1) > ite) .or. & 
          (iv%info(gpspw)%j(1,1) < jts-1) .or.(iv%info(gpspw)%j(1,1) > jte)) then 
         call da_error("da_get_innov_vector_gpspw.inc",94,(/"Should never have obs outside halo by now"/))
         iv%info(gpspw)%proc_domain(:,1) = .false. 
      else 
         iv%info(gpspw)%proc_domain(:,1) = .false.  

         if (iv%info(gpspw)%i(1,1) >= its .and. iv%info(gpspw)%i(1,1) <= ite .and. &  
             iv%info(gpspw)%j(1,1) >= jts .and. iv%info(gpspw)%j(1,1) <= jte) then  
            iv%info(gpspw)%proc_domain(:,1) = .true.  
         end if  
      end if 

      write(unit=message(1),fmt='(a,2f15.5)') pseudo_var, pseudo_val, pseudo_err
      write(unit=message(2),fmt='(3f15.5)')   pseudo_x,  pseudo_y,  pseudo_z
      write(unit=message(3),fmt='(a,f8.2)')   'iv%gpspw: elv=',iv%info(gpspw)%elv(1)
      call da_message(message(1:3))
   end if 

   if (iv%info(gpspw)%nlocal > 0) then


      do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2

         if( iv % gpspw(n) % tpw % qc == fails_error_max .and. it > 1) &
             iv % gpspw(n) % tpw % qc = 0

         ! [1.1] Get horizontal interpolation weights:

         i   = iv%info(gpspw)%i(1,n)
         j   = iv%info(gpspw)%j(1,n)
         dx  = iv%info(gpspw)%dx(1,n)
         dy  = iv%info(gpspw)%dy(1,n)
         dxm = iv%info(gpspw)%dxm(1,n)
         dym = iv%info(gpspw)%dym(1,n)

         model_tpw = dym*(dxm*grid%xb%tpw(i,j)   + dx*grid%xb%tpw(i+1,j)) + &
                     dy *(dxm*grid%xb%tpw(i,j+1) + dx*grid%xb%tpw(i+1,j+1))

            ! To compute the 'inv':
            if (.not.(pseudo_var(1:3) == 'tpw' .and. num_pseudo > 0) ) &
            iv % gpspw(n) % tpw % inv = 0.0
            if (ob % gpspw(n) % tpw > missing_r .AND. &
                iv % gpspw(n) % tpw % qc >= obs_qc_pointer) then
               dpw = 0.0
               obs_terr   = iv%info(gpspw)%elv(n)
               model_terr = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + &
                          dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))
               if (obs_terr <= model_terr) then
                  model_q(1) = dym*(dxm*grid%xb%q(i,j,1)   + dx*grid%xb%q(i+1,j,1)) + &
                             dy *(dxm*grid%xb%q(i,j+1,1) + dx*grid%xb%q(i+1,j+1,1))
                  model_rho(1) = dym*(dxm*grid%xb%rho(i,j,1) + dx*grid%xb%rho(i+1,j,1)) + &
                             dy *(dxm*grid%xb%rho(i,j+1,1) + dx*grid%xb%rho(i+1,j+1,1))
                  dpw = model_rho(1) * model_q(1) *( obs_terr - model_terr)
               else
                  model_z(1) = dym*(dxm*grid%xb%hf(i,j,1)   + dx*grid%xb%hf(i+1,j,1)) + &
                            dy *(dxm*grid%xb%hf(i,j+1,1) + dx*grid%xb%hf(i+1,j+1,1))
                  do k = kts, kte
                     if (model_z(k) >= obs_terr) exit

                     model_z(k+1) = dym*(dxm*grid%xb%hf(i,j,k+1)   + dx*grid%xb%hf(i+1,j,k+1)) + &
                                dy *(dxm*grid%xb%hf(i,j+1,k+1) + dx*grid%xb%hf(i+1,j+1,k+1))
                     model_q(k) = dym*(dxm*grid%xb%q(i,j,k)   + dx*grid%xb%q(i+1,j,k)) + &
                              dy *(dxm*grid%xb%q(i,j+1,k) + dx*grid%xb%q(i+1,j+1,k))
                     model_rho(k) = dym*(dxm*grid%xb%rho(i,j,k)   + dx*grid%xb%rho(i+1,j,k)) + &
                                dy *(dxm*grid%xb%rho(i,j+1,k) + dx*grid%xb%rho(i+1,j+1,k))

                     if (model_z(k+1) <= obs_terr) then
                        ddpw = model_rho(k) * model_q(k) *( model_z(k+1) - model_z(k))
                     else
                        ddpw = model_rho(k) * model_q(k) *( obs_terr - model_z(k))
                     end if

                     dpw = dpw + ddpw
                  end do
               end if
               if ( (pseudo_var(1:3) == 'tpw' .and. num_pseudo > 0) .and. it == 1 ) then
               ! To compute the 'ob':
                 ob % gpspw(n) % tpw = iv % gpspw(n) % tpw % inv + model_tpw - 0.1*dpw
               else
                 iv % gpspw(n) % tpw % inv = ob % gpspw(n) % tpw - model_tpw + 0.1*dpw
               end if
            end if

      end do

      !------------------------------------------------------------------------
      ! [5.0] Perform optional maximum error check:
      !------------------------------------------------------------------------
      if (.not.(pseudo_var(1:3) == 'tpw' .and. num_pseudo > 0) .and. check_max_iv ) &
         call da_check_max_iv_gpspw(iv, it, num_qcstat_conv) 
   end if
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_gpspw")

end subroutine da_get_innov_vector_gpspw


SUBROUTINE da_get_innov_vector_gpsztd ( it, num_qcstat_conv, grid, ob, iv )
!----------------------------------------------------------------
! Innovation for Ground-based ZTD.
!
! Because we can ONLY assimilate either GPS PW or GPS ZTD,
! never assimilate both of them simultaneously, here we 
! used the PW structure for ZTD to avoid declaration of the
! another structure.          
!                                 Y.-R. Guo           05/21/2008
!    Updated for Analysis on Arakawa-C grid
!    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
!----------------------------------------------------------------
   IMPLICIT NONE

!-----
!    INCLUDE 'mpif.h'
!-----

   integer,       intent(in)    :: it      ! External iteration.
   type(domain),  intent(in)    :: grid    ! first guess state
   type(y_type),  intent(inout) :: ob      ! Observation structure.
   type(iv_type), intent(inout) :: iv      ! O-B structure.
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   INTEGER                      :: n        ! Loop counter.
   INTEGER                      :: i, j     ! Index dimension.
   REAL                         :: dx, dxm  ! Interpolation weights.
   REAL                         :: dy, dym  ! Interpolation weights.
   REAL                         :: mdl_ztd  ! Model value u at oblocation.
   INTEGER           :: ittpw,ittpwf

!--------------------------------------------------------------------------
   INTEGER                :: k            ! Index dimension
   REAL                   :: dzd, ddzd    ! adjustment pw [mm]
   REAL                   :: obs_terr     ! real terrain height at GPS site [m]
   REAL                   :: model_terr   ! model terrain height at GPS site[m]
   REAL,DIMENSION(kts:kte):: model_ztd    ! model q at GPS site [kg/kg]
   REAL,DIMENSION(kts:kte):: model_z      ! model z at GPS site [m]
   INTEGER                      :: myrank, ierr, unit_gps

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_gpsztd")

   myrank=0
   unit_gps = myrank + 140
!---------------------------------------------------------------------------

!
! GPS ZTD Pseudo OBS test:
!
   if ( (pseudo_var(1:3) == 'ztd' .and. num_pseudo > 0) .and. &
        it == 1 ) then

! Deallocate:
     if (iv%info(gpspw)%nlocal > 0) then
       write(unit=stdout, fmt='(a,i4)') 'iv%num_gpsztd =', iv%info(gpspw)%nlocal
       deallocate (iv % gpspw)
       deallocate (ob % gpspw)
     endif
     
     use_gpsztdObs = .true.

! Allocate:
     iv %info(gpspw)%nlocal  = num_pseudo
     iv %info(gpspw)%plocal(1)= num_pseudo
     iv %info(gpspw)%ntotal  = num_pseudo
     iv %info(pseudo)%nlocal = 0

     call da_allocate_observations (iv)
     iv%info(gpspw)%n1 = 1
     iv%info(gpspw)%n2 = 1
     allocate (ob % gpspw (1:num_pseudo))
     ob % gpspw(1) % tpw   = 0.0

      write(unit=message(1),fmt='(a,i2)') '==> GPS ZTD pseudo OBS test: num_pseudo=',num_pseudo
      call da_message(message(1:1))
     
     iv % info(gpspw) % x(:,1)  = pseudo_x
     iv % info(gpspw) % y(:,1)  = pseudo_y
     iv % info(gpspw) % i(:,1)  = int(pseudo_x)
     iv % info(gpspw) % j(:,1)  = int(pseudo_y)
     iv % info(gpspw) % dx(:,1) = pseudo_x-real(iv % info(gpspw) % i(1,1))
     iv % info(gpspw) % dy(:,1) = pseudo_y-real(iv % info(gpspw) % j(1,1))
     iv % info(gpspw) % dxm(:,1)= 1.0 - iv % info(gpspw) % dx(1,1)
     iv % info(gpspw) % dym(:,1)= 1.0 - iv % info(gpspw) % dy(1,1)

     iv % gpspw(1) % tpw % inv   = pseudo_val
     iv % gpspw(1) % tpw % qc    = 0
     iv % gpspw(1) % tpw % error = pseudo_err

! To consider the site elevation effect, set the model terrain height
! to elevation for pseudo OBS:

     i   = iv%info(gpspw)%i(1,1)
     j   = iv%info(gpspw)%j(1,1)
     dx  = iv%info(gpspw)%dx(1,1)
     dy  = iv%info(gpspw)%dy(1,1)
     dxm = iv%info(gpspw)%dxm(1,1)
     dym = iv%info(gpspw)%dym(1,1)

     iv%info(gpspw)%elv(1) = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + & 
                            dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))

! Set halo:
     if((iv%info(gpspw)%i(1,1) < its-1) .or. (iv%info(gpspw)%i(1,1) > ite) .or. & 
        (iv%info(gpspw)%j(1,1) < jts-1) .or. (iv%info(gpspw)%j(1,1) > jte)) then 
        call da_error("da_get_innov_vector_gpsztd.inc",106,(/"Should never have obs outside halo by now"/))
        iv%info(gpspw)%proc_domain(:,1) = .false. 
     else 
        iv%info(gpspw)%proc_domain(:,1) = .false.  
     
        if(iv%info(gpspw)%i(1,1) >= its .and. iv%info(gpspw)%i(1,1) <= ite .and. &  
           iv%info(gpspw)%j(1,1) >= jts .and. iv%info(gpspw)%j(1,1) <= jte) then  
           iv%info(gpspw)%proc_domain(:,1) = .true.  
        endif  
     endif 

     write(unit=message(1),fmt='(a,2f15.5)') pseudo_var, pseudo_val, pseudo_err
     write(unit=message(2),fmt='(3f15.5)')   pseudo_x,  pseudo_y,  pseudo_z
     write(unit=message(3),fmt='(a,f8.2)')   'iv%gpsztd: elv=',iv%info(gpspw)%elv(1)
     call da_message(message(1:3))
   end if

! ----------------------------------------------------------------------------

   if ( iv%info(gpspw)%nlocal > 0 ) then

   ittpw   = 0 ; ittpwf  = 0

    write(unit=unit_gps,fmt='(3x,a3,12a10)') ' n ','     lat  ',  &
                       '     lon  ', '  obs ght ', '  mdl ght ',  &
                       ' obsh-mdlh', '   obs ztd', ' model ztd',  &
                       '   O-B ztd', '    Dztd  ', '  O-B+Dztd',  &
                       '   Obs_err', '    qc    ' 

      do n=iv%info(gpspw)%n1,iv%info(gpspw)%n2

         if( iv % gpspw(n) % tpw % qc == fails_error_max .and. it > 1) &
             iv % gpspw(n) % tpw % qc = 0

!        [1.1] Get horizontal interpolation weights:

         i   = iv%info(gpspw)%i(1,n)
         j   = iv%info(gpspw)%j(1,n)
         dx  = iv%info(gpspw)%dx(1,n)
         dy  = iv%info(gpspw)%dy(1,n)
         dxm = iv%info(gpspw)%dxm(1,n)
         dym = iv%info(gpspw)%dym(1,n)

         mdl_ztd   = dym*(dxm*grid%xb%ztd(i,j)   + dx*grid%xb%ztd(i+1,j)) + &
                     dy *(dxm*grid%xb%ztd(i,j+1) + dx*grid%xb%ztd(i+1,j+1))

! To compute the 'inv':
         if ( .not.(pseudo_var(1:3) == 'ztd' .and. num_pseudo > 0) ) &
         iv % gpspw(n) % tpw % inv = 0.0
         if ( ob % gpspw(n) % tpw > missing_r .and. &
                 iv % gpspw(n) % tpw % qc >= obs_qc_pointer ) then

            dzd = 0.0
            obs_terr   = iv%info(gpspw)%elv(n)
            model_terr = dym*(dxm*grid%xb%terr(i,j)   + dx*grid%xb%terr(i+1,j)) + &
                         dy *(dxm*grid%xb%terr(i,j+1) + dx*grid%xb%terr(i+1,j+1))

            if ( obs_terr <= model_terr ) then

               model_ztd(1) = dym*(dxm*grid%xb%ref(i,j,1)   + dx*grid%xb%ref(i+1,j,1)) + &
                              dy *(dxm*grid%xb%ref(i,j+1,1) + dx*grid%xb%ref(i+1,j+1,1))
               dzd = model_ztd(1) * ( obs_terr - model_terr )

            else

               model_z(1) = dym*(dxm*grid%xb%hf(i,j,1)   + dx*grid%xb%hf(i+1,j,1)) + &
                            dy *(dxm*grid%xb%hf(i,j+1,1) + dx*grid%xb%hf(i+1,j+1,1))

               do k = kts, kte
                  
                  if (model_z(k) >= obs_terr ) exit

                  model_z(k+1) = dym*(dxm*grid%xb%hf(i,j,k+1)   + dx*grid%xb%hf(i+1,j,k+1)) + &
                                 dy *(dxm*grid%xb%hf(i,j+1,k+1) + dx*grid%xb%hf(i+1,j+1,k+1))
                  model_ztd(k) = dym*(dxm*grid%xb%ref(i,j,k)   + dx*grid%xb%ref(i+1,j,k)) + &
                                 dy *(dxm*grid%xb%ref(i,j+1,k) + dx*grid%xb%ref(i+1,j+1,k))
                  
                  if ( model_z(k+1) <= obs_terr ) then
                    ddzd = model_ztd(k) * ( model_z(k+1) - model_z(k) )
                  else
                    ddzd = model_ztd(k) * ( obs_terr - model_z(k) )
                  endif

                  dzd = dzd + ddzd
               end do
            end if

            if ( (pseudo_var(1:3) == 'ztd' .and. num_pseudo > 0) .and. it == 1 ) then

! To compute the 'ob':
              ob % gpspw(n) % tpw = iv % gpspw(n) % tpw % inv + mdl_ztd - 1.e-4 * dzd

            else


              iv % gpspw(n) % tpw % inv = ob % gpspw(n) % tpw - mdl_ztd &
                                          + 1.e-4 * dzd
!
! Overwrite the observation error specification (YRG):
!
!              iv % gpspw(n) % tpw % error = 1.0 + 0.02*(ob%gpspw(n)%tpw-200.)   

            end if
         endif
!---   
        write(unit=unit_gps, fmt='(i4,11f10.3,i7)') n, &
              iv%info(gpspw)%lat(1,n), iv%info(gpspw)%lon(1,n), obs_terr, &
              model_terr, obs_terr - model_terr, ob%gpspw(n)%tpw,   &
              mdl_ztd , ob%gpspw(n)%tpw-mdl_ztd, 1.e-4*dzd,       &
              ob%gpspw(n)%tpw-mdl_ztd+1.e-4*dzd, iv%gpspw(n)%tpw%error,&
              iv%gpspw(n)%tpw%qc
!---   
      end do

!------------------------------------------------------------------------
!        [5.0] Perform optional maximum error check:
!------------------------------------------------------------------------
      if ( .not.(pseudo_var(1:3) == 'ztd' .and. num_pseudo > 0) .and. check_max_iv ) &
         call da_check_max_iv_gpspw(iv, it, num_qcstat_conv)
   end if

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_gpsztd")

end subroutine da_get_innov_vector_gpsztd

subroutine da_calculate_grady_gpspw(iv, re, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_gpspw")

   do n=1, iv%info(gpspw)%nlocal
      if (iv%gpspw(n)%tpw%qc < obs_qc_pointer) then
         re%gpspw(n)%tpw = 0.0
      end if

      jo_grad_y%gpspw(n)%tpw = -re%gpspw(n)%tpw / (iv%gpspw(n)%tpw%error * iv%gpspw(n)%tpw%error)
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_gpspw")


end subroutine da_calculate_grady_gpspw




end module da_gpspw

