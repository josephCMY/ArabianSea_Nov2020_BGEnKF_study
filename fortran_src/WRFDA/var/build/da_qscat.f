












module da_qscat 

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, max_error_uv, max_error_t, rootproc, &
      qscat, max_error_p,max_error_q, trace_use_dull, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, anal_type_verify, kms,kme,kts,kte,&
      ob_vars,qcstat_conv_unit, fails_error_max, &
      convert_fd2uv,convert_uv2fd,max_error_spd,max_error_dir,max_omb_spd,max_omb_dir,pi,qc_rej_both,&
      wind_sd_qscat, wind_stats_sd 
   use da_grid_definitions, only : da_ffdduv, da_ffdduv_model, da_ffdduv_diagnose
   use da_physics, only : da_uv_to_sd_lin, da_uv_to_sd_adj
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_to_zk, &
      da_interp_lin_3d,da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_qscat1_type
      real          :: u                        
      real          :: v                        
   end type residual_qscat1_type

   type maxmin_qscat_stats_type
      type (maxmin_type)         :: u, v
   end type maxmin_qscat_stats_type

   type stats_qscat_type
      type (maxmin_qscat_stats_type)  :: maximum, minimum
      type (residual_qscat1_type)     :: average, rms_err
   end type stats_qscat_type

contains

subroutine da_jo_and_grady_qscat(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(in)    :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_qscat")

   jo % qscat_u = 0.0
   jo % qscat_v = 0.0

   do n=1, iv%info(qscat)%nlocal
      jo_grad_y%qscat(n)%u = -re%qscat(n)%u / &
         (iv%qscat(n)%u%error * iv%qscat(n)%u%error)
      jo_grad_y%qscat(n)%v = -re%qscat(n)%v / &
         (iv%qscat(n)%v%error * iv%qscat(n)%v%error)

      if (iv%info(qscat)%proc_domain(1,n)) then
         jo % qscat_u = jo % qscat_u - re%qscat(n)%u * jo_grad_y%qscat(n)%u
         jo % qscat_v = jo % qscat_v - re%qscat(n)%v * jo_grad_y%qscat(n)%v
      end if
   end do

   jo % qscat_u = 0.5 * jo % qscat_u
   jo % qscat_v = 0.5 * jo % qscat_v

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_qscat")

     
end subroutine da_jo_and_grady_qscat


subroutine da_residual_qscat(iv, y, re, np_missing, np_bad_data,np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for qscat obs
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

   if (trace_use_dull) call da_trace_entry("da_residual_qscat")

   n_obs_bad % u % num = number_type(0, 0, 0)
   n_obs_bad % v % num = number_type(0, 0, 0)

   do n=1, iv%info(qscat)%nlocal
      np_available = np_available + 2
      re%qscat(n)%u = da_residual(n, 0, y%qscat(n)%u, iv%qscat(n)%u, n_obs_bad % u)
      re%qscat(n)%v = da_residual(n, 0, y%qscat(n)%v, iv%qscat(n)%v, n_obs_bad % v)
   end do

   np_missing  = np_missing  + n_obs_bad % u % num % miss + n_obs_bad % v % num % miss
   np_bad_data = np_bad_data + n_obs_bad % u % num % bad  + n_obs_bad % v % num % bad
   np_obs_used = np_obs_used + n_obs_bad % u % num % use  + n_obs_bad % v % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_qscat")

end subroutine da_residual_qscat


subroutine da_check_max_iv_qscat(iv, it, num_qcstat_conv)        

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   logical :: failed,failed1,failed2
   integer :: n
   
   if (trace_use_dull) call da_trace_entry("da_check_max_iv_qscat")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(qscat)%n1,iv%info(qscat)%n2
      if(.not. qc_rej_both)then
         if(wind_sd_qscat)then
            failed=.false.
            if( iv%qscat(n)%u%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%u, max_error_spd,failed)
                if( iv%info(qscat)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,qscat,1,1) = num_qcstat_conv(1,qscat,1,1) + 1
                    if(failed) then
                       num_qcstat_conv(2,qscat,1,1) = num_qcstat_conv(2,qscat,1,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                       'qscat',ob_vars(1),iv%info(qscat)%lat(1,n),iv%info(qscat)%lon(1,n),'1013.25'
                    end if
                end if
            end if

            failed=.false.
            if( iv%qscat(n)%v%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%v, max_error_dir,failed)
                if( iv%info(qscat)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,qscat,2,1) = num_qcstat_conv(1,qscat,2,1) + 1
                    if(failed)then
                       num_qcstat_conv(2,qscat,2,1) = num_qcstat_conv(2,qscat,2,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                             'qscat',ob_vars(2),iv%info(qscat)%lat(1,n),iv%info(qscat)%lon(1,n),'1013.25'
                    end if
                end if
            end if

         else
            failed=.false.
            if( iv%qscat(n)%u%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%u, max_error_uv,failed)
                if( iv%info(qscat)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,qscat,1,1) = num_qcstat_conv(1,qscat,1,1) + 1
                    if(failed) then
                       num_qcstat_conv(2,qscat,1,1) = num_qcstat_conv(2,qscat,1,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                             'qscat',ob_vars(1),iv%info(qscat)%lat(1,n),iv%info(qscat)%lon(1,n),'1013.25'
                    end if
                end if
            end if
            failed=.false.
            if( iv%qscat(n)%v%qc >= obs_qc_pointer ) then
                call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%v, max_error_uv,failed)
                if( iv%info(qscat)%proc_domain(1,n) ) then
                    num_qcstat_conv(1,qscat,2,1) = num_qcstat_conv(1,qscat,2,1) + 1
                    if(failed)then
                       num_qcstat_conv(2,qscat,2,1) = num_qcstat_conv(2,qscat,2,1) + 1
                       write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                             'qscat',ob_vars(2),iv%info(qscat)%lat(1,n),iv%info(qscat)%lon(1,n),'1013.25'
                    end if
                end if
             end if
          end if

          if(wind_sd_qscat)then
             if(iv%qscat(n)%u%qc == fails_error_max .or. abs(iv%qscat(n)%u%inv) >= max_omb_spd) then
                iv%qscat(n)%u%qc = fails_error_max
                iv%qscat(n)%u%inv = 0.0
             endif
             if(iv%qscat(n)%v%qc == fails_error_max .or. abs(iv%qscat(n)%v%inv) >= max_omb_dir) then
                iv%qscat(n)%v%qc = fails_error_max
                iv%qscat(n)%v%inv = 0.0
             endif
          endif

       else

          failed1=.false.
          failed2=.false.

          if( iv%qscat(n)%v%qc >= obs_qc_pointer .or. iv%qscat(n)%u%qc >= obs_qc_pointer )  then
              if(wind_sd_qscat)then
                 call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%u, max_error_spd,failed1)
                 call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%v, max_error_dir,failed2)
              else
                 call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%u, max_error_uv,failed1)
                 call da_max_error_qc (it,iv%info(qscat), n, iv%qscat(n)%v, max_error_uv,failed2)
              endif

              if( iv%info(qscat)%proc_domain(1,n) ) then
                  num_qcstat_conv(1,qscat,1,1) = num_qcstat_conv(1,qscat,1,1) + 1
                  num_qcstat_conv(1,qscat,2,1) = num_qcstat_conv(1,qscat,2,1) + 1
 
                  if(failed1 .or. failed2) then
                     num_qcstat_conv(2,qscat,1,1) = num_qcstat_conv(2,qscat,1,1) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                           'qscat',ob_vars(1),iv%info(qscat)%lat(1,n),iv%info(qscat)%lon(1,n),'1013.25'
                     num_qcstat_conv(2,qscat,2,1) = num_qcstat_conv(2,qscat,2,1) + 1
                     write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                           'qscat',ob_vars(2),iv%info(qscat)%lat(1,n),iv%info(qscat)%lon(1,n),'1013.25'
                  end if
              end if
          end if

	  if(wind_sd_qscat)then
             if(iv%qscat(n)%u%qc == fails_error_max .or. iv%qscat(n)%v%qc == fails_error_max .or. &
                abs(iv%qscat(n)%v%inv) >= max_omb_dir .or. abs(iv%qscat(n)%u%inv) >= max_omb_spd )then
                iv%qscat(n)%u%qc = fails_error_max
                iv%qscat(n)%v%qc = fails_error_max
                iv%qscat(n)%u%inv = 0.0
                iv%qscat(n)%v%inv = 0.0
             endif
          else
             if(iv%qscat(n)%u%qc == fails_error_max .or. iv%qscat(n)%v%qc == fails_error_max ) then
                iv%qscat(n)%u%qc = fails_error_max
                iv%qscat(n)%v%qc = fails_error_max
                iv%qscat(n)%u%inv = 0.0
                iv%qscat(n)%v%inv = 0.0
             endif
         endif
      endif

    end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_qscat")

end subroutine da_check_max_iv_qscat


subroutine da_get_innov_vector_qscat (it,num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid    ! first guess state.
   type(y_type),     intent(inout) :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real    :: speed, direction

   real, allocatable :: model_u(:,:)  ! Model value u at ob location.
   real, allocatable :: model_v(:,:)  ! Model value v at ob location.

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_qscat")

   allocate (model_u(iv%info(qscat)%max_lev,iv%info(qscat)%n1:iv%info(qscat)%n2))
   allocate (model_v(iv%info(qscat)%max_lev,iv%info(qscat)%n1:iv%info(qscat)%n2))

   if ( it > 1 ) then
      do n=iv%info(qscat)%n1,iv%info(qscat)%n2
         if (iv%qscat(n)%u%qc == fails_error_max) iv%qscat(n)%u%qc = 0
         if (iv%qscat(n)%v%qc == fails_error_max) iv%qscat(n)%v%qc = 0
      end do
   end if

   do n=iv%info(qscat)%n1,iv%info(qscat)%n2

      ! [1.1] Get horizontal interpolation weights:

      i   = iv%info(qscat)%i(1,n)
      j   = iv%info(qscat)%j(1,n)
      dx  = iv%info(qscat)%dx(1,n)
      dy  = iv%info(qscat)%dy(1,n)
      dxm = iv%info(qscat)%dxm(1,n)
      dym = iv%info(qscat)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j,k)+dx*grid%xb%h(i+1,j,k)) + dy*(dxm*grid%xb%h(i,j+1,k)+dx*grid%xb%h(i+1,j+1,k))
      end do

      if (iv % qscat(n) % h > missing_r) then
         call da_to_zk(iv % qscat(n) % h, v_h, v_interp_h, iv%info(qscat)%zk(1,n))
         if (iv%info(qscat)%zk(1,n) < 1.0) then
            iv%info(qscat)%zk(1,n) = 1.0
         end if
      end if
   end do

   call da_convert_zk (iv%info(qscat))

   if (.not. anal_type_verify) then
      do n=iv%info(qscat)%n1,iv%info(qscat)%n2
         if (iv%info(qscat)%zk(1,n) < 0.0) then
            iv%qscat(n)%u%qc = missing_data
            iv%qscat(n)%v%qc = missing_data
         end if
      end do
   end if

   call da_interp_lin_3d (grid%xb%u, iv%info(qscat), model_u)
   call da_interp_lin_3d (grid%xb%v, iv%info(qscat), model_v)

   do n=iv%info(qscat)%n1,iv%info(qscat)%n2

      !------------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! [3.0] Fast interpolation:
      !------------------------------------------------------------------------
          if (wind_sd_qscat) then
              call da_ffdduv_model (speed,direction,model_u(1,n), model_v(1,n), convert_uv2fd)

              if (ob%qscat(n)%u > missing_r .AND. iv%qscat(n)%u%qc >= obs_qc_pointer) then
                  iv%qscat(n)%u%inv = ob%qscat(n)%u - speed
              end if

              if (ob%qscat(n)%v > missing_r .AND. iv%qscat(n)%v%qc >= obs_qc_pointer) then
                  iv%qscat(n)%v%inv = ob%qscat(n)%v - direction
                  if (iv%qscat(n)%v%inv > 180.0 ) iv%qscat(n)%v%inv = iv%qscat(n)%v%inv - 360.0
                  if (iv%qscat(n)%v%inv < -180.0 ) iv%qscat(n)%v%inv = iv%qscat(n)%v%inv + 360.0
              end if
          else
               if (ob % qscat(n) % u > missing_r .AND. &
                   iv % qscat(n) % u % qc >= obs_qc_pointer) then
                   iv % qscat(n) % u % inv = ob % qscat(n) % u - model_u(1,n)
               end if

               if (ob % qscat(n) % v > missing_r .AND. &
                   iv % qscat(n) % v % qc >= obs_qc_pointer) then
                   iv % qscat(n) % v % inv = ob % qscat(n) % v - model_v(1,n)
               end if
          end if
   end do

   !------------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_qscat(iv, it, num_qcstat_conv)       

   deallocate (model_u)
   deallocate (model_v)

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_qscat")

end subroutine da_get_innov_vector_qscat


subroutine da_ao_stats_qscat (stats_unit, iv, re, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type  (y_type), intent(in)    :: re            ! A - O
   type(y_type),   intent (in)   :: ob            ! Observation structure.

   type (stats_qscat_type)       :: stats
   integer                       :: nu, nv
   integer                       :: n
   real                          :: u_inc, v_inc, u_obs, v_obs
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_qscat")

   nu = 0
   nv = 0

   stats%maximum%u = maxmin_type(-1.0E+20, 0, 0)
   stats%maximum%v = maxmin_type(-1.0E+20, 0, 0)
   stats%minimum%u = maxmin_type (1.0E+20, 0, 0)
   stats%minimum%v = maxmin_type (1.0E+20, 0, 0)
   stats%average = residual_qscat1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(qscat)%nlocal
      if (iv%info(qscat)%proc_domain(1,n)) then

         u_inc = re%qscat(n)%u
         v_inc = re%qscat(n)%v
         u_obs = ob%qscat(n)%u
         v_obs = ob%qscat(n)%v

         if (.not. wind_sd_qscat .and. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%qscat(n)%u%qc, iv%qscat(n)%v%qc, convert_uv2fd)
         if (wind_sd_qscat .and. .not. wind_stats_sd) &
            call da_ffdduv_diagnose(u_obs, u_obs, u_inc, v_obs, v_obs, v_inc, &
                                    iv%qscat(n)%u%qc, iv%qscat(n)%v%qc, convert_fd2uv)

         call da_stats_calculate (n, 0, iv%qscat(n)%u%qc, & 
            u_inc, nu, & 
            stats%minimum%u, stats%maximum%u, &
            stats%average%u, stats%rms_err%u)
         call da_stats_calculate (n, 0, iv%qscat(n)%v%qc, & 
            v_inc, nv, & 
            stats%minimum%v, stats%maximum%v, &
            stats%average%v, stats%rms_err%v)
      end if    ! end if (iv%info(qscat)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   iv%nstats(qscat) = nu + nv
   
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for qscat'
         call da_print_stats_qscat(stats_unit, nu, nv, stats)
      end if
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_qscat")

end subroutine da_ao_stats_qscat


subroutine da_oi_stats_qscat (stats_unit, iv, ob)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI
   type(y_type),   intent (in) :: ob            ! Observation structure.

   type (stats_qscat_type) :: stats
   integer                 :: nu, nv
   integer                 :: n
   real                    :: u_inv, v_inv, u_obs, v_obs

   if (trace_use_dull) call da_trace_entry("da_oi_stats_qscat")

   nu = 0
   nv = 0
   
   stats%maximum%u = maxmin_type(-1.0E+20, 0, 0)
   stats%maximum%v = maxmin_type(-1.0E+20, 0, 0)
   stats%minimum%u = maxmin_type(1.0E+20, 0, 0)
   stats%minimum%v = maxmin_type(1.0E+20, 0, 0)
   stats%average = residual_qscat1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(qscat)%nlocal
      if (iv%info(qscat)%proc_domain(1,n)) then

        u_inv = iv%qscat(n)%u%inv
        v_inv = iv%qscat(n)%v%inv
        u_obs = ob%qscat(n)%u
        v_obs = ob%qscat(n)%v

        if (.not. wind_sd_qscat .and. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%qscat(n)%u%qc, iv%qscat(n)%v%qc, convert_uv2fd)
        if (wind_sd_qscat .and. .not. wind_stats_sd) &
           call da_ffdduv_diagnose(u_obs, u_inv, u_obs, v_obs, v_inv, v_obs, &
                                   iv%qscat(n)%u%qc, iv%qscat(n)%v%qc, convert_fd2uv)

         call da_stats_calculate(n, 0, iv%qscat(n)%u%qc, u_inv, nu, &
            stats%minimum%u, stats%maximum%u, stats%average%u, stats%rms_err%u)
         call da_stats_calculate(n, 0, iv%qscat(n)%v%qc, v_inv, nv, &
            stats%minimum%v, stats%maximum%v, stats%average%v, stats%rms_err%v)
      end if    ! end if (iv%info(qscat)%proc_domain(1,n))
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
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for qscat'
         call da_print_stats_qscat(stats_unit, nu, nv, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_qscat")

end subroutine da_oi_stats_qscat


subroutine da_print_stats_qscat(stats_unit, nu, nv, qscat)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit  
   integer,                 intent(inout) :: nu, nv      
   type (stats_qscat_type), intent(in)    :: qscat

   if (trace_use_dull) call da_trace_entry("da_print_stats_qscat")

      write(unit=stats_unit, fmt='(3a/)') &
        '   var             ', &
        'u (m/s)     n    k    ', &
        'v (m/s)     n    k    '

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
        '  Number: ', nu, nv

   if (nu < 1) nu = 1
   if (nv < 1) nv = 1

   write(unit=stats_unit, fmt='((a,2(f12.4,2i5)))') &
        ' Minimum(n,k): ', qscat%minimum%u, qscat%minimum%v, &
        ' Maximum(n,k): ', qscat%maximum%u, qscat%maximum%v
   
   write(unit=stats_unit, fmt='((a,2(f12.4,10x)))') &
        ' Average     : ', qscat%average%u/real(nu), qscat%average%v/real(nv), &
        '    RMSE     : ', sqrt(qscat%rms_err%u/real(nu)), &
                           sqrt(qscat%rms_err%v/real(nv))

   if (trace_use_dull) call da_trace_exit("da_print_stats_qscat")
   
end subroutine da_print_stats_qscat


subroutine da_transform_xtoy_qscat(grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n        ! Loop counter.

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_qscat")

   allocate (u(1,iv%info(qscat)%n1:iv%info(qscat)%n2))
   allocate (v(1,iv%info(qscat)%n1:iv%info(qscat)%n2))
   allocate (ub(1,iv%info(qscat)%n1:iv%info(qscat)%n2))
   allocate (vb(1,iv%info(qscat)%n1:iv%info(qscat)%n2))

   call da_interp_lin_3d (grid%xa%u, iv%info(qscat), u)
   call da_interp_lin_3d (grid%xa%v, iv%info(qscat), v)
   call da_interp_lin_3d (grid%xb%u, iv%info(qscat), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(qscat), vb)

   do n=iv%info(qscat)%n1,iv%info(qscat)%n2
      if (wind_sd_qscat) then
          call da_uv_to_sd_lin(y%qscat(n)%u,y%qscat(n)%v,u(1,n),v(1,n),ub(1,n),vb(1,n))
      else
          y%qscat(n)%u = u(1,n)
          y%qscat(n)%v = v(1,n)
      end if

   end do

   deallocate (u)
   deallocate (v)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_qscat")

end subroutine da_transform_xtoy_qscat


subroutine da_transform_xtoy_qscat_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none
   type (domain),  intent(in)    :: grid          ! first guess state.
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.

   real, allocatable :: u(:,:)
   real, allocatable :: v(:,:)
   real, allocatable :: ub(:,:)
   real, allocatable :: vb(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_qscat_adj")

   allocate (u(1,iv%info(qscat)%n1:iv%info(qscat)%n2))
   allocate (v(1,iv%info(qscat)%n1:iv%info(qscat)%n2))

   allocate (ub(1,iv%info(qscat)%n1:iv%info(qscat)%n2))
   allocate (vb(1,iv%info(qscat)%n1:iv%info(qscat)%n2))

   call da_interp_lin_3d (grid%xb%u, iv%info(qscat), ub)
   call da_interp_lin_3d (grid%xb%v, iv%info(qscat), vb)

   do n=iv%info(qscat)%n1,iv%info(qscat)%n2
      if (wind_sd_qscat) then
          call da_uv_to_sd_adj(jo_grad_y%qscat(n)%u, &
                               jo_grad_y%qscat(n)%v, u(1,n), v(1,n), ub(1,n), vb(1,n))
      else
          u(1,n) = jo_grad_y%qscat(n)%u
          v(1,n) = jo_grad_y%qscat(n)%v
      end if
   end do

   call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(qscat), u)
   call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(qscat), v)

   deallocate (u)
   deallocate (v)
   deallocate (ub)
   deallocate (vb)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_qscat_adj")

end subroutine da_transform_xtoy_qscat_adj


subroutine da_calculate_grady_qscat(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer :: n
   
   if (trace_use_dull) call da_trace_entry("da_calculate_grady_qscat")

   do n=1, iv%info(qscat)%nlocal
      if (iv%qscat(n)%u%qc < obs_qc_pointer) re%qscat(n)%u = 0.0
      if (iv%qscat(n)%v%qc < obs_qc_pointer) re%qscat(n)%v = 0.0
      jo_grad_y%qscat(n)%u = -re%qscat(n)%u / (iv%qscat(n)%u%error * iv%qscat(n)%u%error)
      jo_grad_y%qscat(n)%v = -re%qscat(n)%v / (iv%qscat(n)%v%error * iv%qscat(n)%v%error)
   end do
   
   if (trace_use_dull) call da_trace_exit("da_calculate_grady_qscat")
        
end subroutine da_calculate_grady_qscat



end module da_qscat
