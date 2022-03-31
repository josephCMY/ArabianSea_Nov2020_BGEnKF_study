












module da_rain 

   use module_domain, only : domain
   
   use module_dm, only : local_communicator, mytask, ntasks, ntasks_x, &
      ntasks_y
   use module_comm_dm, only : halo_em_rain_sub

   use da_control, only : obs_qc_pointer,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing, rootproc, max_error_rain, &
      rain, trace_use,fails_error_max, &
      max_stheight_diff,missing_data,anal_type_verify, &
      anal_type_verify,max_ext_its,qcstat_conv_unit,ob_vars, &
      ids,ide,jds,jde,kds,kde, ims,ime,jms,jme,kms,kme, &
      ips,ipe,jps,jpe,kps,kpe,num_fgat_time

   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_to_zk, &
      da_interp_lin_3d,da_interp_lin_3d_adj, &
      da_interp_lin_2d, da_interp_lin_2d_adj, da_interp_lin_2d_partial
   use da_par_util1, only : da_proc_sum_int
   use da_par_util, only : da_proc_stats_combine
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_obs_sfc_correction, &
       da_convert_zk,map_info,da_llxy_wrf, da_llxy_default
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_rain1_type
      real          :: rain     
   end type residual_rain1_type

   type maxmin_rain_stats_type
      type (maxmin_type)         :: rain
   end type maxmin_rain_stats_type

   type stats_rain_type
      type (maxmin_rain_stats_type)  :: maximum, minimum
      type (residual_rain1_type)     :: average, rms_err
   end type stats_rain_type

contains

subroutine da_ao_stats_rain (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O

   type (stats_rain_type) :: stats
   integer                :: nrain
   integer                :: n
   
   if (trace_use) call da_trace_entry("da_ao_stats_rain")    

   nrain = 0

   stats%maximum%rain = maxmin_type (missing_r, 0, 0)
   stats%minimum%rain = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_rain1_type(0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(rain)%nlocal
      if (iv%info(rain)%proc_domain(1,n)) then
         call da_stats_calculate (n, 0, iv%rain(n)%rain%qc, & 
            re%rain(n)%rain, nrain, & 
            stats%minimum%rain, stats%maximum%rain, &
            stats%average%rain, stats%rms_err%rain)
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nrain)
   iv%nstats(rain) = nrain

   call da_proc_stats_combine(stats%average%rain, stats%rms_err%rain, &
      stats%minimum%rain%value, stats%maximum%rain%value, &
      stats%minimum%rain%n, stats%maximum%rain%n, &
      stats%minimum%rain%l, stats%maximum%rain%l)
   if (rootproc) then
      if (nrain /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for rainfall'
         call da_print_stats_rain(stats_unit, nrain, stats)
      end if
   end if
   
   if (trace_use) call da_trace_exit("da_ao_stats_rain")    

 end subroutine da_ao_stats_rain


subroutine da_jo_and_grady_rain(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type),intent(inout) :: jo          ! Obs cost function.

   integer :: n 

   if (trace_use) call da_trace_entry("da_jo_and_grady_rain")

   jo % rain_r = 0.0
   
   do n=1, iv%info(rain)%nlocal
      jo_grad_y%rain(n)%rain = -re%rain(n)%rain / &
         (iv%rain(n)%rain%error * iv%rain(n)%rain%error)

      if (iv%info(rain)%proc_domain(1,n)) then
         jo % rain_r = jo % rain_r - re%rain(n)%rain * jo_grad_y%rain(n)%rain
      end if
   end do
      
   jo % rain_r = 0.5 * jo % rain_r

   if (trace_use) call da_trace_exit("da_jo_and_grady_rain")

end subroutine da_jo_and_grady_rain


subroutine da_residual_rain(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

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

   type (bad_data_type)          :: n_obs_bad
   integer                       :: n

   if (trace_use) call da_trace_entry("da_residual_rain")

   n_obs_bad % rain % num = number_type(0, 0, 0)

   do n=1, iv%info(rain)%nlocal
      np_available = np_available + 1
      re%rain(n)%rain = da_residual(n, 0, y%rain(n)%rain, iv%rain(n)%rain, n_obs_bad % rain) 
   end do

   np_missing = np_missing + n_obs_bad % rain % num % miss 
   np_bad_data = np_bad_data + n_obs_bad % rain % num % bad  
   np_obs_used = np_obs_used + n_obs_bad % rain % num % use 

   if (trace_use) call da_trace_exit("da_residual_rain")

end subroutine da_residual_rain


subroutine da_oi_stats_rain (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_rain_type) :: stats      
   integer                :: nrain
   integer                :: n

   if (trace_use) call da_trace_entry("da_oi_stats_rain")

   nrain = 0

   stats%maximum%rain = maxmin_type(missing_r, 0, 0)
   stats%minimum%rain = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_rain1_type(0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(rain)%nlocal
      if (iv%info(rain)%proc_domain(1,n)) then
         call da_stats_calculate(iv%info(rain)%obs_global_index(n), &
            0, iv%rain(n)%rain%qc, &
            iv%rain(n)%rain%inv, nrain, &
            stats%minimum%rain, stats%maximum%rain, &
            stats%average%rain, stats%rms_err%rain)
      end if  
   end do

   ! Do inter-processor communication to gather statistics.

   call da_proc_sum_int(nrain)
   call da_proc_stats_combine(stats%average%rain, stats%rms_err%rain, &
      stats%minimum%rain%value, stats%maximum%rain%value, &
      stats%minimum%rain%n, stats%maximum%rain%n, &
      stats%minimum%rain%l, stats%maximum%rain%l)
 
   if (rootproc) then
      if (nrain /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for rainfall'
         call da_print_stats_rain(stats_unit, nrain, stats)
      end if
   end if

   if (trace_use) call da_trace_exit("da_oi_stats_rain")

end subroutine da_oi_stats_rain


subroutine da_print_stats_rain(stats_unit, nrain, rain)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                intent(in)    :: stats_unit
   integer,                intent(inout) :: nrain 
   type (stats_rain_type), intent(in)  :: rain

   if (trace_use) call da_trace_entry("da_print_stats_rain")

   if (nrain > 0) then

      write(unit=stats_unit, fmt='(a/)') '   var        rainfall(mm)     n'

      write(unit=stats_unit, fmt='(a,i14)') '  Number: ', nrain

      write(unit=stats_unit, fmt='(a, f12.4,i5)') &
         ' Minimum(n): ', rain%minimum%rain%value, &
                          rain%minimum%rain%n    , &
         ' Maximum(n): ', rain%maximum%rain%value, &
                          rain%maximum%rain%n
      write(unit=stats_unit, fmt='(a, f12.4,5x)') &
         ' Average   : ', rain%average%rain/real(nrain), &
         '    RMSE   : ', sqrt(rain%rms_err%rain/real(nrain))
   end if

   if (trace_use) call da_trace_exit("da_print_stats_rain")

end subroutine da_print_stats_rain


subroutine da_transform_xtoy_rain(grid, iv, y, hr_rainc, hr_rainnc)
   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !--------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)
   real, dimension(ims:ime,jms:jme), intent(in) :: hr_rainc, hr_rainnc

   integer :: n        ! Loop counter.

   real, allocatable :: model_rainc(:)
   real, allocatable :: model_rainnc(:)

   if (trace_use) call da_trace_entry("da_transform_xtoy_rain")

   allocate (model_rainc(iv%info(rain)%n1:iv%info(rain)%n2))
   allocate (model_rainnc(iv%info(rain)%n1:iv%info(rain)%n2))

   model_rainc=0.0
   model_rainnc=0.0

   call da_interp_lin_2d (hr_rainc, iv%info(rain), 1, model_rainc)
   call da_interp_lin_2d (hr_rainnc, iv%info(rain), 1, model_rainnc)

   do n=iv%info(rain)%n1,iv%info(rain)%n2
      if (iv % rain(n) % rain % qc ==  missing_data) then
         y%rain(n)%rain = 0.0
      else
         y%rain(n)%rain = model_rainc(n) + model_rainnc(n)
      endif
   end do

   deallocate (model_rainc)
   deallocate (model_rainnc)

   if (trace_use) call da_trace_exit("da_transform_xtoy_rain")

end subroutine da_transform_xtoy_rain

subroutine da_transform_xtoy_rain_adj(grid, iv, jo_grad_y, a_hr_rainc, a_hr_rainnc)

   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !--------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   real, dimension(ims:ime,jms:jme), intent(inout) :: a_hr_rainc, a_hr_rainnc

   integer :: n        ! Loop counter.
   integer :: i, j     ! Index dimension.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.

   real, allocatable :: model_rainnc(:)
   real, allocatable :: model_rainc(:)

   if (trace_use) call da_trace_entry("da_transform_xtoy_rain_adj")

      allocate (model_rainnc(iv%info(rain)%n1:iv%info(rain)%n2))
      allocate (model_rainc(iv%info(rain)%n1:iv%info(rain)%n2))

      model_rainnc=0.0
      model_rainc=0.0

      do n=iv%info(rain)%n1,iv%info(rain)%n2
         model_rainnc(n)  = model_rainnc(n) + jo_grad_y%rain(n)%rain
         model_rainc(n)   = model_rainc(n) + jo_grad_y%rain(n)%rain
         jo_grad_y%rain(n)%rain=0.0
      end do
      
      call da_interp_lin_2d_adj (a_hr_rainc, iv%info(rain), 1, model_rainc)
      call da_interp_lin_2d_adj (a_hr_rainnc, iv%info(rain), 1, model_rainnc)

      deallocate (model_rainc)
      deallocate (model_rainnc)

   if (trace_use) call da_trace_exit("da_transform_xtoy_rain_adj")

end subroutine da_transform_xtoy_rain_adj


subroutine da_check_max_iv_rain(iv,ob, it, num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it      ! Outer iteration 
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)
   type(y_type),  intent(in)    :: ob      ! Observation structure.

   logical :: failed 
   integer :: n
   
   if (trace_use) call da_trace_entry("da_check_max_iv_rain")       


   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(rain)%n1,iv%info(rain)%n2
      failed=.false.
      if ( iv%rain(n)%rain%qc >= obs_qc_pointer ) then 
         call da_max_error_qc (it, iv%info(rain), n, iv%rain(n)%rain, max_error_rain, failed)
         if ( iv%info(rain)%proc_domain(1,n) ) then
            num_qcstat_conv(1,rain,10,1)= num_qcstat_conv(1,rain,10,1) + 1
            if (failed) then
               num_qcstat_conv(2,rain,10,1)= num_qcstat_conv(2,rain,10,1) + 1
               write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
                  'Rainfall','Rain',iv%info(rain)%lat(1,n),iv%info(rain)%lon(1,n),'-8888.88'   
            end if
         end if
      end if
   end do

   if (trace_use) call da_trace_exit("da_check_max_iv_rain")       

end subroutine da_check_max_iv_rain
subroutine da_get_hr_rain(k, grid, hr_rainc, hr_rainnc, savegridrainc, savegridrainnc)

   !-----------------------------------------------------------------------
   ! Purpose: TBD    
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: k     
   type(domain),     intent(in)    :: grid    ! first guess state.
   real, dimension(ims:ime,jms:jme,1:num_fgat_time), intent(inout) :: hr_rainc, hr_rainnc
   real, dimension(ims:ime,jms:jme), intent(inout) :: savegridrainc, savegridrainnc
   
   if (trace_use) call da_trace_entry("da_get_hr_rain")

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_RAIN.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_EM_RAIN_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   if (k .lt. num_fgat_time) then
      !hr_rainc(:,:,k+1)  = hr_rainc (:,:,k+1) - grid%rainc(:,:) 
      !hr_rainnc(:,:,k+1) = hr_rainnc(:,:,k+1) - grid%rainnc(:,:)
      hr_rainc(:,:,k+1)  = savegridrainc(:,:)  - grid%rainc(:,:) 
      hr_rainnc(:,:,k+1) = savegridrainnc(:,:) - grid%rainnc(:,:)
      savegridrainc(:,:)  = grid%rainc(:,:)
      savegridrainnc(:,:) = grid%rainnc(:,:)
   else
      hr_rainc(:,:,k)  = grid%rainc(:,:) 
      hr_rainnc(:,:,k) = grid%rainnc(:,:)
      savegridrainc(:,:)  = grid%rainc(:,:)
      savegridrainnc(:,:) = grid%rainnc(:,:)
   endif

 if (trace_use) call da_trace_exit("da_get_hr_rain")

end subroutine da_get_hr_rain
subroutine da_get_innov_vector_rain( it, num_qcstat_conv, grid, ob, iv, hr_rainc, hr_rainnc)
   !-----------------------------------------------------------------------
   ! Purpose: TBD    
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it      ! External iteration.
   type(domain),     intent(in)    :: grid    ! first guess state.
   type(y_type),     intent(inout) :: ob      ! Observation structure.
   type(iv_type),    intent(inout) :: iv      ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)
   real, dimension(ims:ime,jms:jme), intent(inout) :: hr_rainc, hr_rainnc

   integer :: n        ! Loop counter.

   real, allocatable :: model_rainc(:)
   real, allocatable :: model_rainnc(:)

   if (trace_use) call da_trace_entry("da_get_innov_vector_rain")

   if ( it > 1 ) then
      do n=iv%info(rain)%n1,iv%info(rain)%n2
            if (iv%rain(n)%rain%qc == fails_error_max) iv%rain(n)%rain%qc = 0
      end do
   end if

   ! [0.0]  Get hourly rainfall

   allocate (model_rainc(iv%info(rain)%n1:iv%info(rain)%n2))
   allocate (model_rainnc(iv%info(rain)%n1:iv%info(rain)%n2))
   model_rainc  = 0.0
   model_rainnc = 0.0

   ! [1.0] horizontal interpolation:

   call da_interp_lin_2d (hr_rainc,  iv%info(rain), 1, model_rainc)
   call da_interp_lin_2d (hr_rainnc, iv%info(rain), 1, model_rainnc)

   do n=iv%info(rain)%n1,iv%info(rain)%n2

      ! [2.0] Initialise components of innovation vector: 

      iv % rain(n) % rain % inv = 0.0

      ! [3.0] To compute the 'inv':

      if (ob % rain(n) % rain > missing_r .and. iv % rain(n) % rain % qc >=  obs_qc_pointer) then
         iv % rain(n) % rain % inv = ob % rain(n) % rain - model_rainc(n) - model_rainnc(n)
      else
         iv % rain(n) % rain % inv = 0.0
      end if	 
     
   end do

   deallocate(model_rainc)
   deallocate(model_rainnc)
 
   ! -----------------------------------------------------------------------
   ! [4.0] Perform optional maximum error check:
   !-----------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_rain(iv,ob, it, num_qcstat_conv)

   if (trace_use) call da_trace_exit("da_get_innov_vector_rain")

end subroutine da_get_innov_vector_rain


subroutine da_calculate_grady_rain(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n
   
   if (trace_use) call da_trace_entry("da_calculate_grady_rain")       

   do n=1, iv%info(rain)%nlocal
             if (iv%rain(n)%rain%qc < obs_qc_pointer) re%rain(n)%rain = 0.0
             jo_grad_y%rain(n)%rain = -re%rain(n)%rain / (iv%rain(n)%rain%error * iv%rain(n)%rain%error)
   end do
   
   if (trace_use) call da_trace_exit("da_calculate_grady_rain")  
     
end subroutine da_calculate_grady_rain




end module da_rain 

