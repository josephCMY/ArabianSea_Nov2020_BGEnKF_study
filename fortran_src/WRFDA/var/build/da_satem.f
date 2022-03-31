












module da_satem


   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, trace_use_dull, &
      missing, max_error_uv, max_error_t, rootproc, kts,kte, fails_error_max, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv,  &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, satem,max_error_thickness, above_model_lid,&
      ob_vars,qcstat_conv_unit

   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      maxmin_type
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_physics, only : da_tv_profile, da_thickness, da_find_layer, da_thickness_adj, &
      da_find_layer_adj, da_tv_profile_adj, da_find_layer_tl, &
      da_thickness_tl, da_tv_profile_tl
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_satem1_type
      real          :: thickness                
   end type residual_satem1_type

   type maxmin_satem_stats_type
      type (maxmin_type)         :: thickness
   end type maxmin_satem_stats_type

   type stats_satem_type
      type (maxmin_satem_stats_type)  :: maximum, minimum
      type (residual_satem1_type)     :: average, rms_err
   end type stats_satem_type


contains

subroutine da_ao_stats_satem (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type  (y_type), intent (in)    :: re            ! A - O

   type (stats_satem_type) :: stats
   integer                 :: nthickness
   integer                 :: n, k

   if (trace_use_dull) call da_trace_entry("da_ao_stats_satem")

   nthickness = 0

   stats%maximum%thickness = maxmin_type (missing_r, 0, 0)
   stats%minimum%thickness = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_satem1_type(0.0)
   stats%rms_err = stats%average

   nthickness = 0

   stats%maximum%thickness = maxmin_type(0.0, 0, 0)

   stats%minimum = stats%maximum
   stats%average = residual_satem1_type(0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(satem)%nlocal
      if (iv%info(satem)%proc_domain(1,n)) then
         do k=1, iv%info(satem)%levels(n)
            call da_stats_calculate (n, k, iv%satem(n)%thickness(k)%qc, & 
               re%satem(n)%thickness(k), nthickness, &
               stats%minimum%thickness, stats%maximum%thickness, &
               stats%average%thickness, stats%rms_err%thickness)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nthickness)
   iv%nstats(satem) = nthickness

   call da_proc_stats_combine(stats%average%thickness, stats%rms_err%thickness, &
      stats%minimum%thickness%value, stats%maximum%thickness%value, &
      stats%minimum%thickness%n, stats%maximum%thickness%n, &
      stats%minimum%thickness%l, stats%maximum%thickness%l)

   if (rootproc) then
      if (nthickness /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for satem'
         call da_print_stats_satem(stats_unit, nthickness, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_satem")

end subroutine da_ao_stats_satem


subroutine da_jo_and_grady_satem(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector.
   type (y_type), intent(in)    :: re          ! Residual vector.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout):: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_satem")

   jo % satem_thickness = 0.0

   do n=1, iv%info(satem)%nlocal
      do k=1, iv%info(satem)%levels(n)
         jo_grad_y%satem(n)%thickness(k) = -re%satem(n)%thickness(k) / &
            (iv%satem(n)%thickness(k)%error * iv%satem(n)%thickness(k)%error)
      end do

      if (iv%info(satem)%proc_domain(1,n)) then
         do k=1, iv%info(satem)%levels(n)
            jo % satem_thickness = jo % satem_thickness - &
               re%satem(n)%thickness(k) * jo_grad_y%satem(n)%thickness(k)
         end do
      end if
   end do

   jo % satem_thickness = 0.5 * jo % satem_thickness

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_satem")

end subroutine da_jo_and_grady_satem


subroutine da_residual_satem(iv, y, re,np_missing, np_bad_data,np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals for satem obs
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

   if (trace_use_dull) call da_trace_entry("da_residual_satem")

   n_obs_bad % thickness % num = number_type(0, 0, 0)

   do n=1, iv%info(satem)%nlocal
       do k=1, iv%info(satem)%levels(n)
          np_available = np_available + 1
          re%satem(n)%thickness(k) = &
             da_residual(n, k, y%satem(n)%thickness(k), iv%satem(n)%thickness(k), n_obs_bad % thickness)
       end do
   end do

   np_missing  = np_missing  + n_obs_bad % thickness % num % miss
   np_bad_data = np_bad_data + n_obs_bad % thickness % num % bad
   np_obs_used = np_obs_used + n_obs_bad % thickness % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_satem")

end subroutine da_residual_satem


subroutine da_oi_stats_satem (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_satem_type) :: stats
   integer                 :: nthickness
   integer                 :: n, k

   if (trace_use_dull) call da_trace_entry("da_oi_stats_satem")

   nthickness = 0

   stats%maximum%thickness = maxmin_type(missing_r, 0, 0)
   stats%minimum%thickness = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_satem1_type(0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(satem)%nlocal
      if (iv%info(satem)%proc_domain(1,n)) then
         do k=1, iv%info(satem)%levels(n)
            call da_stats_calculate(iv%info(satem)%obs_global_index(n), &
               k, iv%satem(n)%thickness(k)%qc, &
               iv%satem(n)%thickness(k)%inv, nthickness, &
               stats%minimum%thickness, stats%maximum%thickness, &
               stats%average%thickness, stats%rms_err%thickness)
         end do
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nthickness)

   call da_proc_stats_combine(stats%average%thickness, stats%rms_err%thickness, &
      stats%minimum%thickness%value, stats%maximum%thickness%value, &
      stats%minimum%thickness%n, stats%maximum%thickness%n, &
      stats%minimum%thickness%l, stats%maximum%thickness%l)

   if (rootproc) then
      if (nthickness /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for satem'
         call da_print_stats_satem(stats_unit, nthickness, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_satem")

end subroutine da_oi_stats_satem


subroutine da_print_stats_satem(stats_unit, nthickness, satem)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                 intent(in)    :: stats_unit
   integer,                 intent(inout) :: nthickness
   type (stats_satem_type), intent(in)    :: satem

   if (trace_use_dull) call da_trace_entry("da_print_stats_satem")

   write(unit=stats_unit, fmt='(a/)') &
      '   var           thickness(m)  n    k'  

   write(unit=stats_unit, fmt='(a,i16)') &
      '  Number: ', nthickness

   if (nthickness < 1) nthickness = 1

   write(unit=stats_unit, fmt='((a,f12.4,2i5))') &
      ' Minimum(n,k): ', satem%minimum%thickness,    &
      ' Maximum(n,k): ', satem%maximum%thickness
   write(unit=stats_unit, fmt='((a,f12.4,10x))') &
      ' Average     : ', satem%average%thickness/real(nthickness),    &
      '    RMSE     : ', sqrt(satem%rms_err%thickness/real(nthickness))

   if (trace_use_dull) call da_trace_exit("da_print_stats_satem")

end subroutine da_print_stats_satem


subroutine da_transform_xtoy_satem (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer :: n        ! Loop counter.
   integer :: i, j     ! Index dimension.
   real    :: dx, dxm  !
   real    :: dy, dym  !
   integer :: num_levs ! obs vertical levels

   integer :: k
   real    :: pre_ma(kts-1:kte+1)
   real    :: tv_ma(kts-1:kte+1)
   integer :: layer1,layer2
   real    :: tv1,tv2,pres2

   real    :: TGL_pre_ma(kts-1:kte+1)
   real    :: TGL_tv_ma(kts-1:kte+1)
   real    :: TGL_tv1,TGL_tv2

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_satem")

   do n=iv%info(satem)%n1,iv%info(satem)%n2

      num_levs = iv%info(satem)%levels(n)

      ! [1.0] Get horizontal interpolation weights:

      i   = iv%info(satem)%i(1,n)
      dy  = iv%info(satem)%dy(1,n)
      dym = iv%info(satem)%dym(1,n)
      j   = iv%info(satem)%j(1,n)
      dx  = iv%info(satem)%dx(1,n)
      dxm = iv%info(satem)%dxm(1,n)

      ! [2.0] Virtual temperature at obs pt.

      call da_tv_profile_tl(grid,i,j,dx,dxm,dy,dym,                &
         pre_ma,tv_ma,TGL_pre_ma,TGL_tv_ma)

      ! [3.0] Find model vertical position of pressure and do interp.

      call da_find_layer_tl(layer2,tv2,iv%satem(n)%ref_p,              &
         pre_ma,tv_ma,kts,kte,TGL_tv2,TGL_pre_ma,TGL_tv_ma)
      pres2 = iv%satem(n)%ref_p

      ! [4.0] Thickness calculation

      do k=1, num_levs
         if (ABS(iv % satem(n) %p (k) - missing_r) > 1.0) then
            call da_find_layer_tl(layer1,tv1,iv%satem(n)%p(k),            &
               pre_ma,tv_ma,kts,kte,TGL_tv1,TGL_pre_ma,TGL_tv_ma)
            call da_thickness_tl(pre_ma,tv_ma,kts,kte,tv1,tv2,layer1,layer2,&
               iv%satem(n)%p(k),pres2,TGL_pre_ma,TGL_tv_ma,           &
               TGL_tv1,TGL_tv2,y%satem(n)%thickness(k))

            pres2 = iv%satem(n)%p(k)
            layer2 = layer1
            tv2 = tv1
         end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_satem")

end subroutine da_transform_xtoy_satem


subroutine da_transform_xtoy_satem_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.
   integer :: num_levs ! obs vertical levels

   integer :: k
   real    :: pre_ma(kts-1:kte+1)
   real    :: tv_ma(kts-1:kte+1)
   integer :: layer1,layer2
   real    :: tv1,tv2,pres2

   real    :: ADJ_pre_ma(kts-1:kte+1)
   real    :: ADJ_tv_ma(kts-1:kte+1)
   real    :: ADJ_tv1,ADJ_tv2

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_satem_adj")

   ADJ_pre_ma(:) = 0.0
   ADJ_tv_ma(:)  = 0.0
   ADJ_tv1 = 0.0
   ADJ_tv2 = 0.0

   do n=iv%info(satem)%n1,iv%info(satem)%n2
      num_levs = iv%info(satem)%levels(n)

      ! [2.0] Virtual temperature at obs pt.

      call da_tv_profile(grid,iv%info(satem),n,pre_ma,tv_ma)

      ! [3.0] Find model vertical position of pressure and do interp.

      call da_find_layer(layer2,tv2,iv%satem(n)%ref_p,pre_ma,tv_ma,kts,kte)
      pres2 = iv%satem(n)%ref_p

      ! [4.0] Adjoint calculation of satem thickness

      do k=1, num_levs
         if (ABS(iv % satem(n) %p (k) - missing_r) > 1.0) then
            call da_find_layer(layer1,tv1,iv%satem(n)%p(k),pre_ma,tv_ma, &
               kts,kte)

            call da_thickness_adj(pre_ma,tv_ma,kts,kte,tv1,tv2,layer1, &
               layer2, iv%satem(n)%p(k),pres2,ADJ_pre_ma,ADJ_tv_ma, &
               ADJ_tv1,ADJ_tv2,jo_grad_y%satem(n)%thickness(k))

            call da_find_layer_adj(layer1,tv1,iv%satem(n)%p(k),         &
               pre_ma,tv_ma,kts,kte,ADJ_tv1,ADJ_pre_ma,ADJ_tv_ma)

            pres2  = iv%satem(n)%p(k)
            layer2 = layer1
            tv2    = tv1
         end if
      end do

      ! [5.0] Adjoint of layer-finding and vertical interpolation

      call da_find_layer_adj(layer2,tv2,iv%satem(n)%ref_p,              &
         pre_ma,tv_ma,kts,kte,ADJ_tv2,ADJ_pre_ma,ADJ_tv_ma)

      ! [6.0] Adjoint of horizontal interpolation
      call da_tv_profile_adj(grid,jo_grad_x,iv%info(satem),n,         &
         pre_ma,tv_ma,ADJ_pre_ma,ADJ_tv_ma)
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_satem_adj")

end subroutine da_transform_xtoy_satem_adj


subroutine da_check_max_iv_satem(iv, it, num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it       ! External iteration.
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer  :: k,n, ipr
   logical  :: failed

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_satem")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n = iv%info(satem)%n1,iv%info(satem)%n2
      do k = 1, iv%info(satem)%levels(n)
         call da_get_print_lvl(iv%satem(n)%p(k),ipr)
         failed=.false.
         if ( iv%satem(n)%thickness(k)%qc >= obs_qc_pointer ) then 
         ! Thickness
            call da_max_error_qc(it, iv%info(satem), n, iv%satem(n)%thickness(k),& 
                                 max_error_thickness, failed)
            if( iv%info(satem)%proc_domain(k,n) ) then
            num_qcstat_conv(1,satem,9,ipr-1) = num_qcstat_conv(1,satem,9,ipr-1) + 1
            if (failed) then
               num_qcstat_conv(2,satem,9,ipr-1) = num_qcstat_conv(2,satem,9,ipr-1) + 1
               write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
                     'satem',ob_vars(9),iv%info(satem)%lat(k,n),iv%info(satem)%lon(k,n),0.01*iv%satem(n)%p(k)
            end if
            end if
            end if
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_satem")

end subroutine da_check_max_iv_satem
subroutine da_get_innov_vector_satem(it, num_qcstat_conv,grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !
   !  03/10/2008   Y.-R. Guo
   !
   !    To skip the iv % satem(n) % thickness(k) % inv calculation above the
   !    model lid.
   !
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type(domain),  intent(in)    :: grid
   integer,       intent(in)    :: it       ! External iteration.
   type(y_type),  intent(in)    :: ob       ! Observation structure.
   type(iv_type), intent(inout) :: iv       ! O-B structure.
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs ! Number of obs levels.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real    :: model_thickness(1:max_ob_levels) !Model thickness at ob loc

   real    :: pre_ma(kts-1:kte+1)
   real    :: tv_ma(kts-1:kte+1)
   integer :: layer1,layer2,mm
   real    :: tv1,tv2,pres2
   
   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_satem")

   if ( it > 1 ) then
      do n=iv%info(satem)%n1,iv%info(satem)%n2
         do k = 1, iv%info(satem)%levels(n)
            if (iv%satem(n)%thickness(k)%qc == fails_error_max) iv%satem(n)%thickness(k)%qc = 0
         end do
      end do
   end if

   do n=iv%info(satem)%n1,iv%info(satem)%n2
      num_levs = iv%info(satem)%levels(n)

      if (num_levs < 1) cycle

      model_thickness(:) = 0.0

      ! [1.0] Get cross pt. horizontal interpolation weights:

      i   = iv%info(satem)%i(1,n)
      dy  = iv%info(satem)%dy(1,n)
      dym = iv%info(satem)%dym(1,n)
      j   = iv%info(satem)%j(1,n)
      dx  = iv%info(satem)%dx(1,n)
      dxm = iv%info(satem)%dxm(1,n)

      !------------------------------------------------------------------------

      ! [2.0] Calculate vertical profile of virtual temperature at obs pt.

      call da_tv_profile(grid,iv%info(satem),n,pre_ma,tv_ma)

      ! [3.0] Find model vertical position of pressure and do interp.

      call da_find_layer(layer2,tv2,iv%satem(n)%ref_p,pre_ma,tv_ma,kts,kte)
      pres2 = iv%satem(n)%ref_p

      ! [4.0] Thickness innovations calculation

      do k = 1, num_levs
         iv % satem(n) % thickness(k) % inv = 0.0

      if ( iv%satem(n)%p(k) >= iv%ptop) then
         call da_find_layer(layer1,tv1,iv%satem(n)%p(k),pre_ma,tv_ma,kts,kte)
         call da_thickness(pre_ma,tv_ma,kts,kte,tv1,tv2,layer1,layer2,   &
            iv%satem(n)%p(k),pres2,model_thickness(k))

         if (ABS(ob % satem(n) % thickness(k) - missing_r) > 1.0 .and. &
              iv % satem(n) % thickness(k)%qc /= missing_data) then
            iv % satem(n) % thickness(k) % inv = ob % satem(n) % thickness(k) - model_thickness(k)
!                  mm = mm + 1
!                  write(101,'(A,3I6,2x,A, 2F10.3,F10.0,A,F5.0,2(A,F10.3),A,I8)') &
!                    "num, n,k:", mm, n, k, &
!                    "observed, model_thickness, layer = ", &
!                     ob%satem(n)%thickness(k), &
!                     model_thickness(k), 0.01*pres2, " -", &
!                     iv%satem(n)%p(k)*0.01,'hPa  ob_error:',  &
!                     iv%Satem(n)%thickness(k)%error, " inv=", &
!                     iv%Satem(n)%thickness(k)%inv,   " qc =", &
!                     iv%Satem(n)%thickness(k)%qc
         end if

         pres2 = iv%satem(n)%p(k)
         layer2 = layer1
         tv2 = tv1
      else

!   For other type of OBS, such as SOUND, the vertical range check was
!   complted in DA_Interpolation/to_zk.inc, but SATEM never used to_zk.
!   So it is need the check above the model lid (YRG, 02/21/2008):

               iv%Satem(n)%thickness(k)%qc = above_model_lid
      endif

      end do
   end do

   !------------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_satem(iv, it, num_qcstat_conv)        

   !------------------------------------------------------------------------
   ! [6.0] Perform land/ocean check
   !------------------------------------------------------------------------

   do n=iv%info(satem)%n1,iv%info(satem)%n2
      i   = iv%info(satem)%i(1,n)
      j   = iv%info(satem)%j(1,n)
      if (grid%xb%landmask(i,j ) /= 0.0 .or. grid%xb%landmask(i+1,j ) /= 0. .or.  &
          grid%xb%landmask(i,j+1) /= 0.0 .or. grid%xb%landmask(i+1,j+1) /= 0.0) then
         iv % satem(n) % thickness(1) % inv = 0.0
      end if
   end do

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_satem")

end subroutine da_get_innov_vector_satem


subroutine da_calculate_grady_satem(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer  :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_satem")

   do n=1, iv%info(satem)%nlocal
      do k=1, iv%info(satem)%levels(n)
         if (iv%satem(n)%thickness(k)%qc < obs_qc_pointer) then
            re%satem(n)%thickness(k) = 0.0
         end if

         jo_grad_y%satem(n)%thickness(k) = -re%satem(n)%thickness(k) / &
            (iv%satem(n)%thickness(k)%error * iv%satem(n)%thickness(k)%error)
         end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_satem")

end subroutine da_calculate_grady_satem



end module da_satem

