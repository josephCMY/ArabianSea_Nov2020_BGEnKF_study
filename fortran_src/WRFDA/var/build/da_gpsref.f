












module da_gpsref

   use module_domain, only : domain
   use module_dm, only : wrf_dm_sum_real

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, radian, &
      missing, max_error_uv, max_error_t, rootproc,fails_error_max, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv, qcstat_conv_unit, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, ob_vars, &
      max_error_bt, max_error_buv, gpsref,max_error_thickness, &


      pseudo_var, num_pseudo, ims,ime,jms,jme, kms,kme,kts,kte, trace_use_dull, &

      anal_type_verify,fails_error_max,pseudo_err,pseudo_x, pseudo_y, stdout, &
      use_gpsrefobs, gpsref_thinning, pseudo_z,pseudo_val,max_error_ref, pseudo, &
      jts, jte,its,ite, npres_print, pptop
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      maxmin_type, da_allocate_observations
   use da_interpolation, only : da_interp_lin_3d,da_interp_lin_3d_adj, &
      da_to_zk
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type residual_gpsref1_type
      real :: ref                   
      real ::   p                   
      real ::   t                   
      real ::   q                   
   end type residual_gpsref1_type

   type maxmin_gpsref_stats_type
      type (maxmin_type)         :: ref          
   end type maxmin_gpsref_stats_type

   type stats_gpsref_type
      type (maxmin_gpsref_stats_type)  :: maximum, minimum
      type (residual_gpsref1_type)     :: average, rms_err
   end type stats_gpsref_type

contains

subroutine da_ao_stats_gpsref (stats_unit, iv, re)

   !--------------------------------------------------------------------
   ! Purpose: Called by da_minimisation/da_write_diagnostics.inc.
   !--------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type  (y_type), intent(in)    :: re            ! A - O

   type (stats_gpsref_type)         :: stats
   integer                          :: ngpsref
   integer                          :: n, k
   real                             :: o_minus_b, o_minus_a, sigma_o, sigma_b
   real                             :: o_minus_b_0, o_minus_a_0, sigma_o_0, sigma_b_0

   if (trace_use_dull) call da_trace_entry("da_ao_stats_gpsref")

   ngpsref = 0
   o_minus_b_0 = 0.0; o_minus_a_0 = 0.0; sigma_o_0 = 0.0; sigma_b_0 = 0.0
   
   stats%maximum%ref = maxmin_type (missing_r, 0, 0)
   stats%minimum%ref = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_gpsref1_type(0.0,0.0,0.0,0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(gpsref)%nlocal
      if (iv%info(gpsref)%proc_domain(1,n)) then
         do k=1, iv%info(gpsref)%levels(n)
            call da_stats_calculate (n, k, iv%gpsref(n)%ref(k)%qc, & 
               re%gpsref(n)%ref(k), ngpsref, &
               stats%minimum%ref, stats%maximum%ref, &
               stats%average%ref, stats%rms_err%ref)

            if (pseudo_var(1:3) == 'ref' .and. num_pseudo > 0) then
               o_minus_b = iv%GPSRef(n)%ref(k)%inv
               o_minus_a = re%gpsref(n)%ref(k)
               sigma_o   = iv%gpsref(n)%ref(k)%error
            end if
         end do

         if (pseudo_var(1:3) == 'ref' .and. num_pseudo > 0) then
            ! Calculate equivalent sigma_b using
            ! O-A=(O-B)*sigma_o/(sigma_o+sigma_b)

            sigma_b = sqrt ((o_minus_b - o_minus_a) &
               / o_minus_a) * sigma_o
            o_minus_b_0 = wrf_dm_sum_real (o_minus_b)
            o_minus_a_0 = wrf_dm_sum_real (o_minus_a)
            sigma_o_0 = wrf_dm_sum_real (sigma_o)
            sigma_b_0 = wrf_dm_sum_real (sigma_b)
            write (unit=stdout,fmt='(A,F10.2)') &
               'TEST_COVERAGE_da_ao_stats_gpsref:  o_minus_b_0 = ', o_minus_b_0 
            write (unit=stdout,fmt='(A,F10.2)') &
               'TEST_COVERAGE_da_ao_stats_gpsref:  o_minus_a_0 = ', o_minus_a_0
            write (unit=stdout,fmt='(A,F10.2)') & 
               'TEST_COVERAGE_da_ao_stats_gpsref:  sigma_o_0 = ', sigma_o_0
            write (unit=stdout,fmt='(A,F10.2)') &
               'TEST_COVERAGE_da_ao_stats_gpsref:  sigma_b_0 = ', sigma_b_0
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
      end if    ! end if (iv%info(gpsref)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.

   call da_proc_sum_int (ngpsref)
   iv%nstats(gpsref) = ngpsref
    
   call da_proc_stats_combine(stats%average%ref, stats%rms_err%ref, &
      stats%minimum%ref%value, stats%maximum%ref%value, &
      stats%minimum%ref%n, stats%maximum%ref%n, &
      stats%minimum%ref%l, stats%maximum%ref%l)
   
   if (rootproc) then
      if (ngpsref > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for gpsref'
            call da_print_stats_gpsref(stats_unit, ngpsref, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_gpsref")

end subroutine da_ao_stats_gpsref


subroutine da_calculate_grady_gpsref(iv, re, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector         
   !----------------------------------------------------------------------

   implicit none


   type (iv_type), intent(in)     :: iv          ! Innovation vector.
   type (y_type),  intent(inout)  :: re          ! Residual vector.
   type (y_type),  intent(inout)  :: jo_grad_y   ! Grad_y(Jo)

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_gpsref")
   
   do n=1, iv%info(gpsref)%nlocal
      do k=1, iv%info(gpsref)%levels(n)
         if (iv%gpsref(n)%ref(k)%qc < obs_qc_pointer) then
            re%gpsref(n)%ref(k) = 0.0
         end if
         jo_grad_y%gpsref(n)%ref(k) = -re%gpsref(n)%ref(k) / (iv%gpsref(n)%ref(k)%error * iv%gpsref(n)%ref(k)%error)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_gpsref")

end subroutine da_calculate_grady_gpsref


subroutine da_jo_and_grady_gpsref( iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   !----------------------------------------------------------------------
   ! Called by da_minimisation/da_jo_and_grady.inc
   !----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(in)    :: iv          ! Innovation vector.
   type(y_type),  intent(in)    :: re          ! Residual vector.
   type(y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type(jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_gpsref")

   jo % gpsref_ref = 0.0

   do n=1, iv%info(gpsref)%nlocal
      do k=1, iv%info(gpsref)%levels(n)
         jo_grad_y%gpsref(n)%ref(k) = -re%gpsref(n)%ref(k) / &
            ( iv%gpsref(n)%ref(k)%error * iv%gpsref(n)%ref(k)%error)
      end do

      if (iv%info(gpsref)%proc_domain(1,n)) then
         do k=1, iv%info(gpsref)%levels(n)
            jo % gpsref_ref = jo % gpsref_ref - &
                re%gpsref(n)%ref(k) * jo_grad_y%gpsref(n)%ref(k)
         end do
      end if

   end do

   jo % gpsref_ref = 0.5 * jo % gpsref_ref    

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_gpsref")

end subroutine da_jo_and_grady_gpsref


subroutine da_residual_gpsref(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residiual for gpsref obs
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type),intent(inout) :: iv     ! Innovation vector (O-B).
   type (y_type) ,intent(in)    :: y      ! y = H (xa)
   type (y_type) ,intent(inout) :: re     ! Residual vector (O-A).

   integer       ,intent(inout) :: np_available
   integer       ,intent(inout) :: np_obs_used
   integer       ,intent(inout) :: np_missing
   integer       ,intent(inout) :: np_bad_data

   type (bad_data_type) :: n_obs_bad
   integer              :: n, k
!
   real                 :: constant0, g_lat, g_height, weight1, weight2, &
                           gpsref_org

   if (trace_use_dull) call da_trace_entry("da_residual_gpsref")

! Assuming the highest latitude is 90.0 degree:
         constant0= sin(radian * 90.0)

   n_obs_bad % gpsref % num = number_type(0, 0, 0)

   do n=1, iv%info(gpsref)%nlocal
      do k=1, iv%info(gpsref)%levels(n)
         np_available = np_available + 1
!
! Weighted the GPSREF innovation with the latitude: 
         if (iv%gpsref(n)%ref(k)%qc >= obs_qc_pointer ) then

! depend on the height: above 7km, set to 1.0, below 7km, decrease to 0.0:
           g_height = iv%gpsref(n)% h(k)
           weight1 = 1.0 - (7000.0 - g_height) / 7000.0
           if ( weight1 > 1.0 ) weight1 = 1.0
! not depend on height:
           weight1 = 1.0

! depend on the latitude, at 90 degree, weight = 1.0, at 0 degree, weight = 0.0
           g_lat    = iv%info(gpsref)%lat(k,n)
           weight2  = abs(sin(radian * g_lat) / constant0)
! not depend on the latitude:
           weight2   = 1.0

           gpsref_org = iv%gpsref(n)%ref(k)%inv
           iv%gpsref(n)%ref(k)%inv = gpsref_org * weight1 * weight2
         endif
!.............................................................
         re%gpsref(n)%ref(k) = &
            da_residual(n, k, y%gpsref(n)%ref(k), iv%gpsref(n)%ref(k), n_obs_bad%gpsref)
!
         if (iv%gpsref(n)%ref(k)%qc >= obs_qc_pointer ) &
           iv%gpsref(n)%ref(k)%inv = gpsref_org
         
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % gpsref % num % miss
   np_bad_data = np_bad_data + n_obs_bad % gpsref % num % bad
   np_obs_used = np_obs_used + n_obs_bad % gpsref % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_gpsref")

   
end subroutine da_residual_gpsref


subroutine da_oi_stats_gpsref (stats_unit, iv)

   ! -------------------------------------------------------------------
   ! Purpose: TBD
   ! -------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   type (stats_gpsref_type) :: stats
   integer                  :: ngpsref
   integer                  :: n, k

   if (trace_use_dull) call da_trace_entry("da_oi_stats_gpsref")

   ngpsref = 0
   
   stats%maximum%ref = maxmin_type(missing_r, 0, 0)
   stats%minimum%ref = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_gpsref1_type(0.0,0.0,0.0,0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(gpsref)%nlocal
      if (iv%info(gpsref)%proc_domain(1,n)) then
         do k=1, iv%info(gpsref)%levels(n)
            call da_stats_calculate(iv%info(gpsref)%obs_global_index(n), &
               k, iv%gpsref(n)%ref(k)%qc, &
               iv%gpsref(n)%ref(k)%inv, ngpsref, &
               stats%minimum%ref, stats%maximum%ref, &
               stats%average%ref, stats%rms_err%ref)
         end do
      end if
   end do

   ! do inter-processor communication to gather statistics.

   call da_proc_sum_int(ngpsref)
   
   call da_proc_stats_combine(stats%average%ref, stats%rms_err%ref, &
       stats%minimum%ref%value, stats%maximum%ref%value, &
       stats%minimum%ref%n, stats%maximum%ref%n, &
       stats%minimum%ref%l, stats%maximum%ref%l)
   
   if (rootproc .and. (ngpsref > 0)) then
      write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for gpsref'
         call da_print_stats_gpsref(stats_unit, ngpsref, stats)
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_gpsref")

end subroutine da_oi_stats_gpsref


subroutine da_print_stats_gpsref(stats_unit, ngpsref, GPSRef)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                  intent(in)    :: stats_unit
   integer,                  intent(inout) :: ngpsref
   type (stats_gpsref_type), intent(in)    :: gpsref

   if (trace_use_dull) call da_trace_entry("da_print_stats_gpsref")
   
   write (unit=stats_unit, fmt='(a/)') '   var           ref(m)  n    k'  

   write (unit=stats_unit, fmt='(a,i16)') '  Number: ', ngpsref

   if (ngpsref < 1) ngpsref = 1
   
   write(unit=stats_unit, fmt='((a,f12.4,2i5))') &
      ' Minimum(n,k): ', GPSRef%minimum%ref,    &
      ' Maximum(n,k): ', GPSRef%maximum%ref
   write(unit=stats_unit, fmt='((a,f12.4,10x))') &
      ' Average     : ', GPSRef%average%ref/real(ngpsref),    &
      '    RMSE     : ', sqrt(GPSRef%rms_err%ref/real(ngpsref))

   if (trace_use_dull) call da_trace_exit("da_print_stats_gpsref")

end subroutine da_print_stats_gpsref


subroutine da_transform_xtoy_gpsref (grid, iv, y)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer :: n  ! Loop counter.

   real, allocatable :: model_ref(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpsref")

   allocate (model_ref(iv%info(gpsref)%max_lev,iv%info(gpsref)%n1:iv%info(gpsref)%n2))

   call da_interp_lin_3d (grid%xa%ref, iv%info(gpsref), model_ref)

   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
      y%gpsref(n)%ref(:) = model_ref(1:iv%info(gpsref)%levels(n),n)
   end do

   deallocate (model_ref)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpsref")

end subroutine da_transform_xtoy_gpsref


subroutine da_transform_xtoy_gpsref_adj(iv, jo_grad_y, jo_grad_x)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer                       :: n  ! Loop counter.

   real, allocatable :: model_ref(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpsref_adj")

   allocate (model_ref(iv%info(gpsref)%max_lev,iv%info(gpsref)%n1:iv%info(gpsref)%n2))

   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
      model_ref(1:iv%info(gpsref)%levels(n),n) = jo_grad_y%gpsref(n)%ref(:)
   end do

   call da_interp_lin_3d_adj (jo_grad_x%ref, iv%info(gpsref), model_ref)

   deallocate(model_ref)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpsref_adj")

end subroutine da_transform_xtoy_gpsref_adj


subroutine da_check_max_iv_gpsref(iv,it, num_qcstat_conv, opt)        

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !
   !    Added argument: opt = 1, max error check, and count the total num. of obs;
   !                    opt = 2, only count the rejected num. of obs by all of QC
   !                             procedures specified for GPSREF.
   !        
   !                             Shu-Ya Chen and Y.-R. Guo, 12/30/2011  
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it       ! External iteration.
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)  

   integer                           :: k, n, ipr
   logical                           :: failed
   integer, intent(in) :: opt ! 1: only counting, 2: print out  (sychen add)
   integer, parameter :: qc_below = -31, qc_middle = -32, qc_above = -33
   integer, parameter :: qc_step1 = -34, qc_step2  = -35  ! refer to Poli et al. (2009)

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_gpsref")

   !----------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !----------------------------------------------------------------------------
         IF ( opt == 1 ) THEN
   do n = iv%info(gpsref)%n1,iv%info(gpsref)%n2
      do k = 1, iv%info(gpsref)%levels(n)
        call da_get_print_lvl(iv%gpsref(n)%p(k)%inv,ipr)
        if (iv%gpsref(n)%p(k)%inv == missing_r) ipr = 1
        failed=.false.
        if( iv%gpsref(n)%ref(k)%qc >= obs_qc_pointer ) &
        call da_max_error_qc(it, iv%info(gpsref), n, iv%gpsref(n)%ref(k), max_error_ref, failed)  
        if( iv%info(gpsref)%proc_domain(k,n) ) &
                 num_qcstat_conv(1,gpsref,8,ipr) = num_qcstat_conv(1,gpsref,8,ipr) + 1
      end do
   end do
         ENDIF
  
         IF ( opt == 2 ) THEN
   do n = iv%info(gpsref)%n1,iv%info(gpsref)%n2
      do k = 1, iv%info(gpsref)%levels(n)
        call da_get_print_lvl(iv%gpsref(n)%p(k)%inv,ipr)
        if (iv%gpsref(n)%p(k)%inv == missing_r) ipr = 1
        failed=.false.
        if ( ( iv%gpsref(n)%ref(k)%qc == fails_error_max ) .or. &
             ( iv%gpsref(n)%ref(k)%qc == qc_below ) .or. &
             ( iv%gpsref(n)%ref(k)%qc == qc_middle ) .or. &
             ( iv%gpsref(n)%ref(k)%qc == qc_above ) .or. &
             ( iv%gpsref(n)%ref(k)%qc == qc_step1 ) .or. &
             ( iv%gpsref(n)%ref(k)%qc == qc_step2 ) .or. &
             ( iv%gpsref(n)%ref(k)%qc == missing_data ) ) then 
             failed=.true.
        endif
            if(failed) then
            num_qcstat_conv(2,gpsref,8,ipr) = num_qcstat_conv(2,gpsref,8,ipr) + 1
            write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2,I5)')&
             'gpsref',ob_vars(8),iv%info(gpsref)%lat(k,n), &
             iv%info(gpsref)%lon(k,n),0.01*iv%gpsref(n)%p(k)%inv, &
             iv%gpsref(n)%ref(k)%qc 
            end if
      end do
   end do
         ENDIF

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_gpsref")

end subroutine da_check_max_iv_gpsref


subroutine da_get_innov_vector_gpsref(it, num_qcstat_conv, grid, ob, iv)

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
   integer :: i, j, k, kk  ! Index dimension.
   real    :: dx, dxm, dz, dzm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real,allocatable :: model_ref(:,:) !Model gpsref at ob loc
   real    :: v_h(kms:kme), v_p(kms:kme)     ! Model value h at ob

   integer :: Iu_ref, l
   integer,allocatable :: used_lev(:,:) ! for obs. data thinning
                                        ! record the closest level with model
   integer,allocatable :: qc(:)         ! record iv%gpsref(n)%ref(k)%qc
                               ! hor. location.
   real    :: distance_h       ! cal. h-difference between obs and model
   real,allocatable :: min_dis(:)   ! minimum difference
                                               ! hor. location.
! t_iwabuchi 20111216 T for QC 
   real, allocatable :: int_t(:,:,:)     ! for T

   ! For quality control

   real   , parameter :: h_1 = 7000.0, h_2 = 25000.0
   ! Lidia Cucurull values:
   real   , parameter :: pcnt1 = 0.05, pcnt2 = 0.04, pcnt3 = 0.10
   ! testing values:
   ! real   , parameter :: pcnt1 = 0.02, pcnt2 = 0.01, pcnt3 = 0.03
   integer, parameter :: qc_below = -31, qc_middle = -32, qc_above = -33

   integer, parameter :: qc_step1 = -34, qc_step2  = -35  ! refer to Poli et al. (2009)
   integer :: top_level
   real, allocatable :: dndz_obs(:),dndz_mod(:)
   real, allocatable :: dndz2_obs(:),dndz2_mod(:)

   integer :: nn, na, nb, ntotal, nqc0, nqc1, nqc2, nqc3, n_rej
   real    :: percnt
   real    :: height_below(5000)
   character(len=40) :: name_qc(5000)
   character(len=2) :: c_it
   character(len=7) :: c_nt

! t_iwabuchii 20111216 GSI regional QC scheme
   real :: g_height, g_lat, cutoff, stddev     ! QC cutoff GSI
   integer, parameter :: qc_cutoff = -36
   real, allocatable :: model_t(:,:)  ! Model value t at ob location


   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_gpsref")

   n_rej = 0
!
   ! GPS REF Pseudo OBS test:

   if ( (pseudo_var(1:3) == 'ref' .and. num_pseudo > 0) .and. &
        it == 1 ) then

      ! Deallocate:
      if (iv%info(gpsref)%nlocal > 0) then
         do n = 1, iv%info(gpsref)%nlocal
            deallocate(iv % gpsref(n) % ref)
            deallocate(iv % gpsref(n) %  h)
            deallocate(iv % gpsref(n) %  p)
            deallocate(iv % gpsref(n) %  t)
            deallocate(iv % gpsref(n) %  q)
            deallocate(ob % gpsref(n) % ref)
         end do
         deallocate(iv % gpsref)
         deallocate(ob % gpsref)
      end if

      use_gpsrefobs = .true.

      iv%info(gpsref)%nlocal = num_pseudo
      iv%info(gpsref)%plocal(1) = num_pseudo
      iv%info(gpsref)%ntotal   = num_pseudo
      iv%info(gpsref)%max_lev = 1
      iv%info(pseudo)%nlocal = 0

      call da_allocate_observations (iv)
      iv%info(gpsref)%n1 = 1
      iv%info(gpsref)%n2 = 1

      allocate(iv%gpsref(num_pseudo)%ref(1:1))
      allocate(iv%gpsref(num_pseudo)%  h(1:1))
      allocate(iv%gpsref(num_pseudo)%  p(1:1))
      allocate(iv%gpsref(num_pseudo)%  t(1:1))
      allocate(iv%gpsref(num_pseudo)%  q(1:1))
      allocate(ob%gpsref(1:num_pseudo))
      allocate(ob%gpsref(num_pseudo)%ref(1:1))

      write(stdout,'(a,i2)') '==> GPS REF pseudo OBS test: num_pseudo=',num_pseudo

      iv%info(gpsref)%levels(1) = 1

      iv%info(gpsref)%x(:,1) = pseudo_x
      iv%info(gpsref)%y(:,1) = pseudo_y

      iv%info(gpsref)%i(:,1) = int(pseudo_x)
      iv%info(gpsref)%j(:,1) = int(pseudo_y)
      iv % gpsref(1) %  h(1) = pseudo_z

      iv%info(gpsref)%dx(:,1) = pseudo_x-real(iv%info(gpsref)%i(1,1))
      iv%info(gpsref)%dy(:,1) = pseudo_y-real(iv%info(gpsref)%j(1,1))
      iv%info(gpsref)%dxm(:,1)=1.0-iv%info(gpsref)%dx(1,1)
      iv%info(gpsref)%dym(:,1)=1.0-iv%info(gpsref)%dy(1,1)

      iv % gpsref(1) %ref(1) % inv = pseudo_val
      iv % gpsref(1) %ref(1) % error = pseudo_err
      iv % gpsref(1) %ref(1) % qc = 0

      ! Set halo:
      if ((iv%info(gpsref)%i(1,1) < its-1) .or.(iv%info(gpsref)%i(1,1) > ite) .or. &
          (iv%info(gpsref)%j(1,1) < jts-1) .or.(iv%info(gpsref)%j(1,1) > jte)) then
         iv%info(gpsref)%proc_domain(:,1) = .false.
      else
         iv%info(gpsref)%proc_domain(:,1) = .false. 

         if (iv%info(gpsref)%i(1,1) >= its .and. iv%info(gpsref)%i(1,1) <= ite .and. & 
             iv%info(gpsref)%j(1,1) >= jts .and. iv%info(gpsref)%j(1,1) <= jte) then 
            iv%info(gpsref)%proc_domain(:,1) = .true. 
         end if 
      end if

      write(stdout,'(a4,2f15.5)') pseudo_var, pseudo_val, pseudo_err
      write(stdout,'(3f15.5)')    pseudo_x, pseudo_y, pseudo_z
   end if

   if (iv%info(gpsref)%nlocal < 1) then
      if (trace_use_dull) call da_trace_exit("da_get_innov_vector_gpsref")
      return
   end if

   ntotal = 0

   allocate (model_ref(iv%info(gpsref)%max_lev,iv%info(gpsref)%n1:iv%info(gpsref)%n2))

   model_ref(:,:) = 0.0

! t_iwabuchi 20111216 (hmn) allocate model_t for GSI regional QC
   allocate (model_t(iv%info(gpsref)%max_lev,iv%info(gpsref)%n1:iv%info(gpsref)%n2))
   model_t(:,:) = 0.0


   allocate (used_lev(kms:kme,iv%info(gpsref)%n1:iv%info(gpsref)%n2))
   used_lev(:,:) = missing_data 

   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2

      do k=1, iv%info(gpsref)%levels(n)
!sychen        if( iv%gpsref(n)%ref(k)%qc == fails_error_max .and. it > 1 ) &
!sychen            iv%gpsref(n)%ref(k)%qc = 0
            if ( ( iv%gpsref(n)%ref(k)%qc == fails_error_max ) .or. &
                 ( iv%gpsref(n)%ref(k)%qc == qc_below ) .or. &
                 ( iv%gpsref(n)%ref(k)%qc == qc_middle ) .or. &
                 ( iv%gpsref(n)%ref(k)%qc == qc_above ) .or. &
                 ( iv%gpsref(n)%ref(k)%qc == qc_step1 ) .or. &
                 ( iv%gpsref(n)%ref(k)%qc == qc_step2 ) ) then
               if( it > 1 ) iv%gpsref(n)%ref(k)%qc = 0
            endif
      end do

      ! Get cross pt. horizontal interpolation weights:

      i   = iv%info(gpsref)%i(1,n)
      j   = iv%info(gpsref)%j(1,n)
      dx  = iv%info(gpsref)%dx(1,n)
      dy  = iv%info(gpsref)%dy(1,n)
      dxm = iv%info(gpsref)%dxm(1,n)
      dym = iv%info(gpsref)%dym(1,n)

      if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) then

         ! Get the zk from gpsref%h:

         do k=kts,kte
            v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                   + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
!
            v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                   + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
         end do
         do k=1, iv%info(gpsref)%levels(n)
            if (iv%gpsref(n)%h(k) > 0.0) &
               call da_to_zk(iv%gpsref(n)%h(k), v_h, v_interp_h, iv%info(gpsref)%zk(k,n))
            if (iv%info(gpsref)%zk(k,n) < 0.0 .and. .not. anal_type_verify) then
               iv%gpsref(n)%ref(k)%qc = missing_data
            end if
         end do
      else
         iv%info(gpsref)%zk(:,n) = pseudo_z
      end if

!
! To assign the retrieved pressure to GPSREF data (YRG, 06/15/2011)............
   do k=1, iv%info(gpsref)%levels(n)
      kk = int (iv%info(gpsref)%zk(k,n))
      if (kk >= kts .and. kk+1 <= kte) then
        dz  = iv%info(gpsref)%zk(k,n) - real(kk)
        dzm = 1.0 - dz
        iv%gpsref(n)%p(k)%inv = v_p(kk  ) * dzm + v_p(kk+1) * dz 
        ob%gpsref(n)%p(k) = iv%gpsref(n)%p(k)%inv
        iv%gpsref(n)%p(k)%qc = -5
      else
        n_rej = n_rej + 1
      endif
   end do
!..............................................................................

   end do

   call da_convert_zk (iv%info(gpsref))

! t_iwabuchi 20111216 (hmn's update) linear interpolation of log(n) 
   call da_interp_lin_3d (grid%xb%reflog, iv%info(gpsref), model_ref)
   allocate (int_t(ims:ime, jms:jme, kms:kme))
   int_t = grid%xb%t
   call da_interp_lin_3d (int_t, iv%info(gpsref), model_t)
   deallocate(int_t)
! t_iwabuchi END 


   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2

! t_iwabuchi 20111216 compute exp of log (N)
   do k = 1, iv%info(gpsref)%levels(n)
     model_ref(k,n)=exp(model_ref(k,n))
!     model_t(k,n)=exp(model_t(k,n))
   end do
! t_iwabuchi END

      if ( (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) .or. &
           it > 1 ) then
         do k = 1, iv%info(gpsref)%levels(n)
            iv%gpsref(n)%ref(k)%inv = 0.0

            if (ob%gpsref(n)%ref(k) > missing_r .AND. &
                 iv%gpsref(n)%ref(k)%qc >= obs_qc_pointer) then
                 iv%gpsref(n)%ref(k)%inv = ob%gpsref(n)%ref(k) - model_ref(k,n)
            end if
         end do
      else
         ob % gpsref(1)%ref(1) = model_ref(1,n) + iv %gpsref(1)%ref(1)%inv 
      end if
   end do


   ! Quality check 1: Gross error(departure from the background) check 
    if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0) .and. check_max_iv ) &
       call da_check_max_iv_gpsref(iv, it, num_qcstat_conv, 1)

! refer to Poli et al. (2009) ------------------------------------------------
! flag if fit in each of these two qc steps for both of obs. and model
! qc_step1: dN/dz < -50 km^-1
! qc_step2: abs(d^2N/dz^2) > 100 km^-2
!                                  Shu-Ya Chen (2009-07-29)
   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
      if (iv%info(gpsref)%levels(n) <= 2) cycle 
      do k=1,iv%info(gpsref)%levels(n)
         if (model_ref(k,n) > 0.0) top_level=k
      end do
      allocate(dndz_obs(top_level))
      allocate(dndz_mod(top_level))
      allocate(dndz2_obs(top_level))
      allocate(dndz2_mod(top_level))

      ! QC_STEP1

        if (.not. anal_type_verify) then
          if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) then

          ! check for bottom boundary (Forward Difference)
               dndz_obs(1)=(ob%gpsref(n)%ref(2)-ob%gpsref(n)%ref(1))/  &
                           ((iv%gpsref(n)%h(2)-iv%gpsref(n)%h(1))/1000.)
               dndz_mod(1)=(model_ref(2,n)-model_ref(1,n))/  &
                           ((iv%gpsref(n)%h(2)-iv%gpsref(n)%h(1))/1000.)
          ! check for upper boundary (Backward Difference)
               dndz_obs(top_level)= &
                       (ob%gpsref(n)%ref(top_level)-ob%gpsref(n)%ref(top_level-1))/  &
                       ((iv%gpsref(n)%h(top_level)-iv%gpsref(n)%h(top_level-1))/1000.)
               dndz_mod(top_level)= &
                       (model_ref(top_level,n)-model_ref(top_level-1,n))/  &
                       ((iv%gpsref(n)%h(top_level)-iv%gpsref(n)%h(top_level-1))/1000.)
          ! check for middle levels (Central Difference)
            do k=2, top_level-1
               dndz_obs(k)=(ob%gpsref(n)%ref(k+1)-ob%gpsref(n)%ref(k-1))/  &
                           ((iv%gpsref(n)%h(k+1)-iv%gpsref(n)%h(k-1))/1000.)
               dndz_mod(k)=(model_ref(k+1,n)-model_ref(k-1,n))/  &
                           ((iv%gpsref(n)%h(k+1)-iv%gpsref(n)%h(k-1))/1000.)
            end do
            do k=1, top_level
! hmn 20111206
               if (iv%gpsref(n)%ref(k)%qc /= missing_data) then
               if ((dndz_obs(k) < -50.) .or. (dndz_mod(k) < -50.)) then
               iv%gpsref(n)%ref(k)%qc = qc_step1
               end if
! hmn 20111206 
               end if
            end do

      ! QC_STEP2

          ! check for bottom boundary
               dndz2_obs(1)=(dndz_obs(2)-dndz_obs(1))/  &
                            ((iv%gpsref(n)%h(2)-iv%gpsref(n)%h(1))/1000.)
               dndz2_mod(1)=(dndz_mod(2)-dndz_mod(1))/   &
                            ((iv%gpsref(n)%h(2)-iv%gpsref(n)%h(1))/1000.)
          ! check for upper boundary
               dndz2_obs(top_level)=(dndz_obs(top_level)-dndz_obs(top_level-1))/  &
                            ((iv%gpsref(n)%h(top_level)-iv%gpsref(n)%h(top_level-1))/1000.)
               dndz2_mod(top_level)=(dndz_mod(top_level)-dndz_mod(top_level-1))/   &
                            ((iv%gpsref(n)%h(top_level)-iv%gpsref(n)%h(top_level-1))/1000.)
          ! check for middle levels
            do k=2, top_level-1
               dndz2_obs(k)=(dndz_obs(k+1)-dndz_obs(k-1))/  &
                            ((iv%gpsref(n)%h(k+1)-iv%gpsref(n)%h(k-1))/1000.)
               dndz2_mod(k)=(dndz_mod(k+1)-dndz_mod(k-1))/   &
                            ((iv%gpsref(n)%h(k+1)-iv%gpsref(n)%h(k-1))/1000.)
            end do
            do k=1, top_level
! hmn 20111206
               if (iv%gpsref(n)%ref(k)%qc /= missing_data) then
               if ((abs(dndz2_obs(k)) > 100.) .or. (abs(dndz2_mod(k)) > 100.)) then
               iv%gpsref(n)%ref(k)%qc = qc_step2
               end if
               end if

            end do
          end if ! end of if pseudo check

        end if  ! end of if verify check
   deallocate(dndz_obs,dndz_mod)
   deallocate(dndz2_obs,dndz2_mod)
   end do  ! end of do iv%info(gpsref)%n1~n2
!
! End of Poli's check. (2009) -------------------------------------------------
! 


! t_iwabuchi 20111216 GSI's regional QC
! START GSI-QC REGIONAL TI-GSI -------------------------------------------------
! ------------------------------------------------------------------------------
! GSI-QC Implementation, Ted Iwabuchi
!   First release: 2011-03-04
!   Review 3.3.1 : 2011-12-16
!
! HISTORY:
!   2011-05-04 : bug in formula fixed
!   2011-05-11 : add abs, and validated
!   2011-06-06 : reviewd, added comments
!   2011-12-16 : implemented for 3.3.1
!
! SUMMARY:
! Modified version of Cucurull, 2010 for NCEP GSI-Global DA
! NCEP-GSI regional DA  setupref.f90
!
!   O-B cutoff depending on assigned std dev for each height and latitude
!
!  > 30                                           all observation is rejected
!  11 - 30 km    0.25 + 0.5 cos (lambda)   [c1]   3
!   9 - 11 km    (11-h)/2 x c2 + (h-9)/2 x c1     3
!   6 -  9 km    0.5 (if T<= 240)          [c2]   3
!                0.01 x T^2 - 0.455 T + 52.075 (if T>240)
!   4 -  6 km    (6-h)/2 x c3 + (h-4)/2 x c2    3
!   0 -  4 km    1+2.5 cos (lambda)        [c3]   1
!
!   h, lambda, T
! ------------------------------------------------------------------------------
!
! REFERENCES:
!   GSI code: comGSI_v3Beta/src/main/setupref.f90
!
!
! NOTES:
!  See TI-GSI blocks, and don't forget to declare type of variables
!  used in this block:
!
! 1  real :: g_height, g_lat, cutoff, stddev     ! QC cutoff GSI
! 2  real, allocatable :: model_t(:,:)  ! Model value t at ob location
! 3  integer, parameter :: qc_cutoff = -36
!

!  process refracticity data n from n1 to n2
!
if (.not. anal_type_verify) then
   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
       do k=1,iv%info(gpsref)%levels(n)
         if (model_ref(k,n) > 0.0) top_level=k

         g_height = iv%gpsref(n)% h(k)
         g_lat    = cos ( radian * iv%info(gpsref)%lat(k,n) )

         if (     g_height >= 0.0     .and. g_height < 4000.0) then
           cutoff = 1.0
           stddev = 1.0 + 2.5 * g_lat

         else if (g_height >= 4000.0  .and. g_height < 6000.0) then
           cutoff = 3.0
           if (model_t(k,n) <= 240.0) then
             stddev = 0.5
           else
             stddev = 0.001 * model_t(k,n) * model_t(k,n) - 0.455 * model_t(k,n) + 52.075
           endif
           stddev = (6000.0 - g_height)/2000.0 * (1.0 + 2.5 * g_lat)  &
                    + (g_height - 4000.0)/2000.0 * stddev

         else if (g_height >= 6000.0  .and. g_height < 9000.0) then
           cutoff = 3.0
           if (model_t(k,n) <= 240.0) then
             stddev = 0.5
           else
!            0.01 x T^2 - 0.455 T + 52.075 (if T>240)
             stddev = 0.001 * model_t(k,n) * model_t(k,n) - 0.455 * model_t(k,n) + 52.075
           end if

         else if (g_height >= 9000.0  .and. g_height < 11000.0) then
           cutoff = 3.0
           if (model_t(k,n) <= 240.0) then
             stddev = 0.5
           else
             stddev = 0.001 * model_t(k,n) * model_t(k,n) - 0.455 * model_t(k,n) + 52.075
           endif
           stddev = (11000.0 - g_height)/2000.0 * stddev  &
                   + (g_height - 9000)/2000.0 * (0.25 + 0.5 * g_lat)


         else if (g_height >= 11000.0 .and. g_height < 30000.0) then
           cutoff = 3.0
           stddev = 0.25 + 0.5 * g_lat

         else if (g_height >= 30000.0 ) then
           cutoff = 0.0

         endif
! hmn 20111206
               if (iv%gpsref(n)%ref(k)%qc /= missing_data) then
! Check innovation, stddev, and cutoff
         if (abs(iv%gpsref(n)%ref(k)%inv) >= stddev * cutoff) then
           iv%gpsref(n)%ref(k)%qc = qc_cutoff
! hmn 20111202
!           qc_36(k,n)=1
! qc flag
         endif
! hmn 20111206
               end if

       end do
   end do
end if

   deallocate (model_t)
!
! End   GSI-QC REGIONAL TI-GSI ------------------------------------------------
! 
! t_iwabuchi END 




   do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
      ! Quality check 2: Error percentage check.

      if (.not. anal_type_verify) then
         if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) then
            do k=1, iv%info(gpsref)%levels(n)

               ! incremetal refractivity or the relative error:
               !   abs[(O-B)/{(O+B)/2}]              (Lidia Cucurull 2005)

               ntotal = ntotal + 1
               percnt = 2.0 * abs(iv%gpsref(n)%ref(k)%inv / &
                 (ob%gpsref(n)%ref(k) + model_ref(k,n)))

! hmn 20111206
!              if (iv%gpsref(n)%ref(k)%qc >= obs_qc_pointer) then
               if (iv%gpsref(n)%ref(k)%qc /= missing_data) then

                  if (iv%gpsref(n)%h(k) < h_1) then
                     if (percnt > pcnt1) iv%gpsref(n)%ref(k)%qc = qc_below
                  else if (iv%gpsref(n)%h(k) > h_2) then
                     if (percnt > pcnt3) iv%gpsref(n)%ref(k)%qc = qc_above
                  else
                     if (percnt > pcnt2) iv%gpsref(n)%ref(k)%qc = qc_middle
                  end if
               end if
            end do
         end if
      end if  ! end of if verify check
   end do

   ! Quality check 3: Low levels quality control

   if (.not. anal_type_verify) then
      if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) then
         ! Search for the GPS RO's name with the 'qc_below':

       if ( maxval(iv%info(gpsref)%levels(:)) > 1 ) then ! gpsref in profiles
         nn = 0
         height_below = 0.0
         name_qc      = '                                       '

         do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
             nn = nn + 1
             iv%info(gpsref)%levels(n) = iv%info(gpsref)%levels(n)
             do k=1, iv%info(gpsref)%levels(n)
                if (iv%gpsref(n)%ref(k)%qc == qc_below) then
                   name_qc(nn) = iv%info(gpsref)%name(n)
                   height_below(nn) = max(iv%gpsref(n)%h(k),height_below(nn))
                end if
             end do
             if (height_below(nn) == 0.0) nn = nn - 1
         end do

         ! Set the flag qc_below to the levels below percnt < pcnt1::

         ntotal = 0
         nqc0   = 0
         nqc1   = 0
         nqc2   = 0
         nqc3   = 0

         do n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
            do na = 1,nn
               if (iv%info(gpsref)%name(n) == name_qc(na)) then
                  do k=1, iv%info(gpsref)%levels(n)
! hmn 20111202
                     if (iv%gpsref(n)%h(k) < height_below(na) .and. &
!                        iv%gpsref(n)%ref(k)%qc >= 0) then
! hmn 20111206
                         iv%gpsref(n)%ref(k)%qc /= missing_data) then
                       iv%gpsref(n)%ref(k)%qc = qc_below
                     end if

                  end do
                  exit
               end if
            end do

            do k=1, iv%info(gpsref)%levels(n)
               ntotal = ntotal + 1
               if (iv%gpsref(n)%ref(k)%qc == fails_error_max) nqc0 = nqc0 + 1
               if (iv%gpsref(n)%ref(k)%qc == qc_middle) nqc1 = nqc1 + 1
               if (iv%gpsref(n)%ref(k)%qc == qc_below) nqc2 = nqc2 + 1
               if (iv%gpsref(n)%ref(k)%qc == qc_above) nqc3 = nqc3 + 1
            end do
         end do
       else    ! gpsref not in profiles
          do na = iv%info(gpsref)%n1, iv%info(gpsref)%n2
             if ( iv%gpsref(na)%ref(1)%qc == qc_below) then
                do nb = iv%info(gpsref)%n1, iv%info(gpsref)%n2
                   if ( iv%info(gpsref)%id(nb) == iv%info(gpsref)%id(na) .and. &
                        iv%info(gpsref)%name(nb) == iv%info(gpsref)%name(na) .and. &
! hmn 20111206
!                        iv%gpsref(nb)%ref(1)%qc >= obs_qc_pointer .and.        &
                        iv%gpsref(nb)%ref(1)%qc /= missing_data .and.            &

                        iv%gpsref(nb)%h(1) <= iv%gpsref(na)%h(1)   ) then
                      iv%gpsref(nb)%ref(1)%qc = qc_below
                   end if
                end do
             end if
          end do
       end if  ! end of if gpsref profiles
      end if
   end if  ! end of if verify check
12221 continue

! print out the amounts of obs. rejection   

    if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0) .and. check_max_iv ) &
       call da_check_max_iv_gpsref(iv, it, num_qcstat_conv, 2)

! ------------------------------------------------------------------------------
!   GPSRO thinning  (Shu-Ya Chen 20090701)
   if (.not. anal_type_verify) then
   if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) then
    IF ( gpsref_thinning ) THEN
       DO n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
       allocate(min_dis(kms:kme))
       allocate(qc(iv%info(gpsref)%levels(n)))
       i   = iv%info(gpsref)%i(1,n)
       j   = iv%info(gpsref)%j(1,n)
       dx  = iv%info(gpsref)%dx(1,n)
       dy  = iv%info(gpsref)%dy(1,n)
       dxm = iv%info(gpsref)%dxm(1,n)
       dym = iv%info(gpsref)%dym(1,n)

       if (.not.(pseudo_var(1:3) == 'ref' .and. num_pseudo > 0)) then
       ! Get the zk from gpsref%h:
          do k=kts,kte
             v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                    + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
          end do
          do k=kts,kte 
          min_dis(k)=1.0E10
             do l=1, iv%info(gpsref)%levels(n)
                if ( iv%gpsref(n)%ref(l)%qc >= obs_qc_pointer ) then
                distance_h=abs(iv%gpsref(n)%h(l)-v_h(k))
                min_dis(k)=min(min_dis(k),distance_h)
                if ( min_dis(k) == distance_h ) used_lev(k,n)=l
                end if
             end do
          end do

          write(533,*) 'obs_qc_pointer=',obs_qc_pointer,'missing_data=',missing_data
          do k=kts,kte
          write(533,*) n,k,'used_lev=',used_lev(k,n)
          enddo

          do l=1, iv%info(gpsref)%levels(n)
          write(533,*) n,l,iv%gpsref(n)%ref(l)%qc
          enddo

          do l=1, iv%info(gpsref)%levels(n)
          qc(l)=iv%gpsref(n)%ref(l)%qc
          end do
          do k=kts,kte
           qc(used_lev(k,n))=1   ! which level is closest to model level
          end do
       !  data thinning (set thinned levels to be -99)
          do l=1, iv%info(gpsref)%levels(n)
          if ( iv%gpsref(n)%ref(l)%qc >= obs_qc_pointer &
               .and. qc(l) /= 1 ) then
          iv%gpsref(n)%ref(l)%qc = -99
          end if
          end do
       end if
       deallocate(min_dis)
       deallocate(qc)
       END DO
    END IF

    goto 12345

! Write out GPS Ref data:

     DO n=iv%info(gpsref)%n1,iv%info(gpsref)%n2
     Iu_ref = 336
     write(c_it,'(I2.2)') it
     c_nt=iv%info(gpsref)%name(n)(8:11)//iv%info(gpsref)%name(n)(28:30)
     open (unit=Iu_ref, file='RO_Innov_'//iv%info(gpsref)%date_char(n)//'_'//c_nt//'.'//c_it, &
           form='formatted')
           write(unit=Iu_ref, fmt='(/i5,2x,a,2x,a,2x,4f10.3,i5)') n, &
              iv%info(gpsref)%date_char(n), iv%info(gpsref)%id(n),  &
              iv%info(gpsref)%lat(1,n)  , iv%info(gpsref)%lon(1,n), &
              iv%info(gpsref)%x(1,n)  , iv%info(gpsref)%y(1,n), &
              iv%info(gpsref)%levels(n)
           write(unit=Iu_ref, fmt='(a5,3x,6a14)') 'level','     height   ', &
                       '    Obs_ref   ','  model_ref   ','  Innov_ref   ', &
                       '  error_ref   ',' qc_ref       '
           do k = 1, iv%info(gpsref)%levels(n)
!             if ( gpsref_thinning ) then
!               if ( iv%gpsref(n)%ref(l)%qc >= obs_qc_pointer ) then
!                  write(unit=Iu_ref, fmt='(i3,1x,5f14.3,i10)')  k, &
!                  iv%gpsref(n)%h(k), ob%gpsref(n)%ref(k),   &
!                  model_ref(k,n), iv%gpsref(n)%ref(k)%inv,   &
!                  iv%gpsref(n)%ref(k)%error, iv%gpsref(n)%ref(k)%qc
!               end if
!             else
               write(unit=Iu_ref, fmt='(i3,1x,5f14.3,i10)')  k, &
                  iv%gpsref(n)%h(k), ob%gpsref(n)%ref(k),   &
                  model_ref(k,n), iv%gpsref(n)%ref(k)%inv,   &
                  iv%gpsref(n)%ref(k)%error, iv%gpsref(n)%ref(k)%qc
!             end if
           end do
     close(Iu_ref)
     END DO
12345 continue
   ! .........................................................................
   end if  ! end of pseudo test
   end if  ! end of verify check

   deallocate (used_lev)
   deallocate (model_ref)

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_gpsref")

end subroutine da_get_innov_vector_gpsref

end module da_gpsref

