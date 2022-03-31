












module da_statistics
   
   
   
   
   
   use module_domain, only : domain
   use da_control, only : obs_qc_pointer,stdout, missing_r, &
      myproc,rootproc, mjy, mix, mkz, jts,jte,its,ite,kts,kte, trace_use_dull, trace_use,&
      crtm_cloud,use_radar_rf,pptop,ppbot,num_ob_indexes,num_ob_vars,npres_print,&
      obs_names, ob_vars, filename_len 
   use da_define_structures, only : iv_type, maxmin_type, x_type, maxmin_field_type
   use da_par_util1, only : da_proc_sum_real, da_proc_sum_int, da_proc_sum_ints
   use da_par_util, only : da_proc_maxmin_combine
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_tools_serial, only : da_free_unit, da_get_unit
   use da_reporting, only : da_error
   
   implicit none
   
   contains
   
subroutine da_analysis_stats  (grid, stats_unit)
   
   !------------------------------------------------------------------------
   ! Purpose: Calculate min, max, mean and RMS of input 1d field.
   !------------------------------------------------------------------------

   implicit none

   type (domain), intent (in) :: grid
   integer,       intent (in) :: stats_unit ! Output unit for stats.
   
   integer :: i, j, k
   integer :: ij_g, ijk_g   ! ij, ijk for global domain. 
   integer :: kdim  ! k range 

   real    :: um, vm, tm, pm, qm , qcwm, qrnm ! On local domain.
   real    :: rij_g, rijk_g      ! On global domain.

   type (maxmin_field_type) :: max_u(kts:kte), max_v(kts:kte), &
                               max_t(kts:kte), max_p(kts:kte), &
                               max_q(kts:kte), &
                               max_qcw(kts:kte),max_qrn(kts:kte),&
                               min_u(kts:kte), min_v(kts:kte), &
                               min_t(kts:kte), min_p(kts:kte), &
                               min_q(kts:kte), &
                               min_qcw(kts:kte),min_qrn(kts:kte)

 
   real                        uv(kts:kte), vv(kts:kte), &
                               tv(kts:kte), pv(kts:kte), &
                               qv(kts:kte), &
                               qcwv(kts:kte),qrnv(kts:kte)

   call da_trace_entry("da_analysis_stats")

   kdim = kte-kts+1

   ij_g = mix * mjy
   ijk_g = ij_g * mkz

   rij_g  = 1.0/real(ij_g)
   rijk_g = 1.0/real(ijk_g)

   if (rootproc) then
      write(unit=stats_unit, fmt='(/a/)') ' Minimum of gridded analysis increments'
      if (crtm_cloud .or. use_radar_rf) then
         write(unit=stats_unit, fmt='(8a/)') &
              ' Lvl         ', &
              'u     i    j          ', &
              'v     i    j          ', &
              't     i    j          ', &
              'p     i    j          ', &
              'q     i    j          ', &
              'qcw   i    j          ', &
              'qrn   i    j          '
      else
         write(unit=stats_unit, fmt='(6a/)') &
              ' Lvl         ', &
              'u     i    j          ', &
              'v     i    j          ', &
              't     i    j          ', &
              'p     i    j          ', &
              'q     i    j          '
      end if
   end if
   call da_maxmin_in_field(grid%xa%u(its:ite,jts:jte,kts:kte), max_u, min_u)
   call da_proc_maxmin_combine(kdim, max_u, min_u)
   call da_maxmin_in_field(grid%xa%v(its:ite,jts:jte,kts:kte), max_v, min_v)
   call da_proc_maxmin_combine(kdim, max_v, min_v)
   call da_maxmin_in_field(grid%xa%t(its:ite,jts:jte,kts:kte), max_t, min_t)
   call da_proc_maxmin_combine(kdim, max_t, min_t)
   call da_maxmin_in_field(grid%xa%p(its:ite,jts:jte,kts:kte), max_p, min_p)
   call da_proc_maxmin_combine(kdim, max_p, min_p)
   call da_maxmin_in_field(grid%xa%q(its:ite,jts:jte,kts:kte), max_q, min_q)
   call da_proc_maxmin_combine(kdim, max_q, min_q)
   if (crtm_cloud .or. use_radar_rf) then
      call da_maxmin_in_field(grid%xa%qcw(its:ite,jts:jte,kts:kte), max_qcw, min_qcw)
      call da_proc_maxmin_combine(kdim, max_qcw, min_qcw)
      call da_maxmin_in_field(grid%xa%qrn(its:ite,jts:jte,kts:kte), max_qrn, min_qrn)
      call da_proc_maxmin_combine(kdim, max_qrn, min_qrn)
   end if
        

   um = 999999.0
   vm = 999999.0
   tm = 999999.0
   pm = 999999.0
   qm = 999999.0
   if (crtm_cloud .or. use_radar_rf) then
      qcwm = 999999.0
      qrnm = 999999.0
   end if
   do k = kts, kte   
      if (rootproc) then
         if (crtm_cloud .or. use_radar_rf) then
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),3(e12.4,2i5))') k, &
                 min_u(k), min_v(k), min_t(k), min_p(k), min_q(k),min_qcw(k),min_qrn(k)
         else
            if ( abs(min_q(k)%value) < 1.e-30 ) min_q(k)%value = 0.0
            write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),e12.4,2i5)') k, &
                 min_u(k), min_v(k), min_t(k), min_p(k), min_q(k)
         end if
      end if
      um=minval(min_u(:)%value)
      vm=minval(min_v(:)%value)
      tm=minval(min_t(:)%value)
      pm=minval(min_p(:)%value)
      qm=minval(min_q(:)%value)
      if (crtm_cloud .or. use_radar_rf) then
         qcwm=minval(min_qcw(:)%value)
         qrnm=minval(min_qrn(:)%value)
      end if
   end do

   if (rootproc) then
     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),3(e12.4))') ' ALL', &
          um, vm, tm, pm, qm, qcwm, qrnm
     else
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
          um, vm, tm, pm, qm
     end if
  
      write(unit=stats_unit, fmt='(/a/)') &
         ' Maximum of gridded analysis increments'
     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(8a/)') &
           ' Lvl         ', &
           'u     i    j          ', &
           'v     i    j          ', &
           't     i    j          ', &
           'p     i    j          ', &
           'q     i    j          ',&
           'qcw   i    j        ',&
           'qrn   i    j        '
     else
        write(unit=stats_unit, fmt='(6a/)') &
           ' Lvl         ', &
           'u     i    j          ', &
           'v     i    j          ', &
           't     i    j          ', &
           'p     i    j          ', &
           'q     i    j          '
     end if
   end if

   do k = kts, kte
      if (rootproc) then
        if (crtm_cloud .or. use_radar_rf) then
           write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),3(e12.4,2i5))') k, &
             max_u(k), max_v(k), max_t(k), max_p(k), max_q(k),max_qcw(k),max_qrn(k) 
        else
           if ( abs(max_q(k)%value) < 1.e-30 ) max_q(k)%value = 0.0
           write(unit=stats_unit, fmt='(i4,4(f12.4,2i5),e12.4,2i5)') k, &
             max_u(k), max_v(k), max_t(k), max_p(k), max_q(k)
        end if
      end if

      um=maxval(max_u(:)%value)
      vm=maxval(max_v(:)%value)
      tm=maxval(max_t(:)%value)
      pm=maxval(max_p(:)%value)
      qm=maxval(max_q(:)%value)
      if (crtm_cloud .or. use_radar_rf) then
         qcwm=maxval(max_qcw(:)%value)
         qrnm=maxval(max_qrn(:)%value)
      end if
   end do

   if (rootproc) then
     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
           um, vm, tm, pm, qm, qcwm, qrnm
     else
        write(unit=stats_unit, fmt='(a,4(f12.4,10x),e12.4)') ' ALL', &
           um, vm, tm, pm, qm
     end if
      write(unit=stats_unit, fmt='(/a/)') ' Mean of gridded analysis increments'
      
     if (crtm_cloud .or. use_radar_rf) then
        write(unit=stats_unit, fmt='(a/)') &
          ' Lvl        u           v           t           p           q           qcw           qrn'
     else
        write(unit=stats_unit, fmt='(a/)') &
          ' Lvl        u           v           t           p           q'
     end if
   end if

   um = 0.0
   vm = 0.0
   tm = 0.0
   pm = 0.0
   qm = 0.0
   do k = kts, kte
      uv(k) = sum(grid%xa%u(its:ite,jts:jte,k))
      vv(k) = sum(grid%xa%v(its:ite,jts:jte,k))
      tv(k) = sum(grid%xa%t(its:ite,jts:jte,k))
      pv(k) = sum(grid%xa%p(its:ite,jts:jte,k))
      qv(k) = sum(grid%xa%q(its:ite,jts:jte,k))
   end do
   if (crtm_cloud .or. use_radar_rf) then
      qcwm = 0.0
      qrnm = 0.0
       do k = kts, kte
          qcwv(k) = sum(grid%xa%qcw(its:ite,jts:jte,k))
          qrnv(k) = sum(grid%xa%qrn(its:ite,jts:jte,k))
       end do
      call da_proc_sum_real  (qcwv)
      call da_proc_sum_real  (qrnv)
   end if
   call da_proc_sum_real  (uv)
   call da_proc_sum_real  (vv)
   call da_proc_sum_real  (tv)
   call da_proc_sum_real  (pv)
   call da_proc_sum_real  (qv)
   if (rootproc) then
      do k = kts, kte
         write(unit=stats_unit, fmt='(i4,4f12.4,4e12.4)') k, &
            uv(k)*rij_g, vv(k)*rij_g, tv(k)*rij_g, &
            pv(k)*rij_g, qv(k)*rij_g

         um=um+uv(k)
         vm=vm+vv(k)
         tm=tm+tv(k)
         pm=pm+pv(k)
         qm=qm+qv(k)
      end do
   end if

   if (rootproc) then
      write(unit=stats_unit, fmt='(a,4f12.4,4e12.4)') ' ALL', &
         um*rijk_g, vm*rijk_g, tm*rijk_g, pm*rijk_g, qm*rijk_g

      write(unit=stats_unit, fmt='(/a/)') ' RMSE of gridded analysis increments'

      write(unit=stats_unit, fmt='(a/)') &
         ' Lvl        u           v           t           p           q'
   end if

   um = 0.0
   vm = 0.0
   tm = 0.0
   pm = 0.0
   qm = 0.0
   uv = 0.0
   vv = 0.0
   tv = 0.0
   pv = 0.0
   qv = 0.0
   do k = kts, kte
      do j=jts,jte
         do i=its,ite
            uv(k) = uv(k) + grid%xa%u(i,j,k) * grid%xa%u(i,j,k)
            vv(k) = vv(k) + grid%xa%v(i,j,k) * grid%xa%v(i,j,k)
            tv(k) = tv(k) + grid%xa%t(i,j,k) * grid%xa%t(i,j,k)
            pv(k) = pv(k) + grid%xa%p(i,j,k) * grid%xa%p(i,j,k)
            qv(k) = qv(k) + grid%xa%q(i,j,k) * grid%xa%q(i,j,k)
         end do
      end do
   end do

   call da_proc_sum_real  (uv)
   call da_proc_sum_real  (vv)
   call da_proc_sum_real  (tv)
   call da_proc_sum_real  (pv)
   call da_proc_sum_real  (qv)
   if (rootproc) then
      do k = kts, kte
         write(unit=stats_unit, fmt='(i4,4f12.4,4e12.4)') k, &
            sqrt(uv(k)*rij_g), &
            sqrt(vv(k)*rij_g), &
            sqrt(tv(k)*rij_g), &
            sqrt(pv(k)*rij_g), &
            sqrt(qv(k)*rij_g)

         um=um+uv(k)
         vm=vm+vv(k)
         tm=tm+tv(k)
         pm=pm+pv(k)
         qm=qm+qv(k)
      end do
   end if

   if (rootproc) then
      write(unit=stats_unit, fmt='(a,4f12.4,4e12.4)') ' ALL', &
         sqrt(um*rijk_g), sqrt(vm*rijk_g), sqrt(tm*rijk_g), &
         sqrt(pm*rijk_g), sqrt(qm*rijk_g)
   end if

   call da_trace_exit("da_analysis_stats")

contains

subroutine da_maxmin_in_field(field, max, min)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   real,                    intent(in)  :: field(its:ite,jts:jte,kts:kte)
   type(maxmin_field_type), intent(out) :: max(kts:kte)
   type(maxmin_field_type), intent(out) :: min(kts:kte)

   if (trace_use_dull) call da_trace_entry("da_maxmin_in_field")

   do k = kts, kte
      max(k)%value = -1.0e20
      min(k)%value =  1.0e20
      do j = jts, jte
         do i = its, ite
            if (field(i,j,k) > max(k)%value) then
               max(k)%value = field(i,j,k)
               max(k)%i     = i
               max(k)%j     = j
            end if
            if (field(i,j,k) < min(k)%value) then
               min(k)%value = field(i,j,k)
               min(k)%i     = i
               min(k)%j     = j
            end if
         end do
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_maxmin_in_field")

end subroutine da_maxmin_in_field


   
end subroutine da_analysis_stats


subroutine da_correlation_coeff1d(field1, field2, corr_coeff, &
                                   rel_acc)
   
   !---- -------------------------------------------------------------------------
   ! Purpose: Calculate correlation coefficient between two fields.  
   !------------------------------------------------------------------------------
   
   implicit none
   
   real, intent(in)  :: field1(:)       ! input field 1.
   real, intent(in)  :: field2(:)       ! input field 2.
   real, intent(out) :: corr_coeff      ! correlation coefficient
   real, intent(out) :: rel_acc       ! relative error.
   
   integer           :: jx              ! 2nd dimension of field.
   real              :: j_inv           ! 1/(jx)
   real              :: coeff0          ! coefficient.
   real              :: coeff1          ! coefficient.
   real              :: coeff2          ! coefficient.
   real              :: coeff3          ! coefficient.
   real              :: field_mean      ! mean of field.
   real, allocatable :: data1(:)
   real, allocatable :: data2(:)

   if (trace_use) call da_trace_entry("da_correlation_coeff1d")
   
   ! [1.0] Set up scalars:
   
   jx = size(field1(:))
   
   j_inv = 1.0 / real(jx)
   
   ! [2.0] Calculate mean and remove from field:
   
   field_mean = sum(field1(:)) * j_inv
   allocate(data1(1:jx))
   data1(:) = field1(:) - field_mean
   
   field_mean = sum(field2(:)) * j_inv
   allocate(data2(1:jx))
   data2(:) = field2(:) - field_mean
   
   ! [3.0] Calculate correlation coefficient:
   
   coeff0 = sum(data1(:) * data2(:))
   coeff1 = sum(data1(:) * data1(:))
   coeff2 = sum(data2(:) * data2(:))
   
   if (coeff1 /= 0.0 .and. coeff2 /= 0.0) then
      corr_coeff = coeff0 /  sqrt(coeff1 * coeff2)
   else
      corr_coeff = 0.0
   end if
   
   ! [4.0] Calculate accuracy:
   
   coeff3 = sum((data2(:) - data1(:))**2)
   if (coeff2 /= 0.0) then
      rel_acc = 1.0 - min(coeff3/coeff2,1.0)
   else
      rel_acc = 0.0
   end if
   
   deallocate(data1)
   deallocate(data2)

   if (trace_use) call da_trace_exit("da_correlation_coeff1d")
   
end subroutine da_correlation_coeff1d


subroutine da_correlation_coeff2d(field1, field2, corr_coeff, rel_acc)
   
   !--------------------------------------------------------------------------
   ! Purpose: Calculate correlation coefficient between two fields.
   !--------------------------------------------------------------------------

   implicit none
   
   real, intent(in)  :: field1(:,:)     ! input field 1.
   real, intent(in)  :: field2(:,:)     ! input field 2.
   real, intent(out) :: corr_coeff      ! correlation coefficient
   real, intent(out) :: rel_acc      ! relative error.
   
   integer           :: iy              ! size of field dim1.
   integer           :: jx              ! size of field dim2.
   real              :: ij_inv          ! 1/(ij)
   real              :: coeff0          ! coefficient.
   real              :: coeff1          ! coefficient.
   real              :: coeff2          ! coefficient.
   real              :: coeff3          ! coefficient.
   real              :: field_mean      ! mean of field.
   real, allocatable :: data1(:,:)
   real, allocatable :: data2(:,:)

   if (trace_use) call da_trace_entry("da_correlation_coeff2d")
   
   ! [1.0] Set up scalars:
   
   iy = size(field1(:,:),dim=1)
   jx = size(field1(:,:),dim=2)
   
   ij_inv = 1.0 / real(iy*jx)
   
   ! [2.0] Calculate mean and remove from field:
   
   field_mean = sum(field1(:,:)) * ij_inv
   allocate(data1(1:iy,1:jx))
   data1(:,:) = field1(:,:) - field_mean
   
   field_mean = sum(field2(:,:)) * ij_inv
   allocate(data2(1:iy,1:jx))
   data2(:,:) = field2(:,:) - field_mean
   
   ! [3.0] Calculate correlation coefficient:
   
   coeff0 = sum(data1(:,:) * data2(:,:))
   coeff1 = sum(data1(:,:) * data1(:,:))
   coeff2 = sum(data2(:,:) * data2(:,:))
   
   if (coeff1 /= 0.0 .and. coeff2 /= 0.0) then
      corr_coeff = coeff0 /  sqrt(coeff1 * coeff2)
   else
      corr_coeff = 0.0
   end if
   
   ! [4.0] Calculate relative error:
   
   coeff3 = sum((data2(:,:) - data1(:,:))**2)
   if (coeff2 /= 0.0) then
      rel_acc = min(coeff3/coeff2,1.0)
   else
      rel_acc = 0.0
   end if
   
   deallocate(data1)
   deallocate(data2)

   if (trace_use) call da_trace_exit("da_correlation_coeff2d")
   
end subroutine da_correlation_coeff2d


subroutine da_data_distribution(ob_name, num_obs, min_val, max_val, bin_width, ob)
   
   !---------------------------------------------------------------------------
   ! Purpose: Bin ob data to get distribution.
   !---------------------------------------------------------------------------

   implicit none
   
   character (len=*), intent(in) :: ob_name       ! Data description.
   integer,           intent(in) :: num_obs       ! Number of obs.
   real,              intent(in) :: min_val       ! Minimum bin value.
   real,              intent(in) :: max_val       ! Maximum bin value.
   real,              intent(in) :: bin_width     ! Bin width.
   real,              intent(in) :: ob(:)         ! Ob data.
   
   integer              :: num_bins      ! Number of bins
   integer              :: bin           ! Bin counter.
   integer              :: n             ! Data counter.
   integer              :: num_missing   ! Number of missing data.
   
   real, allocatable    :: bin_val(:)    ! Central value of bin.
   real, allocatable    :: bin_min(:)    ! Minimum value of bin.
   integer, allocatable :: num_in_bin(:) ! Number of values in bin. 

   if (trace_use) call da_trace_entry("da_data_distribution")      
   
   !---------------------------------------------------------------------------
   ! [1.0] Initialise bins:
   !---------------------------------------------------------------------------

   num_bins = int((max_val - min_val) / bin_width) + 1

   allocate(bin_val(1:num_bins))
   bin_val(1) = min_val
   do bin = 2, num_bins
      bin_val(bin) = bin_val(bin-1) + bin_width
   end do
   
   allocate(bin_min(1:num_bins+1))
   bin_min(1:num_bins) = bin_val(1:num_bins) - 0.5 * bin_width
   bin_min(num_bins+1) = bin_val(num_bins) + 0.5 * bin_width

   allocate(num_in_bin(0:num_bins+1))
   num_in_bin(0:num_bins+1) = 0
   num_missing = 0
   
   !---------------------------------------------------------------------------
   ! [2.0] Assign data to bins:
   !---------------------------------------------------------------------------
   
   do n = 1, num_obs
      if (ob(n) == missing_r) then
         num_missing = num_missing + 1
      else if (ob(n) < bin_min(1) .AND. ob(n) /= missing_r) then
         num_in_bin(0) = num_in_bin(0) + 1
      else if (ob(n) >= bin_min(num_bins+1)) then
         num_in_bin(num_bins+1) = num_in_bin(num_bins+1) + 1
      else
         do bin = 1, num_bins
            if (ob(n) >= bin_min(bin) .AND. ob(n) < bin_min(bin+1)) then
                 num_in_bin(bin) = num_in_bin(bin) + 1
               exit
            end if
         end do
      end if
   end do

   !---------------------------------------------------------------------------
   ! [3.0] Output statistics:
   !---------------------------------------------------------------------------
   
   write(unit=stdout,fmt='(A,A,A,I8)')' Number of ', trim(ob_name), &
      ' obs = ', num_obs
   write(unit=stdout,fmt='(A,I8)')' Number with missing data indicator = ', &
      num_missing
   write(unit=stdout,fmt='(A,f12.5,A,I8)')' Number below minimum O-B(', &
      min_val-0.5*bin_width, ') = ', num_in_bin(0)
   do bin = 1, num_bins
      write(unit=stdout,fmt='(I4,f12.5,I8)')bin, bin_val(bin), num_in_bin(bin)
   end do
   write(unit=stdout,fmt='(A,f12.5,A,I8)')' Number above maximum O-B(', &
      max_val+0.5*bin_width, ') = ', num_in_bin(num_bins+1)
                               
   !---------------------------------------------------------------------------
   ! [4.0] Tidy up:
   !---------------------------------------------------------------------------

   deallocate(bin_val)
   deallocate(bin_min)
   deallocate(num_in_bin)

   if (trace_use) call da_trace_exit("da_data_distribution")      

end subroutine da_data_distribution


subroutine da_stats_calculate(n, k, qc_flag, x, nn, minimum, maximum, &
   average, rms_err)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,            intent(in)    :: n             ! Sequence number of ob.
   integer,            intent(in)    :: k             ! Level of ob.
   integer,            intent(in)    :: qc_flag       ! QC flag.
   real,               intent(in)    :: x             ! Value.
   integer,            intent(inout) :: nn            ! Number of ok obs.
   type (maxmin_type), intent(inout) :: minimum, maximum
   real,               intent(inout) :: average, rms_err

   if (trace_use_dull) call da_trace_entry("da_stats_calculate")

   if (qc_flag >= obs_qc_pointer) then
      nn = nn + 1

      if (x < minimum%value) then
         minimum%value = x
         minimum%n     = n
         minimum%l     = k
      end if

      if (x > maximum%value) then
         maximum%value = x
         maximum%n     = n
         maximum%l     = k
      end if

      average = average + x
      rms_err = rms_err + x * x
   end if

   if (trace_use_dull) call da_trace_exit("da_stats_calculate")

end subroutine da_stats_calculate


subroutine da_print_qcstat(it, iv, num_qcstat_conv)                       

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv      ! innovation vector.
   integer,           intent(in) :: it
   integer,        intent(inout) ::  num_qcstat_conv(:,:,:,:)

   integer                       :: ii, i, j, k ,n, num, ounit, ios
   character(len=filename_len)   :: file


   logical   :: write_head
   character(1),dimension(115):: cline
   character(4),dimension(2):: typx
   integer,  allocatable :: count(:)
!
   if (trace_use) call da_trace_entry("da_print_qcstat")
   num = num_ob_indexes*num_ob_vars*(npres_print)
   allocate (count(2*num))

  do k=1,115
     cline(k) = '-'
  end do
  typx(1)='used'
  typx(2)='rej '
    count  = 0
    ii = 0    
        do k = 1, npres_print   
        do j = 1, num_ob_vars
        do i = 1, num_ob_indexes
          ii = ii + 1
           count(ii)     = num_qcstat_conv(1,i,j,k)  
           count(num+ii) = num_qcstat_conv(2,i,j,k)  
        end do
        end do
        end do

    call da_proc_sum_ints(count)

  if (rootproc) then
      call da_get_unit(ounit)
      write(unit=file,fmt ='(a,i2.2)')'qcstat_conv_',it
      open(unit=ounit,file=trim(file),form='formatted', status='replace', iostat=ios)
      if (ios /= 0) call da_error("da_print_qcstat.inc",49, &
         (/"Cannot open file "//file/))
    num_qcstat_conv = 0
    ii = 0    
        do k = 1, npres_print    
        do j = 1, num_ob_vars
        do i = 1, num_ob_indexes
          ii = ii + 1
          num_qcstat_conv(1,i,j,k) = count(ii)  
          num_qcstat_conv(2,i,j,k) = count(num+ii)  
        end do
        end do
        end do

        do j = 1, num_ob_vars
        do i = 1, num_ob_indexes
         num_qcstat_conv(1,i,j,npres_print+1) = sum( num_qcstat_conv(1,i,j,1:npres_print) )  
         num_qcstat_conv(2,i,j,npres_print+1) = sum( num_qcstat_conv(2,i,j,1:npres_print) )  
        end do
        end do
      write_head = .false.
   do i = 1, num_ob_indexes
     if (.not. write_head) then
51   format(110a1)
     write(ounit,50)it
50   format(20x,'WRF-Var data utilization statistics for outer iteration ',i3,/)
     write(ounit,510)'ptop',(pptop(k),k=1,npres_print), 0.0
     write(ounit,511)'obs type','var','pbot',(ppbot(k),k=1,npres_print), 2000.0
510  format(15x,a8,1x,13(1x,f6.1))
511  format(1x,a8,1x,a3,6x,a4,1x,13(1x,f6.1))
     write(ounit,500) (cline(j),j=1,115)
500  format(115a1)
     write_head = .true.
     end if
    do j = 1, num_ob_vars
     if( num_qcstat_conv(1,i,j,npres_print+1) > 0 )  then            
      write(ounit,700) obs_names(i),ob_vars(j),typx(1),&
      ((num_qcstat_conv(1,i,j,k) - num_qcstat_conv(2,i,j,k)),k=1,npres_print+1)
      write(ounit,701) typx(2), (num_qcstat_conv(2,i,j,k),k=1,npres_print+1)
700 format(1x,a10,a3,2x,a4,3x, 25(1x,i6) )
701 format(16x,a4,3x,25(1x,i6) )
     end if
    end do
   end do
    write(ounit,500) (cline(j),j=1,115)
  close (ounit) 
  call da_free_unit(ounit)
  end if

   deallocate (count)
   if (trace_use) call da_trace_exit("da_print_qcstat")

end subroutine da_print_qcstat

end module da_statistics

