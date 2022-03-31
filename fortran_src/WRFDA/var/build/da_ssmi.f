












module da_ssmi


   use module_domain, only : xb_type, domain
   use module_ssmi, only : cal_sigma_v,tb,effht,epsalt,spemiss,tbatmos,&
      roughem,effang,tbatmos,filter

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, t_roughem, t_kelvin, &
      missing, max_error_uv, max_error_t, rootproc, pi, trace_use, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv,  &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, max_error_thickness, mkz, &
      max_error_rh,max_error_tb, max_error_pw, trace_use_dull, &
      test_transforms,stdout, use_ssmiretrievalobs, use_ssmitbobs, &
      global, print_detail_obs, max_ssmi_rv_input, max_ssmi_tb_input, &
      its,ite,jts,jte,kts,kte,kms,kme,ids,ide,jds,jde,fails_error_max, &
      ssmi_tb, ssmi_rv, num_ob_indexes, ssmt1, ssmt2, ob_vars,qcstat_conv_unit
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      maxmin_type,residual_ssmi_rv_type, &
      residual_ssmi_tb_type, model_loc_type, info_type, field_type, &
      count_obs_number_type
   use da_interpolation, only : da_to_zk,da_interp_lin_2d, da_interp_lin_2d_adj, &
      da_interp_lin_3d,da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_reporting, only : da_warning, message, da_error
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_llxy, da_convert_zk,da_get_print_lvl
   use da_tools_serial, only : da_get_unit, da_free_unit
   use da_tracing, only : da_trace_entry, da_trace_exit

   

   type maxmin_ssmi_rv_stats_type
      type (maxmin_type)         :: tpw      
      type (maxmin_type)         :: Speed    
   end type maxmin_ssmi_rv_stats_type

   type stats_ssmi_rv_type
      type (maxmin_ssmi_rv_stats_type)      :: maximum, minimum
      type (residual_ssmi_rv_type)   :: average, rms_err
   end type stats_ssmi_rv_type

   

   type maxmin_ssmi_tb_stats_type
      type (maxmin_type)         :: tb19v    
      type (maxmin_type)         :: tb19h    
      type (maxmin_type)         :: tb22v    
      type (maxmin_type)         :: tb37v    
      type (maxmin_type)         :: tb37h    
      type (maxmin_type)         :: tb85v    
      type (maxmin_type)         :: tb85h    
   end type maxmin_ssmi_tb_stats_type

   type stats_ssmi_tb_type
      type (maxmin_ssmi_tb_stats_type)  :: maximum, minimum
      type (residual_ssmi_tb_type)      :: average, rms_err
   end type stats_ssmi_tb_type


contains

subroutine da_ao_stats_ssmi_rv (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type  (y_type), intent (in)    :: re            ! A - O

   integer                          :: ntpw, nspeed, n
   type (stats_ssmi_rv_type) :: stats

   if (trace_use) call da_trace_entry("da_ao_stats_ssmi_rv")

   ntpw = 0
   nspeed = 0

   stats%maximum%tpw   = maxmin_type (missing_r, 0, 0)
   stats%maximum%speed = maxmin_type (missing_r, 0, 0)
   stats%minimum%tpw   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%speed = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_ssmi_rv_type(0.0, 0.0)
   stats%rms_err = stats%average

   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n=1, iv%info(ssmi_rv)%nlocal
         if (iv%info(ssmi_rv)%proc_domain(1,n)) then
            call da_stats_calculate (n, 0, iv%ssmi_rv(n)%tpw%qc, &
               re%ssmi_rv(n)%tpw, ntpw,  &
               stats%minimum%tpw,  stats%maximum%tpw,&
               stats%average%tpw,  stats%rms_err%tpw)

            call da_stats_calculate (n, 0, iv%ssmi_rv(n)%speed%qc, &
               re%ssmi_rv(n)%speed,    nspeed,  &
               stats%minimum%speed,  stats%maximum%speed,&
               stats%average%speed,  stats%rms_err%speed)
         end if
      end do   
   end if

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (ntpw)
   call da_proc_sum_int (nspeed)
   iv%nstats(ssmi_rv) = ntpw + nspeed
   
   call da_proc_stats_combine(stats%average%tpw, stats%rms_err%tpw, &
      stats%minimum%tpw%value, stats%maximum%tpw%value, &
      stats%minimum%tpw%n, stats%maximum%tpw%n, &
      stats%minimum%tpw%l, stats%maximum%tpw%l)

   call da_proc_stats_combine(stats%average%speed, stats%rms_err%speed, &
      stats%minimum%speed%value, stats%maximum%speed%value, &
      stats%minimum%speed%n, stats%maximum%speed%n, &
      stats%minimum%speed%l, stats%maximum%speed%l)

   if (rootproc) then
      if (ntpw /= 0) then

         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_retrieval'
         write(unit=stats_unit, fmt='(a/)') '   var           tpw(cm)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntpw
         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tpw%value, stats%minimum%tpw%n    , &
            ' Maximum(n): ', stats%maximum%tpw%value, stats%maximum%tpw%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tpw/real(ntpw), &
            '    RMSE   : ', sqrt(stats%rms_err%tpw/real(ntpw))
      end if

      if (nspeed /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_retrieval'
         write(unit=stats_unit, fmt='(a/)') '   var           Speed(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', nspeed
         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%speed%value, stats%minimum%speed%n    , &
            ' Maximum(n): ', stats%maximum%speed%value, stats%maximum%speed%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%Speed/real(nspeed), &
            '    RMSE   : ', sqrt(stats%rms_err%Speed/real(nspeed))
      end if
   end if

   if (trace_use) call da_trace_exit("da_ao_stats_ssmi_rv")

end subroutine da_ao_stats_ssmi_rv


subroutine da_ao_stats_ssmi_tb  (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type  (y_type), intent (in)    :: re            ! A - O

   integer                          :: ntb19v,ntb19h,ntb22v,ntb37v,ntb37h, &
                                       ntb85v,ntb85h
   integer                          :: n
   type (stats_ssmi_tb_type)        :: stats

   if (trace_use) call da_trace_entry("da_ao_stats_ssmi_tb")

   ntb19v = 0
   ntb19h = 0
   ntb22v = 0
   ntb37v = 0
   ntb37h = 0
   ntb85v = 0
   ntb85h = 0

   stats%maximum%tb19v = maxmin_type (missing_r, 0, 0)
   stats%maximum%tb19h = maxmin_type (missing_r, 0, 0)
   stats%maximum%tb22v = maxmin_type (missing_r, 0, 0)
   stats%maximum%tb37v = maxmin_type (missing_r, 0, 0)
   stats%maximum%tb37h = maxmin_type (missing_r, 0, 0)
   stats%maximum%tb85v = maxmin_type (missing_r, 0, 0)
   stats%maximum%tb85h = maxmin_type (missing_r, 0, 0)
   stats%minimum%tb19v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb19h = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb22v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb37v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb37h = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb85v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb85h = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_ssmi_tb_type(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   if (iv%info(ssmi_tb)%nlocal .gt. 0) then
      do n=1, iv%info(ssmi_tb)%nlocal 
         if (iv%info(ssmi_tb)%proc_domain(1,n)) then
            call da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb19v%qc, &
               re%ssmi_tb(n)%tb19v, ntb19v,  &
               stats%minimum%tb19v,  stats%maximum%tb19v,&
               stats%average%tb19v,  stats%rms_err%tb19v)

            call da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb19h%qc, &
               re%ssmi_tb(n)%tb19h, ntb19h,  &
               stats%minimum%tb19h,  stats%maximum%tb19h,&
               stats%average%tb19h,  stats%rms_err%tb19h)

            call da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb22v%qc, &
               re%ssmi_tb(n)%tb22v, ntb22v,  &
               stats%minimum%tb22v,  stats%maximum%tb22v,&
               stats%average%tb22v,  stats%rms_err%tb22v)

            call da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb37v%qc, &
               re%ssmi_tb(n)%tb37v, ntb37v,  &
               stats%minimum%tb37v,  stats%maximum%tb37v,&
               stats%average%tb37v,  stats%rms_err%tb37v)

            call da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb37h%qc, &
               re%ssmi_tb(n)%tb37h, ntb37h,  &
               stats%minimum%tb37h,  stats%maximum%tb37h,&
               stats%average%tb37h,  stats%rms_err%tb37h)

            call da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb85v%qc, &
               re%ssmi_tb(n)%tb85v, ntb85v,  &
               stats%minimum%tb85v,  stats%maximum%tb85v,&
               stats%average%tb85v,  stats%rms_err%tb85v)

            call  da_stats_calculate (n, 0, iv%ssmi_tb(n)%tb85h%qc, &
               re%ssmi_tb(n)%tb85h, ntb85h,  &
               stats%minimum%tb85h,  stats%maximum%tb85h,&
               stats%average%tb85h,  stats%rms_err%tb85h)
         end if    ! end if (iv%info(ssmi_tb)%proc_domain(1,n))
      end do
   end if

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (ntb19v)
   call da_proc_sum_int (ntb19h)
   call da_proc_sum_int (ntb22v)
   call da_proc_sum_int (ntb37v)
   call da_proc_sum_int (ntb37h)
   call da_proc_sum_int (ntb85v)
   call da_proc_sum_int (ntb85h)
   iv%nstats(ssmi_tb) = ntb19v + ntb19h + ntb22v + ntb37v + ntb37h + ntb85v + ntb85h

   call da_proc_stats_combine(stats%average%tb19v, stats%rms_err%tb19v, &
      stats%minimum%tb19v%value, stats%maximum%tb19v%value, &
      stats%minimum%tb19v%n, stats%maximum%tb19v%n, &
      stats%minimum%tb19v%l, stats%maximum%tb19v%l)

   call da_proc_stats_combine(stats%average%tb19h, stats%rms_err%tb19h, &
      stats%minimum%tb19h%value, stats%maximum%tb19h%value, &
      stats%minimum%tb19h%n, stats%maximum%tb19h%n, &
      stats%minimum%tb19h%l, stats%maximum%tb19h%l)

   call da_proc_stats_combine(stats%average%tb22v, stats%rms_err%tb22v, &
      stats%minimum%tb22v%value, stats%maximum%tb22v%value, &
      stats%minimum%tb22v%n, stats%maximum%tb22v%n, &
      stats%minimum%tb22v%l, stats%maximum%tb22v%l)

   call da_proc_stats_combine(stats%average%tb37v, stats%rms_err%tb37v, &
      stats%minimum%tb37v%value, stats%maximum%tb37v%value, &
      stats%minimum%tb37v%n, stats%maximum%tb37v%n, &
      stats%minimum%tb37v%l, stats%maximum%tb37v%l)

   call da_proc_stats_combine(stats%average%tb37h, stats%rms_err%tb37h, &
      stats%minimum%tb37h%value, stats%maximum%tb37h%value, &
      stats%minimum%tb37h%n, stats%maximum%tb37h%n, &
      stats%minimum%tb37h%l, stats%maximum%tb37h%l)

   call da_proc_stats_combine(stats%average%tb85v, stats%rms_err%tb85v, &
      stats%minimum%tb85v%value, stats%maximum%tb85v%value, &
      stats%minimum%tb85v%n, stats%maximum%tb85v%n, &
      stats%minimum%tb85v%l, stats%maximum%tb85v%l)

   call da_proc_stats_combine(stats%average%tb85h, stats%rms_err%tb85h, &
      stats%minimum%tb85h%value, stats%maximum%tb85h%value, &
      stats%minimum%tb85h%n, stats%maximum%tb85h%n, &
      stats%minimum%tb85h%l, stats%maximum%tb85h%l)

   if (rootproc) then
      if (ntb19v /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

         write(unit=stats_unit, fmt='(a/)') '   var           tb19v(m/s)     n'
 
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb19v

         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tb19v%value, stats%minimum%tb19v%n,     &
            ' Maximum(n): ', stats%maximum%tb19v%value, stats%maximum%tb19v%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb19v/real(ntb19v), &
            '    RMSE   : ', sqrt(stats%rms_err%tb19v/real(ntb19v))
      end if

      if (ntb19h /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

         write(unit=stats_unit, fmt='(a/)') '   var           tb19h(m/s)     n'

         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb19h

         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tb19h%value, stats%minimum%tb19h%n,     &
            ' Maximum(n): ', stats%maximum%tb19h%value, stats%maximum%tb19h%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb19h/real(ntb19h), &
            '    RMSE   : ', sqrt(stats%rms_err%tb19h/real(ntb19h))
      end if

      if (ntb22v /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

         write(unit=stats_unit, fmt='(a/)') '   var           tb22v(m/s)     n'

         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb22v

         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tb22v%value, stats%minimum%tb22v%n,     &
            ' Maximum(n): ', stats%maximum%tb22v%value, &
                                        stats%maximum%tb22v%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb22v/real(ntb22v), &
            '    RMSE   : ', sqrt(stats%rms_err%tb22v/real(ntb22v))
     end if

     if (ntb37v /= 0) then
        write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

        write(unit=stats_unit, fmt='(a/)') '   var           tb37v(m/s)     n'

        write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb37v

        write(unit=stats_unit, fmt='(a, f12.4,i5)') &
           ' Minimum(n): ', stats%minimum%tb37v%value, stats%minimum%tb37v%n,     &
           ' Maximum(n): ', stats%maximum%tb37v%value, stats%maximum%tb37v%n
        write(unit=stats_unit, fmt='(a, f12.4,5x)') &
           ' Average   : ', stats%average%tb37v/real(ntb37v), &
           '    RMSE   : ', sqrt(stats%rms_err%tb37v/real(ntb37v))
      end if

      if (ntb37h /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

         write(unit=stats_unit, fmt='(a/)') '   var           tb37h(m/s)     n'
  
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb37h

         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tb37h%value, stats%minimum%tb37h%n,     &
            ' Maximum(n): ', stats%maximum%tb37h%value, stats%maximum%tb37h%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb37h/real(ntb37h), &
            '    RMSE   : ', sqrt(stats%rms_err%tb37h/real(ntb37h))
      end if

      if (ntb85v /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

         write(unit=stats_unit, fmt='(a/)') '   var           tb85v(m/s)     n' 

         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb85v

         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tb85v%value, stats%minimum%tb85v%n,     &
            ' Maximum(n): ', stats%maximum%tb85v%value, stats%maximum%tb85v%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb85v/real(ntb85v), &
            '    RMSE   : ', sqrt(stats%rms_err%tb85v/real(ntb85v))
      end if

      if (ntb85h /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for ssmi_tb'

         write(unit=stats_unit, fmt='(a/)') '   var           tb85h(m/s)     n'

         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb85h

         write(unit=stats_unit, fmt='(a, f12.4,i5)') &
            ' Minimum(n): ', stats%minimum%tb85h%value, stats%minimum%tb85h%n,     &
            ' Maximum(n): ', stats%maximum%tb85h%value, stats%maximum%tb85h%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb85h/real(ntb85h), &
            '    RMSE   : ', sqrt(stats%rms_err%tb85h/real(ntb85h))
      end if
   end if

   if (trace_use) call da_trace_exit("da_ao_stats_ssmi_tb")

end subroutine da_ao_stats_ssmi_tb


subroutine da_read_obs_ssmi (iv, filename)

   !---------------------------------------------------------------------------
   ! Purpose: Read a SSMI 2.0 GTS observation file
   !---------------------------------------------------------------------------

   implicit none

   type(iv_type),    intent(inout) :: iv
   character(len=*), intent(in)    :: filename

   character (len =  10)        :: fmt_name

   character (len = 160)        :: fmt_info
   character (len = 160)        :: fmt_loc
   character (len = 120)        :: char_ned

   integer                      :: iost, fm,iunit

   type (model_loc_type)        :: loc
   type (info_type)             :: info
   type (field_type)            :: speed, tpw

   type (field_type)            :: tb19v, tb19h, tb22v
   type (field_type)            :: tb37v, tb37h, tb85v, tb85h

   type (count_obs_number_type) :: count_obs_ssmir
   type (count_obs_number_type) :: count_obs_ssmit

   logical                      :: isfilter,ipass 
   logical                      :: outside, outside_all
   integer                      :: irain, iprecip
   integer                      :: n, ndup
   integer                      :: nlocal(num_ob_indexes)
   integer                      :: ntotal(num_ob_indexes)

   if (trace_use) call da_trace_entry("da_read_obs_ssmi")

   nlocal(:) = iv%info(:)%plocal(iv%time-1)
   ntotal(:) = iv%info(:)%ptotal(iv%time-1)

   count_obs_ssmir = count_obs_number_type(0, 0, 0, 0)
   count_obs_ssmit = count_obs_number_type(0, 0, 0, 0)

   isfilter = .true. ! filter out rain points
   irain = 0

   ! open file

   call da_get_unit(iunit)
   open(unit   = iunit,     &
      FILE   = trim(filename), &
      FORM   = 'FORMATTED',  &
      ACCESS = 'SEQUENTIAL', &
      iostat =  iost,     &
      STATUS = 'OLD')

   if (iost /= 0) then
      call da_warning("da_read_obs_ssmi.inc",59, (/"Cannot open SSMI file "//filename/))
      call da_free_unit(iunit)
      return
   end if

   rewind (unit = iunit)

   ! 2.  read header
   ! ===============

   ! 2.1 read first line
   !     ---------------

   read (unit = iunit, fmt = '(A)', iostat = iost) char_ned
   if (iost /= 0) then
      call da_error("da_read_obs_ssmi.inc",74, (/"Cannot read SSMI file"//filename/))
   end if

   ! 2.3 read number of reports
   !     ---------------------

   do
      read (unit = iunit, fmt = '(A)', iostat = iost) char_ned
      if (iost /= 0) then
         call da_error("da_read_obs_ssmi.inc",83, (/"Cannot read SSMI file"//filename/))
      end if
      if (char_ned(1:6) == 'NESTIX') exit
   end do

   do
     read (unit = iunit, fmt = '(A)', iostat = iost) char_ned
     if (char_ned(1:6) == 'INFO  ') exit
   end do

   read (unit = iunit, fmt = '(A)', iostat = iost) char_ned

   ! read formats
   ! ------------

   read (unit=iunit, fmt = '(A,1X,A)') fmt_name, fmt_info, fmt_name, fmt_loc

   !  skip 1 line
   !  -----------

   read (unit=iunit, fmt = '(A)') fmt_name

   !  loop over records
   !  -----------------

   reports: do
      ! read station general info
      ! =========================

      read (unit=iunit, fmt = fmt_info, iostat = iost) &
         info % platform,    &
         info % date_char,   &
         info % name,        &
         info % levels,      &
         info % lat,         &
         info % lon,         &
         info % elv,         &
         info % id

      read(unit=info % platform (4:6),fmt='(I3)') fm
      if (iost /= 0) exit reports

      select case(fm)
         case (125)    ;
            ! read surface wind speed and precipitable water
            read (unit=iunit, fmt = fmt_loc) speed%inv, speed%qc, speed%error, &
                                             tpw%inv, tpw%qc, tpw%error
         case (126)    ;
            read (unit=iunit, fmt = fmt_loc) &
               tb19v%inv, tb19v%qc, tb19v%error, &
               tb19h%inv, tb19h%qc, tb19h%error, &
               tb22v%inv, tb22v%qc, tb22v%error, &
               tb37v%inv, tb37v%qc, tb37v%error, &
               tb37h%inv, tb37h%qc, tb37h%error, &
               tb85v%inv, tb85v%qc, tb85v%error, &
               tb85h%inv, tb85h%qc, tb85h%error

               tb19v % error = tb19v % error + 2.0
               tb19h % error = tb19h % error + 2.0
               tb22v % error = tb22v % error + 2.0
               tb37h % error = tb37h % error + 2.0
               tb37v % error = tb37v % error + 2.0
               tb85h % error = tb85h % error + 2.0
               tb85v % error = tb85v % error + 2.0

         case default;
            write(unit=message(1), fmt='(a, i6)') 'unsaved ssmi obs found, fm=', fm
            write(unit=message(2), fmt='(a, 2f12.6)') &
               'info%(lon,lat)=', info%lon, info%lat
            call da_warning("da_read_obs_ssmi.inc",152,message(1:2))
      end select

      ! check if obs is in horizontal domain
      ! ====================================

      ! Compute the model horizontal coordinate x, y
      ! Check if obs is wihin horizontal domain

      call da_llxy (info, loc, outside, outside_all)

      if (outside_all) cycle reports

      loc % pw  % inv = missing_r
      loc % pw  % qc  = missing_data
      loc % slp       = loc % pw

      ! Loop over duplicating obs for global
      ndup = 1
      if (global .and. (loc%i < ids .or. loc%i >= ide)) ndup= 2

      ! It is possible that logic for counting obs is incorrect for the
      ! global case with >1 MPI tasks due to obs duplication, halo, etc.
      ! TBH:  20050913

      if (.not.outside) then
         if (print_detail_obs .and. ndup > 1) then
            write(unit=stdout, fmt = fmt_info) &
               info%platform,    &
               info%date_char,   &
               info%name,        &
               info%levels,      &
               info%lat,         &
               info%lon,         &
               info%elv,         &
               info%id

            write(unit=stdout, fmt = '(a,2i5,4e20.10)') &
               ' duplicating obs since loc% i,j,dx,dxm,dy & dym ', &
               loc%i, loc%j, loc%dx, loc%dxm, loc%dy, loc%dym
         end if
      end if

      dup_loop: do n = 1, ndup

      select case(fm)
         case (125) ;
            if (.not. use_ssmiretrievalobs) cycle reports

            if (n==1) ntotal(ssmi_rv) = ntotal(ssmi_rv)
            if (outside) cycle reports

            ! Check if at least one field is present
            if ((tpw % qc == missing_data) .AND. (speed % qc == missing_data)) then
               count_obs_ssmir % num_missing = count_obs_ssmir % num_missing + 1
               cycle reports
            end if

            ! fill permanent structure
            ! ========================

            nlocal(ssmi_rv) = nlocal(ssmi_rv) + 1
            ! Track serial obs index for reassembly of obs during bit-for-bit
            ! tests with different numbers of MPI tasks.  
            loc%obs_global_index = ntotal(ssmi_rv)

            !  One more data used

            count_obs_ssmir % num_used = count_obs_ssmir % num_used + 1
      
            !  Initialize other non read fields

            iv%info(ssmi_rv)%name(nlocal(ssmi_rv))      = info%name
            iv%info(ssmi_rv)%platform(nlocal(ssmi_rv))  = info%platform
            iv%info(ssmi_rv)%id(nlocal(ssmi_rv))        = info%id
            iv%info(ssmi_rv)%date_char(nlocal(ssmi_rv)) = info%date_char
            iv%info(ssmi_rv)%levels(nlocal(ssmi_rv))    = 1
            iv%info(ssmi_rv)%lat(:,nlocal(ssmi_rv))     = info%lat
            iv%info(ssmi_rv)%lon(:,nlocal(ssmi_rv))     = info%lon
            iv%info(ssmi_rv)%elv(nlocal(ssmi_rv))       = info%elv
            iv%info(ssmi_rv)%pstar(nlocal(ssmi_rv))     = info%pstar

            iv%info(ssmi_rv)%slp(nlocal(ssmi_rv))           = loc%slp
            iv%info(ssmi_rv)%pw(nlocal(ssmi_rv))            = loc%pw
            iv%info(ssmi_rv)%x(:,nlocal(ssmi_rv))           = loc%x
            iv%info(ssmi_rv)%y(:,nlocal(ssmi_rv))           = loc%y 
            iv%info(ssmi_rv)%i(:,nlocal(ssmi_rv))           = loc%i 
            iv%info(ssmi_rv)%j(:,nlocal(ssmi_rv))           = loc%j 
            iv%info(ssmi_rv)%dx(:,nlocal(ssmi_rv))          = loc%dx
            iv%info(ssmi_rv)%dxm(:,nlocal(ssmi_rv))         = loc%dxm
            iv%info(ssmi_rv)%dy(:,nlocal(ssmi_rv))          = loc%dy
            iv%info(ssmi_rv)%dym(:,nlocal(ssmi_rv))         = loc%dym
            iv%info(ssmi_rv)%proc_domain(:,nlocal(ssmi_rv)) = loc%proc_domain

            iv%info(ssmi_rv)%obs_global_index(nlocal(ssmi_rv)) = ntotal(ssmi_rv)

            iv % ssmi_rv (nlocal(ssmi_rv)) % speed = speed
            iv % ssmi_rv (nlocal(ssmi_rv)) % tpw   = tpw

         case (126) ;
            if (.not. use_ssmitbobs) cycle reports

            if (n==1) ntotal(ssmi_tb) = ntotal(ssmi_tb) + 1
            if (outside) cycle reports

            ! Check if at least one field is present

            if ((tb19v % qc == missing_data) .AND. (tb19h % qc == missing_data)  .AND. &
                (tb22v % qc == missing_data)                                .AND. &
                (tb37v % qc == missing_data) .AND. (tb37h % qc == missing_data)  .AND. &
                (tb85v % qc == missing_data) .AND. (tb85h % qc == missing_data)) then
               count_obs_ssmit % num_missing = &
               count_obs_ssmit % num_missing + 1
               ! write (unit=stdout,fmt=*) 'missing data'
               cycle reports
            end if

            ! filter rain pixels
            !  ====================================

            if (isfilter) then
               ipass = .false.
               iprecip = 0
               call filter(tb19v%inv, tb19h%inv, tb22v%inv, tb37v%inv, &
                  tb37h%inv, tb85v%inv, tb85h%inv, ipass, iprecip, info%lat)
               if (iprecip .eq. 1) then
                  tb19v % qc    = -88.0
                  tb19h % qc    = -88.0
                  tb22v % qc    = -88.0
                  tb37h % qc    = -88.0
                  tb37v % qc    = -88.0
                  tb85h % qc    = -88.0
                  tb85v % qc    = -88.0
                  irain = irain + 1
                  cycle reports
               end if
            end if

            ! fill permanent structure
            ! ========================

            ! One more data read in

            nlocal(ssmi_tb) = nlocal(ssmi_tb) + 1
            ! Track serial obs index for reassembly of obs during bit-for-bit
            ! tests with different numbers of MPI tasks.  
            loc%obs_global_index = ntotal(ssmi_tb)

            !  One more data used

            count_obs_ssmit % num_used = count_obs_ssmit % num_used + 1

            !  Initialize other non read fields

            iv%info(ssmi_tb)%name(nlocal(ssmi_tb))      = info%name
            iv%info(ssmi_tb)%platform(nlocal(ssmi_tb))  = info%platform
            iv%info(ssmi_tb)%id(nlocal(ssmi_tb))        = info%id
            iv%info(ssmi_tb)%date_char(nlocal(ssmi_tb)) = info%date_char
            iv%info(ssmi_tb)%levels(nlocal(ssmi_tb))    = 1
            iv%info(ssmi_tb)%lat(:,nlocal(ssmi_tb))     = info%lat
            iv%info(ssmi_tb)%lon(:,nlocal(ssmi_tb))     = info%lon
            iv%info(ssmi_tb)%elv(nlocal(ssmi_tb))       = info%elv
            iv%info(ssmi_tb)%pstar(nlocal(ssmi_tb))     = info%pstar

            iv%info(ssmi_tb)%slp(nlocal(ssmi_tb))           = loc%slp
            iv%info(ssmi_tb)%pw(nlocal(ssmi_tb))            = loc%pw
            iv%info(ssmi_tb)%x(:,nlocal(ssmi_tb))           = loc%x
            iv%info(ssmi_tb)%y(:,nlocal(ssmi_tb))           = loc%y 
            iv%info(ssmi_tb)%i(:,nlocal(ssmi_tb))           = loc%i 
            iv%info(ssmi_tb)%j(:,nlocal(ssmi_tb))           = loc%j 
            iv%info(ssmi_tb)%dx(:,nlocal(ssmi_tb))          = loc%dx
            iv%info(ssmi_tb)%dxm(:,nlocal(ssmi_tb))         = loc%dxm
            iv%info(ssmi_tb)%dy(:,nlocal(ssmi_tb))          = loc%dy
            iv%info(ssmi_tb)%dym(:,nlocal(ssmi_tb))         = loc%dym
            iv%info(ssmi_tb)%proc_domain(:,nlocal(ssmi_tb)) = loc%proc_domain

            iv%info(ssmi_tb)%obs_global_index(nlocal(ssmi_tb)) = ntotal(ssmi_tb)

            iv % ssmi_tb (nlocal(ssmi_tb)) % tb19v = tb19v
            iv % ssmi_tb (nlocal(ssmi_tb)) % tb19h = tb19h
            iv % ssmi_tb (nlocal(ssmi_tb)) % tb22v = tb22v
            iv % ssmi_tb (nlocal(ssmi_tb)) % tb37v = tb37v
            iv % ssmi_tb (nlocal(ssmi_tb)) % tb37h = tb37h
            iv % ssmi_tb (nlocal(ssmi_tb)) % tb85v = tb85v
            iv % ssmi_tb (nlocal(ssmi_tb)) % tb85h = tb85h

         case default;
            ! Do nothing.
      end select
      end do dup_loop
   end do reports

   close(iunit)
   call da_free_unit(iunit)

   write(unit=stdout, fmt='(/,25x,a)')   '     used   outdomain  max_er_chk   missing' 
   write(unit=stdout, fmt='(4x,a,4i11)') 'ssmi_rv_diag:        ', count_obs_ssmir
   write(unit=stdout, fmt='(4x,a,4i11)') 'ssmi_tb_diag:        ', count_obs_ssmit

   if (irain > 0) then
      write(unit=stdout, fmt='(/,5x,a,i6/)') '** Rain contaminated SSMI_Tb =', irain
   end if

   write(unit=stdout, fmt='(/,a)') ' '

   if (trace_use) call da_trace_exit("da_read_obs_ssmi")

end subroutine da_read_obs_ssmi


subroutine da_scan_obs_ssmi (iv, filename)

   !---------------------------------------------------------------------------
   ! Purpose: Read the header of a SSMI 2.0 GTS observation file
   !---------------------------------------------------------------------------

   implicit none

   type(iv_type),    intent(inout) :: iv
   character(len=*), intent(in)    :: filename

   character (len =  10)        :: fmt_name
   character (len = 160)        :: fmt_info, fmt_loc
   character (len = 120)        :: char_ned

   integer                      :: iost, fm,iunit

   type (model_loc_type)        :: loc
   type (info_type)             :: info
   type (field_type)            :: speed, tpw

   type (field_type)            :: tb19v, tb19h, tb22v
   type (field_type)            :: tb37v, tb37h, tb85v, tb85h

   logical                      :: isfilter,ipass 
   logical                      :: outside
   integer                      :: irain, iprecip

   if (trace_use) call da_trace_entry("da_scan_obs_ssmi")

   isfilter = .true. ! filter out rain points
   irain = 0

   ! open FILE
   call da_get_unit(iunit)
   open(unit   = iunit,     &
      FILE   = trim(filename), &
      FORM   = 'FORMATTED',  &
      ACCESS = 'SEQUENTIAL', &
      iostat =  iost,     &
      STATUS = 'OLD')

   if (iost /= 0) then
      call da_warning("da_scan_obs_ssmi.inc",44, (/"Cannot open SSMI file "//filename/))
      call da_free_unit(iunit)
      return
   end if

   rewind (unit = iunit)

   ! 2.  read header
   ! ===============

   ! 2.1 read first line
   !     ---------------

   read (unit = iunit, fmt = '(A)', iostat = iost) char_ned
   if (iost /= 0) then
      call da_error("da_scan_obs_ssmi.inc",59, (/"Cannot read SSMI file "//filename/))
   end if

   ! 2.3 read number of reports
   do
      read (unit = iunit, fmt = '(A)', iostat = iost) char_ned
      if (iost /= 0) then
         call da_error("da_scan_obs_ssmi.inc",66, (/"Cannot read SSMI file "//filename/))
      end if
      if (char_ned(1:6) == 'NESTIX') exit
   end do

   do
      read (unit = iunit, fmt = '(A)', iostat = iost) char_ned
      if (char_ned(1:6) == 'INFO  ') exit
   end do

   read (unit = iunit, fmt = '(A)', iostat = iost) char_ned

   ! read formats
   read (unit=iunit, fmt = '(A,1X,A)') &
      fmt_name, fmt_info, &
      fmt_name, fmt_loc

   ! skip 1 line
   read (unit=iunit, fmt = '(A)') fmt_name

   ! loop over records
   reports: do
      ! read station general info
      read (unit=iunit, fmt = fmt_info, iostat = iost) &
         info % platform,    &
         info % date_char,   &
         info % name,        &
         info % levels,      &
         info % lat,         &
         info % lon,         &
         info % elv,         &
         info % id

      read(unit=info % platform (4:6),fmt='(I3)') fm
      if (iost /= 0) exit reports

      select case(fm)
      case (125)    ;
         ! read surface wind speed and precipitable water
         read (unit=iunit, fmt = fmt_loc) &
               speed%inv, speed%qc, speed%error, &
               tpw%inv, tpw%qc, tpw%error
      case (126)    ;
         read (unit=iunit, fmt = fmt_loc)  &
               tb19v%inv, tb19v%qc, tb19v%error, &
               tb19h%inv, tb19h%qc, tb19h%error, &
               tb22v%inv, tb22v%qc, tb22v%error, &
               tb37v%inv, tb37v%qc, tb37v%error, &
               tb37h%inv, tb37h%qc, tb37h%error, &
               tb85v%inv, tb85v%qc, tb85v%error, &
               tb85h%inv, tb85h%qc, tb85h%error
      case default;
         write(unit=message(1), fmt='(a, i6)') 'unsaved ssmi obs found, fm=', fm
         write(unit=message(2), fmt='(a, 2f12.6)') 'info%(lon,lat)=', info%lon, info%lat
         call da_warning("da_scan_obs_ssmi.inc",120,message(1:2))
      end select

      ! check if obs is in horizontal domain

      ! Compute the model horizontal coordinate x, y
      ! Check if obs is wihin horizontal domain

      call da_llxy (info, loc, outside)

      select case(fm)
      case (125) ;
         if (.not. use_ssmiretrievalobs .or. iv%info(ssmi_rv)%ntotal == max_ssmi_rv_input) cycle reports

         ! Check if at least one field is present
         if ((tpw % qc == missing) .AND. (speed % qc == missing)) then
            cycle reports
         end if
         iv%info(ssmi_rv)%ntotal = iv%info(ssmi_rv)%ntotal + 1
         if (outside) cycle reports

         iv%info(ssmi_rv)%nlocal = iv%info(ssmi_rv)%nlocal + 1

      case (126) ;
         if (.not. use_ssmitbobs .or. iv%info(ssmi_tb)%ntotal == max_ssmi_tb_input) cycle reports

         ! Check if at least one field is present

         if ((tb19v % qc == missing) .AND. (tb19h % qc == missing)  .AND. &
             (tb22v % qc == missing)                                .AND. &
             (tb37v % qc == missing) .AND. (tb37h % qc == missing)  .AND. &
             (tb85v % qc == missing) .AND. (tb85h % qc == missing)) then
            cycle reports
         end if

         ! filter rain pixels
         ! ====================================

         if (isfilter) then
            ipass = .false.
            iprecip = 0
            call filter(tb19v%inv, tb19h%inv, tb22v%inv, tb37v%inv, &
               tb37h%inv, tb85v%inv, tb85h%inv, ipass, iprecip, &
               info%lat)
            if (iprecip == 1) then
               irain = irain + 1
               cycle reports
            end if
         end if

         iv%info(ssmi_tb)%ntotal = iv%info(ssmi_tb)%ntotal + 1
         if (outside) cycle reports
         iv%info(ssmi_tb)%nlocal = iv%info(ssmi_tb)%nlocal + 1

      case default;
         ! Do nothing.
      end select

   end do reports

   iv%info(ssmt1)%max_lev   = 1
   iv%info(ssmt2)%max_lev   = 1
   iv%info(ssmi_tb)%max_lev = 1
   iv%info(ssmi_rv)%max_lev = 1

   close(unit=iunit)
   call da_free_unit(iunit)

   if (irain > 0) then
      write(unit=stdout, fmt='(/,5x,a,i6/)') 'Rejected rain contaminated ssmi_tb =', irain
      write(unit=stdout, fmt='(A)') ' '
   end if

   if (trace_use) call da_trace_exit("da_scan_obs_ssmi")

end subroutine da_scan_obs_ssmi


subroutine da_jo_and_grady_ssmi_rv(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Ob Inc. structure.
   type (y_type),  intent(in)    :: re          ! Residual structure.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer  :: n

   if (trace_use) call da_trace_entry("da_jo_and_grady_ssmi_rv")

   jo % ssmir_speed = 0.0
   jo % ssmir_tpw   = 0.0

   do n=1, iv%info(ssmi_rv)%nlocal
      jo_grad_y%ssmi_rv(n)%speed = - re%ssmi_rv(n)%speed / &
         (iv%ssmi_rv(n)%speed%error * iv%ssmi_rv(n)%speed%error)

      jo_grad_y%ssmi_rv(n)%tpw = -re%ssmi_rv(n)%tpw / &
          (iv%ssmi_rv(n)%tpw%error * iv%ssmi_rv(n)%tpw%error)

      if (iv%info(ssmi_rv)%proc_domain(1,n)) then

         jo%ssmir_speed = jo%ssmir_speed - re%ssmi_rv(n)%speed * jo_grad_y%ssmi_rv(n)%speed
         jo%ssmir_tpw   = jo%ssmir_tpw   - re%ssmi_rv(n)%tpw   * jo_grad_y%ssmi_rv(n)%tpw
      end if
   end do
   
   jo % ssmir_speed = 0.5 * jo % ssmir_speed
   jo % ssmir_tpw   = 0.5 * jo % ssmir_tpw

   if (trace_use) call da_trace_exit("da_jo_and_grady_ssmi_rv")

end subroutine da_jo_and_grady_ssmi_rv


subroutine da_jo_and_grady_ssmi_tb(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Ob Inc. structure.
   type (y_type),  intent(in)    :: re          ! Residual structure.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(out)   :: jo          ! Obs cost function.

   integer :: n

   if (trace_use) call da_trace_entry("da_jo_and_grady_ssmi_tb")

   jo % ssmi_tb19v = 0.0
   jo % ssmi_tb19h = 0.0
   jo % ssmi_tb22v = 0.0
   jo % ssmi_tb37v = 0.0
   jo % ssmi_tb37h = 0.0
   jo % ssmi_tb85v = 0.0
   jo % ssmi_tb85h = 0.0

   do n=1, iv%info(ssmi_tb)%nlocal

      jo_grad_y%ssmi_tb(n)%tb19v = - re%ssmi_tb(n)%tb19v / &
         (iv%ssmi_tb(n)%tb19v%error * iv%ssmi_tb(n)%tb19v%error)

      jo_grad_y%ssmi_tb(n)%tb19h = - re%ssmi_tb(n)%tb19h / &
          (iv%ssmi_tb(n)%tb19h%error * iv%ssmi_tb(n)%tb19h%error)

      jo_grad_y%ssmi_tb(n)%tb22v = - re%ssmi_tb(n)%tb22v / &
         (iv%ssmi_tb(n)%tb22v%error * iv%ssmi_tb(n)%tb22v%error)

      jo_grad_y%ssmi_tb(n)%tb37v = - re%ssmi_tb(n)%tb37v / &
         (iv%ssmi_tb(n)%tb37v%error * iv%ssmi_tb(n)%tb37v%error)

      jo_grad_y%ssmi_tb(n)%tb37h = - re%ssmi_tb(n)%tb37h / &
         (iv%ssmi_tb(n)%tb37h%error * iv%ssmi_tb(n)%tb37h%error)

      jo_grad_y%ssmi_tb(n)%tb85v = - re%ssmi_tb(n)%tb85v / &
         (iv%ssmi_tb(n)%tb85v%error * iv%ssmi_tb(n)%tb85v%error)

      jo_grad_y%ssmi_tb(n)%tb85h = - re%ssmi_tb(n)%tb85h / &
         (iv%ssmi_tb(n)%tb85h%error * iv%ssmi_tb(n)%tb85h%error)

      if (iv%info(ssmi_tb)%proc_domain(1,n)) then
         jo % ssmi_tb19v = jo % ssmi_tb19v - re%ssmi_tb(n)%tb19v * jo_grad_y%ssmi_tb(n)%tb19v 
         jo % ssmi_tb19h = jo % ssmi_tb19h - re%ssmi_tb(n)%tb19h * jo_grad_y%ssmi_tb(n)%tb19h
         jo % ssmi_tb22v = jo % ssmi_tb22v - re%ssmi_tb(n)%tb22v * jo_grad_y%ssmi_tb(n)%tb22v
         jo % ssmi_tb37v = jo % ssmi_tb37v - re%ssmi_tb(n)%tb37v * jo_grad_y%ssmi_tb(n)%tb37v
         jo % ssmi_tb37h = jo % ssmi_tb37h - re%ssmi_tb(n)%tb37h * jo_grad_y%ssmi_tb(n)%tb37h
         jo % ssmi_tb85v = jo % ssmi_tb85v - re%ssmi_tb(n)%tb85v * jo_grad_y%ssmi_tb(n)%tb85v
         jo % ssmi_tb85h = jo % ssmi_tb85h - re%ssmi_tb(n)%tb85h * jo_grad_y%ssmi_tb(n)%tb85h
      end if
   end do
   
   jo % ssmi_tb19h = 0.5 * jo % ssmi_tb19h
   jo % ssmi_tb19v = 0.5 * jo % ssmi_tb19v
   jo % ssmi_tb22v = 0.5 * jo % ssmi_tb22v
   jo % ssmi_tb37h = 0.5 * jo % ssmi_tb37h
   jo % ssmi_tb37v = 0.5 * jo % ssmi_tb37v
   jo % ssmi_tb85h = 0.5 * jo % ssmi_tb85h
   jo % ssmi_tb85v = 0.5 * jo % ssmi_tb85v

   if (trace_use) call da_trace_exit("da_jo_and_grady_ssmi_tb")

end subroutine da_jo_and_grady_ssmi_tb


subroutine da_residual_ssmi_rv(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Obs increment vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual structure.

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type) :: n_obs_bad
   integer              :: n

   if (trace_use) call da_trace_entry("da_residual_ssmi_rv")

   n_obs_bad % Speed % num = number_type(0, 0, 0)
   n_obs_bad % q % num     = number_type(0, 0, 0)

   do n=1, iv%info(ssmi_rv)%nlocal
      np_available = np_available + 2

      re%ssmi_rv(n)%Speed = da_residual(n, 0, y%ssmi_rv(n)%Speed, &
         iv%ssmi_rv(n)%Speed, n_obs_bad % Speed)
      re%ssmi_rv(n)%tpw   = da_residual(n, 0, y%ssmi_rv(n)%tpw,   &
         iv%ssmi_rv(n)%tpw, n_obs_bad % q      )
   end do

   np_missing  = np_missing + n_obs_bad % Speed % num % miss + n_obs_bad % q % num % miss
   np_bad_data = np_bad_data + n_obs_bad % Speed % num % bad + n_obs_bad % q % num % bad
   np_obs_used = np_obs_used + n_obs_bad % Speed % num % use + n_obs_bad % q % num % use

   if (trace_use) call da_trace_exit("da_residual_ssmi_rv")

end subroutine da_residual_ssmi_rv


subroutine da_residual_ssmi_tb(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Obs increment vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual structure.

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type) :: n_obs_bad
   integer              :: n

   if (trace_use) call da_trace_entry("da_residual_ssmi_tb")

   n_obs_bad % t % num = number_type(0, 0, 0)

   do n=1, iv%info(ssmi_tb)%nlocal
      np_available = np_available + 7

      re%ssmi_tb(n)%tb19v = da_residual(n, 0, y%ssmi_tb(n)%tb19v, iv%ssmi_tb(n)%tb19v, n_obs_bad % t)
      re%ssmi_tb(n)%tb19h = da_residual(n, 0, y%ssmi_tb(n)%tb19h, iv%ssmi_tb(n)%tb19h, n_obs_bad % t)
      re%ssmi_tb(n)%tb22v = da_residual(n, 0, y%ssmi_tb(n)%tb22v, iv%ssmi_tb(n)%tb22v, n_obs_bad % t)
      re%ssmi_tb(n)%tb37v = da_residual(n, 0, y%ssmi_tb(n)%tb37v, iv%ssmi_tb(n)%tb37v, n_obs_bad % t)
      re%ssmi_tb(n)%tb37h = da_residual(n, 0, y%ssmi_tb(n)%tb37h, iv%ssmi_tb(n)%tb37h, n_obs_bad % t)
      re%ssmi_tb(n)%tb85v = da_residual(n, 0, y%ssmi_tb(n)%tb85v, iv%ssmi_tb(n)%tb85v, n_obs_bad % t)
      re%ssmi_tb(n)%tb85h = da_residual(n, 0, y%ssmi_tb(n)%tb85h, iv%ssmi_tb(n)%tb85h, n_obs_bad % t)
   end do

   np_missing  = np_missing  + n_obs_bad % t % num % miss  
   np_bad_data = np_bad_data + n_obs_bad % t % num % bad  
   np_obs_used = np_obs_used + n_obs_bad % t % num % use    

   if (trace_use) call da_trace_exit("da_residual_ssmi_tb") 

end subroutine da_residual_ssmi_tb


subroutine da_oi_stats_ssmi_rv(stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! iv

   integer                   :: ntpw, nspeed, n
   type (stats_ssmi_rv_type) :: stats

   if (trace_use) call da_trace_entry("da_oi_stats_ssmi_rv")

   ntpw = 0
   nspeed = 0

   stats%maximum%tpw   = maxmin_type(missing_r, 0, 0)
   stats%maximum%speed = maxmin_type(missing_r, 0, 0)
   stats%minimum%tpw   = maxmin_type(-missing_r, 0, 0)
   stats%minimum%speed = maxmin_type(-missing_r, 0, 0)
   stats%average = residual_ssmi_rv_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(ssmi_rv)%nlocal
      if (iv%info(ssmi_rv)%proc_domain(1,n)) then
         call da_stats_calculate(iv%info(ssmi_rv)%obs_global_index(n), &
           0, iv%ssmi_rv(n)%tpw%qc, &
           iv%ssmi_rv(n)%tpw%inv, ntpw, &
           stats%minimum%tpw, &
           stats%maximum%tpw,  &
           stats%average%tpw, &
           stats%rms_err%tpw)
        call da_stats_calculate(iv%info(ssmi_rv)%obs_global_index(n), &
           0, iv%ssmi_rv(n)%speed%qc, &
           iv%ssmi_rv(n)%speed%inv, nspeed, &
           stats%minimum%speed, &
           stats%maximum%speed, &
           stats%average%speed, &
           stats%rms_err%speed)
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(ntpw)
   call da_proc_sum_int(nspeed)

   call da_proc_stats_combine(stats%average%tpw, stats%rms_err%tpw, &
       stats%minimum%tpw%value, stats%maximum%tpw%value, &
       stats%minimum%tpw%n, stats%maximum%tpw%n, &
       stats%minimum%tpw%l, stats%maximum%tpw%l)

   call da_proc_stats_combine(stats%average%speed, stats%rms_err%speed, &
      stats%minimum%speed%value, stats%maximum%speed%value, &
      stats%minimum%speed%n, stats%maximum%speed%n, &
      stats%minimum%speed%l, stats%maximum%speed%l)

   if (rootproc) then
      if (ntpw > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_retrieval'
         write(unit=stats_unit, fmt='(a/)') '   var           tpw(cm)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntpw
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tpw%value, stats%minimum%tpw%n    , &
            ' Maximum(n): ', stats%maximum%tpw%value, stats%maximum%tpw%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tpw/real(ntpw), &
            '    RMSE   : ', sqrt(stats%rms_err%tpw/real(ntpw))
      end if

      if (nspeed > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_retrieval'
         write(unit=stats_unit, fmt='(a/)') '   var           speed(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', nspeed
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
           ' Minimum(n): ', stats%minimum%speed%value, stats%minimum%speed%n    , &
           ' Maximum(n): ', stats%maximum%speed%value, stats%maximum%speed%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
             ' Average   : ', stats%average%speed/real(nspeed), &
             '    RMSE   : ', sqrt(stats%rms_err%speed/real(nspeed))
      end if
   end if

   if (trace_use) call da_trace_exit("da_iv_stats_ssmi_rv")

end subroutine da_oi_stats_ssmi_rv


subroutine da_oi_stats_ssmi_tb (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! OI

   integer                   :: ntb19v,ntb19h,ntb22v,ntb37v,ntb37h, ntb85v,ntb85h
   integer                   :: n
   type (stats_ssmi_tb_type) :: stats

   if (trace_use) call da_trace_entry("da_oi_stats_ssmi_tb")

   ntb19v = 0
   ntb19h = 0
   ntb22v = 0
   ntb37v = 0
   ntb37h = 0
   ntb85v = 0
   ntb85h = 0

   stats%maximum%tb19v = maxmin_type(missing_r, 0, 0)
   stats%maximum%tb19h = maxmin_type(missing_r, 0, 0)
   stats%maximum%tb22v = maxmin_type(missing_r, 0, 0)
   stats%maximum%tb37v = maxmin_type(missing_r, 0, 0)
   stats%maximum%tb37h = maxmin_type(missing_r, 0, 0)
   stats%maximum%tb85v = maxmin_type(missing_r, 0, 0)
   stats%maximum%tb85h = maxmin_type(missing_r, 0, 0)
   stats%minimum%tb19v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb19h = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb22v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb37v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb37h = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb85v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%tb85h = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_ssmi_tb_type(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(ssmi_tb)%nlocal
      if (iv%info(ssmi_tb)%proc_domain(1,n)) then
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb19v%qc, & 
            iv%ssmi_tb(n)%tb19v%inv, ntb19v, &
            stats%minimum%tb19v, stats%maximum%tb19v,&
            stats%average%tb19v, stats%rms_err%tb19v)
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb19h%qc, & 
            iv%ssmi_tb(n)%tb19h%inv, ntb19h, &
            stats%minimum%tb19h, stats%maximum%tb19h,&
            stats%average%tb19h, stats%rms_err%tb19h)
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb22v%qc, & 
            iv%ssmi_tb(n)%tb22v%inv, ntb22v, &
            stats%minimum%tb22v, stats%maximum%tb22v,&
            stats%average%tb22v, stats%rms_err%tb22v)
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb37v%qc, & 
            iv%ssmi_tb(n)%tb37v%inv, ntb37v, &
            stats%minimum%tb37v, stats%maximum%tb37v,&
            stats%average%tb37v, stats%rms_err%tb37v)
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb37h%qc, & 
            iv%ssmi_tb(n)%tb37h%inv, ntb37h, &
            stats%minimum%tb37h, stats%maximum%tb37h,&
            stats%average%tb37h, stats%rms_err%tb37h)
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb85v%qc, & 
            iv%ssmi_tb(n)%tb85v%inv, ntb85v, &
            stats%minimum%tb85v, stats%maximum%tb85v,&
            stats%average%tb85v, stats%rms_err%tb85v)
         call da_stats_calculate(iv%info(ssmi_tb)%obs_global_index(n), &
            0, iv%ssmi_tb(n)%tb85h%qc, & 
            iv%ssmi_tb(n)%tb85h%inv, ntb85h, &
            stats%minimum%tb85h, stats%maximum%tb85h,&
            stats%average%tb85h, stats%rms_err%tb85h)
      end if
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(ntb19v)
   call da_proc_sum_int(ntb19h)
   call da_proc_sum_int(ntb22v)
   call da_proc_sum_int(ntb37v)
   call da_proc_sum_int(ntb37h)
   call da_proc_sum_int(ntb85v)
   call da_proc_sum_int(ntb85h)

   call da_proc_stats_combine(stats%average%tb19v, stats%rms_err%tb19v, &
      stats%minimum%tb19v%value, stats%maximum%tb19v%value, &
      stats%minimum%tb19v%n, stats%maximum%tb19v%n, &
      stats%minimum%tb19v%l, stats%maximum%tb19v%l)

   call da_proc_stats_combine(stats%average%tb19h, stats%rms_err%tb19h, &
      stats%minimum%tb19h%value, stats%maximum%tb19h%value, &
      stats%minimum%tb19h%n, stats%maximum%tb19h%n, &
      stats%minimum%tb19h%l, stats%maximum%tb19h%l)

   call da_proc_stats_combine(stats%average%tb22v, stats%rms_err%tb22v, &
      stats%minimum%tb22v%value, stats%maximum%tb22v%value, &
      stats%minimum%tb22v%n, stats%maximum%tb22v%n, &
      stats%minimum%tb22v%l, stats%maximum%tb22v%l)

   call da_proc_stats_combine(stats%average%tb37v, stats%rms_err%tb37v, &
      stats%minimum%tb37v%value, stats%maximum%tb37v%value, &
      stats%minimum%tb37v%n, stats%maximum%tb37v%n, &
      stats%minimum%tb37v%l, stats%maximum%tb37v%l)

   call da_proc_stats_combine(stats%average%tb37h, stats%rms_err%tb37h, &
      stats%minimum%tb37h%value, stats%maximum%tb37h%value, &
      stats%minimum%tb37h%n, stats%maximum%tb37h%n, &
      stats%minimum%tb37h%l, stats%maximum%tb37h%l)

   call da_proc_stats_combine(stats%average%tb85v, stats%rms_err%tb85v, &
      stats%minimum%tb85v%value, stats%maximum%tb85v%value, &
      stats%minimum%tb85v%n, stats%maximum%tb85v%n, &
      stats%minimum%tb85v%l, stats%maximum%tb85v%l)

   call da_proc_stats_combine(stats%average%tb85h, stats%rms_err%tb85h, &
      stats%minimum%tb85h%value, stats%maximum%tb85h%value, &
      stats%minimum%tb85h%n, stats%maximum%tb85h%n, &
      stats%minimum%tb85h%l, stats%maximum%tb85h%l)

   if (rootproc) then
      if (ntb19v > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb19v(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb19v
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb19v%value, stats%minimum%tb19v%n    , &
            ' Maximum(n): ', stats%maximum%tb19v%value, stats%maximum%tb19v%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb19v/real(ntb19v), &
            '    RMSE   : ', sqrt(stats%rms_err%tb19v/real(ntb19v))
      end if

      if (ntb19h > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb19h(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb19h
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb19h%value, stats%minimum%tb19h%n    , &
            ' Maximum(n): ', stats%maximum%tb19h%value, stats%maximum%tb19h%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb19h/real(ntb19h), &
            '    RMSE   : ', sqrt(stats%rms_err%tb19h/real(ntb19h))
      end if

      if (ntb22v > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb22v(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb22v
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb22v%value, stats%minimum%tb22v%n    , &
            ' Maximum(n): ', stats%maximum%tb22v%value, stats%maximum%tb22v%n
          write(unit=stats_unit, fmt='(a, f12.4,5x)') &
             ' Average   : ', stats%average%tb22v/real(ntb22v), &
             '    RMSE   : ', sqrt(stats%rms_err%tb22v/real(ntb22v))
      end if

      if (ntb37v > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb37v(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb37v
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb37v%value, stats%minimum%tb37v%n    , &
            ' Maximum(n): ', stats%maximum%tb37v%value, stats%maximum%tb37v%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb37v/real(ntb37v), &
            '    RMSE   : ', sqrt(stats%rms_err%tb37v/real(ntb37v))
      end if

      if (ntb37h > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb37h(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb37h
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb37h%value, stats%minimum%tb37h%n    , &
            ' Maximum(n): ', stats%maximum%tb37h%value, stats%maximum%tb37h%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb37h/real(ntb37h), &
            '    RMSE   : ', sqrt(stats%rms_err%tb37h/real(ntb37h))
      end if

      if (ntb85v > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb85v(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb85v
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb85v%value, stats%minimum%tb85v%n    , &
            ' Maximum(n): ', stats%maximum%tb85v%value, stats%maximum%tb85v%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb85v/real(ntb85v), &
            '    RMSE   : ', sqrt(stats%rms_err%tb85v/real(ntb85v))
      end if

      if (ntb85h > 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of OI for ssmi_tb'
         write(unit=stats_unit, fmt='(a/)') '   var           tb85h(m/s)     n'
         write(unit=stats_unit, fmt='(a,i14)') '  Number: ', ntb85h
         write(unit=stats_unit, fmt='(a, f12.4,i9)') &
            ' Minimum(n): ', stats%minimum%tb85h%value, stats%minimum%tb85h%n    , &
            ' Maximum(n): ', stats%maximum%tb85h%value, stats%maximum%tb85h%n
         write(unit=stats_unit, fmt='(a, f12.4,5x)') &
            ' Average   : ', stats%average%tb85h/real(ntb85h), &
            '    RMSE   : ', sqrt(stats%rms_err%tb85h/real(ntb85h))
      end if
   end if

   if (trace_use) call da_trace_exit("da_oi_stats_ssmi_tb")

end subroutine da_oi_stats_ssmi_tb


subroutine da_transform_xtospeed(Speed,U,V)

   !----------------------------------------------------------------------------
   ! Purpose: Convert (U-V in m/s) components into wind speed (Speed in m/s)
   !----------------------------------------------------------------------------

   implicit none

   real, intent (out) :: Speed
   real, intent (in)  :: U, V

   if (trace_use) call da_trace_entry("da_transform_xtospeed")

   Speed = sqrt(U*U + V*V)

   if (trace_use) call da_trace_exit("da_transform_xtospeed")

end subroutine da_transform_xtospeed


subroutine da_transform_xtospeed_lin(TGL_speed,U,V,TGL_u,TGL_v)

   !----------------------------------------------------------------------------
   ! Purpose: Convert (U-V in m/s) components into wind speed (Speed in m/s)
   !----------------------------------------------------------------------------

   implicit none

   real, intent(out) :: TGL_speed
   real, intent(in)  :: U,V
   real, intent(in)  :: TGL_u,TGL_v

   real :: speed

   if (trace_use) call da_trace_entry("da_transform_xtospeed_lin")

   speed = sqrt(U*U + V*V + 1.0e-6)

   TGL_speed = (U*TGL_u + V*TGL_v)/speed

   if (trace_use) call da_trace_exit("da_transform_xtospeed_lin")

end subroutine da_transform_xtospeed_lin


subroutine da_transform_xtospeed_adj(ADJ_speed,U,V,ADJ_u,ADJ_v)

   !----------------------------------------------------------------------------
   ! Purpose: Convert (U-V in m/s) components into wind speed (Speed in m/s)
   !----------------------------------------------------------------------------

   implicit none

   real, intent (in)    :: ADJ_speed
   real, intent (in)    :: U, V
   real, intent (inout) :: ADJ_u, ADJ_v

   real :: speed

   if (trace_use) call da_trace_entry("da_transform_xtospeed_adj")

   speed = sqrt(U*U+V*V+ 1.0e-6)

   ADJ_u = U*ADJ_speed/speed + ADJ_u
   ADJ_v = V*ADJ_speed/speed + ADJ_v

   if (trace_use) call da_trace_exit("da_transform_xtospeed_adj")

end subroutine da_transform_xtospeed_adj


subroutine da_transform_xtoseasfcwind(U,V,Speed,zhmkz)

   !----------------------------------------------------------------------------
   ! Purpose: Convert (U-V in m/s) components into wind speed (Speed in m/s)
   !----------------------------------------------------------------------------

   implicit none

   real, intent (out) :: Speed
   real, intent (in)  :: U, V, zhmkz

   real :: usfc, vsfc

   if (trace_use) call da_trace_entry("da_transform_xtoseasfcwind")

   usfc   = U*log(10./0.0001)/log(zhmkz/0.0001) ! roughness = 0.0001
   vsfc   = V*log(10./0.0001)/log(zhmkz/0.0001) ! roughness = 0.0001
   speed  = sqrt(usfc*usfc + vsfc*vsfc)

   if (trace_use) call da_trace_exit("da_transform_xtoseasfcwind")

end subroutine da_transform_xtoseasfcwind


subroutine da_transform_xtoseasfcwind_lin(grid)

   !----------------------------------------------------------------------------
   ! Purpose: Convert (U-V in m/s) components into wind speed (Speed in m/s)
   !----------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   real    :: const, rgh_fac, height
   integer :: i, j, k

   if (trace_use) call da_trace_entry("da_transform_xtoseasfcwind_lin")

   const = log(10.0/0.0001)
   k = grid%xb%kts
    
   do j=jts,jte
      do i=its,ite

        height = grid%xb%h(i,j,k) - grid%xb%terr(i,j)
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write(unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',grid%xb%h(i,j,k) ,' terr =  ',grid%xb%terr(i,j)
            call da_error("da_transform_xtoseasfcwind_lin.inc",27,message(1:2))
         end if

         rgh_fac = const/log(height/0.0001) ! roughness = 0.0001
         grid%xa%speed(i,j) = (grid%xa%u(i,j,k)*grid%xb%u(i,j,k) &
            + grid%xa%v(i,j,k)*grid%xb%v(i,j,k)) * rgh_fac*rgh_fac / grid%xb%speed(i,j)
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtoseasfcwind_lin")

end subroutine da_transform_xtoseasfcwind_lin


subroutine da_transform_xtoseasfcwind_adj(grid)

   !-------------------------------------------------------------------------
   ! Purpose: Convert (U-V in m/s) components into wind speed (Speed in m/s)
   !-------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   real    :: const, rgh_fac, var, height
   integer :: i, j, is, js, ie, je

   if (trace_use) call da_trace_entry("da_transform_xtoseasfcwind_adj")

   const = log(10./0.0001)

   is = its
   js = jts

   ie = ite
   je = jte

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if


   do j=js, je
      do i=is, ie
        height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write(unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',grid%xb%h(i,j,kts) ,' terr =  ',grid%xb%terr(i,j)
            call da_error("da_transform_xtoseasfcwind_adj.inc",46,message(1:2))
         end if

         rgh_fac = const/log(height/0.0001) ! roughness = 0.0001

         var = rgh_fac*rgh_fac/grid%xb%speed(i,j)

         grid%xa%u(i,j,kts)=grid%xa%u(i,j,kts)+var*grid%xa%speed(i,j)*grid%xb%u(i,j,kts)
         grid%xa%v(i,j,kts)=grid%xa%v(i,j,kts)+var*grid%xa%speed(i,j)*grid%xb%v(i,j,kts)
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtoseasfcwind_adj")

end subroutine da_transform_xtoseasfcwind_adj


subroutine da_transform_xtotb(grid)

   !----------------------------------------------------------------------
   ! Purpose: TBD
   !----------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)   :: grid

   integer            :: i,j,k

   real               :: psfc,ta,gamma,sst,htpw,speed,alw,zcld,tpw,dum1,zrhom

   if (trace_use) call da_trace_entry("da_transform_xtotb")

   do j=jts,jte
      do i=its,ite
         ! surface pressure (mb) (940 -1030)

         psfc          = 0.01*grid%xb%psfc(i,j)

         ! sea surface temperature (k) (273 - 303) (doesnot change) 

         sst           = grid%xb%tgrn(i,j)

         ! effective surface air temperature (263 - 303)

         ta          = grid%xb%tgrn(i,j) + &
                       (grid%xb%t(i,j,kts)-grid%xb%tgrn(i,j))*log(2.0/0.0001)/ &
                       log((grid%xb%h(i,j,kts) - grid%xb%terr(i,j))/0.0001)

         ! gamma is an emperical formula and zcld is given for now

         gamma   = (ta-270)*0.023 + 5.03  ! effective lapse rate(km) (4.0-6.5)

         zcld    = 1                           ! effective cloud height (km)
                                               ! = 1 if no cloud infomation
         ! total precipitable water in cm
         ! total precipitable water in (kg/m**2) (0 - 70)

         tpw          = grid%xb%tpw(i,j)*10.0
         speed        = grid%xb%speed(i,j)

         ! Column liquid water (kg/m**2)  (0-0.5) (no data now. So do it later)

         alw          = 0.0

         ! Column height weighted moisture density on the grid locally 

         zrhom = 0.0
         do k=kts,kte
            zrhom=zrhom+(grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*grid%xb%h(i,j,k)*grid%xb%q(i,j,k)* &
               grid%xb%rho(i,j,k)
         end do

         ! Column moisture density on the grid locally

         htpw          = zrhom/tpw/1000.0

         call tb(1,53.0,psfc,ta,gamma,sst,tpw,htpw, &
            speed,alw,zcld,grid%xb%tb19v(i,j),grid%xb%tb19h(i,j))
         call tb(2,53.0,psfc,ta,gamma,sst,tpw,htpw, &
            speed,alw,zcld,grid%xb%tb22v(i,j),dum1)
         call tb(3,53.0,psfc,ta,gamma,sst,tpw,htpw, &
            speed,alw,zcld,grid%xb%tb37v(i,j),grid%xb%tb37h(i,j))
         call tb(4,53.0,psfc,ta,gamma,sst,tpw,htpw, &
            speed,alw,zcld,grid%xb%tb85v(i,j),grid%xb%tb85h(i,j))
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_xtotb")

end subroutine da_transform_xtotb


subroutine da_transform_xtotb_lin(grid)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: i,j,k
   real    :: dum1, dum2, zrhom, TGL_zrhom

   real    :: psfc,ta,gamma,sst,htpw,speed,alw,zcld,tpw
   real    :: TGL_psfc,TGL_ta,TGL_gamma,TGL_sst,TGL_tpw
   real    :: TGL_htpw,TGL_speed,TGL_alw,TGL_zcld         

   if (trace_use) call da_trace_entry("da_transform_xtotb_lin")

   ! WHY?
   ! call da_transform_xtotpw(grid)

   do j=jts,jte
      do i=its,ite
         ! surface pressure (mb) (940 -1030)

         psfc     =  0.01*grid%xb%psfc(i,j)
         TGL_psfc =  0.01*grid%xa%psfc(i,j)

         ! sea surface temperature (k) (273 - 303) (doesnot change) 

         sst      = grid%xb%tgrn(i,j)
         ! TGL_sst  = grid%xa%tgrn(1,1)
         TGL_sst  = 0.0

         ! effective surface air temperature (263 - 303)

         ta          = grid%xb%tgrn(i,j) + &
                    (grid%xb%t(i,j,kts)-grid%xb%tgrn(i,j))*log(2.0/0.0001)/ &
                    log((grid%xb%h(i,j,kts)-grid%xb%terr(i,j))/0.0001)

         TGL_ta      = (grid%xa%t(i,j,kts)-0.)*log(2.0/0.0001)/ &
                    log((grid%xb%h(i,j,kts)-grid%xb%terr(i,j))/0.0001)

         ! gamma is an emperical formula and zcld is given for now

         gamma = (ta-270)*0.023 + 5.03  ! effective lapse rate (km) (4.0 - 6.5)
         zcld       = 1                           ! effective cloud height (km)
                                                 ! = 1 if no cloud infomation
         TGL_gamma = TGL_ta*0.023
         TGL_zcld  = 0.0
     

         ! total precipitable water in (kg/m**2) (0 - 70)

         tpw     = grid%xb%tpw(i,j)*10.0
         TGL_tpw = grid%xa%tpw(i,j)*10.0

         ! speed, surface wind speed (m/sec)    (0 - 30) , take 10 m wind

         speed     = grid%xb%speed(i,j)
         TGL_speed = grid%xa%speed(i,j)

         ! Column liquid water (kg/m**2)   (0 - 0.5) (no data now. So, do it later.)

         alw     = 0.0
         TGL_alw = 0.0

         ! Column height weighted moisture density on the grid locally 

         zrhom = 0.0
         do k=kts,kte
            zrhom=zrhom+(grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*grid%xb%h(i,j,k)* &
               grid%xb%q(i,j,k)*grid%xb%rho(i,j,k)
         end do

         TGL_zrhom = 0.0
         do k = kts,kte
            TGL_zrhom = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*grid%xb%h(i,j,k)* &
               (grid%xb%q(i,j,k)*grid%xa%rho(i,j,k) + &
               grid%xa%q(i,j,k)*grid%xb%rho(i,j,k)) + TGL_zrhom
         end do

         ! WHY?
         ! call da_transform_xtozrhoq(grid%xb, i, j, zh, zf, zrhom)
         ! call da_transform_xtozrhoq_lin(grid, i, j, zh, zf, TGL_zrhom)

         ! Column moisture density on the grid locally

         htpw     = zrhom/tpw/1000.0
         TGL_htpw = TGL_zrhom/tpw/1000.0 &
                  - TGL_tpw/tpw*htpw

         dum1 = 0.0

         call da_tb_tl(1,53.0,psfc,ta,gamma,sst,tpw,htpw,speed,alw,zcld,  &
            ! grid%xb%tb19v(i,j),grid%xb%tb19h(i,j),                      &
            TGL_psfc,TGL_ta,TGL_gamma,TGL_sst,                &
            TGL_tpw,TGL_htpw,TGL_speed,TGL_alw,               &
            TGL_zcld,grid%xa%tb19v(i,j),grid%xa%tb19h(i,j)             )

         call da_tb_tl(2,53.0,psfc,ta,gamma,sst,tpw,htpw,speed,alw,zcld,  &
            ! grid%xb%tb22v(i,j),dum1,                               &
            TGL_psfc,TGL_ta,TGL_gamma,TGL_sst,                &
            TGL_tpw,TGL_htpw,TGL_speed,TGL_alw,               &
            TGL_zcld,grid%xa%tb22v(i,j),dum2                      )

         call da_tb_tl(3,53.0,psfc,ta,gamma,sst,tpw,htpw,speed,alw,zcld,  &
            ! grid%xb%tb37v(i,j),grid%xb%tb37h(i,j),                      &
            TGL_psfc,TGL_ta,TGL_gamma,TGL_sst,                &
            TGL_tpw,TGL_htpw,TGL_speed,TGL_alw,               &
            TGL_zcld,grid%xa%tb37v(i,j),grid%xa%tb37h(i,j)             )

         call da_tb_tl(4,53.0,psfc,ta,gamma,sst,tpw,htpw,speed,alw,zcld,  &
            ! grid%xb%tb85v(i,j),grid%xb%tb85h(i,j),                      &
            TGL_psfc,TGL_ta,TGL_gamma,TGL_sst,                &
            TGL_tpw,TGL_htpw,TGL_speed,TGL_alw,               &
            TGL_zcld,grid%xa%tb85v(i,j),grid%xa%tb85h(i,j)             )
      end do
   end do       

   if (trace_use) call da_trace_exit("da_transform_xtotb_lin")

end subroutine da_transform_xtotb_lin


subroutine da_transform_xtotb_adj (grid)

   !----------------------------------------------------------------------
   ! Purpose: TBD
   !----------------------------------------------------------------------

   implicit none

   type (domain), intent(inout) :: grid

   integer :: i,j,k
   integer :: is, js, ie, je

   real    :: dx, dy, dxm, dym, zhmkz
   real    :: dum1, dum2, zrhom, ADJ_zrhom

   real    :: psfc,ta,gamma,sst,htpw,speed,alw,zcld,tpw
   real    :: ADJ_psfc,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_tpw
   real    :: ADJ_htpw,ADJ_speed,ADJ_alw,ADJ_zcld        

   if (trace_use) call da_trace_entry("da_transform_xtotb_adj")        

   psfc      = 0.0
   ta        = 0.0
   gamma     = 0.0
   sst       = 0.0
   htpw      = 0.0
   speed     = 0.0
   alw       = 0.0
   zcld      = 0.0
   tpw       = 0.0
   dx        = 0.0
   dy        = 0.0
   dxm       = 0.0
   dym       = 0.0
   zhmkz     = 0.0
   dum1      = 0.0
   dum2      = 0.0
   zrhom     = 0.0
   ADJ_zrhom = 0.0

   is = its
   js = jts

   ie = ite
   je = jte

   if (test_transforms) then
      is = its-1
      js = jts-1

      ie = ite+1
      je = jte+1

      if (is < ids) is = ids
      if (js < jds) js = jds

      if (ie > ide) ie = ide
      if (je > jde) je = jde
   end if

   ! Mean fields

   do j=js, je
      do i=is, ie
         psfc  = 0.01*grid%xb%psfc(i,j)
         ! sst   = grid%xb%tgrn(i,j)
         ta    = grid%xb%tgrn(i,j) + &
                 (grid%xb%t(i,j,kts)-grid%xb%tgrn(i,j))*log(2.0/0.0001)/ &
                 log((grid%xb%h(i,j,kts)- grid%xb%terr(i,j))/0.0001)

         gamma = (ta-270.0)*0.023 + 5.03  ! effective lapse rate   (km) (4.0-6.5)
         zcld  = 1                      ! effective cloud height (km)

         tpw   = grid%xb%tpw(i,j)*10.0
         ! speed = grid%xb%speed(i,j)

         alw   = 0.0

         zrhom = 0.0
         do k=kts,kte
            zrhom=zrhom+(grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*grid%xb%h(i,j,k)*grid%xb%q(i,j,k)* &
               grid%xb%rho(i,j,k)
         end do

         ! call da_transform_xtozrhoq(grid%xb, i, j, zh, zf, zrhom)

         htpw    = zrhom/tpw/1000.0

         dum1=0.0
         dum2=0.0

         ADJ_gamma    = 0.0
         ADJ_speed    = 0.0
         ADJ_psfc     = 0.0
         ADJ_zcld     = 0.0
         ADJ_htpw     = 0.0
         ADJ_sst      = 0.0
         ADJ_alw      = 0.0
         ADJ_tpw      = 0.0
         ADJ_ta       = 0.0
         ADJ_zrhom    = 0.0

         call da_tb_adj(1,53.0,psfc,ta,gamma,grid%xb%tgrn(i,j),tpw,      &
            htpw,grid%xb%speed(i,j),alw,zcld,               &
            ! grid%xb%tb19v(i,j),grid%xb%tb19h(i,j),               &
            ADJ_psfc,ADJ_ta,ADJ_gamma,ADJ_sst,         &
            ADJ_tpw,ADJ_htpw,ADJ_speed,ADJ_alw,        &
            ADJ_zcld,grid%xa%tb19v(i,j),grid%xa%tb19h(i,j)    )

         call da_tb_adj(2,53.0,psfc,ta,gamma,grid%xb%tgrn(i,j),tpw,      &
            htpw,grid%xb%speed(i,j),alw,zcld,               &
            ! grid%xb%tb22v(i,j),dum1,                        &
            ADJ_psfc,ADJ_ta,ADJ_gamma,ADJ_sst,         &
            ADJ_tpw,ADJ_htpw,ADJ_speed,ADJ_alw,        &
            ADJ_zcld,grid%xa%tb22v(i,j),dum2              )

         call da_tb_adj(3,53.0,psfc,ta,gamma,grid%xb%tgrn(i,j),tpw,      &
            htpw,grid%xb%speed(i,j),alw,zcld,               &
            ! grid%xb%tb37v(i,j),grid%xb%tb37h(i,j),               &
            ADJ_psfc,ADJ_ta,ADJ_gamma,ADJ_sst,         &
            ADJ_tpw,ADJ_htpw,ADJ_speed,ADJ_alw,        &
            ADJ_zcld,grid%xa%tb37v(i,j),grid%xa%tb37h(i,j)    )

         call da_tb_adj(4,53.0,psfc,ta,gamma,grid%xb%tgrn(i,j),tpw,      &
            htpw,grid%xb%speed(i,j),alw,zcld,               &
            ! grid%xb%tb85v(i,j),grid%xb%tb85h(i,j),               &
            ADJ_psfc,ADJ_ta,ADJ_gamma,ADJ_sst,         &
            ADJ_tpw,ADJ_htpw,ADJ_speed,ADJ_alw,        &
            ADJ_zcld,grid%xa%tb85v(i,j),grid%xa%tb85h(i,j)    )

         ADJ_zrhom    = ADJ_htpw/tpw/1000.0
         ADJ_tpw      = ADJ_tpw - ADJ_htpw/tpw*htpw

         do k = kts,kte
            grid%xa%rho(i,j,k) = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*grid%xb%h(i,j,k)* &
               grid%xb%q(i,j,k)*ADJ_zrhom + grid%xa%rho(i,j,k)
            grid%xa%q(i,j,k)   = (grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k))*grid%xb%h(i,j,k)* &
               ADJ_zrhom*grid%xb%rho(i,j,k) + grid%xa%q(i,j,k)
         end do

         ! call da_transform_xtozrhoq_adj(grid%xb,grid%xa,i,j,zh,zf,ADJ_zrhom)

         ADJ_alw = 0.0

         grid%xa%speed(i,j)=grid%xa%speed(i,j) + ADJ_speed

         grid%xa%tpw(i,j) = grid%xa%tpw(i,j) + ADJ_tpw*10.0

         ADJ_zcld= 0
         ADJ_ta  = ADJ_ta + ADJ_gamma*0.023

         grid%xa%t(i,j,kts) = grid%xa%t(i,j,kts) + ADJ_ta* &
                   log(2.0/0.0001)/log((grid%xb%h(i,j,kts)-grid%xb%terr(i,j))/0.0001)
         ADJ_sst = 0.0

         grid%xa%psfc(i,j) = grid%xa%psfc(i,j) + ADJ_psfc*0.01 
      end do
   end do   

   if (trace_use) call da_trace_exit("da_transform_xtotb_adj") 

end subroutine da_transform_xtotb_adj


subroutine da_transform_xtoy_ssmi_rv(grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! obs. increment vector (o-b).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer :: n        ! loop counter.

   real, allocatable :: tpw(:)
   real, allocatable :: speed(:)

   if (trace_use) call da_trace_entry("da_transform_xtoy_ssmi_rv")

   ! SSMI observation operator y = H(x):

   allocate (tpw(iv%info(ssmi_rv)%n1:iv%info(ssmi_rv)%n2))
   allocate (speed(iv%info(ssmi_rv)%n1:iv%info(ssmi_rv)%n2))

   call da_interp_lin_2d (grid%xa%tpw,   iv%info(ssmi_rv), 1, tpw)
   call da_interp_lin_2d (grid%xa%speed, iv%info(ssmi_rv), 1, speed)   
   
   do n=iv%info(ssmi_rv)%n1,iv%info(ssmi_rv)%n2
      y%ssmi_rv(n)%tpw   = tpw(n)
      y%ssmi_rv(n)%speed = speed(n)      
   end do

   deallocate (tpw)
   deallocate (speed)

   if (trace_use) call da_trace_exit("da_transform_xtoy_ssmi_rv")

end subroutine da_transform_xtoy_ssmi_rv


subroutine da_transform_xtoy_ssmi_rv_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)       :: grid
   type (iv_type), intent(in)       :: iv          ! obs. inc. vector (o-b).
   type (y_type) , intent(in)       :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout)    :: jo_grad_x   ! grad_x(jo)

   integer :: n

   real, allocatable :: tpw(:)
   real, allocatable :: speed(:)

   if (trace_use) call da_trace_entry("da_transform_xtoy_ssmi_rv_adj")

   allocate (tpw(iv%info(ssmi_rv)%n1:iv%info(ssmi_rv)%n2))
   allocate (speed(iv%info(ssmi_rv)%n1:iv%info(ssmi_rv)%n2))

   do n=iv%info(ssmi_rv)%n1,iv%info(ssmi_rv)%n2
      tpw(n)   = jo_grad_y%ssmi_rv(n)%tpw
      speed(n) = jo_grad_y%ssmi_rv(n)%speed
   end do

   call da_interp_lin_2d_adj(jo_grad_x%tpw,   iv%info(ssmi_rv), 1, tpw)
   call da_interp_lin_2d_adj(jo_grad_x%speed, iv%info(ssmi_rv), 1, speed)
   
   deallocate (tpw)
   deallocate (speed)

   if (trace_use) call da_trace_exit("da_transform_xtoy_ssmi_rv_adj")

end subroutine da_transform_xtoy_ssmi_rv_adj


subroutine da_transform_xtoy_ssmi_tb(grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! obs. increment vector (o-b).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   integer                      :: n        ! loop counter.

   ! real, dimension(mix,mjy)     :: TGL_dotspeed
   ! real, dimension(mix,mjy)     :: psfc,ta,gamma,sst,htpw,speed,alw,zcld,tpw
   ! real, dimension(mix,mjy)     :: TGL_psfc,TGL_ta,TGL_gamma,TGL_sst,TGL_tpw
   ! real, dimension(mix,mjy)     :: TGL_htpw,TGL_speed,TGL_alw,TGL_zcld         

   ! real, dimension(mix,mjy)     :: tb19v,tb19h, &
   !                                  tb22v,       &
   !                                  tb37v,tb37h, &
   !                                  tb85v,tb85h

   ! real, dimension(mix,mjy)     :: TGL_tb19v,TGL_tb19h, &
   !                                  TGL_tb22v,           &
   !                                  TGL_tb37v,TGL_tb37h, &
   !                                  TGL_tb85v,TGL_tb85h

   ! real, dimension(mkz)         :: zh
   ! real, dimension(mkz+1)       :: zf
   !----------------------------------------------------------------------------

   ! Wind at 1st level at dot pivnts

   ! speed, surface wind speed (m/sec)    (0 - 30) , take 10 m wind

   ! total precipitable water in cm

   ! call da_transform_xtotpw(grid%xa, grid%xb)

   ! do j=1,mjy-1
   !    do i=1,mix-1
   !       zf(1) = grid%xb%ztop
   !       do k=1,mkz
   !          zh(k)=grid%xb%h(i,j,k)
   !          zf(k+1)=grid%xb%hf(i,j,k)
   !       end do

   !       ! surface pressure (mb) (940 -1030)

   !       psfc(i,j)     = (grid%xb%p(i,j,mkz)+grid%xb%Psac(i,j)*(1.0-grid%xb % sigmah(mkz)))*0.001
   !       TGL_psfc(i,j) =  grid%xa%p  (i,j,mkz)*0.01

   !       ! sea surface temperature (k) (273 - 303) (doesnot change) 

   !       sst(i,j)      = grid%xb%tgrn(i,j)
   !       TGL_sst(i,j)  = grid%xa%tgrn(1,1)
   !       TGL_sst(i,j)  = 0.0

   !       ! effective surface air temperature (263 - 303)

   !       ta(i,j)     = grid%xb%tgrn(i,j) + &
   !                   (grid%xb%t(i,j,mkz)-grid%xb%tgrn(i,j))*log(2.0/0.0001)/ &
   !                   log(zh(mkz)/0.0001)

   !       TGL_ta(i,j) = grid%xa%tgrn(1,1) + &
   !                   (grid%xa%t(i,j,mkz)-grid%xa%tgrn(1,1))*log(2.0/0.0001)/ &
   !                   log(zh(mkz)/0.0001)

   !       TGL_ta(i,j) = (grid%xa%t(i,j,mkz)-0.)*log(2.0/0.0001)/ &
   !                   log(zh(mkz)/0.0001)

   !       ! gamma is an emperical formula and zcld is given for now

   !       gamma(i,j) = (ta(i,j)-270)*0.023 + 5.03  ! effective lapse rate   (km) (4.0 - 6.5)
   !       zcld(i,j)  = 1                           ! effective cloud height (km)
                                                     ! = 1 if no cloud infomation
   !       TGL_gamma(i,j) = TGL_ta(i,j)*0.023
   !       TGL_zcld(i,j)  = 0.0

   !       ! total precipitable water in (kg/m**2) (0 - 70)

   !       tpw(i,j)     = grid%xb%tpw(i,j)*10.0
   !       TGL_tpw(i,j) = grid%xa%tpw(i,j)*10.0

   !       ! Column liquid water (kg/m**2)   (0 - 0.5) (no data now. So, do it later.)

   !       alw(i,j)     = 0.0
   !       TGL_alw(i,j) = 0.0

   !       ! Column height weighted mivsture density on the grid locally 

   !       call da_transform_xtozrhoq(grid%xb, i, j, zh, zf, zrhom)
   !       call da_transform_xtozrhoq_lin(grid%xb, grid%xa, i, j, zh, zf, TGL_zrhom)

   !       ! Column mivsture density on the grid locally

   !       htpw(i,j)     = zrhom/tpw(i,j)/1000.0
   !       TGL_htpw(i,j) = TGL_zrhom/tpw(i,j)/1000.0 &
   !                   - TGL_tpw(i,j)/tpw(i,j)*htpw(i,j)

   !    end do
   ! end do

   ! do j=1,mjy-1
   !    do i=1,mix-1

   !       call tgl_tb(1,53.,psfc(i,j),ta(i,j),gamma(i,j),sst(i,j),tpw(i,j),  &
   !          htpw(i,j),speed(i,j),alw(i,j),zcld(i,j),               &
   !          grid%xb%tb19v(i,j),grid%xb%tb19h(i,j),                           &
   !          TGL_psfc(i,j),TGL_ta(i,j),TGL_gamma(i,j),TGL_sst(i,j), &
   !          TGL_tpw(i,j),TGL_htpw(i,j),TGL_speed(i,j),TGL_alw(i,j),&
   !          TGL_zcld(i,j),TGL_tb19v(i,j),TGL_tb19h(i,j)           )

   !       call TGL_tb(2,53.,psfc(i,j),ta(i,j),gamma(i,j),sst(i,j),tpw(i,j),  &
   !          htpw(i,j),speed(i,j),alw(i,j),zcld(i,j),               &
   !          grid%xb%tb22v(i,j),dum1,                                    &
   !          TGL_psfc(i,j),TGL_ta(i,j),TGL_gamma(i,j),TGL_sst(i,j), &
   !          TGL_tpw(i,j),TGL_htpw(i,j),TGL_speed(i,j),TGL_alw(i,j),&
   !          TGL_zcld(i,j),TGL_tb22v(i,j),dum2                     )

   !       call TGL_tb(3,53.,psfc(i,j),ta(i,j),gamma(i,j),sst(i,j),tpw(i,j),  &
   !          htpw(i,j),speed(i,j),alw(i,j),zcld(i,j),               &
   !          grid%xb%tb37v(i,j),grid%xb%tb37h(i,j),                           &
   !          TGL_psfc(i,j),TGL_ta(i,j),TGL_gamma(i,j),TGL_sst(i,j), &
   !          TGL_tpw(i,j),TGL_htpw(i,j),TGL_speed(i,j),TGL_alw(i,j),&
   !          TGL_zcld(i,j),TGL_tb37v(i,j),TGL_tb37h(i,j)           )

   !       call TGL_tb(4,53.,psfc(i,j),ta(i,j),gamma(i,j),sst(i,j),tpw(i,j),  &
   !          htpw(i,j),speed(i,j),alw(i,j),zcld(i,j),               &
   !          grid%xb%tb85v(i,j),grid%xb%tb85h(i,j),                           &
   !          TGL_psfc(i,j),TGL_ta(i,j),TGL_gamma(i,j),TGL_sst(i,j), &
   !          TGL_tpw(i,j),TGL_htpw(i,j),TGL_speed(i,j),TGL_alw(i,j),&
   !          TGL_zcld(i,j),TGL_tb85v(i,j),TGL_tb85h(i,j)           )
   !    end do
   ! end do

   real, allocatable :: tb19v(:)
   real, allocatable :: tb19h(:)
   real, allocatable :: tb22v(:)
   real, allocatable :: tb37v(:)
   real, allocatable :: tb37h(:)
   real, allocatable :: tb85v(:)
   real, allocatable :: tb85h(:)

   if (trace_use) call da_trace_entry("da_transform_xtoy_ssmi_tb")

   allocate(tb19v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb19h(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb22v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb37v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb37h(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb85v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb85h(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))

   call da_interp_lin_2d (grid%xa%tb19v, iv%info(ssmi_tb), 1, tb19v)
   call da_interp_lin_2d (grid%xa%tb19h, iv%info(ssmi_tb), 1, tb19h)
   call da_interp_lin_2d (grid%xa%tb22v, iv%info(ssmi_tb), 1, tb22v)
   call da_interp_lin_2d (grid%xa%tb37v, iv%info(ssmi_tb), 1, tb37v)
   call da_interp_lin_2d (grid%xa%tb37h, iv%info(ssmi_tb), 1, tb37h)
   call da_interp_lin_2d (grid%xa%tb85v, iv%info(ssmi_tb), 1, tb85v)
   call da_interp_lin_2d (grid%xa%tb85h, iv%info(ssmi_tb), 1, tb85h)

   do n=iv%info(ssmi_tb)%n1,iv%info(ssmi_tb)%n2
      y%ssmi_tb(n)%tb19v = tb19v(n)
      y%ssmi_tb(n)%tb19h = tb19h(n)
      y%ssmi_tb(n)%tb22v = tb22v(n)
      y%ssmi_tb(n)%tb37v = tb37v(n)
      y%ssmi_tb(n)%tb37h = tb37h(n)
      y%ssmi_tb(n)%tb85v = tb85v(n)
      y%ssmi_tb(n)%tb85h = tb85h(n)
   end do

   deallocate(tb19v)
   deallocate(tb19h)
   deallocate(tb22v)
   deallocate(tb37v)
   deallocate(tb37h)
   deallocate(tb85v)
   deallocate(tb85h)
    
   ! WHY?
   !          y%ssmi_tb(n)%tb19v = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb19v(i,j ),TGL_tb19v(i+1,j ), &
   !              TGL_tb19v(i,j+1),TGL_tb19v(i+1,j+1) )

   !          y%ssmi_tb(n)%tb19h = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb19h(i,j ),TGL_tb19h(i+1,j ), &
   !              TGL_tb19h(i,j+1),TGL_tb19h(i+1,j+1) )

   !          y%ssmi_tb(n)%tb22v = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb22v(i,j ),TGL_tb22v(i+1,j ), &
   !              TGL_tb22v(i,j+1),TGL_tb22v(i+1,j+1) )

   !          y%ssmi_tb(n)%tb37v = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb37v(i,j ),TGL_tb37v(i+1,j ), &
   !              TGL_tb37v(i,j+1),TGL_tb37v(i+1,j+1) )

   !          y%ssmi_tb(n)%tb37h = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb37h(i,j ),TGL_tb37h(i+1,j ), &
   !              TGL_tb37h(i,j+1),TGL_tb37h(i+1,j+1) )

   !          y%ssmi_tb(n)%tb85v = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb85v(i,j ),TGL_tb85v(i+1,j ), &
   !              TGL_tb85v(i,j+1),TGL_tb85v(i+1,j+1) )

   !          y%ssmi_tb(n)%tb85h = hor_interp(dxm,dx,dym,dy, &
   !              TGL_tb85h(i,j ),TGL_tb85h(i+1,j ), &
   !              TGL_tb85h(i,j+1),TGL_tb85h(i+1,j+1) )
   !       end if
   !    end if

   ! end do

   if (trace_use) call da_trace_exit("da_transform_xtoy_ssmi_tb")

end subroutine da_transform_xtoy_ssmi_tb


subroutine da_transform_xtoy_ssmi_tb_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer                       :: n                ! loop counter.

   real, allocatable :: tb19v(:)
   real, allocatable :: tb19h(:)
   real, allocatable :: tb22v(:)
   real, allocatable :: tb37v(:)
   real, allocatable :: tb37h(:)
   real, allocatable :: tb85v(:)
   real, allocatable :: tb85h(:)

   if (trace_use) call da_trace_entry("da_transform_xtoy_ssmi_tb_adj")

   allocate(tb19v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb19h(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb22v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb37v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb37h(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb85v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))
   allocate(tb85v(iv%info(ssmi_tb)%n1:iv%info(ssmi_tb)%n2))

   do n=iv%info(ssmi_tb)%n1,iv%info(ssmi_tb)%n2
      tb19v(n) = jo_grad_y%ssmi_tb(n)%tb19v
      tb19h(n) = jo_grad_y%ssmi_tb(n)%tb19h
      tb22v(n) = jo_grad_y%ssmi_tb(n)%tb22v
      tb37v(n) = jo_grad_y%ssmi_tb(n)%tb37v
      tb37h(n) = jo_grad_y%ssmi_tb(n)%tb37h
      tb85v(n) = jo_grad_y%ssmi_tb(n)%tb85v
      tb85h(n) = jo_grad_y%ssmi_tb(n)%tb85h
   end do

   call da_interp_lin_2d_adj(jo_grad_x%tb19v, iv%info(ssmi_tb), 1, tb19v)
   call da_interp_lin_2d_adj(jo_grad_x%tb19h, iv%info(ssmi_tb), 1, tb19h)
   call da_interp_lin_2d_adj(jo_grad_x%tb22v, iv%info(ssmi_tb), 1, tb22v)
   call da_interp_lin_2d_adj(jo_grad_x%tb37v, iv%info(ssmi_tb), 1, tb37v)
   call da_interp_lin_2d_adj(jo_grad_x%tb37h, iv%info(ssmi_tb), 1, tb37h)
   call da_interp_lin_2d_adj(jo_grad_x%tb85v, iv%info(ssmi_tb), 1, tb85v)
   call da_interp_lin_2d_adj(jo_grad_x%tb85h, iv%info(ssmi_tb), 1, tb85h)
  
   deallocate(tb19v)
   deallocate(tb19h)
   deallocate(tb22v)
   deallocate(tb37v)
   deallocate(tb37h)
   deallocate(tb85v)
   deallocate(tb85v)

   if (trace_use) call da_trace_exit("da_transform_xtoy_ssmi_tb_adj")

end subroutine da_transform_xtoy_ssmi_tb_adj


subroutine da_transform_xtozrhoq(xb, i, j, zh, zf, zrhom)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   integer,        intent(in)      :: i, j
   type (xb_type), intent(in)      :: xb         ! first guess state.
   real,           intent(out)     :: zrhom      ! gridded height weighted moisture
   real,           intent(in)      :: zh(mkz)
   real,           intent(in)      :: zf(mkz+1)

   integer                         :: k

   if (trace_use) call da_trace_entry("da_transform_xtozrhoq")
   
   zrhom = 0.0

   do k = 1,mkz
      zrhom = (zf(k)-zf(k+1))*zh(k)*(xb%q(i,j,k)*xb%rho(i,j,k))+zrhom
   end do

   if (trace_use) call da_trace_exit("da_transform_xtozrhoq")
 
end subroutine da_transform_xtozrhoq


subroutine da_transform_xtozrhoq_lin(xb, xa, i, j, zh, zf, tgl_zrhom)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   integer,        intent(in)      :: i, j
   type (xb_type), intent(in)      :: xb         ! first guess state.
   type (x_type) , intent(in)      :: xa         ! increment
   real,           intent(out)     :: TGL_zrhom  ! gridded height weighted moisture
   real,           intent(in)      :: zh(mkz)
   real,           intent(in)      :: zf(mkz+1)

   integer                         :: k

   if (trace_use) call da_trace_entry("da_transform_xtozrhoq_lin")

   TGL_zrhom = 0.0

   do k = 1,mkz
      TGL_zrhom = (zf(k)-zf(k+1))*zh(k)*(xb%q(i,j,k)*xa%rho(i,j,k) + &
                                         xa%q(i,j,k)*xb%rho(i,j,k)   &
                                       ) + TGL_zrhom
   end do

   if (trace_use) call da_trace_exit("da_transform_xtozrhoq_lin")
 
end subroutine da_transform_xtozrhoq_lin


subroutine da_transform_xtozrhoq_adj(grid, i, j, zh, zf, adj_zrhom)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   type (domain), intent(inout) :: grid
   integer,       intent(in)    :: i, j
   real,          intent(in)    :: ADJ_zrhom  ! gridded height weighted moisture
   real,          intent(in)    :: zh(mkz)
   real,          intent(in)    :: zf(mkz+1)

   integer                      :: k

   if (trace_use) call da_trace_entry("da_transform_xtozrhoq_adj")

   do k = 1,mkz
      grid%xa%rho(i,j,k) = (zf(k)-zf(k+1))*zh(k)*grid%xb%q(i,j,k)*ADJ_zrhom   + grid%xa%rho(i,j,k)
      grid%xa%q(i,j,k)   = (zf(k)-zf(k+1))*zh(k)*ADJ_zrhom*grid%xb%rho(i,j,k) + grid%xa%q(i,j,k)
   end do

   if (trace_use) call da_trace_exit("da_transform_xtozrhoq_adj")
 
end subroutine da_transform_xtozrhoq_adj


subroutine da_jo_and_grady_ssmt1(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(in)    :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_ssmt1")

   jo % ssmt1_t = 0.0

   do n=1, iv%info(ssmt1)%nlocal
      do k=1, iv%info(ssmt1)%levels(n)
         jo_grad_y%ssmt1(n)%t(k) = -re%ssmt1(n)%t(k) / (iv%ssmt1(n)%t(k)%error * iv%ssmt1(n)%t(k)%error)

         jo % ssmt1_t = jo % ssmt1_t - re%ssmt1(n)%t(k) * jo_grad_y%ssmt1(n)%t(k)
      end do
   end do

   jo % ssmt1_t = 0.5 * jo % ssmt1_t

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_ssmt1")

end subroutine da_jo_and_grady_ssmt1


subroutine da_jo_and_grady_ssmt2(iv, re, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(in)    :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   type (jo_type), intent(inout) :: jo          ! Obs cost function.

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_jo_and_grady_ssmt2")

   jo % ssmt2_rh = 0.0
   
   do n=1, iv%info(ssmt2)%nlocal
      do k=1, iv%info(ssmt2)%levels(n)
         jo_grad_y%ssmt2(n)%rh(k) = -re%ssmt2(n)%rh(k) / &
            (iv%ssmt2(n)%rh(k)%error * iv%ssmt2(n)%rh(k)%error)

         jo % ssmt2_rh = jo % ssmt2_rh - re%ssmt2(n)%rh(k) * jo_grad_y%ssmt2(n)%rh(k)
      end do
   end do

   jo % ssmt2_rh = 0.5 * jo % ssmt2_rh

   if (trace_use_dull) call da_trace_exit("da_jo_and_grady_ssmt2")

end subroutine da_jo_and_grady_ssmt2


subroutine da_residual_ssmt1(iv, y, re, np_missing, np_bad_data,np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual structure.

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type) :: n_obs_bad
   integer              :: n, k

   if (trace_use_dull) call da_trace_entry("da_residual_ssmt1")

   n_obs_bad % t % num = number_type(0, 0, 0)

   do n=1, iv%info(ssmt1)%nlocal
      do k=1, iv%info(ssmt1)%levels(n)
         np_available = np_available + 1
         re%ssmt1(n)%t(k) = da_residual(n, k, y%ssmt1(n)%t(k), iv%ssmt1(n)%t(k), n_obs_bad % t)
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % t % num % miss
   np_bad_data = np_bad_data + n_obs_bad % t % num % bad
   np_obs_used = np_obs_used + n_obs_bad % t % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_ssmt1")

end subroutine da_residual_ssmt1


subroutine da_residual_ssmt2(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv     ! Innovation vector (O-B).
   type (y_type) , intent(in)    :: y      ! y = H (xa)
   type (y_type) , intent(inout) :: re     ! Residual structure.

   integer       , intent(inout) :: np_available
   integer       , intent(inout) :: np_obs_used
   integer       , intent(inout) :: np_missing
   integer       , intent(inout) :: np_bad_data

   type (bad_data_type) :: n_obs_bad
   integer              :: n, k

   if (trace_use_dull) call da_trace_entry("da_residual_ssmt2")

   n_obs_bad % rh % num = number_type(0, 0, 0)

   do n=1, iv%info(ssmt2)%nlocal
      do k=1, iv%info(ssmt2)%levels(n)
         np_available = np_available + 1
         re%ssmt2(n)%rh(k) = da_residual(n, k, y%ssmt2(n)%rh(k), iv%ssmt2(n)%rh(k), n_obs_bad % rh)
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % rh % num % miss
   np_bad_data = np_bad_data + n_obs_bad % rh % num % bad
   np_obs_used = np_obs_used + n_obs_bad % rh % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_ssmt2")

end subroutine da_residual_ssmt2


subroutine da_check_max_iv_ssmi_rv(iv, it, num_qcstat_conv)              

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it ! Outer loop 
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)

   logical :: failed
   integer :: n
   
   if (trace_use) call da_trace_entry("da_check_max_iv_ssmi_rv")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(ssmi_rv)%n1,iv%info(ssmi_rv)%n2

      failed=.false.
      if( iv%ssmi_rv(n)%tpw%qc >= obs_qc_pointer ) then 
      call da_max_error_qc (it, iv%info(ssmi_rv), n, iv%ssmi_rv(n)%tpw, max_error_pw, failed)
      if( iv%info(ssmi_rv)%proc_domain(1,n) ) then
                 num_qcstat_conv(1,ssmi_rv,7,1) = num_qcstat_conv(1,ssmi_rv,7,1) + 1
      if(failed) then
        num_qcstat_conv(2,ssmi_rv,7,1) = num_qcstat_conv(2,ssmi_rv,7,1) + 1
        write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
        'ssmi_rv',ob_vars(7),iv%info(ssmi_rv)%lat(1,n),iv%info(ssmi_rv)%lon(1,n),'1013.25'                  
      end if
      end if
      end if

      failed=.false.
      if( iv%ssmi_rv(n)%speed%qc >= obs_qc_pointer ) then
      call da_max_error_qc (it, iv%info(ssmi_rv), n, iv%ssmi_rv(n)%speed, max_error_uv, failed)
      if( iv%info(ssmi_rv)%proc_domain(1,n) ) then
                 num_qcstat_conv(1,ssmi_rv,6,1) = num_qcstat_conv(1,ssmi_rv,6,1) + 1
      if(failed)then
         num_qcstat_conv(2,ssmi_rv,6,1) = num_qcstat_conv(2,ssmi_rv,6,1) + 1
        write(qcstat_conv_unit,'(2x,a10,2x,a4,2f12.2,a12)')&
        'ssmi_rv',ob_vars(6),iv%info(ssmi_rv)%lat(1,n),iv%info(ssmi_rv)%lon(1,n),'1013.25'                  
      endif
      end if
      end if

   end do

   if (trace_use) call da_trace_exit("da_check_max_iv_ssmi_rv")

end subroutine da_check_max_iv_ssmi_rv
subroutine da_check_max_iv_ssmi_tb(iv, it)  

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it

   logical :: failed
   integer :: n

   if (trace_use) call da_trace_entry("da_check_max_iv_ssmi_tb")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(ssmi_tb)%n1, iv%info(ssmi_tb)%n2
      ! Tb19h

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb19h, max_error_tb, failed)

      ! Tb19v

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb19v, max_error_tb, failed)

      ! Tb22v

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb22v, max_error_tb, failed)

      ! Tb37h

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb37h, max_error_tb, failed)

      ! Tb37v

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb37v, max_error_tb, failed)

      ! Tb85h

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb85h, max_error_tb, failed)

      ! Tb85v

      call da_max_error_qc(it, iv%info(ssmi_tb), n, iv%ssmi_tb(n)%tb85v,max_error_tb, failed)
   end do

   if (trace_use) call da_trace_exit("da_check_max_iv_ssmi_tb")

end subroutine da_check_max_iv_ssmi_tb


subroutine da_check_max_iv_ssmt1(iv, it, num_qcstat_conv)              

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    ::  it      ! outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)



   integer ::  k,n, ipr
   logical :: failed

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_ssmt1")

   !---------------------------------------------------------------------------
   ! [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(ssmt1)%n1,iv%info(ssmt1)%n2
      do k = 1, iv%info(ssmt1)%levels(n)
         call da_get_print_lvl(iv%ssmt1(n)%p(k),ipr)

        failed= .false.
        if( iv%ssmt1(n)%t(k)%qc >= obs_qc_pointer ) then 
        call da_max_error_qc(it, iv%info(ssmt1), n, iv%ssmt1(n)%t(k), max_error_t, failed)
        if( iv%info(ssmt1)%proc_domain(k,n) ) then
         num_qcstat_conv(1,ssmt1,3,ipr) = num_qcstat_conv(1,ssmt1,3,ipr) + 1
        if(failed) then
        num_qcstat_conv(2,ssmt1,3,ipr) = num_qcstat_conv(2,ssmt1,3,ipr) + 1
        write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
        'ssmt1',ob_vars(3),iv%info(ssmt1)%lat(k,n),iv%info(ssmt1)%lon(k,n),0.01*iv%ssmt1(n)%p(k)
        endif
        endif
        endif

      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_ssmt1")

end subroutine da_check_max_iv_ssmt1


subroutine da_check_max_iv_ssmt2(iv, it, num_qcstat_conv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Update:
   !    Removed Outerloop check as it is done in da_get_innov
   !    Author: Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 07/12/2009
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: it     ! outer iteration
   integer,       intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: k,n, ipr
   logical :: failed

   if (trace_use_dull) call da_trace_entry("da_check_max_iv_ssmt2")


   !---------------------------------------------------------------------------
   !  [1.0] Perform maximum innovation vector check:
   !---------------------------------------------------------------------------

   do n=iv%info(ssmt2)%n1,iv%info(ssmt2)%n2
      do k = 1, iv%info(ssmt2)%levels(n)
        call da_get_print_lvl(iv%ssmt2(n)%p(k),ipr)

        failed = .false.
        if( iv%ssmt2(n)%rh(k)%qc >= obs_qc_pointer ) then
        call da_max_error_qc (it, iv%info(ssmt2), n, iv%ssmt2(n)%rh(k), max_error_q, failed)
        if( iv%info(ssmt2)%proc_domain(k,n) ) then
        num_qcstat_conv(1,ssmt2,4,ipr) = num_qcstat_conv(1,ssmt2,4,ipr) + 1
        if(failed)then
         num_qcstat_conv(2,ssmt2,4,ipr) = num_qcstat_conv(2,ssmt2,4,ipr) + 1
        write(qcstat_conv_unit,'(2x,a10,2x,a4,3f12.2)')&
        'ssmt2',ob_vars(4),iv%info(ssmt2)%lat(k,n),iv%info(ssmt2)%lon(k,n),0.01*iv%ssmt2(n)%p(k)
        endif
        endif
        endif
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_check_max_iv_ssmt2")

end subroutine da_check_max_iv_ssmt2

subroutine da_get_innov_vector_ssmi_rv (it,num_qcstat_conv, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   integer,       intent(in)    :: it         ! External iteration.
   type(domain),  intent(in)    :: grid       ! first guess state.
   type(y_type),  intent(in)    :: ob         ! Observation structure.
   type(iv_type), intent(inout) :: iv         ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n           ! Loop counter.
   integer :: i, j        ! Index dimension.
   real    :: dx, dxm     ! Interpolation weights.
   real    :: dy, dym     ! Interpolation weights.
   real    :: model_tpw   ! Model value tpw   at oblocation
   real    :: model_speed ! Model value speed at oblocation
   integer :: itpw, itpwf, ispeed, ispeedf 
   
   if (trace_use) call da_trace_entry("da_get_innov_vector_ssmi_rv")

   if ( it > 1 ) then
      do n=iv%info(ssmi_rv)%n1,iv%info(ssmi_rv)%n2
         if ( iv % ssmi_rv(n) % speed % qc == fails_error_max) iv % ssmi_rv(n) % speed % qc = 0
         if ( iv % ssmi_rv(n) % tpw % qc == fails_error_max) iv % ssmi_rv(n) % tpw % qc = 0
      end do
   end if
   
   do n=iv%info(ssmi_rv)%n1,iv%info(ssmi_rv)%n2

      ! compute innovation vector
      ! =========================

      ! Obs coordinates on model grid

      ! TPW

      i   = iv%info(ssmi_rv)%i(1,n)
      j   = iv%info(ssmi_rv)%j(1,n)
      dx  = iv%info(ssmi_rv)%dx(1,n)
      dy  = iv%info(ssmi_rv)%dy(1,n)
      dxm = iv%info(ssmi_rv)%dxm(1,n)
      dym = iv%info(ssmi_rv)%dym(1,n)

      iv % ssmi_rv(n) % tpw % inv  = 0.0
      if (abs(ob % ssmi_rv(n) % tpw - missing_r) > 1.0  .and. iv % ssmi_rv(n) % tpw % qc >= obs_qc_pointer) then
         model_tpw = dym*(dxm*grid%xb%tpw(i,j) + dx*grid%xb%tpw(i+1,j)) + dy *(dxm*grid%xb%tpw(i,j+1) + dx*grid%xb%tpw(i+1,j+1))
         iv % ssmi_rv(n) % tpw % inv = ob % ssmi_rv(n) % tpw - model_tpw
      end if

      ! surface wind speed

      iv % ssmi_rv(n) % speed % inv  = 0.0
      if (abs(ob % ssmi_rv(n) % speed - missing_r) > 1.0 .and.     &
          iv % ssmi_rv(n) % speed % qc >= obs_qc_pointer) then

         model_speed = dym*(dxm*grid%xb%speed(i,j ) + dx*grid%xb%speed(i+1,j )) &
            + dy *(dxm*grid%xb%speed(i,j+1) + dx*grid%xb%speed(i+1,j+1))
         iv % ssmi_rv(n) % speed % inv = ob % ssmi_rv(n) % speed - model_speed
      end if
   end do

   !------------------------------------------------------------------
   ! Perform optional maximum error check:
   !------------------------------------------------------------------

   if ( check_max_iv ) &
      call da_check_max_iv_ssmi_rv(iv, it, num_qcstat_conv)           

   if (trace_use) call da_trace_exit("da_get_innov_vector_ssmi_rv")

end subroutine da_get_innov_vector_ssmi_rv


subroutine da_get_innov_vector_ssmi_tb (it, grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,       intent(in)    :: it         ! External iteration.
   type(domain),  intent(in)    :: grid       ! first guess state.
   type(y_type),  intent(in)    :: ob         ! Observation structure.
   type(iv_type), intent(inout) :: iv         ! O-B structure.

   integer :: n           ! Loop counter.
   integer :: i, j        ! Index dimension.
   real    :: dx, dxm     ! Interpolation weights.
   real    :: dy, dym     ! Interpolation weights.
   real    :: model_tb19h ! Model value tb19h at oblocation.
   real    :: model_tb19v ! Model value tb19v at oblocation.
   real    :: model_tb22v ! Model value tb22v at oblocation.
   real    :: model_tb37h ! Model value tb37h at oblocation.
   real    :: model_tb37v ! Model value tb37v at oblocation.
   real    :: model_tb85h ! Model value tb85h at oblocation.
   real    :: model_tb85v ! Model value tb85v at oblocation.

   if (trace_use) call da_trace_entry("da_get_innov_vector_ssmi_tb")

   if ( it > 1 ) then
      do n=iv%info(ssmi_tb)%n1,iv%info(ssmi_tb)%n2
         if(iv % ssmi_tb(n) % tb19v % qc == fails_error_max ) iv % ssmi_tb(n) % tb19v % qc = 0
         if(iv % ssmi_tb(n) % tb19h % qc == fails_error_max ) iv % ssmi_tb(n) % tb19h % qc = 0
         if(iv % ssmi_tb(n) % tb22v % qc == fails_error_max ) iv % ssmi_tb(n) % tb22v % qc = 0
         if(iv % ssmi_tb(n) % tb37v % qc == fails_error_max ) iv % ssmi_tb(n) % tb37v % qc = 0
         if(iv % ssmi_tb(n) % tb37h % qc == fails_error_max ) iv % ssmi_tb(n) % tb37h % qc = 0
         if(iv % ssmi_tb(n) % tb85v % qc == fails_error_max ) iv % ssmi_tb(n) % tb85v % qc = 0
         if(iv % ssmi_tb(n) % tb85h % qc == fails_error_max ) iv % ssmi_tb(n) % tb85h % qc = 0
      end do
   end if

   do n=iv%info(ssmi_tb)%n1,iv%info(ssmi_tb)%n2
      ! compute innovation vector
      ! =========================

      !  Obs coordinates on model grid

      ! TB

      i   = iv%info(ssmi_tb)%i(1,n)
      j   = iv%info(ssmi_tb)%j(1,n)
      dx  = iv%info(ssmi_tb)%dx(1,n)
      dy  = iv%info(ssmi_tb)%dy(1,n)
      dxm = iv%info(ssmi_tb)%dxm(1,n)
      dym = iv%info(ssmi_tb)%dym(1,n)

      ! Tb19h

      if (abs(ob % ssmi_tb(n) % tb19h - missing_r) > 1.0) then
         model_tb19h = dym*(dxm*grid%xb%tb19h(i,j)   + dx*grid%xb%tb19h(i+1,j)) &
            + dy *(dxm*grid%xb%tb19h(i,j+1) + dx*grid%xb%tb19h(i+1,j+1))
         iv % ssmi_tb(n) % tb19h % inv = ob % ssmi_tb(n) % tb19h - &
            model_tb19h
      else
         iv % ssmi_tb(n) % tb19h % inv = 0.0
      end if

      ! Tb19v

      if (abs(ob % ssmi_tb(n) % tb19v - missing_r) > 1.0) then
         model_tb19v = dym*(dxm*grid%xb%tb19v(i,j)   + dx *grid%xb%tb19v(i+1,j)) &
            + dy *(dxm*grid%xb%tb19v(i,j+1) + dx *grid%xb%tb19v(i+1,j+1))
         iv % ssmi_tb(n) % tb19v % inv = ob % ssmi_tb(n) % tb19v - &
            model_tb19v
      else
         iv % ssmi_tb(n) % tb19v % inv = 0.0
      end if

     ! Tb19v

      if (abs(ob % ssmi_tb(n) % tb22v - missing_r) > 1.0) then
         model_tb22v = dym*(dxm*grid%xb%tb22v(i,j) + dx *grid%xb%tb22v(i+1,j)) &
            + dy *(dxm*grid%xb%tb22v(i,j+1) + dx *grid%xb%tb22v(i+1,j+1))
         iv % ssmi_tb(n) % tb22v % inv = ob % ssmi_tb(n) % tb22v - &
            model_tb22v
      else
         iv % ssmi_tb(n) % tb22v % inv = 0.0
      end if

      ! Tb37h

      if (abs(ob % ssmi_tb(n) % tb37h - missing_r) > 1.0) then
         model_tb37h = dym*(dxm*grid%xb%tb37h(i,j)  + dx *grid%xb%tb37h(i+1,j)) &
            + dy *(dxm*grid%xb%tb37h(i,j+1) + dx *grid%xb%tb37h(i+1,j+1))
         iv % ssmi_tb(n) % tb37h % inv = ob % ssmi_tb(n) % tb37h - &
            model_tb37h
      else
         iv % ssmi_tb(n) % tb37h % inv = 0.0
      end if

      ! Tb37v

      if (abs(ob % ssmi_tb(n) % tb37v - missing_r) > 1.0) then
         model_tb37v = dym*(dxm*grid%xb%tb37v(i,j)  + dx *grid%xb%tb37v(i+1,j)) &
            + dy *(dxm*grid%xb%tb37v(i,j+1) + dx *grid%xb%tb37v(i+1,j+1))
         iv % ssmi_tb(n) % tb37v % inv = ob % ssmi_tb(n) % tb37v - &
            model_tb37v
      else
         iv % ssmi_tb(n) % tb37v % inv = 0.0
      end if

      ! Tb85h

      if (abs(ob % ssmi_tb(n) % tb85h - missing_r) > 1.0) then
         model_tb85h = dym*(dxm*grid%xb%tb85h(i,j) + dx *grid%xb%tb85h(i+1,j)) &
            + dy *(dxm*grid%xb%tb85h(i,j+1) + dx *grid%xb%tb85h(i+1,j+1))
         iv % ssmi_tb(n) % tb85h % inv = ob % ssmi_tb(n) % tb85h - &
            model_tb85h
      else
         iv % ssmi_tb(n) % tb85h % inv = 0.0
      end if

      ! Tb85v

      if (abs(ob % ssmi_tb(n) % tb85v - missing_r) > 1.0) then
         model_tb85v = dym*(dxm*grid%xb%tb85v(i,j) + dx *grid%xb%tb85v(i+1,j)) &
            + dy *(dxm*grid%xb%tb85v(i,j+1) + dx *grid%xb%tb85v(i+1,j+1))
         iv % ssmi_tb(n) % tb85v % inv = ob % ssmi_tb(n) % tb85v -  &
            model_tb85v
      else
         iv % ssmi_tb(n) % tb85v % inv = 0.0
      end if
   end do

   !----------------------------------------------------------------
   !     Perform optional maximum error check:
   !----------------------------------------------------------------

   if (check_max_iv) call da_check_max_iv_ssmi_tb(iv, it)  
   
   if (trace_use) call da_trace_exit("da_get_innov_vector_ssmi_tb")

end subroutine da_get_innov_vector_ssmi_tb


subroutine da_get_innov_vector_ssmt1(it,num_qcstat_conv,grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it       ! External iteration.
   type(domain),     intent(in)    :: grid       ! first guess state.
   type(y_type),     intent(inout) :: ob       ! Observation structure.
   type(iv_type),    intent(inout) :: iv       ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)


   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs ! Number of obs levels.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real, allocatable :: model_t(:,:)  ! Model value t at ob location.

   real   :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real   :: v_p(kms:kme)      ! Model value p at ob hor. location.

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_ssmt1")

   allocate (model_t(1:iv%info(ssmt1)%max_lev,iv%info(ssmt1)%n1:iv%info(ssmt1)%n2))
   model_t(:,:) = 0.0

   if ( it > 1 )then
      do n=iv%info(ssmt1)%n1,iv%info(ssmt1)%n2
         do k = 1, iv%info(ssmt1)%levels(n)
            if(iv % ssmt1(n) % t(k) % qc == fails_error_max)iv % ssmt1(n) % t(k) % qc = 0
         end do
      end do
   end if

   do n=iv%info(ssmt1)%n1,iv%info(ssmt1)%n2

      num_levs = iv%info(ssmt1)%levels(n)

      if (num_levs < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      i   = iv%info(ssmt1)%i(1,n)
      j   = iv%info(ssmt1)%j(1,n)
      dx  = iv%info(ssmt1)%dx(1,n)
      dy  = iv%info(ssmt1)%dy(1,n)
      dxm = iv%info(ssmt1)%dxm(1,n)
      dym = iv%info(ssmt1)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j,k)+dx*grid%xb%h(i+1,j,k)) + dy*(dxm*grid%xb%h(i,j+1,k)+dx*grid%xb%h(i+1,j+1,k))
         v_p(k) = dym*(dxm*grid%xb%p(i,j,k)+dx*grid%xb%p(i+1,j,k)) + dy*(dxm*grid%xb%p(i,j+1,k)+dx*grid%xb%p(i+1,j+1,k))
      end do

      num_levs=0

      do k=1, iv%info(ssmt1)%levels(n)
         if (iv % ssmt1(n) % h(k) > 0.0) then
            call da_to_zk(iv % ssmt1(n) % h(k), v_h, v_interp_h, iv%info(ssmt1)%zk(k,n))
         else if (iv % ssmt1(n) % p(k) > 1.0) then
            call da_to_zk(iv % ssmt1(n) % p(k), v_p, v_interp_p, iv%info(ssmt1)%zk(k,n))
         end if

         if ( iv%info(ssmt1)%zk(k,n) > 0.0) then
            num_levs=num_levs+1
            iv%info(ssmt1)%zk(num_levs,n)=iv%info(ssmt1)%zk(k,n)

            ob % ssmt1(n) % t(num_levs)         = ob % ssmt1(n) % t(k)
            iv % ssmt1(n) % t(num_levs) % qc    = iv % ssmt1(n) % t(k) % qc
            iv % ssmt1(n) % t(num_levs) % error = iv % ssmt1(n) % t(k) % error
         end if
      end do

      iv%info(ssmt1)%levels(n) = num_levs
   end do

   call da_convert_zk (iv%info(ssmt1))

   ! [1.2] Interpolate horizontally to ob:

   call da_interp_lin_3d (grid%xb%t, iv%info(ssmt1), model_t)

   do n=iv%info(ssmt1)%n1,iv%info(ssmt1)%n2

      !---------------------------------------------------------------------
      ! [2.0] Initialise components of innovation vector:
      !---------------------------------------------------------------------

      do k = 1, iv%info(ssmt1)%levels(n)
         iv % ssmt1(n) % t(k) % inv = 0.0

         !----------------------------------------------------------------
         ! [3.0] Interpolation:
         !----------------------------------------------------------------

         if (ob % ssmt1(n) % t(k) > missing_r .AND. &
             iv % ssmt1(n) % t(k) % qc >= obs_qc_pointer) then
            iv % ssmt1(n) % t(k) % inv = ob % ssmt1(n) % t(k) - model_t(k,n)
         end if
      end do
   end do

   !------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------

   if (check_max_iv) call da_check_max_iv_ssmt1(iv, it, num_qcstat_conv)

   deallocate (model_t)
    
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_ssmt1")

end subroutine da_get_innov_vector_ssmt1


subroutine da_get_innov_vector_ssmt2 (it, num_qcstat_conv,grid, ob, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: it       ! External iteration.
   type(domain),     intent(in)    :: grid       ! first guess state.
   type(y_type),     intent(inout) :: ob       ! Observation structure.
   type(iv_type),    intent(inout) :: iv       ! O-B structure.
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   integer :: num_levs ! Number of obs levels.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real, allocatable :: model_rh(:,:) ! Model value rh at ob location.

   real    :: v_h(kms:kme)      ! Model value h at ob hor. location.
   real    :: v_p(kms:kme)      ! Model value p at ob hor. location.

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_ssmt2")

   allocate (model_rh(1:iv%info(ssmt2)%max_lev,iv%info(ssmt2)%n1:iv%info(ssmt2)%n2))
   model_rh(:,:) = 0.0

   if ( it > 1 )then
      do n=iv%info(ssmt2)%n1,iv%info(ssmt2)%n2
         do k = 1, iv%info(ssmt1)%levels(n)
            if(iv % ssmt2(n) % rh(k) % qc == fails_error_max)iv % ssmt2(n) % rh(k) % qc = 0
         end do
      end do
   end if

   do n=iv%info(ssmt2)%n1,iv%info(ssmt2)%n2
      num_levs = iv%info(ssmt2)%levels(n)

      if (num_levs < 1) cycle

      ! [1.1] Get horizontal interpolation weights:

      i   = iv%info(ssmt2)%i(1,n)
      j   = iv%info(ssmt2)%j(1,n)
      dx  = iv%info(ssmt2)%dx(1,n)
      dy  = iv%info(ssmt2)%dy(1,n)
      dxm = iv%info(ssmt2)%dxm(1,n)
      dym = iv%info(ssmt2)%dym(1,n)

      do k=kts,kte
         v_h(k) = dym*(dxm*grid%xb%h(i,j  ,k) + dx*grid%xb%h(i+1,j  ,k)) &
                + dy *(dxm*grid%xb%h(i,j+1,k) + dx*grid%xb%h(i+1,j+1,k))
         v_p(k) = dym*(dxm*grid%xb%p(i,j  ,k) + dx*grid%xb%p(i+1,j  ,k)) &
                + dy *(dxm*grid%xb%p(i,j+1,k) + dx*grid%xb%p(i+1,j+1,k))
      end do

      num_levs=0
      do k=1, iv%info(ssmt2)%levels(n)
         if (iv % ssmt2(n) % h(k) > 0.0) then
            call da_to_zk(iv % ssmt2(n) % h(k), v_h, v_interp_h, iv%info(ssmt2)%zk(k,n))
         else if (iv % ssmt2(n) % p(k) > 1.0) then
            call da_to_zk(iv % ssmt2(n) % p(k), v_p, v_interp_p, iv%info(ssmt2)%zk(k,n))
         end if

         if (iv%info(ssmt2)%zk(k,n) > 0.0) then
            num_levs=num_levs+1
            iv%info(ssmt2)%zk(num_levs,n)=iv%info(ssmt2)%zk(k,n)

            ob % ssmt2(n) % rh(num_levs) = ob % ssmt2(n) % rh(k)

            iv % ssmt2(n) % rh(num_levs) % qc = iv % ssmt2(n) % rh(k) % qc

            iv % ssmt2(n) % rh(num_levs) % error = iv % ssmt2(n) % rh(k) % error
         end if
      end do

      iv%info(ssmt2)%levels(n) = num_levs
   end do

   call da_convert_zk (iv%info(ssmt2))

   call da_interp_lin_3d (grid%xb%rh, iv%info(ssmt2), model_rh)

   do n=iv%info(ssmt2)%n1, iv%info(ssmt2)%n2
      do k = 1, iv%info(ssmt2)%levels(n)
         iv % ssmt2(n) % rh(k) % inv = 0.0

         !-----------------------------------------------------------------
         ! [3.0] Interpolation:
         !-----------------------------------------------------------------

         if (ob % ssmt2(n) % rh(k) > missing_r .AND. &
             iv % ssmt2(n) % rh(k) % qc >= obs_qc_pointer) then
            iv % ssmt2(n) % rh(k) % inv = ob % ssmt2(n) % rh(k) - model_rh(k,n)
         end if

         ! write(122,'(2i4,i8,5f15.5)')n, k, iv%ssmt2(n)%height_qc(k), &
         ! iv%ssmt2(n)%info%lat, iv%ssmt2(n)%info%lon, &
         ! iv%ssmt2(n)%h(k), &
         ! ob % ssmt2(n) % rh(k), model_rh(k,n)

      end do
   end do

   !------------------------------------------------------------------------
   ! [5.0] Perform optional maximum error check:
   !------------------------------------------------------------------------

   if (check_max_iv) call da_check_max_iv_ssmt2(iv, it, num_qcstat_conv)           

   deallocate (model_rh)
   
   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_ssmt2")

end subroutine da_get_innov_vector_ssmt2


subroutine da_ao_stats_ssmt1 (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! OI
   type  (y_type), intent (in)    :: re            ! O-A

   type (maxmin_type) :: minimum
   type (maxmin_type) :: maximum
   integer            :: nt
   integer            :: n, k
   real               :: average, rms_err

   if (trace_use_dull) call da_trace_entry("da_ao_stats_ssmt1")

   nt = 0

   maximum = maxmin_type(-1.0E+20, 0, 0)
   minimum = maxmin_type (1.0E+20, 0, 0)
   average = 0.0
   rms_err = 0.0

   if (iv%info(ssmt1)%nlocal > 0) then
      do n=1, iv%info(ssmt1)%nlocal
         if (iv%info(ssmt1)%proc_domain(1,n)) then
            do k=1, iv%info(ssmt1)%levels(n)
               call da_stats_calculate (n, k, iv%ssmt1(n)%t(k)%qc, re%ssmt1(n)%t(k), nt, &
                  minimum, maximum, average, rms_err)
            end do
         end if    ! end if (iv%info(ssmt1)%proc_domain(1,n))
      end do
   end if

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nt)
   iv%nstats(ssmt1) = nt

   call da_proc_stats_combine(average, rms_err, minimum%value, maximum%value, &
      minimum%n, maximum%n, minimum%l, maximum%l)
   
   if (rootproc) then
      if (nt /= 0) then
          write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of O-A for SSMT1'
          call da_print_stats_ssmt1 (stats_unit, nt, minimum, maximum, average, rms_err)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_ssmt1")

end subroutine da_ao_stats_ssmt1


subroutine da_ao_stats_ssmt2 (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! OI
   type (y_type),  intent(in)    :: re            ! O-A

   type (maxmin_type) :: minimum
   type (maxmin_type) :: maximum
   integer            :: nrh
   integer            :: n, k
   real               :: average, rms_err

   if (trace_use_dull) call da_trace_entry("da_ao_stats_ssmt2")

   nrh = 0
   
   maximum = maxmin_type(-1.0E+20, 0, 0)
   minimum = maxmin_type (1.0E+20, 0, 0)
   average = 0.0
   rms_err = 0.0

   if (iv%info(ssmt2)%nlocal > 0) then
      do n=1, iv%info(ssmt2)%nlocal
         if (iv%info(ssmt2)%proc_domain(1,n)) then
            do k=1, iv%info(ssmt2)%levels(n)
               call da_stats_calculate (n, k, iv%ssmt2(n)%rh(k)%qc, re%ssmt2(n)%rh(k), nrh, &
                  minimum, maximum, average, rms_err)

            end do
         end if    ! end if (iv%info(ssmt2)%proc_domain(1,n))
      end do
   end if

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nrh)
   iv%nstats(ssmt2) = nrh
      
   call da_proc_stats_combine(average, rms_err, minimum%value, maximum%value, &
      minimum%n, maximum%n, minimum%l, maximum%l)
   
   if (rootproc) then 
      if (nrh /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of O-A for SSMT2'
         call da_print_stats_ssmt2 (stats_unit, nrh, minimum, maximum, average, rms_err)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_ssmt2")

end subroutine da_ao_stats_ssmt2


subroutine da_oi_stats_ssmt1 (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in) :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in) :: iv            ! O-B

   type (maxmin_type) :: minimum
   type (maxmin_type) :: maximum
   integer            :: nt
   integer            :: n, k
   real               :: average, rms_err

   if (trace_use_dull) call da_trace_entry("da_oi_stats_ssmt1")

   nt = 0
   
   maximum = maxmin_type(-1.0E+20, 0, 0)
   minimum = maxmin_type(1.0E+20, 0, 0)
   average = 0.0
   rms_err = 0.0

   do n=1, iv%info(ssmt1)%nlocal
      if (iv%info(ssmt1)%proc_domain(1,n)) then
         do k=1, iv%info(ssmt1)%levels(n)   
            call da_stats_calculate(n, k, iv%ssmt1(n)%t(k)%qc, iv%ssmt1(n)%t(k)%inv, nt, &
               minimum, maximum, average, rms_err)
         end do
      end if    ! end if (iv%info(ssmt1)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nt)
   
   call da_proc_stats_combine(average, rms_err, minimum%value, maximum%value, &
      minimum%n, maximum%n, minimum%l, maximum%l)
   
   if (rootproc) then
      if (nt /= 0) then  
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of O-B for SSMT1'
            call da_print_stats_ssmt1(stats_unit, nt, minimum, maximum, average, rms_err)   
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_ssmt1")

end subroutine da_oi_stats_ssmt1


subroutine da_oi_stats_ssmt2 (stats_unit, iv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent (in)      :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (in)      :: iv            ! O-B

   type (maxmin_type) :: minimum
   type (maxmin_type) :: maximum
   integer            :: nrh
   integer            :: n, k
   real               :: average, rms_err

   if (trace_use_dull) call da_trace_entry("da_oi_stats_ssmt2")

   nrh = 0
   
   maximum = maxmin_type(-1.0E+20, 0, 0)
   minimum = maxmin_type(1.0E+20, 0, 0)
   average = 0.0
   rms_err = 0.0

   do n=1, iv%info(ssmt2)%nlocal
      if (iv%info(ssmt2)%proc_domain(1,n)) then
         do k=1, iv%info(ssmt2)%levels(n)   
            call da_stats_calculate(n, k, iv%ssmt2(n)%rh(k)%qc, iv%ssmt2(n)%rh(k)%inv, nrh, &
               minimum, maximum, average, rms_err)
         end do
      end if    ! end if (iv%info(ssmt2)%proc_domain(1,n))
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int(nrh)
   
   call da_proc_stats_combine(average, rms_err, minimum%value, maximum%value, &
      minimum%n, maximum%n, minimum%l, maximum%l)
   
   if (rootproc) then
      if (nrh /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of O-B for SSMT2'
         call da_print_stats_ssmt2(stats_unit, nrh, minimum, maximum, average, rms_err)
      end if 
   end if

   if (trace_use_dull) call da_trace_exit("da_oi_stats_ssmt2")

end subroutine da_oi_stats_ssmt2


subroutine da_print_stats_ssmt1(stats_unit, nt, minimum, maximum, average, rms_err)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,            intent(in)    :: stats_unit
   integer,            intent(inout) :: nt
   type (maxmin_type), intent(in)    :: minimum
   type (maxmin_type), intent(in)    :: maximum
   real,               intent(in)    :: average
   real,               intent(in)    :: rms_err

   if (trace_use_dull) call da_trace_entry("da_print_stats_ssmt1")
   
   write(unit=stats_unit, fmt='(a/)') '   var              T(K)     n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') '  Number: ', nt

   if (nt < 1) nt = 1
   
   write(unit=stats_unit, fmt='((a,f12.4,2i5))') &
      ' Minimum(n,k): ', minimum,    &
      ' Maximum(n,k): ', maximum
   write(unit=stats_unit, fmt='((a,f12.4,10x))') &
      ' Average     : ', average/real(nt),    &
      '    RMSE     : ', sqrt(rms_err/real(nt))

   if (trace_use_dull) call da_trace_exit("da_print_stats_ssmt1")

end subroutine da_print_stats_ssmt1


subroutine da_print_stats_ssmt2(stats_unit, nrh, minimum, maximum, &
                                 average, rms_err)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,            intent(in)    :: stats_unit
   integer,            intent(inout) :: nrh
   type (maxmin_type), intent(in)    :: minimum
   type (maxmin_type), intent(in)    :: maximum
   real,               intent(in)    :: average
   real,               intent(in)    :: rms_err

   if (trace_use_dull) call da_trace_entry("da_print_stats_ssmt2")
   
   write(unit=stats_unit, fmt='(a/)') '   var             rh(K)     n    k'

   write(unit=stats_unit, fmt='(a,i16,4i22)') &
        '  Number: ', nrh

   if (nrh < 1) nrh = 1
   
   write(unit=stats_unit, fmt='((a,f12.4,2i5))') &
        ' Minimum(n,k): ', minimum,    &
        ' Maximum(n,k): ', maximum
   write(unit=stats_unit, fmt='((a,f12.4,10x))') &
        ' Average     : ', average/real(nrh),    &
        '    RMSE     : ', sqrt(rms_err/real(nrh))

end subroutine da_print_stats_ssmt2


subroutine da_transform_xtoy_ssmt1 (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n  ! Loop counter.

   real, allocatable :: t (:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_ssmt1") 

   allocate (t(iv%info(ssmt1)%max_lev,iv%info(ssmt1)%n1:iv%info(ssmt1)%n2))

   call da_interp_lin_3d (grid%xa % t, iv%info(ssmt1), t)

   do n=iv%info(ssmt1)%n1,iv%info(ssmt1)%n2
      y%ssmt1(n)%t(:) = t(1:iv%info(ssmt1)%levels(n),n)
   end do

   deallocate (t)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_ssmt1")

end subroutine da_transform_xtoy_ssmt1


subroutine da_transform_xtoy_ssmt1_adj(iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n  ! Loop counter.

   real, allocatable :: t(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_ssmt1_adj") 

   allocate (t(iv%info(ssmt1)%max_lev,iv%info(ssmt1)%n1:iv%info(ssmt1)%n2))

   do n=iv%info(ssmt1)%n1,iv%info(ssmt1)%n2
      t(1:iv%info(ssmt1)%levels(n),n) = jo_grad_y%ssmt1(n)%t(:)
   end do

   call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(ssmt1), t)

   deallocate (t)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_ssmt1_adj") 

end subroutine da_transform_xtoy_ssmt1_adj


subroutine da_transform_xtoy_ssmt2 (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer  :: n  ! Loop counter.

   real, allocatable :: rh(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_ssmt2") 

   allocate (rh(1:iv%info(ssmt2)%max_lev,iv%info(ssmt2)%n1:iv%info(ssmt2)%n2))

   call da_interp_lin_3d (grid%xa%rh, iv%info(ssmt2), rh)

   do n=iv%info(ssmt2)%n1,iv%info(ssmt2)%n2
      y%ssmt2(n)%rh(:) = rh(1:iv%info(ssmt2)%levels(n),n)
   end do

   deallocate (rh)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_ssmt2") 

end subroutine da_transform_xtoy_ssmt2


subroutine da_transform_xtoy_ssmt2_adj(iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer  :: n  ! Loop counter.

   real, allocatable :: rh(:,:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_ssmt2_adj") 

   allocate (rh(1:iv%info(ssmt2)%max_lev,iv%info(ssmt2)%n1:iv%info(ssmt2)%n2))

   do n=iv%info(ssmt2)%n1,iv%info(ssmt2)%n2
      rh(1:iv%info(ssmt2)%levels(n),n) = jo_grad_y%ssmt2(n)%rh(:)
   end do

   call da_interp_lin_3d_adj (jo_grad_x%rh, iv%info(ssmt2), rh)

   deallocate (rh)

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_ssmt2_adj") 

end subroutine da_transform_xtoy_ssmt2_adj


subroutine da_calculate_grady_ssmi_tb(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Ob Inc. structure.
   type (y_type), intent(inout) :: re          ! Residual structure.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_ssmi_tb")

   do n=1, iv%info(ssmi_tb)%nlocal
      if (iv%ssmi_tb(n)%tb19v%qc < obs_qc_pointer) re%ssmi_tb(n)%tb19v = 0.0
      jo_grad_y%ssmi_tb(n)%tb19v = - re%ssmi_tb(n)%tb19v / (iv%ssmi_tb(n)%tb19v%error * iv%ssmi_tb(n)%tb19v%error)

      if (iv%ssmi_tb(n)%tb19h%qc < obs_qc_pointer) re%ssmi_tb(n)%tb19h = 0.0
      jo_grad_y%ssmi_tb(n)%tb19h = - re%ssmi_tb(n)%tb19h / (iv%ssmi_tb(n)%tb19h%error * iv%ssmi_tb(n)%tb19h%error)

      if (iv%ssmi_tb(n)%tb22v%qc < obs_qc_pointer) re%ssmi_tb(n)%tb22v = 0.0
      jo_grad_y%ssmi_tb(n)%tb22v = - re%ssmi_tb(n)%tb22v / (iv%ssmi_tb(n)%tb22v%error * iv%ssmi_tb(n)%tb22v%error)

      if (iv%ssmi_tb(n)%tb37v%qc < obs_qc_pointer) re%ssmi_tb(n)%tb37v = 0.0
      jo_grad_y%ssmi_tb(n)%tb37v = - re%ssmi_tb(n)%tb37v / (iv%ssmi_tb(n)%tb37v%error * iv%ssmi_tb(n)%tb37v%error)

      if (iv%ssmi_tb(n)%tb37h%qc < obs_qc_pointer) re%ssmi_tb(n)%tb37h = 0.0
      jo_grad_y%ssmi_tb(n)%tb37h = - re%ssmi_tb(n)%tb37h / (iv%ssmi_tb(n)%tb37h%error * iv%ssmi_tb(n)%tb37h%error)

      if (iv%ssmi_tb(n)%tb85v%qc < obs_qc_pointer) re%ssmi_tb(n)%tb85v = 0.0
      jo_grad_y%ssmi_tb(n)%tb85v = - re%ssmi_tb(n)%tb85v / (iv%ssmi_tb(n)%tb85v%error * iv%ssmi_tb(n)%tb85v%error)

      if (iv%ssmi_tb(n)%tb85h%qc < obs_qc_pointer) re%ssmi_tb(n)%tb85h = 0.0
      jo_grad_y%ssmi_tb(n)%tb85h = - re%ssmi_tb(n)%tb85h / (iv%ssmi_tb(n)%tb85h%error * iv%ssmi_tb(n)%tb85h%error)
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_ssmi_tb")

end subroutine da_calculate_grady_ssmi_tb


subroutine da_calculate_grady_ssmi_rv(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Ob Inc. structure.
   type (y_type), intent(inout) :: re          ! Residual structure.
   type (y_type), intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer                      :: n

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_ssmi_rv")

   do n=1, iv%info(ssmi_rv)%nlocal
      if (iv%ssmi_rv(n)%speed%qc < obs_qc_pointer) then
         re%ssmi_rv(n)%speed = 0.0
      end if
      jo_grad_y%ssmi_rv(n)%Speed = - re%ssmi_rv(n)%Speed / (iv%ssmi_rv(n)%Speed%error * iv%ssmi_rv(n)%Speed%error)

      if (iv%ssmi_rv(n)%tpw%qc < obs_qc_pointer) then
         re%ssmi_rv(n)%tpw = 0.0
      end if
      jo_grad_y%ssmi_rv(n)%tpw = -re%ssmi_rv(n)%tpw / (iv%ssmi_rv(n)%tpw%error * iv%ssmi_rv(n)%tpw%error)
   end do  

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_ssmi_rv")

end subroutine da_calculate_grady_ssmi_rv


subroutine da_calculate_grady_ssmt1(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_ssmt1")

   do n=1, iv%info(ssmt1)%nlocal
      do k=1, iv%info(ssmt1)%levels(n)
         if (iv%ssmt1(n)%t(k)%qc < obs_qc_pointer) re%ssmt1(n)%t(k) = 0.0
         jo_grad_y%ssmt1(n)%t(k) = -re%ssmt1(n)%t(k) / &
             (iv%ssmt1(n)%t(k)%error * iv%ssmt1(n)%t(k)%error)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_ssmt1")

end subroutine da_calculate_grady_ssmt1


subroutine da_calculate_grady_ssmt2(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Applies obs inverse on re-vector
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: re          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n, k

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_ssmt2")

   do n=1, iv%info(ssmt2)%nlocal
      do k=1, iv%info(ssmt2)%levels(n)
         if (iv%ssmt2(n)%rh(k)%qc < obs_qc_pointer) re%ssmt2(n)%rh(k) = 0.0
         jo_grad_y%ssmt2(n)%rh(k) = -re%ssmt2(n)%rh(k) / &
            (iv%ssmt2(n)%rh(k)%error * iv%ssmt2(n)%rh(k)%error)
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_ssmt2")

end subroutine da_calculate_grady_ssmt2



subroutine da_tb_adj(ifreq,theta,p0,ta,gamma,sst,wv,hwv,u,alw,zcld,            &
!                  tbv,tbh,                                                  &
                  ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,                   &
                  ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            )

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! Output : ADJ_p0,  ADJ_ta,   ADJ_gamma, ADJ_sst, ADJ_wv, ADJ_hwv, ADJ_u
   !          ADJ_alw, ADJ_zcld
   ! Input  : ADJ_tbv, ADJ_tbh,  tbv,  tbh
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in   ) :: ifreq
   real   , intent(in   ) :: theta,p0,ta,gamma,sst,wv,hwv,u,alw,zcld
   real   , intent(inout) :: ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,     &
                             ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld
   real   , intent(in   ) :: ADJ_tbv,ADJ_tbh 
!   real   , intent(in   ) :: tbv,tbh

   real :: freq(4),ebiasv(4),ebiash(4),cf1(4),cf2(4),cg(4)

   real :: f,costhet,gx2,tbup,tbdn,tauatm,sigma,remv,remh
   real :: effangv,effangh,tbdnv,foam,foam_save,emissv,emissh 
   real :: dum
   real :: refv,refh,semv,semh,scatv,scath,tbdnh
   real :: ADJ_gx2,ADJ_tbup,ADJ_tbdn,ADJ_tauatm,ADJ_sigma
   real :: ADJ_tremv,ADJ_remh,ADJ_effangv,ADJ_effangh
   real :: ADJ_tbdnh,ADJ_dum,ADJ_foam,ADJ_emissv
   real :: ADJ_emissh,ADJ_refv,ADJ_refh,ADJ_semv
   real :: ADJ_semh,ADJ_scatv,ADJ_scath,ADJ_remv,ADJ_tbdnv
   real :: ADJ_theta

   real :: fem

   data freq/19.35,22.235,37.0,85.5/

   ! empirical bias corrections for surface emissivity

   data ebiasv/0.0095, 0.005,-0.014, -0.0010/
   data ebiash/0.004,   0.0,-0.023, 0.043/


   data cf1/0.0015,0.004,0.0027,0.005/
   data cg/4.50e-3, 5.200e-3, 5.5e-3, 7.2e-3 /

   data cf2/0.0023,0.000,0.0002,-0.006/

   ! 'foam' emissivity
   data fem /1.0/

   if (trace_use) call da_trace_entry("da_tb_adj")

   f=0.0;costhet=0.0;gx2=0.0;tbup=0.0;tbdn=0.0;tauatm=0.0
   sigma=0.0;remv=0.0;remh=0.0;effangv=0.0;effangh=0.0
   tbdnv=0.0;foam=0.0;foam_save=0.0;emissv=0.0;emissh=0.0
   dum=0.0;refv=0.0;refh=0.0;semv=0.0;semh=0.0;scatv=0.0
   scath=0.0;tbdnh=0.0;ADJ_gx2=0.0;ADJ_tbup=0.0;ADJ_tbdn=0.0
   ADJ_tauatm=0.0;ADJ_sigma=0.0;ADJ_tremv=0.0;ADJ_remh=0.0
   ADJ_effangv=0.0;ADJ_effangh=0.0;ADJ_tbdnh=0.0;ADJ_dum=0.0
   ADJ_foam=0.0;ADJ_emissv=0.0;ADJ_emissh=0.0;ADJ_refv=0.0
   ADJ_refh=0.0;ADJ_semv=0.0;ADJ_semh=0.0;ADJ_scatv=0.0
   ADJ_scath=0.0;ADJ_remv=0.0;ADJ_tbdnv=0.0
   ADJ_theta =0.0


   ! write (unit=stdout,fmt=*) 'ifreq',ifreq,theta,p0,ta,gamma,sst,wv,hwv,u,alw,zcld,          &
   !              tbv,tbh,                                                  &
   !              ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,                   &
   !               ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            

   f = freq(ifreq)
   costhet = cos(theta*0.017453)

   ! effective surface slope variance

   gx2 = cg(ifreq)*    u

   ! get upwelling atmospheric brightness temperature

   call tbatmos(ifreq,theta,p0,wv,hwv,ta,gamma,alw,zcld, tbup,tbdn,tauatm)

   ! convert transmittance to optical depth

   sigma = -alog(tauatm)*costhet

   ! get rough surface emissivity

   call roughem(ifreq,gx2,sst,theta,remv,remh)

   ! get effective zenith angles for scattered radiation at surface

   call effang(ifreq,theta,gx2,sigma,effangv,effangh)

   ! get effective sky brightness temperatures for scattered radiation

   call tbatmos(ifreq,effangv,p0,wv,hwv,ta,gamma,alw,zcld, dum,tbdnv,dum)

   call tbatmos(ifreq,effangh,p0,wv,hwv,ta,gamma,alw,zcld, dum,tbdnh,dum)

   ! compute 'foam' coverage

   foam = cf1(ifreq)*    u

   if (u .gt. 5.0) then
      foam_save = foam
      foam =     foam + cf2(ifreq)*(   u-5.0)
   end if

   ! compute surface emissivities and reflectivity

   emissv =     foam*fem + (1.0 - foam)*(remv + ebiasv(ifreq))
   emissh =     foam*fem + (1.0 - foam)*(remh + ebiash(ifreq))
   refv =   1.0 - emissv
   refh =   1.0 - emissh

   ! compute surface emission term

   semv = sst*emissv
   semh = sst*emissh

   ! compute surface scattering term

   scatv = refv*tbdnv
   scath = refh*tbdnh

   ! combine to get space-observed brightness temperature

   ! tbv =     tbup + tauatm*(semv + scatv)
   ! tbh =     tbup + tauatm*(semh + scath)


   ! start
   ! write (unit=stdout,fmt=*) 'ifreq 1',ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,  &
   !                 ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            


   ADJ_tbup   = ADJ_tbh                    !!! first
   ADJ_tauatm = ADJ_tbh*(semh + scath)     !!! first
   ADJ_semh   = tauatm*ADJ_tbh             !!! first
   ADJ_scath  = tauatm*ADJ_tbh             !!! first

   ADJ_tbup   = ADJ_tbv                  + ADJ_tbup
   ADJ_tauatm = ADJ_tbv*(semv + scatv)   + ADJ_tauatm
   ADJ_semv   = tauatm*ADJ_tbv             !!! first
   ADJ_scatv  = tauatm*ADJ_tbv             !!! first

   ADJ_refh   = ADJ_scath*tbdnh            !!! first
   ADJ_tbdnh  = refh*ADJ_scath             !!! first
   ADJ_refv   = ADJ_scatv*tbdnv            !!! first
   ADJ_tbdnv  = refv*ADJ_scatv             !!! first
   ADJ_sst    = ADJ_semh*emissh          + ADJ_sst
   ADJ_emissh = sst*ADJ_semh               !!! first
   ADJ_sst    = ADJ_semv*emissv          + ADJ_sst
   ADJ_emissv = sst*ADJ_semv               !!! first

   ADJ_emissh = - ADJ_refh               + ADJ_emissh
   ADJ_emissv = - ADJ_refv               + ADJ_emissv

   ADJ_foam   =   ADJ_emissh*fem                      !!! first
   ADJ_foam   = - ADJ_emissh*(remh + ebiash(ifreq)) + ADJ_foam
   ADJ_remh   =   (1.0 - foam)*ADJ_emissh             !!! first
   ADJ_foam   =   ADJ_emissv*fem                    + ADJ_foam
   ADJ_foam   = - ADJ_emissv*(remv + ebiasv(ifreq)) + ADJ_foam
   ADJ_remv   =   (1.0 - foam)*ADJ_emissv             !!! first

   if (u .gt. 5.0) then
     ADJ_u = cf2(ifreq)*ADJ_foam  + ADJ_u
     foam=foam_save
   end if
   ADJ_u = cf1(ifreq)*ADJ_foam    + ADJ_u
   
   ADJ_dum = 0.0
   dum     = 0.0
   ADJ_effangh = 0.0
   call da_tbatmos_adj(ifreq,effangh,p0,wv,hwv,ta,gamma,alw,    &
                       zcld,dum,tbdnh,dum,                      &
                       ADJ_effangh,ADJ_p0,ADJ_wv,ADJ_hwv,       &
                       ADJ_ta,ADJ_gamma,ADJ_alw,ADJ_zcld,       &
                       ADJ_dum,ADJ_tbdnh,ADJ_dum                )
   dum     = 0.0
   ADJ_dum = 0.0

   ! write (unit=stdout,fmt=*) 'ifreq 2',ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,  &
   !                 ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            

   ADJ_effangv = 0.0
   call da_tbatmos_adj(ifreq,effangv,p0,wv,hwv,ta,gamma,alw,    &
                       zcld,dum,tbdnv,dum,                      &
                       ADJ_effangv,ADJ_p0,ADJ_wv,ADJ_hwv,       &
                       ADJ_ta,ADJ_gamma,ADJ_alw,ADJ_zcld,       & 
                       ADJ_dum,ADJ_tbdnv,ADJ_dum                )

   ADJ_gx2=0.0
   ADJ_sigma=0.0
   ! write (unit=stdout,fmt=*) 'ifreq 3',ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,  &
   !                 ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            

   call da_effang_adj(ifreq,theta,gx2,sigma,effangv,effangh,    &
                      ADJ_gx2,ADJ_sigma,ADJ_effangv,ADJ_effangh )

   call da_roughem_adj(ifreq,gx2,sst,theta,remv,remh,         &
                       ADJ_gx2,ADJ_sst,ADJ_remv,ADJ_remh      )

   ADJ_tauatm = - costhet*ADJ_sigma/tauatm + ADJ_tauatm


   ! write (unit=stdout,fmt=*) 'ifreq 4',ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,  &
   !              ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            

   call da_tbatmos_adj(ifreq,theta,p0,wv,hwv,ta,gamma,alw,zcld, &
                       tbup,tbdn,tauatm,                        &
                       ADJ_theta,ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta,ADJ_gamma,  &
                       ADJ_alw,ADJ_zcld,ADJ_tbup,ADJ_tbdn,      &
                       ADJ_tauatm                               )

   ADJ_theta=0.0   ! first

   ADJ_u = cg(ifreq)*ADJ_gx2 + ADJ_u

   ! write (unit=stdout,fmt=*) 'ifreq 5',ADJ_p0,ADJ_ta,ADJ_gamma,ADJ_sst,ADJ_wv,  &
   !              ADJ_hwv,ADJ_u,ADJ_alw,ADJ_zcld,ADJ_tbv,ADJ_tbh            
   ! end

   if (trace_use) call da_trace_exit("da_tb_adj")

end subroutine da_tb_adj


subroutine da_sigma_v_adj(ifreq,p0,wv,hwv,ta,gamma,sigv,                   &
                           ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta,ADJ_gamma,ADJ_sigma_v)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! output: ADJ_p0, ADJ_wv, ADJ_hwv, ADJ_ta, ADJ_gamma
   ! input: ADJ_sigma_v
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: ifreq
   real,    intent(in)    :: p0,wv,hwv,ta,gamma  ! base field
   real,    intent(inout) :: ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta,ADJ_gamma
   real,    intent(in)    :: ADJ_sigma_v
   real,    intent(out)   :: sigv

   real :: wvc, wvcor(4)
   real :: ADJ_wvc

   real :: voh1,otbar1,pbar1
   real :: term21,term31,term41,term51,term61
   real :: a11,a21,a31,a41,a51,a61
   real :: ADJ_voh1,ADJ_otbar1,ADJ_pbar1
   real :: ADJ_term21,ADJ_term31,ADJ_term41,ADJ_term51,ADJ_term61

   real :: voh2,otbar2,pbar2
   real :: term22,term32,term42,term52,term62
   real :: a12,a22,a32,a42,a52,a62
   real :: ADJ_voh2,ADJ_otbar2,ADJ_pbar2
   real :: ADJ_term22,ADJ_term32,ADJ_term42,ADJ_term52,ADJ_term62

   real :: voh3,otbar3,pbar3
   real :: term23,term33,term43,term53,term63
   real :: a13,a23,a33,a43,a53,a63
   real :: ADJ_voh3,ADJ_otbar3,ADJ_pbar3
   real :: ADJ_term23,ADJ_term33,ADJ_term43,ADJ_term53,ADJ_term63

   real :: voh4,otbar4,pbar4
   real :: term24,term34,term44,term54,term64
   real :: a14,a24,a34,a44,a54,a64
   real :: ADJ_voh4,ADJ_otbar4,ADJ_pbar4
   real :: ADJ_term24,ADJ_term34,ADJ_term44,ADJ_term54,ADJ_term64

   real :: const1,const2,const3,const4
   real :: h1,h2,h3,h4

   real :: ADJ_sigv

   data const1,const2,const3,const4/0.6,2.8,0.2,0.2/
   data h1,h2,h3,h4/5.0,4.9,6.8,6.4/

   data a11,a21,a31,a41,a51,a61/-.13747e-2,-.43061e-4, .14618e+1,  &
     .25101e-3, .14635e-1,-.18588e+3/
   data a12,a22,a32,a42,a52,a62/ .22176e-1,-.32367e-4,-.10840e-4,  &
     -.63578e-1, .16988e-7,-.29824e+2/
   data a13,a23,a33,a43,a53,a63/-.10566e-2,-.12906e-3, .56975e+0,  &
      .10828e-8,-.17551e-7, .48601e-1/
   data a14,a24,a34,a44,a54,a64/-.60808e-2,-.70936e-3, .28721e+1,  &
      .42636e-8,-.82910e-7, .26166e+0/

   ! data wvcor/1.01,0.95,1.06,0.92/
   data wvcor/1.02,0.98,1.02,0.88/

   if (trace_use) call da_trace_entry("da_sigma_v_adj")

   ! use modified water vapor value to correct for errors in theoretical absorption

   wvc        = 0.0
   ADJ_wvc    = 0.0
   voh1       = 0.0
   otbar1     = 0.0
   pbar1      = 0.0
   term21     = 0.0
   term31     = 0.0
   term41     = 0.0
   term51     = 0.0
   term61     = 0.0
   ADJ_voh1   = 0.0
   ADJ_otbar1 = 0.0
   ADJ_pbar1  = 0.0
   ADJ_term21 = 0.0
   ADJ_term31 = 0.0
   ADJ_term41 = 0.0
   ADJ_term51 = 0.0
   ADJ_term61 = 0.0

   voh2       = 0.0
   otbar2     = 0.0
   pbar2      = 0.0
   term22     = 0.0
   term32     = 0.0
   term42     = 0.0
   term52     = 0.0
   term62     = 0.0
   ADJ_voh2   = 0.0
   ADJ_otbar2 = 0.0
   ADJ_pbar2  = 0.0
   ADJ_term22 = 0.0
   ADJ_term32 = 0.0
   ADJ_term42 = 0.0
   ADJ_term52 = 0.0
   ADJ_term62 = 0.0

   voh3       = 0.0
   otbar3     = 0.0
   pbar3      = 0.0
   term23     = 0.0
   term33     = 0.0
   term43     = 0.0
   term53     = 0.0
   term63     = 0.0
   ADJ_voh3   = 0.0
   ADJ_otbar3 = 0.0
   ADJ_pbar3  = 0.0
   ADJ_term23 = 0.0
   ADJ_term33 = 0.0
   ADJ_term43 = 0.0
   ADJ_term53 = 0.0
   ADJ_term63 = 0.0

   voh4       = 0.0
   otbar4     = 0.0
   pbar4      = 0.0
   term24     = 0.0
   term34     = 0.0
   term44     = 0.0
   term54     = 0.0
   term64     = 0.0
   ADJ_voh4   = 0.0
   ADJ_otbar4 = 0.0
   ADJ_pbar4  = 0.0
   ADJ_term24 = 0.0
   ADJ_term34 = 0.0
   ADJ_term44 = 0.0
   ADJ_term54 = 0.0
   ADJ_term64 = 0.0

   sigv     = 0.0
   ADJ_sigv = 0.0

   wvc = wv*wvcor(ifreq)

   if (ifreq == 1) then
      pbar1  = p0/(1.0 + hwv/h1)
      voh1   = wv/hwv
      term21 = a21*voh1
      otbar1 =  1.0/(ta - const1*gamma*hwv)
      term31 = a31*otbar1
      term61 = a61*otbar1*otbar1
      term41 = a41*pbar1*otbar1
      term51 = a51*voh1*otbar1
      sigv   = a11 + term21 + term31 + term41 + term51 + term61
   else if (ifreq == 2) then
      pbar2  = p0/(1.0 + hwv/h2)
      term22 = a22*pbar2
      term52 = a52*pbar2*pbar2
      voh2   = wv/hwv
      term32 = a32*voh2
      otbar2 = 1.0/(ta - const2*gamma*hwv)
      term42 = a42*otbar2
      term62 = a62*otbar2*otbar2
      sigv   = a12 + term22 + term32 + term42 + term52 + term62
   else if (ifreq==3) then
      pbar3  = p0/(1.0 + hwv/h3)
      term43 = a43*pbar3*pbar3
      voh3   = wv/hwv
      term23 = a23*voh3
      otbar3 = 1.0/(ta - const3*gamma*hwv)
      term33 = a33*otbar3
      term53 = a53*pbar3*voh3
      term63 = a63*otbar3*voh3
      sigv   = a13 + term23 + term33 + term43 + term53 + term63
   else if (ifreq == 4) then
      pbar4  = p0/(1.0 + hwv/h4)
      term44 = a44*pbar4*pbar4
      voh4   = wv/hwv
      term24 = a24*voh4
      otbar4 = 1.0/(ta - const4*gamma*hwv)
      term34 = a34*otbar4
      term54 = a54*pbar4*voh4
      term64 = a64*otbar4*voh4
      sigv   = a14 + term24 + term34 + term44 + term54 + term64
   else
      sigv = 0.0
   end if

   ADJ_sigv = ADJ_sigma_v*wvc
   ADJ_wvc  = sigv*ADJ_sigma_v

   if (ifreq == 1) then
      ADJ_term21 = ADJ_sigv 
      ADJ_term31 = ADJ_sigv
      ADJ_term41 = ADJ_sigv
      ADJ_term51 = ADJ_sigv
      ADJ_term61 = ADJ_sigv

      ADJ_voh1   = a51*ADJ_term51*otbar1
      ADJ_otbar1 = a51*voh1*ADJ_term51

      ADJ_pbar1  = a41*ADJ_term41*otbar1
      ADJ_otbar1 = a41*pbar1*ADJ_term41 + ADJ_otbar1
      ADJ_otbar1 = 2.0*a61*otbar1*ADJ_term61 + ADJ_otbar1

      ADJ_otbar1 = a31*ADJ_term31 + ADJ_otbar1

      ADJ_ta    = -otbar1*otbar1*ADJ_otbar1  + ADJ_ta
      ADJ_hwv   = otbar1*otbar1*const1*gamma*ADJ_otbar1  + ADJ_hwv
      ADJ_gamma = otbar1*otbar1*const1*ADJ_otbar1*hwv  + ADJ_gamma      

      ADJ_voh1  = a21*ADJ_term21 + ADJ_voh1

      ADJ_wv    = ADJ_voh1/hwv  + ADJ_wv
      ADJ_hwv   = -voh1*ADJ_voh1/hwv + ADJ_hwv

      ADJ_p0    = ADJ_pbar1/(1.0 + hwv/h1)  + ADJ_p0
      ADJ_hwv   = -pbar1*ADJ_pbar1/(h1*(1.0 + hwv/h1)) + ADJ_hwv
   else if (ifreq==2) then
      ADJ_term22 = ADJ_sigv 
      ADJ_term32 = ADJ_sigv
      ADJ_term42 = ADJ_sigv
      ADJ_term52 = ADJ_sigv
      ADJ_term62 = ADJ_sigv

      ADJ_otbar2 = 2.0*a62*otbar2*ADJ_term62

      ADJ_otbar2 = a42*ADJ_term42 + ADJ_otbar2

      ADJ_ta     = -otbar2*otbar2*ADJ_otbar2  + ADJ_ta
      ADJ_hwv    =  otbar2*otbar2*const2*gamma*ADJ_otbar2 + ADJ_hwv
      ADJ_gamma  =  otbar2*otbar2*const2*ADJ_otbar2*hwv + ADJ_gamma

      ADJ_voh2   = a32*ADJ_term32

      ADJ_wv     = ADJ_voh2/hwv + ADJ_wv
      ADJ_hwv    = -voh2*ADJ_voh2/hwv + ADJ_hwv

      ADJ_pbar2  = 2.0*a52*pbar2*ADJ_term52

      ADJ_pbar2  = a22*ADJ_term22 + ADJ_pbar2

      ADJ_p0     = ADJ_pbar2/(1.0 + hwv/h2) + ADJ_p0
      ADJ_hwv    = -pbar2*ADJ_pbar2/h2/(1.0 + hwv/h2) + ADJ_hwv
   else if (ifreq==3) then
      ADJ_term23 = ADJ_sigv
      ADJ_term33 = ADJ_sigv
      ADJ_term43 = ADJ_sigv
      ADJ_term53 = ADJ_sigv
      ADJ_term63 = ADJ_sigv

      ADJ_otbar3 = a63*ADJ_term63*voh3
      ADJ_voh3   = a63*otbar3*ADJ_term63

      ADJ_pbar3  = a53*ADJ_term53*voh3
      ADJ_voh3   = a53*pbar3*ADJ_term53 + ADJ_voh3

      ADJ_otbar3 = a33*ADJ_term33 + ADJ_otbar3

      ADJ_ta     = -otbar3*otbar3*ADJ_otbar3 + ADJ_ta
      ADJ_hwv    =  otbar3*otbar3*const3*gamma*ADJ_otbar3 + ADJ_hwv
      ADJ_gamma  =  otbar3*otbar3*const3*ADJ_otbar3*hwv + ADJ_gamma

      ADJ_voh3   = a23*ADJ_term23 + ADJ_voh3

      ADJ_wv     = ADJ_voh3/hwv  + ADJ_wv
      ADJ_hwv    = -voh3*ADJ_voh3/hwv + ADJ_hwv

      ADJ_pbar3 = 2.0*a43*pbar3*ADJ_term43 + ADJ_pbar3

      ADJ_p0    = ADJ_pbar3/(1.0 + hwv/h3) + ADJ_p0
      ADJ_hwv   = -pbar3*ADJ_pbar3/h3/(1.0 + hwv/h3) + ADJ_hwv
   else if (ifreq == 4) then
      ADJ_term24 = ADJ_sigv
      ADJ_term34 = ADJ_sigv
      ADJ_term44 = ADJ_sigv
      ADJ_term54 = ADJ_sigv
      ADJ_term64 = ADJ_sigv

      ADJ_otbar4 = a64*ADJ_term64*voh4
      ADJ_voh4   = a64*otbar4*ADJ_term64 

      ADJ_pbar4  = a54*ADJ_term54*voh4
      ADJ_voh4   = a54*pbar4*ADJ_term54 + ADJ_voh4

      ADJ_otbar4 = a34*ADJ_term34 + ADJ_otbar4

      ADJ_ta     = -otbar4*otbar4*ADJ_otbar4  + ADJ_ta
      ADJ_hwv    =  otbar4*otbar4*const4*gamma*ADJ_otbar4 + ADJ_hwv
      ADJ_gamma  =  otbar4*otbar4*const4*ADJ_otbar4*hwv + ADJ_gamma

      ADJ_voh4   = a24*ADJ_term24 + ADJ_voh4

      ADJ_wv     = ADJ_voh4/hwv + ADJ_wv
      ADJ_hwv    = -voh4*ADJ_voh4/hwv + ADJ_hwv

      ADJ_pbar4  = 2.0*a44*pbar4*ADJ_term44 + ADJ_pbar4

      ADJ_p0     = ADJ_pbar4/(1.0 + hwv/h4) + ADJ_p0
      ADJ_hwv    = -pbar4*ADJ_pbar4/h4/(1.0 + hwv/h4) + ADJ_hwv
   end if

   ADJ_wv  = ADJ_wvc*wvcor(ifreq) + ADJ_wv

   if (trace_use) call da_trace_exit("da_sigma_v_adj")

end subroutine da_sigma_v_adj


subroutine da_effang_adj(ifreq,theta,gx2,sigma,effangv,effangh,     &
   ADJ_gx2,ADJ_sigma,ADJ_effangv,ADJ_effangh  )

   implicit none

   !-----------------------------------------------------------------
   ! Purpose: Calculate the effective zenith angle of reflected microwave
   ! radiation at SSM/I frequencies for vertical and horizontal polarization
   !
   ! Output: ADJ_gx2, ADJ_sigma
   ! Input: ADJ_effangv,ADJ_effangh,effangv,effangh
   !-----------------------------------------------------------------
 
   integer, intent(in)    :: ifreq
   real,    intent(in)    :: theta,gx2,sigma
   real,    intent(inout) :: ADJ_gx2,ADJ_sigma
   real,    intent(inout) :: ADJ_effangv,ADJ_effangh
   real,    intent(out)   :: effangv,effangh

   real c19v,c19h,c22v,c22h,c37v,c37h,c85v,c85h
   real s19v(6),s19h(6),s22v(6),s22h(6), &
        s37v(6),s37h(6),s85v(6),s85h(6)

   real :: alnsig,gg,ggg,xd,xx,xd_save,xx_save
   real :: z1,z2,z3,z4,z5,z6
   real :: y,dth,angh,angv,alnthv,alnthh,alnthv_save,alnthh_save
   real :: ADJ_alnsig,ADJ_gg,ADJ_ggg,ADJ_xd
   real :: ADJ_z1,ADJ_z2,ADJ_z3,ADJ_z4,ADJ_z5,ADJ_z6,ADJ_alnthv
   real :: ADJ_y,ADJ_dth,ADJ_angh,ADJ_angv,ADJ_xx,ADJ_alnthh

   data c19v,c19h,c22v,c22h,c37v,c37h,c85v,c85h &
     /-.5108,.5306,-.5108,.5306,-.6931,.1823,-.9163,.3000/
   data s19v /.225E+2,.698E+2,-.238E+2,-.648E+1,.402E+0,.262E+1/
   data s19h /.743E+1,.507E+2,-.206E+2,-.191E+1,.648E-1,.291E+1/
   data s22v /.249E+2,.701E+2,-.240E+2,-.714E+1,.405E+0,.256E+1/
   data s22h /.498E+1,.442E+2,-.190E+2,-.129E+1,.803E-2,.277E+1/
   data s37v /.215E+2,.573E+2,-.211E+2,-.670E+1,.443E+0,.253E+1/
   data s37h /.869E+1,.571E+2,-.257E+2,-.302E+1,.237E+0,.386E+1/
   data s85v /.116E+2,.263E+2,-.101E+2,-.358E+1,.270E+0,.175E+1/
   data s85h /.736E+1,.568E+2,-.254E+2,-.248E+1,.196E+0,.387E+1/

   if (trace_use) call da_trace_entry("da_effang_adj")

   alnsig=0.0;gg=0.0;ggg=0.0;xd=0.0;xx=0.0;xd_save=0.0;xx_save=0.0
   z1=0.0;z2=0.0;z3=0.0;z4=0.0;z5=0.0;z6=0.0;y=0.0;dth=0.0;angh=0.0
   angv=0.0;alnthv=0.0;alnthh=0.0;alnthv_save=0.0;alnthh_save=0.0
   ADJ_alnsig=0.0;ADJ_gg=0.0;ADJ_ggg=0.0;ADJ_xd=0.0
   ADJ_z1=0.0;ADJ_z2=0.0;ADJ_z3=0.0;ADJ_z4=0.0;ADJ_z5=0.0
   ADJ_z6=0.0;ADJ_alnthv=0.0;ADJ_y=0.0;ADJ_dth=0.0;ADJ_angh=0.0
   ADJ_angv=0.0;ADJ_xx=0.0;ADJ_alnthh=0.0


   if (gx2 .le. 0.0 .or. sigma .le. 0.0) then
      effangv = theta
      effangh = theta
      return
   end if
   alnsig = alog(sigma)
   gg  = gx2*gx2
   ggg = gg*gx2

   if (ifreq .eq. 1) then 
      xd =      alnsig - c19v
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthv =  s19v(1)*z1 + s19v(2)*z2 + s19v(3)*z3 + &
                s19v(4)*z4 + s19v(5)*z5 + s19v(6)*z6
      alnthv_save = alnthv
      alnthv =  alnthv + 3.611

      xd_save = xd
      xx_save = xx

      xd =  alnsig - c19h
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2

      alnthh =  s19h(1)*z1 + s19h(2)*z2 + s19h(3)*z3 + &
                s19h(4)*z4 + s19h(5)*z5 + s19h(6)*z6
      alnthh_save = alnthh
      alnthh =  alnthh + 3.611

   else if (ifreq .eq. 2) then 
      xd =      alnsig - c22v
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthv =  s22v(1)*z1 + s22v(2)*z2 + s22v(3)*z3 + &
                s22v(4)*z4 + s22v(5)*z5 + s22v(6)*z6

      alnthv_save = alnthv
      alnthv =  alnthv + 3.611
 
      xd_save = xd
      xx_save = xx

      xd =      alnsig - c22h
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthh =  s22h(1)*z1 + s22h(2)*z2 + s22h(3)*z3 + &
                   s22h(4)*z4 + s22h(5)*z5 + s22h(6)*z6

      alnthh_save = alnthh
      alnthh =  alnthh + 3.611

   else if (ifreq .eq. 3) then 
      xd =      alnsig - c37v
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthv =  s37v(1)*z1 + s37v(2)*z2 + s37v(3)*z3 + &
                s37v(4)*z4 + s37v(5)*z5 + s37v(6)*z6

      alnthv_save = alnthv
      alnthv =  alnthv + 3.611
 
      xd_save = xd
      xx_save = xx

      xd =      alnsig - c37h
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthh =  s37h(1)*z1 + s37h(2)*z2 + s37h(3)*z3 + &
                s37h(4)*z4 + s37h(5)*z5 + s37h(6)*z6
 
      alnthh_save = alnthh
      alnthh =  alnthh + 3.611

   else if (ifreq .eq. 4) then 
      xd =      alnsig - c85v
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthv =  s85v(1)*z1 + s85v(2)*z2 + s85v(3)*z3 + &
                s85v(4)*z4 + s85v(5)*z5 + s85v(6)*z6

      alnthv_save = alnthv
      alnthv =  alnthv + 3.611

      xd_save = xd
      xx_save = xx

      xd =      alnsig - c85h
      xx =  xd*xd
      z1 =  xx*ggg
      z2 =  xd*ggg
      z3 =  xd*gg
      z4 =  xx*gg
      z5 =  xx*gx2
      z6 =  xd*gx2
      alnthh =  s85h(1)*z1 + s85h(2)*z2 + s85h(3)*z3 + &
                s85h(4)*z4 + s85h(5)*z5 + s85h(6)*z6

      alnthh_save = alnthh
      alnthh =  alnthh + 3.611
   end if

   angv =   90.0 - exp(alnthv)
   angh =   90.0 - exp(alnthh)
       y    =   1.0 - 28.0*gx2
   if (y .lt. 0.0) then
       y = 0.0
   end if

   dth     = (theta - 53.0)*y
   effangv = angv + dth
   effangh = angh + dth

   ! start

   if (gx2 .le. 0.0 .or. sigma .le. 0.0) then
      ADJ_effangv = 0.0
      ADJ_effangh = 0.0
      return
   end if

   ADJ_angh  = ADJ_effangh
   ADJ_dth   = ADJ_effangh

   ADJ_angv  = ADJ_effangv
   ADJ_dth   = ADJ_effangv  + ADJ_dth

   ADJ_y     = (theta - 53.0)*ADJ_dth  

   if (y .lt. 0.0) then
      ADJ_y = 0.0
   end if

   ADJ_gx2  = - 28.0*ADJ_y + ADJ_gx2

   ADJ_alnthh = - ADJ_angh*exp(alnthh)
   ADJ_alnthv = - ADJ_angv*exp(alnthv)

   if (ifreq .eq. 1) then 

      alnthh = alnthh_save

      ADJ_z1  = s19h(1)*ADJ_alnthh
      ADJ_z2  = s19h(2)*ADJ_alnthh
      ADJ_z3  = s19h(3)*ADJ_alnthh
      ADJ_z4  = s19h(4)*ADJ_alnthh
      ADJ_z5  = s19h(5)*ADJ_alnthh
      ADJ_z6  = s19h(6)*ADJ_alnthh

      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg

      ADJ_xd  = 2.0*xd*ADJ_xx + ADJ_xd
      ADJ_alnsig = ADJ_xd

      ! vertical
      xd =  xd_save
      xx =  xx_save

      alnthv = alnthv_save

      ADJ_z1   = s19v(1)*ADJ_alnthv  
      ADJ_z2   = s19v(2)*ADJ_alnthv 
      ADJ_z3   = s19v(3)*ADJ_alnthv 
      ADJ_z4   = s19v(4)*ADJ_alnthv
      ADJ_z5   = s19v(5)*ADJ_alnthv
      ADJ_z6   = s19v(6)*ADJ_alnthv

      ADJ_xd  = ADJ_z6*gx2 !+ ADJ_xd
      ADJ_gx2 = xd*ADJ_z6   + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2 !+ ADJ_xx
      ADJ_gx2 = xx*ADJ_z5   + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg   + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg

      ADJ_xd  = 2.0*xd*ADJ_xx + ADJ_xd

      ADJ_alnsig = ADJ_xd + ADJ_alnsig

   else if (ifreq .eq. 2) then 

      ADJ_z1 = s22h(1)*ADJ_alnthh
      ADJ_z2 = s22h(2)*ADJ_alnthh
      ADJ_z3 = s22h(3)*ADJ_alnthh
      ADJ_z4 = s22h(4)*ADJ_alnthh
      ADJ_z5 = s22h(5)*ADJ_alnthh
      ADJ_z6 = s22h(6)*ADJ_alnthh
 
      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx 
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg

      ADJ_xd  =  2.0*xd*ADJ_xx + ADJ_xd
      ADJ_alnsig = ADJ_xd

      ! vertical
      xd =  xd_save
      xx =  xx_save

      alnthv = alnthv_save

      ADJ_z1  = s22v(1)*ADJ_alnthv
      ADJ_z2  = s22v(2)*ADJ_alnthv
      ADJ_z3  = s22v(3)*ADJ_alnthv
      ADJ_z4  = s22v(4)*ADJ_alnthv
      ADJ_z5  = s22v(5)*ADJ_alnthv
      ADJ_z6  = s22v(6)*ADJ_alnthv

      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2 
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg
      ADJ_xd  =  2.0*xd*ADJ_xx + ADJ_xd
      ADJ_alnsig = ADJ_xd     + ADJ_alnsig

   else if (ifreq .eq. 3) then 

      ADJ_z1  = s37h(1)*ADJ_alnthh
      ADJ_z2  = s37h(2)*ADJ_alnthh
      ADJ_z3  = s37h(3)*ADJ_alnthh
      ADJ_z4  = s37h(4)*ADJ_alnthh
      ADJ_z5  = s37h(5)*ADJ_alnthh
      ADJ_z6  = s37h(6)*ADJ_alnthh

      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2 
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg
      ADJ_xd  = 2.0*xd*ADJ_xx + ADJ_xd
      ADJ_alnsig = ADJ_xd

      ! vertical
      xd =  xd_save
      xx =  xx_save

      alnthv = alnthv_save

      ADJ_z1  = s37v(1)*ADJ_alnthv
      ADJ_z2  = s37v(2)*ADJ_alnthv
      ADJ_z3  = s37v(3)*ADJ_alnthv
      ADJ_z4  = s37v(4)*ADJ_alnthv
      ADJ_z5  = s37v(5)*ADJ_alnthv
      ADJ_z6  = s37v(6)*ADJ_alnthv
  
      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  =ADJ_z3*gg   + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd 
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg
      ADJ_xd  = 2.0*xd*ADJ_xx + ADJ_xd
      ADJ_alnsig = ADJ_xd  + ADJ_alnsig

   else if (ifreq .eq. 4) then 

      ADJ_z1  = s85h(1)*ADJ_alnthh
      ADJ_z2  = s85h(2)*ADJ_alnthh
      ADJ_z3  = s85h(3)*ADJ_alnthh
      ADJ_z4  = s85h(4)*ADJ_alnthh
      ADJ_z5  = s85h(5)*ADJ_alnthh
      ADJ_z6  = s85h(6)*ADJ_alnthh

      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg
      ADJ_xd  = 2.0*xd*ADJ_xx + ADJ_xd
      ADJ_alnsig = ADJ_xd

      ! vertical
      xd =  xd_save
      xx =  xx_save

      alnthv = alnthv_save

      ADJ_z1  = s85v(1)*ADJ_alnthv
      ADJ_z2  = s85v(2)*ADJ_alnthv
      ADJ_z3  = s85v(3)*ADJ_alnthv
      ADJ_z4  = s85v(4)*ADJ_alnthv
      ADJ_z5  = s85v(5)*ADJ_alnthv
      ADJ_z6  = s85v(6)*ADJ_alnthv

      ADJ_xd  = ADJ_z6*gx2
      ADJ_gx2 = xd*ADJ_z6  + ADJ_gx2
      ADJ_xx  = ADJ_z5*gx2
      ADJ_gx2 = xx*ADJ_z5  + ADJ_gx2
      ADJ_xx  = ADJ_z4*gg  + ADJ_xx
      ADJ_gg  = xx*ADJ_z4  + ADJ_gg
      ADJ_xd  = ADJ_z3*gg  + ADJ_xd
      ADJ_gg  = xd*ADJ_z3  + ADJ_gg
      ADJ_xd  = ADJ_z2*ggg + ADJ_xd
      ADJ_ggg = xd*ADJ_z2  + ADJ_ggg
      ADJ_xx  = ADJ_z1*ggg + ADJ_xx
      ADJ_ggg = xx*ADJ_z1  + ADJ_ggg
      ADJ_xd  = 2.0*xd*ADJ_xx  + ADJ_xd
      ADJ_alnsig = ADJ_xd  + ADJ_alnsig
   end if

   ADJ_gg  = ADJ_ggg*gx2   + ADJ_gg
   ADJ_gx2 = gg*ADJ_ggg    + ADJ_gx2
   ADJ_gx2 = 2.0*gx2*ADJ_gg + ADJ_gx2

   if (abs(sigma) .gt. 0.0) then
      ADJ_sigma = ADJ_alnsig/sigma + ADJ_sigma
   ! else
   !    ADJ_sigma = 0.0
   end if

   if (trace_use) call da_trace_exit("da_effang_adj")

end subroutine da_effang_adj


subroutine da_effht_adj(ho,hv,sigo,sigv,mu,zcld,hdn,hup,hdninf,hupinf, &
   ADJ_ho,ADJ_hv,ADJ_sigo,ADJ_sigv,ADJ_mu,        &
   ADJ_zcld,ADJ_hdn,ADJ_hup,ADJ_hdninf,ADJ_hupinf )

   implicit none

   !--------------------------------------------------------------------
   ! Purpose: TBD
   ! Output  : ADJ_ho, ADJ_hv, ADJ_sigo, ADJ_sigv, ADJ_zcld, ADJ_mu
   ! Input   : ADJ_hdn, ADJ_hup, ADJ_hdninf, ADJ_hupinf
   !--------------------------------------------------------------------

   real, intent(in)    :: ho,hv,sigo,sigv,mu,zcld
   real, intent(inout) :: ADJ_ho,ADJ_hv,ADJ_sigo,ADJ_sigv,ADJ_zcld, ADJ_mu
   real, intent(inout) :: hdn,hup,hdninf,hupinf
   real, intent(in)    :: ADJ_hdn,ADJ_hup,ADJ_hdninf,ADJ_hupinf

   real :: gint,zgint,hint,zhint
   real :: ginf,zginf,hinf,zhinf
   real :: ADJ_gint,ADJ_zgint,ADJ_hint,ADJ_zhint
   real :: ADJ_ginf,ADJ_zginf,ADJ_hinf,ADJ_zhinf
   real :: ADJ_mu2,ADJ_halfmu,ADJ_sixthmu2,ADJ_etnthmu2
   real :: ADJ_quartmu,ADJ_halfmu2

   real :: hoinv,hvinv,chio,chiv,ezho,ezhv,alpha,alph2,alph3
   real :: chio_save,chiv_save,dplus_save,dmin_save
   real :: beta,beta2,beta3,mu2,mualph,cplus,cmin,dplus,dmin
   real :: chiov,chivv,chioo,chioov,chiovv,chiooo,chivvv
   real :: h11,h21,h12,newh11
   real :: sigoo,sigov,sigvv,sigooo,sigoov,sigovv,sigvvv
   real :: ezhoo,ezhov,ezhvv,ezhooo,ezhoov,ezhovv,ezhvvv
   real :: s,sprim,t,tprim,u,uprim,term1,term2,term3
   real :: halfmu,halfmu2,sixthmu2,etnthmu2,quartmu

   real :: ADJ_hoinv,ADJ_hvinv,ADJ_chio,ADJ_chiv,ADJ_ezho
   real :: ADJ_ezhv,ADJ_alpha,ADJ_alph2,ADJ_alph3
   real :: ADJ_beta,ADJ_beta2,ADJ_beta3,ADJ_mualph
   real :: ADJ_cplus,ADJ_cmin,ADJ_dplus,ADJ_dmin
   real :: ADJ_chiov,ADJ_chivv,ADJ_chioo,ADJ_chioov
   real :: ADJ_chiovv,ADJ_chiooo,ADJ_chivvv
   real :: ADJ_h11,ADJ_h21,ADJ_h12,ADJ_newh11
   real :: ADJ_sigoo,ADJ_sigov,ADJ_sigvv,ADJ_sigooo
   real :: ADJ_sigoov,ADJ_sigovv,ADJ_sigvvv
   real :: ADJ_ezhoo,ADJ_ezhov,ADJ_ezhvv,ADJ_ezhooo
   real :: ADJ_ezhoov,ADJ_ezhovv,ADJ_ezhvvv
   real :: ADJ_s,ADJ_sprim,ADJ_t,ADJ_tprim
   real :: ADJ_u,ADJ_uprim,ADJ_term1,ADJ_term2,ADJ_term3

   if (trace_use) call da_trace_entry("da_effht_adj")

    ! initial zero
   gint=0.0
   zgint=0.0
   hint=0.0
   zhint=0.0
   ginf=0.0
   zginf=0.0
   hinf=0.0
   zhinf=0.0
   hoinv=0.0
   hvinv=0.0
   chio=0.0
   chiv=0.0
   ezho=0.0
   ezhv=0.0
   alpha=0.0
   alph2=0.0
   alph3=0.0
   chio_save=0.0
   chiv_save=0.0
   dplus_save=0.0
   dmin_save=0.0
   beta=0.0
   beta2=0.0
   beta3=0.0
   mu2=0.0
   mualph=0.0
   cplus=0.0
   cmin=0.0
   dplus=0.0
   dmin=0.0
   chiov=0.0
   chivv=0.0
   chioo=0.0
   chioov=0.0
   chiovv=0.0
   chiooo=0.0
   chivvv=0.0
   h11=0.0
   h21=0.0
   h12=0.0
   newh11=0.0
   sigoo=0.0
   sigov=0.0
   sigvv=0.0
   sigooo=0.0
   sigoov=0.0
   sigovv=0.0
   sigvvv=0.0
   ezhoo=0.0
   ezhov=0.0
   ezhvv=0.0
   ezhooo=0.0
   ezhoov=0.0
   ezhovv=0.0
   ezhvvv=0.0
   s=0.0
   sprim=0.0
   t=0.0
   tprim=0.0
   u=0.0
   uprim=0.0
   term1=0.0
   term2=0.0
   term3=0.0
   halfmu=0.0
   halfmu2=0.0
   sixthmu2=0.0
   etnthmu2=0.0
   quartmu=0.0

   ! WHY?
   ! ADJ_ho=0.0
   ! ADJ_hv=0.0
   ! ADJ_sigo=0.0
   ! ADJ_sigv=0.0
   ! ADJ_zcld=0.0
   ! ADJ_mu=0.0

   ADJ_gint=0.0
   ADJ_zgint=0.0
   ADJ_hint=0.0
   ADJ_zhint=0.0
   ADJ_ginf=0.0
   ADJ_zginf=0.0
   ADJ_hinf=0.0
   ADJ_zhinf=0.0
   ADJ_mu2=0.0
   ADJ_halfmu=0.0
   ADJ_sixthmu2=0.0
   ADJ_etnthmu2=0.0
   ADJ_quartmu=0.0
   ADJ_halfmu2=0.0
   ADJ_hoinv=0.0
   ADJ_hvinv=0.0
   ADJ_chio=0.0
   ADJ_chiv=0.0
   ADJ_ezho=0.0
   ADJ_ezhv=0.0
   ADJ_alpha=0.0
   ADJ_alph2=0.0
   ADJ_alph3=0.0
   ADJ_beta=0.0
   ADJ_beta2=0.0
   ADJ_beta3=0.0
   ADJ_mualph=0.0
   ADJ_cplus=0.0
   ADJ_cmin=0.0
   ADJ_dplus=0.0
   ADJ_dmin=0.0
   ADJ_chiov=0.0
   ADJ_chivv=0.0
   ADJ_chioo=0.0
   ADJ_chioov=0.0
   ADJ_chiovv=0.0
   ADJ_chiooo=0.0
   ADJ_chivvv=0.0
   ADJ_h11=0.0
   ADJ_h21=0.0
   ADJ_h12=0.0
   ADJ_newh11=0.0
   ADJ_sigoo=0.0
   ADJ_sigov=0.0
   ADJ_sigvv=0.0
   ADJ_sigooo=0.0
   ADJ_sigoov=0.0
   ADJ_sigovv=0.0
   ADJ_sigvvv=0.0
   ADJ_ezhoo=0.0
   ADJ_ezhov=0.0
   ADJ_ezhvv=0.0
   ADJ_ezhooo=0.0
   ADJ_ezhoov=0.0
   ADJ_ezhovv=0.0
   ADJ_ezhvvv=0.0
   ADJ_s=0.0
   ADJ_sprim=0.0
   ADJ_t=0.0
   ADJ_tprim=0.0
   ADJ_u=0.0
   ADJ_uprim=0.0
   ADJ_term1=0.0
   ADJ_term2=0.0
   ADJ_term3=0.0

   ! base fields

   hoinv =  1.0d0/ho
   hvinv =  1.0d0/hv
   chio = zcld*hoinv
   chiv = zcld*hvinv
   ezho = sigo*exp(-chio)
   ezhv = sigv*exp(-chiv)
   alpha = sigo + sigv
   alph2 = alpha*alpha
   alph3 = alpha*alph2
   beta = ezho + ezhv
   beta2 = beta*beta
   beta3 = beta*beta2

   mu2        = mu*mu
   halfmu     = 0.5d0*    mu
   sixthmu2   =     mu2/6.0d0
   etnthmu2   =     mu2/18.0d0
   quartmu    = 0.25d0*    mu
   halfmu2    = 0.5d0*    mu2

   mualph = mu*alpha
   cplus  = 1.0d0 +     mualph
   cmin   = 1.0d0 -     mualph
   dplus  = halfmu2*alph2
   dmin   =     dplus
  
   dplus_save = dplus
   dplus  =     cplus +     dplus
  
   dmin_save  = dmin
   dmin   =     cmin  +     dmin

   h11    =     hoinv +     hvinv
   h21    =  1.0d0/(h11 + hvinv)
   h12    =  1.0d0/(h11 + hoinv)
   newh11 =  1.0d0/h11
  
   chiov  = 1.0d0 +     chio +     chiv
   chioo  = 1.0d0 +     chio +     chio
   chivv  = 1.0d0 +     chiv +     chiv
   chioov =     chioo +     chiv
   chiovv =     chio  +     chivv
   chiooo =     chioo +     chio
   chivvv =     chivv +     chiv
   chio_save = chio
   chio   = 1.0d0 +     chio
   chiv_save = chiv
   chiv   = 1.0d0 +     chiv
   sigov  = sigo*sigv
   sigoo  = sigo*sigo
   sigvv  = sigv*sigv
   sigooo = sigoo*sigo
   sigoov = sigoo*sigv
   sigovv = sigo*sigvv
   sigvvv = sigvv*sigv
   ezhoo  = ezho*ezho
   ezhov  = ezho*ezhv
   ezhvv  = ezhv*ezhv
   ezhovv = ezho*ezhvv
   ezhoov = ezhoo*ezhv
   ezhooo = ezhoo*ezho
   ezhvvv = ezhvv*ezhv
   s      = sigo*ho + sigv*hv

   sprim  = ezho*ho*chio + ezhv*hv*chiv
   t      = sigoo*ho + 4.0d0*sigov*newh11 + sigvv*hv

   tprim  = ezhoo*ho*chioo + 4.0d0*ezhov*newh11*chiov + ezhvv*hv*chivv
   u      = sigooo*ho + 9.0d0*(sigovv*h21+sigoov*h12) + sigvvv*hv

   uprim  = ezhvvv*hv*chivvv +  &
            9.0d0*(ezhovv*h21*chiovv + ezhoov*h12*chioov) + &
            ezhooo*ho*chiooo

   term1  =     s -     sprim
   term2  = quartmu*(t - tprim)
   term3  = etnthmu2*(    u -     uprim)
   zgint  = dmin*term1 +  cmin*term2 + term3
   zhint  = -dplus*term1 + cplus*term2 - term3
   term2  = quartmu * t
   term3  = etnthmu2*u
   zginf  = dmin*s +  cmin*term2 + term3
   zhinf  = -dplus*s + cplus*term2 - term3

   term1  =     alpha -     beta
   term2  = halfmu*(    alph2 -     beta2)
   term3  = sixthmu2*(    alph3 -     beta3)
   gint   = dmin*term1 +  cmin*term2 + term3
   hint   = -dplus*term1 + cplus*term2 - term3
   term2  = halfmu*alph2
   term3  = sixthmu2*alph3
   ginf   = dmin*alpha +  cmin*term2 + term3
   hinf   = -dplus*alpha + cplus*term2 - term3
   hdn    = zgint/gint
   hup    = zhint/hint
   hdninf = zginf/ginf
   hupinf = zhinf/hinf

   ! start

   ADJ_zhinf  =   ADJ_hupinf/hinf        + ADJ_zhinf
   ADJ_hinf   = - hupinf*ADJ_hupinf/hinf + ADJ_hinf

   ADJ_zginf  =   ADJ_hdninf/ginf        + ADJ_zginf
   ADJ_ginf   = - hdninf*ADJ_hdninf/ginf + ADJ_ginf

   ADJ_zhint  =   ADJ_hup/hint     + ADJ_zhint
   ADJ_hint   = - hup*ADJ_hup/hint + ADJ_hint

   ADJ_zgint  =   ADJ_hdn/gint       + ADJ_zgint
   ADJ_gint   = - hdn * ADJ_hdn/gint + ADJ_gint

   ADJ_dplus  = - ADJ_hinf*alpha  + ADJ_dplus 
   ADJ_alpha = - dplus*ADJ_hinf   + ADJ_alpha
   ADJ_cplus =   ADJ_hinf*term2   + ADJ_cplus
   ADJ_term2 =   cplus*ADJ_hinf 
   ADJ_term3 = - ADJ_hinf         

   ADJ_dmin   =   ADJ_ginf*alpha + ADJ_dmin
   ADJ_alpha  =   dmin*ADJ_ginf  + ADJ_alpha
   ADJ_cmin   =   ADJ_ginf*term2 + ADJ_cmin
   ADJ_term2  =   cmin*ADJ_ginf  + ADJ_term2
   ADJ_term3  =   ADJ_ginf       + ADJ_term3

   ADJ_sixthmu2 = ADJ_term3*alph3    + ADJ_sixthmu2
   ADJ_alph3    = sixthmu2*ADJ_term3 + ADJ_alph3

   ADJ_halfmu =  ADJ_term2*alph2  + ADJ_halfmu
   ADJ_alph2  =  halfmu*ADJ_term2 + ADJ_alph2

   ! new term2,3

   term2 = halfmu*(alph2 - beta2)
   term3 = sixthmu2*(alph3 - beta3)

   ADJ_dplus  = - ADJ_hint*term1 + ADJ_dplus
   ADJ_term1  = - dplus*ADJ_hint 
   ADJ_cplus  =   ADJ_hint*term2 + ADJ_cplus
   ADJ_term2  =   cplus*ADJ_hint 
   ADJ_term3  = - ADJ_hint

   ADJ_dmin   = ADJ_gint*term1 + ADJ_dmin
   ADJ_term1  = dmin*ADJ_gint  + ADJ_term1
   ADJ_cmin   = ADJ_gint*term2 + ADJ_cmin
   ADJ_term2  = cmin*ADJ_gint  + ADJ_term2
   ADJ_term3  = ADJ_gint       + ADJ_term3

   ADJ_sixthmu2 =    ADJ_term3*(alph3 - beta3) + ADJ_sixthmu2
   ADJ_alph3    =    sixthmu2*ADJ_term3        + ADJ_alph3
   ADJ_beta3    =  - sixthmu2*ADJ_term3        + ADJ_beta3
    
   ADJ_halfmu =   ADJ_term2*(alph2 - beta2) + ADJ_halfmu
   ADJ_alph2  =   halfmu*ADJ_term2          + ADJ_alph2
   ADJ_beta2  = - halfmu*ADJ_term2          + ADJ_beta2

   ADJ_alpha  =   ADJ_term1 + ADJ_alpha
   ADJ_beta   = - ADJ_term1 + ADJ_beta

   ! new term2,3

   term2 = quartmu*t
   term3 = etnthmu2*u

   ADJ_dplus = - ADJ_zhinf*s     + ADJ_dplus
   ADJ_s = - dplus*ADJ_zhinf + ADJ_s
   ADJ_cplus =   ADJ_zhinf*term2 + ADJ_cplus
   ADJ_term2 =   cplus*ADJ_zhinf
   ADJ_term3 = - ADJ_zhinf 

   ADJ_dmin   = ADJ_zginf*s     + ADJ_dmin
   ADJ_s      = dmin*ADJ_zginf  + ADJ_s
   ADJ_cmin   = ADJ_zginf*term2 + ADJ_cmin
   ADJ_term2  = cmin*ADJ_zginf  + ADJ_term2
   ADJ_term3  = ADJ_zginf       + ADJ_term3
   
   ADJ_etnthmu2 = ADJ_term3*u        + ADJ_etnthmu2
   ADJ_u        = etnthmu2*ADJ_term3 + ADJ_u
   ADJ_quartmu = ADJ_term2*t         + ADJ_quartmu   
   ADJ_t       = quartmu*ADJ_term2   + ADJ_t

   ! new term1,2,3

   term1 = s - sprim
   term2 = quartmu*(t - tprim)
   term3 = etnthmu2*(u - uprim)

   ADJ_dplus = - ADJ_zhint*term1 + ADJ_dplus
   ADJ_term1 = - dplus*ADJ_zhint
   ADJ_cplus =   ADJ_zhint*term2 + ADJ_cplus
   ADJ_term2 =   cplus*ADJ_zhint
   ADJ_term3 = - ADJ_zhint

   ADJ_dmin  = ADJ_zgint*term1 + ADJ_dmin
   ADJ_term1  = dmin*ADJ_zgint  + ADJ_term1
   ADJ_cmin  = ADJ_zgint*term2 + ADJ_cmin
   ADJ_term2  = cmin*ADJ_zgint  + ADJ_term2
   ADJ_term3  = ADJ_zgint       + ADJ_term3

   ADJ_etnthmu2 =   ADJ_term3*(u - uprim) + ADJ_etnthmu2
   ADJ_u        =   etnthmu2*ADJ_term3    + ADJ_u
   ADJ_uprim    = - etnthmu2*ADJ_term3    + ADJ_uprim

   ADJ_quartmu =   ADJ_term2*(t - tprim) + ADJ_quartmu
   ADJ_t       =   quartmu*ADJ_term2     + ADJ_t
   ADJ_tprim   = - quartmu*ADJ_term2     + ADJ_tprim 

   ADJ_s      =   ADJ_term1 + ADJ_s
   ADJ_sprim  = - ADJ_term1 + ADJ_sprim


   ADJ_ezhvvv = ADJ_uprim*hv*chivvv            + ADJ_ezhvvv
   ADJ_hv     = ezhvvv*ADJ_uprim*chivvv        + ADJ_hv
   ADJ_chivvv = ezhvvv*hv*ADJ_uprim            + ADJ_chivvv
   ADJ_ezhovv = 9.0d0*ADJ_uprim*h21*chiovv     + ADJ_ezhovv
   ADJ_h21    = 9.0d0*ezhovv*ADJ_uprim*chiovv  + ADJ_h21
   ADJ_chiovv = 9.0d0*ezhovv*h21*ADJ_uprim     + ADJ_chiovv
   ADJ_ezhoov = 9.0d0*ADJ_uprim*h12*chioov     + ADJ_ezhoov
   ADJ_h12    = 9.0d0*ezhoov*ADJ_uprim*chioov  + ADJ_h12
   ADJ_chioov = 9.0d0*ezhoov*h12*ADJ_uprim     + ADJ_chioov
   ADJ_ezhooo = ADJ_uprim*ho*chiooo            + ADJ_ezhooo
   ADJ_ho     = ezhooo*ADJ_uprim*chiooo        + ADJ_ho
   ADJ_chiooo = ezhooo*ho*ADJ_uprim            + ADJ_chiooo

   ADJ_sigooo = ADJ_u*ho                       + ADJ_sigooo
   ADJ_ho     = sigooo*ADJ_u                   + ADJ_ho 
   ADJ_sigovv = 9.0d0*ADJ_u*h21                + ADJ_sigovv
   ADJ_h21    = 9.0d0*sigovv*ADJ_u             + ADJ_h21
   ADJ_sigoov = 9.0d0*ADJ_u*h12                + ADJ_sigoov
   ADJ_h12    = 9.0d0*sigoov*ADJ_u             + ADJ_h12
   ADJ_sigvvv = ADJ_u*hv                       + ADJ_sigvvv
   ADJ_hv     = sigvvv*ADJ_u                   + ADJ_hv
 
   ADJ_ezhoo  = ADJ_tprim*ho*chioo             + ADJ_ezhoo
   ADJ_ho     = ezhoo*ADJ_tprim*chioo          + ADJ_ho
   ADJ_chioo  = ezhoo*ho*ADJ_tprim             + ADJ_chioo
   ADJ_ezhov  = 4.0d0*ADJ_tprim*newh11*chiov   + ADJ_ezhov
   ADJ_newh11 = 4.0d0*ezhov*ADJ_tprim*chiov    + ADJ_newh11
   ADJ_chiov  = 4.0d0*ezhov*newh11*ADJ_tprim   + ADJ_chiov
   ADJ_ezhvv  = ADJ_tprim*hv*chivv             + ADJ_ezhvv
   ADJ_hv     = ezhvv*ADJ_tprim*chivv          + ADJ_hv
   ADJ_chivv  = ezhvv*hv*ADJ_tprim             + ADJ_chivv

   ADJ_sigoo  = ADJ_t*ho           + ADJ_sigoo
   ADJ_ho     = sigoo*ADJ_t        + ADJ_ho
   ADJ_sigov  = 4.0d0*ADJ_t*newh11 + ADJ_sigov
   ADJ_newh11 = 4.0d0*sigov*ADJ_t  + ADJ_newh11
   ADJ_sigvv  = ADJ_t*hv           + ADJ_sigvv
   ADJ_hv     = sigvv*ADJ_t        + ADJ_hv

   ADJ_ezho   = ADJ_sprim*ho*chio   + ADJ_ezho
   ADJ_ho     = ezho*ADJ_sprim*chio + ADJ_ho
   ADJ_chio   = ezho*ho*ADJ_sprim   + ADJ_chio
   ADJ_ezhv   = ADJ_sprim*hv*chiv   + ADJ_ezhv
   ADJ_hv     = ezhv*ADJ_sprim*chiv + ADJ_hv
   ADJ_chiv   = ezhv*hv*ADJ_sprim   + ADJ_chiv

   ADJ_sigo   = ADJ_s*ho     + ADJ_sigo
   ADJ_ho     = sigo*ADJ_s   + ADJ_ho
   ADJ_sigv   = ADJ_s*hv     + ADJ_sigv
   ADJ_hv     = sigv*ADJ_s   + ADJ_hv

   ADJ_ezhvv  = ADJ_ezhvvv*ezhv   + ADJ_ezhvv
   ADJ_ezhv   = ezhvv*ADJ_ezhvvv  + ADJ_ezhv
   ADJ_ezhoo  = ADJ_ezhooo*ezho   + ADJ_ezhoo
   ADJ_ezho   = ezhoo*ADJ_ezhooo  + ADJ_ezho
   ADJ_ezhoo  = ADJ_ezhoov*ezhv   + ADJ_ezhoo 
   ADJ_ezhv   = ezhoo*ADJ_ezhoov  + ADJ_ezhv
   ADJ_ezho   = ADJ_ezhovv*ezhvv  + ADJ_ezho
   ADJ_ezhvv  = ezho*ADJ_ezhovv   + ADJ_ezhvv
   ADJ_ezhv   = 2.0*ezhv*ADJ_ezhvv + ADJ_ezhv
   ADJ_ezho   = ADJ_ezhov*ezhv    + ADJ_ezho
   ADJ_ezhv   = ezho*ADJ_ezhov    + ADJ_ezhv
   ADJ_ezho   = 2.0*ezho*ADJ_ezhoo + ADJ_ezho
   ADJ_sigvv  = ADJ_sigvvv*sigv   + ADJ_sigvv
   ADJ_sigv   = sigvv*ADJ_sigvvv  + ADJ_sigv
   ADJ_sigo   = ADJ_sigovv*sigvv  + ADJ_sigo
   ADJ_sigvv  = sigo*ADJ_sigovv   + ADJ_sigvv
   ADJ_sigoo  = ADJ_sigoov*sigv   + ADJ_sigoo
   ADJ_sigv   = sigoo*ADJ_sigoov  + ADJ_sigv
   ADJ_sigoo  = ADJ_sigooo*sigo   + ADJ_sigoo
   ADJ_sigo   = sigoo*ADJ_sigooo  + ADJ_sigo
   ADJ_sigv   = 2.0*sigv*ADJ_sigvv + ADJ_sigv
   ADJ_sigo   = 2.0*sigo*ADJ_sigoo + ADJ_sigo
   ADJ_sigo   = ADJ_sigov*sigv    + ADJ_sigo
   ADJ_sigv   = sigo*ADJ_sigov    + ADJ_sigv

   ! WHY?
   !  ADJ_chiv   =         ADJ_chiv 
   ! ADJ_chio   =         ADJ_chio 

   ! new chio chiv

   chio   = chio_save
   chiv   = chiv_save

   ADJ_chivv  = ADJ_chivvv + ADJ_chivv
   ADJ_chiv   = ADJ_chivvv + ADJ_chiv

   ADJ_chioo  = ADJ_chiooo + ADJ_chioo
   ADJ_chio   = ADJ_chiooo + ADJ_chio

   ADJ_chio   = ADJ_chiovv + ADJ_chio
   ADJ_chivv  = ADJ_chiovv + ADJ_chivv

   ADJ_chioo  = ADJ_chioov + ADJ_chioo
   ADJ_chiv   = ADJ_chioov + ADJ_chiv

   ADJ_chiv   = ADJ_chivv + ADJ_chiv
   ADJ_chiv   = ADJ_chivv + ADJ_chiv

   ADJ_chio   = ADJ_chioo + ADJ_chio
   ADJ_chio   = ADJ_chioo + ADJ_chio

   ADJ_chio   =  ADJ_chiov + ADJ_chio
   ADJ_chiv   =  ADJ_chiov + ADJ_chiv

   ADJ_h11    = -1.0d0*newh11*newh11*ADJ_newh11 + ADJ_h11

   ADJ_h11    = -1.0d0*h12*h12*ADJ_h12 + ADJ_h11
   ADJ_hoinv  = -1.0d0*h12*h12*ADJ_h12 + ADJ_hoinv

   ADJ_h11    = -1.0d0*h21*h21*ADJ_h21 + ADJ_h11
   ADJ_hvinv  = -1.0d0*h21*h21*ADJ_h21 + ADJ_hvinv

   ADJ_hoinv  = ADJ_h11 + ADJ_hoinv
   ADJ_hvinv  = ADJ_h11 + ADJ_hvinv

   ADJ_cmin   = ADJ_dmin  + ADJ_cmin
   ! ADJ_dmin   = ADJ_dmin
   dmin   = dmin_save
 
   ADJ_cplus  = ADJ_dplus + ADJ_cplus
   ! ADJ_dplus  = ADJ_dplus
   dplus      = dplus_save

   ADJ_dplus  = ADJ_dmin           + ADJ_dplus
   ADJ_halfmu2 = ADJ_dplus*alph2   + ADJ_halfmu2
   ADJ_alph2   = halfmu2*ADJ_dplus + ADJ_alph2
   ADJ_mualph = - ADJ_cmin         + ADJ_mualph
   ADJ_mualph = ADJ_cplus          + ADJ_mualph
   ADJ_mu     = ADJ_mualph*alpha   + ADJ_mu
   ADJ_alpha  = mu*ADJ_mualph      + ADJ_alpha

   ADJ_mu2  = 0.5d0*ADJ_halfmu2   + ADJ_mu2
   ADJ_mu   = 0.25d0*ADJ_quartmu  + ADJ_mu
   ADJ_mu2  = ADJ_etnthmu2/18.0d0 + ADJ_mu2
   ADJ_mu2  = ADJ_sixthmu2/6.0d0  + ADJ_mu2
   ADJ_mu   = 0.5d0*ADJ_halfmu    + ADJ_mu
   ADJ_mu   = 2.0*mu*ADJ_mu2       + ADJ_mu

   ADJ_beta = ADJ_beta3*beta2               + ADJ_beta
   ADJ_beta2 = beta*ADJ_beta3               + ADJ_beta2
   ADJ_beta = 2.0*beta*ADJ_beta2             + ADJ_beta
   ADJ_ezho = ADJ_beta                      + ADJ_ezho
   ADJ_ezhv = ADJ_beta                      + ADJ_ezhv
   ADJ_alpha = ADJ_alph3*alph2              + ADJ_alpha
   ADJ_alph2 = alpha*ADJ_alph3              + ADJ_alph2
   ADJ_alpha = 2.0*alpha*ADJ_alph2           + ADJ_alpha
   ADJ_sigo = ADJ_alpha                     + ADJ_sigo
   ADJ_sigv = ADJ_alpha                     + ADJ_sigv
   ADJ_sigv =  ADJ_ezhv*exp(-chiv)          + ADJ_sigv
   ADJ_chiv = -ADJ_ezhv*ezhv                + ADJ_chiv
   ADJ_sigo =  ADJ_ezho*exp(-chio)          + ADJ_sigo
   ADJ_chio = -ADJ_ezho*ezho                + ADJ_chio
   ADJ_zcld = ADJ_chiv*hvinv                + ADJ_zcld
   ADJ_hvinv= zcld*ADJ_chiv                 + ADJ_hvinv
   ADJ_zcld  = ADJ_chio*hoinv               + ADJ_zcld
   ADJ_hoinv = zcld*ADJ_chio                + ADJ_hoinv
   ADJ_hv    = -1.0d0*hvinv*hvinv*ADJ_hvinv + ADJ_hv
   ADJ_ho    = -1.0d0*hoinv*hoinv*ADJ_hoinv + ADJ_ho

   if (trace_use) call da_trace_exit("da_effht_adj")

end subroutine da_effht_adj


subroutine da_epsalt_adj(f,t,ssw,ADJ_t, ADJ_epsr, ADJ_epsi)

   implicit none

   !---------------------------------------------------------------------------
   ! Purpose: TBD
   ! Output: ADJ_t      (ssw is treated as a constant now)
   ! Input: ADJ_epsr, ADJ_epsi, epsr, epsi
   !---------------------------------------------------------------------------

   real, intent(in)    :: f, t
   real, intent(inout) :: ADJ_t
   real, intent(inout) :: ssw
   real, intent(in)    :: ADJ_epsr, ADJ_epsi

   complex :: cdum1,cdum2,cdum3
   complex :: ADJ_cdum1,ADJ_cdum2,ADJ_cdum3
   real    :: ssw2,ssw3,t2,t3,es,a,esnew,tau,b,sig,taunew
   real    :: delt,delt2,beta,signew,om,d1,d2
   real    :: ADJ_t2,ADJ_t3,ADJ_es,ADJ_a,ADJ_esnew,ADJ_tau,ADJ_b,ADJ_taunew
   real    :: ADJ_delt,ADJ_delt2,ADJ_beta,ADJ_signew
   real    :: ADJ_d1,ADJ_d2

   if (trace_use) call da_trace_entry("da_epsalt_adj")

   ssw2=0.0
   ssw3=0.0 
   t2=0.0
   t3=0.0
   es=0.0
   a=0.0
   esnew=0.0
   tau=0.0
   b=0.0
   sig=0.0
   taunew=0.0
   delt=0.0
   delt2=0.0
   beta=0.0 
   signew=0.0
   om=0.0
   d1=0.0
   d2=0.0
   ADJ_t2=0.0
   ADJ_t3=0.0
   ADJ_es=0.0
   ADJ_a=0.0
   ADJ_esnew=0.0
   ADJ_tau=0.0
   ADJ_b=0.0
   ADJ_taunew=0.0
   ADJ_delt=0.0
   ADJ_delt2=0.0
   ADJ_beta=0.0
   ADJ_signew=0.0
   ADJ_d1=0.0
   ADJ_d2=0.0

   if (ssw .lt. 0.0) ssw = 32.54

   ssw2     = ssw*ssw
   ssw3     = ssw2*ssw
   t2     = t*t
   t3     = t2*t
   es     = 87.134 - 1.949e-1*t - 1.276e-2*t2 + 2.491e-4*t3
   a      = 1.0 + 1.613e-5*ssw*t - 3.656e-3*ssw + 3.21e-5*ssw2 - &
            4.232e-7*ssw3
   esnew  = es*a

   tau    = 1.768e-11 - 6.086e-13*t + 1.104e-14*t2 - 8.111e-17*t3
   b      = 1.0 + 2.282e-5*ssw*t - 7.638e-4*ssw - 7.760e-6*ssw2 + &
            1.105e-8*ssw3
   taunew = tau*b

   sig    = ssw*(0.182521 - 1.46192e-3*ssw + 2.09324e-5*ssw2 - &
            1.28205e-7*ssw3)
   delt   = 25.0 - t
   delt2  = delt*delt
   beta   =   2.033e-2 + 1.266e-4*delt      + 2.464e-6*delt2       &
            - ssw*(1.849e-5 - 2.551e-7*delt + 2.551e-8*delt2)
   signew  =   sig*exp(-beta*delt)

   om  = 2.0e9*pi*f
   cdum1  = cmplx(0.0,om*taunew)
   cdum2  = cmplx(0.0,signew/(om*8.854e-12))

   cdum3  = 4.9 + (esnew-4.9)/(1.0 + cdum1) - cdum2

   ADJ_cdum3  = ADJ_epsr + ADJ_epsi *(0.0,1.0)
   ADJ_esnew  =   REAL(ADJ_cdum3/((1.0,0.0) + cdum1))
   ADJ_cdum1  = - ADJ_cdum3*(esnew-4.9)/((1.0 + cdum1)*(1.0 + cdum1))
   ADJ_cdum2  = - ADJ_cdum3


   ADJ_signew = -aimag(ADJ_cdum2/(om*8.854e-12)) 

   ADJ_taunew = om*(-aimag(ADJ_cdum1))

   ADJ_beta   = - signew*ADJ_signew*delt
   ADJ_delt   = - signew*beta*ADJ_signew

   ADJ_delt   =   1.266e-4*ADJ_beta     + ADJ_delt
   ADJ_delt2  =   2.464e-6*ADJ_beta
   ADJ_delt   =   ssw*2.551e-7*ADJ_beta + ADJ_delt
   ADJ_delt2  = - ssw*2.551e-8*ADJ_beta + ADJ_delt2

   ADJ_delt   = 2.0*delt*ADJ_delt2 + ADJ_delt

   ADJ_t      =  - ADJ_delt + ADJ_t

   ADJ_tau    = ADJ_taunew*b 
   ADJ_b      = tau*ADJ_taunew

   ADJ_t      = 2.282e-5*ssw*ADJ_b  + ADJ_t

   ADJ_t      = - 6.086e-13*ADJ_tau + ADJ_t
   ADJ_t2     =   1.104e-14*ADJ_tau
   ADJ_t3     = - 8.111e-17*ADJ_tau

   ADJ_es     = ADJ_esnew*a
   ADJ_a      = es*ADJ_esnew
   ADJ_t      = 1.613e-5*ssw*ADJ_a  + ADJ_t
   ADJ_t      = - 1.949e-1*ADJ_es   + ADJ_t
   ADJ_t2     = - 1.276e-2*ADJ_es   + ADJ_t2
   ADJ_t3     =   2.491e-4*ADJ_es   + ADJ_t3

   ADJ_t2     = ADJ_t3*t     + ADJ_t2
   ADJ_t      = t2*ADJ_t3    + ADJ_t
   ADJ_t      = 2.0*t*ADJ_t2  + ADJ_t

   if (trace_use) call da_trace_exit("da_epsalt_adj")

end subroutine da_epsalt_adj


subroutine da_roughem_adj (ifreq,gx2,tk,theta,remv,remh,ADJ_gx2,ADJ_tk,ADJ_remv,ADJ_remh)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   !-----------------------------------------------------------------
   ! Input  :: ADJ_tk, ADJ_gx2
   ! Output :: ADJ_remv,ADJ_remh, remv,remh
   !
   ! Calculates rough-surface emissivity of ocean surface at SSM/I
   ! frequencies.

   integer, intent(in)    :: ifreq
   real   , intent(in)    :: tk, theta, gx2
   real   , intent(in)    :: ADJ_remv,ADJ_remh
   real   , intent(inout) :: ADJ_tk,ADJ_gx2
   real,    intent(out)   :: remv,remh

   real :: ssw

   real :: a19v(4),a22v(4),a37v(4),a85v(4)
   real :: a19h(4),a22h(4),a37h(4),a85h(4)
   real :: f(4)
   real :: semv,semh,ADJ_semv,ADJ_semh,remv_save,remh_save
   real :: tp,g,x1,x2,x3,x4,dtheta
   real :: ADJ_tp,ADJ_g,ADJ_x1,ADJ_x2,ADJ_x3,ADJ_x4

   data a19v/  -0.111E+01,   0.713E+00,  -0.624E-01,   0.212E-01 /
   data a19h/   0.812E+00,  -0.215E+00,   0.255E-01,   0.305E-02 /
   data a22v/  -0.134E+01,   0.911E+00,  -0.893E-01,   0.463E-01 /
   data a22h/   0.958E+00,  -0.350E+00,   0.566E-01,  -0.262E-01 /
   data a37v/  -0.162E+01,   0.110E+01,  -0.730E-01,   0.298E-01 /
   data a37h/   0.947E+00,  -0.320E+00,   0.624E-01,  -0.300E-01 /
   data a85v/  -0.145E+01,   0.808E+00,  -0.147E-01,  -0.252E-01 /
   data a85h/   0.717E+00,  -0.702E-01,   0.617E-01,  -0.243E-01 /

   data f/ 19.35, 22.235, 37.0, 85.5 /

   if (trace_use) call da_trace_entry("da_roughem_adj")

   semv      = 0.0
   semh      = 0.0
   ADJ_semv  = 0.0
   ADJ_semh  = 0.0
   remv_save = 0.0
   remh_save = 0.0
   tp        = 0.0
   g         = 0.0
   x1        = 0.0
   x2        = 0.0
   x3        = 0.0
   x4        = 0.0
   dtheta    = 0.0
   ADJ_tp    = 0.0
   ADJ_g     = 0.0
   ADJ_x1    = 0.0
   ADJ_x2    = 0.0
   ADJ_x3    = 0.0
   ADJ_x4    = 0.0

   tp     = tk/t_roughem
   dtheta = theta-53.0
   g      = 0.5* gx2
   x1     = g
   x2     = tp*g
   x3     = dtheta* g
   x4     = tp*x3

   if (ifreq == 1) then
      remv = x1*a19v(1) +     x2*a19v(2) +     x3*a19v(3) +     x4*a19v(4)
      remh = x1*a19h(1) +     x2*a19h(2) +     x3*a19h(3) +     x4*a19h(4)
   else if (ifreq == 2) then
      remv = x1*a22v(1) +     x2*a22v(2) +     x3*a22v(3) +     x4*a22v(4)
      remh = x1*a22h(1) +     x2*a22h(2) +     x3*a22h(3) +     x4*a22h(4)
   else if (ifreq == 3) then
      remv = x1*a37v(1) +     x2*a37v(2) +     x3*a37v(3) +     x4*a37v(4)
      remh = x1*a37h(1) +     x2*a37h(2) +     x3*a37h(3) +     x4*a37h(4)
   else if (ifreq == 4) then
      remv = x1*a85v(1) +     x2*a85v(2) +     x3*a85v(3) +     x4*a85v(4)
      remh = x1*a85h(1) +     x2*a85h(2) +     x3*a85h(3) +     x4*a85h(4)
   end if

   ssw=36.5
   call spemiss(f(ifreq),tk,theta,ssw,semv,semh)

   remv_save = remv
   remh_save = remh

   ! start

   ADJ_semh = ADJ_remh
     
   ADJ_semv = ADJ_remv

   remv = remv_save
   remh = remh_save

   call da_spemiss_adj(f(ifreq),tk,theta,ssw,semv,semh, ADJ_tk,ADJ_semv,ADJ_semh          )

 
   if (ifreq == 1) then
      ADJ_x1 = ADJ_remh*a19h(1)
      ADJ_x2 = ADJ_remh*a19h(2)
      ADJ_x3 = ADJ_remh*a19h(3)
      ADJ_x4 = ADJ_remh*a19h(4)
 
      ADJ_x1 = ADJ_remv*a19v(1) + ADJ_x1
      ADJ_x2 = ADJ_remv*a19v(2) + ADJ_x2
      ADJ_x3 = ADJ_remv*a19v(3) + ADJ_x3
      ADJ_x4 = ADJ_remv*a19v(4) + ADJ_x4
   else if (ifreq == 2) then
      ADJ_x1 = ADJ_remh*a22h(1) 
      ADJ_x2 = ADJ_remh*a22h(2)
      ADJ_x3 = ADJ_remh*a22h(3)
      ADJ_x4 = ADJ_remh*a22h(4)

      ADJ_x1 = ADJ_remv*a22v(1) + ADJ_x1
      ADJ_x2 = ADJ_remv*a22v(2) + ADJ_x2
      ADJ_x3 = ADJ_remv*a22v(3) + ADJ_x3
      ADJ_x4 = ADJ_remv*a22v(4) + ADJ_x4
   else if (ifreq == 3) then
      ADJ_x1 = ADJ_remh*a37h(1)
      ADJ_x2 = ADJ_remh*a37h(2)
      ADJ_x3 = ADJ_remh*a37h(3)
      ADJ_x4 = ADJ_remh*a37h(4)         

      ADJ_x1 = ADJ_remv*a37v(1) + ADJ_x1
      ADJ_x2 = ADJ_remv*a37v(2) + ADJ_x2
      ADJ_x3 = ADJ_remv*a37v(3) + ADJ_x3
      ADJ_x4 = ADJ_remv*a37v(4) + ADJ_x4
   else if (ifreq == 4) then
      ADJ_x1 = ADJ_remh*a85h(1) 
      ADJ_x2 = ADJ_remh*a85h(2) 
      ADJ_x3 = ADJ_remh*a85h(3)
      ADJ_x4 = ADJ_remh*a85h(4) 
 
      ADJ_x1 = ADJ_remv*a85v(1) + ADJ_x1
      ADJ_x2 = ADJ_remv*a85v(2) + ADJ_x2
      ADJ_x3 = ADJ_remv*a85v(3) + ADJ_x3
      ADJ_x4 = ADJ_remv*a85v(4) + ADJ_x4
     end if

   ADJ_tp  = ADJ_x4*x3
   ADJ_x3  = tp*ADJ_x4     + ADJ_x3
   ADJ_g   = dtheta*ADJ_x3
   ADJ_tp  = ADJ_x2*g      + ADJ_tp
   ADJ_g   = tp*ADJ_x2     + ADJ_g
   ADJ_g   = ADJ_x1        + ADJ_g
   ADJ_gx2 = 0.5*ADJ_g     + ADJ_gx2
   ADJ_tk  = ADJ_tp/t_roughem  + ADJ_tk

   if (trace_use) call da_trace_exit("da_roughem_adj")
 
end subroutine da_roughem_adj


subroutine da_spemiss_adj(f,tk,theta,ssw,ev,eh, ADJ_tk,ADJ_ev,ADJ_eh)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   !------------------------------------------------------------------------
   ! Output :: ADJ_tk
   ! Input  :: ADJ_ev, ADJ_eh
   !------------------------------------------------------------------------

   real, intent(in ) :: f, tk, theta, ADJ_ev,ADJ_eh
   real, intent(inout) :: ssw
   real, intent(inout) :: ADJ_tk
   real, intent(out)   :: ev, eh

   real   epsr,epsi,ADJ_epsr,ADJ_epsi

   real      tc,costh,sinth,rthet
   complex   etav,etah,eps,cterm1v,cterm1h,cterm2,cterm3v,cterm3h,epsnew
   complex   ADJ_etav,ADJ_eps,ADJ_cterm1v,ADJ_cterm2,ADJ_cterm3v
   complex   ADJ_cterm3h,ADJ_epsnew
   real      tmp1r,tmp1i,tmp2r,tmp2i,tmp0r,tmp0i,rnorm
   real      ADJ_tc,ADJ_tmp0r,ADJ_tmp0i,ADJ_rnorm,ADJ_tmp1r
   real      ADJ_tmp1i,ADJ_tmp2r,ADJ_tmp2i

   if (trace_use) call da_trace_entry("da_spemiss_adj")

   epsr=0.0
   epsi=0.0
   ADJ_epsr=0.0
   ADJ_epsi=0.0
   ev=0.0
   eh=0.0
   tc=0.0
   costh=0.0
   sinth=0.0
   rthet=0.0
   tmp1r=0.0
   tmp1i=0.0
   tmp2r=0.0
   tmp2i=0.0
   tmp0r=0.0
   tmp0i=0.0
   rnorm=0.0
   ADJ_tc=0.0
   ADJ_tmp0r=0.0
   ADJ_tmp0i=0.0
   ADJ_rnorm=0.0
   ADJ_tmp1r=0.0
   ADJ_tmp1i=0.0
   ADJ_tmp2r=0.0
   ADJ_tmp2i=0.0

   tc     =      tk - t_kelvin

   call epsalt(f,tc,ssw,epsr,epsi)

   eps     =  cmplx(epsr,epsi)
   etav    =  eps
   etah    =  (1.0,0.0)
   rthet   =  theta*0.017453292
   costh   =  cos(rthet)
   sinth   =  sin(rthet)
   sinth   =  sinth*sinth
   cterm1v   =  etav*costh
   cterm1h   =  etah*costh
   epsnew   =  eps - sinth
   cterm2   =  csqrt(epsnew)

   cterm3v   =  (cterm1v - cterm2)/(cterm1v + cterm2)
   cterm3h   =  (cterm1h - cterm2)/(cterm1h + cterm2)
   tmp1r   =  real(cterm3v)
   tmp1i   = -aimag(cterm3v)
   ! ev   =  1.0 - (tmp1r*tmp1r+tmp1i*tmp1i)

   tmp2r   =  real(cterm3h)
   tmp2i   = -aimag(cterm3h)
   ! eh   =  1.0 - (tmp2r*tmp2r+tmp2i*tmp2i)

   ADJ_tmp2r   = - 2.0*tmp2r*ADJ_eh
   ADJ_tmp2i   = - 2.0*tmp2i*ADJ_eh

   ADJ_cterm3h =  ADJ_tmp2r + ADJ_tmp2i*(0.0,1.0)

   ADJ_tmp1r   = - 2.0*tmp1r*ADJ_ev
   ADJ_tmp1i   = - 2.0*tmp1i*ADJ_ev

   ADJ_cterm3v =  ADJ_tmp1r + ADJ_tmp1i*(0.0,1.0)

   ADJ_cterm2  = - ADJ_cterm3h/(cterm1h + cterm2)
   ADJ_cterm2  = - cterm3h*ADJ_cterm3h/(cterm1h + cterm2) + ADJ_cterm2

   ADJ_cterm1v =   ADJ_cterm3v/(cterm1v + cterm2)
   ADJ_cterm2  = - ADJ_cterm3v/(cterm1v + cterm2) + ADJ_cterm2
   ADJ_cterm1v = - cterm3v*ADJ_cterm3v/(cterm1v + cterm2) + ADJ_cterm1v
   ADJ_cterm2  = - cterm3v*ADJ_cterm3v/(cterm1v + cterm2) + ADJ_cterm2

   if (cabs(epsnew) .gt. 0.0) then
      ADJ_epsnew  = ADJ_cterm2*0.5/cterm2
   else
      ADJ_epsnew  =  0.0
   end if

   ADJ_eps     =  ADJ_epsnew

   ADJ_etav    =  ADJ_cterm1v*costh

   ADJ_eps     =  ADJ_etav + ADJ_eps

   ADJ_epsr    =  real(ADJ_eps)
   ADJ_epsi    =  -aimag(ADJ_eps) 
   ADJ_tc      =  0.0
   call da_epsalt_adj(f,tc,ssw,ADJ_tc, ADJ_epsr, ADJ_epsi)

   ADJ_tk      =  ADJ_tc + ADJ_tk

   if (trace_use) call da_trace_exit("da_spemiss_adj")

end subroutine da_spemiss_adj


subroutine da_tbatmos_adj(ifreq,theta,p0,wv,hwv,ta,gamma,lw,zcld,   &
   tbup,tbdn,tauatm, ADJ_theta,ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta,ADJ_gamma,   &
   ADJ_lw,ADJ_zcld,ADJ_tbup,ADJ_tbdn, ADJ_tauatm)

   implicit none

   !-----------------------------------------------------------------
   ! Purpose: TBD
   ! Output : ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta,ADJ_gamma,ADJ_lw,ADJ_zcld 
   !          ADJ_theta (somtime theta is a variable)
   ! Input  : ADJ_tbup,ADJ_tbdn,ADJ_tauatm
   ! Output mean fields : tbup,tbdn,tauatm
   !-----------------------------------------------------------------

   integer, intent(in)    :: ifreq
   real,    intent(in)    :: theta,p0,wv,hwv,ta,gamma,lw,zcld
   real,    intent(inout) :: ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta, ADJ_gamma,ADJ_lw,ADJ_zcld,ADJ_theta
   real,    intent(inout) :: ADJ_tbup,ADJ_tbdn,ADJ_tauatm
   real,    intent(out)   :: tbup,tbdn,tauatm

   real :: tbdn_save

   real :: mu,hdn,hup,hdninf,hupinf,ADJ_mu

   real :: b1(4),b2(4),b3(4)
   real :: c(4),d1(4),d2(4),d3(4),zeta(4),kw0(4),kw1(4),kw2(4),kw3(4)
   real :: tau,tau1,tau2,taucld
   real :: tcld,tc,em,em1
   real :: sigv,sigo,sig,sig1,sigcld
   real :: teff1dn,teff1up,teffdn,teffup
   real :: tbcld,tbclrdn,tbclrup,tb1dn,tb1up,tb2dn,tb2up
   real :: otbar,tc2,tc3,hv,ho,alph
   real :: ADJ_sigv,ADJ_otbar,ADJ_sigo,ADJ_tcld,ADJ_tc,ADJ_tc2,ADJ_tc3
   real :: ADJ_sigcld,ADJ_taucld,ADJ_tbcld,ADJ_hv,ADJ_ho
   real :: ADJ_hdn,ADJ_hup,ADJ_hdninf,ADJ_sig,ADJ_sig1,ADJ_tau,ADJ_tau1
   real :: ADJ_tau2,ADJ_em1,ADJ_teff1dn,ADJ_hupinf,ADJ_em,ADJ_teff1up
   real :: ADJ_teffdn,ADJ_teffup,ADJ_tbclrdn,ADJ_tbclrup,ADJ_tb1dn,ADJ_tb1up
   real :: ADJ_tb2dn,ADJ_tb2up,ADJ_alph

   data b1/-.46847e-1,-.57752e-1,-.18885,.10990/
   data b2/.26640e-4,.31662e-4,.9832e-4,.60531e-4/
   data b3/.87560e+1,.10961e+2,.36678e+2,-.37578e+2/
   data c/ .9207,   1.208,     .8253,     .8203/
   data zeta/4.2,4.2,4.2,2.9/
   data d1/-.35908e+1,-.38921e+1,-.43072e+1,-.17020e+0/
   data d2/ .29797e-1, .31054e-1, .32801e-1, .13610e-1/
   data d3/-.23174e-1,-.23543e-1,-.24101e-1,-.15776e+0/
   data kw0/ .786e-1, .103,    .267,    .988/
   data kw1/-.230e-2,-.296e-2,-.673e-2,-.107e-1/
   data kw2/ .448e-4, .557e-4, .975e-4,-.535e-4/
   data kw3/-.464e-6,-.558e-6,-.724e-6, .115e-5/

   if (trace_use) call da_trace_entry("da_tbatmos_adj")

   mu=0.0;hdn=0.0;hup=0.0;hdninf=0.0;hupinf=0.0;ADJ_mu=0.0
   tcld=0.0;tc=0.0;em=0.0;em1=0.0
   sigv=0.0;sigo=0.0;sig=0.0;sig1=0.0;sigcld=0.0
   teff1dn=0.0;teff1up=0.0;teffdn=0.0;teffup=0.0
   tbcld=0.0;tbclrdn=0.0;tbclrup=0.0;tb1dn=0.0;tb1up=0.0;tb2dn=0.0;tb2up=0.0
   otbar=0.0;tc2=0.0;tc3=0.0;hv=0.0;ho=0.0;alph=0.0
   ADJ_sigv=0.0;ADJ_otbar=0.0;ADJ_sigo=0.0;ADJ_tcld=0.0;
   ADJ_tc=0.0;ADJ_tc2=0.0;ADJ_tc3=0.0
   ADJ_sigcld=0.0;ADJ_taucld=0.0;ADJ_tbcld=0.0;ADJ_hv=0.0;ADJ_ho=0.0
   ADJ_hdn=0.0;ADJ_hup=0.0;ADJ_hdninf=0.0;ADJ_sig=0.0;ADJ_sig1=0.0
   ADJ_tau=0.0;ADJ_tau1=0.0
   ADJ_tau2=0.0;ADJ_em1=0.0;ADJ_teff1dn=0.0;ADJ_hupinf=0.0;ADJ_em=0.0
   ADJ_teff1up=0.0;ADJ_teffdn=0;ADJ_teffup=0.0;ADJ_tbclrdn=0.0
   ADJ_tbclrup=0.0;ADJ_tb1dn=0.0;ADJ_tb1up=0.0
   ADJ_tb2dn=0.0;ADJ_tb2up=0.0;ADJ_alph=0.0
   tau=0.0;tau1=0.0;tau2=0.0;taucld=0.0
   tcld=0.0;tc=0.0;em=0.0;em1=0.0
   sigv=0.0;sigo=0.0;sig=0.0;sig1=0.0;sigcld=0.0
   teff1dn=0.0;teff1up=0.0;teffdn=0.0;teffup=0.0


   ! mu = secant(theta)
   ! somtime theta is a variable

   mu     = 1.0/cos(theta*0.0174533)

   ! get water vapor optical depth

   call cal_sigma_v(ifreq,p0,wv,hwv,ta,gamma,sigv)

   ! otbar = one over "mean" temperature

   otbar =   1.0/(ta - gamma*zeta(ifreq))

   ! sigo = dry air optical depth

   sigo = b1(ifreq) + b2(ifreq)*    p0  + b3(ifreq)*    otbar

   ! cloud parameters

   tcld   =     ta - gamma*zcld
   tc     =     tcld - t_kelvin
   tc2    =     tc*tc
   tc3    =     tc2*tc
   sigcld =  ( kw0(ifreq) + tc*kw1(ifreq) + tc2*kw2(ifreq) +  &
               tc3*kw3(ifreq) )*lw
   taucld =   exp(-mu*sigcld)
   tbcld  =   (1.0 - taucld)*tcld

   ! hv, ho = effective absorber scale heights for vapor, dry air

   hv = c(ifreq)*   hwv
   ho = d1(ifreq) + d2(ifreq)* ta + d3(ifreq)* gamma

   ! get effective emission heights for layer 1 and total atmosphere


   call effht(ho,hv,sigo,sigv,mu,zcld,hdn,hup, hdninf,hupinf)

   ! atmospheric transmittances in layer one and two, and combined

   sig =     sigo +     sigv
   sig1 = sigo*(1.0-exp(-zcld/ho)) + sigv*(1.0-exp(-zcld/hv))
   tau  =  exp(-mu*sig)
   tau1 =  exp(-mu*sig1)
   tau2 =  tau/tau1
 
   ! atmospheric "emissivity"

   em1  =   1.0 - tau1
   em   =   1.0 - tau

   ! downwelling and upwelling brightness temperature for each layer

   teff1dn =     ta - gamma*hdn
   teff1up =     ta - gamma*hup
   teffdn =     ta - gamma*hdninf
   teffup =     ta - gamma*hupinf
   tbclrdn = teffdn*em
   tbclrup = teffup*em

   tb1dn = em1*teff1dn
   tb1up = em1*teff1up
   tb2dn = (tbclrdn - tb1dn)/tau1
   tb2up =      tbclrup - tau2*tb1up

   ! total downwelling and upwelling brightness temperature and transmittance

   tbdn  =     tb1dn + tau1*(tbcld + taucld*tb2dn)
   tbup  =     tb2up + tau2*(tbcld + taucld*tb1up)
   tauatm = tau*taucld

   ! the following lines apply an ad hoc correction to improve fit 
   ! at large angles and/or high gaseous opacities 
   !  (downwelling brightness temperatures only)

   alph = (0.636619*atan(mu*sig))**2
   tbdn_save = tbdn
   tbdn = (1.0-alph)*tbdn + em*alph*ta


   ! start

   tbdn = tbdn_save

   ADJ_alph = - ADJ_tbdn*tbdn
   ADJ_em   = ADJ_tbdn*alph*ta
   ADJ_alph = em*ADJ_tbdn*ta + ADJ_alph
   ADJ_ta   = em*alph*ADJ_tbdn + ADJ_ta
   ADJ_tbdn = (1.0-alph)*ADJ_tbdn 

   if (abs(sig) .gt. 0.0) then
      ADJ_mu  = 2.0*0.636619*0.636619*ADJ_alph*sig*atan(mu*sig)/(1.0+mu*mu*sig*sig)
      ADJ_sig = 2.0*0.636619*0.636619*mu*ADJ_alph*atan(mu*sig)/(1.0+mu*mu*sig*sig)
   else
      ADJ_mu  = 0.0
      ADJ_sig = 0.0
   end if

   ADJ_tau    = ADJ_tauatm*taucld
   ADJ_taucld = tau*ADJ_tauatm
   ADJ_tb2up  = ADJ_tbup
   ADJ_tau2   = ADJ_tbup*(tbcld + taucld*tb1up)
   ADJ_tbcld  = tau2*ADJ_tbup
   ADJ_taucld = tau2*ADJ_tbup*tb1up  + ADJ_taucld
   ADJ_tb1up  = tau2*taucld*ADJ_tbup
   ADJ_tb1dn = ADJ_tbdn
   ADJ_tau1  = ADJ_tbdn*(tbcld + taucld*tb2dn)
   ADJ_tbcld = tau1*ADJ_tbdn   + ADJ_tbcld
   ADJ_taucld = tau1*ADJ_tbdn*tb2dn  + ADJ_taucld
   ADJ_tb2dn  = tau1*taucld*ADJ_tbdn

   ADJ_tbclrup =   ADJ_tb2up
   ADJ_tau2    = - ADJ_tb2up*tb1up + ADJ_tau2
   ADJ_tb1up   = - tau2*ADJ_tb2up  + ADJ_tb1up

   ADJ_tbclrdn =   ADJ_tb2dn/tau1
   ADJ_tb1dn   = - ADJ_tb2dn/tau1  + ADJ_tb1dn
   ADJ_tau1    = - tb2dn*ADJ_tb2dn/tau1 + ADJ_tau1

   ADJ_em1     = ADJ_tb1up*teff1up
   ADJ_teff1up = em1*ADJ_tb1up

   ADJ_em1     = ADJ_tb1dn*teff1dn + ADJ_em1
   ADJ_teff1dn = em1*ADJ_tb1dn

   ADJ_teffup  = ADJ_tbclrup*em
   ADJ_em      = teffup*ADJ_tbclrup + ADJ_em

   ADJ_teffdn  = ADJ_tbclrdn*em
   ADJ_em      = teffdn*ADJ_tbclrdn + ADJ_em

   ADJ_ta      =   ADJ_teffup + ADJ_ta
   ADJ_gamma   = - ADJ_teffup*hupinf + ADJ_gamma
   ADJ_hupinf  = - gamma*ADJ_teffup

   ADJ_ta      =   ADJ_teffdn + ADJ_ta
   ADJ_gamma   = - ADJ_teffdn*hdninf + ADJ_gamma
   ADJ_hdninf  = - gamma*ADJ_teffdn

   ADJ_ta    =   ADJ_teff1up + ADJ_ta
   ADJ_gamma = - ADJ_teff1up*hup + ADJ_gamma
   ADJ_hup   = - gamma*ADJ_teff1up

   ADJ_ta    =   ADJ_teff1dn + ADJ_ta
   ADJ_gamma = - ADJ_teff1dn*hdn + ADJ_gamma
   ADJ_hdn   = - gamma*ADJ_teff1dn

   ADJ_tau  = - ADJ_em + ADJ_tau

   ADJ_tau1 = - ADJ_em1 + ADJ_tau1

   ADJ_tau  =   ADJ_tau2/tau1 + ADJ_tau
   ADJ_tau1 = - tau2*ADJ_tau2/tau1 + ADJ_tau1

   ADJ_sig1 = - mu*ADJ_tau1*tau1
   ADJ_mu   = - ADJ_tau1*sig1*tau1 + ADJ_mu
   
   ADJ_mu   = - ADJ_tau*sig*tau + ADJ_mu
   ADJ_sig  = - mu*ADJ_tau*tau + ADJ_sig
 
   ADJ_sigo = ADJ_sig1*(1.0-exp(-zcld/ho))
   ADJ_sigv = ADJ_sig1*(1.0-exp(-zcld/hv)) 
   ADJ_zcld = sigo*ADJ_sig1/ho*exp(-zcld/ho) + ADJ_zcld
   ADJ_ho   = - sigo*zcld*ADJ_sig1/(ho*ho)*exp(-zcld/ho)
   ADJ_zcld = sigv*ADJ_sig1/hv*exp(-zcld/hv) + ADJ_zcld
   ADJ_hv   = - sigv*zcld*ADJ_sig1/(hv*hv)*exp(-zcld/hv)

   ADJ_sigo = ADJ_sig + ADJ_sigo
   ADJ_sigv = ADJ_sig + ADJ_sigv

   call da_effht_adj(ho,hv,sigo,sigv,mu,zcld,hdn,hup,        &
                  hdninf,hupinf,                          &
                  ADJ_ho,ADJ_hv,ADJ_sigo,ADJ_sigv,ADJ_mu, &
                  ADJ_zcld,ADJ_hdn,ADJ_hup,ADJ_hdninf,    &
                  ADJ_hupinf                              )

   ADJ_ta    = d2(ifreq)*ADJ_ho + ADJ_ta
   ADJ_gamma = d3(ifreq)*ADJ_ho + ADJ_gamma
   ADJ_hwv   = c(ifreq)*ADJ_hv  + ADJ_hwv

   ADJ_taucld = - ADJ_tbcld*tcld + ADJ_taucld
   ADJ_tcld = (1.0 - taucld)*ADJ_tbcld

   ADJ_mu     = - ADJ_taucld*sigcld*taucld + ADJ_mu
   ADJ_sigcld = - mu*ADJ_taucld*taucld

   ADJ_tc  = ADJ_sigcld*kw1(ifreq)*lw
   ADJ_tc2 = ADJ_sigcld*kw2(ifreq)*lw
   ADJ_tc3 = ADJ_sigcld*kw3(ifreq)*lw
   ADJ_lw  = (kw0(ifreq)+tc*kw1(ifreq)+tc2*kw2(ifreq)+tc3*kw3(ifreq)) &
             *ADJ_sigcld + ADJ_lw

   ADJ_tc2  = ADJ_tc3*tc + ADJ_tc2
   ADJ_tc   = tc2*ADJ_tc3 + ADJ_tc
   ADJ_tc   = 2.0*tc*ADJ_tc2 + ADJ_tc

   ADJ_tcld = ADJ_tc + ADJ_tcld

   ADJ_ta    =   ADJ_tcld + ADJ_ta
   ADJ_gamma = - ADJ_tcld*zcld + ADJ_gamma
   ADJ_zcld  = - gamma*ADJ_tcld + ADJ_zcld

   ADJ_p0    = b2(ifreq)*ADJ_sigo + ADJ_p0
   ADJ_otbar = b3(ifreq)*ADJ_sigo

   ADJ_ta    = - otbar*otbar*ADJ_otbar + ADJ_ta
   ADJ_gamma =   otbar*otbar*ADJ_otbar*zeta(ifreq) + ADJ_gamma

   call da_sigma_v_adj(ifreq,p0,wv,hwv,ta,gamma,sigv,   &
                        ADJ_p0,ADJ_wv,ADJ_hwv,ADJ_ta,    &
                        ADJ_gamma,ADJ_sigv               )

   ADJ_theta = mu*mu*0.0174533*ADJ_mu*sin(theta*0.0174533) + ADJ_theta

   if (trace_use) call da_trace_exit("da_tbatmos_adj")

end subroutine da_tbatmos_adj



subroutine da_tb_tl(ifreq,theta,p0,ta,gamma,sst,wv,hwv,u,alw,zcld,            &
!  tbv,tbh,                                                  &
  TGL_p0,TGL_ta,TGL_gamma,TGL_sst,TGL_wv,                   &
  TGL_hwv,TGL_u,TGL_alw,TGL_zcld,TGL_tbv,TGL_tbh           )
   !---------------------------------------------------------------------------
   ! Purpose: TBD
   ! Input  : TGL_p0,  TGL_ta,   TGL_gamma, TGL_sst, TGL_wv, TGL_hwv, TGL_u
   !          TGL_alw, TGL_zcld
   ! Output : TGL_tbv, TGL_tbh,  tbv,  tbh
   ! ---------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: ifreq
   real,    intent(in) :: theta,p0,ta,gamma,sst,wv,hwv,u,alw,zcld
   real,    intent(in) :: TGL_p0,TGL_ta,TGL_gamma,TGL_sst,TGL_wv, TGL_hwv,TGL_u,TGL_alw,TGL_zcld
!   real   , intent(in) :: tbv,tbh
   real,    intent(out) :: TGL_tbv,TGL_tbh 

   real :: freq(4),ebiasv(4),ebiash(4),cf1(4),cf2(4),cg(4)

   real :: f,costhet,gx2,tbup,tbdn,tauatm,sigma,remv,remh
   real :: effangv,effangh,tbdnv,dum,foam,emissv,emissh 
   real :: refv,refh,semv,semh,scatv,scath,tbdnh
   real :: TGL_gx2,TGL_tbup,TGL_tbdn,TGL_tauatm,TGL_sigma
   real :: TGL_remh,TGL_effangv,TGL_effangh
   ! real :: TGL_tremv
   real :: TGL_tbdnh,TGL_dum,TGL_foam,TGL_emissv
   real :: TGL_emissh,TGL_refv,TGL_refh,TGL_semv
   real :: TGL_semh,TGL_scatv,TGL_scath,TGL_remv,TGL_tbdnv
   real :: TGL_theta

   real :: fem

   data freq/19.35,22.235,37.0,85.5/

   ! empirical bias corrections for surface emissivity

   data ebiasv/0.0095, 0.005,-0.014, -0.0010/
   data ebiash/0.004,   0.0,-0.023, 0.043/

   data cf1/0.0015,0.004,0.0027,0.005/
   data cg/4.50e-3, 5.200e-3, 5.5e-3, 7.2e-3 /

   data cf2/0.0023,0.000,0.0002,-0.006/

   ! 'foam' emissivity
   data fem /1.0/

   if (trace_use) call da_trace_entry("da_tb_tl")

   f = freq(ifreq)
   costhet = cos(theta*0.017453)

   ! effective surface slope variance

   gx2 = cg(ifreq)*    u
   TGL_gx2 = cg(ifreq)*TGL_u

   ! get upwelling atmospheric brightness temperature

   TGL_theta=0.0
   call da_tbatmos_tl(ifreq,theta,p0,wv,hwv,ta,gamma,alw,zcld, &
                       tbup,tbdn,tauatm,                        &
                       TGL_theta,TGL_p0,TGL_wv,TGL_hwv,TGL_ta,TGL_gamma,  &
                       TGL_alw,TGL_zcld,TGL_tbup,TGL_tbdn,      &
                       TGL_tauatm                              )

   ! convert transmittance to optical depth

   sigma = -alog(tauatm)*costhet
   TGL_sigma = -costhet*TGL_tauatm/tauatm

   ! get rough surface emissivity

   call da_roughem_tl(ifreq,gx2,sst,theta,remv,remh,         &
                       TGL_gx2,TGL_sst,TGL_remv,TGL_remh     )

   ! get effective zenith angles for scattered radiation at surface

   call da_effang_tl(ifreq,theta,gx2,sigma,effangv,effangh,    &
                      TGL_gx2,TGL_sigma,TGL_effangv,TGL_effangh)

   ! get effective sky brightness temperatures for scattered radiation

   call da_tbatmos_tl(ifreq,effangv,p0,wv,hwv,ta,gamma,alw,    &
                       zcld,dum,tbdnv,dum,                      &
                       TGL_effangv,TGL_p0,TGL_wv,TGL_hwv,       &
                       TGL_ta,TGL_gamma,TGL_alw,TGL_zcld,       & 
                       TGL_dum,TGL_tbdnv,TGL_dum               )

   call da_tbatmos_tl(ifreq,effangh,p0,wv,hwv,ta,gamma,alw,    &
                       zcld,dum,tbdnh,dum,                      &
                       TGL_effangh,TGL_p0,TGL_wv,TGL_hwv,       &
                       TGL_ta,TGL_gamma,TGL_alw,TGL_zcld,       &
                       TGL_dum,TGL_tbdnh,TGL_dum               )

   ! compute 'foam' coverage

   foam = cf1(ifreq)*    u
   TGL_foam = cf1(ifreq)*TGL_u

   if (u .gt. 5.0) then
      TGL_foam = TGL_foam + cf2(ifreq)*TGL_u
      foam =     foam + cf2(ifreq)*(  u-5.0)
   end if

   ! compute surface emissivities and reflectivity

   emissv     =     foam*fem + (1.0 - foam)*(remv + ebiasv(ifreq))
   TGL_emissv = TGL_foam*fem - TGL_foam*(remv + ebiasv(ifreq)) &
                                + (1.0 - foam)*TGL_remv
   emissh     =     foam*fem + (1.0 - foam)*(remh + ebiash(ifreq))
   TGL_emissh = TGL_foam*fem - TGL_foam*(remh + ebiash(ifreq)) &
                                + (1.0 - foam)*TGL_remh
   refv     =   1.0 - emissv
   TGL_refv = - TGL_emissv
   refh     =   1.0 - emissh
   TGL_refh = - TGL_emissh

   ! compute surface emission term

   semv     = sst*emissv
   TGL_semv = TGL_sst*emissv + sst*TGL_emissv
   semh     = sst*emissh
   TGL_semh = TGL_sst*emissh + sst*TGL_emissh

   ! compute surface scattering term

   scatv     = refv*tbdnv
   TGL_scatv = TGL_refv*tbdnv + refv*TGL_tbdnv
   scath     = refh*tbdnh
   TGL_scath = TGL_refh*tbdnh + refh*TGL_tbdnh

   ! combine to get space-observed brightness temperature

   ! tbv     =     tbup + tauatm*(semv + scatv)
   TGL_tbv = TGL_tbup + TGL_tauatm*(semv + scatv)   &
                         + tauatm*(TGL_semv + TGL_scatv)
   ! tbh     =     tbup + tauatm*(semh + scath)
   TGL_tbh = TGL_tbup + TGL_tauatm*(semh + scath)   &
                         + tauatm*(TGL_semh + TGL_scath)

   if (trace_use) call da_trace_exit("da_tb_tl")

end subroutine da_tb_tl


subroutine da_tbatmos_tl(ifreq,theta,p0,wv,hwv,ta,gamma,lw,zcld,   &
   tbup,tbdn,tauatm,                         &
   TGL_theta,TGL_p0,TGL_wv,TGL_hwv,TGL_ta,TGL_gamma,   &
   TGL_lw,TGL_zcld,TGL_tbup,TGL_tbdn,        &
   TGL_tauatm                               )

   !-----------------------------------------------------------------
   ! Purpose: TBD
   ! Input  : TGL_p0,TGL_wv,TGL_hwv,TGL_ta,TGL_gamma,TGL_lw,TGL_zcld 
   !          TGL_theta (somtime theta is a variable)
   ! Output : TGL_tbup,TGL_tbdn,TGL_tauatm,tbup,tbdn,tauatm
   ! -----------------------------------------------------------------

   implicit none

   integer,intent(in)  :: ifreq
   real   ,intent(in)  :: theta,p0,wv,hwv,ta,gamma,lw,zcld
   real   ,intent(in)  :: TGL_p0,TGL_wv,TGL_hwv,TGL_ta,TGL_gamma,TGL_lw,TGL_zcld,TGL_theta
   real   ,intent(out) :: tbup,tbdn,tauatm,TGL_tbup,TGL_tbdn, TGL_tauatm

   real :: mu,hdn,hup,hdninf,hupinf,TGL_mu

   real :: b1(4),b2(4),b3(4)
   real :: c(4),d1(4),d2(4),d3(4),zeta(4),kw0(4),kw1(4),kw2(4),kw3(4)
   real :: tau,tau1,tau2,taucld
   real :: tcld,tc,em,em1
   real :: sigv,sigo,sig,sig1,sigcld
   real :: teff1dn,teff1up,teffdn,teffup
   real :: tbcld,tbclrdn,tbclrup,tb1dn,tb1up,tb2dn,tb2up
   real :: otbar,tc2,tc3,hv,ho,alph
   real :: TGL_sigv,TGL_otbar,TGL_sigo,TGL_tcld,TGL_tc,TGL_tc2,TGL_tc3
   real :: TGL_sigcld,TGL_taucld,TGL_tbcld,TGL_hv,TGL_ho
   real :: TGL_hdn,TGL_hup,TGL_hdninf,TGL_sig,TGL_sig1,TGL_tau,TGL_tau1
   real :: TGL_tau2,TGL_em1,TGL_teff1dn,TGL_hupinf,TGL_em,TGL_teff1up
   real :: TGL_teffdn,TGL_teffup,TGL_tbclrdn,TGL_tbclrup,TGL_tb1dn,TGL_tb1up
   real :: TGL_tb2dn,TGL_tb2up,TGL_alph

   data b1/-.46847e-1,-.57752e-1,-.18885,.10990/
   data b2/.26640e-4,.31662e-4,.9832e-4,.60531e-4/
   data b3/.87560e+1,.10961e+2,.36678e+2,-.37578e+2/
   data c/ .9207,   1.208,     .8253,     .8203/
   data zeta/4.2,4.2,4.2,2.9/
   data d1/-.35908e+1,-.38921e+1,-.43072e+1,-.17020e+0/
   data d2/ .29797e-1, .31054e-1, .32801e-1, .13610e-1/
   data d3/-.23174e-1,-.23543e-1,-.24101e-1,-.15776e+0/
   data kw0/ .786e-1, .103,    .267,    .988/
   data kw1/-.230e-2,-.296e-2,-.673e-2,-.107e-1/
   data kw2/ .448e-4, .557e-4, .975e-4,-.535e-4/
   data kw3/-.464e-6,-.558e-6,-.724e-6, .115e-5/

   if (trace_use) call da_trace_entry("da_tbatmos_tl")

   ! mu = secant(theta)
   ! somtime theta is a variable

   mu     = 1.0/cos(theta*0.0174533)
   TGL_mu = mu*mu*0.0174533*TGL_theta*sin(theta*0.0174533)

   ! get water vapor optical depth

   call da_sigma_v_tl(ifreq,p0,wv,hwv,ta,gamma,sigv,   &
                        TGL_p0,TGL_wv,TGL_hwv,TGL_ta,    &
                        TGL_gamma,TGL_sigv              )

   ! otbar = one over "mean" temperature

       otbar =   1.0/(ta - gamma*zeta(ifreq))
   TGL_otbar = - otbar*otbar*(TGL_ta - TGL_gamma*zeta(ifreq))

   ! sigo = dry air optical depth

       sigo = b1(ifreq) + b2(ifreq)*    p0  + b3(ifreq)*    otbar
   TGL_sigo =             b2(ifreq)*TGL_p0  + b3(ifreq)*TGL_otbar

   ! cloud parameters

       tcld   =     ta - gamma*zcld
   TGL_tcld   = TGL_ta - TGL_gamma*zcld - gamma*TGL_zcld
         tc   =     tcld - t_kelvin
     TGL_tc   = TGL_tcld
        tc2   = tc*tc
    TGL_tc2   = 2.0*tc*TGL_tc
        tc3   = tc2*tc
    TGL_tc3   = TGL_tc2*tc + tc2*TGL_tc
       sigcld =  (kw0(ifreq) + tc*kw1(ifreq) + tc2*kw2(ifreq) +  &
                   tc3*kw3(ifreq))*lw
   TGL_sigcld =  (TGL_tc *kw1(ifreq) + TGL_tc2*kw2(ifreq) +  &
                   TGL_tc3*kw3(ifreq))*lw &
               + (kw0(ifreq) + tc*kw1(ifreq) + tc2*kw2(ifreq) +  &
                   tc3*kw3(ifreq))*TGL_lw

       taucld =   exp(-mu*sigcld)
   TGL_taucld = - (TGL_mu*sigcld + mu*TGL_sigcld)*taucld
        tbcld =   (1.0 - taucld)*tcld
    TGL_tbcld = - TGL_taucld*tcld + (1.0 - taucld)*TGL_tcld

   ! hv, ho = effective absorber scale heights for vapor, dry air

       hv = c(ifreq)*    hwv
   TGL_hv = c(ifreq)*TGL_hwv
       ho = d1(ifreq) + d2(ifreq)*    ta + d3(ifreq)*    gamma
   TGL_ho =             d2(ifreq)*TGL_ta + d3(ifreq)*TGL_gamma

   ! get effective emission heights for layer 1 and total atmosphere

   call da_effht_tl(ho,hv,sigo,sigv,mu,zcld,hdn,hup,        &
                  hdninf,hupinf,                          &
                  TGL_ho,TGL_hv,TGL_sigo,TGL_sigv,TGL_mu, &
                  TGL_zcld,TGL_hdn,TGL_hup,TGL_hdninf,    &
                  TGL_hupinf                             )

   ! atmospheric transmittances in layer one and two, and combined

        sig =     sigo +     sigv

    TGL_sig = TGL_sigo + TGL_sigv
       sig1 = sigo*(1.0-exp(-zcld/ho)) + sigv*(1.0-exp(-zcld/hv))
   TGL_sig1 =   TGL_sigo*(1.0-exp(-zcld/ho))  &
              + TGL_sigv*(1.0-exp(-zcld/hv))  &
              + sigo*(TGL_zcld/ho - zcld*TGL_ho/(ho*ho))*exp(-zcld/ho) &
              + sigv*(TGL_zcld/hv - zcld*TGL_hv/(hv*hv))*exp(-zcld/hv)
       tau  =  exp(-mu*sig)
   TGL_tau  = -(TGL_mu*sig + mu*TGL_sig) * tau
       tau1 =  exp(-mu*sig1)
   TGL_tau1 = -(mu*TGL_sig1 + TGL_mu*sig1) *tau1
       tau2 =  tau/tau1
   TGL_tau2 =  TGL_tau/tau1 - tau2*TGL_tau1/tau1

   ! atmospheric "emissivity"

       em1  =   1.0 - tau1
   TGL_em1  = - TGL_tau1
       em   =   1.0 - tau
   TGL_em   = - TGL_tau

   ! downwelling and upwelling brightness temperature for each layer

       teff1dn =     ta - gamma*hdn
   TGL_teff1dn = TGL_ta - TGL_gamma*hdn - gamma*TGL_hdn
       teff1up =     ta - gamma*hup
   TGL_teff1up = TGL_ta - TGL_gamma*hup - gamma*TGL_hup
        teffdn =     ta - gamma*hdninf
    TGL_teffdn = TGL_ta - TGL_gamma*hdninf - gamma*TGL_hdninf
        teffup =     ta - gamma*hupinf
    TGL_teffup = TGL_ta - TGL_gamma*hupinf - gamma*TGL_hupinf
       tbclrdn = teffdn*em
   TGL_tbclrdn = TGL_teffdn*em + teffdn*TGL_em
       tbclrup = teffup*em
   TGL_tbclrup = TGL_teffup*em + teffup*TGL_em

        tb1dn = em1*teff1dn
   TGL_tb1dn = TGL_em1*teff1dn + em1*TGL_teff1dn
        tb1up = em1*teff1up
   TGL_tb1up = TGL_em1*teff1up + em1*TGL_teff1up
        tb2dn = (tbclrdn - tb1dn)/tau1
   TGL_tb2dn = (TGL_tbclrdn - TGL_tb1dn)/tau1 - tb2dn*TGL_tau1/tau1
        tb2up =      tbclrup - tau2*tb1up
   TGL_tb2up =  TGL_tbclrup - TGL_tau2*tb1up - tau2*TGL_tb1up

   ! total downwelling and upwelling brightness temperature and transmittance

       tbdn  =     tb1dn + tau1*(tbcld + taucld*tb2dn)
   TGL_tbdn  = TGL_tb1dn + TGL_tau1*(tbcld + taucld*tb2dn) &
                         + tau1*(TGL_tbcld + TGL_taucld*tb2dn + taucld*TGL_tb2dn)
       tbup  =     tb2up + tau2*(tbcld + taucld*tb1up)
   TGL_tbup  = TGL_tb2up + TGL_tau2*(tbcld + taucld*tb1up) &
                         + tau2*(TGL_tbcld + TGL_taucld*tb1up + taucld*TGL_tb1up)
       tauatm = tau*taucld
   TGL_tauatm = TGL_tau*taucld + tau*TGL_taucld
   !
   ! the following lines apply an ad hoc correction to improve fit 
   ! at large angles and/or high gaseous opacities 
   !  (downwelling brightness temperatures only)

   alph = (0.636619*atan(mu*sig))**2
   if (abs(sig) .gt. 0.0) then
      TGL_alph = 2.0*0.636619*0.636619* &
                 (TGL_mu*sig + mu*TGL_sig)*atan(mu*sig)/(1.0+mu*mu*sig*sig)
   else
      TGL_alph = 0.0
   end if
   TGL_tbdn = - TGL_alph*tbdn + (1.0-alph)*TGL_tbdn &
              + TGL_em*alph*ta + em*TGL_alph*ta + em*alph*TGL_ta 
   tbdn = (1.0-alph)*tbdn + em*alph*ta

   if (trace_use) call da_trace_exit("da_tbatmos_tl")

end subroutine da_tbatmos_tl


subroutine da_effht_tl(ho,hv,sigo,sigv,mu,zcld,hdn,hup,hdninf,hupinf, &
                     TGL_ho,TGL_hv,TGL_sigo,TGL_sigv,TGL_mu,        &
                     TGL_zcld,TGL_hdn,TGL_hup,TGL_hdninf,TGL_hupinf)

   !--------------------------------------------------------------------
   ! Purpose: TBD
   ! Input  : TGL_ho, TGL_hv, TGL_sigo, TGL_sigv, TGL_mu, TGL_zcld
   ! Output : TGL_hdn, hdn, TGL_hup, hup, 
   !         TGL_hdninf, hdninf, TGL_hupinf, hupinf
   !--------------------------------------------------------------------

   implicit none

   real,   intent(in)  :: ho,hv,sigo,sigv,mu,zcld
   real,   intent(in)  :: TGL_ho,TGL_hv,TGL_sigo,TGL_sigv,TGL_zcld, TGL_mu
   real,   intent(out) :: hdn,hup,hdninf,hupinf
   real,   intent(out) :: TGL_hdn,TGL_hup,TGL_hdninf,TGL_hupinf

   real :: gint,zgint,hint,zhint
   real :: ginf,zginf,hinf,zhinf
   real :: TGL_gint,TGL_zgint,TGL_hint,TGL_zhint
   real :: TGL_ginf,TGL_zginf,TGL_hinf,TGL_zhinf
   real :: TGL_mu2,TGL_halfmu,TGL_sixthmu2,TGL_etnthmu2
   real :: TGL_quartmu,TGL_halfmu2

   real :: hoinv,hvinv,chio,chiv,ezho,ezhv,alpha,alph2,alph3
   real :: beta,beta2,beta3,mu2,mualph,cplus,cmin,dplus,dmin
   real :: chiov,chivv,chioo,chioov,chiovv,chiooo,chivvv
   real :: h11,h21,h12,newh11
   real :: sigoo,sigov,sigvv,sigooo,sigoov,sigovv,sigvvv
   real :: ezhoo,ezhov,ezhvv,ezhooo,ezhoov,ezhovv,ezhvvv
   real :: s,sprim,t,tprim,u,uprim,term1,term2,term3
   real :: halfmu,halfmu2,sixthmu2,etnthmu2,quartmu

   real :: TGL_hoinv,TGL_hvinv,TGL_chio,TGL_chiv,TGL_ezho
   real :: TGL_ezhv,TGL_alpha,TGL_alph2,TGL_alph3
   real :: TGL_beta,TGL_beta2,TGL_beta3,TGL_mualph
   real :: TGL_cplus,TGL_cmin,TGL_dplus,TGL_dmin
   real :: TGL_chiov,TGL_chivv,TGL_chioo,TGL_chioov
   real :: TGL_chiovv,TGL_chiooo,TGL_chivvv
   real :: TGL_h11,TGL_h21,TGL_h12,TGL_newh11
   real :: TGL_sigoo,TGL_sigov,TGL_sigvv,TGL_sigooo
   real :: TGL_sigoov,TGL_sigovv,TGL_sigvvv
   real :: TGL_ezhoo,TGL_ezhov,TGL_ezhvv,TGL_ezhooo
   real :: TGL_ezhoov,TGL_ezhovv,TGL_ezhvvv
   real :: TGL_s,TGL_sprim,TGL_t,TGL_tprim
   real :: TGL_u,TGL_uprim,TGL_term1,TGL_term2,TGL_term3

   if (trace_use) call da_trace_entry("da_effht_tl")

       hoinv =  1.0d0/ho
   TGL_hoinv = -1.0d0*hoinv*hoinv*TGL_ho

       hvinv =  1.0d0/hv
   TGL_hvinv = -1.0d0*hvinv*hvinv*TGL_hv

        chio = zcld*hoinv
    TGL_chio = TGL_zcld*hoinv + zcld*TGL_hoinv

           chiv = zcld*hvinv
    TGL_chiv = TGL_zcld*hvinv + zcld*TGL_hvinv

        ezho = sigo*exp(-chio)
    TGL_ezho = TGL_sigo*exp(-chio)-TGL_chio*ezho

        ezhv = sigv*exp(-chiv)
    TGL_ezhv = TGL_sigv*exp(-chiv)-TGL_chiv*ezhv

       alpha = sigo + sigv
   TGL_alpha = TGL_sigo + TGL_sigv

       alph2 = alpha*alpha
   TGL_alph2 = 2.0*alpha*TGL_alpha

       alph3 = alpha*alph2
   TGL_alph3 = TGL_alpha*alph2+alpha*TGL_alph2

        beta = ezho + ezhv
    TGL_beta = TGL_ezho + TGL_ezhv

       beta2 = beta*beta
   TGL_beta2 = 2.0*beta*TGL_beta

       beta3 = beta*beta2
   TGL_beta3 = TGL_beta*beta2+beta*TGL_beta2

       mu2        = mu*mu
   TGL_mu2        = 2.0*mu*TGL_mu
       halfmu     = 0.5d0*    mu
   TGL_halfmu     = 0.5d0*TGL_mu
       sixthmu2   =     mu2/6.0d0
   TGL_sixthmu2   = TGL_mu2/6.0d0
       etnthmu2   =     mu2/18.0d0
   TGL_etnthmu2   = TGL_mu2/18.0d0
       quartmu    = 0.25d0*    mu
   TGL_quartmu    = 0.25d0*TGL_mu
       halfmu2    = 0.5d0*    mu2
   TGL_halfmu2    = 0.5d0*TGL_mu2

          mualph = mu*alpha
   TGL_mualph = TGL_mu*alpha + mu*TGL_alpha

       cplus  = 1.0d0 +     mualph
   TGL_cplus  =         TGL_mualph

       cmin   = 1.0d0 -     mualph
   TGL_cmin   =       - TGL_mualph

       dplus  = halfmu2*alph2
   TGL_dplus  = TGL_halfmu2*alph2 + halfmu2*TGL_alph2

       dmin   =     dplus
   TGL_dmin   = TGL_dplus

   TGL_dplus  = TGL_cplus + TGL_dplus
       dplus  =     cplus +     dplus

   TGL_dmin   = TGL_cmin  + TGL_dmin
       dmin   =     cmin  +     dmin


       h11    =     hoinv +     hvinv
   TGL_h11    = TGL_hoinv + TGL_hvinv

       h21    =  1.0d0/(h11 + hvinv)
   TGL_h21    = -1.0d0*h21*h21*(TGL_h11+TGL_hvinv)

       h12    =  1.0d0/(h11 + hoinv)
   TGL_h12    = -1.0d0*h12*h12*(TGL_h11 + TGL_hoinv)

       newh11 =  1.0d0/h11
   TGL_newh11 = -1.0d0*newh11*newh11*TGL_h11

       chiov  = 1.0d0 +     chio +     chiv
   TGL_chiov  =         TGL_chio + TGL_chiv

       chioo  = 1.0d0 +     chio +     chio
   TGL_chioo  =         TGL_chio + TGL_chio

       chivv  = 1.0d0 +     chiv +     chiv
   TGL_chivv  =         TGL_chiv + TGL_chiv

       chioov =     chioo +     chiv
   TGL_chioov = TGL_chioo + TGL_chiv

          chiovv =     chio  +     chivv
   TGL_chiovv = TGL_chio  + TGL_chivv

       chiooo =     chioo +     chio
   TGL_chiooo = TGL_chioo + TGL_chio

       chivvv =     chivv +     chiv
   TGL_chivvv = TGL_chivv + TGL_chiv

   TGL_chio   =         TGL_chio
       chio   = 1.0d0 +     chio

   TGL_chiv   =         TGL_chiv
       chiv   = 1.0d0 +     chiv

       sigov  = sigo*sigv
   TGL_sigov  = TGL_sigo*sigv + sigo*TGL_sigv

       sigoo  = sigo*sigo
   TGL_sigoo  = 2.0*sigo*TGL_sigo

       sigvv  = sigv*sigv
   TGL_sigvv  = 2.0*sigv*TGL_sigv

       sigooo = sigoo*sigo
   TGL_sigooo = TGL_sigoo*sigo + sigoo*TGL_sigo

       sigoov = sigoo*sigv
   TGL_sigoov = TGL_sigoo*sigv + sigoo*TGL_sigv

       sigovv = sigo*sigvv
   TGL_sigovv = TGL_sigo*sigvv + sigo*TGL_sigvv

       sigvvv = sigvv*sigv
   TGL_sigvvv = TGL_sigvv*sigv + sigvv*TGL_sigv

       ezhoo  = ezho*ezho
   TGL_ezhoo  = 2.0*ezho*TGL_ezho

       ezhov  = ezho*ezhv
   TGL_ezhov  = TGL_ezho*ezhv + ezho*TGL_ezhv

       ezhvv  = ezhv*ezhv
   TGL_ezhvv  = 2.0*ezhv*TGL_ezhv

       ezhovv = ezho*ezhvv
   TGL_ezhovv = TGL_ezho*ezhvv + ezho*TGL_ezhvv

       ezhoov = ezhoo*ezhv
   TGL_ezhoov = TGL_ezhoo*ezhv + ezhoo*TGL_ezhv

       ezhooo = ezhoo*ezho
   TGL_ezhooo = TGL_ezhoo*ezho + ezhoo*TGL_ezho

       ezhvvv = ezhvv*ezhv
   TGL_ezhvvv = TGL_ezhvv*ezhv + ezhvv*TGL_ezhv

       s      = sigo*ho + sigv*hv
   TGL_s      = TGL_sigo*ho + sigo*TGL_ho + TGL_sigv*hv + sigv*TGL_hv

       sprim  = ezho*ho*chio + ezhv*hv*chiv
   TGL_sprim  = TGL_ezho*ho*chio + ezho*TGL_ho*chio + ezho*ho*TGL_chio + &
                TGL_ezhv*hv*chiv + ezhv*TGL_hv*chiv + ezhv*hv*TGL_chiv

       t      = sigoo*ho + 4.0d0*sigov*newh11 + sigvv*hv
   TGL_t      = TGL_sigoo*ho + sigoo*TGL_ho + &
                4.0d0*(TGL_sigov*newh11 + sigov*TGL_newh11) + &
                TGL_sigvv*hv + sigvv*TGL_hv

       tprim  = ezhoo*ho*chioo + 4.0d0*ezhov*newh11*chiov + ezhvv*hv*chivv
   TGL_tprim  = TGL_ezhoo*ho*chioo +ezhoo*TGL_ho*chioo + ezhoo*ho*TGL_chioo + &
                4.0d0*(TGL_ezhov*newh11*chiov +    &
                       ezhov*TGL_newh11*chiov +    &
                       ezhov*newh11*TGL_chiov ) + &
                TGL_ezhvv*hv*chivv + ezhvv*TGL_hv*chivv + ezhvv*hv*TGL_chivv

       u      = sigooo*ho + 9.0d0*(sigovv*h21+sigoov*h12) + sigvvv*hv
   TGL_u      = TGL_sigooo*ho + sigooo*TGL_ho + &
                9.0d0*(TGL_sigovv*h21 + sigovv*TGL_h21 +    &
                       TGL_sigoov*h12 + sigoov*TGL_h12 ) + &
                TGL_sigvvv*hv + sigvvv*TGL_hv

       uprim  = ezhvvv*hv*chivvv +  &
                9.0d0*(ezhovv*h21*chiovv + ezhoov*h12*chioov) + &
                ezhooo*ho*chiooo
   TGL_uprim  = TGL_ezhvvv*hv*chivvv +ezhvvv*TGL_hv*chivvv+ ezhvvv*hv*TGL_chivvv+  &
                 9.0d0*(TGL_ezhovv*h21*chiovv +     &
                        ezhovv*TGL_h21*chiovv +     &
                        ezhovv*h21*TGL_chiovv +     &
                        TGL_ezhoov*h12*chioov +     &
                        ezhoov*TGL_h12*chioov +     &
                        ezhoov*h12*TGL_chioov  ) + &
                 TGL_ezhooo*ho*chiooo + ezhooo*TGL_ho*chiooo + ezhooo*ho*TGL_chiooo

       term1  =     s -     sprim
   TGL_term1  = TGL_s - TGL_sprim

       term2  = quartmu*(t - tprim)
   TGL_term2  = TGL_quartmu*(t - tprim) + quartmu*(TGL_t - TGL_tprim) 

       term3  = etnthmu2*(   u -     uprim)
   TGL_term3  = TGL_etnthmu2*(u - uprim) + etnthmu2*(TGL_u - TGL_uprim)

       zgint  = dmin*term1 +  cmin*term2 + term3
   TGL_zgint  = TGL_dmin*term1 + dmin*TGL_term1 + &
                TGL_cmin*term2 + cmin*TGL_term2 + TGL_term3

    zhint  = -dplus*term1 + cplus*term2 - term3
   TGL_zhint  = -TGL_dplus*term1 - dplus*TGL_term1 + &
                 TGL_cplus*term2 + cplus*TGL_term2 - TGL_term3

       term2  = quartmu * t
   TGL_term2  = TGL_quartmu*t + quartmu*TGL_t

       term3  = etnthmu2*u
   TGL_term3  = TGL_etnthmu2*u + etnthmu2*TGL_u

       zginf  = dmin*s +  cmin*term2 + term3
   TGL_zginf  = TGL_dmin*s + dmin*TGL_s +  &
                TGL_cmin*term2 + cmin*TGL_term2 + TGL_term3

       zhinf  = -dplus*s + cplus*term2 - term3
   TGL_zhinf  = -TGL_dplus*s - dplus*TGL_s + &
                 TGL_cplus*term2 + cplus*TGL_term2 - TGL_term3

       term1  =     alpha -     beta
   TGL_term1  = TGL_alpha - TGL_beta

       term2  = halfmu*(   alph2 -     beta2)
   TGL_term2  = TGL_halfmu*(alph2 - beta2) + halfmu*(TGL_alph2 - TGL_beta2)

       term3  = sixthmu2*(   alph3 -     beta3)
   TGL_term3  = TGL_sixthmu2*(alph3 - beta3) + sixthmu2*(TGL_alph3 - TGL_beta3)

       gint   = dmin*term1 +  cmin*term2 + term3
   TGL_gint   = TGL_dmin*term1 + dmin*TGL_term1 + &
                TGL_cmin*term2 + cmin*TGL_term2 + TGL_term3

       hint   = -dplus*term1 + cplus*term2 - term3
   TGL_hint   = -TGL_dplus*term1 - dplus*TGL_term1 + &
                 TGL_cplus*term2 + cplus*TGL_term2 - TGL_term3

       term2  = halfmu*alph2
   TGL_term2  = TGL_halfmu*alph2 + halfmu*TGL_alph2

       term3  = sixthmu2*alph3
   TGL_term3  = TGL_sixthmu2*alph3 + sixthmu2*TGL_alph3

       ginf   = dmin*alpha +  cmin*term2 + term3
   TGL_ginf   = TGL_dmin*alpha + dmin*TGL_alpha +  &
                TGL_cmin*term2 + cmin*TGL_term2 + TGL_term3

       hinf   = -dplus*alpha + cplus*term2 - term3
   TGL_hinf   = -TGL_dplus*alpha - dplus*TGL_alpha + &
                 TGL_cplus*term2 + cplus*TGL_term2 - TGL_term3

       hdn    = zgint/gint
   TGL_hdn    = TGL_zgint/gint - hdn * TGL_gint/gint

       hup    = zhint/hint
   TGL_hup    = TGL_zhint/hint - hup*TGL_hint/hint

       hdninf = zginf/ginf
   TGL_hdninf = TGL_zginf/ginf - hdninf*TGL_ginf/ginf

       hupinf = zhinf/hinf
   TGL_hupinf = TGL_zhinf/hinf - hupinf*TGL_hinf/hinf

   if (trace_use) call da_trace_exit("da_effht_tl")

end subroutine da_effht_tl


subroutine da_roughem_tl(ifreq,gx2,tk,theta,remv,remh, TGL_gx2,TGL_tk,TGL_remv,TGL_remh         )

   !----------------------------------------------------------------
   ! Purpose: Calculates rough-surface emissivity of ocean surface at SSM/I
   ! frequencies.
   ! Input  : TGL_tk, TGL_gx2
   ! Output : TGL_remv,TGL_remh, remv,remh
   !----------------------------------------------------------------

   implicit none

   integer, intent(in)  :: ifreq
   real,    intent(in)  :: tk, theta, gx2, TGL_tk,TGL_gx2
   real,    intent(out) :: TGL_remv,TGL_remh, remv,remh

   real :: a19v(4),a22v(4),a37v(4),a85v(4)
   real :: a19h(4),a22h(4),a37h(4),a85h(4)
   real :: f(4)
   real :: semv,semh,TGL_semv,TGL_semh,ssw
   real :: tp,g,x1,x2,x3,x4,dtheta
   real :: TGL_tp,TGL_g,TGL_x1,TGL_x2,TGL_x3,TGL_x4

   data a19v/  -0.111E+01,   0.713E+00,  -0.624E-01,   0.212E-01 /
   data a19h/   0.812E+00,  -0.215E+00,   0.255E-01,   0.305E-02 /
   data a22v/  -0.134E+01,   0.911E+00,  -0.893E-01,   0.463E-01 /
   data a22h/   0.958E+00,  -0.350E+00,   0.566E-01,  -0.262E-01 /
   data a37v/  -0.162E+01,   0.110E+01,  -0.730E-01,   0.298E-01 /
   data a37h/   0.947E+00,  -0.320E+00,   0.624E-01,  -0.300E-01 /
   data a85v/  -0.145E+01,   0.808E+00,  -0.147E-01,  -0.252E-01 /
   data a85h/   0.717E+00,  -0.702E-01,   0.617E-01,  -0.243E-01 /

   data f/ 19.35, 22.235, 37.0, 85.5 /

   if (trace_use) call da_trace_entry("da_roughem_tl")

   tp     = tk/t_roughem
   TGL_tp = TGL_tk/t_roughem
   dtheta = theta-53.0
   g      = 0.5*    gx2
   TGL_g  = 0.5*TGL_gx2
   x1     = g
   TGL_x1 = TGL_g
   x2     = tp*g
   TGL_x2 = TGL_tp*g + tp*TGL_g
   x3     = dtheta*    g
   TGL_x3 = dtheta*TGL_g
   x4     = tp*x3
   TGL_x4 = TGL_tp*x3 + tp*TGL_x3

   if (ifreq == 1) then
      remv     =     x1*a19v(1) +     x2*a19v(2) +     x3*a19v(3) +     x4*a19v(4)
      TGL_remv = TGL_x1*a19v(1) + TGL_x2*a19v(2) + TGL_x3*a19v(3) + TGL_x4*a19v(4)
      remh     =     x1*a19h(1) +     x2*a19h(2) +     x3*a19h(3) +     x4*a19h(4)
      TGL_remh = TGL_x1*a19h(1) + TGL_x2*a19h(2) + TGL_x3*a19h(3) + TGL_x4*a19h(4)
   else if (ifreq == 2) then
      remv     =     x1*a22v(1) +     x2*a22v(2) +     x3*a22v(3) +     x4*a22v(4)
      TGL_remv = TGL_x1*a22v(1) + TGL_x2*a22v(2) + TGL_x3*a22v(3) + TGL_x4*a22v(4)
      remh     =     x1*a22h(1) +     x2*a22h(2) +     x3*a22h(3) +     x4*a22h(4)
      TGL_remh = TGL_x1*a22h(1) + TGL_x2*a22h(2) + TGL_x3*a22h(3) + TGL_x4*a22h(4)
   else if (ifreq == 3) then
      remv     =     x1*a37v(1) +     x2*a37v(2) +     x3*a37v(3) +     x4*a37v(4)
      TGL_remv = TGL_x1*a37v(1) + TGL_x2*a37v(2) + TGL_x3*a37v(3) + TGL_x4*a37v(4)
      remh     =     x1*a37h(1) +     x2*a37h(2) +     x3*a37h(3) +     x4*a37h(4)
      TGL_remh = TGL_x1*a37h(1) + TGL_x2*a37h(2) + TGL_x3*a37h(3) + TGL_x4*a37h(4)
   else if (ifreq == 4) then
      remv     =     x1*a85v(1) +     x2*a85v(2) +     x3*a85v(3) +     x4*a85v(4)
      TGL_remv = TGL_x1*a85v(1) + TGL_x2*a85v(2) + TGL_x3*a85v(3) + TGL_x4*a85v(4)
      remh     =     x1*a85h(1) +     x2*a85h(2) +     x3*a85h(3) +     x4*a85h(4)
      TGL_remh = TGL_x1*a85h(1) + TGL_x2*a85h(2) + TGL_x3*a85h(3) + TGL_x4*a85h(4)
   end if

   ssw=36.5
   call da_spemiss_tl(f(ifreq),tk,theta,ssw,semv,semh, TGL_tk,TGL_semv,TGL_semh)

   TGL_remv = TGL_remv + TGL_semv
   remv     = remv + semv
   TGL_remh = TGL_remh + TGL_semh
   remh     = remh + semh

   if (trace_use) call da_trace_exit("da_roughem_tl")
   
end subroutine da_roughem_tl


subroutine da_spemiss_tl(f,tk,theta,ssw,ev,eh,TGL_tk,TGL_ev,TGL_eh                               )

   !-----------------------------------------------------------------------
   ! Purpose: returns the specular emissivity of sea water for given 
   ! freq. (GHz), temperature T (K), incidence angle theta (degrees), 
   ! salinity (permil)
   !     
   ! Returned values verified against data in Klein and Swift (1977) and
   ! against Table 3.8 in Olson (1987, Ph.D. Thesis)
   !
   ! Input  : TGL_tk
   ! Output : TGL_ev, TGL_eh, ev, eh
   !------------------------------------------------------------------------

   implicit none

   real, intent(in)    :: f, tk, theta, TGL_tk
   real, intent(inout) :: ssw
   real, intent(out)   :: TGL_ev, TGL_eh, ev, eh

   real :: epsr,epsi,TGL_epsr,TGL_epsi

   real    ::  tc,costh,sinth,rthet
   complex ::  etav,etah,eps,cterm1v,cterm1h,cterm2,cterm3v,cterm3h,epsnew
   complex ::  TGL_etav,TGL_eps,TGL_cterm1v,TGL_cterm2,TGL_cterm3v
   complex ::  TGL_cterm3h,TGL_epsnew
   ! complex   uniti
   real    ::  tmp1r,tmp1i,tmp2r,tmp2i
   ! real :: rnorm,tmp0i,tmp0r
   real    ::  TGL_tc,TGL_tmp1r
   ! real :: TGL_tmp0r,TGL_tmp0i,TGL_rnorm
   real    ::  TGL_tmp1i,TGL_tmp2r,TGL_tmp2i 

   if (trace_use) call da_trace_entry("da_spemiss_tl")


   tc          =      tk - t_kelvin
   TGL_tc      =  TGL_tk

   call da_epsalt_tl(f,tc,ssw,epsr,epsi,TGL_tc, TGL_epsr, TGL_epsi )

       eps     =  cmplx(epsr,epsi)
   TGL_eps     =  cmplx(TGL_epsr,TGL_epsi)
       etav    =  eps
   TGL_etav    =  TGL_eps
   etah        =  (1.0,0.0)
   rthet       =  theta*0.017453292
   costh       =  cos(rthet)
   sinth       =  sin(rthet)
   sinth       =  sinth*sinth
       cterm1v =  etav*costh
   TGL_cterm1v =  TGL_etav*costh
   cterm1h     =  etah*costh
       epsnew  =      eps - sinth
   TGL_epsnew  =  TGL_eps
   cterm2      =  csqrt(epsnew)

   ! calculate TGL_cterm2

   if (cabs(epsnew) .gt. 0.0) then
      TGL_cterm2      =  TGL_epsnew*0.5/cterm2
   else
      TGL_cterm2      =  0.0
   end if

   ! Wei's Comment
   !     It is not a standard fortran if statement here.

   !     if (0) then
   !               tmp0r   =  real(epsnew)
   !           TGL_tmp0r   =  real(TGL_epsnew)
   !               tmp0i   = -aimag(epsnew)
   !           TGL_tmp0i   = -aimag(TGL_epsnew)
   !               rnorm   =  sqrt(tmp0r*tmp0r+tmp0i*tmp0i)
   !               uniti   =  (0,1)
   !           if (rnorm .gt. 0.0) then
   !             if (abs(tmp0i) .gt. 0.0) then
   !                TGL_rnorm =  (tmp0r*TGL_tmp0r + tmp0i*TGL_tmp0i)/rnorm
   !                TGL_cterm2=  cterm2*0.5*(TGL_rnorm/rnorm  &
   !                                      -uniti*(TGL_tmp0r*rnorm-TGL_rnorm*tmp0r)/(rnorm*tmp0i))
   !             else
   !                TGL_rnorm =  TGL_tmp0r
   !                TGL_cterm2=  TGL_tmp0r*0.5/sqrt(tmp0r)
   !             end if
   !           else 
   !             TGL_rnorm =  0.0
   !             TGL_cterm2=  0.0
   !           end if
   !     end if

   ! End Wei's Comment

       cterm3v =  (cterm1v - cterm2)/(cterm1v + cterm2)
   TGL_cterm3v =  (TGL_cterm1v - TGL_cterm2)/(cterm1v + cterm2) &
                 -cterm3v*(TGL_cterm1v + TGL_cterm2)/(cterm1v + cterm2)
       cterm3h =  (cterm1h - cterm2)/(cterm1h + cterm2)
   TGL_cterm3h = -TGL_cterm2/(cterm1h + cterm2) &
                 -cterm3h*TGL_cterm2/(cterm1h + cterm2)
       tmp1r   =  real(cterm3v)
   TGL_tmp1r   =  real(TGL_cterm3v)
       tmp1i   = -aimag(cterm3v)
   TGL_tmp1i   = -aimag(TGL_cterm3v)
   ! ev      =  1.0 - cabs(cterm3v)**2
       ev      =  1.0 - (tmp1r*tmp1r+tmp1i*tmp1i)
   TGL_ev      = -2.0*tmp1r*TGL_tmp1r - 2.0*tmp1i*TGL_tmp1i

       tmp2r   =  real(cterm3h)
   TGL_tmp2r   =  real(TGL_cterm3h)
       tmp2i   = -aimag(cterm3h)
   TGL_tmp2i   = -aimag(TGL_cterm3h)
   ! eh      =  1.0 - cabs(cterm3h)**2
       eh      =  1.0 - (tmp2r*tmp2r+tmp2i*tmp2i)
   TGL_eh      = -2.0*tmp2r*TGL_tmp2r - 2.0*tmp2i*TGL_tmp2i

   if (trace_use) call da_trace_exit("da_spemiss_tl")

end subroutine da_spemiss_tl


subroutine da_effang_tl(ifreq,theta,gx2,sigma,effangv,effangh,     &
   TGL_gx2,TGL_sigma,TGL_effangv,TGL_effangh )
 
   !--------------------------------------------------------------------------
   ! Purpose : Calculate the effective zenith angle of reflected microwave 
   ! radiation at SSM/I frequencies for vertical and horizontal polarization
   !
   ! Input  :: TGL_gx2, TGL_sigma
   ! Output :: TGL_effangv,TGL_effangh,effangv,effangh
   !--------------------------------------------------------------------------

   implicit none

   integer, intent(in)  :: ifreq
   real,    intent(in)  :: theta,gx2,sigma,TGL_gx2, TGL_sigma
   real,    intent(out) :: TGL_effangv,TGL_effangh,effangv,effangh

   real c19v,c19h,c22v,c22h,c37v,c37h,c85v,c85h
   real s19v(6),s19h(6),s22v(6),s22h(6), &
        s37v(6),s37h(6),s85v(6),s85h(6)

   real :: alnsig,gg,ggg,xd,xx
   real :: z1,z2,z3,z4,z5,z6,alnthv
   real :: y,dth,angh,angv,alnthh
   real :: TGL_alnsig,TGL_gg,TGL_ggg,TGL_xd
   real :: TGL_z1,TGL_z2,TGL_z3,TGL_z4,TGL_z5,TGL_z6,TGL_alnthv
   real :: TGL_y,TGL_dth,TGL_angh,TGL_angv,TGL_xx,TGL_alnthh

   data c19v,c19h,c22v,c22h,c37v,c37h,c85v,c85h &
     /-.5108,.5306,-.5108,.5306,-.6931,.1823,-.9163,.3000/
   data s19v /.225E+2,.698E+2,-.238E+2,-.648E+1,.402E+0,.262E+1/
   data s19h /.743E+1,.507E+2,-.206E+2,-.191E+1,.648E-1,.291E+1/
   data s22v /.249E+2,.701E+2,-.240E+2,-.714E+1,.405E+0,.256E+1/
   data s22h /.498E+1,.442E+2,-.190E+2,-.129E+1,.803E-2,.277E+1/
   data s37v /.215E+2,.573E+2,-.211E+2,-.670E+1,.443E+0,.253E+1/
   data s37h /.869E+1,.571E+2,-.257E+2,-.302E+1,.237E+0,.386E+1/
   data s85v /.116E+2,.263E+2,-.101E+2,-.358E+1,.270E+0,.175E+1/
   data s85h /.736E+1,.568E+2,-.254E+2,-.248E+1,.196E+0,.387E+1/

   if (trace_use) call da_trace_entry("da_effang_tl")

   if (gx2 .le. 0.0 .or. sigma .le. 0.0) then
         effangv = theta
     TGL_effangv = 0.0
         effangh = theta
     TGL_effangh = 0.0
     return
   end if
   alnsig = alog(sigma)
   if (abs(sigma) .gt. 0.0) then
      TGL_alnsig = TGL_sigma/sigma
   else
      TGL_alnsig = 0.0
   end if
       gg  = gx2*gx2
   TGL_gg  = 2.0*gx2*TGL_gx2
       ggg = gg*gx2
   TGL_ggg = TGL_gg*gx2 + gg*TGL_gx2

   if (ifreq .eq. 1) then 

          xd =      alnsig - c19v
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg + xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      alnthv =  s19v(1)*z1 + s19v(2)*z2 + s19v(3)*z3 + &
                s19v(4)*z4 + s19v(5)*z5 + s19v(6)*z6
      TGL_alnthv = s19v(1)*TGL_z1 + s19v(2)*TGL_z2 + s19v(3)*TGL_z3 + &
                   s19v(4)*TGL_z4 + s19v(5)*TGL_z5 + s19v(6)*TGL_z6
      alnthv     =     alnthv + 3.611

          xd =      alnsig - c19h
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg + xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      
      alnthh =  s19h(1)*z1 + s19h(2)*z2 + s19h(3)*z3 + &
                s19h(4)*z4 + s19h(5)*z5 + s19h(6)*z6
      TGL_alnthh = s19h(1)*TGL_z1 + s19h(2)*TGL_z2 + s19h(3)*TGL_z3 + &
                   s19h(4)*TGL_z4 + s19h(5)*TGL_z5 + s19h(6)*TGL_z6

      alnthh     =     alnthh + 3.611

   else if (ifreq .eq. 2) then 
          xd =      alnsig - c22v
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg 
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg + xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      alnthv =  s22v(1)*z1 + s22v(2)*z2 + s22v(3)*z3 + &
                s22v(4)*z4 + s22v(5)*z5 + s22v(6)*z6
      TGL_alnthv = s22v(1)*TGL_z1 + s22v(2)*TGL_z2 + s22v(3)*TGL_z3 + &
                   s22v(4)*TGL_z4 + s22v(5)*TGL_z5 + s22v(6)*TGL_z6
      alnthv     =     alnthv + 3.611
      ! TGL_alnthv = TGL_alnthv

          xd =      alnsig - c22h
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg + xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      alnthh =  s22h(1)*z1 + s22h(2)*z2 + s22h(3)*z3 + &
                s22h(4)*z4 + s22h(5)*z5 + s22h(6)*z6
      TGL_alnthh = s22h(1)*TGL_z1 + s22h(2)*TGL_z2 + s22h(3)*TGL_z3 + &
                   s22h(4)*TGL_z4 + s22h(5)*TGL_z5 + s22h(6)*TGL_z6
      alnthh     =     alnthh + 3.611
      ! TGL_alnthh = TGL_alnthh
   else if (ifreq .eq. 3) then 
          xd =      alnsig - c37v
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg  + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg  + xx*TGL_gg
       z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      alnthv =  s37v(1)*z1 + s37v(2)*z2 + s37v(3)*z3 + &
                s37v(4)*z4 + s37v(5)*z5 + s37v(6)*z6
      TGL_alnthv = s37v(1)*TGL_z1 + s37v(2)*TGL_z2 + s37v(3)*TGL_z3 + &
                   s37v(4)*TGL_z4 + s37v(5)*TGL_z5 + s37v(6)*TGL_z6
      alnthv     =     alnthv + 3.611
      ! TGL_alnthv = TGL_alnthv 

          xd =      alnsig - c37h
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg +  xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg +  xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg  +  xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg  +  xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 +  xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 +  xd*TGL_gx2
      alnthh =  s37h(1)*z1 + s37h(2)*z2 + s37h(3)*z3 + &
                s37h(4)*z4 + s37h(5)*z5 + s37h(6)*z6
      TGL_alnthh = s37h(1)*TGL_z1 + s37h(2)*TGL_z2 + s37h(3)*TGL_z3 + &
                   s37h(4)*TGL_z4 + s37h(5)*TGL_z5 + s37h(6)*TGL_z6
      alnthh     =     alnthh + 3.611
      ! TGL_alnthh = TGL_alnthh
   else if (ifreq .eq. 4) then 
          xd =      alnsig - c85v
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd 
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg  + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg  + xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      alnthv =  s85v(1)*z1 + s85v(2)*z2 + s85v(3)*z3 + &
                   s85v(4)*z4 + s85v(5)*z5 + s85v(6)*z6
      TGL_alnthv = s85v(1)*TGL_z1 + s85v(2)*TGL_z2 + s85v(3)*TGL_z3 + &
                   s85v(4)*TGL_z4 + s85v(5)*TGL_z5 + s85v(6)*TGL_z6
      alnthv     =     alnthv + 3.611
      ! TGL_alnthv = TGL_alnthv 

          xd =      alnsig - c85h
      TGL_xd =  TGL_alnsig
          xx =  xd*xd
      TGL_xx =  2.0*xd*TGL_xd
          z1 =  xx*ggg
      TGL_z1 =  TGL_xx*ggg + xx*TGL_ggg
          z2 =  xd*ggg
      TGL_z2 =  TGL_xd*ggg + xd*TGL_ggg
          z3 =  xd*gg
      TGL_z3 =  TGL_xd*gg  + xd*TGL_gg
          z4 =  xx*gg
      TGL_z4 =  TGL_xx*gg  + xx*TGL_gg
          z5 =  xx*gx2
      TGL_z5 =  TGL_xx*gx2 + xx*TGL_gx2
          z6 =  xd*gx2
      TGL_z6 =  TGL_xd*gx2 + xd*TGL_gx2
      alnthh =  s85h(1)*z1 + s85h(2)*z2 + s85h(3)*z3 + &
                s85h(4)*z4 + s85h(5)*z5 + s85h(6)*z6
      TGL_alnthh = s85h(1)*TGL_z1 + s85h(2)*TGL_z2 + s85h(3)*TGL_z3 + &
                   s85h(4)*TGL_z4 + s85h(5)*TGL_z5 + s85h(6)*TGL_z6
      alnthh     =     alnthh + 3.611
      ! TGL_alnthh = TGL_alnthh
   end if
       angv =   90.0 - exp(alnthv)
   TGL_angv = - TGL_alnthv*exp(alnthv)
       angh =   90.0 - exp(alnthh)
   TGL_angh = - TGL_alnthh*exp(alnthh)
       y    =   1.0 - 28.0*gx2
   TGL_y    = - 28.0*TGL_gx2
   if (y .lt. 0.0) then
          y = 0.0
      TGL_y = 0.0
   end if
       dth     = (theta - 53.0)*y
   TGL_dth     = (theta - 53.0)*TGL_y
       effangv =     angv +     dth
   TGL_effangv = TGL_angv + TGL_dth
       effangh =     angh +     dth
   TGL_effangh = TGL_angh + TGL_dth

   if (trace_use) call da_trace_exit("da_effang_tl")

end subroutine da_effang_tl


subroutine da_epsalt_tl(f,t,ssw,epsr,epsi,TGL_t, TGL_epsr, TGL_epsi)

   !-----------------------------------------------------------------------
   ! Purpose: returns the complex dielectric constant of sea water, using the
   !  model of Klein and Swift (1977)
   !
   ! Input   f = frequency (GHz)
   !         t = temperature (C)
   !         ssw = salinity (permil) (if ssw < 0, ssw = 32.54)
   ! Output  epsr,epsi  = real and imaginary parts of dielectric constant
   ! Input  : TGL_t      (ssw is treated as a constant now)
   ! Output : TGL_epsr, TGL_epsi, epsr, epsi
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: f, t, TGL_t
   real, intent(inout) :: ssw
   real, intent(out)   :: TGL_epsr, TGL_epsi, epsr, epsi

   complex :: cdum1,cdum2,cdum3
   complex :: TGL_cdum1,TGL_cdum2,TGL_cdum3
   real    :: ssw2,ssw3,t2,t3,es,a,esnew,tau,b,sig,taunew
   real    :: delt,delt2,beta,signew,om
   real    :: TGL_t2,TGL_t3,TGL_es,TGL_a,TGL_esnew,TGL_tau,TGL_b,TGL_taunew
   real    :: TGL_delt,TGL_delt2,TGL_beta,TGL_signew

   if (trace_use) call da_trace_entry("da_epsalt_tl")

   if (ssw .lt. 0.0) ssw = 32.54
   ssw2       = ssw*ssw
   ssw3       = ssw2*ssw
       t2     = t*t
   TGL_t2     = 2.0*t*TGL_t
       t3     = t2*t
   TGL_t3     = TGL_t2*t + t2*TGL_t
       es     = 87.134 - 1.949e-1*t     - 1.276e-2*t2     + 2.491e-4*t3
   TGL_es     =        - 1.949e-1*TGL_t - 1.276e-2*TGL_t2 + 2.491e-4*TGL_t3

       a      = 1.0 + 1.613e-5*ssw*t - 3.656e-3*ssw + 3.21e-5*ssw2 - &
                4.232e-7*ssw3
   TGL_a      = 1.613e-5*ssw*TGL_t 
       esnew  = es*a
   TGL_esnew  = TGL_es*a + es*TGL_a

       tau    = 1.768e-11 - 6.086e-13*t     + 1.104e-14*t2     - 8.111e-17*t3
   TGL_tau    =           - 6.086e-13*TGL_t + 1.104e-14*TGL_t2 - 8.111e-17*TGL_t3
       b      = 1.0 + 2.282e-5*ssw*t - 7.638e-4*ssw - 7.760e-6*ssw2 + &
                1.105e-8*ssw3
   TGL_b      = 2.282e-5*ssw*TGL_t 
       taunew = tau*b
   TGL_taunew = TGL_tau*b + tau*TGL_b

   sig        = ssw*(0.182521 - 1.46192e-3*ssw + 2.09324e-5*ssw2 - &
                1.28205e-7*ssw3)
       delt   = 25.0 - t
   TGL_delt   =      - TGL_t
       delt2  = delt*delt
   TGL_delt2  = 2.0*delt*TGL_delt
       beta   =   2.033e-2 + 1.266e-4*delt      + 2.464e-6*delt2       &
                - ssw*(1.849e-5 - 2.551e-7*delt + 2.551e-8*delt2)
   TGL_beta   =              1.266e-4*TGL_delt  + 2.464e-6*TGL_delt2   &
                - ssw*(-2.551e-7*TGL_delt + 2.551e-8*TGL_delt2)
   signew     =   sig*exp(-beta*delt)
   TGL_signew = - signew*(TGL_beta*delt+beta*TGL_delt)

   om         = 2.0e9*pi*f
       cdum1  = cmplx(0.0,om*taunew)
   TGL_cdum1  = cmplx(0.0,om*TGL_taunew)
       cdum2  = cmplx(0.0,signew/(om*8.854e-12))
   TGL_cdum2  = cmplx(0.0,TGL_signew/(om*8.854e-12))

   cdum3      = 4.9 + (esnew-4.9)/(1.0 + cdum1) - cdum2
   TGL_cdum3  =  TGL_esnew/(1.0 + cdum1)  &
               - TGL_cdum1*(esnew-4.9)/((1.0 + cdum1)*(1.0 + cdum1))  &
               - TGL_cdum2
   epsr       = real(cdum3)
   TGL_epsr   = real(TGL_cdum3)
   epsi       = -aimag(cdum3)
   TGL_epsi   = -aimag(TGL_cdum3)

   if (trace_use) call da_trace_exit("da_epsalt_tl")

end subroutine da_epsalt_tl


subroutine da_sigma_v_tl(ifreq,p0,wv,hwv,ta,gamma,sigma_v,                &
                           TGL_p0,TGL_wv,TGL_hwv,TGL_ta,TGL_gamma,TGL_sigma_v)

   !---------------------------------------------------------------------------
   ! Purpose : TBD
   ! Input             : TGL_p0, TGL_wv, TGL_hwv, TGL_ta, TGL_gamma
   ! Output            : TGL_sigma_v
   ! Output base field : sigma_v
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: ifreq
   real, intent(in  ):: p0,wv,hwv,ta,gamma  ! base field
   real, intent(in  ):: TGL_p0,TGL_wv,TGL_hwv,TGL_ta,TGL_gamma
   real, intent(out ):: TGL_sigma_v,sigma_v

   real wvc, wvcor(4)
   real TGL_wvc

   real voh1,otbar1,pbar1
   real term21,term31,term41,term51,term61
   real a11,a21,a31,a41,a51,a61
   real TGL_voh1,TGL_otbar1,TGL_pbar1
   real TGL_term21,TGL_term31,TGL_term41,TGL_term51,TGL_term61

   real voh2,otbar2,pbar2
   real term22,term32,term42,term52,term62
   real a12,a22,a32,a42,a52,a62
   real TGL_voh2,TGL_otbar2,TGL_pbar2
   real TGL_term22,TGL_term32,TGL_term42,TGL_term52,TGL_term62

   real voh3,otbar3,pbar3
   real term23,term33,term43,term53,term63
   real a13,a23,a33,a43,a53,a63
   real TGL_voh3,TGL_otbar3,TGL_pbar3
   real TGL_term23,TGL_term33,TGL_term43,TGL_term53,TGL_term63

   real voh4,otbar4,pbar4
   real term24,term34,term44,term54,term64
   real a14,a24,a34,a44,a54,a64
   real TGL_voh4,TGL_otbar4,TGL_pbar4
   real TGL_term24,TGL_term34,TGL_term44,TGL_term54,TGL_term64

   real const1,const2,const3,const4
   real h1,h2,h3,h4

   real sigv, TGL_sigv

   data const1,const2,const3,const4/0.6,2.8,0.2,0.2/
   data h1,h2,h3,h4/5.0,4.9,6.8,6.4/

   data a11,a21,a31,a41,a51,a61/-.13747e-2,-.43061e-4, .14618e+1,  &
     .25101e-3, .14635e-1,-.18588e+3/
   data a12,a22,a32,a42,a52,a62/ .22176e-1,-.32367e-4,-.10840e-4,  &
     -.63578e-1, .16988e-7,-.29824e+2/
   data a13,a23,a33,a43,a53,a63/-.10566e-2,-.12906e-3, .56975e+0,  &
      .10828e-8,-.17551e-7, .48601e-1/
   data a14,a24,a34,a44,a54,a64/-.60808e-2,-.70936e-3, .28721e+1,  &
      .42636e-8,-.82910e-7, .26166e+0/

   ! data wvcor/1.01,0.95,1.06,0.92/
   data wvcor/1.02,0.98,1.02,0.88/

   if (trace_use) call da_trace_entry("da_sigma_v_tl")

   ! use modified water vapor value to correct for errors in theoretical absorption

   wvc     = wv*wvcor(ifreq)
   TGL_wvc = TGL_wv*wvcor(ifreq)

   if (ifreq==1) then
      pbar1 = p0/(1.0 + hwv/h1)
      TGL_pbar1  = TGL_p0/(1.0 + hwv/h1)-pbar1*TGL_hwv/(h1*(1.0 + hwv/h1))
      voh1       = wv/hwv
      TGL_voh1   = TGL_wv/hwv-voh1*TGL_hwv/hwv
      term21     = a21*voh1
      TGL_term21 = a21*TGL_voh1
      otbar1     =  1.0/(ta - const1*gamma*hwv)
      TGL_otbar1 = -otbar1*otbar1*(TGL_ta-const1*gamma*TGL_hwv - const1*TGL_gamma*hwv)
      term31     = a31*otbar1
      TGL_term31 = a31*TGL_otbar1
      term61     = a61*otbar1*otbar1
      TGL_term61 = 2.0*a61*otbar1*TGL_otbar1
      term41     = a41*pbar1*otbar1
      TGL_term41 = a41*(TGL_pbar1*otbar1+pbar1*TGL_otbar1)
      term51     = a51*voh1*otbar1
      TGL_term51 = a51*(TGL_voh1*otbar1+voh1*TGL_otbar1)
      sigv       = a11 + term21 + term31 + term41 + term51 + term61
      TGL_sigv   = TGL_term21+TGL_term31+TGL_term41+TGL_term51+TGL_term61

   else if (ifreq==2) then
      pbar2      = p0/(1.0 + hwv/h2)
      TGL_pbar2  = TGL_p0/(1.0 + hwv/h2)-pbar2*TGL_hwv/h2/(1.0 + hwv/h2)
      term22     = a22*pbar2
      TGL_term22 = a22*TGL_pbar2
      term52     = a52*pbar2*pbar2
      TGL_term52 = 2.0*a52*pbar2*TGL_pbar2
      voh2       = wv/hwv
      TGL_voh2   = TGL_wv/hwv-voh2*TGL_hwv/hwv
      term32     = a32*voh2
      TGL_term32 = a32*TGL_voh2
      otbar2     = 1.0/(ta - const2*gamma*hwv)
      TGL_otbar2 = -otbar2*otbar2*(TGL_ta-const2*gamma*TGL_hwv &
                                           -const2*TGL_gamma*hwv)
      term42     = a42*otbar2
      TGL_term42 = a42*TGL_otbar2
      term62     = a62*otbar2*otbar2
      TGL_term62 = 2.0*a62*otbar2*TGL_otbar2
      sigv       = a12 + term22 + term32 + term42 + term52 + term62
      TGL_sigv   = TGL_term22 + TGL_term32 + TGL_term42 + TGL_term52 + TGL_term62

   else if (ifreq == 3) then
      pbar3      = p0/(1.0 + hwv/h3)
      TGL_pbar3  = TGL_p0/(1.0 + hwv/h3)-pbar3*TGL_hwv/h3/(1.0 + hwv/h3)
      term43     = a43*pbar3*pbar3
      TGL_term43 = 2.0*a43*pbar3*TGL_pbar3
      voh3       = wv/hwv
      TGL_voh3   = TGL_wv/hwv-voh3*TGL_hwv/hwv
      term23     = a23*voh3
      TGL_term23 = a23*TGL_voh3
      otbar3     = 1.0/(ta - const3*gamma*hwv)
      TGL_otbar3 = -otbar3*otbar3*(TGL_ta-const3*gamma*TGL_hwv &
                                        -const3*TGL_gamma*hwv)
      term33     = a33*otbar3
      TGL_term33 = a33*TGL_otbar3
      term53     = a53*pbar3*voh3
      TGL_term53 = a53*(TGL_pbar3*voh3+pbar3*TGL_voh3)
      term63     = a63*otbar3*voh3
      TGL_term63 = a63*(TGL_otbar3*voh3+otbar3*TGL_voh3)
      sigv       = a13 + term23 + term33 + term43 + term53 + term63
      TGL_sigv   = TGL_term23 + TGL_term33 + TGL_term43 + TGL_term53 + TGL_term63
   else if (ifreq == 4) then
      pbar4      = p0/(1.0 + hwv/h4)
      TGL_pbar4  = TGL_p0/(1.0 + hwv/h4)-pbar4*TGL_hwv/h4/(1.0 + hwv/h4)
      term44     = a44*pbar4*pbar4
      TGL_term44 = 2.0*a44*pbar4*TGL_pbar4
      voh4       = wv/hwv
      TGL_voh4   = TGL_wv/hwv-voh4*TGL_hwv/hwv
      term24     = a24*voh4
      TGL_term24 = a24*TGL_voh4
      otbar4     = 1.0/(ta - const4*gamma*hwv)
      TGL_otbar4 = -otbar4*otbar4*(TGL_ta-const4*gamma*TGL_hwv -const4*TGL_gamma*hwv)
      term34     = a34*otbar4
      TGL_term34 = a34*TGL_otbar4
      term54     = a54*pbar4*voh4
      TGL_term54 = a54*(TGL_pbar4*voh4+pbar4*TGL_voh4)
      term64     = a64*otbar4*voh4
      TGL_term64 = a64*(TGL_otbar4*voh4+otbar4*TGL_voh4)
      sigv       = a14 + term24 + term34 + term44 + term54 + term64
      TGL_sigv   = TGL_term24 + TGL_term34 + TGL_term44 + TGL_term54 + TGL_term64
   else
      sigv     = 0.0
      TGL_sigv = 0.0
   end if

   sigma_v     = sigv*wvc
   TGL_sigma_v = TGL_sigv*wvc+sigv*TGL_wvc

   if (trace_use) call da_trace_exit("da_sigma_v_tl")

end subroutine da_sigma_v_tl


   
end module da_ssmi

