












module da_radiance

   
   
   

   use hdf5


   use module_domain, only : xb_type, domain
   use module_radiance, only : satinfo, &
      i_kind,r_kind, r_double, &
       one, zero, three,deg2rad,rad2deg, &
      q2ppmv, &
      init_constants_derived, &
      rttov_platform_name, rttov_inst_name, crtm_sensor_name  
   use module_radiance, only : crtm_channelinfo_type, crtm_platform_name, crtm_init, &
   CRTM_Planck_Radiance, CRTM_Planck_Temperature   




   use da_control, only : max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, trace_use_dull, &
      missing, max_error_uv, max_error_t, rootproc, &
      max_error_p,max_error_q, radiance, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, &
      max_error_bt, max_error_buv, rtminit_platform,rtminit_satid, &
      rtminit_nsensor,rtminit_sensor,filename_len,read_biascoef,analysis_date, &
      time_window_max,time_window_min,print_detail_obs,use_hsbobs,use_msuobs, &
      use_amsubobs,use_eos_amsuaobs,use_amsuaobs,use_hirs2obs,rtm_option, &
      rtm_option_rttov,rtm_option_crtm,use_airsobs,use_kma1dvar,use_hirs3obs, &
      use_ssmisobs,use_iasiobs,use_seviriobs,use_filtered_rad,print_detail_rad,stderr, mw_emis_sea, &
      rtminit_print, rttov_scatt,comm,root,ierr,biasprep, qc_rad, num_procs, &
      tovs_min_transfer,use_error_factor_rad,num_fgat_time,stdout,trace_use, &
      qc_good, qc_bad,myproc,biascorr,thinning,thinning_mesh, &
      rad_monitoring, monitor_on, kts, kte, kms, kme, calc_weightfunc, &
      use_mwtsobs, use_mwhsobs, use_atmsobs, use_amsr2obs, &
      use_hirs4obs, use_mhsobs,bufr_year, bufr_month,bufr_day,bufr_hour, &
      bufr_minute, bufr_second,bufr_solzen, bufr_station_height, &
      bufr_landsea_mask,bufr_solazi,tovs_end, max_tovs_input, bufr_satzen, nchan_mhs, &
      nchan_msu, nchan_amsua,nchan_hirs2, nchan_hirs3, nchan_hirs4, nchan_airs, &
      bufr_lon, bufr_satellite_id, bufr_ifov,nchan_amsub, tovs_start, bufr_lat, &
      use_pseudo_rad, pseudo_rad_platid,pseudo_rad_satid, pseudo_rad_senid, &
      pseudo_rad_ichan, pseudo_rad_inv, pseudo_rad_lat,pseudo_rad_lon, &
      pseudo_rad_err, use_simulated_rad,use_rttov_kmatrix, use_crtm_kmatrix , &
      use_rad,crtm_cloud, DT_cloud_model, global, use_varbc, freeze_varbc, &
      airs_warmest_fov, time_slots, interp_option, ids, ide, jds, jde, &
      ips, ipe, jps, jpe, simulated_rad_ngrid, obs_qc_pointer, use_blacklist_rad, use_satcv
 
   use da_crtm, only : da_crtm_init, da_get_innov_vector_crtm
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, j_type, &
      bad_data_type, x_type, number_type, bad_data_type, &
      airsr_type,info_type, model_loc_type, varbc_info_type, varbc_type
   use da_interpolation, only : da_to_zk, da_to_zk_new
   use da_tools_serial, only : da_get_unit, da_free_unit
   use da_par_util1, only : da_proc_sum_int,da_proc_sum_ints
   use da_par_util, only :  da_proc_stats_combine, true_mpi_real
   use da_physics, only : da_sfc_pre, da_transform_xtopsfc, &
      da_transform_xtopsfc_adj,da_tpq_to_slp_lin,da_tpq_to_slp_adj
   use da_radiance1, only : num_tovs_before,num_tovs_after,tovs_copy_count, &
      tovs_send_pe, tovs_recv_pe, tovs_send_start, tovs_send_count, &
      tovs_recv_start,con_vars_type,aux_vars_type, datalink_type,da_qc_amsub, &
      da_qc_amsua,da_biascorr, da_detsurtyp,da_biasprep, &
      da_qc_rad, da_cld_eff_radius, da_read_biascoef
   use da_reporting, only : da_message, da_warning, message, da_error
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_residual, da_obs_sfc_correction, &
      da_llxy, da_llxy_new, da_togrid_new, da_get_julian_time, da_get_time_slots, &
      da_xyll, map_info
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace, &
      da_trace_int_sort
   use da_varbc, only : da_varbc_direct,da_varbc_coldstart,da_varbc_precond, &
      da_varbc_pred
   use da_wrf_interfaces, only : wrf_dm_bcast_integer
   use gsi_thinning, only : r999,r360,rlat_min,rlat_max,rlon_min,rlon_max, &
                            dlat_grid,dlon_grid,thinning_grid, &
                            makegrids,map2grids, &
                            destroygrids
                            
   implicit none

   include 'mpif.h'
   
contains


subroutine da_calculate_grady_rad( iv, re, jo_grad_y ) 
!------------------------------------------------------------------------------
!  Purpose: Calculate Gradient_y for radiance data.
!
!  METHOD:  grad_y = -R^-1 (d - H delta_x ).
!
!  HISTORY: 08/2005 - Creation            Zhiquan Liu
!
!------------------------------------------------------------------------------
   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type) , intent(inout) :: re          ! Residual vector.
   type (y_type) , intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer                       :: n, k, i

   if (trace_use_dull) call da_trace_entry("da_calculate_grady_rad")
   do i =1, iv%num_inst
      if ( iv%instid(i)%num_rad < 1 .or. iv%instid(i)%rad_monitoring == monitor_on ) cycle
      do n=1, iv%instid(i)%num_rad
         do k=1, iv%instid(i)%nchan
            if ( iv%instid(i)%tb_qc(k,n) < obs_qc_pointer) re%instid(i)%tb(k,n) = 0.0

            jo_grad_y%instid(i)%tb(k,n) = -re%instid(i)%tb(k,n) / (iv%instid(i)%tb_error(k,n) * iv%instid(i)%tb_error(k,n))


         end do
      end do
   end do

   if (trace_use_dull) call da_trace_exit("da_calculate_grady_rad")

end subroutine da_calculate_grady_rad

subroutine da_read_filtered_rad( iv )

   !---------------------------------------------------------------------------
   ! Purpose: read in QCed radiance data from separated PEs 
   !---------------------------------------------------------------------------

   implicit none


   type (iv_type), intent(inout)  :: iv       ! O-B structure.

   integer              :: n        ! Loop counter.
   integer              :: i        ! Index dimension.

   integer              :: ios, filtered_rad_unit
   character(len=50)    :: filename
   integer              :: ndomain
   logical, allocatable :: outside(:,:)

   if (trace_use) call da_trace_entry("da_read_filtered_rad")

   do i = 1, iv%num_inst

      write(filename, '(a,i4.4)') 'filtered_'//trim(iv%instid(i)%rttovid_string)//'.', myproc

      if (print_detail_rad) then
         write(unit=stdout,fmt='(2a)') 'Reading in ', trim(filename) 
      end if

      call da_get_unit(filtered_rad_unit)
      open(unit=filtered_rad_unit,file=trim(filename),form='unformatted',status='old',iostat=ios)
      if (ios /= 0) then
         write(unit=stdout,fmt='(2a)') 'Cannot open filtered radiance file ', filename
         iv%instid(i)%num_rad = 0
         cycle
      end if

      read(unit=filtered_rad_unit) iv%instid(i)%num_rad
      if (iv%instid(i)%num_rad < 1) cycle
      if (print_detail_rad) then
         write(unit=stdout,fmt='(a,i3,2x,a,3x,i10)')  ' allocating space for radiance innov structure', &
                               n, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad
      end if

      allocate (iv%instid(i)%tb_inv(1:iv%instid(i)%nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_qc(1:iv%instid(i)%nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_error(1:iv%instid(i)%nchan,iv%instid(i)%num_rad))

      do n=1,iv%instid(i)%num_rad
         read(unit=filtered_rad_unit) ndomain, &
            iv%instid(i)%info%date_char(n), &
            iv%instid(i)%scanpos(n)        , &
            iv%instid(i)%ifgat(n)          , &
            iv%instid(i)%landsea_mask(n)   , &
            iv%instid(i)%info%elv(n)   , &
            iv%instid(i)%info%lat(1,n) , &
            iv%instid(i)%info%lon(1,n) , &
            iv%instid(i)%satzen(n)         , &
            iv%instid(i)%satazi(n)         , &
            iv%instid(i)%tb_inv(:,n)       , &
            iv%instid(i)%tb_error(:,n)     , &
            iv%instid(i)%tb_qc(:,n)
         iv%instid(i)%info%lat(2:,n) = iv%instid(i)%info%lat(1,n)
         iv%instid(i)%info%lon(2:,n) = iv%instid(i)%info%lon(1,n)

         allocate (outside(iv%instid(i)%nlevels, iv%instid(i)%num_rad))
         call da_llxy_new (iv%instid(i)%info, outside)
         deallocate (outside)

      end do     ! end do pixels
      close(filtered_rad_unit)
      call da_free_unit(filtered_rad_unit)

      iv%total_rad_pixel   = iv%total_rad_pixel + iv%instid(i)%num_rad
      iv%total_rad_channel = iv%total_rad_channel + iv%instid(i)%num_rad*iv%instid(i)%nchan
   end do      ! end do instruments

   if (trace_use) call da_trace_exit("da_read_filtered_rad")
 
end subroutine da_read_filtered_rad

subroutine da_read_simulated_rad (iv)

   !---------------------------------------------------------------------------
   !  Purpose: Generate simulated radiances for every model grid point and 
   !           every channel
   !
   !  Called from 
   !
   !  HISTORY: 12/12/2008 - Creation                        Tom Auligne
   !---------------------------------------------------------------------------

   use da_control

   implicit none

   type (iv_type),  intent (inout) :: iv

   ! Instrument triplet, follow the convension of RTTOV 
   integer   :: platform_id, satellite_id, sensor_id

   type (datalink_type) :: p

   logical :: outside, outside_all
   integer :: i,j,k,nchan,inst, alloc_stat
   type(info_type)       ::  info
   type(model_loc_type)  ::  loc
 
   call da_trace_entry("da_read_simulated_rad")

   ! Initialize variables

   platform_id  = pseudo_rad_platid
   satellite_id = pseudo_rad_satid
   sensor_id    = pseudo_rad_senid
   if (sensor_id == 0) then
      nchan=19 !nchan_hirs
   else if (sensor_id == 1) then
      nchan=nchan_msu
   else if (sensor_id == 3) then
      nchan=nchan_amsua
   else if (sensor_id == 4)  then
      nchan=nchan_amsub
   else if (sensor_id == 15)  then
      nchan=nchan_mhs
   else if (sensor_id == 10)  then
      nchan=nchan_ssmis
   else if (sensor_id == 11)  then
      nchan=nchan_airs
   end if

   inst = 1    ! single instrument

   iv%info(radiance)%ntotal    = (ide - ids + 1) * (jde - jds + 1)
   iv%info(radiance)%nlocal    = (min(ipe+simulated_rad_ngrid,ide) - max(ips-simulated_rad_ngrid,ids) + 1) * &
                                 (min(jpe+simulated_rad_ngrid,jde) - max(jps-simulated_rad_ngrid,jds) + 1)
   iv%info(radiance)%ptotal(1) = iv%info(radiance)%ntotal
   iv%instid(inst)%num_rad     = iv%info(radiance)%nlocal
   iv%instid(inst)%info%nlocal = iv%info(radiance)%nlocal

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   
   if ( iv%instid(inst)%num_rad > 0 ) then

      call da_allocate_rad_iv(inst,nchan,iv)

      !  6.0 assign sequential structure to innovation structure
      !-------------------------------------------------------------

      allocate (p%tb_inv(1:nchan), stat=alloc_stat)
      if ( alloc_stat /= 0 ) CALL da_error("da_read_simulated_rad.inc",71,(/"error allocating"/))

      p%info%date_char   = "0000-00-00_00:00:00"         
      p%landsea_mask     = 1
      p%scanpos          = 1
      p%satzen           = 0.0
      p%satazi           = 0.0
      p%solzen           = 0.0
      p%tb_inv(1:nchan)  = pseudo_rad_inv
      p%sensor_index     = inst          
      p%ifgat            = 1
     
      k = 0
      do i = ids, ide
         do j = jds, jde
	    if ((i<ips-simulated_rad_ngrid).or.(i>ipe+simulated_rad_ngrid).or. &
	        (j<jps-simulated_rad_ngrid).or.(j>jpe+simulated_rad_ngrid)) cycle
	    k = k + 1
            call da_xyll(map_info, i*1.0, j*1.0, p%info%lat, p%info%lon)
	    p%loc%x   = float(i)
	    p%loc%y   = float(j)
	    p%loc%i   = i
	    p%loc%j   = j
	    p%loc%dx  = 0.0
	    p%loc%dxm = 1.0
	    p%loc%dy  = 0.0
	    p%loc%dym = 1.0	    
            call da_initialize_rad_iv (inst, k, iv, p)
	 end do
      end do

      iv%instid(inst)%tb_error(:,:) = pseudo_rad_err
      iv%instid(inst)%tb_qc(:,:)    = qc_good

      deallocate(p%tb_inv)
   end if

      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad

   call da_trace_exit("da_read_simulated_rad")
   

end subroutine da_read_simulated_rad
subroutine da_write_filtered_rad(ob, iv)

   !---------------------------------------------------------------------------
   ! Purpose: write out QCed radiance data.
   !
   !  METHOD: write out for separate PEs 
   !---------------------------------------------------------------------------

   implicit none

   type (y_type),     intent(in)  :: ob       ! Observation structure.
   type (iv_type),    intent(in)  :: iv       ! O-B structure.

   integer                        :: n        ! Loop counter.
   integer                        :: i        ! Index dimension.

   integer            :: ios, filtered_rad_unit
   character(len=50)  :: filename

   if (trace_use) call da_trace_entry("da_write_filtered_rad")

   do i = 1, iv%num_inst
      if (iv%instid(i)%num_rad < 1) cycle

      write(unit=filename, fmt='(a,i4.4)') &
         'filtered_'//trim(iv%instid(i)%rttovid_string)//'.', myproc
      
      call da_get_unit(filtered_rad_unit)
      open(unit=filtered_rad_unit,file=trim(filename), &
         form='unformatted',iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_filtered_rad.inc",37, &
            (/"Cannot open filtered radiance file"//filename/))
      Endif

      write(unit=filtered_rad_unit) iv%instid(i)%num_rad

      do n =1,iv%instid(i)%num_rad
         write(unit=filtered_rad_unit) n, &
                   iv%instid(i)%info%date_char(n), &
                   iv%instid(i)%scanpos(n)        , &
                   iv%instid(i)%ifgat(n)          , &
                   iv%instid(i)%landsea_mask(n)   , &
                   iv%instid(i)%info%elv(n)   , &
                   iv%instid(i)%info%lat(1,n) , &
                   iv%instid(i)%info%lon(1,n) , &
                   iv%instid(i)%satzen(n)         , &
                   iv%instid(i)%satazi(n)         , &
                   ob%instid(i)%tb(:,n)           , &
                   iv%instid(i)%tb_error(:,n)     , &
                   iv%instid(i)%tb_qc(:,n)

      end do     ! end do pixels
      close(unit=filtered_rad_unit)
      call da_free_unit(filtered_rad_unit)
   end do    !! end do instruments

   if (trace_use) call da_trace_exit("da_write_filtered_rad")
 
end subroutine da_write_filtered_rad


subroutine da_read_obs_bufrtovs (obstype,iv, infile)

   !---------------------------------------------------------------------------
   !  Purpose: read in NCEP bufr tovs 1b data to innovation structure
   !
   !   METHOD: use F90 sequential data structure to avoid reading file twice  
   !            so that da_scan_bufrtovs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !               and deallocate sequential data structure
   !   File History:
   !      Peiming Dong Added NPP atms, 2012/4/17
   !---------------------------------------------------------------------------

   implicit none

   character(5)      ,  intent (in)    :: obstype
   character(20)     ,  intent (in)    :: infile
   type (iv_type)    ,  intent (inout) :: iv


   integer          :: iost
   integer(i_kind), allocatable :: nread(:)

   integer(i_kind),parameter:: n1bhdr=15
   integer(i_kind),parameter:: n2bhdr=2
   integer(i_kind),parameter:: maxinfo=12
   integer(i_kind),parameter:: maxchanl=100

   logical hirs2,hirs3,hirs4,msu,amsua,amsub,mhs,atms
   logical outside, outside_all, iuse
   integer :: inst

   character(10) date
   character(8) subset,subfgn
   character(80) hdr1b
   character(80) hdr2b
   integer(i_kind) ihh,i,j,k,ifov,idd,ireadmg,ireadsb
   integer(i_kind) iret,idate,im,iy,nchan
   integer :: num_bufr(7), numbufr, ibufr
   character(20) :: filename

   ! thinning variables
   integer(i_kind) itt,itx,iobs,iout
   real(r_kind) terrain,timedif,crit,dist
   real(r_kind) dlon_earth,dlat_earth

   real(r_kind) tbmin,tbmax, tbbad
   real(r_kind) panglr,rato
   ! real(r_kind) rmask
   real(r_kind) step,start

   real(r_double),dimension(maxinfo+maxchanl):: data1b8
   real(r_double),dimension(n1bhdr):: bfr1bhdr
   real(r_double),dimension(n2bhdr):: bfr2bhdr

   ! Instrument triplet, follow the convension of RTTOV 
   integer   :: platform_id, satellite_id, sensor_id

   ! pixel information
   integer   ::  year,month,day,hour,minute,second  ! observation time
   real*8    ::  obs_time
   ! real      ::  rlat, rlon                         !  lat/lon in degrees   for Anfovs
   real      ::  satzen, satazi, solzen ,solazi       !  scan angles for Anfovs
   integer   ::  landsea_mask
   real      ::  srf_height
   ! channels' bright temperature
   real , allocatable ::   tb_inv(:)                    !  bright temperatures
   !  end type bright_temperature

   type (datalink_type), pointer    :: head, p, current, prev

   integer                        ::  ifgat
   type(info_type)                ::  info
   type(model_loc_type)           ::  loc

   data hdr1b /'SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA HOLS LSQL SOLAZI'/
   data hdr2b /'CLATH CLONH'/
   !  data hdr1b /'FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA HOLS LSQL SLNM BEARAZ'/

   data tbmin,tbmax,tbbad / 50.0_r_kind, 550.0_r_kind, -9.99e11_r_kind /
   integer :: num_tovs_local, num_tovs_file, num_tovs_global, num_tovs_selected
   integer :: num_tovs_thinned, num_tovs_used, num_tovs_used_tmp
   integer :: lnbufr
   integer :: n
   integer(i_kind), allocatable :: ptotal(:)
   real , allocatable :: in(:), out(:)
   logical :: found, head_found

   call da_trace_entry("da_read_obs_bufrtovs")

   ! Initialize variables

   nchan = 22
   allocate(nread(1:rtminit_nsensor))
   allocate(ptotal(0:num_fgat_time))
   nread(1:rtminit_nsensor) = 0
   ptotal(0:num_fgat_time) = 0

   ! Set various variables depending on type of data to be read

   ! platform_id  = 1                 !! for NOAA series
   ! platform_id  = 10                !! for METOP series

   hirs2=     obstype == 'hirs2'
   hirs3=     obstype == 'hirs3'
   hirs4=     obstype == 'hirs4'
   msu=       obstype == 'msu  '
   amsua=     obstype == 'amsua'
   amsub=     obstype == 'amsub'
   mhs=       obstype == 'mhs  '
   atms=      obstype == 'atms '

   if (hirs2) then
      sensor_id    =  0
      step   = 1.80_r_kind
      start  = -49.5_r_kind
      nchan=nchan_hirs2
      subfgn='NC021021'
      rato=1.1363987_r_kind
   else if (hirs3) then 
      sensor_id    =  0
      step   = 1.80_r_kind
      start  = -49.5_r_kind
      nchan=nchan_hirs3
      subfgn='NC021025'
   else if (hirs4) then 
      sensor_id    =  0
      step   = 1.80_r_kind
      start  = -49.5_r_kind
      nchan=nchan_hirs4
      subfgn='NC021028'
   else if (mhs) then 
      sensor_id    =  15
      step   = 10.0_r_kind/9.0_r_kind
      start  = -445.0_r_kind/9.0_r_kind
      nchan=nchan_mhs
      subfgn='NC021027'
   else if (msu) then
      sensor_id    =  1
      step   = 9.474_r_kind
      start  = -47.37_r_kind
      nchan=nchan_msu
      subfgn='NC021022'
      rato=1.1363987_r_kind
   else if (amsua) then
      sensor_id    =  3
      step   = three + one/three
      start  = -48.33_r_kind
      nchan=nchan_amsua
      subfgn='NC021023'
   else if (amsub)  then
      sensor_id    =  4
      step   = 1.1_r_kind
      start  = -48.95_r_kind
      nchan=nchan_amsub
      subfgn='NC021024'
   else if (atms) then
      sensor_id    =  19
      step   = 1.11_r_kind
      start  = -52.725_r_kind
      nchan=22
      subfgn='NC021203'
      hdr1b='SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA HMSL LSQL SOLAZI'
   end if

   allocate (tb_inv(nchan))

   num_tovs_file     = 0    ! number of obs in file
   num_tovs_global   = 0    ! number of obs within whole domain
   num_tovs_local    = 0    ! number of obs within tile
   num_tovs_thinned  = 0    ! number of obs rejected by thinning
   num_tovs_used     = 0    ! number of obs entered into innovation computation
   num_tovs_selected = 0    ! number of obs limited for debuging
   iobs = 0                 ! for thinning, argument is inout

   ! 0.0  Open unit to satellite bufr file and read file header
   !--------------------------------------------------------------

   num_bufr(:)=0
   numbufr=0
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
	    num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if
  
   if (numbufr==0) numbufr=1

bufrfile:  do ibufr=1,numbufr   
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'   
      end if
   end if

!  We want to use specific unit number for bufr data, so we can control the endian format in environment. 
   lnbufr = 99

   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_warning("da_read_obs_bufrtovs.inc",220, &
         (/"Cannot open file "//infile/))
      call da_trace_exit("da_read_obs_bufrtovs")
      return
   end if

   call openbf(lnbufr,'IN',lnbufr)
   call datelen(10)
   call readmg(lnbufr,subset,idate,iret)
   if (subset /= subfgn) then
      call closbf(lnbufr)
      close(lnbufr)
      message(1)='The file title does not match the data subset'
      write(unit=message(2),fmt=*) &
         'infile=', lnbufr, infile,' subset=', subset, ' subfgn=',subfgn
      call da_error("da_read_obs_bufrtovs.inc",235,message(1:2))
   end if

   iy=0
   im=0
   idd=0
   ihh=0
   write(unit=date,fmt='( i10)') idate
   read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
   write(unit=stdout,fmt=*) &
      'Bufr file date is ',iy,im,idd,ihh,infile

   ! Loop to read bufr file and assign information to a sequential structure
   !-------------------------------------------------------------------------

   if ( ibufr == 1 ) then
      allocate (head)
      !  allocate ( head % tb_inv (1:nchan) )
      nullify  ( head % next )
      p => head
   endif

   if (tovs_start > 1) then
      write (unit=stdout,fmt='(A,I6)') "   Skipping tovs obs before", tovs_start
   end if


   obs: do while (ireadmg(lnbufr,subset,idate)==0 .and. subset==subfgn)
      do while (ireadsb(lnbufr)==0)

         ! 1.0     Read header record and data record
         call ufbint(lnbufr,bfr1bhdr,n1bhdr,1,iret,hdr1b)
         call ufbint(lnbufr,bfr2bhdr,n2bhdr,1,iret,hdr2b)
         call ufbrep(lnbufr,data1b8,1,nchan,iret,'TMBR')
         ! call ufbrep(lnbufr,data1b8,1,1,iret,'BEARAZ')

         ! check if observation outside range

         num_tovs_file = num_tovs_file + 1

         ! 2.0     Extract observation location and other required information
         !     QC1:  judge if data is in the domain, read next record if not
         !------------------------------------------------------------------------
         ! rlat = bfr1bhdr(bufr_lat)
         ! rlon = bfr1bhdr(bufr_lat)
         ! if (rlon < 0.0) rlon = rlon+360.0

         if(abs(bfr2bhdr(1)) <= 90. .and. abs(bfr2bhdr(2)) <= 360.)then
            info%lat = bfr2bhdr(1)
            info%lon = bfr2bhdr(2)
         elseif(abs(bfr1bhdr(9)) <= 90. .and. abs(bfr1bhdr(10)) <= 360.)then
            info%lat = bfr1bhdr(9)
            info%lon = bfr1bhdr(10)
         endif

         call da_llxy (info, loc, outside, outside_all)

         if (outside_all) cycle

         !  3.0     Extract other information
         !------------------------------------------------------
         !  3.1     Extract satellite id and scan position. 
   
         if ( nint(bfr1bhdr(bufr_satellite_id)) >= 200 .and. nint(bfr1bhdr(bufr_satellite_id)) <= 204 ) then
            platform_id = 1
            satellite_id = nint(bfr1bhdr(bufr_satellite_id))-192
         else if ( nint(bfr1bhdr(bufr_satellite_id)) >= 205 .and. nint(bfr1bhdr(bufr_satellite_id)) <= 209 ) then
            platform_id = 1
            satellite_id = nint(bfr1bhdr(bufr_satellite_id))-191
         else if ( nint(bfr1bhdr(bufr_satellite_id)) == 223 ) then ! noaa-19
            platform_id = 1
            satellite_id = nint(bfr1bhdr(bufr_satellite_id))-204
         else if ( nint(bfr1bhdr(bufr_satellite_id)) >= 3 .and. nint(bfr1bhdr(bufr_satellite_id)) <= 5 ) then
            platform_id = 10
            satellite_id = nint(bfr1bhdr(bufr_satellite_id))-2
         else if ( nint(bfr1bhdr(bufr_satellite_id)) == 708 ) then ! tiros-n
            platform_id = 19
            satellite_id = 0
         else if ( nint(bfr1bhdr(bufr_satellite_id)) == 224 ) then ! npp
            platform_id = 17
            satellite_id = 0
         else if ( nint(bfr1bhdr(bufr_satellite_id)) >= 701 .and. nint(bfr1bhdr(bufr_satellite_id)) <= 707 ) then
            platform_id = 1
            satellite_id = nint(bfr1bhdr(bufr_satellite_id))-700
         end if
         ifov = nint(bfr1bhdr(bufr_ifov))    


         !  QC2:  limb pixel rejected (not implemented)

         !  3.2     Extract date information.
    
         year   = bfr1bhdr(bufr_year)   
         month  = bfr1bhdr(bufr_month)  
         day    = bfr1bhdr(bufr_day)    
         hour   = bfr1bhdr(bufr_hour)   
         minute = bfr1bhdr(bufr_minute) 
         second = bfr1bhdr(bufr_second) 
         
         write(unit=info%date_char, fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
            year, '-', month, '-', day, '_', hour, ':', minute, ':', second

         !  QC3: time consistency check with the background date

         if (year <= 99) then
            if (year < 78) then
               year = year + 2000
            else
               year = year + 1900
            end if
         end if

         call da_get_julian_time(year,month,day,hour,minute,obs_time)

         if (obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time)) cycle

         ! 3.2.1 determine FGAT index ifgat
   
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do

         ! 3.3 Find wrfvar instrument index from RTTOV instrument triplet
         !     go to next data if id is not in the lists

         inst = 0
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do
         if (inst == 0) cycle

         ! 3.4 extract satellite and solar angle
   
         panglr=(start+float(ifov-1)*step)*deg2rad
         if (hirs2 .or. msu) then
            satzen = asin(rato*sin(panglr))*rad2deg
            satzen = abs(satzen)
         else
            satzen = bfr1bhdr(bufr_satzen) !*deg2rad   ! local zenith angle
            satzen = abs(satzen)
 
            ! if (amsua .and. ifov .le. 15) satzen=-satzen
            ! if (amsub .and. ifov .le. 45) satzen=-satzen
            ! if (hirs3 .and. ifov .le. 28) satzen=-satzen
         end if
         if ( satzen > 65.0 ) cycle   ! 1 has a satzen > 65.0 check
         satazi = panglr*rad2deg            ! look angle
         ! if (satazi<0.0) satazi = satazi+360.0
         solzen = bfr1bhdr(bufr_solzen)              ! solar zenith angle
         solazi = bfr1bhdr(bufr_solazi)              !RTTOV9_3

         num_tovs_global = num_tovs_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1

         if (num_tovs_global < tovs_start) then
            cycle
         end if

         if (num_tovs_global > tovs_end) then
            write (unit=stdout,fmt='(A,I6)') "   Skipping radiance obs after", tovs_end
            exit obs
         end if

         num_tovs_selected = num_tovs_selected + 1
 
         if (num_tovs_selected > max_tovs_input) then
            write(unit=message(1),fmt='(A,I10,A)') &
               "Max number of tovs",max_tovs_input," reached"
            call da_warning("da_read_obs_bufrtovs.inc",410,message(1:1))
            num_tovs_selected = num_tovs_selected-1
            num_tovs_global   = num_tovs_global-1
            ptotal(ifgat) = ptotal(ifgat) - 1
            exit obs
         end if

         if (outside) cycle ! No good for this PE
         num_tovs_local = num_tovs_local + 1

         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad           
            timedif = 0.0 !2.0_r_kind*abs(tdiff)        ! range:  0 to 6
            terrain = 0.01_r_kind*abs(bfr1bhdr(13))
            crit = 1.0 !0.01_r_kind+terrain + timedif !+ 10.0_r_kind*float(iskip)
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_tovs_thinned=num_tovs_thinned+1
               cycle
            end if
         end if

         num_tovs_used = num_tovs_used + 1
         nread(inst) = nread(inst) + 1

         ! 3.5 extract surface information
         srf_height = bfr1bhdr(bufr_station_height)          ! station height
         if (srf_height < 8888.0 .AND. srf_height > -416.0) then
         else
            srf_height = 0.0
         endif  

         landsea_mask = nint(bfr1bhdr(bufr_landsea_mask))  ! 0:land ; 1:sea (same as RTTOV)
         ! rmask=one                          ! land
         ! if (nint(bfr1bhdr(bufr_landsea_mask))==1) rmask=0.0_r_kind   ! reverse the land/sea mask in bufr
         ! landsea_mask = rmask+.001_r_kind             ! land sea mask

         info%elv = srf_height

         ! 3.6 extract channel bright temperature
   
         tb_inv(1:nchan) = data1b8(1:nchan)
         do k = 1, nchan
            if ( tb_inv(k) < tbmin .or. tb_inv(k) > tbmax) &
               tb_inv(k) = missing_r
         end do
         if ( all(tb_inv<0.0) ) then
            num_tovs_local = num_tovs_local -1
            num_tovs_used = num_tovs_used - 1
            nread(inst) = nread(inst) - 1
            cycle
         end if

         !  4.0   assign information to sequential radiance structure
         !--------------------------------------------------------------------------
         allocate (p % tb_inv (1:nchan))
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = landsea_mask
         p%scanpos          = ifov
         p%satzen           = satzen
         p%satazi           = satazi
         p%solzen           = solzen
         p%tb_inv(1:nchan)  = tb_inv(1:nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat
!RTTOV9_3
         p%solazi           = solazi
 !end of RTTOV9_3
         allocate (p%next)   ! add next data

         p => p%next
         nullify (p%next)
      end do
   end do obs

   call closbf(lnbufr)
   close(lnbufr)

end do bufrfile

   if (thinning .and. num_tovs_global > 0 ) then


      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do
      end do

      allocate ( in  (j) )
      allocate ( out (j) )

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i) 
            end do
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time 
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i)  = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_tovs_used_tmp = num_tovs_used
      do j = 1, num_tovs_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat 
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            endif
         end do 
        
         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            endif
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_tovs_thinned = num_tovs_thinned + 1
            num_tovs_used = num_tovs_used - 1
            nread(n) = nread(n) - 1
            continue
         endif

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         endif

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         endif
        
      end do

   endif  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel   + num_tovs_used
   iv%total_rad_channel = iv%total_rad_channel + num_tovs_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_tovs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_tovs_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1) 
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i) 
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_bufrtovs.inc",607,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'num_tovs_file num_tovs_global num_tovs_local num_tovs_used num_tovs_thinned'
   write(unit=stdout,fmt='(5i10)') num_tovs_file,num_tovs_global, num_tovs_local,num_tovs_used,num_tovs_thinned

   deallocate(tb_inv)  

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   
   do i = 1, iv%num_inst
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad

      call da_allocate_rad_iv(i,nchan,iv)

   end do
   
   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0 
   p => head
   ! do while ( associated(p) )

   do n = 1, num_tovs_used
      i = p%sensor_index
      nread(i) = nread(i) + 1

      call da_initialize_rad_iv (i, nread(i), iv, p)

      current => p
      p => p%next

      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do

   deallocate ( p )

   deallocate (nread)
   deallocate (ptotal)

   ! check if sequential structure has been freed
   !
   ! p => head
   ! do i = 1, num_rad_selected
   !    write (unit=stdout,fmt=*)  i, p%tb_inv(1:nchan)
   !    p => p%next
   ! end do

   call da_trace_exit("da_read_obs_bufrtovs")
  

end subroutine da_read_obs_bufrtovs


subroutine da_read_obs_fy3 (obstype,iv, infile)

   !---------------------------------------------------------------------------
   !  Purpose: read in fy3 1b data to innovation structure
   !
   !   Dong peiming 2012/03/09
   !   METHOD: use F90 sequential data structure to avoid reading file twice  
   !            so that da_scan_bufrtovs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !               and deallocate sequential data structure
   !---------------------------------------------------------------------------

   implicit none

   character(5)      ,  intent (in)    :: obstype
   character(20)     ,  intent (in)    :: infile
   type (iv_type)    ,  intent (inout) :: iv

!
TYPE type_rad_FY3
    INTEGER :: yyyy,mn,dd,hh,mm,ss
    INTEGER :: iscanline,iscanpos
    REAL*4  :: rlat,rlon !lat/lon in degrees for Anfovs
    INTEGER :: isurf_height, isurf_type !height/type for Anfovs
    REAL*4  :: satzen,satazi,solzen,solazi !scan angles for Anfovs
    REAL*4  :: tbb(20) !bright temperatures
!   REAL*4  :: btemps(20)
    INTEGER :: iavhrr(13),ihirsflag
    INTEGER :: iprepro(5) ! values from pre-processing
    REAL*4  :: clfra ! Cloud cover (<1.0)
    REAL*4  :: ts ! Skin temperature
    REAL*4  :: tctop ! Cloud top temperature
END TYPE type_rad_FY3

  TYPE (type_rad_FY3)   :: rad
  integer :: iscan,nscan
!
   integer          :: iost
   integer(i_kind), allocatable :: nread(:)

!dongpm   logical hirs2,hirs3,hirs4,msu,amsua,amsub,mhs
   logical mwts,mwhs
   logical outside, outside_all, iuse
   integer :: inst

   integer(i_kind) i,j,k,ifov
   integer(i_kind) nchan
   integer :: num_bufr(7), numbufr, ibufr
   character(20) :: filename

   ! thinning variables
   integer(i_kind) itt,itx,iobs,iout
   real(r_kind) terrain,timedif,crit,dist
   real(r_kind) dlon_earth,dlat_earth

   real(r_kind) tbmin,tbmax, tbbad
   ! real(r_kind) rmask

   ! Instrument triplet, follow the convension of RTTOV 
   integer   :: platform_id, satellite_id, sensor_id

   ! pixel information
   integer   ::  year,month,day,hour,minute,second  ! observation time
   real*8    ::  obs_time
   ! real      ::  rlat, rlon                         !  lat/lon in degrees   for Anfovs
   real      ::  satzen, satazi, solzen ,solazi       !  scan angles for Anfovs
   integer   ::  landsea_mask
   real      ::  srf_height
   ! channels' bright temperature
   real , allocatable ::   tb_inv(:)                    !  bright temperatures
   !  end type bright_temperature

   type (datalink_type), pointer    :: head, p, current, prev

   integer                        ::  ifgat
   type(info_type)                ::  info
   type(model_loc_type)           ::  loc

   data tbmin,tbmax,tbbad / 50.0_r_kind, 550.0_r_kind, -9.99e11_r_kind /
   integer :: num_tovs_local, num_tovs_file, num_tovs_global, num_tovs_selected
   integer :: num_tovs_thinned, num_tovs_used, num_tovs_used_tmp
   integer :: lnbufr
   integer :: n
   integer(i_kind), allocatable :: ptotal(:)
   real , allocatable :: in(:), out(:)
   logical :: found, head_found

   if (trace_use) call da_trace_entry("da_read_obs_fy3")

   ! Initialize variables

   nchan = 20
   allocate(nread(1:rtminit_nsensor))
   allocate(ptotal(0:num_fgat_time))
   nread(1:rtminit_nsensor) = 0
   ptotal(0:num_fgat_time) = 0

   ! Set various variables depending on type of data to be read

   ! platform_id  = 1                 !! for NOAA series
   ! platform_id  = 10                !! for METOP series

!dongpm   hirs2=     obstype == 'hirs2'
!dongpm   hirs3=     obstype == 'hirs3'
!dongpm   hirs4=     obstype == 'hirs4'
!dongpm   msu=       obstype == 'msu  '
!dongpm   amsua=     obstype == 'amsua'
!dongpm   amsub=     obstype == 'amsub'
!dongpm   mhs=       obstype == 'mhs  '
          mwts=      obstype == 'mwts '
          mwhs=      obstype == 'mwhs '

!dongpm   if (hirs2) then
!dongpm      sensor_id    =  0
!dongpm      step   = 1.80_r_kind
!dongpm      start  = -49.5_r_kind
!dongpm      nchan=nchan_hirs2
!dongpm      subfgn='NC021021'
!dongpm      rato=1.1363987_r_kind
!dongpm   else if (hirs3) then 
!dongpm      sensor_id    =  0
!dongpm      step   = 1.80_r_kind
!dongpm      start  = -49.5_r_kind
!dongpm      nchan=nchan_hirs3
!dongpm      subfgn='NC021025'
!dongpm   else if (hirs4) then 
!dongpm      sensor_id    =  0
!dongpm      step   = 1.80_r_kind
!dongpm      start  = -49.5_r_kind
!dongpm      nchan=nchan_hirs4
!dongpm      subfgn='NC021028'
!dongpm   else if (mhs) then 
!dongpm      sensor_id    =  15
!dongpm      step   = 10.0_r_kind/9.0_r_kind
!dongpm      start  = -445.0_r_kind/9.0_r_kind
!dongpm      nchan=nchan_mhs
!dongpm      subfgn='NC021027'
!dongpm   else if (msu) then
!dongpm      sensor_id    =  1
!dongpm      step   = 9.474_r_kind
!dongpm      start  = -47.37_r_kind
!dongpm      nchan=nchan_msu
!dongpm      subfgn='NC021022'
!dongpm      rato=1.1363987_r_kind
!dongpm   else if (amsua) then
!dongpm      sensor_id    =  3
!dongpm      step   = three + one/three
!dongpm      start  = -48.33_r_kind
!dongpm      nchan=nchan_amsua
!dongpm      subfgn='NC021023'
!dongpm   else if (amsub)  then
!dongpm      sensor_id    =  4
!dongpm      step   = 1.1_r_kind
!dongpm      start  = -48.95_r_kind
!dongpm      nchan=nchan_amsub
!dongpm      subfgn='NC021024'
!dongpm   end if
          if (mwts) then
           sensor_id    =  40
           nchan=4
           nscan=15
          else if(mwhs) then
           sensor_id    =  41
           nchan=5
           nscan=98
          endif

   allocate (tb_inv(nchan))

   num_tovs_file     = 0    ! number of obs in file
   num_tovs_global   = 0    ! number of obs within whole domain
   num_tovs_local    = 0    ! number of obs within tile
   num_tovs_thinned  = 0    ! number of obs rejected by thinning
   num_tovs_used     = 0    ! number of obs entered into innovation computation
   num_tovs_selected = 0    ! number of obs limited for debuging
   iobs = 0                 ! for thinning, argument is inout

   ! 0.0  Open unit to satellite bufr file and read file header
   !--------------------------------------------------------------

   num_bufr(:)=0
   numbufr=0
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.dat'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
            num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if
  
   if (numbufr==0) numbufr=1

bufrfile:  do ibufr=1,numbufr   
   if (num_fgat_time==1) then
      filename=trim(infile)//'.dat'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.dat'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.dat'   
      end if
   end if

!  We want to use specific unit number for bufr data, so we can control the endian format in environment. 
   lnbufr = 99

   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_warning("da_read_obs_fy3.inc",221, &
         (/"Cannot open file "//infile/))
      if (trace_use) call da_trace_exit("da_read_obs_fy3")
      return
   end if


   if ( ibufr == 1 ) then
      allocate (head)
      !  allocate ( head % tb_inv (1:nchan) )
      nullify  ( head % next )
      p => head
   endif

   if (tovs_start > 1) then
      write (unit=stdout,fmt='(A,I6)') "   Skipping tovs obs before", tovs_start
   end if

   obs: do while (.true.)
        do iscan=1,nscan

         ! 1.0     Read fy3 data
         read(lnbufr,end=1000) rad

         num_tovs_file = num_tovs_file + 1

         ! 2.0     Extract observation location and other required information
         !     QC1:  judge if data is in the domain, read next record if not
         !------------------------------------------------------------------------
         ! rlat = bfr1bhdr(bufr_lat)
         ! rlon = bfr1bhdr(bufr_lat)
         ! if (rlon < 0.0) rlon = rlon+360.0

         info%lat  =  rad%rlat

         info%lon  =  rad%rlon
         call da_llxy (info, loc, outside, outside_all)

         if (outside_all) cycle

         !  3.0     Extract other information
         !------------------------------------------------------
         !  3.1     Extract satellite id and scan position. 
   
         platform_id = 23
         if(infile(5:5)=='a') then
            satellite_id = 1
         elseif(infile(5:5)=='b') then
            satellite_id = 2
         else
            call da_warning("da_read_obs_fy3.inc",271,(/"Can not assimilate data from this instrument"/))
            if (trace_use) call da_trace_exit("da_read_obs_fy3")
            return
         endif
         ifov = rad%iscanpos    

         !  QC2:  limb pixel rejected (not implemented)

         !  3.2     Extract date information.
    
         year   = rad%yyyy   
         month  = rad%mn  
         day    = rad%dd    
         hour   = rad%hh   
         minute = rad%mm 
         second = rad%ss 
!dongpm for test
!          year   = 2008
!          month  = 8
!          day    = 5
!          hour   = 18
!          minute = 0
!          second = 0     
         
         write(unit=info%date_char, fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
            year, '-', month, '-', day, '_', hour, ':', minute, ':', second

         !  QC3: time consistency check with the background date

         if (year <= 99) then
            if (year < 78) then
               year = year + 2000
            else
               year = year + 1900
            end if
         end if

         call da_get_julian_time(year,month,day,hour,minute,obs_time)

         if (obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time)) cycle

         ! 3.2.1 determine FGAT index ifgat
   
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do

         ! 3.3 Find wrfvar instrument index from RTTOV instrument triplet
         !     go to next data if id is not in the lists

         inst = 0
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do
         if (inst == 0) cycle

         ! 3.4 extract satellite and solar angle
   
            satzen = rad%satzen !*deg2rad   ! local zenith angle
            satzen = abs(satzen)
 
            ! if (amsua .and. ifov .le. 15) satzen=-satzen
            ! if (amsub .and. ifov .le. 45) satzen=-satzen
            ! if (hirs3 .and. ifov .le. 28) satzen=-satzen
!dongpm         if ( satzen > 65.0 ) cycle   ! 1 has a satzen > 65.0 check
         satazi = rad%satazi*0.01           ! look angle
         ! if (satazi<0.0) satazi = satazi+360.0
         solzen = rad%solzen*0.01              ! solar zenith angle
         solazi = rad%solazi*0.01              !RTTOV9_3

         num_tovs_global = num_tovs_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1

         if (num_tovs_global < tovs_start) then
            cycle
         end if

         if (num_tovs_global > tovs_end) then
            write (unit=stdout,fmt='(A,I6)') "   Skipping radiance obs after", tovs_end
            exit obs
         end if

         num_tovs_selected = num_tovs_selected + 1
 
         if (num_tovs_selected > max_tovs_input) then
            write(unit=message(1),fmt='(A,I10,A)') &
               "Max number of tovs",max_tovs_input," reached"
            call da_warning("da_read_obs_fy3.inc",365,message(1:1))
            num_tovs_selected = num_tovs_selected-1
            num_tovs_global   = num_tovs_global-1
            ptotal(ifgat) = ptotal(ifgat) - 1
            exit obs
         end if

         if (outside) cycle ! No good for this PE
         num_tovs_local = num_tovs_local + 1

         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad           
            timedif = 0.0 !2.0_r_kind*abs(tdiff)        ! range:  0 to 6
!dongpm            terrain = 0.01_r_kind*abs(bfr1bhdr(13))
            terrain = 0.01_r_kind*abs(rad%satzen)
            crit = 1.0 !0.01_r_kind+terrain + timedif !+ 10.0_r_kind*float(iskip)
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_tovs_thinned=num_tovs_thinned+1
               cycle
            end if
         end if

         num_tovs_used = num_tovs_used + 1
         nread(inst) = nread(inst) + 1

         ! 3.5 extract surface information
         srf_height = rad%isurf_height          ! station height
         if (srf_height < 8888.0 .AND. srf_height > -416.0) then
         else
            srf_height = 0.0
         endif  

!dongpm         landsea_mask = rad%isurf_type  ! 0:land ; 1:sea (same as RTTOV)
!fy3 isurf_type is just reversed as RTTOV
         if(rad%isurf_type .eq. 0) then   ! sea
           landsea_mask = 1
         elseif(rad%isurf_type .eq. 1) then   !coast 
           landsea_mask = 0
         elseif(rad%isurf_type .eq. 2) then   !land
           landsea_mask = 0
         else
           landsea_mask = rad%isurf_type
           write(unit=message(1),fmt='(A,I6)') 'Unknown surface type: ', landsea_mask
           call da_warning("da_read_obs_fy3.inc",417,message(1:1))
         endif
         ! rmask=one                          ! land
         ! if (nint(bfr1bhdr(bufr_landsea_mask))==1) rmask=0.0_r_kind   ! reverse the land/sea mask in bufr
         ! landsea_mask = rmask+.001_r_kind             ! land sea mask

         info%elv = srf_height

         ! 3.6 extract channel bright temperature
   
         tb_inv(1:nchan) = rad%tbb(1:nchan)
         do k = 1, nchan
            if ( tb_inv(k) < tbmin .or. tb_inv(k) > tbmax) &
               tb_inv(k) = missing_r
         end do
         if ( all(tb_inv<0.0) ) then
            num_tovs_local = num_tovs_local -1
            num_tovs_used = num_tovs_used - 1
            nread(inst) = nread(inst) - 1
            cycle
         end if

         !  4.0   assign information to sequential radiance structure
         !--------------------------------------------------------------------------
         allocate (p % tb_inv (1:nchan))
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = landsea_mask
         p%scanpos          = ifov
         p%satzen           = satzen
         p%satazi           = satazi
         p%solzen           = solzen
         p%tb_inv(1:nchan)  = tb_inv(1:nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat
!RTTOV9_3
         p%solazi           = solazi
 !end of RTTOV9_3
         allocate (p%next)   ! add next data

         p => p%next
         nullify (p%next)
      end do
   end do obs

   call closbf(lnbufr)
   close(lnbufr)
1000  continue
end do bufrfile

   if (thinning .and. num_tovs_global > 0 ) then


      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do
      end do

      allocate ( in  (j) )
      allocate ( out (j) )

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i) 
            end do
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time 
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i)  = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_tovs_used_tmp = num_tovs_used
      do j = 1, num_tovs_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat 
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            endif
         end do 
        
         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            endif
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_tovs_thinned = num_tovs_thinned + 1
            num_tovs_used = num_tovs_used - 1
            nread(n) = nread(n) - 1
            continue
         endif

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         endif

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         endif
        
      end do

   endif  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel   + num_tovs_used
   iv%total_rad_channel = iv%total_rad_channel + num_tovs_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_tovs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_tovs_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1) 
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i) 
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_fy3.inc",575,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'num_tovs_file num_tovs_global num_tovs_local num_tovs_used num_tovs_thinned'
   write(unit=stdout,fmt='(5i10)') num_tovs_file,num_tovs_global, num_tovs_local,num_tovs_used,num_tovs_thinned

   deallocate(tb_inv)  

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   
   do i = 1, iv%num_inst
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad

      call da_allocate_rad_iv(i,nchan,iv)

   end do
   
   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0 
   p => head
   ! do while ( associated(p) )

   do n = 1, num_tovs_used
      i = p%sensor_index
      nread(i) = nread(i) + 1

      call da_initialize_rad_iv (i, nread(i), iv, p)

      current => p
      p => p%next

      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do

   deallocate ( p )

   deallocate (nread)
   deallocate (ptotal)

   ! check if sequential structure has been freed
   !
   ! p => head
   ! do i = 1, num_rad_selected
   !    write (unit=stdout,fmt=*)  i, p%tb_inv(1:nchan)
   !    p => p%next
   ! end do

   if (trace_use) call da_trace_exit("da_read_obs_fy3")

  

end subroutine da_read_obs_fy3


subroutine da_read_obs_bufratms (obstype,iv, infile)

   !---------------------------------------------------------------------------
   !  Purpose: read in NCEP bufr atms 1b data to innovation structure
   !
   !   METHOD: use F90 sequential data structure to avoid reading file twice  
   !            so that da_scan_bufrtovs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !               and deallocate sequential data structure
   !   Peiming Dong Added NPP atms, 2012/4/17
   !   Peiming Dong seperated the atms from da_read_obs_bufrtovs.inc in that
   !                all the data needs to be read in together first to make 
   !                the spatial average, 2013/1/10
   !---------------------------------------------------------------------------

   implicit none

   character(5)      ,  intent (in)    :: obstype
   character(20)     ,  intent (in)    :: infile
   type (iv_type)    ,  intent (inout) :: iv


   integer          :: iost
   integer(i_kind), allocatable :: nread(:)
!Dongpm for the spatial average
    integer(i_kind) ,parameter :: Num_Obs = 800000
    integer(i_kind) ,parameter :: NChanl  = 22
    integer(i_kind) ,allocatable :: Fov_save(:)
    real(r_kind)    ,allocatable :: Time_save(:)
    real(r_kind)    ,allocatable :: BT_InOut_save(:,:)
    integer(i_kind) ,allocatable :: Scanline_save(:)
    integer(i_kind)  :: Error_Status
    integer(i_kind)  :: nnum, nn
    real(r_kind)    ,allocatable :: lat_save(:)
    real(r_kind)    ,allocatable :: lon_save(:)
    real(r_kind)    ,allocatable :: obs_time_save(:)
    real(r_kind)    ,allocatable :: satzen_save(:)
    real(r_kind)    ,allocatable :: satazi_save(:)
    real(r_kind)    ,allocatable :: solzen_save(:)
    real(r_kind)    ,allocatable :: solazi_save(:)
    real(r_kind)    ,allocatable :: srf_height_save(:)
    integer(i_kind) ,allocatable :: landsea_mask_save(:)
    integer(i_kind) ,allocatable :: satid_save(:)
    character(len=19),allocatable :: date_char_save(:)
!Dongpm
   integer(i_kind),parameter:: n1bhdr=15
   integer(i_kind),parameter:: n2bhdr=2
!Dongpm   integer(i_kind),parameter:: n1bhdr=13
   integer(i_kind),parameter:: maxinfo=12
   integer(i_kind),parameter:: maxchanl=100

   logical atms
   logical outside, outside_all, iuse
   integer :: inst

   character(10) date
   character(8) subset,subfgn
   character(80) hdr1b
!Dongpm
   character(80) hdr2b
   integer(i_kind) ihh,i,j,k,ifov,idd,ireadmg,ireadsb
   integer(i_kind) iret,idate,im,iy,nchan
   integer :: num_bufr(7), numbufr, ibufr
   character(20) :: filename

   ! thinning variables
   integer(i_kind) itt,itx,iobs,iout
   real(r_kind) terrain,timedif,crit,dist
   real(r_kind) dlon_earth,dlat_earth

   real(r_kind) tbmin,tbmax, tbbad
   real(r_kind) panglr,rato
   ! real(r_kind) rmask
   real(r_kind) step,start

   real(r_double),dimension(maxinfo+maxchanl):: data1b8
   real(r_double),dimension(n1bhdr):: bfr1bhdr
!Dongpm
   real(r_double),dimension(n2bhdr):: bfr2bhdr

   ! Instrument triplet, follow the convension of RTTOV 
   integer   :: platform_id, satellite_id, sensor_id

   ! pixel information
   integer   ::  year,month,day,hour,minute,second  ! observation time
   real*8    ::  obs_time
   ! real      ::  rlat, rlon                         !  lat/lon in degrees   for Anfovs
   real      ::  satzen, satazi, solzen ,solazi       !  scan angles for Anfovs
   integer   ::  landsea_mask
   real      ::  srf_height
   ! channels' bright temperature
   real , allocatable ::   tb_inv(:)                    !  bright temperatures
   !  end type bright_temperature

   type (datalink_type), pointer    :: head, p, current, prev

   integer                        ::  ifgat
   type(info_type)                ::  info
   type(model_loc_type)           ::  loc

!Dongpm
   data hdr1b /'SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA HMSL LSQL SOLAZI'/
   data hdr2b /'CLATH CLONH'/
   !  data hdr1b /'FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA HOLS LSQL SLNM BEARAZ'/

   data tbmin,tbmax,tbbad / 50.0_r_kind, 550.0_r_kind, -9.99e11_r_kind /
   integer :: num_tovs_local, num_tovs_file, num_tovs_global, num_tovs_selected
   integer :: num_tovs_thinned, num_tovs_used, num_tovs_used_tmp
   integer :: lnbufr
   integer :: n
   integer(i_kind), allocatable :: ptotal(:)
   real , allocatable :: in(:), out(:)
   logical :: found, head_found

   if (trace_use) call da_trace_entry("da_read_obs_bufratms")

   ! Initialize variables

!Dongpm
!Dongpm   nchan = 20
   nchan = NChanl
   allocate(nread(1:rtminit_nsensor))
   allocate(ptotal(0:num_fgat_time))
   nread(1:rtminit_nsensor) = 0
   ptotal(0:num_fgat_time) = 0

   ! Set various variables depending on type of data to be read

   ! platform_id  = 1                 !! for NOAA series
   ! platform_id  = 10                !! for METOP series

   atms=      obstype == 'atms '

   if (atms) then
      sensor_id    =  19
      step   = 1.11_r_kind
      start  = -52.725_r_kind
      nchan=22
      subfgn='NC021203'
   end if

   allocate (tb_inv(nchan))

   num_tovs_file     = 0    ! number of obs in file
   num_tovs_global   = 0    ! number of obs within whole domain
   num_tovs_local    = 0    ! number of obs within tile
   num_tovs_thinned  = 0    ! number of obs rejected by thinning
   num_tovs_used     = 0    ! number of obs entered into innovation computation
   num_tovs_selected = 0    ! number of obs limited for debuging
   iobs = 0                 ! for thinning, argument is inout

   ! 0.0  Open unit to satellite bufr file and read file header
   !--------------------------------------------------------------
   allocate(Fov_save(1:Num_obs))
   allocate(Time_save(1:Num_Obs))
   allocate(BT_InOut_save(1:NChanl,1:Num_Obs))
   allocate(Scanline_save(1:Num_Obs))
   allocate(lat_save(1:Num_Obs))
   allocate(lon_save(1:Num_Obs))
   allocate(satid_save(1:Num_Obs))
   allocate(obs_time_save(1:Num_Obs))
   allocate(satzen_save(1:Num_Obs))
   allocate(satazi_save(1:Num_Obs))
   allocate(solzen_save(1:Num_Obs))
   allocate(solazi_save(1:Num_Obs))
   allocate(srf_height_save(1:Num_Obs))
   allocate(landsea_mask_save(1:Num_Obs))
   allocate(date_char_save(1:Num_Obs))
!
   num_bufr(:)=0
   numbufr=0
   nnum=1
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
	    num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if
  
   if (numbufr==0) numbufr=1

bufrfile:  do ibufr=1,numbufr   
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'   
      end if
   end if

!  We want to use specific unit number for bufr data, so we can control the endian format in environment. 
   lnbufr = 99

   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_warning("da_read_obs_bufratms.inc",212, &
         (/"Cannot open file "//infile/))
      if (trace_use) call da_trace_exit("da_read_obs_bufratms")
      return
   end if

   call openbf(lnbufr,'IN',lnbufr)
   call datelen(10)
   call readmg(lnbufr,subset,idate,iret)
   if (subset /= subfgn) then
      call closbf(lnbufr)
      close(lnbufr)
      message(1)='The file title does not match the data subset'
      write(unit=message(2),fmt=*) &
         'infile=', lnbufr, infile,' subset=', subset, ' subfgn=',subfgn
      call da_error("da_read_obs_bufratms.inc",227,message(1:2))
   end if

   iy=0
   im=0
   idd=0
   ihh=0
   write(unit=date,fmt='( i10)') idate
   read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
   write(unit=stdout,fmt=*) &
      'Bufr file date is ',iy,im,idd,ihh,infile

   ! Loop to read bufr file and assign information to a sequential structure
   !-------------------------------------------------------------------------

!   if ( ibufr == 1 ) then
!      allocate (head)
!      !  allocate ( head % tb_inv (1:nchan) )
!      nullify  ( head % next )
!      p => head
!   endif

   if (tovs_start > 1) then
      write (unit=stdout,fmt='(A,I6)') "   Skipping tovs obs before", tovs_start
   end if
   bufrobs: do while (ireadmg(lnbufr,subset,idate)==0 .and. subset==subfgn .and. nnum <  Num_Obs)
      do while (ireadsb(lnbufr)==0 .and. nnum <  Num_Obs)

         ! 1.0     Read header record and data record

         call ufbint(lnbufr,bfr1bhdr,n1bhdr,1,iret,hdr1b)
         call ufbint(lnbufr,bfr2bhdr,n2bhdr,1,iret,hdr2b)
         ! Dongpm call ufbrep(lnbufr,data1b8,1,nchan,iret,'TMBR')
         call ufbrep(lnbufr,data1b8,1,nchan,iret,'TMANT')
         ! call ufbrep(lnbufr,data1b8,1,1,iret,'BEARAZ')

         ! check if observation outside range
         ! 1.5     To save the data

         if(abs(bfr2bhdr(1)) <= 90. .and. abs(bfr2bhdr(2)) <= 360.)then
              lat_save(nnum) = bfr2bhdr(1)
              lon_save(nnum) = bfr2bhdr(2)
         elseif(abs(bfr1bhdr(9)) <= 90. .and. abs(bfr1bhdr(10)) <= 360.)then
              lat_save(nnum) = bfr1bhdr(9)
              lon_save(nnum) = bfr1bhdr(10)
         endif
         ifov = nint(bfr1bhdr(bufr_ifov))
         Fov_save(nnum) = ifov
         satid_save(nnum)=nint(bfr1bhdr(bufr_satellite_id))         
         year   = bfr1bhdr(bufr_year)
         month  = bfr1bhdr(bufr_month)
         day    = bfr1bhdr(bufr_day)
         hour   = bfr1bhdr(bufr_hour)
         minute = bfr1bhdr(bufr_minute)
         second = bfr1bhdr(bufr_second)

         write(unit=date_char_save(nnum), fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
            year, '-', month, '-', day, '_', hour, ':', minute, ':', second

         !  QC3: time consistency check with the background date

         if (year <= 99) then
            if (year < 78) then
               year = year + 2000
            else
               year = year + 1900
            end if
         end if

         call da_get_julian_time(year,month,day,hour,minute,obs_time)
         obs_time_save(nnum)=obs_time
         Time_save(nnum)=obs_time_save(nnum)*60.0+second         

         panglr=(start+float(ifov-1)*step)*deg2rad
         satzen_save(nnum) = bfr1bhdr(bufr_satzen) !*deg2rad   ! local zenith angle
         satazi_save(nnum) = panglr*rad2deg            ! look angle
         solzen_save(nnum) = bfr1bhdr(bufr_solzen)              ! solar zenith angle
         solazi_save(nnum) = bfr1bhdr(bufr_solazi)              !RTTOV9_3
         srf_height_save(nnum) = bfr1bhdr(bufr_station_height)          ! station height
         landsea_mask_save(nnum) = nint(bfr1bhdr(bufr_landsea_mask))  ! 0:land ; 1:sea (same as RTTOV)
         BT_InOut_save(1:nchan,nnum)= data1b8(1:nchan)
!
         nnum=nnum+1
         num_tovs_file = num_tovs_file + 1

      end do
   end do bufrobs

!
         call closbf(lnbufr)
         close(lnbufr)

 end do bufrfile
!
         nnum=nnum-1
         if(nnum <= 0) then
            call da_warning("da_read_obs_bufratms.inc",323,(/"No ATMS data were read in"/))
            if (trace_use) call da_trace_exit("da_read_obs_bufratms")
            return
         endif
         write(unit=message(1),fmt='(a,i10)') 'The number of observations is:',nnum-1
         call da_message(message(1:1))
      call ATMS_Spatial_Average(nnum, NChanl, Fov_save(1:nnum), Time_save(1:nnum), BT_InOut_save(1:NChanl,1:nnum), &
                                Scanline_save, Error_Status)
      if(Error_Status==1) then
         WRITE(UNIT=message(1), FMT='(A)')"Error reading ATMS data"
         call da_error("da_read_obs_bufratms.inc",333,message(1:1))
      endif         

 obs: do nn=1, nnum
   if ( nn == 1 ) then
      allocate (head)
!      !  allocate ( head % tb_inv (1:nchan) )
      nullify  ( head % next )
      p => head
   endif

         ! 2.0     Extract observation location and other required information
         !     QC1:  judge if data is in the domain, read next record if not
         !------------------------------------------------------------------------
         ! rlat = bfr1bhdr(bufr_lat)
         ! rlon = bfr1bhdr(bufr_lat)
         ! if (rlon < 0.0) rlon = rlon+360.0
              info%lat = lat_save(nn)
              info%lon = lon_save(nn)

         call da_llxy (info, loc, outside, outside_all)

         if (outside_all) cycle

         !  3.0     Extract other information
         !------------------------------------------------------
         !  3.1     Extract satellite id and scan position. 
   
         if ( satid_save(nn) == 224 ) then ! npp
            platform_id = 17
            satellite_id = 0
         end if
         ifov = Fov_save(nn) 

         !  QC2:  limb pixel rejected (not implemented)

         !  3.2     Extract date information.
    
         info%date_char=date_char_save(nn)
         obs_time=obs_time_save(nn)

         if (obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time)) cycle

         ! 3.2.1 determine FGAT index ifgat
   
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do

         ! 3.3 Find wrfvar instrument index from RTTOV instrument triplet
         !     go to next data if id is not in the lists

         inst = 0
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do
         if (inst == 0) cycle

         ! 3.4 extract satellite and solar angle
   
            satzen = satzen_save(nn)
            satzen = abs(satzen)
            ! if (amsua .and. ifov .le. 15) satzen=-satzen
            ! if (amsub .and. ifov .le. 45) satzen=-satzen
            ! if (hirs3 .and. ifov .le. 28) satzen=-satzen
         if ( satzen > 65.0 ) cycle   ! 1 has a satzen > 65.0 check
         satazi = satazi_save(nn)            ! look angle
         ! if (satazi<0.0) satazi = satazi+360.0
         solzen = solzen_save(nn)              ! solar zenith angle
         solazi = solazi_save(nn)              !RTTOV9_3

         num_tovs_global = num_tovs_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1

         if (num_tovs_global < tovs_start) then
            cycle
         end if

         if (num_tovs_global > tovs_end) then
            write (unit=stdout,fmt='(A,I6)') "   Skipping radiance obs after", tovs_end
            exit obs
         end if

         num_tovs_selected = num_tovs_selected + 1
 
         if (num_tovs_selected > max_tovs_input) then
            write(unit=message(1),fmt='(A,I10,A)') &
               "Max number of tovs",max_tovs_input," reached"
            call da_warning("da_read_obs_bufratms.inc",428,message(1:1))
            num_tovs_selected = num_tovs_selected-1
            num_tovs_global   = num_tovs_global-1
            ptotal(ifgat) = ptotal(ifgat) - 1
            exit obs
         end if

         if (outside) cycle ! No good for this PE
         num_tovs_local = num_tovs_local + 1

         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad           
            timedif = 0.0 !2.0_r_kind*abs(tdiff)        ! range:  0 to 6
            terrain = 0.01_r_kind*abs(bfr1bhdr(13))
            crit = 1.0 !0.01_r_kind+terrain + timedif !+ 10.0_r_kind*float(iskip)
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_tovs_thinned=num_tovs_thinned+1
               cycle
            end if
         end if

        
         num_tovs_used = num_tovs_used + 1
         nread(inst) = nread(inst) + 1

         ! 3.5 extract surface information
         srf_height = srf_height_save(nn)          ! station height
         if (srf_height < 8888.0 .AND. srf_height > -416.0) then
         else
            srf_height = 0.0
         endif  

         landsea_mask = landsea_mask_save(nn)  ! 0:land ; 1:sea (same as RTTOV)
!Dongpm  There is no landsea-mask in atms bufr data
         if (landsea_mask <= 1 .AND. landsea_mask >= 0) then
         else
            landsea_mask = 0
         endif
         
         ! rmask=one                          ! land
         ! if (nint(bfr1bhdr(bufr_landsea_mask))==1) rmask=0.0_r_kind   ! reverse the land/sea mask in bufr
         ! landsea_mask = rmask+.001_r_kind             ! land sea mask

         info%elv = srf_height

         ! 3.6 extract channel bright temperature
   
         tb_inv(1:nchan) = BT_InOut_save(1:nchan,nn)
         do k = 1, nchan
            if ( tb_inv(k) < tbmin .or. tb_inv(k) > tbmax) &
               tb_inv(k) = missing_r
         end do
         if ( all(tb_inv<0.0) ) then
            num_tovs_local = num_tovs_local -1
            num_tovs_used = num_tovs_used - 1
            nread(inst) = nread(inst) - 1
            cycle
         end if
         !  4.0   assign information to sequential radiance structure
         !--------------------------------------------------------------------------
         allocate (p % tb_inv (1:nchan))
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = landsea_mask
         p%scanpos          = ifov
         p%satzen           = satzen
         p%satazi           = satazi
         p%solzen           = solzen
         p%tb_inv(1:nchan)  = tb_inv(1:nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat
!RTTOV9_3
         p%solazi           = solazi
 !end of RTTOV9_3
         allocate (p%next)   ! add next data

         p => p%next
         nullify (p%next)
!      end do
   end do obs

!   call closbf(lnbufr)
!   close(lnbufr)

!end do bufrfile

   if (thinning .and. num_tovs_global > 0 ) then


      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do
      end do

      allocate ( in  (j) )
      allocate ( out (j) )

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i) 
            end do
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time 
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i)  = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_tovs_used_tmp = num_tovs_used
      do j = 1, num_tovs_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat 
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            endif
         end do 
        
         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            endif
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_tovs_thinned = num_tovs_thinned + 1
            num_tovs_used = num_tovs_used - 1
            nread(n) = nread(n) - 1
            continue
         endif

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         endif

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         endif
        
      end do

   endif  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel   + num_tovs_used
   iv%total_rad_channel = iv%total_rad_channel + num_tovs_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_tovs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_tovs_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1) 
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i) 
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_bufratms.inc",631,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'num_tovs_file num_tovs_global num_tovs_local num_tovs_used num_tovs_thinned'
   write(unit=stdout,fmt='(5i10)') num_tovs_file,num_tovs_global, num_tovs_local,num_tovs_used,num_tovs_thinned

   deallocate(tb_inv)  

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   
   do i = 1, iv%num_inst
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad

      call da_allocate_rad_iv(i,nchan,iv)

   end do
   
   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0 
   p => head
   ! do while ( associated(p) )

   do n = 1, num_tovs_used
      i = p%sensor_index
      nread(i) = nread(i) + 1

      call da_initialize_rad_iv (i, nread(i), iv, p)

      current => p
      p => p%next

      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do
!  deallocate the save data
   deallocate(Fov_save)
   deallocate(Time_save)
   deallocate(BT_InOut_save)
   deallocate(Scanline_save)
   deallocate(lat_save)
   deallocate(lon_save)
   deallocate(satid_save)
   deallocate(obs_time_save)
   deallocate(satzen_save)
   deallocate(satazi_save)
   deallocate(solzen_save)
   deallocate(solazi_save)
   deallocate(srf_height_save)
   deallocate(landsea_mask_save)
   deallocate(date_char_save)
   deallocate ( p )

   deallocate (nread)
   deallocate (ptotal)

   ! check if sequential structure has been freed
   !
   ! p => head
   ! do i = 1, num_rad_selected
   !    write (unit=stdout,fmt=*)  i, p%tb_inv(1:nchan)
   !    p => p%next
   ! end do

   if (trace_use) call da_trace_exit("da_read_obs_bufratms")
  

end subroutine da_read_obs_bufratms


!
  SUBROUTINE ATMS_Spatial_Average(Num_Obs, NChanl, FOV, Time, BT_InOut, &
       Scanline, Error_Status)


    IMPLICIT NONE
    
    ! Declare passed variables
    integer(i_kind) ,intent(in   ) :: Num_Obs, NChanl
    integer(i_kind) ,intent(in   ) :: Fov(num_obs)
    real(r_kind)    ,intent(in   ) :: Time(Num_Obs)
    real(r_kind)    ,intent(inout) :: BT_InOut(NChanl,Num_Obs)
    integer(i_kind) ,intent(  out) :: Scanline(Num_Obs)
    integer(i_kind) ,intent(  out) :: Error_Status
    ! Declare local parameters
    integer(i_kind), parameter :: atms1c_h_wmosatid=224
    integer(i_kind), parameter :: lninfile=15
    integer(i_kind), parameter :: max_fov=96
    integer(i_kind), parameter :: max_obs=1000000
    real(r_kind), parameter    :: scan_interval = 8.0_r_kind/3.0_r_kind
    ! Maximum number of channels 
    integer(i_kind), parameter :: MaxChans = 22
    ! Minimum allowed BT as a function of channel number
    real(r_kind), parameter :: MinBT(MaxChans) = &
         (/ 120.0_r_kind, 120.0_r_kind, 190.0_r_kind, 190.0_r_kind, &
            200.0_r_kind, 200.0_r_kind, 200.0_r_kind, 190.0_r_kind, &
            190.0_r_kind, 180.0_r_kind, 180.0_r_kind, 180.0_r_kind, &
            190.0_r_kind, 200.0_r_kind, 200.0_r_kind, 120.0_r_kind, &
            120.0_r_kind, 120.0_r_kind, 150.0_r_kind, 170.0_r_kind, &
            180.0_r_kind, 190.0_r_kind /)
    ! Maximum allowed BT as a function of channel number
    real(r_kind), parameter :: MaxBT(MaxChans) = &
         (/ 320.0_r_kind, 320.0_r_kind, 300.0_r_kind, 300.0_r_kind, &
            300.0_r_kind, 270.0_r_kind, 250.0_r_kind, 240.0_r_kind, &
            240.0_r_kind, 250.0_r_kind, 250.0_r_kind, 270.0_r_kind, &
            280.0_r_kind, 290.0_r_kind, 300.0_r_kind, 320.0_r_kind, &
            320.0_r_kind, 300.0_r_kind, 300.0_r_kind, 300.0_r_kind, &
            300.0_r_kind, 300.0_r_kind /)
    
    ! Declare local variables
    character(30) :: Cline
    integer(i_kind) :: i, iscan, ifov, ichan, nchannels, wmosatid, version
    integer(i_kind) :: ios, max_scan, mintime
    integer(i_kind) :: nxaverage(nchanl), nyaverage(nchanl)
    integer(i_Kind) :: channelnumber(nchanl),qc_dist(nchanl)
    integer(i_kind), ALLOCATABLE ::  scanline_back(:,:)
    real(r_kind) :: sampling_dist, beamwidth(nchanl) 
    real(r_kind) :: newwidth(nchanl), cutoff(nchanl)
    real(r_kind), allocatable, target :: bt_image(:,:,:)
    real(r_kind), pointer :: bt_image1(:,:)
    Error_Status=0
    IF (NChanl > MaxChans) THEN
       WRITE(0,*) 'Unexpected number of ATMS channels: ',nchanl
       Error_Status = 1
       RETURN
    END IF
    ! Read the beamwidth requirements
    OPEN(lninfile,file='radiance_info/atms_beamwidth.txt',form='formatted',status='old', &
         iostat=ios)
    IF (ios /= 0) THEN
       WRITE(*,*) 'Unable to open atms_beamwidth.txt'
       Error_Status=1
       RETURN
    ENDIF
    wmosatid=999
    read(lninfile,'(a30)',iostat=ios) Cline
    DO WHILE (wmosatid /= atms1c_h_wmosatid .AND. ios == 0)
       DO WHILE (Cline(1:1) == '#')
          read(lninfile,'(a30)') Cline
       ENDDO
       READ(Cline,*) wmosatid
       
       read(lninfile,'(a30)') Cline
       DO WHILE (Cline(1:1) == '#')
          read(lninfile,'(a30)') Cline
       ENDDO
       READ(Cline,*) version
       
       read(lninfile,'(a30)') Cline
       DO WHILE (Cline(1:1) == '#')
          read(lninfile,'(a30)') Cline
       ENDDO
       READ(Cline,*) sampling_dist
       
       read(lninfile,'(a30)') Cline
       DO WHILE (Cline(1:1) == '#')
          read(lninfile,'(a30)') Cline
       ENDDO
       READ(Cline,*) nchannels
      
       read(lninfile,'(a30)') Cline
       if (nchannels > 0) then 
          DO ichan=1,nchannels
             read(lninfile,'(a30)') Cline
             DO WHILE (Cline(1:1) == '#')
                read(lninfile,'(a30)') Cline
             ENDDO
             READ(Cline,*) channelnumber(ichan),beamwidth(ichan), &
                  newwidth(ichan),cutoff(ichan),nxaverage(ichan), &
                  nyaverage(ichan), qc_dist(ichan)
          ENDDO
       end if
       read(lninfile,'(a30)',iostat=ios) Cline
    ENDDO
    IF (wmosatid /= atms1c_h_wmosatid) THEN
       WRITE(*,*) 'ATMS_Spatial_Averaging: sat id not matched in atms_beamwidth.dat'
       Error_Status=1
       RETURN
    ENDIF
    CLOSE(lninfile)
  
    ! Determine scanline from time
    MinTime = MINVAL(Time)
    Scanline(:)   = NINT((Time(1:Num_Obs)-MinTime)/Scan_Interval)+1
    Max_Scan=MAXVAL(Scanline)
    write(*,*) 'Max_Scan:',Max_Scan
    
    ALLOCATE(BT_Image(Max_FOV,Max_Scan,nchanl))
    ALLOCATE(Scanline_Back(Max_FOV,Max_Scan))
    BT_Image(:,:,:) = 1000.0_r_kind
    
    ScanLine_Back(:,:) = -1
    DO I=1,Num_Obs
       bt_image(FOV(I),Scanline(I),:)=bt_inout(:,I)
       Scanline_Back(FOV(I),Scanline(I))=I
    END DO
    DO IChan=1,nchanl
    
   
       bt_image1 => bt_image(:,:,ichan)
       ! Set all scan positions to missing in a scanline if one is missing
       do iscan=1,max_scan
          if (ANY(bt_image1(:,iscan) > 500.0_r_kind)) &
	     bt_image1(:,iscan)=1000.0_r_kind
       enddo
       ! If the channel number is present in the channelnumber array we should process it 
       ! (otherwise bt_inout just keeps the same value):
       IF (ANY(channelnumber(1:nchannels) == ichan)) THEN
          CALL MODIFY_BEAMWIDTH ( max_fov, max_scan, bt_image1, &
               sampling_dist, beamwidth(ichan), newwidth(ichan), &
               cutoff(ichan), nxaverage(ichan), nyaverage(ichan), &
               qc_dist(ichan), MinBT(Ichan), MaxBT(IChan), IOS)
          
          IF (IOS == 0) THEN
             do iscan=1,max_scan
                do ifov=1,max_fov
                   IF (Scanline_Back(IFov, IScan) > 0) &
                        bt_inout(channelnumber(ichan),Scanline_Back(IFov, IScan)) = &
                        BT_Image1(ifov,iscan)
                end do
             end do
          ELSE
             Error_Status=1
             RETURN
          END IF
       END IF
    END DO
    DEALLOCATE(BT_Image, Scanline_Back)
    NULLIFY(BT_Image1)
    
END Subroutine ATMS_Spatial_Average

SUBROUTINE MODIFY_BEAMWIDTH ( nx, ny, image, sampling_dist,& 
     beamwidth, newwidth, mtfcutoff, nxaverage, nyaverage, qc_dist, &
     Minval, MaxVal, Error)
     
!-----------------------------------------
! Name: $Id: modify_beamwidth.F 222 2010-08-11 14:39:09Z frna $
!
! Purpose:
!   Manipulate the effective beam width of an image. For example, convert ATMS
!   to AMSU-A-like resolution while reducing the noise.
!
! Method:
!   1) Pad the image to a power of 2 in each dimension.
! If FFT technique is to be used then: 
!   2) Assuming Gaussian beam shapes, calcluate the input and output Modulation
!      Transfer Functions (MTF).
!   3) FFT image to frequency domain (2-D).
!   4) Multiply by output MTF divided by input MTF. If a cut-off is specified
!      (when attempting to make the beam width narrower), attenuate further
!      by an exponential function - factor of 2 at the cutoff. 
!   5) FFT back to image domain 
! Finally,
!   6) Over-write the input image, with averaging if requested.
!
! COPYRIGHT
!    This software was developed within the context of the EUMETSAT Satellite
!    Application Facility on Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 1 December 2006, between EUMETSAT and the
!    Met Office, UK, by one or more partners within the NWP SAF. The partners
!    in the NWP SAF are the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! History:
! Version    Date     Comment
!
!  1.0   22/07/2010   N.C.Atkinson
!  1.1   21/11/2011   Convert to f90. A. Collard
!
! Code Description:
!   FORTRAN 77, following AAPP standards
!
! Declarations:
      IMPLICIT NONE
! Parameters
      real(r_kind), parameter    :: Missing_Value=-888888.0
      INTEGER(I_KIND), PARAMETER :: nxmax=128  !Max number of spots per scan line
      INTEGER(I_KIND), PARAMETER :: nymax=8192 !Max number of lines. Allows 6hrs of ATMS.
! Arguments
      INTEGER(I_KIND), INTENT(IN)  :: nx, ny         !Size of image
      REAL(R_KIND), INTENT(INOUT)  :: image(nx,ny)   !BT or radiance image
      REAL(R_KIND), INTENT(IN)     :: sampling_dist  !typically degrees
      REAL(R_KIND), INTENT(IN)     :: beamwidth      !ditto
      REAL(R_KIND), INTENT(IN)     :: newwidth       !ditto
      REAL(R_KIND), INTENT(IN)     :: mtfcutoff      !0.0 to 1.0
      INTEGER(I_KIND), INTENT(IN)  :: nxaverage      !Number of samples to average (or zero)
      INTEGER(I_KIND), INTENT(IN)  :: nyaverage      !Number of samples to average (or zero)
      INTEGER(I_KIND), INTENT(IN)  :: qc_dist        !Number of samples around missing data to set to 
      REAL(R_KIND), INTENT(IN)     :: maxval         !BTs above this are considered missing
      REAL(R_KIND), INTENT(IN)     :: minval         !BTs below this are considered missing
      INTEGER(I_KIND), INTENT(OUT) :: Error          !Error Status
       
! Local variables
      INTEGER(I_KIND) :: nxpad, nypad, dx, dy
      INTEGER(I_KIND) :: i,j,k,ix,iy, ii, jj
      INTEGER(I_KIND) :: ifirst
      INTEGER(I_KIND) :: xpow2, ypow2
      INTEGER(I_KIND) :: nxav2, nyav2, naverage
      INTEGER(I_KIND) :: deltax, minii, maxii, minjj, maxjj
      REAL(R_KIND), ALLOCATABLE :: mtfxin(:),mtfxout(:)
      REAL(R_KIND), ALLOCATABLE :: mtfyin(:),mtfyout(:)
      REAL(R_KIND) :: mtfin,mtfout,mtf_constant
      REAL(R_KIND), ALLOCATABLE :: mtfpad(:,:)
      REAL(R_KIND), ALLOCATABLE :: imagepad(:,:)
      REAL(R_KIND), ALLOCATABLE :: work(:)
      REAL(R_KIND) :: f,df,factor
      REAL(R_KIND) :: PI, LN2, LNcsquared
      LOGICAL :: missing
      LOGICAL, ALLOCATABLE :: gooddata_map(:,:)
! End of declarations
!-----------------------------------------
      
      PI = 4.0_r_kind*atan(1.0)
      LN2 = LOG(2.0_r_kind)
      MTF_Constant=-(PI/(2*sampling_dist))**2/LN2
      IF (mtfcutoff > 0.0_r_kind) LNcsquared = LOG(mtfcutoff)**2
      nxav2 = nxaverage/2
      nyav2 = nyaverage/2
      naverage = nxaverage*nyaverage
      Error = 0
!1) Pad the image up to the nearest power of 2 in each dimension, by reversing
!the points near the edge.
      xpow2 = INT(LOG(nx*1.0_r_kind)/LN2 + 1.0_r_kind)
      ypow2 = INT(LOG(ny*1.0_r_kind)/LN2 + 1.0_r_kind)
      nxpad = 2**xpow2
      nypad = 2**ypow2
      dx = (nxpad - nx)/2
      dy = (nypad - ny)/2
      IF (nxpad > nxmax) THEN
         write(*,*) 'ATMS_Spatial_Average: nx too large, maximum allowed value is ',nxmax-1
         Error = 1
         RETURN
      END IF
      
      IF (nypad > nymax) THEN
         write(*,*) 'ATMS_Spatial_Average: ny too large, maximum allowed value is ',nymax-1
         Error = 1
         RETURN
      END IF
      ALLOCATE(mtfxin(nxpad),mtfxout(nxpad))
      ALLOCATE(mtfyin(nypad),mtfyout(nypad))
      ALLOCATE(mtfpad(nxpad,nypad))
      ALLOCATE(imagepad(nxpad,nypad))
      ALLOCATE(work(nypad))
      ALLOCATE(gooddata_map(nxmax,nymax))
!Loop over scan positions
      DO j=dy+1,dy+ny
        DO i=dx+1,dx+nx
          if (image(i-dx,j-dy) < minval) &
               image(i-dx,j-dy) = minval - 1.0_r_kind
          if (image(i-dx,j-dy) > maxval ) &
               image(i-dx,j-dy) = maxval + 1.0_r_kind
          imagepad(i,j) = image(i-dx,j-dy)   !Take a copy of the input data
          gooddata_map(i,j) = .TRUE.   ! Initialised for step 6)
        ENDDO
!Interpolate missing points in the along-track direction
        ifirst = -1
        missing = .false.
        
        DO i=dx+1,dx+nx
          IF (.not.missing) THEN
            IF (imagepad(i,j) >= minval .AND. imagepad(i,j) <= maxval) THEN
              ifirst = i
            ELSE
              missing = .true.
            ENDIF
          ELSE
            IF (imagepad(i,j) >= minval .AND. imagepad(i,j) <= maxval) THEN  !First good point 
                                                                             ! after missing
               missing = .false.
               IF (ifirst == -1) THEN
                  DO k=dx+1,i-1
                     imagepad(k,j) = imagepad(i,j)      !Constant
                  ENDDO
               ELSE
                  DO k=ifirst+1,i-1
                     factor = (i-k)*1.0_r_kind/(i-ifirst)      !Interpolate
                     imagepad(k,j) = imagepad(ifirst,j)*factor + &
                          imagepad(i,j)*(1.0_r_kind-factor)
                  ENDDO
               ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF (missing) THEN         !Last scan is missing
          IF (ifirst >= 1) then
            DO k=ifirst+1,dx+nx
              imagepad(k,j) = imagepad(ifirst,j)     !Constant
            ENDDO
          ENDIF
        ENDIF          
!Continue padding the edges
        DO i=1,dx
          imagepad(i,j) = imagepad(dx+dx+2-i,j)
        ENDDO
        DO i=nx+dx+1,nxpad
          imagepad(i,j) = imagepad(nx+dx+nx+dx-i,j)
        ENDDO
     ENDDO
     DO j=1,dy
        DO i=1,nxpad
           imagepad(i,j) = imagepad(i,dy+dy+2-j)
        ENDDO
     ENDDO
     DO j=ny+dy+1,nypad
        DO i=1,nxpad
           imagepad(i,j) = imagepad(i,ny+dy+ny+dy-j)
        ENDDO
     ENDDO
!2) Compute the MTF modifications. Assume beams are Gaussian.
      IF (newwidth > 0) THEN
        df = 1.0_r_kind/nxpad
        DO i=1,nxpad/2+1
          f = df*(i-1)      !DC to Nyquist
          mtfxin(i) = exp(MTF_Constant*(f*beamwidth)**2)
          mtfxout(i) = exp(MTF_Constant*(f*newwidth)**2)
          IF (i > 1 .AND. i < nxpad/2+1) THEN
            mtfxin(nxpad-i+2) = mtfxin(i)
            mtfxout(nxpad-i+2) = mtfxout(i)
          ENDIF
        ENDDO
        df = 1.0_r_kind/nypad
        DO i=1,nypad/2+1
          f = df*(i-1)      !DC to Nyquist
          mtfyin(i) = exp(MTF_Constant*(f*beamwidth)**2)
          mtfyout(i) = exp(MTF_Constant*(f*newwidth)**2)
          IF (i > 1 .AND. i < nypad/2+1) THEN
            mtfyin(nypad-i+2) = mtfyin(i)
            mtfyout(nypad-i+2) = mtfyout(i)
          ENDIF
        ENDDO
        DO i=1,nxpad
          DO j=1,nypad
            mtfin = mtfxin(i)*mtfyin(j)
            mtfout = mtfxout(i)*mtfyout(j)
            if (mtfcutoff > 0.0_r_kind) THEN
              mtfpad(i,j) = (mtfout * &
                exp(-LN2/LNcsquared*(LOG(mtfout))**2))/mtfin
            else
              mtfpad(i,j) = mtfout/mtfin
            endif
          ENDDO
        ENDDO
!3) Fourier transform, line by line then column by column.
!After each FFT, points 1 to nxpad/2+1 contain the real part of the spectrum,
!the rest contain the imaginary part in reverse order.
        DO j=1,nypad
           CALL SFFTCF(imagepad(:,j),nxpad,xpow2)
        ENDDO
        DO i=1,nxpad
           DO j=1,nypad
              work(j) = imagepad(i,j)
           ENDDO
           CALL SFFTCF(work,nypad,ypow2)
           DO j=1,nypad
              imagepad(i,j) = work(j)
           ENDDO
        ENDDO
!4) Multiply the spectrum by the MTF factor
        DO j=1,nypad
           DO i=1,nxpad
            imagepad(i,j) = imagepad(i,j)*mtfpad(i,j)
          ENDDO
        ENDDO
!5) Inverse Fourier transform, column by column then line by line 
        DO i=1,nxpad
          DO j=1,nypad
            work(j) = imagepad(i,j)
          ENDDO
          CALL SFFTCB(work,nypad,ypow2)
          DO j=1,nypad
            imagepad(i,j) = work(j)
          ENDDO
        ENDDO
        DO j=1,nypad
          CALL SFFTCB(imagepad(:,j),nxpad,xpow2)
        ENDDO
     ENDIF   !New width is specified
!6) Reset missing values in gooddata_map, based on qc_dist and the values 
!   in the input image array
     ! Set the ends of the image to missing in the along track direction
     ! (doing the same across track will remove too much data)
     gooddata_map(:,1:qc_dist)=.FALSE.
     gooddata_map(:,ny-qc_dist+1:ny)=.FALSE.
     
     DO j=1,ny
        DO i=1,nx
           IF (image(i,j) <= minval .OR. image(i,j) >= maxval ) THEN
              minjj=max(j+dy-qc_dist,1)
              maxjj=min(j+dy+qc_dist,nymax)
              DO jj=minjj,maxjj
                 deltax=INT(SQRT(REAL(qc_dist**2 - (jj-j-dy)**2 )))
                 minii=max(i+dx-deltax,1)
                 maxii=min(i+dx+deltax,nxmax)
                 DO ii=minii,maxii
                    gooddata_map(ii,jj)=.FALSE.
                 END DO
              END DO
           END IF
        END DO
     END DO
!7) Over-write the input image (points that are not missing)
     DO j=1,ny
        DO i=1,nx
           IF (gooddata_map(i+dx,j+dy)) THEN
              IF (nxav2 == 0.0_r_kind .AND. nyav2 == 0) THEN
                 image(i,j) = imagepad(i+dx,j+dy)
              ELSE
                 image(i,j) = 0.0_r_kind             !Do averaging
                 DO ix = -nxav2,nxav2
                    DO iy = -nyav2,nyav2
                       image(i,j) = image(i,j) + imagepad(i+dx+ix,j+dy+iy)
                    ENDDO
                 ENDDO
                 image(i,j) = image(i,j)/naverage
              ENDIF
           ELSE
              image(i,j) = missing_value
           END IF
        ENDDO
     ENDDO
!8) Deallocate arrays
     DEALLOCATE(mtfxin,mtfxout)
     DEALLOCATE(mtfyin,mtfyout)
     DEALLOCATE(mtfpad)
     DEALLOCATE(imagepad)
     DEALLOCATE(work)
     DEALLOCATE(gooddata_map)
     RETURN
   END SUBROUTINE MODIFY_BEAMWIDTH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  A real-valued, in place, split-radix FFT program
!  Real input and output in data array X
!  Length is N = 2 ** M
!  Decimation-in-time, cos/sin in second loop
!  Output in order:
!         [ Re(0), Re(1), ..., Re(N/2), Im(N/2-1), ..., Im(1) ]
!
!  This FFT computes
!     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)
!
!
!  H.V. Sorensen, Rice University, Oct. 1985
!
!  Reference:  H.V. Sorensen, D.L. Jones, M.T. Heideman, & C.S. Burrus;
!              Real-Valued Fast Fourier Transform Algorithms; IEEE
!              Trans. Acoust., Speech, Signal Process.; Vol ASSP-35,
!              June 1987, pp. 849-863.
!
!  This code was originally named RVFFT.
!
!  History:
!   21/11/2011 Converted to something resembling f90.   A.Collard
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SFFTCF( X, N, M )
      IMPLICIT NONE
! ... Parameters ...
      REAL(R_KIND), PARAMETER :: SQRT2 = 1.4142135623730950488
      REAL(R_KIND), PARAMETER :: TWOPI = 6.2831853071795864769 
! ... Scalar arguments ...
      INTEGER(I_KIND), INTENT(IN) :: N, M
! ... Array arguments ...
      REAL(R_KIND), INTENT(INOUT) ::  X(N)
! ... Local scalars ...
      INTEGER(I_KIND)  J, I, K, IS, ID, I0, I1, I2, I3, I4, I5, I6, I7, I8
      INTEGER(I_KIND)  N1, N2, N4, N8
      REAL(R_KIND)  XT, R1, T1, T2, T3, T4, T5, T6
      REAL(R_KIND)  A, A3, E, CC1, SS1, CC3, SS3
!
! ... Exe. statements ...
!
      IF ( N == 1 ) RETURN
!
 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
         IF ( I >= J ) GOTO 101
         XT = X(J)
         X(J) = X(I)
         X(I) = XT
 101     K = N / 2
 102     IF ( K >= J ) GOTO 103
            J = J - K
            K = K / 2
            GOTO 102
 103     J = J + K
 104  CONTINUE
! 
      IS = 1
      ID = 4
 70   DO 60, I0 = IS, N, ID
         I1 = I0 + 1
         R1 = X(I0)
         X(I0) = R1 + X(I1)
         X(I1) = R1 - X(I1)
 60   CONTINUE
      IS = 2 * ID - 1
      ID = 4 * ID
      IF ( IS < N ) GOTO 70
!
      N2 = 2
      DO 10, K = 2, M
         N2 = N2 * 2
         N4 = N2 / 4
         N8 = N2 / 8
         E = TWOPI / N2
         IS = 0
         ID = N2 * 2
 40      DO 38, I = IS, N-1, ID
            I1 = I + 1
            I2 = I1 + N4
            I3 = I2 + N4
            I4 = I3 + N4
            T1 = X(I4) + X(I3)
            X(I4) = X(I4) - X(I3)
            X(I3) = X(I1) - T1
            X(I1) = X(I1) + T1
            IF ( N4 == 1 ) GOTO 38
            I1 = I1 + N8
            I2 = I2 + N8
            I3 = I3 + N8
            I4 = I4 + N8
            T1 = ( X(I3) + X(I4) ) / SQRT2
            T2 = ( X(I3) - X(I4) ) / SQRT2
            X(I4) = X(I2) - T1
            X(I3) = - X(I2) - T1
            X(I2) = X(I1) - T2
            X(I1) = X(I1) + T2
 38      CONTINUE
         IS = 2 * ID - N2
         ID = 4 * ID
         IF ( IS < N ) GOTO 40
         A = E
         DO 32, J = 2, N8
            A3 = 3 * A
            CC1 = COS(A)
            SS1 = SIN(A)
            CC3 = COS(A3)
            SS3 = SIN(A3)
            A = J * E
            IS = 0
            ID = 2 * N2
 36         DO 30, I = IS, N-1, ID
               I1 = I + J
               I2 = I1 + N4
               I3 = I2 + N4
               I4 = I3 + N4
               I5 = I + N4 - J + 2
               I6 = I5 + N4
               I7 = I6 + N4
               I8 = I7 + N4
               T1 = X(I3) * CC1 + X(I7) * SS1
               T2 = X(I7) * CC1 - X(I3) * SS1
               T3 = X(I4) * CC3 + X(I8) * SS3
               T4 = X(I8) * CC3 - X(I4) * SS3
               T5 = T1 + T3
               T6 = T2 + T4
               T3 = T1 - T3
               T4 = T2 - T4
               T2 = X(I6) + T6
               X(I3) = T6 - X(I6)
               X(I8) = T2
               T2 = X(I2) - T3
               X(I7) = - X(I2) - T3
               X(I4) = T2
               T1 = X(I1) + T5
               X(I6) = X(I1) - T5
               X(I1) = T1
               T1 = X(I5) + T4
               X(I5) = X(I5) - T4
               X(I2) = T1
 30         CONTINUE
            IS = 2 * ID - N2
            ID = 4 * ID
            IF ( IS < N ) GOTO 36
 32      CONTINUE
 10   CONTINUE
      RETURN
!
! ... End of subroutine SFFTCF ...
!
   END SUBROUTINE SFFTCF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  A real-valued, in place, split-radix IFFT program
!  Hermitian symmetric input and real output in array X
!  Length is N = 2 ** M
!  Decimation-in-frequency, cos/sin in second loop
!  Input order:
!         [ Re(0), Re(1), ..., Re(N/2), Im(N/2-1), ..., Im(1) ]
!
!  This FFT computes
!     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)
!
!
!  H.V. Sorensen, Rice University, Nov. 1985
!
!  Reference:  H.V. Sorensen, D.L. Jones, M.T. Heideman, & C.S. Burrus;
!              Real-Valued Fast Fourier Transform Algorithms; IEEE
!              Trans. Acoust., Speech, Signal Process.; Vol ASSP-35,
!              June 1987, pp. 849-863.
!
!  This code was originally named IRVFFT.
!
!  History:
!   21/11/2011 Converted to something resembling f90.   A.Collard
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SFFTCB( X, N, M )
      IMPLICIT NONE
! ... Parameters ...
      REAL(R_KIND), PARAMETER :: SQRT2 = 1.4142135623730950488
      REAL(R_KIND), PARAMETER :: TWOPI = 6.2831853071795864769 
! ... Scalar arguments ...
      INTEGER(I_KIND), INTENT(IN) :: N, M
! ... Array arguments ...
      REAL(R_KIND), INTENT(INOUT) ::  X(N)
! ... Local scalars ...
      INTEGER(I_KIND)  J, I, K, IS, ID, I0, I1, I2, I3, I4, I5, I6, I7, I8
      INTEGER(I_KIND)  N1, N2, N4, N8
      REAL(R_KIND)  XT, R1, T1, T2, T3, T4, T5
      REAL(R_KIND)  A, A3, E, CC1, SS1, CC3, SS3
!
! ... Exe. statements ...
!
      IF ( N == 1 ) RETURN
!
      N2 = 2 * N
      DO 10, K = 1, M-1
         IS = 0
         ID = N2
         N2 = N2 / 2
         N4 = N2 / 4
         N8 = N4 / 2
         E = TWOPI / N2
 17      DO 15, I = IS, N-1, ID
            I1 = I + 1
            I2 = I1 + N4
            I3 = I2 + N4
            I4 = I3 + N4
            T1 = X(I1) - X(I3)
            X(I1) = X(I1) + X(I3)
            X(I2) = 2 * X(I2)
            X(I3) = T1 - 2 * X(I4)
            X(I4) = T1 + 2 * X(I4)
            IF ( N4 == 1 ) GOTO 15
            I1 = I1 + N8
            I2 = I2 + N8
            I3 = I3 + N8
            I4 = I4 + N8
            T1 = ( X(I2) - X(I1) ) / SQRT2
            T2 = ( X(I4) + X(I3) ) / SQRT2
            X(I1) = X(I1) + X(I2)
            X(I2) = X(I4) - X(I3)
            X(I3) = 2 * ( - T2 - T1 )
            X(I4) = 2 * ( -T2 + T1 )
 15      CONTINUE
         IS = 2 * ID - N2
         ID = 4 * ID
         IF ( IS < N-1 ) GOTO 17
         A = E
         DO 20, J = 2, N8
            A3 = 3 * A
            CC1 = COS(A)
            SS1 = SIN(A)
            CC3 = COS(A3)
            SS3 = SIN(A3)
            A = J * E
            IS = 0
            ID = 2 * N2
 40         DO 30, I = IS, N-1, ID
               I1 = I + J
               I2 = I1 + N4
               I3 = I2 + N4
               I4 = I3 + N4
               I5 = I + N4 - J + 2
               I6 = I5 + N4
               I7 = I6 + N4
               I8 = I7 + N4
               T1 = X(I1) - X(I6)
               X(I1) = X(I1) + X(I6)
               T2 = X(I5) - X(I2)
               X(I5) = X(I2) + X(I5)
               T3 = X(I8) + X(I3)
               X(I6) = X(I8) - X(I3)
               T4 = X(I4) + X(I7)
               X(I2) = X(I4) - X(I7)
               T5 = T1 - T4
               T1 = T1 + T4
               T4 = T2 - T3
               T2 = T2 + T3
               X(I3) = T5 * CC1 + T4 * SS1
               X(I7) = - T4 * CC1 + T5 * SS1
               X(I4) = T1 * CC3 - T2 * SS3
               X(I8) = T2 * CC3 + T1 * SS3
 30         CONTINUE
            IS = 2 * ID - N2
            ID = 4 * ID
            IF ( IS < N-1 ) GOTO 40
 20      CONTINUE
 10   CONTINUE
!
      IS = 1
      ID = 4
 70   DO 60, I0 = IS, N, ID
         I1 = I0 + 1
         R1 = X(I0)
         X(I0) = R1 + X(I1)
         X(I1) = R1 - X(I1)
 60   CONTINUE
      IS = 2 * ID - 1
      ID = 4 * ID
      IF ( IS < N ) GOTO 70
!
 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
         IF ( I >= J ) GOTO 101
         XT = X(J)
         X(J) = X(I)
         X(I) = XT
 101     K = N / 2
 102     IF ( K >= J ) GOTO 103
            J = J - K
            K = K / 2
            GOTO 102
 103     J = J + K
 104  CONTINUE
      XT = 1.0_r_kind / FLOAT( N )
      DO 99, I = 1, N
         X(I) = XT * X(I)
 99   CONTINUE
      RETURN
!
! ... End of subroutine SFFTCB ...
! 
      END SUBROUTINE SFFTCB
subroutine da_read_obs_bufrairs(obstype,iv,infile)
   !--------------------------------------------------------
   !  Purpose: read in NCEP bufr eos AIRS/AMSUA/HSB 1b data 
   !            to innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read file twice
   !            so that da_scan_bufrairs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2006/01/03 - Creation            Zhiquan Liu
   !           2008/03/01 - VISNIR Cloud Cover  Tom Auligne
   !           2008/04/01 - Warmest FoV         Tom Auligne
   !           2008/07/15 - 1 format change  Hui-Chuan Lin
   !            NCEP msg type NC021250 (center FOV data) discontinued after Aug 2007
   !            replaced by NC021249 (every FOV data)
   !
   !------------------------------------------------------------------------------




 implicit none

  character(9)      ,  intent (in)  :: obstype
  character(100)    ,  intent (in)  :: infile
  type (iv_type)    ,intent (inout) :: iv


! Number of channels for sensors in 1
  integer(i_kind),parameter :: N_AIRSCHAN = 281  !- 281 subset ch out of 2378 ch for AIRS
  integer(i_kind),parameter :: N_AMSUCHAN =  15  
  integer(i_kind),parameter :: N_HSBCHAN  =   4
  integer(i_kind),parameter :: N_MAXCHAN  = 350
  integer(i_kind),parameter :: maxinfo    =  12

! 1 format for AQUASPOT (SPITSEQN)
  integer(i_kind),parameter :: N_AQUASPOT_LIST = 25
  type aquaspot_list
     sequence
     real(r_double) :: said   ! Satellite identifier
     real(r_double) :: orbn   ! Orbit number
     real(r_double) :: slnm   ! Scan line number 
     real(r_double) :: mjfc   ! Major frame count
     real(r_double) :: selv   ! Height of station
     real(r_double) :: soza   ! Solar zenith angle
     real(r_double) :: solazi ! Solar azimuth angle
     real(r_double) :: intms(2,9) ! SATELLITE inSTRUMENT TEMPERATURES
  end type aquaspot_list
  real(r_double), dimension(1:N_AQUASPOT_LIST) :: aquaspot_list_array


! 1 format for AIRSSPOT (SITPSEQN)
  integer(i_kind),parameter :: N_AIRSSPOT_LIST = 12
  type airsspot_list
     sequence
     real(r_double) :: siid  ! Satellite instruments
     real(r_double) :: year
     real(r_double) :: mnth
     real(r_double) :: days
     real(r_double) :: hour
     real(r_double) :: minu
     real(r_double) :: seco
     real(r_double) :: clath ! Latitude (high accuracy)
     real(r_double) :: clonh ! Longitude (high accuracy)
     real(r_double) :: saza  ! Satellite zenith angle 
     real(r_double) :: bearaz ! Bearing or azimuth
     real(r_double) :: fovn  ! Field of view number
  end type airsspot_list
  real(r_double), dimension(1:N_AIRSSPOT_LIST) :: airsspot_list_array


! 1 format for AIRSCHAN (SCBTSEQN)
  integer(i_kind),parameter :: N_AIRSCHAN_LIST = 4
  type airschan_list
     sequence
     real(r_double) :: chnm    ! Channel number
     real(r_double) :: logrcw  ! Log-10 of (Temperature-radiance central wavenumber
     real(r_double) :: acqf    ! Channel quality flags for ATOVS
     real(r_double) :: tmbrst  ! Brightness temperature
  end type airschan_list
  real(r_double), dimension(1:N_AIRSCHAN_LIST,1:N_MAXCHAN) :: airschan_list_array
  
! 1 talble file sequencial number
  character(len=512)  :: table_file


! Variables for 1 IO    
  type(aquaspot_list) :: aquaspot
  type(airsspot_list) :: airsspot
  type(airschan_list) :: airschan(N_MAXCHAN)
  
  real(r_kind)      :: step, start
  real(r_kind)      :: airdata(N_AIRSCHAN+maxinfo)
  character(len=8)  :: subset
  character(len=4)  :: senname
  character(len=8)  :: spotname
  character(len=8)  :: channame
  integer(i_kind)   :: nchan,nchanr
  integer(i_kind)   :: iret


! Work variables for time
  integer(i_kind)   :: idate, ifgat
  integer(i_kind)   :: idate5(6)
  character(len=10) :: date
  integer(i_kind)   :: nmind
  integer(i_kind)   :: iy, im, idd, ihh
  real*8            :: obs_time


! Other work variables
  integer(i_kind)  :: nreal
  integer(i_kind)  :: k, iobsout
  integer(i_kind),dimension(19)::icount
  real(r_kind)     :: rlat, rlon, dx, dy, dx1, dy1, sstx, dlon, dlat
  real(r_kind)     :: sza, timedif, pred, crit1
  integer(i_kind)  :: klat1, klon1, klatp1, klonp1
  integer(i_kind)  :: ifov,size
  integer(i_kind)  :: num_eos_file, num_eos_global, num_eos_local, num_eos_kept, num_eos_thinned
  integer(i_kind)  :: inst,platform_id,satellite_id,sensor_id
  logical          :: iflag,outside,outside_all
  integer(i_kind)  :: i, l, n, error, airs_table_unit
  integer          :: iost, lnbufr
  real(r_kind),allocatable,dimension(:,:):: airdata_all
  real(kind=8)     :: tocc(1,1)
  integer          ::num_bufr(7),numbufr,ibufr
  character(20)    ::filename
  integer(i_kind), allocatable :: ptotal(:)

! Set standard parameters
  real(r_kind)     :: POinT001 =   0.001_r_kind
  real(r_kind)     :: POinT01  =   0.01_r_kind
  real(r_kind)     :: TEN      =  10.0_r_kind
  real(r_kind)     :: R45      =  45.0_r_kind
  real(r_kind)     :: R60      =  60.0_r_kind
  real(r_kind)     :: R90      =  90.0_r_kind
  real(r_kind)     :: R180     = 180.0_r_kind
  real(r_kind)     :: R360     = 360.0_r_kind

! Thinning variables
  integer(i_kind) itt,itx,iobs,iout,size_tmp, j
  real(r_kind) crit,dist
  real(r_kind) dlon_earth,dlat_earth
  logical luse
  real , allocatable :: in(:), out(:)
  logical :: found, head_found

  logical           :: airs, eos_amsua, hsb, airstab
  type(info_type)            :: info
  type(model_loc_type)       :: loc
  type (datalink_type), pointer   :: head, p, current, prev

  if (trace_use) call da_trace_entry("da_read_obs_bufrairs")

!  0.0  Initialize variables
!-----------------------------------
  platform_id  = 9   ! eos series
  satellite_id = 2   ! eos-2
  nreal  = maxinfo
  allocate(ptotal(0:num_fgat_time))
  ptotal(0:num_fgat_time) = 0
  airs=      obstype == 'airs     '
  eos_amsua= obstype == 'eos_amsua'
  hsb=       obstype == 'hsb      '

  icount=0
  if(airs)then
     sensor_id = 11
     step   = 1.1_r_kind
     start = -48.9_r_kind
     senname = 'AIRS'
     nchan  = N_AIRSCHAN
     nchanr = N_AIRSCHAN
  else if(eos_amsua)then
     sensor_id = 3
     step   = three + one/three
     start  = -48.33_r_kind
     senname = 'AMSU'
     nchan  = N_AMSUCHAN
     nchanr = N_AMSUCHAN
  else if(hsb)then
     sensor_id = 12
     step   = 1.1_r_kind
     start  = -48.95_r_kind
     senname = 'HSB'
     nchan  = N_HSBCHAN
     nchanr = N_HSBCHAN+1
  end if
  spotname = trim(senname)//'SPOT'
  channame = trim(senname)//'CHAN'
 
      do inst = 1, rtminit_nsensor
        if (    platform_id  == rtminit_platform(inst) &
          .and. satellite_id == rtminit_satid(inst)    &
          .and. sensor_id    == rtminit_sensor(inst)    ) then
            exit
        end if
      end do

      if ( inst == rtminit_nsensor .and.           &
           platform_id  /= rtminit_platform(inst)  &
          .or. satellite_id /= rtminit_satid(inst) &
          .or. sensor_id /= rtminit_sensor(inst)  ) then
         if (trace_use) call da_trace_exit("da_read_obs_bufrairs")
         return
      end if

  num_eos_file    = 0  ! number of obs in file
  num_eos_global  = 0  ! number of obs within whole domain
  num_eos_local   = 0  ! number of obs within tile
  num_eos_kept    = 0  ! number of obs kept for innovation computation
  num_eos_thinned = 0  ! number of obs rejected by thinning
  size = 0     !
  iobs = 0     ! for thinning, argument is inout

!    1.0  Open 1 table and 1 file
!--------------------------------------------------------------
  table_file = 'gmao_airs_bufr.tbl'      ! make table file name
  inquire(file=table_file,exist=airstab)
  if (airstab) then
      if (print_detail_rad) then
         write(unit=message(1),fmt=*) &
            'Reading BUFR Table A file: ',trim(table_file)
         call da_message(message(1:1))
      end if
      call da_get_unit(airs_table_unit)
      open(unit=airs_table_unit,file=table_file,iostat = iost)
      if (iost /= 0) then
         call da_error("da_read_obs_bufrairs.inc",233, &
            (/"Cannot open file "//table_file/))
      end if
  end if

! Open 1 file
!  call da_get_unit(lnbufr)

   num_bufr(:)=0
   numbufr=0
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
	    num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if

   if (numbufr==0) numbufr=1

bufrfile:  do ibufr=1,numbufr   
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'   
      end if
   end if

  lnbufr=97

  open(unit=lnbufr,file=trim(filename),form='unformatted',iostat = iost, status = 'old')
  if (iost /= 0) then
     call da_warning("da_read_obs_bufrairs.inc",277, &
        (/"Cannot open file "//infile/))
     close(lnbufr)
     if (airstab) then
        close(airs_table_unit)
        call da_free_unit(airs_table_unit)
     end if
     if (trace_use) call da_trace_exit("da_read_obs_bufrairs")
     return
  end if
  if ( airstab ) then
     call openbf(lnbufr,'IN',airs_table_unit)
  else
     call openbf(lnbufr,'IN',lnbufr)
  end if
  call datelen(10)


!   2.0  Read header
!---------------------------
  call readmg(lnbufr,subset,idate,iret)

  iy = 0
  im = 0
  idd = 0
  ihh = 0
  if( iret /= 0 ) goto 1000     ! no data?

  write(unit=date,fmt='( i10)') idate
  read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
  write(unit=stdout,fmt='(a,4i4,2x,a)') &
      'Bufr file date is ',iy,im,idd,ihh,trim(infile)

!   3.0 Loop over observations
!----------------------------
     if ( ibufr == 1 ) then
        allocate ( head )
      ! allocate ( head % tb (1:nchan) )
        nullify  ( head % next )
        p => head
     endif  

  loop_obspoints: do

!   3.1 Read headder
!-------------------------------
     call readsb(lnbufr,iret)

     if( iret /=0 )then
        call readmg(lnbufr,subset,idate,iret)
        if( iret /= 0 ) exit loop_obspoints     ! end of file
        cycle loop_obspoints
     end if

     num_eos_file = num_eos_file + 1

!   3.2 Read AQUASPOT (SPITSEQN)
!------------------------
     call ufbseq(lnbufr,aquaspot_list_array,N_AQUASPOT_LIST,1,iret,'SPITSEQN')
     aquaspot = aquaspot_list( aquaspot_list_array(1), &
                               aquaspot_list_array(2), &
                               aquaspot_list_array(3), &
                               aquaspot_list_array(4), &
                               aquaspot_list_array(5), &
                               aquaspot_list_array(6), &
                               aquaspot_list_array(7), &
                               RESHAPE(aquaspot_list_array(8:25), (/2,9/)) )

!   3.3 Read AIRSSPOT or AMSUSPOT or HSBSPOT
!-------------------------------------------------
     if ( trim(senname) == 'AIRS' ) then
        call ufbseq(lnbufr,airsspot_list_array,N_AIRSSPOT_LIST,1,iret,'SITPSEQN')
     else
        call ufbseq(lnbufr,airsspot_list_array,N_AIRSSPOT_LIST,1,iret,spotname)
     end if
     airsspot = airsspot_list( airsspot_list_array(1), &
                               airsspot_list_array(2), &
                               airsspot_list_array(3), &
                               airsspot_list_array(4), &
                               airsspot_list_array(5), &
                               airsspot_list_array(6), &
                               airsspot_list_array(7), &
                               airsspot_list_array(8), &
                               airsspot_list_array(9), &
                               airsspot_list_array(10), &
                               airsspot_list_array(11), &
                               airsspot_list_array(12) )

!   3.4 Read AIRSCHAN or AMSUCHAN or HSBCHAN
!-------------------------------------------
     if ( trim(senname) == 'AIRS' ) then
        call ufbseq(lnbufr,airschan_list_array,N_AIRSCHAN_LIST,N_MAXCHAN,iret,'SCBTSEQN')
     else
        call ufbseq(lnbufr,airschan_list_array,N_AIRSCHAN_LIST,N_MAXCHAN,iret,channame)
     end if
     do l = 1 , N_MAXCHAN
        airschan(l) = airschan_list( airschan_list_array(1,l), &
                                     airschan_list_array(2,l), & 
                                     airschan_list_array(3,l), & 
                                     airschan_list_array(4,l)  )
     end do

     if (iret /= nchanr) then
        write(unit=message(1),fmt=*) &
            'Cannot read ', channame, &
            ' bufr data:', &
            iret, ' ch data is read instead of', nchanr
        call da_warning("da_read_obs_bufrairs.inc",384,message(1:1))
        cycle loop_obspoints
     end if

!   3.5 Read Cloud Cover from AIRS/VISNIR
!-------------------------------------------
     call ufbint(lnbufr,tocc,1,1,iret,'TOCC')
     
!   4.0  Check observing position (lat/lon)
!      QC1:  juge if data is in the domain, 
!            read next record if not
!------------------------------------------
     if( abs(airsspot%clath) > R90  .or. &
          abs(airsspot%clonh) > R360 .or. &
          (abs(airsspot%clath) == R90 .and. airsspot%clonh /= ZERO) )then
        cycle loop_obspoints
     end if

!    Retrieve observing position
     if(airsspot%clonh >= R360) then
        airsspot%clonh = airsspot%clonh - R360
!     else if(airsspot%clonh < ZERO) then
!        airsspot%clonh = airsspot%clonh + R360
     end if

        info%lat  = airsspot%clath
        info%lon  = airsspot%clonh 
        call da_llxy (info, loc, outside, outside_all )

        if ( outside_all ) cycle loop_obspoints
	
!  4.1  Check obs time
!-------------------------------------
     idate5(1) = airsspot%year ! year
     idate5(2) = airsspot%mnth ! month
     idate5(3) = airsspot%days ! day
     idate5(4) = airsspot%hour ! hour
     idate5(5) = airsspot%minu ! minute
     idate5(6) = airsspot%seco ! second

     if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
          idate5(2) <    1 .or. idate5(2) >   12 .or. &
          idate5(3) <    1 .or. idate5(3) >   31 .or. &
          idate5(4) <    0 .or. idate5(4) >   24 .or. &
          idate5(5) <    0 .or. idate5(5) >   60 .or. &
          idate5(6) <    0 .or. idate5(6) >   60 ) then
        cycle loop_obspoints
     end if

!  QC3: time consistency check with the background date

      if (idate5(1) .LE. 99) then
        if (idate5(1) .LT. 78) then
          idate5(1) = idate5(1) + 2000
        else
          idate5(1) = idate5(1) + 1900
        end if
      end if

      call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)

      if ( obs_time < time_slots(0) .or.  &
           obs_time >= time_slots(num_fgat_time) ) cycle

!  3.2.1   determine FGAT index ifgat
!
       do ifgat=1,num_fgat_time
           if ( obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat) ) exit
       end do

       num_eos_global = num_eos_global + 1
       ptotal(ifgat) = ptotal(ifgat) + 1

       if (outside) cycle loop_obspoints

       num_eos_local = num_eos_local + 1

       write(unit=info%date_char, &
         fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
         idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
         ':', idate5(5), ':', idate5(6)

       info%elv = 0.0  !aquaspot%selv

!  4.2  Check observational info
!-------------------------------------------------------
     if( airsspot%fovn <    0.0_r_kind .or. airsspot%fovn > 100.0_r_kind .or. &
          airsspot%saza < -360.0_r_kind .or. airsspot%saza > 360.0_r_kind .or. &
          aquaspot%soza < -180.0_r_kind .or. aquaspot%soza > 180.0_r_kind )then
         write(unit=message(1),fmt=*) &
            'Cannot read ', channame, ' bufr data:', &
            ' strange obs info(fov,saza,soza):', &
            airsspot%fovn, airsspot%saza, aquaspot%soza
        call da_warning("da_read_obs_bufrairs.inc",478,message(1:1))
        cycle loop_obspoints
     end if

!    Retrieve observing info
     ifov = int( airsspot%fovn + POinT001 )
     sza  = abs(airsspot%saza)
!     if( ((airs .or. hsb) .and. ifov <= 45) .or. &
!          ( eos_amsua     .and. ifov <= 15) )then
!        sza = - sza
!     end if

!     airdata(6) = (start + float(ifov-1)*step)  ! look angle (deg)
!     airdata(9) = ZERO                          ! surface height
!     airdata(10)= POinT001                      ! land sea mask

!  4.1 Make Thinning
!  Map obs to thinning grid
!-------------------------------------------------------------------
       if (thinning) then
          dlat_earth = info%lat
          dlon_earth = info%lon
          if (dlon_earth<zero) dlon_earth = dlon_earth+r360
          if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
          dlat_earth = dlat_earth*deg2rad
          dlon_earth = dlon_earth*deg2rad
          crit = 1.0 ! 0.01_r_kind+terrain + timedif !+ 10.0_r_kind*float(iskip)
          if (airs_warmest_fov) &
             crit = 1E10 * exp(-(airschan(129)%tmbrst-220.0)/2)   ! warmest bt for window channel (10.36 micron)
          call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,luse)
          if (.not. luse) then
             num_eos_thinned = num_eos_thinned + 1
             cycle loop_obspoints
          end if
       end if


!   4.3 Retrieve Tb
!-----------------------
     iflag = .false.
  
     do l=1,nchan
        airdata(l+nreal) = airschan(l)%tmbrst            ! brightness temperature
        if( airdata(l+nreal) > 0.0_r_kind .and. airdata(l+nreal) < 500.0_r_kind )then
           iflag = .true.
        else
           airdata(l+nreal) = missing_r
        end if
     end do

     if ( .not. iflag )then
        write(unit=message(1),fmt=*) &
          'Error in reading ', channame, ' bufr data:', &
          ' all tb data is missing'
        call da_warning("da_read_obs_bufrairs.inc",532,message(1:1))
        num_eos_local = num_eos_local - 1
        cycle loop_obspoints
     end if

     num_eos_kept = num_eos_kept + 1

!  4.0   assign information to sequential radiance structure
!--------------------------------------------------------------------------
   allocate ( p % tb_inv (1:nchan) )
   p%info             = info
   p%loc              = loc
   p%landsea_mask     = POinT001
   p%scanline         = int(aquaspot%slnm + POinT001)
   p%scanpos          = ifov
   p%satzen           = sza
   p%satazi           = (start + float(ifov-1)*step)  ! look angle (deg) ! airsspot%bearaz
   p%solzen           = aquaspot%soza
   p%solazi           = aquaspot%solazi
   p%tb_inv(1:nchan)  = airdata(nreal+1:nreal+nchan)
   p%sensor_index     = inst
   p%ifgat            = ifgat

   size = size + 1
   allocate ( p%next, stat=error)   ! add next data
   if (error /= 0 ) then
      call da_error("da_read_obs_bufrairs.inc",558, &
          (/"Cannot allocate radiance structure"/))
   end if

   p => p%next
   nullify (p%next)

  end do loop_obspoints

  call closbf(lnbufr)
  close(lnbufr)

end do bufrfile

   if (thinning .and. num_eos_global > 0) then

         
      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do
      end do
         
      allocate ( in  (j) )
      allocate ( out (j) ) 
               
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i)
            end do
         end do             
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i)  = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      if ( size > 0 ) then
         p => head
         prev => head
         head_found = .false.
         size_tmp = size
         do j = 1, size_tmp
            n = p%sensor_index
            ifgat = p%ifgat
            found = .false.

            do i = 1, thinning_grid(n,ifgat)%itxmax
               if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
                  found = .true.
                  exit
               endif
            end do

            ! free current data
            if ( .not. found ) then
               current => p
               p => p%next
               if ( head_found ) then
                  prev%next => p
               else
                  head => p
                  prev => p
               endif
               deallocate ( current % tb_inv )
               deallocate ( current )
               num_eos_thinned = num_eos_thinned + 1
               num_eos_kept = num_eos_kept - 1
               size = size - 1
               continue
            endif

            if ( found .and. head_found ) then
               prev => p
               p => p%next
               continue
            endif

            if ( found .and. .not. head_found ) then
               head_found = .true.
               head => p
               prev => p
               p => p%next
            endif
         end do
      endif

   endif  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + size
   iv%total_rad_channel = iv%total_rad_channel + size*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + size
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_eos_global
   iv%instid(inst)%info%nlocal = size
   
   do ifgat = 1, num_fgat_time
      ptotal(ifgat) = ptotal(ifgat) + ptotal(ifgat-1)
      iv%info(radiance)%ptotal(ifgat) = iv%info(radiance)%ptotal(ifgat) + ptotal(ifgat)
   end do
   
   write(unit=stdout,fmt='(a)') '   num_eos_file num_eos_global  num_eos_local   num_eos_kept num_eos_thinned'
   write(unit=stdout,fmt='(5(5x,i10))') num_eos_file,num_eos_global, num_eos_local,num_eos_kept,num_eos_thinned

!  5.0 allocate innovation radiance structure
!----------------------------------------------------------------
!  do i = 1, iv%num_inst
   if ( size > 0 ) then
      iv%instid(inst)%num_rad = size
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad
   end if
!  end do

!   6.0 assign sequential structure to innovation structure
!-------------------------------------------------------------

  n = 0
  p => head

  call da_allocate_rad_iv (inst, nchan, iv)

  do i = 1, size
!   inst = p%sensor_index
   n = n + 1

   call da_initialize_rad_iv (inst, n, iv, p)
   iv%instid(inst)%rain_flag(n) = tocc(1,1) ! Temporary dumping of AIRS/VISNIR cloud cover
   
   current => p
   p => p%next

! free current data
   deallocate ( current % tb_inv )
   deallocate ( current )

 end do

 deallocate (p)
 deallocate (ptotal)

1000 continue
   call closbf(lnbufr)
   close(lnbufr)
!   call da_free_unit(lnbufr)
   if (airstab) then
      close(airs_table_unit)
      call da_free_unit(airs_table_unit)
   end if

   if (trace_use) call da_trace_exit("da_read_obs_bufrairs")

end subroutine da_read_obs_bufrairs
subroutine da_read_obs_bufrssmis (obstype,iv,infile)

   !---------------------------------------------------------------------------
   !  Purpose: read in NCEP bufr SSM/IS data to innovation structure
   !
   !   METHOD: use F90 sequential data structure to avoid reading file twice  
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !               and deallocate sequential data structure
   !---------------------------------------------------------------------------



   use da_control

   implicit none

   character(5) ,     intent (in)    :: obstype    ! ssmis
   character(20),     intent (in)    :: infile     ! ssmis.bufr
   type (iv_type),    intent (inout) :: iv


   integer(i_kind), parameter :: bufsat_dmsp16 = 249  ! DMSP16 1 identifier
   integer(i_kind), parameter :: bufsat_dmsp17 = 285  ! DMSP17 1 identifier
   integer(i_kind), parameter :: bufsat_dmsp18 = 286  ! DMSP18 1 identifier
   integer(i_kind), parameter :: n1bhdr = 15
   integer(i_kind), parameter :: maxchanl = 24
   real(r_kind),    parameter :: tbmin = 70.0_r_kind
   real(r_kind),    parameter :: tbmax = 320.0_r_kind

   character(80) :: hdr1b
   !data hdr1b /'SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SLNM ORBN SELV SURF RAINF'/
   data hdr1b /'SAID FOVN YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SLNM ORBN SELV SFLG RFLAG'/
   character(10) :: date
   character(8)  :: subset, subfgn
   character(20) :: filename

   logical :: outside, outside_all

   integer(i_kind) :: iost, inst, lnbufr, ifgat
   integer(i_kind) :: num_bufr(7),numbufr,ibufr
   integer(i_kind) :: ihh, i, j, n, k, slnm, ifov, idd, ireadmg, ireadsb
   integer(i_kind) :: iret, idate, im, iy
   integer(i_kind) :: jc, incangl, bch, landsea_mask, rain_flag
   integer(i_kind) :: platform_id, satellite_id, sensor_id, nchan, num_ssmis_file
   integer(i_kind) :: num_ssmis_local, num_ssmis_global, num_ssmis_used, num_ssmis_thinned
   integer(i_kind) :: num_ssmis_used_tmp

   real(r_double), dimension(2,maxchanl) :: bufrtbb
   real(r_double), dimension(n1bhdr)     :: bfr1bhdr

   ! pixel information
   integer(i_kind)   ::  year,month,day,hour,minute,second  ! observation time
   real(kind=8)    ::  obs_time
   real(r_double), allocatable ::  tb_inv(:)          !  bright temperatures

   type (datalink_type), pointer  :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   ! thinning variables
   integer(i_kind) :: itt,itx,iobs,iout
   real(r_kind)    :: terrain,timedif,crit,dist
   real(r_kind)    :: dlon_earth,dlat_earth
   logical         :: iuse
   real, allocatable :: in(:), out(:)
   logical           :: found, head_found

   integer(i_kind), allocatable :: ptotal(:), nread(:)

   call da_trace_entry("da_read_obs_bufrssmis")

   allocate(nread(1:rtminit_nsensor))
   allocate(ptotal(0:num_fgat_time))
   nread(1:rtminit_nsensor) = 0
   ptotal(0:num_fgat_time) = 0

   platform_id  = 2                 ! for DMSP series
   sensor_id    = 10                ! for SSMIS
   nchan        = nchan_ssmis

   allocate (tb_inv(nchan))
   num_ssmis_file    = 0
   num_ssmis_local   = 0
   num_ssmis_global  = 0
   num_ssmis_used    = 0
   num_ssmis_thinned = 0
   iobs = 0                 ! for thinning, argument is inout

   ! 0.0  Open unit to satellite bufr file and read file header
   !--------------------------------------------------------------

!   call da_get_unit(lnbufr)

   num_bufr(:)=0
   numbufr=0
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
            num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if
 
   if (numbufr==0) numbufr=1
 
bufrfile:  do ibufr=1,numbufr
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'
      end if
   end if

   lnbufr=98
   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_warning("da_read_obs_bufrssmis.inc",133, &
         (/"Cannot open file "//infile/))
      call da_trace_exit("da_read_obs_bufrssmis")
      return
   end if

   call openbf(lnbufr,'IN',lnbufr)
   call datelen(10)
   call readmg(lnbufr,subset,idate,iret)

   iy=0
   im=0
   idd=0
   ihh=0
   write(unit=date,fmt='( i10)') idate
   read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
   write(unit=stdout,fmt=*) &
      'Bufr file date is ',iy,im,idd,ihh,infile

   ! Loop to read bufr file and assign information to a sequential structure
   !-------------------------------------------------------------------------
   if ( ibufr == 1 ) then
      allocate (head)
      nullify  ( head % next )
      p => head
   end if

! Set various variables depending on type of data to be read

   !subfgn = 'NC003003'
   subfgn = 'NC021201'
   incangl = 53.2_r_kind

   subset_loop: do while (ireadmg(lnbufr,subset,idate)==0)

      read_loop: do while (ireadsb(lnbufr)==0 .and. subset==subfgn)

         num_ssmis_file = num_ssmis_file + 1
         ! 1.0     Read header record and data record

         call ufbint(lnbufr,bfr1bhdr,n1bhdr,1,iret,hdr1b)
         call ufbrep(lnbufr,bufrtbb,2,maxchanl,iret,"CHNM TMBR" )

         ! check if observation outside range

         ! 2.0     Extract observation location and other required information
         !     QC1:  judge if data is in the domain, read next record if not
         !------------------------------------------------------------------------

         info%lat  =  bfr1bhdr(bufr_lat)
         info%lon  =  bfr1bhdr(bufr_lon)
         call da_llxy (info, loc, outside, outside_all)


         if (outside_all) cycle

         !  3.0     Extract other information

         info%elv  = 0.0
         landsea_mask = nint(bfr1bhdr(bufr_landsea_mask))   ! ssmis surface flag
                                                            ! 0:land, 2:near coast, 3:ice,
                                                            ! 4:possible ice, 5:ocean, 6:coast
         ! RTTOV surftype: 0:land, 1:sea, 2:sea ice
         if ( landsea_mask == 5 ) then
            landsea_mask = 1
         else if ( landsea_mask == 2 .or. landsea_mask == 6 ) then
            landsea_mask = 0
         else if ( landsea_mask == 3 .or. landsea_mask == 4 ) then
            landsea_mask = 2
         end if
         rain_flag = nint(bfr1bhdr(15))    ! 0:no rain, 1:rain

         !------------------------------------------------------
         !  3.1     Extract satellite id and scan position. 
   
         if (nint(bfr1bhdr(bufr_satellite_id)) == bufsat_dmsp16) then
            satellite_id = 16
         else if (nint(bfr1bhdr(bufr_satellite_id)) == bufsat_dmsp17) then
            satellite_id = 17
         else if (nint(bfr1bhdr(bufr_satellite_id)) == bufsat_dmsp18) then
            satellite_id = 18
         end if

         ! 3.3 Find wrfvar instrument index from RTTOV instrument triplet
         !     go to next data if id is not in the lists

         inst = 0
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i)      &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do
         if (inst == 0) cycle read_loop

         !  3.1     Extract scan number and scan position. 

         slnm = nint(bfr1bhdr(11))
         ifov = nint(bfr1bhdr(bufr_ifov))

         !  3.2     Extract date information.
    
         year   = bfr1bhdr(bufr_year)   
         month  = bfr1bhdr(bufr_month)  
         day    = bfr1bhdr(bufr_day)    
         hour   = bfr1bhdr(bufr_hour)   
         minute = bfr1bhdr(bufr_minute) 
         second = bfr1bhdr(bufr_second) 

         write(unit=info%date_char, fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
            year, '-', month, '-', day, '_', hour, ':', minute, ':', second

         !  QC3: time consistency check with the background date

         if (year <= 99) then
            if (year < 78) then
               year = year + 2000
            else
               year = year + 1900
            end if
         end if

         call da_get_julian_time(year,month,day,hour,minute,obs_time)

         if (obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time)) cycle read_loop

         ! 3.2.1 determine FGAT index ifgat
   
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do

         num_ssmis_global = num_ssmis_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1

         if (outside) cycle ! No good for this PE

         num_ssmis_local = num_ssmis_local + 1

         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad
            timedif = 0. !2.0_r_kind*abs(tdiff)        ! range:  0 to 6
            terrain = 0.01_r_kind*abs(bfr1bhdr(13))
            crit = 1. !0.01_r_kind+terrain + timedif !+ 10._r_kind*float(iskip)
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_ssmis_thinned=num_ssmis_thinned+1
               cycle
            end if
         end if

         num_ssmis_used = num_ssmis_used + 1
         nread(inst) = nread(inst) + 1

         if (num_ssmis_used > max_ssmis_input) then
            write(unit=message(1),fmt='(A,I10,A)') &
               "Max number of ssmis",max_ssmis_input," reached"
            call da_warning("da_read_obs_bufrssmis.inc",302,message(1:1))
            num_ssmis_used = num_ssmis_used - 1
            exit read_loop
         end if

         ! 3.4 extract satellite and solar angle
   
         ! 3.5 extract surface information

         ! 3.6 extract channel bright temperature
   
         tb_inv(1:nchan) = missing_r

         do jc = 1, nchan
            bch = nint(bufrtbb(1,jc))     !ch index from bufr
            tb_inv(jc) = bufrtbb(2,jc)
            if (tb_inv(jc) < tbmin .or. tb_inv(jc) > tbmax .or. bch /= jc) then
                tb_inv(jc) = missing_r
            end if
         end do

         if ( maxval(tb_inv(:)) > missing_r ) then

            !  4.0   assign information to sequential radiance structure
            !--------------------------------------------------------------------------
            allocate (p % tb_inv (1:nchan))
            p%info                = info
            p%loc                 = loc
            p%landsea_mask        = landsea_mask
            p%scanline            = slnm
            p%scanpos             = ifov
            p%satzen              = incangl
            p%satazi              = 0.0     ! dummy value
            p%solzen              = 0.0     ! dummy value
            p%tb_inv(1:nchan)     = tb_inv(1:nchan)
            p%sensor_index        = inst
            p%ifgat               = ifgat
            p%rain_flag           = rain_flag

            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)

         else

            num_ssmis_local = num_ssmis_local - 1
            num_ssmis_used  = num_ssmis_used - 1

         end if

      end do read_loop

   end do subset_loop

  call closbf(lnbufr)
  close(lnbufr)

end do bufrfile

   if (thinning .and. num_ssmis_global > 0 ) then

      
      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do 
      end do
   
      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i)
            end do
         end do 
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i) = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_ssmis_used_tmp = num_ssmis_used
      do j = 1, num_ssmis_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_ssmis_thinned = num_ssmis_thinned + 1
            num_ssmis_used = num_ssmis_used - 1
            nread(n) = nread(n) - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_ssmis_used
   iv%total_rad_channel = iv%total_rad_channel + num_ssmis_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_ssmis_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_ssmis_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_bufrssmis.inc",468,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'num_ssmis_file, num_ssmis_global, num_ssmis_local, num_ssmis_used, num_ssmis_thinned'
   write(stdout,*) num_ssmis_file, num_ssmis_global, num_ssmis_local, num_ssmis_used, num_ssmis_thinned

   deallocate(tb_inv)  

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   do i = 1, iv%num_inst 
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad  = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad

      call da_allocate_rad_iv (i, nchan, iv)

   end do

   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0
   p => head

   do n = 1, num_ssmis_used
      i = p%sensor_index
      nread(i) = nread(i) + 1
      call da_initialize_rad_iv (i, nread(i), iv, p)

      iv%instid(i)%rain_flag(n) = p%rain_flag

      current => p
      p => p%next

      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )

   end do

   deallocate ( p )
   deallocate (nread)
   deallocate (ptotal)

   call closbf(lnbufr)
   close(lnbufr)
!   call da_free_unit(lnbufr)

   call da_trace_exit("da_read_obs_bufrssmis")
 

end subroutine da_read_obs_bufrssmis

subroutine da_read_obs_bufriasi (obstype,iv,infile)	 
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_iasi                  read bufr format iasi data

  implicit none


  character(9)      ,  intent (in)  :: obstype
  character(100)    ,  intent (in)  :: infile
  type (iv_type)    ,intent (inout) :: iv
  real(kind=8)    ::  obs_time  
  type (datalink_type), pointer  :: head, p, current, prev
  type(info_type)                :: info
  type(model_loc_type)           :: loc 
  type(model_loc_type)           :: loc_fov  

! Subroutine constants
  real(r_kind)                   :: ten      =  10.0_r_kind
  real(r_kind)                   :: r90      =  90.0_r_kind
  real(r_kind)                   :: r360     = 360.0_r_kind

! Number of channels for sensors in 1
  integer(i_kind),parameter      :: nchan = 616        !--- 616 subset ch out of 8078 ch for AIRS
  integer(i_kind),parameter      :: n_totchan  = 616
  integer(i_kind),parameter      :: maxinfo    =  33
  integer(i_kind)                :: inst,platform_id,satellite_id,sensor_id
  integer(i_kind), allocatable   :: ptotal(:), nread(:)
  real(r_kind)                   :: crit
  integer(i_kind)                :: ifgat, iout, iobs
  logical                        :: outside, outside_all, iuse

! 1 functions
  integer(i_kind)   :: ireadsb,ireadmg

! Variables for 1 IO    
  real(r_double),dimension(5)           :: linele
  real(r_double),dimension(14)          :: allspot
  real(r_double),dimension(2,n_totchan) :: allchan
  real(r_double),dimension(3,10)        :: cscale
  real(r_double),dimension(6)           :: cloud_frac
  integer           :: numbufr,ibufr  
  logical           :: found, head_found 
  real(r_kind)      :: step, start,step_adjust
  character(len=8)  :: subset
  character(len=4)  :: senname
  character(len=80) :: allspotlist
  integer(i_kind)   :: jstart
  integer(i_kind)   :: iret
  character(10)     :: date

! Work variables for time
  integer(i_kind)   :: idate, im, iy, idd, ihh
  integer(i_kind)   :: idate5(6)

! Other work variables
  real(r_kind)                          :: piece
  real(r_kind)                          :: dlon_earth,dlat_earth
  real(r_kind)                          :: lza, lzaest,sat_height_ratio
  real(r_kind)                          :: sat_zenang
  real(r_kind)                          :: radi
  real(r_kind),dimension(10)            :: sscale
  real(r_kind),dimension(n_totchan)     :: temperature
  real(r_kind),allocatable,dimension(:) :: data_all
  real(r_kind)                          :: Effective_Temperature

  logical          :: iasi
  integer(i_kind)  :: nreal, ksatid
  integer(i_kind)  :: ifov, iscn
  integer(i_kind)  :: num_iasi_file
  integer(i_kind)  :: num_iasi_local, num_iasi_global, num_iasi_used, num_iasi_thinned 
  integer(i_kind)  :: num_iasi_used_tmp  
  integer(i_kind)  :: i, j, l, iskip, ifovn, bad_line
  integer(i_kind)  :: itx, k, nele, itt, n
  integer(i_kind)  :: iexponent
  integer(i_kind)  :: num_bufr(7)
  integer          :: iost, lnbufr 
  character(20)    :: filename  
  real,allocatable :: in(:), out(:)

! Set standard parameters
  integer(i_kind),parameter :: ichan=-999  ! fov-based surface code is not channel specific for iasi 
  real(r_kind),parameter    :: expansion=one         ! exansion factor for fov-based location code.                                                 ! use one for ir sensors.
  real(r_kind),parameter    :: tbmin  = 50._r_kind
  real(r_kind),parameter    :: tbmax  = 550._r_kind
  real(r_kind),parameter    :: earth_radius = 6371000._r_kind
  

   if (trace_use) call da_trace_entry("da_read_obs_bufriasi")
!  0.0  Initialize variables
!-----------------------------------
  platform_id  = 10   ! metop series
  sensor_id = 16 !iasi
  nreal  = maxinfo
  iasi=      obstype == 'iasi'
  bad_line=-1
  step  = 3.334
  start = -48.33
  step_adjust = 0.625_r_kind
  senname = 'IASI'
  allspotlist= &
   'SIID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH SAZA BEARAZ SOZA SOLAZI SAID'
  num_bufr(:)=0
  numbufr=0
  allocate(nread(1:rtminit_nsensor))
  allocate(ptotal(0:num_fgat_time))
  nread(1:rtminit_nsensor) = 0
  ptotal(0:num_fgat_time) = 0
  iobs = 0                 ! for thinning, argument is inout  
  num_iasi_file    = 0
  num_iasi_local   = 0
  num_iasi_global  = 0
  num_iasi_used    = 0
  num_iasi_thinned = 0  
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
            num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if
 
   if (numbufr==0) numbufr=1
 
  bufrfile:  do ibufr=1,numbufr
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'
      end if
   end if 
   lnbufr = 95
   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_warning("da_read_obs_bufriasi.inc",149, &
         (/"Cannot open file "//infile/))
      if (trace_use) call da_trace_exit("da_read_obs_bufriasi")
      return
   end if

! Open 1 table
   call openbf(lnbufr,'IN',lnbufr) 
   call datelen(10)
   call readmg(lnbufr,subset,idate,iret)
   iy=0
   im=0
   idd=0
   ihh=0
   write(unit=date,fmt='( i10)') idate
   read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
   write(unit=stdout,fmt='(a,4i4,2x,a)') &
      'Bufr file date is ',iy,im,idd,ihh,trim(infile)

   ! Loop to read bufr file and assign information to a sequential structure
   !-------------------------------------------------------------------------
! Allocate arrays to hold data
  nele=nreal+nchan
  allocate(data_all(nele))
   if ( ibufr == 1 ) then
      allocate (head)
      nullify  ( head % next )
      p => head
   end if
  

! Big loop to read data file

  do while(ireadmg(lnbufr,subset,idate)>=0)  

     read_loop: do while (ireadsb(lnbufr)==0)
         num_iasi_file = num_iasi_file + 1	 

!    Read IASI FOV information
         call ufbint(lnbufr,linele,5,1,iret,'FOVN SLNM QGFQ MJFC SELV')				
         if ( linele(3) /= zero) then
	         cycle read_loop  
	     end if
         if ( bad_line == nint(linele(2))) then
!        zenith angle/scan spot mismatch, reject entire line
           cycle read_loop
         else
           bad_line = -1
         end if
          call ufbint(lnbufr,allspot,14,1,iret,allspotlist)
         if(iret /= 1) cycle read_loop

         ksatid = nint(allspot(14))
         ! SAID 3 is metop-b (metop-1)
         ! SAID 4 is metop-a (metop-2)
         ! SAID 5 is metop-c (metop-3)
         if ( ( ksatid > 5) .or. ( ksatid < 3) ) then 
            write(unit=message(1),fmt='(A,I6)') 'Unknown platform: ', ksatid
            call da_warning("da_read_obs_bufriasi.inc",207,message(1:1))
         end if
         if ( ksatid == 3 ) then
             satellite_id = 1
         else if ( ksatid == 4 ) then
             satellite_id = 2
         else if ( ksatid == 5 ) then
             satellite_id = 3
         end if

!    Check observing position
         info%lat  =  allspot(8)  ! latitude
         info%lon  =  allspot(9)  ! longitude)
         if( abs(info%lat) > r90  .or. abs(info%lon) > r360 .or. &
         (abs(info%lat) == r90 .and. info%lon /= zero) )then
         write(unit=stdout,fmt=*) &
         'READ_IASI:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
               ' STRANGE OBS POINT (LAT,LON):', info%lat, info%lon
            cycle read_loop
         end if		 
	 		 
         call da_llxy (info, loc, outside, outside_all)	
         if (outside_all) cycle
	     inst = 0	 
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do	 
         if (inst == 0) cycle read_loop		 
		 
!    Check obs time
         idate5(1) = nint(allspot(2)) ! year
         idate5(2) = nint(allspot(3)) ! month
         idate5(3) = nint(allspot(4)) ! day
         idate5(4) = nint(allspot(5)) ! hour
         idate5(5) = nint(allspot(6)) ! minute
         idate5(6) = nint(allspot(7)) ! second		
		
         if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
             idate5(2) < 1    .or. idate5(2) >   12 .or. &
             idate5(3) < 1    .or. idate5(3) >   31 .or. &
             idate5(4) <0     .or. idate5(4) >   24 .or. &
             idate5(5) <0     .or. idate5(5) >   60 )then

            write(6,*)'READ_IASI:  ### ERROR IN READING ', 'IASI', ' BUFR DATA:', &
                ' STRANGE OBS TIME (YMDHM):', idate5(1:5)
             cycle read_loop
         end if
         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)		
         if ( obs_time < time_slots(0) .or.  &
           obs_time >= time_slots(num_fgat_time) ) cycle read_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat) ) exit
         end do	
         num_iasi_global = num_iasi_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1  

!    Observational info
         sat_zenang  = allspot(10)            ! satellite zenith angle
         ifov = nint(linele(1))               ! field of view

!    IASI fov ranges from 1 to 120.   Current angle dependent bias
!    correction has a maximum of 90 scan positions.   Geometry
!    of IASI scan allows us to remap 1-120 to 1-60.   Variable
!    ifovn below contains the remapped IASI fov.  This value is
!    passed on to and used in setuprad
         ifovn = (ifov-1)/2 + 1
         iscn = nint(linele(2))               ! scan line

!    Check field of view (FOVN) and satellite zenith angle (SAZA)
         if( ifov <= 0 .or. ifov > 120 .or. sat_zenang > 90._r_kind ) then
         write(unit=stdout,fmt=*) &
         'READ_IASI:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
                ' STRANGE OBS INFO(FOVN,SLNM,SAZA):', ifov, iscn, allspot(10)				
           cycle read_loop
         end if
         if ( ifov <= 60 ) sat_zenang = -sat_zenang

!    Compare IASI satellite scan angle and zenith angle
         piece = -step_adjust
         if ( mod(ifovn,2) == 1) piece = step_adjust
         lza = ((start + float((ifov-1)/4)*step) + piece)*deg2rad
         sat_height_ratio = (earth_radius + linele(5))/earth_radius
         lzaest = asin(sat_height_ratio*sin(lza))*rad2deg
         if (abs(sat_zenang - lzaest) > one) then
		 bad_line = iscn
          write(unit=stdout,fmt=*) 'abs(sat_zenang - lzaest) > one',abs(sat_zenang - lzaest)		   		   	   
           cycle read_loop
         end if
!   Clear Amount  (percent clear)
         call ufbrep(lnbufr,cloud_frac,1,6,iret,'FCPH')
         if (outside) cycle ! No good for this PE		
         num_iasi_local = num_iasi_local + 1
		 
         write(unit=info%date_char, &
         fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
         idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
         ':', idate5(5), ':', idate5(6)
         info%elv = 0.0  !aquaspot%selv		   
		 
         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad
            crit = 1. 
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_iasi_thinned=num_iasi_thinned+1
               cycle
            end if
         end if		
         call ufbrep(lnbufr,cscale,3,10,iret,'STCH ENCH CHSF')
         if(iret /= 10) then
           write(unit=stdout,fmt=*)  'READ_IASI  read scale error ',iret		   
           cycle read_loop
         end if

! The scaling factors are as follows, cscale(1) is the start channel number,
!                                     cscale(2) is the end channel number,
!                                     cscale(3) is the exponent scaling factor
! In our case (616 channels) there are 10 groups of cscale (dimension :: cscale(3,10))
!  The units are W/m2..... you need to convert to mW/m2.... (subtract 5 from cscale(3)
         do i=1,10  ! convert exponent scale factor to int and change units
             iexponent = -(nint(cscale(3,i)) - 5)
             sscale(i)=ten**iexponent
         end do

!    Read IASI channel number(CHNM) and radiance (SCRA)

         call ufbint(lnbufr,allchan,2,n_totchan,iret,'SCRA CHNM')
         if( iret /= n_totchan)then
              write(unit=stdout,fmt=*) &
			  'READ_IASI:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
                iret, ' CH DATA IS READ INSTEAD OF ',n_totchan				
           cycle read_loop
         end if
         num_iasi_used = num_iasi_used + 1		 
         nread(inst) = nread(inst) + 1								
         iskip = 0
         jstart=1
         do i=1,n_totchan
!     check that channel number is within reason
            if (( allchan(1,i) > zero .and. allchan(1,i) < 99999._r_kind) .and. &  ! radiance bounds
               (allchan(2,i) < 8462._r_kind .and. allchan(2,i) > zero )) then     ! chan # bounds
!         radiance to BT calculation
              radi = allchan(1,i)
              scaleloop: do j=jstart,10
                 if(allchan(2,i) >= cscale(1,j) .and. allchan(2,i) <= cscale(2,j))then
                    radi = allchan(1,i)*sscale(j)
                    jstart=j
                    exit scaleloop
                 end if
              end do scaleloop
              if (rtm_option == rtm_option_crtm) then
                 call CRTM_Planck_Temperature(inst,i,radi,temperature(i))
              else if (rtm_option == rtm_option_rttov) then
              end if
              if(temperature(i) < tbmin .or. temperature(i) > tbmax ) then
                 temperature(i) = min(tbmax,max(zero,temperature(i)))

              end if
            else           ! error with channel number or radiance
!             write(6,*)'READ_IASI:  iasi chan error',i,allchan(1,i), allchan(2,i)
              temperature(i) = min(tbmax,max(zero,temperature(i)))

            end if
         end do
         if(iskip > 0)write(6,*) ' READ_IASI : iskip > 0 ',iskip
!         if( iskip >= 10 )cycle read_loop 		
         data_all(5) = sat_zenang             ! satellite zenith angle (deg)
         data_all(6) = allspot(11)            ! satellite azimuth angle (deg)
         data_all(9) = allspot(12)            ! solar zenith angle (deg)
         data_all(10)= allspot(13)            ! solar azimuth angle (deg)
         do l=1,nchan
            data_all(l+nreal) = temperature(l)   ! brightness temerature
         end do
				
!  4.0   assign information to sequential radiance structure
!--------------------------------------------------------------------------
         allocate ( p % tb_inv (1:nchan) )
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = 1
         p%scanpos          = ifovn
         p%satzen           = data_all(5)
         p%satazi           = data_all(6)  ! look angle (deg) ! airsspot%bearaz
         p%solzen           = data_all(9)
         p%solazi           = data_all(10)
         p%tb_inv(1:nchan)  = data_all(nreal+1:nreal+nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat		

         allocate (p%next)   ! add next data
         p => p%next
         nullify (p%next) 		
     end do read_loop
  end do
  call closbf(lnbufr)
end do bufrfile

deallocate(data_all) ! Deallocate data arrays

  
   if (thinning .and. num_iasi_global > 0 ) then

      
      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do 
      end do
   
      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i)
            end do
         end do 
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i) = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_iasi_used_tmp = num_iasi_used
      do j = 1, num_iasi_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_iasi_thinned = num_iasi_thinned + 1
            num_iasi_used = num_iasi_used - 1
            nread(n) = nread(n) - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_iasi_used
   iv%total_rad_channel = iv%total_rad_channel + num_iasi_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_iasi_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_iasi_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_bufriasi.inc",536,message(1:1))
   endif

   write(unit=message(1),fmt='(a)') '   num_iasi_file num_iasi_global  num_iasi_local   num_iasi_used num_iasi_thinned'
   write(unit=message(2),fmt='(5(6x,i10))') num_iasi_file, num_iasi_global, num_iasi_local, num_iasi_used, num_iasi_thinned
   call da_message(message(1:2))

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
 

   do i = 1, iv%num_inst 
    
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad  = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad
      call da_allocate_rad_iv (i, nchan, iv)
	  
   end do

   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0
   p => head

   do n = 1, num_iasi_used
      i = p%sensor_index 
      nread(i) = nread(i) + 1 
  
      call da_initialize_rad_iv (i, nread(i), iv, p)
  
      current => p
      p => p%next
      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do
   deallocate ( p )
   deallocate (nread)
   deallocate (ptotal)

   call closbf(lnbufr)
   close(lnbufr)

!   call da_free_unit(lnbufr)

   if (trace_use) call da_trace_exit("da_read_obs_bufriasi")


end subroutine da_read_obs_bufriasi
subroutine da_read_obs_bufrseviri (obstype,iv,infile)	 

! subprogram:    read_seviri                  read bufr format seviri data
   !--------------------------------------------------------
   !  Purpose: read in NCEP bufr eos AIRS/AMSUA/HSB 1b data
   !            to innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read file twice
   !            so that da_scan_bufrairs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2013/03/26 - Creation            Hongli Wang 
   !
   !------------------------------------------------------------------------------

  implicit none

  real(r_kind)     :: POinT001 =   0.001_r_kind
  real(r_kind)     :: POinT01  =   0.01_r_kind
  real(r_kind)     :: TEN      =  10.0_r_kind
  real(r_kind)     :: R45      =  45.0_r_kind
  real(r_kind)     :: R60      =  60.0_r_kind
  real(r_kind)     :: R90      =  90.0_r_kind
  real(r_kind)     :: R180     = 180.0_r_kind
  real(r_kind)     :: R360     = 360.0_r_kind  

  character(9)      ,  intent (in)  :: obstype
  type (iv_type)    ,intent (inout) :: iv
  real(kind=8)    ::  obs_time  
  type (datalink_type), pointer  :: head, p, current, prev
  type(info_type)                :: info
  type(model_loc_type)           :: loc 
  type(model_loc_type)           :: loc_fov  

!!! for seviri
  character(80),  intent (in) :: infile

  character(8)       :: subset,subcsr,subasr
  character(80)      :: hdrsevi                             ! seviri header

  integer(i_kind)    :: nchanl,ilath,ilonh,ilzah,iszah,irec 
  integer(i_kind)    :: nhdr,nchn,ncld,nbrst !,idate,lnbufr
  integer(i_kind)    :: ireadmg,ireadsb,iret

  real(r_double),allocatable,dimension(:)   :: hdr         
  real(r_double),allocatable,dimension(:,:) :: datasev1,datasev2   

  logical            :: clrsky,allsky,allchnmiss
  real               :: rclrsky
  integer :: kidsat 
  integer(i_kind)    :: idate5(6)
  integer (i_kind), allocatable :: ptotal(:), nread(:) 
  integer(i_kind)    :: idate, im, iy, idd, ihh
!!! end for seviri

  ! Number of channels for sensors in 1
  integer(i_kind),parameter :: nchan = 8       
  integer(i_kind),parameter :: n_totchan  = 12 
  integer(i_kind),parameter :: maxinfo    =  33
  integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id
  real(r_kind)      :: crit
  integer(i_kind)   :: ifgat, iout, iobs
  logical           :: outside, outside_all, iuse,outside_fov

  integer           :: numbufr,ibufr,jj  
  logical :: found, head_found 
  real(r_kind)      :: step, start,step_adjust
  character(len=4)  :: senname
  integer(i_kind)   :: size,size_tmp
  character(10)     :: date
  real(r_kind)      :: sstime, tdiff, t4dv
  integer(i_kind)   :: nmind

  ! Other work variables
  real(r_kind)     :: clr_amt,piece
  real(r_kind)     :: dlon, dlat
  real(r_kind)     :: dlon_earth,dlat_earth,dlon_earth_deg,dlat_earth_deg
  real(r_kind)     :: lza, lzaest,sat_height_ratio
  real(r_kind)     :: pred, crit1, dist1
  real(r_kind)     :: sat_zenang
  real(r_kind)     :: radi
  real(r_kind)     :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr
  real(r_kind),dimension(0:4) :: rlndsea
  real(r_kind),dimension(0:3) :: sfcpct
  real(r_kind),dimension(0:3) :: ts
  real(r_kind),dimension(10) :: sscale
  real(r_kind),dimension(n_totchan) :: temperature
  real(r_kind),allocatable,dimension(:):: data_all
  real(r_kind) disterr,disterrmax,rlon00,rlat00

  logical          :: assim,valid
  logical          :: seviri 
  logical          :: data_on_edges,luse
  integer(i_kind)  :: nreal, ityp,ityp2, isflg
  integer(i_kind)  :: ifov, instr, ioff, ilat, ilon, sensorindex
  integer(i_kind)  :: num_seviri_file
  integer(i_kind)  :: num_seviri_local, num_seviri_global, num_seviri_used, num_seviri_thinned 
  integer(i_kind)  :: num_seviri_used_tmp  
  integer(i_kind)  :: i, j, l, iskip, ifovn, bad_line

  integer(i_kind)  :: itx, k, nele, itt, n
  integer(i_kind)  :: iexponent
  integer(i_kind)  :: idomsfc
  integer(i_kind)  :: ntest
  integer(i_kind)  :: error_status
  integer(i_kind)  :: num_bufr(7)

  integer          :: iost, lnbufr 
  character(20)    ::filename  
  real, allocatable :: in(:), out(:)

  ! Set standard parameters
  integer(i_kind),parameter:: ichan=-999  ! fov-based surface code is not channel specific for seviri 
  real(r_kind),parameter:: expansion=one  ! exansion factor for fov-based location code. 
  real(r_kind),parameter:: tbmin  = 50._r_kind
  real(r_kind),parameter:: tbmax  = 550._r_kind
  real(r_kind),parameter:: earth_radius = 6371000._r_kind

  ilath=8                      ! the position of latitude in the header
  ilonh=9                      ! the position of longitude in the header
  ilzah=10                     ! satellite zenith angle
  iszah=11                     ! solar zenith angle
  subcsr='NC021043'            ! sub message
  subasr='NC021042'            ! sub message

   if(trace_use) call da_trace_entry("da_read_obs_bufrseviri")

  !  0.0  Initialize variables
  !-----------------------------------
  sensor_id = 21 !seviri
  disterrmax=zero
  ntest=0
  nreal  = maxinfo
  seviri= obstype == 'seviri'

  bad_line=-1
  step  = 1.0
  start = 0.0
  step_adjust = 1.0_r_kind
  senname = 'SEVIRI'
  num_bufr(:)=0
  numbufr=0
  allocate(nread(1:rtminit_nsensor))
  allocate(ptotal(0:num_fgat_time))
  nread(1:rtminit_nsensor) = 0
  ptotal(0:num_fgat_time) = 0
  iobs = 0                 ! for thinning, argument is inout  
  num_seviri_file    = 0
  num_seviri_local   = 0
  num_seviri_global  = 0
  num_seviri_used    = 0
  num_seviri_thinned = 0  

  ! 1.0 Assign file names and prepare to read bufr files 
  !-------------------------------------------------------------------------

   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnbufr)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.bufr'
         open(unit   = lnbufr, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numbufr=numbufr+1
            num_bufr(numbufr)=i
         else
            close (lnbufr)
         end if
         call da_free_unit(lnbufr)
      end do
   else
     numbufr=1
   end if

 
   if (numbufr==0) numbufr=1
 
  bufrfile:  do ibufr=1,numbufr
   if (num_fgat_time==1) then
      filename=trim(infile)//'.bufr'
   else
      if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
         filename=trim(infile)//'.bufr'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_bufr(ibufr),'.bufr'
      end if
   end if 

   ! dont change, WRFDA uses specified units to read radiance data
   lnbufr = 92 
   open(unit=lnbufr,file=trim(filename),form='unformatted', &
      iostat = iost, status = 'old' ) !,convert='little_endian')
   if (iost /= 0) then
      call da_warning("da_read_obs_bufrseviri.inc",197, &
         (/"Cannot open file "//infile/))
      if(trace_use) call da_trace_exit("da_read_obs_bufrsevri")
      return
   end if

   ! Open 1 table
   call openbf(lnbufr,'IN',lnbufr) 
   call datelen(10)
   call readmg(lnbufr,subset,idate,iret)

   ! Check the data set
  if( iret/=0) then
     write(message(1),fmt='(A)')'SKIP PROCESSING OF SEVIRI FILE'
     write(message(2),fmt='(A,I2,A)')'infile = ', lnbufr, infile
     call da_warning("da_read_obs_bufrseviri.inc",212,message(1:2))
     if(trace_use) call da_trace_exit("da_read_obs_bufrseviri")
     return 
  endif

  clrsky=.false.
  allsky=.false.
  if(subset == subcsr) then
     clrsky=.true.
     write(message(1),fmt='(A)')'READ SEVIRI FILE'
     write(message(2),fmt='(A,L1,2A)')'clrsky= ', clrsky,' subset= ', subset
     call da_message(message(1:2))
  elseif(subset == subasr) then
     allsky=.true.
     write(message(1),fmt='(A)')'READ SEVIRI FILE'
     write(message(2),fmt='(A,L1,2A)')'allsky= ', allsky,' subset= ', subset
     call da_message(message(1:2))
  else
     write(message(1),fmt='(A)')'SKIP PROCESSING OF UNKNOWN SEVIRI FILE'
     write(message(2),fmt='(A,I2,3A)')'infile=', lnbufr, infile,' subset=', subset
     call da_warning("da_read_obs_bufrseviri.inc",232,message(1:2))
     if(trace_use) call da_trace_exit("da_read_obs_bufrseviri")
     return 
  endif

  ! Set 1 string based on seviri data set
  if (clrsky) then
     hdrsevi='SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH SAZA SOZA'
     nhdr=11
     nchn=12
     ncld=nchn
     nbrst=nchn
  else if (allsky) then
     hdrsevi='SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH'
     nhdr=9
     nchn=11
     ncld=2
     nbrst=nchn*6                ! channel dependent: all, clear, cloudy, low,
                                 !middle and high clouds
  endif
  allocate(datasev1(1,ncld))     ! not channel dependent
  allocate(datasev2(1,nbrst))    ! channel dependent: all, clear, cloudy, low,
                                 !middle and high clouds
  allocate(hdr(nhdr))



   iy=0
   im=0
   idd=0
   ihh=0

   sensorindex=1  
   write(unit=date,fmt='( i10)') idate
   read(unit=date,fmt='(i4,3i2)') iy,im,idd,ihh
   write(unit=stdout,fmt='(a,4i4,2x,a)') &
      'Bufr file date is ',iy,im,idd,ihh,trim(infile)

   ! 2.0 Loop to read bufr file and assign information to a sequential structure
   !-------------------------------------------------------------------------

   ! Allocate arrays to hold data
   nele=nreal+nchan
   allocate(data_all(nele))
   if ( ibufr == 1 ) then
      allocate (head)
      nullify  ( head % next )
      p => head
   end if
  

! Big loop to read data file

  do while(ireadmg(lnbufr,subset,idate)>=0)  

     read_loop: do while (ireadsb(lnbufr)==0)
         num_seviri_file = num_seviri_file + 1	 

    ! Read SEVIRI information
         call ufbint(lnbufr,hdr,nhdr,1,iret,hdrsevi)				

        kidsat = nint(hdr(1))
        ! SAID 55 is meteosat-8  or msg-1
        ! SAID 56 is meteosat-9  or msg-2
        ! SAID 57 is meteosat-10 or msg-3
        if ( ( kidsat > 57) .or. ( kidsat < 55) ) then 
           write(unit=message(1),fmt='(A,I6)') 'Unknown platform: ', kidsat
           call da_warning("da_read_obs_bufrseviri.inc",299,message(1:1))
        end if
        platform_id  = 12  ! MSG - Meteosat Second Generation
        if ( kidsat == 55 ) then
            satellite_id = 1
        else if ( kidsat == 56 ) then
            satellite_id = 2
        else if ( kidsat == 57 ) then
            satellite_id = 3
        end if

        if (clrsky) then     ! asr bufr has no sza
        ! remove the obs whose satellite zenith angles larger than 65 degree
           if ( hdr(ilzah) > 65.0 ) cycle read_loop
        end if

        call ufbint(lnbufr,datasev1,1,ncld,iret,'NCLDMNT')

        if(iret /= 1) cycle read_loop
        do n=1,ncld
           if(datasev1(1,n)>= 0.0 .and. datasev1(1,n) <= 100.0 ) then
              rclrsky=datasev1(1,n)
               ! first QC filter out data with less clear sky fraction
               if ( rclrsky < 70.0 ) cycle read_loop
               !if ( rclrsky < 90.0 ) cycle read_loop
           end if
        end do

        call ufbrep(lnbufr,datasev2,1,nbrst,iret,'TMBRST')
        
        ! Check if data of channel 1-11 are missing

        allchnmiss=.true.
        do n=1,11
           if(datasev2(1,n)<500.)  then
              allchnmiss=.false.
           end if
        end do
        if(allchnmiss) cycle read_loop

        ! Check observing position
         info%lat  =  hdr(ilath)           ! latitude
         info%lon  =  hdr(ilonh)           ! longitude
         if( abs(info%lat) > R90  .or. abs(info%lon) > R360 .or. &
         (abs(info%lat) == R90 .and. info%lon /= ZERO) )then
         write(unit=stdout,fmt=*) &
         'READ_SEVIRI:  ### ERROR IN READING ', senname, ' BUFR DATA:', &
               ' STRANGE OBS POINT (LAT,LON):', info%lat, info%lon
            cycle read_loop
         end if		 
	 		 
         call da_llxy (info, loc, outside, outside_all)	
         if (outside_all) cycle
	     inst = 0	 
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do	 
         if (inst == 0) cycle read_loop		 
		 
        ! Check obs time
         idate5(1) = nint(hdr(2)) ! year
         idate5(2) = nint(hdr(3)) ! month
         idate5(3) = nint(hdr(4)) ! day
         idate5(4) = nint(hdr(5)) ! hour
         idate5(5) = nint(hdr(6)) ! minute
         idate5(6) = nint(hdr(7)) ! second		
		
         if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
             idate5(2) < 1    .or. idate5(2) >   12 .or. &
             idate5(3) < 1    .or. idate5(3) >   31 .or. &
             idate5(4) < 0    .or. idate5(4) >   24 .or. &
             idate5(5) < 0    .or. idate5(5) >   60 ) then

            write(message(1),fmt='(A)')'ERROR IN READING SEVIRI BUFR DATA'
            write(message(2),fmt='(A,5I8)')'STRANGE OBS TIME (YMDHM):', idate5(1:5)
            call da_warning("da_read_obs_bufrseviri.inc",379,message(1:2))
            cycle read_loop
         end if

         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)		
         if ( obs_time < time_slots(0) .or.  &
           obs_time >= time_slots(num_fgat_time) ) cycle read_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat) ) exit
         end do	
         num_seviri_global = num_seviri_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1  

         if (outside) cycle ! No good for this PE		
         num_seviri_local = num_seviri_local + 1
		 
         write(unit=info%date_char, &
         fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
         idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
         ':', idate5(5), ':', idate5(6)
         info%elv = 0.0  !aquaspot%selv		   
		 
         ! 3.0  Make Thinning
         ! Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad
            crit = 1. 
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_seviri_thinned=num_seviri_thinned+1
               cycle
            end if
         end if		

         ! data SEVIRI channel number(CHNM) and radiance (SCRA)

         num_seviri_used = num_seviri_used + 1		 
         nread(inst) = nread(inst) + 1								
         iskip = 0
         do i=1,n_totchan
            ! check that tb is within reasonal bound
            if ( datasev2(1,i) < tbmin .or. datasev2(1,i) > tbmax ) then
               temperature(i) = missing_r
            else 
               temperature(i) = datasev2(1,i)
            end if
         end do

         if(iskip > 0)write(6,*) ' READ_SEVIRI : iskip > 0 ',iskip

         do l=1,nchan
            data_all(l+nreal) = temperature(l+3)   ! brightness temerature
         end do
				
    ! 4.0 assign information to sequential radiance structure
    !--------------------------------------------------------------------------
    !        iscan = nint(hdr(ilzah))+1.001_r_kind 
         allocate ( p % tb_inv (1:nchan ))
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = 1
         p%scanpos          = nint(hdr(ilzah))+1.001_r_kind 
         p%satzen           = hdr(ilzah)               ! satellite zenith angle (deg)
         p%satazi           = 0.0                      ! satellite azimuth angle (deg) 
         p%solzen           = 0.0                      ! solar zenith angle (deg)  
         p%solazi           = 0.0                      ! solar azimuth angle (deg) 
         p%tb_inv(1:nchan)  = data_all(nreal+1:nreal+nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat		

         allocate (p%next)   ! add next data
         p => p%next
         nullify (p%next) 		
     end do read_loop
  end do
  call closbf(lnbufr)
end do bufrfile

deallocate(data_all) ! Deallocate data arrays

  
   if (thinning .and. num_seviri_global > 0 ) then

      
      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do 
      end do
   
      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i)
            end do
         end do 
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i) = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )


      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_seviri_used_tmp = num_seviri_used
      do j = 1, num_seviri_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_seviri_thinned = num_seviri_thinned + 1
            num_seviri_used = num_seviri_used - 1
            nread(n) = nread(n) - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_seviri_used
   iv%total_rad_channel = iv%total_rad_channel + num_seviri_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_seviri_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_seviri_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_bufrseviri.inc",574,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') '   num_seviri_file num_seviri_global  num_seviri_local   num_seviri_used num_seviri_thinned'
   write(stdout,'(5(8x,i10))') num_seviri_file, num_seviri_global, num_seviri_local, num_seviri_used, num_seviri_thinned

  

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
 

   do i = 1, iv%num_inst 
    
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad  = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad
      call da_allocate_rad_iv (i, nchan, iv)
	  
   end do

   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0
   p => head

   do n = 1, num_seviri_used
      i = p%sensor_index 
      nread(i) = nread(i) + 1 
  
      call da_initialize_rad_iv (i, nread(i), iv, p)
  
      current => p
      p => p%next
      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do
   deallocate ( p )
   deallocate (nread)
   deallocate (ptotal)

   call closbf(lnbufr)
   close(lnbufr)

   call da_free_unit(lnbufr)

   if(trace_use) call da_trace_exit("da_read_obs_bufrseviri")


end subroutine da_read_obs_bufrseviri
subroutine da_read_obs_hdf5amsr2 (iv, infile_tb,infile_clw)
   !--------------------------------------------------------
   !  Purpose: read in JAXA AMSR2 Level-1R data in 1 format
   !           and form innovation structure
   !
   !   METHOD: use F90 sequantial data structure to avoid read the file twice
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2013/10/22 - Creation         Syed RH Rizvi, NCAR/NESL/MMM/DAS
   !           2014         Modification     Chun Yang
   !------------------------------------------------------------------------------

   implicit none

   character(len=*), intent(in)    :: infile_tb, infile_clw
   type(iv_type),    intent(inout) :: iv

! fixed parameter values
   integer,parameter::max_scan=2200     ! Maximum allowed NumberOfScans
   integer,parameter::ovr=20            ! Number of OverlapScans
   integer,parameter::hi_rez_fov=486    ! high resolution pixel width
   integer,parameter::lo_rez_fov=243    ! low  resolution pixel width
   integer,parameter::time_dims=6       ! Time dimension
   integer,parameter::nfile_max = 8     ! each hdf file contains ~50min of data
                                        ! at most 8 files for a 6-h time window
! interface variable
   integer iret                         ! return status
   integer(HID_T) fhnd1, fhnd2, fhnd3   ! file handle
   integer(HID_T) ahnd1, ahnd2, ahnd3   ! attribute handle
   integer(HID_T) dhnd1, dhnd2, dhnd3   ! dataset handle
   integer(HID_T) shnd1, shnd2, shnd3   ! dataspace handle
   integer(HSIZE_T) sz1(3)              ! array size 1
   integer(HSIZE_T) sz2(3)              ! array size 2

   integer(4) :: nscan                       ! NumberOfScans
   real(4)    :: sca                         ! ScaleFactor
   real(8)    :: r8d1(max_scan)              ! array buffer for real(8) D1
   integer(4) :: i4d2hi(hi_rez_fov,max_scan) ! array buffer for integer(4) D2 high

   type AM2_COMMON_SCANTIME
      real(8)    tai93sec
      integer(2) year
      integer(2) month
      integer(2) day
      integer(2) hour
      integer(2) minute
      integer(2) second
      integer(2) ms
      integer(2) reserve
   endtype

! array data
   type(AM2_COMMON_SCANTIME) st(max_scan)  ! scantime
   real(4) :: lat89ar(hi_rez_fov,max_scan) ! lat for 89a altitude revised
   real(4) :: latlr(lo_rez_fov,max_scan)   ! lat for low resolution
   real(4) :: lon89ar(hi_rez_fov,max_scan) ! lon for 89a altitude revised
   real(4) :: lonlr(lo_rez_fov,max_scan)   ! lon for low resolution

   real(4) :: tb89ah(hi_rez_fov,max_scan)  ! tb for 89ah
   real(4) :: tb89av(hi_rez_fov,max_scan)  ! tb for 89av

   integer(4) :: loflo(lo_rez_fov,max_scan*4) ! land ocean flag for low
   integer(1) :: lof06(lo_rez_fov,max_scan)   ! land ocean flag for 06
   integer(1) :: lof10(lo_rez_fov,max_scan)   ! land ocean flag for 10
   integer(1) :: lof23(lo_rez_fov,max_scan)   ! land ocean flag for 23
   integer(1) :: lof36(lo_rez_fov,max_scan)   ! land ocean flag for 36

   integer(4) :: lofhi(hi_rez_fov,max_scan*2) ! land ocean flag for high
   integer(1) :: lof89a(hi_rez_fov,max_scan)  ! land ocean flag for 89a
   integer(1) :: lof89b(hi_rez_fov,max_scan)  ! land ocean flag for 89b

   real(4)    :: ear_in(lo_rez_fov,max_scan)  ! earth incidence
   real(4)    :: ear_az(lo_rez_fov,max_scan)  ! earth azimuth
   real(4)    :: clw(lo_rez_fov,max_scan)     ! obs retrieved cloud liquid water

   real(4)    :: sun_az(lo_rez_fov,max_scan)  ! sun_azimuth
   real(4)    :: sun_el(lo_rez_fov,max_scan)  ! sun_elevation
   real(4)    :: sun_zen(lo_rez_fov,max_scan) ! sun_zenith

   real(r_kind)            :: R90    = 90.0_r_kind
   real(r_kind),parameter  :: tbmin  = 50._r_kind
   real(r_kind),parameter  :: tbmax  = 550._r_kind

   real(4)          :: sat_zenith(lo_rez_fov,max_scan)
   real(4)          :: sat_azimuth(lo_rez_fov,max_scan)

   real(kind=8)                   :: obs_time
   type (datalink_type),pointer   :: head, p, current, prev
   type(info_type)                :: info
   type(model_loc_type)           :: loc

   integer(i_kind)    :: idate5(6)

   integer(i_kind)   :: inst,platform_id,satellite_id,sensor_id
   real(r_kind)      :: tb, crit
   integer(i_kind)   :: ifgat, iout, iobs
   logical           :: outside, outside_all, iuse

   integer           :: i,j,k,l,m,n, ifile, landsea_mask
   logical           :: found, head_found, head_allocated

! Other work variables
   real(r_kind)     :: dlon_earth,dlat_earth
   integer(i_kind)  :: num_amsr2_local, num_amsr2_global, num_amsr2_used, num_amsr2_thinned
   integer(i_kind)  :: num_amsr2_used_tmp, num_amsr2_file
   integer(i_kind)  :: num_amsr2_local_local, num_amsr2_global_local, num_amsr2_file_local
   integer(i_kind)  :: itx, itt
   character(80)    :: filename1, filename2
   integer          :: nchan,ifov,iscan,ichannels
   integer          :: nfile
   character(80)    :: fname_tb(nfile_max)
   character(80)    :: fname_clw(nfile_max)
   logical          :: fexist, got_clw_file

! Allocatable arrays
   integer(i_kind),allocatable  :: ptotal(:)
   real,allocatable             :: in(:), out(:)
   real(r_kind),allocatable     :: data_all(:)

   real,allocatable             :: obstime(:,:)

   real(r_kind)    :: sun_zenith, sun_azimuth 

   integer,parameter  :: num_low_freq_chan=12
   real(4)            :: tb_low(lo_rez_fov,max_scan,num_low_freq_chan)
   character(len=80) tb_low_name(num_low_freq_chan)
   data tb_low_name/'Brightness Temperature (res06,6.9GHz,V)','Brightness Temperature (res06,6.9GHz,H)',&
                    'Brightness Temperature (res06,7.3GHz,V)','Brightness Temperature (res06,7.3GHz,H)',&
                    'Brightness Temperature (res10,10.7GHz,V)','Brightness Temperature (res10,10.7GHz,H)',&
                    'Brightness Temperature (res23,18.7GHz,V)','Brightness Temperature (res23,18.7GHz,H)',&
                    'Brightness Temperature (res23,23.8GHz,V)','Brightness Temperature (res23,23.8GHz,H)',&
                    'Brightness Temperature (res36,36.5GHz,V)','Brightness Temperature (res36,36.5GHz,H)'/

   if (trace_use) call da_trace_entry("da_read_obs_hdf5amsr2")

!  0.0  Initialize variables
!-----------------------------------
   head_allocated = .false.
   platform_id  = 29  ! Table-2 Col 1 corresponding to 'gcom-w'
   satellite_id = 1   ! Table-2 Col 3
   sensor_id    = 63  ! Table-3 Col 2 corresponding to 'amsr2'

   allocate(ptotal(0:num_fgat_time))
   ptotal(0:num_fgat_time) = 0
   iobs = 0                 ! for thinning, argument is inout
   num_amsr2_file    = 0
   num_amsr2_local   = 0
   num_amsr2_global  = 0
   num_amsr2_used    = 0
   num_amsr2_thinned = 0

   do i = 1, rtminit_nsensor
      if (platform_id  == rtminit_platform(i) &
          .and. satellite_id == rtminit_satid(i)    &
          .and. sensor_id    == rtminit_sensor(i)) then
         inst = i
         exit
      end if
   end do
   if (inst == 0) then
      call da_warning("da_read_obs_hdf5amsr2.inc",165, &
          (/"The combination of Satellite_Id and Sensor_Id for AMSR-2 is not found"/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5amsr2")
      return
   end if

! Initialize 1 library and Fortran90 interface
   call H5open_f(iret)
   if(iret.lt.0)then
      call da_warning("da_read_obs_hdf5amsr2.inc",174, &
           (/"Problems in Initializing HDF5 library. Can not read AMSR-2 HDF5 data. "/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5amsr2")
      return
   endif

   nchan = iv%instid(inst)%nchan
   write(unit=stdout,fmt=*)'AMSR2 nchan: ',nchan
   allocate(data_all(1:nchan))

! 1.0 Assign file names and prepare to read amsr2 files
!-------------------------------------------------------------------------
   nfile       = 0  !initialize
   fname_tb(:) = '' !initialize
   ! first check if L1SGRTBR.h5 is available
   filename1 = trim(infile_tb)//'.h5'
   filename2 = trim(infile_clw)//'.h5'
   inquire (file=filename1, exist=fexist)
   if ( fexist ) then
      nfile = 1
      fname_tb(nfile)  = filename1
      fname_clw(nfile) = filename2
   else
      ! check if L1SGRTBR-0x.h5 is available for multiple input files
      ! here 0x is the input file sequence number
      ! do not confuse it with fgat time slot index
      do i = 1, nfile_max
         write(filename1,fmt='(A,A,I2.2,A)') trim(infile_tb),'-',i,'.h5'
         write(filename2,fmt='(A,A,I2.2,A)') trim(infile_clw),'-',i,'.h5'
         inquire (file=filename1, exist=fexist)
         if ( fexist ) then
            nfile = nfile + 1
            fname_tb(nfile)  = filename1
            fname_clw(nfile) = filename2
         else
            exit
         end if
      end do
   end if

   if ( nfile == 0 ) then
      call da_warning("da_read_obs_hdf5amsr2.inc",215, &
         (/"No valid AMSR-2 L1SGRTBR.h5 or L1SGRTBR-01.h5 file found."/))
      if (trace_use) call da_trace_exit("da_read_obs_hdf5amsr2")
      return
   end if

   ! Check to see if leap second file exists for graceful failure
      inquire( file='leapsec.dat', exist=fexist )
      if (.not. fexist) call da_error("da_read_obs_hdf5amsr2.inc",223, &
           (/'Can not find leapsec.dat for AMSR2 data: copy or link from WRFDA/var/run'/))

   infile_loop:  do ifile = 1, nfile
      num_amsr2_file_local    = 0
      num_amsr2_local_local   = 0
      num_amsr2_global_local  = 0

   ! open 1 file for read
      call H5Fopen_f(fname_tb(ifile),H5F_ACC_RDONLY_F,fhnd1,iret,H5P_DEFAULT_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",234, &
              (/"Cannot open HDF5 file "//trim(fname_tb(ifile))/))
         cycle infile_loop
      endif
      got_clw_file = .false.
      call H5Fopen_f(fname_clw(ifile),H5F_ACC_RDONLY_F,fhnd2,iret,H5P_DEFAULT_F)
      if ( iret == 0 ) then
         got_clw_file = .true.
      endif
      ! to do: when got_clw_file=.true., need to check GranuleID for consistency
      ! betweee tb and clw files

   ! calculate NumberOfScans from array size and OverlapScans
      call H5Dopen_f(fhnd1,'Scan Time',dhnd1,iret)
      call H5Dget_space_f(dhnd1,shnd1,iret)
      call H5Sget_simple_extent_dims_f(shnd1,sz1,sz2,iret)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",251, &
             (/"HDF5 read problem for: Scan Time"/))
      endif
      call H5Sclose_f(shnd1,iret)
      call H5Dclose_f(dhnd1,iret)

      nscan=sz1(1)-ovr*2

      write(unit=stdout,fmt=*)'NumberOfScans(RETRIEVE BY ARRAY SIZE): ',nscan
      write(unit=stdout,fmt=*)'OverlapScans(FIXED VALUE): ',ovr

   ! check limit
      if(nscan.gt.max_scan)then
         write(unit=stdout,fmt=*)'limit of NumberOfScans = ',max_scan
         call da_warning("da_read_obs_hdf5amsr2.inc",265, &
              (/"HDF5 lmit error for: max_scan"/))
      endif

   ! read array: scantime
   ! read
      call H5Dopen_f(fhnd1,'Scan Time',dhnd1,iret)
      sz1(1)=max_scan
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_DOUBLE,r8d1,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",276, &
             (/"HDF5 read error for: Scan Time"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap
      do j=1,nscan
         r8d1(j)=r8d1(j+ovr)
      enddo
      do j=nscan+1,max_scan
         r8d1(j)=0
      enddo
   ! convert
      call amsr2time(nscan,r8d1,st)
   ! sample display
      allocate  (obstime(1:time_dims,1:nscan))  ! year, month, day, hour, min, sec
      do j = 1, nscan
         obstime(1,j) = st(j)%year
         obstime(2,j) = st(j)%month
         obstime(3,j) = st(j)%day
         obstime(4,j) = st(j)%hour
         obstime(5,j) = st(j)%minute
         obstime(6,j) = st(j)%second
      end do
      write(unit=stdout,fmt=*)'time(scan=1) year: ',st(1)%year,' month:',st(1)%month,' day: ',st(1)%day,&
         ' hour: ',st(1)%hour,' minute: ',st(1)%minute,' second: ',st(1)%second

   ! read array: latlon for 89a altitude revised
   ! read lat
      call H5Dopen_f(fhnd1, &
         'Latitude of Observation Point for 89A',dhnd1,iret)
      sz1(1)=max_scan
      sz1(2)=hi_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,lat89ar,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",311, &
             (/"HDF5 read error for: Latitude of Observation Point for 89A"/))
      endif
      call H5Dclose_f(dhnd1,iret)

   ! read lon
      call H5Dopen_f(fhnd1, &
         'Longitude of Observation Point for 89A',dhnd1,iret)
      sz1(1)=max_scan
      sz1(2)=hi_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,lon89ar,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
          call da_warning("da_read_obs_hdf5amsr2.inc",324, &
              (/"HDF5 read error for: Longitude of Observation Point for 89A"/))
          call da_trace_exit("da_read_obs_hdf5amsr2")
      endif
      call H5Dclose_f(dhnd1,iret)

   ! cutoff overlap
      do j=1,nscan
         lat89ar(:,j)=lat89ar(:,j+ovr)
         lon89ar(:,j)=lon89ar(:,j+ovr)
      enddo
      do j=nscan+1,max_scan
         lat89ar(:,j)=0
         lon89ar(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'latlon89ar(pixel=1,scan=1): ',lat89ar(1,1),lon89ar(1,1)

   ! read array: latlon for low resolution
      do j=1,nscan
         do i=1,lo_rez_fov
            latlr(i,j)=lat89ar(i*2-1,j)
            lonlr(i,j)=lon89ar(i*2-1,j)
         enddo
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)&
      !   'latlonlr(pixel=1,scan=1): ',latlr(1,1),lonlr(1,1)

   ! read array: tb for low frequency channels
      do k=1,num_low_freq_chan
         call H5Dopen_f(fhnd1,tb_low_name(k),dhnd1,iret)
         call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)
         call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
         call H5Aclose_f(ahnd1,iret)
         sz1(1)=max_scan
         sz1(2)=lo_rez_fov
         call H5Dread_f(dhnd1, &
            H5T_NATIVE_REAL,tb_low(:,:,k),sz1,iret,H5S_ALL_F,H5S_ALL_F)
         if(iret.lt.0)then
            call da_warning("da_read_obs_hdf5amsr2.inc",364, &
                  (/"HDF5 read error for: Brightness Temperature"/))
         endif
         call H5Dclose_f(dhnd1,iret)
      ! cutoff overlap & convert to unsignd & change scale
         do j=1,nscan
            do i=1,lo_rez_fov
               tb_low(i,j,k)=tb_low(i,j+ovr,k)
               if(tb_low(i,j,k).lt.65534)tb_low(i,j,k)=tb_low(i,j,k)*sca
            enddo
         enddo
         do j=nscan+1,max_scan
            tb_low(:,j,:)=0
         enddo
      ! sample display
         if (print_detail_rad) then
            write(unit=message(1),fmt='(A,I6,A,F10.4)')&
               'tb_low(pixel=1,scan=1,chan=',k,'): ',tb_low(1,1,k)
            call da_message(message(1:1))
         endif
      enddo

   ! read array: tb for 89ah
   ! read
      call H5Dopen_f(fhnd1, &
         'Brightness Temperature (original,89GHz-A,H)',dhnd1,iret)
      call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)    ! get scale
      call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
      call H5Aclose_f(ahnd1,iret)
      sz1(1)=max_scan
      sz1(2)=hi_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,tb89ah,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",398, &
            (/"HDF5 read error for: Brightness Temperature (original,89GHz-A,H)"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap & convert to unsignd & change scale
      do j=1,nscan
         do i=1,hi_rez_fov
            tb89ah(i,j)=tb89ah(i,j+ovr)
            if(tb89ah(i,j).lt.65534)tb89ah(i,j)=tb89ah(i,j)*sca
         enddo
      enddo
      do j=nscan+1,max_scan
         tb89ah(:,j)=0
      enddo
   ! sample display
         if (print_detail_rad) then
            write(unit=message(1),fmt='(A,F10.4)')&
               'tb89ah(pixel=1,scan=1): ',tb89ah(1,1)
            call da_message(message(1:1))
         endif

   ! read array: tb for 89av
   ! read
      call H5Dopen_f(fhnd1, &
         'Brightness Temperature (original,89GHz-A,V)',dhnd1,iret)
      call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)    ! get scale
      call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
      call H5Aclose_f(ahnd1,iret)
      sz1(1)=max_scan
      sz1(2)=hi_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,tb89av,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",431, &
            (/"HDF5 read error for: Brightness Temperature (original,89GHz-A,V)"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap & convert to unsignd & change scale
      do j=1,nscan
         do i=1,hi_rez_fov
            tb89av(i,j)=tb89av(i,j+ovr)
            if(tb89av(i,j).lt.65534)tb89av(i,j)=tb89av(i,j)*sca
         enddo
      enddo
      do j=nscan+1,max_scan
         tb89av(:,j)=0
      enddo
   ! sample display
         if (print_detail_rad) then
            write(unit=message(1),fmt='(A,F10.4)')&
               'tb89av(pixel=1,scan=1): ',tb89av(1,1)
            call da_message(message(1:1))
         endif

  ! read array: land ocean flag for low
  ! read
      call H5Dopen_f(fhnd1, &
         'Land_Ocean Flag 6 to 36',dhnd1,iret)
      sz1(1)=max_scan*6
      sz1(2)=lo_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_INTEGER,loflo,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",461, &
            (/"HDF5 read error for: Land_Ocean Flag 6 to 36"/))
      endif
      call H5Dclose_f(dhnd1,iret)
  ! separate
      do j=1,nscan+ovr*2
         do i=1,lo_rez_fov
            lof06(i,j)=loflo(i,(nscan+ovr*2)*0+j)
            lof10(i,j)=loflo(i,(nscan+ovr*2)*1+j)
            lof23(i,j)=loflo(i,(nscan+ovr*2)*2+j)
            lof36(i,j)=loflo(i,(nscan+ovr*2)*3+j)
         enddo
      enddo
  ! cutoff overlap
      do j=1,nscan
         do i=1,lo_rez_fov
            lof06(i,j)=lof06(i,j+ovr)
            lof10(i,j)=lof10(i,j+ovr)
            lof23(i,j)=lof23(i,j+ovr)
            lof36(i,j)=lof36(i,j+ovr)
         enddo
      enddo
      do j=nscan+1,max_scan
         lof06(:,j)=0
         lof10(:,j)=0
         lof23(:,j)=0
         lof36(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'lof06(pixel=1,scan=1): ',lof06(1,1)
      !write(unit=stdout,fmt=*)'lof10(pixel=1,scan=1): ',lof10(1,1)
      !write(unit=stdout,fmt=*)'lof23(pixel=1,scan=1): ',lof23(1,1)
      !write(unit=stdout,fmt=*)'lof36(pixel=1,scan=1): ',lof36(1,1)

   ! read array: land ocean flag for high
   ! read
      call H5Dopen_f(fhnd1, &
         'Land_Ocean Flag 89',dhnd1,iret)
      sz1(1)=max_scan*2
      sz1(2)=hi_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_INTEGER,lofhi,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",504, &
            (/"HDF5 read error for: Land_Ocean Flag 89"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! separate
      do j=1,nscan+ovr*2
         do i=1,hi_rez_fov
            lof89a(i,j)=lofhi(i,(nscan+ovr*2)*0+j)
            lof89b(i,j)=lofhi(i,(nscan+ovr*2)*1+j)
         enddo
      enddo
      do j=1,nscan
         do i=1,hi_rez_fov
            lof89a(i,j)=lof89a(i,j+ovr)
            lof89b(i,j)=lof89b(i,j+ovr)
         enddo
      enddo
      do j=nscan+1,max_scan
         lof89a(:,j)=0
         lof89b(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'lof89a(pixel=1,scan=1): ',lof89a(1,1)
      !write(unit=stdout,fmt=*)'lof89b(pixel=1,scan=1): ',lof89b(1,1)

   ! read array: earth incidence
   ! read
      call H5Dopen_f(fhnd1, &
         'Earth Incidence',dhnd1,iret)
      call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)     ! get scale
      call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
      call H5Aclose_f(ahnd1,iret)
      sz1(1)=max_scan
      sz1(2)=lo_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,ear_in,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",541, &
            (/"HDF5 read error for: Earth Incidence"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap & change scale
      do j=1,nscan
         do i=1,lo_rez_fov
            ear_in(i,j)=ear_in(i,j+ovr)
            if(ear_in(i,j).gt.-32767)ear_in(i,j)=ear_in(i,j)*sca
         enddo
      enddo
      do j=nscan+1,max_scan
         ear_in(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'ear_in(pixel=1,scan=1): ',ear_in(1,1)

   ! read array: earth azimuth
   ! read
      call H5Dopen_f(fhnd1, &
         'Earth Azimuth',dhnd1,iret)
      call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)    ! get scale
      call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
      call H5Aclose_f(ahnd1,iret)
      sz1(1)=max_scan
      sz1(2)=lo_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,ear_az,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",570, &
            (/"HDF5 read error for: Earth Azimuth"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap & change scale
      do j=1,nscan
         do i=1,lo_rez_fov
            ear_az(i,j)=ear_az(i,j+ovr)
            if(ear_az(i,j).gt.-32767)ear_az(i,j)=ear_az(i,j)*sca
         enddo
      enddo
      do j=nscan+1,max_scan
         ear_az(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'ear_az(pixel=1,scan=1): ',ear_az(1,1)

   ! read array: sun azimuth
   ! read
      call H5Dopen_f(fhnd1, &
         'Sun Azimuth',dhnd1,iret)
      call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)    ! get scale
      call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
      call H5Aclose_f(ahnd1,iret)
      sz1(1)=max_scan
      sz1(2)=lo_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,sun_az,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",599, &
            (/"HDF5 read error for: Sun Azimuth"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap & change scale
      do j=1,nscan
         do i=1,lo_rez_fov
            sun_az(i,j)=sun_az(i,j+ovr)
           if(sun_az(i,j).gt.-32767)sun_az(i,j)=sun_az(i,j)*sca
         enddo
      enddo
      do j=nscan+1,max_scan
         sun_az(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'sun_az(pixel=1,scan=1): ',sun_az(1,1)

   ! read array: sun elevation
   ! read
      call H5Dopen_f(fhnd1, &
         'Sun Elevation',dhnd1,iret)
      call H5Aopen_f(dhnd1,'SCALE FACTOR',ahnd1,iret)    ! get scale
      call H5Aread_f(ahnd1,H5T_NATIVE_REAL,sca,sz1,iret)
      call H5Aclose_f(ahnd1,iret)
      sz1(1)=max_scan
      sz1(2)=lo_rez_fov
      call H5Dread_f(dhnd1, &
         H5T_NATIVE_REAL,sun_el,sz1,iret,H5S_ALL_F,H5S_ALL_F)
      if(iret.lt.0)then
         call da_warning("da_read_obs_hdf5amsr2.inc",628, &
            (/"HDF5 read error for: Sun Elevation"/))
      endif
      call H5Dclose_f(dhnd1,iret)
   ! cutoff overlap & change scale
      do j=1,nscan
         do i=1,lo_rez_fov
            sun_el(i,j)=sun_el(i,j+ovr)
            if(sun_el(i,j).gt.-32767)sun_el(i,j)=sun_el(i,j)*sca
         enddo
      enddo
      do j=nscan+1,max_scan
         sun_el(:,j)=0
      enddo
   ! sample display
      !write(unit=stdout,fmt=*)'sun_el(pixel=1,scan=1): ',sun_el(1,1)
      sun_zen(:,:)=R90-sun_el(:,:)
      sat_zenith(:,:)=ear_in(:,:)
      sat_azimuth(:,:)=ear_az(:,:)

   ! close file and 1
      call H5Fclose_f(fhnd1,iret)

      if ( got_clw_file ) then
      ! read CLW from infile_clw:
         call H5Dopen_f(fhnd2,'Geophysical Data',dhnd2,iret)
         call H5Aopen_f(dhnd2,'SCALE FACTOR',ahnd2,iret)
         call H5Aread_f(ahnd2,H5T_NATIVE_REAL,sca,sz1,iret)
         call H5Aclose_f(ahnd2,iret)
         sz1(1)=max_scan
         sz1(2)=lo_rez_fov
         call H5Dread_f(dhnd2, &
            H5T_NATIVE_REAL,clw,sz1,iret,H5S_ALL_F,H5S_ALL_F)
         if(iret.lt.0)then
            call da_warning("da_read_obs_hdf5amsr2.inc",662, &
               (/"HDF5 read error for: CLW data"/))
         endif
         call H5Dclose_f(dhnd2,iret)
      ! change scale
         do j=1,nscan
            do i=1,lo_rez_fov
               if(clw(i,j).gt.-32767)clw(i,j)=clw(i,j)*sca
            enddo
         enddo
         do j=nscan+1,max_scan
            clw(:,j)=0
         enddo
      ! sample display
         !write(unit=stdout,fmt=*)'clw(pixel=1,scan=1): ',clw(1,1)
      ! close file and 1
         call H5Fclose_f(fhnd2,iret)
      end if

! 2.0 Loop to read hdf file and assign information to a sequential structure
!-------------------------------------------------------------------------

   ! Allocate arrays to hold data
      if ( .not. head_allocated ) then
         allocate (head)
         nullify  ( head % next )
         p => head
         head_allocated = .true.
      end if
   ! start scan_loop
      scan_loop:     do iscan=1, nscan
         do i = 1, 6
            idate5(i)=obstime(i, iscan)
         end do
         call da_get_julian_time(idate5(1),idate5(2),idate5(3),idate5(4),idate5(5),obs_time)
         if ( obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time) ) cycle scan_loop
         do ifgat=1,num_fgat_time
            if ( obs_time >= time_slots(ifgat-1) .and.  &
               obs_time  < time_slots(ifgat) ) exit
         end do

      ! start fov_loop
         fov_loop:      do ifov=1, lo_rez_fov
            num_amsr2_file       = num_amsr2_file + 1
            num_amsr2_file_local = num_amsr2_file_local + 1
            info%lat  =  latlr(ifov,iscan)
            info%lon  =  lonlr(ifov,iscan)

            call da_llxy (info, loc, outside, outside_all)
            if (outside_all) cycle fov_loop

            num_amsr2_global       = num_amsr2_global + 1
            num_amsr2_global_local = num_amsr2_global_local + 1
            ptotal(ifgat) = ptotal(ifgat) + 1
            if (outside) cycle fov_loop   ! No good for this PE
         ! Discard data over Land (landmask =0 -->Land =1 -->Sea)
            landsea_mask = 0
            if(lof06(ifov,iscan) < 1.0 .and. lof10(ifov,iscan) < 1.0 .and. &
               lof23(ifov,iscan) < 1.0 .and. lof36(ifov,iscan) < 1.0 .and. &
               lof89a(2*ifov-1,iscan)  < 1.0  ) landsea_mask = 1
            if( landsea_mask == 0 ) cycle fov_loop

            num_amsr2_local       = num_amsr2_local + 1
            num_amsr2_local_local = num_amsr2_local_local + 1
            write(unit=info%date_char, &
            fmt='(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
               idate5(1), '-', idate5(2), '-', idate5(3), '_', idate5(4), &
               ':', idate5(5), ':', idate5(6)
            info%elv = 0.0

! 3.0  Make Thinning
! Map obs to thinning grid
!-------------------------------------------------------------------
            if (thinning) then
               dlat_earth = info%lat !degree
               dlon_earth = info%lon
               if (dlon_earth<zero)  dlon_earth = dlon_earth+r360
               if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
               dlat_earth = dlat_earth*deg2rad !radian
               dlon_earth = dlon_earth*deg2rad
               crit = 1.
               call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
               if (.not. iuse) then
                  num_amsr2_thinned = num_amsr2_thinned+1
                  cycle fov_loop
               end if
            end if

            num_amsr2_used = num_amsr2_used + 1
            data_all = missing_r

            do k=1,num_low_freq_chan
               tb = tb_low(ifov,iscan,k)
               if( tb < tbmin .or. tb > tbmax ) tb = missing_r
               data_all(k)= tb
            enddo

            tb = tb89av(2*ifov-1,iscan)
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(13)= tb

            tb = tb89ah(2*ifov-1,iscan)
            if( tb < tbmin .or. tb > tbmax ) tb = missing_r
            data_all(14)= tb

! 4.0 assign information to sequential radiance structure
!--------------------------------------------------------------------------
            allocate ( p % tb_inv (1:nchan ))
            p%info             = info
            p%loc              = loc
            p%landsea_mask     = landsea_mask
            p%scanpos          = ifov
            p%satzen           = sat_zenith(ifov,iscan)
            p%satazi           = sat_azimuth(ifov,iscan)
            p%solzen           = sun_zen(ifov,iscan)
            p%solazi           = sun_az(ifov,iscan)
            p%clw              = clw(ifov,iscan)
            p%tb_inv(1:nchan)  = data_all(1:nchan)
            p%sensor_index     = inst
            p%ifgat            = ifgat

            allocate (p%next)   ! add next data
            p => p%next
            nullify (p%next)
         end do fov_loop
      end do scan_loop

   ! Dellocate arrays
      deallocate  (obstime)

      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_amsr2_file    : ',num_amsr2_file_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_amsr2_global  : ',num_amsr2_global_local
      write(stdout,fmt='(3a,i7)') ' In file: ',trim(fname_tb(ifile)),' got num_amsr2_local   : ',num_amsr2_local_local
   end do infile_loop

   call H5close_f(iret)

   deallocate(data_all) ! Deallocate data arrays

   if (thinning .and. num_amsr2_global > 0 ) then
   ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         j = j + thinning_grid(inst,ifgat)%itxmax
      end do 

      allocate ( in  (j) )
      allocate ( out (j) )
      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            in(j) = thinning_grid(inst,ifgat)%score_crit(i)
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time
         do i = 1, thinning_grid(inst,ifgat)%itxmax
            j = j + 1
            if ( ABS(out(j)-thinning_grid(inst,ifgat)%score_crit(i)) > 1.0E-10 ) &
            thinning_grid(inst,ifgat)%ibest_obs(i) = 0
         end do
      end do

      deallocate( in  )
      deallocate( out )


   ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_amsr2_used_tmp = num_amsr2_used
      do j = 1, num_amsr2_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            end if
         end do

      ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            end if
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_amsr2_thinned = num_amsr2_thinned + 1
            num_amsr2_used = num_amsr2_used - 1
            continue
         end if

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         end if

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         end if

      end do

   end if  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel + num_amsr2_used
   iv%total_rad_channel = iv%total_rad_channel + num_amsr2_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_amsr2_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_amsr2_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1)
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i)
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning("da_read_obs_hdf5amsr2.inc",901,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'AMSR2 data counts: '
   write(stdout,fmt='(a,i7)') ' In file: ',num_amsr2_file
   write(stdout,fmt='(a,i7)') ' Global : ',num_amsr2_global
   write(stdout,fmt='(a,i7)') ' Local  : ',num_amsr2_local
   write(stdout,fmt='(a,i7)') ' Used   : ',num_amsr2_used
   write(stdout,fmt='(a,i7)') ' Thinned: ',num_amsr2_thinned

!  5.0 allocate innovation radiance structure
!----------------------------------------------------------------

   if (num_amsr2_used > 0) then
      iv%instid(inst)%num_rad  = num_amsr2_used
      iv%instid(inst)%info%nlocal = num_amsr2_used
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
         'Allocating space for radiance innov structure', &
         inst, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad
      call da_allocate_rad_iv (inst, nchan, iv)
   end if

!  6.0 assign sequential structure to innovation structure
!-------------------------------------------------------------
   p => head

   do n = 1, num_amsr2_used
      i = p%sensor_index 
      call da_initialize_rad_iv (i, n, iv, p)
      current => p
      p => p%next
   ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do
   deallocate ( p )
   deallocate (ptotal)

   if (trace_use) call da_trace_exit("da_read_obs_hdf5amsr2")
end subroutine da_read_obs_hdf5amsr2
subroutine da_allocate_rad_iv (i, nchan, iv)

   !---------------------------------------------------------------------------
   !  Purpose: allocate radiance innovation structure
   !---------------------------------------------------------------------------

   use da_control

   implicit none

   integer           ,  intent (in)    :: i, nchan
   type (iv_type)    ,  intent (inout) :: iv

   call da_trace_entry("da_allocate_rad_iv")

      allocate (iv%instid(i)%info%date_char(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%name(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%platform(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%id(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%levels(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%lat(kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%lon(kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%elv(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%pstar(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%i  (kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%j  (kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%dx (kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%dy (kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%dxm(kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%dym(kts:kte,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%k  (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%dz (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%dzm(iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%zk (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%info%proc_domain(iv%instid(i)%nlevels,iv%instid(i)%num_rad))

      allocate (iv%instid(i)%t  (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%mr (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tm (kms:kme,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%qm (kms:kme,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%qrn(kms:kme,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%qcw(kms:kme,iv%instid(i)%num_rad))
      if ( crtm_cloud ) then
         allocate (iv%instid(i)%qci(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%qsn(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%qgr(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%qhl(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%rcw(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%rci(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%rrn(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%rsn(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%rgr(kms:kme,iv%instid(i)%num_rad))
         allocate (iv%instid(i)%rhl(kms:kme,iv%instid(i)%num_rad))
      end if
      allocate (iv%instid(i)%pm (kms:kme,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%pf (0:kme,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%u10(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%v10(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%t2m(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%q2m(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%mr2m(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%psfc(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%ts(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%smois(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tslb(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%snowh(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%isflg(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%landsea_mask(iv%instid(i)%num_rad))
      if (rtm_option == rtm_option_rttov) then
         allocate (iv%instid(i)%surftype(iv%instid(i)%num_rad))
         allocate (iv%instid(i)%snow_frac(iv%instid(i)%num_rad))
      end if
      allocate (iv%instid(i)%elevation(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%vegfra(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%vegtyp(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%soiltyp(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%clwp(iv%instid(i)%num_rad))
      if ( index(iv%instid(i)%rttovid_string, 'amsr2') > 0 ) then
         allocate (iv%instid(i)%clw(iv%instid(i)%num_rad))
      end if
      allocate (iv%instid(i)%ps(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_xb(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_qc(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_inv(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_error(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_sens(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%tb_imp(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%rad_xb(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%rad_obs(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%rad_ovc(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%emiss(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%scanpos(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%scanline(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%ifgat(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%cloud_flag(nchan,iv%instid(i)%num_rad))
      allocate (iv%instid(i)%rain_flag(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%satzen(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%satazi(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%solzen(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%solazi(iv%instid(i)%num_rad))
      allocate (iv%instid(i)%gamma_jacobian(nchan,iv%instid(i)%num_rad))
      if ( use_rttov_kmatrix .or. use_crtm_kmatrix ) then
         allocate(iv%instid(i)%ts_jacobian(nchan,iv%instid(i)%num_rad))
         allocate(iv%instid(i)%ps_jacobian(nchan,iv%instid(i)%num_rad))
         allocate(iv%instid(i)%emiss_jacobian(nchan,iv%instid(i)%num_rad))
         allocate(iv%instid(i)%windspeed_jacobian(nchan,iv%instid(i)%num_rad))
         allocate(iv%instid(i)%t_jacobian(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
         allocate(iv%instid(i)%q_jacobian(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      end if
      if (rtm_option == rtm_option_crtm) then
         allocate(iv%instid(i)%crtm_climat(iv%instid(i)%num_rad))
         allocate(iv%instid(i)%water_coverage(iv%instid(i)%num_rad))
         allocate(iv%instid(i)%land_coverage(iv%instid(i)%num_rad))
         allocate(iv%instid(i)%ice_coverage(iv%instid(i)%num_rad))
         allocate(iv%instid(i)%snow_coverage(iv%instid(i)%num_rad))
         if (use_crtm_kmatrix) then
            if ( crtm_cloud ) then
               allocate(iv%instid(i)%water_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%ice_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%rain_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%snow_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%graupel_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%hail_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%water_r_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%ice_r_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%rain_r_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%snow_r_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%graupel_r_jacobian(nchan,kte,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%hail_r_jacobian(nchan,kte,iv%instid(i)%num_rad))
            end if
            if ( calc_weightfunc ) then
               allocate(iv%instid(i)%lod(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%lod_jacobian(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%trans(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%trans_jacobian(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
               allocate(iv%instid(i)%der_trans(nchan,iv%instid(i)%nlevels,iv%instid(i)%num_rad))
            end if
         end if
      end if
	
      call da_trace_exit("da_allocate_rad_iv")

end subroutine da_allocate_rad_iv

subroutine da_initialize_rad_iv (i, n, iv, p)

   !---------------------------------------------------------------------------
   !  Purpose: allocate radiance innovation structure
   !---------------------------------------------------------------------------

   use da_control

   implicit none

   integer,             intent(in)    :: i, n
   type(datalink_type), intent(in)    :: p
   type(iv_type),       intent(inout) :: iv

   call da_trace_entry("da_initialize_rad_iv")

   iv%instid(i)%info%lat(:,n)   = p%info%lat
   iv%instid(i)%info%lon(:,n)   = p%info%lon
   iv%instid(i)%info%elv(n)     = p%info%elv
   iv%instid(i)%info%date_char(n) = p%info%date_char

   iv%instid(i)%info%max_lev    = iv%instid(i)%nlevels
   iv%instid(i)%info%levels(n)  = iv%instid(i)%nlevels
   iv%instid(i)%info%i  (:,n)   = p%loc%i
   iv%instid(i)%info%j  (:,n)   = p%loc%j
   iv%instid(i)%info%k  (:,n)   = 0
   iv%instid(i)%info%dx (:,n)   = p%loc%dx
   iv%instid(i)%info%dy (:,n)   = p%loc%dy
   iv%instid(i)%info%dz (:,n)   = 0.0
   iv%instid(i)%info%dxm(:,n)   = p%loc%dxm
   iv%instid(i)%info%dym(:,n)   = p%loc%dym
   iv%instid(i)%info%dzm(:,n)   = 0.0
   iv%instid(i)%info%proc_domain(:,n) = .false.
   ! z done in da_get_innov_vector_rad
   iv%instid(i)%t(:,n)          = 0.0
   iv%instid(i)%mr(:,n)         = 0.0
   iv%instid(i)%tm(:,n)         = 0.0
   iv%instid(i)%qm(:,n)         = 0.0
   iv%instid(i)%qrn(:,n)        = 0.0
   iv%instid(i)%qcw(:,n)        = 0.0
   if ( crtm_cloud ) then
      iv%instid(i)%qci(:,n)        = 0.0
      iv%instid(i)%qsn(:,n)        = 0.0
      iv%instid(i)%qgr(:,n)        = 0.0
      iv%instid(i)%qhl(:,n)        = 0.0
      iv%instid(i)%rcw(:,n)        = 0.0
      iv%instid(i)%rci(:,n)        = 0.0
      iv%instid(i)%rrn(:,n)        = 0.0
      iv%instid(i)%rsn(:,n)        = 0.0
      iv%instid(i)%rgr(:,n)        = 0.0
      iv%instid(i)%rhl(:,n)        = 0.0
   end if
   iv%instid(i)%pm(:,n)         = 0.0
   iv%instid(i)%pf(:,n)         = 0.0
   iv%instid(i)%u10(n)          = 0.0
   iv%instid(i)%v10(n)          = 0.0
   iv%instid(i)%t2m(n)          = 0.0
   iv%instid(i)%q2m(n)          = 0.0
   iv%instid(i)%mr2m(n)         = 0.0
   iv%instid(i)%psfc(n)         = 0.0
   iv%instid(i)%ts(n)           = 0.0
   iv%instid(i)%smois(n)        = 0.0
   iv%instid(i)%tslb(n)         = 0.0
   iv%instid(i)%snowh(n)        = 0.0
   iv%instid(i)%isflg(n)        = 0
   iv%instid(i)%soiltyp(n)      = 0.0
   iv%instid(i)%landsea_mask(n) = p%landsea_mask
   iv%instid(i)%elevation(n)    = 0.0
   iv%instid(i)%vegfra(n)       = 0.0
   iv%instid(i)%vegtyp(n)       = 0.0
   iv%instid(i)%clwp(n)         = 0.0
   if ( index(iv%instid(i)%rttovid_string, 'amsr2') > 0 ) then
      iv%instid(i)%clw(n)       = p%clw
   end if
   iv%instid(i)%ps(n)           = 0.0
   iv%instid(i)%tb_xb(:,n)      = 0.0
   iv%instid(i)%tb_inv(:,n)     = p%tb_inv(:)
   iv%instid(i)%tb_qc(:,n)      = 0
   iv%instid(i)%tb_error(:,n)   = 500.0
   iv%instid(i)%tb_sens(:,n)    = 0.0
   iv%instid(i)%tb_imp(:,n)     = 0.0
   iv%instid(i)%rad_xb(:,n)     = 0.0
   iv%instid(i)%rad_obs(:,n)    = 0.0
   iv%instid(i)%rad_ovc(:,:,n)  = 0.0
   iv%instid(i)%emiss(:,n)      = 0.0
   iv%instid(i)%scanpos(n)      = p%scanpos
   ! iv%instid(i)%scanline(n)    = p%scanline
   iv%instid(i)%scanline(n)     = 0
   iv%instid(i)%ifgat(n)        = p%ifgat
   iv%instid(i)%cloud_flag(:,n) = qc_good  ! no cloud
   iv%instid(i)%rain_flag(n)    = 0        ! no rain;  1:rain
   iv%instid(i)%satzen(n)       = p%satzen
   iv%instid(i)%satazi(n)       = p%satazi
   iv%instid(i)%solzen(n)       = p%solzen
   iv%instid(i)%solazi(n)       = p%solazi
 !  iv%instid(i)%solazi(n)       = 0.0

   if ( rtm_option == rtm_option_rttov ) then
      iv%instid(i)%surftype(n)     = 0
      iv%instid(i)%snow_frac(n)     = 0.0
   end if

   iv%instid(i)%gamma_jacobian(:,n)=0.0

   if ( use_rttov_kmatrix .or. use_crtm_kmatrix ) then
      iv%instid(i)%ts_jacobian(:,n)=0.0
      iv%instid(i)%ps_jacobian(:,n)=0.0
      iv%instid(i)%emiss_jacobian(:,n)=0.0
      iv%instid(i)%windspeed_jacobian(:,n)=0.0
      iv%instid(i)%t_jacobian(:,:,n)=0.0
      iv%instid(i)%q_jacobian(:,:,n)=0.0
   end if

   if (rtm_option == rtm_option_crtm) then
      iv%instid(i)%crtm_climat(n)=0  ! invalid_model
      iv%instid(i)%water_coverage(n)=1.0
      iv%instid(i)%land_coverage(n)=0.0
      iv%instid(i)%ice_coverage(n)=0.0
      iv%instid(i)%snow_coverage(n)=0.0
      if (use_crtm_kmatrix) then
         if ( crtm_cloud ) then
            iv%instid(i)%water_jacobian(:,:,n)=0.0
            iv%instid(i)%ice_jacobian(:,:,n)=0.0
            iv%instid(i)%rain_jacobian(:,:,n)=0.0
            iv%instid(i)%snow_jacobian(:,:,n)=0.0
            iv%instid(i)%graupel_jacobian(:,:,n)=0.0
            iv%instid(i)%hail_jacobian(:,:,n)=0.0
            iv%instid(i)%water_r_jacobian(:,:,n)=0.0
            iv%instid(i)%ice_r_jacobian(:,:,n)=0.0
            iv%instid(i)%rain_r_jacobian(:,:,n)=0.0
            iv%instid(i)%snow_r_jacobian(:,:,n)=0.0
            iv%instid(i)%graupel_r_jacobian(:,:,n)=0.0
            iv%instid(i)%hail_r_jacobian(:,:,n)=0.0
         end if
         if ( calc_weightfunc ) then
            iv%instid(i)%lod(:,:,n) = 0.0
            iv%instid(i)%lod_jacobian(:,:,n) = 0.0
            iv%instid(i)%trans(:,:,n) = 0.0
            iv%instid(i)%trans_jacobian(:,:,n) = 0.0
            iv%instid(i)%der_trans(:,:,n) = 0.0
         end if
      end if
   end if

   call da_trace_exit("da_initialize_rad_iv")

end subroutine da_initialize_rad_iv

subroutine da_read_kma1dvar (inst,iv, ob, infile)

   !---------------------------------------------------------------------------
   ! Purpose: read in kma 1dvar innovation output to innovation and obs structure
   !
   ! METHOD: use F90 sequantial data structure to avoid read file twice  
   !          so that da_scan_bufrtovs is not necessary any more.
   !          1. read radiance data in sequential data structure
   !          2. assign sequential data structure to innovation structure
   !                and deallocate sequential data structure
   !---------------------------------------------------------------------------

   implicit none

   ! subroutine argument
   integer           ,  intent (in)    :: inst
   character(20)     ,  intent (in)    :: infile
   type (y_type)     ,  intent (inout) :: ob
   type (iv_type)    ,  intent (inout) :: iv

   ! local variables
   integer          :: iost, n, i, j,inunit

   ! Declare local variables
   logical outside,inside_halo
   character(20)  :: instrument

   ! pixel information
   integer   ::  idate, npass, scanpos
   real      ::  rlat, rlon                         !  lat/lon in degrees   for Anfovs
   real      ::  satzen      !  scan angles for Anfovs
   integer   ::  landsea_mask
   real      ::  srf_height
   integer   ::  nchanl,knchan                      !  number of channels

   real      :: fastem(5)
   real , allocatable :: &         !  bright temperatures
                         otb(:),oerr(:),inov(:), emiss(:), &
                         tb(:), inv(:),err(:)
   integer, allocatable :: kchan(:),qc(:)

   type (datalink_type), pointer    :: head, p, current

   integer                        ::  size, error
   type(info_type)                ::  info
   type(model_loc_type)           ::  loc

   if (trace_use) call da_trace_entry("da_read_kma1dvar")

   !**************************************************************************
   ! Initialize variables

   nchanl = iv%instid(inst)%nchan

   allocate ( kchan(nchanl) )
   allocate ( otb(nchanl) )
   allocate ( oerr(nchanl) )
   allocate ( inov(nchanl) )
   allocate ( emiss(nchanl) )

   allocate ( tb(nchanl) )
   allocate ( inv(nchanl) )
   allocate ( err(nchanl) )
   allocate ( qc(nchanl) )

   ! 0.0  Open unit for kma 1dvar file and read file header
   !--------------------------------------------------------------

   call da_get_unit(inunit)
   open(unit=inunit,file=trim(infile),form='formatted', &
       iostat = iost, status = 'old')
   if (iost /= 0) then
      call da_error("da_read_kma1dvar.inc",73,&
         (/"Cannot open file "//infile/))
   end if

   read(unit=inunit,fmt='(A,I10,I12)') instrument,idate,npass

   ! Loop to read pixel and assign information to a sequential structure
   !-------------------------------------------------------------------------

   allocate ( head )
   nullify  ( head % next )
   p => head
   size = 0

   do j=1,npass
      read(inunit,'(2F12.4)')    rlat,rlon
      read(inunit,'(I5,20I5)')   knchan,(kchan(i),i=1,knchan)
      ! landsea_mask: 0:land ; 1:sea (same as RTTOV)
      read(inunit,'(I5,F12.4,I5,F12.4)')  scanpos,satzen,landsea_mask,srf_height
      read(inunit,'(20F12.4)') (oerr(i),i=1,knchan)
      read(inunit,'(20F12.4)') (emiss(i),i=1,knchan)
      read(inunit,'(5F12.4)')  (fastem(i),i=1,5)
      read(inunit,'(20F12.4)') (otb(i),i=1,knchan)
      read(inunit,'(20F12.4)') (inov(i),i=1,knchan)

      ! 1.0 Extract observation location and other required information
      !     judge if data is in the domain, read next record if not
      !-------------------------------------------------------------------

      info%lat  =  rlat
      info%lon  =  rlon
      call da_llxy (info, loc, outside, inside_halo )
      if (outside) cycle

      info%elv = srf_height
      write(unit=info%date_char, fmt='(i10)') idate 

      tb(1:nchanl) = missing_r
      inv(1:nchanl) = missing_r
      err(1:nchanl) = missing_r
      qc(1:nchanl) = qc_bad

      do i=1,knchan
         tb(kchan(i)) = otb(i)
         inv(kchan(i)) = inov(i)
         err(kchan(i)) = oerr(i)
         qc(kchan(i)) = qc_good
      end do 

      !  2.0   assign information to sequential radiance structure
      !--------------------------------------------------------------------

      allocate (p % tb_inv (1:nchanl))
      allocate (p % tb_error (1:nchanl))
      allocate (p % tb_qc (1:nchanl))
      allocate (p % tb_ob (1:nchanl))
      p%info             = info
      p%loc              = loc
      p%landsea_mask     = landsea_mask
      p%scanpos          = scanpos
      p%satzen           = satzen
      p%satazi           = missing_r
      p%solzen           = missing_r
      p%tb_inv(1:nchanl)     = inv(1:nchanl)
      p%tb_error(1:nchanl)   = err(1:nchanl)
      p%tb_qc(1:nchanl)      = qc(1:nchanl)
      p%tb_ob(1:nchanl)      = tb(1:nchanl)
      p%sensor_index         = inst

      size = size + 1
      allocate ( p%next, stat=error)   ! add next data
      if (error /= 0 ) then
         call da_error("da_read_kma1dvar.inc",145, &
            (/"Cannot allocate radiance structure"/))
      end if

      p => p%next
      nullify (p%next)
   end do
   
   iv%total_rad_pixel   = iv%total_rad_pixel + size
   iv%total_rad_channel = iv%total_rad_channel + size*nchanl

   deallocate (kchan)
   deallocate (otb)
   deallocate (oerr)
   deallocate (inov)
   deallocate (emiss)

   deallocate (tb)
   deallocate (inv)
   deallocate (err)
   deallocate (qc)

   close(inunit)
   call da_free_unit(inunit)

   !  3.0 allocate innovation and obs radiance structure
   !----------------------------------------------------------------
   iv%instid(inst)%num_rad = size
   ob%instid(inst)%num_rad = size
      
   write(unit=stdout,fmt='(a,i3,2x,a,3x,i10)')  ' allocating space for radiance innov structure', &
                             i, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad

   !  4.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   p => head
   allocate (iv%instid(inst)%tb_inv(1:nchanl,size))
   allocate (iv%instid(inst)%tb_error(1:nchanl,size))
   allocate (iv%instid(inst)%tb_qc(1:nchanl,size))
   allocate (ob%instid(inst)%tb(1:nchanl,size))
   do n = 1, size

      iv%instid(i)%info%name(n)	 = p%info%name
      iv%instid(i)%info%platform(n)  = p%info%platform
      iv%instid(i)%info%id(n) 	 = p%info%id
      iv%instid(i)%info%date_char(n) = p%info%date_char
      iv%instid(i)%info%levels(n)    = p%info%levels
      iv%instid(i)%info%lat(:,n)	 = p%info%lat
      iv%instid(i)%info%lon(:,n)	 = p%info%lon
      iv%instid(i)%info%elv(n)	 = p%info%elv
      iv%instid(i)%info%pstar(n)     = p%info%pstar

      iv%instid(inst)%info%i(:,n)    = p%loc%i
      iv%instid(inst)%info%j(:,n)    = p%loc%j
      iv%instid(inst)%info%dx(:,n)   = p%loc%dx
      iv%instid(inst)%info%dy(:,n)   = p%loc%dy
      iv%instid(inst)%info%dxm(:,n)  = p%loc%dxm
      iv%instid(inst)%info%dym(:,n)  = p%loc%dym

      iv%instid(inst)%landsea_mask(n) = p%landsea_mask
      iv%instid(inst)%scanpos(n)      = p%scanpos
      iv%instid(inst)%satzen(n) = p%satzen
      iv%instid(inst)%satazi(n) = p%satazi
      iv%instid(inst)%solzen(n) = p%solzen
      iv%instid(inst)%tb_inv(1:nchanl,n)   = p%tb_inv(1:nchanl)
      iv%instid(inst)%tb_error(1:nchanl,n) = p%tb_error(1:nchanl)
      iv%instid(inst)%tb_qc(1:nchanl,n)    = p%tb_qc(1:nchanl)
      ob%instid(inst)%tb(1:nchanl,n)       = p%tb_ob(1:nchanl)

      ! write(unit=stdout,*) inst, nread, iv%instid(inst)%tb_inv(1:nchanl,n)
      current => p
      p => p%next

      ! free current data
      deallocate (current % tb_ob)
      deallocate (current % tb_inv)
      deallocate (current % tb_error)
      deallocate (current % tb_qc)
      deallocate (current)
   end do

   ! check if sequential structure has been freed
   !
   ! p => head
   ! do i = 1, size
   !    write (unit=stdout,fmt=*)  i, p%tb(1:nchanl)%inv
   !    p => p%next
   ! end do

   if (trace_use) call da_trace_exit("da_readkma1dvar")

end subroutine da_read_kma1dvar


subroutine da_sort_rad (iv)

   !---------------------------------------------------------------------------
   ! Purpose: sorting radiance innovation to FGAT time bins
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type), intent (inout) :: iv

   integer                           :: i,j, n,t, error
   integer, allocatable              :: ifgat(:),landsea_mask(:)
   integer, allocatable              :: loc_i(:,:), loc_j(:,:),loc_k(:,:)
   real, allocatable                 :: loc_dx(:,:),loc_dy(:,:),loc_dz(:,:)
   real, allocatable                 :: loc_dxm(:,:),loc_dym(:,:),loc_dzm(:,:)
   character (len = 40), allocatable :: name(:)       
   character (len = 12), allocatable :: platform(:)   
   character (len =  5), allocatable :: id(:)          
   character (len = 19), allocatable :: date_char(:)   
   integer, allocatable              :: levels(:)      
   real, allocatable                 :: lat(:,:)         
   real, allocatable                 :: lon(:,:)         
   real, allocatable                 :: elv(:)        
   real, allocatable                 :: pstar(:)      
   real, allocatable                 :: scanpos(:), satzen(:), satazi(:), solzen(:)
   real, allocatable                 :: tb_inv(:,:)

   if (trace_use) call da_trace_entry("da_sort_rad")

   iv%info(radiance)%plocal(:) = 0
   if (num_fgat_time == 1) then
      do i=1,rtminit_nsensor
         iv%instid(i)%info%plocal(:) = 0
         iv%instid(i)%info%plocal(num_fgat_time) = iv%instid(i)%num_rad
         iv%info(radiance)%plocal(num_fgat_time) = iv%info(radiance)%plocal(num_fgat_time) + iv%instid(i)%num_rad
      end do
      if (trace_use) call da_trace_exit("da_sort_rad")
      return
   end if

   do i=1,rtminit_nsensor
      iv%instid(i)%info%plocal(:) = 0
      if (iv%instid(i)%num_rad < 1) cycle

      allocate (ifgat       (iv%instid(i)%num_rad)) 
      allocate (landsea_mask(iv%instid(i)%num_rad))
      allocate (loc_i       (kts:kte, iv%instid(i)%num_rad))
      allocate (loc_j       (kts:kte, iv%instid(i)%num_rad))
      allocate (loc_k       (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (loc_dx      (kts:kte, iv%instid(i)%num_rad))
      allocate (loc_dy      (kts:kte, iv%instid(i)%num_rad))
      allocate (loc_dz      (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (loc_dxm     (kts:kte, iv%instid(i)%num_rad))
      allocate (loc_dym     (kts:kte, iv%instid(i)%num_rad))
      allocate (loc_dzm     (iv%instid(i)%nlevels,iv%instid(i)%num_rad))
      allocate (name        (iv%instid(i)%num_rad))
      allocate (platform    (iv%instid(i)%num_rad))
      allocate (id          (iv%instid(i)%num_rad))
      allocate (date_char   (iv%instid(i)%num_rad))
      allocate (levels      (iv%instid(i)%num_rad))
      allocate (lat         (kts:kte,iv%instid(i)%num_rad))
      allocate (lon         (kts:kte,iv%instid(i)%num_rad))
      allocate (elv         (iv%instid(i)%num_rad))
      allocate (pstar       (iv%instid(i)%num_rad))
      allocate (scanpos     (iv%instid(i)%num_rad))
      allocate (satzen      (iv%instid(i)%num_rad))
      allocate (satazi      (iv%instid(i)%num_rad))
      allocate (solzen      (iv%instid(i)%num_rad))
      allocate (tb_inv      (iv%instid(i)%nchan,iv%instid(i)%num_rad))

      j = 0
      do t = 1,num_fgat_time
         do n = 1,iv%instid(i)%num_rad
            if (iv%instid(i)%ifgat(n) /= t) cycle
            j = j + 1
            ifgat(j)        = iv%instid(i)%ifgat(n)
            name(j)         = iv%instid(i)%info%name(n)
            platform(j)     = iv%instid(i)%info%platform(n)
            id(j)           = iv%instid(i)%info%id(n)
            date_char(j)    = iv%instid(i)%info%date_char(n)
            levels(j)       = iv%instid(i)%info%levels(n)
            elv(j)          = iv%instid(i)%info%elv(n)
            pstar(j)        = iv%instid(i)%info%pstar(n)

            lat    (:,j) = iv%instid(i)%info%lat(:,n)
            lon    (:,j) = iv%instid(i)%info%lon(:,n)
            loc_i  (:,j) = iv%instid(i)%info%i(:,n)
            loc_j  (:,j) = iv%instid(i)%info%j(:,n)
            loc_k  (:,j) = iv%instid(i)%info%k(:,n)
            loc_dx (:,j) = iv%instid(i)%info%dx(:,n)
            loc_dy (:,j) = iv%instid(i)%info%dy(:,n)
            loc_dxm(:,j) = iv%instid(i)%info%dxm(:,n)
            loc_dym(:,j) = iv%instid(i)%info%dym(:,n)
            loc_dz (:,j) = iv%instid(i)%info%dz(:,n)
            loc_dzm(:,j) = iv%instid(i)%info%dzm(:,n)

            landsea_mask(j) = iv%instid(i)%landsea_mask(n)
            scanpos(j)      = iv%instid(i)%scanpos(n)
            satzen(j)       = iv%instid(i)%satzen(n)
            satazi(j)       = iv%instid(i)%satazi(n)
            solzen(j)       = iv%instid(i)%solzen(n)

            tb_inv(1:iv%instid(i)%nchan,j) = iv%instid(i)%tb_inv(1:iv%instid(i)%nchan,n) 
         end do
         iv%instid(i)%info%plocal(t) = j
         !write (0,*) "da_sort_rad.inc",106,"i,t,iv%instid(i)%info%plocal(t)",i,t,iv%instid(i)%info%plocal(t)
      end do

      iv%info(radiance)%plocal = iv%info(radiance)%plocal + iv%instid(i)%info%plocal

      write(unit=stdout,fmt='(a,2x,a,2x,240i7)') &
         ' FGAT: ',iv%instid(i)%rttovid_string, iv%instid(i)%info%plocal(1:num_fgat_time)

      do n = 1,iv%instid(i)%num_rad
         iv%instid(i)%ifgat(n)        = ifgat(n)

         iv%instid(i)%info%name(n)	    = name(n)
         iv%instid(i)%info%platform(n)  = platform(n)
         iv%instid(i)%info%id(n) 	    = id(n)
         iv%instid(i)%info%date_char(n) = date_char(n)
         iv%instid(i)%info%levels(n)    = levels(n)
         iv%instid(i)%info%lat(:,n)     = lat(:,n)
         iv%instid(i)%info%lon(:,n)     = lon(:,n)
         iv%instid(i)%info%elv(n)	    = elv(n)
         iv%instid(i)%info%pstar(n)     = pstar(n)

         iv%instid(i)%info%i(:,n)    = loc_i(:,n)
         iv%instid(i)%info%j(:,n)    = loc_j(:,n)
         iv%instid(i)%info%k(:,n)    = loc_k(:,n)
         iv%instid(i)%info%dx(:,n)   = loc_dx(:,n)
         iv%instid(i)%info%dy(:,n)   = loc_dy(:,n)
         iv%instid(i)%info%dz(:,n)   = loc_dz(:,n)
         iv%instid(i)%info%dxm(:,n)  = loc_dxm(:,n)
         iv%instid(i)%info%dym(:,n)  = loc_dym(:,n)
         iv%instid(i)%info%dzm(:,n)  = loc_dzm(:,n)
         iv%instid(i)%landsea_mask(n) = landsea_mask(n)
         iv%instid(i)%scanpos(n)      = scanpos(n)
         iv%instid(i)%satzen(n)       = satzen(n)
         iv%instid(i)%satazi(n)       = satazi(n)
         iv%instid(i)%solzen(n)       = solzen(n)

         iv%instid(i)%tb_inv(1:iv%instid(i)%nchan,n) = tb_inv(1:iv%instid(i)%nchan,n)
      end do

      deallocate (ifgat) 
      deallocate (landsea_mask)
      deallocate (name)
      deallocate (platform)
      deallocate (id)
      deallocate (date_char)
      deallocate (levels)
      deallocate (lat)
      deallocate (lon)
      deallocate (elv)
      deallocate (pstar)
      deallocate (loc_i)
      deallocate (loc_j)
      deallocate (loc_k)
      deallocate (loc_dx)
      deallocate (loc_dy)
      deallocate (loc_dz)
      deallocate (loc_dxm)
      deallocate (loc_dym)
      deallocate (loc_dzm)
      deallocate (scanpos)
      deallocate (satzen)
      deallocate (satazi)
      deallocate (solzen)
      deallocate (tb_inv)
   end do

   if (trace_use) call da_trace_exit("da_sort_rad")

end subroutine da_sort_rad


subroutine da_setup_radiance_structures( grid, ob, iv )

   !---------------------------------------------------------------------------
   ! Purpose: Define, allocate and read of tovs raidance observation structure.
   !---------------------------------------------------------------------------

   implicit none

   type (domain) , intent(inout)   :: grid       ! model data
   type ( y_type), intent(inout)   :: ob         ! Observation structure.
   type (iv_type), intent(inout)   :: iv         ! O-B structure.

   character(len=200)          :: filename
   integer                     :: i, j, n, ios, ifgat
   logical                     :: lprinttovs 

   ! thinning variables
   integer  :: istart,iend,jstart,jend
   real     :: rlonlat(4)
   ! crtm_cloud
   integer  :: n1,n2,k,its,ite,jts,jte,kts,kte,inst
   
   if (trace_use) call da_trace_entry("da_setup_radiance_structures")

   !-------------------------------------------------------------------
   ! [1.0] Initialize RTTOV coefs and innovations vector for radiance
   !-------------------------------------------------------------------

    call da_radiance_init(iv, ob)
    
    do n = 1, rtminit_nsensor
       iv%instid(n)%rad_monitoring = rad_monitoring(n)
    enddo

   !-------------------------------
   ! 1.1 Make thinning grids
   !------------------------------
   call init_constants_derived

   if (thinning) then
      rlat_min =  r999
      rlat_max = -r999
      rlon_min =  r999
      rlon_max = -r999

      istart=MINVAL( grid%i_start(1:grid%num_tiles) )
      iend  =MAXVAL( grid%i_end  (1:grid%num_tiles) )
      jstart=MINVAL( grid%j_start(1:grid%num_tiles) )
      jend  =MAXVAL( grid%j_end  (1:grid%num_tiles) )

      do i = istart, iend
         do j = jstart, jend
            rlat_min=min(rlat_min, grid%xb%lat(i,j))
            rlat_max=max(rlat_max, grid%xb%lat(i,j))
            if( grid%xb%lon(i,j) < zero) then
              rlon_min=min(rlon_min, (r360+grid%xb%lon(i,j)))
              rlon_max=max(rlon_max, (r360+grid%xb%lon(i,j)))
            else
              rlon_min=min(rlon_min, grid%xb%lon(i,j))
              rlon_max=max(rlon_max, grid%xb%lon(i,j))
            endif
         enddo
      enddo

      call mpi_reduce(rlat_min, rlonlat(1), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlon_min, rlonlat(2), 1, true_mpi_real, mpi_min, root, comm, ierr)
      call mpi_reduce(rlat_max, rlonlat(3), 1, true_mpi_real, mpi_max, root, comm, ierr)
      call mpi_reduce(rlon_max, rlonlat(4), 1, true_mpi_real, mpi_max, root, comm, ierr)

      CALL mpi_bcast( rlonlat, 4 , true_mpi_real , root , comm, ierr )

      rlat_min = rlonlat(1)
      rlon_min = rlonlat(2)
      rlat_max = rlonlat(3)
      rlon_max = rlonlat(4)

      dlat_grid = rlat_max - rlat_min
      dlon_grid = rlon_max - rlon_min

      allocate(thinning_grid(iv%num_inst,num_fgat_time))
      do ifgat=1,num_fgat_time
         do n=1,iv%num_inst
            call makegrids (n,thinning_mesh(n),ifgat)
         end do
      end do 
   end if

   !-------------------------------------------------------------------
   ! [2.0] Read NCEP bufr tovs data in radiance innovations vector
   !-------------------------------------------------------------------

   if ( (.not. use_filtered_rad) .and. (.not. use_pseudo_rad) .and. (.not. use_simulated_rad) ) then

      if (use_hirs2obs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from hirs2.bufr'
         filename = 'hirs2 '
         call da_read_obs_bufrtovs ('hirs2', iv, filename)
      end if

      if (use_msuobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from msu.bufr'
         filename = 'msu'
         call da_read_obs_bufrtovs ('msu  ', iv, filename)
      end if

      if (use_hirs3obs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from hirs3.bufr'
         filename = 'hirs3'
         call da_read_obs_bufrtovs('hirs3', iv, filename)
      end if

      if (use_amsuaobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from amsua.bufr'
         filename = 'amsua'
         call da_read_obs_bufrtovs ('amsua', iv, filename)
      end if

      if (use_amsubobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from amsub.bufr'
         filename = 'amsub'
         call da_read_obs_bufrtovs ('amsub', iv, filename)
      end if

      if (use_hirs4obs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from hirs4.bufr'
         filename = 'hirs4'
         call da_read_obs_bufrtovs('hirs4', iv, filename)
      end if

      if (use_mhsobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from mhs.bufr'
         filename = 'mhs'
         call da_read_obs_bufrtovs('mhs  ', iv, filename)
      end if

      if (use_mwtsobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from mwtsa.dat and mwtsb.dat'
         filename = 'mwtsa'
         call da_read_obs_fy3('mwts ', iv, filename)
         filename = 'mwtsb'
         call da_read_obs_fy3('mwts ', iv, filename)
      end if

      if (use_mwhsobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from mwhsa.dat and mwhsb.dat'
         filename = 'mwhsa'
         call da_read_obs_fy3('mwhs ', iv, filename)
         filename = 'mwhsb'
         call da_read_obs_fy3('mwhs ', iv, filename)
      end if

      if (use_atmsobs) then
         write(unit=stdout,fmt='(a)') 'Reading radiance 1b data from atms.bufr'
         filename = 'atms'
         call da_read_obs_bufratms('atms ', iv, filename)
      end if

      if (use_airsobs) then
         write(unit=stdout,fmt='(a)') 'Reading airs 1b data from airs.bufr'
         filename = 'airs'
         call da_read_obs_bufrairs ('airs     ',iv, filename)
      end if

      if (use_eos_amsuaobs) then
         write(unit=stdout,fmt='(a)') 'Reading eos_amsua 1b data from airs.bufr'
         filename = 'airs'
         call da_read_obs_bufrairs ('eos_amsua',iv, filename)
      end if

      if (use_hsbobs) then
         write(unit=stdout,fmt='(a)') 'Reading hsb 1b data from airs.bufr'
         filename = 'airs'
         call da_read_obs_bufrairs ('hsb      ',iv, filename)
      end if

      if (use_ssmisobs) then
         write(unit=stdout,fmt='(a)') 'Reading ssmis data from ssmis.bufr'
         filename = 'ssmis'
         call da_read_obs_bufrssmis ('ssmis    ',iv, filename)
      end if
      if (use_iasiobs) then
         write(unit=stdout,fmt='(a)') 'Reading iasi data from iasi.bufr'
         filename = 'iasi'
         call da_read_obs_bufriasi ('iasi     ',iv, filename)
      end if
      if (use_seviriobs) then
         write(unit=stdout,fmt='(a)') 'Reading seviri data from seviri.bufr'
         filename = 'seviri'
         call da_read_obs_bufrseviri ('seviri   ',iv, filename)
      end if
      if (use_amsr2obs) then
         write(unit=stdout,fmt='(a)') 'Reading AMSR2 data in HDF5 format'
         call da_read_obs_hdf5amsr2 (iv, 'L1SGRTBR', 'L2SGCLWLD')
      end if
   end if

   if ( use_filtered_rad ) then
      write(unit=stdout,fmt='(a)') 'Reading filtered radiance'
      call da_read_filtered_rad (iv)
   end if

   if ( use_simulated_rad ) then
      write(unit=stdout,fmt='(a)') 'Reading simulated radiance'
 

      call da_read_simulated_rad (iv)
 

   end if

   if ( use_pseudo_rad ) then
      write(unit=stdout,fmt='(a)') 'Reading pseudo radiance from namelist'    
  
      call da_read_pseudo_rad (iv)
   

   end if

   if (use_kma1dvar) then
      do i=1,rtminit_nsensor
         filename = ' '
         filename='kma1dvar-'//trim(iv%instid(i)%rttovid_string)//'.inv'
         write(unit=stdout,fmt='(a,a)')  ' Reading KMA 1dvar innovation from  ', filename
         call da_read_kma1dvar (i,iv, ob, filename)
      end do
   end if

   if (thinning) then
      do ifgat=1,num_fgat_time
           do n=1,iv%num_inst
              call destroygrids (n,ifgat)
           end do
      end do
      deallocate(thinning_grid)
   end if

   ! sorting obs into FGAT time bins

   call da_sort_rad(iv)


   !-----------------------------------------------------------------------------
   ! [3.0] create (smaller) ob structure:
   !-----------------------------------------------------------------------------

   if (.not. use_kma1dvar) then
      do i = 1,  ob % num_inst
         ob % instid(i) % num_rad = iv % instid(i) % num_rad
         if (ob % instid(i) % num_rad < 1) cycle
         allocate (ob % instid(i) % tb(ob % instid(i) % nchan,ob % instid(i)%num_rad))
         ob % instid(i) % tb(:,:) = iv % instid(i) % tb_inv(:,:)
      end do
   end if

! Calculate DT for Cloudy Radiance DA

   if (use_rad .and. crtm_cloud .and. .not. DT_cloud_model) then
      its = grid%xp % its
      ite = grid%xp % ite
      jts = grid%xp % jts
      jte = grid%xp % jte
      kts = grid%xp % kts
      kte = grid%xp % kte

      grid%xb%delt(its:ite,jts:jte,kts:kte) = 0.0

      do inst= 1, iv % num_inst
         do n=1,iv%instid(inst)%num_rad
	     i = int(iv%instid(inst)%info%i(1,n))
	     j = int(iv%instid(inst)%info%j(1,n))
   	     grid%xb%delt(i:i+1, j:j+1, kts:kte) = 1800.0
         end do
      end do
   endif


   if (trace_use) call da_trace_exit("da_setup_radiance_structures")

end subroutine da_setup_radiance_structures

subroutine da_radiance_init(iv,ob)
!------------------------------------------------------------------------------
!  PURPOSE: subroutine to initialize radiances.
!
!  METHOD:  
!  1.0 Set up from namelist parameter
!  2.0 Set up some common variables for innovation/observation
!  3.0 Initialize RTTOV / 1
!  4.0 Set up bias correction
!  5.0 Read error factor file
!  6.0 Get FGAT time slots
!
!  HISTORY: 10/24/2007 Created from da_crtm_init            Tom Auligne
!  HISTORY: 12/15/2008 getting FGAT time slots is moved to 
!                      da_setup_obs_structures.inc.         Hui-Chuan Lin
!------------------------------------------------------------------------------

 implicit none 

 type (iv_type), intent (inout) :: iv
 type (y_type) , intent (inout) :: ob

!
!  local arguments
!------------------- 
 integer   :: n, j, ichan
 integer :: nsensor, unit_factor_rad
 integer     :: error
 integer, allocatable   ::  nscan(:), nchanl(:)

! local variables
!----------------
 integer             :: idum, wmo_sensor_id, sensor_type, iost
 integer             :: iunit
 character(len=filename_len)  :: filename

! local variables for tuning error factor
!----------------------------------------
 character(len=20)   ::  rttovid_string
 integer             ::  num_tot
 real                ::  joa, jo, trace, factor 

  call da_trace_entry("da_radiance_init")

!--------------------------------------------------------------
!  1.0 setup from namelist parameter
!--------------------------------------------------------------
  nsensor = rtminit_nsensor
  allocate (nscan(nsensor))
  allocate (nchanl(nsensor))

!----------------------------------------------------------------
!  2.0 set up some common variables for innovation/observation structure
!----------------------------------------------------------------
  iv % num_inst = nsensor
  ob % num_inst = nsensor

  allocate (iv%instid(1:nsensor))
  allocate (ob%instid(1:nsensor))
  allocate (satinfo(1:nsensor))

  iv%instid(1:nsensor)%num_rad = 0
  ob%instid(1:nsensor)%num_rad = 0

  loop_sensor: do n = 1, nsensor

   iv%instid(n)%platform_id  = rtminit_platform(n)
   iv%instid(n)%satellite_id = rtminit_satid(n)
   iv%instid(n)%sensor_id    = rtminit_sensor(n)
   if ( rtminit_satid(n) < 10 ) then
      write(iv%instid(n)%rttovid_string, '(a,i1,a)')  &
             trim( rttov_platform_name(rtminit_platform(n)) )//'-',  &
             rtminit_satid(n),     &
             '-'//trim( rttov_inst_name(rtminit_sensor(n)) )
      write(iv%instid(n)%rttovid_string_coef, '(a,i1,a)')  &
             trim( rttov_platform_name(rtminit_platform(n)) )//'_',  &
             rtminit_satid(n),     &
             '_'//trim( rttov_inst_name(rtminit_sensor(n)) )
   else
      write(iv%instid(n)%rttovid_string, '(a,i2.2,a)')  &
             trim( rttov_platform_name(rtminit_platform(n)) )//'-',  &
             rtminit_satid(n),     &
             '-'//trim( rttov_inst_name(rtminit_sensor(n)) )
      write(iv%instid(n)%rttovid_string_coef, '(a,i2.2,a)')  &
             trim( rttov_platform_name(rtminit_platform(n)) )//'_',  &
             rtminit_satid(n),     &
             '_'//trim( rttov_inst_name(rtminit_sensor(n)) )
   end if

   if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'msu' ) then
      nchanl(n)  = 4
      nscan(n)   = 11
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'hirs' ) then
      nchanl(n)  = 19
      nscan(n)   = 56
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'amsua' ) then
      nchanl(n)  = 15
      nscan(n)   = 30
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'amsub' ) then
      nchanl(n)  = 5
      nscan(n)   = 90
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'airs' ) then
      nchanl(n)  = 281
      nscan(n)   = 90
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'hsb' ) then
      nchanl(n)  = 4
      nscan(n)   = 90
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'mhs' ) then
      nchanl(n)  = 5
      nscan(n)   = 90
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'ssmis' ) then
      nchanl(n)  = 24
      nscan(n)   = 60
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'mwts' ) then
      nchanl(n)  = 4
      nscan(n)   = 15
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'mwhs' ) then
      nchanl(n)  = 5
      nscan(n)   = 98
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'atms' ) then
      nchanl(n)  = 22
      nscan(n)   = 96
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'iasi' ) then
     nchanl(n)  = 616
      nscan(n)   = 60	  
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'seviri' ) then
     nchanl(n)  =  8 
     nscan(n)   = 90 
   else if ( trim( crtm_sensor_name(rtminit_sensor(n))) == 'amsr2' ) then
      nchanl(n)  = 14
      nscan(n)   = 486
   else
      call da_error("da_radiance_init.inc",133, &
        (/"Unrecognized instrument"/))
   end if

   iv%instid(n)%nchan  = nchanl(n)
   ob%instid(n)%nchan  = nchanl(n)

   allocate ( iv%instid(n)%ichan(1:nchanl(n)), stat = error )
   if( error /= 0 ) then
      call da_error("da_radiance_init.inc",142, &
         (/"Memory allocation error to iv%instid(n)%ichan"/))
   end if

   allocate ( ob%instid(n)%ichan(1:nchanl(n)), stat = error )
   if( error /= 0 ) then
      call da_error("da_radiance_init.inc",148, &                                                           
         (/"Memory allocation error to ob%instid(n)%ichan"/))
   end if

   call da_get_unit(iunit)
   filename='radiance_info/'//trim(adjustl(iv%instid(n)%rttovid_string))//'.info'
   open(unit=iunit,file=filename, form='formatted',iostat = iost, status='old')

   if (iost /= 0) then
      message(1)="Cannot open radiance info file "//adjustl(filename)
      call da_error("da_radiance_init.inc",158,message(1:1))
   end if

   allocate ( satinfo(n) % ichan(nchanl(n)) )
   allocate ( satinfo(n) % iuse (nchanl(n)) )
   allocate ( satinfo(n) % error(nchanl(n)) )
   allocate ( satinfo(n) % polar(nchanl(n)) )

   read(iunit,*)
   do j = 1, nchanl(n)
     read(iunit,'(1x,5i5,2e18.10)')    &
                     wmo_sensor_id, &
               satinfo(n)%ichan(j), &
                       sensor_type, &
               satinfo(n)%iuse(j) , &
                              idum, &
               satinfo(n)%error(j), &
               satinfo(n)%polar(j)
     iv%instid(n)%ichan(j) = satinfo(n)%ichan(j)
     ob%instid(n)%ichan(j) = satinfo(n)%ichan(j)
   end do
   call da_free_unit(iunit)

   if ( use_blacklist_rad ) then
      call da_blacklist_rad(trim(rttov_platform_name(rtminit_platform(n))), &
                            rtminit_satid(n), &
                            trim(rttov_inst_name(rtminit_sensor(n))), &
                            nchanl(n), &
                            satinfo(n)%iuse )
   end if

  end do loop_sensor

!---------------------------------------------------------------------
! 3.0 Interface to the initialization subroutine of RTTOV and 1
!---------------------------------------------------------------------

    if (rtm_option == rtm_option_rttov) then
       call da_error("da_radiance_init.inc",199, &
          (/"Must compile with $RTTOV option for radiances"/))
    end if

    if (rtm_option == rtm_option_crtm) then
       call da_crtm_init(iv,ob,nsensor)
    end if

!-------------------------------------------------------
!  4.0 read bias correction coefs files
!-------------------------------------------------------

 loop_sensor2: do n = 1, nsensor

   allocate ( satinfo(n) % scanbias  (nchanl(n),nscan(n)) )
   allocate ( satinfo(n) % scanbias_b(nchanl(n),nscan(n),18) )
   allocate ( satinfo(n) % bcoef     (nchanl(n),4) )
   allocate ( satinfo(n) % bcoef0    (nchanl(n)) )
   allocate ( satinfo(n) % error_std (nchanl(n)) )

   satinfo(n) % error_std(:) = 500.0
   satinfo(n) % scanbias(:,:) = 0.0
   satinfo(n) % scanbias_b(:,:,:) = 0.0
   satinfo(n) % bcoef(:,:) = 0.0
   satinfo(n) % bcoef0(:) = 0.0

  if (read_biascoef) then
  !  new bias coefs files
  !  use o-b standard deviation statistics from Harris and Kelly method as obs errors
  !----------------------------------
 
    if ( index(iv%instid(n)%rttovid_string,'eos')  > 0 ) cycle   ! not implemented
    if ( index(iv%instid(n)%rttovid_string,'hirs') > 0 ) cycle   ! not implemented

    call da_read_biascoef(iv%instid(n)%rttovid_string, &
                      nchanl(n),nscan(n),18,4,global, &
                      satinfo(n)%scanbias, &
                      satinfo(n)%scanbias_b, &
                      satinfo(n)%bcoef, &
                      satinfo(n)%bcoef0, &
                      satinfo(n)%error_std)
  else
    ! use values specified in radiance_info files as obs errors
    satinfo(n)%error_std = satinfo(n)%error
  end if
 end do loop_sensor2

!-------------------------------------------------------
!  5.0 read error factor file
!-------------------------------------------------------
 if (use_error_factor_rad) then

    do n=1, rtminit_nsensor
       allocate ( satinfo(n)%error_factor(1:nchanl(n)) )
       satinfo(n)%error_factor(:) = 1.0
    end do

    call da_get_unit(unit_factor_rad)
    open(unit_factor_rad, file='radiance_error.factor', &
         form='formatted',iostat = iost, status='old')

    if (iost /= 0) then
       call da_error("da_radiance_init.inc",267, &
         (/"Cannot open radiance error factor file: radiance_error.factor"/))
    end if

    read(unit_factor_rad, *)
    do
      read(unit_factor_rad,fmt='(a15,i8,i8,3f15.5,f8.3)',iostat=iost)   &
          rttovid_string, ichan, num_tot, joa,jo,trace,factor
      if ( iost == 0 ) then
        do n=1, rtminit_nsensor
          if ( index(rttovid_string,trim(iv%instid(n)%rttovid_string))>0 ) then
             satinfo(n)%error_factor(ichan) = factor
             write(6,'(a,i5,a,f10.3)') trim(rttovid_string)//' Channel ', ichan, '  Error Factor = ', factor
             exit
          end if
        end do
      else
         exit
      end if
    end do
    close(unit_factor_rad)
    call da_free_unit(unit_factor_rad)

 end if

  deallocate(nscan)
  deallocate(nchanl)

  call da_trace_exit("da_radiance_init")
end subroutine da_radiance_init
subroutine da_get_innov_vector_radiance (it, grid, ob, iv)

   !---------------------------------------------------------------------------
   !  PURPOSE: Calculate innovation vector for radiance data.
   !
   !  METHOD:  d = y - H(x) - bc
   !       1. obs BT - simulated BT
   !       2. Bias correction
   !       3. Radiances Quality Control
   !
   !  HISTORY: 10/24/2007 - Creation from da_get_innov_vector_crtm  Tom Auligne
   !---------------------------------------------------------------------------
   
   implicit none
   
   integer, intent(in)            :: it       ! External iteration.
   type(domain),   intent(in)     :: grid
   type (y_type),  intent(inout)  :: ob       ! Observation structure.
   type (iv_type), intent(inout)  :: iv       ! O-B structure.

   integer                        :: inst

   if(trace_use) call  da_trace_entry("da_get_innov_vector_radiance")

   iv%instid(:)%info%n1 = iv%instid(:)%info%plocal(iv%time-1) + 1
   iv%instid(:)%info%n2 = iv%instid(:)%info%plocal(iv%time)

   !------------------------------------------------------------------------
   ! [1.0] calculate components of innovation vector
   !------------------------------------------------------------------------
   if (rtm_option == rtm_option_rttov) then
      call da_error("da_get_innov_vector_radiance.inc",36, &
       (/"Must compile with $RTTOV option for radiances"/))
   elseif (rtm_option == rtm_option_crtm) then
      call da_get_innov_vector_crtm (it, grid, ob, iv )
   else
      call da_warning("da_get_innov_vector_radiance.inc",47,(/"Unknown Radiative Transfer Model"/))
   endif

   !------------------------------------------------------------------------
   ! [2.0] Perform (Variational) bias correction
   !------------------------------------------------------------------------
   if (use_varbc .or. freeze_varbc) then
      call da_varbc_pred(iv)
      !varbc coldstart can not be done here when num_fgat_time>1
      if ( num_fgat_time == 1 ) then
         call da_varbc_coldstart(iv)  
      end if
      call da_varbc_direct(iv)      
   else if (biascorr) then
      do inst = 1, iv%num_inst                 ! loop for sensor
         write(unit=stdout,fmt='(A,A)') 'Performing bias correction for ', &
            trim(iv%instid(inst)%rttovid_string)
         call da_biascorr(inst,ob,iv)
      end do                                   ! end loop for sensor
   end if

   !------------------------------------------------------------------------
   ! [3.0] Perform QC check
   !------------------------------------------------------------------------
   if (qc_rad) then
      call da_qc_rad(it, ob, iv)
   end if

   !------------------------------------------------------------------------
   ! [4.0] Compute preconditioning for Variational bias correction
   !------------------------------------------------------------------------
   !varbc preconditioning shoud be done after get_innov_vector is done for all time slots
   !if (use_varbc .and. it == 1) call da_varbc_precond(iv)  !moved to da_get_innov_vector.inc
   
   !------------------------------------------------------------------------
   ! [5.0] Prepare (QCed) bias statistics files
   !------------------------------------------------------------------------
   if (biasprep) then
      do inst = 1, iv%num_inst
         write(unit=stdout,fmt='(A,A)') 'Preparing bias statistics files for ', &
            trim(iv%instid(inst)%rttovid_string)
         call da_biasprep(inst,ob,iv)
      end do
   end if

   if(trace_use) call  da_trace_exit("da_get_innov_vector_radiance")

end subroutine da_get_innov_vector_radiance
subroutine da_read_pseudo_rad (iv)

   !---------------------------------------------------------------------------
   !  Purpose: read in NCEP bufr tovs 1b data to innovation structure
   !
   !   METHOD: use F90 sequential data structure to avoid reading file twice  
   !            so that da_scan_bufrtovs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !               and deallocate sequential data structure
   !---------------------------------------------------------------------------




   use da_control

   implicit none

   type (iv_type),  intent (inout) :: iv

   ! Instrument triplet, follow the convension of RTTOV 
   integer   :: platform_id, satellite_id, sensor_id

   real, allocatable    :: tb_inv(:)  !  bright temperatures
   type (datalink_type) :: p

   logical :: outside, outside_all
   integer :: i,k,n,nchan,inst, alloc_stat
   type(info_type)       ::  info
   type(model_loc_type)  ::  loc

   call da_trace_entry("da_read_pseudo_rad")

   ! Initialize variables

   platform_id  = pseudo_rad_platid
   satellite_id = pseudo_rad_satid
   sensor_id    = pseudo_rad_senid
   if (sensor_id == 0) then
      nchan=19 !nchan_hirs
   else if (sensor_id == 1) then
      nchan=nchan_msu
   else if (sensor_id == 3) then
      nchan=nchan_amsua
   else if (sensor_id == 4)  then
      nchan=nchan_amsub
   else if (sensor_id == 15)  then
      nchan=nchan_mhs
   else if (sensor_id == 10)  then
      nchan=nchan_ssmis
   else if (sensor_id == 11)  then
      nchan=nchan_airs
   else if (sensor_id == 16)  then
   !iasi
      nchan= 616
   else if (sensor_id == 21)  then
   !seviri
      nchan = 8 
   else if (sensor_id == 63)  then
      !amsr2
      nchan = 14
   else
      call da_error("da_read_pseudo_rad.inc",65, &
           (/"Error in setting up pseudo radiance, unknown sensor_id."/))
   end if

   inst = 1    ! single instrument

   allocate (tb_inv(nchan))

   info%lat  =  pseudo_rad_lat ! in degree
   info%lon  =  pseudo_rad_lon
   info%date_char = "0000-00-00_00:00:00"
   call da_llxy (info, loc, outside, outside_all)

   if (outside) then
      iv%info(radiance)%nlocal          = 0
      iv%info(radiance)%ntotal          = 0
      iv%info(radiance)%ptotal(1)       = 0
      iv%instid(inst)%num_rad           = 0
      iv%instid(inst)%info%nlocal       = 0
   else  
      iv%info(radiance)%nlocal          = 1
      iv%info(radiance)%ntotal          = 1
      iv%info(radiance)%ptotal(1)       = 1
      iv%instid(inst)%num_rad           = 1
      iv%instid(inst)%info%nlocal       = 1
   end if

   do k = 1, nchan
      tb_inv(k) = missing_r
   end do
   tb_inv(pseudo_rad_ichan)=pseudo_rad_inv

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   
   i = 1

   if ( iv%instid(inst)%num_rad > 0 ) then
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(inst)%rttovid_string, iv%instid(inst)%num_rad

      call da_allocate_rad_iv(i,nchan,iv)

      !  6.0 assign sequential structure to innovation structure
      !-------------------------------------------------------------

      n = 1     ! single obs

      allocate (p%tb_inv(1:nchan), stat=alloc_stat)
      if ( alloc_stat /= 0 ) CALL da_error("da_read_pseudo_rad.inc",115,(/"error allocating"/))

      p%info             = info
      p%loc              = loc
      p%landsea_mask     = 1
      p%scanpos          = 1
      p%satzen           = 0.0
      p%satazi           = 0.0
      p%solzen           = 0.0
      p%tb_inv(1:nchan)  = tb_inv(1:nchan)
      p%sensor_index     = inst          
      p%ifgat            = 1
     
      call da_initialize_rad_iv (inst, n, iv, p)

      iv%instid(inst)%tb_qc(:,n)                   = qc_bad
      iv%instid(inst)%tb_error(pseudo_rad_ichan,n) = pseudo_rad_err
      iv%instid(inst)%tb_qc(pseudo_rad_ichan,n)    = qc_good

      deallocate(p%tb_inv)
   end if

   deallocate(tb_inv)

   call da_trace_exit("da_read_pseudo_rad")
  

end subroutine da_read_pseudo_rad

subroutine da_blacklist_rad (platform_name, satid, sensor_name, nchan, iuse)

   !---------------------------------------------------------------------------
   !  PURPOSE: black list for radiance data.
   !  METHOD:  based on black_ds.f90 provided by Paul Poli <Paul.Poli@ecmwf.int>
   !           -1 == CONSTANT
   !           -2 == EXPERIMENTAL
   !  HISTORY: 04/17/2012
   !---------------------------------------------------------------------------

   implicit none

   character(len=*),  intent(in)    :: platform_name, sensor_name
   integer,           intent(in)    :: satid, nchan
   integer,           intent(inout) :: iuse(nchan)

   character(len=8) :: cdate   !ccyymmdd
   character(len=6) :: ctime   !hhmmss

   cdate = analysis_date(1:4)//analysis_date(6:7)//analysis_date(9:10)
   ctime = analysis_date(12:13)//analysis_date(15:16)//analysis_date(18:19)

!============
! from ERA
!============
   select case ( trim(sensor_name) )
   case ( 'amsua' )
      iuse(1)  = -1
      iuse(2)  = -1
      iuse(3)  = -1
      iuse(4)  = -1
      iuse(14) = -1   ! NCEP
      iuse(15) = -1
   case ( 'amsub', 'mhs' )
      iuse(1)  = -1
      iuse(2)  = -1
   case ( 'hirs' )
      iuse(1)     = -1   ! NCEP
      iuse(16:19) = -1
   end select

   if ( trim(platform_name) == 'noaa' ) then
      select case ( satid )
      case ( 15 )
         if ( cdate < '19980801' )           iuse = -1   ! passive period
         if ( trim(sensor_name) == 'amsub' ) iuse = -1   ! bad quality noaa-15-amsub
         if ( trim(sensor_name) == 'hirs' ) then
            ! noisy for several periods, unknown reason
            if ( cdate == '19981006' .and. ctime == '060000' )   iuse = -1
            if ( cdate == '19990302' .and. ctime == '060000' )   iuse = -1
            if ( cdate >= '20000219' .and. cdate <= '20000220' ) iuse = -1
            ! unreliable since 7 June 2000 (filter wheel, ...)
            if ( cdate >= '20000607' ) iuse = -1
         end if   ! end if noaa-15-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            ! NCEP 199811
            ! noisy, unknown reason
            if ( cdate == '19990101' .and. ctime == '060000' ) iuse = -1
            ! channel 14 failure
            if ( cdate >= '20001030' ) iuse(14) = -1
            if ( cdate >= '20020410' ) iuse(11) = -1
         end if   ! end if noaa-15-amsua
      case ( 16 )
         ! data unusable before 26 Oct 2000 due to calibration errors
         if ( cdate < '20001026' ) iuse = -1
         if ( trim(sensor_name) == 'hirs' ) then
            if ( cdate == '20001208' .and. ctime == '000000' ) iuse = -1
            if ( cdate >  '20040524' ) iuse = -1
         end if   ! end if noaa-16-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            if ( cdate == '20010111' ) iuse = -1
            ! channel 8 noisy from the beginning
            iuse(8) = -1
            if ( cdate > '20030512' .and. cdate <  '20040701' ) iuse(9) = -1
            if ( cdate > '20050117' .and. cdate <  '20050125' ) iuse(9:14) = -1
            if ( cdate > '20021121' .and. cdate <= '20021203' ) iuse = -1
            if ( cdate > '20021203' .and. cdate <= '20021216' ) iuse(5:8) = -1
            !to be checked if ( cdate > '20070124' )
            ! most channel noisy due to solar eclipse
            if ( cdate >= '20090123' .and. cdate < '20090202' ) iuse = -1
            ! unexplained instrument anomaly
            if ( cdate >= '20090408' .and. cdate < '20090416' ) iuse = -1
            if ( cdate >  '20090603' ) iuse = -1
         end if   ! end if noaa-16-amsua
         if ( trim(sensor_name) == 'amsub' ) then
            if ( cdate >= '20070918' ) iuse = -1
         end if   ! end if noaa-16-amsub
      case ( 17 )
         if ( trim(sensor_name) == 'hirs' ) then
            ! channels 2, 3, 4, 5 became suddenly much more noisy
            if ( cdate >= '20101023' ) iuse(2:5) = -1
            ! increase noise, some channels are really bad, compromising cloud detection
            if ( cdate >= '20111204' ) iuse = -1
         end if   ! end if noaa-17-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            if ( cdate >  '20030129' ) iuse(7) = -1
            ! unusable due to scan motor problems
            if ( cdate >= '20031028' ) iuse = -1
         end if   ! end if noaa-17-amsua
         if ( trim(sensor_name) == 'amsub' ) then
            ! blackbody/space view anomalies
            if ( cdate > '20091216' ) iuse = -1
         end if   ! end if noaa-17-amsub
      case ( 18 )
         if ( trim(sensor_name) == 'hirs' ) then
            ! noaa-18 hirs data were assimilated until 20090312 in EI, but really
            ! this blacklisting should be done from the beginning
            iuse = -1
         end if   ! end if noaa-18-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            ! channel 9 briefly noisy
            if ( cdate > '20071115' .and. cdate < '20071219' ) iuse(9) = -1
            if ( cdate > '20101215' .and. cdate < '20101226' ) iuse(9) = -1
            ! orbit too close to noaa-19, so data heavily thnned out.
            ! poor sampling can lead to problems with VarBC
            !to be checked if ( cdate >= '20090804' .and. cdate < '20091224' ) iuse = -1
         end if   ! end if noaa-18-amsua
      case ( 19 )
         ! first noaa-19 data to arrive in opearions: 2009040700
         if ( cdate < '20090407' ) iuse = -1
         if ( cdate >= '20090407' .and. cdate < '20090501' ) iuse = -1
         if ( cdate >= '20090501' .and. cdate < '20090602' ) iuse = -1
         if ( trim(sensor_name) == 'amsua' ) then
            ! channel 8 excluded due to temporary noise issues
            if ( cdate >= '20090627' .and. cdate < '20090804' ) iuse(8) = -1
            if ( cdate >= '20091224' ) iuse(8) = -1
         end if   ! end if noaa-19-amsua
         if ( trim(sensor_name) == 'mhs' ) then
            ! allow bias adjustment after update of antenna pattern correction
            if ( cdate > '20090621' .and. cdate < '20090627' ) iuse = -1
            ! noise has slowly increased to above specifications
            if ( cdate >= '20090804' .and. cdate < '20100317' ) iuse = -1
            ! noise stayed high for channel 3
            if ( cdate >= '20100317' ) iuse(3) = -1
         end if   ! end if noaa-19-mhs
      end select
   end if   ! end if noaa

   if ( trim(platform_name) == 'eos' ) then
      select case ( satid )
      case ( 2 )
         if ( trim(sensor_name) == 'amsua' ) then
            iuse(7) = -1
            !to be checked if ( cdate > '20071001' )
            ! channel 5 noise has slowly increased too much
            if ( cdate >= '20100502' ) iuse(5) = -1
         end if   ! end if eos-2-amsua
         if ( trim(sensor_name) == 'airs' ) then
            ! airs data first appear on 2002102200
            if ( cdate < '20030401' ) iuse = -1    ! NCEP 200211
            ! data missing for a few days
            if ( cdate >= '20100101' .and. cdate <= '20100104' ) iuse = -1
            ! instrument anomaly: stand-by mode on 20100109, instrument restart on 20100121,
            ! data flow resumed on 20100128
            if ( cdate >= '20100109' .and. cdate < '20100128' ) iuse = -1
            if ( cdate >  '20070324' .and. cdate < '20070403' ) iuse = -1
            ! preventive blacklist of airs from 31 Oct 2003 onwards (solar storm)
            if ( cdate > '20031030' .and. cdate < '20031202' ) iuse = -1
         end if   ! end if eos-2-airs
      end select
   end if   ! end if eos

   if ( trim(platform_name) == 'metop' ) then
      select case ( satid )
      case ( 2 )
         if ( trim(sensor_name) == 'hirs' ) then
            ! first metop-2 hirs data arrived 20061130 12UTC
            if ( cdate <= '20061129' ) then
               iuse = -1
            else
               ! a few noisy days in Jan 2007
               if ( cdate < '20070201' ) iuse = -1
            end if
         end if   ! end if metop-2-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            ! first metop-2 data arrived 20061102 00UTC for 3 days, then
            ! data gap util 20061129 12UTC
            ! NCEP 200706
            if ( cdate <= '20061128' ) then
               iuse = -1
            else
               !to be checked
               ! calculate varbc coefficients at the beginning and during
               ! a recalibration period
               if ( cdate < '20070111' .or. &
                    (cdate > '20070521' .and. cdate < '20070528') ) then
                  iuse = -1
               end if
            end if
            ! channel 7 gone too noisy
            if ( cdate > '20090105' ) iuse(7) = -1
         end if   ! end if metop-2-amsua
         if ( trim(sensor_name) == 'mhs' ) then
            ! first metop-2 data arrived 20061130 12UTC
            if ( cdate <= '20061129' ) then
               iuse = -1
            else
               !to be checked, calculate varbc coefficients
               if ( cdate < '20070111' ) then
                  iuse = -1
               end if
            end if
         end if   ! end if metop-2-mhs
      end select
   end if   ! end if metop

!============
! from NCEP
!============
   if ( trim(platform_name) == 'noaa' ) then
      select case ( satid )
      case ( 15 )
         if ( trim(sensor_name) == 'hirs' ) then
            if ( cdate <  '19981101' ) iuse = -1
            if ( cdate >= '20000701' ) iuse = -1
         end if   ! end if noaa-15-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            if ( cdate <  '19981101' ) iuse = -1
            if ( cdate >= '20001101' ) iuse(11) = -1
         end if   ! end if noaa-15-amsua
         if ( trim(sensor_name) == 'amsub' ) then
            if ( cdate <  '20000201' ) iuse = -1
            if ( cdate >= '20060301' ) iuse(4) = -1
         end if   ! end if noaa-15-amsub
      case ( 16 )
         if ( trim(sensor_name) == 'hirs' ) then
            if ( cdate <  '20010201' ) iuse = -1
            if ( cdate >= '20040601' ) iuse = -1
         end if   ! end if noaa-16-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            if ( cdate <  '20010201' ) iuse = -1
            if ( cdate >= '20050201' ) iuse(9:13) = -1
         end if   ! end if noaa-16-amsua
         if ( trim(sensor_name) == 'amsub' ) then
            if ( cdate <  '20010201' ) iuse = -1
         end if   ! end if noaa-16-amsub
      case ( 17 )
         if ( trim(sensor_name) == 'hirs' ) then
            if ( cdate <  '20021001' ) iuse = -1
         end if   ! end if noaa-17-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            iuse = -1
         end if   ! end if noaa-17-amsua
         if ( trim(sensor_name) == 'amsub' ) then
            if ( cdate <  '20021001' ) iuse = -1
         end if   ! end if noaa-17-amsub
      case ( 18 )
         if ( trim(sensor_name) == 'hirs' ) then
            iuse = -1
         end if   ! end if noaa-18-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            if ( cdate <  '20051101' ) iuse = -1
            if ( cdate >= '20071201' ) iuse(9) = -1
         end if   ! end if noaa-18-amsua
         if ( trim(sensor_name) == 'mhs' ) then
            if ( cdate < '20051101' ) iuse = -1
         end if   ! end if noaa-18-mhs
      case ( 19 )
         if ( trim(sensor_name) == 'amsua' ) then
         end if   ! end if noaa-19-amsua
         if ( trim(sensor_name) == 'mhs' ) then
         end if   ! end if noaa-19-mhs
      end select
   end if   ! end if noaa

   if ( trim(platform_name) == 'eos' ) then
      select case ( satid )
      case ( 2 )
         if ( trim(sensor_name) == 'amsua' ) then
            iuse(7) = -1
            if ( cdate < '20021101' ) iuse = -1
         end if   ! end if eos-2-amsua
         if ( trim(sensor_name) == 'airs' ) then
            if ( cdate < '20021101' ) iuse = -1
         end if   ! end if eos-2-airs
      end select
   end if   ! end if eos

   if ( trim(platform_name) == 'metop' ) then
      select case ( satid )
      case ( 2 )
         if ( trim(sensor_name) == 'hirs' ) then
            if ( cdate < '20070601' ) iuse = -1
         end if   ! end if metop-2-hirs
         if ( trim(sensor_name) == 'amsua' ) then
            if ( cdate < '20070601' ) iuse = -1
         end if   ! end if metop-2-amsua
         if ( trim(sensor_name) == 'mhs' ) then
            if ( cdate < '20070601' ) iuse = -1
         end if   ! end if metop-2-mhs
      end select
   end if   ! end if metop

end subroutine da_blacklist_rad
  subroutine da_deallocate_radiance ( ob, iv, j)

   !-----------------------------------------------------------------------
   ! Purpose: deallocate radiance related structures/arrays
   ! Extracted from da_solve.inc
   !-----------------------------------------------------------------------

   implicit none

   type (y_type),        intent(inout)  :: ob        ! Observation structure.
   type (iv_type),       intent(inout)  :: iv        ! Obs. increment structure.
   type (j_type),        intent(inout)  :: j         ! Cost function.

   integer                        :: i,n,ichan

   if (trace_use) call da_trace_entry("da_deallocate_radiance")

      do i =1, iv%num_inst
         deallocate (j % jo % rad(i) % jo_ichan)
         deallocate (j % jo % rad(i) % num_ichan)
         deallocate (satinfo(i) % ichan)
         deallocate (satinfo(i) % iuse)
         deallocate (satinfo(i) % error)
         deallocate (satinfo(i) % polar)

         deallocate (satinfo(i) % scanbias) 
         deallocate (satinfo(i) % scanbias_b)
         deallocate (satinfo(i) % bcoef)
         deallocate (satinfo(i) % bcoef0)
         deallocate (satinfo(i) % error_std)

         deallocate (ob%instid(i) % ichan)
         deallocate (iv%instid(i) % ichan)

         if (iv%instid(i)%num_rad > 0) then
 
            deallocate (iv%instid(i)%info%date_char)
            deallocate (iv%instid(i)%info%name)
            deallocate (iv%instid(i)%info%platform)
            deallocate (iv%instid(i)%info%id)
            deallocate (iv%instid(i)%info%levels)     
            deallocate (iv%instid(i)%info%lat)      
            deallocate (iv%instid(i)%info%lon)      
            deallocate (iv%instid(i)%info%elv)   

            deallocate (iv%instid(i)%info%pstar)
            deallocate (iv%instid(i)%info%i)
            deallocate (iv%instid(i)%info%j)
            deallocate (iv%instid(i)%info%k)
            deallocate (iv%instid(i)%info%zk)
            deallocate (iv%instid(i)%info%dx)
            deallocate (iv%instid(i)%info%dy)
            deallocate (iv%instid(i)%info%dz)
            deallocate (iv%instid(i)%info%dxm)
            deallocate (iv%instid(i)%info%dym)
            deallocate (iv%instid(i)%info%dzm)
            deallocate (iv%instid(i)%info%proc_domain)

            deallocate (iv%instid(i)%t)
            deallocate (iv%instid(i)%mr)
            deallocate (iv%instid(i)%tm)
            deallocate (iv%instid(i)%qm)
            deallocate (iv%instid(i)%qrn)
            deallocate (iv%instid(i)%qcw)
            if ( crtm_cloud ) then
               deallocate (iv%instid(i)%qci)
               deallocate (iv%instid(i)%qsn)
               deallocate (iv%instid(i)%qgr)
               deallocate (iv%instid(i)%qhl)
               deallocate (iv%instid(i)%rcw)
               deallocate (iv%instid(i)%rci)
               deallocate (iv%instid(i)%rrn)
               deallocate (iv%instid(i)%rsn)
               deallocate (iv%instid(i)%rgr)
               deallocate (iv%instid(i)%rhl)
            end if
            deallocate (iv%instid(i)%pm)
            deallocate (iv%instid(i)%pf)
            deallocate (iv%instid(i)%u10)
            deallocate (iv%instid(i)%v10)
            deallocate (iv%instid(i)%t2m)
            deallocate (iv%instid(i)%q2m)
            deallocate (iv%instid(i)%mr2m)
            deallocate (iv%instid(i)%psfc)
            deallocate (iv%instid(i)%ts)
            deallocate (iv%instid(i)%smois)
            deallocate (iv%instid(i)%tslb)
            deallocate (iv%instid(i)%snowh)
            deallocate (iv%instid(i)%isflg)
            deallocate (iv%instid(i)%soiltyp)
            deallocate (iv%instid(i)%landsea_mask)
            if (rtm_option == rtm_option_rttov) then
               deallocate (iv%instid(i)%surftype)
               deallocate (iv%instid(i)%snow_frac)
            end if
            deallocate (iv%instid(i)%elevation)
            deallocate (iv%instid(i)%vegfra)
            deallocate (iv%instid(i)%vegtyp)
            deallocate (iv%instid(i)%clwp)
            if ( index(iv%instid(i)%rttovid_string,'amsr2') > 0 ) then
               deallocate (iv%instid(i)%clw)
            end if
            deallocate (iv%instid(i)%ps)
            deallocate (iv%instid(i)%tb_xb)
            deallocate (iv%instid(i)%tb_qc)
            deallocate (iv%instid(i)%tb_inv)
            deallocate (iv%instid(i)%tb_error)
            deallocate (iv%instid(i)%tb_sens)
            deallocate (iv%instid(i)%tb_imp)
            deallocate (iv%instid(i)%rad_xb)
            deallocate (iv%instid(i)%rad_obs)
            deallocate (iv%instid(i)%rad_ovc)
            deallocate (iv%instid(i)%emiss)
            deallocate (iv%instid(i)%scanpos)
            deallocate (iv%instid(i)%scanline)
            deallocate (iv%instid(i)%ifgat)
            deallocate (iv%instid(i)%cloud_flag)
            deallocate (iv%instid(i)%rain_flag)
            deallocate (iv%instid(i)%satzen)
            deallocate (iv%instid(i)%satazi)
            deallocate (iv%instid(i)%solzen)
            deallocate (iv%instid(i)%solazi)
            deallocate (iv%instid(i)%gamma_jacobian)

           if (ANY(use_satcv)) then
	      if (use_satcv(2)) then
                 do n = 1,iv%instid(i)%num_rad
                    deallocate (iv%instid(i)%cv_index(n)%cc)
                    deallocate (iv%instid(i)%cv_index(n)%vtox)
                 end do
	      end if
              deallocate (iv%instid(i)%cv_index)  
           end if
        
           if ( use_rttov_kmatrix .or. use_crtm_kmatrix ) then
              deallocate(iv%instid(i)%ts_jacobian)
              deallocate(iv%instid(i)%ps_jacobian)
              deallocate(iv%instid(i)%emiss_jacobian)
              deallocate(iv%instid(i)%windspeed_jacobian)
              deallocate(iv%instid(i)%t_jacobian)
              deallocate(iv%instid(i)%q_jacobian)
           end if
            if (rtm_option == rtm_option_crtm) then
               deallocate(iv%instid(i)%crtm_climat)
               deallocate(iv%instid(i)%water_coverage)
               deallocate(iv%instid(i)%land_coverage)
               deallocate(iv%instid(i)%ice_coverage)
               deallocate(iv%instid(i)%snow_coverage)
               if (use_crtm_kmatrix) then
                 if ( crtm_cloud ) then
                     deallocate(iv%instid(i)%water_jacobian)
                     deallocate(iv%instid(i)%ice_jacobian)
                     deallocate(iv%instid(i)%rain_jacobian)
                     deallocate(iv%instid(i)%snow_jacobian)
                     deallocate(iv%instid(i)%graupel_jacobian)
                     deallocate(iv%instid(i)%hail_jacobian)
                     deallocate(iv%instid(i)%water_r_jacobian)
                     deallocate(iv%instid(i)%ice_r_jacobian)
                     deallocate(iv%instid(i)%rain_r_jacobian)
                     deallocate(iv%instid(i)%snow_r_jacobian)
                     deallocate(iv%instid(i)%graupel_r_jacobian)
                     deallocate(iv%instid(i)%hail_r_jacobian)
                 end if
                 if ( calc_weightfunc ) then
                    deallocate(iv%instid(i)%lod)
                    deallocate(iv%instid(i)%lod_jacobian)
                    deallocate(iv%instid(i)%trans)
                    deallocate(iv%instid(i)%trans_jacobian)
                    deallocate(iv%instid(i)%der_trans)
                 end if
              end if
            end if

         end if

	  if ( use_rad .and. (use_varbc.or.freeze_varbc) ) then
	     if (iv%instid(i)%varbc_info%npredmax > 0) then
                deallocate (iv%instid(i)%varbc_info%pred)
                deallocate (iv%instid(i)%varbc_info%pred_mean)
                deallocate (iv%instid(i)%varbc_info%pred_std)
                deallocate (iv%instid(i)%varbc_info%nbgerr)
	     end if
	     do ichan = 1, iv%instid(i)%nchan
	        if (iv%instid(i)%varbc(ichan)%npred <= 0) cycle	 
	        deallocate (iv%instid(i)%varbc(ichan)%pred_use)
	        deallocate (iv%instid(i)%varbc(ichan)%ipred)
	        deallocate (iv%instid(i)%varbc(ichan)%index)
	        deallocate (iv%instid(i)%varbc(ichan)%param)
    	        deallocate (iv%instid(i)%varbc(ichan)%bgerr)
	        deallocate (iv%instid(i)%varbc(ichan)%vtox)
	     end do
   	     deallocate (iv%instid(i)%varbc)
	  end if
      end do
      deallocate (iv%instid)
      deallocate (j % jo % rad)     
      deallocate (satinfo)

   if (trace_use) call da_trace_exit ("da_deallocate_radiance")

end subroutine da_deallocate_radiance



end module da_radiance

