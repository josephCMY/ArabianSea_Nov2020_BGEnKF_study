












module da_obs_io

   use module_domain, only : domain

   use da_control, only : xmiss, missing_r, fmt_each, fmt_info, trace_use, &
      fmt_srfc, filtered_obs_unit, num_procs,missing, ierr,comm, rand_unit, &
      obs_qc_pointer, rootproc, omb_unit,omb_add_noise,use_airepobs, &
      use_airepobs,use_bogusobs,use_gpspwobs,use_gpsztdobs,use_gpsrefobs,use_geoamvobs, &
      use_metarobs,use_profilerobs,use_pilotobs,use_buoyobs,use_shipsobs,use_rainobs, &
      use_synopobs,use_soundobs,use_mtgirsobs,use_tamdarobs,use_qscatobs,use_radarobs, &
      test_transforms, use_ssmiretrievalobs, report_start, &
      report_end, global, print_detail_obs, stdout, t_kelvin, stderr, &
      max_ob_levels, missing_data, max_bogus_input, myproc, convert_uv2fd, convert_fd2uv, &
      fails_error_max,standard_atmosphere,zero_t_td,print_detail_f_obs, &
      print_detail_radar,use_satemobs,use_polaramvobs,use_ssmt1obs, &
      use_ssmt2obs, use_airsretobs,convert_fd2uv,anal_type_qcobs,gravity, &
      filename_len, t0, max_airep_input, max_bogus_input, max_ssmi_rv_input, &
      max_buoy_input, max_gpsref_input, max_gpspw_input, max_geoamv_input, &
      max_airsr_input, max_polaramv_input, max_radar_input, &
      max_profiler_input, max_sound_input, max_mtgirs_input, max_tamdar_input, max_ships_input, &
      max_satem_input,max_pilot_input, max_metar_input, max_ssmt1_input, &
      max_synop_input,max_ssmt2_input,  max_qscat_input, &
      obs_names, num_ob_indexes, fm_index, ids,ide, ite, jte, &
      sound, mtgirs,synop, pilot, satem, geoamv, polaramv, airep, gpspw, gpsref, &
      tamdar, tamdar_sfc, metar, ships, ssmi_rv, ssmi_tb, ssmt1, ssmt2, qscat, profiler, buoy, bogus, pseudo, &
      radar, radiance, airsr, sonde_sfc, trace_use_dull, num_fgat_time, time_slots, myproc, &
      qmarker_retain, anal_type_verify, top_km_gpsro, bot_km_gpsro, thin_rainobs, &
      sfc_assi_options, sfc_assi_options_1, sfc_assi_options_2,print_detail_rain,max_rain_input,rain, &
      pi, ob_format_gpsro, ob_format_ascii, analysis_date, kms,kme, v_interp_h,v_interp_p, &
      wind_sd,wind_sd_synop,wind_sd_tamdar,wind_sd_mtgirs,wind_sd_profiler,wind_sd_geoamv,wind_sd_polaramv, &
      wind_sd_airep,wind_sd_sound,wind_sd_metar,wind_sd_ships,wind_sd_qscat,wind_sd_buoy,wind_sd_pilot,wind_stats_sd,&
      thin_conv, thin_conv_ascii

   use da_define_structures, only : iv_type, multi_level_type, multi_level_type_BUFR, &
      radar_multi_level_type, y_type, field_type, each_level_type, &
      radar_each_level_type, info_type, model_loc_type,gpsref_type, rain_single_level_type, rain_each_type
   use da_grid_definitions, only : da_ffdduv,da_ffdduv_model,da_ffdduv_diagnose
   use da_obs, only : da_count_filtered_obs,da_check_missing,da_obs_proc_station, da_set_obs_missing, da_set_3d_obs_missing
   use da_par_util1, only : da_proc_sum_int
   use da_physics, only : da_tp_to_qs
   use da_reporting, only : da_warning, message, da_error
   use da_tools, only : da_llxy, da_get_julian_time, da_geo2msl1, da_msl2geo1
   use da_tools_serial, only : da_free_unit, da_get_unit, da_advance_time
   use da_tracing, only : da_trace_entry, da_trace_exit

   use module_radiance, only : deg2rad, i_kind
   use gsi_thinning, only : map2grids, map2grids_conv, cleangrids_conv, thinning_grid, &
                            map2tgrid, thinning_grid_conv
   use da_control, only : root

   use da_par_util, only : true_mpi_real
   use da_grid_definitions, only : da_earth_2_model_wind
   use da_reporting, only : message, da_message
   use da_interpolation, only : da_to_zk

   implicit none

   include 'mpif.h'

contains

subroutine da_read_obs_ascii (iv, filename, uvq_direct, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Read a GTS observation file
   !
   ! Modifications:      
   !      10/26/2004 - Syed RH Rizvi
   !      Fix Obs-Long as -180 if it is +180
   !
   !      03/04/2005 - Syed RH Rizvi
   !      (a)  AMVs from Geostationary and Polar orbiting satellite are
   !           seperated & used as profile. Currently there is a provision
   !           to use MODIS winds only.
   !      (b)  MODIS obs error are currently assigned through namelist
   !           parameter (modis_cmv_error)
   !
   !      03/30/2005 - Syed RH Rizvi
   !      For global option duplicate obs at East/West boundary
   !                        
   !      08/30/2005 - Syed RH Rizvi
   !      Writing original errors & thicknesses desired for 
   !      outputting QC obs with NTMAX = 0
   !
   !      02/28/2008 - Y.-R. Guo
   !      Satem thickness error should not be divided by gravity because
   !      the unit from obsproc is already meter, not geopotential meter.
   !
   !      03/19/2009 - Y.-R. Guo
   !      Added the time range check when reading in observations.
   !
   !      02/21/2013 - Syed RH Rizvi
   !      Updated with observation-thinning option
   !
   !      03/07/2014 - Feng Gao
   !      Updated wind_sd option for wind obs. types
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(len=*),  intent(in)    :: filename
   type(domain),     intent(in)     :: grid     ! first guess state.

   character (len =  10)        :: fmt_name

   character (len = 160)        :: info_string

!   character (len = 160)        :: fmt_info, fmt_srfc, fmt_each

   integer                      :: i, j, iost, nlevels, old_nlevels, fm,iunit

   type (multi_level_type)      :: platform
   logical                      :: outside
   logical                      :: outside_all
   logical                      :: uvq_direct
   integer                      :: surface_level
   real                         :: height_error, u_comp, v_comp
   integer                      :: nlocal(num_ob_indexes)
   integer                      :: ilocal(num_ob_indexes)
   integer                      :: ntotal(num_ob_indexes)

   integer                      :: ndup, n, report, obs_index

   real*8                       :: obs_time, analysis_time
   integer                      :: iyear, imonth, iday, ihour, imin
   real                         :: tdiff, dlat_earth,dlon_earth,crit
   integer                      :: itt,itx,iout
   real, allocatable            :: in(:), out(:)
   logical                      :: found, iuse, thin_3d, is_surface
   integer                      :: i1,j1,k, levs
   real                         :: dx,dy,dxm,dym,zk
   real                         :: v_p(kms:kme),v_h(kms:kme)


   if (trace_use) call da_trace_entry("da_read_obs_ascii")
 ! Initialize counts
   ntotal(:) = iv%info(:)%ptotal(iv%time-1)
   nlocal(:) = iv%info(:)%plocal(iv%time-1)
   ilocal    = nlocal
   if ( thin_conv_ascii ) then
       do n = 1, num_ob_indexes
          if ( n == radar ) cycle
          call cleangrids_conv(n)
       end do
   end if
   ! open file
   ! =========

   call da_get_unit(iunit)
   open(unit   = iunit,     &
      FILE   = trim(filename), &
      FORM   = 'FORMATTED',  &
      ACCESS = 'SEQUENTIAL', &
      iostat =  iost,     &
      STATUS = 'OLD')

   if (iost /= 0) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iost," opening gts obs file "//trim(filename)
      call da_warning("da_read_obs_ascii.inc",100,message(1:1))
      call da_free_unit(iunit)
      if (trace_use) call da_trace_exit("da_read_obs_ascii")
      return
   end if

   ! read header
   ! ===========

   do
      read (unit = iunit, fmt = '(A)', iostat = iost) info_string
      if (iost /= 0) then
         call da_warning("da_read_obs_ascii.inc",112, &
            (/"Problem reading gts obs header, skipping file"/))
         if (trace_use) call da_trace_exit("da_read_obs_ascii")
         return
      end if

      if (info_string(1:6) == 'EACH  ') exit
   end do

   !  read formats
   !  ------------

   read (iunit, fmt = '(A,1X,A)', iostat = iost) &
      fmt_name, fmt_info, &
      fmt_name, fmt_srfc,  &
      fmt_name, fmt_each

   if (iost /= 0) then
      call da_warning("da_read_obs_ascii.inc",130, &
         (/"Problem reading gts obs file, skipping"/))
         if (trace_use) call da_trace_exit("da_read_obs_ascii")
      return
   end if

   if (print_detail_obs) then
      write(unit=stdout, fmt='(2a)') &
         'fmt_info=', fmt_info, &
         'fmt_srfc =', fmt_srfc,  &
         'fmt_each=', fmt_each
   end if

   ! skip line
   ! ----------

   ! read (iunit, fmt = '(a)') fmt_name
   read (iunit, fmt = '(a)') info_string
   if (info_string(1:21) == 'MODIFIED FILTERED OBS') then
      uvq_direct=.true.
   else
      uvq_direct=.false.
   endif

   ! loop over records
   ! -----------------

   report = 0 ! report number in file

   reports: &
   do

      report = report+1

      ! read station general info
      ! =============================

      read (iunit, fmt = fmt_info, iostat = iost) &
         platform%info%platform,    &
         platform%info%date_char,   &
         platform%info%name,        &
         platform%info%levels,      &
         platform%info%lat,         &
         platform%info%lon,         &
         platform%info%elv,         &
         platform%info%id

      if (iost /= 0) then
         ! FIX? This is expected, but its unclear how we handle failure
         ! here without assuming the fortran2003 convention on
         ! error statuses
         exit reports
      end if

      if (print_detail_obs) then
         write(unit=stdout, fmt=fmt_info) &
            platform%info%platform,    &
            platform%info%date_char,   &
            platform%info%name,        &
            platform%info%levels,      &
            platform%info%lat,         &
            platform%info%lon,         &
            platform%info%elv,         &
            platform%info%id
      end if

      if (platform%info%lon == 180.0  ) platform%info%lon =-180.000 
      ! Fix funny wind direction at Poles
      if (platform%info%lat < -89.9999 .or. platform%info%lat > 89.9999) then
         platform%info%lon = 0.0
      end if

      if (platform%info%platform(6:6) == ' ') then
         read(platform%info%platform(4:5), '(I2)') fm
      else
         read(platform%info%platform(4:6), '(I3)') fm
      end if

      ! read model location
      ! =========================

      read (iunit, fmt = fmt_srfc)  &
         platform%loc%slp%inv, platform%loc%slp%qc, platform%loc%slp%error, &
         platform%loc%pw%inv, platform%loc%pw%qc, platform%loc%pw%error

      ! read each level
      ! ---------------

      do i = 1, platform%info%levels
         platform%each (i) = each_level_type(missing_r, missing, -1.0, & ! height
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! u
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! v
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! p
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! t
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! q
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! rh
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! td
            field_type(missing_r, missing, missing_r, missing, missing_r))  ! speed

         read (unit = iunit, fmt = trim (fmt_each)) &
            platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
            platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
            platform%each(i)%speed%error, &
            ! Here the 'direction' is stored in platform%each (i)%v:
            platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
            platform%each(i)%height,       &
            platform%each(i)%height_qc,    &
            height_error,                  &
            platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
            platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
            platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error

         if (uvq_direct) then
            platform%each (i)%u = platform%each (i)%speed
            if (platform%each(i)%rh%inv .gt. missing_r)  &
               platform%each(i)%rh%inv=platform%each(i)%rh%inv / 1e3    !convert back to kg/kg
            if (platform%each(i)%rh%error .gt. 0.0)  &
               platform%each(i)%rh%error=platform%each(i)%rh%error / 1e3  !convert back to kg/kg
         end if

         if (print_detail_obs) then
            write(unit=stdout, fmt = '(a, i6)') 'i=', i
! because now the field_type included: inv, qc, error, sens, and imp (YRG, 02/23/2009):              
            write(unit=stdout, fmt = trim (fmt_each)) &
               platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
               platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
               platform%each(i)%speed%error,  &
               platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
               platform%each(i)%height,       &
               platform%each(i)%height_qc,    &
               height_error,                  &
               platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
               platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
               platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error
         end if

         ! Uncomment the following if errors in reported heights (test later):             
         ! if (fm == 97 .or. fm == 96 .or. fm == 42) then
         !    platform%each (i)%height_qc = -5 ! Aircraft heights are not above surface.
         ! end if

         ! To convert the wind speed and direction to 
         !      the model wind components: u, v

         if (.not. uvq_direct) then
         if (platform%each (i)%speed%qc /= missing_data .and. &
             platform%each (i)%v%qc /= missing_data) then

             call da_ffdduv (platform%each (i)%speed%inv, platform%each (i)%v%inv,     &
                             U_comp, V_comp, platform%info%lon, convert_fd2uv)
             platform%each (i)%u = platform%each (i)%speed
             platform%each (i)%v = platform%each (i)%v
             platform%each (i)%u%inv = U_comp
             platform%each (i)%v%inv = V_comp
            
         else
            platform%each (i)%u%inv   = missing_r
            platform%each (i)%v%inv   = missing_r
            platform%each (i)%u%error = missing_r
            platform%each (i)%v%error = missing_r
            platform%each (i)%u%qc    = missing_data
            platform%each (i)%v%qc    = missing_data
         end if
         end if
      end do

      ! Check if outside of the time range:

      read (platform%info%date_char,'(i4,4(1x,i2))') &
                                    iyear, imonth, iday, ihour, imin
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
      if ( obs_time < time_slots(0) .or. &
           obs_time >= time_slots(num_fgat_time) ) then
         if (print_detail_obs) then
           write(unit=stdout, fmt='(a)') '*** Outside of the time range:'
           write(unit=stdout, fmt=fmt_info) &
            platform%info%platform,    &
            platform%info%date_char,   &
            platform%info%name,        &
            platform%info%levels,      &
            platform%info%lat,         &
            platform%info%lon,         &
            platform%info%elv,         &
            platform%info%id
         end if
         cycle
      endif

      ! Restrict to a range of reports, useful for debugging

      if (report < report_start) then
         cycle
      end if

      if (report > report_end) then
         exit
      end if

      call da_llxy (platform%info, platform%loc, outside, outside_all)

      if (outside_all) then
         cycle reports
      end if

      if (print_detail_obs) then
         ! Simplistic approach, could be improved to get it all done on PE 0
         if (.NOT. outside) then
            write(unit=stdout,fmt='(A,I5,A,F8.2,A,F8.2,A,I3,A,2F8.2)') &
               "Report",report," at",platform%info%lon,"E",platform%info%lat, &
               "N on processor", myproc,", position", platform%loc%x,platform%loc%y
         end if
      end if

      read (analysis_date,'(i4,4(1x,i2))') &
                                    iyear, imonth, iday, ihour, imin
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,analysis_time)
      tdiff = abs((obs_time - analysis_time)-0.02)
      dlat_earth = platform%info%lat
      dlon_earth = platform%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad

      nlevels = platform%info%levels
      if (nlevels > max_ob_levels) then
         nlevels = max_ob_levels
         write(unit=message(1), fmt='(/4(2a,2x),a,2f8.2,2x,2(a,f9.2,2x),2(a,i4,2x)/)') &
            'Station: ', trim(platform%info%name), &
            'ID: ', trim(platform%info%id), &
            'Platform: ', trim(platform%info%platform), &
            'Date: ', trim(platform%info%date_char), &
            'At Loc(lat, lon): ', platform%info%lat, platform%info%lon, &
            'At elv: ', platform%info%elv, &
            'with pstar: ', platform%info%pstar, &
            'Has level: ', platform%info%levels, &
            'which is great than max_ob_levels: ', max_ob_levels

         write (unit=message(2), fmt = '(A,1X,A,1X,A,1X,I4,2f9.3,f9.1,1X,A)') &
            platform%info%name,        &
            platform%info%platform,    &
            platform%info%date_char,   &
            platform%info%levels,      &
            platform%info%lat,         &
            platform%info%lon,         &
            platform%info%elv,         &
            platform%info%id
         call da_warning("da_read_obs_ascii.inc",377,message(1:2))

         platform%info%levels = nlevels
      else if (nlevels < 1) then
         ! Not GPSPW and GPSZD data:
         if (fm /= 111 .and. fm /= 114) then
            cycle reports
         end if
      end if

      levs = nlevels

      if (fm /= 86) call da_obs_proc_station(platform,fm,uvq_direct)

      ! Loop over duplicating obs for global
      ndup = 1
      if (global .and. (platform%loc%i < ids .or. platform%loc%i >= ide)) ndup= 2

      ! It is possible that logic for counting obs is incorrect for the 
      ! global case with >1 MPI tasks due to obs duplication, halo, etc.  
      ! TBH:  20050913

      if (.not.outside) then
         if (print_detail_obs .and. ndup > 1) then
            write(unit=stdout, fmt = fmt_info) &
               platform%info%platform,    &
               platform%info%date_char,   &
               platform%info%name,        &
               platform%info%levels,      &
               platform%info%lat,         &
               platform%info%lon,         &
               platform%info%elv,         &
               platform%info%id

            write(unit=stdout, fmt = '(a,2i5,4e20.10)') &
               ' duplicating obs since loc% i,j,dx,dxm,dy & dym ', &
               platform%loc%i,  platform%loc%j,   &
               platform%loc%dx, platform%loc%dxm, &
              platform%loc%dy, platform%loc%dym
         end if
      end if
      
      obs_index = fm_index(fm)

      IF ( wind_sd ) THEN
         wind_sd_buoy   = .true.
         wind_sd_synop  = .true.
         wind_sd_ships  = .true.
         wind_sd_metar  = .true.
         wind_sd_sound  = .true.
         wind_sd_pilot  = .true.
         wind_sd_airep  = .true.
         wind_sd_qscat  = .true.
         wind_sd_tamdar = .true.
         wind_sd_geoamv = .true.
         wind_sd_mtgirs = .true.
       wind_sd_polaramv = .true.
       wind_sd_profiler = .true.
      END IF

      dup_loop: do n = 1, ndup
         is_surface=.false.
         select case(fm)

         case (12) ;
            if (.not. use_synopobs .or. ntotal(synop) == max_synop_input  ) cycle reports
            if (n==1) ntotal(synop) = ntotal(synop)+1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(synop,dlat_earth,dlon_earth,crit,nlocal(synop),itx,1,itt,ilocal(synop),iuse)
               if ( .not. iuse ) cycle reports
            else
               nlocal(synop) = nlocal(synop) + 1
               ilocal(synop) = ilocal(synop) + 1
            end if

            if (.not. wind_sd_synop) platform%each(1)%v%error = platform%each(1)%u%error

            iv%synop(ilocal(synop))%h = platform%each(1)%height
            iv%synop(ilocal(synop))%u = platform%each(1)%u
            iv%synop(ilocal(synop))%v = platform%each(1)%v
            iv%synop(ilocal(synop))%t = platform%each(1)%t
            iv%synop(ilocal(synop))%q = platform%each(1)%q
            iv%synop(ilocal(synop))%p = platform%each(1)%p

            if (wind_sd_synop .and. &
                platform%each(1)%u%qc /= missing_data .and. platform%each(1)%v%qc /= missing_data ) &
                call da_ffdduv_model (iv%synop(ilocal(synop))%u%inv, iv%synop(ilocal(synop))%v%inv, &
                                      platform%each(1)%u%inv, platform%each(1)%v%inv, convert_uv2fd)

         case (13, 17) ;                  ! ships 
            if (.not. use_shipsobs .or. ntotal(ships) == max_ships_input  ) cycle reports
            if (n==1) ntotal(ships) = ntotal(ships) + 1 
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(ships,dlat_earth,dlon_earth,crit,nlocal(ships),itx,1,itt,ilocal(ships),iuse)
                if ( .not. iuse ) cycle reports
             else
                nlocal(ships) = nlocal(ships) + 1
                ilocal(ships) = ilocal(ships) + 1
             end if

            if (.not. wind_sd_ships) platform%each(1)%v%error = platform%each(1)%u%error

            iv%ships(ilocal(ships))%h = platform%each(1)%height
            iv%ships(ilocal(ships))%u = platform%each(1)%u
            iv%ships(ilocal(ships))%v = platform%each(1)%v
            iv%ships(ilocal(ships))%t = platform%each(1)%t
            iv%ships(ilocal(ships))%p = platform%each(1)%p
            iv%ships(ilocal(ships))%q = platform%each(1)%q

            if (wind_sd_ships .and. &
                platform%each(1)%u%qc /= missing_data .and. platform%each(1)%v%qc /= missing_data ) & 
                call da_ffdduv_model (iv%ships(ilocal(ships))%u%inv, iv%ships(ilocal(ships))%v%inv, &
                                      platform%each(1)%u%inv, platform%each(1)%v%inv, convert_uv2fd)

         case (15:16) ;
            if (.not. use_metarobs .or. ntotal(metar) == max_metar_input  ) cycle reports
            if (n==1) ntotal(metar) = ntotal(metar) + 1 
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(metar,dlat_earth,dlon_earth,crit,nlocal(metar),itx,1,itt,ilocal(metar),iuse)
                if ( .not. iuse ) cycle reports
             else
                nlocal(metar) = nlocal(metar) + 1
                ilocal(metar) = ilocal(metar) + 1
             end if
 
            if (.not. wind_sd_metar) platform%each(1)%v%error = platform%each(1)%u%error

            iv%metar(ilocal(metar))%h = platform%each(1)%height
            iv%metar(ilocal(metar))%u = platform%each(1)%u
            iv%metar(ilocal(metar))%v = platform%each(1)%v
            iv%metar(ilocal(metar))%t = platform%each(1)%t
            iv%metar(ilocal(metar))%p = platform%each(1)%p
            iv%metar(ilocal(metar))%q = platform%each(1)%q

            if (wind_sd_metar .and. &
                platform%each(1)%u%qc /= missing_data .and. platform%each(1)%v%qc /= missing_data ) & 
                call da_ffdduv_model (iv%metar(ilocal(metar))%u%inv, iv%metar(ilocal(metar))%v%inv, &
                                      platform%each(1)%u%inv, platform%each(1)%v%inv, convert_uv2fd)

         case (32:34) ;
            if (.not. use_pilotobs .or. ntotal(pilot) == max_pilot_input  ) cycle reports
            if (n==1) ntotal(pilot) = ntotal(pilot) + 1 
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(pilot,dlat_earth,dlon_earth,crit,nlocal(pilot),itx,1,itt,ilocal(pilot),iuse)
               if ( .not. iuse ) cycle reports
            else
               nlocal(pilot) = nlocal(pilot) + 1
               ilocal(pilot) = ilocal(pilot) + 1
            end if

            if (nlocal(pilot) == ilocal(pilot)) then
               allocate (iv%pilot(ilocal(pilot))%h(1:iv%info(pilot)%max_lev))
               allocate (iv%pilot(ilocal(pilot))%p(1:iv%info(pilot)%max_lev))
               allocate (iv%pilot(ilocal(pilot))%u(1:iv%info(pilot)%max_lev))
               allocate (iv%pilot(ilocal(pilot))%v(1:iv%info(pilot)%max_lev))
            end if
            do i = 1, nlevels

               if (.not. wind_sd_pilot) platform%each(i)%v%error = platform%each(i)%u%error

               iv%pilot(ilocal(pilot))%h(i) = platform%each(i)%height
               iv%pilot(ilocal(pilot))%p(i) = platform%each(i)%p%inv
               iv%pilot(ilocal(pilot))%u(i) = platform%each(i)%u
               iv%pilot(ilocal(pilot))%v(i) = platform%each(i)%v

               if ( .not. wind_sd_pilot ) platform%each(i)%v%error = platform%each(i)%u%error

               if (wind_sd_pilot .and. &
                   platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) then
                   call da_ffdduv_model (iv%pilot(ilocal(pilot))%u(i)%inv,iv%pilot(ilocal(pilot))%v(i)%inv, &
                                         platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
               end if
            end do

         case (35:38) ;
            if (.not. use_soundobs .or. ntotal(sound) == max_sound_input  ) cycle reports
            if (n==1) ntotal(sound) = ntotal(sound) + 1 
            if (n==1) ntotal(sonde_sfc) = ntotal(sonde_sfc) + 1 
            if (outside) cycle reports

! determine if level corresponds to surface

            levs = nlevels
            surface_level = 0
            do i = 1, nlevels
               ! if (elevation and height are the same, it is surface.)
               if (platform%info%elv.ne.missing_r.and.platform%each(i)%height.ne.missing_r)then
                  if (abs(platform%info%elv - platform%each(i)%height) < 0.1) then
                     is_surface = .true.
                     surface_level = i
                     levs = nlevels - 1
                     exit
                  end if
               end if
            end do

            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(sound,dlat_earth,dlon_earth,crit,nlocal(sound),itx,1,itt,ilocal(sound),iuse)
               crit = tdiff
               call map2grids_conv(sonde_sfc,dlat_earth,dlon_earth,crit,nlocal(sonde_sfc),itx,1,itt,ilocal(sonde_sfc),iuse)
               if ( .not. iuse ) cycle reports
            else
               nlocal(sound) = nlocal(sound) + 1
               ilocal(sound) = ilocal(sound) + 1
               nlocal(sonde_sfc) = nlocal(sonde_sfc) + 1
               ilocal(sonde_sfc) = ilocal(sonde_sfc) + 1
            end if

            ! Search to see if we have surface obs.

            surface_level = 0

            do i = 1, nlevels
               ! if (elevation and height are the same, it is surface.)
               if (platform%info%elv.ne.missing_r.and.platform%each(i)%height.ne.missing_r)then
                  if (abs(platform%info%elv - platform%each(i)%height) < 0.1) then
                     surface_level = i

                     if (.not. wind_sd_sound) platform%each(i)%v%error = platform%each(i)%u%error

                     ! Save surface pressure.
                     iv%sonde_sfc(ilocal(sonde_sfc))%h = platform%each(i)%height
                     iv%sonde_sfc(ilocal(sonde_sfc))%u = platform%each(i)%u
                     iv%sonde_sfc(ilocal(sonde_sfc))%v = platform%each(i)%v
                     iv%sonde_sfc(ilocal(sonde_sfc))%t = platform%each(i)%t
                     iv%sonde_sfc(ilocal(sonde_sfc))%q = platform%each(i)%q
                     iv%sonde_sfc(ilocal(sonde_sfc))%p = platform%each(i)%p

                     if (wind_sd_sound .and. &
                         platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) & 
                         call da_ffdduv_model (iv%sonde_sfc(ilocal(sonde_sfc))%u%inv,iv%sonde_sfc(ilocal(sonde_sfc))%v%inv, &
                                               platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
                     exit
                  end if
               end if
            end do

            ! processing the sonde_sfc data:

            old_nlevels = nlevels

            if (surface_level > 0) then
                nlevels = nlevels - 1
            else
               iv%sonde_sfc(ilocal(sonde_sfc))%h = missing_r
               iv%sonde_sfc(ilocal(sonde_sfc))%u%inv   = missing_r
               iv%sonde_sfc(ilocal(sonde_sfc))%u%qc    = missing
               iv%sonde_sfc(ilocal(sonde_sfc))%u%error = abs(missing_r)
               iv%sonde_sfc(ilocal(sonde_sfc))%v = iv%sonde_sfc(ilocal(sonde_sfc))%u
               iv%sonde_sfc(ilocal(sonde_sfc))%t = iv%sonde_sfc(ilocal(sonde_sfc))%u
               iv%sonde_sfc(ilocal(sonde_sfc))%p = iv%sonde_sfc(ilocal(sonde_sfc))%u
               iv%sonde_sfc(ilocal(sonde_sfc))%q = iv%sonde_sfc(ilocal(sonde_sfc))%u
            end if

            if (nlevels > 0) then
               if( nlocal(sound) == ilocal(sound)) then
                  allocate (iv%sound(ilocal(sound))%h (1:iv%info(sound)%max_lev))   
                  allocate (iv%sound(ilocal(sound))%p (1:iv%info(sound)%max_lev))   
                  allocate (iv%sound(ilocal(sound))%u (1:iv%info(sound)%max_lev))   
                  allocate (iv%sound(ilocal(sound))%v (1:iv%info(sound)%max_lev))   
                  allocate (iv%sound(ilocal(sound))%t (1:iv%info(sound)%max_lev))   
                  allocate (iv%sound(ilocal(sound))%q (1:iv%info(sound)%max_lev))   
               end if

               j = 0
               do i = 1, old_nlevels
                  if (i == surface_level) cycle

                  j=j+1

                  if (.not. wind_sd_sound) platform%each(i)%v%error = platform%each(i)%u%error

                  iv%sound(ilocal(sound))%h(j) = platform%each(i)%height
                  iv%sound(ilocal(sound))%p(j) = platform%each(i)%p%inv
                  iv%sound(ilocal(sound))%u(j) = platform%each(i)%u
                  iv%sound(ilocal(sound))%v(j) = platform%each(i)%v
                  iv%sound(ilocal(sound))%t(j) = platform%each(i)%t
                  iv%sound(ilocal(sound))%q(j) = platform%each(i)%q

                  if (wind_sd_sound .and. &
                      platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                      call da_ffdduv_model (iv%sound(ilocal(sound))%u(j)%inv,iv%sound(ilocal(sound))%v(j)%inv, &
                                            platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
               end do

           end if


         case (101) ;
            if (.not. use_tamdarobs .or. ntotal(tamdar) == max_tamdar_input  ) cycle reports

! determine if level corresponds to surface
            levs = nlevels
            surface_level = 0
            do i = 1, nlevels
               ! if (elevation and height are the same, it is surface.)
               if (platform%info%elv.ne.missing_r.and.platform%each(i)%height.ne.missing_r)then
                  if (abs(platform%info%elv - platform%each(i)%height) < 0.1) then
                     is_surface = .true.
                     surface_level = i
                     levs = nlevels - 1
                    exit
                  end if
               end if
            end do

          if( is_surface) then
            if (n==1) ntotal(tamdar_sfc) = ntotal(tamdar_sfc) + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(tamdar_sfc,dlat_earth,dlon_earth,crit,nlocal(tamdar_sfc),itx,1,itt,ilocal(tamdar_sfc),iuse)
                if ( .not. iuse ) cycle reports
             else
                nlocal(tamdar_sfc) = nlocal(tamdar_sfc) + 1
                ilocal(tamdar_sfc) = ilocal(tamdar_sfc) + 1
             end if

            if (.not. wind_sd_tamdar) platform%each(surface_level)%v%error = platform%each(surface_level)%u%error

            iv%tamdar_sfc(ilocal(tamdar_sfc))%h = platform%each(surface_level)%height
            iv%tamdar_sfc(ilocal(tamdar_sfc))%u = platform%each(surface_level)%u
            iv%tamdar_sfc(ilocal(tamdar_sfc))%v = platform%each(surface_level)%v
            iv%tamdar_sfc(ilocal(tamdar_sfc))%t = platform%each(surface_level)%t
            iv%tamdar_sfc(ilocal(tamdar_sfc))%q = platform%each(surface_level)%q
            iv%tamdar_sfc(ilocal(tamdar_sfc))%p = platform%each(surface_level)%p

            if (wind_sd_tamdar .and. &
                platform%each(surface_level)%u%qc /= missing_data .and. platform%each(surface_level)%v%qc /= missing_data ) &
                call da_ffdduv_model (iv%tamdar_sfc(ilocal(tamdar_sfc))%u%inv,iv%tamdar_sfc(ilocal(tamdar_sfc))%v%inv, &
                                      platform%each(surface_level)%u%inv, platform%each(surface_level)%v%inv, convert_uv2fd)

           else ! not is_surface

            if ( levs > 0 .and. n==1) ntotal(tamdar) = ntotal(tamdar) + 1
            if (outside) cycle reports
            if( .not. thin_conv_ascii) then
                   nlocal(tamdar) = nlocal(tamdar) + 1
                   ilocal(tamdar) = ilocal(tamdar) + 1
                if( nlocal(tamdar) == ilocal(tamdar)) then
                allocate (iv%tamdar(ilocal(tamdar))%h (1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%p (1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%u (1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%v (1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%t (1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%q (1:iv%info(tamdar)%max_lev))
             end if

             do i = 1, nlevels

                if (.not. wind_sd_tamdar) platform%each(i)%v%error = platform%each(i)%u%error

                iv%tamdar(ilocal(tamdar))%h(i) = platform%each(i)%height
                iv%tamdar(ilocal(tamdar))%p(i) = platform%each(i)%p%inv
                iv%tamdar(ilocal(tamdar))%u(i) = platform%each(i)%u
                iv%tamdar(ilocal(tamdar))%v(i) = platform%each(i)%v
                iv%tamdar(ilocal(tamdar))%t(i) = platform%each(i)%t
                iv%tamdar(ilocal(tamdar))%q(i) = platform%each(i)%q

                if (wind_sd_tamdar .and. &
                    platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                    call da_ffdduv_model (iv%tamdar(ilocal(tamdar))%u(i)%inv,iv%tamdar(ilocal(tamdar))%v(i)%inv, &
                                          platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
             end do

            else ! if thin_conv_ascii==true
              thin_3d=.true.
              i1   = platform%loc%i
              j1   = platform%loc%j
              dx   = platform%loc%dx
              dy   = platform%loc%dy
              dxm  = platform%loc%dxm
              dym  = platform%loc%dym

                 do k=kms,kme
                  v_p(k) = dym*(dxm*grid%xb%p(i1,j1,k)+dx*grid%xb%p(i1+1,j1,k)) + &
                            dy*(dxm*grid%xb%p(i1,j1+1,k)+dx*grid%xb%p(i1+1,j1+1,k))
                 end do
                 do k=kms,kme
                  v_h(k) = dym*(dxm*grid%xb%h(i1,j1,k)+dx*grid%xb%h(i1+1,j1,k)) + &
                            dy*(dxm*grid%xb%h(i1,j1+1,k)+dx*grid%xb%h(i1+1,j1+1,k))
                 end do
              do k= 1, nlevels
               zk = missing_r
               if( platform%each(k)%p%qc  >= 0 ) then
                 call da_to_zk(platform%each(k)%p%inv, v_p, v_interp_p, zk)
               else if( platform%each(k)%height_qc  >= 0 ) then
                 call da_to_zk(platform%each(k)%height , v_h, v_interp_h, zk)
               else
                 write(unit=message(1),fmt='(A)')' For tamdar: neither height nor pressure qc is good'
                 call da_error("da_read_obs_ascii.inc",776,message(1:1))
               end if
               if ( zk == missing_r ) cycle
               crit = tdiff
               call map2grids_conv(tamdar,dlat_earth,dlon_earth,crit,nlocal(tamdar),itx,1,itt,ilocal(tamdar),iuse,zk,thin_3d)
               if ( .not. iuse ) cycle
               iv%info(tamdar)%levels(ilocal(tamdar))    = 1
               iv%info(tamdar)%name(ilocal(tamdar))      = platform%info%name
               iv%info(tamdar)%platform(ilocal(tamdar))  = platform%info%platform
               iv%info(tamdar)%id(ilocal(tamdar))        = platform%info%id
               iv%info(tamdar)%date_char(ilocal(tamdar)) = platform%info%date_char
               iv%info(tamdar)%lat(:,ilocal(tamdar))     = platform%info%lat
               iv%info(tamdar)%lon(:,ilocal(tamdar))     = platform%info%lon
               iv%info(tamdar)%elv(ilocal(tamdar))       = platform%info%elv
               iv%info(tamdar)%pstar(ilocal(tamdar))     = platform%info%pstar

               iv%info(tamdar)%slp(ilocal(tamdar))           = platform%loc%slp
               iv%info(tamdar)%pw(ilocal(tamdar))            = platform%loc%pw
               iv%info(tamdar)%x(:,ilocal(tamdar))           = platform%loc%x
               iv%info(tamdar)%y(:,ilocal(tamdar))           = platform%loc%y
               iv%info(tamdar)%i(:,ilocal(tamdar))           = platform%loc%i
               iv%info(tamdar)%j(:,ilocal(tamdar))           = platform%loc%j
               iv%info(tamdar)%dx(:,ilocal(tamdar))          = platform%loc%dx
               iv%info(tamdar)%dxm(:,ilocal(tamdar))         = platform%loc%dxm
               iv%info(tamdar)%dy(:,ilocal(tamdar))          = platform%loc%dy
               iv%info(tamdar)%dym(:,ilocal(tamdar))         = platform%loc%dym
               iv%info(tamdar)%proc_domain(:,ilocal(tamdar)) = platform%loc%proc_domain

               iv%info(tamdar)%obs_global_index(ilocal(tamdar)) = ntotal(tamdar)
               if( nlocal(tamdar) == ilocal(tamdar)) then
                allocate (iv%tamdar(ilocal(tamdar))%h(1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%p(1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%u(1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%v(1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%t(1:iv%info(tamdar)%max_lev))
                allocate (iv%tamdar(ilocal(tamdar))%q(1:iv%info(tamdar)%max_lev))
               end if

               do i = 1, 1

                  if (.not. wind_sd_tamdar) platform%each(k)%v%error = platform%each(k)%u%error

                  iv%tamdar(ilocal(tamdar))%h(i) = platform%each(k)%height
                  iv%tamdar(ilocal(tamdar))%p(i) = platform%each(k)%p%inv
                  iv%tamdar(ilocal(tamdar))%u(i) = platform%each(k)%u
                  iv%tamdar(ilocal(tamdar))%v(i) = platform%each(k)%v
                  iv%tamdar(ilocal(tamdar))%t(i) = platform%each(k)%t
                  iv%tamdar(ilocal(tamdar))%q(i) = platform%each(k)%q

                  if (wind_sd_tamdar .and. &
                      platform%each(k)%u%qc /= missing_data .and. platform%each(k)%v%qc /= missing_data ) &
                      call da_ffdduv_model (iv%tamdar(ilocal(tamdar))%u(i)%inv,iv%tamdar(ilocal(tamdar))%v(i)%inv, &
                                            platform%each(k)%u%inv, platform%each(k)%v%inv, convert_uv2fd)
               end do
              end do ! loop over k levels

            end if ! if over thin_conv_ascii

           end if ! if  is_surface

         case (161) ;
            if (.not. use_mtgirsobs .or. ntotal(mtgirs) == max_mtgirs_input ) cycle reports
            if (n==1) ntotal(mtgirs) = ntotal(mtgirs) + 1
            if (outside) cycle reports
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(mtgirs,dlat_earth,dlon_earth,crit,nlocal(mtgirs),itx,1,itt,ilocal(mtgirs),iuse)
                if ( .not. iuse ) cycle reports
             else
                nlocal(mtgirs) = nlocal(mtgirs) + 1
                ilocal(mtgirs) = ilocal(mtgirs) + 1
             end if

            if (nlevels > 0) then

               if( nlocal(mtgirs) == ilocal(mtgirs)) then
                  allocate (iv%mtgirs(ilocal(mtgirs))%h (1:iv%info(mtgirs)%max_lev))
                  allocate (iv%mtgirs(ilocal(mtgirs))%p (1:iv%info(mtgirs)%max_lev))
                  allocate (iv%mtgirs(ilocal(mtgirs))%u (1:iv%info(mtgirs)%max_lev))
                  allocate (iv%mtgirs(ilocal(mtgirs))%v (1:iv%info(mtgirs)%max_lev))
                  allocate (iv%mtgirs(ilocal(mtgirs))%t (1:iv%info(mtgirs)%max_lev))
                  allocate (iv%mtgirs(ilocal(mtgirs))%q (1:iv%info(mtgirs)%max_lev))
               end if

               j = 0
               do i = 1, nlevels

                  j=j+1

                  if (.not. wind_sd_mtgirs) platform%each(i)%v%error = platform%each(i)%u%error

                  iv%mtgirs(ilocal(mtgirs))%h(j) = platform%each(i)%height
                  iv%mtgirs(ilocal(mtgirs))%p(j) = platform%each(i)%p%inv
                  iv%mtgirs(ilocal(mtgirs))%u(j) = platform%each(i)%u
                  iv%mtgirs(ilocal(mtgirs))%v(j) = platform%each(i)%v
                  iv%mtgirs(ilocal(mtgirs))%t(j) = platform%each(i)%t
                  iv%mtgirs(ilocal(mtgirs))%q(j) = platform%each(i)%q

                  if(iv%mtgirs(ilocal(mtgirs))%q(j)%inv.ne.missing_r)then
                     iv%mtgirs(ilocal(mtgirs))%q(j)%inv = iv%mtgirs(ilocal(mtgirs))%q(j)%inv/1000.0
                  endif
                  if(iv%mtgirs(ilocal(mtgirs))%q(j)%error.ne.missing_r)then
                     iv%mtgirs(ilocal(mtgirs))%q(j)%error = iv%mtgirs(ilocal(mtgirs))%q(j)%error/1000.0/100.0
                  endif

                  if (wind_sd_mtgirs .and. &
                      platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                      call da_ffdduv_model (iv%mtgirs(ilocal(mtgirs))%u(j)%inv,iv%mtgirs(ilocal(mtgirs))%v(j)%inv, &
                                            platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
               end do
            end if

         case (86)    ;
            ! Reject cloudy satem obs.
            if (.not.use_satemobs .or. ntotal(satem) == max_satem_input) cycle reports
            if (platform%loc%pw%inv > 10.0) cycle reports
            if (n==1) ntotal(satem) = ntotal(satem) + 1
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(satem,dlat_earth,dlon_earth,crit,nlocal(satem),itx,1,itt,ilocal(satem),iuse)
                if ( .not. iuse ) cycle reports
             else
                nlocal(satem) = nlocal(satem) + 1
                ilocal(satem) = ilocal(satem) + 1
             end if


            ! The satem ref_p is put in the SLP position in OBSPROC data stream.

            iv%satem(nlocal(satem))%ref_p= platform%loc%slp%inv

            if( nlocal(satem) == ilocal(satem)) then
             allocate (iv%satem(ilocal(satem))%p        (1:iv%info(satem)%max_lev))
             allocate (iv%satem(ilocal(satem))%thickness(1:iv%info(satem)%max_lev))
             allocate (iv%satem(ilocal(satem))%org_thickness(1:iv%info(satem)%max_lev))
            end if

            iv%satem(ilocal(satem))%p(1) = platform%each(1)%p%inv
            iv%satem(ilocal(satem))%thickness(1) = platform%each(1)%t
            !  write original thickness errors for filtered obs
            if (anal_type_qcobs) then
               do i = 1, nlevels
                  iv%satem(ilocal(satem))%org_thickness(i) = platform%each(i)%t
               end do
            end if

            ! Splitting the reported satem data into smaller layers.

            do i = 2, nlevels
               iv%satem(ilocal(satem))%p(i) = platform%each(i)%p%inv
               iv%satem(ilocal(satem))%thickness(i) = platform%each(i)%t
               if (platform%each(i)%t%qc /= missing_data   .and. &
                  platform%each(i-1)%t%qc /= missing_data) then
                  iv%satem(ilocal(satem))%thickness(i)%inv =            &
                  platform%each(i)%t%inv - platform%each(i-1)%t%inv
               else
                  iv%satem(ilocal(satem))%thickness(i)%inv = missing_r
                  iv%satem(ilocal(satem))%thickness(i)%qc  = missing_data
               end if
            end do

            ! Thickness error (m):

            do i = nlevels, 2, -1
               iv%satem(ilocal(satem))%thickness(i)%error = &
               sqrt(iv%satem(ilocal(satem))%thickness(i )%error ** 2 + &
                  iv%satem(ilocal(satem))%thickness(i-1)%error ** 2)
!                  iv%satem(ilocal(satem))%thickness(i-1)%error ** 2) / gravity
            end do

            iv%satem(ilocal(satem))%thickness(1)%error = &
                           sqrt(iv%satem(ilocal(satem))%thickness(1)%error ** 2 + &
                           platform%loc%pw%error ** 2)
!                           platform%loc%pw%error ** 2) / gravity


            ! Geostationary ot Polar orbitting Satellite AMVs:
         case (88)    ;
            if (index(platform%info%name, 'MODIS') > 0 .or. &
                index(platform%info%name, 'modis') > 0 .or. &
                index(platform%info%id, 'AVHRR') > 0)  then
               if (.not.use_polaramvobs .or. ntotal(polaramv) == max_polaramv_input) cycle reports
               if (n==1) ntotal(polaramv) = ntotal(polaramv) + 1
               if (outside) cycle reports
               if (n==1) ntotal(Polaramv) = ntotal(polaramv) + 1
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(polaramv,dlat_earth,dlon_earth,crit,nlocal(polaramv),itx,1,itt,ilocal(polaramv),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(polaramv) = nlocal(polaramv) + 1
                   ilocal(polaramv) = ilocal(polaramv) + 1
                end if

               if (nlocal(polaramv) == ilocal(polaramv) ) then
                allocate (iv%polaramv(ilocal(polaramv))%p (1:iv%info(polaramv)%max_lev))
                allocate (iv%polaramv(ilocal(polaramv))%u (1:iv%info(polaramv)%max_lev))
                allocate (iv%polaramv(ilocal(polaramv))%v (1:iv%info(polaramv)%max_lev))
               end if

               do i = 1, nlevels

                  if (.not. wind_sd_polaramv) platform%each(i)%v%error = platform%each(i)%u%error

                  iv%polaramv(ilocal(polaramv))%p(i) = platform%each(i)%p%inv
                  iv%polaramv(ilocal(polaramv))%u(i) = platform%each(i)%u
                  iv%polaramv(ilocal(polaramv))%v(i) = platform%each(i)%v

                  if (wind_sd_polaramv .and. &
                      platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                      call da_ffdduv_model (iv%polaramv(ilocal(polaramv))%u(i)%inv,iv%polaramv(ilocal(polaramv))%v(i)%inv, &
                                            platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
               end do
               obs_index = polaramv ! geoamv is the fm_index value for 88
            else
               if (.not.use_geoamvobs .or. ntotal(geoamv) == max_geoamv_input) cycle reports
               if (n==1) ntotal(geoamv) = ntotal(geoamv) + 1
               if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(geoamv,dlat_earth,dlon_earth,crit,nlocal(geoamv),itx,1,itt,ilocal(geoamv),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(geoamv) = nlocal(geoamv) + 1
                   ilocal(geoamv) = ilocal(geoamv) + 1
                end if

               if( nlocal(geoamv) == ilocal(geoamv) )then
                allocate (iv%geoamv(ilocal(geoamv))%p (1:iv%info(geoamv)%max_lev))
                allocate (iv%geoamv(ilocal(geoamv))%u (1:iv%info(geoamv)%max_lev))
                allocate (iv%geoamv(ilocal(geoamv))%v (1:iv%info(geoamv)%max_lev))
               end if

               do i = 1, nlevels

                  if (.not. wind_sd_geoamv) platform%each(i)%v%error = platform%each(i)%u%error

                  iv%geoamv(ilocal(geoamv))%p(i) = platform%each(i)%p%inv
                  iv%geoamv(ilocal(geoamv))%u(i) = platform%each(i)%u
                  iv%geoamv(ilocal(geoamv))%v(i) = platform%each(i)%v

                  if (wind_sd_geoamv .and. &
                      platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                      call da_ffdduv_model (iv%geoamv(ilocal(geoamv))%u(i)%inv,iv%geoamv(ilocal(geoamv))%v(i)%inv, &
                                            platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
               end do

            end if

         case (42,96:97) ;

            if (.not.use_airepobs .or. ntotal(airep) == max_airep_input) cycle reports
            if (n==1) ntotal(airep) = ntotal(airep) + 1
            if (outside) cycle reports

            if( .not. thin_conv_ascii) then
                   nlocal(airep) = nlocal(airep) + 1
                   ilocal(airep) = ilocal(airep) + 1
              if( nlocal(airep) == ilocal(airep)) then
              allocate (iv%airep(ilocal(airep))%h        (1:iv%info(airep)%max_lev))
              allocate (iv%airep(ilocal(airep))%p        (1:iv%info(airep)%max_lev))
              allocate (iv%airep(ilocal(airep))%u        (1:iv%info(airep)%max_lev))
              allocate (iv%airep(ilocal(airep))%v        (1:iv%info(airep)%max_lev))
              allocate (iv%airep(ilocal(airep))%t        (1:iv%info(airep)%max_lev))
              allocate (iv%airep(ilocal(airep))%q        (1:iv%info(airep)%max_lev))
             end if

             do i = 1, nlevels

                if (.not. wind_sd_airep) platform%each(i)%v%error = platform%each(i)%u%error

                iv%airep(ilocal(airep))%h(i) = platform%each(i)%height
                iv%airep(ilocal(airep))%p(i) = platform%each(i)%p%inv
                iv%airep(ilocal(airep))%u(i) = platform%each(i)%u
                iv%airep(ilocal(airep))%v(i) = platform%each(i)%v
                iv%airep(ilocal(airep))%t(i) = platform%each(i)%t
                iv%airep(ilocal(airep))%q(i) = platform%each(i)%q

                if (wind_sd_airep .and. &
                    platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                    call da_ffdduv_model (iv%airep(ilocal(airep))%u(i)%inv,iv%airep(ilocal(airep))%v(i)%inv, &
                                          platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
             end do

            else ! if thin_conv_ascii==true
              thin_3d=.true.
              i1   = platform%loc%i
              j1   = platform%loc%j
              dx   = platform%loc%dx
              dy   = platform%loc%dy
              dxm  = platform%loc%dxm
              dym  = platform%loc%dym

                 do k=kms,kme
                  v_p(k) = dym*(dxm*grid%xb%p(i1,j1,k)+dx*grid%xb%p(i1+1,j1,k)) + &
                            dy*(dxm*grid%xb%p(i1,j1+1,k)+dx*grid%xb%p(i1+1,j1+1,k))
                 end do
                 do k=kms,kme
                  v_h(k) = dym*(dxm*grid%xb%h(i1,j1,k)+dx*grid%xb%h(i1+1,j1,k)) + &
                            dy*(dxm*grid%xb%h(i1,j1+1,k)+dx*grid%xb%h(i1+1,j1+1,k))
                 end do
              do k= 1, nlevels
               zk = missing_r
               if( platform%each(k)%p%qc  >= 0 ) then
                 call da_to_zk(platform%each(k)%p%inv, v_p, v_interp_p, zk)
               else if( platform%each(k)%height_qc  >= 0 ) then
                 call da_to_zk(platform%each(k)%height , v_h, v_interp_h, zk)
               else
                 write(unit=message(1),fmt='(A)')' For airep: neither height nor pressure qc is good'
                 call da_error("da_read_obs_ascii.inc",1087,message(1:1))
               end if
               if ( zk == missing_r ) cycle
               crit = tdiff
               call map2grids_conv(airep,dlat_earth,dlon_earth,crit,nlocal(airep),itx,1,itt,ilocal(airep),iuse,zk,thin_3d)
               if ( .not. iuse ) cycle
               iv%info(airep)%levels(ilocal(airep))    = 1
               iv%info(airep)%name(ilocal(airep))      = platform%info%name
               iv%info(airep)%platform(ilocal(airep))  = platform%info%platform
               iv%info(airep)%id(ilocal(airep))        = platform%info%id
               iv%info(airep)%date_char(ilocal(airep)) = platform%info%date_char
               iv%info(airep)%lat(:,ilocal(airep))     = platform%info%lat
               iv%info(airep)%lon(:,ilocal(airep))     = platform%info%lon
               iv%info(airep)%elv(ilocal(airep))       = platform%info%elv
               iv%info(airep)%pstar(ilocal(airep))     = platform%info%pstar

               iv%info(airep)%slp(ilocal(airep))           = platform%loc%slp
               iv%info(airep)%pw(ilocal(airep))            = platform%loc%pw
               iv%info(airep)%x(:,ilocal(airep))           = platform%loc%x
               iv%info(airep)%y(:,ilocal(airep))           = platform%loc%y
               iv%info(airep)%i(:,ilocal(airep))           = platform%loc%i
               iv%info(airep)%j(:,ilocal(airep))           = platform%loc%j
               iv%info(airep)%dx(:,ilocal(airep))          = platform%loc%dx
               iv%info(airep)%dxm(:,ilocal(airep))         = platform%loc%dxm
               iv%info(airep)%dy(:,ilocal(airep))          = platform%loc%dy
               iv%info(airep)%dym(:,ilocal(airep))         = platform%loc%dym
               iv%info(airep)%proc_domain(:,ilocal(airep)) = platform%loc%proc_domain

               iv%info(airep)%obs_global_index(ilocal(airep)) = ntotal(airep)
               if( nlocal(airep) == ilocal(airep)) then
                allocate (iv%airep(ilocal(airep))%h        (1:iv%info(airep)%max_lev))
                allocate (iv%airep(ilocal(airep))%p        (1:iv%info(airep)%max_lev))
                allocate (iv%airep(ilocal(airep))%u        (1:iv%info(airep)%max_lev))
                allocate (iv%airep(ilocal(airep))%v        (1:iv%info(airep)%max_lev))
                allocate (iv%airep(ilocal(airep))%t        (1:iv%info(airep)%max_lev))
                allocate (iv%airep(ilocal(airep))%q        (1:iv%info(airep)%max_lev))
               end if

               do i = 1, 1

                  if (.not. wind_sd_airep) platform%each(k)%v%error = platform%each(k)%u%error

                  iv%airep(ilocal(airep))%h(i) = platform%each(k)%height
                  iv%airep(ilocal(airep))%p(i) = platform%each(k)%p%inv
                  iv%airep(ilocal(airep))%u(i) = platform%each(k)%u
                  iv%airep(ilocal(airep))%v(i) = platform%each(k)%v
                  iv%airep(ilocal(airep))%t(i) = platform%each(k)%t
                  iv%airep(ilocal(airep))%q(i) = platform%each(k)%q

                  if (wind_sd_airep .and. &
                      platform%each(k)%u%qc /= missing_data .and. platform%each(k)%v%qc /= missing_data ) &
                      call da_ffdduv_model (iv%airep(ilocal(airep))%u(i)%inv,iv%airep(ilocal(airep))%v(i)%inv, &
                                            platform%each(k)%u%inv, platform%each(k)%v%inv, convert_uv2fd)

               end do
              end do ! loop over k levels

            end if ! if over thin_conv_ascii

         case (111, 114) ;
            if((.not.use_gpspwobs  .and. fm == 111) .or. ntotal(gpspw) == max_gpspw_input ) cycle reports
            if((.not.use_gpsztdobs .and. fm == 114) .or. ntotal(gpspw) == max_gpspw_input ) cycle reports
            if (n==1) ntotal(gpspw) = ntotal(gpspw) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(gpspw,dlat_earth,dlon_earth,crit,nlocal(gpspw),itx,1,itt,ilocal(gpspw),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(gpspw) = nlocal(gpspw) + 1
                   ilocal(gpspw) = ilocal(gpspw) + 1
                end if

            iv%gpspw(ilocal(gpspw))%tpw  = platform%loc%pw

         case (116) ;
            if(.not.use_gpsrefobs  .or. ntotal(gpsref) == max_gpsref_input ) cycle reports
            if ( ob_format_gpsro /= ob_format_ascii ) cycle reports
            if (n==1) ntotal(gpsref) = ntotal(gpsref) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(gpsref,dlat_earth,dlon_earth,crit,nlocal(gpsref),itx,1,itt,ilocal(gpsref),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(gpsref) = nlocal(gpsref) + 1
                   ilocal(gpsref) = ilocal(gpsref) + 1
                end if
!
! discarded the GPSRO data above 30-km (refer to GSI, Lidia Cucurull, Nov. 2009).
! YRG, 02/01/2010.
            do i = 1, nlevels
               if ( platform%each(i)%height> 30000.0 ) then
                  nlevels = i-1
                  exit
               endif
            enddo
            levs = nlevels
!
            if( nlocal(gpsref) == ilocal(gpsref)) then
               allocate (iv%gpsref(ilocal(gpsref))%h (1:iv%info(gpsref)%max_lev))
               allocate (iv%gpsref(ilocal(gpsref))%ref(1:iv%info(gpsref)%max_lev))
               allocate (iv%gpsref(ilocal(gpsref))%  p(1:iv%info(gpsref)%max_lev))
               allocate (iv%gpsref(ilocal(gpsref))%  t(1:iv%info(gpsref)%max_lev))
               allocate (iv%gpsref(ilocal(gpsref))%  q(1:iv%info(gpsref)%max_lev))
            end if

            do i = 1, nlevels
               iv%gpsref(ilocal(gpsref))%h(i)   = platform%each(i)%height

               ! In OBSPROC, use "td" field to store "gpsref"
               iv%gpsref(ilocal(gpsref))%ref(i) = platform%each(i)%td

               ! check height, only keep data below and above certain height
               if ( iv%gpsref(ilocal(gpsref))%h(i) > top_km_gpsro*1000.0 .or. &
                    iv%gpsref(ilocal(gpsref))%h(i) < bot_km_gpsro*1000.0 ) then
                  iv%gpsref(ilocal(gpsref))%ref(i)%qc = -77
               end if

               ! Keep the retrieved p and t (and q) as "field_type":
               iv%gpsref(ilocal(gpsref))%p(i)   = platform%each(i)%p
               iv%gpsref(ilocal(gpsref))%t(i)   = platform%each(i)%t
               iv%gpsref(ilocal(gpsref))%q(i)   = platform%each(i)%q
            end do

         case (121) ;
            if(.not.use_ssmt1obs  .or. ntotal(ssmt1) == max_ssmt1_input ) cycle reports
            if (n==1) ntotal(ssmt1) = ntotal(ssmt1) + 1
            if (outside) cycle reports
            if (n==1) ntotal(ssmt1) = ntotal(ssmt1) + 1
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(ssmt1,dlat_earth,dlon_earth,crit,nlocal(ssmt1),itx,1,itt,ilocal(ssmt1),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(ssmt1) = nlocal(ssmt1) + 1
                   ilocal(ssmt1) = ilocal(ssmt1) + 1
                end if

            if ( nlocal(ssmt1) == ilocal(ssmt1)) then
             allocate (iv%ssmt1(ilocal(ssmt1))%h (1:iv%info(ssmt1)%max_lev))
             allocate (iv%ssmt1(ilocal(ssmt1))%p (1:iv%info(ssmt1)%max_lev))
             allocate (iv%ssmt1(ilocal(ssmt1))%t (1:iv%info(ssmt1)%max_lev))
            end if

            do i = 1, nlevels
              iv%ssmt1(ilocal(ssmt1))%h(i) = platform%each(i)%height
              iv%ssmt1(ilocal(ssmt1))%p(i) = platform%each(i)%p%inv
              iv%ssmt1(ilocal(ssmt1))%t(i) = platform%each(i)%t
            end do

         case (122) ;
            if(.not.use_ssmt2obs  .or. ntotal(ssmt2) == max_ssmt2_input ) cycle reports
            if (n==1) ntotal(ssmt2) = ntotal(ssmt2) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(ssmt2,dlat_earth,dlon_earth,crit,nlocal(ssmt2),itx,1,itt,ilocal(ssmt2),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(ssmt2) = nlocal(ssmt2) + 1
                   ilocal(ssmt2) = ilocal(ssmt2) + 1
                end if
            if ( nlocal(ssmt2) == ilocal(ssmt2)) then
             allocate (iv%ssmt2(ilocal(ssmt2))%h (1:iv%info(ssmt2)%max_lev))
             allocate (iv%ssmt2(ilocal(ssmt2))%p (1:iv%info(ssmt2)%max_lev))
             allocate (iv%ssmt2(ilocal(ssmt2))%rh(1:iv%info(ssmt2)%max_lev))
            end if

            do i = 1, nlevels
               iv%ssmt2(ilocal(ssmt2))% h(i) = platform%each(i)%height
               iv%ssmt2(ilocal(ssmt2))% p(i) = platform%each(i)%p%inv
               iv%ssmt2(ilocal(ssmt2))%rh(i) = platform%each(i)%rh
            end do

         case (281)    ;
            if(.not.use_qscatobs  .or. ntotal(qscat) == max_qscat_input ) cycle reports
            if (n==1) ntotal(qscat) = ntotal(qscat) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(qscat,dlat_earth,dlon_earth,crit,nlocal(qscat),itx,1,itt,ilocal(qscat),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(qscat) = nlocal(qscat) + 1
                   ilocal(qscat) = ilocal(qscat) + 1
                end if

!            if (nlocal(qscat) == ilocal(qscat)) then

                if (.not. wind_sd_qscat) platform%each(1)%v%error = platform%each(1)%u%error

                iv%qscat(ilocal(qscat))%h = platform%each(1)%height
                iv%qscat(ilocal(qscat))%u = platform%each(1)%u
                iv%qscat(ilocal(qscat))%v = platform%each(1)%v

                if (wind_sd_qscat .and. &
                    platform%each(1)%u%qc /= missing_data .and. platform%each(1)%v%qc /= missing_data ) &
                    call da_ffdduv_model (iv%qscat(ilocal(qscat))%u%inv,iv%qscat(ilocal(qscat))%v%inv, &
                                          platform%each(1)%u%inv, platform%each(1)%v%inv, convert_uv2fd)

!            end if

            ! Impose minimum observation error = 1.0m/s for Quikscat data:
            iv%qscat(ilocal(qscat))%u%error = max(platform%each(1)%u%error,1.0)
            iv%qscat(ilocal(qscat))%v%error = max(platform%each(1)%v%error,1.0)

         case (132) ; ! profiler
            if (.not. use_profilerobs .or. ntotal(profiler) == max_profiler_input ) cycle reports
            if (n==1) ntotal(profiler) = ntotal(profiler) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(profiler,dlat_earth,dlon_earth,crit,nlocal(profiler),itx,1,itt,ilocal(profiler),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(profiler) = nlocal(profiler) + 1
                   ilocal(profiler) = ilocal(profiler) + 1
                end if

            if (nlocal(profiler) == ilocal(profiler) ) then
               allocate (iv%profiler(ilocal(profiler))%h (1:iv%info(profiler)%max_lev))
               allocate (iv%profiler(ilocal(profiler))%p (1:iv%info(profiler)%max_lev))
               allocate (iv%profiler(ilocal(profiler))%u (1:iv%info(profiler)%max_lev))
               allocate (iv%profiler(ilocal(profiler))%v (1:iv%info(profiler)%max_lev))
            end if

            do i = 1, nlevels

               if (.not. wind_sd_profiler) platform%each(i)%v%error = platform%each(i)%u%error

               iv%profiler(ilocal(profiler))%h(i) = platform%each(i)%height
               iv%profiler(ilocal(profiler))%p(i) = platform%each(i)%p%inv
               iv%profiler(ilocal(profiler))%u(i) = platform%each(i)%u
               iv%profiler(ilocal(profiler))%v(i) = platform%each(i)%v

               if (wind_sd_profiler .and. &
                   platform%each(i)%u%qc /= missing_data .and. platform%each(i)%v%qc /= missing_data ) &
                   call da_ffdduv_model (iv%profiler(ilocal(profiler))%u(i)%inv,iv%profiler(ilocal(profiler))%v(i)%inv, &
                                         platform%each(i)%u%inv, platform%each(i)%v%inv, convert_uv2fd)
            end do

         case (135) ; ! Bogus
            if (.not. use_bogusobs .or. ntotal(bogus) == max_bogus_input ) cycle reports
            if (n==1) ntotal(bogus) = ntotal(bogus) + 1
            if (outside) cycle reports
            if (n==1) ntotal(bogus) = ntotal(bogus) + 1
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(bogus,dlat_earth,dlon_earth,crit,nlocal(bogus),itx,1,itt,ilocal(bogus),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(bogus) = nlocal(bogus) + 1
                   ilocal(bogus) = ilocal(bogus) + 1
                end if

            if (ilocal(bogus) > max_bogus_input) then
               write(unit=message(1),fmt='(A,I6,A,I6)') &
                  'Bogus #=', ilocal(bogus), ' > max_bogus_input=', max_bogus_input
               call da_error("da_read_obs_ascii.inc",1346,message(1:1))
            end if

            if (nlocal(bogus) == ilocal(bogus)) then
               allocate (iv%bogus(ilocal(bogus))%h (1:iv%info(bogus)%max_lev))
               allocate (iv%bogus(ilocal(bogus))%p (1:iv%info(bogus)%max_lev))
               allocate (iv%bogus(ilocal(bogus))%u (1:iv%info(bogus)%max_lev))
               allocate (iv%bogus(ilocal(bogus))%v (1:iv%info(bogus)%max_lev))
               allocate (iv%bogus(ilocal(bogus))%t (1:iv%info(bogus)%max_lev))
               allocate (iv%bogus(ilocal(bogus))%q (1:iv%info(bogus)%max_lev))
            end if

            do i = 1, nlevels
               if (.not. wind_sd) platform%each(i)%v%error = platform%each(i)%u%error
               iv%bogus(ilocal(bogus))%h(i) = platform%each(i)%height
               iv%bogus(ilocal(bogus))%p(i) = platform%each(i)%p%inv
               iv%bogus(ilocal(bogus))%u(i) = platform%each(i)%u
               iv%bogus(ilocal(bogus))%v(i) = platform%each(i)%v
               iv%bogus(ilocal(bogus))%t(i) = platform%each(i)%t
               iv%bogus(ilocal(bogus))%q(i) = platform%each(i)%q
            end do

            iv%bogus(ilocal(bogus))%slp    = platform%loc%slp

            if (print_detail_obs) then
               write(unit=stdout,fmt=*) 'nlevels=', nlevels
               write(unit=stdout,fmt=*) 'iv%info(bogus)%nlocal,slp', ilocal(bogus),  &
                  iv % bogus (ilocal(bogus)) % slp
               do i=1,nlevels
                  write(unit=stdout,fmt=*) 'nlocal(bogus), i ', nlocal(bogus),i
                  write(unit=stdout,fmt=*) 'iv%bogus(nlocal(bogus))%u,v=',  &
                     iv%bogus(ilocal(bogus))%u(i),iv%bogus(ilocal(bogus))%v(i)
                  write(unit=stdout,fmt=*) 'iv%bogus(nlocal(bogus))%t,q=',  &
                     iv%bogus(ilocal(bogus))%t(i),iv%bogus(ilocal(bogus))%q(i)
               end do
               write(unit=stdout,fmt='(2(a,i4))') 'nlocal(bogus)=',ilocal(bogus), &
                  'nlevels=',nlevels
            end if

         case (18,19) ;             ! buoy
            if (.not. use_buoyobs .or. ntotal(buoy) == max_buoy_input ) cycle reports
            if (n==1) ntotal(buoy) = ntotal(buoy) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(buoy,dlat_earth,dlon_earth,crit,nlocal(buoy),itx,1,itt,ilocal(buoy),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(buoy) = nlocal(buoy)  + 1
                   ilocal(buoy) = ilocal(buoy)  + 1
                end if

            if (.not. wind_sd_buoy) platform%each(1)%v%error = platform%each(1)%u%error

            iv%buoy(ilocal(buoy))%h = platform%each(1)%height
            iv%buoy(ilocal(buoy))%u = platform%each(1)%u
            iv%buoy(ilocal(buoy))%v = platform%each(1)%v
            iv%buoy(ilocal(buoy))%t = platform%each(1)%t
            iv%buoy(ilocal(buoy))%p = platform%each(1)%p
            iv%buoy(ilocal(buoy))%q = platform%each(1)%q

            if (wind_sd_buoy .and. &
                platform%each(1)%u%qc /= missing_data .and. platform%each(1)%v%qc /= missing_data ) &
                call da_ffdduv_model (iv%buoy(ilocal(buoy))%u%inv,iv%buoy(ilocal(buoy))%v%inv, &
                                      platform%each(1)%u%inv, platform%each(1)%v%inv, convert_uv2fd)


         case (133)    ;         ! AIRS retrievals  
            if (.not. use_airsretobs .or. ntotal(airsr) == max_airsr_input ) cycle reports
            if (n==1) ntotal(airsr) = ntotal(airsr) + 1
            if (outside) cycle reports
                if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(airsr,dlat_earth,dlon_earth,crit,nlocal(airsr),itx,1,itt,ilocal(airsr),iuse)
                   if ( .not. iuse ) cycle reports
                else
                   nlocal(airsr) = nlocal(airsr)  + 1
                   ilocal(airsr) = ilocal(airsr)  + 1
                end if

            if (nlocal(airsr) == ilocal(airsr)) then
               allocate (iv%airsr(ilocal(airsr))%h (1:iv%info(airsr)%max_lev))
               allocate (iv%airsr(ilocal(airsr))%p (1:iv%info(airsr)%max_lev))
               allocate (iv%airsr(ilocal(airsr))%t (1:iv%info(airsr)%max_lev))
               allocate (iv%airsr(ilocal(airsr))%q (1:iv%info(airsr)%max_lev))
            end if
            do i = 1, nlevels
               iv%airsr(ilocal(airsr))%h(i) = platform%each(i)%height
               iv%airsr(ilocal(airsr))%p(i) = platform%each(i)%p%inv
               iv%airsr(ilocal(airsr))%t(i) = platform%each(i)%t
               iv%airsr(ilocal(airsr))%q(i) = platform%each(i)%q
            end do

         case default;

            write(unit=message(1), fmt='(a)') 'unsaved obs found:'
            write(unit=message(2), fmt='(2a)') &
               'platform%info%platform=', platform%info%platform
            write(unit=message(3), fmt='(a, i3)') &
               'platform%info%levels=', platform%info%levels
            call da_warning("da_read_obs_ascii.inc",1446,message(1:3))
            cycle
         end select
         if( is_surface .or. (obs_index == gpspw) .or. (levs > 0 .and. .not. thin_conv_ascii) .or. &
            (levs > 0 .and. (thin_conv_ascii .and. (obs_index /=  airep .and. obs_index /= tamdar))) ) then
            iv%info(obs_index)%name(ilocal(obs_index))      = platform%info%name
            iv%info(obs_index)%platform(ilocal(obs_index))  = platform%info%platform
            iv%info(obs_index)%id(ilocal(obs_index))        = platform%info%id
            iv%info(obs_index)%date_char(ilocal(obs_index)) = platform%info%date_char
            ! nlevels adjusted for some obs types so use that
            iv%info(obs_index)%levels(ilocal(obs_index))    = min(iv%info(obs_index)%max_lev, levs)
            iv%info(obs_index)%lat(:,ilocal(obs_index))     = platform%info%lat
            iv%info(obs_index)%lon(:,ilocal(obs_index))     = platform%info%lon
            iv%info(obs_index)%elv(ilocal(obs_index))       = platform%info%elv
            iv%info(obs_index)%pstar(ilocal(obs_index))     = platform%info%pstar

            iv%info(obs_index)%slp(ilocal(obs_index))           = platform%loc%slp
            iv%info(obs_index)%pw(ilocal(obs_index))            = platform%loc%pw
            iv%info(obs_index)%x(:,ilocal(obs_index))           = platform%loc%x
            iv%info(obs_index)%y(:,ilocal(obs_index))           = platform%loc%y
            iv%info(obs_index)%i(:,ilocal(obs_index))           = platform%loc%i
            iv%info(obs_index)%j(:,ilocal(obs_index))           = platform%loc%j
            iv%info(obs_index)%dx(:,ilocal(obs_index))          = platform%loc%dx
            iv%info(obs_index)%dxm(:,ilocal(obs_index))         = platform%loc%dxm
            iv%info(obs_index)%dy(:,ilocal(obs_index))          = platform%loc%dy
            iv%info(obs_index)%dym(:,ilocal(obs_index))         = platform%loc%dym
            iv%info(obs_index)%proc_domain(:,ilocal(obs_index)) = platform%loc%proc_domain

            iv%info(obs_index)%obs_global_index(nlocal(obs_index)) = ntotal(obs_index)

            ! special case for sonde_sfc, duplicate sound info
            if (obs_index == sound) then
               iv%info(sonde_sfc)%name(ilocal(sonde_sfc))      = platform%info%name
               iv%info(sonde_sfc)%platform(ilocal(sonde_sfc))  = platform%info%platform
               iv%info(sonde_sfc)%id(ilocal(sonde_sfc))        = platform%info%id
               iv%info(sonde_sfc)%date_char(ilocal(sonde_sfc)) = platform%info%date_char
               iv%info(sonde_sfc)%levels(ilocal(sonde_sfc))    = 1
               iv%info(sonde_sfc)%lat(:,ilocal(sonde_sfc))     = platform%info%lat
               iv%info(sonde_sfc)%lon(:,ilocal(sonde_sfc))     = platform%info%lon
               iv%info(sonde_sfc)%elv(ilocal(sonde_sfc))       = platform%info%elv
               iv%info(sonde_sfc)%pstar(ilocal(sonde_sfc))     = platform%info%pstar

               iv%info(sonde_sfc)%slp(ilocal(sonde_sfc))           = platform%loc%slp
               iv%info(sonde_sfc)%pw(ilocal(sonde_sfc))            = platform%loc%pw
               iv%info(sonde_sfc)%x(:,ilocal(sonde_sfc))           = platform%loc%x
               iv%info(sonde_sfc)%y(:,ilocal(sonde_sfc))           = platform%loc%y
               iv%info(sonde_sfc)%i(:,ilocal(sonde_sfc))           = platform%loc%i
               iv%info(sonde_sfc)%j(:,ilocal(sonde_sfc))           = platform%loc%j
               iv%info(sonde_sfc)%dx(:,ilocal(sonde_sfc))          = platform%loc%dx
               iv%info(sonde_sfc)%dxm(:,ilocal(sonde_sfc))         = platform%loc%dxm
               iv%info(sonde_sfc)%dy(:,ilocal(sonde_sfc))          = platform%loc%dy
               iv%info(sonde_sfc)%dym(:,ilocal(sonde_sfc))         = platform%loc%dym
               iv%info(sonde_sfc)%proc_domain(:,ilocal(sonde_sfc)) = platform%loc%proc_domain
               iv%info(sonde_sfc)%obs_global_index(ilocal(sonde_sfc)) = ntotal(obs_index)
            end if

            if (is_surface .and. obs_index == tamdar) then
               iv%info(tamdar_sfc)%name(ilocal(tamdar_sfc))      = platform%info%name
               iv%info(tamdar_sfc)%platform(ilocal(tamdar_sfc))  = platform%info%platform
               iv%info(tamdar_sfc)%id(ilocal(tamdar_sfc))        = platform%info%id
               iv%info(tamdar_sfc)%date_char(ilocal(tamdar_sfc)) = platform%info%date_char
               iv%info(tamdar_sfc)%levels(ilocal(tamdar_sfc))    = 1
               iv%info(tamdar_sfc)%lat(:,ilocal(tamdar_sfc))     = platform%info%lat
               iv%info(tamdar_sfc)%lon(:,ilocal(tamdar_sfc))     = platform%info%lon
               iv%info(tamdar_sfc)%elv(ilocal(tamdar_sfc))       = platform%info%elv
               iv%info(tamdar_sfc)%pstar(ilocal(tamdar_sfc))     = platform%info%pstar

               iv%info(tamdar_sfc)%slp(ilocal(tamdar_sfc))           = platform%loc%slp
               iv%info(tamdar_sfc)%pw(ilocal(tamdar_sfc))            = platform%loc%pw
               iv%info(tamdar_sfc)%x(:,ilocal(tamdar_sfc))           = platform%loc%x
               iv%info(tamdar_sfc)%y(:,ilocal(tamdar_sfc))           = platform%loc%y
               iv%info(tamdar_sfc)%i(:,ilocal(tamdar_sfc))           = platform%loc%i
               iv%info(tamdar_sfc)%j(:,ilocal(tamdar_sfc))           = platform%loc%j
               iv%info(tamdar_sfc)%dx(:,ilocal(tamdar_sfc))          = platform%loc%dx
               iv%info(tamdar_sfc)%dxm(:,ilocal(tamdar_sfc))         = platform%loc%dxm
               iv%info(tamdar_sfc)%dy(:,ilocal(tamdar_sfc))          = platform%loc%dy
               iv%info(tamdar_sfc)%dym(:,ilocal(tamdar_sfc))         = platform%loc%dym
               iv%info(tamdar_sfc)%proc_domain(:,ilocal(tamdar_sfc)) = platform%loc%proc_domain

               iv%info(tamdar_sfc)%obs_global_index(ilocal(tamdar_sfc)) = ntotal(tamdar_sfc)
            end if
         end if  ! for thin_conv_ascii skipping obs_index for which thin_3d is true like airep and tamdir


         if (global .and. n < 2) then
            if (test_transforms) exit dup_loop
            if (platform%loc % i >= ide) then
               platform%loc%i = platform%loc % i - ide
            else if (platform%loc % i < ids) then
               platform%loc%i = platform%loc % i + ide
            end if

            platform%loc%proc_domain = .not. platform%loc%proc_domain
         end if
      end do dup_loop
   end do reports

   close(iunit)

   call da_free_unit(iunit)

   ! thinning check
   if ( thin_conv_ascii ) then
      do n = 1, num_ob_indexes
          if (n==radar) cycle
          allocate ( in(thinning_grid_conv(n)%itxmax) )
          allocate (out(thinning_grid_conv(n)%itxmax) )
            do i = 1, thinning_grid_conv(n)%itxmax
               in(i) = thinning_grid_conv(n)%score_crit(i)
            end do
            ! Get minimum crit and associated processor index.
            call mpi_reduce(in, out, thinning_grid_conv(n)%itxmax, true_mpi_real, mpi_min, root, comm, ierr)
            call wrf_dm_bcast_real (out, thinning_grid_conv(n)%itxmax)
            do i = 1, thinning_grid_conv(n)%itxmax
              if( out(i) < 9.99e6) iv%info(n)%thin_ntotal=  iv%info(n)%thin_ntotal + 1
               if ( abs(out(i)-thinning_grid_conv(n)%score_crit(i)) > 1.0E-10 ) then
                  thinning_grid_conv(n)%ibest_obs(i) = 0
               end if
            end do
!            do j = iv%info(n)%plocal(iv%time -1)+1 , iv%info(n)%plocal(iv%time -1)+nlocal(n)
            do j = iv%info(n)%plocal(iv%time -1)+1 , nlocal(n)
               found = .false.
               do i = 1, thinning_grid_conv(n)%itxmax
                  if ( thinning_grid_conv(n)%ibest_obs(i) == j .and.         &
                       thinning_grid_conv(n)%score_crit(i) < 9.99e6 ) then
                   iv%info(n)%thin_nlocal =  iv%info(n)%thin_nlocal + 1
                     found = .true.

                     exit
                  end if
               end do
               if ( .not. found ) then
                  iv%info(n)%thinned(:,j) = .true.
               end if
            end do
         deallocate( in  )
         deallocate( out )
      end do ! loop over num_ob_indexes
   end if  ! thin_conv_ascii
   if (trace_use) call da_trace_exit("da_read_obs_ascii")

end subroutine da_read_obs_ascii


subroutine da_scan_obs_ascii (iv, filename, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Scan WRFVAR GTS observation file
   ! Updates:
   !       Date: 03/19/2009 -        Y.-R. Guo
   !           Added the time range check when reading in observations.
   !       Syed RH Rizvi NCAR/NESL/MMM/DAS Date:  02/21/2013 
   !          Updated with thinning option
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(*),      intent(in)    :: filename
   type(domain),     intent(in)     :: grid     ! first guess state.


   character (len =  10)   :: fmt_name
   character (len = 160)   :: info_string
!   character (len = 160)   :: fmt_info
!   character (len = 160)   :: fmt_srfc
!   character (len = 160)   :: fmt_each

   integer                 :: i, iost, fm, report, iunit
   type (multi_level_type) :: platform
   logical                 :: outside, outside_all
   real                    :: height_error
   integer                 :: ndup, n, obs_index

   real*8                :: obs_time, analysis_time
   integer               :: iyear, imonth, iday, ihour, imin
   real                  :: tdiff, dlat_earth,dlon_earth,crit
   integer               :: itt,itx,iout
   logical               :: iuse, thin_3d, is_surface
   integer               :: i1,j1,k, nlevels, levs
   real                  :: dx,dy,dxm,dym,zk
   real                  :: v_p(kms:kme),v_h(kms:kme)

   if (trace_use) call da_trace_entry("da_scan_obs_ascii")
! Initialize 
      if ( thin_conv_ascii ) then
          do n = 1, num_ob_indexes
             if ( n == radar ) cycle
             call cleangrids_conv(n)
          end do
      end if
   ! open file
   ! ---------
   call da_get_unit(iunit)
   open(unit   = iunit,     &
      FILE   = trim(filename), &
      FORM   = 'FORMATTED',  &
      ACCESS = 'SEQUENTIAL', &
      iostat =  iost,     &
      STATUS = 'OLD')

   if (iost /= 0) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iost," opening gts obs file "//trim(filename)
      call da_warning("da_scan_obs_ascii.inc",61,message(1:1))
      call da_free_unit(iunit)
      if (trace_use) call da_trace_exit("da_scan_obs_ascii")
      return
   end if

   ! read header

   head_info: do
      read (unit = iunit, fmt = '(A)', iostat = iost) info_string
      if (iost /= 0) then
         write(unit=message(1),fmt='(A,I3,A,I3)') &
            "Error",iost,"reading gts obs header on unit",iunit
         call da_warning("da_scan_obs_ascii.inc",74,message(1:1))
      if (trace_use) call da_trace_exit("da_scan_obs_ascii")
         return
      end if
      if (info_string(1:6) == 'EACH  ') exit
   end do head_info

   ! read formats

   read (iunit, fmt = '(A,1X,A)', iostat = iost) &
       fmt_name, fmt_info, &
       fmt_name, fmt_srfc,  &
       fmt_name, fmt_each

   if (iost /= 0) then
      write(unit=message(1),fmt='(A,I3,A,I3)') &
         "Error",iost,"reading gts obs formats on unit",iunit
         call da_warning("da_scan_obs_ascii.inc",91,message(1:1))
      if (trace_use) call da_trace_exit("da_scan_obs_ascii")
      return
   end if

   ! skip units
   read (iunit, fmt = '(A)') fmt_name

   ! loop over records

   report = 0 ! report number in file

   reports: do
      report = report+1

      ! read station general info

      read (iunit, fmt = fmt_info, iostat = iost) &
         platform%info%platform,    &
         platform%info%date_char,   &
         platform%info%name,        &
         platform%info%levels,      &
         platform%info%lat,         &
         platform%info%lon,         &
         platform%info%elv,         &
         platform%info%id
      if (iost /= 0) then
         ! FIX? This is expected, but its unclear how we handle failure
         ! here without assuming the fortran2003 convention on
         ! error statuses
         exit reports
      end if

      if (print_detail_obs) then
         write(unit=stdout, fmt = fmt_info) &
            platform%info%platform,    &
            platform%info%date_char,   &
            platform%info%name,        &
            platform%info%levels,      &
            platform%info%lat,         &
            platform%info%lon,         &
            platform%info%elv,         &
            platform%info%id
      end if

      if (platform%info%lon == 180.0) platform%info%lon =-180.000
      ! WHY?
      ! Fix funny wind direction at South Poles
      ! if (platform%info%lat < -89.9999 .or. platform%info%lat > 89.9999) then
      !    platform%info%lon = 0.0
      ! end if

      read (platform%info%platform(4:6), '(I3)') fm

      ! read model location
      read (iunit, fmt = fmt_srfc)  &
         platform%loc%slp%inv, platform%loc%slp%qc, platform%loc%slp%error, &
         platform%loc%pw%inv, platform%loc%pw%qc, platform%loc%pw%error

      ! read each level

      nlevels= platform%info%levels 

      do i = 1, platform%info%levels
         platform%each (i) = each_level_type(missing_r, missing, -1.0, & ! height
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! u
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! v
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! p
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! t
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! q
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! rh
            field_type(missing_r, missing, missing_r, missing, missing_r), & ! td
            field_type(missing_r, missing, missing_r, missing, missing_r))  ! speed 

         read (unit = iunit, fmt = trim (fmt_each)) &
            platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
            platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
            platform%each(i)%speed%error, &
            ! Here the 'direction' is stored in platform%each (i)%v:
            platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
            platform%each(i)%height,       &
            platform%each(i)%height_qc,    &
            height_error,                   &
            platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
            platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
            platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error
      end do

      ! Check if outside of the time range:

      read (platform%info%date_char,'(i4,4(1x,i2))') &
                                    iyear, imonth, iday, ihour, imin
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
      if ( obs_time < time_slots(0) .or. &
           obs_time >= time_slots(num_fgat_time) ) then
           cycle
      endif

      ! Restrict to a range of reports, useful for debugging
      if (report < report_start) cycle
      if (report > report_end) exit

      call da_llxy (platform%info, platform%loc, outside, outside_all)
      if( outside_all ) cycle reports

      if (platform%info%levels < 1) then
         if (fm /= 111 .and. fm /= 114) then
            cycle reports
         end if
      end if

      read (analysis_date,'(i4,4(1x,i2))') &
                                    iyear, imonth, iday, ihour, imin
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,analysis_time)
      tdiff = abs((obs_time - analysis_time)-0.02)
      dlat_earth = platform%info%lat
      dlon_earth = platform%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad


      ! Loop over duplicating obs for global
      ndup = 1
      if (global .and. (platform%loc%i < ids .or. platform%loc%i >= ide)) ndup= 2

      if (test_transforms) ndup = 1
      obs_index = fm_index(fm)
      do n = 1, ndup
         select case(fm)

         case (12) ; 
            if (.not.use_synopobs .or. iv%info(synop)%ntotal == max_synop_input) cycle reports
            if (n==1) iv%info(synop)%ntotal = iv%info(synop)%ntotal + 1
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(synop,dlat_earth,dlon_earth,crit,iv%info(synop)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
             else
                iv%info(synop)%nlocal = iv%info(synop)%nlocal + 1
             end if

         case (13, 17) ;    
            if (.not.use_shipsobs .or. iv%info(ships)%ntotal == max_ships_input) cycle reports
            if (n==1) iv%info(ships)%ntotal = iv%info(ships)%ntotal + 1
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(ships,dlat_earth,dlon_earth,crit,iv%info(ships)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
             else
                iv%info(ships)%nlocal = iv%info(ships)%nlocal + 1
             end if


         case (15:16) ;     
            if (.not.use_metarobs .or. iv%info(metar)%ntotal == max_metar_input) cycle reports
            if (n==1) iv%info(metar)%ntotal = iv%info(metar)%ntotal + 1
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(metar,dlat_earth,dlon_earth,crit,iv%info(metar)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
             else
                iv%info(metar)%nlocal = iv%info(metar)%nlocal + 1
             end if

         case (32:34) ;
            if (.not.use_pilotobs .or. iv%info(pilot)%ntotal == max_pilot_input) cycle reports
            if (n==1) iv%info(pilot)%ntotal = iv%info(pilot)%ntotal + 1
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(pilot,dlat_earth,dlon_earth,crit,iv%info(pilot)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
             else
                iv%info(pilot)%nlocal = iv%info(pilot)%nlocal + 1
             end if

         case (35:38) ;
            if (.not.use_soundobs .or. iv%info(sound)%ntotal == max_sound_input) cycle reports
            if (n==1) iv%info(sound)%ntotal     = iv%info(sound)%ntotal + 1
            if (n==1) iv%info(sonde_sfc)%ntotal = iv%info(sonde_sfc)%ntotal + 1
            if (outside) cycle reports
             if ( thin_conv_ascii ) then
                crit = tdiff
                call map2grids_conv(sound,dlat_earth,dlon_earth,crit,iv%info(sound)%nlocal,itx,1,itt,iout,iuse)
                crit = tdiff
                call map2grids_conv(sonde_sfc,dlat_earth,dlon_earth,crit,iv%info(sonde_sfc)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
             else
                iv%info(sound)%nlocal     = iv%info(sound)%nlocal + 1
                iv%info(sonde_sfc)%nlocal = iv%info(sonde_sfc)%nlocal + 1
             end if

         case (101) ;
            if (.not.use_tamdarobs .or. iv%info(tamdar)%ntotal == max_tamdar_input) cycle reports

! determine if level corresponds to surface 
            is_surface=.false.    
            levs = nlevels
            do i = 1, nlevels
               ! if (elevation and height are the same, it is surface.)
               if (platform%info%elv.ne.missing_r.and.platform%each(i)%height.ne.missing_r)then
                  if (abs(platform%info%elv - platform%each(i)%height) < 0.1) then
                     is_surface = .true.
                     levs = nlevels - 1
                    exit
                  end if
               end if
            end do

           if( is_surface) then

             if (n==1) iv%info(tamdar_sfc)%ntotal = iv%info(tamdar_sfc)%ntotal + 1
             if (outside) cycle reports
             if ( thin_conv_ascii ) then
                crit = tdiff
                call map2grids_conv(tamdar_sfc,dlat_earth,dlon_earth,crit,iv%info(tamdar_sfc)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
             else
                iv%info(tamdar_sfc)%nlocal = iv%info(tamdar_sfc)%nlocal + 1
             end if
           else ! not is_surface
             if ( levs > 0 .and. n==1) iv%info(tamdar)%ntotal = iv%info(tamdar)%ntotal + 1
             if (outside) cycle reports
             if( .not. thin_conv_ascii ) then
                iv%info(tamdar)%nlocal         = iv%info(tamdar)%nlocal + 1
             else  ! if thin_conv_ascii
                thin_3d=.true.
                i1   = platform%loc%i
                j1   = platform%loc%j
                dx   = platform%loc%dx
                dy   = platform%loc%dy
                dxm  = platform%loc%dxm
                dym  = platform%loc%dym
                do k=kms,kme
                   v_p(k) = dym*(dxm*grid%xb%p(i1,j1,k)+dx*grid%xb%p(i1+1,j1,k)) + &
                           dy*(dxm*grid%xb%p(i1,j1+1,k)+dx*grid%xb%p(i1+1,j1+1,k))
                end do
                do k=kms,kme
                   v_h(k) = dym*(dxm*grid%xb%h(i1,j1,k)+dx*grid%xb%h(i1+1,j1,k)) + &
                           dy*(dxm*grid%xb%h(i1,j1+1,k)+dx*grid%xb%h(i1+1,j1+1,k))
                end do
                do k=1,nlevels
                   zk = missing_r
                   if( platform%each(k)%p%qc  >= 0 ) then
                      call da_to_zk(platform%each(k)%p%inv, v_p, v_interp_p, zk)
                   else if( platform%each(k)%height_qc  >= 0 ) then
                      call da_to_zk(platform%each(k)%height, v_h, v_interp_h, zk)
                   else
                      write(unit=message(1),fmt='(A)')' For tamdar: neither height nor pressure qc is good'
                      call da_error("da_scan_obs_ascii.inc",345,message(1:1))
                   end if
                   if ( zk == missing_r ) cycle
                   crit = tdiff
                   call map2grids_conv(tamdar,dlat_earth,dlon_earth,crit,iv%info(tamdar)%nlocal,itx,1,itt,iout,iuse,zk,thin_3d)
                   if ( .not. iuse ) cycle
                end do ! loop over k levels
             end if ! if over thin_conv_ascii

           end if ! if is_surface 

         case (161) ;
            if (.not.use_mtgirsobs .or. iv%info(mtgirs)%ntotal == max_mtgirs_input) cycle reports
            if (n==1) iv%info(mtgirs)%ntotal     = iv%info(mtgirs)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(mtgirs,dlat_earth,dlon_earth,crit,iv%info(mtgirs)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
            else
                iv%info(mtgirs)%nlocal     = iv%info(mtgirs)%nlocal + 1
            end if

         case (86) ;
            if (.not.use_satemobs .or. iv%info(satem)%ntotal == max_satem_input) cycle reports
            ! Reject cloudy satem obs.
            if (platform%loc%pw%inv > 10.0) cycle reports
            if (n==1) iv%info(satem)%ntotal = iv%info(satem)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(satem,dlat_earth,dlon_earth,crit,iv%info(satem)%nlocal,itx,1,itt,iout,iuse)
                if ( .not. iuse ) cycle reports
            else
                iv%info(satem)%nlocal = iv%info(satem)%nlocal + 1
            end if

         case (88)    ;
            ! Geostationary or Polar orbitting Satellite AMVs:
            if (index(platform%info%name, 'MODIS') > 0 .or. &
                index(platform%info%name, 'modis') > 0 .or. &
                index(platform%info%id, 'AVHRR') > 0 )  then
               if (.not.use_polaramvobs .or. iv%info(polaramv)%ntotal == max_polaramv_input) cycle reports
               if (n==1) iv%info(polaramv)%ntotal = iv%info(polaramv)%ntotal + 1
               if (outside) cycle reports
               if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(polaramv,dlat_earth,dlon_earth,crit,iv%info(polaramv)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) cycle reports
               else
                   iv%info(polaramv)%nlocal = iv%info(polaramv)%nlocal + 1
               end if
	       obs_index = polaramv ! geoamv is the fm_index value for 88
            else
               if (.not.use_geoamvobs .or. iv%info(geoamv)%ntotal == max_geoamv_input) cycle reports
               if (n==1) iv%info(geoamv)%ntotal = iv%info(geoamv)%ntotal + 1
               if (outside) cycle reports
               if ( thin_conv_ascii ) then
                  crit = tdiff
                  call map2grids_conv(geoamv,dlat_earth,dlon_earth,crit,iv%info(geoamv)%nlocal,itx,1,itt,iout,iuse)
                   if ( .not. iuse ) cycle reports
               else
                   iv%info(geoamv)%nlocal = iv%info(geoamv)%nlocal + 1
               end if
            end if

         case (42,96:97) ;
            if (.not.use_airepobs .or. iv%info(airep)%ntotal == max_airep_input) cycle reports
            if (n==1) iv%info(airep)%ntotal = iv%info(airep)%ntotal + 1
            if (outside) cycle reports

            if( .not. thin_conv_ascii ) then
               iv%info(airep)%nlocal         = iv%info(airep)%nlocal + 1
            else  ! if thin_conv_ascii
               thin_3d=.true.
               i1   = platform%loc%i
               j1   = platform%loc%j
               dx   = platform%loc%dx
               dy   = platform%loc%dy
               dxm  = platform%loc%dxm
               dym  = platform%loc%dym
               do k=kms,kme
                  v_p(k) = dym*(dxm*grid%xb%p(i1,j1,k)+dx*grid%xb%p(i1+1,j1,k)) + &
                           dy*(dxm*grid%xb%p(i1,j1+1,k)+dx*grid%xb%p(i1+1,j1+1,k))
               end do
               do k=kms,kme
                  v_h(k) = dym*(dxm*grid%xb%h(i1,j1,k)+dx*grid%xb%h(i1+1,j1,k)) + &
                           dy*(dxm*grid%xb%h(i1,j1+1,k)+dx*grid%xb%h(i1+1,j1+1,k))
               end do
               do k=1,nlevels
                  zk = missing_r
                  if( platform%each(k)%p%qc  >= 0 ) then
                     call da_to_zk(platform%each(k)%p%inv, v_p, v_interp_p, zk)
                  else if( platform%each(k)%height_qc  >= 0 ) then
                     call da_to_zk(platform%each(k)%height, v_h, v_interp_h, zk)
                  else
                     write(unit=message(1),fmt='(A)')' For airep: neither height nor pressure qc is good'
                     call da_error("da_scan_obs_ascii.inc",442,message(1:1))
                  end if
                  if ( zk == missing_r ) cycle
                  crit = tdiff
                  call map2grids_conv(airep,dlat_earth,dlon_earth,crit,iv%info(airep)%nlocal,itx,1,itt,iout,iuse,zk,thin_3d)
                  if ( .not. iuse ) cycle
               end do ! loop over k levels

            end if ! if over thin_conv_ascii

         case (111, 114) ;       
            if ( (.not.use_gpspwobs  .and. fm == 111) .or. &
                  iv%info(gpspw)%ntotal == max_gpspw_input) cycle reports
            if ( (.not.use_gpsztdobs  .and. fm == 114) .or. &
                  iv%info(gpspw)%ntotal == max_gpspw_input) cycle reports
            if (n==1) iv%info(gpspw)%ntotal = iv%info(gpspw)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(gpspw,dlat_earth,dlon_earth,crit,iv%info(gpspw)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(gpspw)%nlocal = iv%info(gpspw)%nlocal + 1
            end if

         case (116) ;
            if (.not.use_gpsrefobs .or. iv%info(gpsref)%ntotal == max_gpsref_input) cycle reports
            if ( ob_format_gpsro /= ob_format_ascii ) cycle reports
            if (n==1) iv%info(gpsref)%ntotal = iv%info(gpsref)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(gpsref,dlat_earth,dlon_earth,crit,iv%info(gpsref)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(gpsref)%nlocal = iv%info(gpsref)%nlocal + 1
            end if

          case (121) ;
            ! SSM/T1 temperatures
            if (.not.use_ssmt1obs .or. iv%info(ssmt1)%ntotal == max_ssmt1_input) cycle reports
            if (n==1) iv%info(ssmt1)%ntotal = iv%info(ssmt1)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(ssmt1,dlat_earth,dlon_earth,crit,iv%info(ssmt1)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(ssmt1)%nlocal = iv%info(ssmt1)%nlocal + 1
            end if


         case (122) ;
            ! SSM/T2 relative humidities:
            if (.not.use_ssmt2obs .or. iv%info(ssmt2)%ntotal == max_ssmt2_input) cycle reports
            if (n==1) iv%info(ssmt2)%ntotal = iv%info(ssmt2)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(ssmt2,dlat_earth,dlon_earth,crit,iv%info(ssmt2)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(ssmt2)%nlocal = iv%info(ssmt2)%nlocal + 1
            end if

         case (281)    ;
            ! Scatterometer:
            if (.not.use_qscatobs .or. iv%info(qscat)%ntotal == max_qscat_input) cycle reports
            if (n==1) iv%info(qscat)%ntotal = iv%info(qscat)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(qscat,dlat_earth,dlon_earth,crit,iv%info(qscat)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(qscat)%nlocal = iv%info(qscat)%nlocal + 1
            end if
         case (132) ;
            if (.not.use_profilerobs .or. iv%info(profiler)%ntotal == max_profiler_input) cycle reports
            if (n==1) iv%info(profiler)%ntotal = iv%info(profiler)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(profiler,dlat_earth,dlon_earth,crit,iv%info(profiler)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(profiler)%nlocal = iv%info(profiler)%nlocal + 1
            end if

         case (135) ;
            if (.not.use_bogusobs .or. iv%info(bogus)%ntotal == max_bogus_input) cycle reports
            if (n==1) iv%info(bogus)%ntotal = iv%info(bogus)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(bogus,dlat_earth,dlon_earth,crit,iv%info(bogus)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
              iv%info(bogus)%nlocal = iv%info(bogus)%nlocal + 1
            end if

         case (18,19) ;             ! buoy
            if (.not.use_buoyobs .or. iv%info(buoy)%ntotal == max_buoy_input) cycle reports
            if (n==1) iv%info(buoy)%ntotal = iv%info(buoy)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(buoy,dlat_earth,dlon_earth,crit,iv%info(buoy)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(buoy)%nlocal = iv%info(buoy)%nlocal + 1
            end if

         case (133) ;               ! AIRS retrievals
            if (.not.use_airsretobs .or. iv%info(airsr)%ntotal == max_airsr_input) cycle reports
            if (n==1) iv%info(airsr)%ntotal = iv%info(airsr)%ntotal + 1
            if (outside) cycle reports
            if ( thin_conv_ascii ) then
               crit = tdiff
               call map2grids_conv(airsr,dlat_earth,dlon_earth,crit,iv%info(airsr)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) cycle reports
            else
               iv%info(airsr)%nlocal = iv%info(airsr)%nlocal + 1
            end if

         case default;
            write(unit=message(1), fmt='(a)') 'unsaved obs found:'
            write(unit=message(2), fmt='(2a)') &
               'platform%info%platform=', platform%info%platform
            write(unit=message(3), fmt='(a, i3)') &
                 'platform%info%levels=', platform%info%levels
            call da_warning("da_scan_obs_ascii.inc",573,message(1:3))
            cycle
         end select
         iv%info(obs_index)%max_lev = max(iv%info(obs_index)%max_lev, platform%info%levels)
      end do        !  loop over duplicate
   end do reports

   iv%info(sonde_sfc)%max_lev=1
   iv%info(tamdar_sfc)%max_lev=1
   iv%info(synop)%max_lev=1   ! To prevent some bad observations with more than 1 level, should I add more?
   iv%info(ships)%max_lev=1   
   iv%info(qscat)%max_lev=1   
   iv%info(metar)%max_lev=1   

   close(iunit)

   call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_scan_obs_ascii")

end subroutine da_scan_obs_ascii


subroutine da_read_obs_radar (iv, filename, grid)

   !-----------------------------------------------------------------------
   ! Purpose: Read the radar observation file
   !----------------------------------------------------------------------------------------!

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(len=*),  intent(in)    :: filename
   type(domain),     intent(in)     :: grid     ! first guess state.

   integer                       :: i, j, n, nn, iost, nlevels, fm
   integer                       :: total_radar
   integer                       :: iunit

   type (radar_multi_level_type) :: platform

   character (LEN = 120)         :: char_total_radar
   character (LEN = 120)         :: char_ned

   logical                       :: outside, outside_all
   integer                       :: n_dup, ndup

   real*8                        :: obs_time
   integer                       :: iyear, imonth, iday, ihour, imin

   integer                       :: ntotal,nlocal,ilocal
   integer                       :: radar_nlocal
   real, allocatable             :: in(:), out(:)


   if (trace_use) call da_trace_entry("da_read_obs_radar")

   ntotal = iv%info(radar)%ptotal(iv%time-1)
   nlocal = iv%info(radar)%plocal(iv%time-1)
   ilocal = nlocal 
   radar_nlocal = nlocal
 
   ! 1. open file

   call da_get_unit(iunit)
   open(unit   = iunit,     &
        FILE   = trim(filename), &
        FORM   = 'FORMATTED',  &
        ACCESS = 'SEQUENTIAL', &
        iostat =  iost,     &
        STATUS = 'OLD')

   if (iost /= 0) then
      ! Missing file does not matter
      call da_warning("da_read_obs_radar.inc",52, &
         (/"Cannot open radar file "//filename/))
      call da_free_unit(iunit) 
      if (trace_use) call da_trace_exit("da_read_obs_radar")
      return
   end if

   ! 2. read total radar

   !  2.1 read first line

   read (unit=iunit, fmt = '(A)', iostat = iost) char_total_radar

   !  2.2 process error

   if (iost /= 0) then
     call da_error("da_read_obs_radar.inc",68, &
        (/"Cannot read radar file"/))
   end if

   !  2.3 total radar number

   read (unit=char_total_radar (15:17),fmt='(I3)', iostat = iost) total_radar

   if (print_detail_radar) write (unit=stdout,fmt='(/,A,I3,/)') &
       ' TOTAL RADAR: ', total_radar

   !  2.4 skip one line

   read (unit=iunit, fmt = '(A)', iostat = iost)

   ! 3. read radar data

   do nn = 1, total_radar
 
      ! 3.1 skip one blank line

      read (unit=iunit, fmt = '(A)', iostat = iost)

      ! 3.2 read header

      read (unit=iunit, fmt = '(A)', iostat = iost) char_ned

      ! 3.3 read header information

      read (unit=char_ned (1:5),fmt='(A5)', iostat = iost) platform % stn % platform

      if (print_detail_radar) write (unit=stdout,fmt='(A)') 'RADAR Observations information'

      read (unit=char_ned (8:19),fmt='(A12)', iostat = iost) platform % stn % name

!      if (print_detail_radar) write (unit=stdout,fmt='(A,A5,A,A12)')  &
       write (unit=stdout,fmt='(A,A5,A,A12)')  &
                           ' Reading ',platform % stn % platform, &
                           ' data at station:', platform % stn % name

      read (unit=char_ned(20:27),fmt='(F8.3)', iostat = iost) platform % stn % lon

      read (unit=char_ned (30:37),fmt='(F8.3)', iostat = iost) platform % stn % lat

      read (unit=char_ned (40:47),fmt='(F8.1)', iostat = iost) platform % stn % elv

      if (print_detail_radar) write (unit=stdout,fmt='(A,2(F8.3,2X),F8.1)')  &
         'The station longitude, latitude, and altitude are: ', &
         platform % stn % lon, &
         platform % stn % lat, platform % stn % elv

      read (unit=char_ned (50:68),fmt='(A19)', iostat = iost) platform % stn % date_char

      if (print_detail_radar) write (unit=stdout,fmt='(A,A19)')   &
         'The observation time for this data is ',     &
         platform % stn % date_char

      read (unit=char_ned (69:74),fmt='(I6)', iostat = iost) platform % stn % numobs

      if (print_detail_radar) write (unit=stdout,fmt='(A,I6)')   &
         'Total number of Super-observations is ', &
         platform % stn % numobs


      read (unit=char_ned (75:80),fmt='(I6)', iostat = iost) platform % stn % levels

      if (print_detail_radar) write (unit=stdout,fmt='(A,I6)')   &
         'Vertical layers for each Super-observation is ', &
         platform % stn % levels

      ! 3.4 skip two lines

      read (unit=iunit, fmt = '(A)', iostat = iost)
      read (unit=iunit, fmt = '(A)', iostat = iost)

      ! 3.5 loop over records

      reports: do j = 1, platform % stn % numobs

         ! 3.5.1 read station general info

         read (unit = iunit, iostat = iost, &
                      fmt = '(A12,3X,A19,2X,2(F12.3,2X),F8.1,2X,I6)') &
                      platform % info % platform,  &
                      platform % info % date_char, &
                      platform % info % lat,       &
                      platform % info % lon,       &
                      platform % info % elv,       &
                      platform % info % levels
      if (platform%info%lon == 180.0  ) platform%info%lon =-180.000
      ! Fix funny wind direction at Poles
      if (platform%info%lat < -89.9999 .or. platform%info%lat > 89.9999) then
         platform%info%lon = 0.0
      end if

         read(platform % info % platform (4:6), '(I3)') fm

         ! 3.5.2 read each level

         do i = 1, platform % info % levels
            ! height
            platform%each(i) = radar_each_level_type(missing_r, missing, -1.0,       &
               field_type(missing_r, missing, missing_r, missing, missing_r), & ! rv
               field_type(missing_r, missing, missing_r, missing, missing_r))   ! rf

            read (unit = iunit, fmt = '(3X, F12.1, 2(F12.3,I4,F12.3,2X))') &
                             platform % each (i) % height,           &
                             platform % each (i) % rv % inv,         &
                             platform % each (i) % rv % qc,          &
                             platform % each (i) % rv % error,       &
                             platform % each (i) % rf % inv,         &
                             platform % each (i) % rf % qc,          &
                             platform % each (i) % rf % error

            if (platform % each (i) % rv % error == 0.0) then
                 platform % each (i) % rv % error  = 1.0
            end if

            if (platform % each (i) % rf % error == 0.0) then
                 platform % each (i) % rf % error  = 1.0
            end if

            if (platform % each (i) % rv % inv   == missing_r .or. &
                platform % each (i) % rv % error == missing_r) then
                platform % each (i) % rv % qc     = missing_data
            end if

            if (platform % each (i) % rf % inv   == missing_r .or. &
                platform % each (i) % rf % error == missing_r) then
                platform % each (i) % rf % qc     = missing_data
            end if
         end do

         ! Check if outside of the time range:

         read (platform%info%date_char,'(i4,4(1x,i2))') &
               iyear, imonth, iday, ihour, imin
         call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
         ! If you skip this part, you have to skip in da_scan_obs_radar.inc too!
         if ( obs_time < time_slots(0) .or. &
              obs_time >= time_slots(num_fgat_time) ) then
            if (print_detail_radar) then
               write(unit=stdout, fmt='(a)') '*** Outside of the time range:'
               write(unit=stdout, fmt=fmt_info) &
                     platform%info%platform,    &
                     platform%info%date_char,   &
                     platform%stn%name
            end if
            cycle reports
         endif

         call da_llxy (platform%info, platform%loc, outside, outside_all)
         if( outside_all ) then
            if (print_detail_radar) then
               write(unit=stdout, fmt='(a)') '*** Report is outside of domain:'
               write(unit=stdout, fmt='(2x,a,2(2x,f7.3),2x,a)') &
                     platform%info%platform,    &
                     platform%info%lat,   &
                     platform%info%lon,   &
                     platform%stn%name
            end if
            cycle reports
         end if

        read (analysis_date,'(i4,4(1x,i2))') &
                                    iyear, imonth, iday, ihour, imin

         nlevels = platform%info%levels

            iv%info(radar)%max_lev = max(iv%info(radar)%max_lev, platform%info%levels)

         if (nlevels > max_ob_levels) then
            write(unit=message(1),fmt='(A,2I8)') &
               ' radar=> nlevels > max_ob_levels:',nlevels, max_ob_levels
            call da_warning("da_read_obs_radar.inc",242,message(1:1)) 
            nlevels = max_ob_levels
             platform%info%levels = nlevels
         else if (nlevels < 1) then
            cycle reports
         end if

         ! Loop over duplicating obs for global
         n_dup = 1
         if (global .and. &
            (platform%loc%i == ids .or. platform%loc%i == ide)) n_dup= 2
         do ndup = 1, n_dup
            select case (fm)
            case (128)
            if (.not.use_radarobs .or. ntotal == max_radar_input) cycle reports

               if (ndup==1 ) ntotal = ntotal + 1
               if (outside) cycle reports
               nlocal = nlocal + 1
               ilocal = ilocal + 1
               iv % radar (ilocal) % stn_loc % lon = platform % stn % lon
               iv % radar (ilocal) % stn_loc % lat = platform % stn % lat
               iv % radar (ilocal) % stn_loc % elv = platform % stn % elv

               iv%info(radar)%levels(ilocal)    = nlevels
               iv%info(radar)%name(ilocal)      = platform%info%name
               iv%info(radar)%platform(ilocal)  = platform%info%platform
               iv%info(radar)%id(ilocal)        = platform%info%id
               iv%info(radar)%date_char(ilocal) = platform%info%date_char
               iv%info(radar)%lat(:,ilocal)     = platform%info%lat
               iv%info(radar)%lon(:,ilocal)     = platform%info%lon
               iv%info(radar)%elv(ilocal)       = platform%info%elv
               iv%info(radar)%pstar(ilocal)     = platform%info%pstar
               iv%info(radar)%slp(ilocal)           = platform%loc%slp
               iv%info(radar)%pw(ilocal)            = platform%loc%pw
               iv%info(radar)%x(:,ilocal)           = platform%loc%x
               iv%info(radar)%y(:,ilocal)           = platform%loc%y
               iv%info(radar)%i(:,ilocal)           = platform%loc%i
               iv%info(radar)%j(:,ilocal)           = platform%loc%j
               iv%info(radar)%dx(:,ilocal)          = platform%loc%dx
               iv%info(radar)%dxm(:,ilocal)         = platform%loc%dxm
               iv%info(radar)%dy(:,ilocal)          = platform%loc%dy
               iv%info(radar)%dym(:,ilocal)         = platform%loc%dym
               iv%info(radar)%proc_domain(:,ilocal) = platform%loc%proc_domain
               iv%info(radar)%obs_global_index(ilocal) = ntotal
               if( nlocal == ilocal) then
                  allocate (iv % radar (ilocal) % model_p  (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % model_rho(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % model_qrn(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % model_qcl(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % model_qci(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % model_qsn(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % model_qgr(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % height   (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % height_qc(1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rv       (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rf       (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rrn      (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rcl      (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rci      (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rsn      (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rgr      (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rqv      (1:iv%info(radar)%max_lev))

                  allocate (iv % radar (ilocal) % rrno     (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rclo     (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rcio     (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rsno     (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rgro     (1:iv%info(radar)%max_lev))
                  allocate (iv % radar (ilocal) % rqvo     (1:iv%info(radar)%max_lev))
               end if
               do i = 1, nlevels
                  iv % radar (ilocal) % height(i)    = platform % each(i) % height
                  iv % radar (ilocal) % height_qc(i) = platform % each(i) % height_qc
                  iv % radar (ilocal) % rv(i)        = platform % each(i) % rv
                  iv % radar (ilocal) % rf(i)        = platform % each(i) % rf

                  iv % radar (ilocal) % rrn(i) % inv   = missing_r
                  iv % radar (ilocal) % rrn(i) % qc    = missing_data
                  iv % radar (ilocal) % rrn(i) % error = missing_r
                  iv % radar (ilocal) % rrno(i)        = missing_r

                  iv % radar (ilocal) % rcl(i) % inv   = missing_r
                  iv % radar (ilocal) % rcl(i) % qc    = missing_data
                  iv % radar (ilocal) % rcl(i) % error = missing_r
                  iv % radar (ilocal) % rclo(i)        = missing_r

                  iv % radar (ilocal) % rci(i) % inv   = missing_r
                  iv % radar (ilocal) % rci(i) % qc    = missing_data
                  iv % radar (ilocal) % rci(i) % error = missing_r
                  iv % radar (ilocal) % rcio(i)        = missing_r

                  iv % radar (ilocal) % rsn(i) % inv   = missing_r
                  iv % radar (ilocal) % rsn(i) % qc    = missing_data
                  iv % radar (ilocal) % rsn(i) % error = missing_r
                  iv % radar (ilocal) % rsno(i)        = missing_r

                  iv % radar (ilocal) % rgr(i) % inv   = missing_r
                  iv % radar (ilocal) % rgr(i) % qc    = missing_data
                  iv % radar (ilocal) % rgr(i) % error = missing_r
                  iv % radar (ilocal) % rgro(i)        = missing_r

                  iv % radar (ilocal) % rqv(i) % inv   = missing_r
                  iv % radar (ilocal) % rqv(i) % qc    = missing_data
                  iv % radar (ilocal) % rqv(i) % error = missing_r
                  iv % radar (ilocal) % rqvo(i)        = missing_r
               end do

            case default;
               write(unit=message(1), fmt='(a)') 'Unsaved obs found:'
               write(unit=message(2), fmt='(2a)') &
                  'platform % info % platform=', platform % info % platform
               write(unit=message(3), fmt='(a, i3)') &
                  'platform % info % levels=', platform % info % levels
               call da_warning("da_read_obs_radar.inc",356,message(1:3))
            end select
            if (global .and. ndup == 1) then
               if (platform%loc % i >= ide) then
                  platform%loc%i = ids
                  platform%loc%proc_domain = .false.
               else if (platform%loc % i <= ids) then
                  platform%loc%i = ide
                  platform%loc%proc_domain = .false.
               end if
            end if
         end do        !  loop over duplicate
      end do reports

                radar_nlocal = nlocal             


   end do  ! total_radar

   if (print_detail_radar) write (unit=stdout,fmt='(/,A,I3,/)') &
       ' Processed TOTAL RADAR: ', total_radar
   close(iunit)
   call da_free_unit(iunit)
   if (trace_use) call da_trace_exit("da_read_obs_radar")
end subroutine da_read_obs_radar

subroutine da_scan_obs_radar (iv, filename, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Scan the radar observation file
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(len=*),  intent(in)    :: filename
   type(domain),     intent(in)     :: grid     ! first guess state.

   integer                       :: i, j, n, iost, nlevels, fm
   integer                       :: total_radar
   integer                       :: iunit

   type (radar_multi_level_type) :: platform

   character (LEN = 120)         :: char_total_radar
   character (LEN = 120)         :: char_ned

   logical                       :: outside, outside_all
   integer                       :: n_dup, ndup

   real*8                        :: obs_time
   integer                       :: iyear, imonth, iday, ihour, imin


   if (trace_use) call da_trace_entry("da_scan_obs_radar")

   ! 1. open file
   ! ============

   call da_get_unit(iunit)
   open(unit   = iunit,     &
        FILE   = trim(filename), &
        FORM   = 'FORMATTED',  &
        ACCESS = 'SEQUENTIAL', &
        iostat =  iost,     &
        STATUS = 'OLD')

   if (iost /= 0) then
      ! Does not matter of radar file missing
      call da_warning("da_scan_obs_radar.inc",44, &
         (/"Cannot open radar file "//filename/))
      call da_free_unit(iunit) 
      if (trace_use) call da_trace_exit("da_scan_obs_radar")
      return
   end if
   ! 1.1  Initialize
   ! ============


   ! 2. read total radar
   ! ===================

   ! 2.1 read first line
   !     ---------------

   read (unit=iunit, fmt = '(A)', iostat = iost) char_total_radar
   if (iost /= 0) then
      ! Does matter if present and unreadable
      call da_error("da_scan_obs_radar.inc",63, &
         (/"Cannot read radar file"/))
   end if

   ! 2.3 total radar number

   read (unit=char_total_radar (15:17),fmt='(I3)', iostat = iost) total_radar

   ! 2.4 skip one lines

   read (unit=iunit, fmt = '(A)', iostat = iost)

   ! 3. read radar data

   do n = 1, total_radar

      ! 3.1 skip one blank line

      read (unit=iunit, fmt = '(A)', iostat = iost)

      ! 3.2 read header

      read (unit=iunit, fmt = '(A)', iostat = iost) char_ned

      ! 3.3 read header information

      read (unit=char_ned (69:74), fmt='(I6)', iostat = iost) platform % stn % numobs

      ! 3.4 skip two lines

      read (unit=iunit, fmt = '(A)', iostat = iost)
      read (unit=iunit, fmt = '(A)', iostat = iost)

      ! 3.5 loop over records

      reports: do j = 1, platform % stn % numobs

         ! 3.5.1 read station general info

         read (unit = iunit, iostat = iost, &
                      fmt = '(A12,3X,A19,2X,2(F12.3,2X),F8.1,2X,I6)') &
                      platform % info % platform,  &
                      platform % info % date_char, &
                      platform % info % lat,       &
                      platform % info % lon,       &
                      platform % info % elv,       &
                      platform % info % levels

         if (platform%info%lon == 180.0  ) platform%info%lon =-180.000
         ! Fix funny wind direction at Poles
         if (platform%info%lat < -89.9999 .or. platform%info%lat > 89.9999) then
            platform%info%lon = 0.0
         end if

         read(unit=platform % info % platform (4:6), fmt='(I3)') fm

         !     3.5.2 read each level

         do i = 1, platform % info % levels
            ! height
            platform%each (i) = radar_each_level_type(missing_r, missing, -1.0,&
               field_type(missing_r, missing, missing_r, missing, missing_r), & ! rv
               field_type(missing_r, missing, missing_r, missing, missing_r))   ! rf

            read (unit = iunit, fmt = '(3X, F12.1, 2(F12.3,I4,F12.3,2X))') &
                             platform % each (i) % height,           &
                             platform % each (i) % rv % inv,         &
                             platform % each (i) % rv % qc,          &
                             platform % each (i) % rv % error,       &
                             platform % each (i) % rf % inv,         &
                             platform % each (i) % rf % qc,          &
                             platform % each (i) % rf % error

            if (platform % each (i) % rv % error == 0.0) then
                 platform % each (i) % rv % error  = 1.0
            end if

            if (platform % each (i) % rf % error == 0.0) then
                 platform % each (i) % rf % error  = 1.0
            end if

            if (platform % each (i) % rv % inv   == missing_r .or. &
                platform % each (i) % rv % error == missing_r) then
                platform % each (i) % rv % qc     = missing_data
            end if

            if (platform % each (i) % rf % inv   == missing_r .or. &
                platform % each (i) % rf % error == missing_r) then
                platform % each (i) % rf % qc     = missing_data
            end if

         end do

         ! Check if outside of the time range:

         read (platform%info%date_char,'(i4,4(1x,i2))') &
               iyear, imonth, iday, ihour, imin
         call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
         if ( obs_time < time_slots(0) .or. &
              obs_time >= time_slots(num_fgat_time) ) then
            cycle reports
         endif

         call da_llxy (platform%info, platform%loc, outside, outside_all)
         if( outside_all ) cycle reports

         nlevels = platform%info%levels

         if (nlevels > max_ob_levels) then
             write(unit=message(1),fmt='(A,2I8)') &
                ' radar=> nlevels > max_ob_levels:',nlevels, max_ob_levels
             call da_warning("da_scan_obs_radar.inc",174,message(1:1))

             nlevels = max_ob_levels
             platform%info%levels = nlevels
         else if (nlevels < 1) then
            cycle reports
         end if


         ! Loop over duplicating obs for global
         n_dup = 1
         if (global .and. (platform%loc%i == ids .or. platform%loc%i == ide)) n_dup= 2
   
         do ndup = 1, n_dup
            select case (fm)

            case (128)
               if (.not.use_radarobs .or. iv%info(radar)%ntotal == max_radar_input) cycle reports
               if (ndup==1) iv%info(radar)%ntotal = iv%info(radar)%ntotal + 1
               if (outside) cycle reports
               iv%info(radar)%nlocal = iv%info(radar)%nlocal + 1

            case default;
               write(unit=stdout, fmt='(a)') 'Warning: unsaved obs found:'

               write(unit=stdout, fmt='(2a)') &
                  'platform % info % platform=', platform % info % platform

               write(unit=stdout, fmt='(a, i3)') &
                  'platform % info % levels=', platform % info % levels
            end select

            iv%info(radar)%max_lev = max(iv%info(radar)%max_lev, platform%info%levels)
         end do        !  loop over duplicate
      end do reports
      
   end do ! total_radar

   close (iunit)
   call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_scan_obs_radar")

end subroutine da_scan_obs_radar
subroutine da_scan_obs_rain (iv, filename, ifgat)

   !---------------------------------------------------------------------------
   ! Purpose: Scan the rain observation file
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(len=*),  intent(in)    :: filename
   integer,           intent(in)    :: ifgat
   integer                       :: i, j, n, iost, nlevels, fm
   integer                       :: file_rain
   integer                       :: iunit

   type (rain_single_level_type) :: platform

   character (len = 120)         :: char_file_rain
   character (len = 120)         :: fmt_name
   character (len = 160)         :: info_string
   logical                       :: outside, outside_all
   integer                       :: n_dup, ndup

   real*8                        :: obs_time
   integer                       :: iyear, imonth, iday, ihour, imin

   ! for thinning
   real                          :: dlat_earth,dlon_earth,crit
   integer                       :: itt,itx,iout
   logical                       :: iuse
   integer                       :: tp, nlocal
   real, allocatable             :: in(:), out(:)

   integer                      :: num_outside_all, num_outside_time, num_thinned, num_report

   if (trace_use) call da_trace_entry("da_scan_obs_rain")

   num_report       = 0
   num_outside_all  = 0
   num_outside_time = 0
   num_thinned      = 0
   tp               = iv%info(rain)%plocal(ifgat-1)
   nlocal           = 0

   ! 1. open file
   ! ============
   call da_get_unit(iunit)
   open(unit   = iunit,     &
        FILE   = trim(filename), &
        FORM   = 'FORMATTED',  &
        ACCESS = 'SEQUENTIAL', &
        iostat =  iost,     &
        STATUS = 'OLD')

   if (iost /= 0) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iost," opening rainfall obs file "//trim(filename)
      call da_warning("da_scan_obs_rain.inc",58,message(1:1))
      call da_free_unit(iunit) 
      if (trace_use) call da_trace_exit("da_scan_obs_rain")
      return
   end if

   ! 2. read rainfall data
   ! ===================

   ! 2.1 read first line
   !     ---------------

   read (unit=iunit, fmt = '(A)', iostat = iost) char_file_rain

   ! 2.2 process error

   if (iost /= 0) then
      ! Does matter if present and unreadable
      call da_error("da_scan_obs_rain.inc",76, &
         (/"Cannot read rainfall file"/))
   end if

   ! 2.3read header info

   head_info: do

      read (unit=iunit, fmt = '(A)', iostat = iost) info_string

      if (iost /= 0) then
         write(unit=message(1),fmt='(A,I3,A,I3)') &
            "Error",iost,"reading rainfall obs header on unit",iunit
         call da_warning("da_scan_obs_rain.inc",89,message(1:1))
      if (trace_use) call da_trace_exit("da_scan_obs_rain")
         return
      end if

      if (info_string(1:6) == 'EACH  ') exit

   end do head_info

   ! 2.3 total rainfall data info

   read (unit=char_file_rain (8:14),fmt='(I7)', iostat = iost) file_rain

   ! 2.4 skip one lines

   read (unit=iunit, fmt = '(A)', iostat = iost)

   ! 3. read rain data

   reports:   do n = 1, file_rain

      ! 3.1 read station general info

      read (unit = iunit, iostat = iost, &
         fmt = '(A12,1X,A19,1X,I6,2(F12.3,2X),F8.1,1X,A5)') &
         platform % info % platform,  &
         platform % info % date_char, &
         platform % info % levels,    &
         platform % info % lat,       &
         platform % info % lon,       &
         platform % info % elv,       &
         platform % info % id 

      if (print_detail_rain) then
         write(unit=stdout, fmt = '(A12,1X,A19,1X,I6,2(F12.3,2X),F8.1,1X,A5)') &
            platform % info % platform,  &
            platform % info % date_char, &
            platform % info % levels,    &
            platform % info % lat,       &
            platform % info % lon,       &
            platform % info % elv,       &
            platform % info % id
      end if         

      read(unit=platform % info % platform (4:6), fmt='(I3)') fm

      num_report = num_report+1

      ! 3.2 read rainfall data 

      platform%each (1) = rain_each_type(missing_r, missing, -1.0,&
         field_type(missing_r, missing, missing_r, missing, missing_r))

      read (unit = iunit, fmt = '(F12.3,F12.3,I4,F12.3)') &
         platform % each (1) % height,             &
         platform % each (1) % rain % inv,         &
         platform % each (1) % rain % qc,          &
         platform % each (1) % rain % error       

     ! 3.3 Check if outside of the time range:

      read (platform%info%date_char,'(i4,4(1x,i2))') &
            iyear, imonth, iday, ihour, imin
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
      if ( obs_time < time_slots(0) .or. &
           obs_time >= time_slots(num_fgat_time) ) then
           num_outside_time = num_outside_time + 1
         cycle reports
      endif

      call da_llxy (platform%info, platform%loc, outside, outside_all)

      nlevels = platform%info%levels

      if (outside_all) then
         num_outside_all = num_outside_all + 1
         cycle reports
      end if

      dlat_earth = platform%info%lat
      dlon_earth = platform%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad

      ! 3.4  Loop over duplicating obs for global
      n_dup = 1
      if (global .and. &
         (platform%loc%i == ids .or. platform%loc%i == ide)) n_dup= 2
  
      if (test_transforms) ndup = 1
   
      do n_dup = 1, n_dup
         select case (fm)

         case (129)
            if (iv%info(rain)%nlocal > max_rain_input) then
               write(unit=message(1),fmt='(A,I6,A,I6)') &
                  ' rain #= ',iv%info(rain)%nlocal, ' > max_rain_input = ', max_rain_input
               call da_error("da_scan_obs_rain.inc",189,message(1:1))
            end if
            iv%info(rain)%ntotal = iv%info(rain)%ntotal + 1
            if (outside) cycle reports
            if ( thin_rainobs ) then
               crit = 1.0
               call map2grids(rain,ifgat,dlat_earth,dlon_earth,crit,iv%info(rain)%nlocal,itx,1,itt,iout,iuse)
               if ( .not. iuse ) then
                  num_thinned = num_thinned + 1
                  cycle reports
               end if
            else
               iv%info(rain)%nlocal = iv%info(rain)%nlocal + 1
            endif

         case default;
            write(unit=message(1), fmt='(a)') 'unsaved obs found:'
            write(unit=message(2), fmt='(2a)') &
               'platform%info%platform=', platform%info%platform
            write(unit=message(3), fmt='(a, i3)') &
                 'platform%info%levels=', platform%info%levels
            call da_warning("da_scan_obs_rain.inc",210,message(1:3))
            cycle reports
         end select

         iv%info(rain)%max_lev = 1 

      end do        !  loop over duplicate
   end do reports

   ! thinning check
   if ( thin_rainobs ) then
     if ( iv%info(rain)%ntotal > 0 ) then

         ! Get minimum crit and associated processor index.
         allocate ( in  (thinning_grid(rain,ifgat)%itxmax) )
         allocate ( out (thinning_grid(rain,ifgat)%itxmax) )
         do i = 1, thinning_grid(rain,ifgat)%itxmax
            in(i) = thinning_grid(rain,ifgat)%score_crit(i)
         end do

         call mpi_reduce(in, out, thinning_grid(rain,ifgat)%itxmax, true_mpi_real, mpi_min, root, comm, ierr)
         call wrf_dm_bcast_real (out, thinning_grid(rain,ifgat)%itxmax)

         do i = 1, thinning_grid(rain,ifgat)%itxmax
            if ( abs(out(i)-thinning_grid(rain,ifgat)%score_crit(i)) > 1.0E-10 ) then
               thinning_grid(rain,ifgat)%ibest_obs(i) = 0
            end if
         end do

         thinning_grid(rain,ifgat)%score_crit(:) = out(:)
               
         deallocate( in  )
         deallocate( out )
            
         do j = (1+tp), iv%info(rain)%nlocal
            do i = 1, thinning_grid(rain,ifgat)%itxmax
               if ( thinning_grid(rain,ifgat)%ibest_obs(i) == j .and.         &
                    thinning_grid(rain,ifgat)%score_crit(i) < 9.99e6 ) then
                  nlocal = nlocal + 1
                  exit
               end if
            end do
         end do

      num_thinned = num_thinned + iv%info(rain)%nlocal - nlocal
      iv%info(rain)%nlocal = tp + nlocal
      end if 
   end if  ! thin_rainobs

   write(unit=message(1),fmt='(A,4(1x,i7))') &
      'da_scan_obs_rain: num_report, num_outside_all, num_outside_time, num_thinned: ', &
      num_report, num_outside_all, num_outside_time, num_thinned
   call da_message(message(1:1))

   close (iunit)
   call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_scan_obs_rain")


end subroutine da_scan_obs_rain


subroutine da_read_obs_rain (iv, filename, ifgat)

   !-----------------------------------------------------------------------
   ! Purpose: Read the rain observation file
   !----------------------------------------------------------------------------------------!

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(len=*),  intent(in)    :: filename
   integer,           intent(in)    :: ifgat

   character (len = 120)         :: char_total_rain
   character (len = 160)         :: info_string
   
   integer                       :: i, j, n, iost, nlevels, fm

   type (rain_single_level_type) :: platform

   logical                       :: outside, outside_all

   integer                       :: total_rain
   integer                       :: n_dup, ndup, iunit
   integer                       :: nlocal
   integer                       :: ilocal
   integer                       :: ntotal

   real*8                        :: obs_time
   integer                       :: iyear, imonth, iday, ihour, imin

   ! for thinning
   real                          :: dlat_earth,dlon_earth,crit,dist
   integer                       :: itt,itx,iout
   logical                       :: iuse
   integer                       :: tp

   integer                       :: num_outside_all, num_outside_time, num_thinned, num_report

   if (trace_use) call da_trace_entry("da_read_obs_rain")

   nlocal = iv%info(rain)%plocal(ifgat-1)
   ntotal = iv%info(rain)%ptotal(ifgat-1)

   num_report       = 0
   num_outside_all  = 0
   num_outside_time = 0
   num_thinned      = 0
   ilocal           = 0
   tp               = nlocal

   ! 1. open file
   ! ============
   call da_get_unit(iunit)

   open(unit   = iunit,     &
        FILE   = trim(filename), &
        FORM   = 'FORMATTED',  &
        ACCESS = 'SEQUENTIAL', &
        iostat =  iost,     &
        STATUS = 'OLD')

   if (iost /= 0) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iost," opening rainfall obs file "//trim(filename)
      call da_warning("da_read_obs_rain.inc",65,message(1:1))
      call da_free_unit(iunit)
      if (trace_use) call da_trace_exit("da_read_obs_rain")
      return
   end if

   ! 2. read rainfall info 
   ! ==================

   !  2.1 read first line
   !      ---------------

   read (unit=iunit, fmt = '(A)', iostat = iost) char_total_rain

   !  2.2 process error

   if (iost /= 0) then
     call da_error("da_read_obs_rain.inc",82, &
        (/"Cannot read rainfall file"/))
   end if

   !  2.3 read header info

   head_info: do
      read (unit=iunit, fmt = '(A)', iostat = iost) info_string
      if (iost /= 0) then
         write(unit=message(1),fmt='(A,I3,A,I3)') &
            "Error",iost,"reading rainfall obs header on unit",iunit
         call da_warning("da_read_obs_rain.inc",93,message(1:1))
      if (trace_use) call da_trace_exit("da_scan_obs_rain")
      return
      end if
      if (info_string(1:6) == 'EACH  ') exit
   end do head_info

   !  2.4 total rainfall data info

   read (unit=char_total_rain (8:14),fmt='(I7)', iostat = iost) total_rain

   if (print_detail_rain) write (unit=stdout,fmt='(/,A,I7,/)') &
       ' TOTAL RAINFALL DATA: ', total_rain

   !  2.5 skip one line

   read (unit=iunit, fmt = '(A)', iostat = iost)

   ! 3. read rainfall data
   ! =================  

   reports:   do n = 1, total_rain

      if (print_detail_rain) write (unit=stdout,fmt='(A)') 'RAIN Observations information'

      ! 3.1 read station general info

      read (unit = iunit, iostat = iost, &
                   fmt = '(A12,1X,A19,1X,I6,2(F12.3,2X),F8.1,1X,A5)') &
                   platform % info % platform,  &
                   platform % info % date_char, &
                   platform % info % levels,    &
                   platform % info % lat,       &
                   platform % info % lon,       &
                   platform % info % elv,       &
                   platform % info % id 

      if (print_detail_rain) then
         write(unit = stdout, fmt ='(A12,1X,A19,1X,I6,2(F12.3,2X),F8.1,1X,A5)') &
            platform%info%platform,    &
            platform%info%date_char,   &
            platform%info%levels,      &
            platform%info%lat,         &
            platform%info%lon,         &
            platform%info%elv,         &
            platform%info%id
      end if

      read(platform % info % platform (4:6), '(I3)') fm

      ! 3.2 read rainfall data 

      platform%each(1) = rain_each_type(missing_r, missing, -1.0,       &
         field_type(missing_r, missing, missing_r, missing, missing_r))

      read (unit = iunit, fmt = '(F12.3,F12.3,I4,F12.3)') &
         platform % each (1) % height,             &
         platform % each (1) % rain % inv,         &
         platform % each (1) % rain % qc,          &
         platform % each (1) % rain % error

      if (platform % each (1) % rain % error == 0.0) then
         platform % each (1) % rain % error  = 1.0
      end if

      if (platform % each (1) % rain % inv   == missing_r .or. &
          platform % each (1) % rain % error == missing_r) then
          platform % each (1) % rain % qc     = missing_data
      end if

      num_report = num_report+1

      ! 3.3 Check if outside of the time range:

      read (platform%info%date_char,'(i4,4(1x,i2))') &
            iyear, imonth, iday, ihour, imin
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
      if ( obs_time < time_slots(0) .or. &
           obs_time >= time_slots(num_fgat_time) ) then
         num_outside_time = num_outside_time + 1
         if (print_detail_rain) then
            write(unit=stdout, fmt='(a)') '*** Outside of the time range:'
            write(unit=stdout, fmt=fmt_info) &
               platform%info%platform,    &
               platform%info%date_char,   &
               platform%stn%name
         end if
         cycle reports
      endif

      call da_llxy (platform%info, platform%loc, outside, outside_all)

      if (outside_all) then
         num_outside_all = num_outside_all + 1
         cycle reports
      end if

      nlevels = platform%info%levels

      dlat_earth = platform%info%lat
      dlon_earth = platform%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad

      ! 3.4 Loop over duplicating obs for global

      n_dup = 1
      if (global .and. &
         (platform%loc%i == ids .or. platform%loc%i == ide)) n_dup= 2
      do ndup = 1, n_dup
         select case (fm)

         case (129)
            if (ndup==1) ntotal = ntotal + 1
            if (outside) cycle reports
            if ( thin_rainobs ) then
               crit = 1.0
               call map2tgrid(rain,ifgat,dlat_earth,dlon_earth,dist,crit,itx,1,itt,iuse)
               if ( .not. iuse ) then
                  num_thinned = num_thinned + 1
                  cycle reports
               end if
            endif
            nlocal = nlocal + 1
            ilocal = nlocal

            iv % rain (ilocal) % stn_loc % lon = platform % stn % lon
            iv % rain (ilocal) % stn_loc % lat = platform % stn % lat
            iv % rain (ilocal) % stn_loc % elv = platform % stn % elv

            iv%info(rain)%levels(ilocal)    = nlevels
            iv%info(rain)%name(ilocal)      = platform%info%name
            iv%info(rain)%platform(ilocal)  = platform%info%platform
            iv%info(rain)%id(ilocal)        = platform%info%id
            iv%info(rain)%date_char(ilocal) = platform%info%date_char
            iv%info(rain)%lat(:,ilocal)     = platform%info%lat
            iv%info(rain)%lon(:,ilocal)     = platform%info%lon
            iv%info(rain)%elv(ilocal)       = platform%info%elv
            iv%info(rain)%pstar(ilocal)     = platform%info%pstar

            iv%info(rain)%slp(ilocal)           = platform%loc%slp
            iv%info(rain)%pw(ilocal)            = platform%loc%pw
            iv%info(rain)%x(:,ilocal)           = platform%loc%x
            iv%info(rain)%y(:,ilocal)           = platform%loc%y 
            iv%info(rain)%i(:,ilocal)           = platform%loc%i 
            iv%info(rain)%j(:,ilocal)           = platform%loc%j 
            iv%info(rain)%dx(:,ilocal)          = platform%loc%dx
            iv%info(rain)%dxm(:,ilocal)         = platform%loc%dxm
            iv%info(rain)%dy(:,ilocal)          = platform%loc%dy
            iv%info(rain)%dym(:,ilocal)         = platform%loc%dym
            iv%info(rain)%proc_domain(:,ilocal) = platform%loc%proc_domain

            iv%info(rain)%obs_global_index(ilocal) = ntotal

            iv % rain (ilocal) % height    = platform % each(1) % height
            iv % rain (ilocal) % height_qc = platform % each(1) % height_qc
            iv % rain (ilocal) % rain      = platform % each(1) % rain

         case default;
            write(unit=message(1), fmt='(a)') 'Unsaved obs found:'
            write(unit=message(2), fmt='(2a)') &
               'platform % info % platform=', platform % info % platform
            write(unit=message(3), fmt='(a, i3)') &
               'platform % info % levels=', platform % info % levels
            call da_warning("da_read_obs_rain.inc",259,message(1:3))
         end select

         if (global .and. ndup == 1) then
            if (platform%loc % i >= ide) then
               platform%loc%i = ids
               platform%loc%proc_domain = .false.
            else if (platform%loc % i <= ids) then
               platform%loc%i = ide
               platform%loc%proc_domain = .false.
            end if
         end if
      end do        !  loop over duplicate
   end do reports

   iv%info(rain)%ptotal(iv%time)=ntotal
   iv%info(rain)%plocal(iv%time)=nlocal

   write(unit=message(1),fmt='(A,4(1x,i7))') &
      'da_read_obs_rain: num_report, num_outside_all, num_outside_time, num_thinned: ', &
      num_report, num_outside_all, num_outside_time, num_thinned
   call da_message(message(1:1))

   close(iunit)
   call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_read_obs_rain")


end subroutine da_read_obs_rain


subroutine da_read_errfac(ob_name, f1, f2, f3, f4, f5)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   character (len=5), intent(in)  :: ob_name
   real,              intent(out) :: f1
   real,              intent(out) :: f2
   real,              intent(out) :: f3
   real,              intent(out) :: f4
   real,              intent(out) :: f5

   character (len=5)  :: ob_name1
   character (len=21) :: string1
   character (len=91) :: string2
   integer            :: fac_unit
   real               :: d1, d2, d3, d4, d5

   f1 = 1.0
   f2 = 1.0
   f3 = 1.0
   f4 = 1.0
   f5 = 1.0

   if (trace_use_dull) call da_trace_entry("da_read_errfac")

   call da_get_unit(fac_unit)
   open(unit=fac_unit, status='old', file = 'errfac.dat', iostat=ierr)

   if (ierr == 0) then
      do 
         read(unit=fac_unit,fmt='(1x,a5,a21,a91)')ob_name1, string1, string2

         if (ob_name == ob_name1 .and. string1 == ' obs, Error Factor = ') then
            read(unit=string2(17:31),fmt=*)d1
            read(unit=string2(32:46),fmt=*)d2
            read(unit=string2(47:61),fmt=*)d3
            read(unit=string2(62:76),fmt=*)d4
            read(unit=string2(77:91),fmt=*)d5
            if (d1 > 0.0) f1 = d1
            if (d2 > 0.0) f2 = d2
            if (d3 > 0.0) f3 = d3
            if (d4 > 0.0) f4 = d4
            if (d5 > 0.0) f5 = d5

            exit
         else if (ob_name1 == 'Total') then
            write(unit=message(1),fmt='(a,a)') ' No Tuning Error factors for ', ob_name
            write(unit=message(2),fmt='(a)') ' So setting to 1.0 i.e. default errors.'  
            call da_warning("da_read_errfac.inc",53,message(1:2))
            exit
         end if
      end do     
   else   
      call da_warning("da_read_errfac.inc",58, (/"Problem reading errfac.dat - Not tuning ob errors"/))
   end if

   close(fac_unit)
   call da_free_unit(fac_unit)

   if (trace_use_dull) call da_trace_exit("da_read_errfac")

end subroutine da_read_errfac


subroutine da_use_obs_errfac(iv)

   !-------------------------------------------------------------------------
   ! Purpose: Allocates observation structure and fills it from iv.
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv              ! Obs and header structure.

   integer                       :: n, k            ! Loop counters.
   real                          :: d1, d2, d3, d4  ! Dummy values.

   if (trace_use) call da_trace_entry("da_use_obs_errfac")

   !----------------------------------------------------------------------
   ! [2.0] Scale observation errors:
   !-------------------------------------------------------------------

   ! [2.1] Transfer surface obs:

   call da_read_errfac('synop', iv % synop_ef_u, &
                        iv % synop_ef_v, iv % synop_ef_t, &
                        iv % synop_ef_p, iv % synop_ef_q)                           

   if (iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
         iv % synop(n) % u % error = iv % synop(n) % u % error * iv % synop_ef_u
         iv % synop(n) % v % error = iv % synop(n) % v % error * iv % synop_ef_v
         iv % synop(n) % t % error = iv % synop(n) % t % error * iv % synop_ef_t
         iv % synop(n) % p % error = iv % synop(n) % p % error * iv % synop_ef_p
         iv % synop(n) % q % error = iv % synop(n) % q % error * iv % synop_ef_q
      end do
   end if

   ! [2.2] Transfer metar obs:


   call da_read_errfac('metar', iv % metar_ef_u, &
                        iv % metar_ef_v, iv % metar_ef_t, &
                        iv % metar_ef_p, iv % metar_ef_q)
                           
   if (iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
         iv % metar(n) % u % error = iv % metar(n) % u % error * iv % metar_ef_u
         iv % metar(n) % v % error = iv % metar(n) % v % error * iv % metar_ef_v
         iv % metar(n) % t % error = iv % metar(n) % t % error * iv % metar_ef_t
         iv % metar(n) % p % error = iv % metar(n) % p % error * iv % metar_ef_p
         iv % metar(n) % q % error = iv % metar(n) % q % error * iv % metar_ef_q
      end do
   end if

   ! [2.2] Transfer ships obs:

      
   call da_read_errfac('ships', iv % ships_ef_u, &
                        iv % ships_ef_v, iv % ships_ef_t, &
                        iv % ships_ef_p, iv % ships_ef_q)
                           
   if (iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
         iv % ships(n) % u % error = iv % ships(n) % u % error * iv % ships_ef_u
         iv % ships(n) % v % error = iv % ships(n) % v % error * iv % ships_ef_v
         iv % ships(n) % t % error = iv % ships(n) % t % error * iv % ships_ef_t
         iv % ships(n) % p % error = iv % ships(n) % p % error * iv % ships_ef_p
         iv % ships(n) % q % error = iv % ships(n) % q % error * iv % ships_ef_q
      end do
   end if

   ! [2.4.1] Transfer Geo. AMVs Obs:
   
   
   call da_read_errfac('geoamv', iv % geoamv_ef_u, iv % geoamv_ef_v, d1, d2, d3)

   if (iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
        do k = 1, iv%info(geoamv)%levels(n)
         iv % geoamv(n) % u(k) % error = iv % geoamv(n) % u(k) % error * iv % geoamv_ef_u
         iv % geoamv(n) % v(k) % error = iv % geoamv(n) % v(k) % error * iv % geoamv_ef_v
        end do
      end do 
   end if

   ! [2.4.2] Transfer Polar AMVs Obs:


   call da_read_errfac('polaramv', iv % polaramv_ef_u, iv % polaramv_ef_v, d1, d2, d3)

   if (iv%info(polaramv)%nlocal > 0) then 
      do n = 1, iv%info(polaramv)%nlocal
        do k = 1, iv%info(polaramv)%levels(n)
         iv % polaramv(n) % u(k) % error = iv % polaramv(n) % u(k) % error * iv % polaramv_ef_u
         iv % polaramv(n) % v(k) % error = iv % polaramv(n) % v(k) % error * iv % polaramv_ef_v
        end do
      end do
   end if


   ! [2.5] Transfer gpspw obs:


   call da_read_errfac('gpspw', iv % gpspw_ef_tpw, d1, d2, d3, d4)

   if (iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
         iv % gpspw(n) % tpw % error = iv % gpspw(n) % tpw % error * &
                                       iv % gpspw_ef_tpw

      end do
   end if

! [2.5.1] Transfer gpsref obs:

   call da_read_errfac('gpsre', iv % gpsref_ef_ref, d1, d2, d3, d4)

   if (iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
         do k = 1, iv%info(gpsref)%levels(n)
            iv % gpsref(n) % ref(k) % error = iv % gpsref(n) % ref(k) % error * &
                                              iv % gpsref_ef_ref
         enddo
      end do
   end if

   ! [2.6] Transfer sonde obs:


   call da_read_errfac('sound', iv % sound_ef_u, iv % sound_ef_v, &
                        iv % sound_ef_t, iv % sound_ef_q, d1)

   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         do k = 1, iv%info(sound)%levels(n)
            iv % sound(n) % u(k) % error = iv % sound(n) % u(k) % error * &
                                           iv % sound_ef_u
            iv % sound(n) % v(k) % error = iv % sound(n) % v(k) % error * &
                                           iv % sound_ef_v
            iv % sound(n) % t(k) % error = iv % sound(n) % t(k) % error * &
                                           iv % sound_ef_t
            iv % sound(n) % q(k) % error = iv % sound(n) % q(k) % error * &
                                           iv % sound_ef_q
         end do
      end do
   end if

   if (iv%info(sonde_sfc)%nlocal > 0) then
      do n = 1, iv%info(sonde_sfc)%nlocal
         iv % sonde_sfc(n) % u % error = iv % sonde_sfc(n) % u % error * iv % synop_ef_u
         iv % sonde_sfc(n) % v % error = iv % sonde_sfc(n) % v % error * iv % synop_ef_v
         iv % sonde_sfc(n) % t % error = iv % sonde_sfc(n) % t % error * iv % synop_ef_t
         iv % sonde_sfc(n) % p % error = iv % sonde_sfc(n) % p % error * iv % synop_ef_p
         iv % sonde_sfc(n) % q % error = iv % sonde_sfc(n) % q % error * iv % synop_ef_q
      end do
   end if

   ! [2.7] Transfer airep obs:

   
   call da_read_errfac('airep', iv % airep_ef_u, iv % airep_ef_v, &
                        iv % airep_ef_t, iv % airep_ef_q, d1)

   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         do k = 1, iv%info(airep)%levels(n)
            iv % airep(n) % u(k) % error = iv % airep(n) % u(k) % error * &
                                           iv % airep_ef_u
            iv % airep(n) % v(k) % error = iv % airep(n) % v(k) % error * &
                                           iv % airep_ef_v
            iv % airep(n) % t(k) % error = iv % airep(n) % t(k) % error * &
                                           iv % airep_ef_t
            iv % airep(n) % q(k) % error = iv % airep(n) % q(k) % error * &
                                           iv % airep_ef_q
         end do
      end do
   end if

   ! [2.8] Transfer pilot obs:


   call da_read_errfac('pilot', iv % pilot_ef_u, iv % pilot_ef_v, d1, d2, d3)

   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
         do k = 1, iv%info(pilot)%levels(n)
            iv % pilot(n) % u(k) % error = iv % pilot(n) % u(k) % error * &
                                           iv % pilot_ef_u
            iv % pilot(n) % v(k) % error = iv % pilot(n) % v(k) % error * &
                                           iv % pilot_ef_v

         end do
      end do
   end if

   ! [2.9] Transfer SSM/I obs:SSMI:


   call da_read_errfac('ssmir', iv % ssmir_ef_speed, iv % ssmir_ef_tpw, d1, d2, d3)

   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
         iv%ssmi_rv(n)%tpw%error = iv%ssmi_rv(n)%tpw%error * &
                                          iv % ssmir_ef_tpw
         iv%ssmi_rv(n)%speed%error = iv%ssmi_rv(n)%speed%error * &
                                            iv % ssmir_ef_speed
      end do
   end if


   ! iv % ssmit_ef_tb19h = 1.0 ! Tuning not yet coded.
   ! iv % ssmit_ef_tb19v = 1.0 ! Tuning not yet coded.
   ! iv % ssmit_ef_tb22v = 1.0 ! Tuning not yet coded.
   ! iv % ssmit_ef_tb37h = 1.0 ! Tuning not yet coded.
   ! iv % ssmit_ef_tb37v = 1.0 ! Tuning not yet coded.
   ! iv % ssmit_ef_tb85h = 1.0 ! Tuning not yet coded.
   ! iv % ssmit_ef_tb85v = 1.0 ! Tuning not yet coded.

   if (iv%info(ssmi_tb)%nlocal > 0) then
      ! do n = 1, iv%info(ssmi_tb)%nlocal
      !    iv%ssmi_tb(n)%tb19h%error = iv%ssmi_tb(n)%tb19h%error
      !    iv%ssmi_tb(n)%tb19v%error = iv%ssmi_tb(n)%tb19v%error
      !    iv%ssmi_tb(n)%tb22v%error = iv%ssmi_tb(n)%tb22v%error
      !    iv%ssmi_tb(n)%tb37h%error = iv%ssmi_tb(n)%tb37h%error * &
      !                                fac_ssmit_tb37h
      !    iv%ssmi_tb(n)%tb37v%error = iv%ssmi_tb(n)%tb37v%error * &
      !                                fac_ssmit_tb37v
      !    iv%ssmi_tb(n)%tb85h%error = iv%ssmi_tb(n)%tb85h%error * &
      !                                fac_ssmit_tb85h
      !    iv%ssmi_tb(n)%tb85v%error = iv%ssmi_tb(n)%tb85v%error * &
      !                                fac_ssmit_tb85v
      ! end do
   end if

   ! [2.10] Transfer satem obs:

   call da_read_errfac('satem', iv % satem_ef_thickness, d1, d2, d3, d4)

   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         do k = 1, iv%info(satem)%levels(n)
            iv % satem(n) % thickness(k) % error = iv % satem(n) % thickness(k) % error*&
                                                   iv % satem_ef_thickness
         end do
      end do
   end if
   
   ! [2.11] Transfer ssmt1 obs:

   call da_read_errfac('ssmt1', iv % ssmt1_ef_t, d1, d2, d3, d4)
      
   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
         do k = 1, iv%info(ssmt1)%levels(n)   
            iv % ssmt1(n) % t(k) % error = iv % ssmt1(n) % t(k) % error * &
                                        iv % ssmt1_ef_t
         end do
      end do
   end if

   ! [2.12] Transfer ssmt2 obs:

   call da_read_errfac('ssmt2', iv % ssmt2_ef_rh, d1, d2, d3, d4)

   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         do k = 1, iv%info(ssmt2)%levels(n)      
            iv % ssmt2(n) % rh(k) % error = iv % ssmt2(n) % rh(k) % error * &
                                         iv % ssmt2_ef_rh
         end do
      end do
   end if
   
   ! [2.13] Transfer scatterometer obs:

   call da_read_errfac('qscat', iv % qscat_ef_u, &
                        iv % qscat_ef_v, d1, d2, d3)
                           
   if (iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
         iv % qscat(n) % u % error = iv % qscat(n) % u % error * iv % qscat_ef_u
         iv % qscat(n) % v % error = iv % qscat(n) % v % error * iv % qscat_ef_v
      end do
   end if

   ! [2.14] Transfer profiler obs:

   call da_read_errfac('profi', iv % profiler_ef_u, iv % profiler_ef_v, d1, d2, d3)

   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         do k = 1, iv%info(profiler)%levels(n)
            iv % profiler(n) % u(k) % error = iv % profiler(n) % u(k) % error * &
                                           iv % profiler_ef_u
            iv % profiler(n) % v(k) % error = iv % profiler(n) % v(k) % error * &
                                           iv % profiler_ef_v

         end do
      end do
   end if

   ! [2.15] Transfer buoy obs:

   call da_read_errfac('buoy ', iv % buoy_ef_u, &
                        iv % buoy_ef_v, iv % buoy_ef_t, &
                        iv % buoy_ef_p, iv % buoy_ef_q)
                           
   if (iv%info(buoy)%nlocal > 0) then  
      do n = 1, iv%info(buoy)%nlocal
         iv % buoy(n) % u % error = iv % buoy(n) % u % error * iv % buoy_ef_u
         iv % buoy(n) % v % error = iv % buoy(n) % v % error * iv % buoy_ef_v
         iv % buoy(n) % t % error = iv % buoy(n) % t % error * iv % buoy_ef_t
         iv % buoy(n) % p % error = iv % buoy(n) % p % error * iv % buoy_ef_p
         iv % buoy(n) % q % error = iv % buoy(n) % q % error * iv % buoy_ef_q
      end do
   end if

   ! [2.16] Transfer TC bogus obs:

   call da_read_errfac('bogus', iv % bogus_ef_u, iv % bogus_ef_v, &
                        iv % bogus_ef_t, iv % bogus_ef_q, iv % bogus_ef_slp)

   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         do k = 1, iv%info(bogus)%levels(n)
            iv % bogus(n) % u(k) % error = iv % bogus(n) % u(k) % error * &
                                           iv % bogus_ef_u
            iv % bogus(n) % v(k) % error = iv % bogus(n) % v(k) % error * &
                                           iv % bogus_ef_v
            iv % bogus(n) % t(k) % error = iv % bogus(n) % t(k) % error * &
                                           iv % bogus_ef_t
            iv % bogus(n) % q(k) % error = iv % bogus(n) % q(k) % error * &
                                           iv % bogus_ef_q

         end do

         iv % bogus(n) % slp % error = iv % bogus(n) % slp % error * iv % bogus_ef_slp
      end do
   end if

   ! Transfer AIRS retrievals:

   call da_read_errfac('airsr', iv % airsr_ef_t, iv % airsr_ef_q, d1, d3, d3)

   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         do k = 1, iv%info(airsr)%levels(n)
            iv % airsr(n) % t(k) % error = iv % airsr(n) % t(k) % error * &
                                           iv % airsr_ef_t
            iv % airsr(n) % q(k) % error = iv % airsr(n) % q(k) % error * &
                                           iv % airsr_ef_q
         end do
      end do
   end if

   ! Transfer mtgirs obs:

   call da_read_errfac('mtgirs', iv % mtgirs_ef_u, iv % mtgirs_ef_v, &
                        iv % mtgirs_ef_t, iv % mtgirs_ef_q, d1)

   if (iv%info(mtgirs)%nlocal > 0) then
      do n = 1, iv%info(mtgirs)%nlocal
         do k = 1, iv%info(mtgirs)%levels(n)
            iv % mtgirs(n) % u(k) % error = iv % mtgirs(n) % u(k) % error * &
                                           iv % mtgirs_ef_u
            iv % mtgirs(n) % v(k) % error = iv % mtgirs(n) % v(k) % error * &
                                           iv % mtgirs_ef_v
            iv % mtgirs(n) % t(k) % error = iv % mtgirs(n) % t(k) % error * &
                                           iv % mtgirs_ef_t
            iv % mtgirs(n) % q(k) % error = iv % mtgirs(n) % q(k) % error * &
                                           iv % mtgirs_ef_q

         end do

      end do
   end if

   ! Transfer tamdar obs:

   if (iv%info(tamdar)%nlocal > 0) then
      call da_read_errfac('tamdar', iv % tamdar_ef_u, iv % tamdar_ef_v, &
                           iv % tamdar_ef_t, iv % tamdar_ef_q, d1)

      do n = 1, iv%info(tamdar)%nlocal
         do k = 1, iv%info(tamdar)%levels(n)
            iv % tamdar(n) % u(k) % error = iv % tamdar(n) % u(k) % error * &
                                           iv % tamdar_ef_u
            iv % tamdar(n) % v(k) % error = iv % tamdar(n) % v(k) % error * &
                                           iv % tamdar_ef_v
            iv % tamdar(n) % t(k) % error = iv % tamdar(n) % t(k) % error * &
                                           iv % tamdar_ef_t
            iv % tamdar(n) % q(k) % error = iv % tamdar(n) % q(k) % error * &
                                           iv % tamdar_ef_q
         end do
      end do
    end if

   ! Transfer tamdar_sfc obs:

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      do n = 1, iv%info(tamdar_sfc)%nlocal
         iv % tamdar_sfc(n) % u % error = iv % tamdar_sfc(n) % u % error * iv % tamdar_sfc_ef_u
         iv % tamdar_sfc(n) % v % error = iv % tamdar_sfc(n) % v % error * iv % tamdar_sfc_ef_v
         iv % tamdar_sfc(n) % t % error = iv % tamdar_sfc(n) % t % error * iv % tamdar_sfc_ef_t
         iv % tamdar_sfc(n) % p % error = iv % tamdar_sfc(n) % p % error * iv % tamdar_sfc_ef_p
         iv % tamdar_sfc(n) % q % error = iv % tamdar_sfc(n) % q % error * iv % tamdar_sfc_ef_q
      end do
   end if


   if (trace_use) call da_trace_exit("da_use_obs_errfac")

end subroutine da_use_obs_errfac


subroutine da_write_obs(it,ob, iv, re)

   !-------------------------------------------------------------------------
   ! Purpose: Writes out components of iv=O-B structure.
   !-------------------------------------------------------------------------   

   implicit none

   integer,        intent(in)    :: it
   type (y_type),  intent(in)    :: ob      ! Observation structure.
   type (iv_type), intent(in)    :: iv      ! O-B structure.
   type (y_type),  intent(inout) :: re      ! residual vector.
      
   integer                     :: n, k, num_obs, ios
   integer                     :: ounit     ! Output unit           
   character(len=filename_len) :: filename

   if (trace_use) call da_trace_entry("da_write_obs")

   !-------------------------------------------------------------------------
   ! Fix output unit
   !-------------------------------------------------------------------------

   call da_get_unit(ounit)

    write(unit=filename, fmt='(a,i2.2,a,i4.4)') 'gts_omb_oma_',it,'.', myproc

   open (unit=ounit,file=trim(filename),form='formatted',status='replace', &
      iostat=ios)
   if (ios /= 0) then
      call da_error("da_write_obs.inc",35, &
         (/"Cannot open conventional observation omb and oma file"//filename/))
   end if

   num_obs = 0
   do n = 1, iv%info(synop)%nlocal
      if (iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'synop', num_obs  
      num_obs = 0
      do n = 1, iv%info(synop)%nlocal  
         if (iv%info(synop)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1                                 
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
               num_obs , 1, iv%info(synop)%id(n), &  ! Station
               iv%info(synop)%lat(1,n), &       ! Latitude
               iv%info(synop)%lon(1,n), &       ! Longitude
               ob%synop(n)%p,           &       ! Obs Pressure
               ob%synop(n)%u,           & 
               iv%synop(n)%u%inv, iv%synop(n)%u%qc, iv%synop(n)%u%error, &
               re%synop(n)%u,           &
               ob%synop(n)%v,           &
               iv%synop(n)%v%inv, iv%synop(n)%v%qc, iv%synop(n)%v%error, &
               re%synop(n)%v,           &
               ob%synop(n)%t,           &
               iv%synop(n)%t%inv, iv%synop(n)%t%qc, iv%synop(n)%t%error, &
               re%synop(n)%t,           &
               ob%synop(n)%p,           &
               iv%synop(n)%p%inv, iv%synop(n)%p%qc, iv%synop(n)%p%error, &
               re%synop(n)%p,           &
               ob%synop(n)%q,           &
               iv%synop(n)%q%inv, iv%synop(n)%q%qc, iv%synop(n)%q%error, & 
               re%synop(n)%q
         end if
      end do
   end if

  num_obs = 0
  do n = 1, iv%info(metar)%nlocal
     if (iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
  end do
  if (num_obs > 0) then
     write(ounit,'(a20,i8)')'metar', num_obs  
     num_obs = 0
     do n = 1, iv%info(metar)%nlocal  
        if (iv%info(metar)%proc_domain(1,n)) then
           num_obs = num_obs + 1
           write(ounit,'(i8)') 1                                 
           write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
              num_obs  , 1, iv%info(metar)%id(n), &  ! Station
              iv%info(metar)%lat(1,n), &       ! Latitude
              iv%info(metar)%lon(1,n), &       ! Longitude
              ob%metar(n)%p,           &       ! Obs Pressure
              ob%metar(n)%u,           &
              iv%metar(n)%u%inv, iv%metar(n)%u%qc, iv%metar(n)%u%error, &
              re%metar(n)%u,           &
              ob%metar(n)%v,           &
              iv%metar(n)%v%inv, iv%metar(n)%v%qc, iv%metar(n)%v%error, &
              re%metar(n)%v,           &
              ob%metar(n)%t,           &
              iv%metar(n)%t%inv, iv%metar(n)%t%qc, iv%metar(n)%t%error, &
              re%metar(n)%t,           &
              ob%metar(n)%p,           &
              iv%metar(n)%p%inv, iv%metar(n)%p%qc, iv%metar(n)%p%error, &
              re%metar(n)%p,           &
              ob%metar(n)%q,           &
              iv%metar(n)%q%inv, iv%metar(n)%q%qc, iv%metar(n)%q%error, &
              re%metar(n)%q
        end if
     end do
  end if

  num_obs = 0
  do n = 1, iv%info(ships)%nlocal
     if (iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
  end do
  if (num_obs > 0) then
     write(ounit,'(a20,i8)')'ships', num_obs    
     num_obs = 0
     do n = 1, iv%info(ships)%nlocal  
        if (iv%info(ships)%proc_domain(1,n)) then
           write(ounit,'(i8)') 1                                 
           num_obs = num_obs + 1
           write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
              num_obs,1, iv%info(ships)%id(n), &  ! Station
              iv%info(ships)%lat(1,n), &       ! Latitude
              iv%info(ships)%lon(1,n), &       ! Longitude
              ob%ships(n)%p,           &       ! Obs Pressure
              ob%ships(n)%u,           &
              iv%ships(n)%u%inv, iv%ships(n)%u%qc, iv%ships(n)%u%error, &
              re%ships(n)%u,           &
              ob%ships(n)%v,           &
              iv%ships(n)%v%inv, iv%ships(n)%v%qc, iv%ships(n)%v%error, &
              re%ships(n)%v,           &
              ob%ships(n)%t,           &
              iv%ships(n)%t%inv, iv%ships(n)%t%qc, iv%ships(n)%t%error, &
              re%ships(n)%t,           &
              ob%ships(n)%p,           &
              iv%ships(n)%p%inv, iv%ships(n)%p%qc, iv%ships(n)%p%error, &
              re%ships(n)%p,           &
              ob%ships(n)%q,           &
              iv%ships(n)%q%inv, iv%ships(n)%q%qc, iv%ships(n)%q%error, &
              re%ships(n)%q
        end if
     end do
  end if

  num_obs = 0
  do n = 1, iv%info(geoamv)%nlocal
     if (iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
  end do
  if (num_obs > 0) then
     write(ounit,'(a20,i8)')'geoamv', num_obs    
     num_obs = 0
     do n = 1, iv%info(geoamv)%nlocal
        if (iv%info(geoamv)%proc_domain(1,n)) then                  
           num_obs = num_obs + 1
           write(ounit,'(i8)')iv%info(geoamv)%levels(n)
           do k = 1, iv%info(geoamv)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs, 1, iv%info(geoamv)%id(n), &  ! Station
                  iv%info(geoamv)%lat(k,n), &       ! Latitude
                  iv%info(geoamv)%lon(k,n), &       ! Longitude
                  iv%geoamv(n)%p(k),        &       ! Obs Pressure
                  ob%geoamv(n)%u(k),        &
                  iv%geoamv(n)%u(k)%inv, iv%geoamv(n)%u(k)%qc, iv%geoamv(n)%u(k)%error, &
                  re%geoamv(n)%u(k),        &
                  ob%geoamv(n)%v(k),        &
                  iv%geoamv(n)%v(k)%inv, iv%geoamv(n)%v(k)%qc, iv%geoamv(n)%v(k)%error, &
                  re%geoamv(n)%v(k)
           end do
        end if
     end do
  end if

   num_obs = 0
   do n = 1, iv%info(polaramv)%nlocal
      if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'polaramv', num_obs      
      num_obs = 0
      do n = 1, iv%info(polaramv)%nlocal
         if (iv%info(polaramv)%proc_domain(1,n)) then                    
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(polaramv)%levels(n)
            do k = 1, iv%info(polaramv)%levels(n)
                write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                   num_obs, 1, iv%info(polaramv)%id(n), &  ! Station
                   iv%info(polaramv)%lat(k,n), &       ! Latitude
                   iv%info(polaramv)%lon(k,n), &       ! Longitude
                   iv%polaramv(n)%p(k),        &       ! Obs Pressure
                   ob%polaramv(n)%u(k),        &
                   iv%polaramv(n)%u(k)%inv, iv%polaramv(n)%u(k)%qc, iv%polaramv(n)%u(k)%error, &
                   re%polaramv(n)%u(k),        &
                   ob%polaramv(n)%v(k),        &
                   iv%polaramv(n)%v(k)%inv, iv%polaramv(n)%v(k)%qc, iv%polaramv(n)%v(k)%error, &
                   re%polaramv(n)%v(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(gpspw)%nlocal   
      if (iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'gpspw', num_obs    
      num_obs = 0
      do n = 1, iv%info(gpspw)%nlocal
         if (iv%info(gpspw)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1                                 
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
               num_obs, 1, iv%info(gpspw)%id(n), &  ! Station
               iv%info(gpspw)%lat(1,n), &       ! Latitude
               iv%info(gpspw)%lon(1,n), &       ! Longitude
               iv%info(gpspw)%elv(n)  , &
               ob%gpspw(n)%tpw,         &
               iv%gpspw(n)%tpw%inv, iv%gpspw(n)%tpw%qc, iv%gpspw(n)%tpw%error, &
               re%gpspw(n)%tpw
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(sound)%nlocal
     if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'sound', num_obs    
      num_obs = 0
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(sound)%levels(n)
            do k = 1, iv%info(sound)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs,k, iv%info(sound)%id(n), &  ! Station
                  iv%info(sound)%lat(k,n), &       ! Latitude
                  iv%info(sound)%lon(k,n), &       ! Longitude
                  iv%sound(n)%p(k),        &       ! Obs Pressure
                  ob%sound(n)%u(k),        &
                  iv%sound(n)%u(k)%inv, iv%sound(n)%u(k)%qc, iv%sound(n)%u(k)%error, &
                  re%sound(n)%u(k),        &
                  ob%sound(n)%v(k),        &
                  iv%sound(n)%v(k)%inv, iv%sound(n)%v(k)%qc, iv%sound(n)%v(k)%error, &
                  re%sound(n)%v(k),        &
                  ob%sound(n)%t(k),        &
                  iv%sound(n)%t(k)%inv, iv%sound(n)%t(k)%qc, iv%sound(n)%t(k)%error, &
                  re%sound(n)%t(k),        &
                  ob%sound(n)%q(k),        &
                  iv%sound(n)%q(k)%inv, iv%sound(n)%q(k)%qc, iv%sound(n)%q(k)%error, &
                  re%sound(n)%q(k)
            end do
         end if
      end do
   end if

   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'sonde_sfc', num_obs    
      num_obs = 0
      do n = 1, iv%info(sonde_sfc)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) then 
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
               num_obs , 1, iv%info(sonde_sfc)%id(n), &  ! Station
               iv%info(sonde_sfc)%lat(1,n), &       ! Latitude
               iv%info(sonde_sfc)%lon(1,n), &       ! Longitude
               ob%sonde_sfc(n)%p,           &       ! Obs Pressure
               ob%sonde_sfc(n)%u,           &
               iv%sonde_sfc(n)%u%inv, iv%sonde_sfc(n)%u%qc, iv%sonde_sfc(n)%u%error, &
               re%sonde_sfc(n)%u,           &
               ob%sonde_sfc(n)%v,           &
               iv%sonde_sfc(n)%v%inv, iv%sonde_sfc(n)%v%qc, iv%sonde_sfc(n)%v%error, & 
               re%sonde_sfc(n)%v,           &
               ob%sonde_sfc(n)%t,           &
               iv%sonde_sfc(n)%t%inv, iv%sonde_sfc(n)%t%qc, iv%sonde_sfc(n)%t%error, &
               re%sonde_sfc(n)%t,           &
               ob%sonde_sfc(n)%p,           &
               iv%sonde_sfc(n)%p%inv, iv%sonde_sfc(n)%p%qc, iv%sonde_sfc(n)%p%error, & 
               re%sonde_sfc(n)%p,           &
               ob%sonde_sfc(n)%q,           &
               iv%sonde_sfc(n)%q%inv, iv%sonde_sfc(n)%q%qc, iv%sonde_sfc(n)%q%error, &
               re%sonde_sfc(n)%q
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(airep)%nlocal
      if (iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'airep', num_obs  
      num_obs = 0
      do n = 1, iv%info(airep)%nlocal
         if (iv%info(airep)%proc_domain(1,n)) then                  
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(airep)%levels(n)
            do k = 1, iv%info(airep)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,4(2f17.7,i8,2f17.7),f12.2)')&
                  num_obs, k, iv%info(airep)%id(n), &  ! Station
                  iv%info(airep)%lat(k,n), &       ! Latitude
                  iv%info(airep)%lon(k,n), &       ! Longitude
                  iv%airep(n)%p(k),        &       ! Obs pressure
                  ob%airep(n)%u(k),        &
                  iv%airep(n)%u(k)%inv, iv%airep(n)%u(k)%qc, iv%airep(n)%u(k)%error, & 
                  re%airep(n)%u(k),        &
                  ob%airep(n)%v(k),        &
                  iv%airep(n)%v(k)%inv, iv%airep(n)%v(k)%qc, iv%airep(n)%v(k)%error, & 
                  re%airep(n)%v(k),        &
                  ob%airep(n)%t(k),        &
                  iv%airep(n)%t(k)%inv, iv%airep(n)%t(k)%qc, iv%airep(n)%t(k)%error, & 
                  re%airep(n)%t(k),        &
                  ob%airep(n)%q(k),        &
                  iv%airep(n)%q(k)%inv, iv%airep(n)%q(k)%qc, iv%airep(n)%q(k)%error, &
                  re%airep(n)%q(k), iv%info(airep)%zk(k,n)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(pilot)%nlocal
      if (iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'pilot', num_obs   
      num_obs = 0
      do n = 1, iv%info(pilot)%nlocal
         if (iv%info(pilot)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(pilot)%levels(n)
            do k = 1, iv%info(pilot)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs, k, iv%info(pilot)%id(n), &  ! Station
                  iv%info(pilot)%lat(1,n), &       ! Latitude
                  iv%info(pilot)%lon(1,n), &       ! Longitude
                  iv%pilot(n)%h(k),        &       ! Obs Height
                  ob%pilot(n)%u(k),        &
                  iv%pilot(n)%u(k)%inv, iv%pilot(n)%u(k)%qc, iv%pilot(n)%u(k)%error, & 
                  re%pilot(n)%u(k),        &
                  ob%pilot(n)%v(k),        &
                  iv%pilot(n)%v(k)%inv, iv%pilot(n)%v(k)%qc, iv%pilot(n)%v(k)%error, & 
                  re%pilot(n)%v(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(ssmi_rv)%nlocal
      if (iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmir',  num_obs
      num_obs = 0
      do n = 1, iv%info(ssmi_rv)%nlocal
         if (iv%info(ssmi_rv)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
               num_obs, 1, 'SSMIR',              &       ! Station
               iv%info(ssmi_rv)%lat(1,n), &! Latitude
               iv%info(ssmi_rv)%lon(1,n), &! Longitude
               missing_r,                 &       ! Obs height
               ob%ssmi_rv(n)%speed,       &
               iv%ssmi_rv(n)%speed%inv, iv%ssmi_rv(n)%speed%qc, iv%ssmi_rv(n)%speed%error, & 
               re%ssmi_rv(n)%speed,       &
               ob%ssmi_rv(n)%tpw,         &
               iv%ssmi_rv(n)%tpw%inv, iv%ssmi_rv(n)%tpw%qc, iv%ssmi_rv(n)%tpw%error, & 
               re%ssmi_rv(n)%tpw
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(ssmi_tb)%nlocal
      if (iv%info(ssmi_tb)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmiT', num_obs     
      num_obs = 0
      do n = 1, iv%info(ssmi_tb)%nlocal
         ! write(ounit,*)' SSMI radiance output not yet coded.'
         if (iv%info(ssmi_tb)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,7(2f17.7,i8,2f17.7))')&
               num_obs, 1, 'SSMIT',              &        ! Station
               iv%info(ssmi_tb)%lat(1,n), &! Latitude
               iv%info(ssmi_tb)%lon(1,n), &! Longitude
               missing_r,                 &       ! Obs height
               ob%ssmi_tb(n)%tb19h,       &
               iv%ssmi_tb(n)%tb19h%inv, iv%ssmi_tb(n)%tb19h%qc, iv%ssmi_tb(n)%tb19h%error, & 
               re%ssmi_tb(n)%tb19h,       &
               ob%ssmi_tb(n)%tb19v,       &
               iv%ssmi_tb(n)%tb19v%inv, iv%ssmi_tb(n)%tb19v%qc, iv%ssmi_tb(n)%tb19v%error, & 
               re%ssmi_tb(n)%tb19v,       &
               ob%ssmi_tb(n)%tb22v,       &
               iv%ssmi_tb(n)%tb22v%inv, iv%ssmi_tb(n)%tb22v%qc, iv%ssmi_tb(n)%tb22v%error, &
               re%ssmi_tb(n)%tb22v,       &
               ob%ssmi_tb(n)%tb37h,       &
               iv%ssmi_tb(n)%tb37h%inv, iv%ssmi_tb(n)%tb37h%qc, iv%ssmi_tb(n)%tb37h%error, &
               re%ssmi_tb(n)%tb37h,       &
               ob%ssmi_tb(n)%tb37v,       &
               iv%ssmi_tb(n)%tb37v%inv, iv%ssmi_tb(n)%tb37v%qc, iv%ssmi_tb(n)%tb37v%error, & 
               re%ssmi_tb(n)%tb37v,       &
               ob%ssmi_tb(n)%tb85h,       &
               iv%ssmi_tb(n)%tb85h%inv, iv%ssmi_tb(n)%tb85h%qc, iv%ssmi_tb(n)%tb85h%error, & 
               re%ssmi_tb(n)%tb85h,       &
               ob%ssmi_tb(n)%tb85v,       &
               iv%ssmi_tb(n)%tb85v%inv, iv%ssmi_tb(n)%tb85v%qc, iv%ssmi_tb(n)%tb85v%error, & 
               re%ssmi_tb(n)%tb85v
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(satem)%nlocal
      if (iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'satem', num_obs    
      num_obs = 0
      do n = 1, iv%info(satem)%nlocal
         if (iv%info(satem)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(satem)%levels(n)
            do k = 1, iv%info(satem)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs   , k, iv%info(satem)%id(n), &  ! Station
                  iv%info(satem)%lat(1,n), &       ! Latitude
                  iv%info(satem)%lon(1,n), &       ! Longitude
                  iv%satem(n)%p(k),        &       ! Obs Pressure
                  ob%satem(n)%thickness(k),       &
                  iv%satem(n)%thickness(k)%inv,   &
                  iv%satem(n)%thickness(k)%qc,    &
                  iv%satem(n)%thickness(k)%error, &
                  re%satem(n)%thickness(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(ssmt1)%nlocal
      if (iv%info(ssmt1)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmt1', num_obs
      num_obs = 0
      do n = 1, iv%info(ssmt1)%nlocal
         if (iv%info(ssmt1)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(ssmt1)%levels(n)
            do k = 1, iv%info(ssmt1)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs , k, iv%info(ssmt1)%id(n), &  ! Station
                  iv%info(ssmt1)%lat(1,n), &       ! Latitude
                  iv%info(ssmt1)%lon(1,n), &       ! Longitude
                  iv%ssmt1(n)%h(k),       &        ! Obs height
                  ob%ssmt1(n)%t(k),       &
                  iv%ssmt1(n)%t(k)%inv,   &
                  iv%ssmt1(n)%t(k)%qc,    &
                  iv%ssmt1(n)%t(k)%error, &
                  re%ssmt1(n)%t(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(ssmt2)%nlocal
      if (iv%info(ssmt2)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmt2', num_obs    
      num_obs = 0
      do n = 1, iv%info(ssmt2)%nlocal
         if (iv%info(ssmt2)%proc_domain(1,n)) then                   
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(ssmt2)%levels(n)
            do k = 1, iv%info(ssmt2)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs , k, iv%info(ssmt2)%id(n), &  ! Station
                  iv%info(ssmt2)%lat(1,n), &       ! Latitude
                  iv%info(ssmt2)%lon(1,n), &       ! Longitude
                  iv%ssmt2(n)%h(k),        &       ! Obs height
                  ob%ssmt2(n)%rh(k),       &
                  iv%ssmt2(n)%rh(k)%inv,   &
                  iv%ssmt2(n)%rh(k)%qc,    &
                  iv%ssmt2(n)%rh(k)%error, &
                  re%ssmt2(n)%rh(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(qscat)%nlocal  
      if (iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'qscat', num_obs   
      num_obs = 0
      do n = 1, iv%info(qscat)%nlocal  
         if (iv%info(qscat)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                num_obs, 1, iv%info(qscat)%id(n), &  ! Station
                iv%info(qscat)%lat(1,n), &       ! Latitude
                iv%info(qscat)%lon(1,n), &       ! Longitude
                iv%qscat(n)%h,           &       ! Obs height
                ob%qscat(n)%u,           &
                iv%qscat(n)%u%inv, iv%qscat(n)%u%qc, iv%qscat(n)%u%error, &
                re%qscat(n)%u,           &
                ob%qscat(n)%v,           &
                iv%qscat(n)%v%inv, iv%qscat(n)%v%qc, iv%qscat(n)%v%error, & 
                re%qscat(n)%v
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(profiler)%nlocal
      if (iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'profiler',  num_obs
      num_obs = 0
      do n = 1, iv%info(profiler)%nlocal
         if (iv%info(profiler)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(profiler)%levels(n)
            do k = 1, iv%info(profiler)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs, k, iv%info(profiler)%id(n), &  ! Station
                  iv%info(profiler)%lat(1,n), &       ! Latitude
                  iv%info(profiler)%lon(1,n), &       ! Longitude
                  iv%profiler(n)%h(k),        &       ! Obs Height 
                  ob%profiler(n)%u(k),        &
                  iv%profiler(n)%u(k)%inv, iv%profiler(n)%u(k)%qc, iv%profiler(n)%u(k)%error, &
                  re%profiler(n)%u(k),        &
                  ob%profiler(n)%v(k),        &
                  iv%profiler(n)%v(k)%inv, iv%profiler(n)%v(k)%qc, iv%profiler(n)%v(k)%error, &
                  re%profiler(n)%v(k) 
            end do
         end if 
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(buoy)%nlocal  
      if (iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'buoy', num_obs
      num_obs = 0
      do n = 1, iv%info(buoy)%nlocal  
         if (iv%info(buoy)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
               num_obs,1, iv%info(buoy)%id(n), &  ! Station
               iv%info(buoy)%lat(1,n), &       ! Latitude
               iv%info(buoy)%lon(1,n), &       ! Longitude
               ob%buoy(n)%p,           &       ! Obs Pressure
               ob%buoy(n)%u,           &
               iv%buoy(n)%u%inv, iv%buoy(n)%u%qc, iv%buoy(n)%u%error, &
               re%buoy(n)%u,           &
               ob%buoy(n)%v,           &
               iv%buoy(n)%v%inv, iv%buoy(n)%v%qc, iv%buoy(n)%v%error, &
               re%buoy(n)%v,           &
               ob%buoy(n)%t,           &
               iv%buoy(n)%t%inv, iv%buoy(n)%t%qc, iv%buoy(n)%t%error, &
               re%buoy(n)%t,           &
               ob%buoy(n)%p,           &
               iv%buoy(n)%p%inv, iv%buoy(n)%p%qc, iv%buoy(n)%p%error, &
               re%buoy(n)%p,           &
               ob%buoy(n)%q,           &
               iv%buoy(n)%q%inv, iv%buoy(n)%q%qc, iv%buoy(n)%q%error, &
               re%buoy(n)%q
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(bogus)%nlocal
      if (iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'bogus', num_obs
      num_obs = 0
      do n = 1, iv%info(bogus)%nlocal
         if (iv%info(bogus)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                num_obs, 1, iv%info(bogus)%id(n), &  ! Station
                iv%info(bogus)%lat(1,n), &       ! Latitude
                iv%info(bogus)%lon(1,n), &       ! Longitude
                missing_r,               &
                ob%bogus(n)%slp,         &
                iv%bogus(n)%slp%inv, iv%bogus(n)%slp%qc, iv%bogus(n)%slp%error, &
                re%bogus(n)%slp    ! O, O-B, O-A p
            write(ounit,'(i8)')iv%info(bogus)%levels(n)
            do k = 1, iv%info(bogus)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs , k, iv%info(bogus)%id(n), &  ! Station
                  iv%info(bogus)%lat(1,n), &       ! Latitude
                  iv%info(bogus)%lon(1,n), &       ! Longitude
                  iv%bogus(n)%p(k),        &       ! Obs Pressure
                  ob%bogus(n)%u(k),        &
                  iv%bogus(n)%u(k)%inv, iv%bogus(n)%u(k)%qc, iv%bogus(n)%u(k)%error, &
                  re%bogus(n)%u(k),        &
                  ob%bogus(n)%v(k),        &
                  iv%bogus(n)%v(k)%inv, iv%bogus(n)%v(k)%qc, iv%bogus(n)%v(k)%error, &
                  re%bogus(n)%v(k),        &
                  ob%bogus(n)%t(k),        &
                  iv%bogus(n)%t(k)%inv, iv%bogus(n)%t(k)%qc, iv%bogus(n)%t(k)%error, &
                  re%bogus(n)%t(k),        &
                  ob%bogus(n)%q(k),        &
                  iv%bogus(n)%q(k)%inv, iv%bogus(n)%q(k)%qc, iv%bogus(n)%q(k)%error, &
                  re%bogus(n)%q(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(airsr)%nlocal
      if (iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'airsr', num_obs    
      num_obs = 0
      do n = 1, iv%info(airsr)%nlocal
         if (iv%info(airsr)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(airsr)%levels(n)
            do k = 1, iv%info(airsr)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs,k, iv%info(airsr)%id(n), &  ! Station
                  iv%info(airsr)%lat(1,n), &       ! Latitude
                  iv%info(airsr)%lon(1,n), &       ! Longitude
                  iv%airsr(n)%p(k),        &       ! Obs Pressure
                  ob%airsr(n)%t(k),        &
                  iv%airsr(n)%t(k)%inv, iv%airsr(n)%t(k)%qc, iv%airsr(n)%t(k)%error, &
                  re%airsr(n)%t(k),        &
                  ob%airsr(n)%q(k),        &
                  iv%airsr(n)%q(k)%inv, iv%airsr(n)%q(k)%qc, iv%airsr(n)%q(k)%error, &
                  re%airsr(n)%q(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(gpsref)%nlocal
      if (iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'gpsref', num_obs
       num_obs = 0
      do n = 1, iv%info(gpsref)%nlocal
         if (iv%info(gpsref)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(gpsref)%levels(n)
            do k = 1, iv%info(gpsref)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
                  num_obs,k, iv%info(gpsref)%id(n), &  ! Station
                  iv%info(gpsref)%lat(1,n), &       ! Latitude
                  iv%info(gpsref)%lon(1,n), &       ! Longitude
                  iv%gpsref(n)%h(k),        &       ! Obs Height    
                  ob%gpsref(n)%ref(k),      &
                  iv%gpsref(n)%ref(k)%inv, iv%gpsref(n)%ref(k)%qc, iv%gpsref(n)%ref(k)%error, &
                  re%gpsref(n)%ref(k)
            end do
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(tamdar)%nlocal
      if (iv%info(tamdar)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'tamdar', num_obs
      num_obs = 0
      do n = 1, iv%info(tamdar)%nlocal
         if (iv%info(tamdar)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)')iv%info(tamdar)%levels(n)
            do k = 1, iv%info(tamdar)%levels(n)
               write(ounit,'(2i8,a5,2f9.2,f17.7,4(2f17.7,i8,2f17.7),f12.2)')&
                  num_obs,k, iv%info(tamdar)%id(n), &  ! Station
                  iv%info(tamdar)%lat(k,n), &       ! Latitude
                  iv%info(tamdar)%lon(k,n), &       ! Longitude
                  iv%tamdar(n)%p(k),        &       ! Obs Pressure
                  ob%tamdar(n)%u(k),        &
                  iv%tamdar(n)%u(k)%inv, iv%tamdar(n)%u(k)%qc, iv%tamdar(n)%u(k)%error, &
                  re%tamdar(n)%u(k),        &
                  ob%tamdar(n)%v(k),        &
                  iv%tamdar(n)%v(k)%inv, iv%tamdar(n)%v(k)%qc, iv%tamdar(n)%v(k)%error, &
                  re%tamdar(n)%v(k),        &
                  ob%tamdar(n)%t(k),        &
                  iv%tamdar(n)%t(k)%inv, iv%tamdar(n)%t(k)%qc, iv%tamdar(n)%t(k)%error, &
                  re%tamdar(n)%t(k),        &
                  ob%tamdar(n)%q(k),        &
                  iv%tamdar(n)%q(k)%inv, iv%tamdar(n)%q(k)%qc, iv%tamdar(n)%q(k)%error, &
                  re%tamdar(n)%q(k), iv%info(tamdar)%zk(k,n)
            end do
         end if
      end do
   end if

! tamdar_sfc

   num_obs = 0
   do n = 1, iv%info(tamdar_sfc)%nlocal
      if (iv%info(tamdar_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'tamdar_sfc', num_obs  
      num_obs = 0
      do n = 1, iv%info(tamdar_sfc)%nlocal  
         if (iv%info(tamdar_sfc)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1                                 
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7),f12.2)')&
               num_obs , 1, iv%info(tamdar_sfc)%id(n), &  ! Station
               iv%info(tamdar_sfc)%lat(1,n), &       ! Latitude
               iv%info(tamdar_sfc)%lon(1,n), &       ! Longitude
               ob%tamdar_sfc(n)%p,           &       ! Obs Pressure
               ob%tamdar_sfc(n)%u,           & 
               iv%tamdar_sfc(n)%u%inv, iv%tamdar_sfc(n)%u%qc, iv%tamdar_sfc(n)%u%error, &
               re%tamdar_sfc(n)%u,           &
               ob%tamdar_sfc(n)%v,           &
               iv%tamdar_sfc(n)%v%inv, iv%tamdar_sfc(n)%v%qc, iv%tamdar_sfc(n)%v%error, &
               re%tamdar_sfc(n)%v,           &
               ob%tamdar_sfc(n)%t,           &
               iv%tamdar_sfc(n)%t%inv, iv%tamdar_sfc(n)%t%qc, iv%tamdar_sfc(n)%t%error, &
               re%tamdar_sfc(n)%t,           &
               ob%tamdar_sfc(n)%p,           &
               iv%tamdar_sfc(n)%p%inv, iv%tamdar_sfc(n)%p%qc, iv%tamdar_sfc(n)%p%error, &
               re%tamdar_sfc(n)%p,           &
               ob%tamdar_sfc(n)%q,           &
               iv%tamdar_sfc(n)%q%inv, iv%tamdar_sfc(n)%q%qc, iv%tamdar_sfc(n)%q%error, & 
               re%tamdar_sfc(n)%q, iv%info(tamdar_sfc)%zk(1,n)
         end if
      end do
   end if

   num_obs = 0
   do n = 1, iv%info(rain)%nlocal  
      if (iv%info(rain)%proc_domain(1,n)) num_obs = num_obs + 1
   end do
   if (num_obs > 0) then
      write(ounit,'(a20,i8)')'rain', num_obs
      num_obs = 0
      do n = 1, iv%info(rain)%nlocal  
         if (iv%info(rain)%proc_domain(1,n)) then
            num_obs = num_obs + 1
            write(ounit,'(i8)') 1
            write(ounit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))')&
               num_obs,1, iv%info(rain)%id(n), &  ! Station  
               iv%info(rain)%lat(1,n), &          ! Latitude
               iv%info(rain)%lon(1,n), &          ! Longitude
               iv%rain(n)%height,      &          ! Obs Pressure
               ob%rain(n)%rain,        &
               iv%rain(n)%rain%inv, iv%rain(n)%rain%qc, iv%rain(n)%rain%error, &
               re%rain(n)%rain
         end if
      end do
   end if
   
   close (ounit)
   call da_free_unit(ounit)

   if (trace_use) call da_trace_exit("da_write_obs")

end subroutine da_write_obs


subroutine da_write_iv_for_multi_inc(file_index, iv)

   !-------------------------------------------------------------------------
   ! Purpose: Writes out components of iv=O-B structure.
   !-------------------------------------------------------------------------   

   implicit none

   type (iv_type), intent(in)    :: iv      ! O-B structure.
   integer, intent (in)          :: file_index
      
   integer                       :: n, k, ios
   integer                       :: ounit     ! Output unit           
   character(len=filename_len)   :: filename

   if (trace_use) call da_trace_entry("da_write_iv_for_multi_inc")

   !-------------------------------------------------------------------------
   ! Fix output unit
   !-------------------------------------------------------------------------

   call da_get_unit(ounit)
   write(unit=filename, fmt='(a,i3.3,a,i4.4)') 'stub.', file_index, '.', myproc
   ! [1] surface obs:

   if (iv%info(synop)%plocal(iv%time) - iv%info(synop)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.synop',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",35, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'synop',iv%info(synop)%plocal(iv%time) - &
                                     iv%info(synop)%plocal(iv%time-1) 
      do n = iv%info(synop)%plocal(iv%time-1) + 1, &
             iv%info(synop)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
             n , iv%info(synop)%id(n), &  ! Station
             iv%info(synop)%lat(1,n), &       ! Latitude
             iv%info(synop)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,5(E22.13,i8,3E22.13))')&
             iv%synop(n)%h, &
             iv%synop(n)%u, &!  O-B u
             iv%synop(n)%v, &!  O-B v
             iv%synop(n)%t, &!  O-B t
             iv%synop(n)%p, &!  O-B p
             iv%synop(n)%q  !  O-B q
      end do
      close (ounit)
   end if

   ! [2] metar obs:

   if (iv%info(metar)%plocal(iv%time) - iv%info(metar)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.metar',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",65, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'metar', iv%info(metar)%plocal(iv%time) - &
                                       iv%info(metar)%plocal(iv%time-1)
      do n = iv%info(metar)%plocal(iv%time-1) + 1, &
             iv%info(metar)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(metar)%id(n), &  ! Station
                 iv%info(metar)%lat(1,n), &       ! Latitude
                 iv%info(metar)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,5(E22.13,i8,3E22.13))')&
                 iv%metar(n)%h, & 
                 iv%metar(n)%u, &! O-B u
                 iv%metar(n)%v, &! O-B v
                 iv%metar(n)%t, &! O-B t
                 iv%metar(n)%p, &! O-B p
                 iv%metar(n)%q   ! O-B q
      end do
      close (ounit)
   end if

   ! [3] ships obs:

   if (iv%info(ships)%plocal(iv%time) - iv%info(ships)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.ships',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",95, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'ships', iv%info(ships)%plocal(iv%time) - &
                                       iv%info(ships)%plocal(iv%time-1)
      do n = iv%info(ships)%plocal(iv%time-1) + 1, &
             iv%info(ships)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(ships)%id(n), &  ! Station
                 iv%info(ships)%lat(1,n), &       ! Latitude
                 iv%info(ships)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,5(E22.13,i8,3E22.13))')&
                 iv%ships(n)%h, &
                 iv%ships(n)%u, &! O-B u
                 iv%ships(n)%v, &! O-B v
                 iv%ships(n)%t, &! O-B t
                 iv%ships(n)%p, &! O-B p
                 iv%ships(n)%q   ! O-B q
      end do
      close (ounit)
   end if

   ! [4] sonde_sfc obs:

   if (iv%info(sonde_sfc)%plocal(iv%time) - iv%info(sonde_sfc)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.sonde_sfc',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",125, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'sonde_sfc', iv%info(sonde_sfc)%plocal(iv%time) - &
                                       iv%info(sonde_sfc)%plocal(iv%time-1)
      do n = iv%info(sonde_sfc)%plocal(iv%time-1) + 1, &
             iv%info(sonde_sfc)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(sonde_sfc)%id(n), &  ! Station
                 iv%info(sonde_sfc)%lat(1,n), &       ! Latitude
                 iv%info(sonde_sfc)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,5(E22.13,i8,3E22.13))')&
                 iv%sonde_sfc(n)%h, &
                 iv%sonde_sfc(n)%u, &! O-B u
                 iv%sonde_sfc(n)%v, &! O-B v
                 iv%sonde_sfc(n)%t, &! O-B t
                 iv%sonde_sfc(n)%p, &! O-B p
                 iv%sonde_sfc(n)%q   ! O-B q
      end do
      close (ounit)
   end if

   ! [5] sound obs:

   if (iv%info(sound)%plocal(iv%time) - iv%info(sound)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.sound',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",155, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'sound', iv%info(sound)%plocal(iv%time) - &
                                       iv%info(sound)%plocal(iv%time-1)
      do n = iv%info(sound)%plocal(iv%time-1) + 1, &
             iv%info(sound)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(sound)%levels(n), iv%info(sound)%id(n), &  ! Station
                 iv%info(sound)%lat(1,n), &       ! Latitude
                 iv%info(sound)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(sound)%levels(n)
            write(ounit,'(2E22.13,4(E22.13,i8,3E22.13))')&
                     iv%sound(n)%h(k), &
                     iv%sound(n)%p(k), &             ! Obs Pressure
                     iv%sound(n)%u(k), &! O-B u
                     iv%sound(n)%v(k), &! O-B v
                     iv%sound(n)%t(k), &! O-B t
                     iv%sound(n)%q(k)   ! O-B q
         enddo
      end do
      close (ounit)
   end if

   ! [6] mtgirs obs:
   
   if (iv%info(mtgirs)%plocal(iv%time) - iv%info(mtgirs)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.mtgirs',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",187, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'mtgirs', iv%info(mtgirs)%plocal(iv%time) - &
                                       iv%info(mtgirs)%plocal(iv%time-1)
      do n = iv%info(mtgirs)%plocal(iv%time-1) + 1, &
             iv%info(mtgirs)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(mtgirs)%levels(n), iv%info(mtgirs)%id(n), &  ! Station
                 iv%info(mtgirs)%lat(1,n), &       ! Latitude
                 iv%info(mtgirs)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(mtgirs)%levels(n)
            write(ounit,'(2E22.13,4(E22.13,i8,3E22.13))')&
                 iv % mtgirs(n) % h(k), &
                 iv % mtgirs(n) % p(k), &             ! Obs Pressure
                 iv%mtgirs(n)%u(k), &! O-B u
                 iv%mtgirs(n)%v(k), &! O-B v
                 iv%mtgirs(n)%t(k), &! O-B t
                 iv%mtgirs(n)%q(k)   ! O-B q

         enddo
      end do
      close (ounit)
   end if

   ! [7] tamdar

   if (iv%info(tamdar)%plocal(iv%time) - iv%info(tamdar)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.tamdar',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",220, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'tamdar', iv%info(tamdar)%plocal(iv%time) - &
                                       iv%info(tamdar)%plocal(iv%time-1)
      do n = iv%info(tamdar)%plocal(iv%time-1) + 1, &
             iv%info(tamdar)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(tamdar)%levels(n), iv%info(tamdar)%id(n), &  ! Station
                 iv%info(tamdar)%lat(1,n), &       ! Latitude
                 iv%info(tamdar)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(tamdar)%levels(n)
            write(ounit,'(2E22.13,4(E22.13,i8,3E22.13))')&
                     iv%tamdar(n)%h(k), &
                     iv%tamdar(n)%p(k), &             ! Obs Pressure
                     iv%tamdar(n)%u(k), &! O-B u
                     iv%tamdar(n)%v(k), &! O-B v
                     iv%tamdar(n)%t(k), &! O-B t
                     iv%tamdar(n)%q(k)   ! O-B q
         enddo
      end do
      close (ounit)
   end if

   ! [8] tamdar_sfc

   if (iv%info(tamdar_sfc)%plocal(iv%time) - iv%info(tamdar_sfc)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.tamdar_sfc',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",252, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'tamdar_sfc', iv%info(tamdar_sfc)%plocal(iv%time) - &
                                       iv%info(tamdar_sfc)%plocal(iv%time-1)
      do n = iv%info(tamdar_sfc)%plocal(iv%time-1) + 1, &
             iv%info(tamdar_sfc)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(tamdar_sfc)%id(n), &  ! Station
                 iv%info(tamdar_sfc)%lat(1,n), &       ! Latitude
                 iv%info(tamdar_sfc)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,5(E22.13,i8,3E22.13))')&
                 iv%tamdar_sfc(n)%h, &
                 iv%tamdar_sfc(n)%u, &!  O-B u
                 iv%tamdar_sfc(n)%v, &!  O-B v
                 iv%tamdar_sfc(n)%t, &!  O-B t
                 iv%tamdar_sfc(n)%p, &!  O-B p
                 iv%tamdar_sfc(n)%q  !  O-B q
      end do
      close (ounit)
   end if

   ! [9] buoy obs:

   if (iv%info(buoy)%plocal(iv%time) - iv%info(buoy)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.buoy',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",282, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'buoy', iv%info(buoy)%plocal(iv%time) - &
                                       iv%info(buoy)%plocal(iv%time-1)
      do n = iv%info(buoy)%plocal(iv%time-1) + 1, &
             iv%info(buoy)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(buoy)%id(n), &  ! Station
                 iv%info(buoy)%lat(1,n), &       ! Latitude
                 iv%info(buoy)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,5(E22.13,i8,3E22.13))')&
                 iv%buoy(n)%h, &
                 iv%buoy(n)%u, &! O-B u
                 iv%buoy(n)%v, &! O-B v
                 iv%buoy(n)%t, &! O-B t
                 iv%buoy(n)%p, &! O-B p
                 iv%buoy(n)%q   ! O-B q
      end do
      close (ounit)
   end if

   ! [10] Geo AMVs obs:

   if (iv%info(geoamv)%plocal(iv%time) - iv%info(geoamv)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.geoamv',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",312, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'geoamv', iv%info(geoamv)%plocal(iv%time) - &
                                       iv%info(geoamv)%plocal(iv%time-1)
      do n = iv%info(geoamv)%plocal(iv%time-1) + 1, &
             iv%info(geoamv)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(geoamv)%levels(n), iv%info(geoamv)%id(n), &     ! Station
                 iv%info(geoamv)%lat(1,n), &       ! Latitude
                 iv%info(geoamv)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(geoamv)%levels(n)
            write(ounit,'(E22.13,2(E22.13,i8,3E22.13))')&
                      iv%geoamv(n)%p(k), &                ! Obs Pressure
                      iv%geoamv(n)%u(k), &! O-B u
                      iv%geoamv(n)%v(k)
         enddo
      end do
      close (ounit)
   end if

   ! [11] gpspw obs:

   if (iv%info(gpspw)%plocal(iv%time) - iv%info(gpspw)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.gpspw',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",341, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'gpspw', iv%info(gpspw)%plocal(iv%time) - &
                                       iv%info(gpspw)%plocal(iv%time-1)
      do n = iv%info(gpspw)%plocal(iv%time-1) + 1, &
             iv%info(gpspw)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(gpspw)%id(n), &  ! Station
                 iv%info(gpspw)%lat(1,n), &       ! Latitude
                 iv%info(gpspw)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,i8,3E22.13)')&
              iv%gpspw(n)%tpw
      end do
      close (ounit)
   end if

   ! [12] SSM/I obs:

   if (iv%info(ssmi_rv)%plocal(iv%time) - iv%info(ssmi_rv)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.ssmir',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",366, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'ssmir', iv%info(ssmi_rv)%plocal(iv%time) - &
                                       iv%info(ssmi_rv)%plocal(iv%time-1)
      do n = iv%info(ssmi_rv)%plocal(iv%time-1) + 1, &
             iv%info(ssmi_rv)%plocal(iv%time)
         write(ounit,'(i8,2E22.13)')&
                 n, &  ! Station
                 iv%info(ssmi_rv)%lat(1,n), &       ! Latitude
                 iv%info(ssmi_rv)%lon(1,n)          ! Longitude
         write(ounit,'(2(E22.13,i8,3E22.13))')&
                 iv%ssmi_rv(n)%speed, & ! O-B speed
                 iv%ssmi_rv(n)%tpw ! O-BA tpw
      end do
      close (ounit)
   end if

   ! [13] airep obs:

   if (iv%info(airep)%plocal(iv%time) - iv%info(airep)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.airep',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",392, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'airep', iv%info(airep)%plocal(iv%time) - &
                                       iv%info(airep)%plocal(iv%time-1)
      do n = iv%info(airep)%plocal(iv%time-1) + 1, &
             iv%info(airep)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(airep)%levels(n), iv%info(airep)%id(n), &  ! Station
                 iv%info(airep)%lat(1,n), &       ! Latitude
                 iv%info(airep)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(airep)%levels(n)
            write(ounit,'(2E22.13,4(E22.13,i8,3E22.13))')&
                     iv%airep(n)%h(k), &
                     iv%airep(n)%p(k), &             ! Obs pressure
                     iv%airep(n)%u(k), &! O-B u
                     iv%airep(n)%v(k), &! O-B v
                     iv%airep(n)%t(k), &
                     iv%airep(n)%q(k)
         enddo
      end do
      close (ounit)
   end if

   ! [14] Polar AMVs obs:

   if (iv%info(polaramv)%plocal(iv%time) - iv%info(polaramv)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.polaramv',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",424, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'polaramv', iv%info(polaramv)%plocal(iv%time) - &
                                       iv%info(polaramv)%plocal(iv%time-1)
      do n = iv%info(polaramv)%plocal(iv%time-1) + 1, &
             iv%info(polaramv)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(polaramv)%levels(n), iv%info(polaramv)%id(n), &  ! Station
                 iv%info(polaramv)%lat(1,n), &       ! Latitude
                 iv%info(polaramv)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(polaramv)%levels(n)
            write(ounit,'(E22.13,2(E22.13,i8,3E22.13))')&
                      iv%polaramv(n)%p(k), &                ! Obs Pressure
                      iv%polaramv(n)%u(k), &! O-B u
                      iv%polaramv(n)%v(k)
         enddo
      end do
      close (ounit)
   end if

   ! [15] pilot obs:

   if (iv%info(pilot)%plocal(iv%time) - iv%info(pilot)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.pilot',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",453, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'pilot', iv%info(pilot)%plocal(iv%time) - &
                                       iv%info(pilot)%plocal(iv%time-1)
      do n = iv%info(pilot)%plocal(iv%time-1) + 1, &
             iv%info(pilot)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(pilot)%levels(n), iv%info(pilot)%id(n), &  ! Station
                 iv%info(pilot)%lat(1,n), &       ! Latitude
                 iv%info(pilot)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(pilot)%levels(n)
            write(ounit,'(E22.13,2(E22.13,i8,3E22.13))')&
                      iv%pilot(n)%p(k), &                ! Obs Pressure
                      iv%pilot(n)%u(k), &! O-B u
                      iv%pilot(n)%v(k)
         enddo
      end do
      close (ounit)
   end if

   ! [16] ssmi_tb obs:

   if (iv%info(ssmi_tb)%plocal(iv%time) - iv%info(ssmi_tb)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.ssmi_tb',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",482, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'ssmi_tb', iv%info(ssmi_tb)%plocal(iv%time) - &
                                       iv%info(ssmi_tb)%plocal(iv%time-1)
      do n = iv%info(ssmi_tb)%plocal(iv%time-1) + 1, &
             iv%info(ssmi_tb)%plocal(iv%time)
         write(ounit,'(i8,2E22.13)')&
                 n, &  ! Station
                 iv%info(ssmi_tb)%lat(1,n), &       ! Latitude
                 iv%info(ssmi_tb)%lon(1,n)          ! Longitude
         write(ounit,'(7(E22.13,i8,3E22.13))')&
                 iv%ssmi_tb(n)%tb19h, & ! O-B Tb19h
                 iv%ssmi_tb(n)%tb19v, & ! O-B Tb19v
                 iv%ssmi_tb(n)%tb22v, & ! O-B Tb22v
                 iv%ssmi_tb(n)%tb37h, & ! O-B Tb37h
                 iv%ssmi_tb(n)%tb37v, & ! O-B Tb37v
                 iv%ssmi_tb(n)%tb85h, & ! O-B Tb85h
                 iv%ssmi_tb(n)%tb85v    ! O-B Tb85v
      end do
      close (ounit)
   end if

   ! [17] satem obs:

   if (iv%info(satem)%plocal(iv%time) - iv%info(satem)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.satem',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",513, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'satem', iv%info(satem)%plocal(iv%time) - &
                                       iv%info(satem)%plocal(iv%time-1)
      do n = iv%info(satem)%plocal(iv%time-1) + 1, &
             iv%info(satem)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(satem)%levels(n), iv%info(satem)%id(n), &  ! Station
                 iv%info(satem)%lat(1,n), &       ! Latitude
                 iv%info(satem)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(satem)%levels(n)
            write(ounit,'(E22.13,(E22.13,i8,3E22.13))')&
                 iv%satem(n)%p(k), &             ! Obs Pressure
                 iv%satem(n)%thickness(k)
         enddo
      end do
      close (ounit)
   end if

   ! [18] ssmt1 obs:

   if (iv%info(ssmt1)%plocal(iv%time) - iv%info(ssmt1)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.ssmt1',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",541, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'ssmt1', iv%info(ssmt1)%plocal(iv%time) - &
                                       iv%info(ssmt1)%plocal(iv%time-1)
      do n = iv%info(ssmt1)%plocal(iv%time-1) + 1, &
             iv%info(ssmt1)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(ssmt1)%levels(n), iv%info(ssmt1)%id(n), &  ! Station
                 iv%info(ssmt1)%lat(1,n), &       ! Latitude
                 iv%info(ssmt1)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(ssmt1)%levels(n)
            write(ounit,'(E22.13,(E22.13,i8,3E22.13))')&
                  iv%ssmt1(n)%h(k), &             ! Obs height
                  iv%ssmt1(n)%t(k)
         enddo
      end do
      close (ounit)
   end if

   ! [19] ssmt2 obs:

   if (iv%info(ssmt2)%plocal(iv%time) - iv%info(ssmt2)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.ssmt2',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",569, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'ssmt2', iv%info(ssmt2)%plocal(iv%time) - &
                                       iv%info(ssmt2)%plocal(iv%time-1)
      do n = iv%info(ssmt2)%plocal(iv%time-1) + 1, &
             iv%info(ssmt2)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(ssmt2)%levels(n), iv%info(ssmt2)%id(n), &  ! Station
                 iv%info(ssmt2)%lat(1,n), &       ! Latitude
                 iv%info(ssmt2)%lon(1,n)           ! Longitude
         do k = 1 , iv%info(ssmt2)%levels(n)
            write(ounit,'(E22.13,(E22.13,i8,3E22.13))')&
                  iv%ssmt2(n)%h(k), &             ! Obs height
                  iv%ssmt2(n)%rh(k)
         enddo
      end do
      close (ounit)
   end if

   ! [20] scatterometer obs:

   if (iv%info(qscat)%plocal(iv%time) - iv%info(qscat)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.qscat',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",597, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'qscat', iv%info(qscat)%plocal(iv%time) - &
                                       iv%info(qscat)%plocal(iv%time-1)
      do n = iv%info(qscat)%plocal(iv%time-1) + 1, &
             iv%info(qscat)%plocal(iv%time)
         write(ounit,'(i8,a5,2E22.13)')&
                 n, iv%info(qscat)%id(n), &  ! Station
                 iv%info(qscat)%lat(1,n), &       ! Latitude
                 iv%info(qscat)%lon(1,n)          ! Longitude
         write(ounit,'(E22.13,2(E22.13,i8,3E22.13))')&
                   iv%qscat(n)%h, &                ! Obs height
                   iv%qscat(n)%u, &! O-B u
                   iv%qscat(n)%v   ! O-B v
      end do
      close (ounit)
   end if

   ! [21] profiler obs:

   if (iv%info(profiler)%plocal(iv%time) - iv%info(profiler)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.profiler',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",624, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'profiler', iv%info(profiler)%plocal(iv%time) - &
                                       iv%info(profiler)%plocal(iv%time-1)
      do n = iv%info(profiler)%plocal(iv%time-1) + 1, &
             iv%info(profiler)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(profiler)%levels(n), iv%info(profiler)%id(n), &  ! Station
                 iv%info(profiler)%lat(1,n), &       ! Latitude
                 iv%info(profiler)%lon(1,n)          ! Longitude
         do k = 1 , iv%info(profiler)%levels(n)
            write(ounit,'(E22.13,2(E22.13,i8,3E22.13))')&
                     iv%profiler(n)%p(k), &             ! Obs Pressure
                     iv%profiler(n)%u(k), &! O-B u
                     iv%profiler(n)%v(k) ! O-B v
         enddo
      end do
      close (ounit)
   end if

   ! [22] TC bogus obs:

   if (iv%info(bogus)%plocal(iv%time) - iv%info(bogus)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.bogus',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",653, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'bogus', iv%info(bogus)%plocal(iv%time) - &
                                       iv%info(bogus)%plocal(iv%time-1)
      do n = iv%info(bogus)%plocal(iv%time-1) + 1, &
             iv%info(bogus)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(bogus)%levels(n), iv%info(bogus)%id(n), &  ! Station
                 iv%info(bogus)%lat(1,n), &     ! Latitude
                 iv%info(bogus)%lon(1,n)        ! Longitude
         write(ounit,'(E22.13,i8,3E22.13)')&
                 iv%bogus(n)%slp    ! O-B p
         do k = 1 , iv%info(bogus)%levels(n)
            write(ounit,'(2E22.13,4(E22.13,i8,3E22.13))')&
                     iv%bogus(n)%h(k), &
                     iv%bogus(n)%p(k), &             ! Obs Pressure
                     iv%bogus(n)%u(k), &! O-B u
                     iv%bogus(n)%v(k), &! O-B v
                     iv%bogus(n)%t(k), &! O-B t
                     iv%bogus(n)%q(k)   ! O-B q
         enddo
      end do
      close (ounit)
   end if

   ! [23] AIRS retrievals:

   if (iv%info(airsr)%plocal(iv%time) - iv%info(airsr)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.airsr',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",687, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'airsr', iv%info(airsr)%plocal(iv%time) - &
                                       iv%info(airsr)%plocal(iv%time-1)
      do n = iv%info(airsr)%plocal(iv%time-1) + 1, &
             iv%info(airsr)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(airsr)%levels(n), iv%info(airsr)%id(n), &  ! Station
                 iv%info(airsr)%lat(1,n), &    ! Latitude
                 iv%info(airsr)%lon(1,n)       ! Longitude
         do k = 1 , iv%info(airsr)%levels(n)
            write(ounit,'(E22.13,2(E22.13,i8,3E22.13))')&
                     iv%airsr(n)%p(k), &             ! Obs Pressure
                     iv%airsr(n)%t(k), &! O-B t
                     iv%airsr(n)%q(k)   ! O-B q
         enddo
      end do
      close (ounit)
   end if

   ! [24] gpsref obs:

   if (iv%info(gpsref)%plocal(iv%time) - iv%info(gpsref)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.gpsref',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",716, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'gpsref', iv%info(gpsref)%plocal(iv%time) - &
                                       iv%info(gpsref)%plocal(iv%time-1)
      do n = iv%info(gpsref)%plocal(iv%time-1) + 1, &
             iv%info(gpsref)%plocal(iv%time)
         write(ounit,'(2i8,a5,2E22.13)')&
                 n, iv%info(gpsref)%levels(n), iv%info(gpsref)%id(n), &  ! Station
                 iv%info(gpsref)%lat(1,n), &    ! Latitude
                 iv%info(gpsref)%lon(1,n)       ! Longitude
         do k = 1 , iv%info(gpsref)%levels(n)
            write(ounit,'(E22.13,(E22.13,i8,3E22.13))')&
                     iv%gpsref(n)%h(k), &             ! Obs Height
                     iv%gpsref(n)%ref(k) ! O-B ref
         enddo
      end do
      close (ounit)
   end if

   ! [25] radar obs:

   if (iv%info(radar)%plocal(iv%time) - iv%info(radar)%plocal(iv%time-1) > 0) then

      open (unit=ounit,file=trim(filename)//'.radar',form='formatted',status='replace', &
         iostat=ios)
      if (ios /= 0) then
         call da_error("da_write_iv_for_multi_inc.inc",744, &
            (/"Cannot open conventional observation omb file"//filename/))
      end if

      write(ounit,'(a20,i8)')'radar', iv%info(radar)%plocal(iv%time) - &
                                       iv%info(radar)%plocal(iv%time-1)
      do n = iv%info(radar)%plocal(iv%time-1) + 1, &
             iv%info(radar)%plocal(iv%time)
         write(ounit,'(2i8,2E22.13)')&
                 n, iv%info(radar)%levels(n),  &
                 iv%info(radar)%lat(1,n), &      ! Latitude
                 iv%info(radar)%lon(1,n)         ! Longitude
         do k = 1 , iv%info(radar)%levels(n)
            write(ounit,'(E22.13,i8,3E22.13)')&
                     iv%radar(n)%rv(k)           ! radar_rv

         enddo
      end do
      close (ounit)
   end if



   !-------------------------------------------------------------------------------


   call da_free_unit(ounit)

   if (trace_use) call da_trace_exit("da_write_iv_for_multi_inc")

end subroutine da_write_iv_for_multi_inc


subroutine da_read_iv_for_multi_inc(file_index, iv)

   !-----------------------------------------------------------------------
   ! Purpose: Read for Multi-incremental 
   !-----------------------------------------------------------------------

   !-------------------------------------------------------------------------
   ! read iv=O-B structure written by WRFVAR
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout)    :: iv      ! O-B structure.
   integer,  intent(in)             :: file_index
   integer                      :: unit_in
   character(len=filename_len)   :: filename

   integer      :: num_obs, ios 
   character*20 :: ob_type_string               
   
   integer      :: n, gn
   logical      :: found_flag

   if (trace_use) call da_trace_entry("da_read_iv_for_multi_inc")

   !-------------------------------------------------------------------------
   ! Fix input unit
   !-------------------------------------------------------------------------

   call da_get_unit(unit_in)

   write(unit=filename, fmt='(a,i3.3)') 'gts_omb.', file_index

   ! [1] surface obs:

   if (iv%info(synop)%plocal(iv%time)-iv%info(synop)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.synop',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",40, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'synop' ) &
           call da_error("da_read_iv_for_multi_inc.inc",46, &
                           (/"Cannot find synop marker. "/))
       gn = 0
       do n = iv%info(synop)%plocal(iv%time-1) + 1, &
              iv%info(synop)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",53, &
                           (/"Cannot find synop obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(synop)%plocal(iv%time)-iv%info(synop)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",58, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [2] metar obs:

   if (iv%info(metar)%plocal(iv%time)-iv%info(metar)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.metar',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",69, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'metar' ) &
           call da_error("da_read_iv_for_multi_inc.inc",75, &
                           (/"Cannot find metar marker. "/))
       gn = 0
       do n = iv%info(metar)%plocal(iv%time-1) + 1, &
              iv%info(metar)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",82, &
                           (/"Cannot find metar obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(metar)%plocal(iv%time)-iv%info(metar)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",87, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [3] ships obs:

   if (iv%info(ships)%plocal(iv%time)-iv%info(ships)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.ships',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",98, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'ships' ) &
           call da_error("da_read_iv_for_multi_inc.inc",104, &
                           (/"Cannot find ships marker. "/))
       gn = 0
       do n = iv%info(ships)%plocal(iv%time-1) + 1, &
              iv%info(ships)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",111, &
                           (/"Cannot find ships obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(ships)%plocal(iv%time)-iv%info(ships)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",116, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [4] sonde_sfc obs:

   if (iv%info(sonde_sfc)%plocal(iv%time)-iv%info(sonde_sfc)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.sonde_sfc',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",127, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'sonde_sfc' ) &
           call da_error("da_read_iv_for_multi_inc.inc",133, &
                           (/"Cannot find sonde_sfc marker. "/))
       gn = 0
       do n = iv%info(sonde_sfc)%plocal(iv%time-1) + 1, &
              iv%info(sonde_sfc)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",140, &
                           (/"Cannot find sonde_sfc obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(sonde_sfc)%plocal(iv%time)-iv%info(sonde_sfc)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",145, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [5] sound obs:

   if (iv%info(sound)%plocal(iv%time)-iv%info(sound)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.sound',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",156, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'sound' ) &
           call da_error("da_read_iv_for_multi_inc.inc",162, &
                           (/"Cannot find sound marker. "/))
       gn = 0
       do n = iv%info(sound)%plocal(iv%time-1) + 1, &
              iv%info(sound)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",169, &
                           (/"Cannot find sound obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(sound)%plocal(iv%time)-iv%info(sound)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",174, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [6] mtgirs obs:

   if (iv%info(mtgirs)%plocal(iv%time)-iv%info(mtgirs)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.mtgirs',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",185, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'mtgirs' ) &
           call da_error("da_read_iv_for_multi_inc.inc",191, &
                           (/"Cannot find mtgirs marker. "/))
       gn = 0
       do n = iv%info(mtgirs)%plocal(iv%time-1) + 1, &
              iv%info(mtgirs)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",198, &
                           (/"Cannot find mtgirs obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(mtgirs)%plocal(iv%time)-iv%info(mtgirs)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",203, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [7] tamdar obs:

   if (iv%info(tamdar)%plocal(iv%time)-iv%info(tamdar)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.tamdar',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",214, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'tamdar' ) &
           call da_error("da_read_iv_for_multi_inc.inc",220, &
                           (/"Cannot find tamdar marker. "/))
       gn = 0
       do n = iv%info(tamdar)%plocal(iv%time-1) + 1, &
              iv%info(tamdar)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",227, &
                           (/"Cannot find tamdar obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(tamdar)%plocal(iv%time)-iv%info(tamdar)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",232, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [8] tamdar_sfc obs:

   if (iv%info(tamdar_sfc)%plocal(iv%time)-iv%info(tamdar_sfc)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.tamdar_sfc',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",243, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'tamdar_sfc' ) &
           call da_error("da_read_iv_for_multi_inc.inc",249, &
                           (/"Cannot find tamdar_sfc marker. "/))
       gn = 0
       do n = iv%info(tamdar_sfc)%plocal(iv%time-1) + 1, &
              iv%info(tamdar_sfc)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",256, &
                           (/"Cannot find tamdar_sfc obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(tamdar_sfc)%plocal(iv%time)-iv%info(tamdar_sfc)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",261, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [9] buoy obs:

   if (iv%info(buoy)%plocal(iv%time)-iv%info(buoy)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.buoy',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",272, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'buoy' ) &
           call da_error("da_read_iv_for_multi_inc.inc",278, &
                           (/"Cannot find buoy marker. "/))
       gn = 0
       do n = iv%info(buoy)%plocal(iv%time-1) + 1, &
              iv%info(buoy)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",285, &
                           (/"Cannot find buoy obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(buoy)%plocal(iv%time)-iv%info(buoy)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",290, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [10] Geo AMV obs:

   if (iv%info(geoamv)%plocal(iv%time)-iv%info(geoamv)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.geoamv',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",301, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'geoamv' ) &
           call da_error("da_read_iv_for_multi_inc.inc",307, &
                           (/"Cannot find geoamv marker. "/))
       gn = 0
       do n = iv%info(geoamv)%plocal(iv%time-1) + 1, &
              iv%info(geoamv)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",314, &
                           (/"Cannot find geoamv obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(geoamv)%plocal(iv%time)-iv%info(geoamv)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",319, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [11] gpspw obs:

   if (iv%info(gpspw)%plocal(iv%time)-iv%info(gpspw)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.gpspw',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",330, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'gpspw' ) &
           call da_error("da_read_iv_for_multi_inc.inc",336, &
                           (/"Cannot find gpspw marker. "/))
       gn = 0
       do n = iv%info(gpspw)%plocal(iv%time-1) + 1, &
              iv%info(gpspw)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",343, &
                           (/"Cannot find gpspw obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(gpspw)%plocal(iv%time)-iv%info(gpspw)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",348, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [12] SSM/I obs:

   if (iv%info(ssmi_rv)%plocal(iv%time)-iv%info(ssmi_rv)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.ssmir',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",359, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'ssmir' ) &
           call da_error("da_read_iv_for_multi_inc.inc",365, &
                           (/"Cannot find ssmir marker. "/))
       gn = 0
       do n = iv%info(ssmi_rv)%plocal(iv%time-1) + 1, &
              iv%info(ssmi_rv)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",372, &
                           (/"Cannot find ssmir obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(ssmi_rv)%plocal(iv%time)-iv%info(ssmi_rv)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",377, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [13] airep obs:

   if (iv%info(airep)%plocal(iv%time)-iv%info(airep)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.airep',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",388, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'airep' ) &
           call da_error("da_read_iv_for_multi_inc.inc",394, &
                           (/"Cannot find airep marker. "/))
       gn = 0
       do n = iv%info(airep)%plocal(iv%time-1) + 1, &
              iv%info(airep)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",401, &
                           (/"Cannot find airep obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(airep)%plocal(iv%time)-iv%info(airep)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",406, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [14] polaramv obs:

   if (iv%info(polaramv)%plocal(iv%time)-iv%info(polaramv)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.polaramv',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",417, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'polaramv' ) &
           call da_error("da_read_iv_for_multi_inc.inc",423, &
                           (/"Cannot find polaramv marker. "/))
       gn = 0
       do n = iv%info(polaramv)%plocal(iv%time-1) + 1, &
              iv%info(polaramv)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",430, &
                           (/"Cannot find polaramv obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(polaramv)%plocal(iv%time)-iv%info(polaramv)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",435, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [15] pilot obs:

   if (iv%info(pilot)%plocal(iv%time)-iv%info(pilot)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.pilot',form='formatted',status='old',iostat=ios)

       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",447, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'pilot' ) &
           call da_error("da_read_iv_for_multi_inc.inc",453, &
                           (/"Cannot find pilot marker. "/))
       gn = 0
       do n = iv%info(pilot)%plocal(iv%time-1) + 1, &
              iv%info(pilot)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",460, &
                           (/"Cannot find pilot obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(pilot)%plocal(iv%time)-iv%info(pilot)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",465, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [16] ssmi_tb obs:

   if (iv%info(ssmi_tb)%plocal(iv%time)-iv%info(ssmi_tb)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.ssmi_tb',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",476, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'ssmi_tb' ) &
           call da_error("da_read_iv_for_multi_inc.inc",482, &
                           (/"Cannot find ssmi_tb marker. "/))
       gn = 0
       do n = iv%info(ssmi_tb)%plocal(iv%time-1) + 1, &
              iv%info(ssmi_tb)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",489, &
                           (/"Cannot find ssmi_tb obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(ssmi_tb)%plocal(iv%time)-iv%info(ssmi_tb)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",494, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [17] satem obs:

   if (iv%info(satem)%plocal(iv%time)-iv%info(satem)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.satem',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",505, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'satem' ) &
           call da_error("da_read_iv_for_multi_inc.inc",511, &
                           (/"Cannot find satem marker. "/))
       gn = 0
       do n = iv%info(satem)%plocal(iv%time-1) + 1, &
              iv%info(satem)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",518, &
                           (/"Cannot find satem obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(satem)%plocal(iv%time)-iv%info(satem)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",523, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [18] ssmt1 obs:

   if (iv%info(ssmt1)%plocal(iv%time)-iv%info(ssmt1)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.ssmt1',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",534, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'ssmt1' ) &
           call da_error("da_read_iv_for_multi_inc.inc",540, &
                           (/"Cannot find ssmt1 marker. "/))
       gn = 0
       do n = iv%info(ssmt1)%plocal(iv%time-1) + 1, &
              iv%info(ssmt1)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",547, &
                           (/"Cannot find ssmt1 obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(ssmt1)%plocal(iv%time)-iv%info(ssmt1)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",552, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [19] ssmt2 obs:

   if (iv%info(ssmt2)%plocal(iv%time)-iv%info(ssmt2)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.ssmt2',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",563, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'ssmt2' ) &
           call da_error("da_read_iv_for_multi_inc.inc",569, &
                           (/"Cannot find ssmt2 marker. "/))
       gn = 0
       do n = iv%info(ssmt2)%plocal(iv%time-1) + 1, &
              iv%info(ssmt2)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",576, &
                           (/"Cannot find ssmt2 obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(ssmt2)%plocal(iv%time)-iv%info(ssmt2)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",581, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [20] scatterometer obs:

   if (iv%info(qscat)%plocal(iv%time)-iv%info(qscat)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.qscat',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",592, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'qscat' ) &
           call da_error("da_read_iv_for_multi_inc.inc",598, &
                           (/"Cannot find qscat marker. "/))
       gn = 0
       do n = iv%info(qscat)%plocal(iv%time-1) + 1, &
              iv%info(qscat)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",605, &
                           (/"Cannot find qscat obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(qscat)%plocal(iv%time)-iv%info(qscat)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",610, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [21] profiler obs:

   if (iv%info(profiler)%plocal(iv%time)-iv%info(profiler)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.profiler',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",621, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'profiler' ) &
           call da_error("da_read_iv_for_multi_inc.inc",627, &
                           (/"Cannot find profiler marker. "/))
       gn = 0
       do n = iv%info(profiler)%plocal(iv%time-1) + 1, &
              iv%info(profiler)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",634, &
                           (/"Cannot find profiler obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(profiler)%plocal(iv%time)-iv%info(profiler)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",639, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [22] TC bogus obs:

   if (iv%info(bogus)%plocal(iv%time)-iv%info(bogus)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.bogus',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",650, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'bogus' ) &
           call da_error("da_read_iv_for_multi_inc.inc",656, &
                           (/"Cannot find bogus marker. "/))
       gn = 0
       do n = iv%info(bogus)%plocal(iv%time-1) + 1, &
              iv%info(bogus)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",663, &
                           (/"Cannot find bogus obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(bogus)%plocal(iv%time)-iv%info(bogus)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",668, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [23] AIRS retrievals:

   if (iv%info(airsr)%plocal(iv%time)-iv%info(airsr)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.airsr',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",679, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'airsr' ) &
           call da_error("da_read_iv_for_multi_inc.inc",685, &
                           (/"Cannot find airsr marker. "/))
       gn = 0
       do n = iv%info(airsr)%plocal(iv%time-1) + 1, &
              iv%info(airsr)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",692, &
                           (/"Cannot find airsr obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(airsr)%plocal(iv%time)-iv%info(airsr)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",697, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if

   ! [24] gpsref obs:

   if (iv%info(gpsref)%plocal(iv%time)-iv%info(gpsref)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.gpsref',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",708, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'gpsref' ) &
           call da_error("da_read_iv_for_multi_inc.inc",714, &
                           (/"Cannot find gpsref marker. "/))
       gn = 0
       do n = iv%info(gpsref)%plocal(iv%time-1) + 1, &
              iv%info(gpsref)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",721, &
                           (/"Cannot find gpsref obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(gpsref)%plocal(iv%time)-iv%info(gpsref)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",726, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if


   ! [25] radar obs:

   if (iv%info(radar)%plocal(iv%time)-iv%info(radar)%plocal(iv%time-1) > 0) then

       open(unit=unit_in,file=trim(filename)//'.radar',form='formatted',status='old',iostat=ios)
       if (ios /= 0) Then
          call da_error("da_read_iv_for_multi_inc.inc",738, &
             (/"Cannot open file"//filename/))
       end if

       read(unit_in,'(a20,i8)', end = 999, err = 1000) ob_type_string,num_obs
       if ( trim(adjustl(ob_type_string)) .ne. 'radar' ) &
           call da_error("da_read_iv_for_multi_inc.inc",744, &
                           (/"Cannot find radar marker. "/))
       gn = 0
       do n = iv%info(radar)%plocal(iv%time-1) + 1, &
              iv%info(radar)%plocal(iv%time)
          call da_search_obs (ob_type_string, unit_in, num_obs, n, iv, found_flag)
          if (found_flag .eqv. .false.) &
              call da_error("da_read_iv_for_multi_inc.inc",751, &
                           (/"Cannot find radar obs. "/))
          gn = gn + 1
       end do
       if (gn /= iv%info(radar)%plocal(iv%time)-iv%info(radar)%plocal(iv%time-1)) &
           call da_error("da_read_iv_for_multi_inc.inc",756, &
                        (/"Unequal obs. found "/))
       close (unit_in)
   end if


999 continue
   close (unit_in)
   call da_free_unit(unit_in)

   if (trace_use) call da_trace_exit("da_read_iv_for_multi_inc")
   return

1000 continue
   write(unit=message(1), fmt='(a,i3)') &
      'read error on unit: ',unit_in
   call da_warning("da_read_iv_for_multi_inc.inc",772,message(1:1))

end subroutine da_read_iv_for_multi_inc
subroutine da_search_obs (ob_type_string, unit_in, num_obs, nth, iv, found_flag)

   !-----------------------------------------------------------------------
   ! Purpose: Search obs. in gts_omb.000 
   !-----------------------------------------------------------------------

   !-------------------------------------------------------------------------
   ! read iv=O-B structure written by WRFVAR
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout)    :: iv      ! O-B structure.
   integer, intent(in)              :: unit_in, nth, num_obs
   character(len=20), intent(in)    :: ob_type_string
   logical, intent(out)             :: found_flag

   character*5  :: stn_id
   real         :: lat, lon
   integer      :: n, n_dummy, k, levels
   real, parameter :: MIN_ERR=1.0E-6

   if (trace_use) call da_trace_entry("da_search_obs")

   found_flag = .true. 
   SELECT CASE (trim(adjustl(ob_type_string)))

   CASE ('synop')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(synop)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(synop)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(synop)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,5(E22.13,i8,3E22.13))')&
              iv%synop(nth)%h, &
              iv%synop(nth)%u, &!  O-B u
              iv%synop(nth)%v, &!  O-B v
              iv%synop(nth)%t, &!  O-B t
              iv%synop(nth)%p, &!  O-B p
              iv%synop(nth)%q  !  O-B q
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('metar')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(metar)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(metar)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(metar)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,5(E22.13,i8,3E22.13))')&
              iv%metar(nth)%h, &
              iv%metar(nth)%u, &!  O-B u
              iv%metar(nth)%v, &!  O-B v
              iv%metar(nth)%t, &!  O-B t
              iv%metar(nth)%p, &!  O-B p
              iv%metar(nth)%q  !  O-B q
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('ships')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(ships)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(ships)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(ships)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,5(E22.13,i8,3E22.13))')&
              iv%ships(nth)%h, &
              iv%ships(nth)%u, &!  O-B u
              iv%ships(nth)%v, &!  O-B v
              iv%ships(nth)%t, &!  O-B t
              iv%ships(nth)%p, &!  O-B p
              iv%ships(nth)%q  !  O-B q
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('sonde_sfc')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(sonde_sfc)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(sonde_sfc)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(sonde_sfc)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,5(E22.13,i8,3E22.13))')&
              iv%sonde_sfc(nth)%h, &
              iv%sonde_sfc(nth)%u, &!  O-B u
              iv%sonde_sfc(nth)%v, &!  O-B v
              iv%sonde_sfc(nth)%t, &!  O-B t
              iv%sonde_sfc(nth)%p, &!  O-B p
              iv%sonde_sfc(nth)%q  !  O-B q
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('sound')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(sound)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(sound)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(sound)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(2E22.13,4(E22.13,i8,3E22.13))')&
                 iv % sound(nth) % h(k), &
                 iv % sound(nth) % p(k), &             ! Obs Pressure
                 iv%sound(nth)%u(k), &! O-B u
                 iv%sound(nth)%v(k), &! O-B v
                 iv%sound(nth)%t(k), &! O-B t
                 iv%sound(nth)%q(k)   ! O-B q
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('mtgirs')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(mtgirs)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(mtgirs)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(mtgirs)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(2E22.13,4(E22.13,i8,3E22.13))')&
                 iv % mtgirs(nth) % h(k), &
                 iv % mtgirs(nth) % p(k), &             ! Obs Pressure
                 iv%mtgirs(nth)%u(k), &! O-B u
                 iv%mtgirs(nth)%v(k), &! O-B v
                 iv%mtgirs(nth)%t(k), &! O-B t
                 iv%mtgirs(nth)%q(k)   ! O-B q
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('tamdar')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(tamdar)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(tamdar)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(tamdar)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(2E22.13,4(E22.13,i8,3E22.13))')&
                 iv % tamdar(nth) % h(k), &
                 iv % tamdar(nth) % p(k), &             ! Obs Pressure
                 iv%tamdar(nth)%u(k), &! O-B u
                 iv%tamdar(nth)%v(k), &! O-B v
                 iv%tamdar(nth)%t(k), &! O-B t
                 iv%tamdar(nth)%q(k)   ! O-B q
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('tamdar_sfc')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(tamdar_sfc)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(tamdar_sfc)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(tamdar_sfc)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,5(E22.13,i8,3E22.13))')&
              iv%tamdar_sfc(nth)%h, &
              iv%tamdar_sfc(nth)%u, &!  O-B u
              iv%tamdar_sfc(nth)%v, &!  O-B v
              iv%tamdar_sfc(nth)%t, &!  O-B t
              iv%tamdar_sfc(nth)%p, &!  O-B p
              iv%tamdar_sfc(nth)%q  !  O-B q
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('buoy')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(buoy)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(buoy)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(buoy)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,5(E22.13,i8,3E22.13))')&
              iv%buoy(nth)%h, &
              iv%buoy(nth)%u, &!  O-B u
              iv%buoy(nth)%v, &!  O-B v
              iv%buoy(nth)%t, &!  O-B t
              iv%buoy(nth)%p, &!  O-B p
              iv%buoy(nth)%q  !  O-B q
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('geoamv')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(geoamv)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(geoamv)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(geoamv)%lon(1,nth) - lon ) < MIN_ERR ) then
         do k = 1, levels
            read(unit_in,'(E22.13,2(E22.13,i8,3E22.13))')&
                 iv % geoamv(nth) % p(k), &                ! Obs Pressure
                 iv%geoamv(nth)%u(k), &! O-B u
                 iv%geoamv(nth)%v(k)
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('gpspw')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(gpspw)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(gpspw)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(gpspw)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,i8,3E22.13)')&
              iv%gpspw(nth)%tpw
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('radar')

   do n = 1, num_obs
      read(unit_in,'(2i8,2E22.13)') n_dummy, levels, lat, lon

      if ( abs(iv%info(radar)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(radar)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,i8,3E22.13)')&
                 iv%radar(nth)%rv(k) 
         enddo

         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('ssmir')

   do n = 1, num_obs
      read(unit_in,'(i8,2E22.13)') n_dummy, lat, lon
      if ( abs(iv%info(ssmi_rv)%lat(1,nth) - lat ) < MIN_ERR .and. &
           abs(iv%info(ssmi_rv)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(2(E22.13,i8,3E22.13))')&
              iv%ssmi_rv(nth)%speed, & ! O-B speed
              iv%ssmi_rv(nth)%tpw ! O-BA tpw
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('airep')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
 
      if ( trim(iv%info(airep)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(airep)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(airep)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(2E22.13,4(E22.13,i8,3E22.13))')&
                 iv % airep(nth) % h(k), &
                 iv % airep(nth) % p(k), &             ! Obs pressure
                 iv%airep(nth)%u(k), &! O-B u
                 iv%airep(nth)%v(k), &! O-B v
                 iv%airep(nth)%t(k), &! 
                 iv%airep(nth)%q(k)
         enddo

         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('polaramv')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(polaramv)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(polaramv)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(polaramv)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,2(E22.13,i8,3E22.13))')&
                 iv % polaramv(nth) % p(k), &                ! Obs Pressure
                 iv%polaramv(nth)%u(k), &! O-B u
                 iv%polaramv(nth)%v(k)
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('pilot')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(pilot)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(pilot)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(pilot)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,2(E22.13,i8,3E22.13))')&
                 iv % pilot(nth) % p(k), &                ! Obs Pressure
                 iv%pilot(nth)%u(k), &! O-B u
                 iv%pilot(nth)%v(k)
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('ssmi_tb')

   do n = 1, num_obs
      read(unit_in,'(i8,2E22.13)') n_dummy, lat, lon
      if ( abs(iv%info(ssmi_tb)%lat(1,nth) - lat ) < MIN_ERR .and. &
           abs(iv%info(ssmi_tb)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(7(E22.13,i8,3E22.13))')&
              iv%ssmi_tb(nth)%tb19h, & ! O-B Tb19h
              iv%ssmi_tb(nth)%tb19v, & ! O-B Tb19v
              iv%ssmi_tb(nth)%tb22v, & ! O-B Tb22v
              iv%ssmi_tb(nth)%tb37h, & ! O-B Tb37h
              iv%ssmi_tb(nth)%tb37v, & ! O-B Tb37v
              iv%ssmi_tb(nth)%tb85h, & ! O-B Tb85h
              iv%ssmi_tb(nth)%tb85v    ! O-B Tb85v
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('satem')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(satem)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(satem)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(satem)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,(E22.13,i8,3E22.13))')&
                 iv % satem(nth) % p(k), &             ! Obs Pressure
                 iv%satem(nth)%thickness(k)
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('ssmt1')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(ssmt1)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(ssmt1)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(ssmt1)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,(E22.13,i8,3E22.13))')&
                 iv % ssmt1(nth) % h(k), &             ! Obs Pressure
                 iv%ssmt1(nth)%t(k)
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('ssmt2')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(ssmt2)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(ssmt2)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(ssmt2)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,(E22.13,i8,3E22.13))')&
                 iv % ssmt2(nth) % h(k), &             ! Obs Pressure
                 iv%ssmt2(nth)%rh(k)
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('qscat')

   do n = 1, num_obs
      read(unit_in,'(i8,a5,2E22.13)') n_dummy, stn_id, lat, lon
      if ( trim(iv%info(qscat)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(qscat)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(qscat)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,2(E22.13,i8,3E22.13))')&
              iv % qscat(nth) % h, &                ! Obs height
              iv%qscat(nth)%u, &! O-B u
              iv%qscat(nth)%v   ! O-B v
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('profiler')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(profiler)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(profiler)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(profiler)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,2(E22.13,i8,3E22.13))')&
                 iv % profiler(nth) % p(k), &             ! Obs Pressure
                 iv%profiler(nth)%u(k), &! O-B u
                 iv%profiler(nth)%v(k) ! O-B v
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('bogus')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(bogus)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(bogus)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(bogus)%lon(1,nth) - lon ) < MIN_ERR ) then

         read(unit_in,'(E22.13,i8,3E22.13)') iv%bogus(nth)%slp
         do k = 1, levels
            read(unit_in,'(2E22.13,4(E22.13,i8,3E22.13))')&
                 iv % bogus(nth) % h(k), &
                 iv % bogus(nth) % p(k), &             ! Obs Pressure
                 iv%bogus(nth)%u(k), &! O-B u
                 iv%bogus(nth)%v(k), &! O-B v
                 iv%bogus(nth)%t(k), &! O-B t
                 iv%bogus(nth)%q(k)   ! O-B q
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         read(unit_in,*)
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('airsr')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(airsr)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(airsr)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(airsr)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,2(E22.13,i8,3E22.13))')&
                 iv % airsr(nth) % p(k), &             ! Obs Pressure
                 iv%airsr(nth)%t(k), &! O-B t
                 iv%airsr(nth)%q(k)   ! O-B q
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
        do k = 1, levels
          read(unit_in,*)
        enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE ('gpsref')

   do n = 1, num_obs
      read(unit_in,'(2i8,a5,2E22.13)') n_dummy, levels, stn_id, lat, lon
      if ( trim(iv%info(gpsref)%id(nth))  == trim(adjustl(stn_id)) .and.  &
           abs(iv%info(gpsref)%lat(1,nth) - lat ) < MIN_ERR     .and.  &
           abs(iv%info(gpsref)%lon(1,nth) - lon ) < MIN_ERR ) then

         do k = 1, levels
            read(unit_in,'(E22.13,(E22.13,i8,3E22.13))')&
                 iv % gpsref(nth) % h(k), &             ! Obs Height
                 iv%gpsref(nth)%ref(k) ! O-B ref
         enddo
         !found_flag = .true.
         rewind (unit_in)
         read(unit_in,*)
         if (trace_use) call da_trace_exit("da_search_obs")
         return
      else
         do k = 1, levels
            read(unit_in,*)
         enddo
      endif
   enddo
   !found_flag = .false.
   rewind (unit_in)
   read(unit_in,*)

   CASE default;
 
   write(unit=message(1), fmt='(a,a20,a,i3)') &
        'Got unknown obs_type string:', trim(ob_type_string),' on unit ',unit_in
   call da_error("da_search_obs.inc",745,message(1:1))

   END SELECT

   if (trace_use) call da_trace_exit("da_search_obs")
   return
end subroutine da_search_obs
subroutine da_write_obs_etkf(ob, iv, re)

   !-------------------------------------------------------------------------
   ! Purpose: Writes out components of iv=O-B structure.
   !-------------------------------------------------------------------------   
   !
   ! Arthur P. Mizzi (NCAR/MMM) February 2011 Modfied to output the extended ob.etkf file.  

   implicit none

   type (y_type), intent(in)     :: ob      ! Observation structure.
   type (iv_type), intent(in)    :: iv      ! O-B structure.
   type (y_type), intent(inout)  :: re      ! residual vector.
      
   integer                       :: n, k, num_obs, ios
   integer                       :: ounit     ! Output unit           
   character(len=20)             :: filename
   character(len=20)             :: apm_char 
   integer                       :: apm_index, apm_int
   real                          :: apm_plc

   if (trace_use) call da_trace_entry("da_write_obs_etkf")

   !-------------------------------------------------------------------------
   ! Fix output unit
   !-------------------------------------------------------------------------

   apm_index=0
   
   call da_get_unit(ounit)

    write(unit=filename, fmt='(a,i4.4)') 'ob.etkf.', myproc

   open (unit=ounit,file=trim(filename),form='formatted',status='replace', &
      iostat=ios)
   if (ios /= 0) then
      call da_error("da_write_obs_etkf.inc",41, &
         (/"Cannot open ETKF observation file"//filename/))
   end if

   ! [0] Format for extended ob.etkf files (APM 02-10-2011)
 
1000 format(3f17.7,2x,a10,2x,2(f6.0,2x),4(f8.2,2x),i10)

   ! [1] Transfer surface obs:

   if (iv%info(synop)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(synop)%nlocal
         if (iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(synop)%nlocal
            if (iv%info(synop)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(synop)%platform(n)(4:6),*) apm_int 
               if ( iv%synop(n)%u%qc >= 0 .and. ob%synop(n)%u /= missing_r ) then 
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%synop(n)%u, iv%synop(n)%u%inv, iv%synop(n)%u%error, &
                  'WNU', apm_plc, -888.88, iv%info(synop)%lat(1,n), iv%info(synop)%lon(1,n), &
                  (ob%synop(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%synop(n)%v%qc >= 0 .and. ob%synop(n)%v /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%synop(n)%v, iv%synop(n)%v%inv, iv%synop(n)%v%error, &
                  'WNV', apm_plc, -888.88, iv%info(synop)%lat(1,n), iv%info(synop)%lon(1,n), &
                  (ob%synop(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%synop(n)%t%qc >= 0 .and. ob%synop(n)%t /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%synop(n)%t, iv%synop(n)%t%inv, iv%synop(n)%t%error, &
                  'TMP', apm_plc, -888.88, iv%info(synop)%lat(1,n), iv%info(synop)%lon(1,n), &
                  (ob%synop(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%synop(n)%p%qc >= 0 .and. ob%synop(n)%p /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%synop(n)%p, iv%synop(n)%p%inv, iv%synop(n)%p%error, &
                  'PRS', apm_plc, -888.88, iv%info(synop)%lat(1,n), iv%info(synop)%lon(1,n), &
                  (ob%synop(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%synop(n)%q%qc >= 0 .and. ob%synop(n)%q /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%synop(n)%q, iv%synop(n)%q%inv, iv%synop(n)%q%error, &
                  'QVP', apm_plc, -888.88, iv%info(synop)%lat(1,n), iv%info(synop)%lon(1,n), &
                  (ob%synop(n)%p)/100., -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

   ! [2] Transfer metar obs:

   if (iv%info(metar)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(metar)%nlocal
         if (iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(metar)%nlocal
            if (iv%info(metar)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(metar)%platform(n)(4:6),*) apm_int 
               if ( iv%metar(n)%u%qc >= 0 .and. ob%metar(n)%u /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%metar(n)%u, iv%metar(n)%u%inv, iv%metar(n)%u%error, &
                  'WNU', apm_plc, -888.88, iv%info(metar)%lat(1,n), iv%info(metar)%lon(1,n), &
                  (ob%metar(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%metar(n)%v%qc >= 0 .and. ob%metar(n)%v /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%metar(n)%v, iv%metar(n)%v%inv, iv%metar(n)%v%error, &
                  'WNV', apm_plc, -888.88, iv%info(metar)%lat(1,n), iv%info(metar)%lon(1,n), &
                  (ob%metar(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%metar(n)%t%qc >= 0 .and. ob%metar(n)%t /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%metar(n)%t, iv%metar(n)%t%inv, iv%metar(n)%t%error, &
                  'TMP', apm_plc, -888.88, iv%info(metar)%lat(1,n), iv%info(metar)%lon(1,n), &
                  (ob%metar(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%metar(n)%p%qc >= 0 .and. ob%metar(n)%p /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%metar(n)%p, iv%metar(n)%p%inv, iv%metar(n)%p%error, &
                  'PRS', apm_plc, -888.88, iv%info(metar)%lat(1,n), iv%info(metar)%lon(1,n), &
                  (ob%metar(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%metar(n)%q%qc >= 0 .and. ob%metar(n)%q /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%metar(n)%q, iv%metar(n)%q%inv, iv%metar(n)%q%error, &
                  'QVP', apm_plc, -888.88, iv%info(metar)%lat(1,n), iv%info(metar)%lon(1,n), &
                  (ob%metar(n)%p)/100., -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

   ! [3] Transfer ships obs:

   if (iv%info(ships)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ships)%nlocal
         if (iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(ships)%nlocal
            if (iv%info(ships)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(ships)%platform(n)(4:6),*) apm_int 
               if ( iv%ships(n)%u%qc >= 0 .and. ob%ships(n)%u /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ships(n)%u, iv%ships(n)%u%inv, iv%ships(n)%u%error, &
                  'WNU', apm_plc, -888.88, iv%info(ships)%lat(1,n), iv%info(ships)%lon(1,n), &
                  (ob%ships(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%ships(n)%v%qc >= 0 .and. ob%ships(n)%v /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ships(n)%v, iv%ships(n)%v%inv, iv%ships(n)%v%error, &
                  'WNV', apm_plc, -888.88, iv%info(ships)%lat(1,n), iv%info(ships)%lon(1,n), &
                  (ob%ships(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%ships(n)%t%qc >= 0 .and. ob%ships(n)%t /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ships(n)%t, iv%ships(n)%t%inv, iv%ships(n)%t%error, &
                  'TMP', apm_plc, -888.88, iv%info(ships)%lat(1,n), iv%info(ships)%lon(1,n), &
                  (ob%ships(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%ships(n)%p%qc >= 0 .and. ob%ships(n)%p /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ships(n)%p, iv%ships(n)%p%inv, iv%ships(n)%p%error, &
                  'PRS', apm_plc, -888.88, iv%info(ships)%lat(1,n), iv%info(ships)%lon(1,n), &
                  (ob%ships(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%ships(n)%q%qc >= 0 .and. ob%ships(n)%q /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ships(n)%q, iv%ships(n)%q%inv, iv%ships(n)%q%error, &
                  'QVP',apm_plc, -888.88, iv%info(ships)%lat(1,n), iv%info(ships)%lon(1,n), &
                  (ob%ships(n)%p)/100., -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

  ! [4.1] Transfer Geo AMVs Obs:

   if (iv%info(geoamv)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(geoamv)%nlocal
        if (iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(geoamv)%nlocal
            if (iv%info(geoamv)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(geoamv)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(geoamv)%levels(n)
                  if ( iv%geoamv(n)%u(k)%qc >= 0 .and. ob%geoamv(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%geoamv(n)%u(k), iv%geoamv(n)%u(k)%inv, iv%geoamv(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(geoamv)%lat(1,n), iv%info(geoamv)%lon(1,n), &
                     (iv%geoamv(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%geoamv(n)%v(k)%qc >= 0 .and. ob%geoamv(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%geoamv(n)%v(k), iv%geoamv(n)%v(k)%inv, iv%geoamv(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(geoamv)%lat(1,n), iv%info(geoamv)%lon(1,n), &
                     (iv%geoamv(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

  ! [4.2] Transfer Polar AMVs Obs:

   if (iv%info(polaramv)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(polaramv)%nlocal
        if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(polaramv)%nlocal
            if (iv%info(polaramv)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(polaramv)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(polaramv)%levels(n)
                  if ( iv%polaramv(n)%u(k)%qc >= 0 .and. ob%polaramv(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%polaramv(n)%u(k), iv%polaramv(n)%u(k)%inv, iv%polaramv(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(polaramv)%lat(1,n), iv%info(polaramv)%lon(1,n), &
                     (iv%polaramv(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%polaramv(n)%v(k)%qc >= 0 .and. ob%polaramv(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%polaramv(n)%v(k), iv%polaramv(n)%v(k)%inv, iv%polaramv(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(polaramv)%lat(1,n), iv%info(polaramv)%lon(1,n), &
                     (iv%polaramv(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   ! [5] Transfer gpspw obs:

   if (iv%info(gpspw)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(gpspw)%nlocal
         if (iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(gpspw)%nlocal
            if (iv%info(gpspw)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(gpspw)%platform(n)(4:6),*) apm_int 
               if ( iv%gpspw(n)%tpw%qc >= 0 .and. ob%gpspw(n)%tpw /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%gpspw(n)%tpw, iv%gpspw(n)%tpw%inv, iv%gpspw(n)%tpw%error, &
                  'PWT', apm_plc, -888.88, iv%info(gpspw)%lat(1,n), iv%info(gpspw)%lon(1,n), &
                  -888.88, -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

   ! [6] Transfer sonde obs:

   if (iv%info(sound)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(sound)%nlocal
        if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(sound)%nlocal
            if (iv%info(sound)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = 120.
!               read(iv%info(sound)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(sound)%levels(n)
                  if ( iv%sound(n)%u(k)%qc >= 0 .and. ob%sound(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sound(n)%u(k), iv%sound(n)%u(k)%inv, iv%sound(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(sound)%lat(1,n), iv%info(sound)%lon(1,n), &
                     (iv%sound(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%sound(n)%v(k)%qc >= 0 .and. ob%sound(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sound(n)%v(k), iv%sound(n)%v(k)%inv, iv%sound(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(sound)%lat(1,n), iv%info(sound)%lon(1,n), &
                     (iv%sound(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%sound(n)%t(k)%qc >= 0 .and. ob%sound(n)%t(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sound(n)%t(k), iv%sound(n)%t(k)%inv, iv%sound(n)%t(k)%error, &
                     'TMP', apm_plc, -888.88, iv%info(sound)%lat(1,n), iv%info(sound)%lon(1,n), &
                     (iv%sound(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%sound(n)%q(k)%qc >= 0 .and. ob%sound(n)%q(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sound(n)%q(k), iv%sound(n)%q(k)%inv, iv%sound(n)%q(k)%error, &
                     'QVP', apm_plc, -888.88, iv%info(sound)%lat(1,n), iv%info(sound)%lon(1,n), &
                     (iv%sound(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   if (iv%info(sonde_sfc)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(sonde_sfc)%nlocal
        if (iv%info(sonde_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(sonde_sfc)%nlocal
            if (iv%info(sonde_sfc)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = 120.
!               read(iv%info(sonde_sfc)%platform(n)(4:6),*) apm_int 
                  if ( iv%sonde_sfc(n)%u%qc >= 0 .and. ob%sonde_sfc(n)%u /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sonde_sfc(n)%u, iv%sonde_sfc(n)%u%inv, iv%sonde_sfc(n)%u%error, &
                     'WNU', apm_plc, -888.88, iv%info(sonde_sfc)%lat(1,n), iv%info(sonde_sfc)%lon(1,n), &
                     (ob%sonde_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%sonde_sfc(n)%v%qc >= 0 .and. ob%sonde_sfc(n)%v /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sonde_sfc(n)%v, iv%sonde_sfc(n)%v%inv, iv%sonde_sfc(n)%v%error, &
                     'WNV', apm_plc, -888.88, iv%info(sonde_sfc)%lat(1,n), iv%info(sonde_sfc)%lon(1,n), &
                     (ob%sonde_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%sonde_sfc(n)%t%qc >= 0 .and. ob%sonde_sfc(n)%t /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sonde_sfc(n)%t, iv%sonde_sfc(n)%t%inv, iv%sonde_sfc(n)%t%error, &
                     'TMP', apm_plc, -888.88, iv%info(sonde_sfc)%lat(1,n), iv%info(sonde_sfc)%lon(1,n), &
                     (ob%sonde_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%sonde_sfc(n)%p%qc >= 0 .and. ob%sonde_sfc(n)%p /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sonde_sfc(n)%p, iv%sonde_sfc(n)%p%inv, iv%sonde_sfc(n)%p%error, &
                     'PRS', apm_plc, -888.88, iv%info(sonde_sfc)%lat(1,n), iv%info(sonde_sfc)%lon(1,n), &
                     (ob%sonde_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%sonde_sfc(n)%q%qc >= 0 .and. ob%sonde_sfc(n)%q /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%sonde_sfc(n)%q, iv%sonde_sfc(n)%q%inv, iv%sonde_sfc(n)%q%error, &
                     'QVP', apm_plc, -888.88, iv%info(sonde_sfc)%lat(1,n), iv%info(sonde_sfc)%lon(1,n), &
                     (ob%sonde_sfc(n)%p)/100., -888.88, apm_index 
                  end if
            end if
         end do
      end if
   end if

  ! [7] Transfer airep obs:

   if (iv%info(airep)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(airep)%nlocal
        if (iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(airep)%nlocal
            if (iv%info(airep)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(airep)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(airep)%levels(n)
                  if ( iv%airep(n)%u(k)%qc >= 0 .and. ob%airep(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%airep(n)%u(k), iv%airep(n)%u(k)%inv, iv%airep(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(airep)%lat(1,n), iv%info(airep)%lon(1,n), &
                     (iv%airep(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%airep(n)%v(k)%qc >= 0 .and. ob%airep(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%airep(n)%v(k), iv%airep(n)%v(k)%inv, iv%airep(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(airep)%lat(1,n), iv%info(airep)%lon(1,n), &
                     (iv%airep(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%airep(n)%t(k)%qc >= 0 .and. ob%airep(n)%t(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%airep(n)%t(k), iv%airep(n)%t(k)%inv, iv%airep(n)%t(k)%error, &
                     'TMP', apm_plc, -888.88, iv%info(airep)%lat(1,n), iv%info(airep)%lon(1,n), &
                     (iv%airep(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   ! [8] Transfer pilot obs:

   if (iv%info(pilot)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(pilot)%nlocal
        if (iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(pilot)%nlocal
            if (iv%info(pilot)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(pilot)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(pilot)%levels(n)
                  if ( iv%pilot(n)%u(k)%qc >= 0 .and. ob%pilot(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%pilot(n)%u(k), iv%pilot(n)%u(k)%inv, iv%pilot(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(pilot)%lat(1,n), iv%info(pilot)%lon(1,n), &
                     (iv%pilot(n)%p(k))/100, -888.88, apm_index 
                  end if
                  if ( iv%pilot(n)%v(k)%qc >= 0 .and. ob%pilot(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%pilot(n)%v(k), iv%pilot(n)%v(k)%inv, iv%pilot(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(pilot)%lat(1,n), iv%info(pilot)%lon(1,n), &
                     (iv%pilot(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   ! [9] Transfer SSM/I obs:SSMI:

   if (iv%info(ssmi_rv)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ssmi_rv)%nlocal
         if (iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(ssmi_rv)%nlocal
            if (iv%info(ssmi_rv)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(ssmi_rv)%platform(n)(4:6),*) apm_int 
               if ( iv%ssmi_rv(n)%speed%qc >= 0 .and. ob%ssmi_rv(n)%speed /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ssmi_rv(n)%speed, iv%ssmi_rv(n)%speed%inv, &
                                          iv%ssmi_rv(n)%speed%error, &
                  'SPD', apm_plc, -888.88, iv%info(ssmi_rv)%lat(1,n), iv%info(ssmi_rv)%lon(1,n), &
                  -888.88, -888.88, apm_index 
               end if
               if ( iv%ssmi_rv(n)%tpw%qc >= 0 .and. ob%ssmi_rv(n)%tpw /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%ssmi_rv(n)%tpw, iv%ssmi_rv(n)%tpw%inv, &
                                          iv%ssmi_rv(n)%tpw%error, &
                  'PWT', apm_plc, -888.88, iv%info(ssmi_rv)%lat(1,n), iv%info(ssmi_rv)%lon(1,n), &
                  -888.88, -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

! SSM/I TB not coded.

   ! [10] Transfer satem obs:

   if (iv%info(satem)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(satem)%nlocal
        if (iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(satem)%nlocal
            if (iv%info(satem)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(satem)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(satem)%levels(n)
                  if ( iv%satem(n)%thickness(k)%qc >= 0 .and. ob%satem(n)%thickness(k) /= missing_r ) then
                    apm_index = apm_index + 1 
                    write(ounit,1000) ob%satem(n)%thickness(k), iv%satem(n)%thickness(k)%inv, &
                                             iv%satem(n)%thickness(k)%error, &
                     'TCK', apm_plc, -888.88, iv%info(satem)%lat(1,n), iv%info(satem)%lon(1,n), &
                     (iv%satem(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

!  SSMT1 SSMT2 not coded.

  ! [11] Transfer scatterometer obs:

   if (iv%info(qscat)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(qscat)%nlocal
         if (iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(qscat)%nlocal
            if (iv%info(qscat)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(qscat)%platform(n)(4:6),*) apm_int 
               if ( iv%qscat(n)%u%qc >= 0 .and. ob%qscat(n)%u /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%qscat(n)%u, iv%qscat(n)%u%inv, iv%qscat(n)%u%error, &
                  'WNU', apm_plc, -888.88, iv%info(qscat)%lat(1,n), iv%info(qscat)%lon(1,n), &
                  -888.88, -888.88, apm_index 
               end if
               if ( iv%qscat(n)%v%qc >= 0 .and. ob%qscat(n)%v /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%qscat(n)%v, iv%qscat(n)%v%inv, iv%qscat(n)%v%error, &
                  'WNV', apm_plc, -888.88, iv%info(qscat)%lat(1,n), iv%info(qscat)%lon(1,n), &
                  -888.88, -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

  ! [12] Transfer profiler obs:

   if (iv%info(profiler)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(profiler)%nlocal
        if (iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(profiler)%nlocal
            if (iv%info(profiler)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(profiler)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(profiler)%levels(n)
                  if ( iv%profiler(n)%u(k)%qc >= 0 .and. ob%profiler(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%profiler(n)%u(k), iv%profiler(n)%u(k)%inv, iv%profiler(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(profiler)%lat(1,n), iv%info(profiler)%lon(1,n), &
                     (iv%profiler(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%profiler(n)%v(k)%qc >= 0 .and. ob%profiler(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%profiler(n)%v(k), iv%profiler(n)%v(k)%inv, iv%profiler(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(profiler)%lat(1,n), iv%info(profiler)%lon(1,n), &
                     (iv%profiler(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   ! Transfer Buoy obs:

   if (iv%info(buoy)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(buoy)%nlocal
         if (iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(buoy)%nlocal
            if (iv%info(buoy)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(buoy)%platform(n)(4:6),*) apm_int 
               if ( iv%buoy(n)%u%qc >= 0 .and. ob%buoy(n)%u /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%buoy(n)%u, iv%buoy(n)%u%inv, iv%buoy(n)%u%error, &
                  'WNU', apm_plc, -888.88, iv%info(buoy)%lat(1,n), iv%info(buoy)%lon(1,n), &
                  (ob%buoy(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%buoy(n)%v%qc >= 0 .and. ob%buoy(n)%v /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%buoy(n)%v, iv%buoy(n)%v%inv, iv%buoy(n)%v%error, &
                  'WNV', apm_plc, -888.88, iv%info(buoy)%lat(1,n), iv%info(buoy)%lon(1,n), &
                  (ob%buoy(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%buoy(n)%t%qc >= 0 .and. ob%buoy(n)%t /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%buoy(n)%t, iv%buoy(n)%t%inv, iv%buoy(n)%t%error, &
                  'TMP', apm_plc, -888.88, iv%info(buoy)%lat(1,n), iv%info(buoy)%lon(1,n), &
                  (ob%buoy(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%buoy(n)%p%qc >= 0 .and. ob%buoy(n)%p /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%buoy(n)%p, iv%buoy(n)%p%inv, iv%buoy(n)%p%error, &
                  'PRS', apm_plc, -888.88, iv%info(buoy)%lat(1,n), iv%info(buoy)%lon(1,n), &
                  (ob%buoy(n)%p)/100., -888.88, apm_index 
               end if
               if ( iv%buoy(n)%q%qc >= 0 .and. ob%buoy(n)%q /= missing_r ) then
                  apm_index = apm_index + 1
                  write(ounit,1000) ob%buoy(n)%q, iv%buoy(n)%q%inv, iv%buoy(n)%q%error, &
                  'QVP', apm_plc, -888.88, iv%info(buoy)%lat(1,n), iv%info(buoy)%lon(1,n), &
                  (ob%buoy(n)%p)/100., -888.88, apm_index 
               end if
            end if
         end do
      end if
   end if

   ! Transfer TC bogus obs:

   if (iv%info(bogus)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(bogus)%nlocal
        if (iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(bogus)%nlocal
            if (iv%info(bogus)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(bogus)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(bogus)%levels(n)
                  if ( iv%bogus(n)%u(k)%qc >= 0 .and. ob%bogus(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%bogus(n)%u(k), iv%bogus(n)%u(k)%inv, iv%bogus(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(bogus)%lat(1,n), iv%info(bogus)%lon(1,n), &
                     (iv%bogus(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%bogus(n)%v(k)%qc >= 0 .and. ob%bogus(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%bogus(n)%v(k), iv%bogus(n)%v(k)%inv, iv%bogus(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(bogus)%lat(1,n), iv%info(bogus)%lon(1,n), &
                     (iv%bogus(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%bogus(n)%t(k)%qc >= 0 .and.  ob%bogus(n)%t(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%bogus(n)%t(k), iv%bogus(n)%t(k)%inv, iv%bogus(n)%t(k)%error, &
                     'TMP', apm_plc, -888.88, iv%info(bogus)%lat(1,n), iv%info(bogus)%lon(1,n), &
                     (iv%bogus(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%bogus(n)%q(k)%qc >= 0 .and. ob%bogus(n)%q(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%bogus(n)%q(k), iv%bogus(n)%q(k)%inv, iv%bogus(n)%q(k)%error, &
                     'QVP', apm_plc, -888.88, iv%info(bogus)%lat(1,n), iv%info(bogus)%lon(1,n), &
                     (iv%bogus(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   ! Transfer AIRS retrievals:

   if (iv%info(airsr)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(airsr)%nlocal
        if (iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(airsr)%nlocal
            if (iv%info(airsr)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(airsr)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(airsr)%levels(n)
                  if ( iv%airsr(n)%t(k)%qc >= 0 .and. ob%airsr(n)%t(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%airsr(n)%t(k), iv%airsr(n)%t(k)%inv, iv%airsr(n)%t(k)%error, &
                     'TMP', apm_plc, -888.88, iv%info(airsr)%lat(1,n), iv%info(airsr)%lon(1,n), &
                     (iv%airsr(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%airsr(n)%q(k)%qc >= 0 .and. ob%airsr(n)%q(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%airsr(n)%q(k), iv%airsr(n)%q(k)%inv, iv%airsr(n)%q(k)%error, &
                     'QVP', apm_plc, -888.88, iv%info(airsr)%lat(1,n), iv%info(airsr)%lon(1,n), &
                     (iv%airsr(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   ! Transfer gpsref obs:
 
   if (iv%info(gpsref)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(gpsref)%nlocal
        if (iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(gpsref)%nlocal
            if (iv%info(gpsref)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = -888.88
!               read(iv%info(gpsref)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(gpsref)%levels(n)
                  if ( iv%gpsref(n)%ref(k)%qc >= 0 .and. ob%gpsref(n)%ref(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%gpsref(n)%ref(k), iv%gpsref(n)%ref(k)%inv, iv%gpsref(n)%ref(k)%error, &
                     'REF', apm_plc, -888.88, iv%info(gpsref)%lat(1,n), iv%info(gpsref)%lon(1,n), &
                     -888.88, -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if
  
    !  Transfer tamdar obs:

   if (iv%info(tamdar)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(tamdar)%nlocal
        if (iv%info(tamdar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(tamdar)%nlocal
            if (iv%info(tamdar)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = 220.
!               read(iv%info(tamdar)%platform(n)(4:6),*) apm_int 
               do k = 1, iv%info(tamdar)%levels(n)
                  if ( iv%tamdar(n)%u(k)%qc >= 0 .and. ob%tamdar(n)%u(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar(n)%u(k), iv%tamdar(n)%u(k)%inv, iv%tamdar(n)%u(k)%error, &
                     'WNU', apm_plc, -888.88, iv%info(tamdar)%lat(1,n), iv%info(tamdar)%lon(1,n), &
                     (iv%tamdar(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar(n)%v(k)%qc >= 0 .and. ob%tamdar(n)%v(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar(n)%v(k), iv%tamdar(n)%v(k)%inv, iv%tamdar(n)%v(k)%error, &
                     'WNV', apm_plc, -888.88, iv%info(tamdar)%lat(1,n), iv%info(tamdar)%lon(1,n), &
                     (iv%tamdar(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar(n)%t(k)%qc >= 0 .and. ob%tamdar(n)%t(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar(n)%t(k), iv%tamdar(n)%t(k)%inv, iv%tamdar(n)%t(k)%error, &
                     'TMP', apm_plc, -888.88, iv%info(tamdar)%lat(1,n), iv%info(tamdar)%lon(1,n), &
                     (iv%tamdar(n)%p(k))/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar(n)%q(k)%qc >= 0 .and. ob%tamdar(n)%q(k) /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar(n)%q(k), iv%tamdar(n)%q(k)%inv, iv%tamdar(n)%q(k)%error, &
                     'QVP', apm_plc, -888.88, iv%info(tamdar)%lat(1,n), iv%info(tamdar)%lon(1,n), &
                     (iv%tamdar(n)%p(k))/100., -888.88, apm_index 
                  end if
               end do
            end if
         end do
      end if
   end if

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(tamdar_sfc)%nlocal
        if (iv%info(tamdar_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         num_obs = 0
         do n = 1, iv%info(tamdar_sfc)%nlocal
            if (iv%info(tamdar_sfc)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               apm_plc = 220.
!               read(iv%info(tamdar_sfc)%platform(n)(4:6),*) apm_int 
                  if ( iv%tamdar_sfc(n)%u%qc >= 0 .and. ob%tamdar_sfc(n)%u /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar_sfc(n)%u, iv%tamdar_sfc(n)%u%inv, iv%tamdar_sfc(n)%u%error, &
                     'WNU', apm_plc, -888.88, iv%info(tamdar_sfc)%lat(1,n), iv%info(tamdar_sfc)%lon(1,n), &
                     (ob%tamdar_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar_sfc(n)%v%qc >= 0 .and. ob%tamdar_sfc(n)%v /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar_sfc(n)%v, iv%tamdar_sfc(n)%v%inv, iv%tamdar_sfc(n)%v%error, &
                     'WNV', apm_plc, -888.88, iv%info(tamdar_sfc)%lat(1,n), iv%info(tamdar_sfc)%lon(1,n), &
                     (ob%tamdar_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar_sfc(n)%t%qc >= 0 .and. ob%tamdar_sfc(n)%t /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar_sfc(n)%t, iv%tamdar_sfc(n)%t%inv, iv%tamdar_sfc(n)%t%error, &
                     'TMP', apm_plc, -888.88, iv%info(tamdar_sfc)%lat(1,n), iv%info(tamdar_sfc)%lon(1,n), &
                     (ob%tamdar_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar_sfc(n)%p%qc >= 0 .and. ob%tamdar_sfc(n)%p /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar_sfc(n)%p, iv%tamdar_sfc(n)%p%inv, iv%tamdar_sfc(n)%p%error, &
                     'PRS', apm_plc, -888.88, iv%info(tamdar_sfc)%lat(1,n), iv%info(tamdar_sfc)%lon(1,n), &
                     (ob%tamdar_sfc(n)%p)/100., -888.88, apm_index 
                  end if
                  if ( iv%tamdar_sfc(n)%q%qc >= 0 .and. ob%tamdar_sfc(n)%q /= missing_r ) then
                     apm_index = apm_index + 1
                     write(ounit,1000) ob%tamdar_sfc(n)%q, iv%tamdar_sfc(n)%q%inv, iv%tamdar_sfc(n)%q%error, &
                     'QVP', apm_plc, -888.88, iv%info(tamdar_sfc)%lat(1,n), iv%info(tamdar_sfc)%lon(1,n), &
                     (ob%tamdar_sfc(n)%p)/100., -888.88, apm_index 
                  end if
            end if
         end do
      end if
   end if

   close (ounit)
   call da_free_unit(ounit)

   if (trace_use) call da_trace_exit("da_write_obs_etkf")

end subroutine da_write_obs_etkf

subroutine da_write_filtered_obs(it, grid, ob, iv, &
   coarse_ix, coarse_jy, start_x, start_y)

   !------------------------------------------------------------------------
   ! Purpose: Writes filtered observations by WRFVAR
   !------------------------------------------------------------------------   
   !  Note:                                                           
   !   (a) Information just needed for WRFVAR  is only written  
   !   (b) Both U and V should be of good quality (qc > obs_qc_pointer) 
   !   (c) Pressure & Temperature quality should be > fails_error_max (= -3) 
   !   (d) Since there is no check for missing T,P and Q to get Qs in   
   !       "da_tp_to_qs", RH data is recovered from Q if 
   !       q_error = missing_r and both P and T is not missing
   !   (e) Since currently Sfc_Obs_Correction changes quality of
   !       observations to "surface_correction (= 2)" even if it is bad
   !       quality obs (qc < fails_error_max), such obs (element) is 
   !       made missing. Otherwise the bad obs (element) will also get 
   !       assmiilated with "max_iv_check = .false."
   !   (f) AMV's from Geostationary and Polar orbiting satellite are
   !       seperated & used as profile
   !   (g) Currently only MODIS winds from Polar orbiting satellite are used
   !  Modifications:
   !
   !..........................................................................
   ! The Quality Controlled observations (QCed obs) means
   !
   !   i) The check against the Background;
   !
   !  ii) The model bottom/top check, i.e. the observations below the
   !      lowest model level (note, not bottom) and higher than the highest
   !      model level (note, not model top) are filled with missing values
   !      except the Sound data;
   !
   ! iii) For Sound data, the above-top/below-bottom data are not written out.
   !
   ! This can be used for Var with check_max_iv = .false., or for
   ! verification with  analysis_type = 'VERIFY'.
   !---------------------------------------------------------------------------   

   implicit none

   integer, intent(in)            :: it
   type (domain),     intent(in)  :: grid
   type (y_type),     intent(in)  :: ob      ! Observation structure.
   type (iv_type),    intent(inout) :: iv      ! O-B structure.
   integer, intent(in)            :: coarse_ix, coarse_jy
   real, intent(in)               :: start_x, start_y

   integer                      :: j, k, iost, nlevels

   integer            :: n, zero
   real               :: es, qs, speed, dir, rh, rh_error, uu, vv
   real               :: height, height_error, pr_error, slp_error
   real               :: thick, thick_error,dir_error,speed_error
   real               :: td, td_error , elv_error, ships_pr_error

   real               :: slp_inv,speed_err,dir_err, p,p_err
   real               :: t, t_err,ref, ref_err, geom_h, geop_h
   integer            :: thick_qc, t_qc, td_qc, rh_qc, height_qc
   integer            :: slp_qc, speed_qc, dir_qc,p_qc, ref_qc

   integer            :: ounit, nz, mnz
   character(len=filename_len)       :: filename

   if (trace_use) call da_trace_entry("da_write_filtered_obs")

   ! Create noise and output:
   call da_get_unit(ounit)

   write(unit=filename, fmt='(a,i4.4)') 'filtered_obs.', myproc
   open(unit=ounit,file=trim(filename),form='formatted', & 
      status='replace', iostat=iost)
   if (iost /= 0) then
      call da_error("da_write_filtered_obs.inc",78, &
         (/"Cannot open filtered observation file"//filename/))
   end if

   zero = 0
   ! Note: Currently in WRFVAR there is no active role of  
   !       "height" and "TD"        

   height       = missing_r
   height_qc    = missing_data
   height_error = xmiss    
   td           = missing_r 
   td_qc        = missing_data
   td_error     = xmiss        
   ! Currently following errors are fixed 
   rh_error     = 10.0       
   pr_error = 100.
   ships_pr_error = 160.
   slp_error = 200.
   elv_error = 6.0
   dir_error = 5.0
   dir_err   = 5.0

   !--------------------------------------------------------------------------
   !  Write Synop:
   !--------------------------------------------------------------------------

   if (iv%info(synop)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout, fmt='(A,I6)') ' writing filtered obs for synop ', &
            iv%info(synop)%nlocal
      end if

      do n = 1, iv%info(synop)%nlocal
         if (.not.iv%info(synop)%proc_domain(1,n)) cycle
         ! Guo .. not write out the data below the lowest model level:
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(synop)%zk(1,n) == missing_r) cycle

         nlevels = iv%info(synop)%levels(n)
         if (iv%info(synop)%levels(n) > 1) then
            write(unit=stdout, fmt='(3A,I5,A)') &
               ' for SYNOP station ',iv%info(synop)%name(n),' got levels ',&
               iv%info(synop)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if
         write(ounit, fmt = fmt_info)                &
            iv%info(synop)%platform(n),    &
            iv%info(synop)%date_char(n),   &
            iv%info(synop)%name(n),        &
            nlevels,                      &
            iv%info(synop)%lat(1,n),         &
            iv%info(synop)%lon(1,n),         &
            iv%info(synop)%elv(n),         &
            iv%info(synop)%id(n)
         slp_inv= iv%info(synop)%slp(n)%inv
         slp_qc = iv%info(synop)%slp(n)%qc
         if (iv%synop(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(synop)%slp(n)%error, &
            iv%info(synop)%pw(n)%inv, iv%info(synop)%pw(n)%qc, iv%info(synop)%pw(n)%error

         ! get speed and direction
         if (iv%synop(n)%u%qc >= obs_qc_pointer .and. &
              iv%synop(n)%v%qc >= obs_qc_pointer) then
            uu = ob%synop(n)%u
            vv = ob%synop(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(synop)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%synop(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%synop(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
         end if

         ! get RH from Q & T    
         if (iv%synop(n)%q%qc >= obs_qc_pointer .and.  &
              abs(iv%synop(n)%q%error - missing_r) > 1. .and.  &
              abs(ob%synop(n)%t       - missing_r) > 1. .and.  &
              abs(ob%synop(n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%synop(n)%t, ob%synop(n)%p, es, qs)
            rh = (ob%synop(n)%q/qs) * 100.
            rh_qc = iv%synop(n)%q%qc
         else
            rh    = missing_r
            rh_qc = missing_data
         end if

         if (rh_qc < 0) rh_qc = missing_data

         p    = ob%synop(n)%p
         p_qc = iv%synop(n)%p%qc
         p_err = iv%synop(n)%p%error
         if (iv%synop(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error           
         end if

         t    =  ob%synop(n)%t
         t_qc =  iv%synop(n)%t%qc
         if (iv%synop(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
         end if

         write(ounit, fmt = trim (fmt_each))                   &
            p,  p_qc,  p_err,                                  &
            speed,speed_qc,        speed_err  ,                &
            dir  , dir_qc,         dir_err  ,                  &
            iv%info(synop)%elv(n), zero, elv_error,             &
            t, t_qc, iv%synop(n)%t%error,                      &
            td, td_qc, td_error,                               &
            rh, rh_qc,  rh_error                
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Metar:
   !--------------------------------------------------------------------------

   if (iv%info(metar)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for metar ',iv%info(metar)%nlocal
      end if

      do n = 1, iv%info(metar)%nlocal
         if (.not. iv%info(metar)%proc_domain(1,n)) cycle
         ! Do not write out the data below the lowest model level:
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(metar)%zk(1,n) == missing_r) cycle

         nlevels = iv%info(metar)%levels(n)
         if (iv%info(metar)%levels(n)> 1) then
            write(stdout,fmt='(3A,I5,A)') &
               ' for METAR station ',iv%info(metar)%name(n),' got levels ',&
               iv%info(metar)%levels(n),' but wrote only one level'
            nlevels = 1
         end if

         write(ounit, fmt = fmt_info)                &
            iv%info(metar)%platform(n),    &
            iv%info(metar)%date_char(n),   &
            iv%info(metar)%name(n),        &
            nlevels,                      &
            iv%info(metar)%lat(1,n),         &
            iv%info(metar)%lon(1,n),         &
            iv%info(metar)%elv(n),         &
            iv%info(metar)%id(n)
         slp_inv= iv%info(metar)%slp(n)%inv
         slp_qc = iv%info(metar)%slp(n)%qc
         if(iv%metar(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(metar)%slp(n)%error, &
            iv%info(metar)%pw(n)%inv, iv%info(metar)%pw(n)%qc, iv%info(metar)%pw(n)%error

         ! get speed and direction
         if (iv%metar(n)%u%qc >= obs_qc_pointer .and. &
            iv%metar(n)%v%qc >= obs_qc_pointer) then
            uu = ob%metar(n)%u
            vv = ob%metar(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(metar)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%metar(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%metar(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
         end if

         ! get RH from Q & T    
         if (iv%metar(n)%q%qc >= obs_qc_pointer .and.  &
             abs(iv%metar(n)%q%error - missing_r) > 1. .and.  &
             abs(ob%metar(n)%t       - missing_r) > 1. .and.  &
             abs(ob%metar(n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%metar(n)%t, ob%metar(n)%p, es, qs)
            rh = (ob%metar(n)%q/qs) * 100.
            rh_qc = iv%metar(n)%q%qc
         else
            rh    = missing_r
            rh_qc = missing_data
         end if
         if (rh_qc < 0) rh_qc = missing_data

         p    = ob%metar(n)%p
         p_qc = iv%metar(n)%p%qc
         p_err = iv%metar(n)%p%error
         if (iv%metar(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error           
         end if
         t    =  ob%metar(n)%t
         t_qc =  iv%metar(n)%t%qc
         if (iv%metar(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
         end if

         write(ounit, fmt = trim (fmt_each))                   &
            p,  p_qc,  p_err,                                  &
            speed,speed_qc,        speed_err  ,                &
            dir  , dir_qc,         dir_err  ,                  &
            iv%info(metar)%elv(n), zero, elv_error,             &
            t, t_qc, iv%metar(n)%t%error,                      &
            td, td_qc, td_error,                               &
            rh, rh_qc,  rh_error                
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Ships:
   !-------------------------------------------------------------------------- 

   if (iv%info(ships)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for ships ',iv%info(ships)%nlocal
      end if

      do n = 1, iv%info(ships)%nlocal
         if (.not. iv%info(ships)%proc_domain(1,n)) cycle
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(ships)%zk(1,n) == missing_r) cycle
               
         nlevels = iv%info(ships)%levels(n)
         if (iv%info(ships)%levels(n) > 1) then
            write(unit=stdout,fmt='(3A,I5,A)') &
               ' for SHIP station ',iv%info(ships)%name(n),' got levels ',&
               iv%info(ships)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if
         write(ounit, fmt = fmt_info)                &
            iv%info(ships)%platform(n),    &
            iv%info(ships)%date_char(n),   &
            iv%info(ships)%name(n),        &
            nlevels,               &
            iv%info(ships)%lat(1,n),       &
            iv%info(ships)%lon(1,n),       &
            iv%info(ships)%elv(n),         &
            iv%info(ships)%id(n)
            slp_inv= iv%info(ships)%slp(n)%inv
            slp_qc = iv%info(ships)%slp(n)%qc
         if (iv%ships(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(ships)%slp(n)%error, &
            iv%info(ships)%pw(n)%inv, iv%info(ships)%pw(n)%qc, iv%info(ships)%pw(n)%error

         ! get speed and direction
         if (iv%ships(n)%u%qc >= obs_qc_pointer .and. &
             iv%ships(n)%v%qc >= obs_qc_pointer) then
            uu = ob%ships(n)%u
            vv = ob%ships(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(ships)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%ships(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%ships(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
         end if

         ! get RH from Q & T    
         if (iv%ships(n)%q%qc >= obs_qc_pointer .and.  &
             abs(iv%ships(n)%q%error - missing_r) > 1. .and.  &
             abs(ob%ships(n)%t       - missing_r) > 1. .and.  &
             abs(ob%ships(n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%ships(n)%t, ob%ships(n)%p, es, qs)
            rh = (ob%ships(n)%q/qs) * 100.
            rh_qc = iv%ships(n)%q%qc
         else
            rh    = missing_r
            rh_qc = missing_data
         end if
         if (rh_qc < 0) rh_qc = missing_data

         p     = ob%ships(n)%p
         p_qc  = iv%ships(n)%p%qc
         p_err =  ships_pr_error
         if (iv%ships(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
         end if
         t    =  ob%ships(n)%t
         t_qc =  iv%ships(n)%t%qc
         if (iv%ships(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
         end if

         write(ounit, fmt = trim (fmt_each))                   &
            p,  p_qc,  p_err,                                  &
            speed,speed_qc,        speed_err  ,                &
            dir  , dir_qc,         dir_err  ,                  &
            iv%info(ships)%elv(n), zero, elv_error,             &
            t, t_qc, iv%ships(n)%t%error,                      &
            td, td_qc, td_error,                               &
            rh, rh_qc,  rh_error                
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Buoy :
   !--------------------------------------------------------------------------

   if (iv%info(buoy)%nlocal  > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for buoy  ',iv%info(buoy)%nlocal 
      end if

      do n = 1, iv%info(buoy)%nlocal
         if (.not. iv%info(buoy)%proc_domain(1,n)) cycle
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(buoy)%zk(1,n) == missing_r) cycle
                
         nlevels = iv%info(buoy)%levels(n)
         if (iv%info(buoy)%levels(n) > 1) then
            write(unit=stdout, fmt='(3A,I5,A)') &
               ' for BUOY  station ',iv%info(buoy)%name(n),' got levels ',&
               iv%info(buoy)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if

         write(ounit, fmt = fmt_info)                &
            iv%info(buoy)%platform(n),    &
            iv%info(buoy)%date_char(n),   &
            iv%info(buoy)%name(n),        &
            nlevels,         &
            iv%info(buoy)%lat(1,n),         &
            iv%info(buoy)%lon(1,n),         &
            iv%info(buoy)%elv(n),         &
            iv%info(buoy)%id(n)
         slp_inv= iv%info(buoy)%slp(n)%inv
         slp_qc = iv%info(buoy)%slp(n)%qc
         if (iv%buoy (n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(buoy)%slp(n)%error, &
            iv%info(buoy)%pw(n)%inv, iv%info(buoy)%pw(n)%qc, iv%info(buoy)%pw(n)%error

         ! get speed and direction
         if (iv%buoy (n)%u%qc >= obs_qc_pointer .and. &
              iv%buoy (n)%v%qc >= obs_qc_pointer) then
            uu = ob%buoy (n)%u
            vv = ob%buoy (n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(buoy)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%buoy (n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%buoy (n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
         end if

         ! get RH from Q & T    
         if (iv%buoy (n)%q%qc >= obs_qc_pointer .and.  &
             abs(iv%buoy (n)%q%error - missing_r) > 1. .and.  &
             abs(ob%buoy (n)%t       - missing_r) > 1. .and.  &
             abs(ob%buoy (n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%buoy (n)%t, ob%buoy (n)%p, es, qs)
            rh = (ob%buoy (n)%q/qs) * 100.
            rh_qc = iv%buoy (n)%q%qc
         else
            rh    = missing_r
           rh_qc = missing_data
         end if
         if (rh_qc < 0) rh_qc = missing_data

         p     = ob%buoy (n)%p
         p_qc  = iv%buoy (n)%p%qc
         p_err = iv%buoy (n)%p%error
         if (iv%buoy (n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error           
         end if
         t    =  ob%buoy (n)%t
         t_qc =  iv%buoy (n)%t%qc
         if (iv%buoy (n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
         end if

         write(ounit, fmt = trim (fmt_each))                   &
            p,  p_qc,   p_err,                                 &
            speed,speed_qc,        speed_err  ,                &
            dir  , dir_qc,         dir_err  ,                  &
            iv%info(buoy)%elv(n), zero, elv_error,             &
            t, t_qc, iv%buoy (n)%t%error,                      &
            td, td_qc, td_error,                               &
            rh, rh_qc,  rh_error                
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Geo. AMVs Obs:
   !--------------------------------------------------------------------------

   if (iv%info(geoamv)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for geoamv ',iv%info(geoamv)%nlocal
      end if

      do n = 1, iv%info(geoamv)%nlocal               
         if (.not. iv%info(geoamv)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(geoamv)%platform(n),    &
            iv%info(geoamv)%date_char(n),   &
            iv%info(geoamv)%name(n),        &
            iv%info(geoamv)%levels(n),      &
            iv%info(geoamv)%lat(1,n),         &
            iv%info(geoamv)%lon(1,n),         &
            iv%info(geoamv)%elv(n),         &
            iv%info(geoamv)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(geoamv)%slp(n)%inv, iv%info(geoamv)%slp(n)%qc, iv%info(geoamv)%slp(n)%error, &
            iv%info(geoamv)%pw(n)%inv, iv%info(geoamv)%pw(n)%qc, iv%info(geoamv)%pw(n)%error

         do k = 1, iv%info(geoamv)%levels(n)

            ! get speed and direction
            if (iv%geoamv(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%geoamv(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%geoamv(n)%u(k)
               vv = ob%geoamv(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(geoamv)%lon(1,n), convert_uv2fd)                        
               speed_qc  = iv%geoamv(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%geoamv(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc  = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

            write(ounit, fmt = trim (fmt_each))         &
               iv%geoamv(n)%p(k), speed_qc, slp_error,  &
               speed, speed_qc, speed_err,              &
               dir  , dir_qc,             dir_error,    &
               missing_r, standard_atmosphere, xmiss,   &
               missing_r, zero_t_td, xmiss,             & 
               missing_r, zero_t_td, xmiss,             & 
               missing_r, zero_t_td, xmiss   
         end do   
      end do   
   end if

   !--------------------------------------------------------------------------
   ! Write Polar AMVs Obs:
   !--------------------------------------------------------------------------

   if (iv%info(polaramv)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for polaramv ',iv%info(polaramv)%nlocal
      end if

      do n = 1, iv%info(polaramv)%nlocal               
         if (.not. iv%info(polaramv)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(polaramv)%platform(n),    &
            iv%info(polaramv)%date_char(n),   &
            iv%info(polaramv)%name(n),        &
            iv%info(polaramv)%levels(n),      &
            iv%info(polaramv)%lat(1,n),         &
            iv%info(polaramv)%lon(1,n),         &
            iv%info(polaramv)%elv(n),         &
            iv%info(polaramv)%id(n)

         write(ounit, fmt = fmt_srfc) &
            iv%info(polaramv)%slp(n)%inv, iv%info(polaramv)%slp(n)%qc,   &
            iv%info(polaramv)%slp(n)%error,                              &
            iv%info(polaramv)%pw(n)%inv, iv%info(polaramv)%pw(n)%qc,     &
            iv%info(polaramv)%pw(n)%error

         do k = 1, iv%info(polaramv)%levels(n)
            ! get speed and direction
            if (iv%polaramv(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%polaramv(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%polaramv(n)%u(k)
               vv = ob%polaramv(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(polaramv)%lon(k,n), convert_uv2fd)                        
               speed_qc  = iv%polaramv(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%polaramv(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc  = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

            write (ounit, fmt = trim (fmt_each))           &
               iv%polaramv(n)%p(k), speed_qc, slp_error,  &
               speed, speed_qc, speed_err,                &
               dir  , dir_qc,             dir_error,      &
               missing_r, standard_atmosphere, xmiss,     &
               missing_r, zero_t_td, xmiss,               & 
               missing_r, zero_t_td, xmiss,               & 
               missing_r, zero_t_td, xmiss   
         end do   
      end do   
   end if

   !--------------------------------------------------------------------------
   ! Write Sound 
   !--------------------------------------------------------------------------

   if (iv%info(sound)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(stdout,fmt='(A,I5)') &
            ' writing filtered obs for sound  ',iv%info(sound)%nlocal 
      end if
      mnz = 0

      do n = 1, iv%info(sound)%nlocal               
         if (.not. iv%info(sound)%proc_domain(1,n)) cycle
         if (iv%info(sonde_sfc)%platform(n) /= iv%info(sound)%platform(n) .or. &
             iv%info(sonde_sfc)%lat(1,n)    /= iv%info(sound)%lat(1,n)      .or. &
             iv%info(sonde_sfc)%lon(1,n)    /= iv%info(sound)%lon(1,n)      .or. &
             iv%info(sonde_sfc)%elv(n)      /= iv%info(sound)%elv(n)      .or. &
             iv%info(sonde_sfc)%id(n)       /= iv%info(sound)%id(n)      ) then
            write(unit=stderr,fmt='(A)')' Surface Sound Details:            '
            write(unit=stderr, fmt=fmt_info) &
               iv%info(sonde_sfc)%platform(n),    &
               iv%info(sonde_sfc)%date_char(n),   &
               iv%info(sonde_sfc)%name(n),        &
               iv%info(sonde_sfc)%levels(n),      &
               iv%info(sonde_sfc)%lat(1,n),         &
               iv%info(sonde_sfc)%lon(1,n),         &
               iv%info(sonde_sfc)%elv(n),         &
               iv%info(sonde_sfc)%id(n)

            write(unit=stderr,fmt='(A)') ' Upper level  Details: '
            write(unit=stderr, fmt=fmt_info) &
               iv%info(sound)%platform(n),    &
               iv%info(sound)%date_char(n),   &
               iv%info(sound)%name(n),        &
               iv%info(sound)%levels(n),      &
               iv%info(sound)%lat(1,n),         &
               iv%info(sound)%lon(1,n),         &
               iv%info(sound)%elv(n),         &
               iv%info(sound)%id(n)
            call da_error("da_write_filtered_obs.inc",659, &
               (/"Mismatch for Sound surface and upper air report"/))
         end if
         nz = 0
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(sonde_sfc)%zk(1,n) == missing_r) then
            nz = nz + 1
         end if

         do k = 1, iv%info(sound)%levels(n)
            if (iv%info(sound)%zk(k,n) == missing_r) then
               nz = nz + 1
            end if
         end do

         if (nz > 0) then
            mnz = mnz + 1
         end if

         nz = iv%info(sound)%levels(n) + 1 - nz
         if (nz < 1) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(sound)%platform(n),    &
            iv%info(sound)%date_char(n),   &
            iv%info(sound)%name(n),        &
            nz,                           &
            iv%info(sound)%lat(1,n),         &
            iv%info(sound)%lon(1,n),         &
            iv%info(sound)%elv(n),         &
            iv%info(sound)%id(n)
            ! iv%info(sound)%levels(n) + 1,  &
         slp_inv= iv%info(sound)%slp(n)%inv
         slp_qc = iv%info(sound)%slp(n)%qc
         if (iv%sonde_sfc(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(sound)%slp(n)%error,  &
            iv%info(sound)%pw(n)%inv, iv%info(sound)%pw(n)%qc, iv%info(sound)%pw(n)%error
 
         nz = 0
         j  = 0

         ! Not below the first model level
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(sonde_sfc)%zk(1,n) == missing_r) then
            nz = nz + 1
         else
            j = j + 1
            ! First write surface level information
 
            ! get speed and direction
            if (iv%sonde_sfc(n)%u%qc >= obs_qc_pointer .and. &
                iv%sonde_sfc(n)%v%qc >= obs_qc_pointer) then
               uu = ob%sonde_sfc(n)%u
               vv = ob%sonde_sfc(n)%v
               call da_ffdduv(speed, dir, uu, vv, iv%info(sound)%lon(1,n), &
                  convert_uv2fd)             
               speed_qc  = iv%sonde_sfc(n)%u%qc
               dir_qc    = speed_qc
               speed_err = iv%sonde_sfc(n)%u%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

            ! get RH from Q & T    
            if (iv%sonde_sfc(n)%q%qc >= obs_qc_pointer .and.  &
                abs(iv%sonde_sfc(n)%q%error - missing_r) > 1. .and.  &
                abs(ob%sonde_sfc(n)%t       - missing_r) > 1. .and.  &
                abs(ob%sonde_sfc(n)%p       - missing_r) > 1.   ) then 
               call da_tp_to_qs(ob%sonde_sfc(n)%t, ob%sonde_sfc(n)%p, es, qs)
               rh = (ob%sonde_sfc(n)%q/qs) * 100.
               rh_qc = iv%sonde_sfc(n)%q%qc
             else
                rh    = missing_r
                rh_qc = missing_data
             end if
            if (rh_qc < 0) rh_qc = missing_data

            p     = ob%sonde_sfc(n)%p
            p_qc  = iv%sonde_sfc(n)%p%qc
            p_err = iv%sonde_sfc(n)%p%error
            if (iv%sonde_sfc(n)%p%qc <= fails_error_max) then
               p_qc  = missing_data
               p     = missing_r
               p_err   = pr_error           
            end if

            t     =  ob%sonde_sfc(n)%t
            t_qc  =  iv%sonde_sfc(n)%t%qc
            t_err =  iv%sonde_sfc(n)%t%error
            if (iv%sonde_sfc(n)%t%qc <= fails_error_max) then
               t_qc  = missing_data
               t     = missing_r
               t_err = xmiss       
            end if

            write(ounit, fmt = trim (fmt_each))           &
               p,  p_qc,   p_err,                         &
               speed,speed_qc, speed_err ,                &
               dir  , dir_qc,  dir_err ,                  &
               iv%sonde_sfc(n)%h, zero, elv_error,        &
               t, t_qc, t_err,                            &
               td, td_qc, td_error,                       &
               rh, rh_qc,  rh_error
         end if
                
         do k = 1, iv%info(sound)%levels(n)

            if (iv%info(sound)%zk(k,n) == missing_r) then
                nz = nz + 1
            else
               j  = j + 1

               ! get speed and direction
               if (iv%sound(n)%u(k)%qc >= obs_qc_pointer .and. &
                   iv%sound(n)%v(k)%qc >= obs_qc_pointer) then  
                  uu = ob%sound(n)%u(k)
                  vv = ob%sound(n)%v(k)
                  call da_ffdduv(speed, dir, uu, vv, iv%info(sound)%lon(1,n), convert_uv2fd)                        
                  speed_qc  = iv%sound(n)%u(k)%qc
                  dir_qc    = speed_qc
                  speed_err = iv%sound(n)%u(k)%error
               else
                  speed     = missing_r
                  speed_qc  = missing_data
                  dir       = missing_r
                  dir_qc    = missing_data
                  speed_err = xmiss
               end if

               ! get RH from Q & T    
               if (iv%sound(n)%q(k)%qc >= obs_qc_pointer .and.  &
                   abs(iv%sound(n)%q(k)%error - missing_r) > 1. .and.  &
                   abs(ob%sound(n)%t(k)       - missing_r) > 1. .and.  &
                   abs(iv%sound(n)%p(k)       - missing_r) > 1.   ) then 
                  call da_tp_to_qs(ob%sound(n)%t(k), iv%sound(n)%p(k), es, qs)
                  rh = (ob%sound(n)%q(k)/qs) * 100.
                  rh_qc = iv%sound(n)%q(k)%qc
               else
                  rh    = missing_r
                  rh_qc = missing_data
               end if
               if (rh_qc < 0) rh_qc = missing_data

               t    =  ob%sound(n)%t(k)
               t_qc =  iv%sound(n)%t(k)%qc
               t_err= iv%sound(n)%t(k)%error
               if (iv%sound(n)%t(k)%qc <= fails_error_max) then
                  t_qc = missing_data
                  t   = missing_r
               end if

               write(ounit, fmt = trim (fmt_each))           &
                  iv%sound(n)%p(k), zero,    pr_error,       &
                  speed, speed_qc, speed_err ,               &
                  dir  , dir_qc,  dir_err ,                  &
                  iv%sound(n)%h(k), zero, elv_error,         &
                  t,  t_qc, t_err,                           &
                  td, td_qc, td_error,                       &
                  rh, rh_qc,  rh_error                
            end if
         end do
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Airep: 
   !--------------------------------------------------------------------------

   if (iv%info(airep)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for airep  ',iv%info(airep)%nlocal 
      end if

      do n = 1, iv%info(airep)%nlocal               
         if (.not. iv%info(airep)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(airep)%platform(n),    &
            iv%info(airep)%date_char(n),   &
            iv%info(airep)%name(n),        &
            iv%info(airep)%levels(n),      &
            iv%info(airep)%lat(1,n),         &
            iv%info(airep)%lon(1,n),         &
            iv%info(airep)%elv(n),         &
            iv%info(airep)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(airep)%slp(n)%inv, iv%info(airep)%slp(n)%qc, &
            iv%info(airep)%slp(n)%error,                         &
            iv%info(airep)%pw(n)%inv, iv%info(airep)%pw(n)%qc,   &
            iv%info(airep)%pw(n)%error

         do k = 1, iv%info(airep)%levels(n)

            ! get speed and direction
            if (iv%airep(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%airep(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%airep(n)%u(k)
               vv = ob%airep(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(airep)%lon(k,n), convert_uv2fd)             
               speed_qc  = iv%airep(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%airep(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

          ! get RH from Q & T
            if (iv%airep(n)%q(k)%qc >= obs_qc_pointer .and.  &
                abs(iv%airep(n)%q(k)%error - missing_r) > 1. .and.  &
                abs(ob%airep(n)%t(k)       - missing_r) > 1. .and.  &
                abs(iv%airep(n)%p(k)       - missing_r) > 1.   ) then
                call da_tp_to_qs(ob%airep(n)%t(k), iv%airep(n)%p(k), es, qs)
                rh = (ob%airep(n)%q(k)/qs) * 100.
                rh_qc = iv%airep(n)%q(k)%qc
            else
                rh    = missing_r
                rh_qc = missing_data
            end if
            if (rh_qc < 0) rh_qc = missing_data

            t    =  ob%airep(n)%t(k)
            t_qc =  iv%airep(n)%t(k)%qc
            t_err= iv%airep(n)%t(k)%error
            if (iv%airep(n)%t(k)%qc <= fails_error_max) then
               t_qc = missing_data
               t   = missing_r
            end if

            write(ounit, fmt = trim (fmt_each))       &
               iv%airep(n)%p(k),zero , pr_error,      &
               speed, speed_qc, speed_err,            &
               dir  , dir_qc, dir_err,                &
               iv%airep(n)%h(k), zero, elv_error,     &
               t,  t_qc,   t_err,                     &
               td, td_qc, td_error,                   &
               rh, rh_qc,  rh_error
         end do 
      end do     
   end if

   !--------------------------------------------------------------------------
   ! Write Pilot:
   !--------------------------------------------------------------------------

   if (iv%info(pilot)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
           ' writing filtered obs for pilot  ',iv%info(pilot)%nlocal 
      end if

      do n = 1, iv%info(pilot)%nlocal               
         if (.not. iv%info(pilot)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(pilot)%platform(n),    &
            iv%info(pilot)%date_char(n),   &
            iv%info(pilot)%name(n),        &
            iv%info(pilot)%levels(n),      &
            iv%info(pilot)%lat(1,n),         &
            iv%info(pilot)%lon(1,n),         &
            iv%info(pilot)%elv(n),         &
            iv%info(pilot)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(pilot)%slp(n)%inv, iv%info(pilot)%slp(n)%qc, &
            iv%info(pilot)%slp(n)%error,                         &
            iv%info(pilot)%pw(n)%inv, iv%info(pilot)%pw(n)%qc,   &
            iv%info(pilot)%pw(n)%error

         do k = 1, iv%info(pilot)%levels(n)

            ! get speed and direction
            if (iv%pilot(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%pilot(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%pilot(n)%u(k)
               vv = ob%pilot(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(pilot)%lon(k,n), &
                  convert_uv2fd)                        
               speed_qc  = iv%pilot(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%pilot(n)%u(k)%error
            else
               speed     = missing_r
            speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

            write(ounit, fmt = trim (fmt_each))       &
               iv%pilot(n)%p(k),zero, pr_error,       &
               speed, speed_qc, speed_err,            &
               dir  , dir_qc, dir_error,              &
               missing_r, standard_atmosphere, xmiss, &
               missing_r, missing_data, xmiss,        &
               missing_r, missing_data, xmiss,        &
               missing_r, missing_data, xmiss
         end do 
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write satem:          
   ! Note: The satem ref_p is put in the SLP position in OBSPROC data stream.
   !--------------------------------------------------------------------------

   if (iv%info(satem)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for satem  ',iv%info(satem)%nlocal 
      end if

      do n = 1, iv%info(satem)%nlocal           
         if (.not. iv%info(satem)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(satem)%platform(n),    &
            iv%info(satem)%date_char(n),   &
            iv%info(satem)%name(n),        &
            iv%info(satem)%levels(n),      &
            iv%info(satem)%lat(1,n),         &
            iv%info(satem)%lon(1,n),         &
            iv%info(satem)%elv(n),         &
            iv%info(satem)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(satem)%slp(n)%inv, iv%info(satem)%slp(n)%qc, &
            iv%info(satem)%slp(n)%error,                         &
            iv%info(satem)%pw(n)%inv, iv%info(satem)%pw(n)%qc,   &
            iv%info(satem)%pw(n)%error

         do k = 1, iv%info(satem)%levels(n) 
            thick_qc = iv%satem(n)%thickness(k)%qc
            thick = iv%satem(n)%org_thickness(k)%inv   
            thick_error = iv%satem(n)%org_thickness(k)%error
            if (iv%satem(n)%thickness(k)%qc < obs_qc_pointer) &
               thick_qc = missing_data

            write(ounit, fmt = trim (fmt_each))               &
               iv%satem(n)%p(k),zero, slp_error,                 &
               missing_r, missing_data, xmiss,                   &
               missing_r, missing_data, xmiss,                   &
               missing_r, standard_atmosphere, xmiss,            &
               thick, thick_qc, thick_error,                     &
               missing_r, missing_data, xmiss,                   &
               missing_r, missing_data, xmiss
         end do
      end do
   end if

   !--------------------------------------------------------------------------
   ! Write Qscat: 
   !--------------------------------------------------------------------------

   if (iv%info(qscat)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for qscat  ',iv%info(qscat)%nlocal 
      end if
 
      do n = 1, iv%info(qscat)%nlocal               
         if (.not. iv%info(qscat)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(qscat)%platform(n),    &
            iv%info(qscat)%date_char(n),   &
            iv%info(qscat)%name(n),        &
            iv%info(qscat)%levels(n),      &
            iv%info(qscat)%lat(1,n),         &
            iv%info(qscat)%lon(1,n),         &
            iv%info(qscat)%elv(n),         &
            iv%info(qscat)%id(n)

         write(ounit, fmt = fmt_srfc)  &
            iv%info(qscat)%slp(n)%inv, iv%info(qscat)%slp(n)%qc, &
            iv%info(qscat)%slp(n)%error,                         &
            iv%info(qscat)%pw(n)%inv, iv%info(qscat)%pw(n)%qc,   &
            iv%info(qscat)%pw(n)%error

         ! get speed and direction
         if (iv%qscat(n)%u%qc >= obs_qc_pointer .and. &
             iv%qscat(n)%v%qc >= obs_qc_pointer) then  
            uu = ob%qscat(n)%u
            vv = ob%qscat(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(qscat)%lon(1,n), &
               convert_uv2fd)                        
            speed_error = iv%qscat(n)%u%error
            speed_qc  = iv%qscat(n)%u%qc
            dir_qc    = speed_qc
         else
            speed = missing_r
            dir   = missing_r
            speed_error = xmiss
            speed_qc = missing_data
            dir_qc   = missing_data
         end if

         write(ounit, fmt = trim (fmt_each))          &
            iv%qscat(n)%h,iv%qscat(n)%u%qc,slp_error, &  
            speed,speed_qc,        speed_error,       &
            dir  ,dir_qc          ,dir_error,         &
            iv%info(qscat)%elv(n), zero, elv_error,    &
            missing_r, missing_data, xmiss,           & 
            missing_r, missing_data, xmiss,           & 
            missing_r, missing_data, xmiss

      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Profiler: 
   !--------------------------------------------------------------------------
   
   if (iv%info(profiler)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for profiler  ',iv%info(profiler)%nlocal 
      end if

      do n = 1, iv%info(profiler)%nlocal            
         if (.not. iv%info(profiler)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost)    &
            iv%info(profiler)%platform(n),    &
            iv%info(profiler)%date_char(n),   &
            iv%info(profiler)%name(n),        &
            iv%info(profiler)%levels(n),      &
            iv%info(profiler)%lat(1,n),         &
            iv%info(profiler)%lon(1,n),         &
            iv%info(profiler)%elv(n),         &
            iv%info(profiler)%id(n)
         write(ounit, fmt = fmt_srfc) &
            iv%info(profiler)%slp(n)%inv, iv%info(profiler)%slp(n)%qc, &
            iv%info(profiler)%slp(n)%error,                            &
            iv%info(profiler)%pw(n)%inv, iv%info(profiler)%pw(n)%qc,   &
            iv%info(profiler)%pw(n)%error

         do k = 1, iv%info(profiler)%levels(n)

            ! get speed and direction
            if (iv%profiler(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%profiler(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%profiler(n)%u(k)
               vv = ob%profiler(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(profiler)%lon(k,n), &
                  convert_uv2fd)
               speed_error = iv%profiler(n)%u(k)%error
               speed_qc    = iv%profiler(n)%u(k)%qc
               dir_qc      = speed_qc
            else
               speed = missing_r
               dir   = missing_r
               speed_error = 1.1
               speed_qc = missing_data
               dir_qc   = missing_data
           end if

           write(ounit, fmt = trim (fmt_each))          &
              iv%profiler(n)%p(k),zero,pr_error,        &
              speed,speed_qc              ,speed_error, &
              dir  ,dir_qc                ,dir_error,   &
              iv%profiler(n)%h(k), zero, elv_error, &
              missing_r, missing_data, xmiss,           & 
              missing_r, missing_data, xmiss,           & 
              missing_r, missing_data, xmiss

         end do 
      end do 
   end if


   ! [10] Write SSMT1:          
   if (iv%info(ssmt1)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for ssmt1 ', iv%info(ssmt1)%nlocal 
      end if

      do n = 1, iv%info(ssmt1)%nlocal             
         if (.not. iv%info(ssmt1)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(ssmt1)%platform(n),    &
            iv%info(ssmt1)%date_char(n),   &
            iv%info(ssmt1)%name(n),        &
            iv%info(ssmt1)%levels(n),      &
            iv%info(ssmt1)%lat(1,n),         &
            iv%info(ssmt1)%lon(1,n),         &
            iv%info(ssmt1)%elv(n),         &
            iv%info(ssmt1)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(ssmt1)%slp(n)%inv, iv%info(ssmt1)%slp(n)%qc, &
            iv%info(ssmt1)%slp(n)%error,                         &
            iv%info(ssmt1)%pw(n)%inv, iv%info(ssmt1)%pw(n)%qc,   &
            iv%info(ssmt1)%pw(n)%error

         do k = 1, iv%info(ssmt1)%levels(n)
            write(ounit, fmt = trim (fmt_each))                            &
               iv%ssmt1(n)%p(k),zero,  slp_error,                          &
               missing_r, missing, missing_r,                              &
               missing_r, missing, missing_r,                              &
               height,          height_qc,height_error,                    &
               ob%ssmt1(n)%t(k),iv%ssmt1(n)%t(k)%qc,iv%ssmt1(n)%t(k)%error,&
               td, td_qc, td_error,                                        &   
               missing_r, missing, missing_r                                 
         end do 
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write SSMT2:    
   !--------------------------------------------------------------------------      

   if (iv%info(ssmt2)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for ssmt2 ', iv%info(ssmt2)%nlocal 
      end if

      do n = 1, iv%info(ssmt2)%nlocal              
         if (.not. iv%info(ssmt2)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(ssmt2)%platform(n),    &
            iv%info(ssmt2)%date_char(n),   &
            iv%info(ssmt2)%name(n),        &
            iv%info(ssmt2)%levels(n),      &
            iv%info(ssmt2)%lat(1,n),         &
            iv%info(ssmt2)%lon(1,n),         &
            iv%info(ssmt2)%elv(n),         &
            iv%info(ssmt2)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(ssmt2)%slp(n)%inv, iv%info(ssmt2)%slp(n)%qc, &
            iv%info(ssmt2)%slp(n)%error,                         &
            iv%info(ssmt2)%pw(n)%inv, iv%info(ssmt2)%pw(n)%qc,   &
            iv%info(ssmt2)%pw(n)%error

         do k = 1, iv%info(ssmt2)%levels(n)
            write(ounit, fmt = trim (fmt_each))                  &
            iv%ssmt2(n)%p(k),zero,  slp_error,                &
            missing_r, missing, missing_r,                    &
            missing_r, missing, missing_r,                    &
            height,          height_qc,height_error,          &
            missing_r, missing, missing_r,                    &
            td, td_qc, td_error,                              &   
            ob%ssmt2(n)%rh(k),iv%ssmt2(n)%rh(k)%qc,rh_error                 
         end do  
      end do 
   end if

   !--------------------------------------------------------------------------
   !  Write Gpspw:          
   !--------------------------------------------------------------------------

   if (iv%info(gpspw)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for gpspw ', iv%info(gpspw)%nlocal 
      end if

      do n = 1, iv%info(gpspw)%nlocal
         if (.not. iv%info(gpspw)%proc_domain(1,n)) cycle
         ! hardwired the # of levels = 0 for GPSPW (YRG 03/02/2006)
         iv%info(gpspw)%levels(n) = 0
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(gpspw)%platform(n),    &
            iv%info(gpspw)%date_char(n),   &
            iv%info(gpspw)%name(n),        &
            iv%info(gpspw)%levels(n),      &
            iv%info(gpspw)%lat(1,n),         &
            iv%info(gpspw)%lon(1,n),         &
            iv%info(gpspw)%elv(n),         &
            iv%info(gpspw)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(gpspw)%slp(n)%inv, iv%info(gpspw)%slp(n)%qc, &
            iv%info(gpspw)%slp(n)%error,                         &
            iv%info(gpspw)%pw(n)%inv, iv%info(gpspw)%pw(n)%qc,   &
            iv%info(gpspw)%pw(n)%error
      end do 
   end if

   ! Write Gps Refractivity:
   if (iv%info(gpsref)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for gpsref ', iv%info(gpsref)%nlocal
      end if

      do n = 1, iv%info(gpsref)%nlocal
         if (.not. iv%info(gpsref)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(gpsref) %platform(n),    &
            iv%info(gpsref) %date_char(n),   &
            iv%info(gpsref) %name(n),        &
            iv%info(gpsref) %levels(n),      &
            iv%info(gpsref) %lat(1,n),         &
            iv%info(gpsref) %lon(1,n),         &
            iv%info(gpsref) %elv(n),         &
            iv%info(gpsref) %id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(gpsref)%slp(n)%inv, iv%info(gpsref)%slp(n)%qc, &
            iv%info(gpsref)%slp(n)%error,                          &
            iv%info(gpsref)%pw(n)%inv, iv%info(gpsref)%pw(n)%qc,   &
            iv%info(gpsref)%pw(n)%error

         do k = 1, iv%info(gpsref)%levels(n)

            ! get RH from Q & T    
            if (iv%gpsref(n)%q(k)%qc >= obs_qc_pointer .and.  &
                abs(iv%gpsref(n)%q(k)%error - missing_r) > 1. .and.  &
                abs(ob%gpsref(n)%t(k)       - missing_r) > 1. .and.  &
                abs(ob%gpsref(n)%p(k)       - missing_r) > 1.   ) then 
               call da_tp_to_qs(ob%gpsref(n)%t(k), ob%gpsref(n)%p(k), es, qs)
               rh = (ob%gpsref(n)%q(k)/qs) * 100.
               rh_qc = iv%gpsref(n)%q(k)%qc
            else
               rh    = missing_r
               rh_qc = missing_data
            end if

            p    = ob%gpsref(n)%p(k)
            p_qc = iv%gpsref(n)%p(k)%qc 
            if (iv%gpsref(n)%p(k)%qc <= fails_error_max .and. &
                .not.(iv%gpsref(n)%p(k)%qc == -5)) then
               p    = missing_r     
               p_qc = missing_data  
            end if

            t    = ob%gpsref(n)%t(k)
            t_qc  = iv%gpsref(n)%t(k)%qc
            t_err = iv%gpsref(n)%t(k)%error   
            if (iv%gpsref(n)%t(k)%qc <= fails_error_max) then
               t    = missing_r     
               t_qc = missing_data  
            end if

            ref = ob%gpsref(n)%ref(k)
            ref_qc = iv%gpsref(n)%ref(k)%qc
            ref_err = iv%gpsref(n)%ref(k)%error   
            if (iv%gpsref(n)%ref(k)%qc <= fails_error_max) then
               ref    = missing_r     
               ref_qc = missing_data  
            end if
!
!  For GPSREF, convert the geopotential height to geometric height, then write
!  out to filtered_obs (YRG, 06/15/2011):
            geop_h = iv%gpsref(n)%h(k) / 1000.0
            call da_geo2msl1 (geop_h, iv%info(gpsref)%lat(1,n), &
                              iv%info(gpsref)%lon(1,n), geom_h)
! ............................................................................
            write(ounit, fmt = trim (fmt_each))             &
               p,  p_qc, slp_error,                         &
               missing_r, missing_data, xmiss,              &
               missing_r, missing_data, xmiss,              &
               geom_h*1000.0,     zero, height_error,       &
               t, t_qc,  t_err,                             &   
               ref, ref_qc,  ref_err,                       &   
               rh, rh_qc, iv%gpsref(n)%q(k)%error
         end do 
      end do
   end if

   !--------------------------------------------------------------------------
   ! Write AIRS  retrievals:
   !--------------------------------------------------------------------------

   if (iv%info(airsr)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for AIRS retrievals ',iv%info(airsr)%nlocal 
      end if

      do n = 1, iv%info(airsr)%nlocal               
         if (.not. iv%info(airsr)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(airsr)%platform(n),    &
            iv%info(airsr)%date_char(n),   &
            iv%info(airsr)%name(n),        &
            iv%info(airsr)%levels(n),      &
            iv%info(airsr)%lat(1,n),         &
            iv%info(airsr)%lon(1,n),         &
            iv%info(airsr)%elv(n),         &
            iv%info(airsr)%id(n)
         write(ounit, fmt = fmt_srfc) &
            iv%info(airsr)%slp(n)%inv, iv%info(airsr)%slp(n)%qc, &
            iv%info(airsr)%slp(n)%error,                         &
            iv%info(airsr)%pw(n)%inv, iv%info(airsr)%pw(n)%qc,   &
            iv%info(airsr)%pw(n)%error

         do k = 1, iv%info(airsr)%levels(n)
            ! get RH from Q & T    
            if (iv%airsr(n)%q(k)%qc >= obs_qc_pointer .and.  &
                abs(iv%airsr(n)%q(k)%error - missing_r) > 1. .and.  &
                abs(ob%airsr(n)%t(k)       - missing_r) > 1. .and.  &
                abs(iv%airsr(n)%p(k)       - missing_r) > 1.   ) then 
               call da_tp_to_qs(ob%airsr(n)%t(k), iv%airsr(n)%p(k), es, qs)
               rh = (ob%airsr(n)%q(k)/qs) * 100.
               rh_qc = iv%airsr(n)%q(k)%qc
            else
               rh    = missing_r
               rh_qc = missing_data
            end if
            if (rh_qc < 0) rh_qc = missing_data

            t     = ob%airsr(n)%t(k)
            t_qc  = iv%airsr(n)%t(k)%qc
            t_err = iv%airsr(n)%t(k)%error
            if (iv%airsr(n)%t(k)%qc <= fails_error_max) then
               t_qc = missing_data
               t    = missing_r
            end if

            write(ounit, fmt = trim (fmt_each))      &
               iv%airsr(n)%p(k), zero,    pr_error,     &
               missing_r, missing_data, xmiss,          &
               missing_r, missing_data, xmiss,          &
               iv%airsr(n)%h(k), zero, elv_error,       &
               t,  t_qc, t_err,                         &
               td, td_qc, td_error,                     &
               rh, rh_qc,  rh_error                
         end do 
      end do     
   end if
   !--------------------------------------------------------------------------
   ! Write mtgirs
   !--------------------------------------------------------------------------

   if (iv%info(mtgirs)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for mtgirs  ',iv%info(mtgirs)%nlocal
      end if

      do n = 1, iv%info(mtgirs)%nlocal
         if (.not. iv%info(mtgirs)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(mtgirs)%platform(n),    &
            iv%info(mtgirs)%date_char(n),   &
            iv%info(mtgirs)%name(n),        &
            iv%info(mtgirs)%levels(n),      &
            iv%info(mtgirs)%lat(1,n),         &
            iv%info(mtgirs)%lon(1,n),         &
            iv%info(mtgirs)%elv(n),         &
            iv%info(mtgirs)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(mtgirs)%slp(n)%inv, iv%info(mtgirs)%slp(n)%qc, &
            iv%info(mtgirs)%slp(n)%error,                          &
            iv%info(mtgirs)%pw(n)%inv, iv%info(mtgirs)%pw(n)%qc,   &
            iv%info(mtgirs)%pw(n)%error

         do k = 1, iv%info(mtgirs)%levels(n)

            ! get speed and direction
            if (iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer) then
               uu = ob%mtgirs(n)%u(k)
               vv = ob%mtgirs(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(mtgirs)%lon(k,n), convert_uv2fd)
               speed_qc  = iv%mtgirs(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%mtgirs(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

!q from kg/kg to g/kg, for error is to 10e-2g/kg
               if (iv%mtgirs(n)%q(k)%qc >= obs_qc_pointer) then
                  rh = ob%mtgirs(n)%q(k) * 1000.
                  rh_qc = iv%mtgirs(n)%q(k)%qc
                  rh_error=iv%mtgirs(n)%q(k)%error * 1000.*100.
               else
                  rh    = missing_r
                  rh_qc = missing_data
                  rh_error=xmiss
               end if

            write(ounit, fmt = trim (fmt_each))       &
               iv%mtgirs(n)%p(k),zero , pr_error,      &
               speed, speed_qc, speed_err,            &
               dir  , dir_qc, dir_error,              &
               iv%mtgirs(n)%h(k), zero, elv_error,         &
                   t,  t_qc, t_err,                           &
               missing_r, missing_data, xmiss,        &
                  rh, rh_qc,  rh_error
         end do
      end do    
   end if

   !--------------------------------------------------------------------------
   ! Write tamdar
   !--------------------------------------------------------------------------

   if (iv%info(tamdar)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for tamdar  ',iv%info(tamdar)%nlocal
      end if

      do n = 1, iv%info(tamdar)%nlocal
         if (.not. iv%info(tamdar)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(tamdar)%platform(n),    &
            iv%info(tamdar)%date_char(n),   &
            iv%info(tamdar)%name(n),        &
            iv%info(tamdar)%levels(n),      &
            iv%info(tamdar)%lat(1,n),         &
            iv%info(tamdar)%lon(1,n),         &
            iv%info(tamdar)%elv(n),         &
            iv%info(tamdar)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(tamdar)%slp(n)%inv, iv%info(tamdar)%slp(n)%qc, &
            iv%info(tamdar)%slp(n)%error,                         &
            iv%info(tamdar)%pw(n)%inv, iv%info(tamdar)%pw(n)%qc,   &
            iv%info(tamdar)%pw(n)%error

         do k = 1, iv%info(tamdar)%levels(n)

            ! get speed and direction
            if (iv%tamdar(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%tamdar(n)%v(k)%qc >= obs_qc_pointer) then
               uu = ob%tamdar(n)%u(k)
               vv = ob%tamdar(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(tamdar)%lon(k,n), convert_uv2fd)
               speed_qc  = iv%tamdar(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%tamdar(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
            end if

            t    =  ob%tamdar(n)%t(k)
            t_qc =  iv%tamdar(n)%t(k)%qc
            t_err= iv%tamdar(n)%t(k)%error
            if (iv%tamdar(n)%t(k)%qc <= fails_error_max) then
               t_qc = missing_data
               t   = missing_r
            end if
               ! get RH from Q & T
               if (iv%tamdar(n)%q(k)%qc >= obs_qc_pointer .and.  &
                   abs(iv%tamdar(n)%q(k)%error - missing_r) > 1. .and.  &
                   abs(ob%tamdar(n)%t(k)       - missing_r) > 1. .and.  &
                   abs(iv%tamdar(n)%p(k)       - missing_r) > 1.   ) then
                  call da_tp_to_qs(ob%tamdar(n)%t(k), iv%tamdar(n)%p(k), es, qs)
                  rh = (ob%tamdar(n)%q(k)/qs) * 100.
                  rh_qc = iv%tamdar(n)%q(k)%qc
               else
                  rh    = missing_r
                  rh_qc = missing_data
               end if


            write(ounit, fmt = trim (fmt_each))       &
               iv%tamdar(n)%p(k),zero , pr_error,      &
               speed, speed_qc, speed_err,            &
               dir  , dir_qc, dir_error,              &
               missing_r, standard_atmosphere, xmiss, &
               t,  t_qc,   t_err,                     &
               missing_r, missing_data, xmiss,        &
               rh, rh_qc, rh_error               
         end do
      end do
   end if

   !--------------------------------------------------------------------------
   ! Write tamdar_sfc
   !--------------------------------------------------------------------------

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout, fmt='(A,I6)') ' writing filtered obs for tamdar_sfc ', &
            iv%info(tamdar_sfc)%nlocal
      end if

      do n = 1, iv%info(tamdar_sfc)%nlocal
         if (.not.iv%info(tamdar_sfc)%proc_domain(1,n)) cycle
         ! Guo .. not write out the data below the lowest model level:
         if (sfc_assi_options == sfc_assi_options_1 .and.  &
             iv%info(tamdar_sfc)%zk(1,n) == missing_r) cycle

         nlevels = iv%info(tamdar_sfc)%levels(n)
         if (iv%info(tamdar_sfc)%levels(n) > 1) then
            write(unit=stdout, fmt='(3A,I5,A)') &
               ' for SYNOP station ',iv%info(tamdar_sfc)%name(n),' got levels ',&
               iv%info(tamdar_sfc)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if
         write(ounit, fmt = fmt_info)                &
            iv%info(tamdar_sfc)%platform(n),    &
            iv%info(tamdar_sfc)%date_char(n),   &
            iv%info(tamdar_sfc)%name(n),        &
            nlevels,                      &
            iv%info(tamdar_sfc)%lat(1,n),         &
            iv%info(tamdar_sfc)%lon(1,n),         &
            iv%info(tamdar_sfc)%elv(n),         &
            iv%info(tamdar_sfc)%id(n)
         slp_inv= iv%info(tamdar_sfc)%slp(n)%inv
         slp_qc = iv%info(tamdar_sfc)%slp(n)%qc
         if (iv%tamdar_sfc(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(tamdar_sfc)%slp(n)%error, &
            iv%info(tamdar_sfc)%pw(n)%inv, iv%info(tamdar_sfc)%pw(n)%qc, iv%info(tamdar_sfc)%pw(n)%error

         ! get speed and direction
         if (iv%tamdar_sfc(n)%u%qc >= obs_qc_pointer .and. &
              iv%tamdar_sfc(n)%v%qc >= obs_qc_pointer) then
            uu = ob%tamdar_sfc(n)%u
            vv = ob%tamdar_sfc(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(tamdar_sfc)%lon(1,n), &
               convert_uv2fd)
            speed_qc  = iv%tamdar_sfc(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%tamdar_sfc(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
         end if

         ! get RH from Q & T
         if (iv%tamdar_sfc(n)%q%qc >= obs_qc_pointer .and.  &
              abs(iv%tamdar_sfc(n)%q%error - missing_r) > 1. .and.  &
              abs(ob%tamdar_sfc(n)%t       - missing_r) > 1. .and.  &
              abs(ob%tamdar_sfc(n)%p       - missing_r) > 1.   ) then
            call da_tp_to_qs(ob%tamdar_sfc(n)%t, ob%tamdar_sfc(n)%p, es, qs)
            rh = (ob%tamdar_sfc(n)%q/qs) * 100.
            rh_qc = iv%tamdar_sfc(n)%q%qc
         else
            rh    = missing_r
            rh_qc = missing_data
         end if

         if (rh_qc < 0) rh_qc = missing_data

         p    = ob%tamdar_sfc(n)%p
         p_qc = iv%tamdar_sfc(n)%p%qc
         p_err = iv%tamdar_sfc(n)%p%error
         if (iv%tamdar_sfc(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error
         end if

         t    =  ob%tamdar_sfc(n)%t
         t_qc =  iv%tamdar_sfc(n)%t%qc
         if (iv%tamdar_sfc(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
         end if

         write(ounit, fmt = trim (fmt_each))                   &
            p,  p_qc,  p_err,                                  &
            speed,speed_qc,        speed_err  ,                &
            dir  , dir_qc,         dir_err  ,                  &
            iv%info(tamdar_sfc)%elv(n), zero, elv_error,             &
            t, t_qc, iv%tamdar_sfc(n)%t%error,                      &
            td, td_qc, td_error,                               &
            rh, rh_qc,  rh_error
      end do
   end if


   if (iv%info(ssmi_tb)%nlocal > 0) then
      call da_warning("da_write_filtered_obs.inc",1645, &     
        (/"Currently SSMI brightness temperature info is not written for filtered obs"/))
   end if

   close(ounit)
   call da_free_unit(ounit)
   ! Ensure other processors have written their temporary files
   call mpi_barrier(comm, ierr)
   if (rootproc) then
      call da_get_unit(filtered_obs_unit)
      write(unit=filename, fmt='(a,i2.2)') 'filtered_obs_', it
      open (unit=filtered_obs_unit, &
         file= filename ,form='formatted', status='replace', iostat=iost)
      if (iost /= 0) &
         call da_error("da_write_filtered_obs.inc",1661, (/"Cannot open filtered_obs "/))
      call da_count_filtered_obs(&
         grid%xb%ptop, grid%xb%map, grid%xb%ds,  &
         grid%moad_cen_lat, grid%stand_lon, grid%truelat1, grid%truelat2,                 &
         coarse_ix, coarse_jy, start_x, start_y)
      call da_final_write_filtered_obs
      close (filtered_obs_unit)
      call da_free_unit(filtered_obs_unit)
   end if

   if (trace_use) call da_trace_exit("da_write_filtered_obs")

end subroutine da_write_filtered_obs


subroutine da_write_modified_filtered_obs(grid, ob, iv, &
   coarse_ix, coarse_jy, start_x, start_y)

   !------------------------------------------------------------------------
   ! Purpose: Writes observations and innovation by WRFVAR
   !------------------------------------------------------------------------   
   !  Note: [cys] This subroutine is modified from da_write_filtered_obs 
   !  Note:                                                           
   !   (a) Information just needed for WRFVAR  is only written  
   !   (b) Both U and V should be of good quality (qc > obs_qc_pointer) 
   !   (c) Pressure & Temperature quality should be > fails_error_max (= -3) 
   !   (d) Since there is no check for missing T,P and Q to get Qs in   
   !       "da_tp_to_qs", RH data is recovered from Q if 
   !       q_error = missing_r and both P and T is not missing
   !   (e) Since currently Sfc_Obs_Correction changes quality of
   !       observations to "surface_correction (= 2)" even if it is bad
   !       quality obs (qc < fails_error_max), such obs (element) is 
   !       made missing. Otherwise the bad obs (element) will also get 
   !       assmiilated with "max_iv_check = .false."
   !   (f) AMV's from Geostationary and Polar orbiting satellite are
   !       seperated & used as profile
   !   (g) Currently only MODIS winds from Polar orbiting satellite are used
   !  Modifications:
   !
   !..........................................................................
   ! The Quality Controlled observations (QCed obs) means
   !
   !   i) The check against the Background;
   !
   !  ii) The model bottom/top check, i.e. the observations below the
   !      lowest model level (note, not bottom) and higher than the highest
   !      model level (note, not model top) are filled with missing values
   !      except the Sound data;
   !
   ! iii) For Sound data, the above-top/below-bottom data are not written out.
   !
   ! This can be used for Var with check_max_iv = .false., or for
   ! verification with  analysis_type = 'VERIFY'.
   !---------------------------------------------------------------------------   

   implicit none

   type (domain),     intent(in)  :: grid
   type (y_type),     intent(in)  :: ob      ! Observation structure.
   type (iv_type),    intent(inout) :: iv      ! O-B structure.
   integer, intent(in)            :: coarse_ix, coarse_jy
   real, intent(in)               :: start_x, start_y

   integer                      :: j, k, iost, nlevels

   integer            :: n, zero
   real               :: es, qs, speed, dir, rh, rh_error, uu, vv, zz, hh
   real               :: height, height_error, pr_error, slp_error
   real               :: thick, thick_error,dir_error,speed_error
   real               :: td, td_error , elv_error, ships_pr_error

   real               :: slp_inv,speed_err,dir_err, p,p_err
   real               :: t, t_err,ref, ref_err
   integer            :: thick_qc, t_qc, td_qc, rh_qc, height_qc
   integer            :: slp_qc, speed_qc, dir_qc,p_qc, ref_qc

   real               :: q, q_err
   integer            :: q_qc

   integer            :: ounit, nz, mnz
   character(len=filename_len)       :: filename

   character(len=120) :: fmt_each_inv = &
      '(3(f12.3,i4,f7.2),11x,3(f12.3,i4,f7.2),11x,1(f12.3,i4,f7.2),4f17.7)'

   if (trace_use) call da_trace_entry("da_write_modified_filtered_obs")

   ! Create noise and output:
   call da_get_unit(ounit)

   write(unit=filename, fmt='(a,i4.4)') 'filtered_obs.', myproc
   open(unit=ounit,file=trim(filename),form='formatted', & 
      status='replace', iostat=iost)
   if (iost /= 0) then
      call da_error("da_write_modified_filtered_obs.inc",84, &
         (/"Cannot open filtered observation file"//filename/))
   end if

   zero = 0
   ! Note: Currently in WRFVAR there is no active role of  
   !       "height" and "TD"        

   height       = missing_r
   height_qc    = missing_data
   height_error = xmiss    
   td           = missing_r 
   td_qc        = missing_data
   td_error     = xmiss        
   ! Currently following errors are fixed 
   rh_error     = 10.0       
   pr_error = 100.
   ships_pr_error = 160.
   slp_error = 200.
   elv_error = 6.0
   dir_error = 5.0
   dir_err   = 5.0

   !--------------------------------------------------------------------------
   !  Write Synop:
   !--------------------------------------------------------------------------

   if (iv%info(synop)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout, fmt='(A,I6)') ' writing filtered obs for synop ', &
            iv%info(synop)%nlocal
      end if

      do n = 1, iv%info(synop)%nlocal
         if (.not.iv%info(synop)%proc_domain(1,n)) cycle
         ! Guo .. not write out the data below the lowest model level:

         nlevels = iv%info(synop)%levels(n)
         if (iv%info(synop)%levels(n) > 1) then
            write(unit=stdout, fmt='(3A,I5,A)') &
               ' for SYNOP station ',iv%info(synop)%name(n),' got levels ',&
               iv%info(synop)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if
         write(ounit, fmt = fmt_info)                &
            iv%info(synop)%platform(n),    &
            iv%info(synop)%date_char(n),   &
            iv%info(synop)%name(n),        &
            nlevels,                      &
            iv%info(synop)%lat(1,n),         &
            iv%info(synop)%lon(1,n),         &
            iv%info(synop)%elv(n),         &
            iv%info(synop)%id(n)
         slp_inv= iv%info(synop)%slp(n)%inv
         slp_qc = iv%info(synop)%slp(n)%qc
         if (iv%synop(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(synop)%slp(n)%error, &
            iv%info(synop)%pw(n)%inv, iv%info(synop)%pw(n)%qc, iv%info(synop)%pw(n)%error

         ! get speed and direction
         if (iv%synop(n)%u%qc >= obs_qc_pointer .and. &
              iv%synop(n)%v%qc >= obs_qc_pointer) then
            uu = ob%synop(n)%u
            vv = ob%synop(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(synop)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%synop(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%synop(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
            uu        = missing_r
            vv        = missing_r
         end if

         ! get RH from Q & T    
         if (iv%synop(n)%q%qc >= obs_qc_pointer .and.  &
              abs(iv%synop(n)%q%error - missing_r) > 1. .and.  &
              abs(ob%synop(n)%t       - missing_r) > 1. .and.  &
              abs(ob%synop(n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%synop(n)%t, ob%synop(n)%p, es, qs)
            rh = (ob%synop(n)%q/qs) * 100.
            rh_qc = iv%synop(n)%q%qc
            q = ob%synop(n)%q * 1e3    !g/kg
            q_qc = iv%synop(n)%q%qc
            q_err = max(0.01, iv%synop(n)%q%error * 1e3) !g/kg
         else
            rh    = missing_r
            rh_qc = missing_data
            q     = missing_r
            q_qc  = missing_data
            q_err = xmiss
         end if

         if (rh_qc < 0) rh_qc = missing_data
         if ( q_qc < 0)  q_qc = missing_data

         p    = ob%synop(n)%p
         p_qc = iv%synop(n)%p%qc
         p_err = iv%synop(n)%p%error
         if (iv%synop(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error           
         end if

         t    =  ob%synop(n)%t
         t_qc =  iv%synop(n)%t%qc
         t_err =  iv%synop(n)%t%error
         if (iv%synop(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
            t_err = xmiss
         end if

         write(ounit, fmt = trim (fmt_each_inv))       &
            p,  p_qc,   p_err,                         &
            uu, speed_qc, speed_err ,                  &
            vv, speed_qc, speed_err ,                  &
            iv%info(synop)%elv(n), zero, elv_error,        &
            t, t_qc, t_err,                            &
            td, td_qc, td_error,                       &
            q, q_qc, q_err,                            &
            iv%synop(n)%u%inv, iv%synop(n)%v%inv, &
            iv%synop(n)%t%inv, iv%synop(n)%q%inv*1e3            ! q%inv in g/kg
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Metar:
   !--------------------------------------------------------------------------

   if (iv%info(metar)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for metar ',iv%info(metar)%nlocal
      end if

      do n = 1, iv%info(metar)%nlocal
         if (.not. iv%info(metar)%proc_domain(1,n)) cycle
         ! Do not write out the data below the lowest model level:

         nlevels = iv%info(metar)%levels(n)
         if (iv%info(metar)%levels(n)> 1) then
            write(stdout,fmt='(3A,I5,A)') &
               ' for METAR station ',iv%info(metar)%name(n),' got levels ',&
               iv%info(metar)%levels(n),' but wrote only one level'
            nlevels = 1
         end if

         write(ounit, fmt = fmt_info)                &
            iv%info(metar)%platform(n),    &
            iv%info(metar)%date_char(n),   &
            iv%info(metar)%name(n),        &
            nlevels,                      &
            iv%info(metar)%lat(1,n),         &
            iv%info(metar)%lon(1,n),         &
            iv%info(metar)%elv(n),         &
            iv%info(metar)%id(n)
         slp_inv= iv%info(metar)%slp(n)%inv
         slp_qc = iv%info(metar)%slp(n)%qc
         if(iv%metar(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(metar)%slp(n)%error, &
            iv%info(metar)%pw(n)%inv, iv%info(metar)%pw(n)%qc, iv%info(metar)%pw(n)%error

         ! get speed and direction
         if (iv%metar(n)%u%qc >= obs_qc_pointer .and. &
            iv%metar(n)%v%qc >= obs_qc_pointer) then
            uu = ob%metar(n)%u
            vv = ob%metar(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(metar)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%metar(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%metar(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
            uu        = missing_r
            vv        = missing_r
         end if

         ! get RH from Q & T    
         if (iv%metar(n)%q%qc >= obs_qc_pointer .and.  &
             abs(iv%metar(n)%q%error - missing_r) > 1. .and.  &
             abs(ob%metar(n)%t       - missing_r) > 1. .and.  &
             abs(ob%metar(n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%metar(n)%t, ob%metar(n)%p, es, qs)
            rh = (ob%metar(n)%q/qs) * 100.
            rh_qc = iv%metar(n)%q%qc
            q = ob%metar(n)%q * 1e3    !g/kg
            q_qc = iv%metar(n)%q%qc
            q_err = max(0.01, iv%metar(n)%q%error * 1e3) !g/kg
         else
            rh    = missing_r
            rh_qc = missing_data
            q     = missing_r
            q_qc  = missing_data
            q_err = xmiss
         end if
         if (rh_qc < 0) rh_qc = missing_data
         if ( q_qc < 0)  q_qc = missing_data

         p    = ob%metar(n)%p
         p_qc = iv%metar(n)%p%qc
         p_err = iv%metar(n)%p%error
         if (iv%metar(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error           
         end if
         t    =  ob%metar(n)%t
         t_qc =  iv%metar(n)%t%qc
         t_err =  iv%metar(n)%t%error
         if (iv%metar(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
            t_err = xmiss
         end if

         write(ounit, fmt = trim (fmt_each_inv))       &
            p,  p_qc,   p_err,                         &
            uu, speed_qc, speed_err ,                  &
            vv, speed_qc, speed_err ,                  &
            iv%info(metar)%elv(n), zero, elv_error,        &
            t, t_qc, t_err,                            &
            td, td_qc, td_error,                       &
            q, q_qc, q_err,                            &
            iv%metar(n)%u%inv, iv%metar(n)%v%inv, &
            iv%metar(n)%t%inv, iv%metar(n)%q%inv*1e3            ! q%inv in g/kg
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Ships:
   !-------------------------------------------------------------------------- 

   if (iv%info(ships)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for ships ',iv%info(ships)%nlocal
      end if

      do n = 1, iv%info(ships)%nlocal
         if (.not. iv%info(ships)%proc_domain(1,n)) cycle
               
         nlevels = iv%info(ships)%levels(n)
         if (iv%info(ships)%levels(n) > 1) then
            write(unit=stdout,fmt='(3A,I5,A)') &
               ' for SHIP station ',iv%info(ships)%name(n),' got levels ',&
               iv%info(ships)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if
         write(ounit, fmt = fmt_info)                &
            iv%info(ships)%platform(n),    &
            iv%info(ships)%date_char(n),   &
            iv%info(ships)%name(n),        &
            nlevels,               &
            iv%info(ships)%lat(1,n),       &
            iv%info(ships)%lon(1,n),       &
            iv%info(ships)%elv(n),         &
            iv%info(ships)%id(n)
            slp_inv= iv%info(ships)%slp(n)%inv
            slp_qc = iv%info(ships)%slp(n)%qc
         if (iv%ships(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(ships)%slp(n)%error, &
            iv%info(ships)%pw(n)%inv, iv%info(ships)%pw(n)%qc, iv%info(ships)%pw(n)%error

         ! get speed and direction
         if (iv%ships(n)%u%qc >= obs_qc_pointer .and. &
             iv%ships(n)%v%qc >= obs_qc_pointer) then
            uu = ob%ships(n)%u
            vv = ob%ships(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(ships)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%ships(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%ships(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
            uu        = missing_r
            vv        = missing_r
         end if

         ! get RH from Q & T    
         if (iv%ships(n)%q%qc >= obs_qc_pointer .and.  &
             abs(iv%ships(n)%q%error - missing_r) > 1. .and.  &
             abs(ob%ships(n)%t       - missing_r) > 1. .and.  &
             abs(ob%ships(n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%ships(n)%t, ob%ships(n)%p, es, qs)
            rh = (ob%ships(n)%q/qs) * 100.
            rh_qc = iv%ships(n)%q%qc
            q = ob%ships(n)%q * 1e3    !g/kg
            q_qc = iv%ships(n)%q%qc
            q_err = max(0.01, iv%ships(n)%q%error * 1e3) !g/kg
         else
            rh    = missing_r
            rh_qc = missing_data
            q     = missing_r
            q_qc  = missing_data
            q_err = xmiss
         end if
         if (rh_qc < 0) rh_qc = missing_data
         if ( q_qc < 0)  q_qc = missing_data

         p     = ob%ships(n)%p
         p_qc  = iv%ships(n)%p%qc
         p_err =  ships_pr_error
         if (iv%ships(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
         end if
         t    =  ob%ships(n)%t
         t_qc =  iv%ships(n)%t%qc
         t_err =  iv%ships(n)%t%error
         if (iv%ships(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
            t_err = xmiss
         end if

         write(ounit, fmt = trim (fmt_each_inv))       &
            p,  p_qc,   p_err,                         &
            uu, speed_qc, speed_err ,                  &
            vv, speed_qc, speed_err ,                  &
            iv%info(ships)%elv(n), zero, elv_error,        &
            t, t_qc, t_err,                            &
            td, td_qc, td_error,                       &
            q, q_qc, q_err,                            &
            iv%ships(n)%u%inv, iv%ships(n)%v%inv, &
            iv%ships(n)%t%inv, iv%ships(n)%q%inv*1e3            ! q%inv in g/kg
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Buoy :
   !--------------------------------------------------------------------------

   if (iv%info(buoy)%nlocal  > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for buoy  ',iv%info(buoy)%nlocal 
      end if

      do n = 1, iv%info(buoy)%nlocal
         if (.not. iv%info(buoy)%proc_domain(1,n)) cycle
                
         nlevels = iv%info(buoy)%levels(n)
         if (iv%info(buoy)%levels(n) > 1) then
            write(unit=stdout, fmt='(3A,I5,A)') &
               ' for BUOY  station ',iv%info(buoy)%name(n),' got levels ',&
               iv%info(buoy)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if

         write(ounit, fmt = fmt_info)                &
            iv%info(buoy)%platform(n),    &
            iv%info(buoy)%date_char(n),   &
            iv%info(buoy)%name(n),        &
            nlevels,         &
            iv%info(buoy)%lat(1,n),         &
            iv%info(buoy)%lon(1,n),         &
            iv%info(buoy)%elv(n),         &
            iv%info(buoy)%id(n)
         slp_inv= iv%info(buoy)%slp(n)%inv
         slp_qc = iv%info(buoy)%slp(n)%qc
         if (iv%buoy (n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(buoy)%slp(n)%error, &
            iv%info(buoy)%pw(n)%inv, iv%info(buoy)%pw(n)%qc, iv%info(buoy)%pw(n)%error

         ! get speed and direction
         if (iv%buoy (n)%u%qc >= obs_qc_pointer .and. &
              iv%buoy (n)%v%qc >= obs_qc_pointer) then
            uu = ob%buoy (n)%u
            vv = ob%buoy (n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(buoy)%lon(1,n), &
               convert_uv2fd)             
            speed_qc  = iv%buoy (n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%buoy (n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
            uu        = missing_r
            vv        = missing_r
         end if

         ! get RH from Q & T    
         if (iv%buoy (n)%q%qc >= obs_qc_pointer .and.  &
             abs(iv%buoy (n)%q%error - missing_r) > 1. .and.  &
             abs(ob%buoy (n)%t       - missing_r) > 1. .and.  &
             abs(ob%buoy (n)%p       - missing_r) > 1.   ) then 
            call da_tp_to_qs(ob%buoy (n)%t, ob%buoy (n)%p, es, qs)
            rh = (ob%buoy (n)%q/qs) * 100.
            rh_qc = iv%buoy (n)%q%qc
            q = ob%buoy (n)%q * 1e3    !g/kg
            q_qc = iv%buoy (n)%q%qc
            q_err = max(0.01, iv%buoy (n)%q%error * 1e3) !g/kg
         else
            rh    = missing_r
           rh_qc = missing_data
            q     = missing_r
            q_qc  = missing_data
            q_err = xmiss
         end if
         if (rh_qc < 0) rh_qc = missing_data
         if ( q_qc < 0)  q_qc = missing_data

         p     = ob%buoy (n)%p
         p_qc  = iv%buoy (n)%p%qc
         p_err = iv%buoy (n)%p%error
         if (iv%buoy (n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error           
         end if
         t    =  ob%buoy (n)%t
         t_qc =  iv%buoy (n)%t%qc
         t_err =  iv%buoy (n)%t%error
         if (iv%buoy (n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
            t_err = xmiss
         end if

         write(ounit, fmt = trim (fmt_each_inv))       &
            p,  p_qc,   p_err,                         &
            uu, speed_qc, speed_err ,                  &
            vv, speed_qc, speed_err ,                  &
            iv%info(buoy)%elv(n), zero, elv_error,        &
            t, t_qc, t_err,                            &
            td, td_qc, td_error,                       &
            q, q_qc, q_err,                            &
            iv%buoy(n)%u%inv, iv%buoy(n)%v%inv, &
            iv%buoy(n)%t%inv, iv%buoy(n)%q%inv*1e3            ! q%inv in g/kg
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Geo. AMVs Obs:
   !--------------------------------------------------------------------------

   if (iv%info(geoamv)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for geoamv ',iv%info(geoamv)%nlocal
      end if

      do n = 1, iv%info(geoamv)%nlocal               
         if (.not. iv%info(geoamv)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(geoamv)%platform(n),    &
            iv%info(geoamv)%date_char(n),   &
            iv%info(geoamv)%name(n),        &
            iv%info(geoamv)%levels(n),      &
            iv%info(geoamv)%lat(1,n),         &
            iv%info(geoamv)%lon(1,n),         &
            iv%info(geoamv)%elv(n),         &
            iv%info(geoamv)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(geoamv)%slp(n)%inv, iv%info(geoamv)%slp(n)%qc, iv%info(geoamv)%slp(n)%error, &
            iv%info(geoamv)%pw(n)%inv, iv%info(geoamv)%pw(n)%qc, iv%info(geoamv)%pw(n)%error

         do k = 1, iv%info(geoamv)%levels(n)

            ! get speed and direction
            if (iv%geoamv(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%geoamv(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%geoamv(n)%u(k)
               vv = ob%geoamv(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(geoamv)%lon(1,n), convert_uv2fd)                        
               speed_qc  = iv%geoamv(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%geoamv(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc  = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%geoamv(n)%p(k), speed_qc, slp_error,  &
               uu, speed_qc, speed_err,                 &
               vv, speed_qc, speed_err,                 &
               missing_r, standard_atmosphere, xmiss,   &
               missing_r, zero_t_td, xmiss,             & 
               missing_r, zero_t_td, xmiss,             & 
               missing_r, zero_t_td, xmiss,             &
               iv%geoamv(n)%u(k)%inv, iv%geoamv(n)%v(k)%inv, &
               xmiss, xmiss
         end do   
      end do   
   end if

   !--------------------------------------------------------------------------
   ! Write Polar AMVs Obs:
   !--------------------------------------------------------------------------

   if (iv%info(polaramv)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for polaramv ',iv%info(polaramv)%nlocal
      end if

      do n = 1, iv%info(polaramv)%nlocal               
         if (.not. iv%info(polaramv)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(polaramv)%platform(n),    &
            iv%info(polaramv)%date_char(n),   &
            iv%info(polaramv)%name(n),        &
            iv%info(polaramv)%levels(n),      &
            iv%info(polaramv)%lat(1,n),         &
            iv%info(polaramv)%lon(1,n),         &
            iv%info(polaramv)%elv(n),         &
            iv%info(polaramv)%id(n)

         write(ounit, fmt = fmt_srfc) &
            iv%info(polaramv)%slp(n)%inv, iv%info(polaramv)%slp(n)%qc,   &
            iv%info(polaramv)%slp(n)%error,                              &
            iv%info(polaramv)%pw(n)%inv, iv%info(polaramv)%pw(n)%qc,     &
            iv%info(polaramv)%pw(n)%error

         do k = 1, iv%info(polaramv)%levels(n)
            ! get speed and direction
            if (iv%polaramv(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%polaramv(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%polaramv(n)%u(k)
               vv = ob%polaramv(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(polaramv)%lon(k,n), convert_uv2fd)                        
               speed_qc  = iv%polaramv(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%polaramv(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc  = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%polaramv(n)%p(k), speed_qc, slp_error,  &
               uu, speed_qc, speed_err,                 &
               vv, speed_qc, speed_err,                 &
               missing_r, standard_atmosphere, xmiss,   &
               missing_r, zero_t_td, xmiss,             & 
               missing_r, zero_t_td, xmiss,             & 
               missing_r, zero_t_td, xmiss,             &
               iv%polaramv(n)%u(k)%inv, iv%polaramv(n)%v(k)%inv, &
               xmiss, xmiss
         end do   
      end do   
   end if

   !--------------------------------------------------------------------------
   ! Write Sound 
   !--------------------------------------------------------------------------

   if (iv%info(sound)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(stdout,fmt='(A,I5)') &
            ' writing filtered obs for sound  ',iv%info(sound)%nlocal 
      end if
      mnz = 0

      do n = 1, iv%info(sound)%nlocal               
         if (.not. iv%info(sound)%proc_domain(1,n)) cycle
         if (iv%info(sonde_sfc)%platform(n) /= iv%info(sound)%platform(n) .or. &
             iv%info(sonde_sfc)%lat(1,n)    /= iv%info(sound)%lat(1,n)      .or. &
             iv%info(sonde_sfc)%lon(1,n)    /= iv%info(sound)%lon(1,n)      .or. &
             iv%info(sonde_sfc)%elv(n)      /= iv%info(sound)%elv(n)      .or. &
             iv%info(sonde_sfc)%id(n)       /= iv%info(sound)%id(n)      ) then
            write(unit=stderr,fmt='(A)')' Surface Sound Details:            '
            write(unit=stderr, fmt=fmt_info) &
               iv%info(sonde_sfc)%platform(n),    &
               iv%info(sonde_sfc)%date_char(n),   &
               iv%info(sonde_sfc)%name(n),        &
               iv%info(sonde_sfc)%levels(n),      &
               iv%info(sonde_sfc)%lat(1,n),         &
               iv%info(sonde_sfc)%lon(1,n),         &
               iv%info(sonde_sfc)%elv(n),         &
               iv%info(sonde_sfc)%id(n)

            write(unit=stderr,fmt='(A)') ' Upper level  Details: '
            write(unit=stderr, fmt=fmt_info) &
               iv%info(sound)%platform(n),    &
               iv%info(sound)%date_char(n),   &
               iv%info(sound)%name(n),        &
               iv%info(sound)%levels(n),      &
               iv%info(sound)%lat(1,n),         &
               iv%info(sound)%lon(1,n),         &
               iv%info(sound)%elv(n),         &
               iv%info(sound)%id(n)
            call da_error("da_write_modified_filtered_obs.inc",717, &
               (/"Mismatch for Sound surface and upper air report"/))
         end if
         nz = 0

         if (nz > 0) then
            mnz = mnz + 1
         end if

         nz = iv%info(sound)%levels(n) + 1 - nz
         if (nz < 1) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(sound)%platform(n),    &
            iv%info(sound)%date_char(n),   &
            iv%info(sound)%name(n),        &
            nz,                           &
            iv%info(sound)%lat(1,n),         &
            iv%info(sound)%lon(1,n),         &
            iv%info(sound)%elv(n),         &
            iv%info(sound)%id(n)
            ! iv%info(sound)%levels(n) + 1,  &
         slp_inv= iv%info(sound)%slp(n)%inv
         slp_qc = iv%info(sound)%slp(n)%qc
         if (iv%sonde_sfc(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(sound)%slp(n)%error,  &
            iv%info(sound)%pw(n)%inv, iv%info(sound)%pw(n)%qc, iv%info(sound)%pw(n)%error
 
         nz = 0
         j  = 0

         ! Not below the first model level
            j = j + 1
            ! First write surface level information
 
            ! get speed and direction
            if (iv%sonde_sfc(n)%u%qc >= obs_qc_pointer .and. &
                iv%sonde_sfc(n)%v%qc >= obs_qc_pointer) then
               uu = ob%sonde_sfc(n)%u
               vv = ob%sonde_sfc(n)%v
               call da_ffdduv(speed, dir, uu, vv, iv%info(sound)%lon(1,n), &
                  convert_uv2fd)             
               speed_qc  = iv%sonde_sfc(n)%u%qc
               dir_qc    = speed_qc
               speed_err = iv%sonde_sfc(n)%u%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

            ! get RH from Q & T    
            if (iv%sonde_sfc(n)%q%qc >= obs_qc_pointer .and.  &
                abs(iv%sonde_sfc(n)%q%error - missing_r) > 1. .and.  &
                abs(ob%sonde_sfc(n)%t       - missing_r) > 1. .and.  &
                abs(ob%sonde_sfc(n)%p       - missing_r) > 1.   ) then 
               call da_tp_to_qs(ob%sonde_sfc(n)%t, ob%sonde_sfc(n)%p, es, qs)
               rh = (ob%sonde_sfc(n)%q/qs) * 100.
               rh_qc = iv%sonde_sfc(n)%q%qc
               q = ob%sonde_sfc(n)%q * 1e3    !g/kg
               q_qc = iv%sonde_sfc(n)%q%qc
               q_err = max(0.01, iv%sonde_sfc(n)%q%error * 1e3) !g/kg
             else
                rh    = missing_r
                rh_qc = missing_data
                q     = missing_r
                q_qc  = missing_data
                q_err = xmiss
             end if
            if (rh_qc < 0) rh_qc = missing_data
            if ( q_qc < 0)  q_qc = missing_data

            p     = ob%sonde_sfc(n)%p
            p_qc  = iv%sonde_sfc(n)%p%qc
            p_err = iv%sonde_sfc(n)%p%error
            if (iv%sonde_sfc(n)%p%qc <= fails_error_max) then
               p_qc  = missing_data
               p     = missing_r
               p_err   = pr_error           
            end if

            t     =  ob%sonde_sfc(n)%t
            t_qc  =  iv%sonde_sfc(n)%t%qc
            t_err =  iv%sonde_sfc(n)%t%error
            if (iv%sonde_sfc(n)%t%qc <= fails_error_max) then
               t_qc  = missing_data
               t     = missing_r
               t_err = xmiss
            end if

            write(ounit, fmt = trim (fmt_each_inv))       &
               p,  p_qc,   p_err,                         &
               uu, speed_qc, speed_err ,                  &
               vv, speed_qc, speed_err ,                  &
               iv%sonde_sfc(n)%h, zero, elv_error,        &
               t, t_qc, t_err,                            &
               td, td_qc, td_error,                       &
               q, q_qc, q_err,                            &
               iv%sonde_sfc(n)%u%inv, iv%sonde_sfc(n)%v%inv, &
               iv%sonde_sfc(n)%t%inv, iv%sonde_sfc(n)%q%inv*1e3            ! q%inv in g/kg
                
         do k = 1, iv%info(sound)%levels(n)

               j  = j + 1

               ! get speed and direction
               if (iv%sound(n)%u(k)%qc >= obs_qc_pointer .and. &
                   iv%sound(n)%v(k)%qc >= obs_qc_pointer) then  
                  uu = ob%sound(n)%u(k)
                  vv = ob%sound(n)%v(k)
                  call da_ffdduv(speed, dir, uu, vv, iv%info(sound)%lon(1,n), convert_uv2fd)                        
                  speed_qc  = iv%sound(n)%u(k)%qc
                  dir_qc    = speed_qc
                  speed_err = iv%sound(n)%u(k)%error
               else
                  speed     = missing_r
                  speed_qc  = missing_data
                  dir       = missing_r
                  dir_qc    = missing_data
                  speed_err = xmiss
                  uu        = missing_r
                  vv        = missing_r
               end if

               ! get RH from Q & T    
               if (iv%sound(n)%q(k)%qc >= obs_qc_pointer .and.  &
                   abs(iv%sound(n)%q(k)%error - missing_r) > 1. .and.  &
                   abs(ob%sound(n)%t(k)       - missing_r) > 1. .and.  &
                   abs(iv%sound(n)%p(k)       - missing_r) > 1.   ) then 
                  call da_tp_to_qs(ob%sound(n)%t(k), iv%sound(n)%p(k), es, qs)
                  rh = (ob%sound(n)%q(k)/qs) * 100.
                  rh_qc = iv%sound(n)%q(k)%qc
                  q = ob%sound(n)%q(k) * 1e3    !g/kg
                  q_qc = iv%sound(n)%q(k)%qc
                  q_err = max(0.01, iv%sound(n)%q(k)%error * 1e3) !g/kg
               else
                  rh    = missing_r
                  rh_qc = missing_data
                  q     = missing_r
                  q_qc  = missing_data
                  q_err = xmiss
               end if
               if (rh_qc < 0) rh_qc = missing_data
               if ( q_qc < 0)  q_qc = missing_data

               t    =  ob%sound(n)%t(k)
               t_qc =  iv%sound(n)%t(k)%qc
               t_err= iv%sound(n)%t(k)%error
               if (iv%sound(n)%t(k)%qc <= fails_error_max) then
                  t_qc = missing_data
                  t   = missing_r
               end if

               write(ounit, fmt = trim (fmt_each_inv))       &
                  iv%sound(n)%p(k), zero,    pr_error,       &
                  uu, speed_qc, speed_err ,                  &
                  vv, speed_qc, speed_err ,                  &
                  iv%sound(n)%h(k), zero, elv_error,         &
                  t, t_qc, t_err,                            &
                  td, td_qc, td_error,                       &
                  q, q_qc, q_err,                            &
                  iv%sound(n)%u(k)%inv, iv%sound(n)%v(k)%inv, &
                  iv%sound(n)%t(k)%inv, iv%sound(n)%q(k)%inv*1e3   !q%inv g/kg
         end do
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Airep: 
   !--------------------------------------------------------------------------

   if (iv%info(airep)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for airep  ',iv%info(airep)%nlocal 
      end if

      do n = 1, iv%info(airep)%nlocal               
         if (.not. iv%info(airep)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(airep)%platform(n),    &
            iv%info(airep)%date_char(n),   &
            iv%info(airep)%name(n),        &
            iv%info(airep)%levels(n),      &
            iv%info(airep)%lat(1,n),         &
            iv%info(airep)%lon(1,n),         &
            iv%info(airep)%elv(n),         &
            iv%info(airep)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(airep)%slp(n)%inv, iv%info(airep)%slp(n)%qc, &
            iv%info(airep)%slp(n)%error,                         &
            iv%info(airep)%pw(n)%inv, iv%info(airep)%pw(n)%qc,   &
            iv%info(airep)%pw(n)%error

         do k = 1, iv%info(airep)%levels(n)

            ! get speed and direction
            if (iv%airep(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%airep(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%airep(n)%u(k)
               vv = ob%airep(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(airep)%lon(k,n), convert_uv2fd)             
               speed_qc  = iv%airep(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%airep(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

          ! get RH from Q & T
            if (iv%airep(n)%q(k)%qc >= obs_qc_pointer .and.  &
                abs(iv%airep(n)%q(k)%error - missing_r) > 1. .and.  &
                abs(ob%airep(n)%t(k)       - missing_r) > 1. .and.  &
                abs(iv%airep(n)%p(k)       - missing_r) > 1.   ) then
                call da_tp_to_qs(ob%airep(n)%t(k), iv%airep(n)%p(k), es, qs)
                rh = (ob%airep(n)%q(k)/qs) * 100.
                rh_qc = iv%airep(n)%q(k)%qc
                q = ob%airep(n)%q(k) * 1e3    !g/kg
                q_qc = iv%airep(n)%q(k)%qc
                q_err = max(0.01, iv%airep(n)%q(k)%error * 1e3) !g/kg
            else
                rh    = missing_r
                rh_qc = missing_data
                q     = missing_r
                q_qc  = missing_data
                q_err = xmiss
            end if
            if (rh_qc < 0) rh_qc = missing_data
            if ( q_qc < 0)  q_qc = missing_data

            t    =  ob%airep(n)%t(k)
            t_qc =  iv%airep(n)%t(k)%qc
            t_err= iv%airep(n)%t(k)%error
            if (iv%airep(n)%t(k)%qc <= fails_error_max) then
               t_qc = missing_data
               t   = missing_r
               t_err = xmiss
            end if

            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%airep(n)%p(k),zero , pr_error,          &
               uu, speed_qc, speed_err ,                  &
               vv, speed_qc, speed_err ,                  &
               iv%airep(n)%h(k), zero, elv_error,         &
               t, t_qc, t_err,                            &
               td, td_qc, td_error,                       &
               q, q_qc, q_err,                            &
               iv%airep(n)%u(k)%inv, iv%airep(n)%v(k)%inv, &
               iv%airep(n)%t(k)%inv, iv%airep(n)%q(k)%inv*1e3
         end do 
      end do     
   end if

   !--------------------------------------------------------------------------
   ! Write Pilot:
   !--------------------------------------------------------------------------

   if (iv%info(pilot)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
           ' writing filtered obs for pilot  ',iv%info(pilot)%nlocal 
      end if

      do n = 1, iv%info(pilot)%nlocal               
         if (.not. iv%info(pilot)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(pilot)%platform(n),    &
            iv%info(pilot)%date_char(n),   &
            iv%info(pilot)%name(n),        &
            iv%info(pilot)%levels(n),      &
            iv%info(pilot)%lat(1,n),         &
            iv%info(pilot)%lon(1,n),         &
            iv%info(pilot)%elv(n),         &
            iv%info(pilot)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(pilot)%slp(n)%inv, iv%info(pilot)%slp(n)%qc, &
            iv%info(pilot)%slp(n)%error,                         &
            iv%info(pilot)%pw(n)%inv, iv%info(pilot)%pw(n)%qc,   &
            iv%info(pilot)%pw(n)%error

         do k = 1, iv%info(pilot)%levels(n)

            ! get speed and direction
            if (iv%pilot(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%pilot(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%pilot(n)%u(k)
               vv = ob%pilot(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(pilot)%lon(k,n), &
                  convert_uv2fd)                        
               speed_qc  = iv%pilot(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%pilot(n)%u(k)%error
            else
               speed     = missing_r
            speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%pilot(n)%p(k), zero, pr_error,  &
               uu, speed_qc, speed_err,                 &
               vv, speed_qc, speed_err,                 &
               missing_r, standard_atmosphere, xmiss,   &
               missing_r, missing_data, xmiss,             & 
               missing_r, missing_data, xmiss,             & 
               missing_r, missing_data, xmiss,             &
               iv%pilot(n)%u(k)%inv, iv%pilot(n)%v(k)%inv, &
               xmiss, xmiss
         end do 
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write satem:          
   ! Note: The satem ref_p is put in the SLP position in OBSPROC data stream.
   !--------------------------------------------------------------------------

   if (iv%info(satem)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for satem  ',iv%info(satem)%nlocal 
      end if

      do n = 1, iv%info(satem)%nlocal           
         if (.not. iv%info(satem)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(satem)%platform(n),    &
            iv%info(satem)%date_char(n),   &
            iv%info(satem)%name(n),        &
            iv%info(satem)%levels(n),      &
            iv%info(satem)%lat(1,n),         &
            iv%info(satem)%lon(1,n),         &
            iv%info(satem)%elv(n),         &
            iv%info(satem)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(satem)%slp(n)%inv, iv%info(satem)%slp(n)%qc, &
            iv%info(satem)%slp(n)%error,                         &
            iv%info(satem)%pw(n)%inv, iv%info(satem)%pw(n)%qc,   &
            iv%info(satem)%pw(n)%error

         do k = 1, iv%info(satem)%levels(n) 
            thick_qc = iv%satem(n)%thickness(k)%qc
            thick = iv%satem(n)%org_thickness(k)%inv   
            thick_error = iv%satem(n)%org_thickness(k)%error
            if (iv%satem(n)%thickness(k)%qc < obs_qc_pointer) &
               thick_qc = missing_data

            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%satem(n)%p(k),zero, slp_error,                 &
               missing_r, missing_data, xmiss,                   &
               missing_r, missing_data, xmiss,                   &
               missing_r, standard_atmosphere, xmiss,            &
               thick, thick_qc, thick_error,                     &
               missing_r, missing_data, xmiss,                   &
               missing_r, missing_data, xmiss,                   &
               xmiss, xmiss, xmiss, xmiss
         end do
      end do
   end if

   !--------------------------------------------------------------------------
   ! Write Qscat: 
   !--------------------------------------------------------------------------

   if (iv%info(qscat)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for qscat  ',iv%info(qscat)%nlocal 
      end if
 
      do n = 1, iv%info(qscat)%nlocal               
         if (.not. iv%info(qscat)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(qscat)%platform(n),    &
            iv%info(qscat)%date_char(n),   &
            iv%info(qscat)%name(n),        &
            iv%info(qscat)%levels(n),      &
            iv%info(qscat)%lat(1,n),         &
            iv%info(qscat)%lon(1,n),         &
            iv%info(qscat)%elv(n),         &
            iv%info(qscat)%id(n)

         write(ounit, fmt = fmt_srfc)  &
            iv%info(qscat)%slp(n)%inv, iv%info(qscat)%slp(n)%qc, &
            iv%info(qscat)%slp(n)%error,                         &
            iv%info(qscat)%pw(n)%inv, iv%info(qscat)%pw(n)%qc,   &
            iv%info(qscat)%pw(n)%error

         ! get speed and direction
         if (iv%qscat(n)%u%qc >= obs_qc_pointer .and. &
             iv%qscat(n)%v%qc >= obs_qc_pointer) then  
            uu = ob%qscat(n)%u
            vv = ob%qscat(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(qscat)%lon(1,n), &
               convert_uv2fd)                        
            speed_error = iv%qscat(n)%u%error
            speed_qc  = iv%qscat(n)%u%qc
            dir_qc    = speed_qc
         else
            speed = missing_r
            dir   = missing_r
            speed_error = xmiss
            speed_qc = missing_data
            dir_qc   = missing_data
            uu        = missing_r
            vv        = missing_r
         end if

         write(ounit, fmt = trim (fmt_each_inv))       &
            iv%qscat(n)%h,iv%qscat(n)%u%qc,slp_error, &  
            uu, speed_qc, speed_error,               &
            vv, speed_qc, speed_error,               &
            iv%info(qscat)%elv(n), zero, elv_error,    &
            missing_r, zero_t_td, xmiss,             & 
            missing_r, zero_t_td, xmiss,             & 
            missing_r, zero_t_td, xmiss,             &
            iv%qscat(n)%u%inv, iv%qscat(n)%v%inv,    &
            xmiss, xmiss

      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write Profiler: 
   !--------------------------------------------------------------------------
   
   if (iv%info(profiler)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for profiler  ',iv%info(profiler)%nlocal 
      end if

      do n = 1, iv%info(profiler)%nlocal            
         if (.not. iv%info(profiler)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost)    &
            iv%info(profiler)%platform(n),    &
            iv%info(profiler)%date_char(n),   &
            iv%info(profiler)%name(n),        &
            iv%info(profiler)%levels(n),      &
            iv%info(profiler)%lat(1,n),         &
            iv%info(profiler)%lon(1,n),         &
            iv%info(profiler)%elv(n),         &
            iv%info(profiler)%id(n)
         write(ounit, fmt = fmt_srfc) &
            iv%info(profiler)%slp(n)%inv, iv%info(profiler)%slp(n)%qc, &
            iv%info(profiler)%slp(n)%error,                            &
            iv%info(profiler)%pw(n)%inv, iv%info(profiler)%pw(n)%qc,   &
            iv%info(profiler)%pw(n)%error

         do k = 1, iv%info(profiler)%levels(n)

            ! get speed and direction
            if (iv%profiler(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%profiler(n)%v(k)%qc >= obs_qc_pointer) then  
               uu = ob%profiler(n)%u(k)
               vv = ob%profiler(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(profiler)%lon(k,n), &
                  convert_uv2fd)
               speed_error = iv%profiler(n)%u(k)%error
            else
               speed = missing_r
               dir   = missing_r
               speed_error = 1.1
               speed_qc = missing_data
               dir_qc   = missing_data
               uu        = missing_r
               vv        = missing_r
           end if

           write(ounit, fmt = trim (fmt_each_inv))       &
              iv%profiler(n)%p(k),zero,pr_error,        &
              uu, speed_qc, speed_error,               &
              vv, speed_qc, speed_error,               &
              iv%profiler(n)%h(k), zero, elv_error, &
              missing_r, missing_data, xmiss,           & 
              missing_r, missing_data, xmiss,           & 
              missing_r, missing_data, xmiss,           &
              iv%profiler(n)%u(k)%inv, iv%profiler(n)%v(k)%inv, &
              xmiss, xmiss

         end do 
      end do 
   end if


   ! [10] Write SSMT1:          
   if (iv%info(ssmt1)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for ssmt1 ', iv%info(ssmt1)%nlocal 
      end if

      do n = 1, iv%info(ssmt1)%nlocal             
         if (.not. iv%info(ssmt1)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(ssmt1)%platform(n),    &
            iv%info(ssmt1)%date_char(n),   &
            iv%info(ssmt1)%name(n),        &
            iv%info(ssmt1)%levels(n),      &
            iv%info(ssmt1)%lat(1,n),         &
            iv%info(ssmt1)%lon(1,n),         &
            iv%info(ssmt1)%elv(n),         &
            iv%info(ssmt1)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(ssmt1)%slp(n)%inv, iv%info(ssmt1)%slp(n)%qc, &
            iv%info(ssmt1)%slp(n)%error,                         &
            iv%info(ssmt1)%pw(n)%inv, iv%info(ssmt1)%pw(n)%qc,   &
            iv%info(ssmt1)%pw(n)%error

         do k = 1, iv%info(ssmt1)%levels(n)
            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%ssmt1(n)%p(k),zero,  slp_error,                          &
               missing_r, missing, missing_r,                              &
               missing_r, missing, missing_r,                              &
               height,          height_qc,height_error,                    &
               ob%ssmt1(n)%t(k),iv%ssmt1(n)%t(k)%qc,iv%ssmt1(n)%t(k)%error,&
               td, td_qc, td_error,                                        &   
               missing_r, missing, missing_r,                              &
               xmiss, xmiss,                                               &
               iv%ssmt1(n)%t(k)%inv, xmiss
         end do 
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write SSMT2:    
   !--------------------------------------------------------------------------      

   if (iv%info(ssmt2)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for ssmt2 ', iv%info(ssmt2)%nlocal 
      end if

      do n = 1, iv%info(ssmt2)%nlocal              
         if (.not. iv%info(ssmt2)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(ssmt2)%platform(n),    &
            iv%info(ssmt2)%date_char(n),   &
            iv%info(ssmt2)%name(n),        &
            iv%info(ssmt2)%levels(n),      &
            iv%info(ssmt2)%lat(1,n),         &
            iv%info(ssmt2)%lon(1,n),         &
            iv%info(ssmt2)%elv(n),         &
            iv%info(ssmt2)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(ssmt2)%slp(n)%inv, iv%info(ssmt2)%slp(n)%qc, &
            iv%info(ssmt2)%slp(n)%error,                         &
            iv%info(ssmt2)%pw(n)%inv, iv%info(ssmt2)%pw(n)%qc,   &
            iv%info(ssmt2)%pw(n)%error

         do k = 1, iv%info(ssmt2)%levels(n)
            write(ounit, fmt = trim (fmt_each_inv))       &
            iv%ssmt2(n)%p(k),zero,  slp_error,                &
            missing_r, missing, missing_r,                    &
            missing_r, missing, missing_r,                    &
            height,          height_qc,height_error,          &
            missing_r, missing, missing_r,                    &
            td, td_qc, td_error,                              &   
            ob%ssmt2(n)%rh(k),iv%ssmt2(n)%rh(k)%qc,rh_error,  &
            xmiss, xmiss, xmiss, xmiss                           !cys_note: ssmt2 not usable
         end do  
      end do 
   end if

   !--------------------------------------------------------------------------
   !  Write Gpspw:          
   !--------------------------------------------------------------------------

   if (iv%info(gpspw)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for gpspw ', iv%info(gpspw)%nlocal 
      end if

      do n = 1, iv%info(gpspw)%nlocal
         if (.not. iv%info(gpspw)%proc_domain(1,n)) cycle
         ! hardwired the # of levels = 0 for GPSPW (YRG 03/02/2006)
         iv%info(gpspw)%levels(n) = 0
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(gpspw)%platform(n),    &
            iv%info(gpspw)%date_char(n),   &
            iv%info(gpspw)%name(n),        &
            iv%info(gpspw)%levels(n),      &
            iv%info(gpspw)%lat(1,n),         &
            iv%info(gpspw)%lon(1,n),         &
            iv%info(gpspw)%elv(n),         &
            iv%info(gpspw)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(gpspw)%slp(n)%inv, iv%info(gpspw)%slp(n)%qc, &
            iv%info(gpspw)%slp(n)%error,                         &
            iv%info(gpspw)%pw(n)%inv, iv%info(gpspw)%pw(n)%qc,   &
            iv%info(gpspw)%pw(n)%error
      end do 
   end if

   ! Write Gps Refractivity:
   if (iv%info(gpsref)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for gpsref ', iv%info(gpsref)%nlocal
      end if

      do n = 1, iv%info(gpsref)%nlocal
         if (.not. iv%info(gpsref)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(gpsref) %platform(n),    &
            iv%info(gpsref) %date_char(n),   &
            iv%info(gpsref) %name(n),        &
            iv%info(gpsref) %levels(n),      &
            iv%info(gpsref) %lat(1,n),         &
            iv%info(gpsref) %lon(1,n),         &
            iv%info(gpsref) %elv(n),         &
            iv%info(gpsref) %id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(gpsref)%slp(n)%inv, iv%info(gpsref)%slp(n)%qc, &
            iv%info(gpsref)%slp(n)%error,                          &
            iv%info(gpsref)%pw(n)%inv, iv%info(gpsref)%pw(n)%qc,   &
            iv%info(gpsref)%pw(n)%error

         do k = 1, iv%info(gpsref)%levels(n)

            ! get RH from Q & T    
            if (iv%gpsref(n)%q(k)%qc >= obs_qc_pointer .and.  &
                abs(iv%gpsref(n)%q(k)%error - missing_r) > 1. .and.  &
                abs(ob%gpsref(n)%t(k)       - missing_r) > 1. .and.  &
                abs(ob%gpsref(n)%p(k)       - missing_r) > 1.   ) then 
               call da_tp_to_qs(ob%gpsref(n)%t(k), ob%gpsref(n)%p(k), es, qs)
               rh = (ob%gpsref(n)%q(k)/qs) * 100.
               rh_qc = iv%gpsref(n)%q(k)%qc
               q = ob%gpsref(n)%q(k) * 1e3    !g/kg
               q_qc = iv%gpsref(n)%q(k)%qc
               q_err = max(0.01, iv%gpsref(n)%q(k)%error * 1e3) !g/kg
            else
               rh    = missing_r
               rh_qc = missing_data
               q     = missing_r
               q_qc  = missing_data
               q_err = xmiss
            end if
            if ( q_qc < 0)  q_qc = missing_data

            p    = ob%gpsref(n)%p(k)
            p_qc = iv%gpsref(n)%p(k)%qc 
            if (iv%gpsref(n)%p(k)%qc <= fails_error_max) then
               p    = missing_r     
               p_qc = missing_data  
            end if

            t    = ob%gpsref(n)%t(k)
            t_qc  = iv%gpsref(n)%t(k)%qc
            t_err = iv%gpsref(n)%t(k)%error   
            if (iv%gpsref(n)%t(k)%qc <= fails_error_max) then
               t    = missing_r     
               t_qc = missing_data  
               t_err =  xmiss
            end if

            ref = ob%gpsref(n)%ref(k)
            ref_qc = iv%gpsref(n)%ref(k)%qc
            ref_err = iv%gpsref(n)%ref(k)%error   
            if (iv%gpsref(n)%ref(k)%qc <= fails_error_max) then
               ref    = missing_r     
               ref_qc = missing_data  
            end if


!
!  For GPSREF, convert the geopotential height to geometric height, then write
!  out to filtered_obs (YRG, 06/15/2011):
            hh = iv%gpsref(n)%h(k) / 1000.0
            call da_geo2msl1 (hh, iv%info(gpsref)%lat(1,n), &
                              iv%info(gpsref)%lon(1,n), zz)
! ............................................................................

            write(ounit, fmt = trim (fmt_each_inv))         &
               p,  p_qc, slp_error,                         &
               missing_r, missing_data, xmiss,              &
               missing_r, missing_data, xmiss,              &
               zz*1000.0 , zero, height_error,              &
               t, t_qc,  t_err,                             &   
               ref, ref_qc,  ref_err,                       &   
               q, q_qc, q_err,                              &
               xmiss, iv%gpsref(n)%t(k)%inv,                &
               iv%gpsref(n)%ref(k)%inv, missing_r
         end do 
      end do 
   end if

   !--------------------------------------------------------------------------
   ! Write AIRS  retrievals:
   !--------------------------------------------------------------------------

   if (iv%info(airsr)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for AIRS retrievals ',iv%info(airsr)%nlocal 
      end if

      do n = 1, iv%info(airsr)%nlocal               
         if (.not. iv%info(airsr)%proc_domain(1,n)) cycle
         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(airsr)%platform(n),    &
            iv%info(airsr)%date_char(n),   &
            iv%info(airsr)%name(n),        &
            iv%info(airsr)%levels(n),      &
            iv%info(airsr)%lat(1,n),         &
            iv%info(airsr)%lon(1,n),         &
            iv%info(airsr)%elv(n),         &
            iv%info(airsr)%id(n)
         write(ounit, fmt = fmt_srfc) &
            iv%info(airsr)%slp(n)%inv, iv%info(airsr)%slp(n)%qc, &
            iv%info(airsr)%slp(n)%error,                         &
            iv%info(airsr)%pw(n)%inv, iv%info(airsr)%pw(n)%qc,   &
            iv%info(airsr)%pw(n)%error

         do k = 1, iv%info(airsr)%levels(n)
            ! get RH from Q & T    
            if (iv%airsr(n)%q(k)%qc >= obs_qc_pointer .and.  &
                abs(iv%airsr(n)%q(k)%error - missing_r) > 1. .and.  &
                abs(ob%airsr(n)%t(k)       - missing_r) > 1. .and.  &
                abs(iv%airsr(n)%p(k)       - missing_r) > 1.   ) then 
               call da_tp_to_qs(ob%airsr(n)%t(k), iv%airsr(n)%p(k), es, qs)
               rh = (ob%airsr(n)%q(k)/qs) * 100.
               rh_qc = iv%airsr(n)%q(k)%qc
               q = ob%airsr(n)%q(k) * 1e3    !g/kg
               q_qc = iv%airsr(n)%q(k)%qc
               q_err = max(0.01, iv%airsr(n)%q(k)%error * 1e3)  !g/kg
            else
               rh    = missing_r
               rh_qc = missing_data
               q     = missing_r
               q_qc  = missing_data
               q_err = xmiss
            end if
            if (rh_qc < 0) rh_qc = missing_data
            if ( q_qc < 0)  q_qc = missing_data

            t     = ob%airsr(n)%t(k)
            t_qc  = iv%airsr(n)%t(k)%qc
            t_err = iv%airsr(n)%t(k)%error
            if (iv%airsr(n)%t(k)%qc <= fails_error_max) then
               t_qc = missing_data
               t    = missing_r
               t_err = xmiss
            end if

            write(ounit, fmt = trim (fmt_each_inv))     &
               iv%airsr(n)%p(k), zero,    pr_error,     &
               missing_r, missing_data, xmiss,          &
               missing_r, missing_data, xmiss,          &
               iv%airsr(n)%h(k), zero, elv_error,       &
               t,  t_qc, t_err,                         &
               td, td_qc, td_error,                     &
               q, q_qc,  q_err,                         &
               xmiss, xmiss,                            &
               iv%airsr(n)%t(k)%inv, iv%airsr(n)%q(k)%inv*1e3            ! q%inv in g/kg
         end do 
      end do     
   end if
   !--------------------------------------------------------------------------
   ! Write mtgirs
   !--------------------------------------------------------------------------

   if (iv%info(mtgirs)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for mtgirs  ',iv%info(mtgirs)%nlocal
      end if

      do n = 1, iv%info(mtgirs)%nlocal
         if (.not. iv%info(mtgirs)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(mtgirs)%platform(n),    &
            iv%info(mtgirs)%date_char(n),   &
            iv%info(mtgirs)%name(n),        &
            iv%info(mtgirs)%levels(n),      &
            iv%info(mtgirs)%lat(1,n),         &
            iv%info(mtgirs)%lon(1,n),         &
            iv%info(mtgirs)%elv(n),         &
            iv%info(mtgirs)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(mtgirs)%slp(n)%inv, iv%info(mtgirs)%slp(n)%qc, &
            iv%info(mtgirs)%slp(n)%error,                          &
            iv%info(mtgirs)%pw(n)%inv, iv%info(mtgirs)%pw(n)%qc,   &
            iv%info(mtgirs)%pw(n)%error

         do k = 1, iv%info(mtgirs)%levels(n)

            ! get speed and direction
            if (iv%mtgirs(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%mtgirs(n)%v(k)%qc >= obs_qc_pointer) then
               uu = ob%mtgirs(n)%u(k)
               vv = ob%mtgirs(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(mtgirs)%lon(k,n), convert_uv2fd)
               speed_qc  = iv%mtgirs(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%mtgirs(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

!q from kg/kg to g/kg, for error is to 10e-2g/kg
               if (iv%mtgirs(n)%q(k)%qc >= obs_qc_pointer) then
                  rh = ob%mtgirs(n)%q(k) * 1000.
                  rh_qc = iv%mtgirs(n)%q(k)%qc
                  rh_error=iv%mtgirs(n)%q(k)%error * 1000.*100.
                  q = ob%airsr(n)%q(k) * 1e3    !g/kg
                  q_qc = iv%airsr(n)%q(k)%qc
                  q_err = max(0.01, iv%airsr(n)%q(k)%error * 1e3)  !g/kg
               else
                  rh    = missing_r
                  rh_qc = missing_data
                  rh_error=xmiss
                  q     = missing_r
                  q_qc  = missing_data
                  q_err = xmiss
               end if

            write(ounit, fmt = trim (fmt_each_inv))       &
               iv%mtgirs(n)%p(k), zero,    pr_error,       &
               uu, speed_qc, speed_err ,                  &
               vv, speed_qc, speed_err ,                  &
               iv%mtgirs(n)%h(k), zero, elv_error,         &
               t, t_qc, t_err,                            &
               td, td_qc, td_error,                       &
               q, q_qc, q_err,                            &
               iv%mtgirs(n)%u(k)%inv, iv%mtgirs(n)%v(k)%inv, &
               iv%mtgirs(n)%t(k)%inv, iv%mtgirs(n)%q(k)%inv*1e3   !q%inv g/kg
         end do
      end do    
   end if

   !--------------------------------------------------------------------------
   ! Write tamdar
   !--------------------------------------------------------------------------

   if (iv%info(tamdar)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout,fmt='(A,I5)') &
            ' writing filtered obs for tamdar  ',iv%info(tamdar)%nlocal
      end if

      do n = 1, iv%info(tamdar)%nlocal
         if (.not. iv%info(tamdar)%proc_domain(1,n)) cycle

         write(ounit, fmt = fmt_info, iostat = iost) &
            iv%info(tamdar)%platform(n),    &
            iv%info(tamdar)%date_char(n),   &
            iv%info(tamdar)%name(n),        &
            iv%info(tamdar)%levels(n),      &
            iv%info(tamdar)%lat(1,n),         &
            iv%info(tamdar)%lon(1,n),         &
            iv%info(tamdar)%elv(n),         &
            iv%info(tamdar)%id(n)
         write(ounit, fmt = fmt_srfc)  &
            iv%info(tamdar)%slp(n)%inv, iv%info(tamdar)%slp(n)%qc, &
            iv%info(tamdar)%slp(n)%error,                         &
            iv%info(tamdar)%pw(n)%inv, iv%info(tamdar)%pw(n)%qc,   &
            iv%info(tamdar)%pw(n)%error

         do k = 1, iv%info(tamdar)%levels(n)

            ! get speed and direction
            if (iv%tamdar(n)%u(k)%qc >= obs_qc_pointer .and. &
                iv%tamdar(n)%v(k)%qc >= obs_qc_pointer) then
               uu = ob%tamdar(n)%u(k)
               vv = ob%tamdar(n)%v(k)
               call da_ffdduv(speed, dir, uu, vv, iv%info(tamdar)%lon(k,n), convert_uv2fd)
               speed_qc  = iv%tamdar(n)%u(k)%qc
               dir_qc    = speed_qc
               speed_err = iv%tamdar(n)%u(k)%error
            else
               speed     = missing_r
               speed_qc    = missing_data
               dir       = missing_r
               dir_qc    = missing_data
               speed_err = xmiss
               uu        = missing_r
               vv        = missing_r
            end if

            t    =  ob%tamdar(n)%t(k)
            t_qc =  iv%tamdar(n)%t(k)%qc
            t_err= iv%tamdar(n)%t(k)%error
            if (iv%tamdar(n)%t(k)%qc <= fails_error_max) then
               t_qc = missing_data
               t   = missing_r
               t_err = xmiss
            end if
              ! get RH from Q & T
               if (iv%tamdar(n)%q(k)%qc >= obs_qc_pointer .and.  &
                   abs(iv%tamdar(n)%q(k)%error - missing_r) > 1. .and.  &
                   abs(ob%tamdar(n)%t(k)       - missing_r) > 1. .and.  &
                   abs(iv%tamdar(n)%p(k)       - missing_r) > 1.   ) then
                  call da_tp_to_qs(ob%tamdar(n)%t(k), iv%tamdar(n)%p(k), es, qs)
                  rh = (ob%tamdar(n)%q(k)/qs) * 100.
                  rh_qc = iv%tamdar(n)%q(k)%qc
                  q = ob%tamdar(n)%q(k) * 1e3    !g/kg
                  q_qc = iv%tamdar(n)%q(k)%qc
                  q_err = max(0.01, iv%tamdar(n)%q(k)%error * 1e3) !g/kg
               else
                  rh    = missing_r
                  rh_qc = missing_data
                  q     = missing_r
                  q_qc  = missing_data
                  q_err = xmiss
               end if
               if (rh_qc < 0) rh_qc = missing_data
               if ( q_qc < 0)  q_qc = missing_data

               write(ounit, fmt = trim (fmt_each_inv))       &
                  iv%tamdar(n)%p(k), zero,    pr_error,       &
                  uu, speed_qc, speed_err ,                  &
                  vv, speed_qc, speed_err ,                  &
                  iv%tamdar(n)%h(k), zero, elv_error,         &
                  t, t_qc, t_err,                            &
                  td, td_qc, td_error,                       &
                  q, q_qc, q_err,                            &
                  iv%tamdar(n)%u(k)%inv, iv%tamdar(n)%v(k)%inv, &
                  iv%tamdar(n)%t(k)%inv, iv%tamdar(n)%q(k)%inv*1e3   !q%inv g/kg
         end do
      end do
   end if

   !--------------------------------------------------------------------------
   !  Write tamdar_sfc
   !--------------------------------------------------------------------------

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      if (print_detail_f_obs) then
         write(unit=stdout, fmt='(A,I6)') ' writing filtered obs for tamdar_sfc ', &
            iv%info(tamdar_sfc)%nlocal
      end if

      do n = 1, iv%info(tamdar_sfc)%nlocal
         if (.not.iv%info(tamdar_sfc)%proc_domain(1,n)) cycle
         ! Guo .. not write out the data below the lowest model level:

         nlevels = iv%info(tamdar_sfc)%levels(n)
         if (iv%info(tamdar_sfc)%levels(n) > 1) then
            write(unit=stdout, fmt='(3A,I5,A)') &
               ' for SYNOP station ',iv%info(tamdar_sfc)%name(n),' got levels ',&
               iv%info(tamdar_sfc)%levels(n) ,' but wrote only one level'
            nlevels = 1
         end if
         write(ounit, fmt = fmt_info)                &
            iv%info(tamdar_sfc)%platform(n),    &
            iv%info(tamdar_sfc)%date_char(n),   &
            iv%info(tamdar_sfc)%name(n),        &
            nlevels,                      &
            iv%info(tamdar_sfc)%lat(1,n),         &
            iv%info(tamdar_sfc)%lon(1,n),         &
            iv%info(tamdar_sfc)%elv(n),         &
            iv%info(tamdar_sfc)%id(n)
         slp_inv= iv%info(tamdar_sfc)%slp(n)%inv
         slp_qc = iv%info(tamdar_sfc)%slp(n)%qc
         if (iv%tamdar_sfc(n)%p%qc <= fails_error_max) then
            slp_inv = missing_r
            slp_qc  = missing_data
         end if

         write(ounit, fmt = fmt_srfc) &
            slp_inv,slp_qc,iv%info(tamdar_sfc)%slp(n)%error, &
            iv%info(tamdar_sfc)%pw(n)%inv, iv%info(tamdar_sfc)%pw(n)%qc, iv%info(tamdar_sfc)%pw(n)%error

         ! get speed and direction
         if (iv%tamdar_sfc(n)%u%qc >= obs_qc_pointer .and. &
              iv%tamdar_sfc(n)%v%qc >= obs_qc_pointer) then
            uu = ob%tamdar_sfc(n)%u
            vv = ob%tamdar_sfc(n)%v
            call da_ffdduv(speed, dir, uu, vv, iv%info(tamdar_sfc)%lon(1,n), &
               convert_uv2fd)
            speed_qc  = iv%tamdar_sfc(n)%u%qc
            dir_qc    = speed_qc
            speed_err = iv%tamdar_sfc(n)%u%error
         else
            speed     = missing_r
            speed_qc    = missing_data
            dir       = missing_r
            dir_qc    = missing_data
            speed_err = xmiss
            uu        = missing_r
            vv        = missing_r
         end if

         ! get RH from Q & T
         if (iv%tamdar_sfc(n)%q%qc >= obs_qc_pointer .and.  &
              abs(iv%tamdar_sfc(n)%q%error - missing_r) > 1. .and.  &
              abs(ob%tamdar_sfc(n)%t       - missing_r) > 1. .and.  &
              abs(ob%tamdar_sfc(n)%p       - missing_r) > 1.   ) then
            call da_tp_to_qs(ob%tamdar_sfc(n)%t, ob%tamdar_sfc(n)%p, es, qs)
            rh = (ob%tamdar_sfc(n)%q/qs) * 100.
            rh_qc = iv%tamdar_sfc(n)%q%qc
            q = ob%tamdar_sfc(n)%q * 1e3    !g/kg
            q_qc = iv%tamdar_sfc(n)%q%qc
            q_err = max(0.01, iv%tamdar_sfc(n)%q%error * 1e3) !g/kg
         else
            rh    = missing_r
            rh_qc = missing_data
            q     = missing_r
            q_qc  = missing_data
            q_err = xmiss
         end if

         if (rh_qc < 0) rh_qc = missing_data
         if ( q_qc < 0)  q_qc = missing_data

         p    = ob%tamdar_sfc(n)%p
         p_qc = iv%tamdar_sfc(n)%p%qc
         p_err = iv%tamdar_sfc(n)%p%error
         if (iv%tamdar_sfc(n)%p%qc <= fails_error_max) then
            p_qc = missing_data
            p    = missing_r
            p_err  = pr_error
         end if

         t    =  ob%tamdar_sfc(n)%t
         t_qc =  iv%tamdar_sfc(n)%t%qc
         t_err =  iv%tamdar_sfc(n)%t%error
         if (iv%tamdar_sfc(n)%t%qc <= fails_error_max) then
            t_qc = missing_data
            t   = missing_r
            t_err = xmiss
         end if

         write(ounit, fmt = trim (fmt_each_inv))       &
            p,  p_qc,   p_err,                         &
            uu, speed_qc, speed_err ,                  &
            vv, speed_qc, speed_err ,                  &
            iv%info(tamdar_sfc)%elv(n), zero, elv_error,        &
            t, t_qc, t_err,                            &
            td, td_qc, td_error,                       &
            q, q_qc, q_err,                            &
            iv%tamdar_sfc(n)%u%inv, iv%tamdar_sfc(n)%v%inv, &
            iv%tamdar_sfc(n)%t%inv, iv%tamdar_sfc(n)%q%inv*1e3            ! q%inv in g/kg
      end do
   end if



   if (iv%info(ssmi_tb)%nlocal > 0) then
      call da_warning("da_write_modified_filtered_obs.inc",1788, &     
        (/"Currently SSMI brightness temperature info is not written for filtered obs"/))
   end if

   close(ounit)
   call da_free_unit(ounit)
   ! Ensure other processors have written their temporary files
   call mpi_barrier(comm, ierr)
   if (rootproc) then
      call da_get_unit(filtered_obs_unit)
      open (unit=filtered_obs_unit, &
         file= 'filtered_obs' ,form='formatted', status='replace', iostat=iost)
      if (iost /= 0) &
         call da_error("da_write_modified_filtered_obs.inc",1803, (/"Cannot open filtered_obs "/))
      call da_count_filtered_obs(&
         grid%xb%ptop, grid%xb%map, grid%xb%ds,  &
         grid%moad_cen_lat, grid%stand_lon, grid%truelat1, grid%truelat2,                 &
         coarse_ix, coarse_jy, start_x, start_y)
      call da_final_write_modified_filtered_obs
      close (filtered_obs_unit)
      call da_free_unit(filtered_obs_unit)
   end if

   if (trace_use) call da_trace_exit("da_write_modified_filtered_obs")

end subroutine da_write_modified_filtered_obs


subroutine da_write_y (iv, y)

   !-------------------------------------------------------------------------   
   ! Purpose: Writes out components of y=H(x_inc) structure.
   !-------------------------------------------------------------------------   

   implicit none

   type (iv_type), intent(in)    :: iv   ! O-B structure.
   type (y_type), intent(in)     :: y    ! y = H(x_inc) structure.

   integer                       :: ounit ! Output file unit.
   integer                       :: n, k, num_obs, i, ios
   real                          :: f1, f2, f3, f4, f5, f6, f7, dum
   character(len=filename_len)   :: ob_name, filename, file_prefix

   if (trace_use) call da_trace_entry("da_write_y")

   !-------------------------------------------------------------------------   
   ! Fix output unit
   !-------------------------------------------------------------------------   

   if (omb_add_noise) then
      file_prefix='pert_obs.'
   else
      file_prefix='unpert_obs.'
   end if

   dum = -999999.9

    write (unit=filename, fmt='(a,i4.4)') trim(file_prefix), myproc

   call da_get_unit(ounit)
   open (unit=ounit,file=trim(filename),form='formatted', &
         status='replace', iostat=ios )
   if (ios /= 0) then
      call da_error("da_write_y.inc",41, &
         (/"Cannot open (un)perturbed observation file"//filename/))
   end if

   ! [1] Transfer surface obs:

   if (iv%info(synop)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(synop)%nlocal
         if (iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'synop', num_obs
         num_obs = 0
         do n = 1, iv%info(synop)%nlocal
            if (iv%info(synop)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%synop(n)%u%qc, y%synop(n)%u, f1)
               call da_check_missing(iv%synop(n)%v%qc, y%synop(n)%v, f2)
               call da_check_missing(iv%synop(n)%t%qc, y%synop(n)%t, f3)
               call da_check_missing(iv%synop(n)%p%qc, y%synop(n)%p, f4)
               call da_check_missing(iv%synop(n)%q%qc, y%synop(n)%q, f5)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, &
                  dum,dum
            end if
         end do
      end if
   end if

   ! [2] Transfer metar obs:

   if (iv%info(metar)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(metar)%nlocal
         if(iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'metar', num_obs
         num_obs = 0
         do n = 1, iv%info(metar)%nlocal
            if (iv%info(metar)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%metar(n)%u%qc, y%metar(n)%u, f1)
               call da_check_missing(iv%metar(n)%v%qc, y%metar(n)%v, f2)
               call da_check_missing(iv%metar(n)%t%qc, y%metar(n)%t, f3)
               call da_check_missing(iv%metar(n)%p%qc, y%metar(n)%p, f4)
               call da_check_missing(iv%metar(n)%q%qc, y%metar(n)%q, f5)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, &
               dum,dum
            end if
         end do
      end if
   end if

   ! [3] Transfer ships obs:

   if (iv%info(ships)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ships)%nlocal
        if (iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'ships', num_obs
         num_obs = 0
         do n = 1, iv%info(ships)%nlocal
            if (iv%info(ships)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%ships(n)%u%qc, y%ships(n)%u, f1)
               call da_check_missing(iv%ships(n)%v%qc, y%ships(n)%v, f2)
               call da_check_missing(iv%ships(n)%t%qc, y%ships(n)%t, f3)
               call da_check_missing(iv%ships(n)%p%qc, y%ships(n)%p, f4)
               call da_check_missing(iv%ships(n)%q%qc, y%ships(n)%q, f5)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, &
                  dum,dum
            end if
         end do
      end if
   end if

   ! [4.1] Transfer Geo. AMVs Obs:

   if (iv%info(geoamv)%nlocal > 0) then 
      num_obs = 0
      do n = 1, iv%info(geoamv)%nlocal
         if (iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'geoamv', num_obs
         num_obs = 0
         do n = 1, iv%info(geoamv)%nlocal
            if (iv%info(geoamv)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(geoamv)%levels(n)
               do k = 1, iv%info(geoamv)%levels(n)
                  call da_check_missing(iv%geoamv(n)%u(k)%qc, &
                     y%geoamv(n)%u(k), f1)
                  call da_check_missing(iv%geoamv(n)%v(k)%qc, &
                     y%geoamv(n)%v(k), f2)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2 , dum,dum,dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [4.2] Transfer Polar AMVs Obs:

   if (iv%info(polaramv)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(polaramv)%nlocal
         if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'polaramv', num_obs
         num_obs = 0
         do n = 1, iv%info(polaramv)%nlocal
            if (iv%info(polaramv)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)') iv%info(polaramv)%levels(n)
               do k = 1, iv%info(polaramv)%levels(n)
                  call da_check_missing(iv%polaramv(n)%u(k)%qc, &
                     y%polaramv(n)%u(k), f1)
                  call da_check_missing(iv%polaramv(n)%v(k)%qc, &
                     y%polaramv(n)%v(k), f2)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2 , dum,dum,dum,&
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [5] Transfer gpspw obs:

   if (iv%info(gpspw)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(gpspw)%nlocal
         if (iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'gpspw', num_obs
         num_obs = 0
         do n = 1, iv%info(gpspw)%nlocal
            if (iv%info(gpspw)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%gpspw(n)%tpw%qc, y%gpspw(n)%tpw, f1)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, dum,dum,dum,dum, &
                  dum,dum
            end if
         end do
      end if
   end if

   ! [6] Transfer sonde obs:

   if (iv%info(sound)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'sound', num_obs
         num_obs = 0
         do n = 1, iv%info(sound)%nlocal
            if (iv%info(sound)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(sound)%levels(n)
               do k = 1, iv%info(sound)%levels(n)
                  call da_check_missing(iv%sound(n)%u(k)%qc, y%sound(n)%u(k), f1)
                  call da_check_missing(iv%sound(n)%v(k)%qc, y%sound(n)%v(k), f2)
                  call da_check_missing(iv%sound(n)%t(k)%qc, y%sound(n)%t(k), f3)
                  call da_check_missing(iv%sound(n)%q(k)%qc, y%sound(n)%q(k), f4)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2, f3, f4, dum, &
                     dum,dum
               end do
            end if
         end do
      end if

      num_obs = 0
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'sonde_sfc', num_obs
         num_obs = 0
         do n = 1, iv%info(sound)%nlocal
            if (iv%info(sound)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%sonde_sfc(n)%u%qc, y%sonde_sfc(n)%u, f1)
               call da_check_missing(iv%sonde_sfc(n)%v%qc, y%sonde_sfc(n)%v, f2)
               call da_check_missing(iv%sonde_sfc(n)%t%qc, y%sonde_sfc(n)%t, f3)
               call da_check_missing(iv%sonde_sfc(n)%p%qc, y%sonde_sfc(n)%p, f4)
               call da_check_missing(iv%sonde_sfc(n)%q%qc, y%sonde_sfc(n)%q, f5)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, &
                  dum,dum
            end if 
         end do
      end if
   end if

   ! [7] Transfer airep obs:

   if (iv%info(airep)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(airep)%nlocal
         if (iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'airep', num_obs
         num_obs = 0
         do n = 1, iv%info(airep)%nlocal
            if (iv%info(airep)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)') iv%info(airep)%levels(n)
               do k = 1, iv%info(airep)%levels(n)
                  call da_check_missing(iv%airep(n)%u(k)%qc, y%airep(n)%u(k), f1)
                  call da_check_missing(iv%airep(n)%v(k)%qc, y%airep(n)%v(k), f2)
                  call da_check_missing(iv%airep(n)%t(k)%qc, y%airep(n)%t(k), f3)
                  call da_check_missing(iv%airep(n)%q(k)%qc, y%airep(n)%q(k), f4)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2, f3, f4, dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [8] Transfer pilot obs:

   if (iv%info(pilot)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(pilot)%nlocal
         if (iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'pilot', num_obs
         num_obs = 0
         do n = 1, iv%info(pilot)%nlocal
            if (iv%info(pilot)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(pilot)%levels(n)
               do k = 1, iv%info(pilot)%levels(n)
                  call da_check_missing(iv%pilot(n)%u(k)%qc, y%pilot(n)%u(k), f1)
                  call da_check_missing(iv%pilot(n)%v(k)%qc, y%pilot(n)%v(k), f2)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2, dum,dum,dum, &
                     dum,dum
               end do
            end if
         end do
     end if
   end if

   ! [9] Transfer SSM/I obs:SSMI:

   if (iv%info(ssmi_rv)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ssmi_rv)%nlocal
         if (iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'ssmir', num_obs
         num_obs = 0
         do n = 1, iv%info(ssmi_rv)%nlocal
            if (iv%info(ssmi_rv)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%ssmi_rv(n)%speed%qc, &
                                y % ssmi_rv(n) % speed, f1)
               call da_check_missing(iv%ssmi_rv(n)% tpw % qc, &
                                y % ssmi_rv(n) % tpw, f2)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, dum,dum,dum, &
                  dum,dum
            end if 
         end do
      end if
   end if

   if (iv%info(ssmi_tb)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ssmi_tb)%nlocal            
         if (iv%info(ssmi_tb)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'ssmiT', num_obs
         num_obs = 0
         do n = 1, iv%info(ssmi_tb)%nlocal
            if (iv%info(ssmi_tb)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%ssmi_tb(n)%tb19h%qc, &
                                      y %ssmi_tb(n)%tb19h, f1)
               call da_check_missing(iv%ssmi_tb(n)%tb19v%qc, &
                                      y %ssmi_tb(n)%tb19v, f2)
               call da_check_missing(iv%ssmi_tb(n)%tb22v%qc, &
                                      y %ssmi_tb(n)%tb22v, f3)
               call da_check_missing(iv%ssmi_tb(n)%tb37h%qc, &
                                      y %ssmi_tb(n)%tb37h, f4)
               call da_check_missing(iv%ssmi_tb(n)%tb37v%qc, &
                                      y %ssmi_tb(n)%tb37v, f5)
               call da_check_missing(iv%ssmi_tb(n)%tb85h%qc, &
                                      y %ssmi_tb(n)%tb85h, f6)
               call da_check_missing(iv%ssmi_tb(n)%tb85v%qc, &
                                      y %ssmi_tb(n)%tb85v, f7)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, f6, f7
            end if
         end do
      end if
   end if

   ! [10] Transfer satem obs:

   if (iv%info(satem)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(satem)%nlocal            
         if (iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'satem', num_obs
         num_obs = 0
         do n = 1, iv%info(satem)%nlocal
            if (iv%info(satem)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(satem)%levels(n)
               do k = 1, iv%info(satem)%levels(n)
                  call da_check_missing(iv%satem(n)%thickness(k)%qc, &
                     y % satem(n) % thickness(k), f1)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, dum,dum,dum,dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if
   
   ! [11] Transfer ssmt1 obs:

   if (iv%info(ssmt1)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ssmt1)%nlocal            
         if (iv%info(ssmt1)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'ssmt1', num_obs
         num_obs = 0
         do n = 1, iv%info(ssmt1)%nlocal
            if (iv%info(ssmt1)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(ssmt1)%levels(n)
               do k = 1, iv%info(ssmt1)%levels(n)
                  call da_check_missing(iv%ssmt1(n)%t(k)%qc, &
                     y % ssmt1(n) % t(k), f1)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, dum,dum,dum,dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [12] Transfer ssmt2 obs:

   if (iv%info(ssmt2)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(ssmt2)%nlocal            
         if (iv%info(ssmt2)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'ssmt2', num_obs
         num_obs = 0
         do n = 1, iv%info(ssmt2)%nlocal
            if (iv%info(ssmt2)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(ssmt2)%levels(n)
               do k = 1, iv%info(ssmt2)%levels(n)
                  call da_check_missing(iv%ssmt2(n)%rh(k)%qc, &
                  y % ssmt2(n) % rh(k), f1)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, dum,dum,dum,dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [13] Transfer scatterometer obs:

   if (iv%info(qscat)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(qscat)%nlocal            
         if (iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'qscat', num_obs
         num_obs = 0
         do n = 1, iv%info(qscat)%nlocal
            if (iv%info(qscat)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%qscat(n)%u%qc, y%qscat(n)%u, f1)
               call da_check_missing(iv%qscat(n)%v%qc, y%qscat(n)%v, f2)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, dum,dum,dum, &
                  dum,dum
            end if
         end do
      end if
   end if
   
   ! [14] Transfer profiler obs:

   if (iv%info(profiler)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(profiler)%nlocal            
         if (iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'profiler', num_obs
         num_obs = 0
         do n = 1, iv%info(profiler)%nlocal
            if (iv%info(profiler)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(profiler)%levels(n)
               do k = 1, iv%info(profiler)%levels(n)
                  call da_check_missing(iv%profiler(n)%u(k)%qc, &
                     y%profiler(n)%u(k), f1)
                  call da_check_missing(iv%profiler(n)%v(k)%qc, &
                     y%profiler(n)%v(k), f2)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2, dum,dum,dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [15] Transfer buoy  obs:

   if (iv%info(buoy)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(buoy)%nlocal            
         if (iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'buoy', num_obs
         num_obs = 0
         do n = 1, iv%info(buoy)%nlocal
            if (iv%info(buoy)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               call da_check_missing(iv%buoy(n)%u%qc, y%buoy(n)%u, f1)
               call da_check_missing(iv%buoy(n)%v%qc, y%buoy(n)%v, f2)
               call da_check_missing(iv%buoy(n)%t%qc, y%buoy(n)%t, f3)
               call da_check_missing(iv%buoy(n)%p%qc, y%buoy(n)%p, f4)
               call da_check_missing(iv%buoy(n)%q%qc, y%buoy(n)%q, f5)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, &
                  dum,dum
            end if
         end do
      end if
   end if

   ! [16] Transfer TC bogus obs:

   if (iv%info(bogus)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(bogus)%nlocal            
         if (iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'bogus', num_obs
         num_obs = 0
         do n = 1, iv%info(bogus)%nlocal
            if (iv%info(bogus)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1                         
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, dum,dum,dum,dum,dum,dum
               write(ounit,'(i8)')iv%info(bogus)%levels(n)
               do k = 1, iv%info(bogus)%levels(n)
                  call da_check_missing(iv%bogus(n)%u(k)%qc, y%bogus(n)%u(k), f2)
                  call da_check_missing(iv%bogus(n)%v(k)%qc, y%bogus(n)%v(k), f3)
                  call da_check_missing(iv%bogus(n)%t(k)%qc, y%bogus(n)%t(k), f4)
                  call da_check_missing(iv%bogus(n)%q(k)%qc, y%bogus(n)%q(k), f5)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f2, f3, f4, f5, dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [17] Transfer AIRS retrievals:

   if (iv%info(airsr)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(airsr)%nlocal
         if (iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'airsr', num_obs
         num_obs = 0
         do n = 1, iv%info(airsr)%nlocal
            if (iv%info(airsr)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(airsr)%levels(n)
               do k = 1, iv%info(airsr)%levels(n)
                  call da_check_missing(iv%airsr(n)%t(k)%qc, y%airsr(n)%t(k), f1)
                  call da_check_missing(iv%airsr(n)%q(k)%qc, y%airsr(n)%q(k), f2)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2, dum, dum, &
                     dum, dum,dum
               end do
            end if
         end do
      end if
   end if
  
   ! [18] Transfer Radiance obs:

   if (iv%num_inst > 0) then
      do i = 1, iv%num_inst                 ! loop for sensor
         if (iv%instid(i)%num_rad < 1) cycle
         do k = 1,iv%instid(i)%nchan        ! loop for channel
            ! Counting number of obs for channel k
            num_obs = 0
            do n = 1,iv%instid(i)%num_rad      ! loop for pixel
               if (iv%instid(i)%info%proc_domain(1,n) .and. &
                  (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer)) then
                  num_obs = num_obs + 1
               end if
            end do                                ! end loop for pixel
            if (num_obs < 1) cycle

            write(ob_name,'(a,a,i4.4)') trim(iv%instid(i)%rttovid_string),'-', &
	    iv%instid(i)%ichan(k)
            write(ounit,'(a20,i8)')  ob_name,num_obs

            num_obs = 0
            do n= 1, iv%instid(i)%num_rad      ! loop for pixel
              if(iv%instid(i)%info%proc_domain(1,n) .and. &
                 (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer)) then
                    num_obs = num_obs + 1
                    write(ounit,'(2i8,e15.7)')num_obs, 1, y%instid(i)%tb(k,n)
              end if
            end do                                ! end loop for pixel
         end do                                ! end loop for channel
      end do                                   ! end loop for sensor
   end if

   ! [19] Transfer gpsref obs:

   if (iv%info(gpsref)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(gpsref)%nlocal
         if (iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'gpsref', num_obs
         num_obs = 0
         do n = 1, iv%info(gpsref)%nlocal
            if (iv%info(gpsref)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(gpsref)%levels(n)
               do k = 1, iv%info(gpsref)%levels(n)
                  call da_check_missing(iv%gpsref(n)%ref(k)%qc, y%gpsref(n)%ref(k), f1)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, dum, dum, dum, dum, dum,dum
               end do
            end if
         end do
      end if
   end if

   ! [20] Transfer tamdar obs:

   if (iv%info(tamdar)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(tamdar)%nlocal
         if (iv%info(tamdar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'tamdar', num_obs
         num_obs = 0
         do n = 1, iv%info(tamdar)%nlocal
            if (iv%info(tamdar)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')iv%info(tamdar)%levels(n)
               do k = 1, iv%info(tamdar)%levels(n)
                  call da_check_missing(iv%tamdar(n)%u(k)%qc, y%tamdar(n)%u(k), f1)
                  call da_check_missing(iv%tamdar(n)%v(k)%qc, y%tamdar(n)%v(k), f2)
                  call da_check_missing(iv%tamdar(n)%t(k)%qc, y%tamdar(n)%t(k), f3)
                  call da_check_missing(iv%tamdar(n)%q(k)%qc, y%tamdar(n)%q(k), f4)
                  write(ounit,'(2i8,7e15.7)')num_obs, k, f1, f2, f3, f4, dum, &
                     dum,dum
               end do
            end if
         end do
      end if
   end if

! Now tamdar_sfc
   if (iv%info(tamdar_sfc)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(tamdar_sfc)%nlocal
         if (iv%info(tamdar_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'sonde_sfc', num_obs
         num_obs = 0
         do n = 1, iv%info(tamdar_sfc)%nlocal
            if (iv%info(tamdar_sfc)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1
               call da_check_missing(iv%tamdar_sfc(n)%u%qc, y%tamdar_sfc(n)%u, f1)
               call da_check_missing(iv%tamdar_sfc(n)%v%qc, y%tamdar_sfc(n)%v, f2)
               call da_check_missing(iv%tamdar_sfc(n)%t%qc, y%tamdar_sfc(n)%t, f3)
               call da_check_missing(iv%tamdar_sfc(n)%p%qc, y%tamdar_sfc(n)%p, f4)
               call da_check_missing(iv%tamdar_sfc(n)%q%qc, y%tamdar_sfc(n)%q, f5)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, f2, f3, f4, f5, &
                  dum,dum
            end if
         end do
      end if
   end if

   ! [21] Transfer rainfall obs:

   if (iv%info(rain)%nlocal > 0) then
      num_obs = 0
      do n = 1, iv%info(rain)%nlocal
         if (iv%info(rain)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      if (num_obs > 0) then
         write(ounit,'(a20,i8)')'rain', num_obs
         num_obs = 0
         do n = 1, iv%info(rain)%nlocal
            if (iv%info(rain)%proc_domain(1,n)) then
               num_obs = num_obs + 1
               write(ounit,'(i8)')  1
               call da_check_missing(iv%rain(n)%rain%qc, y%rain(n)%rain, f1)
               write(ounit,'(2i8,7e15.7)')num_obs, 1, f1, dum,dum,dum,dum, &
                  dum,dum
            end if
         end do
      end if
   end if

   close (ounit)
   call da_free_unit(ounit)

   if (trace_use) call da_trace_exit("da_write_y")

end subroutine da_write_y


subroutine da_read_obs_bufr (iv)

   !---------------------------------------------------------------------------
   ! Purpose: Read 1 observation file for input to wrfvar
   !---------------------------------------------------------------------------
   !   METHOD: use F90 sequential data structure to avoid read file twice
   !            so that da_scan_obs_bufr is not necessary any more.
   !            1. read prepbufr data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2009/10/13  - F90 sequential structure  Peng Xiu
   ! 
   !----------------------------------------------------------------------------
   
   use da_define_structures, only: da_allocate_observations
 
   implicit none

   type (iv_type),             intent(inout) :: iv
   
   
   character(len=10)    :: filename
   real, parameter    :: r8bfms = 9.0D08  ! 1 missing value threshold
   logical                      :: match, end_of_file
   character(len=8)             :: subst2, csid, csid2
   integer              :: idate2, nlevels2, lv1, lv2
   real               :: hdr_save(7)
   real*8             :: hdr2(7), r8sid, r8sid2
   real               :: pob1, pob2
   real               :: temp(8,255)
   real               :: obs_save(8,255)
   real*8             :: obs2(8,255), qms2(8,255)
   real*8             :: oes2(8,255), pco2(8,255)
   real               :: pmo_save(2,1)
   real*8             :: pmo2(2,1)
   equivalence                     (r8sid, csid), (r8sid2, csid2)
   ! for thinning
   real               :: tdiff                       ! DHR
   real               :: dlat_earth,dlon_earth,crit
   integer              :: itt,itx,iout
   logical                      :: iuse
   
   type (multi_level_type)      :: platform
   logical                      :: outside, outside_all, outside_time
   integer              :: ilocal(num_ob_indexes)
   integer              :: ntotal(num_ob_indexes)
   integer              :: nlocal(num_ob_indexes)
   integer              :: tp(num_ob_indexes)
   character(len=40)            :: obstr,hdstr,qmstr,oestr, pcstr
   character(len=8)             :: subset
   character(len=14)            :: cdate, dmn, obs_date,platform_name
   real*8             :: hdr(7)
   real*8             :: pmo(2,1)
   real*8             :: obs(8,255),qms(8,255),oes(8,255),pco(8,255)
   real               :: woe,toe,qoe,poe,pob,pwe
   real*8             :: obs_time
   integer              :: iyear, imonth, iday, ihour, imin

   integer              :: iost, ndup, n, i, j, k, kk,surface_level, num_report, i1, i2
   integer              :: iret, idate, kx, old_nlevels,nlevels, t29,ifgat,ii
   integer              :: cat,zqm,pqm,qqm,tqm,wqm,pwq,pmq
   integer              :: tpc, wpc,iret2
   integer              :: iunit, fm , obs_index

   integer              :: qflag           ! Flag for retaining data

   real, allocatable  :: in(:), out(:)
   integer :: num_bufr(7)
   logical                      :: found

   logical                      :: use_errtable
   integer              :: junit, itype, ivar
   real                 :: oetab(300,33,6)  ! 300 ob types, 33 levels (rows), 6 variables (columns)
   real                 :: err_uv, err_t, err_p, err_q, err_pw, coef
   integer              :: ibufr
   integer              :: num_outside_all, num_outside_time, num_thinned,num_p,numbufr

!added by yw
   logical :: fexist, ywflag

   type datalink_BUFR              !for PREPBUFR data reading
       type (multi_level_type_BUFR)          :: platform_BUFR
       integer                               :: fm_BUFR
       integer                               :: t29_BUFR
       integer                               :: ifgat_BUFR
       integer                               :: nlevels_BUFR
       integer                               :: kx_BUFR
       real                                  :: pco_BUFR(8,255)
       type(datalink_BUFR), pointer :: next
   end type datalink_BUFR

   type(datalink_BUFR),pointer  :: head=>null(), plink=>null()

   if (trace_use) call da_trace_entry("da_read_obs_bufr")

!  0.0 Initialize variables
!--------------------------------------------------------------
   ilocal(:) = 0
   ntotal(:) = 0
   nlocal(:) = 0
  
   err_uv = 10.0 ! m/s
   err_t  = 5.0  ! degree
   err_p  = 200  ! Pa
   err_q  = 10   ! RH percent
   err_pw = 0.2  ! cm

   ! quality marker 0: Keep (always assimilate)
   !                1: Good
   !                2: Neutral or not checked
   !                3: Suspect
   if ( anal_type_verify ) then
      qflag = min(qmarker_retain,2)
   else
      qflag = qmarker_retain
   end if
   write(unit=message(1),fmt='(a,i1,a)') &
         "PREPBUFR ob with quality marker <= ", qflag, " will be retained."
   call da_message(message(1:1))
   
   num_report       = 0
   num_outside_all  = 0
   num_outside_time = 0
   num_thinned      = 0
   num_p            = 0
   tp(:)            = 0
   

! 1.0  Open file
!---------------------------------------------------------------- 
!
!check if input file exists
num_bufr(:)=0
numbufr=0
if (num_fgat_time>1) then
  do i=1,7
     call da_get_unit(iunit)
     write(filename,fmt='(A,I1,A)') 'ob0',i,'.bufr'
     open(unit   = iunit, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
     if (iost == 0) then
        numbufr=numbufr+1
	num_bufr(numbufr)=i
     else
        close (iunit)
     end if
     call da_free_unit(iunit)
   end do
 else
   numbufr=1
 end if
  
  if (numbufr==0) numbufr=1

!yw added by Wei Yu
  inquire (file="ob1.bufr", exist=fexist)
  ywflag = .false.  
  if (fexist ) then
     numbufr = 2
     ywflag = .true. 
  endif 
!yw end added

     
bufrfile:  do ibufr=1,numbufr   
     if (num_fgat_time==1) then
         filename='ob.bufr'
     else
         if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
	    filename='ob.bufr'
        else
            write(filename,fmt='(A,I1,A)') 'ob0',num_bufr(ibufr),'.bufr'   
        end if
     end if

!yw added by Wei Yu
     if( (ibufr .eq. 2) .and. ywflag)  then
         filename='ob1.bufr'
     endif
!yw end added
       
!
!     We want to use specific unit number to read prepbufr data, which enables us to control its endianness
      iunit = 96 

      open(unit   = iunit, FILE   = trim(filename), &
         iostat =  iost, form = 'unformatted', STATUS = 'OLD')
      if (iost /= 0) then
         write(unit=message(1),fmt='(A,I5,A)') &
            "Error",iost," opening PREPBUFR obs file "//trim(filename)
         call da_warning("da_read_obs_bufr.inc",193,message(1:1))
         if ( num_fgat_time == 1 ) then
            call da_free_unit(iunit)
            if (trace_use) call da_trace_exit("da_read_obs_bufr")
            return
         else
            cycle bufrfile
         end if
      end if
   ! open observation error table if provided.
   call da_get_unit(junit)
   open (unit=junit, file='obs_errtable', form='formatted', status='old', &
         iostat=iost)
   if ( iost /= 0 ) then
      use_errtable = .false.
      call da_free_unit(junit)
   else
      use_errtable = .true.
      write(unit=message(1),fmt='(A)') &
            "obs_errtable file is found. Will use user-provided obs errors."
      call da_message(message(1:1))
   end if
   if ( use_errtable ) then
      read_loop: do
         read (junit,'(1x,i3)',iostat=iost) itype
         if ( iost /=0 ) exit read_loop
         do k = 1, 33
            read (junit,'(1x,6e12.5)',iostat=iost) (oetab(itype,k,ivar),ivar=1,6)
            if ( iost /=0 ) exit read_loop
         end do
      end do read_loop
   end if

   hdstr='SID XOB YOB DHR TYP ELV T29'
   obstr='POB QOB TOB ZOB UOB VOB PWO CAT' ! observation
   qmstr='PQM QQM TQM ZQM WQM NUL PWQ NUL' ! quality marker
   oestr='POE QOE TOE NUL WOE NUL PWE NUL' ! observation error
   pcstr='PPC QPC TPC ZPC WPC NUL PWP NUL' ! program code

   call openbf(iunit,'IN',iunit)
   call datelen(10)

   call readns(iunit,subset,idate,iret)  ! read in the next subset
   if ( iret /= 0 ) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iret," reading PREPBUFR obs file "//trim(filename)
      call da_warning("da_read_obs_bufr.inc",239,message(1:1))
      call closbf(iunit)
      if (trace_use) call da_trace_exit("da_read_obs_bufr")
      return
   end if
   !rewind(iunit)

   write(unit=message(1),fmt='(a,i10)') 'BUFR file date is: ', idate
   call da_message(message(1:1))


!  2.0 read data
!  scan reports first
!--------------------------------------------------------------
  
   match        = .false.
   end_of_file  = .false.
   outside_all  = .false.
   outside_time = .false.
   reports: do while ( .not. end_of_file )

      if ( match .or. outside_all .or. outside_time ) then
         call readns(iunit,subset,idate,iret)  ! read in the next subset
         if ( iret /= 0 ) then 
            write(unit=message(1),fmt='(A,I3,A,I3)') & 
               "return code from readns",iret,       &
               "reach the end of PREPBUFR obs unit",iunit
            !call da_warning("da_read_obs_bufr.inc",266,message(1:1))
            exit reports
         end if
      end if

      num_report = num_report+1

      call ufbint(iunit,hdr,7,1,iret2,hdstr)
      call ufbint(iunit,pmo,2,1,nlevels,'PMO PMQ')
      call ufbint(iunit,qms,8,255,nlevels,qmstr)
      call ufbint(iunit,oes,8,255,nlevels,oestr)
      call ufbint(iunit,pco,8,255,nlevels,pcstr)
      call ufbint(iunit,obs,8,255,nlevels,obstr)
      
      r8sid = hdr(1)
      platform % info % name(1:8)  = subset
      platform % info % name(9:40) = '                                '
      platform % info % id(1:5)    = csid(1:5)
      platform % info % id(6:40)   = '                                   '
      platform % info % dhr        = hdr(4)    ! difference in hour
      platform % info % elv        = hdr(6)
      platform % info % lon        = hdr(2)
      platform % info % lat        = hdr(3)

      ! blacklisted stations should be handled through an external table.
      ! For now, temporary fix is implemented here for known incorrect 
      ! station info in NCEP PREPBUFR file
      if ( trim(platform%info%id) == 'BGQQ' ) then
         platform%info%elv = 19
         platform%info%lon = -69.21
         platform%info%lat = 77.46
      end if
      if ( trim(platform%info%id) == 'UWKE' ) then
         platform%info%elv = 194
         platform%info%lon = 52.09
         platform%info%lat = 55.56
      end if

      ! Put a check on Lon and Lat
      if ( platform%info%lon >= 180.0 ) platform%info%lon = platform%info%lon - 360.0
      ! Fix funny wind direction at Poles
      !if (platform%info%lat < -89.9999 .or. platform%info%lat > 89.9999) then
      !   platform%info%lon = 0.0
      !end if
      platform%info%lat = max(platform%info%lat, -89.95)
      platform%info%lat = min(platform%info%lat,  89.95)

      ! Restrict to a range of reports, useful for debugging

      if (num_report < report_start) cycle reports
      if (num_report > report_end)  exit reports

      call da_llxy (platform%info, platform%loc,outside, outside_all)

      if (outside_all) then
         num_outside_all = num_outside_all + 1
         if ( print_detail_obs ) then
            write(unit=stderr,fmt='(a,1x,a,2(1x,f8.3),a)')  &
               platform%info%name(1:8),platform%info%id(1:5), &
               platform%info%lat, platform%info%lon, '  -> outside_domain'
         end if
         cycle reports
      end if

      ! check date
      write(cdate,'(i10)') idate
      write(dmn,'(i4,a1)') int(platform%info%dhr*60.0), 'm'
      call da_advance_time (cdate(1:10), trim(dmn), obs_date)
      if ( obs_date(13:14) /= '00' ) then
         write(0,*) 'wrong date: ', trim(cdate), trim(dmn), trim(obs_date)
         call da_error("da_read_obs_bufr.inc",336,(/"Wrong date"/))
      else
         read (obs_date(1:12),'(i4,4i2)') iyear, imonth, iday, ihour, imin
      end if
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
      if (obs_time < time_slots(0) .or.  &
         obs_time >= time_slots(num_fgat_time)) then
         outside_time = .true.
         num_outside_time = num_outside_time + 1
         if ( print_detail_obs ) then
            write(unit=stderr,fmt='(a,1x,a,1x,a,a)')  &
               platform%info%name(1:8),platform%info%id(1:5), &
               trim(obs_date), '  -> outside_time'
         end if
         cycle reports       
      else
         outside_time = .false.
      end if
      
!--------  determine FGAT index ifgat
 
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do       
         
      write(unit=platform%info%date_char, fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
         iyear, '-', imonth, '-', iday, '_', ihour, ':', imin, ':', 0

      if (print_detail_obs) then
         ! Simplistic approach, could be improved to get it all done on PE 0
         if (.NOT. outside) then
            write(unit=stdout,fmt='(A,1X,I8,1X,A,F8.2,A,F8.2,A,1X,A,I3,1X,A,2F8.3)') &
               "Report",num_report,"at",platform%info%lon,"E",platform%info%lat,"N", &
               "on processor", myproc,"position", platform%loc%x,platform%loc%y
         end if
      end if

      t29 = int(0.1 + hdr(7))
      kx=int(0.1+hdr(5))

      if ( use_errtable ) then
         do k = 1, nlevels
            pob = obs(1,k)
            do lv1 = 1, 32
               if ( pob >= oetab(kx,lv1+1,1) .and. pob <= oetab(kx,lv1,1) ) then
                  coef = (pob-oetab(kx,lv1,1))/(oetab(kx,lv1,1)-oetab(kx,lv1+1,1))
                  oes(1,k) = (1.0+coef)*oetab(kx,lv1,5)-coef*oetab(kx,lv1+1,5) !p
                  oes(2,k) = (1.0+coef)*oetab(kx,lv1,3)-coef*oetab(kx,lv1+1,3) !q
                  oes(3,k) = (1.0+coef)*oetab(kx,lv1,2)-coef*oetab(kx,lv1+1,2) !t
                  oes(5,k) = (1.0+coef)*oetab(kx,lv1,4)-coef*oetab(kx,lv1+1,4) !uv
                  oes(7,k) = (1.0+coef)*oetab(kx,lv1,6)-coef*oetab(kx,lv1+1,6) !pw
                  exit
               end if
            end do
         end do
      end if

      call readns(iunit,subst2,idate2,iret)

      if ( iret /= 0 ) then
         end_of_file = .true.
      else
         match_check: do
            call ufbint(iunit,hdr2,7,1,iret2,hdstr)
            ! check if this subset and the previous one are matching mass and wind
            match = .true.
            if ( subset /= subst2 ) then
               match = .false.
               exit match_check
            end if
            r8sid2 = hdr2(1)
            if ( csid /= csid2 ) then  ! check SID
               match = .false.
               exit match_check
            end if
            do i = 2, 4   ! check XOB, YOB, DHR
               if ( hdr(i) /= hdr2(i) ) then
                  match = .false.
                  exit match_check
               end if
            end do
            if ( hdr(6) /= hdr2(6) ) then   ! check ELV
               match = .false.
               exit match_check
            end if
            !The two headers match, now read data from the second subset
            call ufbint(iunit,pmo2,2,1,nlevels2,'PMO PMQ')
            call ufbint(iunit,qms2,8,255,nlevels2,qmstr)
            call ufbint(iunit,oes2,8,255,nlevels2,oestr)
            call ufbint(iunit,pco2,8,255,nlevels2,pcstr)
            call ufbint(iunit,obs2,8,255,nlevels2,obstr)

            if ( use_errtable ) then
               kx = nint(hdr2(5))
               do k = 1, nlevels2
                  pob = obs2(1,k)
                  do lv1 = 1, 32
                     if ( pob >= oetab(kx,lv1+1,1) .and. pob <= oetab(kx,lv1,1) ) then
                        coef = (pob-oetab(kx,lv1,1))/(oetab(kx,lv1,1)-oetab(kx,lv1+1,1))
                        oes2(1,k) = (1.0+coef)*oetab(kx,lv1,5)-coef*oetab(kx,lv1+1,5) !p
                        oes2(2,k) = (1.0+coef)*oetab(kx,lv1,3)-coef*oetab(kx,lv1+1,3) !q
                        oes2(3,k) = (1.0+coef)*oetab(kx,lv1,2)-coef*oetab(kx,lv1+1,2) !t
                        oes2(5,k) = (1.0+coef)*oetab(kx,lv1,4)-coef*oetab(kx,lv1+1,4) !uv
                        oes2(7,k) = (1.0+coef)*oetab(kx,lv1,6)-coef*oetab(kx,lv1+1,6) !pw
                        exit
                     end if
                  end do
               end do
            end if

            ! If this is a surface report, the wind subset precedes the
            ! mass subset - switch the subsets around in order to combine
            ! the surface pressure properly
            kx = nint(hdr(5))
            if ( kx == 280 .or. kx == 281 .or. kx == 284  .or.  &
                 kx == 287 .or. kx == 288  ) then
               pmo_save = pmo2
               pmo2 = pmo
               pmo = pmo_save
               temp = obs2
               obs2 = obs
               obs = temp
               hdr_save = hdr2
               hdr2 = hdr
               hdr = hdr_save
               temp = qms2
               qms2 = qms
               qms = temp
               temp = oes2
               oes2 = oes
               oes = temp
               temp = pco2
               pco2 = pco
               pco = temp
            end if

            ! combine the two matching subsets
            do i = 1, 2
               if ( pmo(i,1) > r8bfms ) then
                  pmo(i,1) = pmo2(i,1)
               end if
            end do
            lev_loop: do lv2 = 1, nlevels2
               do lv1 = 1, nlevels
                  pob1 = obs(1,lv1)
                  pob2 = obs2(1,lv2)
                  if ( pob1 == pob2 ) then
                     do i = 1, 7   ! skip the CAT
                        if ( obs(i,lv1) > r8bfms ) then
                           obs(i,lv1) = obs2(i,lv2)
                           if ( obs2(i,lv2) <= r8bfms ) then
                              obs(8,lv1) = obs2(8,lv2)  ! rewrite CAT
                           end if
                        end if
                        if ( oes(i,lv1) > r8bfms ) then
                           oes(i,lv1) = oes2(i,lv2)
                        end if
                        if ( pco(i,lv1) > r8bfms ) then
                           pco(i,lv1) = pco2(i,lv2)
                        end if
                     end do
                     do i = 1, 8
                        if ( qms(i,lv1) > r8bfms ) then
                           qms(i,lv1) = qms2(i,lv2)
                        end if
                     end do
                     cycle lev_loop
                  else if ( (pob2 > pob1) .or. (lv1 .eq. nlevels) ) then
                     nlevels = nlevels + 1
                     obs(:,nlevels) = obs2(:,lv2)
                     qms(:,nlevels) = qms2(:,lv2)
                     oes(:,nlevels) = oes2(:,lv2)
                     pco(:,nlevels) = pco2(:,lv2)
                     cycle lev_loop
                  end if
               end do
            end do lev_loop
            ! sort combined report in descending pressure order
            do i1 = 1, nlevels-1
               do i2 = i1+1, nlevels
                  if ( obs(1,i2) .gt. obs(1,i1) ) then
                     temp(:,1) = obs(:,i1)
                     obs(:,i1) = obs(:,i2)
                     obs(:,i2) = temp(:,1)
                     temp(:,1) = qms(:,i1)
                     qms(:,i1) = qms(:,i2)
                     qms(:,i2) = temp(:,1)
                     temp(:,1) = oes(:,i1)
                     oes(:,i1) = oes(:,i2)
                     oes(:,i2) = temp(:,1)
                     temp(:,1) = pco(:,i1)
                     pco(:,i1) = pco(:,i2)
                     pco(:,i2) = temp(:,1)
                  end if
               end do
            end do
            exit match_check
         end do match_check

         if ( .not. match ) then
            subset = subst2
            idate = idate2
         end if

      end if

      ! skip some types
      !  61: Satellite soundings/retrievals/radiances
      !  66: SSM/I rain rate product
      !  72: NEXTRAD VAD winds
      if ( t29 == 61 .or. t29 == 66 .or. t29 == 72 ) cycle reports 
      
      platform % info % levels    = nlevels

      platform % loc % slp %inv   = missing_r     
      platform % loc % slp %qc    = missing_data
      platform % loc % slp %error = err_p
      platform % loc % pw %inv    = missing_r     
      platform % loc % pw %qc     = missing_data
      platform % loc % pw %error  = err_pw
      pmq=nint(pmo(2,1))
      if ( pmq <= qflag .and. pmq >= 0 .and. &
           pmo(1,1) < r8bfms ) then
         platform % loc % slp % inv =pmo(1,1)*100.0
         platform % loc % slp % qc  =pmq
         platform % loc % slp % error = err_p          ! hardwired
      end if
      pwq=nint(qms(7,1))
      pwe = min(err_pw, oes(7,1))
      if ( pwq <= qflag .and. pwq >= 0 .and. &
           obs(7,1) < r8bfms ) then
         platform % loc % pw % inv = obs(7,1) * 0.1    ! convert to cm
         platform % loc % pw % qc  =  pwq
         platform % loc % pw % error = pwe          ! hardwired
      end if

      !$OMP PARALLEL DO &
      !$OMP PRIVATE (i)
      do i=1,max_ob_levels
         platform % each (i) % height  = missing_r
         platform % each (i) % height_qc = missing_data

         platform % each (i) % zk = missing_r
            
         platform % each (i) % u % inv = missing_r
         platform % each (i) % u % qc  = missing_data
         platform % each (i) % u % error = err_uv

         platform % each (i) % v = platform % each (i) % u

         platform % each (i) % t % inv = missing_r
         platform % each (i) % t % qc  = missing_data
         platform % each (i) % t % error = err_t

         platform % each (i) % p % inv = missing_r      
         platform % each (i) % p % qc  = missing_data
         platform % each (i) % p % error = err_p

         platform % each (i) % q % inv = missing_r
         platform % each (i) % q % qc  = missing_data
         platform % each (i) % q % error = err_q
      end do 
      !$OMP END PARALLEL DO

      !!$OMP PARALLEL DO &
      !!$OMP PRIVATE (k, tpc, wpc, pqm, qqm, tqm, wqm, zqm, cat, toe, poe, qoe, woe) 
      do k = 1, platform % info % levels

         tpc = nint(pco(3,k))
         wpc = nint(pco(5,k))

         ! set t units to Kelvin
         if (obs(3,k) > -200.0 .and. obs(3,k) < 300.0) then
            obs(3,k) = obs(3,k) + t_kelvin
         end if

         ! scale q and compute t from tv, if they aren't missing
         if (obs(2,k) > 0.0 .and. obs(2,k) < r8bfms) then   
            obs(2,k) = obs(2,k)*1e-6
            if (obs(3,k) > -200.0 .and. obs(3,k) < 350.0) then
               if ( tpc >= 8 ) then   ! program code 008 VIRTMP
                  ! 0.61 is used in NCEP prepdata.f to convert T to Tv
                  obs(3,k) = obs(3,k) / (1.0 + 0.61 * obs(2,k))
               end if
            end if
         end if

         pqm=nint(qms(1,k))
         qqm=nint(qms(2,k))
         tqm=nint(qms(3,k))
         zqm=nint(qms(4,k))
         wqm=nint(qms(5,k))
         cat=nint(obs(8,k))

         toe = min(err_t, oes(3,k))
         woe = min(err_uv, oes(5,k))
         qoe = min(err_q, oes(2,k)*10.0)  ! convert to % from PREPBUFR percent divided by 10
         poe = min(err_p, oes(1,k)*100.0) ! convert to Pa

         if ( tqm <= qflag .and. tqm >= 0 .and. &
              obs(3,k) < r8bfms ) then
            platform % each (k) % t % inv =obs(3,k)
            platform % each (k) % t % qc  =tqm
            platform % each (k) % t % error =toe
         end if

         if ( wqm <= qflag .and. wqm >= 0 .and. &
             obs(5,k) < r8bfms .and. obs(6,k) < r8bfms ) then
            platform % each (k) % u % inv =obs(5,k)
            platform % each (k) % v % inv =obs(6,k)
            platform % each (k) % u % qc  =wqm
            platform % each (k) % u % error =woe
            platform % each (k) % v % qc  =wqm
            platform % each (k) % v % error =woe

            ! Convert earth wind to model wind.
            ! note on SPSSMI wind: only wspd available (stored in VOB)
            ! and direction is initially set to to zero and stored in UOB
            ! in wpc = 1 stage. u and v are computed in program wpc = 10 (OIQC).
            if ( kx /= 283 .or. ( kx == 283 .and. wpc == 10 ) ) then
               call da_earth_2_model_wind(obs(5,k), obs(6,k), &
                  platform % each (k) % u % inv, &
                  platform % each (k) % v % inv, &
                  platform%info%lon)
            end if
            if (platform % each (k) % u % inv == 0.0 .and. platform % each (k) % v % inv == 0.0) then
               platform % each (k) % u % inv = missing_r  
               platform % each (k) % v % inv = missing_r  
               platform % each (k) % u % qc  = missing_data
               platform % each (k) % v % qc  = missing_data
            end if
         end if

         if (qqm<=qflag .and. qqm>=0 .and. obs(2,k)>0.0 .and. obs(2,k)<r8bfms) then
            platform % each (k) % q % inv =obs(2,k)
            if ( obs(1,k) >= 300.0 .and. obs(1,k) < r8bfms ) then   ! do not use moisture above 300 hPa
               platform % each (k) % q % qc  =qqm
            end if
            platform % each (k) % q % error = qoe
         end if

         if ( zqm <= qflag .and. zqm >= 0 .and. &
              obs(4,k) < r8bfms )then
            platform % each (k) % height  = obs(4,k)
            platform % each (k) % height_qc =zqm
         end if

         if ( pqm <= qflag .and. pqm >= 0 .and. &
              obs(1,k) > 0.0 .and. obs(1,k) < r8bfms )then
            platform % each (k) % p % inv =obs(1,k)*100.0  ! convert to Pa
            platform % each (k) % p % qc  =pqm
            platform % each (k) % p % error =poe
         end if
      end do
      !!$OMP END PARALLEL DO

      ! assign u,v,t,q obs errors for synop and metar
      if ( t29 == 512 .or. t29 == 511 .or. t29 == 514 ) then
         qqm=nint(qms(2,1))
         tqm=nint(qms(3,1))
         wqm=nint(qms(5,1))
         toe = min(err_t, oes(3,1))
         woe = min(err_uv, oes(5,1))
         qoe = min(err_q, oes(2,1)*10.0)  ! convert to % from PREPBUFR percent divided by 10
         if ( wqm == 8 .or. wqm  == 9 .or. wqm  == 15) then
            platform%each(1)%u%qc  = 88
            platform%each(1)%v%qc  = 88
            if ( use_errtable ) then
               platform%each(1)%u%error = woe
               platform%each(1)%v%error = woe
            else
               platform%each(1)%u%error = 1.1
               platform%each(1)%v%error = 1.1
            end if
            ! Convert earth wind to model wind.
            call da_earth_2_model_wind(obs(5,1), obs(6,1),     &
               platform%each(1)%u%inv, platform%each(1)%v%inv, &
               platform%info%lon)
            if ( platform%each(1)%u%inv == 0.0 .and. platform%each(1)%v%inv == 0.0 ) then
               platform%each(1)%u%inv = missing_r  
               platform%each(1)%v%inv = missing_r  
               platform%each(1)%u%qc  = missing_data
               platform%each(1)%v%qc  = missing_data
            end if
         end if
         if ( tqm == 8 .or. tqm == 9 .or. tqm == 15 ) then
            if ( obs(3,1) < r8bfms ) then
            platform%each(1)%t%inv = obs(3,1)
            platform%each(1)%t%qc  = 88
            if ( use_errtable ) then
               platform%each(1)%t%error = toe
            else
               platform%each(1)%t%error = 2.0
            end if
            end if
         end if
         if ( qqm == 8 .or. qqm == 9 .or. qqm == 15 ) then
            if ( obs(2,1) > 0.0 .and. obs(2,1) < r8bfms ) then
            platform%each(1)%q%inv = obs(2,1)
            platform%each(1)%q%qc  = 88
            if ( use_errtable ) then
               platform%each(1)%q%error = qoe  ! RH percent
            else
               platform%each(1)%q%error = 10   ! RH percent
            end if
            end if
         end if
      end if
      ! assign tpw obs errors for gpspw
      if ( t29 == 74 ) then
         if ( pwq == 8 .or. pwq  == 9 .or. pwq  == 15) then
            if ( obs(7,1) > 0.0 .and. obs(7,1) < r8bfms ) then
            platform%loc%pw%inv   = obs(7,1) * 0.1    ! convert to cm
            platform%loc%pw%qc    = 88
            if ( use_errtable ) then
               platform%loc%pw%error = pwe
            else
               platform%loc%pw%error = 0.2     ! hardwired to 0.2 cm
            end if
            end if
         end if
      end if

      nlevels    = platform%info%levels

      if (nlevels > max_ob_levels) then
         nlevels = max_ob_levels

         write(unit=stderr, fmt='(/a/)') 'Warning: Too many levels.'

         write(unit=stderr, fmt='(/2a/2a/2x,a,2f8.2,2(2x,a,f9.2)/2(2x,a,i4)/)') &
            'Subset:   ', platform%info%name(1:8), &
            'Platfrom: ', trim(platform%info%platform), &
            'Loc(lat, lon): ', platform%info%lat, platform%info%lon, &
            'elv:   ', platform%info%elv, &
            'pstar: ', platform%info%pstar, &
            'level: ', platform%info%levels, &
            'kx:    ', kx
      else if ( nlevels < 1 ) then
         if ( (kx /= 164) .and. (kx /= 174) .and.                   &
              (kx /= 165) .and. (kx /= 175) .and. (kx /= 74) ) then
            write(unit=stderr, fmt='(/a/)') &
               'Warning: Too few levels.'
   
            write(unit=stderr, fmt='(/2a/2a/2x,a,2f8.2,2(2x,a,f9.2)/3(2x,a,i4)/)') &
               'Subset:   ', platform%info%name(1:8), &
               'Platfrom: ', trim(platform%info%platform), &
               'Loc(lat, lon): ', platform%info%lat, platform%info%lon, &
               'elv:   ', platform%info%elv, &
               'pstar: ', platform%info%pstar, &
               'level: ', platform%info%levels, &
               'kx:    ', kx, &
               't29:   ', t29

            cycle reports
         end if
      end if

     if (num_fgat_time > 1 ) then !for 4dvar

      if (.not. associated(head)) then
         nullify ( head )
         allocate ( head )
         allocate ( head%platform_BUFR%each(1:nlevels) )
         nullify ( head%next )
         plink => head
      else
         allocate ( plink%next )
         plink => plink%next
         allocate ( plink%platform_BUFR%each(1:nlevels) )
         nullify ( plink%next )
      end if

      ! Recalculate the tdiff  based on the base time at each slot
      platform%info%dhr = ( obs_time - & 
             time_slots(0) - (ifgat-1)*(time_slots(num_fgat_time)-time_slots(0))/float(num_fgat_time-1) & 
                          ) / 60.0
      plink%platform_BUFR%info = platform%info
      plink%platform_BUFR%loc = platform%loc
      plink%platform_BUFR%each(1:nlevels) = platform%each(1:nlevels)
      plink%ifgat_BUFR=ifgat
      plink%fm_BUFR=0
      plink%nlevels_BUFR=nlevels
      plink%kx_BUFR=kx
      plink%t29_BUFR=t29
      plink%pco_BUFR=pco

      num_p=num_p+1

   else  !3dvar
          
      tdiff = abs(platform%info%dhr-0.02)
      dlat_earth = platform%info%lat
      dlon_earth = platform%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad
      ndup = 1
      if (global .and. &
         (platform%loc%i < ids .or. platform%loc%i >= ide)) ndup= 2
      if (test_transforms) ndup = 1

      ! It is possible that logic for counting obs is incorrect for the
      ! global case with >1 MPI tasks due to obs duplication, halo, etc.
      ! TBH:  20050913
      dup_loop: do n = 1, ndup
         select case(t29)
         case (11, 12, 13, 22, 23, 31)
            select case (kx)
            case (120, 122, 132, 220, 222, 232) ;         ! Sound
               if (.not.use_soundobs) cycle reports
               if (n==1) iv%info(sound)%ntotal     = iv%info(sound)%ntotal + 1
               if (n==1) iv%info(sonde_sfc)%ntotal = iv%info(sonde_sfc)%ntotal + 1
               
               if (outside) cycle reports
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(sound,dlat_earth,dlon_earth,crit,iv%info(sound)%nlocal,itx,1,itt,iout,iuse)
                  call map2grids_conv(sonde_sfc,dlat_earth,dlon_earth,crit,iv%info(sonde_sfc)%nlocal,itx,1,itt,iout,iuse)

                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(sound)%nlocal     = iv%info(sound)%nlocal + 1
                  iv%info(sonde_sfc)%nlocal = iv%info(sound)%nlocal
               end if
               fm = 35
               
            case (221) ;                   ! Pilot
               if (.not.use_pilotobs) cycle reports
               if (n==1) iv%info(pilot)%ntotal = iv%info(pilot)%ntotal + 1
               if (outside) cycle reports
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(pilot,dlat_earth,dlon_earth,crit,iv%info(pilot)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(pilot)%nlocal = iv%info(pilot)%nlocal + 1
               end if
               fm = 32
               
            case default
               exit dup_loop
            end select

         case (41)
            ! case (130:131, 133, 230:231, 233) ; ! Airep
               if (.not.use_airepobs) cycle reports
               if (n==1) iv%info(airep)%ntotal = iv%info(airep)%ntotal + 1
               if (outside)  cycle reports
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(airep,dlat_earth,dlon_earth,crit,iv%info(airep)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(airep)%nlocal = iv%info(airep)%nlocal + 1
               end if
               fm = 42
               
         case (522, 523);        ! Ships
               if (.not.use_shipsobs) cycle reports
               if (n==1) iv%info(ships)%ntotal = iv%info(ships)%ntotal + 1
               if (outside) cycle reports
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(ships,dlat_earth,dlon_earth,crit,iv%info(ships)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(ships)%nlocal = iv%info(ships)%nlocal + 1
               end if
               fm = 13
               
         case (531, 532, 561, 562) ;          ! Buoy  
               if (.not.use_buoyobs) cycle reports
               if (n==1) iv%info(buoy)%ntotal = iv%info(buoy)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(buoy,dlat_earth,dlon_earth,crit,iv%info(buoy)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(buoy)%nlocal = iv%info(buoy)%nlocal + 1
               end if
               fm = 18
               
         case (511, 514)
            ! case (181, 281) ;                   ! Synop
               if (.not.use_synopobs) cycle reports
               
               if (n==1) iv%info(synop)%ntotal = iv%info(synop)%ntotal + 1
               if (outside) cycle reports
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(synop,dlat_earth,dlon_earth,crit,iv%info(synop)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(synop)%nlocal = iv%info(synop)%nlocal + 1
               end if
               fm = 12
               
         case (512)
            ! case (187, 287) ;                        ! Metar
               if (.not.use_metarobs) cycle reports
               
               if (n==1) iv%info(metar)%ntotal = iv%info(metar)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(metar,dlat_earth,dlon_earth,crit,iv%info(metar)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(metar)%nlocal = iv%info(metar)%nlocal + 1
               end if
               fm = 15
               
         case (63)
            ! case (242:246, 252:253, 255) ;         ! Geo. CMVs
               if (.not.use_geoamvobs) cycle reports
               
               if (n==1) iv%info(geoamv)%ntotal = iv%info(geoamv)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(geoamv,dlat_earth,dlon_earth,crit,iv%info(geoamv)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(geoamv)%nlocal = iv%info(geoamv)%nlocal + 1
               end if
               fm = 88
               
         case (581, 582, 583, 584)      ! ERS 581, QuikSCAT 582, WindSat 583, ASCAT 584
               if (.not.use_qscatobs) cycle reports
               
               if (n==1) iv%info(qscat)%ntotal = iv%info(qscat)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(qscat,dlat_earth,dlon_earth,crit,iv%info(qscat)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(qscat)%nlocal = iv%info(qscat)%nlocal + 1
               end if
               fm = 281
               
         case (74)       ! GPS PW
               if (.not.use_gpspwobs) cycle reports
               
               if (n==1) iv%info(gpspw)%ntotal = iv%info(gpspw)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(gpspw,dlat_earth,dlon_earth,crit,iv%info(gpspw)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(gpspw)%nlocal = iv%info(gpspw)%nlocal + 1
               end if
               fm = 111
               
         case (71, 73, 75, 76, 77)    ! Profiler
               if (.not.use_profilerobs) cycle reports

               if (n==1) iv%info(profiler)%ntotal = iv%info(profiler)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(profiler,dlat_earth,dlon_earth,crit,iv%info(profiler)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(profiler)%nlocal = iv%info(profiler)%nlocal + 1
               end if
               fm = 132
               
         case (571, 65)
              if (.not. use_ssmiretrievalobs) cycle reports
               
               if (n==1) iv%info(ssmi_rv)%ntotal = iv%info(ssmi_rv)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(ssmi_rv,dlat_earth,dlon_earth,crit,iv%info(ssmi_rv)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(ssmi_rv)%nlocal = iv%info(ssmi_rv)%nlocal + 1
               end if
               fm = 125      ! ssmi wind speed & tpw
               
         case default 
            select case (kx)
            case (111 , 210)    ;         !  Tropical Cyclone Bogus
               ! Note Tropical cyclone Bougus is given type 135 in Obs-ascii
               if (.not.use_bogusobs) cycle reports
               
               if (n==1) iv%info(bogus)%ntotal = iv%info(bogus)%ntotal + 1
               if (outside) cycle reports
               
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(bogus,dlat_earth,dlon_earth,crit,iv%info(bogus)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     cycle reports
                  end if
               else
                  iv%info(bogus)%nlocal = iv%info(bogus)%nlocal + 1
               end if
               fm = 135
               
            case default
               if ( print_detail_obs ) then
                  write(unit=message(1), fmt='(a, 2i12)') &
                     'unsaved obs found with kx & t29= ',kx,t29
                  call da_warning("da_read_obs_bufr.inc",1108,message(1:1))
               end if
               exit dup_loop
            end select
         end select

      obs_index=fm_index(fm)
      iv%info(obs_index)%max_lev = max(iv%info(obs_index)%max_lev, nlevels)    
      
      if (.not. associated(head)) then
         nullify ( head )
         allocate ( head )
         allocate ( head%platform_BUFR%each(1:nlevels) )
         nullify ( head%next )
         plink => head
      else
         allocate ( plink%next )
         plink => plink%next
         allocate ( plink%platform_BUFR%each(1:nlevels) )
         nullify ( plink%next )
      end if

      plink%platform_BUFR%info = platform%info
      plink%platform_BUFR%loc = platform%loc
      plink%platform_BUFR%each(1:nlevels) = platform%each(1:nlevels)
      plink%fm_BUFR=fm
      plink%ifgat_BUFR=ifgat
      plink%nlevels_BUFR=nlevels
      plink%kx_BUFR=kx
      plink%t29_BUFR=t29
      plink%pco_BUFR=pco
      
      num_p=num_p+1
         end do dup_loop
       end if !3dvar and 4dvar
    end do reports
    
call closbf(iunit)
close(iunit)
if ( use_errtable ) then
  close(junit)
  call da_free_unit(junit)
end if

end do bufrfile   


!  3.0 Thinning based on FGAT
!      Only for 4dvar
!--------------------------------------------------------------

 if (num_fgat_time > 1 ) then

    do kk=1,num_fgat_time
      if ( thin_conv ) then
         do n = 1, num_ob_indexes
            call cleangrids_conv(n)
         end do
      end if

       plink => head    
       reports2: do ii=1,num_p

       if (plink%ifgat_BUFR /= kk) then  !sort iv
            if ( .not. associated(plink%next) ) exit reports2
            plink => plink%next
            cycle reports2
      else
                 ! for thinning
      tdiff = abs(plink%platform_BUFR%info%dhr-0.02)
      dlat_earth = plink%platform_BUFR%info%lat
      dlon_earth = plink%platform_BUFR%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad
      call da_llxy (plink%platform_BUFR%info, plink%platform_BUFR%loc,outside, outside_all)

 ! Loop over duplicating obs for global
      ndup = 1
      if (global .and. &
         (plink%platform_BUFR%loc%i < ids .or. plink%platform_BUFR%loc%i >= ide)) ndup= 2
      if (test_transforms) ndup = 1

      ! It is possible that logic for counting obs is incorrect for the
      ! global case with >1 MPI tasks due to obs duplication, halo, etc.
      ! TBH:  20050913
      dup_loop2: do n = 1, ndup
         select case(plink%t29_BUFR)
         case (11, 12, 13, 22, 23, 31)
            select case (plink%kx_BUFR)
            case (120, 122, 132, 220, 222, 232) ;         ! Sound
               if (.not.use_soundobs) then
                      plink => plink%next
                      cycle reports2
               end if
               if (n==1) iv%info(sound)%ntotal     = iv%info(sound)%ntotal + 1
               if (n==1) iv%info(sonde_sfc)%ntotal = iv%info(sonde_sfc)%ntotal + 1
               
               if (outside) then
                    plink => plink%next
                    cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(sound,dlat_earth,dlon_earth,crit,iv%info(sound)%nlocal,itx,1,itt,iout,iuse)
                  call map2grids_conv(sonde_sfc,dlat_earth,dlon_earth,crit,iv%info(sonde_sfc)%nlocal,itx,1,itt,iout,iuse)

                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(sound)%nlocal     = iv%info(sound)%nlocal + 1
                  iv%info(sonde_sfc)%nlocal = iv%info(sound)%nlocal
               end if
               fm = 35
               
            case (221) ;                   ! Pilot
               if (.not.use_pilotobs) then
                     plink => plink%next
                     cycle reports2
               end if
               if (n==1) iv%info(pilot)%ntotal = iv%info(pilot)%ntotal + 1
               if (outside) then 
                   plink => plink%next
                   cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(pilot,dlat_earth,dlon_earth,crit,iv%info(pilot)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(pilot)%nlocal = iv%info(pilot)%nlocal + 1
               end if
               fm = 32
               
            case default
               exit dup_loop2
            end select

         case (41)
            ! case (130:131, 133, 230:231, 233) ; ! Airep
               if (.not.use_airepobs) then
                      plink => plink%next
                      cycle reports2
               end if
               if (n==1) iv%info(airep)%ntotal = iv%info(airep)%ntotal + 1
               if (outside) then 
                    plink => plink%next
                    cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(airep,dlat_earth,dlon_earth,crit,iv%info(airep)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(airep)%nlocal = iv%info(airep)%nlocal + 1
               end if
               fm = 42
               
         case (522, 523);        ! Ships
               if (.not.use_shipsobs) then
                   plink => plink%next
                   cycle reports2
               end if
               if (n==1) iv%info(ships)%ntotal = iv%info(ships)%ntotal + 1
               if (outside) then
                  plink => plink%next
                  cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(ships,dlat_earth,dlon_earth,crit,iv%info(ships)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(ships)%nlocal = iv%info(ships)%nlocal + 1
               end if
               fm = 13
               
         case (531, 532, 561, 562) ;          ! Buoy  
               if (.not.use_buoyobs) then
                   plink => plink%next
                   cycle reports2
               end if
               if (n==1) iv%info(buoy)%ntotal = iv%info(buoy)%ntotal + 1
               if (outside) then
                     plink => plink%next
                     cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(buoy,dlat_earth,dlon_earth,crit,iv%info(buoy)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(buoy)%nlocal = iv%info(buoy)%nlocal + 1
               end if
               fm = 18
               
         case (511, 514)
            ! case (181, 281) ;                   ! Synop
               if (.not.use_synopobs) then
                   plink => plink%next
                   cycle reports2
               end if
               if (n==1) iv%info(synop)%ntotal = iv%info(synop)%ntotal + 1
               if (outside) then
                     plink => plink%next
                     cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(synop,dlat_earth,dlon_earth,crit,iv%info(synop)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(synop)%nlocal = iv%info(synop)%nlocal + 1
               end if
               fm = 12
               
         case (512)
            ! case (187, 287) ;                        ! Metar
               if (.not.use_metarobs) then
                    plink => plink%next
                    cycle reports2
               end if
               if (n==1) iv%info(metar)%ntotal = iv%info(metar)%ntotal + 1
               if (outside) then 
                    plink => plink%next
                    cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(metar,dlat_earth,dlon_earth,crit,iv%info(metar)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(metar)%nlocal = iv%info(metar)%nlocal + 1
               end if
               fm = 15
               
         case (63)
            ! case (242:246, 252:253, 255) ;         ! Geo. CMVs
               if (.not.use_geoamvobs) then 
                       plink => plink%next
                       cycle reports2
               end if
               if (n==1) iv%info(geoamv)%ntotal = iv%info(geoamv)%ntotal + 1
               if (outside) then 
                  plink => plink%next
                  cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(geoamv,dlat_earth,dlon_earth,crit,iv%info(geoamv)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(geoamv)%nlocal = iv%info(geoamv)%nlocal + 1
               end if
               fm = 88
               
         case (581, 582, 583, 584)      ! ERS 581, QuikSCAT 582, WindSat 583, ASCAT 584
               if (.not.use_qscatobs) then
                      plink => plink%next
                      cycle reports2
               end if
               if (n==1) iv%info(qscat)%ntotal = iv%info(qscat)%ntotal + 1
               if (outside) then
                   plink => plink%next
                   cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(qscat,dlat_earth,dlon_earth,crit,iv%info(qscat)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(qscat)%nlocal = iv%info(qscat)%nlocal + 1
               end if
               fm = 281
               
         case (74)       ! GPS PW
               if (.not.use_gpspwobs) then 
                    plink => plink%next
                    cycle reports2
               end if
               if (n==1) iv%info(gpspw)%ntotal = iv%info(gpspw)%ntotal + 1
               if (outside) then 
                   plink => plink%next
                   cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(gpspw,dlat_earth,dlon_earth,crit,iv%info(gpspw)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(gpspw)%nlocal = iv%info(gpspw)%nlocal + 1
               end if
               fm = 111
               
         case (71, 73, 75, 76, 77)    ! Profiler
               if (.not.use_profilerobs) then
                     plink => plink%next
                     cycle reports2
               end if
               if (n==1) iv%info(profiler)%ntotal = iv%info(profiler)%ntotal + 1
               if (outside) then
                  plink => plink%next
                  cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(profiler,dlat_earth,dlon_earth,crit,iv%info(profiler)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(profiler)%nlocal = iv%info(profiler)%nlocal + 1
               end if
               fm = 132
               
         case (571, 65)
              if (.not. use_ssmiretrievalobs) then
                         plink => plink%next
                         cycle reports2
               end if
               if (n==1) iv%info(ssmi_rv)%ntotal = iv%info(ssmi_rv)%ntotal + 1
               if (outside) then
                    plink => plink%next
                    cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(ssmi_rv,dlat_earth,dlon_earth,crit,iv%info(ssmi_rv)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(ssmi_rv)%nlocal = iv%info(ssmi_rv)%nlocal + 1
               end if
               fm = 125      ! ssmi wind speed & tpw
               
         case default 
            select case (plink%kx_BUFR)
            case (111 , 210)    ;         !  Tropical Cyclone Bogus
               ! Note Tropical cyclone Bougus is given type 135 in Obs-ascii
               if (.not.use_bogusobs) then 
                      plink => plink%next
                      cycle reports2
               end if
               if (n==1) iv%info(bogus)%ntotal = iv%info(bogus)%ntotal + 1
               if (outside) then 
                  plink => plink%next
                  cycle reports2
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(bogus,dlat_earth,dlon_earth,crit,iv%info(bogus)%nlocal,itx,1,itt,iout,iuse)
                  if ( .not. iuse ) then
                     num_thinned = num_thinned + 1
                     plink => plink%next
                     cycle reports2
                  end if
               else
                  iv%info(bogus)%nlocal = iv%info(bogus)%nlocal + 1
               end if
               fm = 135
               
            case default
               if ( print_detail_obs ) then
                  write(unit=message(1), fmt='(a, 2i12)') &
                     'unsaved obs found with kx & t29= ',plink%kx_BUFR,plink%t29_BUFR
                  call da_warning("da_read_obs_bufr.inc",1518,message(1:1))
               end if
               exit dup_loop2
            end select
         end select

      obs_index=fm_index(fm)
      iv%info(obs_index)%max_lev = max(iv%info(obs_index)%max_lev, plink%nlevels_BUFR)    
      
      plink%fm_BUFR=fm

      end do dup_loop2
     end if  !sort iv

     plink => plink%next
     if ( .not. associated(plink) ) exit reports2

   end do reports2

end do !kk

end if 

!  4.0 Allocate iv 
!  
!--------------------------------------------------------------

   iv%info(synop)%max_lev     = 1
   iv%info(metar)%max_lev     = 1
   iv%info(ships)%max_lev     = 1
   iv%info(buoy)%max_lev      = 1
   iv%info(sonde_sfc)%max_lev = 1
  
   call da_allocate_observations (iv) 

!  5.0 Transfer p structure into iv structure
!      Also sort iv structure to FGAT time bins
!--------------------------------------------------------------       
    do kk=1,num_fgat_time

       if ( thin_conv ) then
          do n = 1, num_ob_indexes
             call cleangrids_conv(n)
          end do
       end if

       plink => head 
   
      reports3: do ii=1,num_p
         
     if (plink%ifgat_BUFR /=  kk) then  !sort iv
            if ( .not. associated(plink%next) ) exit reports3
            plink => plink%next
            cycle reports3
      else
      ! for thinning
      tdiff = abs(plink%platform_BUFR%info%dhr-0.02)
      dlat_earth = plink%platform_BUFR%info%lat
      dlon_earth = plink%platform_BUFR%info%lon
      if (dlon_earth < 0.0) dlon_earth = dlon_earth + 360.0
      if (dlon_earth >= 360.0) dlon_earth = dlon_earth - 360.0
      dlat_earth = dlat_earth * deg2rad
      dlon_earth = dlon_earth * deg2rad
      call da_llxy (plink%platform_BUFR%info, plink%platform_BUFR%loc,outside, outside_all)
      ndup = 1
      if (global .and. &
         (plink%platform_BUFR%loc%i < ids .or. plink%platform_BUFR%loc%i >= ide)) ndup= 2
      if (test_transforms) ndup = 1
      dup_loop3: do n = 1, ndup

         select case(plink%t29_BUFR)
           case (11, 12, 13, 22, 23, 31)
            select case (plink%kx_BUFR)
              case (120, 122, 132, 220, 222, 232) ;         ! Sound
               if (.not.use_soundobs) then
                      plink => plink%next
                      cycle reports3
               end if
                 if (n==1) ntotal(sound) = ntotal(sound) + 1 
                 if (n==1) ntotal(sonde_sfc) = ntotal(sonde_sfc) + 1
                 if (outside) then
                    plink => plink%next
                    cycle reports3
               end if
                if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(sound,dlat_earth,dlon_earth,crit,nlocal(sound),itx,1,itt,ilocal(sound),iuse)
                  call map2grids_conv(sonde_sfc,dlat_earth,dlon_earth,crit,nlocal(sonde_sfc),itx,1,itt,ilocal(sonde_sfc),iuse)
                   if ( .not. iuse ) then
                     plink => plink%next
                     cycle reports3
                  end if
               else
                  nlocal(sound) = nlocal(sound) + 1
                  nlocal(sonde_sfc) = nlocal(sound)
                  ilocal(sound) = nlocal(sound)
                  ilocal(sonde_sfc) = ilocal(sound)
               end if

               platform_name ='FM-35 TEMP  '
               
               old_nlevels = plink%nlevels_BUFR

               ! Search to see if we have surface obs.

               surface_level = 0

               do i = 1, plink%nlevels_BUFR
                  ! if (elevation and height are the same, it is surface)
                  if (abs(plink%platform_BUFR%info%elv - &
                     plink%platform_BUFR%each(i)%height) < 0.1) then
                     surface_level = i

                     ! Save surface pressure.
                     iv%sonde_sfc(ilocal(sonde_sfc))%h = plink%platform_BUFR%each(i)%height
                     iv%sonde_sfc(ilocal(sonde_sfc))%u = plink%platform_BUFR%each(i)%u
                     iv%sonde_sfc(ilocal(sonde_sfc))%v = plink%platform_BUFR%each(i)%v
                     iv%sonde_sfc(ilocal(sonde_sfc))%t = plink%platform_BUFR%each(i)%t
                     iv%sonde_sfc(ilocal(sonde_sfc))%q = plink%platform_BUFR%each(i)%q
                     iv%sonde_sfc(ilocal(sonde_sfc))%p = plink%platform_BUFR%each(i)%p
                    
                     exit
                  end if
               end do

               ! processing the sound_sfc data:

               if (surface_level > 0) then
                  plink%nlevels_BUFR = plink%nlevels_BUFR - 1
               else   !missing surface data
                  iv%sonde_sfc(ilocal(sonde_sfc))%h = missing_r
                  iv%sonde_sfc(ilocal(sonde_sfc))%u%inv   = missing_r
                  iv%sonde_sfc(ilocal(sonde_sfc))%u%qc    = missing_data
                  iv%sonde_sfc(ilocal(sonde_sfc))%u%error = abs(missing_r)
                  iv%sonde_sfc(ilocal(sonde_sfc))%v = iv%sonde_sfc(ilocal(sonde_sfc))%u
                  iv%sonde_sfc(ilocal(sonde_sfc))%t = iv%sonde_sfc(ilocal(sonde_sfc))%u
                  iv%sonde_sfc(ilocal(sonde_sfc))%p = iv%sonde_sfc(ilocal(sonde_sfc))%u
                  iv%sonde_sfc(ilocal(sonde_sfc))%q = iv%sonde_sfc(ilocal(sonde_sfc))%u
                  
               end if

               if (plink%nlevels_BUFR > 0) then

                  if ( ilocal(sound) == nlocal(sound) .or. &
                       .not. associated(iv%sound(ilocal(sound))%h) ) then
                     allocate (iv%sound(ilocal(sound))%h(1:iv%info(sound)%max_lev))
                     allocate (iv%sound(ilocal(sound))%p(1:iv%info(sound)%max_lev))
                     allocate (iv%sound(ilocal(sound))%u(1:iv%info(sound)%max_lev))
                     allocate (iv%sound(ilocal(sound))%v(1:iv%info(sound)%max_lev))
                     allocate (iv%sound(ilocal(sound))%t(1:iv%info(sound)%max_lev))
                     allocate (iv%sound(ilocal(sound))%q(1:iv%info(sound)%max_lev))
                   endif  

                  j = 0
                  do i = 1, old_nlevels
                     if (i == surface_level) cycle
                     j=j+1
                     iv%sound(ilocal(sound))%h(j) = plink%platform_BUFR%each(i)%height
                     iv%sound(ilocal(sound))%p(j) = plink%platform_BUFR%each(i)%p%inv
                     iv%sound(ilocal(sound))%u(j) = plink%platform_BUFR%each(i)%u
                     iv%sound(ilocal(sound))%v(j) = plink%platform_BUFR%each(i)%v
                     iv%sound(ilocal(sound))%t(j) = plink%platform_BUFR%each(i)%t
                     iv%sound(ilocal(sound))%q(j) = plink%platform_BUFR%each(i)%q
                     
                  end do
               end if

            case (221) ;           ! Pilot 
                if (.not.use_pilotobs) then
                     plink => plink%next
                     cycle reports3
               end if
               if (n==1) ntotal(pilot) = ntotal(pilot) + 1
               if (outside) then 
                   plink => plink%next
                   cycle reports3
               end if
               if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(pilot,dlat_earth,dlon_earth,crit,nlocal(pilot),itx,1,itt,ilocal(pilot),iuse)
                  if ( .not. iuse ) then
                     plink => plink%next
                     cycle reports3
                  end if
               else
                  nlocal(pilot) = nlocal(pilot) + 1
                  ilocal(pilot) = nlocal(pilot)
               end if
                platform_name='FM-32 PILOT '

               if ( ilocal(pilot) == nlocal(pilot) ) then
                  allocate (iv%pilot(ilocal(pilot))%h(1:iv%info(pilot)%max_lev))
                  allocate (iv%pilot(ilocal(pilot))%p(1:iv%info(pilot)%max_lev))
                  allocate (iv%pilot(ilocal(pilot))%u(1:iv%info(pilot)%max_lev))
                  allocate (iv%pilot(ilocal(pilot))%v(1:iv%info(pilot)%max_lev))
               endif   

               do i = 1, plink%nlevels_BUFR
                  iv%pilot(ilocal(pilot))%h(i) = plink%platform_BUFR%each(i)%height
                  iv%pilot(ilocal(pilot))%p(i) = plink%platform_BUFR%each(i)%p%inv
                  iv%pilot(ilocal(pilot))%u(i) = plink%platform_BUFR%each(i)%u
                  iv%pilot(ilocal(pilot))%v(i) = plink%platform_BUFR%each(i)%v
                  
               end do
                case default
               exit dup_loop3
            end select 
         case (41)
             if (.not.use_airepobs) then
                      plink => plink%next
                      cycle reports3
               end if

            if (n==1) ntotal(airep) = ntotal(airep) + 1
             if (outside) then 
                    plink => plink%next
                    cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(airep,dlat_earth,dlon_earth,crit,nlocal(airep),itx,1,itt,ilocal(airep),iuse)
               if ( .not. iuse ) then
                    plink => plink%next
                    cycle reports3
                  end if
            else
               nlocal(airep) = nlocal(airep) + 1
               ilocal(airep) = nlocal(airep)
            end if
            platform_name ='FM-97 AIREP '
            if ( ilocal(airep) == nlocal(airep) ) then
               allocate (iv%airep(ilocal(airep))%h(1:iv%info(airep)%max_lev))
               allocate (iv%airep(ilocal(airep))%p(1:iv%info(airep)%max_lev))
               allocate (iv%airep(ilocal(airep))%u(1:iv%info(airep)%max_lev))
               allocate (iv%airep(ilocal(airep))%v(1:iv%info(airep)%max_lev))
               allocate (iv%airep(ilocal(airep))%t(1:iv%info(airep)%max_lev))
               allocate (iv%airep(ilocal(airep))%q(1:iv%info(airep)%max_lev))
            end if
            do i = 1, plink%nlevels_BUFR
               iv % airep (ilocal(airep)) % h(i) = plink%platform_BUFR % each(i) % height
               iv % airep (ilocal(airep)) % p(i) = plink%platform_BUFR % each(i) % p % inv
               iv % airep (ilocal(airep)) % u(i) = plink%platform_BUFR % each(i) % u
               iv % airep (ilocal(airep)) % v(i) = plink%platform_BUFR % each(i) % v
               iv % airep (ilocal(airep)) % t(i) = plink%platform_BUFR % each(i) % t
               iv % airep (ilocal(airep)) % q(i) = plink%platform_BUFR % each(i) % q
            end do

         case (522,523);        ! Ships
          if (.not.use_shipsobs) then
                   plink => plink%next
                   cycle reports3
           end if

            if (n==1) ntotal(ships) = ntotal(ships) + 1
             if (outside) then
                  plink => plink%next
                  cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(ships,dlat_earth,dlon_earth,crit,nlocal(ships),itx,1,itt,ilocal(ships),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if
            else
               nlocal(ships) = nlocal(ships) + 1
               ilocal(ships) = nlocal(ships)
            end if
            platform_name ='FM-13 SHIP  '
            
            iv % ships (ilocal(ships)) % h = plink%platform_BUFR % each(1) % height
            iv % ships (ilocal(ships)) % u = plink%platform_BUFR % each(1) % u
            iv % ships (ilocal(ships)) % v = plink%platform_BUFR % each(1) % v
            iv % ships (ilocal(ships)) % t = plink%platform_BUFR % each(1) % t
            iv % ships (ilocal(ships)) % p = plink%platform_BUFR % each(1) % p
            iv % ships (ilocal(ships)) % q = plink%platform_BUFR % each(1) % q
            

          case (531, 532, 561, 562) ;          ! Buoy
           if (.not.use_buoyobs) then
                   plink => plink%next
                   cycle reports3
               end if
            if (n==1)  ntotal(buoy) = ntotal(buoy) + 1
             if (outside) then
                plink => plink%next
                cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(buoy,dlat_earth,dlon_earth,crit,nlocal(buoy),itx,1,itt,ilocal(buoy),iuse)
               if ( .not. iuse ) then
                  plink => plink%next
                  cycle reports3
               end if
            else
               nlocal(buoy) = nlocal(buoy) + 1
               ilocal(buoy) = nlocal(buoy)
            end if
            platform_name ='FM-18 BUOY  '
            
            iv%buoy(ilocal(buoy))%h = plink%platform_BUFR%each(1)%height
            iv%buoy(ilocal(buoy))%u = plink%platform_BUFR%each(1)%u
            iv%buoy(ilocal(buoy))%v = plink%platform_BUFR%each(1)%v
            iv%buoy(ilocal(buoy))%t = plink%platform_BUFR%each(1)%t
            

            if ( plink%kx_BUFR == 282 ) then  ! ATLAS BUOY, reported Pstn and Pmslp are both
                                   ! missing. Pstn is set to 1013hPa in PREPBUFR
               iv%buoy(ilocal(buoy))%p%inv   = missing_r
               iv%buoy(ilocal(buoy))%p%qc    = missing_data
               iv%buoy(ilocal(buoy))%p%error = err_p
            else
               iv%buoy(ilocal(buoy))%p = plink%platform_BUFR%each(1)%p
            end if
            iv%buoy(ilocal(buoy))%q = plink%platform_BUFR%each(1)%q

         case (511, 514)
              if (.not.use_synopobs) then
                 plink => plink%next
                 cycle reports3
               end if
            if (n==1) ntotal(synop) = ntotal(synop) + 1
            if (outside) then
               plink => plink%next
               cycle reports3
            end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(synop,dlat_earth,dlon_earth,crit,nlocal(synop),itx,1,itt,ilocal(synop),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                end if
            else
               nlocal(synop) = nlocal(synop) + 1
               ilocal(synop) = nlocal(synop)
            end if
            platform_name ='FM-12 SYNOP '
            
            iv % synop (ilocal(synop)) % h = plink%platform_BUFR % each(1) % height
            iv % synop (ilocal(synop)) % u = plink%platform_BUFR % each(1) % u
            iv % synop (ilocal(synop)) % v = plink%platform_BUFR % each(1) % v
            iv % synop (ilocal(synop)) % t = plink%platform_BUFR % each(1) % t
            iv % synop (ilocal(synop)) % p = plink%platform_BUFR % each(1) % p
            iv % synop (ilocal(synop)) % q = plink%platform_BUFR % each(1) % q
            
            if (iv % synop(ilocal(synop)) % h < plink%platform_BUFR % info % elv) then
               iv % synop(ilocal(synop)) % h = plink%platform_BUFR % info % elv
            end if

         case (512)
          if (.not.use_metarobs) then
                   plink => plink%next
                   cycle reports3
             end if
            if (n==1) ntotal(metar) = ntotal(metar) + 1
             if (outside) then 
                   plink => plink%next
                   cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(metar,dlat_earth,dlon_earth,crit,nlocal(metar),itx,1,itt,ilocal(metar),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if

            else
               nlocal(metar) = nlocal(metar) + 1
               ilocal(metar) = nlocal(metar)
            end if
               
             platform_name='FM-15 METAR '
           
            iv % metar (ilocal(metar)) % h = plink%platform_BUFR % each(1) % height
            iv % metar (ilocal(metar)) % u = plink%platform_BUFR % each(1) % u
            iv % metar (ilocal(metar)) % v = plink%platform_BUFR % each(1) % v
            iv % metar (ilocal(metar)) % t = plink%platform_BUFR % each(1) % t
            iv % metar (ilocal(metar)) % p = plink%platform_BUFR % each(1) % p
            iv % metar (ilocal(metar)) % q = plink%platform_BUFR % each(1) % q
          

         case (63)
              if (.not.use_geoamvobs) then 
                   plink => plink%next
                   cycle reports3
               end if
               if (n==1) ntotal(geoamv) = ntotal(geoamv) + 1
                if (outside) then 
                   plink => plink%next
                   cycle reports3
               end if
               if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(geoamv,dlat_earth,dlon_earth,crit,nlocal(geoamv),itx,1,itt,ilocal(geoamv),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if
            else
               nlocal(geoamv) = nlocal(geoamv) + 1
               ilocal(geoamv) = nlocal(geoamv)
            end if
               
               platform_name ='FM-88 SATOB '
            if ( ilocal(geoamv) == nlocal(geoamv) ) then
               allocate (iv%geoamv(ilocal(geoamv))%p(1:iv%info(geoamv)%max_lev))
               allocate (iv%geoamv(ilocal(geoamv))%u(1:iv%info(geoamv)%max_lev))
               allocate (iv%geoamv(ilocal(geoamv))%v(1:iv%info(geoamv)%max_lev))
               
            end if

            do i = 1, plink%nlevels_BUFR
               iv % geoamv (ilocal(geoamv)) % p(i)  = plink%platform_BUFR % each(i) % p % inv
               iv % geoamv (ilocal(geoamv)) % u(i)  = plink%platform_BUFR % each(i) % u
               iv % geoamv (ilocal(geoamv)) % v(i)  = plink%platform_BUFR % each(i) % v
               
            end do

          case (581, 582, 583, 584)      ! ERS 581, QuikSCAT 582, WindSat 583, ASCAT 584
          if (.not.use_qscatobs) then
                   plink => plink%next
                   cycle reports3
               end if
            if (n==1) ntotal(qscat) = ntotal(qscat) + 1
            if (outside) then
                   plink => plink%next
                   cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(qscat,dlat_earth,dlon_earth,crit,nlocal(qscat),itx,1,itt,ilocal(qscat),iuse)
               if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if
            else
               nlocal(qscat) = nlocal(qscat) + 1
               ilocal(qscat) = nlocal(qscat)
            end if
               
            platform_name ='FM-281 Quiks'
 
            ! prepbufr uses pressure not height, so hardwire height to 
            ! 0 (sea-level)
            iv%qscat(ilocal(qscat))%h = 0.0
            iv%qscat(ilocal(qscat))%u = plink%platform_BUFR%each(1)%u
            iv%qscat(ilocal(qscat))%v = plink%platform_BUFR%each(1)%v
            iv%qscat(ilocal(qscat))%u%error = max(plink%platform_BUFR%each(1)%u%error,1.0)
            iv%qscat(ilocal(qscat))%v%error = max(plink%platform_BUFR%each(1)%v%error,1.0)
            

         case (74)       ! GPS PW
          if (.not.use_gpspwobs) then 
                   plink => plink%next
                   cycle reports3
               end if
            if (n==1)  ntotal(gpspw) = ntotal(gpspw) + 1
             if (outside) then 
                   plink => plink%next
                   cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(gpspw,dlat_earth,dlon_earth,crit,nlocal(gpspw),itx,1,itt,ilocal(gpspw),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if
            else
               nlocal(gpspw) = nlocal(gpspw) + 1
               ilocal(gpspw) = nlocal(gpspw)
            end if
            platform_name ='FM-111 GPSPW'
 
            iv%gpspw(ilocal(gpspw))%tpw  = plink%platform_BUFR%loc%pw
            

          case (71, 73, 75, 76, 77)
                if (.not.use_profilerobs) then
                   plink => plink%next
                   cycle reports3
               end if
               if (n==1) ntotal(profiler) = ntotal(profiler) + 1
                if (outside) then
                   plink => plink%next
                   cycle reports3
               end if
               if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(profiler,dlat_earth,dlon_earth,crit,nlocal(profiler),itx,1,itt,ilocal(profiler),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if
            else
               nlocal(profiler) = nlocal(profiler) + 1
               ilocal(profiler) = nlocal(profiler)
            end if
               
               platform_name ='FM-132 PRFLR'
            if ( ilocal(profiler) == nlocal(profiler) ) then
               allocate (iv%profiler(ilocal(profiler))%h(1:iv%info(profiler)%max_lev))
               allocate (iv%profiler(ilocal(profiler))%p(1:iv%info(profiler)%max_lev))
               allocate (iv%profiler(ilocal(profiler))%u(1:iv%info(profiler)%max_lev))
               allocate (iv%profiler(ilocal(profiler))%v(1:iv%info(profiler)%max_lev))
           end if     

            do i = 1, plink%nlevels_BUFR 
               iv%profiler(ilocal(profiler))%h(i) = plink%platform_BUFR%each(i)%height
               iv%profiler(ilocal(profiler))%p(i) = plink%platform_BUFR%each(i)%p%inv
               iv%profiler(ilocal(profiler))%u(i) = plink%platform_BUFR%each(i)%u
               iv%profiler(ilocal(profiler))%v(i) = plink%platform_BUFR%each(i)%v
               
            end do

          case (571, 65)   ! SSM/I wind speed & TPW
          if (.not. use_ssmiretrievalobs) then
                   plink => plink%next
                   cycle reports3
               end if
           if (n==1)  ntotal(ssmi_rv) = ntotal(ssmi_rv) + 1
           if (outside) then
                   plink => plink%next
                   cycle reports3
               end if
            if ( thin_conv ) then
               crit = tdiff
               call map2grids_conv(ssmi_rv,dlat_earth,dlon_earth,crit,nlocal(ssmi_rv),itx,1,itt,ilocal(ssmi_rv),iuse)
                if ( .not. iuse ) then
                   plink => plink%next
                   cycle reports3
                  end if
            else
               nlocal(ssmi_rv) = nlocal(ssmi_rv) + 1
               ilocal(ssmi_rv) = nlocal(ssmi_rv)
            end if
             
            platform_name ='FM-125 SSMI '
  
            select case (plink%kx_BUFR)
            case ( 283)  ! wind speed
               do i = 1, plink%nlevels_BUFR
                  wpc = nint(plink%pco_BUFR(5,i))
                  ! if wpc == 1, UOB is set to zero, VOB is speed
                  iv%ssmi_rv(ilocal(ssmi_rv))%speed  = plink%platform_BUFR%each(i)%v
                  if ( wpc == 10 ) then
                     iv%ssmi_rv(ilocal(ssmi_rv))%speed%inv  = sqrt ( &
                        plink%platform_BUFR%each(i)%u%inv * plink%platform_BUFR%each(i)%u%inv + & 
                        plink%platform_BUFR%each(i)%v%inv * plink%platform_BUFR%each(i)%v%inv  )
                  end if
                  iv%ssmi_rv(ilocal(ssmi_rv))%tpw    = plink%platform_BUFR%loc%pw
                  
               end do
           case ( 152 )  ! tpw
              do i = 1, plink%nlevels_BUFR
                 iv%ssmi_rv(ilocal(ssmi_rv))%speed  = plink%platform_BUFR%each(i)%u
                 iv%ssmi_rv(ilocal(ssmi_rv))%tpw    = plink%platform_BUFR%loc%pw
                 
               end do
          case default
              exit dup_loop3

           end select  

          case default 
            select case (plink%kx_BUFR)
            case (111 , 210) 
          if (.not.use_bogusobs) then 
                   plink => plink%next
                   cycle reports3
               end if
               if (n==1)  ntotal(bogus) = ntotal(bogus) + 1
                if (outside) then 
                   plink => plink%next
                   cycle reports3
               end if
                if ( thin_conv ) then
                  crit = tdiff
                  call map2grids_conv(bogus,dlat_earth,dlon_earth,crit,nlocal(bogus),itx,1,itt,ilocal(bogus),iuse)
                  if ( .not. iuse ) then
                     plink => plink%next
                     cycle reports3
                  end if
               else
                  nlocal(bogus) = nlocal(bogus) + 1
                  ilocal(bogus) = nlocal(bogus)
               end if
          
               platform_name ='FM-135 TCBOG'
               if ( ilocal(bogus) == nlocal(bogus) ) then
                  allocate (iv%bogus(ilocal(bogus))%h(1:iv%info(bogus)%max_lev))
                  allocate (iv%bogus(ilocal(bogus))%p(1:iv%info(bogus)%max_lev))
                  allocate (iv%bogus(ilocal(bogus))%u(1:iv%info(bogus)%max_lev))
                  allocate (iv%bogus(ilocal(bogus))%v(1:iv%info(bogus)%max_lev))
                  allocate (iv%bogus(ilocal(bogus))%t(1:iv%info(bogus)%max_lev))
                  allocate (iv%bogus(ilocal(bogus))%q(1:iv%info(bogus)%max_lev))
               end if   

               do i = 1, plink%nlevels_BUFR
                  iv%bogus(ilocal(bogus))%h(i) = plink%platform_BUFR%each(i)%height
                  iv%bogus(ilocal(bogus))%p(i) = plink%platform_BUFR%each(i)%p%inv
                  iv%bogus(ilocal(bogus))%u(i) = plink%platform_BUFR%each(i)%u
                  iv%bogus(ilocal(bogus))%v(i) = plink%platform_BUFR%each(i)%v
                  iv%bogus(ilocal(bogus))%t(i) = plink%platform_BUFR%each(i)%t
                  iv%bogus(ilocal(bogus))%q(i) = plink%platform_BUFR%each(i)%q
                  
               end do

               iv%bogus(ilocal(bogus))%slp    = plink%platform_BUFR%loc%slp
              case default
              exit dup_loop3
             end select
         end select

         obs_index = fm_index(plink%fm_BUFR)
         iv%info(obs_index)%name(ilocal(obs_index))      = plink%platform_BUFR%info%name
         iv%info(obs_index)%platform(ilocal(obs_index))  = platform_name
         iv%info(obs_index)%id(ilocal(obs_index))        = plink%platform_BUFR%info%id
         iv%info(obs_index)%date_char(ilocal(obs_index)) = plink%platform_BUFR%info%date_char
         ! nlevels adjusted for some obs types so use that
         iv%info(obs_index)%levels(ilocal(obs_index))    = min(iv%info(obs_index)%max_lev, plink%nlevels_BUFR)
         iv%info(obs_index)%lat(:,ilocal(obs_index))     = plink%platform_BUFR%info%lat
         iv%info(obs_index)%lon(:,ilocal(obs_index))     = plink%platform_BUFR%info%lon
         iv%info(obs_index)%elv(ilocal(obs_index))       = plink%platform_BUFR%info%elv
         iv%info(obs_index)%pstar(ilocal(obs_index))     = plink%platform_BUFR%info%pstar

         iv%info(obs_index)%slp(ilocal(obs_index))           = plink%platform_BUFR%loc%slp
         iv%info(obs_index)%pw(ilocal(obs_index))            = plink%platform_BUFR%loc%pw
         iv%info(obs_index)%x(:,ilocal(obs_index))           = plink%platform_BUFR%loc%x
         iv%info(obs_index)%y(:,ilocal(obs_index))           = plink%platform_BUFR%loc%y
         iv%info(obs_index)%i(:,ilocal(obs_index))           = plink%platform_BUFR%loc%i
         iv%info(obs_index)%j(:,ilocal(obs_index))           = plink%platform_BUFR%loc%j
         iv%info(obs_index)%dx(:,ilocal(obs_index))          = plink%platform_BUFR%loc%dx
         iv%info(obs_index)%dxm(:,ilocal(obs_index))         = plink%platform_BUFR%loc%dxm
         iv%info(obs_index)%dy(:,ilocal(obs_index))          = plink%platform_BUFR%loc%dy
         iv%info(obs_index)%dym(:,ilocal(obs_index))         = plink%platform_BUFR%loc%dym
         iv%info(obs_index)%proc_domain(:,ilocal(obs_index)) = plink%platform_BUFR%loc%proc_domain

         iv%info(obs_index)%obs_global_index(ilocal(obs_index)) = iv%info(obs_index)%ntotal
         ! special case for sonde_sfc, duplicate sound info
         if (obs_index == sound) then

            iv%info(sonde_sfc)%name(ilocal(sonde_sfc))      = plink%platform_BUFR%info%name
            iv%info(sonde_sfc)%platform(ilocal(sonde_sfc))  = platform_name
            iv%info(sonde_sfc)%id(ilocal(sonde_sfc))        = plink%platform_BUFR%info%id
            iv%info(sonde_sfc)%date_char(ilocal(sonde_sfc)) = plink%platform_BUFR%info%date_char
            iv%info(sonde_sfc)%levels(ilocal(sonde_sfc))    = 1
            iv%info(sonde_sfc)%lat(:,ilocal(sonde_sfc))     = plink%platform_BUFR%info%lat
            iv%info(sonde_sfc)%lon(:,ilocal(sonde_sfc))     = plink%platform_BUFR%info%lon
            iv%info(sonde_sfc)%elv(ilocal(sonde_sfc))       = plink%platform_BUFR%info%elv
            iv%info(sonde_sfc)%pstar(ilocal(sonde_sfc))     = plink%platform_BUFR%info%pstar

            iv%info(sonde_sfc)%slp(ilocal(sonde_sfc))           = plink%platform_BUFR%loc%slp
            iv%info(sonde_sfc)%pw(ilocal(sonde_sfc))            = plink%platform_BUFR%loc%pw
            iv%info(sonde_sfc)%x(:,ilocal(sonde_sfc))           = plink%platform_BUFR%loc%x
            iv%info(sonde_sfc)%y(:,ilocal(sonde_sfc))           = plink%platform_BUFR%loc%y
            iv%info(sonde_sfc)%i(:,ilocal(sonde_sfc))           = plink%platform_BUFR%loc%i
            iv%info(sonde_sfc)%j(:,ilocal(sonde_sfc))           = plink%platform_BUFR%loc%j
            iv%info(sonde_sfc)%dx(:,ilocal(sonde_sfc))          = plink%platform_BUFR%loc%dx
            iv%info(sonde_sfc)%dxm(:,ilocal(sonde_sfc))         = plink%platform_BUFR%loc%dxm
            iv%info(sonde_sfc)%dy(:,ilocal(sonde_sfc))          = plink%platform_BUFR%loc%dy
            iv%info(sonde_sfc)%dym(:,ilocal(sonde_sfc))         = plink%platform_BUFR%loc%dym
            iv%info(sonde_sfc)%proc_domain(:,ilocal(sonde_sfc)) = plink%platform_BUFR%loc%proc_domain

            iv%info(sonde_sfc)%obs_global_index(ilocal(sonde_sfc)) =iv%info(sonde_sfc)%ntotal

         end if
 
       end do dup_loop3
     end if !sort iv

     plink => plink%next
     if ( .not. associated(plink) ) exit reports3

   end do reports3
   
   if (num_fgat_time >1 ) then
      do n = 1, num_ob_indexes
          iv%info(n)%ptotal(kk)=ntotal(n)
          iv%info(n)%plocal(kk)=nlocal(n)
      end do
   else
      do n = 1, num_ob_indexes
         ntotal(n)=iv%info(n)%ntotal
	 nlocal(n)=iv%info(n)%nlocal
         iv%info(n)%ptotal(1)=iv%info(n)%ntotal
         iv%info(n)%plocal(1)=iv%info(n)%nlocal
       end do
    end if

   ! thinning check
   if ( thin_conv ) then
      do n = 1, num_ob_indexes

        if ( ntotal(n)>0 ) then

            allocate ( in  (thinning_grid_conv(n)%itxmax) )
            allocate ( out (thinning_grid_conv(n)%itxmax) )
            do i = 1, thinning_grid_conv(n)%itxmax
               in(i) = thinning_grid_conv(n)%score_crit(i)
            end do
            ! Get minimum crit and associated processor index.

            call mpi_reduce(in, out, thinning_grid_conv(n)%itxmax, true_mpi_real, mpi_min, root, comm, ierr)
            call wrf_dm_bcast_real (out, thinning_grid_conv(n)%itxmax)
            do i = 1, thinning_grid_conv(n)%itxmax
               if ( abs(out(i)-thinning_grid_conv(n)%score_crit(i)) > 1.0E-10 ) then
                  thinning_grid_conv(n)%ibest_obs(i) = 0
               else
                  if(thinning_grid_conv(n)%score_crit(i) < 9.99e6 ) &
                  iv%info(n)%thin_ntotal=  iv%info(n)%thin_ntotal + 1
               end if
            end do

            deallocate( in  )
            deallocate( out )
         
            do j = (1+tp(n)), nlocal(n)
               found = .false.
               do i = 1, thinning_grid_conv(n)%itxmax
                  if ( thinning_grid_conv(n)%ibest_obs(i) == j .and.         &
                       thinning_grid_conv(n)%score_crit(i) < 9.99e6 ) then
                     found = .true.
                     iv%info(n)%thin_nlocal =  iv%info(n)%thin_nlocal + 1
                     exit
                  end if
               end do
               if ( .not. found ) then
                  iv%info(n)%thinned(:,j) = .true.
               end if
            end do
            
         tp(n)=nlocal(n) 
         end if
      end do
       if (num_fgat_time >1 ) then
          do n = 1, num_ob_indexes
              iv%info(n)%thin_plocal(kk)=iv%info(n)%thin_nlocal
              iv%info(n)%thin_ptotal(kk)=iv%info(n)%thin_ntotal
          end do
       else
          do n = 1, num_ob_indexes
             iv%info(n)%thin_plocal(1)=iv%info(n)%thin_nlocal
             iv%info(n)%thin_ptotal(1)=iv%info(n)%thin_ntotal
           end do
        end if

   end if  ! thin_conv

 end do  !kk cycle

if ( thin_conv ) then
    do n = 1, num_ob_indexes
        if ( ntotal(n)>0 ) then
            if ( nlocal(n) > 0 ) then
               if ( ANY(iv%info(n)%thinned(:,:)) ) then
                  call da_set_obs_missing(iv,n)  ! assign missing values to those thinned=true data
               end if
            end if
        end if
    end do
end if

write(unit=message(1),fmt='(A,4(1x,i7))') &   
      'da_read_obs_bufr: num_report, num_outside_all, num_outside_time, num_thinned: ', &
      num_report, num_outside_all, num_outside_time, num_thinned
   call da_message(message(1:1))


! Release the linked list
plink => head

DO WHILE ( ASSOCIATED (plink) )
   head => plink%next
   DEALLOCATE ( plink%platform_BUFR%each )
   DEALLOCATE ( plink )
   plink => head
ENDDO

NULLIFY (head)

   if (trace_use) call da_trace_exit("da_read_obs_bufr")

end subroutine da_read_obs_bufr



subroutine da_read_obs_bufrgpsro (iv)

   !---------------------------------------------------------------------------
   ! Purpose: Read NCEP GPSRO 1 observation file for input to wrfvar
   !---------------------------------------------------------------------------
   !   METHOD: use F90 sequantial data structure to avoid read file twice
   !            so that da_scan_obs_bufr is not necessary any more.
   !            1. read gpsro data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !                  and deallocate sequential data structure
   !
   !  HISTORY: 2009/12/17  - F90 sequantial structure  Peng Xiu
   ! 
   !----------------------------------------------------------------------------
   
   use da_control

   implicit none

   type (iv_type),             intent(inout) :: iv


   real,    parameter   :: r8bfms = 9.0D08  ! 1 missing value threshold
   integer, parameter   :: maxlevs = 500
   integer              :: iunit, iost, idate, iret, nlev1, nlev2,k,i,ii
   integer              :: num_report, num_outside_all, num_outside_time
   integer              :: iyear,imonth,iday,ihour,imin
   integer              :: ntotal, nlocal, nlev, ref_qc
   real*8               :: obs_time
   real*8               :: hdr(10)
   real                 :: ntotal_ifgat(0:num_fgat_time)
   real*8               :: rdata1(25,maxlevs), rdata2(25,maxlevs)
   real                 :: height, ref_data, ref_error
   character(len=8)     :: subset
   character(len=80)    :: hdstr
   character(len=12)    :: filename
   logical              :: outside, outside_all
   type(info_type)      :: info
   type(model_loc_type) :: loc
   character(len=5)     :: id
   character(len=19)    :: date_char
   integer              :: ifgat, kk, num_p
   integer              :: err_opt   ! 1: WRF-Var/obsproc, 2: GSI
   real                 :: erh90, erh0, err90, err0
   integer(i_kind)      :: ireadns
   integer              :: num_bufr(7), numbufr, ibufr
   type datalink_gpsro              !for gpsro data reading
       type (info_type)         :: info
       type(model_loc_type)     :: loc
       type (gpsref_type)       :: gpsref
       integer                  :: ifgat
       type (field_type)        :: slp         
       type (field_type)        :: pw
       integer                  :: obs_global_index
       type(datalink_gpsro), pointer :: next
   end type datalink_gpsro
   type(datalink_gpsro),pointer  :: head=>null(), plink=>null()
   integer(i_kind)      :: iprof

   if (trace_use) call da_trace_entry("da_read_obs_bufrgosro")

   iprof  = 0
   nlocal = 0
   ntotal = 0
   num_report       = 0
   num_outside_all  = 0
   num_outside_time = 0
   num_p=0
   ntotal_ifgat(0:num_fgat_time)=0
   
   ! open file
   !  ---------
   num_bufr(:)=0
   numbufr=0
   if (num_fgat_time>1) then
      do i=1,7
        call da_get_unit(iunit)
        write(filename,fmt='(A,I1,A)') 'gpsro0',i,'.bufr'
        open(unit   = iunit, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
        if (iost == 0) then
           numbufr=numbufr+1
           num_bufr(numbufr)=i
        else
           close (iunit)
        end if
        call da_free_unit(iunit)
      end do
   else
      numbufr=1
   end if

   if (numbufr==0) numbufr=1

bufrfile:  do ibufr=1,numbufr
   if (num_fgat_time==1) then
       filename='gpsro.bufr'
   else
       if ((numbufr ==1) .and. (num_bufr(ibufr) == 0)) then
          filename='gpsro.bufr'
      else
          write(filename,fmt='(A,I1,A)') 'gpsro0',num_bufr(ibufr),'.bufr'
      end if
   end if

   ! Use specified unit number, so we can control its endian format in environment.
   iunit = 96
   call closbf(iunit)
   open(unit   = iunit, FILE   = trim(filename), &
      iostat =  iost, form = 'unformatted', STATUS = 'OLD')
   if (iost /= 0) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iost," opening PREPBUFR obs file "//trim(filename)
      call da_warning("da_read_obs_bufrgpsro.inc",115,message(1:1))
      if (trace_use) call da_trace_exit("da_read_obs_bufrgpsro")
      return
   end if

   !--------------------------------
   ! open bufr file then check date
   !--------------------------------
   call openbf(iunit,'IN',iunit)
   call datelen(10)
   call readmg(iunit,subset,idate,iret)
   if ( iret /= 0 ) then
      write(unit=message(1),fmt='(A,I5,A)') &
         "Error",iret," reading GPSRO BUFR obs file "//trim(filename)
      call da_warning("da_read_obs_bufrgpsro.inc",129,message(1:1))
      call closbf(iunit)
      close(iunit)
      if (trace_use) call da_trace_exit("da_read_obs_bufrgpsro")
      return
   end if
   write(unit=message(1),fmt='(a,i10)') 'GPSRO BUFR file date is: ', idate
   call da_message(message(1:1))
   rewind(iunit)

   hdstr = 'YEAR MNTH DAYS HOUR MINU PCCF ELRC SAID PTID GEODU'

   reports: do while ( ireadns(iunit,subset,idate) == 0 )

      num_report = num_report + 1

      call ufbint(iunit,hdr,10,1,iret,hdstr)

      iyear  = int(hdr(1))
      imonth = int(hdr(2))
      iday   = int(hdr(3))
      ihour  = int(hdr(4))
      imin   = int(hdr(5))

      write(id, '(i3.3,i2.2)') int(hdr(8)), int(hdr(9)) ! construct id using SAID and PTID
      write(date_char, fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
         iyear, '-', imonth, '-', iday, '_', ihour, ':', imin, ':', 0

      ! check date
      call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
      if (obs_time < time_slots(0) .or.  &
          obs_time >= time_slots(num_fgat_time)) then
         num_outside_time = num_outside_time + 1
         if ( print_detail_obs ) then
            write(unit=stderr,fmt='(a,1x,i4.4,4i2.2,a)')  &
               info%id(1:5),iyear,imonth,iday,ihour,imin, '  -> outside_time'
         end if
         cycle reports
      end if
         
      if ( hdr(6) < 100.0 ) then   ! check percentage of confidence PCCF
         cycle reports
      end if
         
      call ufbseq(iunit,rdata1,25,maxlevs,nlev1,'ROSEQ1')  ! RAOC PROFILE LOCATIONS SEQUENCE
      call ufbseq(iunit,rdata2,25,maxlevs,nlev2,'ROSEQ3')  ! RAOC HEIGHT/REFRACTIVITY SEQUENCE

      if ( nlev1 /= nlev2 ) then
         cycle reports
      end if

 !--------  determine FGAT index ifgat
 
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do  

      iprof = iprof + 1

      lev_loop: do k = 1, nlev1

         info%lat  = rdata1(1,k)
         info%lon  = rdata1(2,k)

         ! gpsro.bufr contains missing longitude, occasionally
         if ( info%lat > r8bfms .or. info%lon > r8bfms ) then
            cycle lev_loop
         end if

         height    = rdata2(1,k)
         ref_data  = rdata2(2,k)

         ! check for missing data
         if ( height > r8bfms .or. ref_data > r8bfms ) then
            cycle lev_loop
         end if

         ref_qc    = 0                ! initialized to be good
         ref_error = ref_data * 0.01

         ! check loc
         info%lat = max(info%lat, -89.95)
         info%lat = min(info%lat,  89.95)
         call da_llxy(info, loc, outside, outside_all)
         if ( outside_all ) then
            num_outside_all = num_outside_all + 1
            if ( print_detail_obs ) then
               write(unit=stderr,fmt='(a,2(1x,f8.3),a)')  &
                  id(1:5), info%lat, info%lon, '  -> outside_domain'
            end if
            cycle lev_loop
         end if
         ntotal = ntotal + 1
         ntotal_ifgat(ifgat)=ntotal_ifgat(ifgat)+1
         if ( outside ) then
            cycle lev_loop
         end if
         nlocal = nlocal + 1

         ! check height, only keep data below certain height (default is 30km)
         if ( height > top_km_gpsro*1000.0 .or. &
              height < bot_km_gpsro*1000.0 ) then
            ref_qc = -77
         end if
         
         err_opt = 1
         ! observation errors  WRF-Var/obsproc
         if ( err_opt == 1 ) then
            if ( height >= 12000.0 ) then
               ref_error = ref_error * 0.3
            else
               erh90 = (0.3-1.5)*(height-12000.0)/(12000.0-0.0) + 0.3
               if ( height >= 5500.0 ) then
                  erh0 = (0.3-1.3)*(height-12000.0)/(12000.0-5500.0) + 0.3
               else if ( height >= 2500.0) then
                  erh0 = (1.3-2.5)*(height-5500.0)/(5500.0-2500.0) + 1.3
               else
                  erh0 = 2.5
               end if
               err90 = ref_error * erh90
               err0  = ref_error * erh0
               ref_error = err90 - (1.0-abs(info%lat)/90.0)*(err90-err0)
            end if
         end if

         ! observation errors  GSI_Q1FY09,  Kuo et al. 2003
         if ( err_opt == 2 ) then
            if ( (info%lat >= -30.0) .and. (info%lat <= 30.0) ) then   ! tropics
               if ( (height >= 7000.0) .and. (height <= 31000.0) ) then
                  ref_error = ref_error*(0.1125+(1.25e-5*height))
               else if ( height > 31000.0 ) then
                  ref_error = ref_error*0.5
               else if ( height < 7000.0  ) then
                  ref_error = ref_error*(3.0-(4.0e-4*height))
               else
                  write(unit=message(1),fmt='(a,f8.1,a,f8.2)') 'unable to process with height = ', &
                     height, ' at lat = ', info%lat
                  call da_error("da_read_obs_bufrgpsro.inc",267,message(1:1))
               end if
            else   ! mid-latitudes
               if ( (height >= 5000.0) .and. (height <= 25000.0) ) then
                  ref_error = ref_error*0.3
               else if ( (height >= 25000.0) .and. (height <= 31000.0) ) then
                  ref_error = ref_error*(-3.45+(1.5e-4*height))
               else if ( height > 31000.0 ) then
                  ref_error = ref_error*1.2
               else if ( height < 5000.0 ) then
                  ref_error = ref_error*(0.75-(9.0e-5*height))
               else
                  write(unit=message(1),fmt='(a,f8.1,a,f8.2)') 'unable to process with height = ', &
                     height, ' at lat = ', info%lat
                  call da_error("da_read_obs_bufrgpsro.inc",281,message(1:1))
               end if
            end if
         end if

         !write(info%name, '(a,i6.6,a,a)') 'NCEP_GPSRO_', nlocal, '_', date_char
         write(info%name, '(a,i6.6,a,a)') 'NCEP_GPSRO_', iprof, '_', date_char

         if ( print_detail_obs ) then
            write(unit=stdout,fmt='(a,1x,a,1x,i4.4,4i2.2,2f8.2,f8.1,f8.2,i3,f9.5)')  &
               info%name,id(1:5),iyear,imonth,iday,ihour,imin, &
               info%lat,info%lon,height,ref_data,ref_qc,ref_error 
         end if
        
         if (.not. associated(head)) then
             nullify ( head )
             allocate ( head )
             nullify ( head%next )
             plink => head
         else
             allocate ( plink%next )
             plink => plink%next
             nullify ( plink%next )
         end if

         plink%info%name      = info%name
         plink%info%platform  = 'FM-116 GPSRF'
         plink%info%id        = id
         plink%info%date_char = date_char
         plink%info%levels    = 1              ! each level is treated as separate obs
         plink%info%elv       = 0.0            ! not used
         plink%info%lat     = info%lat
         plink%info%lon     = info%lon

         plink%loc%x       = loc%x
         plink%loc%y       = loc%y
         plink%loc%i       = loc%i
         plink%loc%j       = loc%j
         plink%loc%dx      = loc%dx
         plink%loc%dxm     = loc%dxm
         plink%loc%dy      = loc%dy
         plink%loc%dym     = loc%dym

         plink%slp%inv   = missing_r
         plink%slp%qc    = missing_data
         plink%slp%error = xmiss
         plink%pw%inv    = missing_r
         plink%pw%qc     = missing_data
         plink%pw%error  = xmiss

         plink%obs_global_index = ntotal

         nlev = 1
         allocate (plink%gpsref%h  (1:nlev))
         allocate (plink%gpsref%ref(1:nlev))
         allocate (plink%gpsref%p  (1:nlev))
         allocate (plink%gpsref%t  (1:nlev))
         allocate (plink%gpsref%q  (1:nlev))
         do i = 1, nlev
            plink%gpsref%h(i)         = height
            plink%gpsref%ref(i)%inv   = ref_data
            plink%gpsref%ref(i)%qc    = ref_qc
            plink%gpsref%ref(i)%error = ref_error
            plink%gpsref%p(i)%inv     = missing_r
            plink%gpsref%p(i)%qc      = missing_data
            plink%gpsref%p(i)%error   = xmiss
            plink%gpsref%t(i)%inv     = missing_r
            plink%gpsref%t(i)%qc      = missing_data
            plink%gpsref%t(i)%error   = xmiss
            plink%gpsref%q(i)%inv     = missing_r
            plink%gpsref%q(i)%qc      = missing_data
            plink%gpsref%q(i)%error   = xmiss
         end do
         plink%ifgat=ifgat
         num_p=num_p+1
         
         end do lev_loop
        end do reports  

        call closbf(iunit)
        close(iunit)

end do bufrfile
         
         iv%info(gpsref)%ntotal=ntotal
         iv%info(gpsref)%nlocal=nlocal
         iv%info(gpsref)%max_lev=1
         
!allocate iv 
          
       if (iv%info(gpsref)%nlocal > 0) allocate(iv%gpsref(1:iv%info(gpsref)%nlocal))
       
       if (iv%info(gpsref)%nlocal > 0) then
         allocate (iv%info(gpsref)%name(iv%info(gpsref)%nlocal))     
         allocate (iv%info(gpsref)%platform(iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%id(iv%info(gpsref)%nlocal))       
         allocate (iv%info(gpsref)%date_char(iv%info(gpsref)%nlocal))
         allocate (iv%info(gpsref)%levels(iv%info(gpsref)%nlocal))   
         allocate (iv%info(gpsref)%lat(iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal))    
         allocate (iv%info(gpsref)%lon(iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal))    
         allocate (iv%info(gpsref)%elv(iv%info(gpsref)%nlocal))      
         allocate (iv%info(gpsref)%pstar(iv%info(gpsref)%nlocal))    
         allocate (iv%info(gpsref)%slp(iv%info(gpsref)%nlocal))   
         allocate (iv%info(gpsref)%pw(iv%info(gpsref)%nlocal))    
         allocate (iv%info(gpsref)%x  (kms:kme,iv%info(gpsref)%nlocal))   
         allocate (iv%info(gpsref)%y  (kms:kme,iv%info(gpsref)%nlocal))   
         allocate (iv%info(gpsref)%i  (kms:kme,iv%info(gpsref)%nlocal))   
         allocate (iv%info(gpsref)%j  (kms:kme,iv%info(gpsref)%nlocal))      
         allocate (iv%info(gpsref)%dx (kms:kme,iv%info(gpsref)%nlocal))  
         allocate (iv%info(gpsref)%dxm(kms:kme,iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%dy (kms:kme,iv%info(gpsref)%nlocal))  
         allocate (iv%info(gpsref)%dym(kms:kme,iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%k  (iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal))
         allocate (iv%info(gpsref)%dz (iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal))  
         allocate (iv%info(gpsref)%dzm(iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%zk (iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%proc_domain(iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%thinned(iv%info(gpsref)%max_lev,iv%info(gpsref)%nlocal)) 
         allocate (iv%info(gpsref)%obs_global_index(iv%info(gpsref)%nlocal)) 
         iv%info(gpsref)%proc_domain(:,:)  = .false.
         iv%info(gpsref)%thinned(:,:)      = .false.
         iv%info(gpsref)%zk(:,:)           = missing_r
        end if
!
!sort for 4d var
         nlocal=0
         do kk=1,num_fgat_time
         plink => head 
         iv%info(gpsref)%ptotal(kk)=0
         
         reports2: do ii=1,num_p
         
         if (plink%ifgat /= kk) then  !sort iv
            plink => plink%next
            cycle reports2
         else
         nlocal=nlocal+1
         iv%info(gpsref)%name(nlocal)      = plink%info%name
         iv%info(gpsref)%platform(nlocal)  = plink%info%platform
         iv%info(gpsref)%id(nlocal)        = plink%info%id
         iv%info(gpsref)%date_char(nlocal) = plink%info%date_char
         iv%info(gpsref)%levels(nlocal)    = plink%info%levels              ! each level is treated as separate obs
         iv%info(gpsref)%elv(nlocal)       = plink%info%elv            ! not used
         iv%info(gpsref)%lat(:,nlocal)     = plink%info%lat
         iv%info(gpsref)%lon(:,nlocal)     = plink%info%lon

         iv%info(gpsref)%x(:,nlocal)       = plink%loc%x
         iv%info(gpsref)%y(:,nlocal)       = plink%loc%y
         iv%info(gpsref)%i(:,nlocal)       = plink%loc%i
         iv%info(gpsref)%j(:,nlocal)       = plink%loc%j
         iv%info(gpsref)%dx(:,nlocal)      = plink%loc%dx
         iv%info(gpsref)%dxm(:,nlocal)     = plink%loc%dxm
         iv%info(gpsref)%dy(:,nlocal)      = plink%loc%dy
         iv%info(gpsref)%dym(:,nlocal)     = plink%loc%dym

         iv%info(gpsref)%slp(nlocal)%inv   = plink%slp%inv
         iv%info(gpsref)%slp(nlocal)%qc    = plink%slp%qc
         iv%info(gpsref)%slp(nlocal)%error = plink%slp%error
         iv%info(gpsref)%pw(nlocal)%inv    = plink%pw%inv
         iv%info(gpsref)%pw(nlocal)%qc     = plink%pw%qc
         iv%info(gpsref)%pw(nlocal)%error  = plink%pw%error

         iv%info(gpsref)%obs_global_index(nlocal) = plink%obs_global_index

         nlev = 1
         allocate (iv%gpsref(nlocal)%h  (1:nlev))
         allocate (iv%gpsref(nlocal)%ref(1:nlev))
         allocate (iv%gpsref(nlocal)%p  (1:nlev))
         allocate (iv%gpsref(nlocal)%t  (1:nlev))
         allocate (iv%gpsref(nlocal)%q  (1:nlev))

         do i = 1, nlev
            iv%gpsref(nlocal)%h(i)         = plink%gpsref%h(i)
            iv%gpsref(nlocal)%ref(i)%inv   = plink%gpsref%ref(i)%inv
            iv%gpsref(nlocal)%ref(i)%qc    = plink%gpsref%ref(i)%qc
            iv%gpsref(nlocal)%ref(i)%error = plink%gpsref%ref(i)%error

            iv%gpsref(nlocal)%p(i)%inv     = plink%gpsref%p(i)%inv
            iv%gpsref(nlocal)%p(i)%qc      = plink%gpsref%p(i)%qc
            iv%gpsref(nlocal)%p(i)%error   = plink%gpsref%p(i)%error
            iv%gpsref(nlocal)%t(i)%inv     = plink%gpsref%t(i)%inv
            iv%gpsref(nlocal)%t(i)%qc      = plink%gpsref%t(i)%qc
            iv%gpsref(nlocal)%t(i)%error   = plink%gpsref%t(i)%error
            iv%gpsref(nlocal)%q(i)%inv     = plink%gpsref%q(i)%inv
            iv%gpsref(nlocal)%q(i)%qc      = plink%gpsref%q(i)%qc
            iv%gpsref(nlocal)%q(i)%error   = plink%gpsref%q(i)%error
         end do
       end if
       plink => plink%next
      end do reports2

      ntotal_ifgat(kk)=ntotal_ifgat(kk)+ntotal_ifgat(kk-1)
      iv%info(gpsref)%ptotal(kk)=ntotal_ifgat(kk)
      iv%info(gpsref)%plocal(kk)=nlocal
      ! thinning is not applied to GPSRO 1, simply make a copy from total number
      iv%info(gpsref)%thin_ptotal(kk)=ntotal_ifgat(kk)
      iv%info(gpsref)%thin_plocal(kk)=nlocal
   end do 

   write(unit=message(1),fmt='(A,3(1x,i7))') &
      'da_read_obs_bufrgpsro: num_report, num_outside_all, num_outside_time: ', &
      num_report, num_outside_all, num_outside_time
   call da_message(message(1:1))

   if ( nlocal /= iv%info(gpsref)%nlocal ) then
      call da_error("da_read_obs_bufrgpsro.inc",486,(/"numbers mismatch between scanning and reading NCEP GSPRO BUFR file"/))
   end if

   ! Release the linked list
   plink => head

   DO WHILE ( ASSOCIATED (plink) )
      head => plink%next
      deallocate (plink%gpsref%h)
      deallocate (plink%gpsref%ref)
      deallocate (plink%gpsref%p)
      deallocate (plink%gpsref%t)
      deallocate (plink%gpsref%q)
      deallocate ( plink )
      plink => head
   ENDDO

   NULLIFY (head)

   if (trace_use) call da_trace_exit("da_read_obs_bufrgpsro")

end subroutine da_read_obs_bufrgpsro
subroutine da_final_write_obs(it,iv)

   !-------------------------------------------------------------------------
   ! Purpose: Writes full diagnostics for O, (O-B) & OMA together
   !-------------------------------------------------------------------------   

   implicit none
 
   integer,        intent(in)    :: it
   type (iv_type), intent(in)    :: iv      ! O-B structure.
   integer                       :: n, k, iunit
   integer                       :: ios  ! Error code from MPI routines.
   integer                       :: num_obs
   logical                       :: if_wind_sd
   character(len=filename_len), allocatable     :: filename(:)
   character(len=filename_len)                  :: file


   if (trace_use) call da_trace_entry("da_final_write_obs")

   ! Wait to ensure all temporary files have been written
   call mpi_barrier(comm, ierr)

   if (rootproc) then
      call da_get_unit(iunit)
      allocate (filename(0:num_procs-1))
      do k = 0,num_procs-1
         write(unit=filename(k),fmt ='(a,i2.2,a,i4.4)')'gts_omb_oma_',it,'.',k
      end do 
      call da_get_unit(omb_unit)
       write(unit=file,fmt ='(a,i2.2)')'gts_omb_oma_',it
      open(unit=omb_unit,file=trim(file),form='formatted', status='replace', iostat=ios) 
      if (ios /= 0) call da_error("da_final_write_obs.inc",35, &
         (/"Cannot open file "//file/))
   end if

   num_obs = 0
   if (iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
         if(iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_synop) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'synop', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'synop',5,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [2] writing Metar
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
         if (iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_metar) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,20i8)')'metar', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'metar',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [3] writing Ships
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
         if(iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_ships) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'ships', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'ships',5,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [4] writing GeoAMV
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(geoamv)%nlocal > 0) then
       do n = 1, iv%info(geoamv)%nlocal
          if (iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_geoamv) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'geoamv', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'geoamv',6,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [5] writing PolarAMV
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
         if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_polaramv) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'polaramv', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'polaramv',8,if_wind_sd) 
      end do
   end if

   !------------------------------------------------------------------
   ! [5] writing GPSPW  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
         if(iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'gpspw', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'gpspw',5,if_wind_sd)     
      end do
   end if

   !------------------------------------------------------------------
   ! [6] writing Sonde  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_sound) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'sound', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'sound',5,if_wind_sd) 
      end do
   end if

!  Now sonde_sfc
   num_obs = 0
   if (iv%info(sonde_sfc)%nlocal > 0) then
      do n = 1, iv%info(sonde_sfc)%nlocal
         if(iv%info(sonde_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'sonde_sfc', num_obs
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'sonde_sfc',9,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [7] writing Airep  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         if(iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_airep) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'airep', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'airep',5,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [8] writing Pilot  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
         if(iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_pilot) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'pilot', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'pilot',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [9] writing ssmi_rv
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
         if(iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'ssmir', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'ssmir',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [10] writing SSMITB
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
         if (iv%info(ssmi_tb)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'ssmiT', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'ssmiT',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [11] writing SATEM  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         if(iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'satem', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'satem',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [12] writing SSMT1  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
         if(iv%info(ssmt1)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'ssmt1', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'ssmt1',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [13] writing SSMT2  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         if(iv%info(ssmt2)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'ssmt2', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'ssmt2',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [14] writing QSCAT  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
         if(iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_qscat) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'qscat', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'qscat',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [15] writing Profiler
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         if(iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_profiler) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'profiler', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'profiler',8,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [16] writing Buoy 
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
         if(iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_buoy) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'buoy', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'buoy',4,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [17] writing Bogus 
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         if(iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'bogus', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'bogus',5,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [18] writing Tamdar
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(tamdar)%nlocal > 0) then
      do n = 1, iv%info(tamdar)%nlocal
         if (iv%info(tamdar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (wind_sd_tamdar) if_wind_sd = .true.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'tamdar', num_obs
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'tamdar',6,if_wind_sd)
      end do

   end if

!  Now tamdar_sfc
   num_obs = 0
   if (iv%info(tamdar_sfc)%nlocal > 0) then
      do n = 1, iv%info(tamdar_sfc)%nlocal
         if(iv%info(tamdar_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'tamdar_sfc', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'tamdar_sfc',10,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [19] writing AIRS retrievals:
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         if(iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'airsr', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'airsr',5,if_wind_sd)
      end do
   end if

   !------------------------------------------------------------------
   ! [20] writing GPS refractivity
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
         if(iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'gpsref', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'gpsref',6,if_wind_sd)    
      end do
   end if

   !------------------------------------------------------------------
   ! [21] writing rain 
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(rain)%nlocal > 0) then
      do n = 1, iv%info(rain)%nlocal
         if(iv%info(rain)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if_wind_sd = .false.
   if (num_obs > 0 .and. rootproc) then
      write(omb_unit,'(a20,i8)')'rain', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_omb_tmp(filename(k),iunit,num_obs,'rain',4,if_wind_sd)    
      end do
   end if


   if (rootproc) then
      close(iunit)
      close(omb_unit)
      call da_free_unit(iunit)
      call da_free_unit(omb_unit)
      deallocate (filename)
   end if

   if (trace_use) call da_trace_exit("da_final_write_obs")
   
end subroutine da_final_write_obs


subroutine da_final_write_y(iv)

   !-------------------------------------------------------------------------   
   ! Purpose: Writes full diagnostics for y=H(x_inc)    
   !-------------------------------------------------------------------------   

   implicit none

   type (iv_type), intent(in)    :: iv   ! O-B structure.
      
   integer                       :: n, k,kk,i, ounit
   integer                       :: sound_num_obs, num_obs, ios
   character(len=filename_len), allocatable     :: filename(:)
   character(len=filename_len)                  :: ob_name, file_prefix  

   if (trace_use_dull) call da_trace_entry("da_final_write_y")

   ! Ensure other processors have written their temporary files
   call mpi_barrier(comm, ierr)

   if (omb_add_noise) then
      ! perturbed ob run.
      file_prefix='pert_obs'
   else
      ! unperturbed ob run.
      file_prefix='unpert_obs'
   end if

   if (rootproc) then
      allocate (filename(0:num_procs-1))
      do k = 0,num_procs-1
         write(unit=filename(k),fmt ='(a,a,i4.4)')trim(file_prefix),'.',k
      end do
      call da_get_unit(ounit)
      open(unit=ounit,file=trim(file_prefix),form='formatted', &
         status='replace' , iostat=ios)
      if (ios /= 0) call da_error("da_final_write_y.inc",39, &
         (/"Cannot open random observation error file"//file_prefix/))
   end if

   !------------------------------------------------------------------
   ! [1] writing Surface
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
         if (iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'synop', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'synop',5)
      end do
   end if

   !------------------------------------------------------------------
   ! [2] writing Metar
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
        if (iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,20i8)')'metar', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'metar',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [3] writing Ships
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
        if (iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'ships', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'ships',5)
      end do
   end if

   !---------------------------------------------------------------
   ! [4] writing GeoAMV
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
         if (iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'geoamv', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'geoamv',6)
      end do
   end if

   !------------------------------------------------------------------
   ! [5] writing PolarAMV
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
        if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'polaramv', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'polaramv',8) 
      end do
   end if

   !------------------------------------------------------------------
   ! [5] writing GPSPW  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
         if (iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'gpspw', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'gpspw',5)     
      end do
   end if

   !------------------------------------------------------------------
   ! [6] writing Sonde  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   sound_num_obs = num_obs
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'sound', sound_num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'sound',5)     
      end do  
      !  writing Sonde_sfc  
      write(ounit,'(a20,i8)')'sonde_sfc', sound_num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'sonde_sfc',9)
      end do
   end if

   !------------------------------------------------------------------
   ! [7] writing Airep  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         if (iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'airep', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'airep',5)   
      end do
   end if

   !------------------------------------------------------------------
   ! [8] writing Pilot  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
        if (iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'pilot', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'pilot',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [9] writing ssmi_rv
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
         if (iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmir', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'ssmir',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [10] writing SSMITB
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
         if (iv%info(ssmi_tb)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmiT', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'ssmiT',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [11] writing SATEM  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         if (iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'satem', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'satem',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [12] writing SSMT1  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
        if (iv%info(ssmt1)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmt1', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'ssmt1',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [13] writing SSMT2  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         if (iv%info(ssmt2)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'ssmt2', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'ssmt2',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [14] writing QSCAT  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
         if (iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'qscat', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'qscat',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [15] writing Profiler
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         if (iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'profiler', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'profiler',8)    
      end do
   end if

   !---------------------------------------------------------------
   ! [16] writing Buoy 
   !---------------------------------------------------------------

   num_obs = 0
   if (iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
         if (iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'buoy', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'buoy',4)    
      end do
   end if

   !---------------------------------------------------------------
   ! [17] writing  Bogus  
   !---------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         if (iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'bogus', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'bogus',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [18] writing AIRS retrievals:
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         if (iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'airsr', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'airsr',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [19] writing Radiance data:
   !------------------------------------------------------------------

   if (iv%num_inst > 0) then
      do i = 1, iv%num_inst                 ! loop for sensor
         do k = 1,iv%instid(i)%nchan        ! loop for channel
            !  Counting number of obs for channel k
            num_obs = 0
            do n = 1,iv%instid(i)%num_rad      ! loop for pixel
               if (iv%instid(i)%info%proc_domain(1,n) .and. &
                 (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer)) then
                  num_obs = num_obs + 1
               end if
            end do                                ! end loop for pixel
            call da_proc_sum_int(num_obs)
            if (rootproc .and. num_obs > 0) then
               write(ob_name,'(a,a,i4.4)') &
                  trim(iv%instid(i)%rttovid_string),'-', &
		  iv%instid(i)%ichan(k)
               write(ounit,'(a20,i8)')  ob_name,num_obs
               num_obs = 0
               do kk = 0,num_procs-1
                  call da_read_y_unit(trim(filename(kk)),ounit,num_obs, &
                     trim(ob_name),len(trim(ob_name)))    
               end do
            end if
         end do                                  ! end loop for Channel 
      end do                                   ! end loop for sensor
   end if
   !------------------------------------------------------------------
   ! [20] writing gpsref:
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
         if (iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'gpsref', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'gpsref',6)    
      end do
   end if

   !------------------------------------------------------------------
   ! [21] writing RAIN  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(rain)%nlocal > 0) then      
      do n = 1, iv%info(rain)%nlocal       
         if (iv%info(rain)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if   
   call da_proc_sum_int(num_obs)
   if (rootproc .and. num_obs > 0) then
      write(ounit,'(a20,i8)')'rain', num_obs
      num_obs = 0
      do k = 0,num_procs-1                        
         call da_read_y_unit(trim(filename(k)),ounit,num_obs,'rain',4)
      end do
   end if      


   !------------------------------------------------------------------

   if (rootproc) then
      close(ounit)
      call da_free_unit(ounit)
      deallocate (filename)
   end if

   if (trace_use_dull) call da_trace_exit("da_final_write_y")
   
end subroutine da_final_write_y


subroutine da_read_y_unit(filename,unit_out,num,obs_type_in, nc)

   !-------------------------------------------------------------------------  
   ! Purpose: read diagnostics written on yp_unit or y_unit by WRF-Var
   !-------------------------------------------------------------------------   

   implicit none

   integer      ,intent (in)    :: unit_out
   integer      ,intent (inout) :: num
   character*(*),intent (in)    :: obs_type_in, filename                 
   integer      ,intent (in)    :: nc      

   integer      :: num_obs , unit_in, ios
   character*20 :: iv_type               
   logical      :: if_write
   
   real         :: fld(7), fld_rad                          
   integer      :: n, k, n1,n2, levels

   if (trace_use_dull) call da_trace_entry("da_read_y_unit")

   iv_type="Unknown"

   call da_get_unit(unit_in)
   open(unit=unit_in,file=trim(filename),form='formatted',iostat=ios,status='old')
   if (ios /= 0) then
      call da_error("da_read_y_unit.inc",28, &
         (/"Cannot open random observation error file"//filename/))
   end if

   reports: do
   read(unit_in,'(a20,i8)', end = 999, err = 1000) iv_type,num_obs
   
   if_write = .false.
   if (index(iv_type,OBS_type_in(1:nc)) > 0) if_write = .true.

   ! If radiance data treat differently
   if ( (index(iv_type,'noaa') > 0) .or. (index(iv_type,'eos') > 0) .or.   &
        (index(iv_type,'dmsp') > 0) .or. (index(iv_type,'metop') > 0) .or. &
        (index(iv_type,'tiros') > 0) .or. (index(iv_type,'msg') > 0) .or. &
        (index(iv_type,'jpss') > 0)  .or. (index(iv_type,'gcom-w') >0) ) then

      do n = 1, num_obs
         if (if_write) num = num + 1
         read(unit_in,'(2i8,e15.7)')n1, n2, fld_rad
         if (if_write)write(unit_out,'(2i8,e15.7)')num,n2, fld_rad
      end do
   else
      do n = 1, num_obs
         if (if_write) num = num + 1
         if (index(iv_type,'bogus') > 0) then
            read(unit_in,'(i8)', err=1000)levels
            if (if_write) write(unit_out,'(i8)')levels
            read(unit_in,'(2i8,7e15.7)', err= 1000) n1, n2, fld
            if (if_write) write(unit_out,'(2i8,7e15.7)') num, levels, fld  
         end if
         read(unit_in,'(i8)', err=1000)levels
         if (if_write) write(unit_out,'(i8)')levels
         do k = 1, levels
            read(unit_in,'(2i8,7e15.7)', err= 1000) n1, n2, fld  
            if (if_write) write(unit_out,'(2i8,7e15.7)') num, k, fld
         end do
      end do
   end if
   if (if_write) exit reports
   cycle reports
1000  continue 
   write(unit=message(1), fmt='(a,i3,a,a)') &
      'read error on unit: ',unit_in, ' for iv_type', trim(iv_type)
   ! call da_warning("da_read_y_unit.inc",71,message(1:1))
   end do reports
 999  continue 
   close (unit_in)
   call da_free_unit(unit_in)

   if (trace_use_dull) call da_trace_exit("da_read_y_unit")

end subroutine da_read_y_unit


subroutine da_read_rand_unit(filename, unit_in,num,obs_type_in, nc)

   !----------------------------------------------------------------------------   
   ! Purpose: Program to read diagnostics written on rand_unit by WRF-Var
   !----------------------------------------------------------------------------   

   implicit none

   integer      ,intent (in)    :: unit_in
   integer      ,intent (inout) :: num
   character*(*),intent (in)    :: obs_type_in 
   character*(*),intent (inout) :: filename                    
   integer      ,intent (in)    :: nc      

   integer      :: num_obs 
   character*20 :: iv_type               
   logical   :: if_write
   
   real         :: fld(10), fld1_rad , fld2_rad                         
   integer      :: n, k, n1,n2, levels, ios

   if (trace_use) call da_trace_entry("da_read_rand_unit")

   open(unit=unit_in,   file=trim(filename), status='old',iostat=ios)
   if (ios /= 0) then
      call da_error("da_read_rand_unit.inc",26, &
         (/"Cannot open file"//filename/))
   end if
1  continue
  

   read(unit_in,'(a20,i8)', end = 999, err = 1000)iv_type,num_obs                    
   
   if_write = .false.
   if (index(iv_type,OBS_type_in(1:nc)) > 0) if_write = .true.

   ! If radiance data treat differently
   if ( (index(iv_type,'noaa') > 0)  .or. (index(iv_type,'eos') > 0) .or.  &
        (index(iv_type,'dmsp') > 0) .or. (index(iv_type,'metop') > 0) ) then
      do n = 1, num_obs
         if (if_write) num = num + 1
         read(unit_in,'(2i8,f10.3,e15.7)')n1, n2, fld1_rad,fld2_rad
         if (if_write)write(rand_unit,'(2i8,f10.3,e15.7)')num,n2, fld1_rad,fld2_rad
      end do
   else
      do n = 1, num_obs        
         if (if_write) num = num + 1
         if (index(iv_type,'bogus')     > 0)  then
            read(unit_in,'(i8)', err=1000)levels
            if (if_write) write(rand_unit,'(i8)')levels
            read(unit_in,'(2i8,10e15.7)', err= 1000) n1, n2, fld
            if (if_write) write(rand_unit,'(2i8,10e15.7)') num, levels, fld  
         end if
         read(unit_in,'(i8)', err=1000)levels
         if (if_write) write(rand_unit,'(i8)')levels
         do k = 1, levels
            read(unit_in,'(2i8,10e15.7)', err= 1000) n1, n2, fld  
            if (if_write) write(rand_unit,'(2i8,10e15.7)') num, k, fld
         end do
      end do
   end if

   if (if_write)  goto 999
   goto 1

999 continue
   close (unit_in)
   if (trace_use) call da_trace_exit("da_read_rand_unit")
   return

1000  continue 
   write(unit=message(1), fmt='(/a,i3/a/)') &
              'read error on unit: ',unit_in, &
              'in da_read_rand_unit'
   call da_warning("da_read_rand_unit.inc",75,message(1:1))

   if (trace_use) call da_trace_exit("da_read_rand_unit")

end subroutine da_read_rand_unit


subroutine da_read_omb_tmp(filename,unit_in,num,obs_type_in,nc,if_wind_sd)

   !-------------------------------------------------------------------------
   ! read diagnostics written to temporary file by WRFVAR
   !
   ! 07 MAR 2014 -- Variables of OMB/OMA for wind obs. in diagnostic files 
   !                are optional, i.e. SPD/DIR or U/V        -- Feng Gao
   !-------------------------------------------------------------------------

   implicit none

   integer      ,intent (in)    :: unit_in
   integer      ,intent (inout) :: num
   character*(*),intent (in)    :: obs_type_in, filename 
   integer      ,intent (in)    :: nc

   integer      :: num_obs, ios 
   character*20 :: iv_type
   logical      :: if_write, if_wind_sd
   
   character*5  :: stn_id
   integer      :: n, k, kk, l, levels, dummy_i
   real         :: lat, lon, press, height, dummy
   real         :: tpw_obs, tpw_inv, tpw_err, tpw_inc
   real         :: u_obs, u_inv, u_error, u_inc, & 
                   v_obs, v_inv, v_error, v_inc, &
                   t_obs, t_inv, t_error, t_inc, &
                   p_obs, p_inv, p_error, p_inc, &
                   q_obs, q_inv, q_error, q_inc, &
                   spd_obs, spd_inv, spd_err, spd_inc,   &
                   ref_obs, ref_inv, ref_error, ref_inc, &
		   rain_obs, rain_inv, rain_error, rain_inc, zk
   integer     :: u_qc, v_qc, t_qc, p_qc, q_qc, tpw_qc, spd_qc, ref_qc, rain_qc

   if (trace_use_dull) call da_trace_entry("da_read_omb_tmp")

   open(unit=unit_in,file=trim(filename),form='formatted',status='old',iostat=ios)
   if (ios /= 0) then
      call da_error("da_read_omb_tmp.inc",39, (/"Cannot open file"//filename/))
   end if

   reports: do

      read(unit_in,'(a20,i8)', end = 999, err = 1000) iv_type,num_obs
      if_write = .false.
      if (index(iv_type,OBS_type_in(1:nc)) > 0) if_write = .true.

      select case (trim(adjustl(iv_type)))

      case ('synop', 'ships', 'buoy', 'metar', 'sonde_sfc', 'tamdar_sfc')
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)')levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
               num = num + 1
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, press, &       ! Lat/lon, pressure
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc, &
                     t_obs, t_inv, t_qc, t_error, t_inc, &
                     p_obs, p_inv, p_qc, p_error, p_inc, &
                     q_obs, q_inv, q_qc, q_error, q_inc

                  if (.not. if_wind_sd .and. wind_stats_sd) & 
                     call da_ffdduv_diagnose(u_obs, u_inv, u_inc, v_obs, v_inv, v_inc, u_qc, v_qc, convert_uv2fd)
                  if (if_wind_sd .and. .not. wind_stats_sd) &
                     call da_ffdduv_diagnose(u_obs, u_inv, u_inc, v_obs, v_inv, v_inc, u_qc, v_qc, convert_fd2uv)

                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num, k, stn_id, &          ! Station
                        lat, lon, press, &       ! Lat/lon, pressure
                        u_obs, u_inv, u_qc, u_error, u_inc, & 
                        v_obs, v_inv, v_qc, v_error, v_inc, &
                        t_obs, t_inv, t_qc, t_error, t_inc, &
                        p_obs, p_inv, p_qc, p_error, p_inc, &
                        q_obs, q_inv, q_qc, q_error, q_inc
               end do
            end do
         end if

         if (if_write) exit reports
         cycle reports

      case ('pilot', 'profiler', 'geoamv', 'qscat', 'polaramv')
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)')levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
               num = num + 1
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                      kk, l, stn_id, &          ! Station
                      lat, lon, press, &        ! Lat/lon, pressure
                      u_obs, u_inv, u_qc, u_error, u_inc, & 
                      v_obs, v_inv, v_qc, v_error, v_inc

                  if (.not. if_wind_sd .and. wind_stats_sd) &
                     call da_ffdduv_diagnose(u_obs, u_inv, u_inc, v_obs, v_inv, v_inc, u_qc, v_qc, convert_uv2fd)
                  if (if_wind_sd .and. .not. wind_stats_sd) &
                     call da_ffdduv_diagnose(u_obs, u_inv, u_inc, v_obs, v_inv, v_inc, u_qc, v_qc, convert_fd2uv)

                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num, k, stn_id, &          ! Station
                        lat, lon, press, &         ! Lat/lon, pressure
                        u_obs, u_inv, u_qc, u_error, u_inc, & 
                        v_obs, v_inv, v_qc, v_error, v_inc

               end do 
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('gpspw' )
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)')levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
               num = num + 1
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, dummy, &       ! Lat/lon, dummy    
                     tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num, k, stn_id,  &       ! Station
                        lat, lon, dummy, &       ! Lat/lon, dummy    
                        tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('sound', 'tamdar', 'airep')
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)')levels
               if (if_write) then
                   write(omb_unit,'(i8)')levels
                   num = num + 1 
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, press, &       ! Lat/lon, dummy    
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc, &
                     t_obs, t_inv, t_qc, t_error, t_inc, &
                     q_obs, q_inv, q_qc, q_error, q_inc

                  if (.not. if_wind_sd .and. wind_stats_sd) &
                     call da_ffdduv_diagnose(u_obs, u_inv, u_inc, v_obs, v_inv, v_inc, u_qc, v_qc, convert_uv2fd)  
                  if (if_wind_sd .and. .not. wind_stats_sd) &
                     call da_ffdduv_diagnose(u_obs, u_inv, u_inc, v_obs, v_inv, v_inc, u_qc, v_qc, convert_fd2uv) 

                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num, k, stn_id,  &       ! Station
                        lat, lon, press, &       ! Lat/lon, dummy    
                        u_obs, u_inv, u_qc, u_error, u_inc, & 
                        v_obs, v_inv, v_qc, v_error, v_inc, &
                        t_obs, t_inv, t_qc, t_error, t_inc, &
                        q_obs, q_inv, q_qc, q_error, q_inc
               end do 
            end do
         end if
     if (if_write) exit reports
     cycle reports

      case ('ssmir' )
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)')levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
                  num = num + 1 
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, dummy, &       ! Lat/lon, dummy    
                     spd_obs, spd_inv, spd_qc, spd_err, spd_inc, &
                     tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num, k, stn_id,  &       ! Station
                        lat, lon, dummy, &       ! Lat/lon, dummy    
                        spd_obs, spd_inv, spd_qc, spd_err, spd_inc, &
                        tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports
   
      case ('ssmit' )
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)')levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
                  num = num + 1 
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,7(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, dummy, &       ! Lat/lon, dummy    
                     dummy, dummy, dummy_i, dummy, dummy, &    
                     dummy, dummy, dummy_i, dummy, dummy, &    
                     dummy, dummy, dummy_i, dummy, dummy, &    
                     dummy, dummy, dummy_i, dummy, dummy, &    
                     dummy, dummy, dummy_i, dummy, dummy, &    
                     dummy, dummy, dummy_i, dummy, dummy, &    
                     dummy, dummy, dummy_i, dummy, dummy
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,7(2f17.7,i8,2f17.7))', err= 1000)&
                        num,k,stn_id, &          ! Station
                        lat, lon, dummy, &       ! Lat/lon, dummy    
                        dummy, dummy, dummy_i, dummy, dummy, &    
                        dummy, dummy, dummy_i, dummy, dummy, &    
                        dummy, dummy, dummy_i, dummy, dummy, &    
                        dummy, dummy, dummy_i, dummy, dummy, &    
                        dummy, dummy, dummy_i, dummy, dummy, &    
                        dummy, dummy, dummy_i, dummy, dummy, &    
                        dummy, dummy, dummy_i, dummy, dummy
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('satem' )
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)') levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
                  num = num + 1 
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, press, &       ! Lat/lon, dummy    
                     tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num,k,stn_id, &          ! Station
                        lat, lon, press, &       ! Lat/lon, dummy    
                        tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
               end do  
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('ssmt1' , 'ssmt2' )
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)') levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
                  num = num + 1 
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, dummy, &       ! Lat/lon, dummy    
                     dummy,dummy, dummy_i, dummy, dummy
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num,k,stn_id, &          ! Station
                        lat, lon, dummy, &       ! Lat/lon, dummy    
                        dummy,dummy, dummy_i, dummy, dummy
               end do 
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('bogus' )          
         ! TC Bogus data is written in two records
         ! 1st record holds info about surface level
         ! 2nd is for upper air

         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)') levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
                  num = num + 1 
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                      kk,l, stn_id, &          ! Station
                      lat, lon, press, &       ! Lat/lon, dummy    
                      u_obs, u_inv, u_qc, u_error, u_inc, & 
                      v_obs, v_inv, v_qc, v_error, v_inc
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                         num,l,stn_id, &          ! Station
                         lat, lon, press, &       ! Lat/lon, dummy    
                         u_obs, u_inv, u_qc, u_error, u_inc, & 
                         v_obs, v_inv, v_qc, v_error, v_inc
               end do
               read(unit_in,'(i8)') levels
               if (if_write) then
                  write(omb_unit,'(i8)')levels
               end if
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, press, &       ! Lat/lon, dummy    
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                         num,l,stn_id, &          ! Station
                         lat, lon, press, &       ! Lat/lon, dummy    
                         u_obs, u_inv, u_qc, u_error, u_inc, & 
                         v_obs, v_inv, v_qc, v_error, v_inc
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('airsr' )          
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)') levels
               if (if_write) write(omb_unit,'(i8)')levels
               num = num + 1
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, press, &       ! Lat/lon, dummy    
                     t_obs, t_inv, t_qc, t_error, t_inc, & 
                     q_obs, q_inv, q_qc, q_error, q_inc
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                         num,k,stn_id, &          ! Station
                         lat, lon, press, &       ! Lat/lon, dummy    
                         t_obs, t_inv, t_qc, t_error, t_inc, & 
                         q_obs, q_inv, q_qc, q_error, q_inc
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('gpsref' )          
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)') levels
               if (if_write) write(omb_unit,'(i8)')levels
               num = num + 1
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, height, &       ! Lat/lon, height   
                     ref_obs, ref_inv, ref_qc, ref_error, ref_inc 
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num,k,stn_id, &          ! Station
                        lat, lon, height, &       ! Lat/lon, height   
                        ref_obs, ref_inv, ref_qc, ref_error, ref_inc 
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case ('rain' )          
         if (num_obs > 0) then
            do n = 1, num_obs    
               read(unit_in,'(i8)') levels
               if (if_write) write(omb_unit,'(i8)')levels
               num = num + 1
               do k = 1, levels
                  read(unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          ! Station
                     lat, lon, height, &       ! Lat/lon, height   
                     rain_obs, rain_inv, rain_qc, rain_error, rain_inc 
                  if (if_write) &
                     write(omb_unit,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                        num,k,stn_id, &          ! Station
                        lat, lon, height, &       ! Lat/lon, height   
                        rain_obs, rain_inv, rain_qc, rain_error, rain_inc 
               end do
            end do
         end if
         if (if_write) exit reports
         cycle reports

      case default;
                  

         write(unit=message(1), fmt='(a,a20,a,i3)') &
            'Got unknown obs_type string:', trim(iv_type),' on unit ',unit_in
         call da_error("da_read_omb_tmp.inc",412,message(1:1))
      end select
   end do reports 

999 continue
   close (unit_in)

   if (trace_use_dull) call da_trace_exit("da_read_omb_tmp")
   return

1000 continue
   write(unit=message(1), fmt='(a,i3)') &
      'read error on unit: ',unit_in
   call da_warning("da_read_omb_tmp.inc",425,message(1:1))

end subroutine da_read_omb_tmp

subroutine da_write_noise_to_ob (iv) 

   !-------------------------------------------------------------------------
   ! Purpose: Write consolidated obs noise created in da_add_noise_to_ob
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv   ! Obs and header structure.
   integer                       :: n, k,kk,i, iunit
   integer                       :: num_obs
   character(len=filename_len), allocatable     :: filename(:)   
   character*20                  :: ob_name   

   call da_trace_entry("da_write_noise_to_ob")

   ! Ensure other processors have written their temporary files
   call mpi_barrier(comm, ierr)
   call da_get_unit (iunit)
   allocate (filename(0:num_procs-1))
   do k = 0,num_procs-1
      write(unit=filename(k), fmt='(a,i4.4)')'rand_obs_error.',k
   end do
   if (rootproc) then
      call da_get_unit (rand_unit)
      open(unit=rand_unit,file='rand_obs_error',form='formatted', &
         iostat=ierr,status='unknown')
!         iostat=ierr,status='new')
      if (ierr /= 0) &
         call da_error("da_write_noise_to_ob.inc",32, (/"Cannot open file rand_obs_error"/))
   end if

   num_obs = 0
   if (iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
         if (iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'synop', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'synop',5)
      end do
   end if

   !------------------------------------------------------------------
   ! [2] writing Metar
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
         if (iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,20i8)')'metar', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'metar',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [3] writing Ships
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
         if (iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'ships', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'ships',5)
      end do
   end if

   !------------------------------------------------------------------
   ! [4] writing GeoAMV
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
         if (iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'geoamv', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'geoamv',6)
      end do
   end if

   !------------------------------------------------------------------
   ! [5] writing PolarAMV
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
         if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      if (rootproc) write(rand_unit,'(a20,i8)')'polaramv', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'polaramv',8) 
      end do
   end if

   !------------------------------------------------------------------
   ! [5] writing GPSPW  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
         if (iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'gpspw', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'gpspw',5)     
      end do
   end if

   !------------------------------------------------------------------
   ! [6] writing Sonde  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         if (iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'sound', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'sound',5)     
      end do
      !  writing Sonde_sfc  
      if (rootproc) write(rand_unit,'(a20,i8)')'sonde_sfc', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'sonde_sfc',9)
      end do
   end if

   !------------------------------------------------------------------
   ! [7] writing Airep  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         if (iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'airep', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'airep',5)   
      end do
   end if

   !------------------------------------------------------------------
   ! [8] writing   
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
         if (iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'pilot', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'pilot',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [9] writing ssmi_rv
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
         if (iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'ssmir', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'ssmir',5)    
     end do
   end if

   !------------------------------------------------------------------
   ! [10] writing SSMITB
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
         if (iv%info(ssmi_tb)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'ssmiT', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'ssmiT',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [11] writing SATEM  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         if (iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'satem', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'satem',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [12] writing SSMT1  
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
         if (iv%info(ssmt1)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'ssmt1', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'ssmt1',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [13] writing SSMT2  
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         if (iv%info(ssmt2)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'ssmt2', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'ssmt2',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [14] writing QSCAT  
   !------------------------------------------------------------------
    
   num_obs = 0
   if (iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
         if (iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'qscat', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'qscat',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [15] writing Profiler
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         if (iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'profiler', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'profiler',8)    
      end do
   end if

   !------------------------------------------------------------------
   ! [16] writing Buoy 
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
         if (iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'buoy', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'buoy',4)    
      end do
   end if

   !------------------------------------------------------------------
   ! [17] writing  Bogus 
   !------------------------------------------------------------------
  
   num_obs = 0
   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         if (iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'bogus', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'bogus',5)    
      end do
   end if

   !------------------------------------------------------------------
   ! [18] writing AIRS retrievals
   !------------------------------------------------------------------
   
   num_obs = 0
   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         if (iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
       write(rand_unit,'(a20,i8)')'airsr', num_obs  
       num_obs = 0
       do k = 0,num_procs-1
          call da_read_rand_unit(filename(k),iunit,num_obs,'airsr',5)    
       end do
    end if

   !------------------------------------------------------------------
   ! [19] writing gpsref
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
         if (iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,20i8)')'gpsref', num_obs  
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'gpsref',6)    
      end do
   end if


   !------------------------------------------------------------------
   ! [20] writing Radiance data:  
   !------------------------------------------------------------------

   if (iv%num_inst > 0) then
      do i = 1, iv%num_inst                 ! loop for sensor
         !if (iv%instid(i)%num_rad < 1) cycle ! may broken da_proc_sum_int(num_obs) if some PE num_rad=0
         do k = 1,iv%instid(i)%nchan        ! loop for channel
            ! Counting number of obs for channle k
            num_obs = 0
            do n = 1,iv%instid(i)%num_rad      ! loop for pixel
               if (iv%instid(i)%info%proc_domain(1,n) .and. &
                  (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer)) then
                  num_obs = num_obs + 1
               end if
            end do                                ! end loop for pixel
            call da_proc_sum_int(num_obs)
            if (rootproc .and. num_obs > 0) then
               write (ob_name,'(a,a,i4.4)') &
                  trim(iv%instid(i)%rttovid_string),'-',iv%instid(i)%ichan(k)
               write (rand_unit,'(a20,i8)')  ob_name,num_obs
               num_obs = 0
               do kk = 0,num_procs-1
                  call da_read_rand_unit(filename(kk),iunit,num_obs, &
                     trim(ob_name),len(trim(ob_name)))
               end do
            end if
         end do                           ! end loop for channel
      end do                            ! end loop for sensor
   end if

   !------------------------------------------------------------------
   ! [21] writing RAIN 
   !------------------------------------------------------------------

   num_obs = 0
   if (iv%info(rain)%nlocal > 0) then
      do n = 1, iv%info(rain)%nlocal
         if (iv%info(rain)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
   end if
   call da_proc_sum_int(num_obs)
   if (num_obs > 0 .and. rootproc) then
      write(rand_unit,'(a20,i8)')'rain', num_obs
      num_obs = 0
      do k = 0,num_procs-1
         call da_read_rand_unit(filename(k),iunit,num_obs,'rain',4)
      end do
   end if


   call da_free_unit (iunit)
!rizvi   call da_free_unit (rand_unit)
   if (rootproc ) call da_free_unit (rand_unit)
   deallocate (filename)
   call da_trace_exit("da_write_noise_to_ob")
   
end subroutine da_write_noise_to_ob  


subroutine da_final_write_filtered_obs

   !---------------------------------------------------------------------------
   !  Purpose: Scans intermediate Filtered Obs files  
   !           and writes the data part on filtered_obs_unit
   !---------------------------------------------------------------------------
   implicit none        

   integer                    :: i,iost,iunit, files     
   type (multi_level_type)    :: platform
   real                       :: height_error
   character(len=filename_len) :: filename

   if (trace_use) call da_trace_entry("da_final_write_filtered_obs")

   call da_get_unit(iunit)

   ! Loop over all data files
   !--------------------------

   do files = 0, num_procs-1

      write(unit=filename, fmt='(a,i4.4)') 'filtered_obs.',files      
      open(unit=iunit, file= trim(filename),form='formatted', &
       status='old', iostat=iost)
      if(iost /= 0) &
      call da_error("da_final_write_filtered_obs.inc",27, (/"Cannot open "//filename/))

      !  loop over records
      !  -----------------
      reports: do
         !     read station general info
         !     =============================
         read (iunit, fmt = fmt_info, iostat = iost) &
                      platform%info%platform,    &
                      platform%info%date_char,   &
                      platform%info%name,        &
                      platform%info%levels,      &
                      platform%info%lat,         &
                      platform%info%lon,         &
                      platform%info%elv,         &
                      platform%info%id
         if (iost /= 0) then
            !write (0,'(/A,I9)') ' end OF OBS unit: ',iunit
            !write (0,'(/A,I9)') ' iostat:          ',iost
            exit reports
         end if
         write(filtered_obs_unit, fmt = fmt_info)           &
                      platform%info%platform,    &
                      platform%info%date_char,   &
                      platform%info%name,        &
                      platform%info%levels,      &
                      platform%info%lat,         &
                      platform%info%lon,         &
                      platform%info%elv,         &
                      platform%info%id

         !  Read surface Info
         !  -------------------
         read (iunit, fmt = fmt_srfc)  &
            platform%loc%slp%inv, platform%loc%slp%qc, platform%loc%slp%error, &
            platform%loc%pw%inv, platform%loc%pw%qc, platform%loc%pw%error
         write(filtered_obs_unit, fmt = fmt_srfc)  &
            platform%loc%slp%inv, platform%loc%slp%qc, platform%loc%slp%error, &
            platform%loc%pw%inv, platform%loc%pw%qc, platform%loc%pw%error

         ! == levels < 1 and not GPSPW, go back to reports

         if ((platform%info%levels < 1) .AND.            &
             (index(platform%info%platform, 'GPSPW') <= 0)) then
              cycle reports
         end if

         !     read EACH LEVELS
         !     ----------------
         do i = 1, platform%info%levels

            platform%each (i) = each_level_type(missing_r, missing, -1.0, & ! height
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! u
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! v
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! p
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! t
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! q
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! rh
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! td
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r))  ! speed 

            read (unit = iunit, fmt = trim (fmt_each)) &
               platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
               platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
               platform%each(i)%speed%error, &
               platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
               platform%each(i)%height,       &
               platform%each(i)%height_qc,    &
               height_error,                   &
               platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
               platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
               platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error

            write (unit = filtered_obs_unit, fmt = trim (fmt_each)) &
               platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
               platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
               platform%each(i)%speed%error, &
               platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
               platform%each(i)%height,       &
               platform%each(i)%height_qc,    &
               height_error,                   &
               platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
               platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
               platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error
         end do
      end do reports                  !  Loop over reports              
      close (iunit)
   end do !  Loop over all data files

   call da_free_unit (iunit)

   if (trace_use) call da_trace_exit("da_final_write_filtered_obs")

end subroutine da_final_write_filtered_obs 


subroutine da_final_write_modified_filtered_obs

   !---------------------------------------------------------------------------
   !  Purpose: Scans intermediate Filtered Obs files  
   !           and writes the data part on filtered_obs_unit
   !---------------------------------------------------------------------------
   implicit none        

   integer                    :: i,iost,iunit, files     
   type (multi_level_type)    :: platform
   real                       :: height_error
   character(len=filename_len) :: filename

   character(len=120) :: fmt_each_inv = &
      '(3(f12.3,i4,f7.2),11x,3(f12.3,i4,f7.2),11x,1(f12.3,i4,f7.2),4f17.7)'
   real                       :: uinv, vinv, tinv, qinv

   if (trace_use) call da_trace_entry("da_final_write_modified_filtered_obs")

   call da_get_unit(iunit)

   ! Loop over all data files
   !--------------------------

   do files = 0, num_procs-1

      write(unit=filename, fmt='(a,i4.4)') 'filtered_obs.',files      
      open(unit=iunit, file= trim(filename),form='formatted', &
       status='old', iostat=iost)
      if(iost /= 0) &
      call da_error("da_final_write_modified_filtered_obs.inc",31, (/"Cannot open "//filename/))

      !  loop over records
      !  -----------------
      reports: do
         !     read station general info
         !     =============================
         read (iunit, fmt = fmt_info, iostat = iost) &
                      platform%info%platform,    &
                      platform%info%date_char,   &
                      platform%info%name,        &
                      platform%info%levels,      &
                      platform%info%lat,         &
                      platform%info%lon,         &
                      platform%info%elv,         &
                      platform%info%id
         if (iost /= 0) then
             write (0,'(/A,I9)') ' end OF OBS unit: ',iunit
             write (0,'(/A,I9)') ' iostat:          ',iost
            exit reports
         end if
         write(filtered_obs_unit, fmt = fmt_info)           &
                      platform%info%platform,    &
                      platform%info%date_char,   &
                      platform%info%name,        &
                      platform%info%levels,      &
                      platform%info%lat,         &
                      platform%info%lon,         &
                      platform%info%elv,         &
                      platform%info%id

         !  Read surface Info
         !  -------------------
         read (iunit, fmt = fmt_srfc)  &
            platform%loc%slp%inv, platform%loc%slp%qc, platform%loc%slp%error, &
            platform%loc%pw%inv, platform%loc%pw%qc, platform%loc%pw%error
         write(filtered_obs_unit, fmt = fmt_srfc)  &
            platform%loc%slp%inv, platform%loc%slp%qc, platform%loc%slp%error, &
            platform%loc%pw%inv, platform%loc%pw%qc, platform%loc%pw%error

         ! == levels < 1 and not GPSPW, go back to reports

         if ((platform%info%levels < 1) .AND.            &
             (index(platform%info%platform, 'GPSPW') <= 0)) then
              cycle reports
         end if

         !     read EACH LEVELS
         !     ----------------
         do i = 1, platform%info%levels

            platform%each (i) = each_level_type(missing_r, missing, -1.0, & ! height
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! u
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! v
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! p
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! t
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! q
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! rh
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r), & ! td
                          field_type(missing_r, missing_data, missing_r, missing_r, missing_r))  ! speed 
            uinv = missing_r
            vinv = missing_r
            tinv = missing_r
            qinv = missing_r

            read (unit = iunit, fmt = trim (fmt_each_inv)) &
               platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
               platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
               platform%each(i)%speed%error, &
               platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
               platform%each(i)%height,       &
               platform%each(i)%height_qc,    &
               height_error,                   &
               platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
               platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
               platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error, &
               uinv, vinv, tinv, qinv

            write (unit = filtered_obs_unit, fmt = trim (fmt_each_inv)) &
               platform%each(i)%p%inv, platform%each(i)%p%qc, platform%each(i)%p%error, &
               platform%each(i)%speed%inv, platform%each(i)%speed%qc, &
               platform%each(i)%speed%error, &
               platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
               platform%each(i)%height,       &
               platform%each(i)%height_qc,    &
               height_error,                   &
               platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
               platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
               platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error, &
               uinv, vinv, tinv, qinv
         end do
      end do reports                  !  Loop over reports              
      close (iunit)
   end do !  Loop over all data files

   call da_free_unit (iunit)

   if (trace_use) call da_trace_exit("da_final_write_modified_filtered_obs")

end subroutine da_final_write_modified_filtered_obs 



end module da_obs_io
