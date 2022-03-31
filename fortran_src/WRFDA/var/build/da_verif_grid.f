












program da_verif_grid 
   
   
   
   
   
   
   
   
   
   
   
   
   

   use da_verif_grid_control, only : control_main, control_times, control_vars, sub_domain, &
       max_3d_variables, max_2d_variables,num_vert_levels,verification_file_string,&
       missing,namelist_unit,time_series_unit, time_series_2d, profile_time_series_3d,&
       filename, stime, etime, hstart, hend, hdate, date, pdate, vert_levels, &
       nx, ny, nz, io_status,  debug1, debug2, verify_its_own_analysis, &
       num_verifying_experiments, verify_forecast_hour, domain, control_exp_dir, verif_dirs, &
       out_dirs,start_year, end_year, start_month, end_month, start_day, end_day, &
       start_hour, end_hour,start_minutes, end_minutes, start_seconds, end_seconds,interval_hour, &
       num3dvar, num2dvar, var3d, var2d, num_scores, score_names, vertical_type, &
       istart, iend, jstart, jend, unit_all, unit_land, unit_water

   use da_netcdf_interface, only : da_get_dims_cdf, da_get_gl_att_int_cdf, da_get_gl_att_real_cdf, &
      da_get_var_3d_real_cdf, da_get_var_2d_real_cdf, da_get_var_2d_int_cdf

   implicit none

   character (len=512) :: control_file, verif_file
   character (len=512) :: out_dir, outname_all, outname_land, outname_water

   integer, parameter                   :: imiss = -99
   real, parameter                      :: rmiss = -99.99
   integer                              :: time_loop_count
   integer                              :: time(6), ptime(6)
   integer                              :: nx1, ny1, nz1
   integer                              :: nx2, ny2, nz2
   integer                              :: i,k
   integer                              :: ivar, iexp, iscore
   character (len=10)                   :: sdate
   character (len=20)                   :: file_string, domain_string, out_hr
   logical, allocatable,dimension(:)    :: first_score
   logical                              :: exist_file1, exist_file2

   real, allocatable, dimension(:,:,:)  :: data_out1, data_out2
   real, allocatable, dimension(:,:,:)  :: u1, u2, v1, v2       
 
   real, allocatable, dimension(:,:,:)  :: sum3d, asum3d, sqr3d, diff, absdiff, sqdiff
   real, allocatable, dimension(:,:)    :: score_avg_prof
   real, allocatable, dimension(:)      :: avg_prof
   real, allocatable, dimension(:,:)    :: landmask  
   integer, parameter                   :: island = 1, iswater = 0
   integer, allocatable, dimension(:)   :: num_counter
   integer, allocatable, dimension(:,:,:)  :: mask_count
   real, allocatable, dimension(:,:,:)  :: mask_stats 
                                                      
                                                      
   

   verify_forecast_hour = 0

   debug1 = .false.
   debug2 = .false.
   vertical_type = 'p'


   num3dvar=max_3d_variables
   var3d(1)='U'
   var3d(2)='V'
   var3d(3)='T'
   var3d(4)='QVAPOR'
   var3d(5)='Z'
   var3d(6)='WV'


   num2dvar=max_2d_variables
   var2d(1)='SLP'
   var2d(2)='PSFC'
   var2d(3)='U10M'
   var2d(4)='V10M'
   var2d(5)='T2M'
   var2d(6)='Q2M'
   var2d(7)='MU'


   num_scores = 3
   score_names(1) = 'BIAS'
   score_names(2) = 'RMSE'
   score_names(3) = 'ABIAS'
                                     
   domain = 1
   verification_file_string = 'wrfout'
   out_hr ='_00'



   io_status = 0

                                                                      
   open(unit = namelist_unit, file = 'namelist.in', &
          status = 'old' , access = 'sequential', &
          form   = 'formatted', action = 'read', &
          iostat = io_status )

                                                                   
   if(io_status /= 0) then
      print *, 'Error to open namelist.in file: '
   else
      read(unit=namelist_unit, nml = control_main , iostat = io_status)

      if(io_status /= 0) then
         print *, 'Error to read control_main. Stopped.'
         stop
      endif

                                                                   
      read(unit=namelist_unit, nml = control_times , iostat = io_status )
      
      if(io_status /= 0) then
         print *, 'Error to read control_times Stopped.'
         stop
      endif
   
      read(unit=namelist_unit, nml = control_vars , iostat = io_status )
     
      if(io_status /= 0) then
         print *, 'Error to read control_vars Stopped.'
         stop
      endif
   
      read(unit=namelist_unit, nml = sub_domain , iostat = io_status )
     
      if(io_status /= 0) then
         print *, 'Error reading sub_domain'
         istart = 1
         iend   = 10000
         jstart = 1
         jend   = 10000
         
      endif
   
      close(unit=namelist_unit)
   endif



    write(domain_string, fmt ='("_d",i2.2,"_")') domain
    file_string = trim(verification_file_string)//trim(domain_string) 
    write(out_hr, fmt ='("_",i2.2)') verify_forecast_hour
    allocate(first_score(num_scores))


   stime(1) = start_year
   stime(2) = start_month
   stime(3) = start_day  
   stime(4) = start_hour 
   stime(5) = start_minutes
   stime(6) = start_seconds


   call build_hdate(hstart, stime )

   etime(1) = end_year
   etime(2) = end_month
   etime(3) = end_day  
   etime(4) = end_hour 
   etime(5) = end_minutes
   etime(6) = end_seconds

   call build_hdate(hend, etime )

                                                                   
   if ( hend < hstart ) then
      print*, '****************************************************************'
      print*, 'End time is before the start time'
      print*, ' Start time = ', hstart ,'  & End time = ', hend
      print*, '****************************************************************'
      stop
   endif

   hdate = hstart 
   call build_hdate(sdate, stime )

   filename = trim(file_string)//hdate
   

   loop_verif_exp : do iexp = 1, num_verifying_experiments

      first_score  = .true.
      out_dir = trim(out_dirs(iexp))//'/'



      allocate( avg_prof(num_vert_levels) )
      allocate( num_counter(num_vert_levels) )
      allocate( score_avg_prof(num_scores, num_vert_levels) )

      loop_3d : do ivar=1,num3dvar



        profile_time_series_3d = trim(out_dir)//trim(var3d(ivar))//'_time_series'//trim(out_hr) 
        open(time_series_unit, file=profile_time_series_3d,form='formatted', &
                               status='unknown')
       
        outname_all   = trim(profile_time_series_3d)//'_sub'
        outname_land  = trim(profile_time_series_3d)//'_sub_land'
        outname_water = trim(profile_time_series_3d)//'_sub_water'
        open(unit_all,   file=trim(outname_all),  form='formatted', status='unknown')
        open(unit_land,  file=trim(outname_land), form='formatted', status='unknown')
        open(unit_water, file=trim(outname_water),form='formatted', status='unknown')



        time_loop_count = 0
        hdate = hstart
        time = stime
     
        time_loop_3d : do 
              

           call build_hdate(hdate, time)
           print*,' processing exp: ',iexp,' 3d var: ',trim(var3d(ivar)),' for time: ',hdate

           if ( hdate > hend ) exit time_loop_3d

           call build_hdate(date, time )
           ptime = time
           call advance_date(ptime,-verify_forecast_hour) 
           call build_hdate(pdate, ptime)


           filename = trim(file_string)//hdate

           if( verify_its_own_analysis) then
           control_file = trim(verif_dirs(iexp))//'/'//date//'/'//trim(filename) 
           else
           control_file = trim(control_exp_dir)//'/'//date//'/'//trim(filename) 
           endif

           verif_file = trim(verif_dirs(iexp))//'/'//pdate//'/'//trim(filename)

           inquire(file=trim(control_file), exist=exist_file1)
           inquire(file=trim(verif_file), exist=exist_file2)
           if ( (.not. exist_file1) .or. (.not. exist_file2) ) then
              avg_prof       = rmiss
              num_counter    = imiss
              score_avg_prof = rmiss
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, &
                                 time_series_unit) 
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, unit_all)
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, unit_land)
              call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, unit_water)
              call advance_date(time,interval_hour) 
              if ( .not. exist_file1 ) print*, trim(control_file), ' file not found'
              if ( .not. exist_file2 ) print*, trim(verif_file), ' file not found'
              print*, 'skipping date ', date
              cycle time_loop_3d
           end if





           call get_dimensions(control_file,nx1,ny1,nz1)
           call get_dimensions(verif_file,nx2,ny2,nz2)

           if ( nx1 /= nx2 .or. ny1 /= ny2 .or. nz1 /= nz2 ) then
              print*, '********************************************************'
              print*, 'Dimension mismatch between files of the experiments ....' 
              print*, '********************************************************'
              stop
           else
             nx = nx1
             ny = ny1
             nz = nz1
             if (time_loop_count == 0 ) then
               allocate( sum3d(nx,ny,num_vert_levels))
               allocate( asum3d(nx,ny,num_vert_levels))
               allocate( sqr3d(nx,ny,num_vert_levels))
               sum3d = 0.0
               asum3d = 0.0
               sqr3d = 0.0

               if ( ivar == 1 ) then
                  allocate(landmask(nx,ny))
                  call da_get_var_2d_real_cdf(verif_file,"LANDMASK", landmask, nx, ny, 1, .false.) 
                  istart = max(1, istart)
                  iend   = min(nx, iend)
                  jstart = max(1, jstart)
                  jend   = min(ny, jend)
               end if

             endif
           endif
           allocate(diff(nx,ny,num_vert_levels))
           allocate(absdiff(nx,ny,num_vert_levels))
           allocate(sqdiff(nx,ny,num_vert_levels))

           if( trim(var3d(ivar)) .eq. "WV" ) then
           allocate ( u1 (nx, ny, num_vert_levels) )
           allocate ( u2 (nx, ny, num_vert_levels) )
           allocate ( v1 (nx, ny, num_vert_levels) )
           allocate ( v2 (nx, ny, num_vert_levels) )
           call compute_wind_3d( control_file, nx, ny, nz, num_vert_levels, 1, &
                                 vert_levels, vertical_type, missing, u1, v1, debug1 )
           call compute_wind_3d( verif_file, nx, ny, nz, num_vert_levels, 1, &
                                 vert_levels, vertical_type, missing, u2, v2, debug1 )
           call get_diffs_wind(u1, v1, u2, v2, diff, absdiff, sqdiff, nx, ny, &
                          num_vert_levels, missing)
           deallocate(u1, v1, u2, v2)
           else
           allocate ( data_out1 (nx, ny, num_vert_levels) )
           allocate ( data_out2 (nx, ny, num_vert_levels) )
           call compute_data_3d( control_file, trim(var3d(ivar)), len_trim(var3d(ivar)),    &
                                nx, ny, nz, num_vert_levels, 1, vert_levels, &
                                vertical_type, missing, data_out1, debug1 )

           call compute_data_3d( verif_file, trim(var3d(ivar)), len_trim(var3d(ivar)),    &
                                nx, ny, nz, num_vert_levels, 1, vert_levels, &
                                vertical_type, missing, data_out2, debug2 )

           call get_diffs(data_out1, data_out2, diff, absdiff, sqdiff, nx, ny,  &
                          num_vert_levels, missing)
           deallocate(data_out1)
           deallocate(data_out2)
           end if

           do iscore = 1, num_scores
             if ( trim(score_names(iscore)) == 'BIAS' ) then
                call domain_average( diff, avg_prof, num_counter, nx, ny, num_vert_levels, missing,0)
             elseif ( trim(score_names(iscore)) == 'RMSE' ) then
                call domain_average( sqdiff, avg_prof, num_counter, nx, ny, num_vert_levels, missing,1)
             elseif ( trim(score_names(iscore)) == 'ABIAS' ) then
                call domain_average( absdiff, avg_prof, num_counter, nx, ny, num_vert_levels, missing,0)
             endif
             score_avg_prof(iscore,:) = avg_prof(:)
           enddo
           call write_profile(date, score_avg_prof, num_counter, num_vert_levels, num_scores, &
                              time_series_unit) 

           allocate(mask_stats(3,3,nz))
           allocate(mask_count(3,3,nz))
           call mask_domain_average( diff,    mask_stats(1,:,:), mask_count(1,:,:), nx, ny, num_vert_levels, missing,0)
           call mask_domain_average( sqdiff,  mask_stats(2,:,:), mask_count(2,:,:), nx, ny, num_vert_levels, missing,1)
           call mask_domain_average( absdiff, mask_stats(3,:,:), mask_count(3,:,:), nx, ny, num_vert_levels, missing,0)
           call write_profile(date, mask_stats(:,1,:), mask_count(1,1,:), num_vert_levels, 3, unit_all) 
           call write_profile(date, mask_stats(:,2,:), mask_count(1,2,:), num_vert_levels, 3, unit_land) 
           call write_profile(date, mask_stats(:,3,:), mask_count(1,3,:), num_vert_levels, 3, unit_water) 
           deallocate(mask_stats)
           deallocate(mask_count)
           
           call get_sum(sum3d,diff,nx,ny,num_vert_levels,missing)
           call get_sum(asum3d,absdiff,nx,ny,num_vert_levels,missing)
           call get_sum(sqr3d,sqdiff,nx,ny,num_vert_levels,missing)

           deallocate(diff)
           deallocate(absdiff)
           deallocate(sqdiff)

           time_loop_count = time_loop_count + 1

           call advance_date(time,interval_hour) 

        enddo  time_loop_3d       

        close(time_series_unit)
        close(unit_all)
        close(unit_land)
        close(unit_water)
        deallocate(sum3d)
        deallocate(asum3d)
        deallocate(sqr3d)

     enddo loop_3d
     print*, ' successful completion of loop_3d '

     deallocate( avg_prof )
     deallocate( num_counter )
     deallocate( score_avg_prof )




     allocate( avg_prof(1) )
     allocate( num_counter(1) )
     allocate( score_avg_prof(num_scores, 1) )

     loop_2d : do ivar = 1, num2dvar



        time_series_2d  = trim(out_dir)//trim(var2d(ivar))//'_time_series'//trim(out_hr) 
        open(time_series_unit, file=time_series_2d,form='formatted', &
                               status='unknown')
       
        outname_all   = trim(time_series_2d)//'_sub'
        outname_land  = trim(time_series_2d)//'_sub_land'
        outname_water = trim(time_series_2d)//'_sub_water'
        open(unit_all,   file=trim(outname_all),  form='formatted', status='unknown')
        open(unit_land,  file=trim(outname_land), form='formatted', status='unknown')
        open(unit_water, file=trim(outname_water),form='formatted', status='unknown')



        
        time_loop_count = 0
        hdate = hstart
        time = stime
     
        time_loop_2d : do 
              

           call build_hdate(hdate, time )
           print*,' processing exp: ',iexp,' 2d var: ',trim(var2d(ivar)),' for time: ',hdate

           if ( hdate > hend ) exit time_loop_2d

           call build_hdate(date, time )
           ptime = time
           call advance_date(ptime,-verify_forecast_hour) 
           call build_hdate(pdate, ptime)

           filename = trim(file_string)//hdate
           if( verify_its_own_analysis) then
           control_file = trim(verif_dirs(iexp))//'/'//date//'/'//trim(filename) 
           else
           control_file = trim(control_exp_dir)//'/'//date//'/'//trim(filename) 
           endif

           verif_file = trim(verif_dirs(iexp))//'/'//pdate//'/'//trim(filename)

           inquire(file=trim(control_file), exist=exist_file1)
           inquire(file=trim(verif_file), exist=exist_file2)
           if ( (.not. exist_file1) .or. (.not. exist_file2) ) then
              avg_prof       = rmiss
              num_counter    = imiss
              score_avg_prof = rmiss
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, &
                                 time_series_unit) 
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, unit_all)
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, unit_land)
              call write_profile(date, score_avg_prof, num_counter, 1, num_scores, unit_water)
              call advance_date(time,interval_hour) 
              if ( .not. exist_file1 ) print*, trim(control_file), ' file not found'
              if ( .not. exist_file2 ) print*, trim(verif_file), ' file not found'
              print*, 'skipping date ', date
              cycle time_loop_2d
           end if



           call get_dimensions(control_file,nx1,ny1,nz1)
           call get_dimensions(verif_file,nx2,ny2,nz2)

           if ( nx1 /= nx2 .or. ny1 /= ny2 .or. nz1 /= nz2 ) then
              print*, '********************************************************'
              print*, 'Dimension mismatch between files of the experiments ....' 
              print*, '********************************************************'
              stop
           else
             nx = nx1
             ny = ny1
             nz = nz1
             if (time_loop_count == 0 ) then
               allocate( sum3d(nx,ny,1))
               allocate( asum3d(nx,ny,1))
               allocate( sqr3d(nx,ny,1))
               sum3d = 0.0
               asum3d = 0.0
               sqr3d = 0.0
             endif
           endif

           allocate(data_out1(nx, ny, 1))
           allocate(data_out2(nx, ny, 1))
 
           call g_output_2d (control_file, 1, trim(var2d(ivar)), len_trim(var2d(ivar)), &
                             nx, ny, nz, data_out1, debug1)

           call g_output_2d (verif_file, 1, trim(var2d(ivar)), len_trim(var2d(ivar)), &
                             nx, ny, nz, data_out2, debug2)

           allocate(diff(nx,ny,1))
           allocate(absdiff(nx,ny,1))
           allocate(sqdiff(nx,ny,1))
           call get_diffs(data_out1, data_out2, diff, absdiff, sqdiff, nx, ny, 1, missing)
           deallocate(data_out1)
           deallocate(data_out2)

           do iscore = 1, num_scores
             if ( trim(score_names(iscore)) == 'BIAS' ) then
                call domain_average( diff, avg_prof, num_counter, nx, ny, 1, missing,0)
             elseif ( trim(score_names(iscore)) == 'RMSE' ) then
                call domain_average( sqdiff, avg_prof, num_counter, nx, ny, 1, missing,1)
             elseif ( trim(score_names(iscore)) == 'ABIAS' ) then
                call domain_average( absdiff, avg_prof, num_counter, nx, ny, 1, missing,0)
             endif
             score_avg_prof(iscore,:) = avg_prof(:)
           enddo
           call write_profile(date, score_avg_prof, num_counter, 1, num_scores, &
                              time_series_unit)

           allocate(mask_stats(3,3,1))
           allocate(mask_count(3,3,1))
           call mask_domain_average( diff,    mask_stats(1,:,:), mask_count(1,:,:), nx, ny, 1, missing,0)
           call mask_domain_average( sqdiff,  mask_stats(2,:,:), mask_count(2,:,:), nx, ny, 1, missing,1)
           call mask_domain_average( absdiff, mask_stats(3,:,:), mask_count(3,:,:), nx, ny, 1, missing,0)
           call write_profile(date, mask_stats(:,1,:), mask_count(1,1,:), 1, 3, unit_all) 
           call write_profile(date, mask_stats(:,2,:), mask_count(1,2,:), 1, 3, unit_land) 
           call write_profile(date, mask_stats(:,3,:), mask_count(1,3,:), 1, 3, unit_water) 
           deallocate(mask_stats)
           deallocate(mask_count)
           
           call get_sum(sum3d,diff,nx,ny,1,missing)
           call get_sum(asum3d,absdiff,nx,ny,1,missing)
           call get_sum(sqr3d,sqdiff,nx,ny,1,missing)

           deallocate(diff)
           deallocate(absdiff)
           deallocate(sqdiff)
           time_loop_count = time_loop_count + 1

           call advance_date(time,interval_hour) 


        enddo  time_loop_2d
        
        close(time_series_unit)
        close(unit_all)
        close(unit_land)
        close(unit_water)
        deallocate(sum3d)
        deallocate(asum3d)
        deallocate(sqr3d)

     enddo loop_2d

     deallocate( avg_prof )
     deallocate( num_counter )
     deallocate( score_avg_prof )
     deallocate( landmask )

     print*, ' successful completion of loop_2d '
   print*,' Finished Experiment : ', trim(verif_dirs(iexp))
   enddo loop_verif_exp



   contains  

   subroutine advance_date( time, delta )

      implicit none

      integer, intent(inout) :: time(6)         
      integer, intent(in)    :: delta

      integer                :: ccyy, mm, dd, hh
      integer, dimension(12) :: mmday



      mmday = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      mmday(2) = 28

    ccyy = time(1)
    mm   = time(2)
    dd   = time(3)
    hh   = time(4)


    hh = hh + delta

    do while (hh < 0)
       hh = hh + 24
       call change_date ( ccyy, mm, dd, -1 )
    end do

    do while (hh > 23)
       hh = hh - 24
       call change_date ( ccyy, mm, dd, 1 )
    end do


    time(1) = ccyy
    time(2) = mm
    time(3) = dd
    time(4) = hh

    end subroutine advance_date

    subroutine change_date ( ccyy, mm, dd, delta)
      integer, intent(inout) :: ccyy, mm, dd
      integer, intent(in) :: delta

      integer, dimension(12) :: mmday
      mmday = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      mmday(2) = 28

      if (mod(ccyy,4) == 0) then
         mmday(2) = 29

         if ( mod(ccyy,100) == 0) then
            mmday(2) = 28
         endif

         if(mod(ccyy,400) == 0) then
            mmday(2) = 29
         end if
      endif

      dd = dd + delta

      if(dd == 0) then
         mm = mm - 1

         if(mm == 0) then
            mm = 12         
            ccyy = ccyy - 1
         endif

         dd = mmday(mm)
      elseif ( dd .gt. mmday(mm) ) then
         dd = 1
         mm = mm + 1
         if(mm > 12 ) then
            mm = 1
            ccyy = ccyy + 1
         end if
      end if
   end subroutine change_date

   subroutine build_hdate(hdate, time)







      integer, intent(in) :: time(6) 

      character*(*), intent(out) :: hdate 


      integer iyr     
      integer imo     
      integer idy     
      integer ihr     
      integer imi     
      integer isc     


      integer hlen 

      hlen = len(hdate)
      iyr = time(1)
      imo = time(2)
      idy = time(3)
      ihr = time(4)
      imi = time(5)
      isc = time(6)

      if (hlen.eq.19) then
         write(hdate,19) iyr, imo, idy, ihr, imi, isc
 19      format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2,':',i2.2)

      elseif (hlen.eq.16) then
         write(hdate,16) iyr, imo, idy, ihr, imi
 16      format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2)

      elseif (hlen.eq.13) then
         write(hdate,13) iyr, imo, idy, ihr
 13      format(i4,'-',i2.2,'-',i2.2,'_',i2.2)

      elseif (hlen.eq.10) then
         write(hdate,10) iyr, imo, idy, ihr
 10      format(i4,i2.2,i2.2,i2.2)
      endif

      return
      end subroutine build_hdate


  subroutine write_profile(date, profile, counter, nlevel, nscore, out_unit)
  
  integer, intent(in)                 :: nlevel, nscore, out_unit
  real, intent(in), dimension(:,:)    :: profile
  integer, intent(in), dimension(:)   :: counter
  character (len=10), intent(in)      :: date
  write(out_unit,fmt='(a10,1x,100(i8,1x,3(f14.8,1x)))') date,  &
                              (counter(k), (profile(i,k),i=1,nscore),k=1,nlevel)

  end subroutine write_profile
  

  subroutine time_calc( time, timestamp, datestamp, debug , tdef,it)

  implicit none

  character (len=19), intent(in) :: time
  character (len=35), intent(inout)             :: tdef
  integer, intent(out)           :: timestamp, datestamp
  logical, intent(in)            :: debug

   integer, intent(in) :: it
  integer :: hours, minutes, seconds, year, month, day,hour1,hourint
  integer :: mins1,minsint

  save hourint 
  save minsint

  read(time(18:19),*) seconds
  read(time(15:16),*) minutes
  read(time(12:13),*) hours
  read(time(1:4),*)   year
  read(time(6:7),*)   month
  read(time(9:10),*)  day

  if(debug) write(6,*) ' day, month, year, hours, minutes, seconds '
  if(debug) write(6,*) day, month, year, hours, minutes, seconds 

  if ( it == 1) then
    write (tdef(19:20),'(i2)') hours
    if ( day < 10 ) then
      write (tdef(23:23),'(i1)') day
    else
      write (tdef(22:23),'(i2)') day
    endif
    write (tdef(27:30),'(i4)') year
    if (month == 1) write (tdef(24:26),'(a3)') 'jan'
    if (month == 2) write (tdef(24:26),'(a3)') 'feb'
    if (month == 3) write (tdef(24:26),'(a3)') 'mar'
    if (month == 4) write (tdef(24:26),'(a3)') 'apr'
    if (month == 5) write (tdef(24:26),'(a3)') 'may'
    if (month == 6) write (tdef(24:26),'(a3)') 'jun'
    if (month == 7) write (tdef(24:26),'(a3)') 'jul'
    if (month == 8) write (tdef(24:26),'(a3)') 'aug'
    if (month == 9) write (tdef(24:26),'(a3)') 'sep'
    if (month ==10) write (tdef(24:26),'(a3)') 'oct'
    if (month ==11) write (tdef(24:26),'(a3)') 'nov'
    if (month ==12) write (tdef(24:26),'(a3)') 'dec'
    hour1=hours
    mins1=minutes
  elseif ( it == 2) then
    hourint = abs(hours-hour1)
    minsint = abs(minutes-mins1)
    if (hourint == 0 ) then
      if (minsint == 0 ) minsint = 1
      if(debug) write(6,*) "interval is",minsint
      write (tdef(34:35),'(a2)') "mn"
      write (tdef(32:33),'(i2)') minsint
      if(debug) write(6,*) "TDEF is",tdef
    else
      if(debug) write(6,*) "Interval is",hourint
      write (tdef(32:33),'(i2)') hourint
      if(debug) write(6,*) "TDEF is",tdef
    endif
  endif

  timestamp = seconds+100*minutes+10000*hours

  if((year > 1800) .and. (year < 2000)) year = year-1900
  if((year >= 2000)) year = year-2000

  if(month >= 2) day = day+31  
  if(month >= 3) day = day+28  
  if(month >= 4) day = day+31  
  if(month >= 5) day = day+30  
  if(month >= 6) day = day+31  
  if(month >= 7) day = day+30  
  if(month >= 8) day = day+31  
  if(month >= 9) day = day+31  
  if(month >= 10) day = day+30 
  if(month >= 11) day = day+31 
  if(month >= 12) day = day+30 
  if((month > 2) .and. (mod(year,4) == 0)) day = day+1  

  datestamp = (year)*1000 + day


  if(debug) then
    write(6,*) ' time, timestamp, datestamp ',time(1:19),timestamp,datestamp
  endif

  end subroutine time_calc



  subroutine g_output_3d (file, file_time_index, var, length_var,            &
                          nx, ny, nz, data_out, debug)
  implicit none

  character (len=*), intent(in)             ::   file
  integer, intent(in)                       ::   file_time_index
  character (len=*), intent(in)             ::   var
  integer, intent(in)                       ::   length_var
  integer , intent(in)                      ::   nx, ny, nz           
  real, intent(out), dimension(:,:,:)       ::   data_out
  logical, intent(in)                                   ::   debug
  real,    allocatable, dimension(:,:,:)    ::   data_tmp, data_tmp2
  real,    allocatable, dimension(:,:,:)    ::   u, v
  real,    allocatable, dimension(:,:)      ::   xlat, xlon

  real,    allocatable, dimension(:,:,:)    ::   ph, phb  
  real,    allocatable, dimension(:,:,:)    ::   p, pb  
  real,    allocatable, dimension(:,:,:)    ::   t, qv 
  integer                                   ::   map_proj
  real                                      ::   cen_lon, truelat1, truelat2


   REAL    , PARAMETER :: g            = 9.81  
   REAL    , PARAMETER :: r_d          = 287.
   REAL    , PARAMETER :: r_v          = 461.6
   REAL    , PARAMETER :: cp           = 7.*r_d/2.
   REAL    , PARAMETER :: cv           = cp-r_d
   REAL    , PARAMETER :: cliq         = 4190.
   REAL    , PARAMETER :: cice         = 2106.
   REAL    , PARAMETER :: psat         = 610.78
   REAL    , PARAMETER :: rcv          = r_d/cv
   REAL    , PARAMETER :: rcp          = r_d/cp
   REAL    , PARAMETER :: c2           = cp * rcv
   REAL    , PARAMETER :: T0           = 273.16

   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL    , PARAMETER :: cvovcp       = 1./cpovcv


  if(debug) then
    write(6,*) ' calculations for variable ',var
  end if

       if(var == 'U' ) then

          allocate ( data_tmp(nx+1,ny,nz) )
          call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,            &
                                file_time_index, debug  )
          data_out = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))

          deallocate ( data_tmp )

  else if(var == 'V' ) then

          allocate ( data_tmp(nx,ny+1,nz) )
          call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,            &
                                file_time_index, debug  )
          data_out = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
          deallocate ( data_tmp )

  else if(var == 'UMET' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate (        u(nx,ny,nz)   )
              allocate (        v(nx,ny,nz)   )
              allocate (     xlat(nx,ny)      )             
              allocate (     xlon(nx,ny)      )             
    
              allocate ( data_tmp(nx+1,ny,nz) )
              call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,        &
                                    file_time_index, debug  )
              u = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))
              deallocate ( data_tmp )
    
              allocate ( data_tmp(nx,ny+1,nz) )
              call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,        &
                                    file_time_index, debug  )
              v = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
              deallocate ( data_tmp )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u,v,nx,ny,nz,var,                                &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u    )
              deallocate ( v    )

          ELSE

              allocate ( data_tmp(nx+1,ny,nz) )
              call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,        &
                                    file_time_index, debug  )
              data_out = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))
              deallocate ( data_tmp )
    
          ENDIF

  else if(var == 'VMET' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate (        u(nx,ny,nz)   )
              allocate (        v(nx,ny,nz)   )
              allocate (     xlat(nx,ny)      )             
              allocate (     xlon(nx,ny)      )             
    
              allocate ( data_tmp(nx+1,ny,nz) )
              call da_get_var_3d_real_cdf( file,"U", data_tmp, nx+1, ny, nz,        &
                                    file_time_index, debug  )
              u = 0.5*(data_tmp(1:nx,:,:)+data_tmp(2:nx+1,:,:))
              deallocate ( data_tmp )
    
              allocate ( data_tmp(nx,ny+1,nz) )
              call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,        &
                                    file_time_index, debug  )
              v = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
              deallocate ( data_tmp )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u,v,nx,ny,nz,var,                                &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u    )
              deallocate ( v    )

          ELSE

              allocate ( data_tmp(nx,ny+1,nz) )
              call da_get_var_3d_real_cdf( file,"V", data_tmp, nx, ny+1, nz,        &
                                    file_time_index, debug  )
              data_out = 0.5*(data_tmp(:,1:ny,:)+data_tmp(:,2:ny+1,:))
              deallocate ( data_tmp )
 
          ENDIF

  else if(var == 'W' ) then

          allocate ( data_tmp(nx,ny,nz+1) )
          call da_get_var_3d_real_cdf( file,"W", data_tmp, nx, ny, nz+1,            &
                                file_time_index, debug  )
          data_out = 0.5*(data_tmp(:,:,1:nz)+data_tmp(:,:,2:nz+1))
          deallocate ( data_tmp )
 
  else if(var == 'P' ) then

          allocate (  p(nx,ny,nz) )
          allocate ( pb(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          data_out = (p+pb)*.01

          deallocate (  p )
          deallocate ( pb )             
 
  else if(var == 'Z' ) then

          allocate (  ph(nx,ny,nz+1) )
          allocate ( phb(nx,ny,nz+1) )             

          call da_get_var_3d_real_cdf( file,"PH", ph, nx, ny, nz+1,                 &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PHB", phb, nx, ny, nz+1,               &
                                file_time_index, debug  ) 
          ph = (ph+phb)/9.81
          data_out = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))

          deallocate (  ph )
          deallocate ( phb )             

  else if(var == 'THETA' ) then

          call da_get_var_3d_real_cdf( file,"T", data_out, nx, ny, nz,              &
                                file_time_index, debug  )
          data_out = data_out + 300.

  else if(var == 'T' ) then

          allocate (        p(nx,ny,nz) )
          allocate (       pb(nx,ny,nz) )             
          allocate ( data_tmp(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", data_tmp, nx, ny, nz,              &
                                file_time_index, debug  )
          data_out = (data_tmp+300.)*(p/p1000mb)**rcp

          deallocate (  p       )
          deallocate ( pb       )             
          deallocate ( data_tmp )

  else if(var == 'TC' ) then

          allocate (        p(nx,ny,nz) )
          allocate (       pb(nx,ny,nz) )             
          allocate ( data_tmp(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", data_tmp, nx, ny, nz,              &
                                file_time_index, debug  )
          data_out = (data_tmp+300.)*(p/p1000mb)**rcp -T0

          deallocate (  p       )
          deallocate ( pb       )             
          deallocate ( data_tmp )

  else if(var == 'TD' ) then

          allocate (        p(nx,ny,nz) )
          allocate (       pb(nx,ny,nz) )             
          allocate (       qv(nx,ny,nz) )             
          allocate ( data_tmp(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"QVAPOR", qv, nx, ny, nz,               &
                                file_time_index, debug  )

          data_tmp = qv*(p/100.)/(0.622+qv)
          data_tmp = AMAX1(data_tmp,0.001)
          data_out = (243.5*log(data_tmp)-440.8)/(19.48-log(data_tmp))

          deallocate (  p       )
          deallocate ( pb       )             
          deallocate ( qv       )             
          deallocate ( data_tmp )

  else if(var == 'RH' ) then

          allocate (         p(nx,ny,nz) )
          allocate (        pb(nx,ny,nz) )             
          allocate (        qv(nx,ny,nz) )             
          allocate (         t(nx,ny,nz) )             
          allocate (  data_tmp(nx,ny,nz) )             
          allocate ( data_tmp2(nx,ny,nz) )             

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny, nz,                   &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", t, nx, ny, nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"QVAPOR", qv, nx, ny, nz,               &
                                file_time_index, debug  )

          t = (t+300.)*(p/p1000mb)**rcp
          data_tmp2 = 10.*0.6112*exp(17.67*(t-T0)/(t-29.65))
          data_tmp  = 0.622*data_tmp2/(0.01 * p -  (1.-0.622)*data_tmp2)
          data_out  = 100.*AMAX1(AMIN1(qv/data_tmp,1.0),0.0)


          deallocate (  p        )
          deallocate ( pb        )             
          deallocate ( qv        )             
          deallocate ( t         )             
          deallocate ( data_tmp  )
          deallocate ( data_tmp2 )

  else 
          call da_get_var_3d_real_cdf( file,var(1:length_var),                      &
                                    data_out, nx,ny,nz,                          &
                                    file_time_index, debug  )
  endif


  end subroutine g_output_3d



  subroutine g_output_2d (file, file_time_index, var, length_var,            &
                          nx, ny, nz, data_out, debug)
  implicit none

  character (len=*), intent(in)             ::   file
  integer, intent(in)                       ::   file_time_index
  character (len=*), intent(in)             ::   var
  integer, intent(in)                       ::   length_var
  integer, intent(in)                       ::   nx, ny, nz           
  real, intent(out), dimension(:,:,:)       ::   data_out
  logical, intent(in)                                   ::   debug
  integer, allocatable, dimension(:,:,:)    ::   data_int
  real,    allocatable, dimension(:,:,:)    ::   u10, v10
  real,    allocatable, dimension(:,:)      ::   psfc,t2m,q2m,mu
  real,    allocatable, dimension(:,:)      ::   xlat, xlon
  real,    allocatable, dimension(:,:,:)    ::   z,ph,phb  
  real,    allocatable, dimension(:,:,:)    ::   p,pb  
  real,    allocatable, dimension(:,:,:)    ::   ts,qv 
  integer                                   ::   map_proj
  real                                      ::   cen_lon, truelat1, truelat2

  if(debug) then
     write(6,*) ' calculations for variable ',var
  end if

       if(var == 'SLP') then

          allocate (   z(nx,ny,nz)   )
          allocate (  ph(nx,ny,nz+1) )
          allocate ( phb(nx,ny,nz+1) )             
          allocate (   p(nx,ny,nz)   )             
          allocate (  pb(nx,ny,nz)   )             
          allocate (  ts(nx,ny,nz)   )             
          allocate (  qv(nx,ny,nz)   )             

          call da_get_var_3d_real_cdf( file,"PH", ph, nx, ny,nz+1,                  &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PHB", phb, nx, ny,nz+1,                &
                                file_time_index, debug  ) 
          ph = (ph+phb)/9.81
          z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))

          call da_get_var_3d_real_cdf( file,"P", p, nx, ny,nz,                      &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"PB", pb, nx, ny,nz,                    &
                                file_time_index, debug  ) 
          p = p+pb

          call da_get_var_3d_real_cdf( file,"T", ts, nx, ny,nz,                     &
                                file_time_index, debug  )
          call da_get_var_3d_real_cdf( file,"QVAPOR", qv, nx, ny,nz,                &
                                file_time_index, debug  )

          call compute_seaprs (nx, ny, nz, z, ts, p, qv, data_out, debug)


          deallocate (   z )              
          deallocate (  ph )              
          deallocate ( phb )                         
          deallocate (   p )                       
          deallocate (  pb )                       
          deallocate (  ts )                       
          deallocate (  qv )                       

  else if(var == 'MU' ) then
          allocate ( mu(nx,ny) )
          call da_get_var_2d_real_cdf( file,"MU", mu, nx, ny,                 &
                                file_time_index, debug  )
          data_out(:,:,1) = mu(:,:)
          deallocate ( mu )

  else if(var == 'PSFC' ) then
          allocate ( psfc(nx,ny) )
          call da_get_var_2d_real_cdf( file,"PSFC", psfc, nx, ny,                 &
                                file_time_index, debug  )
          data_out(:,:,1) = psfc(:,:)
          deallocate ( psfc )

  else if(var == 'T2M'  ) then
          allocate ( t2m(nx,ny) )
          call da_get_var_2d_real_cdf( file,"T2", t2m, nx, ny,                  &
                                file_time_index, debug  )
          data_out(:,:,1) = t2m(:,:)
          deallocate ( t2m )

  else if(var == 'Q2M'  ) then
          allocate ( q2m(nx,ny) )
          call da_get_var_2d_real_cdf( file,"Q2", q2m, nx, ny,                  &
                                file_time_index, debug  )
          data_out(:,:,1) = q2m(:,:)
          deallocate ( q2m )

  else if(var == 'U10M' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate ( u10(nx,ny,1) )
              allocate ( v10(nx,ny,1) )
              allocate ( xlat(nx, ny) )             
              allocate ( xlon(nx, ny) )             
              call da_get_var_2d_real_cdf( file,"U10", u10, nx, ny,                 &
                                    file_time_index, debug  )
              call da_get_var_2d_real_cdf( file,"V10", v10, nx, ny,                 &
                                    file_time_index, debug  )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u10,v10,nx,ny,1,var,                             &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u10  )
              deallocate ( v10  )

          ELSE

              call da_get_var_2d_real_cdf( file,"U10", data_out, nx, ny,            &
                                    file_time_index, debug  )

          ENDIF

  else if(var == 'V10M' ) then

          call da_get_gl_att_int_cdf ( file, 'MAP_PROJ', map_proj, debug )

          IF ( map_proj == 1  .OR.  map_proj == 2 ) THEN

              allocate ( u10(nx,ny,1) )
              allocate ( v10(nx,ny,1) )
              allocate ( xlat(nx, ny) )             
              allocate ( xlon(nx, ny) )             
              call da_get_var_2d_real_cdf( file,"U10", u10, nx, ny,                 &
                                    file_time_index, debug  )
              call da_get_var_2d_real_cdf( file,"V10", v10, nx, ny,                 &
                                    file_time_index, debug  )
 
              call da_get_gl_att_real_cdf( file, 'STAND_LON', cen_lon, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT1', truelat1, debug )
              call da_get_gl_att_real_cdf( file, 'TRUELAT2', truelat2, debug )
              call da_get_var_2d_real_cdf( file, 'XLAT', xlat,nx,ny, 1,debug )
              call da_get_var_2d_real_cdf( file, 'XLONG',xlon,nx,ny, 1,debug )

              call rotate_wind (u10,v10,nx,ny,1,var,                             &
                                map_proj,cen_lon,xlat,xlon,                      &
                                truelat1,truelat2,data_out)

              deallocate ( xlat )             
              deallocate ( xlon )             
              deallocate ( u10  )
              deallocate ( v10  )

          ELSE

              call da_get_var_2d_real_cdf( file,"V10", data_out, nx, ny,            &
                                    file_time_index, debug  )

          ENDIF

  else if(var == 'XLONG' ) then
          call da_get_var_2d_real_cdf( file,var(1:length_var),                      &
                                    data_out, nx,ny,                             &
                                    file_time_index, debug  )
          WHERE ( data_out < 0.0 )
             data_out = data_out + 360.0
          ENDWHERE

  else if(var == 'IVGTYP' .or. var == 'ISLTYP') then

          allocate (data_int(nx,ny,1))
          call da_get_var_2d_int_cdf( file,var(1:length_var),                       &
                                    data_int, nx,ny,                             &
                                    file_time_index, debug  )
          data_out = data_int
          deallocate (data_int)

  else 
          call da_get_var_2d_real_cdf( file,var(1:length_var),                      &
                                    data_out, nx,ny,                             &
                                    file_time_index, debug  )
  endif


  end subroutine g_output_2d



  subroutine interp_to_z( data_in , nx_in , ny_in , nz_in , &
                              data_out, nx_out, ny_out, nz_out, &
                              z_in, z_out, missing_value, &
                              vertical_type, debug   )
  implicit none
  integer, intent(in)                                  :: nx_in , ny_in , nz_in 
  integer, intent(in)                                  :: nx_out, ny_out, nz_out
  real, intent(in)                                     :: missing_value
  real, dimension(nx_in , ny_in , nz_in ), intent(in ) :: data_in, z_in
  real, dimension(nx_out, ny_out, nz_out), intent(out) :: data_out
  real, dimension(nz_out), intent(in)                  :: z_out
  logical, intent(in)                                  :: debug
  character (len=1)                , intent(in)                     :: vertical_type

  real, dimension(nz_in)                               :: data_in_z, zz_in
  real, dimension(nz_out)                              :: data_out_z

  integer :: i,j,k

    do i=1,nx_in
    do j=1,ny_in

      do k=1,nz_in
        data_in_z(k) = data_in(i,j,k)
        zz_in(k) = z_in(i,j,k)
      enddo





      call interp_1d( data_in_z, zz_in, nz_in, &
                      data_out_z, z_out, nz_out, &
                      vertical_type, missing_value )

      do k=1,nz_out
        data_out(i,j,k) = data_out_z(k)
      enddo


    enddo
    enddo

  end subroutine interp_to_z



  subroutine interp_1d( a, xa, na, &
                        b, xb, nb, vertical_type, missing_value )
  implicit none
  integer, intent(in)              ::  na, nb
  real, intent(in), dimension(na)  :: a, xa
  real, intent(in), dimension(nb)  :: xb
  real, intent(out), dimension(nb) :: b
  real, intent(in)                 :: missing_value

  integer                          :: n_in, n_out
  logical                          :: interp
  real                             :: w1, w2
  character (len=1) ,intent(in)               :: vertical_type


  if ( vertical_type == 'p' ) then

  do n_out = 1, nb

    b(n_out) = missing_value
    interp = .false.
    n_in = 1

    do while ( (.not.interp) .and. (n_in < na) )

      if( (xa(n_in)   >= xb(n_out)) .and. &
          (xa(n_in+1) <= xb(n_out))        ) then
        interp = .true.
        w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
        w2 = 1. - w1
        b(n_out) = w1*a(n_in) + w2*a(n_in+1)
      end if
      n_in = n_in +1

    enddo

  enddo
  
  else

  do n_out = 1, nb

    b(n_out) = missing_value
    interp = .false.
    n_in = 1

    do while ( (.not.interp) .and. (n_in < na) )

      if( (xa(n_in)   <= xb(n_out)) .and. &
          (xa(n_in+1) >= xb(n_out))        ) then
        interp = .true.
        w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
        w2 = 1. - w1
        b(n_out) = w1*a(n_in) + w2*a(n_in+1)
      end if
      n_in = n_in +1

    enddo

  enddo
  
  endif

  end subroutine interp_1d
















  subroutine compute_seaprs ( nx , ny , nz  ,         &
                                  z, t , p , q ,          &
                                  sea_level_pressure,debug)

      IMPLICIT NONE

      INTEGER, intent(in) :: nx , ny , nz
      REAL, intent(in) ::    z(nx,ny,nz)
      REAL, intent(in)  ::  p(nx,ny,nz) , q(nx,ny,nz)
      REAL, intent(inout)  ::  t(nx,ny,nz)

      REAL, intent(out)   :: sea_level_pressure(nx,ny)
      INTEGER level(nx,ny)
      REAL t_surf(nx,ny) , t_sea_level(nx,ny)
      LOGICAL, intent(in) :: debug



      REAL R, G, GAMMA
      PARAMETER (R=287.04, G=9.81, GAMMA=0.0065)



      REAL    TC, PCONST
      PARAMETER (TC=273.16+17.5, PCONST = 10000)
      LOGICAL ridiculous_mm5_test
      PARAMETER (ridiculous_mm5_test = .TRUE.)




      INTEGER i , j , k
      INTEGER klo , khi


      REAL plo , phi , tlo, thi , zlo , zhi
      REAL p_at_pconst , t_at_pconst , z_at_pconst
      REAL z_half_lowest

      REAL    , PARAMETER :: cp           = 7.*R/2.
      REAL    , PARAMETER :: rcp          = R/cp
      REAL    , PARAMETER :: p1000mb      = 100000.

      LOGICAL  l1 , l2 , l3, found





      t(:,:,:) = (t(:,:,:)+300.)*(p(:,:,:)/p1000mb)**rcp

      DO j = 1 , ny
         DO i = 1 , nx
            level(i,j) = -1

            k = 1
            found = .false.
            do while( (.not. found) .and. (k.le.nz))
               IF ( p(i,j,k) .LT. p(i,j,1)-PCONST ) THEN
                  level(i,j) = k
                  found = .true.
               END IF
               k = k+1
            END DO

            IF ( level(i,j) .EQ. -1 ) THEN
            PRINT '(A,I4,A)','Troubles finding level ',   &
                        NINT(PCONST)/100,' above ground.'
            PRINT '(A,I4,A,I4,A)',                        &
                  'Problems first occur at (',i,',',j,')'
            PRINT '(A,F6.1,A)',                           &
                  'Surface pressure = ',p(i,j,1)/100,' hPa.'
            STOP 'Error_in_finding_100_hPa_up'
         END IF


         END DO
      END DO




      DO j = 1 , ny
         DO i = 1 , nx

            klo = MAX ( level(i,j) - 1 , 1      )
            khi = MIN ( klo + 1        , nz - 1 )

            IF ( klo .EQ. khi ) THEN
               PRINT '(A)','Trapping levels are weird.'
               PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
                            ': and they should not be equal.'
               STOP 'Error_trapping_levels'
            END IF

         plo = p(i,j,klo)
         phi = p(i,j,khi)
         tlo = t(i,j,klo)*(1. + 0.608 * q(i,j,klo) )
         thi = t(i,j,khi)*(1. + 0.608 * q(i,j,khi) )


         zlo = z(i,j,klo)
         zhi = z(i,j,khi)

         p_at_pconst = p(i,j,1) - pconst
         t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
         z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

         t_surf(i,j) = t_at_pconst*(p(i,j,1)/p_at_pconst)**(gamma*R/g)
         t_sea_level(i,j) = t_at_pconst+gamma*z_at_pconst

         END DO
      END DO




      IF ( ridiculous_mm5_test ) THEN
         DO j = 1 , ny
            DO i = 1 , nx
               l1 = t_sea_level(i,j) .LT. TC
               l2 = t_surf     (i,j) .LE. TC
               l3 = .NOT. l1
               IF ( l2 .AND. l3 ) THEN
                  t_sea_level(i,j) = TC
               ELSE
                  t_sea_level(i,j) = TC - 0.005*(t_surf(i,j)-TC)**2
               END IF
            END DO
         END DO
      END IF



      DO j = 1 , ny
      DO i = 1 , nx

         z_half_lowest=z(i,j,1)
         sea_level_pressure(i,j) = p(i,j,1) *              &
                               EXP((2.*g*z_half_lowest)/   &
                                   (R*(t_sea_level(i,j)+t_surf(i,j))))
      END DO
      END DO

    if (debug) then
      print *,'sea pres input at weird location i=20,j=1,k=1'
      print *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
      print *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
      print *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
      print *,'slp=',sea_level_pressure(20,1),     &
               sea_level_pressure(20,2),sea_level_pressure(20,3)
    endif






  end subroutine compute_seaprs




  subroutine rotate_wind (u,v,d1,d2,d3,var,                 &
                          map_proj,cen_lon,xlat,xlon,       &
                          truelat1,truelat2,data_out)

  implicit none

  integer, intent(in)            ::  d1, d2, d3

  real, dimension(d1,d2,d3), intent(out)      :: data_out
  integer, intent(in)                        :: map_proj
  integer                        ::i,j,k
  real, intent(in)                           :: cen_lon, truelat1, truelat2
  real                          :: cone
  real, dimension(d1,d2,d3), intent(in)      :: u,v 
  real, dimension(d1,d2), intent(in)         :: xlat, xlon
  real, dimension(d1,d2)         :: diff, alpha

  character (len=10), intent(in) :: var

   REAL    , PARAMETER           :: pii = 3.14159265
   REAL    , PARAMETER           :: radians_per_degree = pii/180.




       cone = 1.                                          
       if( map_proj .eq. 1) then                          
         IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
            cone=(ALOG(COS(truelat1*radians_per_degree))-            &
                  ALOG(COS(truelat2*radians_per_degree))) /          &
            (ALOG(TAN((90.-ABS(truelat1))*radians_per_degree*0.5 ))- &
             ALOG(TAN((90.-ABS(truelat2))*radians_per_degree*0.5 )) )
         ELSE
            cone = SIN(ABS(truelat1)*radians_per_degree )  
         ENDIF
       end if


       diff = xlon - cen_lon

       do i = 1, d1
       do j = 1, d2
         if(diff(i,j) .gt. 180.) then
           diff(i,j) = diff(i,j) - 360.
         end if
         if(diff(i,j) .lt. -180.) then
           diff(i,j) = diff(i,j) + 360.
         end if
       end do
       end do

       do i = 1, d1
       do j = 1, d2
          if(xlat(i,j) .lt. 0.) then
            alpha(i,j) = - diff(i,j) * cone * radians_per_degree
          else
            alpha(i,j) = diff(i,j) * cone * radians_per_degree
          end if
       end do
       end do


       if(var(1:1) .eq. "U") then
         do k=1,d3
           data_out(:,:,k) = v(:,:,k)*sin(alpha) + u(:,:,k)*cos(alpha)
         end do
       else if(var(1:1) .eq. "V") then    
         do k=1,d3
           data_out(:,:,k) = v(:,:,k)*cos(alpha) - u(:,:,k)*sin(alpha)
         end do
       end if


  end subroutine rotate_wind



  subroutine handle_err(rmarker,nf_status)

!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at 1 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with 1 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)
      integer, intent(in) :: nf_status
      character*(*), intent(in)        :: rmarker
      if (nf_status .ne. nf_noerr) then
         write(*,*)  'NetCDF error : ',rmarker
         write(*,*)  '  ',nf_strerror(nf_status)
         stop 
      endif
      
  end subroutine handle_err


  subroutine get_dimensions(infile,nx,ny,nz)
    
!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at 1 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with 1 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)
     character (len=512), intent(in)                          :: infile
     integer                                                  :: ncid, dimid, nf_status
     integer, intent(inout)                                   :: nx, ny, nz
     integer                                                  :: nlgen


      nf_status = nf_open (infile, nf_nowrite, ncid)
      call handle_err('Error opening file: '//trim(infile),nf_status)

        nf_status = nf_inq_dimid (ncid, 'west_east_stag', dimid)
        call handle_err('west_east_stag',nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, nx)
        call handle_err('Get NX',nf_status)
        nx = nx-1

        nf_status = nf_inq_dimid (ncid, 'south_north_stag', dimid)
        call handle_err('south_north_stag',nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, ny)
        call handle_err('Get NY',nf_status)
        ny = ny-1

        nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid)
        call handle_err('bottom_top',nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, nz)
        call handle_err('Get NZ',nf_status)
        nlgen = nz

      nf_status = nf_close (ncid)
      call handle_err('Error closing file',nf_status)

  end subroutine get_dimensions

  
  subroutine get_diffs(var1, var2, diff, absdiff, sqdiff, nx, ny, nz, missing)
    
    real, intent(in), dimension(:,:,:)             :: var1, var2
    real, intent(out), dimension(:,:,:)            :: diff, absdiff, sqdiff
    integer, intent(in)                            :: nx, ny, nz
    integer                                        :: i,j,k
    real, intent(in)                               :: missing

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( var1(i,j,k) /= missing .and. var2(i,j,k) /= missing ) then
         diff(i,j,k)   = var2(i,j,k) - var1(i,j,k) 
         absdiff(i,j,k)   = abs(var2(i,j,k) - var1(i,j,k)) 
         sqdiff(i,j,k) = (var2(i,j,k) - var1(i,j,k) )*(var2(i,j,k) - var1(i,j,k) ) 
       else
         diff(i,j,k)   = missing 
         absdiff(i,j,k)   = missing 
         sqdiff(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine get_diffs 


  subroutine get_diffs_wind(u1, v1, u2,v2,diff, absdiff, sqdiff, nx, ny, nz, missing)
    
    real, intent(in), dimension(:,:,:)             :: u1, u2, v1, v2
    real, intent(out), dimension(:,:,:)            :: diff, absdiff, sqdiff
    integer, intent(in)                            :: nx, ny, nz
    integer                                        :: i,j,k
    real, intent(in)                               :: missing
    real                                           :: u, v

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( u1(i,j,k) /= missing .and. v1(i,j,k) /= missing .and. &
            u2(i,j,k) /= missing .and. v2(i,j,k) /= missing ) then
         u  = u1(i,j,k) - u2(i,j,k) 
         v  = v1(i,j,k) - v2(i,j,k) 
         diff(i,j,k)   = u + v
         absdiff(i,j,k)   = abs(u) + abs (v)                
         sqdiff(i,j,k) = u*u + v*v
       else
         diff(i,j,k)   = missing 
         absdiff(i,j,k)   = missing 
         sqdiff(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine get_diffs_wind 


  subroutine domain_average(var, avg_prof, counter, nx, ny, nz, missing,isq)

  integer, intent(in)                                 :: nx,ny,nz,isq
  real, intent(in), dimension(:,:,:)                  :: var
  integer, intent(out), dimension(:)                  :: counter
  real, intent(out), dimension(:)                     :: avg_prof
  real, intent(in)                                    :: missing

  integer                                             :: i,j,k,icount
  real                                                :: dsum, dmiss
  integer                                             :: imiss




  dmiss = -99.99
  imiss = -99
  do k = 1, nz
     icount = 0
     dsum = 0
     do j = 1, ny
     do i = 1, nx
         if ( var(i,j,k) /= missing ) then
            icount = icount + 1
            dsum = dsum + var(i,j,k)
         endif
     enddo
     enddo   
     avg_prof(k) = dmiss

     counter(k) = imiss
     if (icount /= 0 ) then
        counter(k) = icount
        if ( isq .eq. 0 ) then
          avg_prof(k) = dsum /float(icount)
        else
          avg_prof(k) = sqrt(dsum /float(icount))
        endif
      endif
  enddo

  end subroutine domain_average

  
  subroutine mask_domain_average(var, avg_prof, counter, nx, ny, nz, missing,isq)

  integer, intent(in)                                 :: nx,ny,nz,isq
  real, intent(in), dimension(nx,ny,nz)               :: var
  integer, intent(out), dimension(3,nz)               :: counter
  real, intent(out), dimension(3,nz)                  :: avg_prof
  real, intent(in)                                    :: missing

  integer                                             :: i,j,k,imask
  integer                                             :: icount(3)
  real                                                :: dsum(3)
  real                                                :: dmiss
  integer                                             :: imiss

  

  dmiss = -99.99
  imiss = -99
  do k = 1, nz
     icount = 0
     dsum = 0.0
     do j = 1, ny
        do i = 1, nx
           if ( var(i,j,k) /= missing ) then
              if ( i>=istart .and. i<=iend .and. &
                   j>=jstart .and. j<=jend ) then
                  icount(1) = icount(1) + 1
                  dsum(1) = dsum(1) + var(i,j,k)
                  if ( int(landmask(i,j)) == island ) then
                     icount(2) = icount(2) + 1
                     dsum(2) = dsum(2) + var(i,j,k)
                  else if ( int(landmask(i,j)) == iswater ) then
                     icount(3) = icount(3) + 1
                     dsum(3) = dsum(3) + var(i,j,k)
                  end if
              end if
           end if
        end do
     end do
     avg_prof(:,k) = dmiss
     counter(:,k) = imiss
     do imask = 1, 3
        if ( icount(imask) /= 0 ) then
           counter(imask,k) = icount(imask)
           if ( isq .eq. 0 ) then
             avg_prof(imask,k) = dsum(imask)/float(icount(imask))
           else
             avg_prof(imask,k) = sqrt(dsum(imask)/float(icount(imask)))
           end if
        end if
     end do
  enddo

  end subroutine mask_domain_average


  subroutine get_sum(dsum, dvar, nx, ny, nz, missing)
    
    integer, intent(in)                               :: nx, ny, nz
    real, intent(in)                                  :: missing
    real, intent(in),dimension(:,:,:)                 :: dvar
    real, intent(inout),dimension(:,:,:)              :: dsum

    integer                                           :: i,j,k

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( dvar(i,j,k) /= missing .and. dsum(i,j,k) /= missing ) then
         dsum(i,j,k)   = dsum(i,j,k) + dvar(i,j,k) 
       else
         dsum(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine get_sum

  subroutine time_average(dvar, nx, ny, nz, time_count,missing, isqr)
    
    integer, intent(in)                               :: nx, ny, nz,time_count,isqr
    real, intent(in)                                  :: missing
    real, intent(inout), dimension(:,:,:)             :: dvar 

    integer                                           :: i,j,k

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
       if ( dvar(i,j,k) /= missing ) then
         if (isqr .eq. 1 ) then
           dvar(i,j,k) = sqrt(dvar(i,j,k)/float(time_count))
         else
           dvar(i,j,k)   =  dvar(i,j,k)/float(time_count)
         endif
       else
         dvar(i,j,k) = missing
       endif
    enddo
    enddo
    enddo

  end subroutine time_average

  subroutine compute_data_3d( infile, var, length, nx, ny, nz, levels, time_idx,   &
                 vert_levels, vertical_type, missing, data_out_z, debug )

  integer, intent(in)                      :: time_idx
  integer, intent(in)                      :: nx, ny, nz, levels
  integer, intent(in)                      :: length
  real, intent(in)                         :: missing
  real, intent(in)                         :: vert_levels(:)   
  character (len=1), intent(in)            :: vertical_type
  character (len=*), intent(in)            :: var
  character (len=*), intent(in)            :: infile
  logical, intent(in)                      :: debug

  real, intent(out), dimension(:,:,:)      :: data_out_z
  real, allocatable, dimension(:,:,:)      :: data_out
  real, allocatable, dimension(:,:,:)      :: z, ph, phb
  real, allocatable, dimension(:,:,:)      :: p, pb


  

      if ( vertical_type == 'p' ) then
        allocate( p (nx, ny, nz)  )
        allocate( pb(nx, ny, nz)  )
        call da_get_var_3d_real_cdf( infile,'PB',pb,nx,ny,nz,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'P',p, nx, ny, nz, time_idx,debug )
        p = (p+pb)/100.0 
      endif
      if ( vertical_type == 'z' ) then
        allocate( z (nx, ny, nz)  )
        allocate( ph (nx, ny, nz+1)  )
        allocate( phb(nx, ny, nz+1)  )
        call da_get_var_3d_real_cdf( infile,'PHB',phb,nx,ny,nz+1,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'PH',ph, nx, ny, nz+1, time_idx,debug )
        ph = (ph+phb)/9.81
        z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))
        z = z/1000.  
      endif

      allocate ( data_out (nx, ny, nz) )

      call g_output_3d (infile, time_idx, var, length, nx, ny, nz, data_out, debug)

      if ( vertical_type == 'p' ) then
         call interp_to_z( data_out, nx, ny, nz, data_out_z, nx, ny, levels,     &
                         p, vert_levels, missing, vertical_type, debug )

      else if ( vertical_type == 'z' ) then
         call interp_to_z( data_out, nx, ny, nz, data_out_z, nx, ny, levels,     &
                         z, vert_levels, missing, vertical_type, debug )
      else
        data_out_z = data_out
      endif
      deallocate ( data_out )
      if ( vertical_type == 'p' ) then
              deallocate( p )
              deallocate( pb )
      endif
      if ( vertical_type == 'z' ) then
              deallocate( z )
              deallocate( ph )
              deallocate( phb )
      endif

  end subroutine compute_data_3d
  
  subroutine compute_wind_3d( infile, nx, ny, nz, levels, time_idx,   &
                 vert_levels, vertical_type, missing, data_out_z1, data_out_z2, debug )

  integer, intent(in)                      :: time_idx
  integer, intent(in)                      :: nx, ny, nz, levels
  real, intent(in)                         :: missing
  real, intent(in)                         :: vert_levels(:)   
  character (len=1), intent(in)            :: vertical_type
  character (len=*), intent(in)            :: infile
  logical, intent(in)                      :: debug

  real, intent(out), dimension(:,:,:)      :: data_out_z1
  real, intent(out), dimension(:,:,:)      :: data_out_z2
  real, allocatable, dimension(:,:,:)      :: data_out
  real, allocatable, dimension(:,:,:)      :: z, ph, phb
  real, allocatable, dimension(:,:,:)      :: p, pb
  character (len=1)                        :: var


  

      if ( vertical_type == 'p' ) then
        allocate( p (nx, ny, nz)  )
        allocate( pb(nx, ny, nz)  )
        call da_get_var_3d_real_cdf( infile,'PB',pb,nx,ny,nz,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'P',p, nx, ny, nz, time_idx,debug )
        p = (p+pb)/100.0 
      endif
      if ( vertical_type == 'z' ) then
        allocate( z (nx, ny, nz)  )
        allocate( ph (nx, ny, nz+1)  )
        allocate( phb(nx, ny, nz+1)  )
        call da_get_var_3d_real_cdf( infile,'PHB',phb,nx,ny,nz+1,time_idx,debug )
        call da_get_var_3d_real_cdf( infile,'PH',ph, nx, ny, nz+1, time_idx,debug )
        ph = (ph+phb)/9.81
        z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))
        z = z/1000.  
      endif

      allocate ( data_out2 (nx, ny, nz) )
      allocate ( data_out1 (nx, ny, nz) )
      var='U'
      call g_output_3d (infile, time_idx, var, 1, nx, ny, nz, data_out1, debug)
      var='V'
      call g_output_3d (infile, time_idx, var, 1, nx, ny, nz, data_out2, debug)

      if ( vertical_type == 'p' ) then
         call interp_to_z( data_out1, nx, ny, nz, data_out_z1, nx, ny, levels,     &
                         p, vert_levels, missing, vertical_type, debug )
         call interp_to_z( data_out2, nx, ny, nz, data_out_z2, nx, ny, levels,     &
                         p, vert_levels, missing, vertical_type, debug )

      else if ( vertical_type == 'z' ) then
         call interp_to_z( data_out1, nx, ny, nz, data_out_z1, nx, ny, levels,     &
                         z, vert_levels, missing, vertical_type, debug )
         call interp_to_z( data_out2, nx, ny, nz, data_out_z2, nx, ny, levels,     &
                         z, vert_levels, missing, vertical_type, debug )
      else
        data_out_z1 = data_out1
        data_out_z2 = data_out2
      endif
      deallocate ( data_out1 )
      deallocate ( data_out2 )
      if ( vertical_type == 'p' ) then
              deallocate( p )
              deallocate( pb )
      endif
      if ( vertical_type == 'z' ) then
              deallocate( z )
              deallocate( ph )
              deallocate( phb )
      endif

  end subroutine compute_wind_3d
  

end program da_verif_grid
