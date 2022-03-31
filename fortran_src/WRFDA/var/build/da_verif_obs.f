












program da_verif_obs  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   use da_verif_obs_control, only : surface_type, upr_type, gpspw_type, &
      gpsref_type, record1, record2, record3, &
      record4, record5, record6, stats_value, exp_dirs, out_dirs, nstd,nstdh, &
      rmiss, diag_unit_out, nml_unit, alpha, &
      diag_unit_in, info_unit, exp_num, end_date, file_path_string, &
      if_plot_bias, if_plot_airsret, if_plot_airep,if_plot_abias, &
      if_plot_buoy, if_plot_gpspw, if_plot_gpsref, if_plot_pilot, &
      if_plot_profiler, if_plot_polaramv, if_plot_qscat, if_plot_rmse, &
      if_plot_sound, if_plot_sonde_sfc, if_plot_synop, if_plot_surface, &
      if_plot_upr, if_plot_ships, if_plot_metar, if_plot_tamdar,interval, stdp, start_date, &
      if_plot_geoamv, stdh, num_miss, &
      wrf_file, istart, iend, jstart, jend
   use da_verif_obs_init, only : initialize_surface_type, initialize_upr_type, &
      initialize_gpspw_type, initialize_gpsref_type, da_advance_cymdh , &
      initialize_t_tab      
   use da_verif_tools, only : map_info, proj_merc, proj_ps,proj_lc,proj_latlon, &
      da_llxy_wrf,da_xyll,da_map_set
  
   implicit none

   integer      :: num_obs 
   character*20 :: obs_type, dummy_c
   
   character*5  :: stn_id               
   integer      :: n, k, kk, l, levels, dummy_i
   real         :: lat, lon, press, height, dummy           
   real         :: u_obs, u_inv, u_error, u_inc, & 
                   v_obs, v_inv, v_error, v_inc, &
                   t_obs, t_inv, t_error, t_inc, &
                   p_obs, p_inv, p_error, p_inc, &
                   q_obs, q_inv, q_error, q_inc, &
                   spd_obs, spd_inv, spd_err, spd_inc
   real         :: tpw_obs, tpw_inv, tpw_err, tpw_inc
   real         :: ref_obs, ref_inv, ref_err, ref_inc
   integer      :: u_qc, v_qc, t_qc, p_qc, q_qc, tpw_qc, spd_qc, ref_qc
   integer      :: npr, nht, ier, iexp
   character*10 :: date, new_date             
   integer      :: sdate, cdate, edate        
   logical      :: if_write, is_file
   logical, allocatable :: l_skip(:)

   character(len=512)     :: out_dir,filename
   type (surface_type)    :: surface
   type (upr_type)        :: upr, gupr  
   type (gpspw_type)      :: gpspw
   type (gpsref_type)     :: gpsref, ggpsref

   integer :: nx, ny, nz, num_date, idate
   real    :: dx, cen_lat, cen_lon, truelat1, truelat2, stand_lon
   integer :: map_proj_wrf
   logical :: l_got_info, inside

   inside = .true.

   nml_unit      = 10
   diag_unit_in  = 50
   diag_unit_out = 20
   info_unit     = 30

   exp_num   = 0
   exp_dirs = ''
   out_dirs = ''

   if_plot_rmse  = .false.
   if_plot_bias  = .false.
   if_plot_abias = .false.

   if_plot_synop     = .false.
   if_plot_sonde_sfc = .false.
   if_plot_metar     = .false.
   if_plot_ships     = .false.
   if_plot_qscat     = .false.
   if_plot_buoy      = .false.

   if_plot_sound     = .false.
   if_plot_geoamv    = .false.
   if_plot_polaramv  = .false.
   if_plot_profiler  = .false.
   if_plot_airep     = .false.
   if_plot_pilot     = .false.

   if_plot_gpspw     = .false.
   if_plot_gpsref    = .false.
   if_plot_airsret   = .false.

   if_plot_tamdar    = .false.

   file_path_string = 'wrfvar/gts_omb_oma_01'
   wrf_file = 'foo.nc'

   istart = 1
   iend   = 10000
   jstart = 1
   jend   = 10000

   

   open ( unit=nml_unit, file='namelist.plot_diag', STATUS='OLD',  &
         form='formatted' )
   read ( unit=nml_unit, nml=record1, IOSTAT=ier )
   write ( unit=*, nml = record1 )
   if ( ier /= 0 ) then
      write (*,*) 'error in reading namelist record1'
      stop
   end if

   read ( unit=nml_unit, nml=record2, iostat=ier )
   write ( unit=*, nml = record2 )
   if ( ier /= 0 ) then
      write (*,*) 'error in reading namelist record2'
      stop
   end if
   read ( unit=nml_unit, nml=record3, iostat=ier )
   write ( unit=*, nml = record3 )
   if ( ier /= 0 ) then
      write (*,*) 'error in reading namelist record3'
      stop
   end if
   read ( unit=nml_unit, nml=record4, iostat=ier )
   write ( unit=*, nml = record4 )
   if ( ier /= 0 ) then
      write (*,*) 'error in reading namelist record4'
      stop
   end if
   read ( unit=nml_unit, nml=record5, iostat=ier )
   write ( unit=*, nml = record5 )
   if ( ier /= 0 ) then
      write (*,*) 'error in reading namelist record5'
      stop
   end if
   read ( unit=nml_unit, nml=record6, iostat=ier )
   if ( ier /= 0 ) then
      write (*,*) 'error in reading namelist record6'
      
   else
      write ( unit=*, nml = record6 )
   end if
   close(nml_unit)
   call initialize_t_tab

   call get_fileinfo
   if ( l_got_info ) then
      call set_mapinfo
      istart = max(1, istart)
      iend   = min(nx, iend)
      jstart = max(1, jstart)
      jend   = min(ny, jend)
   end if

   if_plot_surface = .false.
   if (if_plot_synop .or. if_plot_metar .or. if_plot_ships .or. if_plot_buoy .or. &
        if_plot_sonde_sfc .or. if_plot_qscat   ) if_plot_surface = .true.

   if_plot_upr = .false.
   if (if_plot_sound .or. if_plot_pilot .or. if_plot_profiler   .or.    &
       if_plot_geoamv .or. if_plot_polaramv  .or. if_plot_airep .or.    &
       if_plot_airsret .or. if_plot_tamdar ) if_plot_upr= .true.

   read(start_date(1:10), fmt='(i10)')sdate
   read(end_date(1:10), fmt='(i10)')edate
   write(6,'(4a)')' Diag Starting date ', start_date, ' Ending date ', end_date
   write(6,'(a,i8,a)')' Interval between dates = ', interval, ' hours.'

   num_date = 0
   date = start_date
   do while ( date <= end_date )
      num_date = num_date + 1
      call da_advance_cymdh(date, interval, date)
   end do
   allocate(l_skip(num_date))
   l_skip(:) = .false.

   
   idate = 0  
   date = start_date
   do while ( date <= end_date )
      idate = idate + 1
      do iexp = 1, exp_num
         filename = TRIM(exp_dirs(iexp))//'/'//date//'/'//trim(file_path_string)
         inquire ( file=trim(filename), exist=is_file)
         if ( .not. is_file ) then
            l_skip(idate) = .true.
         end if
         if ( l_skip(idate) ) exit
      end do
      call da_advance_cymdh(date, interval, date)
   end do

   
   

   do iexp  =1,exp_num

      idate = 0
      date = start_date
      cdate = sdate
      call initialize_upr_type(gupr)
      call initialize_gpsref_type(ggpsref)

      do while ( cdate <= edate )         
         
         call initialize_surface_type(surface)
         call initialize_upr_type(upr)
         call initialize_gpspw_type(gpspw)
         call initialize_gpsref_type(gpsref)

         idate = idate + 1

         
         filename = TRIM(exp_dirs(iexp))//'/'//date//'/'//trim(file_path_string)

         inquire ( file=trim(filename), exist=is_file) 
         if ( l_skip(idate) .or. .not. is_file ) then
            print*, 'skipping file ', trim(filename)
            
            
            out_dir=trim(out_dirs(iexp))
            if (if_plot_surface  )  then
               call write_diag_single_level(out_dir,diag_unit_out,date,'surface_u',surface%uomb,surface%uoma)     
               call write_diag_single_level(out_dir,diag_unit_out,date,'surface_v',surface%vomb,surface%voma)     
               call write_diag_single_level(out_dir,diag_unit_out,date,'surface_t',surface%tomb,surface%toma)     
               call write_diag_single_level(out_dir,diag_unit_out,date,'surface_p',surface%pomb,surface%poma)     
               call write_diag_single_level(out_dir,diag_unit_out,date,'surface_q',surface%qomb,surface%qoma)     
            end if      

            if (if_plot_gpspw )  then
               call write_diag_single_level(out_dir,diag_unit_out,date,'gpspw_tpw',gpspw%tpwomb,gpspw%tpwoma)     
            end if

            if (if_plot_gpsref  )  then
             call write_diag_multi_level_h(out_dir,diag_unit_out,date,'gps_ref',gpsref%refomb,gpsref%refoma)
            end if

            if (if_plot_upr ) then
               call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_u',upr%uomb,upr%uoma)
               call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_v',upr%vomb,upr%voma)
               call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_t',upr%tomb,upr%toma)
               call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_q',upr%qomb,upr%qoma)
            end if
            
            call da_advance_cymdh( date, interval, new_date )
            date = new_date
            read(date(1:10), fmt='(i10)')cdate
            cycle
         end if

         open (unit=diag_unit_in, file=trim(filename), form='formatted',   &
                  status='old', iostat=ier)
1           continue

         if_write = .false.
         read(diag_unit_in,'(a20,i8)', end=2000, err = 1000)obs_type,num_obs                    
         if (index( obs_type,'synop') > 0 ) then
            if (if_plot_synop ) if_write = .true.
            goto 10

         elseif (index( obs_type,'metar') > 0 ) then 
            if (if_plot_metar ) if_write = .true.
            goto 10

         elseif (index( obs_type,'ships') > 0 )  then
            if (if_plot_ships ) if_write = .true.
            goto 10

         elseif (index( obs_type,'buoy' ) > 0 )  then
            if (if_plot_buoy ) if_write = .true.
            goto 10

         elseif (index( obs_type,'sonde_sfc') > 0 )  then
            if (if_plot_sonde_sfc ) if_write = .true.
            goto 10

         elseif (index( obs_type,'polaramv') > 0)  then
            if (if_plot_polaramv ) if_write = .true.
            goto 20

         elseif (index( obs_type,'geoamv'  ) > 0)  then
            if (if_plot_geoamv ) if_write = .true.
            goto 20

         elseif (index( obs_type,'gpspw') > 0)  then
            if ( if_plot_gpspw ) if_write = .true.
            goto 30

         elseif (index( obs_type,'sound') > 0)  then
            if (if_plot_sound ) if_write = .true.
            goto 40

         elseif (index( obs_type,'airep') > 0)  then
            if (if_plot_airep ) if_write = .true.
            goto 50

         elseif (index( obs_type,'pilot')    > 0)  then
            if (if_plot_pilot ) if_write = .true.
            goto 60

         elseif (index( obs_type,'profiler') > 0)  then
            if (if_plot_profiler ) if_write = .true.
            goto 60

         elseif (index( obs_type,'ssmir') > 0)  then
            goto 70

         elseif (index( obs_type,'ssmiT') > 0)  then
            goto 80

         elseif (index( obs_type,'satem') > 0)  then
            goto 90

         elseif (index( obs_type,'ssmt1') > 0)  then
            goto 100

         elseif (index( obs_type,'ssmt2') > 0)  then
            goto 100

         elseif (index( obs_type,'qscat') > 0)  then
            if (if_plot_qscat ) if_write = .true.
            goto 110
         elseif (index( obs_type,'gpsref' ) > 0) then
            if (if_plot_gpsref ) if_write = .true.
               goto 120
         elseif (index( obs_type,'airsr') > 0)  then
               if (if_plot_airsret ) if_write = .true.
               goto 130

         elseif (index( obs_type,'tamdar') > 0)  then
            if (if_plot_tamdar ) if_write = .true.
            goto 140

         else
         print*,' Got unknown OBS_TYPE ',obs_type(1:20),' on unit ',diag_unit_in
         stop    
         end if

10       continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                                  kk, l, stn_id, &          
                                  lat, lon, press, &       
                                  u_obs, u_inv, u_qc, u_error, u_inc, & 
                                  v_obs, v_inv, v_qc, v_error, v_inc, &
                                  t_obs, t_inv, t_qc, t_error, t_inc, &
                                  p_obs, p_inv, p_qc, p_error, p_inc, &
                                  q_obs, q_inv, q_qc, q_error, q_inc
                  if (if_write) then
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if (u_qc >=  0) call update_stats(surface%uomb, surface%uoma, u_inv, u_inc)
                        if (v_qc >=  0) call update_stats(surface%vomb, surface%voma, v_inv, v_inc)
                        if (t_qc >=  0) call update_stats(surface%tomb, surface%toma, t_inv, t_inc)
                        if (p_qc >=  0) call update_stats(surface%pomb, surface%poma, p_inv, p_inc)
                        if (q_qc >=  0) call update_stats(surface%qomb, surface%qoma, q_inv, q_inc)
                     end if
                  end if
               end do      
            end do      
         end if
         goto 1

20       continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk, l, stn_id, &          
                     lat, lon, press, &        
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc

                  if (if_write .and. press > 0 ) then
                     call get_std_pr_level(press, npr, stdp, nstd) 
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if( u_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%uomb(npr),upr%uoma(npr),u_inv,u_inc)
                           call update_stats(gupr%uomb(npr),gupr%uoma(npr),u_inv,u_inc)
                        endif
                        if( v_qc >= 0 .and. npr > 0 ) then
                          call update_stats(upr%vomb(npr),upr%voma(npr),v_inv,v_inc)
                          call update_stats(gupr%vomb(npr),gupr%voma(npr),v_inv,v_inc)
                        endif
                     end if
                  end if
               end do      
            end do      
         end if

         goto 1

30       continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk, l, stn_id, &          
                     lat, lon, dummy, &       
                     tpw_obs, tpw_inv, tpw_qc, tpw_err, tpw_inc
                  if (if_write) then
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if (tpw_qc >=  0) call update_stats(gpspw%tpwomb,gpspw%tpwoma,tpw_inv,tpw_inc)
                     end if
                  end if
               end do      
            end do      
         end if

         goto 1

40       continue      

         

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, press, &       
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc, &
                     t_obs, t_inv, t_qc, t_error, t_inc, &
                     q_obs, q_inv, q_qc, q_error, q_inc
                  if (if_write .and. press > 0 ) then
                     call get_std_pr_level(press, npr, stdp, nstd) 
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if( u_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%uomb(npr),upr%uoma(npr),u_inv,u_inc)
                           call update_stats(gupr%uomb(npr),gupr%uoma(npr),u_inv,u_inc)
                        endif
                        if( v_qc >= 0 .and. npr > 0 ) then
                          call update_stats(upr%vomb(npr),upr%voma(npr),v_inv,v_inc)
                          call update_stats(gupr%vomb(npr),gupr%voma(npr),v_inv,v_inc)
                        endif
                        if( t_qc >= 0 .and. npr > 0 )  then
                          call update_stats(upr%tomb(npr),upr%toma(npr),t_inv,t_inc)
                          call update_stats(gupr%tomb(npr),gupr%toma(npr),t_inv,t_inc)
                        endif
                        if( q_qc >= 0 .and. npr > 0 )  then
                          call update_stats(upr%qomb(npr),upr%qoma(npr),q_inv,q_inc)
                          call update_stats(gupr%qomb(npr),gupr%qoma(npr),q_inv,q_inc)
                        endif
                     end if

                  end if
                end do      
             end do      
         end if
         goto 1

50       continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)') levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, press, &       
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc, &
                     t_obs, t_inv, t_qc, t_error, t_inc    
                  if (if_write .and. press > 0 ) then
                     call get_std_pr_level(press, npr, stdp, nstd) 
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if( u_qc >= 0 .and. npr > 0 ) then
                          call update_stats(upr%uomb(npr),upr%uoma(npr),u_inv,u_inc)
                          call update_stats(gupr%uomb(npr),gupr%uoma(npr),u_inv,u_inc)
                        endif
                        if( v_qc >= 0 .and. npr > 0 ) then
                          call update_stats(upr%vomb(npr),upr%voma(npr),v_inv,v_inc)
                          call update_stats(gupr%vomb(npr),gupr%voma(npr),v_inv,v_inc)
                        endif
                        if( t_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%tomb(npr),upr%toma(npr),t_inv,t_inc)
                           call update_stats(gupr%tomb(npr),gupr%toma(npr),t_inv,t_inc)
                        endif
                     end if

                  end if
               end do      
            end do     
         end if

         goto 1

60       continue      

         
         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, press, &       
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc
                  if (if_write .and. press > 0 ) then
                     call get_std_pr_level(press, npr, stdp, nstd) 
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if( u_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%uomb(npr),upr%uoma(npr),u_inv,u_inc)
                           call update_stats(gupr%uomb(npr),gupr%uoma(npr),u_inv,u_inc)
                        endif
                        if( v_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%vomb(npr),upr%voma(npr),v_inv,v_inc)
                           call update_stats(gupr%vomb(npr),gupr%voma(npr),v_inv,v_inc)
                        endif
                     end if

                  end if
               end do      
            end do      
         end if
         goto 1

70       continue      

         if ( num_obs > 0 ) then
          do n = 1, num_obs    
            read(diag_unit_in,'(i8)')levels
            do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, dummy, &       
                     spd_obs, spd_inv, spd_qc, spd_err, spd_inc 
            end do
          end do
         end if

         goto 1

80       continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,*)dummy_c                                 
               read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,7(2f17.7,i8,2f17.7))', err= 1000)&
                  k, l, stn_id, &          
                  lat, lon, dummy, &       
                  dummy, dummy, dummy_i, dummy, dummy, &    
                  dummy, dummy, dummy_i, dummy, dummy, &    
                  dummy, dummy, dummy_i, dummy, dummy, &    
                  dummy, dummy, dummy_i, dummy, dummy, &    
                  dummy, dummy, dummy_i, dummy, dummy, &    
                  dummy, dummy, dummy_i, dummy, dummy, &    
                  dummy, dummy, dummy_i, dummy, dummy
            end do
         end if
         goto 1

90       continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)') levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, dummy, &       
                     dummy,dummy, dummy_i, dummy, dummy
               end do      
            end do     
         end if

         goto 1

100      continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)') levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, dummy, &       
                     dummy,dummy, dummy_i, dummy, dummy
               end do      
            end do      
         end if

         goto 1

110      continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)') levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk, l, stn_id, &          
                     lat, lon, press, &       
                     u_obs, u_inv, u_qc, u_error, u_inc, & 
                     v_obs, v_inv, v_qc, v_error, v_inc
                  if (if_write) then
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if (u_qc >=  0) call update_stats(surface%uomb,surface%uoma,u_inv,u_inc)
                        if (v_qc >=  0) call update_stats(surface%vomb,surface%voma,v_inv,v_inc)
                     end if
                  end if
               end do      
            end do      
         end if
         goto 1

120      continue      

   IF ( num_obs > 0 ) THEN
      DO n = 1, num_obs
         read(diag_unit_in,'(i8)') levels
         DO k = 1, levels
            read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                         kk, l, stn_id, &          
                         lat, lon, height, &       
                         ref_obs, ref_inv, ref_qc, ref_err, ref_inc
            if (if_write .and. height > 0.0) then
               call get_std_ht_level(height, nht, stdh, nstdh)
               if ( l_got_info ) call check_domain(lat, lon, inside)
               if ( inside ) then
                  if ( ref_qc >=  0) then
                     call update_stats(gpsref%refomb(nht),gpsref%refoma(nht),ref_inv,ref_inc)
                     call update_stats(ggpsref%refomb(nht),ggpsref%refoma(nht),ref_inv,ref_inc)
                  end if
               end if
            end if
         END DO      
      END DO      
   ENDIF
   go to 1


130      continue      

         if ( num_obs > 0 ) then
            do n = 1, num_obs    
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &          
                     lat, lon, press, &       
                     t_obs, t_inv, t_qc, t_error, t_inc, &
                     q_obs, q_inv, q_qc, q_error, q_inc
                  if (if_write .and. press > 0 ) then
                     call get_std_pr_level(press, npr, stdp, nstd) 
                     if ( l_got_info ) call check_domain(lat, lon, inside)
                     if ( inside ) then
                        if( t_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%tomb(npr),upr%toma(npr),t_inv,t_inc)
                           call update_stats(gupr%tomb(npr),gupr%toma(npr),t_inv,t_inc)
                        endif
                        if( q_qc >= 0 .and. npr > 0 ) then
                           call update_stats(upr%qomb(npr),upr%qoma(npr),q_inv,q_inc)
                           call update_stats(gupr%qomb(npr),gupr%qoma(npr),q_inv,q_inc)
                        endif
                     end if
                  end if
               end do      
            end do      
         end if
         goto 1

140      continue

         if ( num_obs > 0 ) then
            do n = 1, num_obs
               read(diag_unit_in,'(i8)')levels
               do k = 1, levels
                  read(diag_unit_in,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', err= 1000)&
                     kk,l, stn_id, &
                     lat, lon, press, &
                     u_obs, u_inv, u_qc, u_error, u_inc, &
                     v_obs, v_inv, v_qc, v_error, v_inc, &
                     t_obs, t_inv, t_qc, t_error, t_inc, &
                     q_obs, q_inv, q_qc, q_error, q_inc
                  if (if_write .and. press > 0 ) then
                     call get_std_pr_level(press, npr, stdp, nstd)

                   if( u_qc >=  0) then
                     call update_stats(upr%uomb(npr),upr%uoma(npr),u_inv,u_inc)
                     call update_stats(gupr%uomb(npr),gupr%uoma(npr),u_inv,u_inc)
                   endif
                   if( v_qc >=  0) then
                    call update_stats(upr%vomb(npr),upr%voma(npr),v_inv,v_inc)
                    call update_stats(gupr%vomb(npr),gupr%voma(npr),v_inv,v_inc)
                   endif
                   if( t_qc >=  0)  then
                    call update_stats(upr%tomb(npr),upr%toma(npr),t_inv,t_inc)
                    call update_stats(gupr%tomb(npr),gupr%toma(npr),t_inv,t_inc)
                   endif
                   if( q_qc >=  0)  then
                    call update_stats(upr%qomb(npr),upr%qoma(npr),q_inv,q_inc)
                    call update_stats(gupr%qomb(npr),gupr%qoma(npr),q_inv,q_inc)
                   endif

                  end if
                end do
             end do
         end if
         goto 1


         
2000     continue

         close (diag_unit_in)
         
         out_dir=trim(out_dirs(iexp))
         if (if_plot_surface  )  then
            call write_diag_single_level(out_dir,diag_unit_out,date,'surface_u',surface%uomb,surface%uoma)     
            call write_diag_single_level(out_dir,diag_unit_out,date,'surface_v',surface%vomb,surface%voma)     
            call write_diag_single_level(out_dir,diag_unit_out,date,'surface_t',surface%tomb,surface%toma)     
            call write_diag_single_level(out_dir,diag_unit_out,date,'surface_p',surface%pomb,surface%poma)     
            call write_diag_single_level(out_dir,diag_unit_out,date,'surface_q',surface%qomb,surface%qoma)     
         end if      

         if (if_plot_gpspw )  then
            call write_diag_single_level(out_dir,diag_unit_out,date,'gpspw_tpw',gpspw%tpwomb,gpspw%tpwoma)     
         end if

         if (if_plot_gpsref  )  then
          call write_diag_multi_level_h(out_dir,diag_unit_out,date,'gps_ref',gpsref%refomb,gpsref%refoma)

         end if

         if (if_plot_upr ) then
            call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_u',upr%uomb,upr%uoma)
            call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_v',upr%vomb,upr%voma)
            call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_t',upr%tomb,upr%toma)
            call write_diag_multi_level(out_dir,diag_unit_out,date,'upr_q',upr%qomb,upr%qoma)
         end if
         
         call da_advance_cymdh( date, interval, new_date )
         date = new_date
         read(date(1:10), fmt='(i10)')cdate
      end do     
       if( if_plot_upr ) then
        call write_diag_multi_level(out_dir,diag_unit_out,date,'gupr_u',gupr%uomb,gupr%uoma)
        call write_diag_multi_level(out_dir,diag_unit_out,date,'gupr_v',gupr%vomb,gupr%voma)
        call write_diag_multi_level(out_dir,diag_unit_out,date,'gupr_t',gupr%tomb,gupr%toma)
        call write_diag_multi_level(out_dir,diag_unit_out,date,'gupr_q',gupr%qomb,gupr%qoma)
       endif
       if (if_plot_gpsref  )  then
        call write_diag_multi_level_h(out_dir,diag_unit_out,date,'ggps_ref',ggpsref%refomb,ggpsref%refoma)
       endif

   end do   
   stop
  
1000 print*,' Error while reading unit ',diag_unit_in,' for experiment ',exp_dirs(iexp)   
   stop
      
contains

subroutine get_std_pr_level(prs, npr, stdp, nstd) 

   implicit none

   integer, intent(in )      :: nstd
   real,    intent(in)       :: stdp(nstd)    
   integer, intent(out)      :: npr
   real,    intent(in)       :: prs           

   real             :: pr
   integer               :: k   

   npr = num_miss  
   pr = prs/100.0
   if        ( pr >= stdp(1)    ) then
       npr = 1
       return
   else if ( pr < stdp(nstd-1) .and. pr >= stdp(nstd) ) then
      npr = nstd
      return
   else
      do k = 2,nstd - 1
         if (pr   >= stdp(k) ) then
             npr = k 
            return
         end if
      end do
   end if
     
end subroutine get_std_pr_level

subroutine get_std_ht_level(height, nht, stdh, nstdh) 

   implicit none

   integer, intent(in )      :: nstdh
   real,    intent(in)       :: stdh(nstdh)    
   integer, intent(out)      :: nht
   real,    intent(in)       :: height

   real    :: ht
   integer :: k   

   ht = height*0.001  
   if ( ht <= stdh(1)    ) then
      nht = 1
      return
   else if ( ht > stdh(nstdh-1) ) then
      nht = nstdh
      return
   else
      do k = 2,nstdh - 1
         if ( ht <= stdh(k) ) then
            nht = k 
            return
         end if
      end do
   end if
     
end subroutine get_std_ht_level

subroutine update_stats(stats_omb, stats_oma, omb, oma) 

   implicit none
   type(stats_value),   intent(inout)   :: stats_omb, stats_oma  
   real, intent (in)                    :: omb, oma

   real      :: x1, x2

   stats_omb%num  = stats_omb%num + 1
   stats_oma%num  = stats_omb%num 

   x1 = 1.0/ stats_omb%num  
   x2 =  (stats_omb%num-1)*x1

   stats_omb%bias  = x2*stats_omb%bias  + omb  * x1   
   stats_oma%bias  = x2*stats_oma%bias  + oma  * x1   

   stats_omb%abias = x2*stats_omb%abias + abs(omb) * x1   
   stats_oma%abias = x2*stats_oma%abias + abs(oma) * x1   

   stats_omb%rmse  = x2*stats_omb%rmse  + omb*omb  * x1   
   stats_oma%rmse  = x2*stats_oma%rmse  + oma*oma  * x1   
     
end subroutine update_stats 

subroutine write_diag_single_level(out_dir,ounit,ldate,obs_type,omb,oma)     

   implicit none

   integer, intent(in)            :: ounit
   character*512,intent(in)       :: out_dir          
   character*10,intent(in)        :: ldate       
   character*(*),intent(in)       :: obs_type
   type (stats_value),intent(in)  :: omb
   type (stats_value),intent(in)  :: oma
 
   character*512                  :: filename         
   integer                        :: ounit1, ounit2
   real                           :: sigt,bar


   ounit1 = ounit
   ounit2 = ounit + 1

   filename = trim(out_dir)//'/'//trim(obs_type)//'_omb.diag'
   open (ounit1, file = trim(filename), form='formatted',status='unknown',position='append')                         
   filename = trim(out_dir)//'/'//trim(obs_type)//'_oma.diag'
   open (ounit2, file = trim(filename), form='formatted',status='unknown',position='append')                         
   if ( omb%num <= 1 ) then
     sigt=1.       ; bar =  rmiss
     write(ounit1,'(1x,a10,1x,6(f6.2,1x))') ldate,rmiss, rmiss, rmiss, rmiss,bar,sigt
     write(ounit2,'(1x,a10,1x,6(f6.2,1x))') ldate,rmiss, rmiss, rmiss, rmiss,bar,sigt
   else
      
     if (index(obs_type,'_q') > 0 ) then
     call sig_test(omb%num, omb%bias, omb%rmse, sigt,bar)
     bar=bar*1000.0
     write(ounit1,'(1x,a10,1x,i9,1x,5(f6.2,1x))') ldate,omb%num, 1000.0*omb%bias, 1000.0*omb%abias, 1000.0*sqrt(omb%rmse),bar,sigt
     call sig_test(oma%num, oma%bias, oma%rmse, sigt,bar)
     bar=bar*1000.0
     write(ounit2,'(1x,a10,1x,i9,1x,5(f6.2,1x))') ldate,oma%num, 1000.0*oma%bias, 1000.0*oma%abias, 1000.0*sqrt(oma%rmse),bar,sigt
     else if( index(obs_type,'_p') > 0 ) then
     call sig_test(omb%num, omb%bias, omb%rmse, sigt,bar)
     bar=bar/100.0
     write(ounit1,'(1x,a10,1x,i9,1x,5(f6.2,1x))')ldate,omb%num, omb%bias/100.0, omb%abias/100.0, sqrt(omb%rmse)/100.0,bar,sigt
     call sig_test(oma%num, oma%bias, oma%rmse, sigt,bar)
     bar=bar/100.0
     write(ounit2,'(1x,a10,1x,i9,5(1x,f6.2))') ldate,oma%num, oma%bias/100.0, oma%abias/100.0, sqrt(oma%rmse)/100.0,bar,sigt
     else
     call sig_test(omb%num, omb%bias, omb%rmse, sigt,bar)
     write(ounit1,'(1x,a10,1x,i9,1x,5(f6.2,1x))') ldate,omb%num, omb%bias, omb%abias, sqrt(omb%rmse),bar,sigt
     call sig_test(oma%num, oma%bias, oma%rmse, sigt,bar)
     write(ounit2,'(1x,a10,1x,i9,5(1x,f6.2))') ldate,oma%num, oma%bias, oma%abias, sqrt(oma%rmse),bar,sigt
     endif

   end if
   close(ounit1)
   close(ounit2)
     
end subroutine write_diag_single_level     

subroutine write_diag_multi_level(out_dir,ounit,ldate,obs_type,omb,oma)     

   implicit none

   integer, intent(in)            :: ounit
   character*512,intent(in)       :: out_dir         
   character*10,intent(in)        :: ldate       
   character*(*),intent(in)       :: obs_type
   type (stats_value),intent(in)  :: omb(nstd)
   type (stats_value),intent(in)  :: oma(nstd)
 
   character*512                  :: filename         
   integer                        :: k
   real                           :: xnum(nstd)
   integer                        :: num(nstd)
   real, dimension(nstd)          :: rmse, bias, abias,sigt,bar
   integer                        :: ounit1, ounit2

   ounit1 = ounit
   ounit2 = ounit + 1

   filename = trim(out_dir)//'/'//trim(obs_type)//'_omb.diag'
   open (ounit1, file = trim(filename), form='formatted',status='unknown',position='append')                         
   filename = trim(out_dir)//'/'//trim(obs_type)//'_oma.diag'
   open (ounit2, file = trim(filename), form='formatted',status='unknown',position='append')                         

   do k = 1, nstd
      num(k) = omb(k)%num
      if (num(k) <= 1 ) then
         xnum(k)  = rmiss     
         rmse(k)  = rmiss       
         bias(k)  = rmiss       
         abias(k) = rmiss                 
         bar(k)   = rmiss
         sigt(k)  = 1.0
      else
         if (index(obs_type,'_q') > 0 ) then
            rmse(k) = sqrt(omb(k)%rmse) * 1000
            bias(k) = omb(k)%bias * 1000
            abias(k) = omb(k)%abias * 1000
            call sig_test(num(k), omb(k)%bias, omb(k)%rmse, sigt(k),bar(k))
            bar(k) = bar(k)*1000.
         else
            rmse(k) = sqrt(omb(k)%rmse)
            bias(k) = omb(k)%bias
            abias(k) = omb(k)%abias
            call sig_test(num(k), omb(k)%bias, omb(k)%rmse, sigt(k),bar(k))
         end if
      xnum(k) = num(k)
      end if
   end do

    write(ounit1,'(1x,a10,1x,16(6(1x,f12.2)))')ldate, (xnum(k), bias(k), abias(k),&
         rmse(k),bar(k),sigt(k),k=1,nstd)

   do k = 1, nstd   
      num(k) = oma(k)%num
      if( num(k) <= 1 ) then
         xnum(k)  = rmiss     
         rmse(k)  = rmiss       
         bias(k)  = rmiss       
         abias(k) = rmiss                 
      else
         if (index(obs_type,'_q') > 0 ) then
            rmse(k) = sqrt(oma(k)%rmse) * 1000
            bias(k) = oma(k)%bias * 1000
            abias(k) = oma(k)%abias * 1000
            call sig_test(num(k), oma(k)%bias, oma(k)%rmse, sigt(k),bar(k))
            bar(k) = bar(k)*1000.
         else
            rmse(k) = sqrt(oma(k)%rmse)
            bias(k) = oma(k)%bias
            abias(k) = oma(k)%abias
            call sig_test(num(k), oma(k)%bias, oma(k)%rmse, sigt(k),bar(k))
         end if
       xnum(k) = num(k)
      end if
   end do
   write(ounit2,'(1x,a10,1x,16(6(1x,f12.2)))')ldate, (xnum(k), bias(k), abias(k),&
             rmse(k),bar(k),sigt(k),k=1,nstd)
       
   close(ounit1)
   close(ounit2)
     
end subroutine write_diag_multi_level     
     subroutine write_diag_multi_level_h(out_dir,ounit,date,obs_type,omb,oma)
     implicit none
     integer, intent(in)            :: ounit
     character*512,intent(in)       :: out_dir
     character*10,intent(in)        :: date
     character*(*),intent(in)       :: obs_type
     type (stats_value),intent(in)  :: omb(nstdh)
     type (stats_value),intent(in)  :: oma(nstdh)

     character*512                  :: filename
     integer                        :: k
     real                           :: xnum(nstdh)
     integer                        :: num(nstdh)
     real, dimension(nstdh)         :: rmse, bias, abias, sigt, bar
  
     integer                        :: ounit1, ounit2

     ounit1 = ounit
     ounit2 = ounit + 1


     filename = trim(out_dir)//'/'//trim(obs_type)//'_omb.diag'
     open (ounit1, file = trim(filename), form='formatted',status='unknown',position='append')
     filename = trim(out_dir)//'/'//trim(obs_type)//'_oma.diag'
     open (ounit2, file = trim(filename), form='formatted',status='unknown',position='append')

     do k = 1, nstdh
     num(k) = omb(k)%num
     if( num(k) <= 1 ) then
     xnum(k)  = rmiss
     rmse(k)  = rmiss
     bias(k)  = rmiss
     abias(k) = rmiss
     bar(k)   = rmiss
     sigt(k)  = 1.0
     else
        if( index(obs_type,'_q') > 0 ) then

         rmse(k) = sqrt(omb(k)%rmse) * 1000
         bias(k) = omb(k)%bias * 1000
         abias(k) = omb(k)%abias * 1000
         call sig_test(num(k), omb(k)%bias, omb(k)%rmse, sigt(k),bar(k))
         bar(k) = bar(k)*1000.
       else
        rmse(k) = sqrt(omb(k)%rmse)
        bias(k) = omb(k)%bias
        abias(k) = omb(k)%abias
        call sig_test(num(k), omb(k)%bias, omb(k)%rmse, sigt(k),bar(k))
       endif
      xnum(k) = num(k)
     endif
     enddo
     write(ounit1,'(1x,a10,1x,150(6(1x,f12.2)))')date, (xnum(k), bias(k), abias(k), &
           rmse(k),bar(k),sigt(k), k=1,nstdh)

     do k = 1, nstdh
     num(k) = oma(k)%num
     if( num(k) <= 1 ) then
     xnum(k)  = rmiss
     rmse(k)  = rmiss
     bias(k)  = rmiss
     abias(k) = rmiss
     bar(k)   = rmiss
     sigt(k)  = 1.0
     else
        if( index(obs_type,'_q') > 0 ) then

         rmse(k) = sqrt(oma(k)%rmse) * 1000
         bias(k) = oma(k)%bias * 1000
         abias(k) = oma(k)%abias * 1000
         call sig_test(num(k), oma(k)%bias, oma(k)%rmse, sigt(k),bar(k))
         bar(k) = bar(k)*1000.
       else
        rmse(k) = sqrt(oma(k)%rmse)
        bias(k) = oma(k)%bias
        abias(k) = oma(k)%abias
        call sig_test(num(k), oma(k)%bias, oma(k)%rmse, sigt(k),bar(k))
       endif
      xnum(k) = num(k)
     endif
     enddo
     write(ounit2,'(1x,a10,1x,150(6(1x,f12.2)))')date, (xnum(k), bias(k), abias(k), &                  
           rmse(k),bar(k),sigt(k), k=1,nstdh)


     close(ounit1)
     close(ounit2)

     end subroutine write_diag_multi_level_h

     subroutine sig_test(num, bias, rmse, sigt,bar)
     implicit none
     integer, intent(in)      :: num
     real,    intent(in)      :: bias, rmse
     real,    intent(out)     :: sigt, bar

     real                     :: t_val, sd, tmp

     sigt=0.
     tmp = num/real(num-1)

     sd = sqrt ( tmp*( rmse - bias*bias )  )
     do k=2,34
       if (real(num-1) < alpha(k,1)) exit
     end do

     t_val = bias*sqrt( real(num) ) /sd
     bar = alpha(k-1,2) * sd /sqrt( real(num) )

     if (abs(t_val) >= alpha(k-1,2)) sigt=1.

    end subroutine sig_test


   subroutine get_fileinfo
      implicit none
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
      integer :: iost(12)
      integer :: ncid
      l_got_info = .false.
      iost(1) = nf_open(trim(wrf_file), NF_NOWRITE, ncid)
      if ( iost(1) /= NF_NOERR ) then
         print*, 'INFO: wrf_file: '//trim(wrf_file)//' does not exist for retrieving mapping info'
         return
      else
         print*, 'Retrieving mapping info from wrf_file: ',trim(wrf_file)
      end if
      iost(2)  = nf_get_att_int(ncid, NF_GLOBAL, 'WEST-EAST_GRID_DIMENSION', nx)
      iost(3)  = nf_get_att_int(ncid, NF_GLOBAL, 'SOUTH-NORTH_GRID_DIMENSION', ny)
      iost(4)  = nf_get_att_int(ncid, NF_GLOBAL, 'BOTTOM-TOP_GRID_DIMENSION', nz)
      iost(5)  = nf_get_att_double(ncid, NF_GLOBAL, 'DX', dx)
      iost(6)  = nf_get_att_double(ncid, NF_GLOBAL, 'CEN_LAT', cen_lat)
      iost(7)  = nf_get_att_double(ncid, NF_GLOBAL, 'CEN_LON', cen_lon)
      iost(8)  = nf_get_att_double(ncid, NF_GLOBAL, 'TRUELAT1', truelat1)
      iost(9)  = nf_get_att_double(ncid, NF_GLOBAL, 'TRUELAT2', truelat2)
      iost(10) = nf_get_att_double(ncid, NF_GLOBAL, 'STAND_LON', stand_lon)
      iost(11) = nf_get_att_int(ncid, NF_GLOBAL, 'MAP_PROJ', map_proj_wrf)
      iost(12) = nf_close(ncid)
      if ( .not. any(iost/=NF_NOERR) ) then
         l_got_info = .true.
      end if
      print*, 'nx, ny, nz, dx, map_proj, cen_lat, cen_lon, truelat1, truelat2, stand_lon = ', &
         nx, ny, nz, dx, map_proj_wrf, cen_lat, cen_lon, truelat1, truelat2, stand_lon
   end subroutine get_fileinfo

   subroutine set_mapinfo
      implicit none
      integer :: map_proj_util
      real    :: xref, yref
      real :: start_x, start_y
      xref = nx/2.0
      yref = ny/2.0
      if ( map_proj_wrf == 0 .or. map_proj_wrf == 6 ) then
         map_proj_util = proj_latlon
      else if ( map_proj_wrf == 3 ) then
         map_proj_util = proj_merc
      else if ( map_proj_wrf == 1 ) then
         map_proj_util = proj_lc
      else if ( map_proj_wrf == 2 ) then
         map_proj_util = proj_ps
      end if
      call da_map_set(map_proj_util, cen_lat,cen_lon, xref, yref, dx, &
         stand_lon, truelat1, truelat2, truelat1, stand_lon, map_info)
      
   end subroutine set_mapinfo

   subroutine check_domain(lat, lon, inside)
      implicit none
      real,    intent(in)  :: lat, lon
      logical, intent(out) :: inside
      real                 :: xx, yy
      call da_llxy_wrf(map_info, lat, lon, xx, yy)
      inside = .false.
      if ( xx >= istart .and. xx <= iend .and.  &
           yy >= jstart .and. yy <= jend ) then
         inside = .true.
      end if
   end subroutine check_domain

end program da_verif_obs
