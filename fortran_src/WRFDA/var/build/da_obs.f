












module da_obs

   use da_define_structures, only : multi_level_type, y_type, iv_type, infa_type, &
      field_type, each_level_type,da_allocate_y, da_random_seed,da_allocate_y_rain, &
      da_allocate_y_radar
   use module_domain, only : domain, x_type

   use da_airep, only : da_transform_xtoy_airep, da_transform_xtoy_airep_adj 
   use da_airsr, only : da_transform_xtoy_airsr, da_transform_xtoy_airsr_adj 
   use da_bogus, only : da_transform_xtoy_bogus, da_transform_xtoy_bogus_adj
   use da_buoy, only : da_transform_xtoy_buoy,da_transform_xtoy_buoy_adj
   use da_control, only : use_shipsobs, use_synopobs, use_ssmt2obs, &
      use_soundobs,use_mtgirsobs,use_satemobs, use_profilerobs, use_pilotobs, &
      use_qscatobs,use_metarobs, use_polaramvobs, use_geoamvobs, &
      use_bogusobs,use_buoyobs, use_airsretobs, use_tamdarobs, trace_use, num_procs, &
      xmiss, missing_r, missing, use_airepobs,use_gpspwobs,use_gpsztdobs,use_gpsrefobs, &
      use_ssmt1obs,filtered_obs_unit,fmt_each,fmt_info,fmt_srfc, ide, jde, &
      pseudo_x, fg_format, fg_format_kma_global, fg_format_wrf_arw_regional,fg_format_wrf_nmm_regional, &
      missing_data, pseudo_var, pseudo_val,stdout, num_pseudo, pseudo_y, pseudo_z, &
      pseudo_err,obs_qc_pointer,myproc,rtm_option,rtm_option_rttov, &
      rtm_option_crtm,use_rad, base_temp, base_lapse, base_pres, &
      ob_format,ob_format_ascii,filename_len, trace_use_dull, &
      sound, mtgirs, synop, profiler, gpsref, gpspw, polaramv, geoamv, ships, metar, &
      satem, radar, ssmi_rv, ssmi_tb, ssmt1, ssmt2, airsr, pilot, airep, sonde_sfc,rain, &
      bogus, buoy, qscat, tamdar, tamdar_sfc, pseudo, num_ob_indexes, its,ite,jds,jts,jte,ids, &
      write_mod_filtered_obs, radiance, use_varbc, obs_names
   
   use da_crtm, only : da_transform_xtoy_crtm, da_transform_xtoy_crtm_adj
      
      
   use da_geoamv,    only : da_transform_xtoy_geoamv, da_transform_xtoy_geoamv_adj
   use da_gpspw,     only : da_transform_xtoy_gpspw,da_transform_xtoy_gpspw_adj, &
                            da_transform_xtoy_gpsztd,da_transform_xtoy_gpsztd_adj
   use da_gpsref,    only : da_transform_xtoy_gpsref,da_transform_xtoy_gpsref_adj
   use da_metar,     only : da_transform_xtoy_metar, da_transform_xtoy_metar_adj
   use da_physics,   only : da_tp_to_qs,da_get_q_error
   use da_pilot,     only : da_transform_xtoy_pilot,da_transform_xtoy_pilot_adj
   use da_polaramv,  only : da_transform_xtoy_polaramv, da_transform_xtoy_polaramv_adj
   use da_profiler,  only : da_transform_xtoy_profiler, da_transform_xtoy_profiler_adj
   use da_pseudo,    only : da_transform_xtoy_pseudo, da_transform_xtoy_pseudo_adj
   use da_qscat,     only : da_transform_xtoy_qscat,da_transform_xtoy_qscat_adj
   use da_radar,     only : da_transform_xtoy_radar,da_transform_xtoy_radar_adj
   use da_rain,      only : da_transform_xtoy_rain,da_transform_xtoy_rain_adj
   use da_reporting, only : da_error, message, da_warning, da_message
   use da_satem,     only : da_transform_xtoy_satem, da_transform_xtoy_satem_adj
   use da_ships,     only : da_transform_xtoy_ships, da_transform_xtoy_ships_adj
   use da_sound,     only : da_transform_xtoy_sound, da_transform_xtoy_sonde_sfc, &
      da_transform_xtoy_sound_adj, da_transform_xtoy_sonde_sfc_adj
   use da_mtgirs,    only : da_transform_xtoy_mtgirs, da_transform_xtoy_mtgirs_adj
  use da_tamdar,    only : da_transform_xtoy_tamdar, da_transform_xtoy_tamdar_adj, &
                            da_transform_xtoy_tamdar_sfc, da_transform_xtoy_tamdar_sfc_adj
   use da_ssmi,      only : da_transform_xtoy_ssmt1, da_transform_xtoy_ssmt2, &
      da_transform_xtoy_ssmi_tb, da_transform_xtoy_ssmi_rv, &
      da_transform_xtoy_ssmi_tb_adj, da_transform_xtoy_ssmi_rv_adj, &
      da_transform_xtoy_ssmt1_adj, da_transform_xtoy_ssmt2_adj
   use da_synop,     only : da_transform_xtoy_synop,da_transform_xtoy_synop_adj
   use da_tools_serial,    only : da_free_unit, da_get_unit
   use da_tools,     only : da_add_noise, da_add_noise_new,da_random_omb, &
                            da_geo2msl1, da_msl2geo1
   use da_tracing,   only : da_trace_entry, da_trace_exit 
   use module_dm,    only : wrf_dm_sum_real, wrf_dm_sum_reals

   implicit none

contains

subroutine da_obs_proc_station(obs,fm,uvq_direct)

   !-----------------------------------------------------------------------
   ! Purpose: Processing obs data read in at a station: 
   !
   ! (1) Surface correction
   ! (2) Vertical checks
   ! (3) Missing data check
   !-----------------------------------------------------------------------

   implicit none
                             
   type (multi_level_type), intent(inout) :: obs

   logical,optional :: uvq_direct
   real    :: po, to, rho, es, qs, qvo
   integer :: i,fm

   if (trace_use_dull) call da_trace_entry("da_obs_proc_station")

   do i = 1, obs % info % levels
      !-----------------------------------------------------------------------
      ! Gross check for t, p & q
      !-----------------------------------------------------------------------

      if (obs%each(i)%t %inv <  75.0 .or. obs%each(i)%t %inv > 350.0  ) then
         obs%each(i)%t %inv = missing_r
         obs%each(i)%t %qc  = missing
      end if
      if (obs%each(i)%p %inv <=  0.0 .or. obs%each(i)%p %inv > 110000.0) then
         obs%each(i)%p %inv = missing_r
         obs%each(i)%p %qc  = missing
      end if
      if (obs%each(i)%rh%inv < 0.0)  then
         obs%each(i)%rh%inv = missing_r
         obs%each(i)%rh%qc  = missing
      end if
     
     if(.not. (present(uvq_direct) .and. uvq_direct) .and. fm.ne.161) then
      po   = obs % each(i) % p  % inv
      to   = obs % each(i) % t  % inv
      rho  = obs % each(i) % rh  % inv

      if (ob_format == ob_format_ascii) then ! Calculate q if possible:
         if (abs(po -  missing_r) > 1.0 .and. abs(to -  missing_r) > 1.0 .and. &
             abs(rho-  missing_r) > 1.0) then

            call da_tp_to_qs(to, po, es, qs)

            if (rho > 100.0) then
               qvo = qs
            else
               qvo  = rho * qs / 100.0
            end if
         else
            qvo       = missing_r
         end if

         obs % each(i) % q  % inv = qvo
         obs % each(i) % q  % qc = obs % each(i) % rh % qc
         obs % each(i) % q  % error = obs % each(i) % rh % error
      end if
    else
         obs % each(i) % q  % inv = obs % each(i) % rh  % inv
         obs % each(i) % q  % qc = obs % each(i) % rh % qc
         obs % each(i) % q  % error = obs % each(i) % rh % error
    end if
   end do
   if (trace_use_dull) call da_trace_exit("da_obs_proc_station")

end subroutine da_obs_proc_station


subroutine da_transform_xtoy(cv_size, cv, grid, iv, y)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------------

   implicit none
   
   integer, intent(in)           :: cv_size         ! Size of cv array.
   real, intent(in)              :: cv(1:cv_size)   ! control variables.
   type (domain),  intent(inout) :: grid
   type (iv_type), intent(inout) :: iv       ! obs. increment vector (o-b).
   type (y_type),  intent(inout) :: y        ! y = h (grid%xa)

   if (trace_use) call da_trace_entry("da_transform_xtoy")
   
   !--------------------------------------------------------------------------
   ! [1.0] observation operator y = H(x):
   !--------------------------------------------------------------------------
  
   if (iv%info(sound)%nlocal          > 0) call da_transform_xtoy_sound    (grid, iv, y)
   if (iv%info(sonde_sfc)%nlocal      > 0) call da_transform_xtoy_sonde_sfc(grid, iv, y)
   if (iv%info(mtgirs)%nlocal         > 0) call da_transform_xtoy_mtgirs   (grid, iv, y)
   if (iv%info(tamdar)%nlocal         > 0) call da_transform_xtoy_tamdar   (grid, iv, y)
   if (iv%info(tamdar_sfc)%nlocal     > 0) call da_transform_xtoy_tamdar_sfc(grid, iv, y)
   if (iv%info(synop)%nlocal          > 0) call da_transform_xtoy_synop    (grid, iv, y)
   if (iv%info(geoamv)%nlocal         > 0) call da_transform_xtoy_geoamv   (grid, iv, y)
   if (iv%info(polaramv)%nlocal       > 0) call da_transform_xtoy_polaramv (grid, iv, y)
   if (iv%info(airep)%nlocal          > 0) call da_transform_xtoy_airep    (grid, iv, y)
   if (iv%info(metar)%nlocal          > 0) call da_transform_xtoy_metar    (grid, iv, y)
   if (iv%info(ships)%nlocal          > 0) call da_transform_xtoy_ships    (grid, iv, y)
   if (iv%info(gpspw)%nlocal          > 0) then
      if (use_gpspwobs) then
         call da_transform_xtoy_gpspw    (grid, iv, y)
      else if (use_gpsztdobs) then
         call da_transform_xtoy_gpsztd   (grid, iv, y)
      endif
   end if
   if (iv%info(ssmi_tb)%nlocal        > 0) call da_transform_xtoy_ssmi_tb  (grid, iv, y)
   if (iv%info(ssmi_rv)%nlocal        > 0) call da_transform_xtoy_ssmi_rv  (grid, iv, y)
   if (iv%info(pilot)%nlocal          > 0) call da_transform_xtoy_pilot    (grid, iv, y)
   if (iv%info(satem)%nlocal          > 0) call da_transform_xtoy_satem    (grid, iv, y)
   if (iv%info(ssmt1)%nlocal          > 0) call da_transform_xtoy_ssmt1    (grid, iv, y)
   if (iv%info(ssmt2)%nlocal          > 0) call da_transform_xtoy_ssmt2    (grid, iv, y)
   if (iv%info(qscat)%nlocal          > 0) call da_transform_xtoy_qscat    (grid, iv, y)
   if (iv%info(profiler)%nlocal       > 0) call da_transform_xtoy_profiler (grid, iv, y)
   if (iv%info(buoy)%nlocal           > 0) call da_transform_xtoy_buoy     (grid, iv, y)
   if (iv%info(gpsref)%nlocal         > 0) call da_transform_xtoy_gpsref   (grid, iv, y)
   if (iv%info(radar)%nlocal          > 0) call da_transform_xtoy_radar    (grid, iv, y)
   if (iv%info(bogus)%nlocal          > 0) call da_transform_xtoy_bogus    (grid, iv, y)
   if (iv%info(airsr)%nlocal          > 0) call da_transform_xtoy_airsr    (grid, iv, y)
   if (iv%info(pseudo)%nlocal         > 0) call da_transform_xtoy_pseudo   (grid, iv, y)

   if (use_rad) then
      if (rtm_option == rtm_option_rttov) then
      elseif (rtm_option == rtm_option_crtm) then
         !if (use_crtm_kmatrix) then
         !   call da_transform_xtoy_crtmk (grid, iv, y)
         !else if (use_crtm_kmatrix_fast) then
         !   call da_transform_xtoy_crtmk_f (grid, iv, y)
         !else
            call da_transform_xtoy_crtm (cv_size, cv, grid, iv, y)
         !end if
       else
          call da_warning("da_transform_xtoy.inc",70,(/"Unknown radiative transfer model"/))
       end if
   end if

   if (trace_use) call da_trace_exit("da_transform_xtoy")

end subroutine da_transform_xtoy


subroutine da_transform_xtoy_adj(cv_size, cv, grid, iv, jo_grad_y, jo_grad_x)

   !--------------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !--------------------------------------------------------------------------
   
   implicit none
   
   integer, intent(in)           :: cv_size         ! Size of cv array.
   real, intent(inout)           :: cv(1:cv_size)   ! control variables.
   type (domain),  intent(inout) :: grid
   type (iv_type), intent(inout) :: iv          ! obs. inc vector (o-b).
   type (y_type),  intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type),  intent(inout) :: jo_grad_x   ! grad_x(jo)

   if (trace_use) call da_trace_entry("da_transform_xtoy_adj")
  
   !--------------------------------------------------------------------------
   ! [1.0] observation operator y = H(x):
   !--------------------------------------------------------------------------
  
   if (iv%info(sound)%nlocal    > 0) call da_transform_xtoy_sound_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(sonde_sfc)%nlocal  > 0) call da_transform_xtoy_sonde_sfc_adj(grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(mtgirs)%nlocal   > 0) call da_transform_xtoy_mtgirs_adj   (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(tamdar)%nlocal   > 0) call da_transform_xtoy_tamdar_adj   (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(tamdar_sfc)%nlocal > 0) call da_transform_xtoy_tamdar_sfc_adj(grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(synop)%nlocal    > 0) call da_transform_xtoy_synop_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(geoamv)%nlocal   > 0) call da_transform_xtoy_geoamv_adj   (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(polaramv)%nlocal > 0) call da_transform_xtoy_polaramv_adj (grid, iv, jo_grad_y, jo_grad_x)   
   if (iv%info(airep)%nlocal    > 0) call da_transform_xtoy_airep_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(metar)%nlocal    > 0) call da_transform_xtoy_metar_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(ships)%nlocal    > 0) call da_transform_xtoy_ships_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(gpspw)%nlocal    > 0) then
      if (use_gpspwobs) then
         call da_transform_xtoy_gpspw_adj    (grid, iv, jo_grad_y, jo_grad_x)
      else if (use_gpsztdobs) then
         call da_transform_xtoy_gpsztd_adj   (grid, iv, jo_grad_y, jo_grad_x)
      endif
   end if
   if (iv%info(ssmi_tb)%nlocal  > 0) call da_transform_xtoy_ssmi_tb_adj  (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(ssmi_rv)%nlocal  > 0) call da_transform_xtoy_ssmi_rv_adj  (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(pilot)%nlocal    > 0) call da_transform_xtoy_pilot_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(satem)%nlocal    > 0) call da_transform_xtoy_satem_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(ssmt1)%nlocal    > 0) call da_transform_xtoy_ssmt1_adj    (iv, jo_grad_y, jo_grad_x)
   if (iv%info(ssmt2)%nlocal    > 0) call da_transform_xtoy_ssmt2_adj    (iv, jo_grad_y, jo_grad_x)
   if (iv%info(qscat)%nlocal    > 0) call da_transform_xtoy_qscat_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(profiler)%nlocal > 0) call da_transform_xtoy_profiler_adj (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(buoy)%nlocal     > 0) call da_transform_xtoy_buoy_adj     (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(gpsref)%nlocal   > 0) call da_transform_xtoy_gpsref_adj   (iv, jo_grad_y, jo_grad_x)
   if (iv%info(radar)%nlocal    > 0) call da_transform_xtoy_radar_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(bogus)%nlocal    > 0) call da_transform_xtoy_bogus_adj    (grid, iv, jo_grad_y, jo_grad_x)
   if (iv%info(airsr)%nlocal    > 0) call da_transform_xtoy_airsr_adj    (iv, jo_grad_y, jo_grad_x)
   if (iv%info(pseudo)%nlocal   > 0) call da_transform_xtoy_pseudo_adj   (iv, jo_grad_y, jo_grad_x)

   if (use_rad) then
      if (rtm_option == rtm_option_rttov) then
      elseif (rtm_option == rtm_option_crtm) then
         call da_transform_xtoy_crtm_adj (cv_size, cv, iv, jo_grad_y, jo_grad_x)
      else
         call da_warning("da_transform_xtoy_adj.inc",68,(/"Unknown radiative transfer model"/))
      end if
   end if

   if (trace_use) call da_trace_exit("da_transform_xtoy_adj")

end subroutine da_transform_xtoy_adj


subroutine da_add_noise_to_ob( iv, ob )
!----------------------------------------------------------------------------   
!  History:
!
!  Additions:
!             07/08/2003  -   Profiler and Buoy Obs         Syed RH Rizvi    
!             03/08/2006      Add radiance part               Zhiquan Liu
!             06/23/2006  -   MPI update                    Syed RH Rizvi    
!             07/03/2006  -   update for AIRS retrievals    Syed RH Rizvi    
!
!  Purpose: Allocates observation structure and fills it fro iv.
!---------------------------------------------------------------------------- 
  
   implicit none


   type (iv_type), intent(inout) :: iv   ! Obs and header structure.
   type (y_type), intent(inout)  :: ob   ! (Smaller) observation structure.

   real                          :: z1, z2, z3, z4, z5, z6, z7, dum ! Random numbers.
   integer                       :: n, k, i     ! Loop counters.
   integer                       :: ounit     ! Output unit
   integer                       :: num_obs, ios                 
   character(len=20)             :: ob_name, filename

   if (trace_use_dull) call da_trace_entry("da_add_noise_to_ob")

!----------------------------------------------------------------------------
!  Fix output unit
!----------------------------------------------------------------------------
   call da_get_unit(ounit)

      dum = -999999.9
!----------------------------------------------------------------------
!  [1.0] Initiate random number sequence:
!----------------------------------------------------------------------

   call da_random_seed
   
!----------------------------------------------------------------------
!  [2.0] Create noise and output:
!----------------------------------------------------------------------
      write(unit=filename, fmt='(a,i4.4)') 'rand_obs_error.', myproc

   open(unit=ounit,file=trim(filename),form='formatted',iostat=ios)
   if (ios /= 0 ) then
      call da_error("da_add_noise_to_ob.inc",51, &
         (/"Cannot open random observation error file"//filename/))
   Endif

!  [2.1] Transfer surface obs:

   if ( iv%info(synop)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(synop)%nlocal
       if(iv%info(synop)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'synop', num_obs   
      num_obs = 0 

      do n = 1, iv%info(synop)%nlocal
       if(iv%info(synop)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')  1
!        Add random perturbation:
         call da_add_noise( iv % synop(n) % u, ob % synop(n) % u, z1 )
         call da_add_noise( iv % synop(n) % v, ob % synop(n) % v, z2 )
         call da_add_noise( iv % synop(n) % t, ob % synop(n) % t, z3 )
         call da_add_noise( iv % synop(n) % p, ob % synop(n) % p, z4 )
         call da_add_noise( iv % synop(n) % q, ob % synop(n) % q, z5 )

!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, iv % synop(n) % u % error, z1, &
                                  iv % synop(n) % v % error, z2, &
                                  iv % synop(n) % t % error, z3, &
                                  iv % synop(n) % p % error, z4, &
                                  iv % synop(n) % q % error, z5
       end if
      end do
   end if

!  [2.2] Transfer metar obs:

   if ( iv%info(metar)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(metar)%nlocal
       if(iv%info(metar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'metar', num_obs   
      num_obs = 0 
      do n = 1, iv%info(metar)%nlocal
       if(iv%info(metar)%proc_domain(1,n)) then
         num_obs = num_obs + 1 
         write(ounit,'(i8)')  1
!        Add random perturbation:
         call da_add_noise( iv % metar(n) % u, ob % metar(n) % u, z1 )
         call da_add_noise( iv % metar(n) % v, ob % metar(n) % v, z2 )
         call da_add_noise( iv % metar(n) % t, ob % metar(n) % t, z3 )
         call da_add_noise( iv % metar(n) % p, ob % metar(n) % p, z4 )
         call da_add_noise( iv % metar(n) % q, ob % metar(n) % q, z5 )

!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, &
                                  iv % metar(n) % u % error, z1, &
                                  iv % metar(n) % v % error, z2, &
                                  iv % metar(n) % t % error, z3, &
                                  iv % metar(n) % p % error, z4, &
                                  iv % metar(n) % q % error, z5
       end if
      end do
   end if

!  [2.3] Transfer ships obs:

   if ( iv%info(ships)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(ships)%nlocal
       if(iv%info(ships)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'ships', num_obs   
      num_obs = 0 
      do n = 1, iv%info(ships)%nlocal 
       if(iv%info(ships)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')  1
!        Add random perturbation:
         call da_add_noise( iv % ships(n) % u, ob % ships(n) % u, z1 )
         call da_add_noise( iv % ships(n) % v, ob % ships(n) % v, z2 )
         call da_add_noise( iv % ships(n) % t, ob % ships(n) % t, z3 )
         call da_add_noise( iv % ships(n) % p, ob % ships(n) % p, z4 )
         call da_add_noise( iv % ships(n) % q, ob % ships(n) % q, z5 )
!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, &
                                  iv % ships(n) % u % error, z1, &
                                  iv % ships(n) % v % error, z2, &
                                  iv % ships(n) % t % error, z3, &
                                  iv % ships(n) % p % error, z4, &
                                  iv % ships(n) % q % error, z5
       end if
      end do
   end if


!  [2.4.1] Transfer Geostationary AMVs obs:

   if ( iv%info(geoamv)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(geoamv)%nlocal
       if(iv%info(geoamv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'geoamv', num_obs   
      num_obs = 0 
      do n = 1, iv%info(geoamv)%nlocal
       if(iv%info(geoamv)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(geoamv)%levels(n)
         do k = 1, iv%info(geoamv)%levels(n)
!        Add random perturbation:
            call da_add_noise( iv % geoamv(n) % u(k), ob % geoamv(n) % u(k), z1)
            call da_add_noise( iv % geoamv(n) % v(k), ob % geoamv(n) % v(k), z2)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % geoamv(n) % u(k) % error, z1, &
                               iv % geoamv(n) % v(k) % error, z2, &
                               dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  [2.4.2] Transfer Polar AMVs obs:

   if ( iv%info(polaramv)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(polaramv)%nlocal
       if (iv%info(polaramv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'polaramv', num_obs   
      num_obs = 0 
      do n = 1, iv%info(polaramv)%nlocal
       if (iv%info(polaramv)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(polaramv)%levels(n)
         do k = 1, iv%info(polaramv)%levels(n)
!        Add random perturbation:
            call da_add_noise( iv % polaramv(n) % u(k), ob % polaramv(n) % u(k), z1)
            call da_add_noise( iv % polaramv(n) % v(k), ob % polaramv(n) % v(k), z2)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % polaramv(n) % u(k) % error, z1, &
                               iv % polaramv(n) % v(k) % error, z2, &
                               dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  [2.5] Transfer gpspw obs:

   if ( iv%info(gpspw)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(gpspw)%nlocal
       if(iv%info(gpspw)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'gpspw', num_obs   
      num_obs = 0 
      do n = 1, iv%info(gpspw)%nlocal
       if(iv%info(gpspw)%proc_domain(1,n)) then
        num_obs = num_obs + 1
        write(ounit,'(i8)')  1
!        Add random perturbation:
         call da_add_noise( iv % gpspw(n) % tpw, ob % gpspw(n) % tpw, z1 )
         
!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, iv % gpspw(n) % tpw % error, z1, &
                               dum, dum, dum, dum, dum, dum, dum, dum
       end if
      end do
   end if

!  [2.6] Transfer sonde obs:

   if ( iv%info(sound)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(sound)%nlocal
       if(iv%info(sound)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'sound', num_obs   
      num_obs = 0 
      do n = 1, iv%info(sound)%nlocal
       if(iv%info(sound)%proc_domain(1,n)) then
          num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(sound)%levels(n)
         do k = 1, iv%info(sound)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % sound(n) % u(k), ob % sound(n) % u(k), z1)
            call da_add_noise( iv % sound(n) % v(k), ob % sound(n) % v(k), z2)
            call da_add_noise( iv % sound(n) % t(k), ob % sound(n) % t(k), z3)
            call da_add_noise( iv % sound(n) % q(k), ob % sound(n) % q(k), z4)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % sound(n) % u(k) % error, z1, &
                               iv % sound(n) % v(k) % error, z2, &
                               iv % sound(n) % t(k) % error, z3, &
                               iv % sound(n) % q(k) % error, z4, &
                               dum, dum
         end do
       end if
      end do
   end if

! Transfer sonde_sfc obs:
   if ( iv%info(sonde_sfc)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(sonde_sfc)%nlocal
       if(iv%info(sonde_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'sonde_sfc', num_obs   
      num_obs = 0 
      do n = 1, iv%info(sonde_sfc)%nlocal
       if(iv%info(sonde_sfc)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') 1
!           Add random perturbation:
         call da_add_noise( iv % sonde_sfc(n) % u, ob % sonde_sfc(n) % u, z1 )
         call da_add_noise( iv % sonde_sfc(n) % v, ob % sonde_sfc(n) % v, z2 )
         call da_add_noise( iv % sonde_sfc(n) % t, ob % sonde_sfc(n) % t, z3 )
         call da_add_noise( iv % sonde_sfc(n) % p, ob % sonde_sfc(n) % p, z4 )
         call da_add_noise( iv % sonde_sfc(n) % q, ob % sonde_sfc(n) % q, z5 )

!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, iv % sonde_sfc(n) % u % error, z1, &
                                  iv % sonde_sfc(n) % v % error, z2, &
                                  iv % sonde_sfc(n) % t % error, z3, &
                                  iv % sonde_sfc(n) % p % error, z4, &
                                  iv % sonde_sfc(n) % q % error, z5
       end if
      end do
   end if

!  [2.7] Transfer airep obs:

   if ( iv%info(airep)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(airep)%nlocal
         if (iv%info(airep)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'airep', num_obs   
      num_obs = 0 
      do n = 1, iv%info(airep)%nlocal
       if (iv%info(airep)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(airep)%levels(n)
         do k = 1, iv%info(airep)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % airep(n) % u(k), ob % airep(n) % u(k), z1)
            call da_add_noise( iv % airep(n) % v(k), ob % airep(n) % v(k), z2)
            call da_add_noise( iv % airep(n) % t(k), ob % airep(n) % t(k), z3)
            call da_add_noise( iv % airep(n) % q(k), ob % airep(n) % q(k), z4)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % airep(n) % u(k) % error, z1, &
                               iv % airep(n) % v(k) % error, z2, &
                               iv % airep(n) % t(k) % error, z3, &
                               iv % airep(n) % q(k) % error, z4, &
                               dum, dum
         end do
       end if
      end do
   end if

!  [2.8] Transfer pilot obs:

   if ( iv%info(pilot)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(pilot)%nlocal
       if(iv%info(pilot)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'pilot', num_obs   
      num_obs = 0 
      do n = 1, iv%info(pilot)%nlocal
       if(iv%info(pilot)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') iv%info(pilot)%levels(n)
         do k = 1, iv%info(pilot)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % pilot(n) % u(k), ob % pilot(n) % u(k), z1)
            call da_add_noise( iv % pilot(n) % v(k), ob % pilot(n) % v(k), z2)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % pilot(n) % u(k) % error, z1, &
                               iv % pilot(n) % v(k) % error, z2, &
                               dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  [2.9] Transfer SSM/I obs:SSMI:

   if ( iv%info(ssmi_rv)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(ssmi_rv)%nlocal
       if(iv%info(ssmi_rv)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'ssmir', num_obs   
      num_obs = 0 
      do n = 1, iv%info(ssmi_rv)%nlocal
       if(iv%info(ssmi_rv)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') 1 
 
!        Add random perturbation:
         call da_add_noise( iv % ssmi_rv(n) % speed, &
                            ob % ssmi_rv(n) % speed, z1 )
         call da_add_noise( iv % ssmi_rv(n) % tpw, &
                            ob % ssmi_rv(n) % tpw, z2 )
!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, &
                                  iv % ssmi_rv(n) % speed % error, z1, &
                                  iv % ssmi_rv(n) % tpw % error, z2,   & 
                                  dum, dum, dum, dum, dum, dum
       end if
      end do
   end if

   if ( iv%info(ssmi_tb)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(ssmi_tb)%nlocal
       if(iv%info(ssmi_tb)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'ssmiT', num_obs   
      num_obs = 0 
      do n = 1, iv%info(ssmi_tb)%nlocal
       if(iv%info(ssmi_tb)%proc_domain(1,n)) then
         num_obs = num_obs + 1
!        Add random perturbation:
         call da_add_noise( iv % ssmi_tb(n) % tb19h, &
                            ob % ssmi_tb(n) % tb19h, z1)
         call da_add_noise( iv % ssmi_tb(n) % tb19v, &
                            ob % ssmi_tb(n) % tb19v, z2)
         call da_add_noise( iv % ssmi_tb(n) % tb22v, &
                            ob % ssmi_tb(n) % tb22v, z3)
         call da_add_noise( iv % ssmi_tb(n) % tb37h, &
                            ob % ssmi_tb(n) % tb37h, z4)
         call da_add_noise( iv % ssmi_tb(n) % tb37v, &
                            ob % ssmi_tb(n) % tb37v, z5)
         call da_add_noise( iv % ssmi_tb(n) % tb85h, &
                            ob % ssmi_tb(n) % tb85h, z6)
         call da_add_noise( iv % ssmi_tb(n) % tb85v, &
                            ob % ssmi_tb(n) % tb85v, z7)

!        Write out data:
         write(ounit,'(i8)') 1 
         write(ounit,'(2i8,14e15.7)')num_obs, 1, &
                                  iv % ssmi_tb(n) % tb19h % error, z1, &
                                  iv % ssmi_tb(n) % tb19v % error, z2, &
                                  iv % ssmi_tb(n) % tb22v % error, z3, &
                                  iv % ssmi_tb(n) % tb37h % error, z4, &
                                  iv % ssmi_tb(n) % tb37v % error, z5, &
                                  iv % ssmi_tb(n) % tb85h % error, z6, &
                                  iv % ssmi_tb(n) % tb85v % error, z7
       end if
      end do
   end if

!  [2.10] Transfer satem obs:

   if ( iv%info(satem)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(satem)%nlocal
       if(iv%info(satem)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'satem', num_obs   
      num_obs = 0 
      do n = 1, iv%info(satem)%nlocal
       if(iv%info(satem)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(satem)%levels(n)
         do k = 1, iv%info(satem)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % satem(n) % thickness(k), &
                               ob % satem(n) % thickness(k), z1 )
!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                                     iv % satem(n) % thickness(k) % error, z1, &
                                     dum, dum, dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if
   
!  [2.11] Transfer ssmt1 obs:

   if ( iv%info(ssmt1)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(ssmt1)%nlocal
       if(iv%info(ssmt1)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'ssmt1', num_obs   
      num_obs = 0 

      do n = 1, iv%info(ssmt1)%nlocal
       if(iv%info(ssmt1)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(ssmt1)%levels(n)
         
         do k = 1, iv%info(ssmt1)%levels(n)

!           Add random perturbation:
            call da_add_noise( iv % ssmt1(n) % t(k), &
                               ob % ssmt1(n) % t(k), z1 )
!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, iv % ssmt1(n) % t(k) % error, z1, &
                                     dum, dum, dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  [2.12] Transfer ssmt2 obs:

   if ( iv%info(ssmt2)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(ssmt2)%nlocal
       if(iv%info(ssmt2)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'ssmt2', num_obs   
      num_obs = 0 

      do n = 1, iv%info(ssmt2)%nlocal
       if(iv%info(ssmt2)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(ssmt2)%levels(n)
         
         do k = 1, iv%info(ssmt2)%levels(n)

!           Add random perturbation:
            call da_add_noise( iv % ssmt2(n) % rh(k), &
                               ob % ssmt2(n) % rh(k), z1 )
!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, iv % ssmt2(n) % rh(k) % error, z1, &
                                     dum, dum, dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if
   
!  [2.13] Transfer scatterometer obs:

   if ( iv%info(qscat)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(qscat)%nlocal
       if(iv%info(qscat)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'qscat', num_obs   
      num_obs = 0 
      do n = 1, iv%info(qscat)%nlocal
       if(iv%info(qscat)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') 1
!        Add random perturbation:
        call da_add_noise( iv % qscat(n) % u, ob % qscat(n) % u, z1 )
        call da_add_noise( iv % qscat(n) % v, ob % qscat(n) % v, z2 )

!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, &
                                  iv % qscat(n) % u % error, z1, &
                                  iv % qscat(n) % v % error, z2, &
                                  dum, dum, dum, dum, dum, dum
       end if
      end do
   end if

!  [2.14] Transfer buoy obs:

   if ( iv%info(buoy)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(buoy)%nlocal
       if(iv%info(buoy)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'buoy', num_obs   
      num_obs = 0 
      do n = 1, iv%info(buoy)%nlocal
       if(iv%info(buoy)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') 1
!        Add random perturbation:
         call da_add_noise( iv % buoy(n) % u, ob % buoy(n) % u, z1 )
         call da_add_noise( iv % buoy(n) % v, ob % buoy(n) % v, z2 )
         call da_add_noise( iv % buoy(n) % t, ob % buoy(n) % t, z3 )
         call da_add_noise( iv % buoy(n) % p, ob % buoy(n) % p, z4 )
         call da_add_noise( iv % buoy(n) % q, ob % buoy(n) % q, z5 )

!        Write out data:
        write(ounit,'(2i8,10e15.7)')num_obs, 1, &
                                  iv % buoy(n) % u % error, z1, &
                                  iv % buoy(n) % v % error, z2, &
                                  iv % buoy(n) % t % error, z3, &
                                  iv % buoy(n) % p % error, z4, &
                                  iv % buoy(n) % q % error, z5
      end if
     end do
   end if

!  [2.15] Transfer profiler obs:

   if ( iv%info(profiler)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(profiler)%nlocal
         if(iv%info(profiler)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'profiler', num_obs   
      num_obs = 0 
      do n = 1, iv%info(profiler)%nlocal
       if(iv%info(profiler)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(profiler)%levels(n)
         do k = 1, iv%info(profiler)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % profiler(n) % u(k), ob % profiler(n) % u(k), z1)
            call da_add_noise( iv % profiler(n) % v(k), ob % profiler(n) % v(k), z2)
!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % profiler(n) % u(k) % error, z1, &
                               iv % profiler(n) % v(k) % error, z2, &
                               dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  [2.16] Transfer TC bogus obs:

   if ( iv%info(bogus)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(bogus)%nlocal
       if(iv%info(bogus)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'bogus', num_obs   
      num_obs = 0 

      do n = 1, iv%info(bogus)%nlocal
       if(iv%info(bogus)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') 1
         call da_add_noise( iv % bogus(n) % slp, ob % bogus(n) % slp, z1 )
         write(ounit,'(2i8,10e15.7)')num_obs, 1, &
                                  iv % bogus(n) % slp % error, z1, &
                                  dum, dum, dum, dum, dum, dum, dum, dum

         write(ounit,'(i8)')iv%info(bogus)%levels(n)
         do k = 1, iv%info(bogus)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % bogus(n) % u(k), ob % bogus(n) % u(k), z1)
            call da_add_noise( iv % bogus(n) % v(k), ob % bogus(n) % v(k), z2)
            call da_add_noise( iv % bogus(n) % t(k), ob % bogus(n) % t(k), z3)
            call da_add_noise( iv % bogus(n) % q(k), ob % bogus(n) % q(k), z4)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % bogus(n) % u(k) % error, z1, &
                               iv % bogus(n) % v(k) % error, z2, &
                               iv % bogus(n) % t(k) % error, z3, &
                               iv % bogus(n) % q(k) % error, z4, &
                               dum, dum
         end do
       end if
      end do
   end if
!
!  Transfer AIRS retrievals:
!
   if ( iv%info(airsr)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(airsr)%nlocal
       if(iv%info(airsr)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'airsr', num_obs   
      num_obs = 0 
      do n = 1, iv%info(airsr)%nlocal
       if(iv%info(airsr)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(airsr)%levels(n)
         do k = 1, iv%info(airsr)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % airsr(n) % t(k), ob % airsr(n) % t(k), z1)
            call da_add_noise( iv % airsr(n) % q(k), ob % airsr(n) % q(k), z2)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % airsr(n) % t(k) % error, z1, &
                               iv % airsr(n) % q(k) % error, z2, &
                               dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  Transfer gpsref obs:

   if ( iv%info(gpsref)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(gpsref)%nlocal
       if(iv%info(gpsref)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'gpsref', num_obs   
      num_obs = 0 
      do n = 1, iv%info(gpsref)%nlocal
       if(iv%info(gpsref)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(gpsref)%levels(n)
         do k = 1, iv%info(gpsref)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % gpsref(n) % ref(k), ob % gpsref(n) % ref(k), z1)
!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % gpsref(n) % ref(k) % error, z1, &
                               dum, dum, dum, dum, dum, dum, dum, dum
         end do
       end if
      end do
   end if

!  Transfer mtgirs obs:

   if ( iv%info(mtgirs)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(mtgirs)%nlocal
       if(iv%info(mtgirs)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'mtgirs', num_obs

      num_obs = 0
      do n = 1, iv%info(mtgirs)%nlocal
       if(iv%info(mtgirs)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(mtgirs)%levels(n)
         do k = 1, iv%info(mtgirs)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % mtgirs(n) % u(k), ob % mtgirs(n) % u(k), z1)
            call da_add_noise( iv % mtgirs(n) % v(k), ob % mtgirs(n) % v(k), z2)
            call da_add_noise( iv % mtgirs(n) % t(k), ob % mtgirs(n) % t(k), z3)
            call da_add_noise( iv % mtgirs(n) % q(k), ob % mtgirs(n) % q(k), z4)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % mtgirs(n) % u(k) % error, z1, &
                               iv % mtgirs(n) % v(k) % error, z2, &
                               iv % mtgirs(n) % t(k) % error, z3, &
                               iv % mtgirs(n) % q(k) % error, z4, &
                               dum, dum
         end do
       end if
      end do
   end if

!  Transfer tamdar obs:

   if ( iv%info(tamdar)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(tamdar)%nlocal
       if(iv%info(tamdar)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'tamdar', num_obs

      num_obs = 0
      do n = 1, iv%info(tamdar)%nlocal
       if(iv%info(tamdar)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)')iv%info(tamdar)%levels(n)
         do k = 1, iv%info(tamdar)%levels(n)
!           Add random perturbation:
            call da_add_noise( iv % tamdar(n) % u(k), ob % tamdar(n) % u(k), z1)
            call da_add_noise( iv % tamdar(n) % v(k), ob % tamdar(n) % v(k), z2)
            call da_add_noise( iv % tamdar(n) % t(k), ob % tamdar(n) % t(k), z3)
            call da_add_noise( iv % tamdar(n) % q(k), ob % tamdar(n) % q(k), z4)

!           Write out data:
            write(ounit,'(2i8,10e15.7)')num_obs, k, &
                               iv % tamdar(n) % u(k) % error, z1, &
                               iv % tamdar(n) % v(k) % error, z2, &
                               iv % tamdar(n) % t(k) % error, z3, &
                               iv % tamdar(n) % q(k) % error, z4, &
                               dum, dum
         end do
       end if
      end do
   end if

!  Transfer tamdar_sfc obs:

   if ( iv%info(tamdar_sfc)%nlocal > 0 ) then
      num_obs = 0
      do n = 1, iv%info(tamdar_sfc)%nlocal
       if(iv%info(tamdar_sfc)%proc_domain(1,n)) num_obs = num_obs + 1
      end do
      write(ounit,'(a20,i8)')'tamdar_sfc', num_obs
      num_obs = 0
      do n = 1, iv%info(tamdar_sfc)%nlocal
       if(iv%info(tamdar_sfc)%proc_domain(1,n)) then
         num_obs = num_obs + 1
         write(ounit,'(i8)') 1
!           Add random perturbation:
         call da_add_noise( iv % tamdar_sfc(n) % u, ob % tamdar_sfc(n) % u, z1 )
         call da_add_noise( iv % tamdar_sfc(n) % v, ob % tamdar_sfc(n) % v, z2 )
         call da_add_noise( iv % tamdar_sfc(n) % t, ob % tamdar_sfc(n) % t, z3 )
         call da_add_noise( iv % tamdar_sfc(n) % p, ob % tamdar_sfc(n) % p, z4 )
         call da_add_noise( iv % tamdar_sfc(n) % q, ob % tamdar_sfc(n) % q, z5 )

!        Write out data:
         write(ounit,'(2i8,10e15.7)')num_obs, 1, iv % tamdar_sfc(n) % u % error, z1, &
                                  iv % tamdar_sfc(n) % v % error, z2, &
                                  iv % tamdar_sfc(n) % t % error, z3, &
                                  iv % tamdar_sfc(n) % p % error, z4, &
                                  iv % tamdar_sfc(n) % q % error, z5
       end if
      end do

   end if

!
!  Transfer Radiance obs:
!

   if ( iv%num_inst > 0 ) then
      do i = 1, iv%num_inst                 ! loop for sensor
         if ( iv%instid(i)%num_rad < 1 ) cycle
         do k = 1,iv%instid(i)%nchan        ! loop for channel
!  Counting number of obs for channle k
         num_obs = 0
         do n = 1,iv%instid(i)%num_rad      ! loop for pixel
           if(iv%instid(i)%info%proc_domain(1,n) .and. &
              (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer)) then
                num_obs = num_obs + 1
              end if
         end do                                ! end loop for pixel
         if (num_obs < 1) cycle

         write(ob_name,'(a,a,i4.4)') trim(iv%instid(i)%rttovid_string),'-',k
         write(ounit,'(a20,i8)')  ob_name,num_obs

         num_obs = 0
         do n= 1, iv%instid(i)%num_rad      ! loop for pixel
               if(iv%instid(i)%info%proc_domain(1,n) .and. &
                  (iv%instid(i)%tb_qc(k,n) >= obs_qc_pointer)) then
                     num_obs = num_obs + 1
                     call da_add_noise_new( iv%instid(i)%tb_qc(k,n), &
                                        iv%instid(i)%tb_error(k,n),  &
                                        iv%instid(i)%tb_inv(k,n),  &
                                        ob%instid(i)%tb(k,n), z1)

                     write(ounit,'(2i8,f10.3,e15.7)') num_obs, 1,  &
                              iv%instid(i)%tb_error(k,n), z1
               end if
         end do                                ! end loop for pixel
         end do                                ! end loop for channel
      end do                                   ! end loop for sensor
   end if

  close (ounit)
  call da_free_unit(ounit)

   if (trace_use_dull) call da_trace_exit("da_add_noise_to_ob")

end subroutine da_add_noise_to_ob


subroutine da_check_missing(qc_flag, y_in, y_out)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none
   
   integer, intent(in)   :: qc_flag
   real, intent(in)      :: y_in
   real, intent(out)     :: y_out

   if (qc_flag < obs_qc_pointer) then
      y_out = missing_r
   else
      y_out = y_in
   end if

end subroutine da_check_missing


subroutine da_fill_obs_structures(iv, ob, uvq_direct)

   !----------------------------------------------------------------------------   
   ! Purpose: Allocates observation structure and fills it from iv.
   !----------------------------------------------------------------------------   

   implicit none

   type (iv_type), intent(inout) :: iv   ! Obs and header structure.
   type (y_type), intent(out)    :: ob   ! (Smaller) observation structure.
   logical, optional             :: uvq_direct  !flag for having direct u,v,q obs

   integer :: n, k     ! Loop counters.
   real    :: rh_error ! RH obs. error.
   real    :: q_error  ! q obs. error.
   real    :: geometric_h, geopotential_h
   integer :: i,j
   logical :: outside

   if (trace_use) call da_trace_entry("da_fill_obs_structures")

   !---------------------------------------------------------------------------
   ! Initialise obs error factors (which will be overwritten in use_obs_errfac)
   !---------------------------------------------------------------------------

   iv % synop_ef_u = 1.0
   iv % synop_ef_v = 1.0
   iv % synop_ef_t = 1.0
   iv % synop_ef_p = 1.0
   iv % synop_ef_q = 1.0

   iv % metar_ef_u = 1.0
   iv % metar_ef_v = 1.0
   iv % metar_ef_t = 1.0
   iv % metar_ef_p = 1.0
   iv % metar_ef_q = 1.0

   iv % ships_ef_u = 1.0
   iv % ships_ef_v = 1.0
   iv % ships_ef_t = 1.0
   iv % ships_ef_p = 1.0
   iv % ships_ef_q = 1.0

   iv % geoamv_ef_u = 1.0
   iv % geoamv_ef_v = 1.0

   iv % polaramv_ef_u = 1.0
   iv % polaramv_ef_v = 1.0

   iv % gpspw_ef_tpw = 1.0

   iv % gpsref_ef_ref = 1.0
   iv % gpsref_ef_p = 1.0
   iv % gpsref_ef_t = 1.0
   iv % gpsref_ef_q = 1.0

   iv % sound_ef_u = 1.0
   iv % sound_ef_v = 1.0
   iv % sound_ef_t = 1.0
   iv % sound_ef_q = 1.0

   iv % airep_ef_u = 1.0
   iv % airep_ef_v = 1.0
   iv % airep_ef_t = 1.0
   iv % airep_ef_q = 1.0

   iv % pilot_ef_u = 1.0
   iv % pilot_ef_v = 1.0

   iv % ssmir_ef_speed = 1.0
   iv % ssmir_ef_tpw = 1.0

   iv % satem_ef_thickness = 1.0

   iv % ssmt1_ef_t = 1.0

   iv % ssmt2_ef_rh = 1.0

   iv % qscat_ef_u = 1.0
   iv % qscat_ef_v = 1.0

   iv % profiler_ef_u = 1.0
   iv % profiler_ef_v = 1.0
   
   iv % buoy_ef_u = 1.0
   iv % buoy_ef_v = 1.0
   iv % buoy_ef_t = 1.0
   iv % buoy_ef_p = 1.0
   iv % buoy_ef_q = 1.0

   iv % radar_ef_rv = 1.0
   iv % radar_ef_rf = 1.0

   iv % rain_ef_r  = 1.0

   iv % bogus_ef_u = 1.0
   iv % bogus_ef_v = 1.0
   iv % bogus_ef_t = 1.0
   iv % bogus_ef_p = 1.0
   iv % bogus_ef_q = 1.0
   iv % bogus_ef_slp = 1.0

   iv % airsr_ef_t = 1.0
   iv % airsr_ef_q = 1.0

   iv % mtgirs_ef_u = 1.0
   iv % mtgirs_ef_v = 1.0
   iv % mtgirs_ef_t = 1.0
   iv % mtgirs_ef_q = 1.0

   iv % tamdar_ef_u = 1.0
   iv % tamdar_ef_v = 1.0
   iv % tamdar_ef_t = 1.0
   iv % tamdar_ef_q = 1.0

   iv % tamdar_sfc_ef_u = 1.0
   iv % tamdar_sfc_ef_v = 1.0
   iv % tamdar_sfc_ef_t = 1.0
   iv % tamdar_sfc_ef_p = 1.0
   iv % tamdar_sfc_ef_q = 1.0

   !----------------------------------------------------------------------
   ! [1.0] Allocate innovation vector and observation structures:
   !----------------------------------------------------------------------
   call da_allocate_y(iv, ob)

   !----------------------------------------------------------------------
   ! [2.0] Transfer observations:
   !----------------------------------------------------------------------

   ! [2.1] Transfer surface obs:

   if (iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
         ob % synop(n) % u = iv % synop(n) % u % inv
         ob % synop(n) % v = iv % synop(n) % v % inv
         ob % synop(n) % t = iv % synop(n) % t % inv
         ob % synop(n) % q = iv % synop(n) % q % inv
         ob % synop(n) % p = iv % synop(n) % p % inv

         ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
         rh_error = iv%synop(n)%q%error ! q error is rh at this stage!

         ! if((ob % synop(n) % p > iv%ptop) .AND. &
         !    (ob % synop(n) % t > 100.0) .AND. &
         !    (ob % synop(n) % q > 0.0) .AND. &
         !    (iv % synop(n) % p % qc >= obs_qc_pointer) .and. &
         !    (iv % synop(n) % t % qc >= obs_qc_pointer) .and. &
         !    (iv % synop(n) % q % qc >= obs_qc_pointer)) then
         call da_get_q_error(ob % synop(n) % p, &
                              ob % synop(n) % t, &
                              ob % synop(n) % q, &
                              iv % synop(n) % t % error, &
                              rh_error, iv % synop(n) % q % error)
         if (iv%synop(n)% q % error == missing_r) iv%synop(n)% q % qc = missing_data

         ! end if
         end if
      end do      
   end if

   ! [2.2] Transfer metar obs:

   if (iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
         ob % metar(n) % u = iv % metar(n) % u % inv
         ob % metar(n) % v = iv % metar(n) % v % inv
         ob % metar(n) % t = iv % metar(n) % t % inv
         ob % metar(n) % q = iv % metar(n) % q % inv
         ob % metar(n) % p = iv % metar(n) % p % inv

         ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
         rh_error = iv%metar(n)%q%error ! q error is rh at this stage!
         call da_get_q_error(iv % metar(n) % p % inv, &
                              ob % metar(n) % t, &
                              ob % metar(n) % q, &
                              iv % metar(n) % t % error, &
                              rh_error, q_error)
         iv % metar(n) % q % error = q_error
         if (iv%metar(n)% q % error == missing_r) &
            iv%metar(n)% q % qc = missing_data
         end if
      end do
   end if

   ! [2.2] Transfer ships obs:

   if (iv%info(ships)%nlocal > 0) then   
      do n = 1, iv%info(ships)%nlocal
         ob % ships(n) % u = iv % ships(n) % u % inv
         ob % ships(n) % v = iv % ships(n) % v % inv
         ob % ships(n) % t = iv % ships(n) % t % inv
         ob % ships(n) % q = iv % ships(n) % q % inv
         ob % ships(n) % p = iv % ships(n) % p % inv

         ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
         rh_error = iv%ships(n)%q%error ! q error is rh at this stage!
         call da_get_q_error(iv % ships(n) % p % inv, &
                              ob % ships(n) % t, &
                              ob % ships(n) % q, &
                              iv % ships(n) % t % error, &
                              rh_error, q_error)
         iv % ships(n) % q % error = q_error

         if(iv%ships(n)% q % error == missing_r) iv%ships(n)% q % qc = missing_data
         end if
      end do
      
   end if

   ! [2.4.1] Transfer Geo. AMVs Obs:

   if (iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
         do k = 1, iv%info(geoamv)%levels(n)
            ob % geoamv(n) % u(k) = iv % geoamv(n) % u(k) % inv
            ob % geoamv(n) % v(k) = iv % geoamv(n) % v(k) % inv
         end do
      end do
   end if

   ! [2.4.2] Transfer  Polar AMVs Obs:

   if (iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
         do k = 1,iv%info(polaramv)%levels(n)
            ob % polaramv(n) % u(k) = iv % polaramv(n) % u(k) % inv
            ob % polaramv(n) % v(k) = iv % polaramv(n) % v(k) % inv
         end do
      end do
   end if

   ! [2.5] Transfer gpspw obs:

   if (iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
         ob % gpspw(n) % tpw = iv % gpspw(n) % tpw % inv
      end do

   end if

   ! [2.6] Transfer GPS REF obs:

   if (iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
         do k = 1, iv%info(gpsref)%levels(n)
! ............................................................................
! Convert the geometric height to the geopotential height for GPSREF data
! because GPSRO used geometric height for impact parameter, and the WRF model
! always uses the geopotential height as the vertical coordinate for P, T, q.
! YRG (10/05/2010):
            geometric_h = iv%gpsref(n)%h(k) / 1000.0
            call da_msl2geo1 (geometric_h, iv%info(gpsref)%lat(1,n),  &
                           iv%info(gpsref)%lon(1,n), geopotential_h)
            iv % gpsref(n) % h(k) = geopotential_h * 1000.0
!            write (999,'("n=",i6," k=",i3,2x,"gmh, gph, lat, lon:",4f15.5)') &
!                          n, k, geometric_h*1000.0, iv % gpsref(n) % h(k), &
!                          iv%info(gpsref)%lat(1,n), iv%info(gpsref)%lon(1,n)
! ............................................................................
            ob % gpsref(n) % ref(k) = iv % gpsref(n) % ref(k) % inv
            ob % gpsref(n) %   p(k) = iv % gpsref(n) %   p(k) % inv
            ob % gpsref(n) %   t(k) = iv % gpsref(n) %   t(k) % inv
            ob % gpsref(n) %   q(k) = iv % gpsref(n) %   q(k) % inv
         end do
      end do
   end if

   ! [2.7] Transfer sonde obs:

   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         do k = 1, iv%info(sound)%levels(n)
            ob % sound(n) % u(k) = iv % sound(n) % u(k) % inv
            ob % sound(n) % v(k) = iv % sound(n) % v(k) % inv
            ob % sound(n) % t(k) = iv % sound(n) % t(k) % inv
            ob % sound(n) % q(k) = iv % sound(n) % q(k) % inv

            ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
            rh_error = iv%sound(n)%q(k)%error ! q error is rh at this stage!
            call da_get_q_error(iv % sound(n) % p(k), &
                                 ob % sound(n) % t(k), &
                                 ob % sound(n) % q(k), &
                                 iv % sound(n) % t(k) % error, &
                                 rh_error, q_error)

            iv % sound(n) % q(k) % error = q_error
         if (iv%sound(n)% q(k) % error == missing_r) &
            iv%sound(n)% q(k) % qc = missing_data
         end if
         end do
      end do
   end if

   if (iv%info(sonde_sfc)%nlocal > 0) then
      do n = 1, iv%info(sonde_sfc)%nlocal
         ob % sonde_sfc(n) % u = iv % sonde_sfc(n) % u % inv
         ob % sonde_sfc(n) % v = iv % sonde_sfc(n) % v % inv
         ob % sonde_sfc(n) % t = iv % sonde_sfc(n) % t % inv
         ob % sonde_sfc(n) % q = iv % sonde_sfc(n) % q % inv
         ob % sonde_sfc(n) % p = iv % sonde_sfc(n) % p % inv

         ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
         rh_error = iv%sonde_sfc(n)%q%error ! q error is rh at this stage!
         call da_get_q_error(iv % sonde_sfc(n) % p % inv, &
                              ob % sonde_sfc(n) % t, &
                              ob % sonde_sfc(n) % q, &
                              iv % sonde_sfc(n) % t % error, &
                              rh_error, iv % sonde_sfc(n) % q % error)
         if (iv%sonde_sfc(n)% q % error == missing_r) &
            iv%sonde_sfc(n)% q % qc = missing_data
         end if
      end do
   end if

   ! [2.8] Transfer airep obs:

   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         do k = 1, iv%info(airep)%levels(n)
            ob % airep(n) % u(k) = iv % airep(n) % u(k) % inv
            ob % airep(n) % v(k) = iv % airep(n) % v(k) % inv
            ob % airep(n) % t(k) = iv % airep(n) % t(k) % inv
            ob % airep(n) % q(k) = iv % airep(n) % q(k) % inv

            if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
                rh_error = iv%airep(n)%q(k)%error ! q error is rh at this stage!
                call da_get_q_error(iv % airep(n) % p(k), &
                                    ob % airep(n) % t(k), &
                                    ob % airep(n) % q(k), &
                                    iv % airep(n) % t(k) % error, &
                                    rh_error, q_error)

                iv % airep(n) % q(k) % error = q_error
                if (iv%airep(n)% q(k) % error == missing_r) &
                    iv%airep(n)% q(k) % qc = missing_data
            end if
         end do
      end do
   end if

   ! [2.9] Transfer pilot obs:

   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
         do k = 1, iv%info(pilot)%levels(n)
            ob % pilot(n) % u(k) = iv % pilot(n) % u(k) % inv
            ob % pilot(n) % v(k) = iv % pilot(n) % v(k) % inv
         end do
      end do
   end if

   ! [2.10] Transfer SSM/I obs:SSMI:

   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
         ob % ssmi_rv(n) % speed = iv % ssmi_rv(n) % speed % inv
         ob % ssmi_rv(n) % tpw   = iv % ssmi_rv(n) % tpw % inv
      end do
   end if

   if (iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
         ob % ssmi_tb(n) % tb19v = iv % ssmi_tb(n) % tb19v % inv
         ob % ssmi_tb(n) % tb19h = iv % ssmi_tb(n) % tb19h % inv
         ob % ssmi_tb(n) % tb22v = iv % ssmi_tb(n) % tb22v % inv
         ob % ssmi_tb(n) % tb37v = iv % ssmi_tb(n) % tb37v % inv
         ob % ssmi_tb(n) % tb37h = iv % ssmi_tb(n) % tb37h % inv
         ob % ssmi_tb(n) % tb85v = iv % ssmi_tb(n) % tb85v % inv
         ob % ssmi_tb(n) % tb85h = iv % ssmi_tb(n) % tb85h % inv
      end do
   end if

   ! [2.11] Transfer satem obs:

   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         do k = 1, iv%info(satem)%levels(n)
            ob % satem(n) % thickness(k) = iv % satem(n) % thickness(k) % inv
         end do
      end do
   end if
   
   ! [2.12] Transfer ssmt1 obs:

   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
         do k = 1, iv%info(ssmt1)%levels(n)
            ob % ssmt1(n) % t(k) = iv % ssmt1(n) % t(k) % inv
         end do
      end do

   end if

   ! [2.13] Transfer ssmt2 obs:

   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         do k = 1, iv%info(ssmt2)%levels(n)
            ob % ssmt2(n) % rh(k) = iv % ssmt2(n) % rh(k) % inv
         end do
      end do
   end if
   
   ! [2.14] Setup pseudo observations:

   if (iv%info(pseudo)%nlocal > 0) call da_setup_pseudo_obs(iv, ob)

   ! [2.15] Transfer scatterometer obs:

   if (iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
         ob % qscat(n) % u = iv % qscat(n) % u % inv
         ob % qscat(n) % v = iv % qscat(n) % v % inv
      end do     
   end if

   ! [2.16] Transfer profiler obs:

   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         do k = 1, iv%info(profiler)%levels(n)
            ob % profiler(n) % u(k) = iv % profiler(n) % u(k) % inv
            ob % profiler(n) % v(k) = iv % profiler(n) % v(k) % inv
         end do
      end do
   end if

   ! [2.17] Transfer buoy obs:

   if (iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
         ob % buoy(n) % p = iv % buoy(n) % p % inv
      end do
      do n = 1, iv%info(buoy)%nlocal
         ob % buoy(n) % u = iv % buoy(n) % u % inv
         ob % buoy(n) % v = iv % buoy(n) % v % inv
         ob % buoy(n) % t = iv % buoy(n) % t % inv
         ob % buoy(n) % q = iv % buoy(n) % q % inv

         ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
         rh_error = iv%buoy(n)%q%error ! q error is rh at this stage!
         call da_get_q_error(iv % buoy(n) % p % inv, &
                              ob % buoy(n) % t, &
                              ob % buoy(n) % q, &
                              iv % buoy(n) % t % error, &
                              rh_error, q_error)
         iv % buoy(n) % q % error = q_error

         if(iv%buoy (n)% q % error == missing_r) iv%buoy (n)% q % qc = missing_data
         end if
      end do
   end if

   ! [2.18] Transfer radar obs: this section has been moved to da_fill_obs_structures_radar.inc

!   if (iv%info(radar)%nlocal > 0) then
!      do n = 1, iv%info(radar)%nlocal
!         do k = 1, iv%info(radar)%levels(n)
!            ! Copy observation variables:
!            ob % radar(n) % rv(k) = iv % radar(n) % rv(k) % inv
!           ob % radar(n) % rf(k) = iv % radar(n) % rf(k) % inv
!         end do
!      end do
!   end if

   ! [2.19] Transfer TC bogus:

   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         do k = 1, iv%info(bogus)%levels(n)

            ! Copy observation variables:

            ob % bogus(n) % u(k) = iv % bogus(n) % u(k) % inv
            ob % bogus(n) % v(k) = iv % bogus(n) % v(k) % inv
            ob % bogus(n) % t(k) = iv % bogus(n) % t(k) % inv
            ob % bogus(n) % q(k) = iv % bogus(n) % q(k) % inv

            ! Calculate q error from rh error:

            rh_error = iv%bogus(n)%q(k)%error ! q error is rh at this stage!
            call da_get_q_error(iv % bogus(n) % p(k), &
                                 ob % bogus(n) % t(k), &
                                 ob % bogus(n) % q(k), &
                                 iv % bogus(n) % t(k) % error, &
                                 rh_error, q_error)

            iv % bogus(n) % q(k) % error = q_error
            if (iv%bogus(n)% q(k) % error == missing_r) &
               iv%bogus(n)% q(k) % qc = missing_data
         end do
         ob % bogus(n) % slp = iv % bogus(n) % slp % inv
      end do
   end if

   ! [2.20] Transfer rain obs:

   if (iv%info(rain)%nlocal > 0) then
      do n = 1, iv%info(rain)%nlocal
            ob % rain(n) % rain = iv % rain(n) % rain % inv
      end do
   end if

   ! Transfer AIRS  retrievals:

   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         do k = 1, iv%info(airsr)%levels(n)

            ! Copy observation variables:

            ob % airsr(n) % t(k) = iv % airsr(n) % t(k) % inv
            ob % airsr(n) % q(k) = iv % airsr(n) % q(k) % inv

            ! Calculate q error from rh error:

         if (.not. present(uvq_direct) .or. (present(uvq_direct) .and. (.not. uvq_direct))) then
            rh_error = iv%airsr(n)%q(k)%error ! q error is rh at this stage!
            call da_get_q_error(iv % airsr(n) % p(k), &
                                 ob % airsr(n) % t(k), &
                                 ob % airsr(n) % q(k), &
                                 iv % airsr(n) % t(k) % error, &
                                 rh_error, q_error)

            iv % airsr(n) % q(k) % error = q_error
            if (iv%airsr(n)% q(k) % error == missing_r) &
               iv%airsr(n)% q(k) % qc = missing_data
         end if
         end do
      end do
   end if
   if (iv%info(mtgirs)%nlocal > 0) then
      do n = 1, iv%info(mtgirs)%nlocal
         do k = 1, iv%info(mtgirs)%levels(n)
            ob % mtgirs(n) % u(k) = iv % mtgirs(n) % u(k) % inv
            ob % mtgirs(n) % v(k) = iv % mtgirs(n) % v(k) % inv
            ob % mtgirs(n) % t(k) = iv % mtgirs(n) % t(k) % inv
            ob % mtgirs(n) % q(k) = iv % mtgirs(n) % q(k) % inv
         if (iv%mtgirs(n)% q(k) % error == missing_r) &
            iv%mtgirs(n)% q(k) % qc = missing_data
         end do
      end do
   end if

   if (iv%info(tamdar)%nlocal > 0) then
      do n = 1, iv%info(tamdar)%nlocal
         do k = 1, iv%info(tamdar)%levels(n)
            ob % tamdar(n) % u(k) = iv % tamdar(n) % u(k) % inv
            ob % tamdar(n) % v(k) = iv % tamdar(n) % v(k) % inv
            ob % tamdar(n) % t(k) = iv % tamdar(n) % t(k) % inv
            ob % tamdar(n) % q(k) = iv % tamdar(n) % q(k) % inv

         if (iv%tamdar(n)% u(k) % error == missing_r) &
            iv%tamdar(n)% u(k) % qc = missing_data
         if (iv%tamdar(n)% v(k) % error == missing_r) &
            iv%tamdar(n)% v(k) % qc = missing_data
         if (iv%tamdar(n)% t(k) % error == missing_r) &
            iv%tamdar(n)% t(k) % qc = missing_data

            ! Calculate q error from rh error:

            rh_error = iv%tamdar(n)%q(k)%error ! q error is rh at this stage!
            call da_get_q_error(iv % tamdar(n) % p(k), &
                                ob % tamdar(n) % t(k), &
                                ob % tamdar(n) % q(k), &
                                iv % tamdar(n) % t(k) % error, &
                                rh_error, q_error)

            iv % tamdar(n) % q(k) % error = q_error

         if (iv%tamdar(n)% q(k) % error == missing_r) &
            iv%tamdar(n)% q(k) % qc = missing_data
         end do
      end do
   end if

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      do n = 1, iv%info(tamdar_sfc)%nlocal

         ob % tamdar_sfc(n) % u = iv % tamdar_sfc(n) % u % inv
         ob % tamdar_sfc(n) % v = iv % tamdar_sfc(n) % v % inv
         ob % tamdar_sfc(n) % t = iv % tamdar_sfc(n) % t % inv
         ob % tamdar_sfc(n) % q = iv % tamdar_sfc(n) % q % inv
         ob % tamdar_sfc(n) % p = iv % tamdar_sfc(n) % p % inv

         if (iv%tamdar_sfc(n)% u % error == missing_r) &
            iv%tamdar_sfc(n)% u % qc = missing_data
         if (iv%tamdar_sfc(n)% v % error == missing_r) &
            iv%tamdar_sfc(n)% v % qc = missing_data
         if (iv%tamdar_sfc(n)% t % error == missing_r) &
            iv%tamdar_sfc(n)% t % qc = missing_data

         ! Calculate q error from rh error:

         rh_error = iv%tamdar_sfc(n)%q%error ! q error is rh at this stage!
         call da_get_q_error(iv % tamdar_sfc(n) % p % inv, &
                              ob % tamdar_sfc(n) % t, &
                              ob % tamdar_sfc(n) % q, &
                              iv % tamdar_sfc(n) % t % error, &
                              rh_error, iv % tamdar_sfc(n) % q % error)
         if (iv%tamdar_sfc(n)% q % error == missing_r) &
            iv%tamdar_sfc(n)% q % qc = missing_data
      end do
   end if

   if (trace_use) call da_trace_exit("da_fill_obs_structures")

end subroutine da_fill_obs_structures


subroutine da_fill_obs_structures_radar(iv, ob)

   !----------------------------------------------------------------------------   
   ! Purpose: Allocates observation structure and fills it from iv.
   !----------------------------------------------------------------------------   

   implicit none

   type (iv_type), intent(inout) :: iv   ! Obs and header structure.
   type (y_type), intent(out)    :: ob   ! (Smaller) observation structure.

   integer :: n, k     ! Loop counters.
   real    :: rh_error ! RH obs. error.
   real    :: q_error  ! q obs. error.
   real    :: geometric_h, geopotential_h
   integer :: i,j
   logical :: outside

   if (trace_use) call da_trace_entry("da_fill_obs_structures_radar")

   !---------------------------------------------------------------------------
   ! Initialise obs error factors (which will be overwritten in use_obs_errfac)
   !---------------------------------------------------------------------------

   iv % radar_ef_rv = 1.0
   iv % radar_ef_rf = 1.0

   !----------------------------------------------------------------------
   ! [1.0] Allocate innovation vector and observation structures:
   !----------------------------------------------------------------------
   call da_allocate_y_radar(iv, ob)

   !----------------------------------------------------------------------
   ! [2.0] Transfer observations:
   !----------------------------------------------------------------------

   ! [2.18] Transfer radar obs:

   if (iv%info(radar)%nlocal > 0) then
      do n = 1, iv%info(radar)%nlocal
         do k = 1, iv%info(radar)%levels(n)
            ! Copy observation variables:
            ob % radar(n) % rv(k) = iv % radar(n) % rv(k) % inv
            ob % radar(n) % rf(k) = iv % radar(n) % rf(k) % inv
         end do
      end do
   end if

   if (trace_use) call da_trace_exit("da_fill_obs_structures_radar")

end subroutine da_fill_obs_structures_radar


subroutine da_fill_obs_structures_rain(iv, ob)

   !----------------------------------------------------------------------------   
   ! Purpose: Allocates observation structure and fills it from iv.
   !----------------------------------------------------------------------------   

   implicit none

   type (iv_type), intent(inout) :: iv   ! Obs and header structure.
   type (y_type), intent(out)    :: ob   ! (Smaller) observation structure.

   integer :: n, k     ! Loop counters.
   real    :: rh_error ! RH obs. error.
   real    :: q_error  ! q obs. error.
   real    :: geometric_h, geopotential_h
   integer :: i,j
   logical :: outside

   if (trace_use) call da_trace_entry("da_fill_obs_structures_rain")

   !---------------------------------------------------------------------------
   ! Initialise obs error factors (which will be overwritten in use_obs_errfac)
   !---------------------------------------------------------------------------

   iv % rain_ef_r  = 1.0

   !----------------------------------------------------------------------
   ! [1.0] Allocate innovation vector and observation structures:
   !----------------------------------------------------------------------
   call da_allocate_y_rain(iv, ob)

   !----------------------------------------------------------------------
   ! [2.0] Transfer observations:
   !----------------------------------------------------------------------

   ! [2.20] Transfer rain obs:

   if (iv%info(rain)%nlocal > 0) then
      do n = 1, iv%info(rain)%nlocal
            ob % rain(n) % rain = iv % rain(n) % rain % inv
      end do
   end if

   if (trace_use) call da_trace_exit("da_fill_obs_structures_rain")

end subroutine da_fill_obs_structures_rain


subroutine da_random_omb_all(iv, ob)

   !-------------------------------------------------------------------------
   ! Purpose: Allocates observation structure and fills it fro iv.
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(inout) :: iv   ! Obs and header structure.
   type (y_type), intent(inout)  :: ob   ! (Smaller) observation structure.

   integer                       :: n, k ! Loop counters.

   if (trace_use) call da_trace_entry("da_random_omb_all")

   !----------------------------------------------------------------------
   ! [1.0] Initialise random number sequence:
   !----------------------------------------------------------------------

   call da_random_seed
   
   !----------------------------------------------------------------------
   !  [2.0] Randomize each ob in turn:
   !----------------------------------------------------------------------

   ! [2.1] Transfer surface obs:

   if (iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
         call da_random_omb(iv % synop(n) % u % error, ob % synop(n) % u, &
                           iv % synop(n) % u % qc, iv % synop(n) % u % inv)
         call da_random_omb(iv % synop(n) % v % error, ob % synop(n) % v, &
                           iv % synop(n) % v % qc, iv % synop(n) % v % inv)
         call da_random_omb(iv % synop(n) % t % error, ob % synop(n) % t, &
                           iv % synop(n) % t % qc, iv % synop(n) % t % inv )
         call da_random_omb(iv % synop(n) % p % error, ob % synop(n) % p, &
                           iv % synop(n) % p % qc, iv % synop(n) % p % inv)
         call da_random_omb(iv % synop(n) % q % error, ob % synop(n) % q, &
                           iv % synop(n) % q % qc, iv % synop(n) % q % inv)
      end do
   end if

   ! [2.2] Transfer metar obs:

   if (iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
         call da_random_omb(iv % metar(n) % u % error, ob % metar(n) % u, &
                           iv % metar(n) % u % qc, iv % metar(n) % u % inv)
         call da_random_omb(iv % metar(n) % v % error, ob % metar(n) % v, &
                           iv % metar(n) % v % qc, iv % metar(n) % v % inv)
         call da_random_omb(iv % metar(n) % t % error, ob % metar(n) % t, &
                           iv % metar(n) % t % qc, iv % metar(n) % t % inv)
         call da_random_omb(iv % metar(n) % p % error, ob % metar(n) % p, &
                           iv % metar(n) % p % qc, iv % metar(n) % p % inv)
         call da_random_omb(iv % metar(n) % q % error, ob % metar(n) % q, &
                           iv % metar(n) % q % qc, iv % metar(n) % q % inv)
      end do
   end if

   ! [2.3] Transfer ships obs:

   if (iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
         call da_random_omb(iv % ships(n) % u % error, ob % ships(n) % u, &
                           iv % ships(n) % u % qc, iv % ships(n) % u % inv)
         call da_random_omb(iv % ships(n) % v % error, ob % ships(n) % v, &
                           iv % ships(n) % v % qc, iv % ships(n) % v % inv)
         call da_random_omb(iv % ships(n) % t % error, ob % ships(n) % t, &
                           iv % ships(n) % t % qc, iv % ships(n) % t % inv)
         call da_random_omb(iv % ships(n) % p % error, ob % ships(n) % p, &
                           iv % ships(n) % p % qc, iv % ships(n) % p % inv)
         call da_random_omb(iv % ships(n) % q % error, ob % ships(n) % q, &
                           iv % ships(n) % q % qc, iv % ships(n) % q % inv)
      end do
   end if

   ! [2.4.1] Transfer Geo. AMVs Obs:

   if (iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal 
        do k = 1, iv%info(geoamv)%levels(n)
         call da_random_omb(iv % geoamv(n) % u(k) % error, ob % geoamv(n) % u(k), &
                           iv % geoamv(n) % u(k) % qc, iv % geoamv(n) % u(k) % inv)
         call da_random_omb(iv % geoamv(n) % v(k) % error, ob % geoamv(n) % v(k), &
                           iv % geoamv(n) % v(k) % qc, iv % geoamv(n) % v(k) % inv)
        end do
      end do 
   end if 

   ! [2.4.2] Transfer Polar  AMVs Obs:

   if (iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
        do k = 1, iv%info(polaramv)%levels(n)
         call da_random_omb(iv % polaramv(n) % u(k) % error, ob % polaramv(n) % u(k), &
                           iv % polaramv(n) % u(k) % qc, iv % polaramv(n) % u(k) % inv)
         call da_random_omb(iv % polaramv(n) % v(k) % error, ob % polaramv(n) % v(k), &
                           iv % polaramv(n) % v(k) % qc, iv % polaramv(n) % v(k) % inv)
        end do
      end do
   end if

   ! [2.5] Transfer gpspw obs:

   if (iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
         call da_random_omb(iv % gpspw(n) % tpw % error, ob % gpspw(n) % tpw, &
                           iv % gpspw(n) % tpw % qc, iv % gpspw(n) % tpw % inv)
      end do
   end if

   ! [2.6] Transfer sonde obs:

   if (iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
         do k = 1, iv%info(sound)%levels(n)
            call da_random_omb(iv % sound(n) % u(k) % error, ob % sound(n) % u(k), &
                              iv % sound(n) % u(k) % qc, iv % sound(n) % u(k) % inv)
            call da_random_omb(iv % sound(n) % v(k) % error, ob % sound(n) % v(k), &
                              iv % sound(n) % v(k) % qc, iv % sound(n) % v(k) % inv)
            call da_random_omb(iv % sound(n) % t(k) % error, ob % sound(n) % t(k), &
                              iv % sound(n) % t(k) % qc, iv % sound(n) % t(k) % inv)
            call da_random_omb(iv % sound(n) % q(k) % error, ob % sound(n) % q(k), &
                              iv % sound(n) % q(k) % qc, iv % sound(n) % q(k) % inv)
         end do
      end do
   end if

   if (iv%info(sonde_sfc)%nlocal > 0) then
      do n = 1, iv%info(sonde_sfc)%nlocal
         call da_random_omb(iv % sonde_sfc(n) % u % error, ob % sonde_sfc(n) % u, &
                             iv % sonde_sfc(n) % u % qc, iv % sonde_sfc(n) % u % inv)
         call da_random_omb(iv % sonde_sfc(n) % v % error, ob % sonde_sfc(n) % v, &
                             iv % sonde_sfc(n) % v % qc, iv % sonde_sfc(n) % v % inv)
         call da_random_omb(iv % sonde_sfc(n) % t % error, ob % sonde_sfc(n) % t, &
                             iv % sonde_sfc(n) % t % qc, iv % sonde_sfc(n) % t % inv )
         call da_random_omb(iv % sonde_sfc(n) % p % error, ob % sonde_sfc(n) % p, &
                             iv % sonde_sfc(n) % p % qc, iv % sonde_sfc(n) % p % inv)
         call da_random_omb(iv % sonde_sfc(n) % q % error, ob % sonde_sfc(n) % q, &
                             iv % sonde_sfc(n) % q % qc, iv % sonde_sfc(n) % q % inv)
      end do
   end if

   ! [2.7] Transfer airep obs:

   if (iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
         do k = 1, iv%info(airep)%levels(n)
            call da_random_omb(iv % airep(n) % u(k) % error, ob % airep(n) % u(k), &
                              iv % airep(n) % u(k) % qc, iv % airep(n) % u(k) % inv)
            call da_random_omb(iv % airep(n) % v(k) % error, ob % airep(n) % v(k), &
                              iv % airep(n) % v(k) % qc, iv % airep(n) % v(k) % inv)
            call da_random_omb(iv % airep(n) % t(k) % error, ob % airep(n) % t(k), &
                              iv % airep(n) % t(k) % qc, iv % airep(n) % t(k) % inv)
            call da_random_omb(iv % airep(n) % q(k) % error, ob % airep(n) % q(k), &
                              iv % airep(n) % q(k) % qc, iv % airep(n) % q(k) % inv)
         end do
      end do
   end if

   ! [2.8] Transfer pilot obs:

   if (iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
         do k = 1, iv%info(pilot)%levels(n)
            call da_random_omb(iv % pilot(n) % u(k) % error, ob % pilot(n) % u(k), &
                              iv % pilot(n) % u(k) % qc, iv % pilot(n) % u(k) % inv)
            call da_random_omb(iv % pilot(n) % v(k) % error, ob % pilot(n) % v(k), &
                              iv % pilot(n) % v(k) % qc, iv % pilot(n) % v(k) % inv)
         end do
      end do
   end if

   ! [2.9] Transfer SSM/I obs:SSMI:

   if (iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
         call da_random_omb(iv % ssmi_rv(n) % speed % error, &
                           ob % ssmi_rv(n) % speed, &
                           iv % ssmi_rv(n) % speed % qc, &
                           iv % ssmi_rv(n) % speed % inv)
         call da_random_omb(iv % ssmi_rv(n) % tpw % error, &
                           ob % ssmi_rv(n) % tpw, &
                           iv % ssmi_rv(n) % tpw % qc, &
                           iv % ssmi_rv(n) % tpw % inv)
      end do
   end if

   if (iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
         call da_random_omb(iv % ssmi_tb(n) % tb19h % error, &
                           ob % ssmi_tb(n) % tb19h, &
                           iv % ssmi_tb(n) % tb19h % qc, &
                           iv % ssmi_tb(n) % tb19h % inv)
         call da_random_omb(iv % ssmi_tb(n) % tb19v % error, &
                           ob % ssmi_tb(n) % tb19v, &
                           iv % ssmi_tb(n) % tb19v % qc, &
                           iv % ssmi_tb(n) % tb19v % inv)
         call da_random_omb(iv % ssmi_tb(n) % tb22v % error, &
                           ob % ssmi_tb(n) % tb22v, &
                           iv % ssmi_tb(n) % tb22v % qc, &
                           iv % ssmi_tb(n) % tb22v % inv)
         call da_random_omb(iv % ssmi_tb(n) % tb37h % error, &
                           ob % ssmi_tb(n) % tb37h, &
                           iv % ssmi_tb(n) % tb37h % qc, &
                           iv % ssmi_tb(n) % tb37h % inv)
         call da_random_omb(iv % ssmi_tb(n) % tb37v % error, &
                           ob % ssmi_tb(n) % tb37v, &
                           iv % ssmi_tb(n) % tb37v % qc, &
                           iv % ssmi_tb(n) % tb37v % inv)
         call da_random_omb(iv % ssmi_tb(n) % tb85h % error, &
                           ob % ssmi_tb(n) % tb85h, &
                           iv % ssmi_tb(n) % tb85h % qc, &
                           iv % ssmi_tb(n) % tb85h % inv)
         call da_random_omb(iv % ssmi_tb(n) % tb85v % error, &
                           ob % ssmi_tb(n) % tb85v, &
                           iv % ssmi_tb(n) % tb85v % qc, &
                           iv % ssmi_tb(n) % tb85v % inv)
      end do
   end if

   ! [2.10] Transfer satem obs:

   if (iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
         do k = 1, iv%info(satem)%levels(n)
            call da_random_omb(iv % satem(n) % thickness(k) % error, &
                              ob % satem(n) % thickness(k), &
                              iv % satem(n) % thickness(k) % qc, &
                              iv % satem(n) % thickness(k) % inv)
         end do
      end do
   end if
   
   ! [2.11] Transfer ssmt1 obs:

   if (iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
         do k = 1, iv%info(ssmt1)%levels(n)
            call da_random_omb(iv % ssmt1(n) % t(k) % error, &
                              ob % ssmt1(n) % t(k), &
                              iv % ssmt1(n) % t(k) % qc, &
                              iv % ssmt1(n) % t(k) % inv)
         end do
      end do
   end if

   ! [2.12] Transfer ssmt2 obs:

   if (iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
         do k = 1, iv%info(ssmt2)%levels(n)
            call da_random_omb(iv % ssmt2(n) % rh(k) % error, &
                              ob % ssmt2(n) % rh(k), &
                              iv % ssmt2(n) % rh(k) % qc, &
                              iv % ssmt2(n) % rh(k) % inv)
         end do
      end do
   end if
   
   ! [2.13] Transfer scatterometer obs:

   if (iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
         call da_random_omb(iv % qscat(n) % u % error, ob % qscat(n) % u, &
                           iv % qscat(n) % u % qc, iv % qscat(n) % u % inv)
         call da_random_omb(iv % qscat(n) % v % error, ob % qscat(n) % v, &
                           iv % qscat(n) % v % qc, iv % qscat(n) % v % inv)
      end do
   end if

   ! [2.14] Transfer buoy obs:

   if (iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
         call da_random_omb(iv % buoy(n) % u % error, ob % buoy(n) % u, &
                           iv % buoy(n) % u % qc, iv % buoy(n) % u % inv)
         call da_random_omb(iv % buoy(n) % v % error, ob % buoy(n) % v, &
                           iv % buoy(n) % v % qc, iv % buoy(n) % v % inv)
         call da_random_omb(iv % buoy(n) % t % error, ob % buoy(n) % t, &
                           iv % buoy(n) % t % qc, iv % buoy(n) % t % inv)
         call da_random_omb(iv % buoy(n) % p % error, ob % buoy(n) % p, &
                           iv % buoy(n) % p % qc, iv % buoy(n) % p % inv)
         call da_random_omb(iv % buoy(n) % q % error, ob % buoy(n) % q, &
                           iv % buoy(n) % q % qc, iv % buoy(n) % q % inv)
      end do
   end if

   ! [2.15] Transfer profiler obs:

   if (iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
         do k = 1, iv%info(profiler)%levels(n)
            call da_random_omb(iv % profiler(n) % u(k) % error, &
               ob % profiler(n) % u(k), &
               iv % profiler(n) % u(k) % qc, iv % profiler(n) % u(k) % inv)
            call da_random_omb(iv % profiler(n) % v(k) % error, &
               ob % profiler(n) % v(k), &
               iv % profiler(n) % v(k) % qc, iv % profiler(n) % v(k) % inv)
         end do
      end do
   end if

   ! [2.16] Transfer TC bogus obs:

   if (iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
         do k = 1, iv%info(bogus)%levels(n)
            call da_random_omb(iv % bogus(n) % u(k) % error, &
               ob % bogus(n) % u(k), &
               iv % bogus(n) % u(k) % qc, iv % bogus(n) % u(k) % inv)
            call da_random_omb(iv % bogus(n) % v(k) % error, &
               ob % bogus(n) % v(k), &
               iv % bogus(n) % v(k) % qc, iv % bogus(n) % v(k) % inv)
            call da_random_omb(iv % bogus(n) % t(k) % error, &
               ob % bogus(n) % t(k), &
               iv % bogus(n) % t(k) % qc, iv % bogus(n) % t(k) % inv)
            call da_random_omb(iv % bogus(n) % q(k) % error, &
              ob % bogus(n) % q(k), &
              iv % bogus(n) % q(k) % qc, iv % bogus(n) % q(k) % inv)
         end do

         call da_random_omb(iv % bogus(n) % slp % error, ob % bogus(n) % slp, &
                            iv % bogus(n) % slp % qc, iv % bogus(n) % slp % inv)
      end do
   end if

   ! Transfer AIRS retrievals:

   if (iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
         do k = 1, iv%info(airsr)%levels(n)
            call da_random_omb(iv % airsr(n) % t(k) % error, &
               ob % airsr(n) % t(k), &
               iv % airsr(n) % t(k) % qc, iv % airsr(n) % t(k) % inv)
            call da_random_omb(iv % airsr(n) % q(k) % error, &
               ob % airsr(n) % q(k), &
               iv % airsr(n) % q(k) % qc, iv % airsr(n) % q(k) % inv)
         end do
      end do
   end if

   if (iv%info(mtgirs)%nlocal > 0) then
      do n = 1, iv%info(mtgirs)%nlocal
         do k = 1, iv%info(mtgirs)%levels(n)
            call da_random_omb(iv % mtgirs(n) % u(k) % error, ob %  mtgirs(n) % u(k), &
                              iv %  mtgirs(n) % u(k) % qc, iv %  mtgirs(n) % u(k) % inv)
            call da_random_omb(iv %  mtgirs(n) % v(k) % error, ob %  mtgirs(n) % v(k), &
                              iv %  mtgirs(n) % v(k) % qc, iv %  mtgirs(n) % v(k) % inv)
            call da_random_omb(iv %  mtgirs(n) % t(k) % error, ob %  mtgirs(n) % t(k), &
                              iv %  mtgirs(n) % t(k) % qc, iv %  mtgirs(n) % t(k) % inv)
            call da_random_omb(iv %  mtgirs(n) % q(k) % error, ob %  mtgirs(n) % q(k), &
                              iv %  mtgirs(n) % q(k) % qc, iv %  mtgirs(n) % q(k) % inv)
         end do
      end do
    end if

   if (iv%info(tamdar)%nlocal > 0) then
      do n = 1, iv%info(tamdar)%nlocal
         do k = 1, iv%info(tamdar)%levels(n)
            call da_random_omb(iv % tamdar(n) % u(k) % error, ob %  tamdar(n) % u(k), &
                              iv %  tamdar(n) % u(k) % qc, iv %  tamdar(n) % u(k) % inv)
            call da_random_omb(iv %  tamdar(n) % v(k) % error, ob %  tamdar(n) % v(k), &
                              iv %  tamdar(n) % v(k) % qc, iv %  tamdar(n) % v(k) % inv)
            call da_random_omb(iv %  tamdar(n) % t(k) % error, ob %  tamdar(n) % t(k), &
                              iv %  tamdar(n) % t(k) % qc, iv %  tamdar(n) % t(k) % inv)
            call da_random_omb(iv %  tamdar(n) % q(k) % error, ob %  tamdar(n) % q(k), &
                              iv %  tamdar(n) % q(k) % qc, iv %  tamdar(n) % q(k) % inv)
         end do
      end do
   end if

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      do n = 1, iv%info(tamdar_sfc)%nlocal
         call da_random_omb(iv % tamdar_sfc(n) % u % error, ob % tamdar_sfc(n) % u, &
                             iv % tamdar_sfc(n) % u % qc, iv % tamdar_sfc(n) % u % inv)
         call da_random_omb(iv % tamdar_sfc(n) % v % error, ob % tamdar_sfc(n) % v, &
                             iv % tamdar_sfc(n) % v % qc, iv % tamdar_sfc(n) % v % inv)
         call da_random_omb(iv % tamdar_sfc(n) % t % error, ob % tamdar_sfc(n) % t, &
                             iv % tamdar_sfc(n) % t % qc, iv % tamdar_sfc(n) % t % inv )
         call da_random_omb(iv % tamdar_sfc(n) % p % error, ob % tamdar_sfc(n) % p, &
                             iv % tamdar_sfc(n) % p % qc, iv % tamdar_sfc(n) % p % inv)
         call da_random_omb(iv % tamdar_sfc(n) % q % error, ob % tamdar_sfc(n) % q, &
                             iv % tamdar_sfc(n) % q % qc, iv % tamdar_sfc(n) % q % inv)

      end do
    end if

   if (trace_use) call da_trace_exit("da_random_omb_all")

end subroutine da_random_omb_all


subroutine da_setup_pseudo_obs(iv, ob)

   !-------------------------------------------------------------------------
   ! Purpose: Sets up pseudo ob part of observation structure.
   !-------------------------------------------------------------------------

   implicit none

   type(iv_type),    intent(inout) :: iv   ! Obs and header structure.
   type(y_type),     intent(inout) :: ob   ! (Smaller) observation structure.

   integer                       :: n    ! Loop counters.

   if (trace_use_dull) call da_trace_entry("da_setup_pseudo_obs")

   allocate (iv%pseudo(1:iv%info(pseudo)%nlocal))
   !do n=1, iv%info(pseudo)%nlocal
      iv%pseudo(:) % u % inv = missing_r
      iv%pseudo(:) % v % inv = missing_r
      iv%pseudo(:) % t % inv = missing_r
      iv%pseudo(:) % p % inv = missing_r
      iv%pseudo(:) % q % inv = missing_r

      iv%pseudo(:) % u % error = missing_r
      iv%pseudo(:) % v % error = missing_r
      iv%pseudo(:) % t % error = missing_r
      iv%pseudo(:) % p % error = missing_r
      iv%pseudo(:) % q % error = missing_r

      iv%pseudo(:) % u % qc  = missing_data
      iv%pseudo(:) % v % qc  = missing_data
      iv%pseudo(:) % t % qc  = missing_data
      iv%pseudo(:) % p % qc  = missing_data
      iv%pseudo(:) % q % qc  = missing_data

      ob%pseudo(:) % u = missing_r
      ob%pseudo(:) % v = missing_r
      ob%pseudo(:) % t = missing_r
      ob%pseudo(:) % p = missing_r
      ob%pseudo(:) % q = missing_r

      !---------------------------------------------------------------
      ! [1.0] Initialise components of innovation vector:
      !---------------------------------------------------------------

      iv%info(pseudo)%x(:,:)  = pseudo_x
      iv%info(pseudo)%y(:,:)  = pseudo_y
      iv%info(pseudo)%zk(:,:) = pseudo_z

      iv%info(pseudo)%i(:,:) = int(pseudo_x)
      iv%info(pseudo)%j(:,:) = int(pseudo_y)
      iv%info(pseudo)%k(:,:) = int(pseudo_z)

      iv%info(pseudo)%dx(:,:) = pseudo_x-real(iv%info(pseudo)%i(:,:))
      iv%info(pseudo)%dy(:,:) = pseudo_y-real(iv%info(pseudo)%j(:,:))
      iv%info(pseudo)%dxm(:,:)=1.0-iv%info(pseudo)%dx(:,:)
      iv%info(pseudo)%dym(:,:)=1.0-iv%info(pseudo)%dy(:,:)
      iv%info(pseudo)%levels(:) = 1


      if (pseudo_var(1:1) == 'u' .or. pseudo_var(1:1) == 'U') then
         iv%pseudo(:) % u % inv = pseudo_val
         iv%pseudo(:) % u % error = pseudo_err
         iv%pseudo(:) % u % qc = 0
      else if (pseudo_var(1:1) == 'v' .or. pseudo_var(1:1) == 'V') then
         iv%pseudo(:) % v % inv = pseudo_val
         iv%pseudo(:) % v % error = pseudo_err
         iv%pseudo(:) % v % qc = 0
      else if (pseudo_var(1:1) == 't' .or. pseudo_var(1:1) == 'T') then
         iv%pseudo(:) % t % inv = pseudo_val
         iv%pseudo(:) % t % error = pseudo_err
         iv%pseudo(:) % t % qc = 0
      else if (pseudo_var(1:1) == 'p' .or. pseudo_var(1:1) == 'P') then
         iv%pseudo(:) % p % inv = pseudo_val
         iv%pseudo(:) % p % error = pseudo_err
         iv%pseudo(:) % p % qc = 0
      else if (pseudo_var(1:1) == 'q' .or. pseudo_var(1:1) == 'Q') then
         iv%pseudo(:) % q % inv = pseudo_val
         iv%pseudo(:) % q % error = pseudo_err
         iv%pseudo(:) % q % qc = 0
      end if 
      
!      write(unit=stdout,fmt='(a4,2f15.5)')pseudo_var, pseudo_val, pseudo_err
!      write(unit=stdout,fmt='(3f15.5)')pseudo_x, pseudo_y, pseudo_z
   !end do

   if (trace_use_dull) call da_trace_exit("da_setup_pseudo_obs")

end subroutine da_setup_pseudo_obs


subroutine da_store_obs_grid_info (info)

   !-----------------------------------------------------------------------
   ! Purpose: this is in parallel of da_store_obs_grid_info_rad but with
   !          an extra thinned check to decide proc_domain.
   !-----------------------------------------------------------------------

   implicit none

   type(infa_type), intent(inout) :: info

   integer :: n

   if (trace_use) call da_trace_entry("da_store_obs_grid_info")

   info%proc_domain(:,:) = .false.

   do n=1,info%nlocal
      if (info%i(1,n) >= its .and. info%i(1,n) <= ite .and. info%j(1,n) >= jts .and. info%j(1,n) <= jte) then
         if ( .not. info%thinned(1,n) ) then
            info%proc_domain(:,n) = .true.
         end if
      end if
   end do

   if (trace_use) call da_trace_exit("da_store_obs_grid_info")

end subroutine da_store_obs_grid_info


subroutine da_store_obs_grid_info_rad (info)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(infa_type), intent(inout) :: info

   integer :: n

   if (trace_use) call da_trace_entry("da_store_obs_grid_info_rad")

   info%proc_domain(:,:) = .false.
   
   do n=1,info%nlocal
      if (info%i(1,n) >= its .and. info%i(1,n) <= ite .and. info%j(1,n) >= jts .and. info%j(1,n) <= jte) then
         info%proc_domain(:,n) = .true.
      end if
   end do

   if (trace_use) call da_trace_exit("da_store_obs_grid_info_rad")

end subroutine da_store_obs_grid_info_rad


subroutine da_count_filtered_obs (ptop, map, ds, phic, xlonc, truelat1, truelat2, &
   coarse_ix, coarse_jy, start_x, start_y)

   !---------------------------------------------------------------------------
   ! Purpose: Scans intermediate Filtered Obs file, 
   !          counts various obs type and writes on filtered_obs_unit
   !---------------------------------------------------------------------------

   implicit none

   real,              intent(in) :: ptop, ds
   real,              intent(in) :: phic, xlonc, truelat1, truelat2
   integer,           intent(in) :: coarse_ix, coarse_jy
   real,              intent(in) :: start_x, start_y
   integer,           intent(in) :: map

   integer                      :: i, iost, fm
   type (multi_level_type)      :: platform
   real                         :: height_error
   integer                      :: nlocal(num_ob_indexes)
   integer                      :: num_others

   integer                        :: iunit, files, total_obs     
   integer                        :: maxnes, numc, idd
   character(len=filename_len)    :: filename


   if (trace_use) call da_trace_entry("da_count_filtered_obs")

   nlocal(:) = 0
   num_others = 0

   call da_get_unit(iunit)

   ! Loop over all data files
   do files = 0, num_procs-1
      write(unit=filename, fmt='(a,i4.4)') 'filtered_obs.',files  
      open(unit=iunit, file= trim(filename), form='formatted',iostat=iost)
      if (iost /= 0) call da_error("da_count_filtered_obs.inc",39, (/"Cannot open "//filename/))

      !  loop over records

      reports: do
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
            !write (0,'(/A,I9)') ' end OF OBS unit: ',iunit
            !write (0,'(/A,I9)') ' iostat:          ',iost
            exit reports
         end if

         ! Read surface Info

         read (iunit, fmt = fmt_srfc)  &
            platform%loc%slp%inv, platform%loc%slp%qc, &
            platform%loc%slp%error,                    &
            platform%loc%pw%inv, platform%loc%pw%qc,   &
            platform%loc%pw%error
         ! levels < 1 and not GPSPW, go back to reports

         if ((platform%info%levels < 1) .AND.            &
             (index(platform%info%platform, 'GPSPW') <= 0)) then
              cycle reports
         end if
         read(platform%info%platform(4:6), '(I3)') fm

         ! read each level
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
               platform%each(i)%speed%inv, platform%each(i)%speed%qc,                   &
               platform%each(i)%speed%error,                                            &
               platform%each(i)%v%inv, platform%each(i)%v%qc, platform%each(i)%v%error, &
               platform%each(i)%height,                                                 &
               platform%each(i)%height_qc,                                              &
               height_error,                                                            &
               platform%each(i)%t%inv, platform%each(i)%t%qc, platform%each(i)%t%error, &
               platform%each(i)%td%inv, platform%each(i)%td%qc, platform%each(i)%td%error, &
               platform%each(i)%rh%inv, platform%each(i)%rh%qc, platform%each(i)%rh%error
         end do

         if (platform%info%levels < 1) then
            if (fm /= 111) then
               cycle reports
            end if
         end if
         select case(fm)

         case (12) ;

            if (.not.use_synopobs) cycle reports
            nlocal(synop) = nlocal(synop) + 1

         case (13, 17) ;                  ! ships          

            if (.not.use_shipsobs) cycle reports
            nlocal(ships)  = nlocal(ships)  + 1

         case (15:16) ;

            if (.not.use_metarobs) cycle reports
            nlocal(metar) = nlocal(metar) + 1

         case (32:34) ;

            if (.not.use_pilotobs) cycle reports
            nlocal(pilot) = nlocal(pilot) + 1

         case (35:38) ;
            if (.not.use_soundobs) cycle reports
            nlocal(sound) = nlocal(sound) + 1
            nlocal(sonde_sfc) = nlocal(sonde_sfc) + 1

         case (161) ;
            if (.not.use_mtgirsobs) cycle reports
            nlocal(mtgirs) = nlocal(mtgirs) + 1

         case (86) ;

            if (.not.use_satemobs) cycle reports

            ! Reject cloudy satem obs.
 
            if (platform%loc%pw%inv > 10.0) then
               cycle reports
            end if

            nlocal(satem) = nlocal(satem) + 1

         case (88)    ;
            ! Geostationary or Polar orbitting Satellite AMVs:

            if (index(platform%info%name, 'MODIS') > 0 .or. &
                index(platform%info%name, 'modis') > 0 .or. &
                index(platform%info%id, 'AVHRR') > 0)  then
               if (.not.use_polaramvobs) cycle reports
               nlocal(polaramv) = nlocal(polaramv) + 1
            else
               if (.not.use_geoamvobs) cycle reports 
               nlocal(geoamv) = nlocal(geoamv) + 1
            end if

         case (42,96:97) ;

            if (.not.use_airepobs) cycle reports
            nlocal(airep) = nlocal(airep) + 1

         case (101) ;
            if (.not.use_tamdarobs) cycle reports
            nlocal(tamdar) = nlocal(tamdar) + 1
            nlocal(tamdar_sfc) = nlocal(tamdar_sfc) + 1

         case (111) ;
         
            if (.not.use_gpspwobs) cycle reports
            nlocal(gpspw) = nlocal(gpspw) + 1

         case (116) ;
         
            if (.not.use_gpsrefobs) cycle reports
            nlocal(gpsref) = nlocal(gpsref) + 1

         case (121) ;
            ! SSM/T1 temperatures:

            if (.not.use_ssmt1obs) cycle reports
            nlocal(ssmt1) = nlocal(ssmt1) + 1

         case (122) ;
            ! SSM/T2 relative humidities:

            if (.not.use_ssmt2obs) cycle reports
            nlocal(ssmt2) = nlocal(ssmt2) + 1

         case (281)    ;
            ! Scatterometer:

            if (.not.use_qscatobs) cycle reports
            nlocal(qscat)  = nlocal(qscat)  + 1

         case (132) ;

            if (.not.use_profilerobs) cycle reports
            nlocal(profiler) = nlocal(profiler) + 1

         case (135) ;

            if (.not.use_bogusobs) cycle reports
            nlocal(bogus) = nlocal(bogus) + 1

         case (18,19) ;             ! bouy

            if (.not.use_buoyobs) cycle reports
            nlocal(buoy)  = nlocal(buoy)  + 1

         case (133) ;      !  AIRS retrievals
            if (.not.use_airsretobs) cycle reports
            nlocal(airsr) = nlocal(airsr) + 1

         case default;
            num_others = num_others + 1
            write(unit=message(1), fmt='(a)') 'unsaved obs found:'
            write(unit=message(2), fmt='(2a)') &
               'platform%info%platform=', platform%info%platform
            write(unit=message(3), fmt='(a, i3)') &
               'platform%info%levels=', platform%info%levels
            call da_warning("da_count_filtered_obs.inc",227,message(1:3))
         end select
      end do reports                  !  Loop over reports              
      close (iunit)
   end do               !  Loop over all data files
   call da_free_unit (iunit)

   ! write counts in the header 

   total_obs = nlocal(synop) + nlocal(metar) + nlocal(ships) + &
      nlocal(buoy) + nlocal(sound) + nlocal(sonde_sfc) + nlocal(airep) + nlocal(pilot) + &
      nlocal(geoamv) + nlocal(polaramv) + nlocal(gpspw) + nlocal(gpsref) + &
      nlocal(profiler) + nlocal(qscat) + nlocal(ssmt1) + nlocal(ssmt2) +  &
      nlocal(satem)  + nlocal(bogus) +  nlocal(airsr) + nlocal(mtgirs) + nlocal(tamdar) + nlocal(tamdar_sfc) + num_others

   write  (unit = filtered_obs_unit, fmt = '(A,I7,A,F8.0,A)')    &
      "TOTAL =",total_obs, ", MISS. =",missing_r,","  

   ! Write other counts       

   write  (unit = filtered_obs_unit, fmt = '(6(A,I7,A))')    &
      "SYNOP =",nlocal(synop),", ", &
      "METAR =",nlocal(metar),", ", &
      "SHIP  =",nlocal(ships),", ", &
      "BUOY  =",nlocal(buoy),", ", &
      "TEMP  =",nlocal(sound),", ", &
      "TEMP_SFC  =",nlocal(sonde_sfc),", ", &
      "AIREP =",nlocal(airep),", ", &
      "PILOT =",nlocal(pilot),", ", &
      "GeAMV =",nlocal(geoamv),", ", &
      "PoAMV =",nlocal(polaramv),", ", &
      "GPSPW =",nlocal(gpspw),", ", &
      "GPSRF =",nlocal(gpsref),", ", &
      "PROFL =",nlocal(profiler),", ", &
      "QSCAT =",nlocal(qscat),", ", &
      "SSMT1 =",nlocal(ssmt1),", ", &
      "SSMT2 =",nlocal(ssmt2),", ", &
      "SATEM =",nlocal(satem), ", ", &
      "BOGUS =",nlocal(bogus),", ", &
      "AIRSR =",nlocal(airsr),", ", &
      "MTGIRS=",nlocal(mtgirs),", ", &
      "TAMDAR=",nlocal(tamdar),", ", &
      "TAMDAR_SFC=",nlocal(tamdar_sfc),", ", &
      "OTHER =",num_others,", "

   ! write reference state info

   write (unit = filtered_obs_unit, fmt = '(A,F7.2,A,F7.2,4(A,F7.2),A)') &
      "PHIC  =", phic,", XLONC =", xlonc,", TRUE1 =", truelat1,&
      ", TRUE2 =",truelat2, ", XIM11 =", start_y, ", XJM11 =", start_x, ","
   write (unit = filtered_obs_unit, fmt = '(2(A,F7.2),A,F7.0,A,F7.0,A)') &
      "TS0   =",  base_temp, ", TLP   =", base_lapse, &
      ", PTOP  =",  ptop,", PS0   =",  base_pres,"," 

   ! write domain info 

   !  hardwire following variables
      maxnes = 2; numc = 1
   if ( coarse_ix == ide+1 .and. coarse_jy == jde+1 ) then
      idd = 2
   else
      idd = 1
   end if
   write (unit = filtered_obs_unit, fmt = '(5(A,I7),A)') &
      "IXC   =", coarse_jy   ,", JXC   =", coarse_ix   ,", IPROJ =", map,&
      ", IDD   =", idd,", MAXNES=",maxnes,","
   write (unit = filtered_obs_unit, fmt = '(A,10(:,I7,A))')  &
      "NESTIX=", coarse_jy, ", ", jde+1, ", "
   write (unit = filtered_obs_unit, fmt = '(A,10(:,I7,A))')  &
      "NESTJX=", coarse_ix, ", ", ide+1, ", "
   write (unit = filtered_obs_unit, fmt = '(A,10(:,I7,A))')  &
      "NUMC  =",(numc ,", ", i = 1, maxnes)
   write (unit = filtered_obs_unit, fmt = '(A,10(:,F7.2,A))')&
      "DIS   =",(ds/1000.0,", ",i = 1, maxnes)
   write (unit = filtered_obs_unit, fmt = '(A,10(:,I7,A))')  &
      "NESTI=", 1, ", ", int(start_y), ", "
   write (unit = filtered_obs_unit, fmt = '(A,10(:,I7,A))')  &
      "NESTJ=", 1, ", ", int(start_x), ", "

   ! write variable names and unit

   write (unit = filtered_obs_unit, fmt = '(A)') &
      "INFO  = PLATFORM, DATE, NAME, LEVELS, LATITUDE, LONGITUDE, ELEVATION, ID."
   write (unit = filtered_obs_unit, fmt = '(A)') &
      "SRFC  = SLP, PW (DATA,QC,ERROR)."
   write (unit = filtered_obs_unit, fmt = '(A)') &
      "EACH  = PRES, SPEED, DIR, HEIGHT, TEMP, DEW PT, HUMID (DATA,QC,ERROR)*LEVELS."

   ! write format info
 
   write (unit = filtered_obs_unit, fmt = '(2A)') 'INFO_fmt = ', trim (fmt_info)
   write (unit = filtered_obs_unit, fmt = '(2A)') 'SRFC_fmt = ', trim (fmt_srfc)
   write (unit = filtered_obs_unit, fmt = '(2A)') 'EACH_fmt = ', trim (fmt_each)

   ! write end of header record

   if (write_mod_filtered_obs) then
      write (unit = filtered_obs_unit, fmt = '(A)') &
         "MODIFIED FILTERED OBS #--------------------------------------------------------#"
   else
      write (unit = filtered_obs_unit, fmt = '(A)') &
         "FILTERED OBS #-----------------------------------------------------------------#"
   end if     

   if (trace_use) call da_trace_exit("da_count_filtered_obs")

end subroutine da_count_filtered_obs 


subroutine da_obs_sensitivity(ktr, iv)

   !-------------------------------------------------------------------------
   ! Purpose:        Apply R-1 to KT.dF/dx to obtain Observation sensitivity
   !                 stored in "iv" for every observation
   !
   ! Called from da_minimise_lz
   !
   ! History: 12/12/2008  Creation (Tom Auligne)
   !
   !-------------------------------------------------------------------------

   implicit none

   type (y_type), intent(in)         :: ktr         
   type (iv_type), intent(inout)     :: iv      
   
   integer, parameter                :: nbvar = 7
   integer, parameter                :: npredmax = 8

   integer, parameter                :: u     = 1
   integer, parameter                :: v     = 2
   integer, parameter                :: t     = 3
   integer, parameter                :: p     = 4
   integer, parameter                :: q     = 5
   integer, parameter                :: rv    = 6
   integer, parameter                :: mx    = 7
   character(len=6), parameter       :: var_names(nbvar) = (/ 'U     ', &
                                                              'V     ', &
                                                              'T     ', &
                                                              'P     ', &
                                                              'Q     ', &
                                                              'GPSREF', &
                                                              'SATEM ' /)

   integer                           :: imsg
   integer                           :: i, j, k, n, inst
   integer                           :: npred, ipred
   real                              :: ktd(num_ob_indexes,nbvar)
   real                              :: ktdrad(iv%num_inst)
   real                              :: ktbrad(iv%num_inst,npredmax) 
   real                              :: ktd_global(num_ob_indexes,nbvar)
   real                              :: ktdrad_global(iv%num_inst)
   real                              :: ktbrad_global(iv%num_inst, npredmax) 

   integer                           :: iunit
   if (trace_use) call da_trace_entry("da_obs_sensitivity")

   ktd    = 0.0
   ktdrad = 0.0
   ktbrad = 0.0

   ktdrad_global = 0.0
   ktbrad_global = 0.0

   call da_get_unit(iunit)
       
        if (num_pseudo > 0) then
            do n=1, iv%info(pseudo)%nlocal
	       if (.not. iv%info(pseudo)%proc_domain(1,n)) cycle
	       iv%pseudo(n)%u%sens = ktr%pseudo(n)%u / (iv%pseudo(n)%u%error **2)
	       iv%pseudo(n)%v%sens = ktr%pseudo(n)%v / (iv%pseudo(n)%v%error **2)
	       iv%pseudo(n)%t%sens = ktr%pseudo(n)%t / (iv%pseudo(n)%t%error **2)
	       iv%pseudo(n)%p%sens = ktr%pseudo(n)%p / (iv%pseudo(n)%p%error **2)
	       iv%pseudo(n)%q%sens = ktr%pseudo(n)%q / (iv%pseudo(n)%q%error **2)	   
	       
	       iv%pseudo(n)%u%inv  = iv%pseudo(n)%u%inv * iv%pseudo(n)%u%sens
	       iv%pseudo(n)%v%inv  = iv%pseudo(n)%v%inv * iv%pseudo(n)%v%sens
	       iv%pseudo(n)%t%inv  = iv%pseudo(n)%t%inv * iv%pseudo(n)%t%sens
	       iv%pseudo(n)%p%inv  = iv%pseudo(n)%p%inv * iv%pseudo(n)%p%sens
	       iv%pseudo(n)%q%inv  = iv%pseudo(n)%q%inv * iv%pseudo(n)%q%sens	  
	       
	       ktd(pseudo,u) = ktd(pseudo,u) + iv%pseudo(n)%u%inv
               ktd(pseudo,v) = ktd(pseudo,v) + iv%pseudo(n)%v%inv 
	       ktd(pseudo,t) = ktd(pseudo,t) + iv%pseudo(n)%t%inv 
               ktd(pseudo,p) = ktd(pseudo,p) + iv%pseudo(n)%p%inv 
	       ktd(pseudo,q) = ktd(pseudo,q) + iv%pseudo(n)%q%inv	       	             	       	            	       
	    end do
         end if

         if (iv%info(synop)%nlocal > 0) then
            do n=1, iv%info(synop)%nlocal
	       if (.not. iv%info(synop)%proc_domain(1,n)) cycle
	       iv%synop(n)%u%sens = ktr%synop(n)%u / (iv%synop(n)%u%error **2)
	       iv%synop(n)%v%sens = ktr%synop(n)%v / (iv%synop(n)%v%error **2)
	       iv%synop(n)%t%sens = ktr%synop(n)%t / (iv%synop(n)%t%error **2)
	       iv%synop(n)%p%sens = ktr%synop(n)%p / (iv%synop(n)%p%error **2)
	       iv%synop(n)%q%sens = ktr%synop(n)%q / (iv%synop(n)%q%error **2)	        	       
	       
	       iv%synop(n)%u%inv  = iv%synop(n)%u%inv * iv%synop(n)%u%sens
	       iv%synop(n)%v%inv  = iv%synop(n)%v%inv * iv%synop(n)%v%sens
	       iv%synop(n)%t%inv  = iv%synop(n)%t%inv * iv%synop(n)%t%sens
	       iv%synop(n)%p%inv  = iv%synop(n)%p%inv * iv%synop(n)%p%sens
	       iv%synop(n)%q%inv  = iv%synop(n)%q%inv * iv%synop(n)%q%sens	  
	       
	       ktd(synop,u) = ktd(synop,u) + iv%synop(n)%u%inv
               ktd(synop,v) = ktd(synop,v) + iv%synop(n)%v%inv 
	       ktd(synop,t) = ktd(synop,t) + iv%synop(n)%t%inv 
               ktd(synop,p) = ktd(synop,p) + iv%synop(n)%p%inv 
	       ktd(synop,q) = ktd(synop,q) + iv%synop(n)%q%inv	       	             	       	            	       
	    end do
	 end if   
	       
         if (iv%info(ships)%nlocal > 0) then
            do n=1, iv%info(ships)%nlocal
	       if (.not. iv%info(ships)%proc_domain(1,n)) cycle
	       iv%ships(n)%u%sens = ktr%ships(n)%u / (iv%ships(n)%u%error **2)
	       iv%ships(n)%v%sens = ktr%ships(n)%v / (iv%ships(n)%v%error **2)
	       iv%ships(n)%t%sens = ktr%ships(n)%t / (iv%ships(n)%t%error **2)
	       iv%ships(n)%p%sens = ktr%ships(n)%p / (iv%ships(n)%p%error **2)
	       iv%ships(n)%q%sens = ktr%ships(n)%q / (iv%ships(n)%q%error **2)	        	       
	       
	       iv%ships(n)%u%inv  = iv%ships(n)%u%inv * iv%ships(n)%u%sens
	       iv%ships(n)%v%inv  = iv%ships(n)%v%inv * iv%ships(n)%v%sens
	       iv%ships(n)%t%inv  = iv%ships(n)%t%inv * iv%ships(n)%t%sens
	       iv%ships(n)%p%inv  = iv%ships(n)%p%inv * iv%ships(n)%p%sens
	       iv%ships(n)%q%inv  = iv%ships(n)%q%inv * iv%ships(n)%q%sens	  
	       
	       ktd(ships,u) = ktd(ships,u) + iv%ships(n)%u%inv
               ktd(ships,v) = ktd(ships,v) + iv%ships(n)%v%inv 
	       ktd(ships,t) = ktd(ships,t) + iv%ships(n)%t%inv 
               ktd(ships,p) = ktd(ships,p) + iv%ships(n)%p%inv 
	       ktd(ships,q) = ktd(ships,q) + iv%ships(n)%q%inv	       	             	       	            	       
	    end do
	 end if   
	       
         if (iv%info(metar)%nlocal > 0) then
            do n=1, iv%info(metar)%nlocal
	       if (.not. iv%info(metar)%proc_domain(1,n)) cycle
	       iv%metar(n)%u%sens = ktr%metar(n)%u / (iv%metar(n)%u%error **2)
	       iv%metar(n)%v%sens = ktr%metar(n)%v / (iv%metar(n)%v%error **2)
	       iv%metar(n)%t%sens = ktr%metar(n)%t / (iv%metar(n)%t%error **2)
	       iv%metar(n)%p%sens = ktr%metar(n)%p / (iv%metar(n)%p%error **2)
	       iv%metar(n)%q%sens = ktr%metar(n)%q / (iv%metar(n)%q%error **2)	        	       
	       
	       iv%metar(n)%u%inv  = iv%metar(n)%u%inv * iv%metar(n)%u%sens
	       iv%metar(n)%v%inv  = iv%metar(n)%v%inv * iv%metar(n)%v%sens
	       iv%metar(n)%t%inv  = iv%metar(n)%t%inv * iv%metar(n)%t%sens
	       iv%metar(n)%p%inv  = iv%metar(n)%p%inv * iv%metar(n)%p%sens
	       iv%metar(n)%q%inv  = iv%metar(n)%q%inv * iv%metar(n)%q%sens	  
	       
	       ktd(metar,u) = ktd(metar,u) + iv%metar(n)%u%inv
               ktd(metar,v) = ktd(metar,v) + iv%metar(n)%v%inv 
	       ktd(metar,t) = ktd(metar,t) + iv%metar(n)%t%inv 
               ktd(metar,p) = ktd(metar,p) + iv%metar(n)%p%inv 
	       ktd(metar,q) = ktd(metar,q) + iv%metar(n)%q%inv	       	             	       	            	       
	    end do
	 end if   
	 
         if (iv%info(buoy)%nlocal > 0) then
            do n=1, iv%info(buoy)%nlocal
	       if (.not. iv%info(buoy)%proc_domain(1,n)) cycle
	       iv%buoy(n)%u%sens = ktr%buoy(n)%u / (iv%buoy(n)%u%error **2)
	       iv%buoy(n)%v%sens = ktr%buoy(n)%v / (iv%buoy(n)%v%error **2)
	       iv%buoy(n)%t%sens = ktr%buoy(n)%t / (iv%buoy(n)%t%error **2)
	       iv%buoy(n)%p%sens = ktr%buoy(n)%p / (iv%buoy(n)%p%error **2)
	       iv%buoy(n)%q%sens = ktr%buoy(n)%q / (iv%buoy(n)%q%error **2)	        	       
	       
	       iv%buoy(n)%u%inv  = iv%buoy(n)%u%inv * iv%buoy(n)%u%sens
	       iv%buoy(n)%v%inv  = iv%buoy(n)%v%inv * iv%buoy(n)%v%sens
	       iv%buoy(n)%t%inv  = iv%buoy(n)%t%inv * iv%buoy(n)%t%sens
	       iv%buoy(n)%p%inv  = iv%buoy(n)%p%inv * iv%buoy(n)%p%sens
	       iv%buoy(n)%q%inv  = iv%buoy(n)%q%inv * iv%buoy(n)%q%sens	  
	       
	       ktd(buoy,u) = ktd(buoy,u) + iv%buoy(n)%u%inv
               ktd(buoy,v) = ktd(buoy,v) + iv%buoy(n)%v%inv 
	       ktd(buoy,t) = ktd(buoy,t) + iv%buoy(n)%t%inv 
               ktd(buoy,p) = ktd(buoy,p) + iv%buoy(n)%p%inv 
	       ktd(buoy,q) = ktd(buoy,q) + iv%buoy(n)%q%inv	       	             	       	            	       
	    end do
	 end if   
	 
         if (iv%info(sound)%nlocal > 0) then
            do n=1, iv%info(sound)%nlocal
	       if (.not. iv%info(sound)%proc_domain(1,n)) cycle
	       iv%sonde_sfc(n)%u%sens = ktr%sonde_sfc(n)%u / (iv%sonde_sfc(n)%u%error **2)
	       iv%sonde_sfc(n)%v%sens = ktr%sonde_sfc(n)%v / (iv%sonde_sfc(n)%v%error **2)
	       iv%sonde_sfc(n)%t%sens = ktr%sonde_sfc(n)%t / (iv%sonde_sfc(n)%t%error **2)
	       iv%sonde_sfc(n)%p%sens = ktr%sonde_sfc(n)%p / (iv%sonde_sfc(n)%p%error **2)
	       iv%sonde_sfc(n)%q%sens = ktr%sonde_sfc(n)%q / (iv%sonde_sfc(n)%q%error **2)	        	       
	       
	       iv%sonde_sfc(n)%u%inv  = iv%sonde_sfc(n)%u%inv * iv%sonde_sfc(n)%u%sens
	       iv%sonde_sfc(n)%v%inv  = iv%sonde_sfc(n)%v%inv * iv%sonde_sfc(n)%v%sens
	       iv%sonde_sfc(n)%t%inv  = iv%sonde_sfc(n)%t%inv * iv%sonde_sfc(n)%t%sens
	       iv%sonde_sfc(n)%p%inv  = iv%sonde_sfc(n)%p%inv * iv%sonde_sfc(n)%p%sens
	       iv%sonde_sfc(n)%q%inv  = iv%sonde_sfc(n)%q%inv * iv%sonde_sfc(n)%q%sens	  
	       
	       ktd(sound,u) = ktd(sound,u) + iv%sonde_sfc(n)%u%inv
               ktd(sound,v) = ktd(sound,v) + iv%sonde_sfc(n)%v%inv 
	       ktd(sound,t) = ktd(sound,t) + iv%sonde_sfc(n)%t%inv 
               ktd(sound,p) = ktd(sound,p) + iv%sonde_sfc(n)%p%inv 
	       ktd(sound,q) = ktd(sound,q) + iv%sonde_sfc(n)%q%inv	       	             	       	            	       
	    end do
	 end if   
	 
         if (iv%info(qscat)%nlocal > 0) then
            do n=1, iv%info(qscat)%nlocal
	       if (.not. iv%info(qscat)%proc_domain(1,n)) cycle
	       iv%qscat(n)%u%sens = ktr%qscat(n)%u / (iv%qscat(n)%u%error **2)
	       iv%qscat(n)%v%sens = ktr%qscat(n)%v / (iv%qscat(n)%v%error **2)
	       
	       iv%qscat(n)%u%inv  = iv%qscat(n)%u%inv * iv%qscat(n)%u%sens
	       iv%qscat(n)%v%inv  = iv%qscat(n)%v%inv * iv%qscat(n)%v%sens
	       
	       ktd(qscat,u) = ktd(qscat,u) + iv%qscat(n)%u%inv
               ktd(qscat,v) = ktd(qscat,v) + iv%qscat(n)%v%inv 
	    end do
	 end if   
	 
         if (iv%info(sound)%nlocal > 0) then
            do n=1, iv%info(sound)%nlocal
               do k=1, iv%info(sound)%levels(n)
	          if (.not. iv%info(sound)%proc_domain(k,n)) cycle
	          iv%sound(n)%u(k)%sens = ktr%sound(n)%u(k) / (iv%sound(n)%u(k)%error **2)
	          iv%sound(n)%v(k)%sens = ktr%sound(n)%v(k) / (iv%sound(n)%v(k)%error **2)
	          iv%sound(n)%t(k)%sens = ktr%sound(n)%t(k) / (iv%sound(n)%t(k)%error **2)
	          iv%sound(n)%q(k)%sens = ktr%sound(n)%q(k) / (iv%sound(n)%q(k)%error **2)		  
	       
	          iv%sound(n)%u(k)%inv  = iv%sound(n)%u(k)%inv * iv%sound(n)%u(k)%sens
	          iv%sound(n)%v(k)%inv  = iv%sound(n)%v(k)%inv * iv%sound(n)%v(k)%sens
	          iv%sound(n)%t(k)%inv  = iv%sound(n)%t(k)%inv * iv%sound(n)%t(k)%sens
	          iv%sound(n)%q(k)%inv  = iv%sound(n)%q(k)%inv * iv%sound(n)%q(k)%sens	  
	       
	          ktd(sound,u) = ktd(sound,u) + iv%sound(n)%u(k)%inv
                  ktd(sound,v) = ktd(sound,v) + iv%sound(n)%v(k)%inv 
	          ktd(sound,t) = ktd(sound,t) + iv%sound(n)%t(k)%inv 
	          ktd(sound,q) = ktd(sound,q) + iv%sound(n)%q(k)%inv	       	             	       	            	       
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(mtgirs)%nlocal > 0) then
            do n=1, iv%info(mtgirs)%nlocal
               do k=1, iv%info(mtgirs)%levels(n)
	          if (.not. iv%info(mtgirs)%proc_domain(k,n)) cycle
	          iv%mtgirs(n)%u(k)%sens = ktr%mtgirs(n)%u(k) / (iv%mtgirs(n)%u(k)%error **2)
	          iv%mtgirs(n)%v(k)%sens = ktr%mtgirs(n)%v(k) / (iv%mtgirs(n)%v(k)%error **2)
	          iv%mtgirs(n)%t(k)%sens = ktr%mtgirs(n)%t(k) / (iv%mtgirs(n)%t(k)%error **2)
	          iv%mtgirs(n)%q(k)%sens = ktr%mtgirs(n)%q(k) / (iv%mtgirs(n)%q(k)%error **2)
	       
	          iv%mtgirs(n)%u(k)%inv  = iv%mtgirs(n)%u(k)%inv * iv%mtgirs(n)%u(k)%sens
	          iv%mtgirs(n)%v(k)%inv  = iv%mtgirs(n)%v(k)%inv * iv%mtgirs(n)%v(k)%sens
	          iv%mtgirs(n)%t(k)%inv  = iv%mtgirs(n)%t(k)%inv * iv%mtgirs(n)%t(k)%sens
	          iv%mtgirs(n)%q(k)%inv  = iv%mtgirs(n)%q(k)%inv * iv%mtgirs(n)%q(k)%sens	  
	       
	          ktd(mtgirs,u) = ktd(mtgirs,u) + iv%mtgirs(n)%u(k)%inv
                  ktd(mtgirs,v) = ktd(mtgirs,v) + iv%mtgirs(n)%v(k)%inv 
	          ktd(mtgirs,t) = ktd(mtgirs,t) + iv%mtgirs(n)%t(k)%inv 
	          ktd(mtgirs,q) = ktd(mtgirs,q) + iv%mtgirs(n)%q(k)%inv	       	             	       	            	       
	       end do
	    end do
	 end if   
	       
         if (iv%info(bogus)%nlocal > 0) then
            do n=1, iv%info(bogus)%nlocal
               do k=1, iv%info(bogus)%levels(n)
	          if (.not. iv%info(bogus)%proc_domain(k,n)) cycle
	          iv%bogus(n)%u(k)%sens = ktr%bogus(n)%u(k) / (iv%bogus(n)%u(k)%error **2)
	          iv%bogus(n)%v(k)%sens = ktr%bogus(n)%v(k) / (iv%bogus(n)%v(k)%error **2)
	          iv%bogus(n)%t(k)%sens = ktr%bogus(n)%t(k) / (iv%bogus(n)%t(k)%error **2)
	          iv%bogus(n)%q(k)%sens = ktr%bogus(n)%q(k) / (iv%bogus(n)%q(k)%error **2)
	       
	          iv%bogus(n)%u(k)%inv  = iv%bogus(n)%u(k)%inv * iv%bogus(n)%u(k)%sens
	          iv%bogus(n)%v(k)%inv  = iv%bogus(n)%v(k)%inv * iv%bogus(n)%v(k)%sens
	          iv%bogus(n)%t(k)%inv  = iv%bogus(n)%t(k)%inv * iv%bogus(n)%t(k)%sens
	          iv%bogus(n)%q(k)%inv  = iv%bogus(n)%q(k)%inv * iv%bogus(n)%q(k)%sens	  
	       
	          ktd(bogus,u) = ktd(bogus,u) + iv%bogus(n)%u(k)%inv
                  ktd(bogus,v) = ktd(bogus,v) + iv%bogus(n)%v(k)%inv 
	          ktd(bogus,t) = ktd(bogus,t) + iv%bogus(n)%t(k)%inv 
	          ktd(bogus,q) = ktd(bogus,q) + iv%bogus(n)%q(k)%inv	       	             	       	            	       
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(pilot)%nlocal > 0) then
            do n=1, iv%info(pilot)%nlocal
               do k=1, iv%info(pilot)%levels(n)
	          if (.not. iv%info(pilot)%proc_domain(k,n)) cycle
	          iv%pilot(n)%u(k)%sens = ktr%pilot(n)%u(k) / (iv%pilot(n)%u(k)%error **2)
	          iv%pilot(n)%v(k)%sens = ktr%pilot(n)%v(k) / (iv%pilot(n)%v(k)%error **2)
	       
	          iv%pilot(n)%u(k)%inv  = iv%pilot(n)%u(k)%inv * iv%pilot(n)%u(k)%sens
	          iv%pilot(n)%v(k)%inv  = iv%pilot(n)%v(k)%inv * iv%pilot(n)%v(k)%sens
	       
	          ktd(pilot,u) = ktd(pilot,u) + iv%pilot(n)%u(k)%inv
                  ktd(pilot,v) = ktd(pilot,v) + iv%pilot(n)%v(k)%inv 
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(airep)%nlocal > 0) then
            do n=1, iv%info(airep)%nlocal
               do k=1, iv%info(airep)%levels(n)
	          if (.not. iv%info(airep)%proc_domain(k,n)) cycle
	          iv%airep(n)%u(k)%sens = ktr%airep(n)%u(k) / (iv%airep(n)%u(k)%error **2)
	          iv%airep(n)%v(k)%sens = ktr%airep(n)%v(k) / (iv%airep(n)%v(k)%error **2)
	          iv%airep(n)%t(k)%sens = ktr%airep(n)%t(k) / (iv%airep(n)%t(k)%error **2)
	          iv%airep(n)%q(k)%sens = ktr%airep(n)%q(k) / (iv%airep(n)%q(k)%error **2)
	       
	          iv%airep(n)%u(k)%inv  = iv%airep(n)%u(k)%inv * iv%airep(n)%u(k)%sens
	          iv%airep(n)%v(k)%inv  = iv%airep(n)%v(k)%inv * iv%airep(n)%v(k)%sens
	          iv%airep(n)%t(k)%inv  = iv%airep(n)%t(k)%inv * iv%airep(n)%t(k)%sens
	          iv%airep(n)%q(k)%inv  = iv%airep(n)%q(k)%inv * iv%airep(n)%q(k)%sens
	       
	          ktd(airep,u) = ktd(airep,u) + iv%airep(n)%u(k)%inv
                  ktd(airep,v) = ktd(airep,v) + iv%airep(n)%v(k)%inv 
	          ktd(airep,t) = ktd(airep,t) + iv%airep(n)%t(k)%inv 
	          ktd(airep,q) = ktd(airep,q) + iv%airep(n)%q(k)%inv 
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(geoamv)%nlocal > 0) then
            do n=1, iv%info(geoamv)%nlocal
               do k=1, iv%info(geoamv)%levels(n)
	          if (.not. iv%info(geoamv)%proc_domain(k,n)) cycle
	          iv%geoamv(n)%u(k)%sens = ktr%geoamv(n)%u(k) / (iv%geoamv(n)%u(k)%error **2)
	          iv%geoamv(n)%v(k)%sens = ktr%geoamv(n)%v(k) / (iv%geoamv(n)%v(k)%error **2)
	       
	          iv%geoamv(n)%u(k)%inv  = iv%geoamv(n)%u(k)%inv * iv%geoamv(n)%u(k)%sens
	          iv%geoamv(n)%v(k)%inv  = iv%geoamv(n)%v(k)%inv * iv%geoamv(n)%v(k)%sens
	       
	          ktd(geoamv,u) = ktd(geoamv,u) + iv%geoamv(n)%u(k)%inv
                  ktd(geoamv,v) = ktd(geoamv,v) + iv%geoamv(n)%v(k)%inv 
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(polaramv)%nlocal > 0) then
            do n=1, iv%info(polaramv)%nlocal
               do k=1, iv%info(polaramv)%levels(n)
	          if (.not. iv%info(polaramv)%proc_domain(k,n)) cycle
	          iv%polaramv(n)%u(k)%sens = ktr%polaramv(n)%u(k) / (iv%polaramv(n)%u(k)%error **2)
	          iv%polaramv(n)%v(k)%sens = ktr%polaramv(n)%v(k) / (iv%polaramv(n)%v(k)%error **2)
	       
	          iv%polaramv(n)%u(k)%inv  = iv%polaramv(n)%u(k)%inv * iv%polaramv(n)%u(k)%sens
	          iv%polaramv(n)%v(k)%inv  = iv%polaramv(n)%v(k)%inv * iv%polaramv(n)%v(k)%sens
	       
	          ktd(polaramv,u) = ktd(polaramv,u) + iv%polaramv(n)%u(k)%inv
                  ktd(polaramv,v) = ktd(polaramv,v) + iv%polaramv(n)%v(k)%inv 
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(profiler)%nlocal > 0) then
            do n=1, iv%info(profiler)%nlocal
               if (.not. iv%info(profiler)%proc_domain(1,n)) cycle
	       do k=1, iv%info(profiler)%levels(n)
	          iv%profiler(n)%u(k)%sens = ktr%profiler(n)%u(k) / (iv%profiler(n)%u(k)%error **2)
	          iv%profiler(n)%v(k)%sens = ktr%profiler(n)%v(k) / (iv%profiler(n)%v(k)%error **2)
	       
	          iv%profiler(n)%u(k)%inv  = iv%profiler(n)%u(k)%inv * iv%profiler(n)%u(k)%sens
	          iv%profiler(n)%v(k)%inv  = iv%profiler(n)%v(k)%inv * iv%profiler(n)%v(k)%sens
	       
	          ktd(profiler,u) = ktd(profiler,u) + iv%profiler(n)%u(k)%inv
                  ktd(profiler,v) = ktd(profiler,v) + iv%profiler(n)%v(k)%inv 
	       end do
	    end do
	 end if   
	       	       
         if (iv%info(satem)%nlocal > 0) then
            do n=1, iv%info(satem)%nlocal
               if (.not. iv%info(satem)%proc_domain(1,n)) cycle
	       do k=1, iv%info(satem)%levels(n)	    
	          iv%satem(n)%thickness(k)%sens = ktr%satem(n)%thickness(k) / (iv%satem(n)%thickness(k)%error **2)
	       
	          iv%satem(n)%thickness(k)%inv  = iv%satem(n)%thickness(k)%inv * iv%satem(n)%thickness(k)%sens
	       
	          ktd(satem,mx) = ktd(satem,mx) + iv%satem(n)%thickness(k)%inv
               end do
	    end do
	 end if   
	       
         if (iv%info(gpspw)%nlocal > 0) then
            do n=1, iv%info(gpspw)%nlocal
	       if (.not. iv%info(gpspw)%proc_domain(1,n)) cycle
	       iv%gpspw(n)%tpw%sens = ktr%gpspw(n)%tpw / (iv%gpspw(n)%tpw%error **2)
	       
	       iv%gpspw(n)%tpw%inv  = iv%gpspw(n)%tpw%inv * iv%gpspw(n)%tpw%sens
	       
	       ktd(gpspw,q) = ktd(gpspw,q) + iv%gpspw(n)%tpw%inv
	    end do
	 end if   
	       
          if (iv%info(gpsref)%nlocal > 0) then
            do n=1, iv%info(gpsref)%nlocal
               if (.not. iv%info(gpsref)%proc_domain(1,n)) cycle
               do k=1, iv%info(gpsref)%levels(n)
	          if (iv%gpsref(n)%ref(k)%qc < obs_qc_pointer) cycle	    
	          iv%gpsref(n)%ref(k)%sens = ktr%gpsref(n)%ref(k) / (iv%gpsref(n)%ref(k)%error **2)
 	       
	          iv%gpsref(n)%ref(k)%inv  = iv%gpsref(n)%ref(k)%inv * iv%gpsref(n)%ref(k)%sens
	       
	          ktd(gpsref,rv) = ktd(gpsref,rv) + iv%gpsref(n)%ref(k)%inv
              end do
	    end do
	 end if   
	       
         if (iv%info(ssmi_rv)%nlocal > 0) then
            do n=1, iv%info(ssmi_rv)%nlocal
	       if (.not. iv%info(ssmi_rv)%proc_domain(1,n)) cycle
	       iv%ssmi_rv(n)%Speed%sens = ktr%ssmi_rv(n)%Speed / (iv%ssmi_rv(n)%Speed%error **2)
 	       iv%ssmi_rv(n)%tpw%sens = ktr%ssmi_rv(n)%tpw / (iv%ssmi_rv(n)%tpw%error **2)
	       
	       iv%ssmi_rv(n)%Speed%inv  = iv%ssmi_rv(n)%Speed%inv * iv%ssmi_rv(n)%Speed%sens
 	       iv%ssmi_rv(n)%tpw%inv  = iv%ssmi_rv(n)%tpw%inv * iv%ssmi_rv(n)%tpw%sens
	       
 	       ktd(ssmi_rv,mx) = ktd(ssmi_rv,mx) + iv%ssmi_rv(n)%Speed%inv + iv%ssmi_rv(n)%tpw%inv
	    end do
	 end if   

         if (iv%info(tamdar)%nlocal > 0) then
            do n=1, iv%info(tamdar)%nlocal
               do k=1, iv%info(tamdar)%levels(n)
                  if (.not. iv%info(tamdar)%proc_domain(k,n)) cycle
                  iv%tamdar(n)%u(k)%sens = ktr%tamdar(n)%u(k) / (iv%tamdar(n)%u(k)%error **2)
                  iv%tamdar(n)%v(k)%sens = ktr%tamdar(n)%v(k) / (iv%tamdar(n)%v(k)%error **2)
                  iv%tamdar(n)%t(k)%sens = ktr%tamdar(n)%t(k) / (iv%tamdar(n)%t(k)%error **2)
                  iv%tamdar(n)%q(k)%sens = ktr%tamdar(n)%q(k) / (iv%tamdar(n)%q(k)%error **2)

                  iv%tamdar(n)%u(k)%inv  = iv%tamdar(n)%u(k)%inv * iv%tamdar(n)%u(k)%sens
                  iv%tamdar(n)%v(k)%inv  = iv%tamdar(n)%v(k)%inv * iv%tamdar(n)%v(k)%sens
                  iv%tamdar(n)%t(k)%inv  = iv%tamdar(n)%t(k)%inv * iv%tamdar(n)%t(k)%sens
                  iv%tamdar(n)%q(k)%inv  = iv%tamdar(n)%q(k)%inv * iv%tamdar(n)%q(k)%sens

                  ktd(tamdar,u) = ktd(tamdar,u) + iv%tamdar(n)%u(k)%inv
                  ktd(tamdar,v) = ktd(tamdar,v) + iv%tamdar(n)%v(k)%inv
                  ktd(tamdar,t) = ktd(tamdar,t) + iv%tamdar(n)%t(k)%inv
                  ktd(tamdar,q) = ktd(tamdar,q) + iv%tamdar(n)%q(k)%inv
               end do
            end do
         end if

         if (iv%info(tamdar_sfc)%nlocal > 0) then
            do n=1, iv%info(tamdar_sfc)%nlocal
               if (.not. iv%info(tamdar_sfc)%proc_domain(1,n)) cycle
               iv%tamdar_sfc(n)%u%sens = ktr%tamdar_sfc(n)%u / (iv%tamdar_sfc(n)%u%error **2)
               iv%tamdar_sfc(n)%v%sens = ktr%tamdar_sfc(n)%v / (iv%tamdar_sfc(n)%v%error **2)
               iv%tamdar_sfc(n)%t%sens = ktr%tamdar_sfc(n)%t / (iv%tamdar_sfc(n)%t%error **2)
               iv%tamdar_sfc(n)%q%sens = ktr%tamdar_sfc(n)%q / (iv%tamdar_sfc(n)%q%error **2)

               iv%tamdar_sfc(n)%u%inv  = iv%tamdar_sfc(n)%u%inv * iv%tamdar_sfc(n)%u%sens
               iv%tamdar_sfc(n)%v%inv  = iv%tamdar_sfc(n)%v%inv * iv%tamdar_sfc(n)%v%sens
               iv%tamdar_sfc(n)%t%inv  = iv%tamdar_sfc(n)%t%inv * iv%tamdar_sfc(n)%t%sens
               iv%tamdar_sfc(n)%q%inv  = iv%tamdar_sfc(n)%q%inv * iv%tamdar_sfc(n)%q%sens

               ktd(tamdar,u) = ktd(tamdar,u) + iv%tamdar_sfc(n)%u%inv
               ktd(tamdar,v) = ktd(tamdar,v) + iv%tamdar_sfc(n)%v%inv
               ktd(tamdar,t) = ktd(tamdar,t) + iv%tamdar_sfc(n)%t%inv
               ktd(tamdar,q) = ktd(tamdar,q) + iv%tamdar_sfc(n)%q%inv
            end do
         end if

	       
         if (iv%num_inst > 0) then
            do inst = 1, iv%num_inst                                       ! loop for sensor
               if (iv%instid(inst)%num_rad < 1) cycle
	       do n= 1, iv%instid(inst)%num_rad                            ! loop for pixel
                  if (.not. iv%instid(inst)%info%proc_domain(1,n)) cycle
                  do k=1, iv%instid(inst)%nchan	                           ! loop for channel
		     if ( iv%instid(inst)%tb_qc(k,n) < obs_qc_pointer ) cycle
		     iv%instid(inst)%tb_sens(k,n) = ktr%instid(inst)%tb(k,n) / (iv%instid(inst)%tb_error(k,n) **2)
		     
		     iv%instid(inst)%tb_inv(k,n) = iv%instid(inst)%tb_inv(k,n) * iv%instid(inst)%tb_sens(k,n)
		     
	             ktdrad(inst) = ktdrad(inst) + iv%instid(inst)%tb_inv(k,n)
		     
                    if ( use_varbc ) then
                   ! Impact of Bias Predictors 
		     npred = iv%instid(inst)%varbc(k)%npred
	             do i = 1, npred
		        ipred = iv%instid(inst)%varbc(k)%ipred(i)
	                ktbrad(inst,ipred) = ktbrad(inst,ipred) - &
			                     iv%instid(inst)%varbc(k)%param(i) * &
					     iv%instid(inst)%varbc_info%pred(ipred,n) * &
					     iv%instid(inst)%tb_inv(k,n)
		     end do
                    end if
		     
                  end do                                                   ! loop for channel
	       end do                                                      ! loop for pixel
	    end do                                                         ! loop for sensor
	 end if   
	       
	! Sum across processors
	 do i = 1, nbvar 
            call wrf_dm_sum_reals(ktd(:,i), ktd_global(:,i))
	 end do

         if ( iv%num_inst > 0 ) &
	    call wrf_dm_sum_reals(ktdrad, ktdrad_global)

         if ( iv%num_inst > 0 .and. use_varbc ) then
	    do i = 1, npredmax 
               call wrf_dm_sum_reals(ktbrad(:,i), ktbrad_global(:,i))
	    end do
         endif
         
         write(unit=message(1),fmt='(A)') 'Impact of Conventional Observations for each variable type: '
         do i = 1, nbvar
            write(unit=message(1+i),fmt='(3x,a,2x,e15.5)') var_names(i), SUM(ktd_global(:,i))
         end do
         call da_message(message(1:nbvar+1))
         imsg = 1
         write(unit=message(imsg),fmt='(A)') 'Impact of Conventional Observations for each observation type: '
         do i = 1, num_ob_indexes
            if ( (i == ssmi_tb) .or. (i == ssmt1) .or. (i == ssmt2) .or. &
                 (i == radar ) .or. (i == radiance) .or. (i == airsr) .or. &
                 (i == sonde_sfc) .or. (i == tamdar_sfc) .or. (i == rain) ) cycle
            imsg = imsg + 1
            write(unit=message(imsg),fmt='(3x,a,e15.5)') obs_names(i), SUM(ktd_global(i,:))
         end do
         call da_message(message(1:imsg))

         if ( iv%num_inst > 0 ) then
            imsg = 1
            write(unit=message(imsg),fmt='(A)') 'Impact of Satellite Radiances for each instrument: '
            do i = 1, iv%num_inst
               imsg = imsg + 1
               write(unit=message(imsg),fmt='(3x,a,e15.5)') iv%instid(i)%rttovid_string, ktdrad_global(i)
            end do
            call da_message(message(1:imsg))
            if ( use_varbc ) then
               imsg = 1
               write(unit=message(imsg),fmt='(A)') 'Impact of Satellite Bias Correction for each predictor: '
               do i = 1, npredmax
                  imsg = imsg + 1
                  write(unit=message(imsg),fmt='(3x,e15.5)') SUM(ktbrad_global(:,i))
               end do
               call da_message(message(1:imsg))
               imsg = 1
               write(unit=message(imsg),fmt='(A)') 'Impact of Satellite Bias Correction for each instrument:'
               do i = 1, iv%num_inst
                  imsg = imsg + 1
                  write(unit=message(imsg),fmt='(3x,a,e15.5)') iv%instid(i)%rttovid_string, SUM(ktbrad_global(i,:))
               end do
               call da_message(message(1:imsg))
            end if
         endif

         ! output the impact for ploting
         open(unit=iunit,file='obs_impact',form='formatted')

         write(iunit,'(A)') 'Impact of Conventional Observations: '
         write(iunit,'(7(A5,2x))') 'U', 'V', 'T', 'P', 'Q', 'GPS', 'SATEM'
         write(iunit,'(7(e15.5,2x))') SUM(ktd_global,dim=1)
         write(iunit,'(26(A10,2x))') 'Sound', 'Synop', 'Pilot', 'Satem', 'GeoAMV', 'PolarAMV', 'AIREP',&
                              'GPSZTD', 'GPSRF', 'METAR', 'Ships', 'SSMI_RV', 'SSMI_TB', &
                              'SSMT1', 'SSMT2', 'QSCAT', 'Profiler', 'Buoy', &
                              'Bogus', 'Pseudo', 'Radar', 'Radiance', 'AIRSR', 'Sonde_sfc', 'MTGIRS', 'TAMDAR'
         write(iunit,'(26(e15.5,2x))') SUM(ktd_global,dim=2)
         if ( iv%num_inst > 0 ) then
         write(iunit,'(A)') 'Impact of Satellite Radiances for each instrument: '
         write(iunit,'(100e15.5)') ktdrad_global
         write(iunit,'(A)') 'Impact of Satellite Bias Correction for each predictor: '
         write(iunit,'(100e15.5)') SUM(ktbrad_global,dim=1)
         write(iunit,'(A)') 'Impact of Satellite Bias Correction for each instrument:'
         write(iunit,'(100e15.5)') SUM(ktbrad_global,dim=2)
         endif

         close(iunit)
         call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_obs_sensitivity")

end subroutine da_obs_sensitivity

subroutine da_set_obs_missing (iv, n)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: n    ! obs index

   integer :: i, k

   if (trace_use) call da_trace_entry("da_set_obs_missing")

   do i = 1, iv%info(n)%nlocal
      if ( iv%info(n)%thinned(1,i) ) then
         iv%info(n)%slp(i) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         iv%info(n)%pw(i)  = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
      end if
   end do

   select case (n)
   case (sound)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%sound(i)%h(:) = missing_r
            iv%sound(i)%p(:) = missing_r
            iv%sound(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sound(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sound(i)%t(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sound(i)%q(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case (synop)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%synop(i)%h = missing_r
            iv%synop(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%synop(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%synop(i)%t = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%synop(i)%p = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%synop(i)%q = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case (pilot)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%pilot(i)%h(:) = missing_r
            iv%pilot(i)%p(:) = missing_r
            iv%pilot(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%pilot(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(satem)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%satem(i)%ref_p        = missing_r
            iv%satem(i)%p(:)         = missing_r
            iv%satem(i)%thickness(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(geoamv)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%geoamv(i)%p(:) = missing_r
            iv%geoamv(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%geoamv(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(polaramv)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%polaramv(i)%p(:) = missing_r
            iv%polaramv(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%polaramv(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(airep)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%airep(i)%h(:) = missing_r
            iv%airep(i)%p(:) = missing_r
            iv%airep(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%airep(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%airep(i)%t(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%airep(i)%q(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(gpspw)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%gpspw(i)%tpw = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(gpsref)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%gpsref(i)%h(:)   = missing_r
            iv%gpsref(i)%ref(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%gpsref(i)%p(:)   = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%gpsref(i)%t(:)   = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%gpsref(i)%q(:)   = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(metar)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%metar(i)%h = missing_r
            iv%metar(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%metar(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%metar(i)%t = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%metar(i)%p = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%metar(i)%q = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(ships)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%ships(i)%h = missing_r
            iv%ships(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%ships(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%ships(i)%t = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%ships(i)%p = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%ships(i)%q = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(ssmi_rv)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%ssmi_rv(i)%speed = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%ssmi_rv(i)%tpw   = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(qscat)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%qscat(i)%h = missing_r
            iv%qscat(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%qscat(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(profiler)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%profiler(i)%h(:) = missing_r
            iv%profiler(i)%p(:) = missing_r
            iv%profiler(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%profiler(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(buoy)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%buoy(i)%h = missing_r
            iv%buoy(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%buoy(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%buoy(i)%t = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%buoy(i)%p = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%buoy(i)%q = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(bogus)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%bogus(i)%h(:) = missing_r
            iv%bogus(i)%p(:) = missing_r
            iv%bogus(i)%u(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%bogus(i)%v(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%bogus(i)%t(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%bogus(i)%q(:) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%bogus(i)%slp  = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(sonde_sfc)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%sonde_sfc(i)%h = missing_r
            iv%sonde_sfc(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sonde_sfc(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sonde_sfc(i)%t = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sonde_sfc(i)%p = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%sonde_sfc(i)%q = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case(rain)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%rain(i)%rain = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case default
      if (trace_use) call da_trace_exit("da_set_obs_missing")
      return
   end select

   if (trace_use) call da_trace_exit("da_set_obs_missing")
 
end subroutine da_set_obs_missing


subroutine da_set_3d_obs_missing (iv, n)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(iv_type), intent(inout) :: iv
   integer,       intent(in)    :: n    ! obs index

   integer :: i, k
   real    :: xmiss

   if (trace_use) call da_trace_entry("da_set_obs_missing")

   xmiss = -888.0

   select case (n)
   case (radar)
      do i = 1, iv%info(n)%nlocal
        do k=1,iv%info(n)%levels(i)
         if ( iv%info(n)%thinned(k,i) ) then
            iv%radar(i)%height(k) = missing_r
            iv%radar(i)%height_qc(k) = missing_data
            iv%radar(i)%rv(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%radar(i)%rf(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
        end do
      end do
   case (airep)
      do i = 1, iv%info(n)%nlocal
        do k=1,iv%info(n)%levels(i)
         if ( iv%info(n)%thinned(k,i) ) then
            iv%airep(i)%h(k) = missing_r
            iv%airep(i)%p(k) = missing_r
            iv%airep(i)%u(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%airep(i)%v(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%airep(i)%t(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%airep(i)%q(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
        end do
      end do
   case (tamdar)
      do i = 1, iv%info(n)%nlocal
        do k=1,iv%info(n)%levels(i)
         if ( iv%info(n)%thinned(k,i) ) then
            iv%tamdar(i)%h(k) = missing_r
            iv%tamdar(i)%p(k) = missing_r
            iv%tamdar(i)%u(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar(i)%v(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar(i)%t(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar(i)%q(k) = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
        end do
      end do
   case (tamdar_sfc)
      do i = 1, iv%info(n)%nlocal
         if ( iv%info(n)%thinned(1,i) ) then
            iv%tamdar_sfc(i)%h = missing_r
            iv%tamdar_sfc(i)%u = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar_sfc(i)%v = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar_sfc(i)%t = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar_sfc(i)%p = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
            iv%tamdar_sfc(i)%q = field_type(missing_r, missing_data, xmiss, missing_r, missing_r)
         end if
      end do
   case default
      write(unit=message(1),fmt='(A,I4)') 'Wrong obs_index= ',n
      call da_error("da_set_3d_obs_missing.inc",70,message(1:1))
   end select

   if (trace_use) call da_trace_exit("da_set_3d_obs_missing")
 
end subroutine da_set_3d_obs_missing



end module da_obs
