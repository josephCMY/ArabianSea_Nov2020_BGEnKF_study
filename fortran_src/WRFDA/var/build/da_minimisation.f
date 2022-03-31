












module da_minimisation

   
   
   

   use module_configure, only : grid_config_rec_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use module_dm, only : local_communicator, mytask, ntasks, ntasks_x, &
      ntasks_y, data_order_xy, data_order_xyz
   use module_comm_dm, only : halo_wpec_sub, halo_wpec_adj_sub

   use module_domain, only : domain, ep_type, vp_type, x_type, domain_clockprint, &
                             domain_clockadvance, domain_clock_get, domain_clock_set
   use module_state_description, only : dyn_em,dyn_em_tl,dyn_em_ad,p_g_qv, &
       p_g_qc, p_g_qr, num_moist, PARAM_FIRST_SCALAR





   use da_airep, only : da_calculate_grady_airep, da_ao_stats_airep, &
      da_oi_stats_airep, da_get_innov_vector_airep, da_residual_airep, &
      da_jo_and_grady_airep
   use da_airsr , only : da_calculate_grady_airsr, da_ao_stats_airsr, &
      da_oi_stats_airsr, da_get_innov_vector_airsr, da_residual_airsr, &
      da_jo_and_grady_airsr
   use da_bogus, only : da_calculate_grady_bogus, da_ao_stats_bogus, &
      da_oi_stats_bogus, da_get_innov_vector_bogus, da_residual_bogus, &
      da_jo_and_grady_bogus
   use da_buoy , only : da_calculate_grady_buoy, da_ao_stats_buoy, &
      da_oi_stats_buoy,da_get_innov_vector_buoy, da_residual_buoy, &
      da_jo_and_grady_buoy
   use da_control, only : trace_use, var4d_bin, trajectory_io, analysis_date, &
      var4d, rootproc,jcdfi_use,jcdfi_diag,ierr,comm,num_fgat_time, &
      var4d_lbc, stdout, eps, stats_unit, test_dm_exact, global, multi_inc, &
      calculate_cg_cost_fn,anal_type_randomcv,cv_size_domain,je_factor, &
      jb_factor,ntmax,omb_add_noise,write_iv_rad_ascii,use_obs_errfac, &
      rtm_option,rtm_option_rttov, rtm_option_crtm, anal_type_verify, &
      write_filtered_rad,omb_set_rand,use_rad,var_scaling2,var_scaling1, &
      var_scaling4,var_scaling5,var_scaling3, jo_unit, test_gradient, &
      print_detail_grad,omb_set_rand,grad_unit,cost_unit, num_pseudo, cv_options, &
      cv_size_domain_je,cv_size_domain_jb, cv_size_domain_jp, cv_size_domain_js, cv_size_domain_jl, &
      sound, mtgirs, sonde_sfc, synop, profiler, gpsref, gpspw, polaramv, geoamv, ships, metar, &
      satem, radar, ssmi_rv, ssmi_tb, ssmt1, ssmt2, airsr, pilot, airep,tamdar, tamdar_sfc, rain, &
      bogus, buoy, qscat,pseudo, radiance, monitor_on, max_ext_its, use_rttov_kmatrix,&
      use_crtm_kmatrix,precondition_cg, precondition_factor, use_varbc, varbc_factor, &
      num_procs, myproc, use_gpspwobs, use_rainobs, use_gpsztdobs, &
      use_radar_rf, use_radar_rhv,use_radar_rqv,pseudo_var, num_pseudo, &
      num_ob_indexes, num_ob_vars, npres_print, pptop, ppbot, qcstat_conv_unit, gas_constant, &
      orthonorm_gradient, its, ite, jts, jte, kts, kte, ids, ide, jds, jde, kds, kde, cp, &
      use_satcv, sensitivity_option, print_detail_outerloop, adj_sens, filename_len, &
      ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe, fgat_rain_flags, var4d_bin_rain, freeze_varbc, &
      use_wpec, wpec_factor
   use da_define_structures, only : iv_type, y_type,  j_type, be_type, &
      xbx_type, jo_type, da_allocate_y,da_zero_x,da_zero_y,da_deallocate_y, &
      da_zero_vp_type, qhat_type
   use da_dynamics, only : da_wpec_constraint_lin,da_wpec_constraint_adj
   use da_obs, only : da_transform_xtoy_adj,da_transform_xtoy, &
      da_add_noise_to_ob,da_random_omb_all, da_obs_sensitivity
   use da_geoamv, only : da_calculate_grady_geoamv, da_ao_stats_geoamv, &
      da_oi_stats_geoamv, da_get_innov_vector_geoamv,da_residual_geoamv, &
      da_jo_and_grady_geoamv
   use da_gpspw, only : da_calculate_grady_gpspw, da_ao_stats_gpspw, &
      da_oi_stats_gpspw, da_get_innov_vector_gpspw, da_residual_gpspw, &
      da_jo_and_grady_gpspw, da_get_innov_vector_gpsztd
   use da_gpsref, only : da_calculate_grady_gpsref, da_ao_stats_gpsref, &
      da_oi_stats_gpsref, da_get_innov_vector_gpsref, da_residual_gpsref, &
      da_jo_and_grady_gpsref
   use da_obs_io, only : da_final_write_y, da_write_y, da_final_write_obs, &
      da_write_obs,da_write_obs_etkf,da_write_noise_to_ob, da_use_obs_errfac, &
      da_write_iv_for_multi_inc, da_read_iv_for_multi_inc
   use da_metar, only : da_calculate_grady_metar, da_ao_stats_metar, &
      da_oi_stats_metar, da_get_innov_vector_metar, da_residual_metar, &
      da_jo_and_grady_metar
   use da_pilot, only : da_calculate_grady_pilot, da_ao_stats_pilot, &
      da_oi_stats_pilot, da_get_innov_vector_pilot, da_residual_pilot, &
      da_jo_and_grady_pilot
   use da_par_util, only : da_system,da_cv_to_global
   use da_par_util1, only : da_proc_sum_real,da_proc_sum_ints
   use da_polaramv, only : da_calculate_grady_polaramv, da_ao_stats_polaramv, &
      da_oi_stats_polaramv, da_get_innov_vector_polaramv, da_residual_polaramv, &
      da_jo_and_grady_polaramv
   use da_profiler, only : da_calculate_grady_profiler, da_ao_stats_profiler, &
      da_oi_stats_profiler,da_get_innov_vector_profiler, da_residual_profiler, &
      da_jo_and_grady_profiler
   use da_pseudo, only : da_calculate_grady_pseudo, da_ao_stats_pseudo, &
      da_oi_stats_pseudo, da_get_innov_vector_pseudo, da_residual_pseudo, &
      da_jo_and_grady_pseudo
   use da_qscat, only : da_calculate_grady_qscat, da_ao_stats_qscat, &
      da_oi_stats_qscat, da_get_innov_vector_qscat, da_residual_qscat, &
      da_jo_and_grady_qscat
   use da_mtgirs, only : da_calculate_grady_mtgirs, &
      da_ao_stats_mtgirs, da_oi_stats_mtgirs,da_oi_stats_mtgirs, &
      da_get_innov_vector_mtgirs, &
      da_jo_and_grady_mtgirs, da_residual_mtgirs
   use da_tamdar, only : da_calculate_grady_tamdar, &
      da_ao_stats_tamdar, da_oi_stats_tamdar,da_oi_stats_tamdar, &
      da_get_innov_vector_tamdar, &
      da_jo_and_grady_tamdar, da_residual_tamdar, &
      da_calculate_grady_tamdar_sfc, &
      da_ao_stats_tamdar_sfc, da_oi_stats_tamdar_sfc,da_oi_stats_tamdar_sfc, &
      da_get_innov_vector_tamdar_sfc, &
      da_jo_and_grady_tamdar_sfc, da_residual_tamdar_sfc

   use da_radiance, only : da_calculate_grady_rad, da_write_filtered_rad, &
      da_get_innov_vector_radiance, satinfo
   use da_radiance1, only : da_ao_stats_rad,da_oi_stats_rad, &
      da_write_iv_rad_ascii,da_residual_rad,da_jo_and_grady_rad
   use da_radar, only :  da_calculate_grady_radar, da_ao_stats_radar, &
      da_oi_stats_radar, da_get_innov_vector_radar, da_residual_radar, &
      da_jo_and_grady_radar

   use da_rain, only :  da_calculate_grady_rain, da_ao_stats_rain, &
      da_oi_stats_rain, da_get_innov_vector_rain, da_residual_rain, &
      da_jo_and_grady_rain, da_get_hr_rain, da_transform_xtoy_rain, &
      da_transform_xtoy_rain_adj

   use da_reporting, only : da_message, da_warning, da_error
   use da_satem, only : da_calculate_grady_satem, da_ao_stats_satem, &
      da_oi_stats_satem, da_get_innov_vector_satem, da_residual_satem, &
      da_jo_and_grady_satem
   use da_ships, only : da_calculate_grady_ships, da_ao_stats_ships, &
      da_oi_stats_ships, da_get_innov_vector_ships, da_residual_ships, &
      da_jo_and_grady_ships
   use da_sound, only : da_calculate_grady_sound,da_calculate_grady_sonde_sfc, &
      da_ao_stats_sound, da_oi_stats_sound,da_oi_stats_sound, &
      da_oi_stats_sonde_sfc,da_ao_stats_sonde_sfc,da_get_innov_vector_sound, &
      da_get_innov_vector_sonde_sfc,da_jo_and_grady_sound, da_residual_sound, &
      da_jo_and_grady_sound,da_jo_and_grady_sonde_sfc,da_residual_sonde_sfc
   use da_ssmi, only : da_calculate_grady_ssmi_tb,da_calculate_grady_ssmi_rv,da_calculate_grady_ssmt1, &
      da_calculate_grady_ssmt2, da_ao_stats_ssmi_tb ,da_ao_stats_ssmt2, &
      da_ao_stats_ssmt2, da_oi_stats_ssmt1, da_oi_stats_ssmt2, &
      da_oi_stats_ssmi_tb,da_oi_stats_ssmi_rv,da_ao_stats_ssmt1,da_get_innov_vector_ssmi_tb, &
      da_get_innov_vector_ssmi_rv, da_residual_ssmi_rv, da_residual_ssmi_tb, &
      da_get_innov_vector_ssmt1,da_get_innov_vector_ssmt2, &
      da_jo_and_grady_ssmt1, da_jo_and_grady_ssmt2,da_jo_and_grady_ssmi_tb, &
      da_jo_and_grady_ssmi_rv, &
      da_residual_ssmt1,da_residual_ssmt2, da_ao_stats_ssmi_rv
   use da_synop, only : da_calculate_grady_synop, da_ao_stats_synop, &
      da_oi_stats_synop, da_get_innov_vector_synop, da_residual_synop, &
      da_jo_and_grady_synop
   use da_statistics, only : da_analysis_stats, da_print_qcstat
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tracing, only : da_trace_entry, da_trace_exit,da_trace
   use da_transfer_model, only : da_transfer_wrftltoxa,da_transfer_xatowrftl, &
      da_transfer_xatowrftl_adj,da_transfer_wrftltoxa_adj
   use da_varbc, only : da_varbc_tl,da_varbc_adj,da_varbc_precond,da_varbc_coldstart
   use da_vtox_transforms, only : da_transform_vtox,da_transform_vtox_adj,da_transform_xtoxa,da_transform_xtoxa_adj, &
       da_transform_xtoxa_all,da_transform_xtoxa_adj_all
   use da_wrf_interfaces, only : wrf_dm_bcast_real, wrf_get_dm_communicator
   use module_symbols_util, only : wrfu_finalize
   use da_lapack, only : dsteqr
   use da_wrfvar_io, only : da_med_initialdata_input
   use da_transfer_model, only : da_transfer_wrftoxb

   implicit none

    include 'mpif.h'

   private :: da_dot, da_dot_cv

contains
      
subroutine da_calculate_j(it, iter, cv_size, cv_size_jb, cv_size_je, cv_size_jp, &
                           cv_size_jl, xbx, be, iv, xhat, cv, &
                           re, y, j, grad, grid, config_flags                     )

   !---------------------------------------------------------------------------
   ! Purpose: Initialises the Y-array
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)                :: it     ! external iteration #.
   integer, intent(in)                :: iter   ! internal iteration #.
   integer, intent(in)                :: cv_size    ! Total cv size.
   integer, intent(in)                :: cv_size_jb ! Jb cv size.
   integer, intent(in)                :: cv_size_je ! Je cv size.
   integer, intent(in)                :: cv_size_jp ! Jp cv size.
   integer, intent(in)                :: cv_size_jl ! Jl cv size.
   type (xbx_type),intent(inout)      :: xbx    ! For header & non-grid arrays.
   type (be_type), intent(in)         :: be     ! background error structure.
   type (iv_type), intent(inout)      :: iv     ! innovation vector (o-b).
   real, intent(in)                   :: xhat(1:cv_size) ! control variables.
   real, intent(in)                   :: cv(1:cv_size)   ! control variables.
   type (y_type) , intent(inout)      :: re     ! residual vector (o-a).
   type (y_type) , intent(inout)      :: y      ! y = H(x_inc).
   type (j_type) , intent(out)        :: j      ! cost function j
   real, intent(out)                  :: grad(cv_size)        ! gradient of cost function

   type(domain), intent(inout)  :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   integer          :: je_start, je_end             ! Start/end indices of Je.
   integer          :: jl_start, jl_end             ! Start/end indices of Je.
   real             :: jo_partial                   ! jo for this processor
   type (y_type)    :: jo_grad_y ! Grad_y(jo)
   real             :: cv_xhat_jb(cv_size_jb), cv_xhat_je(cv_size_je), cv_xhat_jl(cv_size_jl)
   integer          :: mz(7)
   integer          :: ndynopt
   real             :: dtemp1x
   integer          :: i, jj, k
   real             :: subarea, whole_area

   ! Variables for VarBC background constraint
   real                              :: cv_xhat_jp(cv_size_jp) ! Jp control variable.
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                           :: inst, ichan, npred, ipred, id
   real                              :: bgerr, gnorm_jp  
    
   integer                           :: n, cldtoplevel(1), icld, nclouds, ncv, minlev_cld
   real                              :: jd_local
   real                              :: js_local
   real, allocatable                 :: cc(:)
   
   if (trace_use) call da_trace_entry("da_calculate_j")

   !-------------------------------------------------------------------------
   ! [0.0] initialization:
   !-------------------------------------------------------------------------
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be % ne /)
   je_start = cv_size_jb + 1
   je_end   = cv_size_jb + cv_size_je
   jp_start = cv_size_jb + cv_size_je + 1
   jp_end   = cv_size_jb + cv_size_je + cv_size_jp
   jl_start = cv_size_jb + cv_size_je + cv_size_jp + 1
   jl_end =   cv_size_jb + cv_size_je + cv_size_jp + cv_size_jl

   call da_allocate_y(iv, jo_grad_y)

   !-------------------------------------------------------------------------
   ! [1.0] calculate jo:
   !-------------------------------------------------------------------------

   ! [1.1] transform from control variable to model grid space:

   if (iter > 0) &
      call da_transform_vtoy(cv_size, be, grid%ep, xhat, iv, grid%vp, grid%vv,&
                              grid%vp6, grid%vv6, xbx, y, &
                              grid, config_flags                      )

   ! [1.2] compute residual (o-a) = (o-b) - h x~

   call da_calculate_residual(iv, y, re)

   ! [1.3] calculate jo:

   call da_jo_and_grady(iv, re, jo_partial, j % jo, jo_grad_y)

   if (test_dm_exact) then
      ! jo_partial has been already summed at lower level
      j % jo % total = jo_partial
   else
      j % jo % total = wrf_dm_sum_real(jo_partial)
   end if

   ! [1.4] calculate jc-dfi:

   j % jc = 0.0

   if ( var4d .and. (grid%jcdfi_use .or. grid%jcdfi_diag == 1) .and. iter > 0 ) then

   end if

   !-------------------------------------------------------------------------
   ! [2.0] calculate jb:
   !-------------------------------------------------------------------------

   j % jb = 0.0
   if (cv_size_jb > 0) then
      cv_xhat_jb(1:cv_size_jb) = cv(1:cv_size_jb) + xhat(1:cv_size_jb)
      j % jb = jb_factor * 0.5 * da_dot_cv(cv_size_jb,  cv_xhat_jb, cv_xhat_jb, grid, mz)
   end if

   !-------------------------------------------------------------------------
   ! [3.0] calculate je:
   !-------------------------------------------------------------------------

   j % je = 0.0
   if (be % ne > 0) then
      cv_xhat_je(1:cv_size_je) = cv(je_start:je_end) + xhat(je_start:je_end)
      j % je = je_factor * 0.5 * da_dot_cv(cv_size_je, cv_xhat_je, cv_xhat_je, grid, mz)
   end if


   !----------------------------------------------------------------------
   ![1.0.1] calculate grad_v (jd):
   !----------------------------------------------------------------------

   j % jd = 0.0

   if (use_wpec) then

      if (var4d) call da_error("da_calculate_j.inc",167,(/'Cannot use 4dvar with dynamic constraint'/))
      if (wpec_factor <= 0) call da_error("da_calculate_j.inc",168,(/'"wpec_factor" for dynamic constraint must be greater than zero'/))

      grid%xa%grad_p_x(:,:,:)=0.0
      grid%xa%grad_p_y(:,:,:)=0.0

      call da_transform_vtod_wpec(cv_size, be, grid%ep, xhat+cv, grid%vp, grid%vv, xbx, grid)

      do i=its,ite
         do jj=jts,jte
            do k=kts,kte
               j % jd = j % jd + 0.5*(grid%xa%grad_p_x(i,jj,k)**2+grid%xa%grad_p_y(i,jj,k)**2)/wpec_factor
            end do
         end do
      end do

      jd_local = j % jd
      ! summation across processors:
      j % jd  = wrf_dm_sum_real(jd_local)

   end if

   !-------------------------------------------------------------------------
   ! [4.0] calculate jl:
   !-------------------------------------------------------------------------
   j % jl = 0.0
   if ( var4d ) then
      cv_xhat_jl(1:cv_size_jl) = cv (jl_start:jl_end) + xhat(jl_start:jl_end)

      j % jl = 0.5 * da_dot_cv(cv_size_jl, cv_xhat_jl, cv_xhat_jl, grid, mz)

   endif

   !-------------------------------------------------------------------------
   ! [5.0] calculate jp:
   !-------------------------------------------------------------------------
   j % jp = 0.0
   if (use_varbc .and. cv_size_jp > 0) then
      cv_xhat_jp = 0.0
      do inst = 1, iv % num_inst   
         do ichan = 1, iv%instid(inst)%nchan
            npred    = iv%instid(inst)%varbc(ichan)%npred
            if (npred <= 0) cycle               !! VarBC channels only	 
            do ipred = 1, npred
               id     = iv%instid(inst)%varbc(ichan)%index(ipred)
	       bgerr  = iv%instid(inst)%varbc(ichan)%bgerr(ipred)
	       if (bgerr > 0.0) &
    	          cv_xhat_jp(id-jp_start+1) = (1/sqrt(bgerr)) * &
	             SUM((cv(id)+xhat(id)) * iv%instid(inst)%varbc(ichan)%vtox(ipred,1:npred))            
	    end do
         end do
      end do
      j % jp = 0.5 * da_dot(cv_size_jp, cv_xhat_jp, cv_xhat_jp)
   end if

   !-------------------------------------------------------------------------
   ! [6.0] calculate js:
   !-------------------------------------------------------------------------
   j % js = 0.0
   if (ANY(use_satcv)) then
      do inst = 1, iv % num_inst   
         do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
         ! Skin Temperature
         !-----------------
	    if (use_satcv(1)) then
               j % js = j % js + 0.5 * xhat(iv%instid(inst)%cv_index(n)%ts) **2
	       
!	       !!! Super-TMP dump of Tskin increment for plotting purposes
!               if (iter > 0) iv%instid(inst)%tb_xb(1,n)  = xhat(iv%instid(inst)%cv_index(n)%ts) 
	    end if	 
	    
         ! Cloud cover(s)
         !---------------
	    if (use_satcv(2)) then
	    j % js = j % js + 0.5 * SUM( xhat(iv%instid(inst)%cv_index(n)%cc) **2)

	    j % js = j % js + 0.5 * SUM( (10.0 * xhat(iv%instid(inst)%cv_index(n)%cc)) **2,      &
	                                  MASK = xhat(iv%instid(inst)%cv_index(n)%cc) < 0.0 .or. &
				                 xhat(iv%instid(inst)%cv_index(n)%cc) > 1.0 )

	       if (iter > 0) then
	          nclouds = iv%instid(inst)%cv_index(n)%nclouds
     	          ncv     = iv%instid(inst)%cv_index(n)%ncv
		  allocate(cc(nclouds))

		  cc = xhat(iv%instid(inst)%cv_index(n)%cc)
	       !---------------------------------------------------------------
               ! Change of variable (preconditioning) 
               !---------------------------------------------------------------
!		  do icld = 1, nclouds
!    	             cc(icld) = SUM( xhat(iv%instid(inst)%cv_index(n)%cc) * &
!	                                iv%instid(inst)%cv_index(n)%vtox(icld,1:ncv) )
!	          end do
		  
	          if (use_satcv(1)) then
		     write (*, '(i6,100F8.2)')n,xhat(iv%instid(inst)%cv_index(n)%ts), SUM(cc)*100, cc*100
		  else
		     write (*, '(i6,100F8.2)')n,SUM(cc)*100, cc*100						  
                  end if
		  
!		  !!! Super-TMP dump of Cloud Cover increment for plotting purposes	 
!                  iv%instid(inst)%tb_inv(1,n) = SUM(cc)*100.0 
!                  
!		  !!! Super-TMP dump of Cloud Top Pressure for plotting purposes
!		  minlev_cld = 5
!		  if (ANY(cc(minlev_cld:nclouds) > 0.01)) then
!		     cldtoplevel = MINLOC(cc(minlev_cld:nclouds), MASK = cc(minlev_cld:nclouds) > 0.01)
!		  else
!		     cldtoplevel = nclouds
!		  end if   
!		  cldtoplevel = cldtoplevel + kte - nclouds !!!+ minlev_cld
!!                  if (rtm_option == rtm_option_rttov) then
!!                     re%instid(inst)%tb(1,n) = coefs(inst)%ref_prfl_p(cldtoplevel(1))
!!                  elseif (rtm_option == rtm_option_crtm) then
!                     re%instid(inst)%tb(1,n) = iv%instid(inst)%pm(cldtoplevel(1),n)
!!                  end if  	    
		  
		  deallocate(cc)
	       end if    
	    end if
	 end do
      end do	      
      js_local = j % js
      ! summation across processors:
      j % js = wrf_dm_sum_real(js_local)
   end if

   !-------------------------------------------------------------------------
   ! [7.0] calculate total cost function j = jo + jb + jc + je + jd + jp + js:
   !-------------------------------------------------------------------------

   j % total = j % jb + j % jo % total + j % je + j % jd + j % jp + j % js
   if (grid%jcdfi_use) j % total = j % total  + j % jc
   if (var4d) j % total = j % total  + j % jl

   !-------------------------------------------------------------------------
   ! [8.0] write cost function:
   !-------------------------------------------------------------------------
   if (rootproc) then
      if (it == 1 .and. iter == 0) then
         write(unit=cost_unit,fmt='(a)')'Outer    EPS     Inner      J           Jb       Jo           Jc         Je         Jd         Jp         Js        jl'
         write(unit=cost_unit,fmt='(a)')'Iter             Iter                            '
         write(unit=grad_unit,fmt='(a)')'Outer    EPS     Inner      G           Gb       Go           Ge         Gd         Gp         Gs        Gl'
         write(unit=grad_unit,fmt='(a)')'Iter             Iter                            '
      end if

      write(unit=cost_unit,fmt='(2x,i2,1x,e10.3,2x,i4,9(1x,f10.3))') &
         it, EPS(it), iter, j % total, j % jb, j % jo % total, j % jc, j % je, j % jd, j % jp, j%js, j%jl
   end if
         
   !-------------------------------------------------------------------------
   ! [9.0] Calculate Gradient:
   !-------------------------------------------------------------------------
   call da_calculate_gradj(it,iter,cv_size,cv_size_jb,cv_size_je,cv_size_jp, &
                           cv_size_jl,xbx, be, iv, xhat+cv, y, grad, grid, config_flags, re)

   call da_deallocate_y (jo_grad_y)
if (trace_use) call da_trace_exit("da_calculate_j")

end subroutine da_calculate_j

subroutine da_calculate_gradj(it, iter, cv_size, cv_size_jb, cv_size_je, cv_size_jp, &
                              cv_size_jl, xbx, be, iv, cv, y, grad, grid, config_flags, re )

   !---------------------------------------------------------------------------
   ! Purpose: Calculates the gradient of the cost function w/r to cv
   !
   ! Called from da_minimise_cg (or da_minimise_lz)
   !
   ! History: 12/12/08 - Creation from da_calculate_j (Tom Auligne)
   !
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)                :: it     ! external iteration #.
   integer, intent(in)                :: iter   ! internal iteration #.
   integer, intent(in)                :: cv_size    ! Total cv size.
   integer, intent(in)                :: cv_size_jb, cv_size_je, cv_size_jp, cv_size_jl
   type (xbx_type),intent(inout)      :: xbx    ! For header & non-grid arrays.
   type (be_type), intent(in)         :: be     ! background error structure.
   type (iv_type), intent(inout)      :: iv     ! innovation vector (o-b).
   real, intent(in)                   :: cv     (1:cv_size)   ! control variables.
   type (y_type), intent(inout)       :: y
   real, intent(out)                  :: grad(cv_size)        ! gradient of cost function
   type (y_type), optional, intent(inout) :: re     ! residual vector (o-a).
   
   type(domain), intent(inout)  :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   integer          :: je_start, je_end             ! Start/end indices of Je.
   integer          :: jl_start, jl_end             ! Start/end indices of Jl.
   real             :: jo_partial                   ! jo for this processor
   type (y_type)    :: jo_grad_y                    ! Grad_y(jo)
   integer          :: mz(7)
   real             :: grad_jo(cv_size)
   real             :: grad_jb(cv_size)
   real             :: grad_je(cv_size)
   real             :: grad_jd(cv_size)
   real             :: grad_jp(cv_size)
   real             :: grad_js(cv_size)
   real             :: grad_jl(cv_size)
   real             :: gnorm_j, gnorm_jo, gnorm_jb, gnorm_je, gnorm_jd, gnorm_jp, gnorm_js, gnorm_jl
   logical          :: jcdf_flag

   ! Variables for VarBC background constraint
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                           :: inst, ichan, npred, ipred, id
   real                              :: bgerr
   integer                           :: n

   if (trace_use) call da_trace_entry("da_calculate_gradj")

   !-------------------------------------------------------------------------
   ! [0.0] initialization:
   !-------------------------------------------------------------------------
  mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz,be%alpha%mz, be % ne /)
   je_start   = cv_size_jb + 1
   je_end     = cv_size_jb + cv_size_je
   jp_start   = cv_size_jb + cv_size_je + 1
   jp_end     = cv_size_jb + cv_size_je + cv_size_jp
   jl_start   = cv_size_jb + cv_size_je + cv_size_jp + 1
   jl_end     = cv_size_jb + cv_size_je + cv_size_jp + cv_size_jl
   
   grad_jo = 0.0
   grad_jb = 0.0
   grad_je = 0.0
   grad_jd = 0.0
   grad_jp = 0.0
   grad_js = 0.0
   grad_jl = 0.0

   jcdf_flag = .false.

   !-------------------------------------------------------------------------
   ! [1.0] calculate grad_v (jo):
   !-------------------------------------------------------------------------
   call da_allocate_y(iv, jo_grad_y)
   
   if (present(re)) then
      call da_calculate_grady(iv, re, jo_grad_y)   
      if ( iter > 0 .and. test_gradient ) jcdf_flag = .true.
      call da_transform_vtoy_adj(cv_size, be, grid%ep, grad_jo, iv, &
              grid%vp, grid%vv, grid%vp6, grid%vv6, xbx, jo_grad_y, grid, config_flags, jcdf_flag)
   else
      call da_transform_vtoy(cv_size, be, grid%ep, cv, iv, grid%vp, &
              grid%vv, grid%vp6, grid%vv6, xbx, y, grid, config_flags)
      call da_calculate_grady(iv, y, jo_grad_y)   
      call da_transform_vtoy_adj(cv_size, be, grid%ep, grad_jo, iv, &
              grid%vp, grid%vv, grid%vp6, grid%vv6, xbx, jo_grad_y, grid, config_flags, .true.)
      grad_jo = - grad_jo    !! Compensate for sign in calculation of grad_v (Jo)
   end if
      
   call da_deallocate_y(jo_grad_y)

   !-------------------------------------------------------------------------
   ! [2.0] calculate grad_v (jb):
   !-------------------------------------------------------------------------
   if (cv_size_jb > 0) grad_jb(1:cv_size_jb) = jb_factor * cv(1:cv_size_jb)

   !-------------------------------------------------------------------------
   ! [3.0] calculate grad_v (je):
   !-------------------------------------------------------------------------
   if (cv_size_je > 0) grad_je(je_start:je_end) = je_factor * cv(je_start:je_end)
   
   !----------------------------------------------------------------------
   ! [3.1] calculate grad_v (jd):
   !----------------------------------------------------------------------
   if (use_wpec) then

      if (var4d) call da_error("da_calculate_gradj.inc",118,(/'Cannot use 4dvar with dynamic constraint'/))
      if (wpec_factor <= 0) call da_error("da_calculate_gradj.inc",119,(/'"wpec_factor" for dynamic constraint must be greater than zero'/))

      grid%xa%grad_p_x(:,:,:)=0.0
      grid%xa%grad_p_y(:,:,:)=0.0

      call da_transform_vtod_wpec(cv_size, be, grid%ep, cv, grid%vp, grid%vv, xbx, grid)

      grid%xa%grad_p_x=(grid%xa%grad_p_x)/wpec_factor
      grid%xa%grad_p_y=(grid%xa%grad_p_y)/wpec_factor

      call da_transform_vtod_wpec_adj(cv_size, be, grid%ep, grad_jd, grid%vp, grid%vv, xbx, grid)

   end if

   !-------------------------------------------------------------------------
   ! [4.0] calculate grad_v (jp):
   !-------------------------------------------------------------------------
   if (use_varbc .and. cv_size_jp > 0) then
      do inst = 1, iv % num_inst   
         do ichan = 1, iv%instid(inst)%nchan
            npred    = iv%instid(inst)%varbc(ichan)%npred
            if (npred <= 0) cycle               !! VarBC channels only	 
            do ipred = 1, npred
               id     = iv%instid(inst)%varbc(ichan)%index(ipred)
	       bgerr  = iv%instid(inst)%varbc(ichan)%bgerr(ipred)
	       if (bgerr > 0.0) &
                  grad_jp(id) = (1/sqrt(bgerr)) * &
                     SUM(cv(id) * iv%instid(inst)%varbc(ichan)%vtox(ipred,1:npred))
	    end do
         end do
      end do
   end if
      
   !-------------------------------------------------------------------------
   ! [5.0] calculate grad_v (js):
   !-------------------------------------------------------------------------
   if (ANY(use_satcv)) then
      do inst = 1, iv % num_inst   
         do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
         ! Skin Temperature
         !-----------------
	    if (use_satcv(1)) &
            grad_js(iv%instid(inst)%cv_index(n)%ts) = cv(iv%instid(inst)%cv_index(n)%ts)
	    
         ! Cloud cover(s)
         !---------------
	    if (use_satcv(2)) then
	    grad_js(iv%instid(inst)%cv_index(n)%cc) = cv(iv%instid(inst)%cv_index(n)%cc)

	    WHERE (cv(iv%instid(inst)%cv_index(n)%cc) < 0.0 .or.                                &
	           cv(iv%instid(inst)%cv_index(n)%cc) > 1.0 )                                   &
	    grad_js(iv%instid(inst)%cv_index(n)%cc) = grad_js(iv%instid(inst)%cv_index(n)%cc) + &
                                                       10.0 * cv(iv%instid(inst)%cv_index(n)%cc)
            end if
	 end do
      end do	      
   end if

   !-------------------------------------------------------------------------
   ! [6.0] calculate grad_v (jl):
   !-------------------------------------------------------------------------
   if (cv_size_jl > 0) grad_jl(jl_start:jl_end) = cv(jl_start:jl_end)

   !--------------------------------------------------------------------------------------------------
   ! [7.0] calculate grad_v (j) = grad_v (jb) + grad_v (jo) + grad_v (je) + grad_v (jd) + grad_v (jp) + grad_v (js) + grad_v (jl)
   !--------------------------------------------------------------------------------------------------   
   grad = grad_jo + grad_jb + grad_je + grad_jd + grad_jp + grad_js + grad_jl

   !-------------------------------------------------------------------------
   ! [8.0] write Gradient:
   !-------------------------------------------------------------------------
   if (present(re)) then
      gnorm_j  = sqrt(da_dot_cv(cv_size, grad,    grad,    grid, mz, jp_start, jp_end))
      gnorm_jo = sqrt(da_dot_cv(cv_size, grad_jo, grad_jo, grid, mz))
      gnorm_jb = sqrt(da_dot_cv(cv_size, grad_jb, grad_jb, grid, mz))
      gnorm_je = sqrt(da_dot_cv(cv_size, grad_je, grad_je, grid, mz))
      gnorm_jd = sqrt(da_dot_cv(cv_size, grad_jd, grad_jd, grid, mz))
      gnorm_jp = sqrt(da_dot_cv(cv_size, grad_jp, grad_jp, grid, mz, jp_start, jp_end))
      gnorm_js = sqrt(da_dot_cv(cv_size, grad_js, grad_js, grid, mz))
      gnorm_jl = sqrt(da_dot_cv(cv_size, grad_jl, grad_jl, grid, mz))
 
      if (rootproc) &
         write(grad_unit,fmt='(2x,i2,1x,e10.3,2x,i4,8(1x,f10.3))') & 
               it, eps(it), iter, gnorm_j, gnorm_jb, gnorm_jo, gnorm_je, gnorm_jd, gnorm_jp, gnorm_js, gnorm_jl
   end if
   if (trace_use) call da_trace_exit("da_calculate_gradj")

end subroutine da_calculate_gradj


subroutine da_jo_and_grady(iv, re, jot, jo, jo_grad_y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)   :: iv          ! Innovation vector (O-B).
   type (y_type),  intent(in)   :: re          ! Residual vector (O-A).
   real,           intent(out)  :: jot         ! Obs cost function.
   type (jo_type), intent(out)  :: jo          ! Obs cost function.
   type (y_type),  intent(out)  :: jo_grad_y   ! Grad_y(Jo)


   real    :: jo_sound, jo_sonde_sfc,jo_synop, jo_geoamv, jo_polaramv, &
              jo_airep, jo_pilot, jo_satem, &
              jo_metar, jo_ships, jo_gpspw, &
              jo_ssmi_tb, jo_ssmi_rv, jo_ssmt1, jo_ssmt2, &
              jo_pseudo, jo_qscat, jo_buoy, &
              jo_profiler, jo_radar, jo_gpsref, jo_bogus, jo_rain, &
              jo_radiance, jo_airsr, jo_mtgirs, jo_tamdar, jo_tamdar_sfc
   integer :: i,k

   if (trace_use) call da_trace_entry("da_jo_and_grady")

   !-------------------------------------------------------------------------
   ! [1.0] Compute components of Grad_y(Jo):
   !-------------------------------------------------------------------------


   if (iv%info(sound)%nlocal > 0) then
      call da_jo_and_grady_sound(iv, re, jo, jo_grad_y)
      jo_sound = jo%sound_u + jo%sound_v + jo%sound_t + jo%sound_q

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_sound          ', jo_sound, &
            '   jo%sound_u      ', jo%sound_u, &
            '   jo%sound_v      ', jo%sound_v, &
            '   jo%sound_t      ', jo%sound_t, &
            '   jo%sound_q      ', jo%sound_q
      end if

   else
      jo % sound_u = 0.0
      jo % sound_v = 0.0
      jo % sound_t = 0.0
      jo % sound_q = 0.0
      jo_sound     = 0.0
   end if

   if (iv%info(sonde_sfc)%nlocal > 0) then
      call da_jo_and_grady_sonde_sfc(iv, re, jo, jo_grad_y)
      jo_sonde_sfc = jo%sonde_sfc_u + jo%sonde_sfc_v + jo%sonde_sfc_t + &
         jo%sonde_sfc_q + jo%sonde_sfc_p

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_sonde_sfc      ', jo_sonde_sfc,     &
            '   jo%sonde_sfc_u  ', jo%sonde_sfc_u, &
            '   jo%sonde_sfc_v  ', jo%sonde_sfc_v, &
            '   jo%sonde_sfc_t  ', jo%sonde_sfc_t, &
            '   jo%sonde_sfc_p  ', jo%sonde_sfc_p, &
            '   jo%sonde_sfc_q  ', jo%sonde_sfc_q
      end if
   else
      jo % sonde_sfc_u = 0.0
      jo % sonde_sfc_v = 0.0
      jo % sonde_sfc_t = 0.0
      jo % sonde_sfc_p = 0.0
      jo % sonde_sfc_q = 0.0
      jo_sonde_sfc = 0.0
   end if



   if (iv%info(mtgirs)%nlocal > 0) then
      call da_jo_and_grady_mtgirs(iv, re, jo, jo_grad_y)
      jo_mtgirs = jo%mtgirs_u + jo%mtgirs_v + jo%mtgirs_t + jo%mtgirs_q

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_mtgirs          ', jo_mtgirs, &
            '   jo%mtgirs_u      ', jo%mtgirs_u, &
            '   jo%mtgirs_v      ', jo%mtgirs_v, &
            '   jo%mtgirs_t      ', jo%mtgirs_t, &
            '   jo%mtgirs_q      ', jo%mtgirs_q
      end if

   else
      jo % mtgirs_u = 0.0
      jo % mtgirs_v = 0.0
      jo % mtgirs_t = 0.0
      jo % mtgirs_q = 0.0
      jo_mtgirs     = 0.0
   end if

   if (iv%info(tamdar)%nlocal > 0) then
      call da_jo_and_grady_tamdar(iv, re, jo, jo_grad_y)
      jo_tamdar = jo%tamdar_u + jo%tamdar_v + jo%tamdar_t + jo%tamdar_q

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_tamdar          ', jo_tamdar, &
            '   jo%tamdar_u      ', jo%tamdar_u, &
            '   jo%tamdar_v      ', jo%tamdar_v, &
            '   jo%tamdar_t      ', jo%tamdar_t, &
            '   jo%tamdar_q      ', jo%tamdar_q
      end if
   else
      jo % tamdar_u = 0.0
      jo % tamdar_v = 0.0
      jo % tamdar_t = 0.0
      jo % tamdar_q = 0.0
      jo_tamdar     = 0.0
   end if

   if (iv%info(tamdar_sfc)%nlocal > 0) then
      call da_jo_and_grady_tamdar_sfc(iv, re, jo, jo_grad_y)
      jo_tamdar_sfc = jo%tamdar_sfc_u + jo%tamdar_sfc_v + jo%tamdar_sfc_t + &
         jo%tamdar_sfc_q + jo%tamdar_sfc_p

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_tamdar_sfc      ', jo_tamdar_sfc,     &
            '   jo%tamdar_sfc_u  ', jo%tamdar_sfc_u, &
            '   jo%tamdar_sfc_v  ', jo%tamdar_sfc_v, &
            '   jo%tamdar_sfc_t  ', jo%tamdar_sfc_t, &
            '   jo%tamdar_sfc_p  ', jo%tamdar_sfc_p, &
            '   jo%tamdar_sfc_q  ', jo%tamdar_sfc_q
      end if

   else
      jo % tamdar_sfc_u = 0.0
      jo % tamdar_sfc_v = 0.0
      jo % tamdar_sfc_t = 0.0
      jo % tamdar_sfc_p = 0.0
      jo % tamdar_sfc_q = 0.0
      jo_tamdar_sfc     = 0.0

   end if

   if (iv%info(synop)%nlocal > 0) then
      call da_jo_and_grady_synop(iv, re, jo, jo_grad_y)
      jo_synop = jo%synop_u + jo%synop_v + jo%synop_t + jo%synop_p + jo%synop_q

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_synop          ', jo_synop, &
            '   jo%synop_u      ', jo%synop_u, &
            '   jo%synop_v      ', jo%synop_v, &
            '   jo%synop_t      ', jo%synop_t, &
            '   jo%synop_p      ', jo%synop_p, &
            '   jo%synop_q      ', jo%synop_q
      end if
   else
      jo % synop_u = 0.0
      jo % synop_v = 0.0
      jo % synop_t = 0.0
      jo % synop_p = 0.0
      jo % synop_q = 0.0
      jo_synop = 0.0
   end if

   if (iv%info(geoamv)%nlocal > 0) then
      call da_jo_and_grady_geoamv(iv, re, jo, jo_grad_y)
      jo_geoamv = jo%geoamv_u + jo%geoamv_v
      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_geoamv         ', jo_geoamv, &
            '   jo%geoamv_u     ', jo%geoamv_u, &
            '   jo%geoamv_v     ', jo%geoamv_v
      end if
   else
      jo % geoamv_u = 0.0
      jo % geoamv_v = 0.0
      jo_geoamv = 0.0
   end if

   if (iv%info(polaramv)%nlocal > 0) then
      call da_jo_and_grady_polaramv(iv, re, jo, jo_grad_y)
      jo_polaramv = jo%polaramv_u + jo%polaramv_v
      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_polaramv       ', jo_polaramv, &
            '   jo%polaramv_u   ', jo%polaramv_u, &
            '   jo%polaramv_v   ', jo%polaramv_v
      end if
   else
      jo % polaramv_u = 0.0
      jo % polaramv_v = 0.0
      jo_polaramv = 0.0
   end if

   if (iv%info(airep)%nlocal > 0) then
      call da_jo_and_grady_airep(iv, re, jo, jo_grad_y)
      jo_airep = jo%airep_u + jo%airep_v + jo%airep_t + jo%airep_q

      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_airep          ', jo_airep, &
               '   jo%airep_u      ', jo%airep_u, &
               '   jo%airep_v      ', jo%airep_v, &
               '   jo%airep_t      ', jo%airep_t, &
               '   jo%airep_q      ', jo%airep_q
      end if
   else
      jo % airep_u = 0.0
      jo % airep_v = 0.0
      jo % airep_t = 0.0
      jo % airep_q = 0.0
      jo_airep = 0.0
   end if

   if (iv%info(pilot)%nlocal > 0) then
      call da_jo_and_grady_pilot(iv, re, jo, jo_grad_y)
      jo_pilot = jo%pilot_u + jo%pilot_v
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_pilot          ', jo_pilot, &
               '   jo%pilot_u      ', jo%pilot_u, &
               '   jo%pilot_v      ', jo%pilot_v
      end if
   else
      jo % pilot_u = 0.0
      jo % pilot_v = 0.0
      jo_pilot = 0.0
   end if

   if (iv%info(satem)%nlocal > 0) then
      call da_jo_and_grady_satem(iv, re, jo, jo_grad_y)
      jo_satem = jo%satem_thickness
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_satem          ', jo_satem, &
               '   jo%satem_thckns ', jo%satem_thickness
      end if
   else
      jo % satem_thickness = 0.0
      jo_satem = 0.0
   end if

   if (iv%info(metar)%nlocal > 0) then
      call da_jo_and_grady_metar(iv, re, jo, jo_grad_y)
      jo_metar = jo%metar_u + jo%metar_v + jo%metar_t + jo%metar_p + jo%metar_q
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_metar          ', jo_metar, &
               '   jo%metar_u      ', jo%metar_u, &
               '   jo%metar_v      ', jo%metar_v, &
               '   jo%metar_t      ', jo%metar_t, &
               '   jo%metar_p      ', jo%metar_p, &
               '   jo%metar_q      ', jo%metar_q
      end if
   else
      jo % metar_u = 0.0
      jo % metar_v = 0.0
      jo % metar_t = 0.0
      jo % metar_p = 0.0
      jo % metar_q = 0.0
      jo_metar = 0.0
   end if

   if (iv%info(ships)%nlocal > 0) then
      call da_jo_and_grady_ships(iv, re, jo, jo_grad_y)
      jo_ships = jo%ships_u + jo%ships_v + jo%ships_t + jo%ships_p + jo%ships_q
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_ships          ', jo_ships, &
               '   jo%ships_u      ', jo%ships_u, &
               '   jo%ships_v      ', jo%ships_v, &
               '   jo%ships_t      ', jo%ships_t, &
               '   jo%ships_p      ', jo%ships_p, &
               '   jo%ships_q      ', jo%ships_q
      end if
   else
      jo % ships_u = 0.0
      jo % ships_v = 0.0
      jo % ships_t = 0.0
      jo % ships_p = 0.0
      jo % ships_q = 0.0
      jo_ships = 0.0
   end if

   if (iv%info(gpspw)%nlocal > 0) then
      call da_jo_and_grady_gpspw(iv, re, jo, jo_grad_y)
      jo_gpspw = jo%gpspw_tpw
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_gpspw          ', jo_gpspw, &
               '   jo%gpspw_tpw    ', jo%gpspw_tpw
      end if
   else
      jo % gpspw_tpw = 0.0
      jo_gpspw = 0.0
   end if

   if (iv%info(gpsref)%nlocal > 0) then
      call da_jo_and_grady_gpsref(iv, re, jo, jo_grad_y)
      jo_gpsref = jo%gpsref_ref
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_gpsref         ', jo_gpsref, &
               '   jo%gpsref_ref   ', jo%gpsref_ref
      end if
   else
      jo % gpsref_ref = 0.0
      jo_gpsref = 0.0
   end if

   if (iv%info(ssmi_tb)%nlocal > 0) then
      call da_jo_and_grady_ssmi_tb (iv, re, jo, jo_grad_y)
      jo_ssmi_tb = jo % ssmi_tb19v + jo % ssmi_tb19h + jo % ssmi_tb22v + &
         jo % ssmi_tb37v + jo % ssmi_tb37h + jo % ssmi_tb85v + &
         jo % ssmi_tb85h 
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_ssmi_tb        ', jo_ssmi_tb, &
               '   jo%ssmi_tb19v   ', jo%ssmi_tb19v, &
               '   jo%ssmi_tb19h   ', jo%ssmi_tb19h, &
               '   jo%ssmi_tb22v   ', jo%ssmi_tb22v, &
               '   jo%ssmi_tb37v   ', jo%ssmi_tb37v, &
               '   jo%ssmi_tb37h   ', jo%ssmi_tb37h, &
               '   jo%ssmi_tb85v   ', jo%ssmi_tb85v, &
               '   jo%ssmi_tb85h   ', jo%ssmi_tb85h
      end if
   else
      jo % ssmi_tb19v = 0.0
      jo % ssmi_tb19h = 0.0
      jo % ssmi_tb22v = 0.0
      jo % ssmi_tb37v = 0.0
      jo % ssmi_tb37h = 0.0
      jo % ssmi_tb85v = 0.0
      jo % ssmi_tb85h = 0.0
      jo_ssmi_tb = 0.0
   end if

   if (iv%info(ssmi_rv)%nlocal > 0) then
      call da_jo_and_grady_ssmi_rv(iv, re, jo, jo_grad_y)
      jo_ssmi_rv = jo % ssmir_speed + jo % ssmir_tpw
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_ssmi_rv        ', jo_ssmi_rv, &
               '   jo%ssmir_speed  ', jo%ssmir_speed, &
               '   jo%ssmir_tpw    ', jo%ssmir_tpw
      end if
   else
      jo % ssmir_speed = 0.0
      jo % ssmir_tpw   = 0.0
      jo_ssmi_rv = 0.0
   end if

   if (iv%info(ssmt1)%nlocal > 0) then
      call da_jo_and_grady_ssmt1(iv, re, jo, jo_grad_y)
      jo_ssmt1 = jo%ssmt1_t
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_ssmt1          ', jo_ssmt1, &
               '   jo%ssmt1_t      ', jo%ssmt1_t
      end if
   else
      jo % ssmt1_t = 0.0
      jo_ssmt1 = 0.0
   end if

   if (iv%info(ssmt2)%nlocal > 0) then
      call da_jo_and_grady_ssmt2(iv, re, jo, jo_grad_y)  
      jo_ssmt2 = jo%ssmt2_rh
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_ssmt2          ', jo_ssmt2, &
               '   jo%ssmt2_rh     ', jo%ssmt2_rh
      end if
   else
      jo % ssmt2_rh = 0.0
      jo_ssmt2 = 0.0
   end if

   if (iv%info(radar)%nlocal > 0) then
      call da_jo_and_grady_radar(iv, re, jo, jo_grad_y)
      jo_radar = jo%radar_rv + jo%radar_rf + jo%radar_rrn + jo%radar_rsn + jo%radar_rgr + jo%radar_rqv
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_radar          ', jo_radar, &
               '   jo%radar_rv     ', jo%radar_rv, &
               '   jo%radar_rf     ', jo%radar_rf, &
               '   jo%radar_rrn    ', jo%radar_rrn, &
               '   jo%radar_rsn    ', jo%radar_rsn, &
               '   jo%radar_rgr    ', jo%radar_rgr, &
               '   jo%radar_rqv    ', jo%radar_rqv
      end if
   else
      jo % radar_rv  = 0.0
      jo % radar_rf  = 0.0
      jo % radar_rrn = 0.0
      jo % radar_rsn = 0.0
      jo % radar_rgr = 0.0
      jo % radar_rqv = 0.0
      jo_radar = 0.0
   end if

   if (iv%info(rain)%nlocal > 0) then
       call da_jo_and_grady_rain(iv, re, jo, jo_grad_y)
       jo_rain = jo%rain_r
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_rain          ', jo_rain, &
               '   jo%rain_r      ', jo%rain_r
      end if
   else
      jo % rain_r = 0.0
      jo_rain = 0.0
   end if

   if (iv%info(pseudo)%nlocal > 0) then
      call da_jo_and_grady_pseudo(iv, re, jo, jo_grad_y)    
      jo_pseudo = jo%pseudo_u + jo%pseudo_v + jo%pseudo_t + jo%pseudo_p + jo%pseudo_q
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_pseudo          ', jo_pseudo, &
               '   jo%pseudo_u      ', jo%pseudo_u, &
               '   jo%pseudo_v      ', jo%pseudo_v, &
               '   jo%pseudo_t      ', jo%pseudo_t, &
               '   jo%pseudo_p      ', jo%pseudo_p, &
               '   jo%pseudo_q      ', jo%pseudo_q
      end if
   else
      jo % pseudo_u = 0.0
      jo % pseudo_v = 0.0
      jo % pseudo_t = 0.0
      jo % pseudo_p = 0.0
      jo % pseudo_q = 0.0
      jo_pseudo = 0.0
   end if
   
   if (iv%info(qscat)%nlocal > 0) then
      call da_jo_and_grady_qscat(iv, re, jo, jo_grad_y)
      jo_qscat = jo%qscat_u + jo%qscat_v
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_qscat           ', jo_qscat, &
               '   jo%qscat_u       ', jo%qscat_u, &
               '   jo%qscat_v       ', jo%qscat_v
      end if
   else
      jo % qscat_u = 0.0
      jo % qscat_v = 0.0
      jo_qscat = 0.0
   end if

   if (iv%info(profiler)%nlocal > 0) then
      call da_jo_and_grady_profiler (iv, re, jo, jo_grad_y)
      jo_profiler = jo%profiler_u + jo%profiler_v
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_profiler        ', jo_profiler, &
               '   jo%profiler_u    ', jo%profiler_u, &
               '   jo%profiler_v    ', jo%profiler_v
      end if
   else
      jo % profiler_u = 0.0
      jo % profiler_v = 0.0
      jo_profiler = 0.0
   end if

   if (iv%info(bogus)%nlocal > 0) then
      call da_jo_and_grady_bogus (iv, re, jo, jo_grad_y)
      jo_bogus = jo%bogus_u + jo%bogus_v + jo%bogus_slp + jo%bogus_t + jo%bogus_q
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_bogus           ', jo_bogus, &
               '   jo%bogus_u       ', jo%bogus_u, &
               '   jo%bogus_v       ', jo%bogus_v, &
               '   jo%bogus_t       ', jo%bogus_t, &
               '   jo%bogus_slp     ', jo%bogus_slp, &
               '   jo%bogus_q       ', jo%bogus_q
      end if
   else
      jo % bogus_u   = 0.0
      jo % bogus_v   = 0.0
      jo % bogus_t   = 0.0
      jo % bogus_q   = 0.0
      jo % bogus_slp = 0.0
      jo_bogus = 0.0
   end if

   if (iv%info(buoy)%nlocal > 0) then
      call da_jo_and_grady_buoy (iv, re, jo, jo_grad_y)
      jo_buoy = jo%buoy_u + jo%buoy_v + jo%buoy_t + jo%buoy_p + jo%buoy_q
      if (print_detail_grad) then
          write(unit=stdout, fmt='(a, e24.12)') &
               ' jo_buoy            ', jo_buoy, &
               '   jo%buoy_u        ', jo%buoy_u, &
               '   jo%buoy_v        ', jo%buoy_v, &
               '   jo%buoy_t        ', jo%buoy_t, &
               '   jo%buoy_p        ', jo%buoy_p, &
               '   jo%buoy_q        ', jo%buoy_q
      end if
   else
      jo % buoy_u = 0.0
      jo % buoy_v = 0.0
      jo % buoy_t = 0.0
      jo % buoy_p = 0.0
      jo % buoy_q = 0.0
      jo_buoy = 0.0
   end if

   jo_radiance = 0.0
   if (iv%num_inst > 0) then
      call da_jo_and_grady_rad (iv, re, jo, jo_grad_y)

      if (use_rad) then
         do i=1,iv%num_inst
            do k=1,iv%instid(i)%nchan
               jo_radiance = jo_radiance + jo%rad(i)%jo_ichan(k)
            end do
         end do
      end if
      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_radiance       ', jo_radiance
         do i = 1, iv%num_inst
            write(unit=stdout, fmt='(a, e24.12)') &
               trim('   jo_'//iv%instid(i)%rttovid_string), sum(jo%rad(i)%jo_ichan(:))
         end do
      end if
   end if

   if (iv%info(airsr)%nlocal > 0) then
      call da_jo_and_grady_airsr(iv, re, jo, jo_grad_y)
      jo_airsr = jo%airsr_t + jo%airsr_q

      if (print_detail_grad) then
         write(unit=stdout, fmt='(a, e24.12)') &
            ' jo_airsr          ', jo_airsr, &
            '   jo%airsr_t      ', jo%airsr_t, &
            '   jo%airsr_q      ', jo%airsr_q
      end if
   else
      jo%airsr_t = 0.0
      jo%airsr_q = 0.0
      jo_airsr = 0.0
   end if

   !-------------------------------------------------------------------------
   ! [2.0] Jo = 1/2 * (yo-y)**2/ob_err_variance:
   !-------------------------------------------------------------------------

   jo%total = jo_sound + jo_sonde_sfc+jo_geoamv + jo_polaramv + jo_synop + jo_satem + &
      jo_pilot + jo_airep + jo_metar + jo_ships + &
      jo_gpspw + jo_ssmi_tb + jo_ssmi_rv + jo_ssmt1 + jo_ssmt2 + &
      jo_pseudo + jo_qscat + jo_profiler + jo_buoy + &
      jo_radar + jo_gpsref + jo_bogus + jo_radiance + jo_airsr + jo_mtgirs + &
      jo_tamdar + jo_tamdar_sfc + jo_rain 

   jot = jo%total

   if (print_detail_grad) then
      write(unit=stdout, fmt='(a, e24.12)') &
         '   jo%total      ', jot

      write(unit=stdout, fmt='(a, e24.12)') &
         '   jo_sound        ', jo_sound, &
         '   jo_sonde_sfc    ', jo_sonde_sfc, &
         '   jo_geoamv       ', jo_geoamv, &
         '   jo_polaramv     ', jo_polaramv, &
         '   jo_synop        ', jo_synop, &
         '   jo_satem        ', jo_satem, &
         '   jo_pilot        ', jo_pilot, &
         '   jo_airep        ', jo_airep, &
         '   jo_metar        ', jo_metar, &
         '   jo_ships        ', jo_ships, &
         '   jo_gpspw        ', jo_gpspw, &
         '   jo_ssmi_tb      ', jo_ssmi_tb, &
         '   jo_ssmi_rv      ', jo_ssmi_rv, &
         '   jo_ssmt1        ', jo_ssmt1, &
         '   jo_ssmt2        ', jo_ssmt2, &
         '   jo_pseudo       ', jo_pseudo, &
         '   jo_qscat        ', jo_qscat, &
         '   jo_profiler     ', jo_profiler, &
         '   jo_buoy         ', jo_buoy, &
         '   jo_radar        ', jo_radar, &
         '   jo_gpsref       ', jo_gpsref, &
         '   jo_bogus        ', jo_bogus,  &
         '   jo_radiance     ', jo_radiance, &
         '   jo_airsr        ', jo_airsr,&
         '   jo_mtgirs       ', jo_mtgirs, &
         '   jo_tamdar       ', jo_tamdar, &
         '   jo_tamdar_sfc   ', jo_tamdar_sfc, &
         '   jo_rain         ', jo_rain

   end if

   if (trace_use) call da_trace_exit("da_jo_and_grady")

end subroutine da_jo_and_grady


subroutine da_calculate_residual(iv, y, re)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residuals
   !-----------------------------------------------------------------------

   implicit none
      
   type (iv_type), intent(inout)     :: iv     ! Innovation vector (O-B).
   type (y_type),  intent(in)        :: y      ! y = H (xa)
   type (y_type),  intent(inout)     :: re     ! Residual (O-A).

   integer    :: np_available, np_obs_used, np_missing, np_bad_data 

   if (trace_use) call da_trace_entry("da_calculate_residual")
      
   np_available = 0
   np_obs_used  = 0
   np_missing   = 0
   np_bad_data  = 0

   !-------------------------------------------------------------------------
   ! [1.0] (O-A) = (O-B) - H x~:
   !-------------------------------------------------------------------------

   if (iv%info(sound)%nlocal > 0) &
      call da_residual_sound(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(sonde_sfc)%nlocal > 0) &
      call da_residual_sonde_sfc(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(mtgirs)%nlocal > 0) & 
      call da_residual_mtgirs(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(tamdar)%nlocal > 0)  &   
      call da_residual_tamdar(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(tamdar_sfc)%nlocal > 0)  &   
      call da_residual_tamdar_sfc(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(synop)%nlocal > 0) &
      call da_residual_synop(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(geoamv)%nlocal > 0) &
      call da_residual_geoamv(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(polaramv)%nlocal > 0) &
      call da_residual_polaramv(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(airep)%nlocal > 0) &
      call da_residual_airep(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(metar)%nlocal > 0) &
      call da_residual_metar(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(ships)%nlocal > 0) &
      call da_residual_ships(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(gpspw)%nlocal > 0) &
      call da_residual_gpspw(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(gpsref)%nlocal > 0) &
      call da_residual_gpsref(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(ssmi_tb)%nlocal > 0) &
      call da_residual_ssmi_tb(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(ssmi_rv)%nlocal > 0) &
      call da_residual_ssmi_rv(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if ( iv%info(ssmt2)%nlocal > 0) &
      call da_residual_ssmt1(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(ssmt2)%nlocal > 0) &
      call da_residual_ssmt2(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(pilot)%nlocal > 0) &
      call da_residual_pilot(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(bogus)%nlocal > 0) &
      call da_residual_bogus(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(satem)%nlocal > 0) &
      call da_residual_satem(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (num_pseudo > 0) &
      call da_residual_pseudo(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(qscat)%nlocal > 0) &
      call da_residual_qscat(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(radar)%nlocal > 0) &
      call da_residual_radar(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(profiler)%nlocal > 0) &
      call da_residual_profiler(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(buoy)%nlocal > 0) &
      call da_residual_buoy(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

    if (iv%info(rain)%nlocal > 0) &
      call da_residual_rain(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%num_inst > 0) &
      call da_residual_rad(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (iv%info(airsr)%nlocal > 0) &
      call da_residual_airsr(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   if (trace_use) call da_trace_exit("da_calculate_residual")

end subroutine da_calculate_residual


subroutine da_get_var_diagnostics(it, iv, j)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)  :: it
   type(iv_type), intent(inout):: iv      ! innovation vector.
   type(j_type), intent(inout) :: j       ! Cost function.

   integer                      :: num_stats_tot
   integer                      :: i,k
   real                         :: jo_radiance
   real                         :: temp(78)

   if (trace_use) call da_trace_entry("da_get_var_diagnostics")

   !--------------------------------------------------------------------------
   ! [1.0] Sum up Jo across processors:
   !--------------------------------------------------------------------------

   num_stats_tot = sum(iv%nstats(:))

   temp(1)  = j % jo % synop_u
   temp(2)  = j % jo % synop_v
   temp(3)  = j % jo % synop_t
   temp(4)  = j % jo % synop_p
   temp(5)  = j % jo % synop_q
   temp(6)  = j % jo % metar_u
   temp(7)  = j % jo % metar_v
   temp(8)  = j % jo % metar_t
   temp(9)  = j % jo % metar_p
   temp(10) = j % jo % metar_q
   temp(11) = j % jo % ships_u
   temp(12) = j % jo % ships_v
   temp(13) = j % jo % ships_t
   temp(14) = j % jo % ships_p
   temp(15) = j % jo % ships_q
   temp(16) = j % jo % geoamv_u
   temp(17) = j % jo % geoamv_v
   temp(18) = j % jo % polaramv_u
   temp(19) = j % jo % polaramv_v      
   temp(20) = j % jo % gpspw_tpw       
   temp(21) = j % jo % gpsref_ref      
   temp(22) = j % jo % sound_u         
   temp(23) = j % jo % sound_v         
   temp(24) = j % jo % sound_t         
   temp(25) = j % jo % sound_q         
   temp(26) = j % jo % sonde_sfc_u     
   temp(27) = j % jo % sonde_sfc_v     
   temp(28) = j % jo % sonde_sfc_t     
   temp(29) = j % jo % sonde_sfc_p     
   temp(30) = j % jo % sonde_sfc_q     
   temp(31) = j % jo % airep_u         
   temp(32) = j % jo % airep_v         
   temp(33) = j % jo % airep_t         
   temp(78) = j % jo % airep_q
   temp(34) = j % jo % pilot_u         
   temp(35) = j % jo % pilot_v         
   temp(36) = j % jo % bogus_u         
   temp(37) = j % jo % bogus_v         
   temp(38) = j % jo % bogus_t         
   temp(39) = j % jo % bogus_q         
   temp(40) = j % jo % bogus_slp       
   temp(41) = j % jo % ssmir_speed     
   temp(42) = j % jo % ssmir_tpw       
   temp(43) = j % jo % ssmi_tb19v      
   temp(44) = j % jo % ssmi_tb19h      
   temp(45) = j % jo % ssmi_tb22v      
   temp(46) = j % jo % ssmi_tb37v      
   temp(47) = j % jo % ssmi_tb37h      
   temp(48) = j % jo % ssmi_tb85v      
   temp(49) = j % jo % ssmi_tb85h      
   temp(50) = j % jo % satem_thickness 
   temp(51) = j % jo % ssmt1_t         
   temp(52) = j % jo % ssmt2_rh        
   temp(53) = j % jo % qscat_u         
   temp(54) = j % jo % qscat_v         
   temp(55) = j % jo % profiler_u      
   temp(56) = j % jo % profiler_v      
   temp(57) = j % jo % buoy_u          
   temp(58) = j % jo % buoy_v          
   temp(59) = j % jo % buoy_t          
   temp(60) = j % jo % buoy_p          
   temp(61) = j % jo % buoy_q          
   temp(62) = j % jo % airsr_t         
   temp(63) = j % jo % airsr_q         
   temp(64) = j % jo % mtgirs_t
   temp(65) = j % jo % mtgirs_q
   temp(66) = j % jo % mtgirs_u
   temp(67) = j % jo % mtgirs_v
   temp(68) = j % jo % tamdar_t
   temp(69) = j % jo % tamdar_q
   temp(70) = j % jo % tamdar_u
   temp(71) = j % jo % tamdar_v
   temp(72) = j % jo % tamdar_sfc_u
   temp(73) = j % jo % tamdar_sfc_v
   temp(74) = j % jo % tamdar_sfc_t
   temp(75) = j % jo % tamdar_sfc_p
   temp(76) = j % jo % tamdar_sfc_q
   temp(77) = j % jo % rain_r

   call da_proc_sum_real(temp(:))

   j % jo % synop_u         = temp(1)  
   j % jo % synop_v         = temp(2)  
   j % jo % synop_t         = temp(3)  
   j % jo % synop_p         = temp(4)  
   j % jo % synop_q         = temp(5)  
   j % jo % metar_u         = temp(6)  
   j % jo % metar_v         = temp(7)  
   j % jo % metar_t         = temp(8)  
   j % jo % metar_p         = temp(9)  
   j % jo % metar_q         = temp(10) 
   j % jo % ships_u         = temp(11) 
   j % jo % ships_v         = temp(12) 
   j % jo % ships_t         = temp(13) 
   j % jo % ships_p         = temp(14) 
   j % jo % ships_q         = temp(15) 
   j % jo % geoamv_u        = temp(16) 
   j % jo % geoamv_v        = temp(17) 
   j % jo % polaramv_u      = temp(18) 
   j % jo % polaramv_v      = temp(19) 
   j % jo % gpspw_tpw       = temp(20) 
   j % jo % gpsref_ref      = temp(21) 
   j % jo % sound_u         = temp(22) 
   j % jo % sound_v         = temp(23) 
   j % jo % sound_t         = temp(24) 
   j % jo % sound_q         = temp(25) 
   j % jo % sonde_sfc_u     = temp(26) 
   j % jo % sonde_sfc_v     = temp(27) 
   j % jo % sonde_sfc_t     = temp(28) 
   j % jo % sonde_sfc_p     = temp(29) 
   j % jo % sonde_sfc_q     = temp(30) 
   j % jo % airep_u         = temp(31) 
   j % jo % airep_v         = temp(32) 
   j % jo % airep_t         = temp(33) 
   j % jo % airep_q         = temp(78) 
   j % jo % pilot_u         = temp(34) 
   j % jo % pilot_v         = temp(35) 
   j % jo % bogus_u         = temp(36) 
   j % jo % bogus_v         = temp(37) 
   j % jo % bogus_t         = temp(38) 
   j % jo % bogus_q         = temp(39) 
   j % jo % bogus_slp       = temp(40) 
   j % jo % ssmir_speed     = temp(41) 
   j % jo % ssmir_tpw       = temp(42) 
   j % jo % ssmi_tb19v      = temp(43) 
   j % jo % ssmi_tb19h      = temp(44) 
   j % jo % ssmi_tb22v      = temp(45) 
   j % jo % ssmi_tb37v      = temp(46) 
   j % jo % ssmi_tb37h      = temp(47) 
   j % jo % ssmi_tb85v      = temp(48) 
   j % jo % ssmi_tb85h      = temp(49) 
   j % jo % satem_thickness = temp(50) 
   j % jo % ssmt1_t         = temp(51) 
   j % jo % ssmt2_rh        = temp(52) 
   j % jo % qscat_u         = temp(53) 
   j % jo % qscat_v         = temp(54) 
   j % jo % profiler_u      = temp(55) 
   j % jo % profiler_v      = temp(56) 
   j % jo % buoy_u          = temp(57) 
   j % jo % buoy_v          = temp(58) 
   j % jo % buoy_t          = temp(59) 
   j % jo % buoy_p          = temp(60) 
   j % jo % buoy_q          = temp(61) 
   j % jo % airsr_t         = temp(62) 
   j % jo % airsr_q         = temp(63) 

   j % jo % mtgirs_t        = temp(64)
   j % jo % mtgirs_q        = temp(65)
   j % jo % mtgirs_u        = temp(66)
   j % jo % mtgirs_v        = temp(67)

   j % jo % tamdar_t        = temp(68)
   j % jo % tamdar_q        = temp(69)
   j % jo % tamdar_u        = temp(70)
   j % jo % tamdar_v        = temp(71)
   j % jo % tamdar_sfc_u    = temp(72)
   j % jo % tamdar_sfc_v    = temp(73)
   j % jo % tamdar_sfc_t    = temp(74)
   j % jo % tamdar_sfc_p    = temp(75)
   j % jo % tamdar_sfc_q    = temp(76)
   j % jo % rain_r          = temp(77)

   if (use_rad) then
      jo_radiance = 0.0
      do i = 1, iv%num_inst                 ! loop for sensor
         call da_proc_sum_ints(j % jo % rad(i)% num_ichan(:))
         call da_proc_sum_real(j % jo % rad(i) % jo_ichan(:))
         jo_radiance = jo_radiance + sum(j % jo % rad(i) % jo_ichan(:))
      end do
   end if

   !-----------------------------------------------------------------------------
   ! [2.0] Print out VAR diagnostics:
   !-----------------------------------------------------------------------------

   if (rootproc) then

      write(unit=stdout,fmt=*) ' '
      write(unit=stdout,fmt='(A)') 'Diagnostics'
      write(unit=stdout,fmt='(A,F12.2)')   '   Final cost function J       = ', j % total
      write(unit=stdout,fmt=*) ' '

      write(unit=stdout,fmt='(a,i8)')    '   Total number of obs.        = ', num_stats_tot
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of J            = ', j % total
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Jo           = ', j % jo % total
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Jd           = ', j % jd
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Jb           = ', j % jb
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Jc           = ', j % jc
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Je           = ', j % je
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Jp           = ', j % jp
      write(unit=stdout,fmt='(a,f15.5)') '   Final value of Jl           = ', j % jl
      if (num_stats_tot > 0) &
         write(unit=stdout,fmt='(a,f15.5)') '   Final J / total num_obs     = ', j % total / &
                                                          real(num_stats_tot)
      if (cv_options /= 3) then
        write(unit=stdout,fmt='(a,(5f15.5))') '   Jb factor used(1)           = ', var_scaling1(it)
        write(unit=stdout,fmt='(a,(5f15.5))') '   Jb factor used(2)           = ', var_scaling2(it)
        write(unit=stdout,fmt='(a,(5f15.5))') '   Jb factor used(3)           = ', var_scaling3(it)
        write(unit=stdout,fmt='(a,(5f15.5))') '   Jb factor used(4)           = ', var_scaling4(it)
        write(unit=stdout,fmt='(a,(5f15.5))') '   Jb factor used(5)           = ', var_scaling5(it)
      endif

      write(unit=stdout,fmt='(a, f15.5)') '   Jb factor used              = ', jb_factor
      write(unit=stdout,fmt='(a, f15.5)') '   Je factor used              = ', je_factor
      write(unit=stdout,fmt='(a, f15.5)') '   VarBC factor used           = ', varbc_factor
      write(unit=stdout,fmt=*) ' '

      if (use_rad) then
         write(unit=stdout,fmt='(a,i8)')    '   Total number of radiances    = ', iv%nstats(radiance)
         write(unit=stdout,fmt='(a,f15.5)') '   Cost function for radiances  = ', jo_radiance
         write(unit=stdout,fmt=*) ' '
      end if

      ! [4.2] Output components of Jo:

      if (iv%info(synop)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    synop obs, Jo(actual)  = ', &
                                  iv%info(synop)%ntotal, iv%nstats(synop), &
                                  j % jo % synop_u, iv % synop_ef_u, &
                                  j % jo % synop_v, iv % synop_ef_v, &
                                  j % jo % synop_t, iv % synop_ef_t, &
                                  j % jo % synop_p, iv % synop_ef_p, &
                                  j % jo % synop_q, iv % synop_ef_q

      end if

      if (trace_use) call da_trace("da_get_var_diagnostics", &
         message="Memory increase from internal write")

      if (iv%info(metar)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    metar obs, Jo(actual)  = ', &
                               iv%info(metar)%ntotal, iv%nstats(metar), &
                               j % jo % metar_u, iv % metar_ef_u, &
                               j % jo % metar_v, iv % metar_ef_v, &
                               j % jo % metar_t, iv % metar_ef_t, &
                               j % jo % metar_p, iv % metar_ef_p, &
                               j % jo % metar_q, iv % metar_ef_q    
      end if

      if (iv%info(ships)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    ships obs, Jo(actual)  = ', &
                               iv%info(ships)%ntotal, iv%nstats(ships), &
                               j % jo % ships_u, iv % ships_ef_u, &
                               j % jo % ships_v, iv % ships_ef_v, &
                               j % jo % ships_t, iv % ships_ef_t, &
                               j % jo % ships_p, iv % ships_ef_p, &
                               j % jo % ships_q, iv % ships_ef_q                                
      end if


      if (iv%info(geoamv)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    geoamv ob, Jo(actual)  = ', &
                               iv%info(geoamv)%ntotal, iv%nstats(geoamv), &
                               j % jo % geoamv_u, iv % geoamv_ef_u, &
                               j % jo % geoamv_v, iv % geoamv_ef_v, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(polaramv)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    polaramv,  Jo(actual)  = ', &
                               iv%info(polaramv)%ntotal, iv%nstats(polaramv), &
                               j % jo % polaramv_u, iv % polaramv_ef_u, &
                               j % jo % polaramv_v, iv % polaramv_ef_v, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if


      if (iv%info(gpspw)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    gpspw obs, Jo(actual)  = ', &
                               iv%info(gpspw)%ntotal, iv%nstats(gpspw), &
                               j % jo % gpspw_tpw, iv % gpspw_ef_tpw, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(gpsref)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    gpsref obs, Jo(actual)  = ', &
                               iv%info(gpsref)%ntotal, iv%nstats(gpsref), &
                               j % jo % gpsref_ref, iv % gpsref_ef_ref, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(sound)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    sound obs, Jo(actual)  = ', &
                               iv%info(sound)%ntotal, iv%nstats(sound), &
                               j % jo % sound_u, iv % sound_ef_u, &
                               j % jo % sound_v, iv % sound_ef_v, &
                               j % jo % sound_t, iv % sound_ef_t, &
                               j % jo % sound_q, iv % sound_ef_q, 0.0, 1.0 
      end if

      if (iv%info(sonde_sfc)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'sonde sfc obs, Jo(actual)  = ', &
                               iv%info(sonde_sfc)%ntotal, iv%nstats(sonde_sfc), &
                               j % jo % sonde_sfc_u, iv % synop_ef_u, &
                               j % jo % sonde_sfc_v, iv % synop_ef_v, &
                               j % jo % sonde_sfc_t, iv % synop_ef_t, &
                               j % jo % sonde_sfc_p, iv % synop_ef_p, &
                               j % jo % sonde_sfc_q, iv % synop_ef_q
      end if

      if (iv%info(mtgirs)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'   mtgirs obs, Jo(actual)  = ', &
                               iv%info(mtgirs)%ntotal, iv%nstats(mtgirs), &
                               j % jo % mtgirs_u, iv % mtgirs_ef_u, &
                               j % jo % mtgirs_v, iv % mtgirs_ef_v, &
                               j % jo % mtgirs_t, iv % mtgirs_ef_t, &
                               j % jo % mtgirs_q, iv % mtgirs_ef_q, 0.0, 1.0
      end if

      if (iv%info(tamdar)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'   tamdar obs, Jo(actual)  = ', &
                               iv%info(tamdar)%ntotal, iv%nstats(tamdar), &
                               j % jo % tamdar_u, iv % tamdar_ef_u, &
                               j % jo % tamdar_v, iv % tamdar_ef_v, &
                               j % jo % tamdar_t, iv % tamdar_ef_t, &
                               j % jo % tamdar_q, iv % tamdar_ef_q, 0.0, 1.0
      end if

      if (iv%info(tamdar_sfc)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'tamdar sfc obs,Jo(actual)  = ', &
                               iv%info(tamdar_sfc)%ntotal, iv%nstats(tamdar_sfc), &
                               j % jo % tamdar_sfc_u, iv % tamdar_sfc_ef_u, &
                               j % jo % tamdar_sfc_v, iv % tamdar_sfc_ef_v, &
                               j % jo % tamdar_sfc_t, iv % tamdar_sfc_ef_t, &
                               j % jo % tamdar_sfc_p, iv % tamdar_sfc_ef_p, &
                               j % jo % tamdar_sfc_q, iv % tamdar_sfc_ef_q
      end if

      if (iv%info(airep)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    airep obs, Jo(actual)  = ', &
                               iv%info(airep)%ntotal, iv%nstats(airep), &
                               j % jo % airep_u, iv % airep_ef_u, &
                               j % jo % airep_v, iv % airep_ef_v, &
                               j % jo % airep_t, iv % airep_ef_t, &
                               j % jo % airep_q, iv % airep_ef_q, &
                               0.0, 1.0
      end if

      if (iv%info(bogus)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    bogus obs, Jo(actual)  = ', &
                               iv%info(bogus)%ntotal, iv%nstats(bogus), &
                               j % jo % bogus_u, iv % bogus_ef_u, &
                               j % jo % bogus_v, iv % bogus_ef_v, &
                               j % jo % bogus_t, iv % bogus_ef_t, &
                               j % jo % bogus_q, iv % bogus_ef_q, &
                               j % jo % bogus_slp, iv % bogus_ef_slp
      end if

      if (iv%info(pilot)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    pilot obs, Jo(actual)  = ', &
                               iv%info(pilot)%ntotal, iv%nstats(pilot), &
                               j % jo % pilot_u, iv % pilot_ef_u, &
                               j % jo % pilot_v, iv % pilot_ef_v, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(ssmi_rv)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    ssmir obs, Jo(actual) = ', &
                                  iv%info(ssmi_rv)%ntotal, iv%nstats(ssmi_rv), &
                                  j % jo % ssmir_speed, iv % ssmir_ef_speed, &
                                  j % jo % ssmir_tpw, iv % ssmir_ef_tpw, &
                                  0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(satem)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    satem obs, Jo(actual)  = ', &
                               iv%info(satem)%ntotal, iv%nstats(satem), &
                               j % jo % satem_thickness, iv % satem_ef_thickness, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(ssmt1)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    ssmt1 obs, Jo(actual)  = ', &
                               iv%info(ssmt1)%ntotal, iv%nstats(ssmt1), &
                               j % jo % ssmt1_t, iv % ssmt1_ef_t, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(ssmt2)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    ssmt2 obs, Jo(actual)  = ', &
                               iv%info(ssmt2)%ntotal, iv%nstats(ssmt2), &
                               j % jo % ssmt2_rh, iv % ssmt2_ef_rh, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(qscat)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    qscat obs, Jo(actual)  = ', &
                               iv%info(qscat)%ntotal, iv%nstats(qscat), &
                               j % jo % qscat_u, iv % qscat_ef_u, &
                               j % jo % qscat_v, iv % qscat_ef_v, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(buoy)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    buoy  obs, Jo(actual)  = ', &
                               iv%info(buoy)%ntotal, iv%nstats(buoy), &
                               j % jo % buoy_u, iv % buoy_ef_u, &
                               j % jo % buoy_v, iv % buoy_ef_v, &
                               j % jo % buoy_t, iv % buoy_ef_t, &
                               j % jo % buoy_p, iv % buoy_ef_p, &
                               j % jo % buoy_q, iv % buoy_ef_q                                
      end if

      if (iv%info(rain)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')' rainfall  obs, Jo(actual)  = ', &
                               iv%info(rain)%ntotal, iv%nstats(rain), &
                               j % jo % rain_r, iv % rain_ef_r, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if

      if (iv%info(profiler)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    profiler,  Jo(actual)  = ', &
                               iv%info(profiler)%ntotal, iv%nstats(profiler), &
                               j % jo % profiler_u, iv % profiler_ef_u, &
                               j % jo % profiler_v, iv % profiler_ef_v, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if
      if (iv%info(airsr)%ntotal > 0) then
         write(unit=jo_unit,fmt='(a30,2i8,10f15.5)')'    airsr obs, Jo(actual)  = ', &
                               iv%info(airsr)%ntotal, iv%nstats(airsr), &
                               j % jo % airsr_t, iv % airsr_ef_t, &
                               j % jo % airsr_q, iv % airsr_ef_q, &
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0
      end if
      do i = 1, iv%num_inst                 ! loop for sensor
         do k = 1, iv%instid(i)%nchan
            if (j % jo % rad(i) % num_ichan(k) > 0) then
               write(unit=jo_unit,fmt='(a30,a16,i5,i10,8f15.5)')'    radiance,  Jo(actual)  = ', &
                  iv%instid(i)%rttovid_string, iv%instid(i)%ichan(k) , &
                  j % jo % rad(i) % num_ichan(k), &
                  j % jo % rad(i) % jo_ichan(k)
            end if
         end do
      end do
   end if

   if (trace_use) call da_trace_exit("da_get_var_diagnostics")
      
end subroutine da_get_var_diagnostics


subroutine da_get_innov_vector (it, num_qcstat_conv, ob, iv, grid, config_flags)


   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: it
   integer,                    intent(inout) :: num_qcstat_conv(:,:,:,:)
   type(y_type),               intent(inout) :: ob ! Observations.
   type(iv_type),              intent(inout) :: iv ! Innovation vector(O-B).
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   type(xbx_type)     :: xbx          ! Header & non-gridded vars.

   character(len=120) :: filename, filename1

   integer            :: n, ios, time_step_seconds, ierr
   character*256 :: timestr, timestr1

   real, dimension(:,:,:), allocatable :: hr_rainc, hr_rainnc
   real, dimension(:,:),   allocatable :: savegridrainc, savegridrainnc
   integer            :: fgat_rain

   if (trace_use) call da_trace_entry("da_get_innov_vector") 

   call da_message((/"Calculate innovation vector(iv)"/))
   call da_get_unit(qcstat_conv_unit)

    write(unit=filename, fmt='(a,i2.2,a,i4.4)') 'rej_obs_conv_',it,'.', myproc

   open (unit=qcstat_conv_unit,file=trim(filename),form='formatted',status='replace', &
      iostat=ios)
   if (ios /= 0) then
      call da_error("da_get_innov_vector.inc",44, &
         (/"Cannot open qc observation file"//filename/))
   end if


   iv%ptop = grid%xb%ptop

   filename = ' '
   
   if ( use_rainobs .and. num_fgat_time > 1 .and. var4d) then
      allocate (hr_rainc (ims:ime,jms:jme,1:num_fgat_time))
      allocate (hr_rainnc(ims:ime,jms:jme,1:num_fgat_time))
      hr_rainc =0.0
      hr_rainnc=0.0
      allocate (savegridrainc (ims:ime,jms:jme))
      allocate (savegridrainnc (ims:ime,jms:jme))
      savegridrainc  =0.0
      savegridrainnc =0.0

      call domain_clock_get (grid, stop_timestr=timestr1)
      call domain_clock_set( grid, current_timestr=timestr1 )
      call domain_clock_set (grid, time_step_seconds=-1*var4d_bin)
      call domain_clockprint(150, grid, 'get CurrTime from clock,')

      fgat_rain = num_fgat_time
      do n= num_fgat_time , 1, -1
         call domain_clock_get( grid, current_timestr=timestr )
         call da_read_basicstates ( xbx, grid, config_flags, timestr ) 
         if ( fgat_rain_flags(n) ) then
            call da_get_hr_rain (fgat_rain, grid, hr_rainc, hr_rainnc, savegridrainc, savegridrainnc)
            fgat_rain = fgat_rain - 1
         end if
         if (n > 1) call domain_clockadvance (grid)
         call domain_clockprint(150, grid, 'get CurrTime from clock,')
      enddo
   endif

   if (num_fgat_time > 1) then
      call domain_clock_get (grid, stop_timestr=timestr1)
      call domain_clock_set( grid, current_timestr=timestr1 )
      call domain_clock_set (grid, time_step_seconds=-1*var4d_bin)
      call domain_clockprint(150, grid, 'get CurrTime from clock,')
   endif

   do n= num_fgat_time , 1, -1
      iv%time = n
      iv%info(:)%n1 = iv%info(:)%plocal(iv%time-1) + 1
      iv%info(:)%n2 = iv%info(:)%plocal(iv%time)

      if (num_fgat_time > 1) then
         if (var4d) then
            call domain_clock_get( grid, current_timestr=timestr )
            call da_read_basicstates ( xbx, grid, config_flags, timestr ) 
         else
            write(unit=filename1, fmt='(a,i2.2)') 'fg',n
            write(unit=stdout,fmt='(A,A)') 'Reading first guess ',trim(filename1)
            call da_read_basicstates ( xbx, grid, config_flags, timestr, filename1 ) 
         endif 
      end if

      ! Radiosonde:
      if (iv%info(sound)%nlocal > 0) then
         call da_get_innov_vector_sound     (it, num_qcstat_conv, grid, ob, iv)
         call da_get_innov_vector_sonde_sfc (it, num_qcstat_conv, grid, ob, iv)
      end if
      if (iv%info(mtgirs)%nlocal    > 0) &
         call da_get_innov_vector_mtgirs    (it, num_qcstat_conv, grid, ob, iv)
      if (iv%info(tamdar)%nlocal    > 0) &
         call da_get_innov_vector_tamdar    (it, num_qcstat_conv, grid, ob, iv)
      if (iv%info(tamdar_sfc)%nlocal> 0) &
         call da_get_innov_vector_tamdar_sfc(it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(synop)%nlocal     > 0) &
         call da_get_innov_vector_synop     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(geoamv)%nlocal    > 0) &
         call da_get_innov_vector_geoamv    (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(polaramv)%nlocal  > 0) &
         call da_get_innov_vector_polaramv  (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(airep)%nlocal     > 0) &
         call da_get_innov_vector_airep     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(pilot)%nlocal     > 0) &
         call da_get_innov_vector_pilot     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(bogus)%nlocal     > 0) &
         call da_get_innov_vector_bogus     (it, num_qcstat_conv, grid, ob, iv)
      if (iv%info(metar)%nlocal     > 0) &
         call da_get_innov_vector_metar     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(ships)%nlocal     > 0) &
         call da_get_innov_vector_ships     (it, num_qcstat_conv,grid, ob, iv)

      if ( (use_gpspwObs .and. iv%info(gpspw)%nlocal > 0) .or. &
           (pseudo_var(1:3) == 'tpw' .and. iv%info(pseudo)%nlocal > 0) ) then
         call da_get_innov_vector_gpspw     (it, num_qcstat_conv,grid, ob, iv)
      else if ( (use_gpsztdObs .and. iv%info(gpspw)%nlocal > 0) .or. &
           (pseudo_var(1:3) == 'ztd' .and. iv%info(pseudo)%nlocal > 0) ) then
         call da_get_innov_vector_gpsztd    ( it, num_qcstat_conv, grid, ob, iv )
      endif 
      if ( (iv%info(gpsref)%nlocal  > 0  ) .or. &
           (pseudo_var(1:3) == 'ref' .and. iv%info(pseudo)%nlocal > 0) ) &
         call da_get_innov_vector_gpsref    (it, num_qcstat_conv, grid, ob, iv)
      if (iv%info(ssmi_tb)%nlocal   > 0) &
         call da_get_innov_vector_ssmi_tb   (it, grid, ob, iv)
      if (iv%info(ssmi_rv)%nlocal   > 0) &
         call da_get_innov_vector_ssmi_rv   (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(ssmt1)%nlocal     > 0) &
         call da_get_innov_vector_ssmt1     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(ssmt2)%nlocal     > 0) &
         call da_get_innov_vector_ssmt2     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(satem)%nlocal     > 0) &
         call da_get_innov_vector_satem     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(radar)%nlocal     > 0) &
         call da_get_innov_vector_radar     (it, grid, ob, iv)
      if (iv%info(qscat)%nlocal     > 0) &
         call da_get_innov_vector_qscat     (it, num_qcstat_conv,grid, ob, iv)
      if (iv%info(profiler)%nlocal  > 0) &
         call da_get_innov_vector_profiler  (it,num_qcstat_conv, grid, ob, iv)
      if (iv%info(buoy)%nlocal      > 0) &
         call da_get_innov_vector_buoy      (it,num_qcstat_conv, grid, ob, iv)
      if (iv%info(rain)%nlocal      > 0) &
         call da_get_innov_vector_rain      (it,num_qcstat_conv, grid, ob, iv, &
                                             hr_rainc(:,:,n),hr_rainnc(:,:,n))

      if (use_rad                      ) call da_get_innov_vector_radiance (it, grid, ob, iv)
      if (iv%info(pseudo)%nlocal    ==1) call da_get_innov_vector_pseudo   (grid, ob, iv)
      if (iv%info(airsr)%nlocal     > 0) &
      call da_get_innov_vector_airsr    (it,num_qcstat_conv, grid, ob, iv)

   !----------------------------------------------
   ! [5]  write out iv in ascii format
   !-----------------------------------------------

      if ( multi_inc == 1 ) then

          call da_write_iv_for_multi_inc(n, iv)

      elseif ( multi_inc == 2 ) then

          call da_read_iv_for_multi_inc(n, iv)

      endif

      if (n > 1 .and. var4d) call domain_clockadvance (grid)
      call domain_clockprint(150, grid, 'DEBUG Adjoint Forcing:  get CurrTime from clock,')

   end do

   if (use_rad) then
      if ( use_varbc .or. freeze_varbc ) then
         if ( num_fgat_time > 1 ) call da_varbc_coldstart(iv)
      end if
      if ( use_varbc .and. it == 1 ) call da_varbc_precond(iv)
   end if

   if ( use_rainobs .and. num_fgat_time > 1 .and. var4d) then
      deallocate (hr_rainc )
      deallocate (hr_rainnc)
      deallocate (savegridrainc)
      deallocate (savegridrainnc)
   endif

   ! For FGAT , we need to restore the the field of the analysis time at this point
   ! At this point, the current time is the beginning of the time window
   if ( num_fgat_time > 1 .and. .not. var4d) then
      ! Set the current time to the analysis_date
      call domain_clock_set( grid, current_timestr=trim(analysis_date(1:19)) )
      ! set the basic states to the analysis time.
      write(unit=filename1, fmt='(a)') 'fg'
      write(unit=stdout,fmt='(A,A)') 'Restore to first guess :fg at ',trim(analysis_date(1:19))
      call da_read_basicstates ( xbx, grid, config_flags, timestr, filename1)
   end if

   if (num_fgat_time > 1) then
      call nl_get_time_step ( grid%id, time_step_seconds)
      call domain_clock_set (grid, time_step_seconds=time_step_seconds)
      call domain_clockprint(150, grid, 'get CurrTime from clock,')
   end if

   if ( multi_inc == 1 ) then
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if ( myproc == 0 )  call da_join_iv_for_multi_inc()
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       call wrf_message("*** WRF-Var multi-increment stage 1 completed successfully ***")
       !call wrfu_finalize
       !call wrf_shutdown
       !stop
       close(qcstat_conv_unit)
       call da_free_unit(qcstat_conv_unit)
       if (trace_use) call da_trace_exit("da_get_innov_vector")
       return
   endif

   !-----------------------------------------------------------------------
   ! [2] Having calculated the real O-Bs, optionally overwrite with scaled,
   !    random values:
   !----------------------------------------------------------------------- 
   
   if (omb_set_rand) call da_random_omb_all( iv, ob)
   
   !------------------------------------------------------------------------  
   ! [3] Optionally rescale observation errors:
   !------------------------------------------------------------------------ 
   
   if (use_obs_errfac) call da_use_obs_errfac( iv)

   !------------------------------------------------------------------------  
   ! [4] Optionally add Gaussian noise to O, O-B:
   !------------------------------------------------------------------------ 

   if (omb_add_noise) then
      call da_add_noise_to_ob( iv, ob)
   !#ifdef 1
   !      if ((num_procs > 1) .and.(.not. use_rad)) call da_write_noise_to_ob(iv)
   !      if ((.not. use_rad)) call da_write_noise_to_ob(iv)
      call da_write_noise_to_ob(iv)
   !#endif
   end if

   !----------------------------------------------
   ! [5]  write out radiance iv in ascii format
   !-----------------------------------------------
   if (use_rad .and. write_iv_rad_ascii) then
      call da_write_iv_rad_ascii(it,ob,iv)
   end if

   !----------------------------------------------------------
   ! [6]  write out filtered radiance obs in binary format
   !----------------------------------------------------------

   if (write_filtered_rad) then
      write(unit=stdout,fmt='(A)') 'Writing filtered radiance'
      call da_write_filtered_rad(ob,iv)
   end if

   close(qcstat_conv_unit)
   call da_free_unit(qcstat_conv_unit)

   if (trace_use) call da_trace_exit("da_get_innov_vector")

end subroutine da_get_innov_vector


real function da_dot(n,x,y)

   !-----------------------------------------------------------------------
   ! Purpose: forms the dot product of two vectors.
   ! uses unrolled loops for increments equal to one.
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in) :: n
   real,    intent(in) :: x(n)
   real,    intent(in) :: y(n)

   real    :: dtemp1
   integer :: i,m,mp1

   da_dot = 0.0
   if (n <= 0) return

   if (trace_use) call da_trace_entry("da_dot")    

   dtemp1 = 0.0

   ! code for both increments equal to 1

   if (n > 0) then
      m = mod(n,5)
      if (m /= 0) then
         do i = 1,m
            dtemp1 = dtemp1 + x(i)*y(i)
         end do
      end if
      if (n >= 5) then
         mp1 = m + 1
         do i = mp1,n,5
            dtemp1 = dtemp1 + x(i   )*y(i   ) + x(i + 1)*y(i + 1) + &
                              x(i + 2)*y(i + 2) + x(i + 3)*y(i + 3) + &
                              x(i + 4)*y(i + 4)
         end do
      end if
   end if

   da_dot = dtemp1

   if (trace_use) call da_trace_exit("da_dot")    


end function da_dot


real function da_dot_cv(cv_size, x, y, grid, mzs, jp_start, jp_end)

   !-----------------------------------------------------------------------
   ! Purpose: Forms the dot product of two vectors that are organized in the 
   ! format of a "cv_type".  
   !
   ! Capable of producing bitwise-exact results for distributed-memory 
   ! parallel runs for testing.  This feature is very slow and consumes 
   ! lots of memory. 
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in) :: cv_size           ! Size of array (tile).
   real,             intent(in) :: x(cv_size)        ! 1st vector.
   real,             intent(in) :: y(cv_size)        ! 1st vector.
   type(domain),     intent(in) :: grid              ! decomposed dimensions
   integer,          intent(in) :: mzs(7)            ! mz for each variable
                                                     ! (to identify 2D arrays)
   integer, optional,intent(in) :: jp_start, jp_end

   logical       :: lvarbc
   real, pointer :: xx(:), yy(:)                     ! Temporary vectors.
   real, pointer :: xg(:), yg(:)                     ! Serial data arrays.
   real          :: dtemp1(1), dtmp                  ! Temporary.

   if (trace_use) call da_trace_entry("da_dot_cv")

   allocate(xx(1:cv_size))
   allocate(yy(1:cv_size))
   xx = x
   yy = y

 ! VarBC parameters are global (no summation across processors)
 !-------------------------------------------------------------
   lvarbc = present(jp_start) .and. present(jp_end)
   if (lvarbc) lvarbc = lvarbc .and. (jp_end >= jp_start)   
   if (lvarbc .and. .not. rootproc) then
      xx(jp_start:jp_end) = 0.
      yy(jp_start:jp_end) = 0.
   end if
      
   ! Bitwise-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not
   ! testing.  This will always be .false. for a serial run or 
   ! one-processor 1 run.

   if (test_dm_exact) then

      ! Collect local cv arrays x and y to globally-sized serially-ordered 
      ! arrays xg and yg.  Note that xg and yg will only exist on the 
      ! monitor task.  

      if (rootproc) then
!         cv_size_domain = cv_size_domain_jb+cv_size_domain_je+cv_size_domain_jp+cv_size_domain_js   
         cv_size_domain = wrf_dm_sum_integer(cv_size)
         allocate(xg(1:cv_size_domain))
         allocate(yg(1:cv_size_domain))
      end if

      call da_cv_to_global(cv_size, cv_size_domain, xx, grid, mzs, xg)
      call da_cv_to_global(cv_size, cv_size_domain, yy, grid, mzs, yg)

      if (rootproc) then
         dtemp1(1) = da_dot(cv_size_domain, xg, yg)
         deallocate(xg, yg)
      end if

      ! Broadcast result from monitor to other tasks.  
      call wrf_dm_bcast_real(dtemp1, 1)

   else

      dtemp1(1) = da_dot(cv_size, xx, yy)
      
      !if (.not. global) then
         dtmp = dtemp1(1)
         ! summation across processors:
         dtemp1(1) = wrf_dm_sum_real(dtmp)
      !end if

   end if
  
   deallocate(xx,yy)

   da_dot_cv = dtemp1(1)

   if (trace_use) call da_trace_exit("da_dot_cv")

end function da_dot_cv
subroutine da_write_diagnostics(it, grid,num_qcstat_conv, ob, iv, re, y, j)

   !---------------------------------------------------------------------------
   ! Purpose: Output data assmiilation diagnostics
   !---------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: it
   type (domain),  intent(in)    :: grid
   integer,        intent(inout) :: num_qcstat_conv(:,:,:,:)
   type (y_type),  intent(in)    :: ob      ! Observation structure.
   type (iv_type), intent(inout) :: iv      ! innovation vector.
   type (y_type),  intent(inout) :: re      ! residual vector.
   type (y_type),  intent(in)    :: y       ! y = H(x_inc) structure.
   type (j_type),  intent(inout) :: j       ! Cost function.

   character(len=filename_len)   :: filename

   if (trace_use) call da_trace_entry("da_write_diagnostics")

   if (rootproc .and. print_detail_outerloop) then
      write(filename,'(a,i2.2)') "statistics_",it
      call da_get_unit (stats_unit)
      open(unit=stats_unit,file=trim(filename),status="replace")
   endif

   iv%nstats(:) = 0

   !---------------------------------------------------------------------------
   ! [1.0] Calculate innovation vector (O-B) statistics:
   !---------------------------------------------------------------------------

   if (iv%info(sound)%ntotal    > 0) call da_oi_stats_sound    (stats_unit, iv, ob)
   if (iv%info(sonde_sfc)%ntotal    > 0) call da_oi_stats_sonde_sfc(stats_unit, iv, ob)
   if (iv%info(mtgirs)%ntotal   > 0) call da_oi_stats_mtgirs   (stats_unit, iv, ob)
   if (iv%info(tamdar)%ntotal   > 0) call da_oi_stats_tamdar   (stats_unit, iv, ob)
   if (iv%info(tamdar_sfc)%ntotal   > 0) call da_oi_stats_tamdar_sfc(stats_unit, iv, ob)
   if (iv%info(synop)%ntotal    > 0) call da_oi_stats_synop    (stats_unit, iv, ob)
   if (iv%info(geoamv)%ntotal   > 0) call da_oi_stats_geoamv   (stats_unit, iv, ob)
   if (iv%info(polaramv)%ntotal > 0) call da_oi_stats_polaramv (stats_unit, iv, ob)
   if (iv%info(airep)%ntotal    > 0) call da_oi_stats_airep    (stats_unit, iv, ob)
   if (iv%info(pilot)%ntotal    > 0) call da_oi_stats_pilot    (stats_unit, iv, ob)
   if (iv%info(metar)%ntotal    > 0) call da_oi_stats_metar    (stats_unit, iv, ob)
   if (iv%info(ships)%ntotal    > 0) call da_oi_stats_ships    (stats_unit, iv, ob)
   if (iv%info(gpspw)%ntotal    > 0) call da_oi_stats_gpspw    (stats_unit, iv)
   if (iv%info(gpsref)%ntotal   > 0) call da_oi_stats_gpsref   (stats_unit, iv)
   if (iv%info(ssmi_tb)%ntotal  > 0) call da_oi_stats_ssmi_tb  (stats_unit, iv)
   if (iv%info(ssmi_rv)%ntotal  > 0) call da_oi_stats_ssmi_rv  (stats_unit, iv)
   if (iv%info(satem)%ntotal    > 0) call da_oi_stats_satem    (stats_unit, iv)
   if (iv%info(ssmt1)%ntotal    > 0) call da_oi_stats_ssmt1    (stats_unit, iv)
   if (iv%info(ssmt2)%ntotal    > 0) call da_oi_stats_ssmt2    (stats_unit, iv)
   if (iv%info(qscat)%ntotal    > 0) call da_oi_stats_qscat    (stats_unit, iv, ob)
   if (iv%info(profiler)%ntotal > 0) call da_oi_stats_profiler (stats_unit, iv, ob)
   if (iv%info(buoy)%ntotal     > 0) call da_oi_stats_buoy     (stats_unit, iv, ob)
   if (iv%info(radar)%ntotal    > 0) call da_oi_stats_radar    (stats_unit, iv)
   if (iv%info(bogus)%ntotal    > 0) call da_oi_stats_bogus    (stats_unit, iv)
   if (iv%info(airsr)%ntotal    > 0) call da_oi_stats_airsr    (stats_unit, iv)
   if (iv%info(rain)%ntotal     > 0) call da_oi_stats_rain     (stats_unit, iv)

   if (iv%num_inst > 0) call da_oi_stats_rad (stats_unit, iv)

   if (num_pseudo > 0) call da_oi_stats_pseudo (stats_unit, iv)

   !---- ----------------------------------------------------------------------
   ! [2.0] Calculate residual vector (O-A) statistics:
   !---------------------------------------------------------------------------
if (.not. anal_type_verify) then

   if (iv%info(sound)%ntotal    > 0) call da_ao_stats_sound    (stats_unit, iv, re, ob)
   if (iv%info(sonde_sfc)%ntotal    > 0) call da_ao_stats_sonde_sfc (stats_unit, iv, re, ob)
   if (iv%info(mtgirs)%ntotal   > 0) call da_ao_stats_mtgirs   (stats_unit, iv, re, ob)
   if (iv%info(tamdar)%ntotal   > 0) call da_ao_stats_tamdar   (stats_unit, iv, re, ob)
   if (iv%info(tamdar_sfc)%ntotal   > 0) call da_ao_stats_tamdar_sfc(stats_unit, iv, re, ob)
   if (iv%info(synop)%ntotal    > 0) call da_ao_stats_synop    (stats_unit, iv, re, ob)
   if (iv%info(geoamv)%ntotal   > 0) call da_ao_stats_geoamv   (stats_unit, iv, re, ob)
   if (iv%info(polaramv)%ntotal > 0) call da_ao_stats_polaramv (stats_unit, iv, re, ob)
   if (iv%info(airep)%ntotal    > 0) call da_ao_stats_airep    (stats_unit, iv, re, ob)
   if (iv%info(pilot)%ntotal    > 0) call da_ao_stats_pilot    (stats_unit, iv, re, ob)
   if (iv%info(metar)%ntotal    > 0) call da_ao_stats_metar    (stats_unit, iv, re, ob)
   if (iv%info(ships)%ntotal    > 0) call da_ao_stats_ships    (stats_unit, iv, re, ob)
   if (iv%info(gpspw)%ntotal    > 0) call da_ao_stats_gpspw    (stats_unit, iv, re)
   if (iv%info(gpsref)%ntotal   > 0) call da_ao_stats_gpsref   (stats_unit, iv, re)
   if (iv%info(ssmi_tb)%ntotal  > 0) call da_ao_stats_ssmi_tb  (stats_unit, iv, re)
   if (iv%info(ssmi_rv)%ntotal  > 0) call da_ao_stats_ssmi_rv  (stats_unit, iv, re)
   if (iv%info(satem)%ntotal    > 0) call da_ao_stats_satem    (stats_unit, iv, re)
   if (iv%info(ssmt1)%ntotal    > 0) call da_ao_stats_ssmt1    (stats_unit, iv, re)
   if (iv%info(ssmt2)%ntotal    > 0) call da_ao_stats_ssmt2    (stats_unit, iv, re)
   if (iv%info(qscat)%ntotal    > 0) call da_ao_stats_qscat    (stats_unit, iv, re, ob)
   if (iv%info(profiler)%ntotal > 0) call da_ao_stats_profiler (stats_unit, iv, re, ob)
   if (iv%info(buoy)%ntotal     > 0) call da_ao_stats_buoy     (stats_unit, iv, re, ob)
   if (iv%info(radar)%ntotal    > 0) call da_ao_stats_radar    (stats_unit, iv, re)
   if (iv%info(bogus)%ntotal    > 0) call da_ao_stats_bogus    (stats_unit, iv, re)
   if (iv%info(airsr)%ntotal    > 0) call da_ao_stats_airsr    (stats_unit, iv, re)
   if (iv%info(rain)%ntotal     > 0) call da_ao_stats_rain     (stats_unit, iv, re)  

   if (iv%num_inst > 0) call da_ao_stats_rad (stats_unit, iv, re)

   if (num_pseudo > 0) call da_ao_stats_pseudo (stats_unit, iv, re)

   !---------------------------------------------------------------------------
   ! [3.0] Calculate analysis increment (A-B) statistics:
   !---------------------------------------------------------------------------

   call da_analysis_stats (grid, stats_unit)

   !---------------------------------------------------------------------------
   ! [4.0] Write VAR diagnostic :
   !------------------------------------------------------------------------------

   call da_get_var_diagnostics (it, iv, j)

end if

   !---- -------------------------------------------------------------------------
   !  [5.0] Write observation data (O, O-B, O-A, y=hx'):
   !------------------------------------------------------------------------------

   call da_write_obs(it, ob, iv, re)

   ! Write ETKF observation files if required (note - 1PE only at present):
   if (anal_type_verify) then
      call da_write_obs_etkf(ob, iv, re)
   end if

   call da_final_write_obs(it, iv)

if (.not. anal_type_verify) then
   call da_write_y(iv, y)

   call da_final_write_y(iv)
   call da_print_qcstat(it, iv, num_qcstat_conv)
end if

   if (rootproc .and. print_detail_outerloop) then
       close(stats_unit)
       call da_free_unit (stats_unit)
   end if

   if (trace_use) call da_trace_exit("da_write_diagnostics")

end subroutine da_write_diagnostics


subroutine da_minimise_cg(grid, config_flags,            &
                           it, cv_size, xbx, be, iv, &
                           j_grad_norm_target, xhat, cv, &
                           re, y, j_cost)

   !-------------------------------------------------------------------------
   ! Purpose:         Main Conjugate Gradient minimisation routine 
   !
   ! Here 
   !    cv   is updated in outer-loop.
   !    xhat is the control variable in inner-loop.
   !
   ! Called from da_solve
   !
   ! History: 12/12/08 - Split J and GradJ calculations (Tom Auligne)
   !          12/12/08 - Re-orthonormalization option   (Tom Auligne)
   !
   !          Sep. 2010 - Add Cloud control variables   (hongli Wang)
   !          09/06/12 - Allow for variable ntmax in each outerloop (Mike Kavulich)
   !-------------------------------------------------------------------------

   implicit none

   integer, intent(in)               :: it    ! external iteration.
   integer, intent(in)               :: cv_size          ! Total cv size
   type (xbx_type),intent(inout)     :: xbx   ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be    ! background error structure.
   type (iv_type), intent(inout)     :: iv    ! ob. increment vector.
   real, intent(inout)               :: j_grad_norm_target ! Target norm.
   real, intent(inout)               :: xhat(1:cv_size)  ! control variable (local).
   real, intent(inout)               :: cv(1:cv_size)    ! control variable (local).
   type (y_type), intent(inout)      :: re    ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y     ! y = H(x_inc) structure.

   type (j_type), intent(out)        :: j_cost                 ! cost function

   type(domain), intent(inout)       :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   integer                           :: iter            
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                           :: mz(7)
   real                              :: fhat(1:cv_size)        ! cv copy.
   real                              :: ghat(1:cv_size)        ! cv copy.
   real                              :: ghat0(1:cv_size)       ! cv copy.
   real                              :: phat(1:cv_size)        ! cv copy.
   type(qhat_type), allocatable      :: qhat(:)                ! cv copy.
   real                              :: apdotp,step,rrmold,rrmnew,ratio 
   real                              :: ob_grad, rrmnew_norm, gdot
   real                              :: j_total, j0_total
 
   ! Variables for Conjugate Gradient preconditioning
   real                              :: precon(1:cv_size)      ! cv copy.
   real                              :: g_total, g_partial, jo_partial                          
   integer                           :: i, ii, nv, nn, istart, iend
   integer                           ::  sz(5)

   if (trace_use) call da_trace_entry("da_minimise_cg")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using CG method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be % ne /)
   sz = (/ be%cv%size1, be%cv%size2, be%cv%size3, be%cv%size4, be%cv%size5 /)
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

   call da_calculate_j(it, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)

   j0_total = j_cost%total
   if (j0_total == 0.0) then
      if (trace_use) call da_trace_exit("da_minimise_cg")
      return
   end if
   ghat0 = ghat
   
   ! [1.1] Preconditioning:
   !-----------------------
   precon  = 1.0
   
   if (precondition_cg) then
      g_total = da_dot(cv_size,ghat,ghat)
      
      iend    = 0
      do nv = 1, 5
         nn = sz(nv) / mz(nv)
	 do ii = 1, mz(nv)
            istart     = iend + 1
            iend       = istart + nn - 1
	    g_partial  = da_dot(nn, ghat(istart:iend), ghat(istart:iend))
            jo_partial = j0_total / SUM(mz(1:5))
	    precon(istart:iend)=  1 / &
	       (1 + precondition_factor*(g_partial/g_total)/(jo_partial/j0_total)) 
	 end do
      end do
   end if
   
   phat  = - precon * ghat

   rrmold = da_dot_cv(cv_size, -phat, ghat, grid, mz, jp_start, jp_end)
   j_grad_norm_target = sqrt (rrmold)

   if (orthonorm_gradient) then
      allocate(qhat(0:ntmax(it)))
      allocate(qhat(0)%values(1:cv_size))
      qhat(0)%values = ghat / sqrt(rrmold)
   end if

   write(unit=stdout,fmt='("Starting outer iteration : ",i3)') it
   write(unit=stdout,fmt=11) j0_total, sqrt(rrmold), eps(it)*j_grad_norm_target
11 format('Starting cost function: ' ,1PD22.15,', Gradient= ',1PD22.15,/,&
          'For this outer iteration gradient target is:       ',1PD22.15)
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------------------'
   write(unit=stdout,fmt='(A)') &
      '              Loop Iter     Cost Function              Gradient                   Step'

   write(unit=stdout,fmt=12) " minimize_cg ", it, 0, j0_total, sqrt(rrmold), 0.0

   !-------------------------------------------------------------------------
   ! [2.0] iteratively solve for minimum of cost function:
   !-------------------------------------------------------------------------

   do iter=1, ntmax(it)
      if (rrmold == 0.0) exit

      call da_calculate_gradj(it,iter,cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                              be%cv%size_jl,xbx,be,iv,phat,y,fhat,grid,config_flags)				 
      
      apdotp = da_dot_cv(cv_size, fhat, phat, grid, mz, jp_start, jp_end)

      step = 0.0
      if (apdotp .gt. 0.0) step = rrmold/apdotp
      
      ghat = ghat + step * fhat
      xhat = xhat + step * phat
      
    ! Orthonormalize new gradient (using modified Gramm-Schmidt algorithm)
      if (orthonorm_gradient) then
         do i = iter-1, 0, -1
            gdot = da_dot_cv(cv_size, ghat, qhat(i)%values, grid, mz, jp_start, jp_end)
            ghat = ghat - gdot * qhat(i)%values
         end do
      end if
      
      rrmnew = da_dot_cv (cv_size, precon*ghat, ghat, grid, mz, jp_start, jp_end)
      rrmnew_norm = sqrt(rrmnew)

      ratio = 0.0
      if (rrmold .gt. 0.0) ratio = rrmnew/rrmold

      if (orthonorm_gradient) then
         allocate(qhat(iter)%values(1:cv_size))
         qhat(iter)%values = ghat / rrmnew_norm
      end if

      phat         = - precon * ghat       + ratio * phat

      rrmold=rrmnew

    ! Print Gradient (and Cost Function)
    !-----------------------------------

      if (calculate_cg_cost_fn) then
         call da_calculate_j(it, iter, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
	                     be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, fhat, grid, config_flags)
         j_total = j_cost%total
      else
         j_total = j0_total + 0.5 * da_dot_cv(cv_size,ghat0,xhat,grid,mz,jp_start,jp_end)
      endif

      write(unit=stdout,fmt=12) " minimize_cg ", it, iter, j_total, rrmnew_norm, step         	 
      if (rrmnew_norm  < eps(it) * j_grad_norm_target) exit

12    format(a13,1x,i3,1x,i4,5x,1PD22.15,5x,1PD22.15,5x,1PD22.15)

   end do

   !-------------------------------------------------------------------------
   ! End of the minimization of cost function
   !-------------------------------------------------------------------------
   iter = MIN(iter, ntmax(it))

   ! Free memory used for reorthonormalization
   !------------------------------------------
   if (orthonorm_gradient) then
      do i = iter-1, 0, -1
         if (allocated(qhat(i)%values)) deallocate(qhat(i)%values)
      end do
      deallocate(qhat)
   end if
   
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------------------'
   write(unit=stdout,fmt='(A)') " "
   write(unit=stdout, &
      fmt='("Inner iteration stopped after ",i4," iterations")') iter
   write(unit=stdout,fmt='(A)') " "

   call da_calculate_j(it, iter, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)

   rrmnew_norm = SQRT(da_dot_cv(cv_size,ghat,ghat,grid,mz,jp_start,jp_end))

    write(unit=stdout,fmt=15) iter, j_cost%total , rrmnew_norm
15  format('Final: ',I3,' iter, J=',1PD22.15,', g=',1PD22.15)
    write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------------------'

   if (trace_use) call da_trace_exit("da_minimise_cg")

end subroutine da_minimise_cg
subroutine da_minimise_lz(grid, config_flags,            &
                           it, cv_size, xbx, be, iv,     &
                           j_grad_norm_target, xhat, qhat, cv, &
                           re, y, j_cost, eignvec, eignval, neign)

   !-------------------------------------------------------------------------
   ! Purpose:         Main Lanczos minimisation routine 
   !
   ! Here 
   !    cv   is updated in outer-loop.
   !    xhat is the control variable in inner-loop.
   !
   ! Called from da_solve
   !
   ! History: 07/30/2008  Creation (Tom Auligne)
   !
   ! Reference: Golub and Van Loan 1996 (p493)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: it                            ! external iteration.
   integer, intent(in)               :: cv_size                       ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                           ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                            ! background error structure.
   type (iv_type), intent(inout)     :: iv                            ! ob. increment vector.
   real, intent(inout)               :: j_grad_norm_target            ! Target norm.
   real, intent(out)                 :: xhat(1:cv_size)               ! control variable (local).
   real, intent(out)                 :: qhat(1:cv_size, 0:ntmax(it))  ! cv copy.
   real, intent(inout)               :: cv(1:cv_size)                 ! control variable (local).
   type (y_type), intent(inout)      :: re                            ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                             ! y = H(x_inc) structure.
   type (j_type), intent(out)        :: j_cost                        ! cost function
   real*8, intent(out)               :: eignvec(ntmax(it), ntmax(it))
   real*8, intent(out)               :: eignval(ntmax(it))
   integer, intent(out)              :: neign
   
   integer                           :: i, j, k, n, inst, iter
   integer                           :: mz(7), info, nsstwrk
   integer                           :: npred, ipred
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   real                              :: fhat(1:cv_size)        ! cv copy.
   real                              :: ghat(1:cv_size)        ! cv copy.
   real                              :: ghat0(1:cv_size)       ! cv copy.
   real                              :: shat(1:cv_size)        ! cv copy.
   real                              :: gdot, rho, mu, bndlm
   real                              :: alpha(ntmax(it)), beta(0:ntmax(it))
   real*8                            :: subdiag(ntmax(it))
   real                              :: c(cv_size), d(ntmax(it))
   real*8                            :: sstwrk(2*ntmax(it)-2)
   real                              :: bnds(ntmax(it))
   real                              :: adj_rhs, adj_sum_rhs ! < cv, cv >
   real                              :: j_total, j0_total
   real                              :: dfs, dfsmax
   
   if (trace_use) call da_trace_entry("da_minimise_lz")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using Lanczos method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz,be%alpha%mz, be % ne /)
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
   
   call da_calculate_j(it, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat0, grid, config_flags)
 
   j0_total = j_cost%total
   if (j0_total == 0.0) then
      if (trace_use) call da_trace_exit("da_minimise_lz")
      return
   end if
   ghat      = - ghat0
   beta(0)   = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
   qhat(:,0) = 0.0
   j_grad_norm_target = beta(0)

   write(unit=stdout,fmt='("Starting outer iteration : ",i3)') it
   write(unit=stdout,fmt=11) j0_total, beta(0),eps(it)
11 format('Starting cost function: ' ,1PD15.8,', Gradient: ',1PD15.8,/,&
          'For this outer iteration, the targeted reduction is:       ',1PD15.8)
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------'
   write(unit=stdout,fmt='(A)') 'Iter     Cost Function       Info Content'

   !-------------------------------------------------------------------------
   ! [2.0] Iteratively solve for minimum of cost function:
   !-------------------------------------------------------------------------
   do iter=1, ntmax(it)
      qhat(:,iter) = ghat / beta(iter-1)
      
      call da_calculate_gradj(it,iter,cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                              be%cv%size_jl, xbx,be,iv,qhat(:,iter),y,fhat,grid,config_flags)
       
    ! Apply Lanczos recurrence and orthonormalize new gradient (using modified Gramm-Schmidt)
    !----------------------------------------------------------------------------------------
      alpha(iter) = da_dot_cv(cv_size, qhat(:,iter), fhat, grid, mz, jp_start, jp_end)

      ghat        = fhat - alpha(iter)*qhat(:,iter) - beta(iter-1)*qhat(:,iter-1)
      do i = iter, 1, -1
         gdot = da_dot_cv(cv_size, ghat, qhat(:,i), grid, mz, jp_start, jp_end)
         ghat = ghat - gdot * qhat(:,i)
      end do

      beta(iter)  = SQRT(da_dot_cv (cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
      
    ! Lanczos iteration  
    !------------------
      if (iter == 1) then
         d(1) = alpha(1)
	 c    = qhat(:,1)
	 rho  = beta(0) / alpha(1)
	 xhat = rho * qhat(:,1)
         dfs  = da_dot_cv(cv_size, xhat, xhat, grid, mz, jp_start, jp_end)
         dfsmax = dfs
      else
         mu      = beta(iter-1) / d(iter-1)
	 d(iter) = alpha(iter) - beta(iter-1)*mu
	 c       = qhat(:,iter) - mu*c
	 rho     = - mu*d(iter-1)*rho / d(iter)
	 xhat    = xhat + rho*c
         dfs     = da_dot_cv(cv_size, rho*c, rho*c, grid, mz, jp_start, jp_end)
         if (dfs > dfsmax) dfsmax = dfs
      end if
      
    ! Determine eigenvalues and eigenvectors of the Lanczos tri-diagonal matrix
    !--------------------------------------------------------------------------
      eignval(1:iter)   = alpha(1:iter)
      subdiag(1:iter-1) = beta(1:iter-1)
      nsstwrk           = MAX(2*iter-2,1)
      info              = 0
      call DSTEQR('I',iter,eignval(1:iter),subdiag(1:iter-1),eignvec(:,1:iter),ntmax(it),&
	           sstwrk(1:nsstwrk),info)
      if (info /=0) write(stdout,*) 'Error in Lanczos minimization: SSTEQR returned ',info 
            
    ! Count converged eigenpairs (using Arnoldi relation)
    !----------------------------------------------------
      bndlm        = eps(it) * eignval(iter)  
      bnds(1:iter) = abs(beta(iter) * eignvec(iter,1:iter))
      neign        = COUNT(bnds(1:iter) <= bndlm)
               
    ! Print Cost Function and Infomation Content
    !-----------------------------------
      j_total = j0_total + 0.5 * da_dot_cv(cv_size,ghat0,xhat,grid,mz,jp_start,jp_end)
      write(unit=stdout,fmt=12)iter, j_total, dfs      	 
12    format(i3,5x,1PD15.8,5x,1PD15.8)

    ! Stop minimization based on Information Content
    ! ----------------------------------------------
      if (sqrt(dfs) < eps(it) * sqrt(dfsmax)) exit

    ! Option to calculate J Cost Function explicitly
    ! ----------------------------------------------
      if (calculate_cg_cost_fn) &
         call da_calculate_j(it, iter, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
	                     be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, fhat, grid, config_flags)
   end do

   !-------------------------------------------------------------------------
   ! End of the minimization of cost function
   !-------------------------------------------------------------------------
   iter = MIN(iter, ntmax(it))
   
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------'
   write(unit=stdout,fmt='(A)') " "
   write(unit=stdout, fmt='("Inner iteration stopped after ",i4," iterations")') iter
   write(unit=stdout,fmt='(A)') " "

   call da_calculate_j(it, iter, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)

   write(unit=stdout,fmt=15) iter, j_cost%total , &
                             SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
15  format('Final: ',I3,' iter, J=',1PD15.8,', g=',1PD15.8)
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------'

   write(stdout,*) 'Ritz eigenvalues: ',eignval(iter:1:-1)     
   write(stdout,*) 'Number of converged eigenpairs: ', neign
   neign = iter

   if (trace_use) call da_trace_exit("da_minimise_lz")

end subroutine da_minimise_lz


subroutine da_calculate_grady(iv, re, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Does part of the obs gradient operation   
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector (O-B).
   type (y_type),  intent(inout) :: re          ! Residual vector (O-A).
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   if (trace_use) call da_trace_entry("da_calculate_grady")

   !-------------------------------------------------------------------------
   ! [1.0] Compute components of Grad_y(Jo):
   !-------------------------------------------------------------------------

   if (iv%info(sound)%nlocal     > 0) call da_calculate_grady_sound    (iv, re, jo_grad_y)
   if (iv%info(sonde_sfc)%nlocal > 0) call da_calculate_grady_sonde_sfc(iv, re, jo_grad_y)
   if (iv%info(synop)%nlocal     > 0) call da_calculate_grady_synop    (iv, re, jo_grad_y)
   if (iv%info(geoamv)%nlocal    > 0) call da_calculate_grady_geoamv   (iv, re, jo_grad_y)
   if (iv%info(polaramv)%nlocal  > 0) call da_calculate_grady_polaramv (iv, re, jo_grad_y)
   if (iv%info(airep)%nlocal     > 0) call da_calculate_grady_airep    (iv, re, jo_grad_y)
   if (iv%info(pilot)%nlocal     > 0) call da_calculate_grady_pilot    (iv, re, jo_grad_y)
   if (iv%info(profiler)%nlocal  > 0) call da_calculate_grady_profiler (iv, re, jo_grad_y)
   if (iv%info(satem)%nlocal     > 0) call da_calculate_grady_satem    (iv, re, jo_grad_y)
   if (iv%info(metar)%nlocal     > 0) call da_calculate_grady_metar    (iv, re, jo_grad_y)
   if (iv%info(ships)%nlocal     > 0) call da_calculate_grady_ships    (iv, re, jo_grad_y)
   if (iv%info(buoy)%nlocal      > 0) call da_calculate_grady_buoy     (iv, re, jo_grad_y)
   if (iv%info(gpspw)%nlocal     > 0) call da_calculate_grady_gpspw    (iv, re, jo_grad_y)
   if (iv%info(gpsref)%nlocal    > 0) call da_calculate_grady_gpsref   (iv, re, jo_grad_y)
   if (iv%info(ssmi_tb)%nlocal   > 0) call da_calculate_grady_ssmi_tb  (iv, re, jo_grad_y) 
   if (iv%info(ssmi_rv)%nlocal   > 0) call da_calculate_grady_ssmi_rv  (iv, re, jo_grad_y) 
   if (iv%info(ssmt1)%nlocal     > 0) call da_calculate_grady_ssmt1    (iv, re, jo_grad_y)
   if (iv%info(ssmt2)%nlocal     > 0) call da_calculate_grady_ssmt2    (iv, re, jo_grad_y)
   if (iv%info(pseudo)%nlocal    > 0) call da_calculate_grady_pseudo   (iv, re, jo_grad_y)
   if (iv%info(bogus)%nlocal     > 0) call da_calculate_grady_bogus    (iv, re, jo_grad_y)  
   if (iv%info(qscat)%nlocal     > 0) call da_calculate_grady_qscat    (iv, re, jo_grad_y)
   if (iv%info(radar)%nlocal     > 0) call da_calculate_grady_radar    (iv, re, jo_grad_y)
   if (iv%info(mtgirs)%nlocal    > 0) call da_calculate_grady_mtgirs   (iv, re, jo_grad_y)
   if (iv%info(tamdar)%nlocal    > 0) call da_calculate_grady_tamdar   (iv, re, jo_grad_y)
   if (iv%info(tamdar_sfc)%nlocal> 0) call da_calculate_grady_tamdar_sfc(iv, re, jo_grad_y)
   if (iv%info(rain)%nlocal      > 0) call da_calculate_grady_rain     (iv, re, jo_grad_y)

   if (iv%num_inst               > 0) call da_calculate_grady_rad      (iv, re, jo_grad_y)
   if (iv%info(airsr)%nlocal     > 0) call da_calculate_grady_airsr    (iv, re, jo_grad_y)

   if (trace_use) call da_trace_exit("da_calculate_grady")

end subroutine da_calculate_grady


subroutine da_transform_vtoy(cv_size, be, ep, cv, iv, vp, vv, vp6, vv6, xbx, &
                              y, grid, config_flags)

   !----------------------------------------------------------------------
   ! Purpose:  Does transform of control variable (V) to Obs-space (Y)
   !----------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! Ensemble perturbation structure.
   real,                       intent(in)    :: cv(1:cv_size)     ! control variables.
   type(iv_type),              intent(inout) :: iv     ! innovation vector (o-b).
   type(vp_type),              intent(inout) :: vp     ! Grdipt/level CV.
   type(vp_type),              intent(inout) :: vp6     ! Grdipt/level CV for 6h.
   type(vp_type),              intent(inout) :: vv     ! Grdipt/EOF CV.
   type(vp_type),              intent(inout) :: vv6     ! Grdipt/EOF CV for 6h.
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(y_type),               intent(inout) :: y      ! y = H(x_inc).
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   type(x_type) :: shuffle
   integer :: nobwin, jl_start, jl_end, time_step_seconds

   character(len=4) :: filnam
   character(len=256) :: timestr, timestr1
   real, dimension(:,:,:), allocatable :: hr_rainc, hr_rainnc
   real, dimension(:,:,:), allocatable :: g_rainc, g_rainnc

   if (trace_use) call da_trace_entry("da_transform_vtoy")

   if (var4d) then
   else  ! not var4d
      call da_transform_vtox(grid, cv_size, xbx, be, ep, cv, vv, vp)
      iv%info(:)%n1 = 1
      iv%info(:)%n2 = iv%info(:)%nlocal
      if ( use_rad ) then
         iv%instid(:)%info%n1 = 1
         iv%instid(:)%info%n2 = iv%instid(:)%num_rad
      end if
      call da_transform_xtoxa(grid)
      call da_transform_xtoy(cv_size, cv, grid, iv, y)
   end if ! var4d

   !--------------------------------------------------------------
   ! TL of Variational Bias Correction
   !--------------------------------------------------------------
   if (use_varbc) call da_varbc_tl(cv_size, cv, iv, y)
   if (trace_use) call da_trace_exit("da_transform_vtoy")

end subroutine da_transform_vtoy


subroutine da_transform_vtoy_adj(cv_size, be, ep, cv, iv, vp, vv, vp6, vv6, xbx, y, &
   grid, config_flags, jcdf_flag)

   !-------------------------------------------------------------------------
   ! Purpose:  Does Adjoint of control variable (V) transform to Obs-space(Y)
   !-------------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! ensemble perturbation structure.
   real,                       intent(out)   :: cv(1:cv_size) ! control variables.
   type(iv_type),              intent(inout) :: iv     ! innovation vector (o-b).
   type(vp_type),              intent(inout) :: vp     ! Grdipt/level CV.
   type(vp_type),              intent(inout) :: vp6    ! Grdipt/level CV for 6h.
   type(vp_type),              intent(inout) :: vv     ! Grdipt/EOF CV.
   type(vp_type),              intent(inout) :: vv6    ! Grdipt/EOF CV for 6h.
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(y_type),               intent(inout) :: y      ! y = H(x_inc).
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   logical,                    intent(in)    :: jcdf_flag       ! additional flag to switch off JcDF, used to switch off JcDF in calculation of J.

   type(x_type) :: shuffle
   integer :: nobwin,ndynopt, jl_start, jl_end, time_step_seconds
   character*256 :: timestr, timestr1
   integer       :: i, j, k
   real          :: subarea, whole_area
   real          :: f_flag
   real, dimension(:,:,:), allocatable :: a_hr_rainc, a_hr_rainnc
   real, dimension(:,:,:), allocatable :: a_rainc, a_rainnc

   character(len=4) :: filnam

   call da_trace_entry("da_transform_vtoy_adj")

   cv = 0.0

   if (var4d) then

   else  ! not var4d
      iv%info(:)%n1 = 1
      iv%info(:)%n2 = iv%info(:)%nlocal
      if ( use_rad ) then
         iv%instid(:)%info%n1 = 1
         iv%instid(:)%info%n2 = iv%instid(:)%num_rad
      end if
      call da_zero_x(grid%xa)

      call da_transform_xtoy_adj(cv_size, cv, grid, iv, y,grid%xa)
      call da_transform_xtoxa_adj(grid)
      call da_transform_vtox_adj(grid, cv_size, xbx, be, ep, vp, vv, cv)
   end if ! var4d

   !--------------------------------------------------------------
   ! ADJ of Variational Bias Correction
   !--------------------------------------------------------------
   if (use_varbc) call da_varbc_adj(cv_size, cv, iv, y)
   if (trace_use) call da_trace_exit("da_transform_vtoy_adj")

end subroutine da_transform_vtoy_adj


subroutine da_transform_vtod_wpec(cv_size, be, ep, cv, vp, vv, xbx, grid)

   !----------------------------------------------------------------------
   ! Purpose:  Transform control variable (V) to observation-space (Y)
   !----------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! Ensemble perturbation structure.
   real,                       intent(in)    :: cv(1:cv_size)     ! control variables.
   type(vp_type),              intent(inout) :: vp     ! Grdipt/level CV.
   type(vp_type),              intent(inout) :: vv     ! Grdipt/EOF CV.
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(domain),               intent(inout) :: grid


   if (trace_use) call da_trace_entry("da_transform_vtod_wpec")

      call da_transform_vtox(grid, cv_size, xbx, be, ep, cv, vv, vp)
      call da_transform_xtoxa_all(grid)

      call da_wpec_constraint_lin(grid, xbx)


   if (trace_use) call da_trace_exit("da_transform_vtod_wpec")

end subroutine da_transform_vtod_wpec


subroutine da_transform_vtod_wpec_adj(cv_size, be, ep, cv, vp, vv, xbx, grid)

   !-------------------------------------------------------------------------
   ! Purpose:  Does Adjoint of control variable (V) transform to Obs-space(Y)
   !-------------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! ensemble perturbation structure.
   real,                       intent(out)   :: cv(1:cv_size) ! control variables.
   type(vp_type),              intent(inout) :: vp     ! Grdipt/level CV.
   type(vp_type),              intent(inout) :: vv     ! Grdipt/EOF CV.
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(domain),               intent(inout) :: grid

   call da_trace_entry("da_transform_vtod_wpec_adj")

   cv = 0.0

      call da_zero_x(grid%xa)

      call da_wpec_constraint_adj(grid, xbx)

      call da_transform_xtoxa_adj_all(grid)

      call da_transform_vtox_adj(grid, cv_size, xbx, be, ep, vp, vv, cv)

   if (trace_use) call da_trace_exit("da_transform_vtod_wpec_adj")

end subroutine da_transform_vtod_wpec_adj


subroutine da_adjoint_sensitivity(grid, config_flags, cv_size, xbx, be, iv, xhat, cv, y, shat)

   !-------------------------------------------------------------------------
   ! Purpose:        Read from WRF adjoint run or define derivative of
   !                 initial sensitivity for observation sensitivity studies
   !                 Shat = dF/dx (e.g. F = 0.5<x',B-1.x'> --> dF/dx=xhat)
   !
   ! Called from da_minimise_lz
   !
   ! History: 12/12/2008  Creation (Tom Auligne)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout)     :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer,        intent(in)        :: cv_size                ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                    ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                     ! background error structure.
   type (iv_type), intent(inout)     :: iv                     ! ob. increment vector.
   real,           intent(in)        :: xhat(1:cv_size)        ! control variable (local).
   real,           intent(in)        :: cv(1:cv_size)          ! control variable (local).
   type (y_type),  intent(inout)     :: y                      ! y = H(x_inc) structure.
   real,           intent(out)       :: shat(1:cv_size)        ! control variable (local).
   
   type (y_type)                     :: re
   type (j_type)                     :: j_cost                 ! cost function
   integer                           :: ix, jy, cvs, cve, var
   integer                           :: ibs, ibe, jbs, jbe
   integer                           :: i, j, mz(6)
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   real                              :: value

   if (trace_use) call da_trace_entry("da_adjoint_sensitivity")
   
   shat = 0.0

   select case(sensitivity_option) 
!--------------------------------------------------------
! Read Initial Sensitivity from WRF-Adjoint (NetCDF file)
!--------------------------------------------------------        
   case(0); 
         call da_zero_x (grid%xa)
         call da_transfer_xatowrftl_adj(grid, config_flags, 'gr01')
         call da_transform_xtoxa_adj(grid)
         call da_transform_vtox_adj (grid, cv_size, xbx, be, grid%ep, grid%vp, grid%vv, shat)

!---------------------------
! Define Initial Sensitivity
!---------------------------
   case(1); 
   ! Shat = dF/dx (e.g. F = 0.5<x',B-1.x'> --> dF/dx=xhat)
   !------------------------------------------------------
      shat = xhat

   case(2); 
    ! VarBC parameters 
    !-----------------
      jp_start = be % cv % size_jb + be % cv % size_je + 1
      jp_end   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
      if (jp_end >= jp_start) shat(jp_start:jp_end) = xhat(jp_start:jp_end)

   case(3); 
    ! Jo for all observations
    !------------------------
      call da_allocate_y (iv, re)
      call da_calculate_j(1,1,cv_size,0,0,0,0,xbx,be,iv,xhat,cv,re,y,j_cost,shat,grid,config_flags)
      call da_deallocate_y(re)
			    
   case default;
   end select
   if (trace_use) call da_trace_exit("da_adjoint_sensitivity")

end subroutine da_adjoint_sensitivity

subroutine da_sensitivity(grid, config_flags, it, cv_size, xbx, be, iv, xhat, qhat, cv, y, &
                          eignvec, eignval, neign)

   !-------------------------------------------------------------------------
   ! Purpose:        Compute observation sensitivity and impact
   !
   ! Called either from da_minimise_lz or da_solve
   !
   ! History: 03/05/2009  Creation (Tom Auligne)
   !          09/06/2012  Modified to allow variable ntmax in each outerloop (Mike Kavulich)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout)     :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: it                            ! external iteration.
   integer,        intent(in)        :: cv_size                       ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                           ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                            ! background error structure.
   type (iv_type), intent(inout)     :: iv                            ! ob. increment vector.
   real,           intent(in)        :: xhat(1:cv_size)               ! control variable (local).
   real,           intent(in)        :: qhat(1:cv_size, 0:ntmax(it))  
   real,           intent(in)        :: cv(1:cv_size)                 ! control variable (local).
   type (y_type),  intent(inout)     :: y                             ! y = H(x_inc) structure.
   real*8,         intent(in)        :: eignvec(ntmax(it), ntmax(it))
   real*8,         intent(in)        :: eignval(ntmax(it))
   integer,        intent(in)        :: neign                         ! Number of eigenpairs
   
   type (y_type)                     :: ktr         
   integer                           :: i, j, mz(7)
   integer                           :: jp_start, jp_end              ! Start/end indices of Jp.
   real                              :: shat(1:cv_size)               ! control variable (local).
   real                              :: amat(1:cv_size)               ! cv copy.
   real                              :: ritz(ntmax(it), ntmax(it))
    mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz,be%alpha%mz, be % ne /)
    jp_start = be % cv % size_jb + be % cv % size_je + 1
    jp_end   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

    ! Define Shat = dF/dx (e.g. F=1/2<x',B-1.x'> --> dF/dx=xhat)
    !-----------------------------------------------------------
    call da_adjoint_sensitivity(grid, config_flags, cv_size, xbx, be, iv, xhat, cv, y, shat)

    ! Apply Analysis Error Covariance Matrix estimation (A) to dF/dx 
    !---------------------------------------------------------------
      amat = 0.0
      do i = 1, neign
         do j = 1, neign
            ritz(i,j) = SUM( eignvec(i,1:neign) * (1.0/eignval(1:neign)) * eignvec(j,1:neign) )
	    amat      = amat + qhat(:,i) * ritz(i,j) * &
	                       da_dot_cv(cv_size, qhat(:,j), shat, grid, mz, jp_start, jp_end)
	 end do
      end do

    ! Calculate observation sensitivity: Kt = R-1.H.A
    !------------------------------------------------
      call da_allocate_y  (iv, ktr)

    ! Apply observation operator H 
      call da_transform_vtoy(cv_size, be, grid%ep, amat, iv, grid%vp, grid%vv, &
           grid%vp6, grid%vv6, xbx, ktr, grid, config_flags)
				
    ! Apply R-1 (for Observation Sensitivity) and then Dot Product with initial innovations (for Observation Impact)
      call da_obs_sensitivity(ktr, iv)
             
      call da_deallocate_y(ktr)

    ! Adjoint test    
!      write(stdout,*) 'ADJOINT_TEST2:', da_dot_cv(cv_size, xhat, xhat, grid, mz, jp_start, jp_end)

end subroutine da_sensitivity
subroutine da_amat_mul(be, grid, cv_size, ntmaxit, neign, eignval, eignvec, qhat, shat, xhat)

   !-------------------------------------------------------------------------
   ! Purpose:  Multiply a control vector by the Analysis Error Cov Matrix A 
   !
   ! Called from da_solve
   !
   ! History: 08/16/2010  Creation (Tom Auligne)
   !
   !-------------------------------------------------------------------------

   implicit none

   type (be_type), intent(in) :: be                     ! Background error structure.
   type (domain),  intent(in) :: grid
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: ntmaxit
   integer,        intent(in) :: neign
   real*8,         intent(in) :: eignvec(ntmaxit, ntmaxit)  
   real*8,         intent(in) :: eignval(ntmaxit)
   real,           intent(in) :: qhat(cv_size, 0:ntmaxit) ! Ritz vectors
   real,           intent(in) :: shat(cv_size)          ! Input vector to multiply by A
   real,           intent(out):: xhat(cv_size)          ! Output vector: xhat = A.shat

   integer                    :: mz(7)
   integer                    :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                    :: i, j
   real                       :: dot_cv

   if (trace_use) call da_trace_entry("da_amat_mul")

   mz       = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be % ne /)
   jp_start = be % cv % size_jb + be % cv % size_je + 1
   jp_end   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

   xhat     = 0.0
   do j = 1, neign
      dot_cv = da_dot_cv(cv_size, qhat(:,j), shat, grid, mz, jp_start, jp_end)
      do i = 1, neign
         xhat = xhat + qhat(:,i) * &
                SUM(eignvec(i,1:neign)*eignval(1:neign)*eignvec(j,1:neign)) * dot_cv
      end do
   end do
 
   if (trace_use) call da_trace_exit ("da_amat_mul")

end subroutine da_amat_mul
subroutine da_kmat_mul(grid, config_flags,            &
                        it, cv_size, xbx, be, iv,     &
                        xhat, qhat, cv, &
                        re, y, j_cost, eignvec, eignval, neign)

   !-------------------------------------------------------------------------
   ! Purpose:  Multiply the innovation vector by the Kalman Gain Matrix K
   !              using the precomputed Lanczos eigenpairs 
   !
   ! Called from da_solve
   !
   ! History: 05/04/2011  Creation (Tom Auligne)
   !          09/06/2012  Modified to allow variable ntmax size for each outerloop (Mike Kavulich)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: it                           ! external iteration.
   integer, intent(in)               :: cv_size                      ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                          ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                           ! background error structure.
   type (iv_type), intent(inout)     :: iv                           ! ob. increment vector.
   real, intent(out)                 :: xhat(1:cv_size)              ! Output vector: xhat=K.d
   real, intent(in)                  :: qhat(1:cv_size, 0:ntmax(it)) ! Ritz vectors
   real, intent(in)                  :: cv(1:cv_size)                ! control variable (local).
   type (y_type), intent(inout)      :: re                           ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                            ! y = H(x_inc) structure.
   type (j_type), intent(out)        :: j_cost                       ! cost function
   real*8, intent(in)                :: eignvec(ntmax(it), ntmax(it))
   real*8, intent(in)                :: eignval(ntmax(it))
   integer, intent(in)               :: neign

   real                              :: shat(1:cv_size)          ! cv copy.

   if (trace_use) call da_trace_entry("da_kmat_mul")

   write(*,*) 'Computing Analysis Increment: v = K.d = A.H^T.R-1.[y^o-H(x_b)]'

 ! Transfer [y^o-H(x_b)] information from iv(iv_type) into re(y_type)  
   call da_zero_y(iv,y)
   call da_calculate_residual(iv,y,re)
	    
 ! H^T.R^-1.[y^o-H(x_b)]
   call da_calculate_gradj(1,1,cv_size,0,0,0,0,xbx,be,iv,cv,y,shat,grid,config_flags,re)
   shat = - shat    !! Compensate for sign in calculation of grad_v (Jo)
	    
 ! A.H^T.R^-1.[y^o-H(x_b)]
   call da_amat_mul(be, grid, cv_size, ntmax(it), neign, 1.0/eignval, eignvec, qhat, shat, xhat)

   if (trace_use) call da_trace_exit ("da_kmat_mul")

end subroutine da_kmat_mul
subroutine da_lanczos_io (io_config, cv_size, ntmaxit, neign, eignvec, eignval, qhat)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write Lanczos eigenpairs
   !
   ! Called from da_solve
   !
   ! History: 08/16/2010  Creation (Tom Auligne)
   !
   !-------------------------------------------------------------------------

   implicit none
   
   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   integer,          intent(inout) :: cv_size
   integer,          intent(in)    :: ntmaxit
   integer,          intent(inout) :: neign
   real*8,           intent(inout) :: eignvec(ntmaxit, ntmaxit)
   real*8,           intent(inout) :: eignval(ntmaxit)
   real,             intent(inout) :: qhat(1:cv_size, 0:ntmaxit)

   character(len=filename_len)     :: filename               ! I/O filename
   character*10                    :: cproc
   integer                         :: ep_unit
   integer                         :: i

   if (trace_use) call da_trace_entry("da_lanczos_io")
   
   write(cproc,fmt='(i4.4)') myproc
   filename = '../lanczos_eigenpairs.'//trim(adjustl(cproc))

   call da_get_unit (ep_unit)

   if (io_config == 'r') then
      write(*,*) 'Reading Lanczos eigenpairs'
      open (unit=ep_unit, file = filename, form = 'unformatted', status = 'old')
      read(unit=ep_unit) neign, cv_size
      do i = 1, neign
         read(unit=ep_unit) eignval(i)
         read(unit=ep_unit) eignvec(1:neign,i)
         read(unit=ep_unit) qhat(1:cv_size,i)
      end do
      close(unit=ep_unit)
   else if (io_config == 'w') then
      write(*,*) 'Writing Lanczos eigenpairs'
      open (unit=ep_unit, file = filename, form = 'unformatted', status = 'replace')
      write(unit=ep_unit) neign, cv_size
      do i = 1, neign
         write(unit=ep_unit) eignval(i)
         write(unit=ep_unit) eignvec(1:neign,i)
         write(unit=ep_unit) qhat(1:cv_size,i)
      end do
      close(unit=ep_unit)
   else
      write(*,*) 'Unknow configuration for Lanczos I/O routine'
   end if

   call da_free_unit (ep_unit)

   if (trace_use) call da_trace_exit ("da_lanczos_io")

end subroutine da_lanczos_io
subroutine da_swap_xtraj ( grid )

   !-------------------------------------------------------------------------
   ! Purpose:        Swap KJ dimensions of fields from WRF to fit the fields of WRFDA
   !
   ! History: 07/16/2010  Creation (Xin Zhang )
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout)     :: grid

end subroutine da_swap_xtraj

subroutine da_read_basicstates ( xbx, grid, config_flags, timestr, fgat )

   !-------------------------------------------------------------------------
   ! Purpose: Read basic state at time = timestr
   !
   ! History: 10/01/2010  Creation (Xin Zhang )
   !
   !-------------------------------------------------------------------------

   implicit none

   character(len=256),         intent(in)    ::   timestr
   type(domain),               intent(inout) ::   grid
   type(xbx_type),             intent(inout) ::   xbx
   type(grid_config_rec_type), intent(inout) ::   config_flags
   character(len=120), intent(in), optional  ::   fgat

! Following codes are for FGAT
   if ( .not. grid%var4d .and. grid%num_fgat_time > 1 .and. present(fgat) ) then
      call da_med_initialdata_input (grid, config_flags, trim(fgat))
      call da_transfer_wrftoxb(xbx, grid, config_flags)
   end if
   return

end subroutine da_read_basicstates

end module da_minimisation
