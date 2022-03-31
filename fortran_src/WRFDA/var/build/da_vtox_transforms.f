












module da_vtox_transforms

   
   
   
   
  
   use module_dm, only : wrf_dm_sum_real,wrf_dm_sum_reals
   use module_dm, only : local_communicator, mytask, ntasks, ntasks_x, & 
      ntasks_y, data_order_xy, data_order_xyz, data_order_xzy, get_dm_max_halo_width
   use module_comm_dm, only : halo_psichi_uv_adj_sub,halo_ssmi_xa_sub,halo_sfc_xa_sub, &
      halo_radar_xa_w_sub, halo_xa_sub, halo_psichi_uv_sub, &
      halo_psichi_uv_sub, halo_xa_all_sub,halo_xb_all_sub
   use module_domain, only : xb_type, xpose_type, ep_type, vp_type, x_type, domain, get_ijk_from_grid

   use da_control, only : trace_use, var4d, ims,ime,jms,jme,kms,kme,jds,jde,kds,kde, &
      its,ite,jts,jte,kts,kte, cos_xls, cos_xle, sin_xle, sin_xls, pi, global, &
      vertical_ip,alphacv_method,use_ssmitbobs,use_rainobs, &
      use_radarobs,use_radar_rf,use_radar_rhv,use_radar_rqv,use_3dvar_phy, &
      use_ssmiretrievalobs, use_ssmt2obs, use_ssmt1obs, use_gpspwobs, use_gpsztdobs, &
      use_gpsrefobs,sfc_assi_options, test_transforms, vert_corr, fg_format, &
      fg_format_kma_global, fg_format_wrf_arw_regional,fg_format_wrf_nmm_regional, &
      ids, ide, stdout, use_rad, crtm_cloud, vert_corr_2, fg_format_wrf_arw_global, &
      alphacv_method_vp, alphacv_method_xa, vertical_ip_0, trace_use_dull, &
      ips,ipe,jps,jpe,kps,kpe, cv_size, cv_options, cv_options_hum, cloud_cv_options, &
      use_background_errors,do_normalize,use_rf,len_scaling1, len_scaling2, len_scaling3, len_scaling4, &
      len_scaling5, len_scaling6, len_scaling7, len_scaling8, len_scaling9, len_scaling10, len_scaling11

   use da_control, only : anal_type_hybrid_dual_res, myproc, num_procs,dual_res_upscale_opt
   use da_control, only : its_int,ite_int,jts_int,jte_int,kts_int,kte_int,shw, &
                          ims_int,ime_int,jms_int,jme_int,kms_int,kme_int, &
                          ids_int,ide_int,jds_int,jde_int,kds_int,kde_int, &
                          ips_int,ipe_int,jps_int,jpe_int,kps_int,kpe_int
   use da_control, only : dual_res_type, ob_locs, total_here


   use da_define_structures, only : be_type, xbx_type,da_zero_vp_type,da_zero_x
   use da_dynamics, only : da_psichi_to_uv,da_psichi_to_uv_adj
   use da_physics, only : da_uvprho_to_w_lin,da_uvprho_to_w_adj, &
      da_pt_to_rho_adj, da_pt_to_rho_lin,da_moist_phys_lin, &
      da_tprh_to_q_lin, da_tprh_to_q_adj,       &
      da_moist_phys_adj, da_transform_xtogpsref_lin, da_transform_xtotpw, &
      da_transform_xtowtq, da_transform_xtotpw_adj, &
      da_transform_xtogpsref_adj, da_transform_xtowtq_adj, &
      da_transform_xtoztd_lin, da_transform_xtoztd_adj
   use da_par_util, only : da_vv_to_cv, da_cv_to_vv

   use da_recursive_filter, only : da_transform_through_rf, &
      da_transform_through_rf_adj, da_apply_rf, da_apply_rf_adj, &
      da_transform_through_rf_dual_res, da_transform_through_rf_adj_dual_res
   use da_reporting, only : da_error, message, da_warning, da_message
   use da_spectral, only : da_vtovv_spectral,da_vtovv_spectral_adj
   use da_ssmi, only : da_transform_xtoseasfcwind_lin,da_transform_xtotb_adj, &
      da_transform_xtoseasfcwind_adj, da_transform_xtotb_lin
   use da_tools, only : da_set_boundary_xa
   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_wrf_interfaces, only : wrf_debug
   use da_wavelet, only: da_transform_through_wavelet, da_transform_through_wavelet_adj, nij

   implicit none


   contains

subroutine da_add_flow_dependence_vp (ne, ep, vp, &
                                      its,ite, jts,jte, kts,kte)

   !-----------------------------------------------------------------------
   ! Purpose: Add flow-dependent increments in control variable space (vp).
   !-----------------------------------------------------------------------
                                      
   implicit none
   
   integer, intent(in)                  :: ne  ! Ensemble size.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbations.
   type (vp_type), intent(inout)        :: vp  ! CV on grid structure.
   integer, intent(in)                  :: its, ite, jts, jte, kts, kte
   
   integer                              :: n, k! Loop counters.

   do n = 1, ne
      do k = kts, kte

         ! psi:
         vp % v1(its:ite,jts:jte,k) = vp % v1(its:ite,jts:jte,k) + &
                                      vp % alpha(its:ite,jts:jte,k,n) * &
                                      ep % v1(its:ite,jts:jte,k,n)
         ! chi_u:
         vp % v2(its:ite,jts:jte,k) = vp % v2(its:ite,jts:jte,k) + &
                                      vp % alpha(its:ite,jts:jte,k,n) * &
                                      ep % v2(its:ite,jts:jte,k,n)
         ! t_u:
         vp % v3(its:ite,jts:jte,k) = vp % v3(its:ite,jts:jte,k) + &
                                      vp % alpha(its:ite,jts:jte,k,n) * &
                                      ep % v3(its:ite,jts:jte,k,n)
         ! rh:
         vp % v4(its:ite,jts:jte,k) = vp % v4(its:ite,jts:jte,k) + &
                                      vp % alpha(its:ite,jts:jte,k,n) * &
                                      ep % v4(its:ite,jts:jte,k,n)
        ! ps_u
        vp % v5(its:ite,jts:jte,k) = vp % v5(its:ite,jts:jte,k) + &
                                     vp % alpha(its:ite,jts:jte,k,n) * &
                                     ep % v5(its:ite,jts:jte,k,n)
      end do
   end do

   if (trace_use) call da_trace_exit("da_add_flow_dependence_vp")

end subroutine da_add_flow_dependence_vp


subroutine da_add_flow_dependence_vp_adj (ne, ep, vp)

   !-----------------------------------------------------------------------
   ! Purpose: Add flow-dependent increments in control variable space (vp).
   !-----------------------------------------------------------------------
                                      
   implicit none
   
   integer, intent(in)                  :: ne  ! Ensemble size.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbations.
   type (vp_type), intent(inout)        :: vp  ! CV on grid structure.
   
   integer                              :: n, k! Loop counters.

   vp % alpha(:,:,:,:) = 0.0

   do n = ne, 1, -1
      do k = kte, kts, -1

         vp % alpha(its:ite,jts:jte,k,n) = vp % alpha(its:ite,jts:jte,k,n) + &
                                         ep % v5(its:ite,jts:jte,k,n) * &
                                         vp % v5(its:ite,jts:jte,k)

         ! rh:
         vp % alpha(its:ite,jts:jte,k,n) = vp % alpha(its:ite,jts:jte,k,n) + &
                                         ep % v4(its:ite,jts:jte,k,n) * &
                                         vp % v4(its:ite,jts:jte,k)

         ! t_u:
         vp % alpha(its:ite,jts:jte,k,n) = vp % alpha(its:ite,jts:jte,k,n) + &
                                         ep % v3(its:ite,jts:jte,k,n) * &
                                         vp % v3(its:ite,jts:jte,k)

         ! chi_u:
         vp % alpha(its:ite,jts:jte,k,n) = vp % alpha(its:ite,jts:jte,k,n) + &
                                         ep % v2(its:ite,jts:jte,k,n) * &
                                         vp % v2(its:ite,jts:jte,k)

         ! psi:
         vp % alpha(its:ite,jts:jte,k,n) = vp % alpha(its:ite,jts:jte,k,n) + &
                                         ep % v1(its:ite,jts:jte,k,n) * &
                                         vp % v1(its:ite,jts:jte,k)

      end do
   end do

   if (trace_use) call da_trace_exit("da_add_flow_dependence_vp_adj")

end subroutine da_add_flow_dependence_vp_adj


subroutine da_add_flow_dependence_xa (grid, ne,  ep, vp)

   !-----------------------------------------------------------------------
   ! Purpose: Add flow-dependent increments in model space (grid%xa).
   !-----------------------------------------------------------------------
                                      
   implicit none

   type (domain), intent(inout)         :: grid
   integer, intent(in)                  :: ne  ! Ensemble size.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbations.
   type (vp_type), intent(in)           :: vp  ! CV on grid structure.
   
   integer                              :: i, j, k, n  ! Loop counters.
   real                                 :: alpha       ! Local alpha copy.

   do k = kts, kte
      do j = jts, jte
         do i = its, ite

            do n = 1, ne
               alpha = vp % alpha(i,j,k,n)

               ! u:
               grid%xa % u(i,j,k) = grid%xa % u(i,j,k) + alpha * ep % v1(i,j,k,n) ! v1 = u

               ! v:
               grid%xa % v(i,j,k) = grid%xa % v(i,j,k) + alpha * ep % v2(i,j,k,n) ! v2 = v

               ! t:
               grid%xa % t(i,j,k) = grid%xa % t(i,j,k) + alpha * ep % v3(i,j,k,n) ! v3 = t

               ! q:
               grid%xa % q(i,j,k) = grid%xa % q(i,j,k) + alpha * ep % v4(i,j,k,n) ! v4 = q

            end do
         end do
      end do
   end do

   ! ps:
   do n = 1, ne
      grid%xa % psfc(its:ite,jts:jte) = grid%xa % psfc(its:ite,jts:jte) + &
                                        vp % alpha(its:ite,jts:jte,1,n) * &
                                        ep % v5(its:ite,jts:jte,1,n) ! v5 = ps
   end do

   if (trace_use) call da_trace_exit("da_add_flow_dependence_xa")

end subroutine da_add_flow_dependence_xa


subroutine da_add_flow_dependence_xa_dual_res (grid, ne,  ep, vp)

   !-----------------------------------------------------------------------
   ! Purpose: Add flow-dependent increments in model space (grid%xa).
   !-----------------------------------------------------------------------
                                      
   implicit none

   type (domain), intent(inout)         :: grid
   integer, intent(in)                  :: ne  ! Ensemble size.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbations.
   type (vp_type), intent(in)           :: vp  ! CV on grid structure.
   
   integer                              :: i, j, k, n  ! Loop counters.
   real                                 :: alpha       ! Local alpha copy.

   real, allocatable, dimension(:,:,:)  :: ens_contrib_u, ens_contrib_v, ens_contrib_t, ens_contrib_q, ens_contrib_p
   real, allocatable, dimension(:,:,:)  :: output_u,output_v,output_t,output_q, output_p

   integer  :: thisdomain_max_halo_width

   integer  :: cids, cide, ckds, ckde, cjds, cjde, &
               cims, cime, ckms, ckme, cjms, cjme, &
               cips, cipe, ckps, ckpe, cjps, cjpe, &
               nids, nide, nkds, nkde, njds, njde, &
               nims, nime, nkms, nkme, njms, njme, &
               nips, nipe, nkps, nkpe, njps, njpe

   integer :: nj, cj, nk, ck, ni, ci ! for testing

   ! HALO STUFF
   integer :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
   integer :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
   logical, external :: rsl_comm_iter


   ! Get coarse (ensemble) grid dimensions ( grid%intermediate_grid)
   CALL get_ijk_from_grid (  grid%intermediate_grid ,               &
                             cids, cide, cjds, cjde, ckds, ckde,    &
                             cims, cime, cjms, cjme, ckms, ckme,    &
                             cips, cipe, cjps, cjpe, ckps, ckpe    )

   ! Get fine (analysis) grid dimensions ( grid)
   CALL get_ijk_from_grid (  grid,                                  &
                             nids, nide, njds, njde, nkds, nkde,    &
                             nims, nime, njms, njme, nkms, nkme,    &
                             nips, nipe, njps, njpe, nkps, nkpe )
   CALL get_dm_max_halo_width ( grid%id , thisdomain_max_halo_width ) ! Can omit?

   ! Input: Ensemble contribution to increment -- low-res domain (x,z,y) order
   allocate( ens_contrib_u(cims:cime,ckms:ckme,cjms:cjme) ) 
   allocate( ens_contrib_v(cims:cime,ckms:ckme,cjms:cjme) ) 
   allocate( ens_contrib_t(cims:cime,ckms:ckme,cjms:cjme) ) 
   allocate( ens_contrib_q(cims:cime,ckms:ckme,cjms:cjme) ) 
   allocate( ens_contrib_p(cims:cime,1:1,cjms:cjme) )

   ! Output: Ensemble contribution to increment interpolated to hi-res domain, (x,z,y) order
   allocate( output_u(nims:nime,nkms:nkme,njms:njme) )
   allocate( output_v(nims:nime,nkms:nkme,njms:njme) )
   allocate( output_t(nims:nime,nkms:nkme,njms:njme) )
   allocate( output_q(nims:nime,nkms:nkme,njms:njme) )
   allocate( output_p(nims:nime,1:1,njms:njme) )

   ens_contrib_u = 0.
   ens_contrib_v = 0.
   ens_contrib_t = 0.
   ens_contrib_q = 0.
   ens_contrib_p = 0.

   output_u = 0.
   output_v = 0.
   output_t = 0.
   output_q = 0.
   output_p = 0.

   !
   ! Determine the ensemble contribution to the increment (low-res) and put in (x,z,y) order for interpolation
   !

   do j = jts_int, jte_int
      do k = kts_int, kte_int
         do i = its_int, ite_int
            do n = 1, ne

               alpha = vp % alpha(i,j,k,n)

               ens_contrib_u(i,k,j) = ens_contrib_u(i,k,j) + alpha * ep % v1(i,j,k,n) ! v1 = u
               ens_contrib_v(i,k,j) = ens_contrib_v(i,k,j) + alpha * ep % v2(i,j,k,n) ! v2 = v
               ens_contrib_t(i,k,j) = ens_contrib_t(i,k,j) + alpha * ep % v3(i,j,k,n) ! v3 = t
               ens_contrib_q(i,k,j) = ens_contrib_q(i,k,j) + alpha * ep % v4(i,j,k,n) ! v4 = q

            end do
         end do
      end do
   end do

   do n = 1,ne
      ens_contrib_p(its_int:ite_int,1,jts_int:jte_int) = ens_contrib_p(its_int:ite_int,1,jts_int:jte_int) + & 
                                                     vp % alpha(its_int:ite_int,jts_int:jte_int,1,n)  * &
                                                     ep % v5   (its_int:ite_int,jts_int:jte_int,1,n)  ! v5 = ps
   end do


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!! DO HALO STUFF !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL rsl_comm_iter_init(4,cjps,cjpe)
   DO WHILE ( rsl_comm_iter( grid%intermediate_grid%id , grid%intermediate_grid%is_intermediate, 4 , &
			    0 , cjds,cjde,cjps,cjpe, grid%intermediate_grid%njds, grid%intermediate_grid%njde, &
	                    rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
	                    rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))

      CALL RSL_LITE_INIT_EXCH ( local_communicator, 4, 0, &
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   &
	   4, 1, 8, &
	   0, 0, 4, &
	   0, 0, 8, &
	    0,  0, 4, &
	    myproc, ntasks, ntasks_x, ntasks_y,   &
	    cips, cipe, cjps, cjpe, ckps, MAX(1,1&
      ,ckpe &
      ))

      IF ( SIZE(ens_contrib_u,1)*SIZE(ens_contrib_u,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_u, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_v,1)*SIZE(ens_contrib_v,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_v, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_t,1)*SIZE(ens_contrib_t,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_t, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_q,1)*SIZE(ens_contrib_q,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_q, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_p,1)*SIZE(ens_contrib_p,2) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_p(:,1,:), 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 0, DATA_ORDER_XY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, 1, 1,    &
         cims, cime, cjms, cjme, 1, 1,    &
         cips, cipe, cjps, cjpe, 1, 1    )
      ENDIF


     CALL RSL_LITE_EXCH_Y ( local_communicator , myproc, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )

      IF ( SIZE(ens_contrib_u,1)*SIZE(ens_contrib_u,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_u, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_v,1)*SIZE(ens_contrib_v,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_v, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_t,1)*SIZE(ens_contrib_t,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_t, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_q,1)*SIZE(ens_contrib_q,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_q, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_p,1)*SIZE(ens_contrib_p,2) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_p(:,1,:), 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 0, 1, DATA_ORDER_XY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, 1, 1,    &
         cims, cime, cjms, cjme, 1, 1,    &
         cips, cipe, cjps, cjpe, 1, 1    )
      ENDIF
   END DO


   CALL rsl_comm_iter_init(4,cips,cipe)
   DO WHILE ( rsl_comm_iter( grid%intermediate_grid%id , grid%intermediate_grid%is_intermediate, 4 , &
			    1 , cids,cide,cips,cipe, grid%intermediate_grid%nids, grid%intermediate_grid%nide, &
	                    rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
	                    rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))

      CALL RSL_LITE_INIT_EXCH ( local_communicator, 4, 1, &
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   &
	   4, 1, 8, &
	   0, 0, 4, &
	   0, 0, 8, &
	    0,  0, 4, &
	    myproc, ntasks, ntasks_x, ntasks_y,   &
	    cips, cipe, cjps, cjpe, ckps, MAX(1,1&
      ,ckpe &
      ))

      IF ( SIZE(ens_contrib_u,1)*SIZE(ens_contrib_u,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_u, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_v,1)*SIZE(ens_contrib_v,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_v, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_t,1)*SIZE(ens_contrib_t,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_t, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_q,1)*SIZE(ens_contrib_q,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_q, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 0, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_p,1)*SIZE(ens_contrib_p,2) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_p(:,1,:), 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 0, DATA_ORDER_XY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, 1, 1,    &
         cims, cime, cjms, cjme, 1, 1,    &
         cips, cipe, cjps, cjpe, 1, 1    )
      ENDIF

     CALL RSL_LITE_EXCH_X ( local_communicator , myproc, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )

      IF ( SIZE(ens_contrib_u,1)*SIZE(ens_contrib_u,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_u, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_v,1)*SIZE(ens_contrib_v,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_v, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_t,1)*SIZE(ens_contrib_t,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_t, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_q,1)*SIZE(ens_contrib_q,3) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_q, 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 1, DATA_ORDER_XZY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, ckds, ckde,    &
         cims, cime, cjms, cjme, ckms, ckme,    &
         cips, cipe, cjps, cjpe, ckps, ckpe    )
      ENDIF
      IF ( SIZE(ens_contrib_p,1)*SIZE(ens_contrib_p,2) .GT. 1 ) THEN
	 CALL RSL_LITE_PACK ( local_communicator,&
	  ens_contrib_p(:,1,:), 4,&
	 rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	 rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	 8, 1, 1, DATA_ORDER_XY, 0, &
	 myproc, ntasks, ntasks_x, ntasks_y,       &
         cids, cide, cjds, cjde, 1, 1,    &
         cims, cime, cjms, cjme, 1, 1,    &
         cips, cipe, cjps, cjpe, 1, 1    )
      ENDIF
   ENDDO

   !!!!!! END HALO STUFF !!!!!!!!!!!!!

   !------------------------------------------------------------------------------
   ! Now, interpolate the  ensemble contributions to in increment to the high-res grid
   !------------------------------------------------------------------------------

   ! see .../frame/module_dm.f90 and look for interp_fcn, patch_rsl_lite
   ! see .../share/mediation_interp_domain.F
   ! Note, grid%intermediate_grid%domain = 2
 
   ! Input is first entry in interp_fcn (low-res ensemble contribution to increment in (x,z,y) order)
   ! Output is the ensemble contribution to the increment on the hi-res grid in (x,z,y) order

   if ( dual_res_upscale_opt .le. 2 ) then
      ! U
      call interp_fcn (  &
                  ens_contrib_u,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  output_u,  &   ! ND field...nested grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


     ! For testing
     DO nj = njps, njpe
        cj =  grid%j_parent_start + (nj-1) / grid%parent_grid_ratio     ! j coord of CD point 
           ck = nk
           DO ni = nips, nipe
              ci = grid%i_parent_start + (ni-1) / grid%parent_grid_ratio      ! j coord of CD point 
           ENDDO
     ENDDO
     ! For testing

      ! V
      call interp_fcn (  &
                  ens_contrib_v,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  output_v,  &   ! ND field...nested grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )

      ! T
      call interp_fcn (  &
                  ens_contrib_t,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  output_t,  &   ! ND field...nested grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )

      ! Q
      call interp_fcn (  &
                  ens_contrib_q,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         ! CD dims
                  output_q,  &   ! ND field...nested grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         ! ND dims
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )

      ! P
      call interp_fcn (  &
                  ens_contrib_p,   &       ! CD field ... intermediate grid
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         ! CD dims
                  output_p,  &   ! ND field...nested grid
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         ! ND dims
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )

   else

      !call da_message((/"Using adjoint-based interpolation"/))

      !$OMP PARALLEL DO &
      !$OMP PRIVATE (n,k)
      do n=1,total_here
	 do k = kts, kte

	    output_u(ob_locs(n)%xx,k,ob_locs(n)%yy) = &
			ob_locs(n)%dym*(ob_locs(n)%dxm*ens_contrib_u(ob_locs(n)%i,k,ob_locs(n)%j) &
		      + ob_locs(n)%dx*ens_contrib_u(ob_locs(n)%i+1,k,ob_locs(n)%j)) &
		      + ob_locs(n)%dy *(ob_locs(n)%dxm*ens_contrib_u(ob_locs(n)%i,k,ob_locs(n)%j+1)  &
		      + ob_locs(n)%dx*ens_contrib_u(ob_locs(n)%i+1,k,ob_locs(n)%j+1))

	    output_v(ob_locs(n)%xx,k,ob_locs(n)%yy) = &
			ob_locs(n)%dym*(ob_locs(n)%dxm*ens_contrib_v(ob_locs(n)%i,k,ob_locs(n)%j) &
		      + ob_locs(n)%dx*ens_contrib_v(ob_locs(n)%i+1,k,ob_locs(n)%j)) &
		      + ob_locs(n)%dy *(ob_locs(n)%dxm*ens_contrib_v(ob_locs(n)%i,k,ob_locs(n)%j+1)  &
		      + ob_locs(n)%dx*ens_contrib_v(ob_locs(n)%i+1,k,ob_locs(n)%j+1))

	    output_t(ob_locs(n)%xx,k,ob_locs(n)%yy) = &
			ob_locs(n)%dym*(ob_locs(n)%dxm*ens_contrib_t(ob_locs(n)%i,k,ob_locs(n)%j) &
		      + ob_locs(n)%dx*ens_contrib_t(ob_locs(n)%i+1,k,ob_locs(n)%j)) &
		      + ob_locs(n)%dy *(ob_locs(n)%dxm*ens_contrib_t(ob_locs(n)%i,k,ob_locs(n)%j+1)  &
		      + ob_locs(n)%dx*ens_contrib_t(ob_locs(n)%i+1,k,ob_locs(n)%j+1))

	    output_q(ob_locs(n)%xx,k,ob_locs(n)%yy) = &
			ob_locs(n)%dym*(ob_locs(n)%dxm*ens_contrib_q(ob_locs(n)%i,k,ob_locs(n)%j) &
		      + ob_locs(n)%dx*ens_contrib_q(ob_locs(n)%i+1,k,ob_locs(n)%j)) &
		      + ob_locs(n)%dy *(ob_locs(n)%dxm*ens_contrib_q(ob_locs(n)%i,k,ob_locs(n)%j+1)  &
		      + ob_locs(n)%dx*ens_contrib_q(ob_locs(n)%i+1,k,ob_locs(n)%j+1))

	 end do

	 output_p(ob_locs(n)%xx,1,ob_locs(n)%yy) = &
	    ob_locs(n)%dym*(ob_locs(n)%dxm*ens_contrib_p(ob_locs(n)%i,1,ob_locs(n)%j) &
	    + ob_locs(n)%dx*ens_contrib_p(ob_locs(n)%i+1,1,ob_locs(n)%j)) &
	    + ob_locs(n)%dy *(ob_locs(n)%dxm*ens_contrib_p(ob_locs(n)%i,1,ob_locs(n)%j+1)  &
	    + ob_locs(n)%dx*ens_contrib_p(ob_locs(n)%i+1,1,ob_locs(n)%j+1))

      end do
      !$OMP END PARALLEL DO

   endif


   ! 
   ! Now add the hi-res ensemble contribution to the increment to the static increment.
   !  This forms the total hi-res increment
   ! 

   do k = kts, kte
      do j = jts, jte
         do i = its, ite

               ! u:
               grid%xa % u(i,j,k) = grid%xa % u(i,j,k) + output_u(i,k,j) ! u

               ! v:
               grid%xa % v(i,j,k) = grid%xa % v(i,j,k) + output_v(i,k,j) ! v

               ! t:
               grid%xa % t(i,j,k) = grid%xa % t(i,j,k) + output_t(i,k,j) ! t

               ! q:
               grid%xa % q(i,j,k) = grid%xa % q(i,j,k) + output_q(i,k,j) ! q

         end do
      end do
   end do

   ! ps:
   grid%xa % psfc(its:ite,jts:jte) = grid%xa % psfc(its:ite,jts:jte) + &
                                      output_p(its:ite,1,jts:jte) ! ps

   !
   ! Clean-up
   !

   if (trace_use) call da_trace_exit("da_add_flow_dependence_xa_dual_res")

   deallocate(ens_contrib_u,ens_contrib_v,ens_contrib_t,ens_contrib_q,ens_contrib_p)
   deallocate(output_u,output_v,output_t,output_q,output_p)

end subroutine da_add_flow_dependence_xa_dual_res


subroutine da_add_flow_dependence_xa_adj (ne, ep, xa, vp)

   !-----------------------------------------------------------------------
   ! Purpose: Add flow-dependent increments in model space (xa).
   !-----------------------------------------------------------------------
                                      
   implicit none
   
   integer, intent(in)                  :: ne  ! Ensemble size.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbations.
   type (x_type), intent(in)            :: xa  ! Analysis increments.
   type (vp_type), intent(inout)        :: vp  ! CV on grid structure.
   
   integer                              :: i, j, k, n ! Loop counters.
   real                                 :: alpha       ! Local alpha copy.

   vp % alpha = 0.0

   do n = ne, 1, -1
      ! ps:
      vp % alpha(its:ite,jts:jte,1,n) = vp % alpha(its:ite,jts:jte,1,n) + &
                                        ep % v5(its:ite,jts:jte,1,n) * & ! v5 = ps
                                        xa % psfc(its:ite,jts:jte)
   end do

   do k = kte, kts, -1
      do j = jte, jts, -1
         do i = ite, its, -1

            do n = ne, 1, -1
               alpha = 0.0
               alpha = alpha + ep % v4(i,j,k,n) * xa % q(i,j,k)
               alpha = alpha + ep % v3(i,j,k,n) * xa % t(i,j,k)
               alpha = alpha + ep % v2(i,j,k,n) * xa % v(i,j,k)
               alpha = alpha + ep % v1(i,j,k,n) * xa % u(i,j,k)
               vp % alpha(i,j,k,n) = vp % alpha(i,j,k,n) + alpha
            end do
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_add_flow_dependence_xa_adj")

end subroutine da_add_flow_dependence_xa_adj

subroutine da_add_flow_dependence_xa_adj_dual_res (ne, ep, grid, vp)

   !-----------------------------------------------------------------------
   ! Purpose: Add flow-dependent increments in model space (xa).
   !-----------------------------------------------------------------------
                                      
   implicit none
   
   integer, intent(in)                  :: ne  ! Ensemble size.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbations.
   type (domain), intent(in)            :: grid  ! Analysis increments.
   type (vp_type), intent(inout)        :: vp  ! CV on grid structure.
   
   integer                              :: i, j, k, n ! Loop counters.
   real                                 :: alpha       ! Local alpha copy.

   real, allocatable, dimension(:,:,:)  :: output_u,output_v,output_t,output_q, output_p
   real, allocatable, dimension(:,:,:)  :: input_u,input_v,input_t,input_q,input_p
   
   integer  :: cids, cide, ckds, ckde, cjds, cjde, &
               cims, cime, ckms, ckme, cjms, cjme, &
               cips, cipe, ckps, ckpe, cjps, cjpe, &
               nids, nide, nkds, nkde, njds, njde, &
               nims, nime, nkms, nkme, njms, njme, &
               nips, nipe, nkps, nkpe, njps, njpe

   ! HALO STUFF
   integer :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
   integer :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
   logical, external :: rsl_comm_iter


   ! Get coarse (ensemble) grid dimensions ( grid%intermediate_grid)
   CALL get_ijk_from_grid (  grid%intermediate_grid ,               &
                             cids, cide, cjds, cjde, ckds, ckde,    &
                             cims, cime, cjms, cjme, ckms, ckme,    &
                             cips, cipe, cjps, cjpe, ckps, ckpe    )

   ! Get fine (analysis) grid dimensions (grid)
   CALL get_ijk_from_grid (  grid,                                  &
                             nids, nide, njds, njde, nkds, nkde,    &
                             nims, nime, njms, njme, nkms, nkme,    &
                             nips, nipe, njps, njpe, nkps, nkpe   )
 
   !
   ! Allocate and initialize arrays
   !

   ! Input is hi-res domain
   allocate( input_u(nims:nime,nkms:nkme,njms:njme) )
   allocate( input_v(nims:nime,nkms:nkme,njms:njme) )
   allocate( input_t(nims:nime,nkms:nkme,njms:njme) )
   allocate( input_q(nims:nime,nkms:nkme,njms:njme) )
   allocate( input_p(nims:nime,1:1,njms:njme) )

   ! Output is low-res domain
   allocate( output_u(cims:cime,ckms:ckme,cjms:cjme) )
   allocate( output_v(cims:cime,ckms:ckme,cjms:cjme) )
   allocate( output_t(cims:cime,ckms:ckme,cjms:cjme) )
   allocate( output_q(cims:cime,ckms:ckme,cjms:cjme) )
   allocate( output_p(cims:cime,1:1,cjms:cjme) )

   output_u = 0. ; input_u = 0.
   output_v = 0. ; input_v = 0.
   output_t = 0. ; input_t = 0.
   output_q = 0. ; input_q = 0.
   output_p = 0. ; input_p = 0.


   !
   ! Get input (hi-res) data into (x,z,y) order for interpolation
   !

   do j = jte, jts, -1
      do k = kte, kts, -1
         do i = ite, its, -1
            input_u(i,k,j) = grid%xa%u(i,j,k)
            input_v(i,k,j) = grid%xa%v(i,j,k)
            input_t(i,k,j) = grid%xa%t(i,j,k)
            input_q(i,k,j) = grid%xa%q(i,j,k)
         end do
      end do
   end do

   input_p(:,1,:) = grid%xa%psfc(:,:)

   write (unit=message(1),fmt='(A,2I8)')' istart, jstart = ',grid%i_parent_start, grid%j_parent_start
   call wrf_debug(2, message(1))
   write (unit=message(2),fmt='(A,4I8)')' grid%j_start(1), grid%j_end(1), jts_int, jte_int', grid%j_start(1), grid%j_end(1),jts_int,jte_int
   call wrf_debug(2, message(2))
   write (unit=message(3),fmt='(A,2F12.5)')' min/max  Input U = ',minval(input_u),maxval(input_u)
   call wrf_debug(2, message(3))
   !call da_message(message(1:3))



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!! DO HALO STUFF !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL rsl_comm_iter_init(4,njps,njpe)
   DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 4 , &
                         0 , njds,njde,njps,njpe, grid%njds, grid%njde, &
                         rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
                         rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))

      CALL RSL_LITE_INIT_EXCH ( local_communicator, 4, 0, &
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   &
	   4, 1, 8, &
	   0, 0, 4, &
	   0, 0, 8, &
	    0,  0, 4, &
	    myproc, ntasks, ntasks_x, ntasks_y,   &
	    nips, nipe, njps, njpe, nkps, MAX(1,1&
      ,nkpe&
      ))

	IF ( SIZE(input_u,1)*SIZE(input_u,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_u, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_v,1)*SIZE(input_v,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_v, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_t,1)*SIZE(input_t,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_t, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_q,1)*SIZE(input_q,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_q, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_p,1)*SIZE(input_p,2) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_p(:,1,:), 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 0, DATA_ORDER_XY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, 1, 1,    &
           nims, nime, njms, njme, 1, 1,    &
           nips, nipe, njps, njpe, 1, 1   )
	ENDIF

        CALL RSL_LITE_EXCH_Y ( local_communicator , myproc, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )

	IF ( SIZE(input_u,1)*SIZE(input_u,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_u, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_v,1)*SIZE(input_v,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_v, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_t,1)*SIZE(input_t,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_t, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_q,1)*SIZE(input_q,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_q, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_p,1)*SIZE(input_p,2) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_p(:,1,:), 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 0, 1, DATA_ORDER_XY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, 1, 1,    &
           nims, nime, njms, njme, 1, 1,    &
           nips, nipe, njps, njpe, 1, 1   )
	ENDIF

      ENDDO

 
   CALL rsl_comm_iter_init(4,nips,nipe)
   DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 4 , &
                         1 , nids,nide,nips,nipe, grid%nids, grid%nide, &
                         rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
                         rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))

      CALL RSL_LITE_INIT_EXCH ( local_communicator, 4, 1, &
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   &
	   4, 1, 8, &
	   0, 0, 4, &
	   0, 0, 8, &
	    0,  0, 4, &
	    myproc, ntasks, ntasks_x, ntasks_y,   &
	    nips, nipe, njps, njpe, nkps, MAX(1,1&
      ,nkpe&
      ))

	IF ( SIZE(input_u,1)*SIZE(input_u,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_u, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_v,1)*SIZE(input_v,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_v, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_t,1)*SIZE(input_t,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_t, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_q,1)*SIZE(input_q,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_q, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 0, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_p,1)*SIZE(input_p,2) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_p(:,1,:), 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 0, DATA_ORDER_XY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, 1, 1,    &
           nims, nime, njms, njme, 1, 1,    &
           nips, nipe, njps, njpe, 1, 1   )
	ENDIF


      CALL RSL_LITE_EXCH_X ( local_communicator , myproc, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )


	IF ( SIZE(input_u,1)*SIZE(input_u,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_u, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_v,1)*SIZE(input_v,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_v, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_t,1)*SIZE(input_t,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_t, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_q,1)*SIZE(input_q,3) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_q, 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 1, DATA_ORDER_XZY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, nkds, nkde,    &
           nims, nime, njms, njme, nkms, nkme,    &
           nips, nipe, njps, njpe, nkps, nkpe   )
	ENDIF
	IF ( SIZE(input_p,1)*SIZE(input_p,2) .GT. 1 ) THEN
	   CALL RSL_LITE_PACK ( local_communicator,&
	   input_p(:,1,:), 4,&
	   rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
	   rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
	   8, 1, 1, DATA_ORDER_XY, 0, &
	   myproc, ntasks, ntasks_x, ntasks_y,       &
           nids, nide, njds, njde, 1, 1,    &
           nims, nime, njms, njme, 1, 1,    &
           nips, nipe, njps, njpe, 1, 1   )
	ENDIF

   ENDDO

   !!!!!!!!! END HALO STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!


   !
   ! Now interpolate grid%xa from  hi-res to low-res
   !  Output is first entry in copy_fcn/copy_fcnm subroutine, and this is the low-resolution grid%xa in (x,z,y) order
   !

   if ( dual_res_upscale_opt .eq. 1 ) then

      ! U 
      CALL copy_fcnm(  &
                  output_u,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_u,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


      ! V 
      CALL copy_fcnm (  &
                  output_v,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_v,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


      ! T 
      CALL copy_fcnm (  &
                  output_t,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_t,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


      ! Q 
      CALL copy_fcnm (  &
                  output_q,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_q,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )



      ! P 
      CALL copy_fcnm (  &
                  output_p,   &       ! CD field ... intermediate grid
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  input_p,  &   ! ND field ... fine grid
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )
 
   else if ( dual_res_upscale_opt .eq. 2 ) then
      ! U 

      CALL copy_fcn(  &
                  output_u,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_u,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


      ! V 
      CALL copy_fcn (  &
                  output_v,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_v,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


      ! T 
      CALL copy_fcn (  &
                  output_t,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_t,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )


      ! Q 
      CALL copy_fcn (  &
                  output_q,   &       ! CD field ... intermediate grid
                 cids, cide, ckds, ckde, cjds, cjde,   &         ! CD dims
                 cims, cime, ckms, ckme, cjms, cjme,   &         ! CD dims
                 cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe,   &         
                  input_q,  &   ! ND field ... fine grid
                 nids, nide, nkds, nkde, njds, njde,   &         ! ND dims
                 nims, nime, nkms, nkme, njms, njme,   &         ! ND dims
                 nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )



      ! P 
      CALL copy_fcn (  &
                  output_p,   &       ! CD field ... intermediate grid
                 cids, cide, 1, 1, cjds, cjde,   &         ! CD dims
                 cims, cime, 1, 1, cjms, cjme,   &         ! CD dims
                 cips, cipe, 1, 1, cjps, cjpe,   &         
                  input_p,  &   ! ND field ... fine grid
                 nids, nide, 1, 1, njds, njde,   &         ! ND dims
                 nims, nime, 1, 1, njms, njme,   &         ! ND dims
                 nips, nipe, 1, 1, njps, njpe,   &         
                  shw(1), grid%imask_nostag,         &         ! stencil half width
                  .FALSE., .FALSE.,                                                &         ! xstag, ystag
                  grid%i_parent_start, grid%j_parent_start,                     &
                  grid%parent_grid_ratio, grid%parent_grid_ratio                &
                  )
   else if ( dual_res_upscale_opt .eq. 3 ) then

      !call da_message((/'Using adjoint of bilinear interpolation'/))
      do n=1,total_here
         do k = kts,kte
	    output_u(ob_locs(n)%i  ,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dxm  * input_u( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_u(ob_locs(n)%i  ,k,ob_locs(n)%j)   
	    output_u(ob_locs(n)%i+1,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dx   * input_u( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_u(ob_locs(n)%i+1,k,ob_locs(n)%j)   
	    output_u(ob_locs(n)%i  ,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dxm  * input_u( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_u(ob_locs(n)%i  ,k,ob_locs(n)%j+1) 
	    output_u(ob_locs(n)%i+1,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dx   * input_u( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_u(ob_locs(n)%i+1,k,ob_locs(n)%j+1) 

	    output_v(ob_locs(n)%i  ,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dxm  * input_v( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_v(ob_locs(n)%i  ,k,ob_locs(n)%j)   
	    output_v(ob_locs(n)%i+1,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dx   * input_v( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_v(ob_locs(n)%i+1,k,ob_locs(n)%j)   
	    output_v(ob_locs(n)%i  ,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dxm  * input_v( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_v(ob_locs(n)%i  ,k,ob_locs(n)%j+1) 
	    output_v(ob_locs(n)%i+1,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dx   * input_v( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_v(ob_locs(n)%i+1,k,ob_locs(n)%j+1) 
 
	    output_t(ob_locs(n)%i  ,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dxm  * input_t( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_t(ob_locs(n)%i  ,k,ob_locs(n)%j)   
	    output_t(ob_locs(n)%i+1,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dx   * input_t( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_t(ob_locs(n)%i+1,k,ob_locs(n)%j)   
	    output_t(ob_locs(n)%i  ,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dxm  * input_t( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_t(ob_locs(n)%i  ,k,ob_locs(n)%j+1) 
	    output_t(ob_locs(n)%i+1,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dx   * input_t( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_t(ob_locs(n)%i+1,k,ob_locs(n)%j+1) 

	    output_q(ob_locs(n)%i  ,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dxm  * input_q( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_q(ob_locs(n)%i  ,k,ob_locs(n)%j)   
	    output_q(ob_locs(n)%i+1,k,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dx   * input_q( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_q(ob_locs(n)%i+1,k,ob_locs(n)%j)   
	    output_q(ob_locs(n)%i  ,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dxm  * input_q( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_q(ob_locs(n)%i  ,k,ob_locs(n)%j+1) 
	    output_q(ob_locs(n)%i+1,k,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dx   * input_q( ob_locs(n)%xx,k,ob_locs(n)%yy) + output_q(ob_locs(n)%i+1,k,ob_locs(n)%j+1) 
	 end do

	 output_p(ob_locs(n)%i  ,1,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dxm  * input_p( ob_locs(n)%xx,1,ob_locs(n)%yy) + output_p(ob_locs(n)%i  ,1,ob_locs(n)%j)   
	 output_p(ob_locs(n)%i+1,1,ob_locs(n)%j)   = ob_locs(n)%dym * ob_locs(n)%dx   * input_p( ob_locs(n)%xx,1,ob_locs(n)%yy) + output_p(ob_locs(n)%i+1,1,ob_locs(n)%j)   
	 output_p(ob_locs(n)%i  ,1,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dxm  * input_p( ob_locs(n)%xx,1,ob_locs(n)%yy) + output_p(ob_locs(n)%i  ,1,ob_locs(n)%j+1) 
	 output_p(ob_locs(n)%i+1,1,ob_locs(n)%j+1) = ob_locs(n)%dy  * ob_locs(n)%dx   * input_p( ob_locs(n)%xx,1,ob_locs(n)%yy) + output_p(ob_locs(n)%i+1,1,ob_locs(n)%j+1) 

      end do

   else 

      write(unit=message(1),fmt='(A,I4)') 'Invalid value of dual_res_upscale_opt: ', dual_res_upscale_opt
      write(unit=message(2),fmt='(A)') 'Valid values are 1, 2, or 3'
      call da_error("da_add_flow_dependence_xa_adj_dual_res.inc",596,message(1:2))

   endif

   write (unit=message(1), fmt='(A,2F12.5)') ' min/max U  Input = ',minval(input_u),maxval(input_u)
   call wrf_debug(2, message(1))
   write (unit=message(2), fmt='(A,2F12.5)') ' min/max U Output = ',minval(output_u),maxval(output_u)
   call wrf_debug(2, message(2))
   write (unit=message(3), fmt='(A,2F12.5)') ' min/max V  Input = ',minval(input_v),maxval(input_v)
   call wrf_debug(2, message(3))
   write (unit=message(4), fmt='(A,2F12.5)') ' min/max V Output = ',minval(output_v),maxval(output_v)
   call wrf_debug(2, message(4))
   write (unit=message(5), fmt='(A,2F12.5)') ' min/max T  Input = ',minval(input_t),maxval(input_t)
   call wrf_debug(2, message(5))
   write (unit=message(6), fmt='(A,2F12.5)') ' min/max T Output = ',minval(output_t),maxval(output_t)
   call wrf_debug(2, message(6))
   write (unit=message(7), fmt='(A,2F12.5)') ' min/max Q  Input = ',minval(input_q),maxval(input_q)
   call wrf_debug(2, message(7))
   write (unit=message(8), fmt='(A,2F12.5)') ' min/max Q Output = ',minval(output_q),maxval(output_q)
   call wrf_debug(2, message(8))
   write (unit=message(9), fmt='(A,2F12.5)') ' min/max P  Input = ',minval(input_p),maxval(input_p)
   call wrf_debug(2, message(9))
   write (unit=message(10),fmt='(A,2F12.5)') ' min/max P Output = ',minval(output_p),maxval(output_p)
   call wrf_debug(2, message(10))
   !call da_message(message(1:10))

    !
    ! Smooth things here if you want...call sm121 from interp_fcn.F. Input is similar to copy_fcn
    !    For dual-res application, smoothing is probably not necessary
    !
    !Input and output is first entry

! SUBROUTINE sm121 ( output_p , &
!                     ids_int, ide_int + 1, kds_int, kde_int + 1, jds_int, jde_int + 1,   &         ! CD dims
!                     ims_int, ime_int, kms_int, kme_int, jms_int, jme_int,   &         ! CD dims
!                     its_int, ite_int, kts_int, kte_int, jts_int, jte_int,   &         ! CD dims
!                     .FALSE., .FALSE.,                     &  ! staggering of field
!                    ids, ide, kds, kde, jds, jde,   &
!                     ims, ime, kms, kme, jms, jme,   &
!                     its, ite, kts, kte, jts, jte,   &
!                     grid%i_parent_start, grid%j_parent_start, &
!                     grid%parent_grid_ratio, grid%parent_grid_ratio &  ! Position of lower left of nest in 
!                     )



   ! 
   ! Now compute alpha on the low-res domain
   !

   vp % alpha = 0.0

   do n = ne, 1, -1
      ! ps:
      vp % alpha(its_int:ite_int,jts_int:jte_int,1,n) = vp % alpha(its_int:ite_int,jts_int:jte_int,1,n) + &
                                        ep % v5(its_int:ite_int,jts_int:jte_int,1,n) * & ! v5 = ps
                                       output_p(its_int:ite_int,1,jts_int:jte_int)
   end do


   do j = jte_int, jts_int, -1
      do k = kte_int, kts_int, -1
         do i = ite_int, its_int, -1

            do n = ne, 1, -1
               alpha = 0.0
               alpha = alpha + ep % v4(i,j,k,n) * output_q(i,k,j)
               alpha = alpha + ep % v3(i,j,k,n) * output_t(i,k,j)
               alpha = alpha + ep % v2(i,j,k,n) * output_v(i,k,j)
               alpha = alpha + ep % v1(i,j,k,n) * output_u(i,k,j)
               vp % alpha(i,j,k,n) = vp % alpha(i,j,k,n) + alpha
            end do
         end do
      end do
   end do

   ! Clean-up
   deallocate(output_u,output_v,output_t,output_q,output_p)
   deallocate(input_u,input_v,input_t,input_q,input_p)

   if (trace_use) call da_trace_exit("da_add_flow_dependence_xa_adj_dual_res")

end subroutine da_add_flow_dependence_xa_adj_dual_res

subroutine da_check_eof_decomposition(be_eigenval, be_eigenvec, name)

   !---------------------------------------------------------------------------
   ! Purpose: Check eigenvectors E of vertical covariance matrix B_{x} which 
   ! have been read in from NMC-statistics file.
   !
   ! Method:  E and L (eigenvalues) are computed using LAPACK/BLAS software in 
   ! the NMC code using the definition E^{T} B_{x} E = L. This routine checks 
   ! for eigenvector orthogonality and completeness as defined below.
   !---------------------------------------------------------------------------

   implicit none
      
   real*8, intent(in)            :: be_eigenval(:)   ! Back. error eigenvalues.
   real*8, intent(in)            :: be_eigenvec(:,:) ! Back. error eigenvector
   character(len=*), intent(in)  :: name             ! Variable name.

   integer                       :: kz               ! Size of 3rd dimension.
   integer                       :: k, k1, k2, m     ! Loop counters
   real                          :: tot_variance     ! Total variance.
   real                          :: cum_variance     ! Cumulative variance.
   real                          :: max_off_diag     ! Maximum off-diagonal.
      
   real, allocatable             :: matrix2(:,:)     ! 2D Work matrix.
   logical, allocatable          :: array_mask(:)    ! Array mask for MAXVAL.

   if (trace_use) call da_trace_entry("da_check_eof_decomposition")

   !----------------------------------------------------------------------
   ! [1.0]: Initialise:
   !----------------------------------------------------------------------  

   kz = size(be_eigenval)
                         
   !----------------------------------------------------------------------
   ! [2.0]: Print out global eigenvalues (used for truncation):
   !----------------------------------------------------------------------  

   cum_variance = 0.0
   tot_variance = sum(be_eigenval(1:kz))

   write(unit=stdout,fmt='(A)') 'Name Mode  Eigenvalue Normalised Variance'
   do k = 1, kz
      cum_variance =  be_eigenval(k) + cum_variance
      write(unit=stdout,fmt='(A,I4,e14.4,f10.4)') &
         trim(name), k, be_eigenval(k), cum_variance / tot_variance
   end do
   write(unit=stdout,fmt=*) ' '
   write(unit=stdout,fmt='(A,e13.5)') 'Total variance = Tr(Bv) = ', tot_variance
   write(unit=stdout,fmt=*) ' '

   !--------------------------------------------------------------------------
   ! [2.0]: Test global eigenvectors:
   !--------------------------------------------------------------------------

   ! [2.1] Print out global eigenvectors:

   write(unit=stdout,fmt='(2A)') 'Domain eigenvectors for ', trim(name)

   write(unit=stdout,fmt='(50i13)')(m, m=1,kz)
   do k = 1, kz      
      write(unit=stdout,fmt='(I3,50e13.5)')k, (be_eigenvec(k,m), m=1,kz)
   end do
   write(unit=stdout,fmt='(A)') " "

   ! [2.2]: Test eigenvector orthogonality: sum_k (e_m(k) e_n(k)) = delta_mn:

   allocate(array_mask(1:kz))
   allocate(matrix2(1:kz,1:kz))
      
   write(unit=stdout,fmt='(2A)') &
      'Eigenvector orthogonality check for ', trim(name)
   write(unit=stdout,fmt='(A)')' Mode     Diagonal         Maximum off-diagonal'
   do k1 = 1, kz
      do k2 = 1, kz
         matrix2(k1,k2) = sum(be_eigenvec(:,k1) * be_eigenvec(:,k2))
      end do
         
      array_mask(1:kz) =.true.
      array_mask(k1) = .false.
      max_off_diag = MAXVAL(ABS(matrix2(k1,:)),mask=array_mask(:))
      write(unit=stdout,fmt='(I4,4x,1pe12.4,10x,1pe12.4)')k1, matrix2(k1,k1), &
         max_off_diag
   end do
   write(unit=stdout,fmt=*) ' '

   ! [2.3] Test eigenvectors completeness - sum_m (e_m(k1) e_m(k2)) = delta_k1k2

   write(unit=stdout,fmt='(2A)') &
      'Eigenvector completeness check for ', trim(name)
   write(unit=stdout,fmt='(A)')' Level    Diagonal         Maximum off-diagonal'
   do k1 = 1, kz
      do k2 = 1, kz
         matrix2(k1,k2) = sum(be_eigenvec(k1,:) * be_eigenvec(k2,:))
      end do
         
      array_mask(1:kz) =.true.
      array_mask(k1) = .false.
      max_off_diag = MAXVAL(ABS(matrix2(k1,:)),mask=array_mask(:))
      write(unit=stdout,fmt='(I4,4x,1pe12.4,10x,1pe12.4)')k1, matrix2(k1,k1), &
         max_off_diag
   end do
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [3.0]: Tidy up:
   !-------------------------------------------------------------------------  

   deallocate(matrix2)
   deallocate(array_mask)

   if (trace_use) call da_trace_exit("da_check_eof_decomposition")
       
end subroutine da_check_eof_decomposition


subroutine da_transform_vtovv(grid, cv_size, be, cv, vv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain),  intent(inout) :: grid
   integer,       intent(in)    :: cv_size ! Size of cv array.
   type(be_type), intent(in)    :: be   ! Background error structure.
   real,          intent(in)    :: cv(cv_size)   ! control variables.
   type(vp_type), intent(inout) :: vv   ! Grid point/EOF equivalent.

   integer :: s(4)   ! Index bounds into arrays.
   integer :: n      ! Loop counter.
   integer :: mz     ! Vertical truncations.
   integer :: ne     ! Ensemble size.

   logical :: scaling

   scaling = .false.

   if (trace_use) call da_trace_entry("da_transform_vtovv")

   if( use_rf )then
      !-------------------------------------------------------------------------
      ! [1.0] Fill vv arrays from 1-dimensional cv array.
      !-------------------------------------------------------------------------
      call da_cv_to_vv(cv_size, cv, (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /), vv)
   endif					! use wavelets:
   if( .not. use_rf .or. do_normalize ) s(1:2)=1

   !-------------------------------------------------------------------------
   ! [2.0] Perform VToVV Transform:
   !-------------------------------------------------------------------------

   ! [2.1] Transform 1st control variable:

   mz = be % v1 % mz
   s(3)=s(1)+mz-1
   if( use_rf .and. mz > 0 .and. len_scaling1(1) /= 0.0) then
      call da_transform_through_rf(grid, mz, be % v1 % rf_alpha, be % v1 % val, vv % v1)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v1)
      s(2)=s(4)+1
   else
!      print'(a,": be%v1%mz=",I0)',"da_transform_vtovv.inc",mz
   endif
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v1)
   s(1)=s(3)+1

   ! [2.2] Transform 2nd control variable:

   mz = be % v2 % mz
   s(3)=s(1)+mz-1
   if( use_rf .and. mz > 0 .and. len_scaling2(1) /= 0.0) then
      call da_transform_through_rf(grid, mz, be % v2 % rf_alpha, be % v2 % val, vv % v2)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v2)
      s(2)=s(4)+1
   else
!      print'(a,": be%v2%mz=",I0)',"da_transform_vtovv.inc",mz
   endif
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v2)
   s(1)=s(3)+1

   ! [2.3] Transform 3rd control variable

   mz = be % v3 % mz
   s(3)=s(1)+mz-1
   if( use_rf .and. mz > 0 .and. len_scaling3(1) /= 0.0) then
      call da_transform_through_rf(grid, mz, be % v3 % rf_alpha,be % v3 % val, vv % v3)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v3)
      s(2)=s(4)+1
   else
!      print'(a,": be%v3%mz=",I0)',"da_transform_vtovv.inc",mz
   endif
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v3)
   s(1)=s(3)+1

   ! [2.4] Transform 4th control variable
      
   mz = be % v4 % mz
   s(3)=s(1)+mz-1
   if( use_rf .and. mz > 0 .and. len_scaling4(1) /= 0.0) then
      call da_transform_through_rf(grid, mz, be % v4 % rf_alpha, be % v4 % val, vv % v4)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v4)
      s(2)=s(4)+1
   else
!      print'(a,": be%v4%mz=",I0)',"da_transform_vtovv.inc",mz
   endif
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v4)
   s(1)=s(3)+1

   ! [2.5] Transform 5th control variable

   mz = be % v5 % mz
   s(3)=s(1)+mz-1
   if( use_rf .and. mz > 0 .and. len_scaling5(1) /= 0.0) then
      call da_transform_through_rf(grid, mz, be % v5 % rf_alpha, be % v5 % val, vv % v5)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v5)
      s(2)=s(4)+1
   else
!      print'(a,": be%v5%mz=",I0)',"da_transform_vtovv.inc",mz
   endif
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v5)
   s(1)=s(3)+1


   ! [2.12] Transform alpha control variable

   ne = be % ne
   if (ne > 0) then
      mz = be % alpha % mz
      if( use_rf )then
         do n = 1, ne
            if ( anal_type_hybrid_dual_res ) then
               call da_transform_through_rf_dual_res(grid%intermediate_grid, mz, be % alpha % rf_alpha, & 
                             be % alpha % val, vv % alpha(:,:,:,n) )
            else
               call da_transform_through_rf(grid, mz, be % alpha % rf_alpha, be % alpha % val, vv % alpha(:,:,:,n) )
            endif
         end do
      else
         do n = 1, ne
            s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
            call da_transform_through_wavelet(grid,mz,be%alpha%wsd,cv(s(2):s(4)),vv%alpha(:,:,:,n))
            s(2)=s(4)+1
         end do
      endif
      if( do_normalize )then
         do n = 1, ne
            call da_transform_rescale(mz,be%alpha%sd,vv%alpha(:,:,:,n))
         end do
      endif
   endif

   if (trace_use) call da_trace_exit("da_transform_vtovv")

endsubroutine da_transform_vtovv
subroutine da_transform_vtovv_adj(grid, cv_size, be, cv, vv)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type(domain),  intent(inout) :: grid
   integer,       intent(in)    :: cv_size ! Size of cv array.
   type(be_type), intent(in)    :: be   ! Background error structure.
   real,          intent(inout) :: cv(cv_size)   ! control variables.
   type(vp_type), intent(inout) :: vv   ! Grid point/EOF control var.

   integer :: s(4)   ! Index bounds into arrays.
   integer :: n      ! Loop counter.
   integer :: mz     ! Vertical truncation.
   integer :: ne     ! Ensemble size.

   logical :: scaling
 
   if (trace_use) call da_trace_entry("da_transform_vtovv_adj")

   if( .not. use_rf .or. do_normalize ) s(1:2)=1


   !-------------------------------------------------------------------------
   ! [2.0] Perform VToVV Transform:
   !-------------------------------------------------------------------------

   ! [2.1] Transform 1st control variable:
   mz = be % v1 % mz
   s(3)=s(1)+mz-1
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v1)
   if( use_rf .and. mz > 0 .and. len_scaling1(1) /= 0.0) then
      call da_transform_through_rf_adj(grid, mz, be % v1 % rf_alpha, be % v1 % val, vv % v1)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet_adj(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v1)
      s(2)=s(4)+1
   else
!      print'(a,": be%v1%mz=",I0)',"da_transform_vtovv_adj.inc",mz
   endif
   s(1)=s(3)+1

   ! [2.2] Transform 2nd control variable:

   mz = be % v2 % mz
   s(3)=s(1)+mz-1
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v2)
   if( use_rf .and. mz > 0 .and. len_scaling2(1) /= 0.0) then
      call da_transform_through_rf_adj(grid, mz, be % v2 % rf_alpha, be % v2 % val, vv % v2)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet_adj(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v2)
      s(2)=s(4)+1
   else
!      print'(a,": be%v2%mz=",I0)',"da_transform_vtovv_adj.inc",mz
   endif
   s(1)=s(3)+1

   ! [2.3] Transform 3rd control variable

   mz = be % v3 % mz
   s(3)=s(1)+mz-1
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v3)
   if( use_rf .and. mz > 0 .and. len_scaling3(1) /= 0.0) then
      call da_transform_through_rf_adj(grid, mz, be % v3 % rf_alpha, be % v3 % val, vv % v3)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet_adj(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v3)
      s(2)=s(4)+1
   else
!      print'(a,": be%v3%mz=",I0)',"da_transform_vtovv_adj.inc",mz
   endif
   s(1)=s(3)+1
   
   ! [2.4] Transform 4th control variable
      
   mz = be % v4 % mz
   s(3)=s(1)+mz-1
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v4)
   if( use_rf .and. mz > 0 .and. len_scaling4(1) /= 0.0) then
      call da_transform_through_rf_adj(grid, mz, be % v4 % rf_alpha, be % v4 % val, vv % v4)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet_adj(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v4)
      s(2)=s(4)+1
   else
!      print'(a,": be%v4%mz=",I0)',"da_transform_vtovv_adj.inc",mz
   endif
   s(1)=s(3)+1

   ! [2.5] Transform 5th control variable

   mz = be % v5 % mz
   s(3)=s(1)+mz-1
   if( do_normalize )call da_transform_rescale(mz,be%sd(:,:,s(1):s(3)),vv%v5)
   if( use_rf .and. mz > 0 .and. len_scaling5(1) /= 0.0) then
      call da_transform_through_rf_adj(grid, mz, be % v5 % rf_alpha, be % v5 % val, vv % v5)
   elseif( mz > 0 ) then
      s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
      call da_transform_through_wavelet_adj(grid,mz,be%wsd(:,:,s(1):s(3)),cv(s(2):s(4)),vv%v5)
      s(2)=s(4)+1
   else
!      print'(a,": be%v5%mz=",I0)',"da_transform_vtovv_adj.inc",mz
   endif
   s(1)=s(3)+1


   ! [2.12] Transform alpha control variable

   ne = be % ne
   if (ne > 0) then
      mz = be % alpha % mz
      if( do_normalize )then
         do n = 1, ne
            call da_transform_rescale(mz,be%alpha%sd,vv%alpha(:,:,:,n))
         end do
      endif
      if( use_rf )then
         do n = 1, ne
            if ( anal_type_hybrid_dual_res ) then
               call da_transform_through_rf_adj_dual_res(grid % intermediate_grid, mz, be % alpha % rf_alpha, &
                                                         be % alpha % val, vv % alpha(:,:,:,n))
            else
               call da_transform_through_rf_adj(grid, mz, be % alpha % rf_alpha, be % alpha % val, vv % alpha(:,:,:,n))
            endif
         end do
      else
         do n = 1, ne
            s(4)=s(2)+nij(0,0,2)*nij(0,1,2)*mz-1
            call da_transform_through_wavelet_adj(grid,mz,be%alpha%wsd,cv(s(2):s(4)),vv%alpha(:,:,:,n))
            s(2)=s(4)+1
         end do
      endif
   endif

   if( use_rf )then
      !-------------------------------------------------------------------------
      ! [1.0] Fill 1D cv array from 3-dimensional vv arrays.
      !-------------------------------------------------------------------------
      call da_vv_to_cv( vv, grid%xp, (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /), cv_size, cv)
   endif

   if (trace_use) call da_trace_exit("da_transform_vtovv_adj")

endsubroutine da_transform_vtovv_adj
 subroutine da_transform_rescale(mz, scaling_factor, field)

 !----------------------------------------------------------------------
 ! Purpose: Control-variable transform through rescaling
 ! Author: Yann Michel 2009/9/29, Aime' Fournier 2011/1/2
 !----------------------------------------------------------------------

 implicit none
 integer,          intent(in)    :: mz
 real,             intent(in)    :: scaling_factor(:,:,:)        ! Scaling factor
 real,             intent(inout) :: field(:,:,:)                 ! Field to be transformed. 
 integer :: kt                                                   ! Loop counters.

 if (trace_use_dull) call da_trace_entry("da_transform_rescale")

 do kt = max(1,kts), min(mz,kte)
   field(its:ite,jts:jte,kt)=field(its:ite,jts:jte,kt)*scaling_factor(its:ite,jts:jte,kt)
 enddo
  
 if (trace_use_dull) call da_trace_exit("da_transform_rescale")

endsubroutine da_transform_rescale
subroutine da_transform_vtox(grid, cv_size, xbx, be, ep, cv, vv, vp)

   !--------------------------------------------------------------------------
   ! Purpose: Control variable transform x' = Uv. 
   !--------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   integer,        intent(in)    :: cv_size ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx  ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be   ! background errors.
   type(ep_type),  intent(in)    :: ep   ! Ensemble perturbations.
   real,           intent(in)    :: cv(1:cv_size)   ! control variables.
   type(vp_type),  intent(out)   :: vv   ! grdipt/eof cv (local).
   type(vp_type),  intent(inout) :: vp   ! grdipt/level cv (local).

   if (trace_use) call da_trace_entry("da_transform_vtox")

   call da_zero_x (grid%xa)
   
   if (.not. use_background_errors) then
      if (trace_use) call da_trace_exit("da_transform_vtox")
      return
   end if

   !----------------------------------------------------------------------
   ! [1.0]: Perform vv = u_h cv transform:
   !----------------------------------------------------------------------

   if (global) then
      call da_transform_vtovv_global(cv_size, xbx, be, cv, vv)
   else if ( (fg_format == fg_format_wrf_arw_regional .or.   &
              fg_format == fg_format_wrf_nmm_regional) .and. &
              (.not. cv_options == 3) )then
      call da_transform_vtovv(grid, cv_size, be, cv, vv)
   end if
   
   !----------------------------------------------------------------------
   ! [2.0]: Perform vp = u_v vv transform:
   !----------------------------------------------------------------------

   if (  cv_options == 3 ) then

      call da_apply_be( be, cv, vp, grid)
      call da_transform_bal( vp, be, grid)
   
   else

   if (vert_corr == vert_corr_2) then      
      call da_vertical_transform(grid, 'u', be, grid%xb % vertical_inner_product, vv, vp)
   else
      vp % v1(its:ite,jts:jte,kts:kte) = vv % v1(its:ite,jts:jte,kts:kte)
      vp % v2(its:ite,jts:jte,kts:kte) = vv % v2(its:ite,jts:jte,kts:kte)
      vp % v3(its:ite,jts:jte,kts:kte) = vv % v3(its:ite,jts:jte,kts:kte)
      vp % v4(its:ite,jts:jte,kts:kte) = vv % v4(its:ite,jts:jte,kts:kte)
      vp % v5(its:ite,jts:jte,kts:kte) = vv % v5(its:ite,jts:jte,kts:kte)
      if (be % ne > 0) then
!        vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne) = vv%alpha(its:ite,jts:jte,kts:kte,1:be%ne)
         vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) =  &
             vv%alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)
      end if
   end if

   !----------------------------------------------------------------------  
   ! [3.0]: Perform x = u_p vp transform::
   !----------------------------------------------------------------------

   call da_transform_vptox(grid, vp, be, ep)

   end if
   
   if (trace_use) call da_trace_exit("da_transform_vtox")

end subroutine da_transform_vtox

subroutine da_transform_xtoxa(grid)

   !--------------------------------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF TL fields to dignostic fields needed by observational operators
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !--------------------------------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid

   integer :: i, j, k
   real    :: sdmd, s1md, mu
   real    :: p(kms:kme)
   real    :: mr_a(kms:kme)
   real    :: mr_b(kms:kme)
   real    :: PU, PD, coeff

   if (trace_use) call da_trace_entry("da_transform_xtoxa")

   !----------------------------------------------------------------------
   ! [4.0]: Move the following:
   !----------------------------------------------------------------------

   if ( .not. var4d ) then ! for 4dvar the xa%p is come from the model

   do j=jts,jte
      do i=its,ite
        if ((fg_format==fg_format_wrf_arw_regional) .or. &
            (fg_format==fg_format_wrf_arw_global  ) ) then

            sdmd=0.0
            s1md=0.0
            do k=kts,kte
               mr_a(k) = grid%xa%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))**2
               mr_b(k) = grid%xb%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))

               sdmd=sdmd+mr_a(k)*grid%xb%dnw(k)
               s1md=s1md+(1.0+mr_b(k))*grid%xb%dnw(k)
            end do

            mu=-(grid%xa%psfc(i,j)+grid%xb%psac(i,j)*sdmd)/s1md

            p(kte+1)=0.0

            do k=kte,kts,-1
               p(k)=p(k+1)-(mu*(1.0+mr_b(k)) &
                       + grid%xb%psac(i,j)*mr_a(k))*grid%xb%dnw(k)

               grid%xa%p(i,j,k)=0.5*(p(k)+p(k+1))
            end do
         else if (fg_format == fg_format_kma_global) then
            do k=kts,kte
               if (k == kte) then
                  coeff=grid%xb%KMA_B(K)/(grid%xb%KMA_A(K)+grid%xb%KMA_B(K)*grid%xb%psfc(I,J)/100.0)
               else
                  PU = grid%xb%KMA_A(K+1) + grid%xb%KMA_B(K+1)*grid%xb%psfc(I,J)/100.0
                  PD = grid%xb%KMA_A(K ) + grid%xb%KMA_B(K )*grid%xb%psfc(I,J)/100.0
                  coeff=grid%xb%KMA_B(K)  *1.0/(PD-PU)**2*(-PU*(LOG(PD)-LOG(PU)) &
                    + PD-PU)&
                    + grid%xb%KMA_B(K+1)*1.0/(PD-PU)**2*(PD*(LOG(PD)-LOG(PU))-PD+PU)
               end if
               ! Here since grid%xa%psfc holds value in Pa. dlnp -> dp
               grid%xa%p(i,j,k) =  grid%xb%p(i,j,k) * grid%xa%psfc(I,J)/100.0 * coeff

            end do
         end if
      end do
   end do

   endif ! only for 3dvar

   call da_pt_to_rho_lin(grid)
 
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE


   ! If test_transforms = .true., not "XToY" transform needed to do here:

   if (.not.test_transforms) then
      ! Exchange grid%xa halo region.

      if (sfc_assi_options == 2) then
         call da_transform_xtowtq (grid)
         ! Exchange grid%xa (surface variable) halo region.
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_SFC_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_SFC_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
      end if

      if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. use_gpsztdobs .or. &
          use_ssmitbobs .or. use_ssmiretrievalobs .or. use_gpsrefobs) then

         ! Now do something for PW
         call da_transform_xtotpw(grid)

         ! Space-based GPS Refractivity and Ground-based GPS ZTD: 
         if ( use_gpsrefObs .or. use_gpsztdObs ) then  
            call da_transform_xtogpsref_lin(grid)  
            if (use_GpsztdObs) call da_transform_xtoztd_lin(grid) 
         endif 

         if (use_ssmt1obs .or. use_ssmt2obs .or. &
              use_ssmitbobs .or. use_ssmiretrievalobs) then
            if (global) then
               call da_error("da_transform_xtoxa.inc",135, &
                  (/"grid%xb%speed is not available, see da_transfer_kmatoxb.inc"/))
            end if
            call da_transform_xtoseasfcwind_lin(grid)
         end if
         if (use_ssmitbobs) call da_transform_xtotb_lin (grid)

         ! Exchange grid%xa halo region.
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_SSMI_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_SSMI_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
      end if
   end if

   ! Compute w increments using Richardson's eqn.

   if ( Use_RadarObs ) then
      if ( .not. var4d .and. cloud_cv_options == 1) call da_uvprho_to_w_lin(grid)

      do k=kts,kte
         do j=jts,jte
            do i=its,ite
               grid%xa%wh(i,j,k)=0.5*(grid%xa%w(i,j,k)+grid%xa%w(i,j,k+1))
            end do
         end do
      end do

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_RADAR_XA_W.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_RADAR_XA_W_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
   end if

   if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) ) then
      if ( cloud_cv_options == 1 .and. use_3dvar_phy)then
         ! Partition of hydrometeor increments via warm rain process
         call da_moist_phys_lin(grid)
      end if
   end if

   !---------------------------------------------------------------
   ! Polar treatment for Global 
   !---------------------------------------------------------------

   if (global)  then   
      call da_get_vpoles(grid%xa%u,grid%xa%v, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%t, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%p, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%q, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%psfc, &
         ids, ide, jds, jde, ims, ime, jms, jme,   1,   1, its, ite, jts, jte,   1,   1)
      call da_set_boundary_xa(grid)
   end if   

   if (trace_use) call da_trace_exit("da_transform_xtoxa")

end subroutine da_transform_xtoxa


subroutine da_transform_xtoxa_all(grid)

   !--------------------------------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF TL fields to dignostic fields needed by observational operators
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !--------------------------------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid

   integer :: i, j, k
   real    :: sdmd, s1md, mu
   real    :: p(kms:kme),geoh(kms:kme)
   real    :: mr_a(kms:kme)
   real    :: mr_b(kms:kme)
   real    :: PU, PD, coeff

   if (trace_use) call da_trace_entry("da_transform_xtoxa_all")

   !----------------------------------------------------------------------
   ! [4.0]: Move the following:
   !----------------------------------------------------------------------

   if ( .not. var4d ) then ! for 4dvar the xa%p is come from the model

   do j=jts,jte
      do i=its,ite
        if ((fg_format==fg_format_wrf_arw_regional) .or. &
            (fg_format==fg_format_wrf_arw_global  ) ) then

            sdmd=0.0
            s1md=0.0
            do k=kts,kte
               mr_a(k) = grid%xa%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))**2
               mr_b(k) = grid%xb%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))

               sdmd=sdmd+mr_a(k)*grid%xb%dnw(k)
               s1md=s1md+(1.0+mr_b(k))*grid%xb%dnw(k)
            end do

            grid%xa%mu(i,j)=-(grid%xa%psfc(i,j)+grid%xb%psac(i,j)*sdmd)/s1md

            p(kte+1)=0.0

            do k=kte,kts,-1
               p(k)=p(k+1)-(grid%xa%mu(i,j)*(1.0+mr_b(k)) &
                       + grid%xb%psac(i,j)*mr_a(k))*grid%xb%dnw(k)

               grid%xa%p(i,j,k)=0.5*(p(k)+p(k+1))
            end do
         else if (fg_format == fg_format_kma_global) then
            do k=kts,kte
               if (k == kte) then
                  coeff=grid%xb%KMA_B(K)/(grid%xb%KMA_A(K)+grid%xb%KMA_B(K)*grid%xb%psfc(I,J)/100.0)
               else
                  PU = grid%xb%KMA_A(K+1) + grid%xb%KMA_B(K+1)*grid%xb%psfc(I,J)/100.0
                  PD = grid%xb%KMA_A(K ) + grid%xb%KMA_B(K )*grid%xb%psfc(I,J)/100.0
                  coeff=grid%xb%KMA_B(K)  *1.0/(PD-PU)**2*(-PU*(LOG(PD)-LOG(PU)) &
                    + PD-PU)&
                    + grid%xb%KMA_B(K+1)*1.0/(PD-PU)**2*(PD*(LOG(PD)-LOG(PU))-PD+PU)
               end if
               ! Here since grid%xa%psfc holds value in Pa. dlnp -> dp
               grid%xa%p(i,j,k) =  grid%xb%p(i,j,k) * grid%xa%psfc(I,J)/100.0 * coeff

            end do
         end if
      end do
   end do

   call da_pt_to_rho_lin(grid)

!by lixin update perturbated geoh

   do j=jts,jte
      do i=its,ite
         if ((fg_format==fg_format_wrf_arw_regional) .or. &
            (fg_format==fg_format_wrf_arw_global  ) ) then
            geoh(kts)=0.0
            do k=kts,kte
               geoh(k+1)=geoh(k)+(-grid%xa%mu(i,j) &
                       + grid%xa%rho(i,j,k)*grid%xb%psac(i,j)/grid%xb%rho(i,j,k))*grid%xb%dnw(k)/grid%xb%rho(i,j,k)
               grid%xa%geoh(i,j,k)=0.5*(geoh(k)+geoh(k+1))
            end do
         end if
      end do
   end do

   endif ! only for 3dvar
 
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XA_ALL.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XA_ALL_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XB_ALL.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XB_ALL_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE


   ! If test_transforms = .true., not "XToY" transform needed to do here:

   if (.not.test_transforms) then
      ! Exchange grid%xa halo region.

      if (sfc_assi_options == 2) then
         call da_transform_xtowtq (grid)
         ! Exchange grid%xa (surface variable) halo region.
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_SFC_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_SFC_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
      end if

      if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. use_gpsztdobs .or. &
          use_ssmitbobs .or. use_ssmiretrievalobs .or. use_gpsrefobs) then

         ! Now do something for PW
         call da_transform_xtotpw(grid)

         ! Space-based GPS Refractivity and Ground-based GPS ZTD: 
         if ( use_gpsrefObs .or. use_gpsztdObs ) then  
            call da_transform_xtogpsref_lin(grid)  
            if (use_GpsztdObs) call da_transform_xtoztd_lin(grid) 
         endif 

         if (use_ssmt1obs .or. use_ssmt2obs .or. &
              use_ssmitbobs .or. use_ssmiretrievalobs) then
            if (global) then
               call da_error("da_transform_xtoxa_all.inc",152, &
                  (/"grid%xb%speed is not available, see da_transfer_kmatoxb.inc"/))
            end if
            call da_transform_xtoseasfcwind_lin(grid)
         end if
         if (use_ssmitbobs) call da_transform_xtotb_lin (grid)

         ! Exchange grid%xa halo region.
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_SSMI_XA.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_SSMI_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
      end if
   end if

   ! Compute w increments using Richardson's eqn.

   if ( Use_RadarObs ) then
      if ( .not. var4d) call da_uvprho_to_w_lin(grid)

      do k=kts,kte
         do j=jts,jte
            do i=its,ite
               grid%xa%wh(i,j,k)=0.5*(grid%xa%w(i,j,k)+grid%xa%w(i,j,k+1))
            end do
         end do
      end do

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_RADAR_XA_W.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_RADAR_XA_W_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
   end if

   if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) ) then
      if ( cloud_cv_options == 1 )then
         ! Partition of hydrometeor increments via warm rain process
         call da_moist_phys_lin(grid)
      end if
   end if

   !---------------------------------------------------------------
   ! Polar treatment for Global 
   !---------------------------------------------------------------

   if (global)  then   
      call da_get_vpoles(grid%xa%u,grid%xa%v, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%t, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%p, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%q, &
         ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call da_get_spoles(grid%xa%psfc, &
         ids, ide, jds, jde, ims, ime, jms, jme,   1,   1, its, ite, jts, jte,   1,   1)
      call da_set_boundary_xa(grid)
   end if   

   if (trace_use) call da_trace_exit("da_transform_xtoxa_all")

end subroutine da_transform_xtoxa_all


subroutine da_transform_vtox_adj(grid, cv_size, xbx, be, ep, vp, vv, cv)

   !--------------------------------------------------------------------------
   ! Purpose: Control variable transform Adjoint
   !--------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   integer,        intent(in)    :: cv_size ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx  ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be   ! background errors.
   type(ep_type),  intent(in)    :: ep   ! ensemble perturbation structure.
   type(vp_type),  intent(out)   :: vp   ! grdipt/level cv (local).
   type(vp_type),  intent(out)   :: vv   ! grdipt/eof cv (local).
   real,           intent(inout) :: cv(1:cv_size) ! control variables.


   if (.not. use_background_errors) return

   if (trace_use) call da_trace_entry("da_transform_vtox_adj")

   if (  cv_options == 3 ) then

      call da_transform_bal_adj( vp, be, grid)

      call da_apply_be_adj( be, cv, vp, grid)

   else
  
   !-------------------------------------------------------------------------
   ! [3.0]: Perform x = u_p vp transform::
   !-------------------------------------------------------------------------

   call da_zero_vp_type (vp)
   call da_transform_vptox_adj(grid, vp, be, ep)


   !-------------------------------------------------------------------------
   ! [2.0]: Perform vp = u_v vv transform:
   !-------------------------------------------------------------------------
   call da_zero_vp_type (vv)

   if (vert_corr == 2) then      
      call da_vertical_transform(grid, 'u_adj', be, &
         grid%xb % vertical_inner_product, vv, vp)
   else
      vv % v1(its:ite,jts:jte,kts:kte) = vp % v1(its:ite,jts:jte,kts:kte)
      vv % v2(its:ite,jts:jte,kts:kte) = vp % v2(its:ite,jts:jte,kts:kte)
      vv % v3(its:ite,jts:jte,kts:kte) = vp % v3(its:ite,jts:jte,kts:kte)
      vv % v4(its:ite,jts:jte,kts:kte) = vp % v4(its:ite,jts:jte,kts:kte)
      vv % v5(its:ite,jts:jte,kts:kte) = vp % v5(its:ite,jts:jte,kts:kte)
   !  Uv for alpha is an identity transform:
      if (be % ne > 0) then
!        vv % alpha(its:ite,jts:jte,kts:kte,1:be%ne) = &
!        vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne) 
         vv % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne) = & 
         vp % alpha(its_int:ite_int,jts_int:jte_int,kts_int:kte_int,1:be%ne)
      end if
   end if
   
   end if

   !-------------------------------------------------------------------------
   ! [1.0]: perform vv = u_h cv transform:
   !-------------------------------------------------------------------------

   if (global) then
      call da_transform_vtovv_global_adj(cv_size, xbx, be, cv, vv)
   else if ( (fg_format == fg_format_wrf_arw_regional .or. &
              fg_format == fg_format_wrf_nmm_regional) .and. &
              (.not. cv_options == 3) ) then
      call da_transform_vtovv_adj(grid, cv_size, be, cv, vv)
   end if

   if (trace_use) call da_trace_exit("da_transform_vtox_adj")

end subroutine da_transform_vtox_adj


subroutine da_transform_xtoxa_adj(grid)

   !--------------------------------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF TL fields to dignostic fields needed by observational operators
   !          Adjoint   
   !--------------------------------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid

   integer :: i, j, k
   real    :: sdmd, s1md, mu
   real    :: p(kms:kme), mr_a(kms:kme), mr_b(kms:kme)
   real    :: PU, PD, coeff

   if (trace_use) call da_trace_entry("da_transform_xtoxa_adj")

   !-------------------------------------------------------------------------
   ! Polar treatment for Global
   !-------------------------------------------------------------------------

   if (global) then     
      ! Poles treatment for global WRFVAR
      call da_get_avpoles(grid%xa%u,grid%xa%v, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%t, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%p, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%q, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%psfc, &
         ids,ide,jds,jde,ims,ime,jms,jme,1,1,its,ite,jts,jte,1,1)
   end if     

   ! Compute w increments using Richardson's eqn.
   if ( Use_RadarObs)  then
      do k=kts,kte
         do j=jts,jte
            do i=its,ite
               grid%xa%w(i,j,k)=grid%xa%w(i,j,k)+0.5*grid%xa%wh(i,j,k)
               grid%xa%w(i,j,k+1)=grid%xa%w(i,j,k+1)+0.5*grid%xa%wh(i,j,k)
               grid%xa%wh(i,j,k)=0.0
            end do
         end do
      end do
      if ( .not. var4d .and. cloud_cv_options == 1 ) call da_uvprho_to_w_adj(grid)
   end if

   if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) ) then
      if ( cloud_cv_options == 1 .and. use_3dvar_phy) then
         ! Partition of hydrometeor increments via warm rain process
         call da_moist_phys_adj(grid)
      end if
   end if

   !-------------------------------------------------------------------------
   ! If test_transforms = .true., not "XToY" transform needed to do here (YRG):

   if (.not.test_transforms) then
      if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. use_gpsztdobs .or. &
          use_ssmitbobs .or. use_ssmiretrievalobs .or. use_gpsrefobs) then
         if (use_ssmitbobs) call da_transform_xtotb_adj(grid)
         if (use_ssmt1obs .or. use_ssmt2obs .or. &
             use_ssmitbobs .or. use_ssmiretrievalobs) then
            if (global) then
               call da_error("da_transform_xtoxa_adj.inc",68, &
                  (/"grid%xb%speed is not available, see da_transfer_kmatoxb.inc"/))
            end if
            call da_transform_xtoseasfcwind_adj(grid)
         end if

         ! GPS Refractivity:
         if ( use_gpsrefObs .or. use_gpsztdObs) then 
            if (use_gpsztdObs) call da_transform_xtoztd_adj(grid)  
            call da_transform_XTogpsref_adj(grid)  
         endif 

         ! Now for PW.
         call da_transform_xtotpw_adj(grid)
      end if

      if (sfc_assi_options == 2) call da_transform_xtowtq_adj(grid)
   end if
 
   call da_pt_to_rho_adj(grid)

   if ( .not. var4d ) then

   do j=jts,jte
      do i=its,ite
         if ((fg_format==fg_format_wrf_arw_regional) .or. &
             (fg_format==fg_format_wrf_arw_global  )  ) then
            mu=0.0
            s1md=0.0

            p(:)=0.0

            do k=kts,kte
               mr_b(k) = grid%xb%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))
               s1md=s1md+(1.0+mr_b(k))*grid%xb%dnw(k)

               p(k) = p(k) + 0.5*grid%xa%p(i,j,k)
               p(k+1) = p(k+1) + 0.5*grid%xa%p(i,j,k)

               mu = mu - p(k)*(1.0+mr_b(k))*grid%xb%dnw(k)

               mr_a(k) = - p(k)*grid%xb%psac(i,j)*grid%xb%dnw(k)
               p(k+1) = p(k+1) + p(k)
            end do

            grid%xa%psfc(i,j) = grid%xa%psfc(i,j) - mu/s1md
            sdmd=-mu*grid%xb%psac(i,j)/s1md

            do k=kts,kte
               mr_a(k) = mr_a(k) + sdmd*grid%xb%dnw(k)
               grid%xa%q(i,j,k) = grid%xa%q(i,j,k) + mr_a(k)/(1.0 - grid%xb%q(i,j,k))**2
            end do
         else if (fg_format == fg_format_kma_global)then
            do k=kts,kte
               if (k == kte) then
                  coeff = grid%xb%KMA_B(K)/(grid%xb%KMA_A(K)+grid%xb%KMA_B(K)* &
                     grid%xb%psfc(I,J)/100.0)
               else
                  PU = grid%xb%KMA_A(K+1) + grid%xb%KMA_B(K+1)*grid%xb%psfc(I,J)/100.0
                  PD = grid%xb%KMA_A(K ) + grid%xb%KMA_B(K )*grid%xb%psfc(I,J)/100.0
                  coeff=grid%xb%KMA_B(K)*1.0/(PD-PU)**2*(-PU*(LOG(PD)-LOG(PU))+PD-PU)&
                     + grid%xb%KMA_B(K+1)*1.0/(PD-PU)**2*(PD*(LOG(PD)-LOG(PU))-PD+PU)
               end if
      
               grid%xa%psfc(i,j) = grid%xa % psfc(i,j) + &
                  grid%xb%p(i,j,k) * grid%xa % p(i,j,k)/100.0 * coeff 
            end do
         end if
      end do
   end do

   endif ! only for 3dvar

   if (global) then     
      call da_set_boundary_xa(grid)
   end if

   if (trace_use) call da_trace_exit("da_transform_xtoxa_adj")

end subroutine da_transform_xtoxa_adj


subroutine da_transform_xtoxa_adj_all(grid)

   !--------------------------------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF TL fields to dignostic fields needed by observational operators
   !          Adjoint   
   !--------------------------------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid

   integer :: i, j, k
   real    :: sdmd, s1md, mu
   real    :: p(kms:kme),geoh(kms:kme), mr_a(kms:kme), mr_b(kms:kme)
   real    :: PU, PD, coeff

   if (trace_use) call da_trace_entry("da_transform_xtoxa_adj_all")

   !-------------------------------------------------------------------------
   ! Polar treatment for Global
   !-------------------------------------------------------------------------

   if (global) then     
      ! Poles treatment for global WRFVAR
      call da_get_avpoles(grid%xa%u,grid%xa%v, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%t, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%p, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%q, &
         ids,ide,jds,jde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      call da_get_aspoles(grid%xa%psfc, &
         ids,ide,jds,jde,ims,ime,jms,jme,1,1,its,ite,jts,jte,1,1)
   end if     

   ! Compute w increments using Richardson's eqn.
   if ( Use_RadarObs)  then
      do k=kts,kte
         do j=jts,jte
            do i=its,ite
               grid%xa%w(i,j,k)=grid%xa%w(i,j,k)+0.5*grid%xa%wh(i,j,k)
               grid%xa%w(i,j,k+1)=grid%xa%w(i,j,k+1)+0.5*grid%xa%wh(i,j,k)
               grid%xa%wh(i,j,k)=0.0
            end do
         end do
      end do

      if ( .not. var4d ) call da_uvprho_to_w_adj(grid)
   end if

   if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) ) then
     if ( cloud_cv_options == 1) then
      ! Partition of hydrometeor increments via warm rain process
      call da_moist_phys_adj(grid)
     end if
   end if

   !-------------------------------------------------------------------------
   ! If test_transforms = .true., not "XToY" transform needed to do here (YRG):

   if (.not.test_transforms) then
      if (use_ssmt1obs .or. use_ssmt2obs .or. use_gpspwobs .or. use_gpsztdobs .or. &
          use_ssmitbobs .or. use_ssmiretrievalobs .or. use_gpsrefobs) then
         if (use_ssmitbobs) call da_transform_xtotb_adj(grid)
         if (use_ssmt1obs .or. use_ssmt2obs .or. &
             use_ssmitbobs .or. use_ssmiretrievalobs) then
            if (global) then
               call da_error("da_transform_xtoxa_adj_all.inc",69, &
                  (/"grid%xb%speed is not available, see da_transfer_kmatoxb.inc"/))
            end if
            call da_transform_xtoseasfcwind_adj(grid)
         end if

         ! GPS Refractivity:
         if ( use_gpsrefObs .or. use_gpsztdObs) then 
            if (use_gpsztdObs) call da_transform_xtoztd_adj(grid)  
            call da_transform_XTogpsref_adj(grid)  
         endif 

         ! Now for PW.
         call da_transform_xtotpw_adj(grid)
      end if

      if (sfc_assi_options == 2) call da_transform_xtowtq_adj(grid)
   end if



   if ( .not. var4d ) then

   grid%xa%mu=0.0
   grid%xa%rho=0.0

   do j=jts,jte
      do i=its,ite
         if ((fg_format==fg_format_wrf_arw_regional) .or. &
             (fg_format==fg_format_wrf_arw_global  )  ) then
            geoh(:)=0.0
            do k=kte,kts,-1
               geoh(k) = geoh(k) + 0.5*grid%xa%geoh(i,j,k)
               geoh(k+1) = geoh(k+1) + 0.5*grid%xa%geoh(i,j,k)
               grid%xa%mu(i,j)=grid%xa%mu(i,j)-geoh(k+1)*grid%xb%dnw(k)/grid%xb%rho(i,j,k)
               grid%xa%rho(i,j,k)=grid%xa%rho(i,j,k)+geoh(k+1)*grid%xb%psac(i,j)*grid%xb%dnw(k)/(grid%xb%rho(i,j,k)**2)
               geoh(k) = geoh(k) + geoh(k+1)
            end do
         end if
      end do
   end do

   call da_pt_to_rho_adj(grid)

   do j=jts,jte
      do i=its,ite
         if ((fg_format==fg_format_wrf_arw_regional) .or. &
             (fg_format==fg_format_wrf_arw_global  )  ) then
            
            s1md=0.0
            p(:)=0.0

            do k=kts,kte
               mr_b(k) = grid%xb%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))
               s1md=s1md+(1.0+mr_b(k))*grid%xb%dnw(k)

               p(k) = p(k) + 0.5*grid%xa%p(i,j,k)
               p(k+1) = p(k+1) + 0.5*grid%xa%p(i,j,k)

               grid%xa%mu(i,j) = grid%xa%mu(i,j) - p(k)*(1.0+mr_b(k))*grid%xb%dnw(k)

               mr_a(k) = - p(k)*grid%xb%psac(i,j)*grid%xb%dnw(k)
               p(k+1) = p(k+1) + p(k)
            end do

            grid%xa%psfc(i,j) = grid%xa%psfc(i,j) - grid%xa%mu(i,j)/s1md
            sdmd=-grid%xa%mu(i,j)*grid%xb%psac(i,j)/s1md

            do k=kts,kte
               mr_a(k) = mr_a(k) + sdmd*grid%xb%dnw(k)
               grid%xa%q(i,j,k) = grid%xa%q(i,j,k) + mr_a(k)/(1.0 - grid%xb%q(i,j,k))**2
            end do
         else if (fg_format == fg_format_kma_global)then
            do k=kts,kte
               if (k == kte) then
                  coeff = grid%xb%KMA_B(K)/(grid%xb%KMA_A(K)+grid%xb%KMA_B(K)* &
                     grid%xb%psfc(I,J)/100.0)
               else
                  PU = grid%xb%KMA_A(K+1) + grid%xb%KMA_B(K+1)*grid%xb%psfc(I,J)/100.0
                  PD = grid%xb%KMA_A(K ) + grid%xb%KMA_B(K )*grid%xb%psfc(I,J)/100.0
                  coeff=grid%xb%KMA_B(K)*1.0/(PD-PU)**2*(-PU*(LOG(PD)-LOG(PU))+PD-PU)&
                     + grid%xb%KMA_B(K+1)*1.0/(PD-PU)**2*(PD*(LOG(PD)-LOG(PU))-PD+PU)
               end if
      
               grid%xa%psfc(i,j) = grid%xa % psfc(i,j) + &
                  grid%xb%p(i,j,k) * grid%xa % p(i,j,k)/100.0 * coeff 
            end do
         end if
      end do
   end do

   endif ! only for 3dvar

   if (global) then     
      call da_set_boundary_xa(grid)
   end if

   if (trace_use) call da_trace_exit("da_transform_xtoxa_adj_all")

end subroutine da_transform_xtoxa_adj_all


subroutine da_transform_vptox(grid, vp, be, ep)

   !-----------------------------------------------------------------------
   ! Purpose: Physical transform of analysis increment variables.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   ! Updates:
   !
   !       Implementation of multi-variate BE for cv_options=6
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !-----------------------------------------------------------------------

   implicit none

   type (domain), intent(inout)         :: grid
   
   type (vp_type), intent(inout)        :: vp  ! CV on grid structure.
   type (be_type), intent(in), optional :: be  ! Background errors.
   type (ep_type), intent(in), optional :: ep  ! Ensemble perturbations.

   integer :: i, k, j, k1, ij            ! Loop counters.
   real, allocatable          :: chi_u(:,:,:)  ! Unbalanced chi

   if (trace_use) call da_trace_entry("da_transform_vptox") 

   !---------------------------------------------------------------------------
   !  [1] Add flow-dependent increments in control variable space (vp):
   !---------------------------------------------------------------------------

   if (be % ne > 0 .and. alphacv_method == alphacv_method_vp) then
      call da_add_flow_dependence_vp(be % ne, ep, vp, its,ite, jts,jte, kts,kte)
   end if

   !--------------------------------------------------------------------------
   ! [2] Impose statistical balance constraints:
   !--------------------------------------------------------------------------

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k1, k, j, i)
   do ij = 1 , grid%num_tiles

   if ( cv_options == 6 ) then
      allocate (chi_u(its:ite,grid%j_start(ij):grid%j_end(ij),kts:kte) )
      do k = kts, kte
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               chi_u(i,j,k) = vp%v2(i,j,k)
            end do
         end do
      end do
   end if

   ! Chi:
   if (cv_options /= 7) then
      do k = kts, kte
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               vp%v2(i,j,k) = vp%v2(i,j,k) + be%reg_psi_chi(j,k)* vp%v1(i,j,k)
            end do
         end do
      end do
   end if
  
   ! Temperature:
   do k = kts, kte
      do j = grid%j_start(ij), grid%j_end(ij)
         do i = its, ite
            grid%xa%t(i,j,k) = vp%v3(i,j,k)
         end do
      end do
   end do

   if (cv_options /= 7) then
      do k1 = kts, kte
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                  grid%xa%t(i,j,k) = grid%xa%t(i,j,k) + be%reg_psi_t(j,k,k1)*vp%v1(i,j,k1)
               end do
            end do
         end do
      end do
   end if

   if ( cv_options == 6 ) then
      do k1 = kts, kte
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                  grid%xa%t(i,j,k) = grid%xa%t(i,j,k) + be%reg_chi_u_t(j,k,k1)*chi_u(i,j,k1)
               end do
            end do
         end do
      end do
   end if

   ! Surface Pressure
   do j = grid%j_start(ij), grid%j_end(ij)
      do i = its, ite
         grid%xa%psfc(i,j) = vp%v5(i,j,1) 
      end do
   end do

   if (cv_options /= 7) then
      do k = kts,kte
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               grid%xa%psfc(i,j) = grid%xa%psfc(i,j) + be%reg_psi_ps(j,k)*vp%v1(i,j,k)
            end do
         end do
      end do
   end if

   if ( cv_options == 6 ) then
      do k = kts,kte
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               grid%xa%psfc(i,j) = grid%xa%psfc(i,j) + be%reg_chi_u_ps(j,k)*chi_u(i,j,k)
            end do
         end do
      end do
   end if

   ! Moisture
   if ( cv_options == 6 ) then
      do k1 = kts, kte
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                  vp%v4(i,j,k1) = vp%v4(i,j,k1) + be%reg_psi_rh(j,k1,k)*vp%v1(i,j,k) + &
                  be%reg_chi_u_rh(j,k1,k)*chi_u(i,j,k) + be%reg_t_u_rh(j,k1,k)*vp%v3(i,j,k)
               end do
            end do
         end do
      end do
!
      do k = kts, kte
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               vp%v4(i,j,k) = vp%v4(i,j,k) + be%reg_ps_u_rh(j,k)*vp%v5(i,j,1)
            end do
         end do
      end do
   end if

!
   if ( cv_options == 6 ) deallocate (chi_u )

   end do
   !$OMP END PARALLEL DO
   !--------------------------------------------------------------------------
   ! [3] Transform to model variable space:
   !--------------------------------------------------------------------------
  

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_PSICHI_UV.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_PSICHI_UV_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE


   ! Psi and chi to u and v:
   if ( cv_options == 5 .or. cv_options == 6 ) then
      call da_psichi_to_uv(vp % v1, vp % v2, grid%xb % coefx, &
           grid%xb % coefy , grid%xa % u, grid%xa % v)
   else if ( cv_options == 7 ) then
      grid%xa%u = vp%v1
      grid%xa%v = vp%v2
   end if

  if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud).or. &
       (use_radarobs .and. use_radar_rhv) .or. (use_radarobs .and. use_radar_rqv) .or. cloud_cv_options .ge. 2 .or. & 
       (grid%pseudo_var(1:1).eq.'q' .and. grid%pseudo_var(2:2).ne.' ') .or.  &
       (grid%pseudo_var(1:1).eq.'Q' .and. grid%pseudo_var(2:2).ne.' ') ) then

     if ( cloud_cv_options == 1 .and. use_3dvar_phy) then
      ! Pseudo RH --> Total water mixing ratio:
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( ij, i, j, k )
      do ij = 1 , grid%num_tiles
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                 grid%xa % qt(i,j,k) = vp%v4(i,j,k) * grid%xb%qs(i,j,k)
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
     end if
     if ( cloud_cv_options .ge. 2 ) then 
      ! Pseudo RH --> Water vapor mixing ratio:
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( ij, i, j, k )
      do ij = 1 , grid%num_tiles
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                  grid%xa % q(i,j,k) =  vp%v4(i,j,k) * grid%xb%qs(i,j,k)
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
     end if
  else  ! no rf or cloud radiance
      ! Pseudo RH --> Water vapor mixing ratio:
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( ij, i, j, k )
      do ij = 1 , grid%num_tiles
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                  grid%xa % q(i,j,k) =  vp%v4(i,j,k) * grid%xb%qs(i,j,k)
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
  end if ! RF or Radiance
   !---------------------------------------------------------------------------
   !  [4] Add flow-dependent increments in model space (grid%xa):
   !---------------------------------------------------------------------------

!  if (be % ne > 0 .and. alphacv_method == alphacv_method_xa) then
!     call da_add_flow_dependence_xa(grid, be % ne, ep, vp)
!  end if
   if (be % ne > 0 .and. alphacv_method == alphacv_method_xa) then
      if ( anal_type_hybrid_dual_res ) then
         call da_add_flow_dependence_xa_dual_res(grid, be % ne, ep, vp)
      else
         call da_add_flow_dependence_xa(grid, be % ne, ep, vp)
      endif
   end if

   if (trace_use) call da_trace_exit("da_transform_vptox") 
 
end subroutine da_transform_vptox

subroutine da_transform_vptox_adj(grid, vp, be, ep)

   !--------------------------------------------------------------------------
   ! Purpose: Adjoint for Physical transform of variables 
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !
   ! Updates:
   !
   !       Implementation of multi-variate BE for cv_options=6
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !--------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout)        :: grid
   type (vp_type), intent(inout)        :: vp  ! CV on grid structure.
   type (ep_type), intent(in)           :: ep  ! Ensemble perturbation.
   type (be_type), intent(in), optional :: be  ! Background errors.

   integer :: i, k, j, ij, k1              ! Loop counters.
   real, allocatable                    :: chi_u(:,:,:)  ! Unbalanced chi

   if (trace_use) call da_trace_entry("da_transform_vptox_adj")

   !---------------------------------------------------------------------------
   !  [4] Add flow-dependent increments in model space (grid%xa):
   !---------------------------------------------------------------------------
      
!  if (be % ne > 0 .and. alphacv_method == alphacv_method_xa) then
!     call da_add_flow_dependence_xa_adj(be % ne, ep, grid%xa, vp)
!  end if
   if (be % ne > 0 .and. alphacv_method == alphacv_method_xa) then
      if ( anal_type_hybrid_dual_res ) then
         call da_add_flow_dependence_xa_adj_dual_res(be % ne, ep, grid, vp)
      else
         call da_add_flow_dependence_xa_adj(be % ne, ep, grid%xa, vp)
      endif
   endif

   !--------------------------------------------------------------------------
   ! [3] Transform to model variable space:
   !--------------------------------------------------------------------------

  if ( (use_radarobs .and. use_radar_rf) .or. (use_rad .and. crtm_cloud) .or. &
       (use_radarobs .and. use_radar_rhv) .or. (use_radarobs .and. use_radar_rqv) .or. cloud_cv_options .ge. 2 .or. &
       (grid%pseudo_var(1:1).eq.'q' .and. grid%pseudo_var(2:2).ne.' ') .or.  &
       (grid%pseudo_var(1:1).eq.'Q' .and. grid%pseudo_var(2:2).ne.' ') ) then

   if ( cloud_cv_options == 1 .and. use_3dvar_phy ) then
      ! Pseudo RH --> Total water mixing ratio:
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( ij )
      do ij = 1 , grid%num_tiles
         do k = kts,kte
            do j = grid%j_start(ij),grid%j_end(ij)
               do i =  its, ite
                  vp%v4(i,j,k)  = vp%v4(i,j,k) + grid%xa%qt(i,j,k) * grid%xb%qs(i,j,k)
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO
   end if

   if ( cloud_cv_options .ge. 2) then 
      ! Pseudo RH --> Water vapor mixing ratio:
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( ij )
      do ij = 1 , grid%num_tiles
         do k = kts,kte
            do j = grid%j_start(ij),grid%j_end(ij)
               do i =  its, ite
                  vp%v4(i,j,k)  = vp%v4(i,j,k) + grid%xa%q(i,j,k) * grid%xb%qs(i,j,k)   
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO
   end if
  else ! no rf or cloud radiance
      ! Pseudo RH --> Water vapor mixing ratio:
      !$OMP PARALLEL DO &
      !$OMP PRIVATE ( ij )

      do ij = 1 , grid%num_tiles
         do k = kts,kte
            do j = grid%j_start(ij),grid%j_end(ij)
               do i =  its, ite
                  vp%v4(i,j,k)  = vp%v4(i,j,k) + grid%xa%q(i,j,k) * grid%xb%qs(i,j,k)
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO
  end if
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_PSICHI_UV_ADJ.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_PSICHI_UV_ADJ_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE



   ! Transform psi and chi to u and v:
   if ( cv_options == 5 .or. cv_options == 6 ) then
      call da_psichi_to_uv_adj(grid%xa % u, grid%xa % v, grid%xb % coefx, grid%xb % coefy, vp % v1, vp % v2)
   else if ( cv_options == 7) then
      vp%v1 = grid%xa%u + vp%v1
      vp%v2 = grid%xa%v + vp%v2
   end if

   !--------------------------------------------------------------------------
   ! [2] Impose statistical balance constraints:
   !--------------------------------------------------------------------------

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k, j, k1, i )
   do ij = 1 , grid%num_tiles

   if ( cv_options == 6 ) then
      allocate (chi_u(its:ite,grid%j_start(ij):grid%j_end(ij),kts:kte) )
      chi_u = 0
   end if

   ! Moisture
!
   if ( cv_options == 6 ) then
      do k = kts, kte
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = its, ite
               vp%v5(i,j,1) = vp%v5(i,j,1) + be%reg_ps_u_rh(j,k) * vp%v4(i,j,k)
            end do
         end do
      end do

      do k1 = kts, kte
         do k = kts, kte
            do j = grid%j_start(ij), grid%j_end(ij)
               do i = its, ite
                  vp%v3(i,j,k1) = vp%v3(i,j,k1) + be%reg_t_u_rh(j,k,k1)*vp%v4(i,j,k)
                  vp%v1(i,j,k1) = vp%v1(i,j,k1) + be%reg_psi_rh(j,k,k1)*vp%v4(i,j,k)
                  chi_u(i,j,k1) = chi_u(i,j,k1) + be%reg_chi_u_rh(j,k,k1)*vp%v4(i,j,k)
               end do
            end do
         end do
      end do
   end if

   ! Surface Pressure
   if ( cv_options == 6 ) then
      do k= kts,kte
         do j= grid%j_start(ij),grid%j_end(ij)
            do i =  its, ite
               chi_u(i,j,k) = chi_u(i,j,k) + be%reg_chi_u_ps(j,k)*grid%xa%psfc(i,j)
            end do
         end do
      end do
   end if
   if (cv_options /= 7) then
      do k= kts,kte
         do j= grid%j_start(ij),grid%j_end(ij)
            do i =  its, ite
               vp%v1(i,j,k) = vp%v1(i,j,k) + be%reg_psi_ps(j,k)*grid%xa%psfc(i,j)
            end do
         end do
      end do
   end if
   do j= grid%j_start(ij),grid%j_end(ij)
      do i =  its, ite
         vp%v5(i,j,1) = vp%v5(i,j,1) + grid%xa%psfc(i,j) 
      end do
   end do

   ! Temperature
   if ( cv_options == 6 ) then
      do k1 = kts,kte
         do k = kts,kte
            do j = grid%j_start(ij),grid%j_end(ij)
               do i =  its, ite
                  chi_u(i,j,k1) = chi_u(i,j,k1) + be%reg_chi_u_t(j,k,k1)*grid%xa%t(i,j,k)
               end do
            end do
         end do
      end do
   end if
   if (cv_options /= 7) then
      do k1 = kts,kte
         do k = kts,kte
            do j = grid%j_start(ij),grid%j_end(ij)
               do i =  its, ite
                  vp%v1(i,j,k1) = vp%v1(i,j,k1) + be%reg_psi_t(j,k,k1)*grid%xa%t(i,j,k)
               end do
            end do
         end do
      end do
   end if
   do k = kts,kte
      do j = grid%j_start(ij),grid%j_end(ij)
         do i =  its, ite
            vp%v3(i,j,k) = vp%v3(i,j,k) + grid%xa%t(i,j,k)
         end do
      end do
   end do

   ! Chi
   if (cv_options /= 7) then
      do k = kts,kte
         do j = grid%j_start(ij),grid%j_end(ij)
            do i =  its, ite
               vp%v1(i,j,k) = vp%v1(i,j,k) + be%reg_psi_chi(j,k)*vp%v2(i,j,k)
            end do
          end do
      end do
   end if
   if ( cv_options == 6 ) then
      do k = kts,kte
         do j = grid%j_start(ij),grid%j_end(ij)
            do i =  its, ite
               vp%v2(i,j,k) = vp%v2(i,j,k) + chi_u(i,j,k)
            end do
          end do
      end do
      deallocate (chi_u)
   end if

   enddo
   !$OMP END PARALLEL DO

   !---------------------------------------------------------------------------
   !  [1] Add flow-dependent increments in control variable space (vp):
   !---------------------------------------------------------------------------
   
   if (be % ne > 0 .and. alphacv_method == alphacv_method_vp) then
      call da_add_flow_dependence_vp_adj(be % ne, ep, vp)
   end if

   if (trace_use) call da_trace_exit("da_transform_vptox_adj")

end subroutine da_transform_vptox_adj


subroutine da_transform_vvtovp(grid, evec, eval, vertical_wgt, vv, vp, mz, levels)

   !---------------------------------------------------------------------------
   ! Purpose: Transform from fields on vertical EOFS to fields on vertical 
   ! levels.
   !
   ! Method:  Perform vp(i,j,k) = P E L^{1/2} vv(i,j,m) transform.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(in)  :: grid
   integer, intent(in)  :: mz                         ! # vertical modes.
   integer, intent(in)  :: levels                     ! # no. of levels  

   real*8,  intent(in)  :: evec(jds:jde,kds:kde,1:mz) ! Eigenvectors.
   real*8,  intent(in)  :: eval(jds:jde,1:mz)         ! Eigenvalues.
   real,    intent(in)  :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   real,    intent(in)  :: vv(ims:ime,jms:jme,kms:kme)   ! CV in EOF space.
   real,    intent(out) :: vp(ims:ime,jms:jme,kms:kme)! CV in level space.
   
   integer :: i, j, k, m, ij             ! Loop counters.
   real    :: temp

   if (trace_use_dull) call da_trace_entry("da_transform_vvtovp")

   !-------------------------------------------------------------------
   ! [1.0] Perform vp(i,j,k) = E L^{1/2} vv(i,j,m) transform:
   !------------------------------------------------------------------- 

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k, m, j, i, temp )
   do ij = 1 , grid%num_tiles
      vp(:,grid%j_start(ij):grid%j_end(ij),:) = 0.0
      do k = kts, levels
         do m = 1, mz
            do j = grid%j_start(ij), grid%j_end(ij)
               temp = evec(j,k,m) * eval(j,m)
   
               do i = its, ite
                  vp(i,j,k) = vp(i,j,k) + temp*vv(i,j,m)
               end do
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO
   
   !-------------------------------------------------------------------
   ! [2.0] Apply inner-product weighting if vertical_ip /= vertical_ip_0:
   !------------------------------------------------------------------- 

   if (vertical_ip /= vertical_ip_0) then
      vp(its:ite,jts:jte,kts:levels) = vp(its:ite,jts:jte,kts:levels) / &
         vertical_wgt(its:ite,jts:jte,kts:levels)                          
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_vvtovp")

end subroutine da_transform_vvtovp


subroutine da_transform_vvtovp_adj(grid, evec, eval, vertical_wgt, vp, vv, mz, levels)

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_transform_vvtovp.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(in) :: grid
   integer, intent(in) :: mz                         ! # vertical modes.
   integer, intent(in) :: levels                     ! no. of vertical levels

   real*8, intent(in)  :: evec(jds:jde,kds:kde,1:mz) ! Eigenvectors.
   real*8, intent(in)  :: eval(jds:jde,1:mz)         ! Eigenvalues.
   real, intent(in)    :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   real, intent(inout) :: vp(ims:ime,jms:jme,kms:kme)! CV in level space.
   real, intent(out)   :: vv(ims:ime,jms:jme,kms:kme)! CV in EOF space.
 
   integer :: i, j, m, k, ij             ! Loop counters.
   real    :: temp

   if (trace_use_dull) call da_trace_entry("da_transform_vvtovp_adj")

   !-------------------------------------------------------------------
   ! [1.0] Apply inner-product weighting if vertical_ip /= vertical_ip_0:
   !------------------------------------------------------------------- 

   if (vertical_ip /= vertical_ip_0) then
      vp(its:ite,jts:jte,kts:levels) = vp(its:ite,jts:jte,kts:levels) / &
         vertical_wgt(its:ite,jts:jte,kts:levels)
   end if

   !-------------------------------------------------------------------
   ! [2.0] Perform vp(i,j,k) = E L^{1/2} vv(i,j,m) transform:
   !------------------------------------------------------------------- 

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, k, j, i, temp )
   do ij = 1 , grid%num_tiles
      vv(:,grid%j_start(ij):grid%j_end(ij),:) = 0.0
      do m = 1, mz
         do k = kts, levels
            do j = grid%j_start(ij), grid%j_end(ij)
               temp = evec(j,k,m) * eval(j,m)
   
               do i = its, ite
                  vv(i,j,m) = vv(i,j,m) + temp*vp(i,j,k)
               end do
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   if (trace_use_dull) call da_trace_exit("da_transform_vvtovp_adj")

end subroutine da_transform_vvtovp_adj


subroutine da_transform_vptovv (evec, eval, vertical_wgt, vp, vv, mz, &
   kds,kde, ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte)

   !---------------------------------------------------------------------------
   ! Purpose: Transform from fields on vertical levels to fields on vertical 
   ! EOFS.
   !
   ! Method: Perform vv(i,j,m) = L^{-1/2} E^{T} vp(i,j,k) transform.
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)    :: mz                         ! # vertical modes.
   integer, intent(in)    :: kds,kde  ! domain dims.
   integer, intent(in)    :: ims,ime, jms,jme, kms,kme  ! memory dims.
   integer, intent(in)    :: its,ite, jts,jte, kts,kte  ! tile   dims
   real*8,  intent(in)    :: evec(jts:jte,kds:kde,1:mz) ! Eigenvectors.
   real*8,  intent(in)    :: eval(jts:jte,1:mz)         ! Eigenvalues.
   real,    intent(in)    :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   real,    intent(inout) :: vp(ims:ime,jms:jme,kms:kme)! CV in level space.
   real,    intent(out)   :: vv(ims:ime,jms:jme,1:mz)   ! CV in EOF space.
   
   integer :: i, j, m                    ! Loop counters.
   real    :: ETVp                       ! E(k,m)^{T}*vp(i,j,k)

   ! Do not add trace, as routine used by gen_be
   ! if (trace_use) call da_trace_entry("da_transform_vptovv")
   
   !-------------------------------------------------------------------
   ! [1.0] Apply inner-product weighting if vertical_ip /= vertical_ip_0
   !------------------------------------------------------------------- 

   if (vertical_ip /= vertical_ip_0) then
      vp(its:ite,jts:jte,kts:kte) = vp(its:ite,jts:jte,kts:kte) * &
                                    vertical_wgt(its:ite,jts:jte,kts:kte)
   end if

   !-------------------------------------------------------------------
   ! [2.0] Perform vv(i,j,m) = L^{-1/2} E^T vp(i,j,k) transform:
   !-------------------------------------------------------------------

   do m = 1, mz
      do j = jts, jte
         do i = its, ite
            ETVp = sum(evec(j,kts:kte,m) * vp(i,j,kts:kte))
            vv(i,j,m) = ETVp / eval(j,m)
         end do
      end do
   end do

end subroutine da_transform_vptovv


subroutine da_vertical_transform(grid, string, be, vertical_wgt, vv, vp)

   !---------------------------------------------------------------------
   ! Purpose: TBD
   !---------------------------------------------------------------------

   implicit none   

   type (domain),    intent(in)    :: grid
   character(len=*), intent(in)    :: string      ! Character operation
   type (be_type),   intent(in)    :: be          ! Background error structure.
   real,             intent(in)    :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   type (vp_type),   intent(inout) :: vv          ! CV in gridpt/EOF space.
   type (vp_type),   intent(inout) :: vp          ! CV in gridpt/level space.

   integer                         :: j, m, n     ! Loop counters.
   real                            :: alpha_stddev_inv ! 1/ sigma_alpha
   real                            :: size_inv    ! 1 / size.
   real                            :: alpha_me, alpha_ms, alpha_sd ! Alpha statistics.

   if (trace_use) call da_trace_entry("da_vertical_transform")

   select case(string)
      
   case ('u');
      
      !-------------------------------------------------------------------
      ! [1.0] Perform vp(i,j,k) = E L^{1/2} vv(i,j,m) transform:
      !------------------------------------------------------------------- 

      if (be % v1 % mz > 0) then
         call da_transform_vvtovp (grid, be % v1 % evec, be % v1 % val, vertical_wgt, &
            vv % v1, vp % v1, be % v1 % mz, kte)
      else
         vp % v1(its:ite,jts:jte,kts:kte) = 0.0
      end if

      if (be % v2 % mz > 0) then
         call da_transform_vvtovp (grid, be % v2 % evec, be % v2 % val, vertical_wgt, &
            vv % v2, vp % v2, be % v2 % mz, kte)
      else
         vp % v2(its:ite,jts:jte,kts:kte) = 0.0
      end if

      if (be % v3 % mz > 0) then
         call da_transform_vvtovp (grid, be % v3 % evec, be % v3 % val, vertical_wgt, &
            vv % v3, vp % v3, be % v3 % mz, kte)
      else
         vp % v3(its:ite,jts:jte,kts:kte) = 0.0
      end if

      if (be % v4 % mz > 0) then
         call da_transform_vvtovp (grid, be % v4 % evec, be % v4 % val, vertical_wgt, &
            vv % v4, vp % v4, be % v4 % mz, kte)
      else
         vp % v4(its:ite,jts:jte,kts:kte) = 0.0
      end if

      if (be % v5 % mz > 0) then
         if (global) then
            vp % v5(its:ite,jts:jte,1) = vv % v5(its:ite,jts:jte,1)
         else 
         call da_transform_vvtovp (grid, be % v5 % evec, be % v5 % val, vertical_wgt, & 
            vv % v5, vp % v5, be % v5 % mz, kts)    
         end if
      else
         vp % v5(its:ite,jts:jte,kts:kts) = 0.0
      end if
      if ( be % ne > 0 .and. be % alpha % mz > 0 ) then
         do n = 1, be % ne
            if ( anal_type_hybrid_dual_res ) then
               call da_transform_vvtovp_dual_res (grid%intermediate_grid, be % alpha % evec, be % alpha % val, vertical_wgt, &
                                       vv % alpha(:,:,:,n), vp % alpha(:,:,:,n), be % alpha % mz, kte_int)
            else
               call da_transform_vvtovp (grid, be % alpha % evec, be % alpha % val, vertical_wgt, &
                                       vv % alpha(:,:,:,n), vp % alpha(:,:,:,n), be % alpha % mz, kte_int)
            endif
         end do

!        Calculate alpha standard deviation diagnostic:
!         size_inv = 1.0 / ( (ite-its+1) * ( jte-jts+1) * be % ne * be % alpha % mz )
!         alpha_me = sum(vp % alpha(its:ite,jts:jte,:)) * size_inv
!         alpha_ms = sum(vp % alpha(its:ite,jts:jte,:) * vp % alpha(its:ite,jts:jte,:)) * &
!                    size_inv
!         alpha_sd = sqrt( alpha_ms - alpha_me * alpha_me )
!         write(6,'(a,f15.5)')' Alpha std. dev = ', alpha_sd
      end if

   case ('u_inv');
     
      !------------------------------------------------------------------- 
      ! [2.0] Perform vv(i,j,m) = L^{-1/2} E^T vp(i,j,k) transform:
      !------------------------------------------------------------------- 

      if (be % v1 % mz > 0) then
         call da_transform_vptovv (be % v1 % evec, be % v1 % val, vertical_wgt, &
            vp % v1, vv % v1, be % v1 % mz, kds,kde, ims,ime, jms,jme, kms,kme, &
            its,ite, jts,jte, kts,kte)
      end if

      if (be % v2 % mz > 0) then
         call da_transform_vptovv (be % v2 % evec, be % v2 % val, vertical_wgt, &
            vp % v2, vv % v2, be % v2 % mz, kds,kde, ims,ime, jms,jme, kms,kme, &
            its,ite, jts,jte, kts,kte)
      end if

      if (be % v3 % mz > 0) then
         call da_transform_vptovv (be % v3 % evec, be % v3 % val, vertical_wgt, &
            vp % v3, vv % v3, be % v3 % mz, kds,kde, ims,ime, jms,jme, kms,kme, &
            its,ite, jts,jte, kts,kte)
      end if

      if (be % v4 % mz > 0) then
         call da_transform_vptovv (be % v4 % evec, be % v4 % val, vertical_wgt, &
            vp % v4, vv % v4, be % v4 % mz, kds,kde, ims,ime, jms,jme, kms,kme, &
            its,ite, jts,jte, kts,kte)
      end if

      if (be % v5 % mz > 0) then
         if (global) then
            vv % v5(its:ite,jts:jte,1) = vp % v5(its:ite,jts:jte,1)
         else
            call da_transform_vptovv (be % v5 % evec, be % v5 % val, vertical_wgt, &
               vp % v5, vv % v5, be % v5 % mz, kds,kde, ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte)
         end if
      end if

      if ( be % ne > 0 .and. be % alpha % mz > 0 ) then
         do n = 1, be % ne
!           call da_transform_vptovv (be % alpha % evec, be % alpha % val, vertical_wgt, &
!                                     vp % alpha(:,:,:,n), vv % alpha(:,:,:,n), be % alpha % mz, kds,kde, &
!                                     ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte)
            call da_transform_vptovv (be % alpha % evec, be % alpha % val, vertical_wgt, &
                                      vp % alpha(:,:,:,n), vv % alpha(:,:,:,n), be % alpha % mz, kds_int,kde_int, &
                                      ims_int,ime_int, jms_int,jme_int, kms_int,kme_int, its_int,ite_int, &
                                      jts_int,jte_int, kts_int,kte_int)
         end do
      end if

   case ('u_adj');
    
      !------------------------------------------------------------------- 
      ! [3.0] Perform vv_adj = U_{v}^{T} vp_adj transform:
      !------------------------------------------------------------------- 

      if (be % v1 % mz > 0) then
         call da_transform_vvtovp_adj (grid, be % v1 % evec, be % v1 % val, vertical_wgt, &
            vp % v1, vv % v1, be % v1 % mz, kte)
      end if

      if (be % v2 % mz > 0) then
         call da_transform_vvtovp_adj (grid, be % v2 % evec, be % v2 % val, vertical_wgt, &
            vp % v2, vv % v2, be % v2 % mz, kte)
      end if

      if (be % v3 % mz > 0) then
         call da_transform_vvtovp_adj (grid, be % v3 % evec, be % v3 % val, vertical_wgt, &
            vp % v3, vv % v3, be % v3 % mz, kte)
      end if

      if (be % v4 % mz > 0) then
         call da_transform_vvtovp_adj (grid, be % v4 % evec, be % v4 % val, vertical_wgt, &
            vp % v4, vv % v4, be % v4 % mz, kte)
      end if

      if (be % v5 % mz > 0) then
         if (global) then
            vv % v5(its:ite,jts:jte,1) = vp % v5(its:ite,jts:jte,1)
         else
            call da_transform_vvtovp_adj (grid, be % v5 % evec, be % v5 % val, vertical_wgt, &
               vp % v5, vv % v5, be % v5 % mz, kts)
         end if
      end if

      if ( be % ne > 0 .and. be % alpha % mz > 0 ) then
         do n = 1, be % ne
            if ( anal_type_hybrid_dual_res ) then
               call da_transform_vvtovp_adj_dual_res (grid%intermediate_grid, be % alpha % evec, be % alpha % val, vertical_wgt, &
                                           vp % alpha(:,:,:,n), vv % alpha(:,:,:,n), be % alpha % mz, kte_int)
            else
               call da_transform_vvtovp_adj (grid, be % alpha % evec, be % alpha % val, vertical_wgt, &
                                           vp % alpha(:,:,:,n), vv % alpha(:,:,:,n), be % alpha % mz, kte_int)
            endif
         end do
      end if

   case default;
   
      call da_error("da_vertical_transform.inc",399, &
         (/"Invalid da_vertical_transform option "//string/))

   end select

   if (trace_use) call da_trace_exit("da_vertical_transform")

end subroutine da_vertical_transform


subroutine da_get_vpoles(u,v,        &
          ids, ide, jds, jde, &
          ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte  )

   !---------------------------------------------------------------------------
   !  Purpose: Treatment for Polar winds  
   !---------------------------------------------------------------------------
   
   implicit none
   
   integer, intent(in)    :: ids, ide, jds, jde
   integer, intent(in)    :: ims, ime, jms, jme, kms, kme
   integer, intent(in)    :: its, ite, jts, jte, kts, kte
   real, intent(inout)    :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real, intent(inout)    :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.
 
   real                   :: tmpvar                                         
   real                   :: tmpu,tmp_u,tmpv,tmp_v
   integer                :: k

   if (trace_use) call da_trace_entry("da_get_vpoles")

   ! cos_xls etc in da_control, calculated in da_setup_firstguess

   tmpvar       = 1.0/real(ide-ids+1)

   do k = kts,kte
      tmp_u =0.0
      tmp_v =0.0
      tmpu = 0.0
      tmpv = 0.0

      if (jts == jds) then 
         tmp_u = tmpvar*sum(-u(its:ite,jts+1,k)*cos_xls(its:ite)& 
                            +v(its:ite,jts+1,k)*sin_xls(its:ite))
         tmp_v = tmpvar*sum(-u(its:ite,jts+1,k)*sin_xls(its:ite)& 
                            -v(its:ite,jts+1,k)*cos_xls(its:ite))
      end if

      tmpu = wrf_dm_sum_real( tmp_u)
      tmpv = wrf_dm_sum_real( tmp_v)

      if (jts == jds) then 
        u(its:ite,jts,k) = -tmpu*cos_xls(its:ite) -tmpv*sin_xls(its:ite)
        v(its:ite,jts,k) =  tmpu*sin_xls(its:ite) -tmpv*cos_xls(its:ite)
      end if
 
      tmp_u =0.0
      tmp_v =0.0
      tmpu = 0.0
      tmpv = 0.0

      if (jte == jde) then 
         tmp_u = tmpvar*sum(-u(its:ite,jte-1,k)*cos_xle(its:ite)& 
                            -v(its:ite,jte-1,k)*sin_xle(its:ite))
         tmp_v = tmpvar*sum( u(its:ite,jte-1,k)*sin_xle(its:ite)& 
                            -v(its:ite,jte-1,k)*cos_xle(its:ite))
      end if
      tmpu = wrf_dm_sum_real( tmp_u)
      tmpv = wrf_dm_sum_real( tmp_v)
      if (jte == jde) then 
         u(its:ite,jte,k) = -tmpu*cos_xle(its:ite) +tmpv*sin_xle(its:ite)
         v(its:ite,jte,k) = -tmpu*sin_xle(its:ite) -tmpv*cos_xle(its:ite)
      end if
   end do

   if (trace_use) call da_trace_exit("da_get_vpoles")

end subroutine da_get_vpoles


subroutine da_get_spoles(x,              &
          ids, ide, jds, jde, &
          ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte  )

   !---------------------------------------------------------------------------
   ! Purpose: Treatment for Scalar field at Poles
   !---------------------------------------------------------------------------

   implicit none
   

   integer, intent(in)    :: ids, ide, jds, jde
   integer, intent(in)    :: ims, ime, jms, jme, kms, kme
   integer, intent(in)    :: its, ite, jts, jte, kts, kte
   real,    intent(inout) :: x(ims:ime,jms:jme,kms:kme)   

   integer                :: k
   real                   :: tmpvar,tmps,tmp_s

   if (trace_use) call da_trace_entry("da_get_spoles")

   tmpvar      = 1.0/real(ide-ids+1)

   do k = kts, kte
      tmps =0.0  ; tmp_s =0.0
      if (jts == jds)  tmp_s = tmpvar*sum(x(its:ite,jts+1,k))

      tmps = wrf_dm_sum_real( tmp_s)
      if (jts == jds) x(its:ite,jts,k) = tmps
 
      tmps =0.0  ; tmp_s =0.0
      if (jte == jde) tmp_s = tmpvar*sum(x(its:ite,jte-1,k))

      tmps = wrf_dm_sum_real( tmp_s)
      if (jte == jde) x(its:ite,jte,k) = tmps
   end do

   if (trace_use) call da_trace_exit("da_get_spoles")

end subroutine da_get_spoles


subroutine da_get_avpoles(u,v,       &
          ids, ide, jds, jde, &
          ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte  )

   !--------------------------------------------------------------------------- 
   ! Purpose: Treatment for Adjoint of Polar winds 
   !---------------------------------------------------------------------------

   implicit none
   
   integer, intent(in)    :: ids, ide, jds, jde
   integer, intent(in)    :: ims, ime, jms, jme, kms, kme
   integer, intent(in)    :: its, ite, jts, jte, kts, kte
   real,    intent(inout) :: u(ims:ime,jms:jme,kms:kme)   ! u wind comp.
   real,    intent(inout) :: v(ims:ime,jms:jme,kms:kme)   ! v wind comp.
 
   real                   :: tmpvar                                         
   real                   :: tmpu(kts:kte)
   real                   :: tmp_u(kts:kte)
   real                   :: tmpv(kts:kte)
   real                   :: tmp_v(kts:kte)
   integer                :: k

   if (trace_use) call da_trace_entry("da_get_avpoles")

   tmpvar      = 1.0/real(ide-ids+1)

   ! cos_xls etc in da_control, calculated in da_setup_firstguess

   tmp_u(:) =0.0
   tmp_v(:) =0.0

   if (jts == jds) then 
      do k = kts,kte
         tmp_u(k) = tmpvar*sum(-u(its:ite,jts,k)*cos_xls(its:ite) & 
                            +v(its:ite,jts,k)*sin_xls(its:ite))
         tmp_v(k) = tmpvar*sum(-u(its:ite,jts,k)*sin_xls(its:ite) & 
                            -v(its:ite,jts,k)*cos_xls(its:ite))
      end do
   end if

   call wrf_dm_sum_reals(tmp_u(:), tmpu(:))
   call wrf_dm_sum_reals(tmp_v(:), tmpv(:))

   if (jts == jds) then 
      do k = kts,kte
         u(its:ite,jts+1,k) = u(its:ite,jts+1,k) -tmpu(k)*cos_xls(its:ite) &
                              - tmpv(k)*sin_xls(its:ite)
         v(its:ite,jts+1,k) = v(its:ite,jts+1,k) +tmpu(k)*sin_xls(its:ite) &
                              - tmpv(k)*cos_xls(its:ite)
         u(its:ite,jts,k) = 0.0
         v(its:ite,jts,k) = 0.0
      end do
   end if

   tmp_u(:) =0.0
   tmp_v(:) =0.0

   if (jte == jde) then 
      do k = kts,kte
         tmp_u(k) = tmpvar*sum(-u(its:ite,jte,k)*cos_xle(its:ite) & 
                            -v(its:ite,jte,k)*sin_xle(its:ite))
         tmp_v(k) = tmpvar*sum( u(its:ite,jte,k)*sin_xle(its:ite) & 
                            -v(its:ite,jte,k)*cos_xle(its:ite))
      end do
   end if

   call wrf_dm_sum_reals(tmp_u(:), tmpu(:))
   call wrf_dm_sum_reals(tmp_v(:), tmpv(:))

   if (jte == jde) then 
      do k = kts,kte
         u(its:ite,jte-1,k) = u(its:ite,jte-1,k) -tmpu(k)*cos_xle(its:ite) &
                              + tmpv(k)*sin_xle(its:ite)
         v(its:ite,jte-1,k) = v(its:ite,jte-1,k) -tmpu(k)*sin_xle(its:ite) &
                              - tmpv(k)*cos_xle(its:ite)
         u(its:ite,jte,k) = 0.0
         v(its:ite,jte,k) = 0.0
      end do
   end if

   if (trace_use) call da_trace_exit("da_get_avpoles")

end subroutine da_get_avpoles


subroutine da_get_aspoles(x,              &
          ids, ide, jds, jde, &
          ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte  )

   !---------------------------------------------------------------------------
   ! Purpose: Treatment for Adjoint of Scalar field at Poles
   !---------------------------------------------------------------------------

   implicit none
   
   integer, intent(in)    :: ids, ide, jds, jde
   integer, intent(in)    :: ims, ime, jms, jme, kms, kme
   integer, intent(in)    :: its, ite, jts, jte, kts, kte
   real,    intent(inout) :: x(ims:ime,jms:jme,kms:kme)   

   integer                :: k
   real                   :: tmpvar
   real                   :: tmps(kts:kte)
   real                   :: tmp_s(kts:kte)

   tmpvar      = 1.0/real(ide-ids+1)

   if (trace_use) call da_trace_entry("da_get_aspoles")

   tmp_s(:) = 0.0

   if (jts == jds) then
      do k = kts, kte
         tmp_s(k) = tmpvar*sum(x(its:ite,jts,k))
      end do
   end if

   call wrf_dm_sum_reals(tmp_s(:),tmps(:))

   if (jts == jds) then
      do k = kts, kte
         x(its:ite,jts+1,k) = x(its:ite,jts+1,k) + tmps(k)
         x(its:ite,jts,k) = 0.0
      end do
   end if

   tmp_s(:) = 0.0

   if (jte == jde) then
      do k = kts, kte
         tmp_s(k) = tmpvar*sum(x(its:ite,jte,k))
      end do
   end if

   call wrf_dm_sum_reals(tmp_s(:),tmps(:))

   if (jte == jde) then
      do k = kts, kte
         x(its:ite,jte-1,k) = x(its:ite,jte-1,k) + tmps(k)
         x(its:ite,jte,k) = 0.0
      end do
   end if

   if (trace_use) call da_trace_exit("da_get_aspoles")

end subroutine da_get_aspoles


subroutine da_transform_vtovv_global (cv_size, xbx, be, cv, vv)

   !--------------------------------------------------------------------------
   ! Purpose: Control variable transform for global WRF-Var 
   !--------------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: cv_size       ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx           ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be            ! background errors.
   real,           intent(in)    :: cv(1:cv_size) ! control variables.
   type(vp_type),  intent(inout) :: vv            ! grdipt/eof cv (local).

   integer :: k, m, n ! Loop counters.
   integer :: cv_s ! Counter.
   integer :: cv_e ! Counter.

   if (trace_use) call da_trace_entry("da_transform_vtovv_global")

   !-------------------------------------------------------------------------
   ! [1] Spectral to grid transform for standard control variables:
   !-------------------------------------------------------------------------
  
   cv_s = 1        
   do k = 1, be%v1%mz
      cv_e = cv_s + 2 * be % cv % size1c - 1 
      call da_vtovv_spectral(be % v1 % max_wave, be % cv % size1c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v1%power(0:be % v1 % max_wave,k), &
         cv(cv_s:cv_e), vv%v1(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v2%mz
      cv_e = cv_s + 2 * be % cv % size2c - 1
      call da_vtovv_spectral(be % v2 % max_wave, be % cv % size2c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v2%power(0:be % v2 % max_wave,k), &
         cv(cv_s:cv_e), vv%v2(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v3%mz
      cv_e = cv_s + 2 * be % cv % size3c - 1
      call da_vtovv_spectral(be % v3 % max_wave, be % cv % size3c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v3%power(0:be % v3 % max_wave,k), &
         cv(cv_s:cv_e), vv%v3(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v4%mz
      cv_e = cv_s + 2 * be % cv % size4c - 1
      call da_vtovv_spectral(be % v4 % max_wave, be % cv % size4c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v4%power(0:be % v4 % max_wave,k), &
         cv(cv_s:cv_e), vv%v4(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v5%mz
     cv_e = cv_s + 2 * be % cv % size5c - 1
     call da_vtovv_spectral(be % v5 % max_wave, be % cv % size5c, &
        xbx % lenr, xbx % lenwrk, xbx % lensav, &
        xbx % inc, xbx % alp_size, xbx % alp, &
        xbx % wsave, be%v5%power(0:be % v5 % max_wave,k), &
        cv(cv_s:cv_e), vv%v5(its:ite,jts:jte,kts:kte))
     cv_s = cv_e + 1
   end do

   !-------------------------------------------------------------------------
   ! [2] Spectral to grid transform for flow-dependent control variables:
   !-------------------------------------------------------------------------
  
   do n = 1, be % ne
      do m = 1, be % alpha % mz
         cv_e = cv_s + 2 * be % cv % size_alphac - 1
         call da_vtovv_spectral(be % alpha % max_wave, be % cv % size_alphac, &
            xbx % lenr, xbx % lenwrk, xbx % lensav, &
            xbx % inc, xbx % alp_size, xbx % alp, &
            xbx % wsave, be % alpha % power(0:be%alpha%max_wave,1), &
            cv(cv_s:cv_e), vv%alpha(its:ite,jts:jte,m,n))
         cv_s = cv_e + 1
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_vtovv_global")

end subroutine da_transform_vtovv_global


subroutine da_transform_vtovv_global_adj (cv_size, xbx, be, cv, vv)

   !--------------------------------------------------------------------------
   ! Purpose: Control variable transform for global WRF-Var 
   !--------------------------------------------------------------------------

   implicit none

   integer,        intent(in)  :: cv_size       ! Size of cv array.
   type(xbx_type), intent(in)  :: xbx           ! For header & non-grid arrays.
   type(be_type),  intent(in)  :: be            ! background errors.
   real,           intent(out) :: cv(1:cv_size) ! control variables.
   type(vp_type),  intent(in)  :: vv            ! grdipt/eof cv (local).

   integer :: k, m, n ! Loop counters.
   integer :: cv_s ! Counter.
   integer :: cv_e ! Counter.

   if (trace_use) call da_trace_entry("da_transform_vtovv_global_adj")

   !-------------------------------------------------------------------------
   ! [1] Spectral to grid transform for standard control variables:
   !-------------------------------------------------------------------------

   cv_s = 1        
   do k = 1, be%v1%mz
      cv_e = cv_s + 2 * be % cv % size1c - 1
      call da_vtovv_spectral_adj(be % v1 % max_wave, be % cv % size1c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v1%power(0:be % v1 % max_wave,k), &
         cv(cv_s:cv_e), vv%v1(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v2%mz
      cv_e = cv_s + 2 * be % cv % size2c - 1
      call da_vtovv_spectral_adj(be % v2 % max_wave, be % cv % size2c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v2%power(0:be % v2 % max_wave,k), &
         cv(cv_s:cv_e), vv%v2(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v3%mz
      cv_e = cv_s + 2 * be % cv % size3c - 1
      call da_vtovv_spectral_adj(be % v3 % max_wave, be % cv % size3c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v3%power(0:be % v3 % max_wave,k), &
         cv(cv_s:cv_e), vv%v3(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v4%mz
      cv_e = cv_s + 2 * be % cv % size4c - 1
      call da_vtovv_spectral_adj(be % v4 % max_wave, be % cv % size4c, &
         xbx % lenr, xbx % lenwrk, xbx % lensav, &
         xbx % inc, xbx % alp_size, xbx % alp, &
         xbx % wsave, be%v4%power(0:be % v4 % max_wave,k), &
         cv(cv_s:cv_e), vv%v4(its:ite,jts:jte,k))
      cv_s = cv_e + 1
   end do

   do k = 1, be%v5%mz
     cv_e = cv_s + 2 * be % cv % size5c - 1 
     call da_vtovv_spectral_adj(be % v5 % max_wave, be % cv % size5c, &
        xbx % lenr, xbx % lenwrk, xbx % lensav, &
        xbx % inc, xbx % alp_size, xbx % alp, &
        xbx % wsave, be%v5%power(0:be % v5 % max_wave,k), &
        cv(cv_s:cv_e), vv%v5(its:ite,jts:jte,kts:kte))
     cv_s = cv_e + 1
   end do

   !-------------------------------------------------------------------------
   ! [2] Spectral to grid transform for flow-dependent control variables:
   !-------------------------------------------------------------------------

   do n = 1, be % ne
      do m = 1, be % alpha % mz
         cv_e = cv_s + 2 * be % cv % size_alphac - 1
         call da_vtovv_spectral_adj(be % alpha % max_wave, be % cv % size_alphac, &
            xbx % lenr, xbx % lenwrk, xbx % lensav, &
            xbx % inc, xbx % alp_size, xbx % alp, &
            xbx % wsave, be % alpha % power(0:be%alpha%max_wave,1), &
            cv(cv_s:cv_e), vv%alpha(its:ite,jts:jte,m,n))
         cv_s = cv_e + 1
      end do
   end do

   if (trace_use) call da_trace_exit("da_transform_vtovv_global_adj")

end subroutine da_transform_vtovv_global_adj



SUBROUTINE da_transform_bal( vp, be, grid )

   IMPLICIT NONE

   TYPE (vp_type), INTENT(INOUT)        :: vp ! work array.
   TYPE (be_type), INTENT(IN)           :: be ! Background errors.
   type (domain) , intent(inout)        :: grid   ! Domain variables.

   INTEGER                              :: i, j, k, kk, ij  ! Loop counters.
   
!-------------------------------------------------------------------
!  [1.0] Initialise:
!-------------------------------------------------------------------
!
!  linear balance btw psi and t-b, Psfc_b and chi_b 
!  [3.1] Calculate t_b from psi

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, k, ij)
   DO ij = 1, grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               grid%xa%t(i,j,k)=vp%v3(i,j,k)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, k, kk, ij)
   DO ij = 1, grid%num_tiles
      DO kk = kts,kte
         DO k = kts,kte
            DO j = grid%j_start(ij), grid%j_end(ij)
               DO i= its,ite
                  grid%xa%t(i,j,k) = grid%xa%t(i,j,k) + &
                                     be%agvz(i,j,k,kk) * vp%v1(i,j,kk)
               END DO
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!  [3.2] Calculate chi_b from psi

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, k, ij)
   DO ij = 1, grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               vp%v2(i,j,k) = vp%v2(i,j,k) + &
                              be%bvz(i,j,k) * vp%v1(i,j,k)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!  [3.3] Calculate Psfc_b from psi

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, ij)
   DO ij = 1, grid%num_tiles
      DO j = grid%j_start(ij), grid%j_end(ij)
         DO i= its,ite
            grid%xa%psfc(i,j)=vp%v5(i,j,1)
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO


   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, k, ij)
   DO ij = 1, grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               grid%xa%psfc(i,j) = grid%xa%psfc(i,j) + &
                                   be%wgvz(i,j,k) * vp%v1(i,j,k)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!--convert from delt.ln(ps) to delt.ps
   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, ij)
   DO ij = 1, grid%num_tiles
      DO j = grid%j_start(ij), grid%j_end(ij)
         DO i= its,ite
            grid%xa%psfc(i,j) = grid%xa%psfc(i,j) * grid%xb%psfc(i,j) 
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!  [3.4] Transform psi and chi to u and v:

!  Communicate halo region.
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_PSICHI_UV.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_PSICHI_UV_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   call da_psichi_to_uv( vp%v1, vp%v2, grid%xb%coefx, &
                         grid%xb%coefy, grid%xa%u, grid%xa%v    )

!  [3.5] treat humidity                         


   IF ( cv_options == 3 ) THEN
   !$OMP PARALLEL DO &
   !$OMP PRIVATE (i, j, k, ij)
   DO ij = 1, grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               grid%xa%q(i,j,k) = vp%v4(i,j,k) * grid%xb%qs(i,j,k)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO
   ELSE IF ( cv_options_hum == 1 ) THEN

      grid%xa%q(its:ite,jts:jte,kts:kte) = vp%v4(its:ite,jts:jte,kts:kte)

   ELSE IF ( cv_options_hum == 2 ) THEN

      grid%xa%rh(its:ite,jts:jte,kts:kte) = vp%v4(its:ite,jts:jte,kts:kte)

      CALL DA_TPRH_To_Q_Lin( grid )

   END IF

END SUBROUTINE da_transform_bal   

SUBROUTINE da_transform_bal_adj( vp, be, grid )

   IMPLICIT NONE

   TYPE (vp_type), INTENT(INOUT)        :: vp ! CV on grid structure.out
   TYPE (be_type), INTENT(IN)           :: be ! Background errors.
   type (domain) , intent(inout) :: grid   ! Dimensions and xpose buffers.

   INTEGER                              :: i, j, k, kk, ij, iunit  ! Loop counters.
!-------------------------------------------------------------------
!  [1.0] Initialise:
!------------------------------------------------------------------- 
  
!  linear balance btw psi and t-b, Psfc_b and chi_b 


!  [3.4]  adj of Transform psi and chi to u and v:        

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_PSICHI_UV_ADJ.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_PSICHI_UV_ADJ_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, ij )
   DO ij = 1 , grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               vp%v1(i,j,k)=0.
               vp%v2(i,j,k)=0.
            END DO
         END DO
      END DO
   ENDDO
   !$OMP END PARALLEL DO

   call da_psichi_to_uv_adj( grid%xa % u, grid%xa % v, grid%xb % coefx,   &
                            grid%xb % coefy,  vp % v1, vp % v2  )

!  [3.3]  adj Calculate Psfc_b from psi

!--convert from delt.ps to delt.ln(ps)
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, ij )
   DO ij = 1 , grid%num_tiles
      DO j = grid%j_start(ij), grid%j_end(ij)
         DO i= its,ite
            grid%xa%psfc(i,j) = grid%xa%psfc(i,j) * grid%xb%psfc(i,j)
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, ij )
   DO ij = 1 , grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               vp % v1(i,j,k)= vp % v1(i,j,k)+ &
                                be % wgvz(i,j,k) * grid%xa % psfc(i,j)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, ij )
   DO ij = 1 , grid%num_tiles
      DO j = grid%j_start(ij), grid%j_end(ij)
         DO i= its,ite
            vp % v5(i,j,1)=grid%xa % psfc(i,j)
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!  [3.2] adj Calculate chi_b from psi

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, ij )
   DO ij = 1 , grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               vp % v1(i,j,k)= vp % v1(i,j,k)+ &
                         be % bvz(i,j,k) * vp % v2(i,j,k) 
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!  [3.1] Calculate t_b from psi

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, kk, ij )
   DO ij = 1 , grid%num_tiles
      DO kk = kts,kte
         DO k = kts,kte
            DO j = grid%j_start(ij), grid%j_end(ij)
               DO i= its,ite
                  vp % v1(i,j,kk)= vp % v1(i,j,kk)+ &
                          be % agvz(i,j,k,kk) * grid%xa % t(i,j,k)     
               END DO
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, ij )
   DO ij = 1 , grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               vp % v3(i,j,k)= grid%xa % t(i,j,k)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

!  [3.5] treat humidity                         

   IF ( cv_options == 3 ) THEN
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( i, j, k, ij )
   DO ij = 1 , grid%num_tiles
      DO k = kts,kte
         DO j = grid%j_start(ij), grid%j_end(ij)
            DO i= its,ite
               vp % v4(i,j,k) = grid%xa % q(i,j,k) * grid%xb % qs(i,j,k)
            END DO
         END DO
      END DO
   END DO
   !$OMP END PARALLEL DO

   ELSE IF ( cv_options_hum == 1 ) THEN

      vp % v4(its:ite,jts:jte,kts:kte) = vp % v4(its:ite,jts:jte,kts:kte)+&
                                         grid%xa % q(its:ite,jts:jte,kts:kte)

   ELSE IF ( cv_options_hum == 2 ) THEN

      CALL DA_TPRH_To_Q_Adj( grid )

      vp % v4(its:ite,jts:jte,kts:kte) = vp % v4(its:ite,jts:jte,kts:kte)+&
                                         grid%xa % rh(its:ite,jts:jte,kts:kte)

   END IF

END SUBROUTINE da_transform_bal_adj   

SUBROUTINE da_apply_be( be, cv, vp, grid )

   IMPLICIT NONE

   TYPE (be_type), INTENT(IN)           :: be   ! Background error structure.
   REAL, INTENT(IN)                     :: cv(:)! Control variable.
   TYPE (vp_type), INTENT(INOUT)        :: vp   ! Grid point/EOF equivalent.
   type (domain) , intent(inout) :: grid   ! Dimensions and xpose buffers.

   INTEGER                              :: i,j,k,ij

!-------------------------------------------------------------------------
!  [1.0] Make local-grid copy of vp from 1-dimensional global-grid cv.
!-------------------------------------------------------------------------

   call da_cv_to_vv( cv_size, cv,&
           (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /), vp )

!  [2.0] Transform control variable:

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k, j, i )
   do ij = 1, grid%num_tiles
      do k=grid%xp%kts,grid%xp%kte
      do j=grid%j_start(ij), grid%j_end(ij)
      do i=grid%xp%its,grid%xp%ite
         vp % v1(i,j,k)=vp % v1(i,j,k)*be % corz(i,j,k,1)
         vp % v2(i,j,k)=vp % v2(i,j,k)*be % corz(i,j,k,2)
         vp % v3(i,j,k)=vp % v3(i,j,k)*be % corz(i,j,k,3)
         vp % v4(i,j,k)=vp % v4(i,j,k)*be % corz(i,j,k,4)
      enddo
      enddo
      enddo
   enddo
   !$OMP END PARALLEL DO

!-----Transform 5th control variable
      k=1
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j )
   do ij = 1, grid%num_tiles
      do j=grid%j_start(ij),grid%j_end(ij)
      do i=grid%xp%its,grid%xp%ite
         vp % v5(i,j,k)=vp % v5(i,j,k)*be % corp(i,j)
      enddo
      enddo
   enddo
   !$OMP END PARALLEL DO

   CALL da_apply_rf( be, vp , grid )

END SUBROUTINE da_apply_be          

SUBROUTINE da_apply_be_adj( be, cv, vp, grid )

   IMPLICIT NONE

   TYPE (be_type), INTENT(IN)           :: be     ! Background error structure.
   REAL, INTENT(INOUT)                  :: cv(:)  ! Control variable.
   TYPE (vp_type), INTENT(INOUT)        :: vp     ! Grid point/EOF equivalent.
   type (domain), intent(inout)         :: grid   ! Dimensions and xpose buffers.

   INTEGER                              :: i,j,k,ijk,ij,iunit

!-------------------------------------------------------------------------
!  [2.0] Transform from cv to vp:
!-------------------------------------------------------------------------

   CALL da_apply_rf_adj( be, vp, grid )

!  [2.1] Transform control variable:
!!!!!!!!!!!!!!!!!!!!!!!

   do k=kts,kte
   do j=jts,jte
   do i=its,ite
      vp % v1(i,j,k)=vp % v1(i,j,k)*be % corz(i,j,k,1)
      vp % v2(i,j,k)=vp % v2(i,j,k)*be % corz(i,j,k,2)
      vp % v3(i,j,k)=vp % v3(i,j,k)*be % corz(i,j,k,3)
      vp % v4(i,j,k)=vp % v4(i,j,k)*be % corz(i,j,k,4)
   enddo
   enddo         
   enddo         

!-----Transform 5th control variable
      k=1
   do j=jts,jte
      do i=its,ite
         vp % v5(i,j,k)=vp % v5(i,j,k)*be % corp(i,j)
      enddo
   enddo         

!-------------------------------------------------------------------------
!  [1.0] Make global-grid copy of cc from 3-dimensional local-grid vv.
!-------------------------------------------------------------------------

   call da_vv_to_cv( vp, grid%xp,&
          (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /), &
          cv_size, cv )

END SUBROUTINE da_apply_be_adj


subroutine da_transform_vvtovp_dual_res(grid, evec, eval, vertical_wgt, vv, vp, mz, levels)

   !---------------------------------------------------------------------------
   ! Purpose: Transform from fields on vertical EOFS to fields on vertical 
   ! levels.
   !
   ! Method:  Perform vp(i,j,k) = P E L^{1/2} vv(i,j,m) transform.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(in)  :: grid
   integer, intent(in)  :: mz                         ! # vertical modes.
   integer, intent(in)  :: levels                     ! # no. of levels  

   real*8,  intent(in)  :: evec(jds_int:jde_int,kds_int:kde_int,1:mz) ! Eigenvectors.
   real*8,  intent(in)  :: eval(jds_int:jde_int,1:mz)         ! Eigenvalues.
   real,    intent(in)  :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   real,    intent(in)  :: vv(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int)   ! CV in EOF space.
   real,    intent(out) :: vp(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int)! CV in level space.
   
   integer :: i, j, k, m, ij             ! Loop counters.
   real    :: temp

   if (trace_use_dull) call da_trace_entry("da_transform_vvtovp")

   !-------------------------------------------------------------------
   ! [1.0] Perform vp(i,j,k) = E L^{1/2} vv(i,j,m) transform:
   !------------------------------------------------------------------- 

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k, m, j, i, temp )
   do ij = 1 , grid%num_tiles
      vp(:,grid%j_start(ij):grid%j_end(ij),:) = 0.0
      do k = kts_int, levels
         do m = 1, mz
            do j = grid%j_start(ij), grid%j_end(ij)
               temp = evec(j,k,m) * eval(j,m)
   
               do i = its_int, ite_int
                  vp(i,j,k) = vp(i,j,k) + temp*vv(i,j,m)
               end do
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO
   
   !-------------------------------------------------------------------
   ! [2.0] Apply inner-product weighting if vertical_ip /= vertical_ip_0:
   !------------------------------------------------------------------- 

   if (vertical_ip /= vertical_ip_0) then
      vp(its:ite,jts:jte,kts:levels) = vp(its:ite,jts:jte,kts:levels) / &
         vertical_wgt(its:ite,jts:jte,kts:levels)                          
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_vvtovp")

end subroutine da_transform_vvtovp_dual_res


subroutine da_transform_vvtovp_adj_dual_res(grid, evec, eval, vertical_wgt, vp, vv, mz, levels)

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_transform_vvtovp.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain), intent(in) :: grid
   integer, intent(in) :: mz                         ! # vertical modes.
   integer, intent(in) :: levels                     ! no. of vertical levels

   real*8, intent(in)  :: evec(jds_int:jde_int,kds_int:kde_int,1:mz) ! Eigenvectors.
   real*8, intent(in)  :: eval(jds_int:jde_int,1:mz)         ! Eigenvalues.

!  real*8, intent(in)  :: evec(1:jde_int-jds_int+1,kds_int:kde_int,1:mz) ! Eigenvectors.
!  real*8, intent(in)  :: eval(1:jde_int-jds_int+1,1:mz)         ! Eigenvalues.
   real, intent(in)    :: vertical_wgt(ims:ime,jms:jme,kms:kme) ! Weighting.
   real, intent(inout) :: vp(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int)! CV in level space.
   real, intent(out)   :: vv(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int)! CV in EOF space.
 
   integer :: i, j, m, k, ij             ! Loop counters.
   real    :: temp

   if (trace_use_dull) call da_trace_entry("da_transform_vvtovp_adj")

   !-------------------------------------------------------------------
   ! [1.0] Apply inner-product weighting if vertical_ip /= vertical_ip_0:
   !------------------------------------------------------------------- 

   if (vertical_ip /= vertical_ip_0) then
      vp(its:ite,jts:jte,kts:levels) = vp(its:ite,jts:jte,kts:levels) / &
         vertical_wgt(its:ite,jts:jte,kts:levels)
   end if

   !-------------------------------------------------------------------
   ! [2.0] Perform vp(i,j,k) = E L^{1/2} vv(i,j,m) transform:
   !------------------------------------------------------------------- 

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, m, k, j, i, temp )
   do ij = 1 , grid%num_tiles
      vv(:,grid%j_start(ij):grid%j_end(ij),:) = 0.0
      do m = 1, mz
         do k = kts_int, levels
            do j = grid%j_start(ij), grid%j_end(ij)
               temp = evec(j,k,m) * eval(j,m)
   
               do i = its_int, ite_int
                  vv(i,j,m) = vv(i,j,m) + temp*vp(i,j,k)
               end do
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   if (trace_use_dull) call da_trace_exit("da_transform_vvtovp_adj")

end subroutine da_transform_vvtovp_adj_dual_res



end module da_vtox_transforms
