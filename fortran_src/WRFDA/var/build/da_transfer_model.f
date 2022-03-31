












module da_transfer_model

   
   
   

   use module_configure, only : grid_config_rec_type, model_config_rec
   use module_date_time, only : geth_julgmt, current_date, start_date
   use module_domain, only : domain, domain_clock_get, x_type, vp_type, ep_type
   use module_io_domain, only : open_r_dataset, close_dataset, input_auxinput17, &
      output_auxinput7, open_w_dataset
   use module_state_description, only : dyn_em_ad, dyn_em, dyn_em_tl, &
      p_qv, p_qh, p_qr, p_qi, p_qs, p_qg, p_qc, param_first_scalar, num_moist, &
      p_g_qv, p_g_qh, p_g_qr, p_g_qi, p_g_qs, p_g_qg, p_g_qc, &
      p_a_qv, p_a_qh, p_a_qr, p_a_qi, p_a_qs, p_a_qg, p_a_qc
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_reals
   use module_dm, only : local_communicator, &
      ntasks_x, ntasks_y, data_order_xyz, mytask, &
      ntasks, data_order_xy
   use module_comm_dm, only : halo_xa_sub, halo_init_sub, halo_psichi_uv_adj_sub, &
                              halo_xb_sub, halo_xb_uv_sub, halo_em_c_sub, halo_em_c_tl_sub, &
                              halo_xa_a_sub, halo_x6a_a_sub, halo_em_bdy_sub, halo_em_e_tl_sub, &
                              halo_em_e_sub

   use da_control, only : cos_xls, sin_xls, cos_xle, sin_xle, trace_use, &
      coarse_jy, coarse_ix, cone_factor, delt_lon, delt_lat, gas_constant, &
      map_projection,earth_omega,mix,pi,phic,mkz,start_lon,start_lat, &
      start_x,xlonc,start_y,mjy, global, rad_to_deg, deg_to_rad, earth_radius, &
      var4d,var4d_lbc,analysis_date,coarse_ds,analysis_accu,dsm,pole, fg_format_kma_global, &
      fg_format, fg_format_wrf_arw_regional, fg_format_wrf_nmm_regional, &  
      print_detail_map,stdout,truelat1_3dv, base_pres, fg_format_wrf_arw_global, &
      truelat2_3dv, periodic_x,write_increments,max_ext_its, gravity, &
      kappa, print_detail_xa,rd_over_rv,t0, print_detail_xa, check_rh, adj_sens,&
      print_detail_xb,test_dm_exact,base_lapse,base_temp,vertical_ip,ptop, &
      use_gpsztdobs, use_ssmitbobs, use_radarobs, use_radar_rf, use_radar_rhv,&
      dt_cloud_model, cp, use_ssmiretrievalobs, var4d_detail_out, &
      vertical_ip_sqrt_delta_p, vertical_ip_delta_p,check_rh_simple, check_rh_tpw, &
      t_kelvin, num_fgat_time, num_pseudo, iso_temp, interval_seconds, trajectory_io, &
      ids,ide,jds,jde,kds,kde, ims,ime,jms,jme,kms,kme, num_fft_factors, &
      its,ite,jts,jte,kts,kte, ips,ipe,jps,jpe,kps,kpe, qlimit, &
      update_sfcdiags, use_wrf_sfcinfo
   use da_control, only: base_pres_strat, base_lapse_strat
   use da_define_structures, only : xbx_type, be_type
   use da_par_util, only : da_patch_to_global
   use da_physics, only : da_check_rh_simple,da_roughness_from_lanu, &
      da_sfc_wtq,da_tpq_to_rh,da_trh_to_td,da_wrf_tpq_2_slp,da_integrat_dz, &
      da_tp_to_qs, da_check_rh,da_transform_xtogpsref, da_transform_xtoztd, &
      sfclayinit
   use da_reporting, only : da_error,message, da_message, da_warning
   use da_setup_structures, only : da_setup_runconstants,da_write_increments, &
      da_write_kma_increments,da_cloud_model, da_write_increments_for_wrf_nmm_regional
   use da_ssmi, only : da_transform_xtotb
   use da_tools, only : map_info, proj_merc, proj_ps,proj_lc,proj_latlon, &
      da_llxy_default,da_llxy_wrf,da_xyll,da_diff_seconds,da_map_set, &
      da_set_boundary_xb
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace
   use da_vtox_transforms, only : da_get_vpoles
   
   
   

   implicit none

   contains

subroutine da_transfer_wrftoxb(xbx, grid, config_flags)

   !---------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF to first guess structure.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------------

   implicit none
   
   type (xbx_type), intent(inout)     :: xbx        ! Header & non-gridded vars.

   type(domain), intent(inout)        :: grid
   type(grid_config_rec_type), intent(in) :: config_flags

   integer :: i, j, k, ij

   real    :: theta, tmpvar, alt

   real, dimension(ims:ime,jms:jme) :: rgh_fac

   character(len=19) :: current_date

   real :: loc_psac_mean

   real, dimension(jds:jde) :: loc_latc_mean

   integer :: size2d

   real, dimension(kms:kme) :: DDT

   real   :: qvf1, cvpm, cpovcv, p_surf, pfd, pfm, phm, pfu, ppb, ttb, albn, aln, height, temp
   real, allocatable :: arrayglobal(:,:)

   logical:: no_ppb
   logical :: has_znt, has_regime

   real, dimension(ims:ime,kms:kme) :: pf, pp



   call sfclayinit !for da_sfc_wtq

   ! If grid%pb does not existed in FG (YRG, 08/26/2010):
     ppb = sum(grid%pb*grid%pb)
     no_ppb = ppb == 0.0

   !---------------------------------------------------------------------------
   ! Set xb array range indices for processor subdomain.
   !---------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_transfer_wrftoxb")

   grid%xb % map  = grid%map_proj
   grid%xb % ds   = grid%dx

   grid%xb % mix = grid%xp % ide - grid%xp % ids + 1
   grid%xb % mjy = grid%xp % jde - grid%xp % jds + 1
   grid%xb % mkz = grid%xp % kde - grid%xp % kds + 1

   ! WHY?
   !rizvi   xbx%big_header%bhi(5,5) = grid%start_year
   !rizvi   xbx%big_header%bhi(6,5) = grid%start_month
   !rizvi   xbx%big_header%bhi(7,5) = grid%start_day
   !rizvi   xbx%big_header%bhi(8,5) = grid%start_hour

   !---------------------------------------------------------------------------
   ! WRF-specific fitelds:
   !---------------------------------------------------------------------------

   ptop = grid%p_top

   grid%xb%sigmaf(kte+1) = grid%znw(kte+1)

   grid%xb%znw(kte+1) = grid%znw(kte+1)
   grid%xb%znu(kte+1) = 0.0
 
   do k=kts,kte
      grid%xb%sigmah(k) = grid%znu(k)
      grid%xb%sigmaf(k) = grid%znw(k)

      grid%xb%znu(k) = grid%znu(k)
      grid%xb%znw(k) = grid%znw(k)
      grid%xb%dn(k)  = grid%dn(k)
      grid%xb%dnw(k) = grid%dnw(k)
   end do

   grid%xb % ptop = ptop
      
   !---------------------------------------------------------------------------
   ! Convert WRF fields to xb:
   !---------------------------------------------------------------------------

   if (print_detail_xb) then
      write(unit=stdout, fmt='(3a, i8)') &
         'file:', "da_transfer_wrftoxb.inc", ', line:', 99

      write(unit=stdout, fmt=*) 'its,ite=', its,ite
      write(unit=stdout, fmt=*) 'jts,jte=', jts,jte
      write(unit=stdout, fmt=*) 'kts,kte=', kts,kte

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         dnw                dn             rdnw       rdn'

      do k=kts,kte+1
         write(unit=stdout, fmt='(i3,8f16.8)') k, &
            grid%dnw(k), grid%dn(k), grid%rdnw(k), grid%rdn(k)
      end do

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         znu                 znw           rdnw       rdn'

      do k=kts,kte+1
         write(unit=stdout, fmt='(i3,8f16.8)') k, &
            grid%xb%sigmah(k), grid%xb%sigmaf(k), grid%rdnw(k), &
            grid%rdn(k)
      end do

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         phb                 ph_2'

      do k=kts,kte
         write(unit=stdout, fmt='(i3,8e20.12)') k, &
               grid%phb(its,jts,k), grid%ph_2(its,jts,k)
      end do

      write(unit=stdout, fmt=*) 'simple variables:'

      if (jte == jde) then
         write(unit=stdout, fmt=*) ' '

         do k=kts+5,kte,10
            do i=its,ite,10
               write(unit=stdout, fmt=*) &
                    '  grid%v_2(', i, ',', jde+1, ',', k, ')=', &
                       grid%v_2(i, jde+1,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      if (ite == ide) then
         write(unit=stdout, fmt=*) ' '

         do k=kts+5,kte,10
            do j=jts,jte,10
               write(unit=stdout, fmt=*) &
                  '  grid%u_2(', ide+1, ',', j, ',', k, ')=', &
                  grid%u_2(ide+1,j,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      write(unit=stdout, fmt=*) 'simple variables:'

      write(unit=stdout,fmt=*) &
         '  grid%u_1(its,jts,kts)=', grid%u_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%v_1(its,jts,kts)=', grid%v_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%w_1(its,jts,kts)=', grid%w_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%t_1(its,jts,kts)=', grid%t_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%ph_1(its,jts,kts)=',grid%ph_1(its,jts,kts)


      write(unit=stdout,fmt=*) &
         '  grid%u_2(its,jte,kts)=', grid%u_2(its,jte,kts)
      write(unit=stdout,fmt=*) &
         '  grid%v_2(ite,jts,kts)=', grid%v_2(ite,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%w_2(its,jts,kts)=', grid%w_2(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%t_2(its,jts,kts)=', grid%t_2(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%ph_2(its,jts,kts)=',grid%ph_2(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%phb(its,jts,kts)=', grid%phb(its,jts,kts)

      write(unit=stdout, fmt=*) &
         '  grid%sm31,grid%em31,grid%sm32,grid%em32, grid%sm33,grid%em33=', &
         grid%sm31,grid%em31,grid%sm32,grid%em32,grid%sm33,grid%em33

      write(unit=stdout, fmt=*) '  grid%p_top=', grid%p_top
      write(unit=stdout, fmt=*) '  grid%znu(kts)=', grid%znu(kts)
      write(unit=stdout, fmt=*) '  grid%mub(its,jts)=', grid%mub(its,jts)
      write(unit=stdout, fmt=*) '  grid%mu_2(its,jts)=', &
         grid%mu_2(its,jts)

      write(unit=stdout, fmt=*) '  hbot(its,jts)=', grid%hbot(its,jts)
      write(unit=stdout, fmt=*) '  htop(its,jts)=', grid%htop(its,jts)

      write(unit=stdout, fmt=*) '  grid%p_top=', grid%p_top
      write(unit=stdout, fmt=*) '  num_moist=', num_moist
      write(unit=stdout, fmt=*) '  P_QV=', P_QV

      write(unit=stdout, fmt=*) '  moist(its,jts,kts,2)=', &
         grid%moist(its,jts,kts,2)
      write(unit=stdout, fmt=*) ' '
   end if

   !---------------------------------------------------------------
   ! Need this to exchange values in the halo region.
   ! grid%xa%u and grid%xa%v are used as temporary arrays and so
   ! it is easy to use the existing exchange scheme.
   !
   ! Note, this is needed as u_2 and v_2 has no guarantee
   ! the most east column, and the most north row are
   ! properly initailized for each tile.
   !---------------------------------------------------------------

   ! Fill the halo region for u and v.

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_C.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_EM_C_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE


   if (print_detail_xb) then
      write(unit=stdout, fmt=*) &
         ' ids,ide,jds,jde,kds,kde=', ids,ide,jds,jde,kds,kde
      write(unit=stdout, fmt=*) &
         ' its,ite,jts,jte,kts,kte=', its,ite,jts,jte,kts,kte
      write(unit=stdout, fmt=*) &
          ' ims,ime,jms,jme,kms,kme=', ims,ime,jms,jme,kms,kme
         
      write(unit=stdout, fmt=*) &
         ' lbound(grid%xb%u)=',   lbound(grid%xb%u)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%xb%v)=',   lbound(grid%xb%v)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%u_2)=', lbound(grid%u_2)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%v_2)=', lbound(grid%v_2)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%xb%u)=',   ubound(grid%xb%u)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%xb%v)=',   ubound(grid%xb%v)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%u_2)=', ubound(grid%u_2)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%v_2)=', ubound(grid%v_2)
   end if

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j, k, cvpm, cpovcv, ppb, temp, ttb ) &
   !$OMP PRIVATE ( albn, qvf1, aln, theta ) &
   !$OMP PRIVATE ( p_surf, pfu, pfd, phm )
   do ij = 1 , grid%num_tiles

   do j=grid%j_start(ij), grid%j_end(ij)
      k = kte+1

      do i=its,ite
         grid%p(i,j,k) = 0.0
         grid%xb%map_factor(i,j) = grid%msft(i,j)
         grid%xb%cori(i,j) = grid%f(i,j)
         grid%xb%tgrn(i,j) = grid%sst(i,j)
         if (grid%xb%tgrn(i,j) < 100.0) &
            grid%xb%tgrn(i,j) = grid%tmn(i,j)
         grid%xb%lat(i,j) = grid%xlat(i,j)
         grid%xb%lon(i,j) = grid%xlong(i,j)
         grid%xb%terr(i,j) = grid%ht(i,j)
         grid%xb%snow(i,j) = grid%snowc(i,j)
         grid%xb%lanu(i,j) = grid%lu_index(i,j)
         grid%xb%landmask(i,j) = grid%landmask(i,j)
         grid%xb%xland(i,j) = grid%xland(i,j)
         ! Z. Liu below are variables used by RTTOV
         grid%xb%tsk(i,j) = grid%tsk(i,j)
         grid%xb%smois(i,j) = grid%smois(i,j,1)
         grid%xb%tslb(i,j) = grid%tslb(i,j,1)
         grid%xb%xice(i,j) = grid%xice(i,j)
         grid%xb%ivgtyp(i,j) = grid%ivgtyp(i,j)
         grid%xb%isltyp(i,j) = grid%isltyp(i,j)
         grid%xb%vegfra(i,j) = grid%vegfra(i,j)
         grid%xb%snowh(i,j) = grid%snowh(i,j)*1000.0 ! meter to mm    
         if ( grid%xb%ivgtyp(i,j) == grid%iswater .and. &
              grid%xb%snow(i,j) > 0.0 )then
            grid%xb%snow(i,j) = 0.0
            grid%xb%snowh(i,j) = 0.0
         end if
      end do
   end do

      ! WHY?
      ! Adapted the code from "real.init.code" by Y.-R. Guo 05/13/2004:

      ! do i=its,ite
      !    k = kte
      !    qvf1 = 0.5*(grid%moist(i,j,k,P_QV)+grid%moist(i,j,k,P_QV))
      !    qvf2 = 1.0/(1.0+qvf1)
      !    qvf1 = qvf1*qvf2
      !    grid%xb%p(i,j,k) = -0.5*(grid%mu_2(i,j)+qvf1* &
      !       grid%mub(i,j))/grid%rdnw(k)/qvf2

      !    do k = kte-1,1,-1
      !       qvf1 = 0.5*(grid%moist(i,j,k,P_QV)+grid%moist(i,j,k+1,P_QV))
      !       qvf2 = 1.0/(1.0+qvf1)
      !       qvf1 = qvf1*qvf2
      !       grid%p(i,j,k) = grid%p(i,j,k+1) - &
      !          (grid%mu_2(i,j)+qvf1*grid%mub(i,j))/qvf2/rdn(k+1)
      !    end do
      ! end do

      ! Adapted the code from WRF module_big_step_utilitites_em.F ----
      !         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

      ! NOTE: as of V3.1, P (pressure perturbation) and PB (base state pressure)
      ! are included in the wrfinput file. However, P and PB are still
      ! re-calculated here.

      ! NOTE: as of 7/01/2010, P and PB are directly from the first guess
      ! and no longer re-calculated here.

      cvpm =  - (1.0 - gas_constant/cp)
      cpovcv = cp / (cp - gas_constant)

      ! In case of var4d, pb etc. will be recalculated in start_em with realsize=8, 
      ! However, the originals are computed with realsize=4.
      if ( no_ppb ) then
         do k=kts,kte
          do j=grid%j_start(ij), grid%j_end(ij)
            do i=its,ite
               ! The base specific volume (from real.init.code)
               ppb  = grid%znu(k) * grid%mub(i,j) + ptop
               grid%pb(i,j,k) = ppb
               temp = MAX ( iso_temp, base_temp + base_lapse*log(ppb/base_pres) )
               if ( grid%pb(i,j,k) < base_pres_strat ) then
                  temp = iso_temp + base_lapse_strat*log(grid%pb(i,j,k)/base_pres_strat)
               end if
               ttb  = temp * (base_pres/ppb)**kappa
               ! ttb  = (base_temp + base_lapse*log(ppb/base_pres)) * &
               !   (base_pres/ppb)**kappa
               albn = (gas_constant/base_pres) * ttb * (ppb/base_pres)**cvpm

               qvf1 = 1.0 + grid%moist(i,j,k,P_QV) / rd_over_rv
               aln  = -1.0 / (grid%mub(i,j)+grid%mu_2(i,j)) * &
                      (albn*grid%mu_2(i,j) + grid%rdnw(k) * &
                      (grid%ph_2(i,j,k+1) - grid%ph_2(i,j,k)))
               ! total pressure:
               grid%xb%p(i,j,k) = base_pres * &
                                 ((gas_constant*(t0+grid%t_2(i,j,k))*qvf1) / &
                                 (base_pres*(aln+albn)))**cpovcv
               ! total density
               grid%xb%rho(i,j,k)= 1.0 / (albn+aln)
               ! pressure purtubation:
               grid%p(i,j,k) = grid%xb%p(i,j,k) - ppb
            end do
           end do
         end do
      else
         do j=grid%j_start(ij), grid%j_end(ij)
          do i=its,ite
            p_surf = base_pres * EXP ( -base_temp/base_lapse + ( (base_temp/base_lapse)**2 - 2.*gravity*grid%ht(i,j)/base_lapse/gas_constant ) **0.5 )
            do k = kts, kte
               grid%pb(i,j,k) = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
               temp = MAX ( iso_temp, base_temp + base_lapse*log(grid%pb(i,j,k)/base_pres) )
               if ( grid%pb(i,j,k) < base_pres_strat ) then
                  temp = iso_temp + base_lapse_strat*log(grid%pb(i,j,k)/base_pres_strat)
               end if
               grid%t_init(i,j,k) = temp*(base_pres/grid%pb(i,j,k))**(gas_constant/cp) - t0
               grid%alb(i,j,k) = (gas_constant/base_pres)*(grid%t_init(i,j,k)+t0)*(grid%pb(i,j,k)/base_pres)**cvpm
            end do
            grid%mub(i,j) = p_surf - grid%p_top
            grid%phb(i,j,1) = grid%ht(i,j) * gravity
            if (grid%hypsometric_opt == 1) then
               do k = kts,kte
                  grid%phb(i,j,k+1) = grid%phb(i,j,k) - grid%dnw(k)*grid%mub(i,j)*grid%alb(i,j,k)
               end do
            else if (grid%hypsometric_opt == 2) then
               do k = kts,kte
                  pfu = grid%mub(i,j)*grid%znw(k+1) + grid%p_top
                  pfd = grid%mub(i,j)*grid%znw(k  ) + grid%p_top
                  phm = grid%mub(i,j)*grid%znu(k  ) + grid%p_top
                  grid%phb(i,j,k+1) = grid%phb(i,j,k) + grid%alb(i,j,k)*phm*LOG(pfd/pfu)
               end do
            end if
          end do
         end do

         if (grid%hypsometric_opt == 1) then
           do k=kts,kte
             do j=grid%j_start(ij), grid%j_end(ij)
               do i=its,ite
                grid%al(i,j,k)=-1./(grid%mub(i,j)+grid%mu_2(i,j))*(grid%alb(i,j,k)*grid%mu_2(i,j)  &
                     +grid%rdnw(k)*(grid%ph_2(i,j,k+1)-grid%ph_2(i,j,k)))
                ! total density
                grid%xb%rho(i,j,k)= 1.0 / (grid%alb(i,j,k)+grid%al(i,j,k))
              end do
            end do
           end do
         elseif (grid%hypsometric_opt == 2) then
           do k=kts,kte
             do j=grid%j_start(ij), grid%j_end(ij)
              do i=its,ite
               pfu = (grid%mub(i,j)+grid%mu_2(i,j))*grid%znw(k+1)+grid%p_top
               pfd = (grid%mub(i,j)+grid%mu_2(i,j))*grid%znw(k)  +grid%p_top
               phm = (grid%mub(i,j)+grid%mu_2(i,j))*grid%znu(k)  +grid%p_top
               grid%al(i,j,k) = (grid%ph_2(i,j,k+1)-grid%ph_2(i,j,k)+grid%phb(i,j,k+1)-grid%phb(i,j,k)) &
                                 /phm/LOG(pfd/pfu)-grid%alb(i,j,k)
               ! total density
               grid%xb%rho(i,j,k)= 1.0 / (grid%alb(i,j,k)+grid%al(i,j,k))
            end do
           end do
          end do
         endif
         do k=kts,kte
           do j=grid%j_start(ij), grid%j_end(ij)
            do i=its,ite
              qvf1 = 1.+grid%moist(i,j,k,P_QV) / rd_over_rv
              grid%xb%p(i,j,k)=base_pres*( (gas_constant*(t0+grid%t_2(i,j,k))*qvf1)/       &
                         (base_pres*(grid%al(i,j,k)+grid%alb(i,j,k))) )**cpovcv  
            end do
           end do
         end do
      endif

      do k=kts,kte+1
        do j=grid%j_start(ij), grid%j_end(ij)
         do i=its,ite
            grid%xb%hf(i,j,k) = (grid%phb(i,j,k)+grid%ph_2(i,j,k))/gravity
            grid%xb%w (i,j,k) = grid%w_2(i,j,k)
         end do
        end do
      end do

      do j=grid%j_start(ij), grid%j_end(ij)

      if (grid%hypsometric_opt == 2) then
        ! Compute full-level pressure
        ! The full and half level dry hydrostatic pressure is easily computed (YRG, 12/20/2011):
        do i=its,ite
           do k=kts, kte+1
              pf(i,k) = (grid%mub(i,j)+grid%mu_2(i,j)) * grid%znw(k) + ptop
              pp(i,k) = (grid%mub(i,j)+grid%mu_2(i,j)) * grid%znu(k) + ptop
           enddo
        enddo
      endif

      do k=kts,kte
         do i=its,ite
            grid%xb%u(i,j,k) = 0.5*(grid%u_2(i,j,k)+grid%u_2(i+1,j,k))
            grid%xb%v(i,j,k) = 0.5*(grid%v_2(i,j,k)+grid%v_2(i,j+1,k))
            grid%xb%wh(i,j,k)= 0.5*(grid%xb%w(i,j,k)+grid%xb%w(i,j,k+1))
            if (grid%hypsometric_opt == 1) then
               grid%xb%h(i,j,k) = 0.5*(grid%xb%hf(i,j,k)+grid%xb%hf(i,j,k+1))
            elseif (grid%hypsometric_opt == 2) then
            !  This is almost correct for pressure, but not for height. 
            !  Arithmetic mean of full-level heights is always higher than the actual,
            !  leading to a biased model for height-based observations (e.g., GPS RO)
            !  and surface variables (2-meter TQ and 10-meter wind).
            !  wee 11/22/2011

               grid%xb%h(i,j,k) = grid%xb%hf(i,j,k)+(grid%xb%hf(i,j,k+1)-grid%xb%hf(i,j,k)) &
                    *LOG(pf(i,k)/pp(i,k))/LOG(pf(i,k)/pf(i,k+1))
            endif

            if ( num_pseudo == 0 ) then
               grid%moist(i,j,k,P_QV) = max(grid%moist(i,j,k,P_QV), qlimit)
            end if
            grid%xb%q(i,j,k) = grid%moist(i,j,k,P_QV)

            theta = t0 + grid%t_2(i,j,k)
            grid%xb%t(i,j,k) = theta*(grid%xb%p(i,j,k)/base_pres)**kappa

            ! Convert to specific humidity from mixing ratio of water vapor:
            grid%xb%q(i,j,k)=grid%xb%q(i,j,k)/(1.0+grid%xb%q(i,j,k))
   
            ! Background qrn needed for radar radial velocity assmiilation:

            if (size(grid%moist,dim=4) >= 4) then
               grid%xb%qcw(i,j,k) = max(grid%moist(i,j,k,p_qc), 0.0)
               grid%xb%qrn(i,j,k) = max(grid%moist(i,j,k,p_qr), 0.0)
!rizvi For doing single obs test with radiance one need to ensure non-zero hydrometeor
!  For doing so uncomment following two cards below
!               grid%xb%qcw(i,j,k) = 0.0001                             
!               grid%xb%qrn(i,j,k) = 0.0001                             

               grid%xb%qt (i,j,k) = grid%xb%q(i,j,k) + grid%xb%qcw(i,j,k) + &
                  grid%xb%qrn(i,j,k)
            end if

            if (size(grid%moist,dim=4) >= 6) then
               grid%xb%qci(i,j,k) = max(grid%moist(i,j,k,p_qi), 0.0)
               grid%xb%qsn(i,j,k) = max(grid%moist(i,j,k,p_qs), 0.0)
!rizvi For doing single obs test with radiance one need to ensure non-zero hydrometeor
!  For doing so uncomment following two cards below
!               grid%xb%qci(i,j,k) = 0.0001                             
!               grid%xb%qsn(i,j,k) = 0.0001                             
            end if

            if (size(grid%moist,dim=4) >= 7) then
               grid%xb%qgr(i,j,k) = max(grid%moist(i,j,k,p_qg), 0.0)
            end if

            if ( config_flags%mp_physics == 3 ) then   ! WSM3-class scheme
               if ( grid%xb%t(i,j,k) <= t_kelvin ) then
                  grid%xb%qci(i,j,k) = grid%xb%qcw(i,j,k)
                  grid%xb%qcw(i,j,k) = 0.0
                  grid%xb%qsn(i,j,k) = grid%xb%qrn(i,j,k)
                  grid%xb%qrn(i,j,k) = 0.0
               end if
            end if

         end do
      end do
   end do

   do j=grid%j_start(ij), grid%j_end(ij)
      do i=its,ite
         grid%xb%psac(i,j) = grid%mub(i,j)+grid%mu_2(i,j)
! To make the Psfc consistent with WRF (YRG, 04/06/2010):
         grid%xb%psfc(i,j) = grid%psfc(i,j)

         if (grid%xb%tgrn(i,j) < 100.0) &    
            grid%xb%tgrn(i,j) = grid%xb%t(i,j,kts)+ &
            0.0065*(grid%xb%h(i,j,kts)-grid%xb%hf(i,j,kts))
      end do
   end do

   end do
   !$OMP END PARALLEL DO

!   write (999,'("MABS=",e13.5)') sum(abs(grid%xb%psfc(:,:)-grid%psfc(:,:))) / &
!                                       float((ite-its+1)*(jte-jts+1))
      
   grid%xb%ztop = grid%xb%hf(its,jts,kte+1)

   if (print_detail_xb) then
      write(unit=stdout, fmt=*) ' '
      if (print_detail_xb) then
         write(unit=stdout, fmt='(/5a/)') &
            'lvl         h                 p                t'

         do k=kts,kte
            write(unit=stdout, fmt='(i3,8e20.12)') k, &
               grid%xb%h(its,jts,k), grid%xb%p(its,jts,k), grid%xb%t(its,jts,k)
         end do
      end if

      write(unit=stdout,fmt=*) ' '
      write(unit=stdout,fmt=*) 'grid%xb%u(its,jte,kte)=', grid%xb%u(its,jte,kte)
      write(unit=stdout,fmt=*) 'grid%xb%v(ite,jts,kte)=', grid%xb%v(ite,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%w(its,jts,kte)=', grid%xb%w(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%t(its,jts,kte)=', grid%xb%t(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%p(its,jts,kte)=', grid%xb%p(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%q(its,jts,kte)=', grid%xb%q(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%h(its,jts,kte)=', grid%xb%h(its,jts,kte)
      write(unit=stdout,fmt=*) &
         'grid%xb%hf(its,jts,kte)=', grid%xb%hf(its,jts,kte)
      write(unit=stdout,fmt=*) &
         'grid%xb%map_factor(its,jts)=', grid%xb%map_factor(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%cori(its,jts)=', grid%xb%cori(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%tgrn(its,jts)=', grid%xb%tgrn(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lat(its,jts)=', grid%xb%lat(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lon(its,jts)=', grid%xb%lon(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%terr(its,jts)=', grid%xb%terr(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%snow(its,jts)=', grid%xb%snow(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lanu(its,jts)=', grid%xb%lanu(its,jts)
      write(unit=stdout,fmt=*) &
         'grid%xb%landmask(its,jts)=', grid%xb%landmask(its,jts)
      write(unit=stdout,fmt=*) '(ite,jte)=', ite,jte                   
      write(unit=stdout,fmt=*) 'grid%xb%lat(ite,jte)=', grid%xb%lat(ite,jte)
      write(unit=stdout,fmt=*) 'grid%xb%lon(ite,jte)=', grid%xb%lon(ite,jte)
      write(unit=stdout,fmt=*) ' '
   end if

   !---------------------------------------------------------------------------
   ! [3.0] Calculate vertical inner product for use in vertical transform:
   !---------------------------------------------------------------------------
      
   if (vertical_ip == vertical_ip_sqrt_delta_p) then
      ! Vertical inner product is sqrt(Delta p):
      do k=kts,kte
         grid%xb % vertical_inner_product(its:ite,jts:jte,k) = &
            sqrt(grid%xb % psac(its:ite,jts:jte) * grid%xb%sigmah(k))
      end do 
   else if (vertical_ip == vertical_ip_delta_p) then

      ! Vertical inner product is Delta p:
      do k=1,grid%xb%mkz
         grid % xb % vertical_inner_product(its:ite,jts:jte,k) = &
         grid % xb % psac(its:ite,jts:jte) * grid%xb%sigmah(k)
      end do
   end if

   !---------------------------------------------------------------------------
   ! Roughness
   !---------------------------------------------------------------------------

   current_date = 'yyyy-mm-dd_hh:mm:ss'

   write(current_date(1:19), fmt='(i4.4, 5(a1, i2.2))') &
      grid%start_year, '-', &
      grid%start_month, '-', &
      grid%start_day, '_', &
      grid%start_hour, ':', &
      grid%start_minute, ':', &
      grid%start_second

   !xbx % mminlu = 'USGS'
   xbx % mminlu = trim(grid%mminlu)

   has_regime = sum(grid%regime*grid%regime) > 0.0
   has_znt    = sum(grid%znt*grid%znt)       > 0.0
   if ( has_regime .neqv. has_znt ) then
      ! make sure they are consistent
      ! either both from model, or both calculated here
      has_regime = .false.
      has_znt    = .false.
   end if
   if ( has_znt ) then
      grid%xb%rough = grid%znt
   else
      ! calculate rough only when it is not available from fg
      call da_roughness_from_lanu(19, xbx % mminlu, current_date, &
         grid%xb % lanu, grid%xb % rough)
   end if

   !---------------------------------------------------------------------------
   ! Calculate 1/grid box areas:
   !---------------------------------------------------------------------------

   if (print_detail_xb) then
      write(unit=stdout, fmt='(/a, e24.16)') &
         'grid%xb % ds=', grid%xb % ds

      write(unit=stdout, fmt='(a, e24.16/)') &
           'grid%xb % map_factor(its,jts)=', grid%xb % map_factor(its,jts)
   end if

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j, tmpvar, height, message )
   do ij = 1 , grid%num_tiles

   do j=grid%j_start(ij),grid%j_end(ij)
      do i=its,ite
         if (grid%xb%ztop < grid%xb%hf(i,j,kte+1)) &
             grid%xb%ztop = grid%xb%hf(i,j,kte+1)

         tmpvar = grid%xb%ds / grid%xb%map_factor(i,j)

         grid%xb % grid_box_area(i,j) = tmpvar*tmpvar

         ! Calculate surface variable(wind, moisture, temperature)
         ! sfc variables: 10-m wind, and 2-m T, Q, at cross points

         height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write (unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',grid%xb%h(i,j,kts) ,' terr =  ',grid%xb%terr(i,j)
            call da_error("da_transfer_wrftoxb.inc",701, message(1:2))
         end if

         if ( has_regime ) then
            ! if fg contains valid regime info, use it
            grid%xb%regime(i,j) = grid%regime(i,j)
         else
            ! xb t2/q2/u10/v10 will be coming directly from fg.
            ! but still call da_sfc_wtq here in order to get regimes.
            ! regimes can not be zero for sfc_assi_options=2.
            ! regimes calculated here could be very different from WRF values
            ! due to the lack of some input sfc info
            call da_sfc_wtq(grid%xb%psfc(i,j), grid%xb%tgrn(i,j), &
               grid%xb%p(i,j,kts), grid%xb%t(i,j,kts), grid%xb%q(i,j,kts), &
               grid%xb%u(i,j,kts), grid%xb%v(i,j,kts), &
               height,  grid%xb%rough(i,j),grid%xb%xland(i,j), grid%xb%ds, &
               grid%xb%u10(i,j), grid%xb%v10(i,j), grid%xb%t2(i,j), &
               grid%xb%q2(i,j), grid%xb%regime(i,j))
         end if !if has_regime

         ! use t2/q2/u10/v10 from fg
         grid%xb%u10(i,j) = grid%u10(i,j)
         grid%xb%v10(i,j) = grid%v10(i,j)
         grid%xb%t2(i,j) = grid%t2(i,j)
         grid%xb%q2(i,j) = grid%q2(i,j)

      end do
   end do

   end do
   !$OMP END PARALLEL DO


   !---------------------------------------------------------------------------
   ! Calculate saturation vapour pressure and relative humidity:
   !---------------------------------------------------------------------------

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k, j, i )
   do ij = 1 , grid%num_tiles
      do k=kts,kte
         do j=grid%j_start(ij),grid%j_end(ij)
            do i=its,ite
               call da_tpq_to_rh(grid%xb % t(i,j,k), grid%xb % p(i,j,k), &
                  grid%xb % q(i,j,k), grid%xb %es(i,j,k), grid%xb %qs(i,j,k), &
                  grid%xb %rh(i,j,k))
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO 

   !---------------------------------------------------------------------------
   ! Calculate dew point temperature:
   !---------------------------------------------------------------------------

   call da_trh_to_td (grid)

   if (print_detail_xb) then
      i=its; j=jts; k=kts

      write(unit=stdout, fmt=*) 'i,j,k=', i,j,k
      write(unit=stdout, fmt=*) 'grid%xb % td(i,j,k)=', grid%xb % td(i,j,k)
      write(unit=stdout, fmt=*) 'grid%xb % es(i,j,k)=', grid%xb % es(i,j,k)
      write(unit=stdout, fmt=*) 'grid%xb % rh(i,j,k)=', grid%xb % rh(i,j,k)
      write(unit=stdout, fmt=*) 'grid%xb % qs(i,j,k)=', grid%xb % qs(i,j,k)
      write(unit=stdout, fmt=*) ' '
   end if

   !---------------------------------------------------------------------------
   ! Sea level pressure and total precipitable water
   !---------------------------------------------------------------------------

   call da_wrf_tpq_2_slp (grid)

   ! WHY?
   ! do j = jts,jte
   !    do i = its,ite
   !       call da_tpq_to_slp(grid%xb%t(i,j,:), grid%xb%q(i,j,:), &
   !          grid%xb%p(i,j,:), grid%xb%terr(i,j), &
   !          grid%xb%psfc(i,j), grid%xb%slp(i,j), grid%xp)
   !    end do
   ! end do

   call da_integrat_dz(grid)

   !---------------------------------------------------------------------------
   ! Surface wind speed
   !---------------------------------------------------------------------------

   tmpvar = log(10.0/0.0001)

   !$OMP PARALLEL DO &
   !$OMP PRIVATE (ij, i, j, height)
   do ij = 1, grid%num_tiles

   do j=grid%j_start(ij), grid%j_end(ij)
      do i=its,ite
         height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)
         rgh_fac(i,j) = 1.0/log(height/0.0001)
         grid%xb%speed(i,j) = sqrt(grid%xb%u(i,j,kts)*grid%xb%u(i,j,kts) &
                         + grid%xb%v(i,j,kts)*grid%xb%v(i,j,kts) + 1.0e-6) &
                    *tmpvar*rgh_fac(i,j)
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   !---------------------------------------------------------------------------
   ! Brightness temperature SH Chen
   !---------------------------------------------------------------------------

   if (use_ssmitbobs)   &
      call da_transform_xtotb(grid)

   !---------------------------------------------------------------------------
   ! GPS Refractivity linked by Y.-R. Guo 05/28/2004
   !---------------------------------------------------------------------------

   call da_transform_xtogpsref(grid)

   !---------------------------------------------------------------------------
   ! Ground-based GPS ZTD must follow the GPS Refractivity calculation.
   !---------------------------------------------------------------------------

   ! WHY? For certain computation method, not current one.
   if (use_gpsztdobs) then
      call da_transform_xtoztd(grid)
      if (print_detail_xb) then
        i=its; j=jts
        write(unit=stdout, fmt=*) 'grid%xb % tpw(i,j)=', grid%xb % tpw(i,j)
        write(unit=stdout, fmt=*) 'grid%xb % ztd(i,j)=', grid%xb % ztd(i,j)
        write(unit=stdout, fmt=*) ' '
      end if
   end if

   !---------------------------------------------------------------------------
   ! Calculate means for later use in setting up background errors.
   !---------------------------------------------------------------------------

   ! WHY?
   ! if (.not. associated(xbx % latc_mean)) then
   allocate (xbx % latc_mean(jds:jde))
   if (trace_use) call da_trace("da_transfer_wrftoxb",&
      message="allocated xbx%latc_mean")
   ! end if

   size2d = (ide-ids+1)*(jde-jds+1)

   tmpvar = 1.0/real(size2d)

   ! Bitwitse-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not
   ! testing.  Thits will always be .false. for a serial or 1-process MPI run.
   if (test_dm_exact) then
      allocate(arrayglobal(ids:ide, jds:jde))
      call da_patch_to_global(grid,grid%xb%psac, arrayglobal)
      loc_psac_mean = tmpvar*sum(arrayglobal(ids:ide,jds:jde))
      deallocate(arrayglobal)
   else
      loc_latc_mean = 0.0
      loc_psac_mean = tmpvar*sum(grid%xb % psac(its:ite,jts:jte))
   end if

   tmpvar = 1.0/real(ide-ids+1)

   if (test_dm_exact) then
      allocate(arrayglobal(ids:ide, jds:jde))
      call da_patch_to_global(grid,grid%xb%lat, arrayglobal)
      do j=jds,jde
         loc_latc_mean(j) = tmpvar*sum(arrayglobal(ids:ide, j))
      end do
      deallocate(arrayglobal)
   else
      loc_latc_mean = 0.0
      do j=jts,jte
         loc_latc_mean(j) = tmpvar*sum(grid%xb % lat(its:ite, j))
      end do
   end if

   if (test_dm_exact) then
      ! Broadcast result from monitor to other tasks.
      call wrf_dm_bcast_real(loc_psac_mean, 1)
      xbx % psac_mean = loc_psac_mean
      ! Broadcast result from monitor to other tasks.
      call wrf_dm_bcast_real(loc_latc_mean, (jde-jds+1))
      xbx % latc_mean = loc_latc_mean
   else
      xbx % psac_mean = wrf_dm_sum_real(loc_psac_mean)
      call wrf_dm_sum_reals(loc_latc_mean, xbx % latc_mean)
   end if

   if (print_detail_xb) then
      ! write(unit=stdout, fmt=*) 'loc_psac_mean  =', loc_psac_mean
      write(unit=stdout, fmt=*) 'xbx % psac_mean=', xbx % psac_mean

      ! write(unit=stdout, fmt=*) 'loc_latc_mean  =', loc_latc_mean(jts)
      write(unit=stdout, fmt=*) 'xbx % latc_mean=', xbx % latc_mean(jts)
   end if


   ! Fill the halo region for xb        

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XB.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XB_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   ! Calculate time step from one dimensional cloud model parameterization

   if (dt_cloud_model) then
      do j = jts, jte
         do i = its, ite
            call da_cloud_model (grid%xb%t(I,J,:),  grid%xb%p(I,J,:), &
               grid%xb%q(I,J,:), grid%xb%qcw(I,J,:), grid%xb%qrn(I,J,:), &
               grid%xb%h(I,J,:), grid%xb%hf(I,J,:), ddt, kts, kte)

            do k = kts, kte
               grid%xb%delt(i,j,k) = DDT(k)
            end do
         end do
      end do
   end if

   deallocate (xbx % latc_mean)

   if (trace_use) call da_trace_exit("da_transfer_wrftoxb")

end subroutine da_transfer_wrftoxb

subroutine da_transfer_wrf_nmm_regional_toxb(xbx, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF-nmm model to first guess structure.
   ! Author :  Syed RH Rizvi,     MMM/NCAR    
   !           06/01/2008
   !---------------------------------------------------------------------------

   implicit none
   
   type (xbx_type), intent(inout) :: xbx        
   type(domain), intent(inout)    :: grid
   integer                        :: i, j, k
   character*19                   :: current_date
   real, allocatable              :: eta1(:), eta2(:)
   real                           :: tv, tmpvar

   if (trace_use) call da_trace_entry("da_transfer_wrf_nmm_regional_toxb")
               
   allocate (eta1(kts:kte+1))
   allocate (eta2(kts:kte+1))
   grid%xb % map  = grid%map_proj
   grid%xb % ds   = grid%dx  

   grid%xb % mix = grid%xp % ide - grid%xp % ids + 1
   grid%xb % mjy = grid%xp % jde - grid%xp % jds + 1
   grid%xb % mkz = grid%xp % kde - grid%xp % kds + 1


   ptop  = grid%p_top
   grid%xb % ptop = ptop
      
   eta1 = grid%znw
   eta2(kte+1) = 0.
   eta2(kts:kte) = grid%znu(kts:kte)
   if (print_detail_xb) then
      write(unit=stdout, fmt='(3a, i8)') &
         'file:', "da_transfer_wrf_nmm_regional_toxb.inc", ', line:', 38

      write(unit=stdout, fmt=*) 'its,ite=', its,ite
      write(unit=stdout, fmt=*) 'jts,jte=', jts,jte
      write(unit=stdout, fmt=*) 'kts,kte=', kts,kte

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         eta1                eta2    '
      do k=kts,kte+1
         write(unit=stdout, fmt='(i3,8f16.8)') k, &
            eta1(k), eta2(k)    
      end do

      write(unit=stdout,fmt=*) &
         '  grid%u(its,jte,kts)=', grid%u_2(its,jte,kts)
      write(unit=stdout,fmt=*) &
         '  grid%v(ite,jts,kts)=', grid%v_2(ite,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%t(its,jts,kts)=', grid%t_2(its,jts,kts)
      write(unit=stdout, fmt=*) &
         '  grid%sm31,grid%em31,grid%sm32,grid%em32, grid%sm33,grid%em33=', &
         grid%sm31,grid%em31,grid%sm32,grid%em32,grid%sm33,grid%em33

      write(unit=stdout, fmt=*) '  grid%p_top=', grid%p_top

      write(unit=stdout, fmt=*) '  num_moist=', num_moist
      write(unit=stdout, fmt=*) '  P_QV=', P_QV
      write(unit=stdout, fmt=*) '  P_QC=', P_QC
      write(unit=stdout, fmt=*) '  P_QR=', P_QR
      write(unit=stdout, fmt=*) '  P_QI=', P_QI
      write(unit=stdout, fmt=*) '  P_QS=', P_QS
      write(unit=stdout, fmt=*) '  P_QG=', P_QG

      write(unit=stdout, fmt=*) '  moist(its,jts,kts,p_qv)=', &
         grid%moist(its,jts,kts,p_qv)
      write(unit=stdout, fmt=*) ' '
   end if

   !---------------------------------------------------------------

   do j=jts,jte
      do k=kts,kte
         do i=its,ite+1
            grid%xa%u(i,j,k) = grid%u_2(i,j,k)
         end do
      end do
   end do

   do j=jts,jte+1
      do k=kts,kte
         do i=its,ite
            grid%xa%v(i,j,k) = grid%v_2(i,j,k)
         end do
      end do
   end do

   ! Fill the halo region for u and v.

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


   if (print_detail_xb) then
      write(unit=stdout, fmt=*) &
         ' ids,ide,jds,jde,kds,kde=', ids,ide,jds,jde,kds,kde
      write(unit=stdout, fmt=*) &
         ' its,ite,jts,jte,kts,kte=', its,ite,jts,jte,kts,kte
      write(unit=stdout, fmt=*) &
          ' ims,ime,jms,jme,kms,kme=', ims,ime,jms,jme,kms,kme
         
      write(unit=stdout, fmt=*) &
         ' lbound(grid%xb%u)=',   lbound(grid%xb%u)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%xb%v)=',   lbound(grid%xb%v)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%u)=', lbound(grid%u_2)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%v)=', lbound(grid%v_2)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%xb%u)=',   ubound(grid%xb%u)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%xb%v)=',   ubound(grid%xb%v)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%u)=', ubound(grid%u_2)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%v)=', ubound(grid%v_2)
   end if

   do j=jts,jte

      do i=its,ite
         grid%xb%map_factor(i,j) = grid%msftx(i,j)
         grid%xb%cori(i,j) = grid%f(i,j)
         grid%xb%tgrn(i,j) = grid%sst(i,j)
         if (grid%xb%tgrn(i,j) < 100.0) &
            grid%xb%tgrn(i,j) = grid%tmn(i,j)
         grid%xb%lat(i,j) = grid%xlat(i,j)
         grid%xb%lon(i,j) = grid%xlong(i,j)
         grid%xb%terr(i,j) = grid%ht(i,j) 
         grid%xb%hf(i,j,kts) = grid%ht(i,j) 
!                               PD          +   PDTOP           + PTOP
         grid%xb%psfc(i,j) = grid%psfc(i,j) + grid%mu_2(i,j) + grid%p_top
      end do

      do k=kts,kte+1
        do i=its,ite
          grid%xb%p(i,j,k) = eta1(k)* grid%mu_2(i,j)  + eta2(k)*grid%psfc(i,j) + grid%p_top
          grid%xa%w (i,j,k) = grid%w_2(i,j,k)
          grid%xb%w (i,j,k) = grid%w_2(i,j,k)
        end do
      end do
      do k=kts+1,kte+1
        do i=its,ite
         tv = grid%t_2(i,j,k)*(1.0 + 0.61 * grid%moist(i,j,k,p_qv) )
         grid%xb%hf(i,j,k) =   grid%xb%hf(i,j,k-1)+ &
            ( gas_constant* tv * log(grid%xb%p(i,j,k-1)/grid%xb%p(i,j,k)) )/gravity 
        end do
      end do


      do k=kts,kte
        do i=its,ite
          grid%xb%u(i,j,k) = 0.5*(grid%xa%u(i,j,k)+grid%xa%u(i+1,j,k))
          grid%xb%v(i,j,k) = 0.5*(grid%xa%v(i,j,k)+grid%xa%v(i,j+1,k))
          grid%xb%wh(i,j,k)= 0.5*(grid%xb%w(i,j,k)+grid%xb%w(i,j,k+1))
          grid%xb%h(i,j,k) = 0.5*(grid%xb%hf(i,j,k)+grid%xb%hf(i,j,k+1))
          grid%xb%t(i,j,k) = grid%t_2(i,j,k)
          grid%xb%q(i,j,k) = grid%moist(i,j,k,P_QV)
          if( num_pseudo == 0 .and. grid%xb%q(i,j,k) < 1.0e-6) &
             grid%xb%q(i,j,k) = 1.0e-6
        end do
      end do

      do i=its,ite
        if (grid%xb%tgrn(i,j) < 100.0) &    
         grid%xb%tgrn(i,j) = grid%xb%t(i,j,kts)+ &
         0.0065*(grid%xb%h(i,j,kts)-grid%xb%hf(i,j,kts))
      end do
   end do     ! Loop over Latitudes (j)

   grid%xb%ztop = grid%xb%hf(its,jts,kte+1)

   do j=jts,jte
      do i=its,ite
         if (grid%xb%ztop < grid%xb%hf(i,j,kte+1)) &
             grid%xb%ztop = grid%xb%hf(i,j,kte+1)

         tmpvar = grid%xb%ds / grid%xb%map_factor(i,j)

         grid%xb % grid_box_area(i,j) = tmpvar*tmpvar
      end do
   end do
   if (print_detail_xb) then
      write(unit=stdout, fmt=*) ' '
      if (print_detail_xb) then
         write(unit=stdout, fmt='(/5a/)') &
            'lvl         h                 p                t'
         do k=kts,kte
            write(unit=stdout, fmt='(i3,8e20.12)') k, &
               grid%xb%h(its,jts,k), grid%xb%p(its,jts,k), grid%xb%t(its,jts,k)
         end do
      end if

      write(unit=stdout,fmt=*) ' '
      write(unit=stdout,fmt=*) 'grid%xb%u(its,jte,kte)=', grid%xb%u(its,jte,kte)
      write(unit=stdout,fmt=*) 'grid%xb%v(ite,jts,kte)=', grid%xb%v(ite,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%t(its,jts,kte)=', grid%xb%t(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%p(its,jts,kte)=', grid%xb%p(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%q(its,jts,kte)=', grid%xb%q(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%h(its,jts,kte)=', grid%xb%h(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%w(its,jts,kte)=', grid%xb%w(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%cori(its,jts)=', grid%xb%cori(its,jts)
      write(unit=stdout,fmt=*) &
         'grid%xb%hf(its,jts,kte)=', grid%xb%hf(its,jts,kte)
      write(unit=stdout,fmt=*) &
         'grid%xb%map_factor(its,jts)=', grid%xb%map_factor(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%tgrn(its,jts)=', grid%xb%tgrn(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lat(its,jts)=', grid%xb%lat(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lat(ite,jte)=', grid%xb%lat(ite,jte)
      write(unit=stdout,fmt=*) 'grid%xb%lon(its,jts)=', grid%xb%lon(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lon(ite,jte)=', grid%xb%lon(ite,jte)
      write(unit=stdout,fmt=*) 'grid%xb%terr(its,jts)=', grid%xb%terr(its,jts)
      write(unit=stdout,fmt=*) ' '
   end if

   !---------------------------------------------------------------------------
   ! [3.0] Calculate vertical inner product for use in vertical transform:
   !---------------------------------------------------------------------------
      
   if (vertical_ip == vertical_ip_sqrt_delta_p) then
      ! Vertical inner product is sqrt(Delta p):
      do k=kts,kte
       grid%xb % vertical_inner_product(its:ite,jts:jte,k) = &
       sqrt(grid%xb%p(its:ite,jts:jte,k) -  grid%xb%p(its:ite,jts:jte,k+1) )
      end do 
   else if (vertical_ip == vertical_ip_delta_p) then

      ! Vertical inner product is Delta p:
      do k=1,grid%xb%mkz
       grid % xb % vertical_inner_product(its:ite,jts:jte,k) = &
       grid%xb%p(its:ite,jts:jte,k) -  grid%xb%p(its:ite,jts:jte,k+1) 
      end do
   end if


   current_date = 'yyyy-mm-dd_hh:mm:ss'

   write(current_date(1:19), fmt='(i4.4, 5(a1, i2.2))') &
      grid%start_year, '-', &
      grid%start_month, '-', &
      grid%start_day, '_', &
      grid%start_hour, ':', &
      grid%start_minute, ':', &
      grid%start_second


   if (print_detail_xb) then
      write(unit=stdout, fmt='(/a, e24.16)') &
         'grid%xb % ds=', grid%xb % ds

      write(unit=stdout, fmt='(a, e24.16/)') &
           'grid%xb % map_factor(its,jts)=', grid%xb % map_factor(its,jts)
   end if
   !---------------------------------------------------------------------------
   ! Calculate saturation vapour pressure and relative humidity:
   !---------------------------------------------------------------------------

   do j=jts,jte
      do k=kts,kte
         do i=its,ite
            call da_tpq_to_rh(grid%xb % t(i,j,k), grid%xb % p(i,j,k), &
               grid%xb % q(i,j,k), grid%xb %es(i,j,k), grid%xb %qs(i,j,k), &
               grid%xb %rh(i,j,k))
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_transfer_wrf_nmm_regional_toxb")

end subroutine da_transfer_wrf_nmm_regional_toxb

subroutine da_transfer_kmatoxb(xbx, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Transfers fields from KMA to first guess (xb) structure.
   !---------------------------------------------------------------------------

   implicit none
   
   type (xbx_type), intent(inout)     :: xbx          ! Header & non-gridded vars.

   type(domain), intent(inout)        :: grid

   integer :: i, j, k
   integer :: is, ie, js, je, ks, ke
   real    :: tmpvar, earth_radius_sq,conv
   real    :: tmpp,tmp_p,tmpps,tmp_ps
   ! real    :: rgh_fac(grid%xp%ims:grid%xp%ime,grid%xp%jms:grid%xp%jme)
   character(len=19) :: current_date
 
   real :: TV(grid%xp%kms:grid%xp%kme)
   real :: ALPHA(grid%xp%kms:grid%xp%kme)
   real :: PHALF(grid%xp%kms:grid%xp%kme)
   real :: PHALFL(grid%xp%kms:grid%xp%kme)                  

   real :: pu, pd, coef, delp, hydro, rgasg, shgt

   real, dimension(jds:jde) :: loc_latc_mean

   real, allocatable :: arrayglobal(:,:) 

   if (trace_use) call da_trace_entry("da_transfer_kmatoxb")

   if (use_ssmiretrievalobs) then
      call da_error("da_transfer_kmatoxb.inc",34, &
         (/"Cannot use ssmi obs with kma global runs"/))
   end if

   conv = pi/180.0
   earth_radius_sq = earth_radius*1000.0 * earth_radius*1000.0 * &
                     conv*(360.0/(coarse_ix-1))*(180.0/(coarse_jy-2))*conv
   COEF=0.6080                                               
   RGASG = gas_constant/gravity                                            

   !---------------------------------------------------------------------------
   ! Set xb array range indices for processor subdomain.
   !---------------------------------------------------------------------------

   is = grid%xp % its
   ie = grid%xp % ite
   js = grid%xp % jts
   je = grid%xp % jte
   ks = grid%xp % kts
   ke = grid%xp % kte

   grid%xb % ds    = grid%dx
   grid%xb % kma_a = grid%kma_a
   grid%xb % kma_b = grid%kma_b

   if (print_detail_xb) then
       write(unit=stdout, fmt='(/a/)') &
          'lvl         kma_a                 kma_b'

       do k=ks,ke
          write(unit=stdout, fmt='(i3,8e20.12)') k, grid%xb%kma_a(k), grid%xb%kma_b(k)
      end do
   end if

   grid%xb % mix = grid%xp % ide - grid%xp % ids + 1
   grid%xb % mjy = grid%xp % jde - grid%xp % jds + 1
   grid%xb % mkz = grid%xp % kde - grid%xp % kds + 1

   !---------------------------------------------------------------------------
   ! KMA-specific fields:
   !---------------------------------------------------------------------------

   ! Fix ptop as 0.4 hPa.
   ptop = 40.0     

   grid%xb % ptop = ptop
      
   !---------------------------------------------------------------------------
   ! Convert KMA fields to xb:
   !---------------------------------------------------------------------------

   grid%xb%lat(is:ie,js:je) =  grid%xlat(is:ie,js:je)
   grid%xb%lon(is:ie,js:je) = grid%xlong(is:ie,js:je)

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

   tmpvar  = 1.0/real(ide-ids+1)

   !---------------------------------------------------------------
   ! transfer u,v,t & q(specific humidity in g/g)
   !---------------------------------------------------------------

   do k=ks,ke
      do j=js,je
         do i=is,ie
            grid%xb%u(i,j,k) = grid%u_2(i,j,k)
            grid%xb%v(i,j,k) = grid%v_2(i,j,k)
            grid%xb%t(i,j,k) = grid%t_2(i,j,k)
            grid%xb%q(i,j,k) = grid%moist(i,j,k,P_QV)
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! Fix wind at Poles
   !---------------------------------------------------------------------------

    call da_get_vpoles(grid%xb%u,grid%xb%v,          &
       ids,ide,jds,jde, &
       grid%xp%ims,grid%xp%ime,grid%xp%jms,grid%xp%jme,grid%xp%kms,grid%xp%kme, &
       grid%xp%its,grid%xp%ite,grid%xp%jts,grid%xp%jte,grid%xp%kts,grid%xp%kte )

   if (print_detail_xb) then
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%its, grid%xp%ite=', grid%xp%its, grid%xp%ite
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%jts, grid%xp%jte=', grid%xp%jts, grid%xp%jte
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%kts, grid%xp%kte=', grid%xp%kts, grid%xp%kte
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%ims, grid%xp%ime=', grid%xp%ims, grid%xp%ime
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%jms, grid%xp%jme=', grid%xp%jms, grid%xp%jme
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%kms, grid%xp%kme=', grid%xp%kms, grid%xp%kme
      write(unit=stdout,fmt='(A,2I5)') 'ids, ide=', ids, ide
      write(unit=stdout,fmt='(A,2I5)') 'jds, jde=', jds, jde
      write(unit=stdout,fmt='(A,2I5)') 'grid%xp%kds, grid%xp%kde=', grid%xp%kds, grid%xp%kde

      write(unit=stdout,fmt='(A,I5)') 'size(grid%xb%u, dim=1)=', size(grid%xb%u, dim=1)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xb%u, dim=2)=', size(grid%xb%u, dim=2)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xb%u, dim=3)=', size(grid%xb%u, dim=3)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xb%v, dim=1)=', size(grid%xb%v, dim=1)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xb%v, dim=2)=', size(grid%xb%v, dim=2)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xb%v, dim=3)=', size(grid%xb%v, dim=3)

      write(unit=stdout,fmt='(A,I5)') 'size(grid%xa%u, dim=1)=', size(grid%xa%u, dim=1)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xa%u, dim=2)=', size(grid%xa%u, dim=2)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xa%u, dim=3)=', size(grid%xa%u, dim=3)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xa%v, dim=1)=', size(grid%xa%v, dim=1)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xa%v, dim=2)=', size(grid%xa%v, dim=2)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%xa%v, dim=3)=', size(grid%xa%v, dim=3)

      write(unit=stdout,fmt='(A,I5)') 'size(grid%u_2, dim=1)=', size(grid%u_2, dim=1)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%u_2, dim=2)=', size(grid%u_2, dim=2)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%u_2, dim=3)=', size(grid%u_2, dim=3)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%v_2, dim=1)=', size(grid%v_2, dim=1)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%v_2, dim=2)=', size(grid%v_2, dim=2)
      write(unit=stdout,fmt='(A,I5)') 'size(grid%v_2, dim=3)=', size(grid%v_2, dim=3)
   end if

   ALPHA(ke) = LOG(2.0)                                           

   do j=js,je
      do i=is,ie
         grid%xb%cori(i,j) = grid%f(i,j)
         grid%xb%psfc(i,j) = grid%psfc(i,j)
         grid%xb%terr(i,j) = grid%ht(i,j)

         ! Zhiquan Liu add some RTTOV variables
         !--------------------------------------
         grid%xb%t2(i,j) = grid%t2(i,j)
         grid%xb%q2(i,j) = grid%q2(i,j)
         grid%xb%u10(i,j) = grid%u10(i,j)
         grid%xb%v10(i,j) = grid%v10(i,j)
         grid%xb%tsk(i,j) = grid%tsk(i,j)
         ! currently KMA guess have no SST, replaced by TSK
         ! grid%xb%tgrn(i,j) = grid%sst(i,j)
         grid%xb%tgrn(i,j) = grid%tsk(i,j)
         grid%xb%landmask(i,j) = grid%landmask(i,j)
         grid%xb%xland(i,j) = grid%xland(i,j)

         grid%xb%smois(i,j) = grid%smois(i,j,1)
         grid%xb%tslb(i,j) = grid%tslb(i,j,1)
         grid%xb%xice(i,j) = grid%xice(i,j)
         grid%xb%ivgtyp(i,j) = grid%ivgtyp(i,j)
         grid%xb%isltyp(i,j) = grid%isltyp(i,j)
         grid%xb%vegfra(i,j) = grid%vegfra(i,j)
         ! convert snow water equivalent kg/m2 to snow depth(mm)
         grid%xb%snowh(i,j) = grid%snow(i,j)*10.0 

         !---------------------------------------------------------------------
         !  Syed RH Rizvi
         ! 
         ! The following code is to be activated after
         ! getting sst, land use, land mask and snow information for KMA grid
         ! This infor needs to be packed in "KMA2NETCDF" convertor
         !                              
         !---------------------------------------------------------------------

         ! grid%xb%tgrn(i,j) = grid%sst(i,j)
         ! if (grid%xb%tgrn(i,j) < 100.0) grid%xb%tgrn(i,j) = tmn(i,j)
         ! grid%xb%lanu(i,j)     = grid%lu_index(i,j)
         ! grid%xb%landmask(i,j) = grid%landmask(i,j)
         ! grid%xb%xland(i,j) = grid%xland(i,j)
         ! grid%xb%snow(i,j)     = grid%snowc(i,j)

         ! Since data is on full levels get full level pr & ht.
         ! Note p = A + B * Psfc formula needs Psfc to be in hPa 

         do K = ks,ke-1    
            PU  = grid%xb%KMA_A(K+1) + grid%xb%KMA_B(K+1)*grid%xb%psfc(i,j)/100.0
            PD  = grid%xb%KMA_A(K ) + grid%xb%KMA_B(K )*grid%xb%psfc(i,j)/100.0
            if (PU == PD)then
               message(1)='PU, PD equal and so denominator will be zero'
               write(unit=message(2),fmt='(A,3I5)') ' i, j, k = ',i,j,k
               write(unit=message(3),fmt='(A,2F10.2)') &
                  ' KMA_A ',grid%KMA_A(k),grid%KMA_A(K+1)
               write(unit=message(4),fmt='(A,2F10.2)') &
                  ' KMA_B ',grid%KMA_B(k),grid%KMA_B(K+1)
               write(unit=message(5),fmt='(A,2F10.2)') &
                  ' psfc(Pa) = ',grid%xb%psfc(i,j)
               call da_error("da_transfer_kmatoxb.inc",212,message(1:5))
            end if
            grid%xb%p(i,j,k)= 100.0 * exp ((PD*LOG(PD)-PU*LOG(PU))/(PD-PU) -1.0)
         end do 
 
         grid%xb%p(i,j,ke) = 100.0*(grid%xb%KMA_A(ke)+grid%xb%KMA_B(ke)*grid%xb%psfc(i,j)/100.0)/2.0

         do k=ks,ke
            if (grid%xb%t(i,j,k) <= 0.0) then
               write(unit=message(1),fmt='(A,3I5,F10.2)') &
                  'Problem in  temp = ',i,j,k,grid%xb%t(i,j,k)   
               call da_error("da_transfer_kmatoxb.inc",223,message(1:1))
            end if

            grid%xb%rho(i,j,k)=grid%xb%p(i,j,k)/(gas_constant*  &
               (1.0+COEF*grid%xb%q(I,J,K))*grid%xb%t(I,J,K))   
         end do

         ! compute full level height

         do K = ks,ke    
            PHALF(K) = grid%xb%KMA_A(K) + grid%xb%KMA_B(K)*grid%xb%psfc(I,J)/100.0                             
            TV   (K) = (1.0+COEF*grid%xb%q(I,J,K))*grid%xb%t(I,J,K)*rgasg    
         end do                                                           

         do K = ks,ke-1                                              
            DELP     = PHALF(K) - PHALF(K+1)                              
            PHALFL(K)= LOG(PHALF(K)/PHALF(K+1))                           
            ALPHA(K) = 1.0-PHALF(K+1)*PHALFL(K)/DELP                     
         end do  

         SHGT = grid%xb%terr(i,j)
         do K = ks, ke                                               
            grid%xb%h(I,J,K) = SHGT + ALPHA(K)*TV(K)                           
         end do 

         HYDRO = 0.0                                                    
         do K = ks+1, ke
            HYDRO = HYDRO + TV(K-1)*PHALFL(K-1)                           
            grid%xb%h(I,J,K) = grid%xb%h(I,J,K) + HYDRO                             
         end do                                                        
      end do
   end do

   !---------------------------------------------------------------------------
   ! Sigma values are needed for interpolating 
   ! background error statistics in vertical
   ! Fix sigmah values based on global average surface pressure 
   !                    and level wise pressure
   !---------------------------------------------------------------------------

   tmp_ps = sum(grid%xb%psfc(is:ie,js:je)) /real((ide-ids+1)*(jde-jds+1))

   tmpps = wrf_dm_sum_real(tmp_ps)

   do k=ks,ke
      tmp_p = sum(grid%xb%p(is:ie,js:je,k)) /real((ide-ids+1)*(jde-jds+1))
      tmpp = wrf_dm_sum_real(tmp_p)

      ! 0.01 is added to see that sigmah should not become negative
      grid%xb%sigmah(ke+ks-k) = (tmpp - ptop+0.01) / (tmpps- ptop+0.01) 
      if (grid%xb%sigmah(ke+ks-k) < 0.0) then
         write(unit=message(1),fmt='(A,F10.2)') &
            ' average surface pressure = ',tmpps  
         write(unit=message(2),fmt='(A,I5,A,F10.2)') &
            ' average pressure at this level= ',k,' = ',tmpp  
         write(unit=message(3),fmt='(A,I5,A,F10.2)') &
            ' negative sigmah(',ke+ks-k,') = ',grid%xb%sigmah(ke+ks-k) 
         call da_error("da_transfer_kmatoxb.inc",280,message(1:3))
      end if
   end do

   ! Fix ztop
   
   grid%xb%ztop = grid%xb%h(is,js,ke)

   if (print_detail_xb) then
      write(unit=stdout, fmt='(/5a/)') &
         'lvl         h                 p                t'

      do k=ks,ke
         write(unit=stdout, fmt='(i3,8e20.12)') k, &
            grid%xb%h(is,js,k), grid%xb%p(is,js,k), grid%xb%t(is,js,k)
      end do

      write (unit=stdout,fmt=*) ' '
      write (unit=stdout,fmt=*) 'grid%xb%u(is,je,ke)=', grid%xb%u(is,je,ke)
      write (unit=stdout,fmt=*) 'grid%xb%v(ie,js,ke)=', grid%xb%v(ie,js,ke)
      write (unit=stdout,fmt=*) 'grid%xb%w(is,js,ke)=', grid%xb%w(is,js,ke)
      write (unit=stdout,fmt=*) 'grid%xb%t(is,js,ke)=', grid%xb%t(is,js,ke)
      write (unit=stdout,fmt=*) 'grid%xb%p(is,js,ke)=', grid%xb%p(is,js,ke)
      write (unit=stdout,fmt=*) 'grid%xb%q(is,js,ke)=', grid%xb%q(is,js,ke)
      write (unit=stdout,fmt=*) 'grid%xb%hf(is,js,ke)=', grid%xb%hf(is,js,ke)
      write (unit=stdout,fmt=*) 'grid%xb%map_factor(is,js)=', grid%xb%map_factor(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%cori(is,js)=', grid%xb%cori(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%tgrn(is,js)=', grid%xb%tgrn(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%lat(is,js)=', grid%xb%lat(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%lon(is,js)=', grid%xb%lon(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%terr(is,js)=', grid%xb%terr(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%snow(is,js)=', grid%xb%snow(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%lanu(is,js)=', grid%xb%lanu(is,js)
      write (unit=stdout,fmt=*) 'grid%xb%landmask(is,js)=', grid%xb%landmask(is,js)
      write (unit=stdout,fmt=*) ' '
   end if

   write(current_date(1:19), fmt='(i4.4, 5(a1, i2.2))') &
      grid%start_year, '-', &
      grid%start_month, '-', &
      grid%start_day,   '_', &
      grid%start_hour,  ':', &
      grid%start_minute,':', &
      grid%start_second

   write(unit=stdout,fmt='(2A)') 'Current date is ',current_date 

   !---------------------------------------------------------------------------
   ! Syed RH Rizvi
   ! 
   ! Following code for calculating roughness length needs to be activated  
   ! after getting land use data for KMA grid 
   !---------------------------------------------------------------------------

   !
   ! Set proper value for "mminlu" if using KMA info

   !xbx % mminlu = 'USGS'
   !call da_roughness_from_lanu(19, xbx % mminlu, current_date, grid%xp, &
   !   grid%xb % lanu(grid%xp%ims,grid%xp%jms), grid%xb % rough(grid%xp%ims,grid%xp%jms))

   do j=js,je
      do i=is,ie
         if (grid%xb%ztop < grid%xb%hf(i,j,ke+1)) grid%xb%ztop = grid%xb%hf(i,j,ke+1)
         ! Calculate grid_box area and vertical inner product for use in 
         ! vertical transform:
         grid%xb % grid_box_area(i,j) = earth_radius_sq * cos(grid%xlat(i,j)*conv)

         if (vertical_ip == vertical_ip_sqrt_delta_p) then
            ! Vertical inner product is sqrt(Delta p):
            do k=1,grid%xb%mkz
               grid%xb % vertical_inner_product(i,j,k) = &
                  sqrt(grid%xb%p(i,j,k)-grid%xb%p(i,j,k+1))
            end do
         else if (vertical_ip == vertical_ip_delta_p) then
            ! Vertical inner product is Delta p:
            do k=1,grid%xb%mkz
            grid%xb % vertical_inner_product(i,j,k) = grid%xb%p(i,j,k)-grid%xb%p(i,j,k+1)
            end do
         end if

         !------------------------------------------------------------------------------
         !  Syed RH Rizvi
         ! 
         ! Following code for calculating roughness length needs to be activated  
         ! to calculate surface variables (10-m wind, and 2-m T, Q) , 
         ! After testing KMA-WRFVAR for SFC_ASSIM_OPTIONS = 2
         !
         !------------------------------------------------------------------------------

         ! call da_sfc_wtq(grid%xb%psfc(i,j), grid%xb%tgrn(i,j), &
         !    grid%xb%p(i,j,ks), grid%xb%t(i,j,ks), grid%xb%q(i,j,ks), &
         !    grid%xb%u(i,j,ks), grid%xb%v(i,j,ks), &
         !    grid%xb%p(i,j,ks+1), grid%xb%t(i,j,ks+1), grid%xb%q(i,j,ks+1), &
         !    grid%xb%h(i,j,ks), grid%xb%rough(i,j),grid%xb%landmask(i,j), &
         !    grid%xb%u10(i,j), grid%xb%v10(i,j), grid%xb%t2(i,j), grid%xb%q2(i,j), &
         !    grid%xb%regime(i,j))

         ! write (unit=stdout,fmt='(7I5)') &
         !    i,j,grid%xb%psfc(i,j),grid%xb%t2(i,j),grid%xb%q2(i,j), &
         !    grid%xb%u10(i,j),grid%xb%v10(i,j)
         ! ------------------------------------------------------------------------------

      end do
   end do

   !---------------------------------------------------------------------------
   ! Calculate saturation vapour pressure and relative humidity:
   !---------------------------------------------------------------------------

   do j=js,je
      do k=ks,ke
         do i=is,ie
            call da_tpq_to_rh(grid%xb % t(i,j,k), grid%xb % p(i,j,k), &
               grid%xb % q(i,j,k), &
               grid%xb %es(i,j,k), grid%xb %qs(i,j,k), grid%xb %rh(i,j,k))
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! Calculate dew point temperature:
   !---------------------------------------------------------------------------

   call da_trh_to_td (grid)

   if (print_detail_xb) then
      i=is; j=js; k=ks

      write (unit=stdout,fmt='(A,3I5)') 'i,j,k=', i,j,k
      write (unit=stdout,fmt='(A,F10.2)') 'grid%xb % td(i,j,k)=', grid%xb % td(i,j,k)
      write (unit=stdout,fmt='(A,F10.2)') 'grid%xb % es(i,j,k)=', grid%xb % es(i,j,k)
      write (unit=stdout,fmt='(A,F10.2)') 'grid%xb % rh(i,j,k)=', grid%xb % rh(i,j,k)
      write (unit=stdout,fmt='(A,F10.2)') 'grid%xb % qs(i,j,k)=', grid%xb % qs(i,j,k)
      write (unit=stdout,fmt='(A)') ' '
   end if

   !---------------------------------------------------------------------------
   ! Sea level pressure and total precipitable water
   !---------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! Following code for calculating roughness length needs to be activated  
   ! for calculating roughness length if sea level pressure is desired
   !---------------------------------------------------------------------------

   ! call da_wrf_tpq_2_slp (grid)

   call da_integrat_dz(grid)
   !---------------------------------------------------------------------------
   ! Following code for calculating roughness length needs to be activated  
   ! for surface wind speed for SFC_ASSIM_OPTIONS = 2 
   ! 
   ! It is not working because KMA terrain height output from wave to grid 
   ! transform is negative at many grid points
   !---------------------------------------------------------------------------

   ! tmpvar = log(10.0/0.0001)
   ! do j=js,je
   !    do i=is,ie
   !       if (grid%xb%hf(i,j,ks) <= 0) then
   !          write(unit=message(1), fmt=*) ' zero hf at i/j ',i,j,' ks= ',ks
   !          call da_error("da_transfer_kmatoxb.inc",442,message(1:1))
   !       end if
   !       rgh_fac(i,j) = 1.0/log(grid%xb%hf(i,j,ks)/0.0001)

   !       grid%xb%speed(i,j) = sqrt(grid%xb%u(i,j,ks)*grid%xb%u(i,j,ks) &
   !          + grid%xb%v(i,j,ks)*grid%xb%v(i,j,ks) + 1.0e-6) &
   !          *tmpvar*rgh_fac(i,j)
   !    end do
   ! end do

   !---------------------------------------------------------------------------
   ! Following code for calculating roughness length needs to be activated  
   ! if SSMI brightness temperature are used
   !---------------------------------------------------------------------------

   ! call da_transform_xtotb(xa)

   !---------------------------------------------------------------------------
   ! Calculate means for later use in setting up background errors.
   !---------------------------------------------------------------------------

   allocate (xbx % latc_mean(jds:jde))

   tmpvar = 1.0/real(ide-ids+1)
   loc_latc_mean(:) = 0.0

   ! Bitwise-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not     
   ! testing.  This will always be .false. for a serial or 1-process MPI run.  
   if (test_dm_exact) then
      allocate(arrayglobal(ids:ide, jds:jde))
      call da_patch_to_global(grid, grid%xb%lat, arrayglobal)
      do j=jds,jde
         loc_latc_mean(j) = tmpvar*sum(arrayglobal(ids:ide, j))
      end do
      deallocate(arrayglobal)
      ! Broadcast result from monitor to other tasks.
      call wrf_dm_bcast_real(loc_latc_mean, (jde-jds+1))
      xbx % latc_mean = loc_latc_mean
   else
      do j=js,je
         loc_latc_mean(j) = tmpvar*sum(grid%xb % lat(is:ie, j))
      end do
      call wrf_dm_sum_reals (loc_latc_mean, xbx % latc_mean)
   end if

   ! WHY?
   !  do j=js,je
   !    do i=is,ie
   !       write(unit=stdout,fmt='(2i4,3f12.2,e12.4,2f12.2)') &
   !          j,i,grid%xb%psfc(i,j),grid%xb%tsk(i,j),grid%xb%t2(i,j), &
   !          grid%xb%q2(i,j),grid%xb%u10(i,j),grid%xb%v10(i,j)
   !    end do
   ! end do

   if (trace_use) call da_trace_exit("da_transfer_kmatoxb")

end subroutine da_transfer_kmatoxb


subroutine da_transfer_xatowrf(grid, config_flags)

   !---------------------------------------------------------------------------
   !  Purpose: Convert analysis increments into WRF increments
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !
   !  The following WRF fields are modified:  
   !    grid%u_2
   !    grid%v_2
   !    grid%w_2
   !    grid%mu_2
   !    grid%ph_2
   !    grid%t_2
   !    grid%moist
   !    grid%p
   !    grid%psfc
   !    grid%t2, grid%q2, grid%u10, grid%v10, grid%th2
   !
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)            :: grid
   type(grid_config_rec_type), intent(in) :: config_flags

   integer :: i, j, k

   real    :: sdmd, s1md

   ! arrays to hold wrf increments on the c-grid 

   real, dimension(ims:ime,jms:jme, kms:kme) :: &
      u_cgrid, v_cgrid, q_cgrid, ph_cgrid

   real, dimension(ims:ime,jms:jme) :: mu_cgrid

   real :: t_full, p_full, rho_dry, q_full, ph_full, ph_xb_hd, &
           qvf1, qvf2, qvf1_b, qvf2_b

   real :: uu, vv, ps1, ts1, qv1, height
   real :: zs1, zs2, pf2, theta, thetam, mu_full, pfu, pfd, phm, cvpm, p2m
   real, dimension(kms:kme) :: ald, ph
   logical :: has_lsm_info

   if (trace_use) call da_trace_entry("da_transfer_xatowrf")

   has_lsm_info = .false.
   if ( config_flags%sf_surface_physics == 2 ) then
       if ( sum(grid%hfx*grid%hfx)   > 0.0 .and. &
            sum(grid%qfx*grid%qfx)   > 0.0 ) then
          has_lsm_info = .true.
       end if
   end if

   ! To keep the background PH perturbation:

   do j=jts,jte
      do i=its,ite
         do k=kts, kte+1
            ph_cgrid(i,j,k) = grid%ph_2(i,j,k)
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! [1.0] Get the mixing ratio of moisture first, as it its easy.
   !---------------------------------------------------------------------------

   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            if ((grid%xb%q(i,j,k)+grid%xa%q(i,j,k)) < 0.0) then
               q_cgrid(i,j,k) =-grid%xb%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))**2
            else
               q_cgrid(i,j,k) = grid%xa%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))**2
            end if
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! [2.0] compute increments of dry-column air mass per unit area
   !---------------------------------------------------------------------------

   do j=jts,jte
      do i=its,ite
         sdmd=0.0
         s1md=0.0
         do k=kts,kte
            sdmd=sdmd+q_cgrid(i,j,k)*grid%dnw(k)
            s1md=s1md+(1.0+grid%moist(i,j,k,P_QV))*grid%dnw(k)
         end do

         mu_cgrid(i,j)=-(grid%xa%psfc(i,j)+grid%xb%psac(i,j)*sdmd)/s1md
      end do
   end do

   !---------------------------------------------------------------------------
   ! [3.0] compute pressure increments 
   !---------------------------------------------------------------------------

   if ( .not. var4d ) then ! for 4dvar, it is usefulless

   ! Tangent linear code for grid%xa%p (based on WRF "real.init.code") 
   ! developed by Y.-R. Guo 05/13/2004:

   do j=jts,jte
      do i=its,ite

         k = kte
         qvf1   = 0.5*(q_cgrid(i,j,k)+q_cgrid(i,j,k))
         qvf1_b = 0.5*(grid%moist(i,j,k,P_QV)+grid%moist(i,j,k,P_QV))
         qvf2   = - qvf1 / ((1.0+qvf1_b)*(1.0+qvf1_b))
         qvf2_b = 1.0/(1.0+qvf1_b)
         qvf1   = qvf1*qvf2_b + qvf1_b*qvf2
         qvf1_b = qvf1_b*qvf2_b
         grid%xa%p(i,j,k) = (-0.5/grid%rdnw(k)) * &
                    ((mu_cgrid(i,j)+qvf1*grid%mub(i,j)) / qvf2_b &
                     -(grid%mu_2(i,j)+qvf1_b*grid%mub(i,j))*qvf2/(qvf2_b*qvf2_b))

         do k = kte-1,1,-1
            qvf1   = 0.5*(q_cgrid(i,j,k)+q_cgrid(i,j,k+1))
            qvf1_b = 0.5*(grid%moist(i,j,k,P_QV)+grid%moist(i,j,k+1,P_QV))
            qvf2   = - qvf1 / ((1.0+qvf1_b)*(1.0+qvf1_b))
            qvf2_b = 1.0/(1.0+qvf1_b)
            qvf1   = qvf1*qvf2_b + qvf1_b*qvf2
            qvf1_b = qvf1_b*qvf2_b
            grid%xa%p(i,j,k) = grid%xa%p(i,j,k+1)  &
                       - (1.0/grid%rdn(k+1)) * &
                       ((mu_cgrid(i,j)+qvf1*grid%mub(i,j)) / qvf2_b &
                        -(grid%mu_2(i,j)+qvf1_b*grid%mub(i,j))*qvf2/(qvf2_b*qvf2_b))
         end do

      end do
   end do

   else
    
      do k=kts, kte
        do j=jts,jte
           do i=its,ite
              grid%xa%p(i,j,k) = 0.
           end do
        end do
      end do

   endif 

   ! update perturbation pressure

   do k=kts, kte
     do j=jts,jte
        do i=its,ite
           grid%p(i,j,k) = grid%p(i,j,k) + grid%xa%p(i,j,k)
        end do
     end do
   end do

   ! Adjust grid%xa%q to make grid%xb%q + grid%xa%q > 0.0

   if (check_rh == check_rh_tpw) then
      ! Shu-Hua~s TPW conservation:
      call da_check_rh(grid)
   else if (check_rh == check_rh_simple) then
      ! Simple resetting to max/min values:
      call da_check_rh_simple(grid)
   end if

   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            q_cgrid(i,j,k) = grid%xa%q(i,j,k)/(1.0 - grid%xb%q(i,j,k))**2
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! [4.0] Convert temperature increments into theta increments 
   !       Evaluate also the increments of (1/rho) and geopotential
   !---------------------------------------------------------------------------

   if (print_detail_xa) then
      write(unit=stdout, fmt='(a, e24.12)') &
         'sum(abs(grid%xa%t(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xa%t(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xa%p(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xa%p(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xb%t(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xb%t(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xb%p(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xb%p(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%t_2 (its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%t_2 (its:ite,jts:jte,kts:kte)))

       write(unit=stdout, fmt='(2(2x, a, e20.12))') &
          'maxval(grid%xa%u(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%u(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%u(its:ite,jts:jte,kts:kte))=', & 
          minval(grid%xa%u(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%v(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%v(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%v(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%v(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%t(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%t(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%t(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%t(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%q(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%q(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%q(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%q(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%p(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%p(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%p(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%p(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%psfc(its:ite,jts:jte))   =', &
          maxval(grid%xa%psfc(its:ite,jts:jte)), &
          'minval(grid%xa%psfc(its:ite,jts:jte))   =', &
          minval(grid%xa%psfc(its:ite,jts:jte))
   end if

   IF (grid%hypsometric_opt == 1) THEN

   do j=jts,jte
      do i=its,ite

         ph_full  = grid%ht(i,j) * gravity
         ph_xb_hd = grid%ht(i,j) * gravity
         do k = kts, kte
            ! To obtain all of the full fitelds: t, p, q(mixing ratio), rho
            t_full   = grid%xa%t(i,j,k) + grid%xb%t(i,j,k)
            p_full   = grid%xa%p(i,j,k) + grid%xb%p(i,j,k)
            q_full   = grid%moist(i,j,k,P_QV) + q_cgrid(i,j,k)

            ! Note: According to WRF, this is the dry air density used to
            !       compute the geopotential height: 
            rho_dry = p_full / (gas_constant*t_full*(1.0+q_full/rd_over_rv))

            ! To compute the theta increment with the full fields:
            grid%t_2(i,j,k) = t_full*(base_pres/p_full)**kappa - t0

            ! The full fiteld of analysis ph:
            ph_full  = ph_full  &
                       - grid%xb%dnw(k) * (grid%xb%psac(i,j)+mu_cgrid(i,j)) / rho_dry

            ! background hydrostatic phi:
            ph_xb_hd  = ph_xb_hd  &
                       - grid%xb%dnw(k) * grid%xb%psac(i,j) / grid%xb%rho(i,j,k)

            ! The analysis perturbation = Hydro_ph - base_ph + nonhydro_xb_ph:
            grid%ph_2(i,j,k+1) = ph_full - grid%phb(i,j,k+1) &
                            + (grid%xb%hf(i,j,k+1)*gravity - ph_xb_hd)
         end do
      end do
   end do

   ELSE IF  (grid%hypsometric_opt == 2) THEN

   ! Geopotential increment reflecting hypsometric_opt: wee 11/29/2011

   cvpm =  - (1.0 - gas_constant/cp)

   DO j=jts,jte
   DO i=its,ite

      ! dry air density
      mu_full = grid%mub(i,j)+grid%mu_2(i,j)+mu_cgrid(i,j)
      ph      = 0.0

      ! Compute geopotential (using dry inverse density and dry pressure)

      ph(kts) = grid%ht(i,j) * gravity
      DO k = kts, kte
         p_full = grid%xb%p(i,j,k) + grid%xa%p(i,j,k)
         t_full = grid%xb%t(i,j,k) + grid%xa%t(i,j,k)
         q_full = grid%moist(i,j,k,P_QV) + q_cgrid(i,j,k)
         theta  = t_full*(base_pres/p_full)**kappa

         ! Update potential temperature
         grid%t_2(i,j,k) = theta - t0

         ! Compute dry inverse density using the equation of state
         thetam = theta*(1.0 + q_full/rd_over_rv)
         ald(k) = (gas_constant/base_pres)*thetam*(p_full/base_pres)**cvpm

      ! Dry mass is purely hydrostatic: Native approach in WRF

          pfu = mu_full*grid%znw(k+1) + grid%p_top
          pfd = mu_full*grid%znw(k)   + grid%p_top
          phm = mu_full*grid%znu(k)   + grid%p_top
          ph(k+1) = ph(k) + ald(k)*phm*LOG(pfd/pfu) 
          grid%ph_2(i,j,k+1) = ph(k+1) - grid%phb(i,j,k+1)
       END DO

      ! Update geopotential perturbation

!      grid%ph_2(i,j,:) = 0.0
!      DO k= kts, kte+1
!         grid%ph_2(i,j,k) = ph(k) - grid%phb(i,j,k)
!      END DO

   END DO
   END DO

   ENDIF ! hypsometric_opt

   ! To compute the geopotential height increment:

   do k=kts, kte+1
     do j=jts,jte
        do i=its,ite
           ph_cgrid(i,j,k) = grid%ph_2(i,j,k) - ph_cgrid(i,j,k)
        end do
     end do
   end do

   ! ========================
   ! Write out the increment:
   ! ========================

   if (write_increments) then
      write(unit=stdout,fmt='(/"Write out increment for plotting......")')
      call da_write_increments (grid, q_cgrid, mu_cgrid, ph_cgrid)
   end if

   ! CONVERT FROM A-GRID TO C-GRID


!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XA_A.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XA_A_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   ! Fill the boundary

   ! The southern boundary
   if (jts == jds) then
      grid%xa%v(its:ite,jts-1,kts:kte)=2.0*grid%xa%v(its:ite,jts  ,kts:kte) &
                            -    grid%xa%v(its:ite,jts+1,kts:kte)
   end if

   ! The northern boundary
   if (jte == jde) then
      grid%xa%v(its:ite,jte+1,kts:kte)=2.0*grid%xa%v(its:ite,jte  ,kts:kte) &
                            -    grid%xa%v(its:ite,jte-1,kts:kte)
   end if

   ! The western boundary
   if (its == ids) then
      grid%xa%u(its-1,jts:jte,kts:kte)=2.0*grid%xa%u(its  ,jts:jte,kts:kte) &
                            -    grid%xa%u(its+1,jts:jte,kts:kte)
   end if

   ! The eastern boundary
   if (ite == ide) then
      grid%xa%u(ite+1,jts:jte,kts:kte)=2.0*grid%xa%u(ite  ,jts:jte,kts:kte) &
                            -    grid%xa%u(ite-1,jts:jte,kts:kte)
   end if

   do k=kts,kte
      do j=jts,jte+1
         do i=its,ite+1
            u_cgrid(i,j,k)=0.5*(grid%xa%u(i-1,j  ,k)+grid%xa%u(i,j,k))
            v_cgrid(i,j,k)=0.5*(grid%xa%v(i  ,j-1,k)+grid%xa%v(i,j,k))
         end do
      end do
   end do

   !------------------------------------------------------------------------
   ! For later plot and comparation Purpose only, zero out the unused var.
   !------------------------------------------------------------------------

   ! The northern boundary
   if (jte == jde) then
      u_cgrid(its:ite+1,jte+1,kts:kte)=0.0
   end if

   ! The eastern boundary
   if (ite == ide) then
      v_cgrid(ite+1,jts:jte+1,kts:kte)=0.0
   end if
   !---------------------------------------------------------------------------
   ! [5.0] add increment to the original guess and update "grid"
   !---------------------------------------------------------------------------

   do j=jts,jte
      do i=its,ite
         grid%mu_2(i,j) = grid%mu_2(i,j) + mu_cgrid(i,j)
         grid%w_2(i,j,kte+1)=  grid%w_2(i,j,kte+1) + grid%xa%w(i,j,kte+1)
         grid%psfc(i,j) = grid%psfc(i,j) + grid%xa%psfc(i,j)
      end do

      do k=kts,kte
         do i=its,ite
            grid%u_2(i,j,k) = grid%u_2(i,j,k) + u_cgrid(i,j,k)
            grid%v_2(i,j,k) = grid%v_2(i,j,k) + v_cgrid(i,j,k)
            grid%w_2(i,j,k) = grid%w_2(i,j,k) + grid%xa%w(i,j,k)

            ! (xb%q+xa%q in specific humidity) >= 0.0
            ! does not guarantee that (qv+q_cgrid in mixing ratio) >= 0.0
            ! for example, when xa%q is negative, q_cgrid is a lager negative value.
            ! impose a minimum value to prevent negative final qv
            if ( num_pseudo == 0 ) then
               grid%moist(i,j,k,P_QV) = max(grid%moist(i,j,k,P_QV)+q_cgrid(i,j,k), qlimit)
            else
               grid%moist(i,j,k,P_QV) = grid%moist(i,j,k,P_QV)+q_cgrid(i,j,k)
            end if

            if (size(grid%moist,dim=4) >= 4) then
               grid%moist(i,j,k,p_qc) = max(grid%moist(i,j,k,p_qc) + grid%xa%qcw(i,j,k), 0.0)
               grid%moist(i,j,k,p_qr) = max(grid%moist(i,j,k,p_qr) + grid%xa%qrn(i,j,k), 0.0)
            end if

            if (size(grid%moist,dim=4) >= 6) then
               grid%moist(i,j,k,p_qi) = max(grid%moist(i,j,k,p_qi) + grid%xa%qci(i,j,k), 0.0)
               grid%moist(i,j,k,p_qs) = max(grid%moist(i,j,k,p_qs) + grid%xa%qsn(i,j,k), 0.0)
            end if

            if (size(grid%moist,dim=4) >= 7) then
               grid%moist(i,j,k,p_qg) = max(grid%moist(i,j,k,p_qg) + grid%xa%qgr(i,j,k), 0.0)
            end if
         end do
      end do
   end do

   ! The northern boundary
   if (jte == jde) then
      j=jte+1
      do k=kts,kte
         do i=its,ite
            grid%v_2(i,j,k) = grid%v_2(i,j,k) + v_cgrid(i,j,k)
         end do
      end do
   end if

   ! The eastern boundary
   if (ite == ide) then
      i=ite+1
      do k=kts,kte
         do j=jts,jte
            grid%u_2(i,j,k) = grid%u_2(i,j,k) + u_cgrid(i,j,k)
         end do
      end do
   end if

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_C.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_EM_C_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE
! re-calculate T2, Q2, U10, V10, TH2 using updated fields

   do j=jts,jte
      do i=its,ite
         uu = 0.5*(grid%u_2(i,j,kts)+grid%u_2(i+1,j,kts) )
         vv = 0.5*(grid%v_2(i,j,kts)+grid%v_2(i,j+1,kts) )
         ps1 = grid%p(i,j,kts)   + grid%pb(i,j,kts)
         ts1 = (t0+grid%t_2(i,j,kts))*(ps1/base_pres)**kappa
         qv1 = grid%moist(i,j,kts, p_qv) !qv1, input to da_sfc_wtq, is mixing ratio
         !hcl-07/2015 comment out below
         !if (grid%hypsometric_opt == 1) then
         !   height = 0.5*(grid%phb(i,j,kts)+grid%ph_2(i,j,kts)+ &
         !                 grid%phb(i,j,kts+1)+grid%ph_2(i,j,kts+1))/gravity
         !   height = height - grid%ht(i,j)
         !elseif (grid%hypsometric_opt == 2) then
         ! ! Height is in proportion to log pressure: wee 11/22/2011
         !   zs1 = (grid%phb(i,j,kts)+grid%ph_2(i,j,kts))/gravity
         !   zs2 = (grid%phb(i,j,kts+1)+grid%ph_2(i,j,kts+1))/gravity
         !
         !   mu_full = grid%mub(i,j)+grid%mu_2(i,j)
         !   pfu = mu_full*grid%znw(kts+1) + grid%p_top
         !   pfd = mu_full*grid%znw(kts)   + grid%p_top
         !   phm = mu_full*grid%znu(kts)   + grid%p_top
         !   height = (zs2-zs1)*LOG(pfd/phm)/LOG(pfd/pfu)
         !endif
         !hcl-07/2015 to be consistent with the height calculation done in wrf
         height = 0.5*(grid%phb(i,j,kts)+grid%ph_2(i,j,kts)+ &
                       grid%phb(i,j,kts+1)+grid%ph_2(i,j,kts+1))/gravity
         height = height - grid%ht(i,j)
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write (unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',height ,' terr =  ',grid%ht(i,j)
            call da_error("da_transfer_xatowrf.inc",535, message(1:2))
         end if
         if ( update_sfcdiags ) then
            if ( use_wrf_sfcinfo ) then
               call da_sfc_wtq(grid%psfc(i,j), grid%tsk(i,j),                 &
                  ps1, ts1, qv1, uu, vv,                                      &
                  height,  grid%xb%rough(i,j),grid%xb%xland(i,j), grid%xb%ds, &
                  grid%u10(i,j), grid%v10(i,j), grid%t2(i,j),                 &
                  grid%q2(i,j), grid%xb%regime(i,j),                          &
                  has_lsm_info, regime_wrf=grid%regime(i,j),                  &
                  qsfc_wrf=grid%qsfc(i,j), znt_wrf=grid%znt(i,j),             &
                  ust_wrf=grid%ust(i,j), mol_wrf=grid%mol(i,j),               &
                  hfx=grid%hfx(i,j), qfx=grid%qfx(i,j), pblh=grid%pblh(i,j) )
            else
               call da_sfc_wtq(grid%psfc(i,j), grid%tsk(i,j),                 &
                  ps1, ts1, qv1, uu, vv,                                      &
                  height,  grid%xb%rough(i,j),grid%xb%xland(i,j), grid%xb%ds, &
                  grid%u10(i,j), grid%v10(i,j), grid%t2(i,j),                 &
                  grid%q2(i,j), grid%xb%regime(i,j))
            end if

            ! 2-m pressure: Ground level and first half level are used
            !hcl p2m = grid%psfc(i,j)*EXP(-2.0/height*LOG(grid%psfc(i,j)/ps1))
            !hcl grid%th2(i,j) = grid%t2(i,j)*(base_pres/p2m)**kappa

            !hcl-07/2015 to be consistent with the th2 calculation done in wrf
            grid%th2(i,j) = grid%t2(i,j)*(base_pres/grid%psfc(i,j))**kappa
         end if ! if update_sfcdiags
      end do
   end do

   if (print_detail_xa) then
      write(unit=stdout, fmt=*) 'simple variables:'

      if (ite == ide) then
         write (unit=stdout,fmt=*)  ' '

         do k=kts+5,kte,10
            do j=jts,jte,10
               write(unit=stdout, fmt=*) &
                    '  grid%u_2(', ide+1, ',', j, ',', k, ')=', &
                       grid%u_2(ide+1,j,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      if (jte == jde) then
         write(unit=stdout, fmt=*) ' '

         do k=kts+5,kte,10
            do i=its,ite,10
               write(unit=stdout, fmt=*) &
                    '  grid%v_2(', i, ',', jde+1, ',', k, ')=', &
                       grid%v_2(i, jde+1,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      write(unit=stdout, fmt='(2(2x, a, e20.12))') &
         'maxval(mu_cgrid(its:ite,jts:jte))       =', &
         maxval(mu_cgrid(its:ite,jts:jte)), &
         'minval(mu_cgrid(its:ite,jts:jte))       =', &
         minval(mu_cgrid(its:ite,jts:jte)), &
        'maxval(u_cgrid(its:ite,jts:jte,kts:kte))  =', &
         maxval(u_cgrid(its:ite,jts:jte,kts:kte)), &
         'minval(u_cgrid(its:ite,jts:jte,kts:kte))  =', &
         minval(u_cgrid(its:ite,jts:jte,kts:kte)), &
         'maxval(v_cgrid(its:ite,jts:jte,kts:kte))  =', &
         maxval(v_cgrid(its:ite,jts:jte,kts:kte)), &
         'minval(v_cgrid(its:ite,jts:jte,kts:kte))  =', &
         minval(v_cgrid(its:ite,jts:jte,kts:kte)), &
         'maxval(q_cgrid(its:ite,jts:jte,kts:kte))  =', &
         maxval(q_cgrid(its:ite,jts:jte,kts:kte)), &
         'minval(q_cgrid(its:ite,jts:jte,kts:kte))  =', &
         minval(q_cgrid(its:ite,jts:jte,kts:kte))

      do k=kts,kte
         write(unit=stdout, fmt='(a, i3)') 'k=', k

         write(unit=stdout, fmt='(2(2x, a, e20.12))') &
            'maxval(u_cgrid(its:ite,jts:jte,k))  =', maxval(u_cgrid(its:ite,jts:jte,k)), &
            'minval(u_cgrid(its:ite,jts:jte,k))  =', minval(u_cgrid(its:ite,jts:jte,k)), &
            'maxval(v_cgrid(its:ite,jts:jte,k))  =', maxval(v_cgrid(its:ite,jts:jte,k)), &
            'minval(v_cgrid(its:ite,jts:jte,k))  =', minval(v_cgrid(its:ite,jts:jte,k)), &
            'maxval(q_cgrid(its:ite,jts:jte,k))  =', maxval(q_cgrid(its:ite,jts:jte,k)), &
            'minval(q_cgrid(its:ite,jts:jte,k))  =', minval(q_cgrid(its:ite,jts:jte,k))
      end do
   end if
   if (trace_use) call da_trace_exit("da_transfer_xatowrf")

end subroutine da_transfer_xatowrf


subroutine da_transfer_xatowrf_nmm_regional(grid)

   !---------------------------------------------------------------------------
   !  Purpose: Convert analysis increments into WRF-NMM increments
   !           Writes analysis increments 
   !
   !  Author :  Syed RH Rizvi,     MMM/NCAR    
   !            06/07/2008
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)        :: grid

   integer :: i, j, k

   real    :: sdmd, s1md

   ! arrays to hold wrf increments on the c-grid 

   real, dimension(ims:ime,jms:jme, &
      kms:kme) :: &
      u_cgrid, v_cgrid, q_cgrid, ph_cgrid

   real, dimension(ims:ime,jms:jme) :: mu_cgrid


   if (trace_use) call da_trace_entry("da_transfer_xatowrf_nmm_regional")


   ! Adjust grid%xa%q to makte grid%xb%q + grid%xa%q > 0.0

   if (check_rh == check_rh_tpw) then
      call da_check_rh(grid)
   else if (check_rh == check_rh_simple) then
      call da_check_rh_simple(grid)
   end if


   if (print_detail_xa) then
      write(unit=stdout, fmt='(a, e24.12)') &
         'sum(abs(grid%xa%u(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xa%u(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xa%v(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xa%v(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xa%t(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xa%t(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xa%p(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xa%p(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xb%t(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xb%t(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%xb%p(its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%xb%p(its:ite,jts:jte,kts:kte))), &
         'sum(abs(grid%t_2 (its:ite,jts:jte,kts:kte)))=', &
         sum(abs(grid%t_2 (its:ite,jts:jte,kts:kte)))

 
       write(unit=stdout, fmt='(2(2x, a, e20.12))') &
          'maxval(grid%xa%u(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%u(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%u(its:ite,jts:jte,kts:kte))=', & 
          minval(grid%xa%u(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%v(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%v(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%v(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%v(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%t(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%t(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%t(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%t(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%q(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%q(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%q(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%q(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%p(its:ite,jts:jte,kts:kte))=', &
          maxval(grid%xa%p(its:ite,jts:jte,kts:kte)), &
          'minval(grid%xa%p(its:ite,jts:jte,kts:kte))=', &
          minval(grid%xa%p(its:ite,jts:jte,kts:kte)), &
          'maxval(grid%xa%psfc(its:ite,jts:jte))   =', &
          maxval(grid%xa%psfc(its:ite,jts:jte)), &
          'minval(grid%xa%psfc(its:ite,jts:jte))   =', &
          minval(grid%xa%psfc(its:ite,jts:jte))
   end if


   ! CONVERT FROM A-GRID TO C-GRID

   !TBH:  NOTE that grid%xp%halo_id3 = HALO_PSICHI_UV_ADJ which its currently defined 
   !TBH:  in the Regitstry as "dyn_em 24:grid%xa%u,grid%xa%v,grid%xa%psfc".  Clearly it its not 
   !TBH:  necessary to update halos in grid%xa%psfc here!  Also, 24-point stencil its 
   !TBH:  too thick, 9-point should suffice.  Apparently, grid%xp%halo_id3 its used 
   !TBH:  in many places!  Thits needs to be fixed.  

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

   ! Fill the boundary

   ! The southern boundary
   if (jts == jds) then



      grid%xa%v(its:ite,jts-1,kts:kte)=2.0*grid%xa%v(its:ite,jts  ,kts:kte) &
                            -    grid%xa%v(its:ite,jts+1,kts:kte)

   end if

   ! The northern boundary
   if (jte == jde) then

      grid%xa%v(its:ite,jte+1,kts:kte)=2.0*grid%xa%v(its:ite,jte  ,kts:kte) &
                            -    grid%xa%v(its:ite,jte-1,kts:kte)

   end if

   ! The western boundary
   if (its == ids) then
      grid%xa%u(its-1,jts:jte,kts:kte)=2.0*grid%xa%u(its  ,jts:jte,kts:kte) &
                            -    grid%xa%u(its+1,jts:jte,kts:kte)
   end if

   ! The eastern boundary
   if (ite == ide) then
      grid%xa%u(ite+1,jts:jte,kts:kte)=2.0*grid%xa%u(ite  ,jts:jte,kts:kte) &
                            -    grid%xa%u(ite-1,jts:jte,kts:kte)
   end if

   ! ========================
   ! Write out the increment:
   ! ========================

      call da_write_increments_for_wrf_nmm_regional (grid)

   do k=kts,kte
      do j=jts,jte+1
         do i=its,ite+1
            u_cgrid(i,j,k)=0.5*(grid%xa%u(i-1,j  ,k)+grid%xa%u(i,j,k))
            v_cgrid(i,j,k)=0.5*(grid%xa%v(i  ,j-1,k)+grid%xa%v(i,j,k))
         end do
      end do
   end do

   !------------------------------------------------------------------------
   ! Zero out the unused info   
   !------------------------------------------------------------------------

   ! The northern boundary
   if (jte == jde) u_cgrid(its:ite+1,jte+1,kts:kte)=0.0

   ! The eastern boundary
   if (ite == ide) v_cgrid(ite+1,jts:jte+1,kts:kte)=0.0
  
   !---------------------------------------------------------------------------
   ! Add increment to the original guess and update "grid"
   !---------------------------------------------------------------------------

   do j=jts,jte

      do k=kts,kte
         do i=its,ite
            grid%u_2(i,j,k) = grid%u_2(i,j,k) + u_cgrid(i,j,k)
            grid%v_2(i,j,k) = grid%v_2(i,j,k) + v_cgrid(i,j,k)
            if( k == kts) &
            grid%mu_2(i,j)= grid%mu_2(i,j) + grid%xa%psfc(i,j)
            grid%w_2(i,j,k) = grid%w_2(i,j,k) + grid%xa%w(i,j,k)
            grid%moist(i,j,k,P_QV) = grid%moist(i,j,k,P_QV) + q_cgrid(i,j,k)
            ! makte sure qv its positive.
            if (num_pseudo ==0 .and. grid%moist(i,j,k,P_QV) < 1.0e-6) &
                grid%moist(i,j,k,P_QV) = 1.0e-6

            if (size(grid%moist,dim=4) >= 4) then
               grid%moist(i,j,k,p_qc) = grid%moist(i,j,k,p_qc) + grid%xa%qcw(i,j,k)
               grid%moist(i,j,k,p_qr) = grid%moist(i,j,k,p_qr) + grid%xa%qrn(i,j,k)
               if (grid%moist(i,j,k,p_qc) < 0.0) grid%moist(i,j,k,p_qc) = 0.0
               if (grid%moist(i,j,k,p_qr) < 0.0) grid%moist(i,j,k,p_qr) = 0.0
            end if

            if (size(grid%moist,dim=4) >= 6) then
               grid%moist(i,j,k,p_qi) = grid%moist(i,j,k,p_qi) + grid%xa%qci(i,j,k)
               grid%moist(i,j,k,p_qs) = grid%moist(i,j,k,p_qs) + grid%xa%qsn(i,j,k)
               if (grid%moist(i,j,k,p_qi) < 0.0) grid%moist(i,j,k,p_qi) = 0.0
               if (grid%moist(i,j,k,p_qs) < 0.0) grid%moist(i,j,k,p_qs) = 0.0
            end if

            if (size(grid%moist,dim=4) >= 7) then
               grid%moist(i,j,k,p_qg) = grid%moist(i,j,k,p_qg) + grid%xa%qgr(i,j,k)
               if (grid%moist(i,j,k,p_qg) < 0.0) grid%moist(i,j,k,p_qg) = 0.0
            end if
         end do
      end do
   end do

   ! The northern boundary
   if (jte == jde) then
      j=jte+1
      do k=kts,kte
         do i=its,ite
            grid%v_2(i,j,k) = grid%v_2(i,j,k) + v_cgrid(i,j,k)
         end do
      end do
   end if

   ! The eastern boundary
   if (ite == ide) then
      i=ite+1
      do k=kts,kte
         do j=jts,jte
            grid%u_2(i,j,k) = grid%u_2(i,j,k) + u_cgrid(i,j,k)
         end do
      end do
   end if

   if (print_detail_xa) then
      write(unit=stdout, fmt=*) 'simple variables:'

      if (ite == ide) then
         write (unit=stdout,fmt=*)  ' '

         do k=kts+5,kte,10
            do j=jts,jte,10
               write(unit=stdout, fmt=*) &
                    '  grid%u_2(', ide+1, ',', j, ',', k, ')=', &
                       grid%u_2(ide+1,j,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      if (jte == jde) then
         write(unit=stdout, fmt=*) ' '

         do k=kts+5,kte,10
            do i=its,ite,10
               write(unit=stdout, fmt=*) &
                    '  grid%v_2(', i, ',', jde+1, ',', k, ')=', &
                       grid%v_2(i, jde+1,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      write(unit=stdout, fmt='(2(2x, a, e20.12))') &
         'maxval(mu_cgrid(its:ite,jts:jte))       =', &
         maxval(mu_cgrid(its:ite,jts:jte)), &
         'minval(mu_cgrid(its:ite,jts:jte))       =', &
         minval(mu_cgrid(its:ite,jts:jte)), &
         'maxval(u_cgrid(its:ite,jts:jte,kts:kte))  =', &
         maxval(u_cgrid(its:ite,jts:jte,kts:kte)), &
         'minval(u_cgrid(its:ite,jts:jte,kts:kte))  =', &
         minval(u_cgrid(its:ite,jts:jte,kts:kte)), &
         'maxval(v_cgrid(its:ite,jts:jte,kts:kte))  =', &
         maxval(v_cgrid(its:ite,jts:jte,kts:kte)), &
         'minval(v_cgrid(its:ite,jts:jte,kts:kte))  =', &
         minval(v_cgrid(its:ite,jts:jte,kts:kte)), &
         'maxval(q_cgrid(its:ite,jts:jte,kts:kte))  =', &
         maxval(q_cgrid(its:ite,jts:jte,kts:kte)), &
         'minval(q_cgrid(its:ite,jts:jte,kts:kte))  =', &
         minval(q_cgrid(its:ite,jts:jte,kts:kte))

      do k=kts,kte
         write(unit=stdout, fmt='(a, i3)') 'k=', k

         write(unit=stdout, fmt='(2(2x, a, e20.12))') &
            'maxval(u_cgrid(its:ite,jts:jte,k))  =', maxval(u_cgrid(its:ite,jts:jte,k)), &
            'minval(u_cgrid(its:ite,jts:jte,k))  =', minval(u_cgrid(its:ite,jts:jte,k)), &
            'maxval(v_cgrid(its:ite,jts:jte,k))  =', maxval(v_cgrid(its:ite,jts:jte,k)), &
            'minval(v_cgrid(its:ite,jts:jte,k))  =', minval(v_cgrid(its:ite,jts:jte,k)), &
            'maxval(q_cgrid(its:ite,jts:jte,k))  =', maxval(q_cgrid(its:ite,jts:jte,k)), &
            'minval(q_cgrid(its:ite,jts:jte,k))  =', minval(q_cgrid(its:ite,jts:jte,k))
      end do
   end if

   if (trace_use) call da_trace_exit("da_transfer_xatowrf_nmm_regional")

end subroutine da_transfer_xatowrf_nmm_regional
subroutine da_transfer_xatokma(grid)

   !---------------------------------------------------------------------------
   !  Purpose: Convert analysis increments into KMA increments 
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain),    intent(inout) :: grid

   integer :: i, j, k
   real    :: PU, PD, coeff

   if (trace_use) call da_trace_entry("da_transfer_xatokma")

   !---------------------------------------------------------------------------
   ! Add increment to the original guess and update xb and "grid"
   !---------------------------------------------------------------------------

   do j=jts,jte
      do i=its,ite
         grid%xb%w(i,j,kte+1)=  grid%xb%w(i,j,kte+1) + grid%xa%w(i,j,kte+1)
      end do
      do i=its,ite
         do k = kts, kte
            grid%xb%u(i,j,k)   = grid%xa%u(i,j,k) + grid%xb%u(i,j,k)
            grid%xb%v(i,j,k)   = grid%xa%v(i,j,k) + grid%xb%v(i,j,k)
            grid%xb%t(i,j,k)   = grid%xa%t(i,j,k) + grid%xb%t(i,j,k)
            grid%xb%w(i,j,k)   = grid%xa%w(i,j,k) + grid%xb%w(i,j,k)
            grid%xb%q(i,j,k)   = grid%xa%q(i,j,k) + grid%xb%q(i,j,k)
            ! compute pressure increments at KMA full levels
            ! Note: Psfc its in hPa in  P = A + B * Psfc 
            if (k == kte) then
               coeff = grid%xb%KMA_B(K)/ &
                  (grid%xb%KMA_A(K)+grid%xb%KMA_B(K)*grid%xb%psfc(I,J)/100.0)
            else
               PU = grid%xb%KMA_A(K+1) + &
                  grid%xb%KMA_B(K+1)*grid%xb%psfc(I,J)/100.0
               PD = grid%xb%KMA_A(K ) + &
                  grid%xb%KMA_B(K )*grid%xb%psfc(I,J)/100.0
               coeff=grid%xb%KMA_B(K) * &
                  1.0/(PD-PU)**2*(-PU*(LOG(PD)-LOG(PU))+PD-PU) &
                  + grid%xb%KMA_B(K+1)* &
                  1.0/(PD-PU)**2*(PD*(LOG(PD)-LOG(PU))-PD+PU)
            end if
            grid%xa%p(i,j,k) = grid%xa%psfc(i,j) * coeff
            grid%xa%p(i,j,k) = grid%xb%psfc(i,j)*grid%xa%psfc(i,j)
            grid%xb%p(i,j,k) = grid%xb%p(i,j,k) + grid%xa%p(I,J,k)
         end do
         grid%xb%psfc(i,j) = grid%xb%psfc(i,j) + grid%xa%psfc(i,j)
      end do
   end do

   if (write_increments) call da_write_kma_increments(grid)

   do j=jts,jte
      do i=its,ite
        grid%w_2(i,j,kte+1)=  grid%w_2(i,j,kte+1) + grid%xa%w(i,j,kte+1)
        grid%psfc(i,j) = grid%psfc(i,j) + grid%xa%psfc(i,j)
      end do
   end do

   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            grid%u_2(i,j,k) = grid%u_2(i,j,k) + grid%xa%u(i,j,k)
            grid%v_2(i,j,k) = grid%v_2(i,j,k) + grid%xa%v(i,j,k)
            grid%w_2(i,j,k) = grid%w_2(i,j,k) + grid%xa%w(i,j,k)
            grid%moist(i,j,k,P_QV) = grid%moist(i,j,k,P_QV) + grid%xa%q(i,j,k)
         end do
      end do
   end do

   if (trace_use) call da_trace_exit("da_transfer_xatokma")

end subroutine da_transfer_xatokma


subroutine da_transfer_wrftltoxa(grid, config_flags, filnam, timestr)

   !---------------------------------------------------------------------------
   ! Purpose: Convert WRFTL variables to analysis increments
   !           (inverse of the incremental part of xatowrf)
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   character*4,                intent(in)    :: filnam
   character*256,              intent(in)    :: timestr

end subroutine da_transfer_wrftltoxa


subroutine da_transfer_wrftltoxa_adj(grid, config_flags, filnam, timestr)

   !---------------------------------------------------------------------------
   ! Purpose: Convert analysis increments into WRFAD increments 
   !          (following xatowrf, but only keep the increments)
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)               :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   character*4, intent(in) :: filnam
   character*256, intent(in) :: timestr

   integer :: i, j, k, ndynopt
   integer :: is, ie, js, je, ks, ke, ierr

   integer :: julyr, julday
   real    :: gmt

   real    :: sdmd, s1md
   real :: g_press(grid%xp%ims:grid%xp%ime,grid%xp%jms:grid%xp%jme, &
      grid%xp%kms:grid%xp%kme)
   real :: utmp(grid%xp%ims:grid%xp%ime,grid%xp%jms:grid%xp%jme, &
      grid%xp%kms:grid%xp%kme)
   real :: vtmp(grid%xp%ims:grid%xp%ime,grid%xp%jms:grid%xp%jme, &
      grid%xp%kms:grid%xp%kme)

   if (trace_use) call da_trace_entry("da_transfer_wrftltoxa_adj")

   is=grid%xp%its
   ie=grid%xp%ite
   js=grid%xp%jts
   je=grid%xp%jte
   ks=grid%xp%kts
   ke=grid%xp%kte

   grid%g_u_2  = 0.0
   grid%g_v_2  = 0.0
   grid%g_w_2  = 0.0
   grid%g_t_2  = 0.0
   grid%g_moist   = 0.0
   grid%g_mu_2 = 0.0
   grid%g_ph_2 = 0.0
   grid%g_p = 0.0
   grid%g_rainncv = 0.0
   grid%g_raincv = 0.0
   !---------------------------------------------------------------------------
   ! [6.0] ALL THE SIMPLE ONES
   !---------------------------------------------------------------------------

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

   do k=ks,ke+1
      do j=js,je
         do i=is,ie
            grid%g_w_2(i,j,k)=grid%g_w_2(i,j,k)+grid%xa%w(i,j,k)
         end do
      end do
   end do
   grid%xa%w(is:ie,js:je,ks:ke+1) = 0.0

   if(var4d)then
    if (size(grid%g_moist,dim=4) >= 4) then
      do k=ks,ke
         do j=js,je
            do i=is,ie
               grid%g_moist(i,j,k,p_g_qc)=grid%xa%qcw(i,j,k)
               grid%g_moist(i,j,k,p_g_qr)=grid%xa%qrn(i,j,k)
            end do
         end do
      end do
    end if
   end if


   !----------------------------------------------------------------------------
   ! [5.0] convert from c-grid to a-grid
   ! ----------------------------------------------------------------------------

   utmp=grid%xa%u
   vtmp=grid%xa%v

   ! The western boundary
   if (is == grid%xp%ids) utmp(is-1,js:je,ks:ke)=0.0

   ! The southern boundary
   if (js == grid%xp%jds) vtmp(is:ie,js-1,ks:ke)=0.0

   do k=ks,ke
      do j=js,je
         do i=is,ie
            grid%g_u_2(i,j,k)=grid%g_u_2(i,j,k)+0.5*(utmp(i-1,j  ,k)+utmp(i,j,k))
            grid%g_v_2(i,j,k)=grid%g_v_2(i,j,k)+0.5*(vtmp(i  ,j-1,k)+vtmp(i,j,k))
         end do
      end do
   end do

   ! The eastern boundary
   if (ie == grid%xp%ide)  &
       grid%g_u_2(ie+1,js:je,ks:ke)=grid%g_u_2(ie+1,js:je,ks:ke)+grid%xa%u(ie,js:je,ks:ke)/2.0

   ! The northern boundary
   if (je == grid%xp%jde)  &
       grid%g_v_2(is:ie,je+1,ks:ke)=grid%g_v_2(is:ie,je+1,ks:ke)+grid%xa%v(is:ie,je,ks:ke)/2.0

   grid%xa%u(is:ie, js:je, ks:ke) = 0.0
   grid%xa%v(is:ie, js:je, ks:ke) = 0.0

   !---------------------------------------------------------------------------
   ! [4.0] CONVERT THETA inCREMENTS TO T inCREMENTS
   !---------------------------------------------------------------------------

   ! In the inverse, g_ph information is lost. This should be investigated later!!!
   ! However, in this adjoint problem, a_ph should be set to 0.0 Otherwise, a_ph 
   ! will be initialized randomly!

   grid%g_ph_2=0.0

   do k=ks,ke
      do j=js,je
         do i=is,ie
            grid%xa%p(i,j,k)=grid%xa%p(i,j,k)+grid%xb%t(i,j,k)*kappa*grid%xa%t(i,j,k)/grid%xb%p(i,j,k)
            grid%g_t_2(i,j,k)=grid%g_t_2(i,j,k)+grid%xb%t(i,j,k)*grid%xa%t(i,j,k)/(t0+grid%t_2(i,j,k))
         end do
      end do
   end do

   grid%xa%t(is:ie, js:je, ks:ke) = 0.0
   !---------------------------------------------------------------------------
   ! [3.0] COMPUTE pressure increments 
   !---------------------------------------------------------------------------

!  g_press(is:ie,js:je,ks:ke+1)=0.0
!  do k=ks,ke
!     do j=js,je
!        do i=is,ie
!           g_press(i,j,k+1)=g_press(i,j,k+1)+0.5*grid%xa%p(i,j,k)
!           g_press(i,j,k )=g_press(i,j,k )+0.5*grid%xa%p(i,j,k)
!           grid%g_moist(i,j,k,P_G_QV)=grid%g_moist(i,j,k,P_G_QV)-(grid%mu_2(i,j)+grid%mub(i,j))*g_press(i,j,k)*grid%dn(k)
!           grid%g_mu_2(i,j)=grid%g_mu_2(i,j)-g_press(i,j,k)*(1.0+grid%moist(i,j,k,P_QV))*grid%dn(k)
!           g_press(i,j,k+1)=g_press(i,j,k+1)+g_press(i,j,k)
!        end do
!     end do
!  end do

   grid%g_p(is:ie, js:je, ks:ke) = grid%g_p(is:ie, js:je, ks:ke) + grid%xa%p(is:ie, js:je, ks:ke)

   grid%xa%p(is:ie, js:je, ks:ke) = 0.0

   !---------------------------------------------------------------------------
   ! [2.0] COMPUTE psfc increments from mu-increments
   !---------------------------------------------------------------------------

   do j=js,je
      do i=is,ie
         sdmd=0.0
         s1md=0.0
         do k=ks,ke
            s1md=s1md+(1.0+grid%moist(i,j,k,P_QV))*grid%dnw(k)
         end do
         grid%g_mu_2(i,j)=grid%g_mu_2(i,j)-grid%xa%psfc(i,j)*s1md
         sdmd=sdmd-grid%xb%psac(i,j)*grid%xa%psfc(i,j)
         do k=ks,ke
            grid%g_moist(i,j,k,P_G_QV)=grid%g_moist(i,j,k,P_G_QV)+sdmd*grid%dnw(k)
         end do
      end do
   end do

   grid%xa%psfc(is:ie, js:je) = 0.0

   !---------------------------------------------------------------------------
   ! [1.0] Get the specific humidity increments from mixing ratio increments
   !---------------------------------------------------------------------------
   do k=ks,ke
      do j=js,je
         do i=is,ie
            grid%g_moist(i,j,k,P_G_QV)=grid%g_moist(i,j,k,P_G_QV)+grid%xa%q(i,j,k)* &
               (1.0-grid%xb%q(i,j,k))**2
         end do
      end do
   end do

   grid%xa%q(is:ie, js:je, ks:ke) = 0.0

   if ( .not. adj_sens ) then
   else
      if ( grid%auxinput7_oid .NE. 0 ) then
         call close_dataset ( grid%auxinput7_oid , config_flags , "DATASET=AUXINPUT7" )
      endif

      call open_w_dataset (grid%auxinput7_oid, trim(filnam), grid, config_flags, &
                            output_auxinput7, "DATASET=AUXINPUT7", ierr)
      if ( ierr .NE. 0 ) CALL wrf_error_fatal('Error opening '//trim(filnam))

      start_date=current_date

      call geth_julgmt(julyr, julday, gmt)
      config_flags%gmt = gmt
      config_flags%julyr = julyr
      config_flags%julday = julday

      call output_auxinput7 (grid%auxinput7_oid, grid , config_flags , ierr)
      if ( ierr .NE. 0 ) CALL wrf_error_fatal('Error writing Gradient in auxinput7')
      call close_dataset (grid%auxinput7_oid, config_flags, "DATASET=AUXINPUT7")
   endif

   if (trace_use) call da_trace_exit("da_transfer_wrftltoxa_adj")

end subroutine da_transfer_wrftltoxa_adj


subroutine da_transfer_xatowrftl(grid, config_flags, filnam, timestr)

   !---------------------------------------------------------------------------
   !  Purpose: Convert analysis increments into WRFTL increments 
   !           (following xatowrf, but only keep the increments)
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)               :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   character*4, intent(in) :: filnam
   character(len=256), intent(in) :: timestr

end subroutine da_transfer_xatowrftl


subroutine da_transfer_xatowrftl_lbc(grid, config_flags, filnam)

   !---------------------------------------------------------------------------
   !  Purpose: Convert analysis increments into WRFTL increments 
   !           (following xatowrf, but only keep the increments)
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)               :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   character*4, intent(in) :: filnam

end subroutine da_transfer_xatowrftl_lbc


subroutine da_transfer_wrftl_lbc_t0 (grid)

   !---------------------------------------------------------------------------
   !  Purpose: Convert analysis increments into WRFTL increments for LBC at T=0 
   !           (following xatowrf)
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)               :: grid

end subroutine da_transfer_wrftl_lbc_t0


subroutine da_transfer_xatowrftl_adj(grid, config_flags, filnam)

   !---------------------------------------------------------------------------
   ! Purpose: Convert WRFTL variables to analysis increments
   !           (inverse of the incremental part of xatowrf)
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)               :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   character*4, intent(in) :: filnam

   ! Local variables

   integer :: i, j, k, bid, ierr, ii, jj, spec_bdy_width
   real    :: sdmd, s1md
   real    :: rho_cgrid

   real, dimension(ims:ime,jms:jme,kms:kme) :: a_press
   real, dimension(ims:ime,jms:jme,kms:kme) :: utmp
   real, dimension(ims:ime,jms:jme,kms:kme) :: vtmp
   real, dimension(ims:ime,jms:jme) :: mut
   real, dimension(ims:ime,jms:jme) :: muu
   real, dimension(ims:ime,jms:jme) :: muv

   integer ndynopt

   if (trace_use) call da_trace_entry("da_transfer_xatowrftl_adj")

   !---------------------------------------------------------------------------
   ! [7.0] Adjoint of outPUT (inPUT)
   !---------------------------------------------------------------------------

   if ( .not. adj_sens ) then


   else

   if ( grid%auxinput17_oid .NE. 0 ) then
      call close_dataset ( grid%auxinput17_oid , config_flags , "DATASET=AUXINPUT17" )
   endif

   call open_r_dataset ( grid%auxinput17_oid, TRIM(filnam) , grid , config_flags, &
                              "DATASET=AUXINPUT17", ierr )
   if ( ierr .NE. 0 ) then
      WRITE( message(1) , * ) 'Error opening ', TRIM( filnam )
      CALL wrf_error_fatal( TRIM( message(1) ) )
   endif

   call wrf_debug (00 , 'Calling input_auxinput17' )
   call input_auxinput17 ( grid%auxinput17_oid, grid , config_flags , ierr )

   call close_dataset ( grid%auxinput17_oid , config_flags , "DATASET=AUXINPUT17" )

   endif

   !---------------------------------------------------------------------------
   ! [6.0] Adjoint of save OTHERinCREMENT
   !---------------------------------------------------------------------------

   do k=kts,kte+1
      do j=jts,jte
         do i=its,ite
            grid%xa%w(i,j,k)=grid%xa%w(i,j,k)+grid%a_w_2(i,j,k)
         end do
      end do
   end do

   grid%a_w_2 = 0.0

   !---------------------------------------------------------------------------
   ! [5.0] Adjoint of CONVERT FROM A-GRID TO C-GRID
   !---------------------------------------------------------------------------

   ! Fill the halo region for a_u and a_v.
   utmp=grid%xa%u
   vtmp=grid%xa%v
   grid%xa%u=grid%a_u_2
   grid%xa%v=grid%a_v_2

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_XA_A.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_XA_A_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE

   grid%a_u_2=grid%xa%u
   grid%a_v_2=grid%xa%v
   grid%xa%u=utmp
   grid%xa%v=vtmp
   utmp=0.0
   vtmp=0.0

   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            utmp(i,j,k)=utmp(i,j,k)+0.5*(grid%a_u_2(i+1,j  ,k)+grid%a_u_2(i,j,k))
            vtmp(i,j,k)=vtmp(i,j,k)+0.5*(grid%a_v_2(i  ,j+1,k)+grid%a_v_2(i,j,k))
         end do
      end do
   end do

   utmp(its-1,jts:jte,kts:kte)=utmp(its-1,jts:jte,kts:kte)+0.5*grid%a_u_2(its,jts:jte,kts:kte)
   utmp(ite+1,jts:jte,kts:kte)=utmp(ite+1,jts:jte,kts:kte)+0.5*grid%a_u_2(ite+1,jts:jte,kts:kte)
   vtmp(its:ite,jts-1,kts:kte)=vtmp(its:ite,jts-1,kts:kte)+0.5*grid%a_v_2(its:ite,jts,kts:kte)
   vtmp(its:ite,jte+1,kts:kte)=vtmp(its:ite,jte+1,kts:kte)+0.5*grid%a_v_2(its:ite,jte+1,kts:kte)

   ! The western boundary
   if (its == grid%xp%ids) then
      grid%xa%u(its  ,jts:jte,kts:kte)=grid%xa%u(its  ,jts:jte,kts:kte)+2.0*utmp(its-1,jts:jte,kts:kte)
      grid%xa%u(its+1,jts:jte,kts:kte)=grid%xa%u(its+1,jts:jte,kts:kte)-utmp(its-1,jts:jte,kts:kte)
   end if

   ! The eastern boundary
   if (ite == grid%xp%ide) then
      grid%xa%u(ite  ,jts:jte,kts:kte)=grid%xa%u(ite  ,jts:jte,kts:kte)+2.0*utmp(ite+1,jts:jte,kts:kte)
      grid%xa%u(ite-1,jts:jte,kts:kte)=grid%xa%u(ite-1,jts:jte,kts:kte)-utmp(ite+1,jts:jte,kts:kte)
   end if

   grid%xa%u=grid%xa%u+utmp

   ! The southern boundary
   if (jts == grid%xp%jds) then
      grid%xa%v(its:ite,jts  ,kts:kte)=grid%xa%v(its:ite,jts  ,kts:kte)+2.0*vtmp(its:ite,jts-1,kts:kte)
      grid%xa%v(its:ite,jts+1,kts:kte)=grid%xa%v(its:ite,jts+1,kts:kte)-vtmp(its:ite,jts-1,kts:kte)
   end if

   ! The northern boundary
   if (jte == grid%xp%jde) then
      grid%xa%v(its:ite,jte  ,kts:kte)=grid%xa%v(its:ite,jte  ,kts:kte)+2.0*vtmp(its:ite,jte+1,kts:kte)
      grid%xa%v(its:ite,jte-1,kts:kte)=grid%xa%v(its:ite,jte-1,kts:kte)-vtmp(its:ite,jte+1,kts:kte)
   end if

   grid%xa%v=grid%xa%v+vtmp

   grid%a_u_2 = 0.0
   grid%a_v_2 = 0.0

   !---------------------------------------------------------------------------
   ! [4.0] Adjoint of CONVERT TEMPERATURE inCREMENTS inTO THETA inCREMENTS
   !       EVALUATE ALSO THE inCREMENTS OF (1/rho) AND GEOPOTENTIAL
   !---------------------------------------------------------------------------

   a_press(its:ite,jts:jte,kts:kte+1)=0.0

   do k=kte,kts,-1
      do j=jts,jte
         do i=its,ite
            grid%xa%p(i,j,k)= grid%xa%p(i,j,k)+grid%a_p(i,j,k)
         end do
      end do
   end do

   grid%a_p = 0.0

!  do k=kte,kts,-1
!     do j=jts,jte
!        do i=its,ite
!           rho_cgrid=-(grid%ph_2(i,j,k+1)-grid%ph_2(i,j,k))*grid%a_ph_2(i,j,k+1)/grid%xb%rho(i,j,k)
!           a_press(i,j,k )=a_press(i,j,k )+grid%a_ph_2(i,j,k+1)/grid%xb%rho(i,j,k)
!           a_press(i,j,k+1)=a_press(i,j,k+1)-grid%a_ph_2(i,j,k+1)/grid%xb%rho(i,j,k)
!           grid%a_ph_2(i,j,k ) =grid%a_ph_2(i,j,k)   +grid%a_ph_2(i,j,k+1)
!           grid%xa%q(i,j,k)=grid%xa%q(i,j,k)-grid%xb%rho(i,j,k)*0.61*rho_cgrid/(1.0+0.61*grid%xb%q(i,j,k))
!           grid%xa%t(i,j,k)=grid%xa%t(i,j,k)-grid%xb%rho(i,j,k)*rho_cgrid/grid%xb%t(i,j,k)
!           grid%xa%p(i,j,k)= grid%xa%p(i,j,k)+grid%xb%rho(i,j,k)*rho_cgrid/grid%xb%p(i,j,k)
!        end do
!     end do
!  end do

   do k=kts,kte
      do j=jts,jte
         do i=its,ite 
            grid%xa%p(i,j,k)=grid%xa%p(i,j,k)-(t0+grid%t_2(i,j,k))*kappa*grid%a_t_2(i,j,k)/grid%xb%p(i,j,k)
            grid%xa%t(i,j,k)=grid%xa%t(i,j,k)+(t0+grid%t_2(i,j,k))*grid%a_t_2(i,j,k)/grid%xb%t(i,j,k)
         end do
      end do
   end do

   grid%a_t_2 = 0.0
   grid%a_ph_2 = 0.0
 
   !---------------------------------------------------------------------------
   ! [3.0] Adjoint of COMPUTE pressure increments (for computing theta increments)
   !---------------------------------------------------------------------------

   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            a_press(i,j,k+1)=a_press(i,j,k+1)+0.5*grid%xa%p(i,j,k)
            a_press(i,j,k )=a_press(i,j,k )+0.5*grid%xa%p(i,j,k)
            grid%xa%p(i,j,k)=0.0
            grid%a_moist(i,j,k,P_A_QV)=grid%a_moist(i,j,k,P_A_QV)-(grid%mu_2(i,j)+grid%mub(i,j))*a_press(i,j,k)*grid%dn(k)
            grid%a_mu_2(i,j)=grid%a_mu_2(i,j)-a_press(i,j,k)*(1.0+grid%moist(i,j,k,P_QV))*grid%dn(k)
            a_press(i,j,k+1)=a_press(i,j,k+1)+a_press(i,j,k)
         end do
      end do
   end do

   !---------------------------------------------------------------------------
   ! [2.0] Adjoint of COMPUTE increments of dry-column air mass per unit area
   !---------------------------------------------------------------------------

   do j=jts,jte
      do i=its,ite
         sdmd=0.0
         s1md=0.0
         do k=kts,kte
            s1md=s1md+(1.0+grid%moist(i,j,k,P_QV))*grid%dnw(k)
         end do
         sdmd=sdmd-grid%xb%psac(i,j)*grid%a_mu_2(i,j)/s1md
         grid%xa%psfc(i,j)=grid%xa%psfc(i,j)-grid%a_mu_2(i,j)/s1md
         do k=kts,kte
            grid%a_moist(i,j,k,P_A_QV)=grid%a_moist(i,j,k,P_A_QV)+sdmd*grid%dnw(k)
         end do
      end do
   end do

   grid%a_mu_2 = 0.0
   !---------------------------------------------------------------------------
   ! [1.0] Adjoint of Get the mixing ratio of moisture 
   !---------------------------------------------------------------------------
   do k=kts,kte
      do j=jts,jte
         do i=its,ite
            grid%xa%q(i,j,k)=grid%xa%q(i,j,k)+grid%a_moist(i,j,k,P_A_QV)/(1.0-grid%xb%q(i,j,k))**2
         end do
      end do
   end do

   if (size(grid%a_moist,dim=4) >= 4) then
      do k=kts,kte
         do j=jts,jte
            do i=its,ite
             grid%xa%qcw(i,j,k)=grid%xa%qcw(i,j,k) + grid%a_moist(i,j,k,P_A_QC)
             grid%xa%qrn(i,j,k)=grid%xa%qrn(i,j,k) + grid%a_moist(i,j,k,P_A_QR)
            end do
         end do
      end do
   end if

   grid%a_moist = 0.0

   grid%a_rainnc = 0.0
   grid%a_rainncv = 0.0
   grid%a_rainc = 0.0
   grid%a_raincv = 0.0

   if (trace_use) call da_trace_exit("da_transfer_xatowrftl_adj")

end subroutine da_transfer_xatowrftl_adj


subroutine da_transfer_xatowrftl_adj_lbc(grid, config_flags, filnam)

   !---------------------------------------------------------------------------
   ! Purpose: Convert WRFTL variables to analysis increments
   !           (inverse of the incremental part of xatowrf)
   !---------------------------------------------------------------------------

   implicit none
   

   type(domain), intent(inout)               :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   character*4, intent(in) :: filnam

end subroutine da_transfer_xatowrftl_adj_lbc


subroutine da_transfer_wrftl_lbc_t0_adj (grid)

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of WRF boundary for LBC at T=0
   ! All adjoint of tendencies variable are kept , because they will be used in T6 calculation. 
   !---------------------------------------------------------------------------

   implicit none
   
   type(domain), intent(inout)               :: grid

end subroutine da_transfer_wrftl_lbc_t0_adj


subroutine da_transfer_xatoanalysis(it, xbx, grid, config_flags)

   !---------------------------------------------------------------------------
   ! Purpose: Transfer xb and xa (increments) to analysis.
   !---------------------------------------------------------------------------

   implicit none

   integer,         intent(in)    :: it    ! outer-loop index
   type (xbx_type), intent(out)   :: xbx    ! Header & non-gridded vars.
   type(domain),    intent(inout) :: grid

   type (grid_config_rec_type), intent(inout) :: config_flags

   character*4 filnam
   character(len=256) ::  timestr

   if (trace_use) call da_trace_entry("da_transfer_xatoanalysis")

   !---------------------------------------------------------------------------
   ! Write out analysis in differing formats:
   !---------------------------------------------------------------------------      

   if (fg_format == fg_format_wrf_arw_regional) then
      if (write_increments .and. var4d) then
         write(unit=filnam,fmt='(a3,i1)') 'inc',it
         call domain_clock_get( grid, current_timestr=timestr )
         call da_transfer_xatowrftl(grid, config_flags, filnam, timestr)
      end if

      call da_transfer_xatowrf(grid, config_flags)

      if (it < max_ext_its) then
         call da_transfer_wrftoxb(xbx, grid, config_flags)
      end if
   else if (fg_format == fg_format_wrf_arw_global) then
      if( var4d) then
      write(unit=message(1),fmt='(A,I5)') &
         "var4d is not possible with Global WRF fg_format = ",fg_format
      call da_error("da_transfer_xatoanalysis.inc",40,message(1:1))
      else
       call da_transfer_xatowrf(grid, config_flags)
       if (it < max_ext_its)call da_transfer_wrftoxb(xbx, grid, config_flags)
      end if
   else if (fg_format == fg_format_wrf_nmm_regional) then
      call da_transfer_xatowrf_nmm_regional(grid)
      if (it < max_ext_its) then
         if (var4d) then
      write(unit=message(1),fmt='(A,I5)') &
         "var4d is not possible for fg_format = ",fg_format
      call da_error("da_transfer_xatoanalysis.inc",51,message(1:1))
         end if

         call da_transfer_wrf_nmm_regional_toxb(xbx, grid)
      end if
   else if (fg_format == fg_format_kma_global) then
      call da_transfer_xatokma(grid)
      if (it < max_ext_its) then
         call da_transfer_kmatoxb(xbx, grid)
      end if
   end if

   if (trace_use) call da_trace_exit("da_transfer_xatoanalysis")

end subroutine da_transfer_xatoanalysis


subroutine da_setup_firstguess(xbx, grid, config_flags, ens )

   !---------------------------------------------------------------------------
   ! Purpose: Allocate and read in components of first guess state.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------------

   implicit none

   type (xbx_type),intent(out)  :: xbx   ! Header & non-gridded vars.

   type(domain),intent(inout)   :: grid
   type(grid_config_rec_type), intent(in) :: config_flags
   logical, intent(in) :: ens

   integer :: is, ie, js, je, ij, i, j
   real    :: ddx , ddy    

   if (trace_use) call da_trace_entry("da_setup_firstguess")

   is = grid%xp % its
   ie = grid%xp % ite
   js = grid%xp % jts
   je = grid%xp % jte

   ! Calculate sin and cosine values used in da_get_avpoles

   if (global) then
      if (grid%xp%jts == grid%xp%jds) then

         allocate(cos_xls(grid%xp%its:grid%xp%ite))
         allocate(sin_xls(grid%xp%its:grid%xp%ite))
         cos_xls(grid%xp%its:grid%xp%ite) = & 
            cos(deg_to_rad*grid%xlong(grid%xp%its:grid%xp%ite,grid%xp%jts))
         sin_xls(grid%xp%its:grid%xp%ite) = &
            sin(deg_to_rad*grid%xlong(grid%xp%its:grid%xp%ite,grid%xp%jts))
      end if

      if (grid%xp%jte == grid%xp%jde) then 
         allocate(cos_xle(grid%xp%its:grid%xp%ite))
         allocate(sin_xle(grid%xp%its:grid%xp%ite))
         cos_xle(grid%xp%its:grid%xp%ite) = &
            cos(deg_to_rad*grid%xlong(grid%xp%its:grid%xp%ite,grid%xp%jte))
         sin_xle(grid%xp%its:grid%xp%ite) = &
            sin(deg_to_rad*grid%xlong(grid%xp%its:grid%xp%ite,grid%xp%jte))
      end if
   end if

   !---------------------------------------------------------------------------      
   ! [1.0] Setup and read in fields from first guess:
   !---------------------------------------------------------------------------      

   if ((fg_format==fg_format_wrf_arw_regional) .or. &
      (fg_format==fg_format_wrf_arw_global  ) ) then
      call da_setup_firstguess_wrf(xbx, grid, config_flags,ens)
      ! when ens=.true., da_setup_firstguess(_wrf) is called solely for map_info,
      ! the rest of the code should be skipped
      if ( ens ) then
         if (trace_use) call da_trace_exit("da_setup_firstguess")
         return
      end if
   else if (fg_format == fg_format_wrf_nmm_regional ) then
      call da_setup_firstguess_wrf_nmm_regional(xbx, grid)
   else if (fg_format == fg_format_kma_global) then
      ! First guess is an KMA format file:
      call da_setup_firstguess_kma(xbx, grid)
   end if

   !---------------------------------------------------------------------------
   ! Exchange halo region for XB arrays.
   !---------------------------------------------------------------------------

   if ((fg_format==fg_format_wrf_arw_regional) .or. &
      (fg_format==fg_format_wrf_arw_global  ) ) then
      ! Calculate multiplicative constants for PsiChi_TO_UV 
      !$OMP PARALLEL DO &
      !$OMP PRIVATE (ij, i, j)
      do ij = 1, grid%num_tiles
         do j = grid%j_start(ij), grid%j_end(ij)
            do i = is, ie
               grid%xb%coefx(i,j) = 0.5 * grid%xb%map_factor(i,j)/grid%xb%ds
               grid%xb%coefy(i,j) = grid%xb%coefx(i,j)
               grid%xb%coefz(i,j) = 0.5 / (grid%xb%map_factor(i,j)*grid%xb%ds)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
   else if (fg_format == fg_format_wrf_nmm_regional) then
      grid%xb%coefx(is:ie,js:je) = 0.5/grid%mu0(is:ie,js:je)
      grid%xb%coefy(is:ie,js:je) = 0.5/grid%xb%ds         
   else if (fg_format == fg_format_kma_global) then
      ! Calculate multiplicative constants for PsiChi_TO_UV 
      ddx =  earth_radius*1000 * 2.0 * pi / (grid%xb%ide-grid%xb%ids+1)
      ddy =  earth_radius*1000       * pi / (grid%xb%jde-grid%xb%jds)
      grid%xb% coefx(is:ie,js:je) = 0.5 / (ddx * cos(grid%xlat(is:ie,js:je)*pi/180.))
      grid%xb% coefy(is:ie,js:je) = 0.5 /  ddy
   else
      write(unit=message(1),fmt='(A,I5)') &
         "Wrong choice for fg_format = ",fg_format
      call da_error("da_setup_firstguess.inc",101,message(1:1))
   end if

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_INIT.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL HALO_INIT_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
!ENDOFREGISTRYGENERATEDINCLUDE


   periodic_x = grid%periodic_x

   if (global) then     
      ! Set East-West boundary for Xb-array 
      call da_set_boundary_xb(grid)
   end if

   !---------------------------------------------------------------------------      
   ! [2.0] Setup grid-dependent constants used:
   !---------------------------------------------------------------------------

   ! [2.1] Set up fast Fourier & Legendre transform constants:

   if (SIZE(xbx%fft_factors_x) /=num_fft_factors) &
      call da_setup_runconstants(grid, xbx)

   if (trace_use) call da_trace_exit("da_setup_firstguess")

end subroutine da_setup_firstguess


subroutine da_setup_firstguess_wrf(xbx, grid, config_flags, ens)

   !---------------------------------------------------------------------------
   ! Purpose: Define/allocate components of WRF model state.
   !---------------------------------------------------------------------------

   implicit none

   type (xbx_type), intent(out)         :: xbx    ! Header & non-gridded vars.

   type (domain), intent(inout)         :: grid
   type(grid_config_rec_type), intent(in) :: config_flags
   logical, intent(in) :: ens

   integer           :: map_util_project
   real              :: x, y, lat_cen, lon_cen
  
   real              :: buf(2)

   character(len=24) :: xb_date, an_date
   integer           :: len, seconds, i_grid,  j_grid, m_expand


   if (trace_use) call da_trace_entry("da_setup_firstguess_wrf")

   !-----------------------------------------------------------------------
   ! [0.0] check the xb_date for 3DVAR
   !-----------------------------------------------------------------------

   if ( num_fgat_time == 1 ) then
      write(unit=xb_date,fmt='(i4.4,2("-",i2.2),"_",i2.2,2(":",i2.2),".0000")')  &
           grid%start_year, grid%start_month, grid%start_day, &
           grid%start_hour, grid%start_minute,grid%start_second

      len = len_trim(ANALYSIS_DATE)

      write(unit=an_date(1:len), fmt='(a)') trim(ANALYSIS_DATE)

      seconds = int(da_diff_seconds(an_date, xb_date))

      if (seconds > analysis_accu) then
         message(1)="The time of your first guess file is different from the analysis date"
         write(unit=message(2),fmt='(A,A)') &
             "xb_date    = ",xb_date
         write(unit=message(3),fmt='(A,A)') &
             "an_date    = ",an_date
         write(unit=message(4),fmt='(A,I10,A)') &
            "Difference between analysis_date and first guess date is ",seconds," seconds."
         message(5)="Correct your choice of fg file, or change analysis_date or analysis_accu in namelist.input"
         call da_error("da_setup_firstguess_wrf.inc",50,message(1:5))
      end if
   end if

   !------------------------------------------------------------------------
   ! [1.0] Read original WRF format first guess:
   !------------------------------------------------------------------------
   
   !------------------------------------------------------------------------
   ! [2.0] Copy header info:
   !------------------------------------------------------------------------

   if ((grid%xp%its == grid%xp%ids) .and. (grid%xp%jts == grid%xp%jds)) then
      buf(1) = grid%xlat(grid%xp%its, grid%xp%jts)
      buf(2) = grid%xlong(grid%xp%its, grid%xp%jts)
   end if
   
   call wrf_dm_bcast_real(buf, 2)
   start_lat=buf(1)
   start_lon=buf(2)

   !------------------------------------------------------------------------
   ! Setup map utility
   !------------------------------------------------------------------------

   call nl_get_map_proj     (grid%id , grid%map_proj)
   call nl_get_truelat1     (grid%id , grid%truelat1)
   call nl_get_truelat2     (grid%id , grid%truelat2)
   call nl_get_dx           (grid%id , grid%dx)
   call nl_get_cen_lat      (grid%id , grid%cen_lat)
   call nl_get_cen_lon      (grid%id , grid%cen_lon)
   call nl_get_pole_lat     (grid%id , grid%pole_lat)
   call nl_get_moad_cen_lat (grid%id , grid%moad_cen_lat)
   call nl_get_stand_lon    (grid%id , grid%stand_lon)

   phic = grid%moad_cen_lat
   xlonc = grid%stand_lon

   truelat1_3dv = grid%truelat1
   truelat2_3dv = grid%truelat2
   pole = 90.0
   dsm = 0.001 * grid%dx

   map_util_project = grid%map_proj

   ! Print mapping info to log file(s)
   write(unit=message(1), fmt='(a)') 'Domain mapping info:'
   write(unit=message(2), fmt='(a, i6)') &
        'map_proj =', grid%map_proj
   write(unit=message(3), fmt='(a, e16.6)') &
        'cen_lat   =', grid%cen_lat
   write(unit=message(4), fmt='(a, e16.6)') &
        'cen_lon   =', grid%cen_lon
   write(unit=message(5), fmt='(a, e16.6)') &
        'truelat1  =', grid%truelat1
   write(unit=message(6), fmt='(a, e16.6)') &
        'truelat2  =', grid%truelat2
   write(unit=message(7), fmt='(a, e16.6)') &
        'start_lat =', start_lat
   write(unit=message(8), fmt='(a, e16.6)') &
        'start_lon =', start_lon
   write(unit=message(9), fmt='(a, e16.6)') &
        'pole_lat  =', grid%pole_lat
   write(unit=message(10), fmt='(a, e16.6)') &
        'dsm       =', dsm
   call da_message(message(1:10))

   ! Set map projection in WRFSI world.
   map_util_project = PROJ_LC

   if (grid%map_proj == 0 .or. grid%map_proj == 6 ) then
      map_util_project = PROJ_LATLON

      if (grid%pole_lat < 89.9) then
         write(unit=message(1),fmt='(A,E10.6)')'POLE_LAT = ',grid%pole_lat
         write(unit=message(2),fmt='(A)')"WRFDA does not support rotated cylindrical equidistant projection"
         write(unit=message(3),fmt='(A)')"Choose a first guess file with a valid projection"
         call da_error("da_setup_firstguess_wrf.inc",127,message(1:3))
      end if
   else if (grid%map_proj == 1) then
      map_util_project = PROJ_LC
   else if (grid%map_proj == 2) then
      map_util_project = PROJ_PS
   else if (grid%map_proj == 3) then
      map_util_project = PROJ_MERC
   else
      write(unit=message(1),fmt='(A,I6)')'map_proj = ',grid%map_proj
      write(unit=message(2),fmt='(A)')"WRFDA does not support the selected map projection"
      write(unit=message(3),fmt='(A)')"Choose a first guess file with a valid projection"
      call da_error("da_setup_firstguess_wrf.inc",139,message(1:3))
   end if

   call da_map_set(map_util_project,grid%cen_lat,grid%cen_lon,   &
                real(grid%xp%ide-grid%xp%ids+2)/2.0, real(grid%xp%jde-grid%xp%jds+2)/2.0, &
                grid%dx,grid%stand_lon,grid%truelat1,grid%truelat2,grid%truelat1,grid%stand_lon,map_info)

   ! Need to set map projection in WRF world.
   map_projection = grid%map_proj

   cone_factor = map_info%cone

   if (.not. global .and. print_detail_map) then
     
      !----------------------------------------------------------------------
      ! Check the ll_to_ij:
      !----------------------------------------------------------------------

      message(1)="Check the map_set correctness::::::::::::::::::::::::"

      ! Domain center:
      call  da_llxy_wrf(map_info, grid%cen_lat, grid%cen_lon, start_x, start_y)
      write(unit=message(2),fmt='("Center: latc,lonc,x,y, Xc, Yc:",6f10.3)') &
                  grid%cen_lat, grid%cen_lon, start_x, start_y, &
                  real(grid%xp%ide-grid%xp%ids+2)/2.0, real(grid%xp%jde-grid%xp%jds+2)/2.0

      start_x = real(grid%xp%ide-grid%xp%ids+2)/2.0
      start_y = real(grid%xp%jde-grid%xp%jds+2)/2.0
      lat_cen = -999.9
      lon_cen = -999.9
      call  da_xyll(map_info, start_x, start_y, lat_cen, lon_cen)
      write(unit=message(3), &
         fmt='("Center: X, Y, latc, lonc, phic, xlonc:",6f10.3)') &
         start_x, start_y, lat_cen, lon_cen,   &
         grid%cen_lat, grid%cen_lon
      call da_message(message(1:3))
   end if

   ! Setup the domain definition for use of the GRAPH:

   coarse_ds = 0.001 * grid%dx
   coarse_ix = grid%e_we - grid%s_we + 1
   coarse_jy = grid%e_sn - grid%s_sn + 1
   start_x = 1.0
   start_y = 1.0

   if( fg_format==fg_format_kma_global) then
   delt_lat = 180.0/real(grid%e_sn - grid%s_sn - 1)
   delt_lon = 360.0/real(grid%e_we - grid%s_we)
   else if( fg_format==fg_format_wrf_arw_global) then
   delt_lat = 180.0/real(grid%e_sn - grid%s_sn)
   delt_lon = 360.0/real(grid%e_we - grid%s_we)
   end if

   !--------------------------------------------------------------------------
   ! [3.0] Interpolate WRF C-grid winds to p points of WRFVAR grid (interpolate 
   ! u to west, v to south?
   !---------------------------------------------------------------------------

   grid%xb % mix = grid%xp%ide - grid%xp%ids + 1
   grid%xb % mjy = grid%xp%jde - grid% xp%jds + 1
   grid%xb % mkz = grid%xp%kde - grid%xp%kds + 1

   grid%xb % ds  = 0.001 * grid%dx

   mix = grid%xb % mix
   mjy = grid%xb % mjy
   mkz = grid%xb % mkz
   
   if ( .not. ens ) then
      call da_transfer_wrftoxb(xbx, grid, config_flags)
   endif

   if (trace_use) call da_trace_exit("da_setup_firstguess_wrf")

end subroutine da_setup_firstguess_wrf


subroutine da_setup_firstguess_wrf_nmm_regional(xbx, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Define/allocate components of WRF model state.
   !---------------------------------------------------------------------------

   implicit none

   type (xbx_type), intent(out)         :: xbx    ! Header & non-gridded vars.

   type (domain), intent(inout)         :: grid

   integer           :: map_util_project
   real              :: x, y, xxc, yyc, lat_cen, lon_cen
  
   real              :: buf(2)

   character(len=24) :: xb_date, an_date
   integer           :: len, seconds, i_grid,  j_grid, m_expand
   real              :: latinc, loninc

   if (trace_use) call da_trace_entry("da_setup_firstguess_wrf_nmm_regional")

   !-----------------------------------------------------------------------
   ! [0.0] check the xb_date for 3DVAR
   !-----------------------------------------------------------------------

   write(unit=xb_date,fmt='(i4.4,2("-",i2.2),"_",i2.2,2(":",i2.2),".0000")')  &
        grid%start_year, grid%start_month, grid%start_day, &
        grid%start_hour, grid%start_minute,grid%start_second

   len = len_trim(ANALYSIS_DATE)

   write(unit=an_date(1:len), fmt='(a)') trim(ANALYSIS_DATE)

   seconds = int(da_diff_seconds(an_date, xb_date))

   if (seconds > ANALYSIS_ACCU) then
      write(unit=message(1),fmt='(A,A,A,A)') &
         "xb_date=",xb_date," an_date=", an_date
      write(unit=message(2),fmt='(A,I6,A,I6)') &
         "diff=",seconds,"   ANALYSIS_ACCU=",ANALYSIS_ACCU
      message(3)="=======> Wrong xb time found???"
      call da_warning("da_setup_firstguess_wrf_nmm_regional.inc",44,message(1:3))
   end if

   !------------------------------------------------------------------------
   ! [1.0] Read original WRF format first guess:
   !------------------------------------------------------------------------
   
   !------------------------------------------------------------------------
   ! [2.0] Copy header info:
   !------------------------------------------------------------------------

   if ((grid%xp%its == grid%xp%ids) .and. (grid%xp%jts == grid%xp%jds)) then
      buf(1) = grid%xlat(grid%xp%its, grid%xp%jts)
      buf(2) = grid%xlong(grid%xp%its, grid%xp%jts)
   end if
   
   call wrf_dm_bcast_real(buf, 2)
   start_lat=buf(1)
   start_lon=buf(2)

   !------------------------------------------------------------------------
   ! Setup map utility
   !------------------------------------------------------------------------

   call nl_get_map_proj     (grid%id , grid%map_proj)
   call nl_get_truelat1     (grid%id , grid%truelat1)
   call nl_get_truelat2     (grid%id , grid%truelat2)
   call nl_get_dx           (grid%id , grid%dx)
   call nl_get_cen_lat      (grid%id , grid%cen_lat)
   call nl_get_cen_lon      (grid%id , grid%cen_lon)
   call nl_get_moad_cen_lat (grid%id , grid%moad_cen_lat)
   call nl_get_stand_lon    (grid%id , grid%stand_lon)

   phic   = grid%moad_cen_lat
   xlonc  = grid%cen_lon   
   loninc = grid%cf1       
   latinc = grid%cf2       

   truelat1_3dv = grid%truelat1
   truelat2_3dv = grid%truelat2
   pole = 90.0
   dsm = 0.001 * grid%dx

   map_util_project = grid%map_proj
   ! Set map projection in WRFSI world.

   if (grid%map_proj == 0 .or. grid%map_proj == 6) then
      map_util_project = PROJ_LATLON
   else if (grid%map_proj == 1) then
      map_util_project = PROJ_LC
   else if (grid%map_proj == 2) then
      map_util_project = PROJ_PS
   else if (grid%map_proj == 3) then
      map_util_project = PROJ_MERC
   end if

   call da_map_set(map_util_project,grid%cen_lat,grid%cen_lon,   &
                real(grid%xp%ide-grid%xp%ids+2)/2.0, real(grid%xp%jde-grid%xp%jds+2)/2.0, &
                grid%dx,grid%stand_lon,grid%truelat1,grid%truelat2,latinc,loninc,map_info)
   ! Need to set map projection in WRF world.
   map_projection = grid%map_proj
   cone_factor = map_info%cone

   if (print_detail_map) then
      write(unit=stdout, fmt='(a, i6)') &
           'map_proj =', grid%map_proj

      write(unit=stdout, fmt='(a, e16.6)') &
           'cen_lat  =', grid%cen_lat,  &
           'cen_lon  =', grid%cen_lon,  &
           'truelat1 =', grid%truelat1, &
           'truelat2 =', grid%truelat2, &
           'start_lat =', start_lat, &
           'start_lon =', start_lon, &
           'latinc    =', latinc   , &
           'loninc    =', loninc   , &
           'cone_fact =', cone_factor, &
           'dsm      =', dsm
   end if


    mix = grid%xp%ide - grid%xp%ids  + 1
    mjy = grid%xp%jde - grid% xp%jds + 1
    mkz = grid%xp%kde - grid%xp%kds  + 1

   call da_transfer_wrf_nmm_regional_toxb(xbx, grid)

   if (trace_use) call da_trace_exit("da_setup_firstguess_wrf_nmm_regional")

end subroutine da_setup_firstguess_wrf_nmm_regional


subroutine da_setup_firstguess_kma(xbx, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Define/allocate components of WRF model state.
   !---------------------------------------------------------------------------

   implicit none

   type (xbx_type),intent(out)         :: xbx    ! Header & non-gridded vars.

   type (domain), intent(inout)        :: grid

   integer           :: i, j
   integer           :: is, ie, js, je
   integer           :: max_wavenumber

   if (trace_use) call da_trace_entry("da_setup_firstguess_kma")

   is = grid%xp % its
   ie = grid%xp % ite
   js = grid%xp % jts
   je = grid%xp % jte

   !------------------------------------------------------------------------
   ! [2.0] Copy header info:
   !------------------------------------------------------------------------

   ! rizvi set it to 1 . Actually it should be decided by KMA 
   grid%map_proj = 0
   map_projection = grid%map_proj
   coarse_ix = grid%e_we - grid%s_we + 1
   coarse_jy = grid%e_sn - grid%s_sn + 1

   grid%xb % mix = grid%xp%ide - grid%xp%ids + 1
   grid%xb % mjy = grid%xp%jde - grid%xp%jds + 1
   grid%xb % mkz = grid%xp%kde - grid%xp%kds + 1

   mix = grid%xb % mix
   mjy = grid%xb % mjy
   mkz = grid%xb % mkz

   grid%xb % ds  = 0.001 * grid%dx

   start_x = 1.0
   start_y = 1.0
   start_lat = -90.0
   start_lon = -180.0
   delt_lat = 180.0/real(grid%e_sn - grid%s_sn - 1)
   delt_lon = 360.0/real(grid%e_we - grid%s_we)

   phic        = 0.0
   xlonc       = 0.0
   cone_factor = 0.0


   do j = js,je
      do i = is,ie
         grid%xlat(i,j) = start_lat + real(j-1)*delt_lat
         grid%xlong(i,j) = start_lon + real(i-1)*delt_lon
      end do
   end do

   ! Avoid assigning -90,90 value                                           
   if (grid%xb%jts == grid%xb%jds) then
      grid%xlat(is:ie,j) = -89.9                               
   end if

   if (grid%xb%jte == grid%xb%jde) then
      grid%xlat(is:ie,j) = 89.9                               
   end if

   ! fix map factor and coriolis parameter
      
   grid%f(is:ie,js:je) = 2.0 *earth_omega*sin(pi*grid%xlat(is:ie,js:je)/180.0)

   xbx%inc = 1
   xbx%ni  = grid%e_we - grid%s_we
   xbx%nj  = grid%e_sn - grid%s_sn
   xbx%nk  = grid%e_vert - 1
   xbx% lenwrk    = xbx%ni
   xbx% lenr           = xbx%inc * (xbx%ni - 1) + 1
   max_wavenumber = xbx%ni/2-1
   xbx % alp_size = (xbx%nj+1)*(max_wavenumber+1)*(max_wavenumber+2)/4

   call da_transfer_kmatoxb(xbx, grid)

   if (trace_use) call da_trace_exit("da_setup_firstguess_kma")

end subroutine da_setup_firstguess_kma


 SUBROUTINE da_get_2nd_firstguess ( grid )

!-------------------------------------------------------------------------
!  Calculate the first guess at the end of thr time window
!  The original grid%u_2,v_2 etc. will be overwrittedm but we need the boundary area only
!  Author: Xin Zhang, 10/7/2010
!-------------------------------------------------------------------------

   IMPLICIT NONE

   TYPE(domain), INTENT(INOUT)        :: grid


 END SUBROUTINE da_get_2nd_firstguess

end module da_transfer_model
