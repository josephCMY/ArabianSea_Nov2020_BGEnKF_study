!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_feedbackup_pack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
msize = (14 + ((num_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_scalar - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_a_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_g_moist - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_a_scalar - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_g_scalar - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_chem - PARAM_FIRST_SCALAR + 1)) & 
 + ((num_tracer - PARAM_FIRST_SCALAR + 1)) )* nlev + 48
CALL rsl_lite_to_parent_info( msize*RWORDSIZE                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,nids,nide,njds,njde                               &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(parent_grid%id)                            &
                        ,nest_pes_y(parent_grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
DO while ( retval .eq. 1 )
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%lu_index(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%u_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%u_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%v_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%v_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%w_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k)= intermediate_grid%w_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%ph_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k)= intermediate_grid%ph_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%phb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k)= intermediate_grid%phb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%t_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%t_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%mu_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%mu_2(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%mub) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%mub(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%alb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%alb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%pb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%pb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%q2(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%t2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%t2(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%th2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%th2(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%psfc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%psfc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%u10(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%v10(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_moist
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%moist(pig,pjg,k,itrace)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDDO
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%scalar(pig,pjg,k,itrace)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDDO
IF ( SIZE(grid%smois) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= intermediate_grid%smois(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= intermediate_grid%sh2o(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%smcrel) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= intermediate_grid%smcrel(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xice(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%ivgtyp(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%isltyp(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%vegfra(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%snow) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%snow(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%snowh(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%tslb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= intermediate_grid%tslb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msft) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msft(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfu) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfu(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfv(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msftx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msftx(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfty) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfty(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfux) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfux(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfuy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfuy(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfvx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfvx(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfvx_inv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfvx_inv(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfvy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%msfvy(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%f) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%f(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%e) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%e(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%sina) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%sina(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%cosa) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%cosa(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%ht) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%ht(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%tsk) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%tsk(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%rainc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%rainnc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xlat(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlong) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xlong(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlat_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xlat_u(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlong_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xlong_u(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlat_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xlat_v(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlong_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xlong_v(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%clat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%clat(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%tmn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%tmn(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xland) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%xland(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%snowc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%snowc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%g_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%g_h_diabatic(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%a_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%a_h_diabatic(pig,pjg,k)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%g_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%g_rainc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%a_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%a_rainc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%g_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%g_rainnc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%a_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%a_rainnc(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)= intermediate_grid%landmask(pig,pjg)
CALL rsl_lite_to_parent_msg(RWORDSIZE,xv)
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_tracer
DO k = ckds,(ckde-1)
xv(k)= intermediate_grid%tracer(pig,k,pjg,itrace)
ENDDO
CALL rsl_lite_to_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDDO
CALL rsl_lite_to_parent_info( msize*RWORDSIZE                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,nids,nide,njds,njde                               &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(parent_grid%id)                            &
                        ,nest_pes_y(parent_grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
