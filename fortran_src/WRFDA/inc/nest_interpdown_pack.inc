!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_interpdown_pack.inc'
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
 + ((num_tracer - PARAM_FIRST_SCALAR + 1)) )* nlev + 49
CALL rsl_lite_to_child_info( msize*RWORDSIZE                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,iids,iide,ijds,ijde                               &
                        ,nids,nide,njds,njde                               &
                        ,pgr , sw                                          &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(grid%id)                            &
                        ,nest_pes_y(grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
DO while ( retval .eq. 1 )
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%lu_index(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%u_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%u_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%v_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%v_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%w_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k)= grid%w_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%ph_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k)= grid%ph_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%phb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,ckde
xv(k)= grid%phb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%t_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%t_2(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%mu_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%mu_2(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%mub) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%mub(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%alb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%alb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%pb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%pb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%q2(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%t2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%t2(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%th2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%th2(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%psfc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%psfc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%u10(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%v10(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_moist
DO k = ckds,(ckde-1)
xv(k)= moist(pig,pjg,k,itrace)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDDO
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
DO k = ckds,(ckde-1)
xv(k)= scalar(pig,pjg,k,itrace)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDDO
IF ( SIZE(grid%smois) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= grid%smois(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= grid%sh2o(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%smcrel) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= grid%smcrel(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xice(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%ivgtyp(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%isltyp(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%vegfra(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%snow) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%snow(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%snowh(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%tslb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = 1,config_flags%num_soil_layers
xv(k)= grid%tslb(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msft) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msft(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfu) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfu(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfv(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msftx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msftx(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfty) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfty(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfux) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfux(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfuy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfuy(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfvx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfvx(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfvx_inv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfvx_inv(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%msfvy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%msfvy(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%f) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%f(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%e) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%e(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%sina) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%sina(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%cosa) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%cosa(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%ht) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%ht(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%tsk) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%tsk(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%rainc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%rainnc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xlat(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlong) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xlong(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlat_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xlat_u(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlong_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xlong_u(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlat_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xlat_v(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xlong_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xlong_v(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%clat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%clat(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%tmn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%tmn(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%xland) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%xland(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%snowc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%snowc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%g_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%g_h_diabatic(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%a_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
DO k = ckds,(ckde-1)
xv(k)= grid%a_h_diabatic(pig,pjg,k)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%g_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%g_rainc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%a_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%a_rainc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%g_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%g_rainnc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%a_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%a_rainnc(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%landmask(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
xv(1)=grid%sst(pig,pjg)
CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_tracer
DO k = ckds,(ckde-1)
xv(k)= tracer(pig,k,pjg,itrace)
ENDDO
CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
ENDDO
CALL rsl_lite_to_child_info( msize*RWORDSIZE                               &
                        ,cips,cipe,cjps,cjpe                               &
                        ,iids,iide,ijds,ijde                               &
                        ,nids,nide,njds,njde                               &
                        ,pgr , sw                                          &
                        ,nest_task_offsets(ngrid%id)                      &
                        ,nest_pes_x(grid%id)                            &
                        ,nest_pes_y(grid%id)                            &
                        ,nest_pes_x(intermediate_grid%id)                 &
                        ,nest_pes_y(intermediate_grid%id)                 &
                        ,thisdomain_max_halo_width                         &
                        ,icoord,jcoord                                     &
                        ,idim_cd,jdim_cd                                   &
                        ,pig,pjg,retval )
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
