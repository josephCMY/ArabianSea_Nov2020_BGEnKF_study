!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_feedbackup_unpack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL rsl_lite_from_child_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%lu_index(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%u_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .TRUE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%u_2(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%v_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%v_2(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%w_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,ckde
NEST_INFLUENCE(grid%w_2(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%ph_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,ckde
NEST_INFLUENCE(grid%ph_2(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%phb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,ckde
NEST_INFLUENCE(grid%phb(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%t_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%t_2(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%mu_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%mu_2(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%mub) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%mub(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%alb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%alb(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%pb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%pb(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%q2(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%t2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%t2(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%th2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%th2(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%psfc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%psfc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%u10(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%v10(pig,pjg),xv(1))
ENDIF
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_moist
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(moist(pig,pjg,k,itrace),xv(k))
ENDDO
ENDIF
ENDDO
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(scalar(pig,pjg,k,itrace),xv(k))
ENDDO
ENDIF
ENDDO
IF ( SIZE(grid%smois) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = 1,config_flags%num_soil_layers
NEST_INFLUENCE(grid%smois(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = 1,config_flags%num_soil_layers
NEST_INFLUENCE(grid%sh2o(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%smcrel) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = 1,config_flags%num_soil_layers
NEST_INFLUENCE(grid%smcrel(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%xice(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%ivgtyp(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%isltyp(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%vegfra(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%snow) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%snow(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%snowh(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%tslb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = 1,config_flags%num_soil_layers
NEST_INFLUENCE(grid%tslb(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%msft) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%msft(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfu) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .TRUE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%msfu(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
NEST_INFLUENCE(grid%msfv(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msftx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%msftx(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfty) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%msfty(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfux) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .TRUE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%msfux(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfuy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .TRUE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%msfuy(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfvx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
NEST_INFLUENCE(grid%msfvx(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfvx_inv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
NEST_INFLUENCE(grid%msfvx_inv(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%msfvy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
NEST_INFLUENCE(grid%msfvy(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%f) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%f(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%e) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%e(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%sina) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%sina(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%cosa) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%cosa(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%ht) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%ht(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%tsk) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%tsk(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%rainc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%rainnc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xlat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%xlat(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xlong) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%xlong(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xlat_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .TRUE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%xlat_u(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xlong_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .TRUE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%xlong_u(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xlat_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
NEST_INFLUENCE(grid%xlat_v(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xlong_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .TRUE. ) ) THEN
NEST_INFLUENCE(grid%xlong_v(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%clat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%clat(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%tmn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%tmn(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%xland) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%xland(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%snowc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%snowc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%g_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%g_h_diabatic(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%a_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(grid%a_h_diabatic(pig,pjg,k),xv(k))
ENDDO
ENDIF
ENDIF
IF ( SIZE(grid%g_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%g_rainc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%a_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%a_rainc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%g_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%g_rainnc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%a_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%a_rainnc(pig,pjg),xv(1))
ENDIF
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_child_msg(RWORDSIZE,xv)
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
NEST_INFLUENCE(grid%landmask(pig,pjg),xv(1))
ENDIF
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_tracer
CALL rsl_lite_from_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv) ;
IF ( cd_feedback_mask( pig, ips_save, ipe_save , pjg, jps_save, jpe_save, .FALSE., .FALSE. ) ) THEN
DO k = ckds,(ckde-1)
NEST_INFLUENCE(tracer(pig,k,pjg,itrace),xv(k))
ENDDO
ENDIF
ENDDO
CALL rsl_lite_from_child_info(pig,pjg,retval)
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
