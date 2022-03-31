!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_interpdown_unpack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
IF ( SIZE(grid%lu_index) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%lu_index(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%u_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%u_2(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%v_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%v_2(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%w_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,ckde
grid%w_2(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%ph_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,ckde
grid%ph_2(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%phb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,ckde
grid%phb(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%t_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%t_2(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%mu_2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%mu_2(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%mub) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%mub(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%alb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%alb(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%pb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%pb(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%q2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%q2(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%t2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%t2(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%th2) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%th2(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%psfc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%psfc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%u10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%u10(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%v10) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%v10(pig,pjg) = xv(1)
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_moist
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
moist(pig,pjg,k,itrace) = xv(k)
ENDDO
ENDDO
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
scalar(pig,pjg,k,itrace) = xv(k)
ENDDO
ENDDO
IF ( SIZE(grid%smois) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
DO k = 1,config_flags%num_soil_layers
grid%smois(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%sh2o) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
DO k = 1,config_flags%num_soil_layers
grid%sh2o(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%smcrel) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
DO k = 1,config_flags%num_soil_layers
grid%smcrel(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%xice) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xice(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%ivgtyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%ivgtyp(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%isltyp) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%isltyp(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%vegfra) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%vegfra(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%snow) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%snow(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%snowh) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%snowh(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%tslb) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(((config_flags%num_soil_layers)-(1)+1)*RWORDSIZE,xv)
DO k = 1,config_flags%num_soil_layers
grid%tslb(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%msft) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msft(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfu) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfu(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfv(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msftx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msftx(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfty) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfty(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfux) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfux(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfuy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfuy(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfvx) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfvx(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfvx_inv) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfvx_inv(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%msfvy) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%msfvy(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%f) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%f(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%e) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%e(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%sina) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%sina(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%cosa) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%cosa(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%ht) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%ht(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%tsk) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%tsk(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%rainc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%rainnc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xlat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xlat(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xlong) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xlong(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xlat_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xlat_u(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xlong_u) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xlong_u(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xlat_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xlat_v(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xlong_v) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xlong_v(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%clat) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%clat(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%tmn) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%tmn(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%xland) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%xland(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%snowc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%snowc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%g_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%g_h_diabatic(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%a_h_diabatic) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
grid%a_h_diabatic(pig,pjg,k) = xv(k)
ENDDO
ENDIF
IF ( SIZE(grid%g_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%g_rainc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%a_rainc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%a_rainc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%g_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%g_rainnc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%a_rainnc) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%a_rainnc(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%landmask) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%landmask(pig,pjg) = xv(1)
ENDIF
IF ( SIZE(grid%sst) .GT. 1 ) THEN ! okay for intermediate_grid too. see comment in gen_comms.c
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%sst(pig,pjg) = xv(1)
ENDIF
DO itrace =  PARAM_FIRST_SCALAR, num_tracer
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
tracer(pig,k,pjg,itrace) = xv(k)
ENDDO
ENDDO
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
