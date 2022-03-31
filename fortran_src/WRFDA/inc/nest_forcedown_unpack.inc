!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_forcedown_unpack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
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
DO itrace =  PARAM_FIRST_SCALAR, num_tracer
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
tracer(pig,k,pjg,itrace) = xv(k)
ENDDO
ENDDO
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
