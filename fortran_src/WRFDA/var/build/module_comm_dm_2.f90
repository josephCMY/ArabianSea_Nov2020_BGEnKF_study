


























MODULE module_comm_dm_2

   IMPLICIT NONE

   PRIVATE module_comm_dm_dummy_2

   INTEGER, PRIVATE :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
   INTEGER, PRIVATE :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
   LOGICAL, EXTERNAL :: rsl_comm_iter

   INTEGER, PRIVATE :: idim1, idim2, idim3, idim4, idim5, idim6, idim7


CONTAINS

   
   SUBROUTINE module_comm_dm_dummy_2
     USE module_domain, ONLY:domain
     USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
     USE module_state_description, ONLY:PARAM_FIRST_SCALAR
     USE module_driver_constants
     RETURN
   END SUBROUTINE module_comm_dm_dummy_2








SUBROUTINE HALO_XA_CLOUD_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_XA_CLOUD_inline.inc')
CALL rsl_comm_iter_init(2,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     7, 0, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xa%q,1)*SIZE(grid%xa%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%t,1)*SIZE(grid%xa%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qcw,1)*SIZE(grid%xa%qcw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qcw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qrn,1)*SIZE(grid%xa%qrn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qrn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qci,1)*SIZE(grid%xa%qci,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qci, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qsn,1)*SIZE(grid%xa%qsn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qsn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qgr,1)*SIZE(grid%xa%qgr,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qgr, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%q,1)*SIZE(grid%xa%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%t,1)*SIZE(grid%xa%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qcw,1)*SIZE(grid%xa%qcw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qcw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qrn,1)*SIZE(grid%xa%qrn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qrn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qci,1)*SIZE(grid%xa%qci,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qci, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qsn,1)*SIZE(grid%xa%qsn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qsn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qgr,1)*SIZE(grid%xa%qgr,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qgr, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
ENDDO
CALL rsl_comm_iter_init(2,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     7, 0, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xa%q,1)*SIZE(grid%xa%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%t,1)*SIZE(grid%xa%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qcw,1)*SIZE(grid%xa%qcw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qcw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qrn,1)*SIZE(grid%xa%qrn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qrn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qci,1)*SIZE(grid%xa%qci,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qci, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qsn,1)*SIZE(grid%xa%qsn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qsn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qgr,1)*SIZE(grid%xa%qgr,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qgr, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%q,1)*SIZE(grid%xa%q,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%t,1)*SIZE(grid%xa%t,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qcw,1)*SIZE(grid%xa%qcw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qcw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qrn,1)*SIZE(grid%xa%qrn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qrn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qci,1)*SIZE(grid%xa%qci,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qci, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qsn,1)*SIZE(grid%xa%qsn,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qsn, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%qgr,1)*SIZE(grid%xa%qgr,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%qgr, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_XA_CLOUD_sub







SUBROUTINE HALO_SFC_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_SFC_XA_inline.inc')
CALL rsl_comm_iter_init(2,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     0, 4, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
))
IF ( SIZE(grid%xa%u10,1)*SIZE(grid%xa%u10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%u10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%v10,1)*SIZE(grid%xa%v10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%v10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%t2,1)*SIZE(grid%xa%t2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%q2,1)*SIZE(grid%xa%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%u10,1)*SIZE(grid%xa%u10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%u10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%v10,1)*SIZE(grid%xa%v10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%v10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%t2,1)*SIZE(grid%xa%t2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%q2,1)*SIZE(grid%xa%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
ENDDO
CALL rsl_comm_iter_init(2,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     0, 4, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
))
IF ( SIZE(grid%xa%u10,1)*SIZE(grid%xa%u10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%u10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%v10,1)*SIZE(grid%xa%v10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%v10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%t2,1)*SIZE(grid%xa%t2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%q2,1)*SIZE(grid%xa%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%u10,1)*SIZE(grid%xa%u10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%u10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%v10,1)*SIZE(grid%xa%v10,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%v10, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%t2,1)*SIZE(grid%xa%t2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%t2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%q2,1)*SIZE(grid%xa%q2,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%q2, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_SFC_XA_sub







SUBROUTINE HALO_SSMI_XA_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_SSMI_XA_inline.inc')
CALL rsl_comm_iter_init(2,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     1, 10, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xa%ztd,1)*SIZE(grid%xa%ztd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ztd, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tpw,1)*SIZE(grid%xa%tpw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tpw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%speed,1)*SIZE(grid%xa%speed,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%speed, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19v,1)*SIZE(grid%xa%tb19v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19h,1)*SIZE(grid%xa%tb19h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb22v,1)*SIZE(grid%xa%tb22v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb22v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37v,1)*SIZE(grid%xa%tb37v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37h,1)*SIZE(grid%xa%tb37h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85v,1)*SIZE(grid%xa%tb85v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85h,1)*SIZE(grid%xa%tb85h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%ref,1)*SIZE(grid%xa%ref,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ref, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%ztd,1)*SIZE(grid%xa%ztd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ztd, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tpw,1)*SIZE(grid%xa%tpw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tpw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%speed,1)*SIZE(grid%xa%speed,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%speed, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19v,1)*SIZE(grid%xa%tb19v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19h,1)*SIZE(grid%xa%tb19h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb22v,1)*SIZE(grid%xa%tb22v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb22v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37v,1)*SIZE(grid%xa%tb37v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37h,1)*SIZE(grid%xa%tb37h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85v,1)*SIZE(grid%xa%tb85v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85h,1)*SIZE(grid%xa%tb85h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%ref,1)*SIZE(grid%xa%ref,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ref, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
ENDDO
CALL rsl_comm_iter_init(2,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     1, 10, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xa%ztd,1)*SIZE(grid%xa%ztd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ztd, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tpw,1)*SIZE(grid%xa%tpw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tpw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%speed,1)*SIZE(grid%xa%speed,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%speed, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19v,1)*SIZE(grid%xa%tb19v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19h,1)*SIZE(grid%xa%tb19h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb22v,1)*SIZE(grid%xa%tb22v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb22v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37v,1)*SIZE(grid%xa%tb37v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37h,1)*SIZE(grid%xa%tb37h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85v,1)*SIZE(grid%xa%tb85v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85h,1)*SIZE(grid%xa%tb85h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%ref,1)*SIZE(grid%xa%ref,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ref, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%ztd,1)*SIZE(grid%xa%ztd,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ztd, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tpw,1)*SIZE(grid%xa%tpw,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tpw, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%speed,1)*SIZE(grid%xa%speed,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%speed, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19v,1)*SIZE(grid%xa%tb19v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb19h,1)*SIZE(grid%xa%tb19h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb19h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb22v,1)*SIZE(grid%xa%tb22v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb22v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37v,1)*SIZE(grid%xa%tb37v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb37h,1)*SIZE(grid%xa%tb37h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb37h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85v,1)*SIZE(grid%xa%tb85v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%tb85h,1)*SIZE(grid%xa%tb85h,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%tb85h, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%xa%ref,1)*SIZE(grid%xa%ref,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%ref, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_SSMI_XA_sub







SUBROUTINE HALO_2D_WORK_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_2D_WORK_inline.inc')
CALL rsl_comm_iter_init(2,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     0, 1, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
))
IF ( SIZE(grid%xp%vxy,1)*SIZE(grid%xp%vxy,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xp%vxy, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xp%vxy,1)*SIZE(grid%xp%vxy,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xp%vxy, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
ENDDO
CALL rsl_comm_iter_init(2,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     0, 1, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
))
IF ( SIZE(grid%xp%vxy,1)*SIZE(grid%xp%vxy,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xp%vxy, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xp%vxy,1)*SIZE(grid%xp%vxy,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xp%vxy, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_2D_WORK_sub







SUBROUTINE HALO_RADAR_XA_W_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_RADAR_XA_W_inline.inc')
CALL rsl_comm_iter_init(2,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     2, 0, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xa%w,1)*SIZE(grid%xa%w,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%w, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%wh,1)*SIZE(grid%xa%wh,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%wh, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%w,1)*SIZE(grid%xa%w,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%w, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%wh,1)*SIZE(grid%xa%wh,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%wh, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
ENDDO
CALL rsl_comm_iter_init(2,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     2, 0, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xa%w,1)*SIZE(grid%xa%w,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%w, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%wh,1)*SIZE(grid%xa%wh,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%wh, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xa%w,1)*SIZE(grid%xa%w,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%w, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xa%wh,1)*SIZE(grid%xa%wh,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xa%wh, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_RADAR_XA_W_sub







SUBROUTINE HALO_XB_UV_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_XB_UV_inline.inc')
CALL rsl_comm_iter_init(2,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     2, 0, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xb%u,1)*SIZE(grid%xb%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%u, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xb%v,1)*SIZE(grid%xb%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xb%u,1)*SIZE(grid%xb%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%u, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xb%v,1)*SIZE(grid%xb%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
ENDDO
CALL rsl_comm_iter_init(2,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 2 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 2, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     2, 0, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%xb%u,1)*SIZE(grid%xb%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%u, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xb%v,1)*SIZE(grid%xb%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%xb%u,1)*SIZE(grid%xb%u,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%u, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%xb%v,1)*SIZE(grid%xb%v,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%xb%v, 2,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_XB_UV_sub







SUBROUTINE HALO_EM_BDY_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )
  USE module_domain, ONLY:domain
  USE module_configure, ONLY:grid_config_rec_type,in_use_for_config
  USE module_state_description, ONLY:PARAM_FIRST_SCALAR
  USE module_driver_constants
  TYPE(domain) ,               INTENT(IN) :: grid
  INTEGER ,                    INTENT(IN) :: local_communicator
  INTEGER ,                    INTENT(IN) :: mytask, ntasks, ntasks_x, ntasks_y
  INTEGER ,                    INTENT(IN) :: ids, ide, jds, jde, kds, kde
  INTEGER ,                    INTENT(IN) :: ims, ime, jms, jme, kms, kme
  INTEGER ,                    INTENT(IN) :: ips, ipe, jps, jpe, kps, kpe
  INTEGER :: itrace
  INTEGER :: rsl_sendw_p, rsl_sendbeg_p, rsl_recvw_p, rsl_recvbeg_p
  INTEGER :: rsl_sendw_m, rsl_sendbeg_m, rsl_recvw_m, rsl_recvbeg_m
  LOGICAL, EXTERNAL :: rsl_comm_iter
  INTEGER :: idim1, idim2, idim3, idim4, idim5, idim6, idim7
  
CALL push_communicators_for_domain( grid%id )






CALL wrf_debug(2,'calling inc/HALO_EM_BDY_inline.inc')
CALL rsl_comm_iter_init(1,jps,jpe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 1 , &
                         0 , jds,jde,jps,jpe, grid%njds, grid%njde , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 1, 0, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     5, 1, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%a_u_1,1)*SIZE(grid%a_u_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_u_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_v_1,1)*SIZE(grid%a_v_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_v_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_ph_1,1)*SIZE(grid%a_ph_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_ph_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_mu_1,1)*SIZE(grid%a_mu_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_mu_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%a_t_1,1)*SIZE(grid%a_t_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_t_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_w_1,1)*SIZE(grid%a_w_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_w_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%a_u_1,1)*SIZE(grid%a_u_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_u_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_v_1,1)*SIZE(grid%a_v_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_v_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_ph_1,1)*SIZE(grid%a_ph_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_ph_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_mu_1,1)*SIZE(grid%a_mu_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_mu_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%a_t_1,1)*SIZE(grid%a_t_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_t_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_w_1,1)*SIZE(grid%a_w_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_w_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 0, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
ENDDO
CALL rsl_comm_iter_init(1,ips,ipe)
DO WHILE ( rsl_comm_iter( grid%id , grid%is_intermediate, 1 , &
                         1 , ids,ide,ips,ipe, grid%nids, grid%nide , & 
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p    ))
 CALL RSL_LITE_INIT_EXCH ( local_communicator, 1, 1, &
     rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p,   & 
     rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p,   & 
     5, 1, 8, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
IF ( SIZE(grid%a_u_1,1)*SIZE(grid%a_u_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_u_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_v_1,1)*SIZE(grid%a_v_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_v_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_ph_1,1)*SIZE(grid%a_ph_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_ph_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_mu_1,1)*SIZE(grid%a_mu_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_mu_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%a_t_1,1)*SIZE(grid%a_t_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_t_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_w_1,1)*SIZE(grid%a_w_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_w_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 0, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y, &
                          rsl_sendw_m,  rsl_sendw_p, rsl_recvw_m,  rsl_recvw_p    )
IF ( SIZE(grid%a_u_1,1)*SIZE(grid%a_u_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_u_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_v_1,1)*SIZE(grid%a_v_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_v_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_ph_1,1)*SIZE(grid%a_ph_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_ph_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_mu_1,1)*SIZE(grid%a_mu_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_mu_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XY, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
ENDIF
IF ( SIZE(grid%a_t_1,1)*SIZE(grid%a_t_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_t_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
IF ( SIZE(grid%a_w_1,1)*SIZE(grid%a_w_1,2) .GT. 1 ) THEN
CALL RSL_LITE_PACK ( local_communicator,&
 grid%a_w_1, 1,&
rsl_sendbeg_m, rsl_sendw_m, rsl_sendbeg_p, rsl_sendw_p, &
rsl_recvbeg_m, rsl_recvw_m, rsl_recvbeg_p, rsl_recvw_p, &
8, 1, 1, DATA_ORDER_XYZ, 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDIF
    ENDDO

CALL pop_communicators_for_domain
  
  END SUBROUTINE HALO_EM_BDY_sub


END MODULE module_comm_dm_2

