!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nest_interpdown_interp_vert.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF ( SIZE( grid%u_2, 1 ) * SIZE( grid%u_2, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%u_2, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%v_2, 1 ) * SIZE( grid%v_2, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%v_2, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%w_2, 1 ) * SIZE( grid%w_2, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting_w( &
                                  grid%w_2, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  alt_w_c, alt_w_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%ph_2, 1 ) * SIZE( grid%ph_2, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting_w( &
                                  grid%ph_2, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  alt_w_c, alt_w_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%phb, 1 ) * SIZE( grid%phb, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting_w( &
                                  grid%phb, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  alt_w_c, alt_w_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%t_2, 1 ) * SIZE( grid%t_2, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%t_2, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%alb, 1 ) * SIZE( grid%alb, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%alb, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%pb, 1 ) * SIZE( grid%pb, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%pb, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_moist
IF ( SIZE( moist, 1 ) * SIZE( moist, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  moist(grid%sm31,grid%sm32,grid%sm33,itrace), & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
IF ( SIZE( scalar, 1 ) * SIZE( scalar, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  scalar(grid%sm31,grid%sm32,grid%sm33,itrace), & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
ENDDO
IF ( SIZE( grid%g_h_diabatic, 1 ) * SIZE( grid%g_h_diabatic, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%g_h_diabatic, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
IF ( SIZE( grid%a_h_diabatic, 1 ) * SIZE( grid%a_h_diabatic, 2 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  grid%a_h_diabatic, & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
DO itrace = PARAM_FIRST_SCALAR, num_tracer
IF ( SIZE( tracer, 1 ) * SIZE( tracer, 3 ) .GT. 1 ) THEN 
    CALL vert_interp_vert_nesting( &
                                  tracer(grid%sm31,grid%sm32,grid%sm33,itrace), & !CD field
                                  ids, ide, kds, kde, jds, jde, & !CD dims
                                  ims, ime, kms, kme, jms, jme, & !CD dims
                                  ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                  pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                  pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                  alt_u_c, alt_u_n ) !coordinates for parent and nest
ENDIF
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
