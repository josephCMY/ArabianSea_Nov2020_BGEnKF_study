!STARTOFREGISTRYGENERATEDINCLUDE 'inc/scalar_indices.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  P_qv = 1 ; F_qv = .FALSE. 
  P_qc = 1 ; F_qc = .FALSE. 
  P_qr = 1 ; F_qr = .FALSE. 
  P_qi = 1 ; F_qi = .FALSE. 
  P_qs = 1 ; F_qs = .FALSE. 
  P_qg = 1 ; F_qg = .FALSE. 
  P_qh = 1 ; F_qh = .FALSE. 
  P_qndrop = 1 ; F_qndrop = .FALSE. 
  P_qni = 1 ; F_qni = .FALSE. 
  P_qt = 1 ; F_qt = .FALSE. 
  P_qns = 1 ; F_qns = .FALSE. 
  P_qnr = 1 ; F_qnr = .FALSE. 
  P_qng = 1 ; F_qng = .FALSE. 
  P_qnh = 1 ; F_qnh = .FALSE. 
  P_qnn = 1 ; F_qnn = .FALSE. 
  P_qnc = 1 ; F_qnc = .FALSE. 
  P_qvolg = 1 ; F_qvolg = .FALSE. 
  P_a_qv = 1 ; F_a_qv = .FALSE. 
  P_a_qc = 1 ; F_a_qc = .FALSE. 
  P_a_qr = 1 ; F_a_qr = .FALSE. 
  P_a_qi = 1 ; F_a_qi = .FALSE. 
  P_a_qs = 1 ; F_a_qs = .FALSE. 
  P_a_qg = 1 ; F_a_qg = .FALSE. 
  P_a_qh = 1 ; F_a_qh = .FALSE. 
  P_g_qv = 1 ; F_g_qv = .FALSE. 
  P_g_qc = 1 ; F_g_qc = .FALSE. 
  P_g_qr = 1 ; F_g_qr = .FALSE. 
  P_g_qi = 1 ; F_g_qi = .FALSE. 
  P_g_qs = 1 ; F_g_qs = .FALSE. 
  P_g_qg = 1 ; F_g_qg = .FALSE. 
  P_g_qh = 1 ; F_g_qh = .FALSE. 
  P_tr17_1 = 1 ; F_tr17_1 = .FALSE. 
  P_tr17_2 = 1 ; F_tr17_2 = .FALSE. 
  P_tr17_3 = 1 ; F_tr17_3 = .FALSE. 
  P_tr17_4 = 1 ; F_tr17_4 = .FALSE. 
  P_tr17_5 = 1 ; F_tr17_5 = .FALSE. 
  P_tr17_6 = 1 ; F_tr17_6 = .FALSE. 
  P_tr17_7 = 1 ; F_tr17_7 = .FALSE. 
  P_tr17_8 = 1 ; F_tr17_8 = .FALSE. 
  IF (model_config_rec%sf_surface_physics(idomain)==1)THEN
  END IF
  IF (model_config_rec%sf_surface_physics(idomain)==2)THEN
  END IF
  IF (model_config_rec%sf_surface_physics(idomain)==3)THEN
  END IF
  IF (model_config_rec%sf_surface_physics(idomain)==4)THEN
  END IF
  IF (model_config_rec%dfi_stage(idomain)==0)THEN
  END IF
  IF (model_config_rec%dfi_stage(idomain)==1)THEN
  END IF
  IF (model_config_rec%dfi_stage(idomain)==2)THEN
  END IF
  IF (model_config_rec%dfi_stage(idomain)==3)THEN
  END IF
  IF (model_config_rec%dfi_stage(idomain)==4)THEN
  END IF
  IF (model_config_rec%dfi_stage(idomain)==5)THEN
  END IF
  IF (model_config_rec%dfi_opt==0)THEN
  END IF
  IF (model_config_rec%dfi_opt==1)THEN
  END IF
  IF (model_config_rec%dfi_opt==2)THEN
  END IF
  IF (model_config_rec%dfi_opt==3)THEN
  END IF
  IF (model_config_rec%io_form_restart==1)THEN
  END IF
  IF (model_config_rec%io_form_restart==2)THEN
  END IF
  IF (model_config_rec%io_form_restart==3)THEN
  END IF
  IF (model_config_rec%io_form_restart==4)THEN
  END IF
  IF (model_config_rec%io_form_restart==5)THEN
  END IF
  IF (model_config_rec%io_form_restart==6)THEN
  END IF
  IF (model_config_rec%io_form_restart==7)THEN
  END IF
  IF (model_config_rec%io_form_restart==8)THEN
  END IF
  IF (model_config_rec%io_form_restart==9)THEN
  END IF
  IF (model_config_rec%io_form_restart==10)THEN
  END IF
  IF (model_config_rec%io_form_restart==11)THEN
  END IF
  IF (model_config_rec%io_form_restart==12)THEN
  END IF
  IF (model_config_rec%wrf_hydro==0)THEN
  END IF
  IF (model_config_rec%wrf_hydro==1)THEN
  END IF
  IF (model_config_rec%dyn_opt==0)THEN
  END IF
  IF (model_config_rec%dyn_opt==2)THEN
  END IF
  IF (model_config_rec%dyn_opt==102)THEN
  END IF
  IF (model_config_rec%dyn_opt==202)THEN
  END IF
  IF (model_config_rec%dyn_opt==302)THEN
  END IF
  IF (model_config_rec%dyn_opt==402)THEN
  END IF
  IF (model_config_rec%dyn_opt==502)THEN
  END IF
  IF (model_config_rec%mp_physics(idomain)==0)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==1)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==2)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==3)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==4)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==5)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( scalar_index_table( PARAM_qt , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qt = scalar_num_table(idomain)
     scalar_index_table( PARAM_qt , idomain ) = P_qt
   ELSE
     P_qt = scalar_index_table( PARAM_qt , idomain )
   END IF
   scalar_boundary_table( idomain, P_qt ) = .TRUE.
   scalar_dname_table( idomain, P_qt ) = 'CWM'
   scalar_desc_table( idomain, P_qt ) = 'Total condensate mixing ratio'
   scalar_units_table( idomain, P_qt ) = 'kg kg-1'
   scalar_streams_table( idomain, P_qt )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qt )%stream(2) = 2097152 ! 00200000 
   F_qt = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==15)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==6)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==7)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==8)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==9)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( moist_index_table( PARAM_qh , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qh = moist_num_table(idomain)
     moist_index_table( PARAM_qh , idomain ) = P_qh
   ELSE
     P_qh = moist_index_table( PARAM_qh , idomain )
   END IF
   moist_boundary_table( idomain, P_qh ) = .TRUE.
   moist_dname_table( idomain, P_qh ) = 'QHAIL'
   moist_desc_table( idomain, P_qh ) = 'Hail mixing ratio'
   moist_units_table( idomain, P_qh ) = 'kg kg-1'
   moist_streams_table( idomain, P_qh )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qh )%stream(2) = 2097152 ! 00200000 
   F_qh = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qh , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qh = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qh , idomain ) = P_g_qh
   ELSE
     P_g_qh = g_moist_index_table( PARAM_g_qh , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qh ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qh ) = 'G_QHAIL'
   g_moist_desc_table( idomain, P_g_qh ) = 'Hail mixing ratio'
   g_moist_units_table( idomain, P_g_qh ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qh )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qh )%stream(2) = 2097152 ! 00200000 
   F_g_qh = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qh , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qh = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qh , idomain ) = P_a_qh
   ELSE
     P_a_qh = a_moist_index_table( PARAM_a_qh , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qh ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qh ) = 'A_QHAIL'
   a_moist_desc_table( idomain, P_a_qh ) = 'Hail mixing ratio'
   a_moist_units_table( idomain, P_a_qh ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qh )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qh )%stream(2) = 2097152 ! 00200000 
   F_a_qh = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qns , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qns = scalar_num_table(idomain)
     scalar_index_table( PARAM_qns , idomain ) = P_qns
   ELSE
     P_qns = scalar_index_table( PARAM_qns , idomain )
   END IF
   scalar_boundary_table( idomain, P_qns ) = .TRUE.
   scalar_dname_table( idomain, P_qns ) = 'QNSNOW'
   scalar_desc_table( idomain, P_qns ) = 'Snow Number concentration'
   scalar_units_table( idomain, P_qns ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qns )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qns )%stream(2) = 2097152 ! 00200000 
   F_qns = .TRUE.
   IF ( scalar_index_table( PARAM_qng , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qng = scalar_num_table(idomain)
     scalar_index_table( PARAM_qng , idomain ) = P_qng
   ELSE
     P_qng = scalar_index_table( PARAM_qng , idomain )
   END IF
   scalar_boundary_table( idomain, P_qng ) = .TRUE.
   scalar_dname_table( idomain, P_qng ) = 'QNGRAUPEL'
   scalar_desc_table( idomain, P_qng ) = 'Graupel Number concentration'
   scalar_units_table( idomain, P_qng ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qng )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qng )%stream(2) = 2097152 ! 00200000 
   F_qng = .TRUE.
   IF ( scalar_index_table( PARAM_qnh , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnh = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnh , idomain ) = P_qnh
   ELSE
     P_qnh = scalar_index_table( PARAM_qnh , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnh ) = .TRUE.
   scalar_dname_table( idomain, P_qnh ) = 'QNHAIL'
   scalar_desc_table( idomain, P_qnh ) = 'Hail Number concentration'
   scalar_units_table( idomain, P_qnh ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnh )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnh )%stream(2) = 2097152 ! 00200000 
   F_qnh = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==10)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qns , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qns = scalar_num_table(idomain)
     scalar_index_table( PARAM_qns , idomain ) = P_qns
   ELSE
     P_qns = scalar_index_table( PARAM_qns , idomain )
   END IF
   scalar_boundary_table( idomain, P_qns ) = .TRUE.
   scalar_dname_table( idomain, P_qns ) = 'QNSNOW'
   scalar_desc_table( idomain, P_qns ) = 'Snow Number concentration'
   scalar_units_table( idomain, P_qns ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qns )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qns )%stream(2) = 2097152 ! 00200000 
   F_qns = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qng , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qng = scalar_num_table(idomain)
     scalar_index_table( PARAM_qng , idomain ) = P_qng
   ELSE
     P_qng = scalar_index_table( PARAM_qng , idomain )
   END IF
   scalar_boundary_table( idomain, P_qng ) = .TRUE.
   scalar_dname_table( idomain, P_qng ) = 'QNGRAUPEL'
   scalar_desc_table( idomain, P_qng ) = 'Graupel Number concentration'
   scalar_units_table( idomain, P_qng ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qng )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qng )%stream(2) = 2097152 ! 00200000 
   F_qng = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==11)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qns , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qns = scalar_num_table(idomain)
     scalar_index_table( PARAM_qns , idomain ) = P_qns
   ELSE
     P_qns = scalar_index_table( PARAM_qns , idomain )
   END IF
   scalar_boundary_table( idomain, P_qns ) = .TRUE.
   scalar_dname_table( idomain, P_qns ) = 'QNSNOW'
   scalar_desc_table( idomain, P_qns ) = 'Snow Number concentration'
   scalar_units_table( idomain, P_qns ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qns )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qns )%stream(2) = 2097152 ! 00200000 
   F_qns = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==13)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==14)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( scalar_index_table( PARAM_qnn , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnn = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnn , idomain ) = P_qnn
   ELSE
     P_qnn = scalar_index_table( PARAM_qnn , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnn ) = .TRUE.
   scalar_dname_table( idomain, P_qnn ) = 'QNCCN'
   scalar_desc_table( idomain, P_qnn ) = 'CCN Number concentration'
   scalar_units_table( idomain, P_qnn ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnn )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnn )%stream(2) = 2097152 ! 00200000 
   F_qnn = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==16)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( scalar_index_table( PARAM_qnn , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnn = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnn , idomain ) = P_qnn
   ELSE
     P_qnn = scalar_index_table( PARAM_qnn , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnn ) = .TRUE.
   scalar_dname_table( idomain, P_qnn ) = 'QNCCN'
   scalar_desc_table( idomain, P_qnn ) = 'CCN Number concentration'
   scalar_units_table( idomain, P_qnn ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnn )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnn )%stream(2) = 2097152 ! 00200000 
   F_qnn = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==17)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( moist_index_table( PARAM_qh , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qh = moist_num_table(idomain)
     moist_index_table( PARAM_qh , idomain ) = P_qh
   ELSE
     P_qh = moist_index_table( PARAM_qh , idomain )
   END IF
   moist_boundary_table( idomain, P_qh ) = .TRUE.
   moist_dname_table( idomain, P_qh ) = 'QHAIL'
   moist_desc_table( idomain, P_qh ) = 'Hail mixing ratio'
   moist_units_table( idomain, P_qh ) = 'kg kg-1'
   moist_streams_table( idomain, P_qh )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qh )%stream(2) = 2097152 ! 00200000 
   F_qh = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qh , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qh = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qh , idomain ) = P_g_qh
   ELSE
     P_g_qh = g_moist_index_table( PARAM_g_qh , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qh ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qh ) = 'G_QHAIL'
   g_moist_desc_table( idomain, P_g_qh ) = 'Hail mixing ratio'
   g_moist_units_table( idomain, P_g_qh ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qh )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qh )%stream(2) = 2097152 ! 00200000 
   F_g_qh = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qh , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qh = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qh , idomain ) = P_a_qh
   ELSE
     P_a_qh = a_moist_index_table( PARAM_a_qh , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qh ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qh ) = 'A_QHAIL'
   a_moist_desc_table( idomain, P_a_qh ) = 'Hail mixing ratio'
   a_moist_units_table( idomain, P_a_qh ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qh )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qh )%stream(2) = 2097152 ! 00200000 
   F_a_qh = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qns , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qns = scalar_num_table(idomain)
     scalar_index_table( PARAM_qns , idomain ) = P_qns
   ELSE
     P_qns = scalar_index_table( PARAM_qns , idomain )
   END IF
   scalar_boundary_table( idomain, P_qns ) = .TRUE.
   scalar_dname_table( idomain, P_qns ) = 'QNSNOW'
   scalar_desc_table( idomain, P_qns ) = 'Snow Number concentration'
   scalar_units_table( idomain, P_qns ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qns )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qns )%stream(2) = 2097152 ! 00200000 
   F_qns = .TRUE.
   IF ( scalar_index_table( PARAM_qng , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qng = scalar_num_table(idomain)
     scalar_index_table( PARAM_qng , idomain ) = P_qng
   ELSE
     P_qng = scalar_index_table( PARAM_qng , idomain )
   END IF
   scalar_boundary_table( idomain, P_qng ) = .TRUE.
   scalar_dname_table( idomain, P_qng ) = 'QNGRAUPEL'
   scalar_desc_table( idomain, P_qng ) = 'Graupel Number concentration'
   scalar_units_table( idomain, P_qng ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qng )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qng )%stream(2) = 2097152 ! 00200000 
   F_qng = .TRUE.
   IF ( scalar_index_table( PARAM_qnh , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnh = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnh , idomain ) = P_qnh
   ELSE
     P_qnh = scalar_index_table( PARAM_qnh , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnh ) = .TRUE.
   scalar_dname_table( idomain, P_qnh ) = 'QNHAIL'
   scalar_desc_table( idomain, P_qnh ) = 'Hail Number concentration'
   scalar_units_table( idomain, P_qnh ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnh )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnh )%stream(2) = 2097152 ! 00200000 
   F_qnh = .TRUE.
   IF ( scalar_index_table( PARAM_qvolg , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qvolg = scalar_num_table(idomain)
     scalar_index_table( PARAM_qvolg , idomain ) = P_qvolg
   ELSE
     P_qvolg = scalar_index_table( PARAM_qvolg , idomain )
   END IF
   scalar_boundary_table( idomain, P_qvolg ) = .TRUE.
   scalar_dname_table( idomain, P_qvolg ) = 'QVGRAUPEL'
   scalar_desc_table( idomain, P_qvolg ) = 'Graupel Particle Volume'
   scalar_units_table( idomain, P_qvolg ) = 'm(3) kg(-1)'
   scalar_streams_table( idomain, P_qvolg )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qvolg )%stream(2) = 2097152 ! 00200000 
   F_qvolg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==18)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( moist_index_table( PARAM_qh , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qh = moist_num_table(idomain)
     moist_index_table( PARAM_qh , idomain ) = P_qh
   ELSE
     P_qh = moist_index_table( PARAM_qh , idomain )
   END IF
   moist_boundary_table( idomain, P_qh ) = .TRUE.
   moist_dname_table( idomain, P_qh ) = 'QHAIL'
   moist_desc_table( idomain, P_qh ) = 'Hail mixing ratio'
   moist_units_table( idomain, P_qh ) = 'kg kg-1'
   moist_streams_table( idomain, P_qh )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qh )%stream(2) = 2097152 ! 00200000 
   F_qh = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qh , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qh = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qh , idomain ) = P_g_qh
   ELSE
     P_g_qh = g_moist_index_table( PARAM_g_qh , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qh ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qh ) = 'G_QHAIL'
   g_moist_desc_table( idomain, P_g_qh ) = 'Hail mixing ratio'
   g_moist_units_table( idomain, P_g_qh ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qh )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qh )%stream(2) = 2097152 ! 00200000 
   F_g_qh = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qh , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qh = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qh , idomain ) = P_a_qh
   ELSE
     P_a_qh = a_moist_index_table( PARAM_a_qh , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qh ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qh ) = 'A_QHAIL'
   a_moist_desc_table( idomain, P_a_qh ) = 'Hail mixing ratio'
   a_moist_units_table( idomain, P_a_qh ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qh )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qh )%stream(2) = 2097152 ! 00200000 
   F_a_qh = .TRUE.
   IF ( scalar_index_table( PARAM_qnn , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnn = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnn , idomain ) = P_qnn
   ELSE
     P_qnn = scalar_index_table( PARAM_qnn , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnn ) = .TRUE.
   scalar_dname_table( idomain, P_qnn ) = 'QNCCN'
   scalar_desc_table( idomain, P_qnn ) = 'CCN Number concentration'
   scalar_units_table( idomain, P_qnn ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnn )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnn )%stream(2) = 2097152 ! 00200000 
   F_qnn = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qns , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qns = scalar_num_table(idomain)
     scalar_index_table( PARAM_qns , idomain ) = P_qns
   ELSE
     P_qns = scalar_index_table( PARAM_qns , idomain )
   END IF
   scalar_boundary_table( idomain, P_qns ) = .TRUE.
   scalar_dname_table( idomain, P_qns ) = 'QNSNOW'
   scalar_desc_table( idomain, P_qns ) = 'Snow Number concentration'
   scalar_units_table( idomain, P_qns ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qns )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qns )%stream(2) = 2097152 ! 00200000 
   F_qns = .TRUE.
   IF ( scalar_index_table( PARAM_qng , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qng = scalar_num_table(idomain)
     scalar_index_table( PARAM_qng , idomain ) = P_qng
   ELSE
     P_qng = scalar_index_table( PARAM_qng , idomain )
   END IF
   scalar_boundary_table( idomain, P_qng ) = .TRUE.
   scalar_dname_table( idomain, P_qng ) = 'QNGRAUPEL'
   scalar_desc_table( idomain, P_qng ) = 'Graupel Number concentration'
   scalar_units_table( idomain, P_qng ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qng )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qng )%stream(2) = 2097152 ! 00200000 
   F_qng = .TRUE.
   IF ( scalar_index_table( PARAM_qnh , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnh = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnh , idomain ) = P_qnh
   ELSE
     P_qnh = scalar_index_table( PARAM_qnh , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnh ) = .TRUE.
   scalar_dname_table( idomain, P_qnh ) = 'QNHAIL'
   scalar_desc_table( idomain, P_qnh ) = 'Hail Number concentration'
   scalar_units_table( idomain, P_qnh ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnh )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnh )%stream(2) = 2097152 ! 00200000 
   F_qnh = .TRUE.
   IF ( scalar_index_table( PARAM_qvolg , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qvolg = scalar_num_table(idomain)
     scalar_index_table( PARAM_qvolg , idomain ) = P_qvolg
   ELSE
     P_qvolg = scalar_index_table( PARAM_qvolg , idomain )
   END IF
   scalar_boundary_table( idomain, P_qvolg ) = .TRUE.
   scalar_dname_table( idomain, P_qvolg ) = 'QVGRAUPEL'
   scalar_desc_table( idomain, P_qvolg ) = 'Graupel Particle Volume'
   scalar_units_table( idomain, P_qvolg ) = 'm(3) kg(-1)'
   scalar_streams_table( idomain, P_qvolg )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qvolg )%stream(2) = 2097152 ! 00200000 
   F_qvolg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==19)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( moist_index_table( PARAM_qh , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qh = moist_num_table(idomain)
     moist_index_table( PARAM_qh , idomain ) = P_qh
   ELSE
     P_qh = moist_index_table( PARAM_qh , idomain )
   END IF
   moist_boundary_table( idomain, P_qh ) = .TRUE.
   moist_dname_table( idomain, P_qh ) = 'QHAIL'
   moist_desc_table( idomain, P_qh ) = 'Hail mixing ratio'
   moist_units_table( idomain, P_qh ) = 'kg kg-1'
   moist_streams_table( idomain, P_qh )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qh )%stream(2) = 2097152 ! 00200000 
   F_qh = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qh , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qh = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qh , idomain ) = P_g_qh
   ELSE
     P_g_qh = g_moist_index_table( PARAM_g_qh , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qh ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qh ) = 'G_QHAIL'
   g_moist_desc_table( idomain, P_g_qh ) = 'Hail mixing ratio'
   g_moist_units_table( idomain, P_g_qh ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qh )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qh )%stream(2) = 2097152 ! 00200000 
   F_g_qh = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qh , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qh = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qh , idomain ) = P_a_qh
   ELSE
     P_a_qh = a_moist_index_table( PARAM_a_qh , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qh ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qh ) = 'A_QHAIL'
   a_moist_desc_table( idomain, P_a_qh ) = 'Hail mixing ratio'
   a_moist_units_table( idomain, P_a_qh ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qh )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qh )%stream(2) = 2097152 ! 00200000 
   F_a_qh = .TRUE.
   IF ( scalar_index_table( PARAM_qvolg , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qvolg = scalar_num_table(idomain)
     scalar_index_table( PARAM_qvolg , idomain ) = P_qvolg
   ELSE
     P_qvolg = scalar_index_table( PARAM_qvolg , idomain )
   END IF
   scalar_boundary_table( idomain, P_qvolg ) = .TRUE.
   scalar_dname_table( idomain, P_qvolg ) = 'QVGRAUPEL'
   scalar_desc_table( idomain, P_qvolg ) = 'Graupel Particle Volume'
   scalar_units_table( idomain, P_qvolg ) = 'm(3) kg(-1)'
   scalar_streams_table( idomain, P_qvolg )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qvolg )%stream(2) = 2097152 ! 00200000 
   F_qvolg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==21)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==22)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( scalar_index_table( PARAM_qndrop , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qndrop = scalar_num_table(idomain)
     scalar_index_table( PARAM_qndrop , idomain ) = P_qndrop
   ELSE
     P_qndrop = scalar_index_table( PARAM_qndrop , idomain )
   END IF
   scalar_boundary_table( idomain, P_qndrop ) = .TRUE.
   scalar_dname_table( idomain, P_qndrop ) = 'QNDROP'
   scalar_desc_table( idomain, P_qndrop ) = 'Droplet number mixing ratio'
   scalar_units_table( idomain, P_qndrop ) = '  kg-1'
   scalar_streams_table( idomain, P_qndrop )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qndrop )%stream(2) = 2097152 ! 00200000 
   F_qndrop = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qns , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qns = scalar_num_table(idomain)
     scalar_index_table( PARAM_qns , idomain ) = P_qns
   ELSE
     P_qns = scalar_index_table( PARAM_qns , idomain )
   END IF
   scalar_boundary_table( idomain, P_qns ) = .TRUE.
   scalar_dname_table( idomain, P_qns ) = 'QNSNOW'
   scalar_desc_table( idomain, P_qns ) = 'Snow Number concentration'
   scalar_units_table( idomain, P_qns ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qns )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qns )%stream(2) = 2097152 ! 00200000 
   F_qns = .TRUE.
   IF ( scalar_index_table( PARAM_qng , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qng = scalar_num_table(idomain)
     scalar_index_table( PARAM_qng , idomain ) = P_qng
   ELSE
     P_qng = scalar_index_table( PARAM_qng , idomain )
   END IF
   scalar_boundary_table( idomain, P_qng ) = .TRUE.
   scalar_dname_table( idomain, P_qng ) = 'QNGRAUPEL'
   scalar_desc_table( idomain, P_qng ) = 'Graupel Number concentration'
   scalar_units_table( idomain, P_qng ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qng )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qng )%stream(2) = 2097152 ! 00200000 
   F_qng = .TRUE.
   IF ( scalar_index_table( PARAM_qvolg , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qvolg = scalar_num_table(idomain)
     scalar_index_table( PARAM_qvolg , idomain ) = P_qvolg
   ELSE
     P_qvolg = scalar_index_table( PARAM_qvolg , idomain )
   END IF
   scalar_boundary_table( idomain, P_qvolg ) = .TRUE.
   scalar_dname_table( idomain, P_qvolg ) = 'QVGRAUPEL'
   scalar_desc_table( idomain, P_qvolg ) = 'Graupel Particle Volume'
   scalar_units_table( idomain, P_qvolg ) = 'm(3) kg(-1)'
   scalar_streams_table( idomain, P_qvolg )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qvolg )%stream(2) = 2097152 ! 00200000 
   F_qvolg = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==28)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qi , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qi = moist_num_table(idomain)
     moist_index_table( PARAM_qi , idomain ) = P_qi
   ELSE
     P_qi = moist_index_table( PARAM_qi , idomain )
   END IF
   moist_boundary_table( idomain, P_qi ) = .TRUE.
   moist_dname_table( idomain, P_qi ) = 'QICE'
   moist_desc_table( idomain, P_qi ) = 'Ice mixing ratio'
   moist_units_table( idomain, P_qi ) = 'kg kg-1'
   moist_streams_table( idomain, P_qi )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qi )%stream(2) = 2097152 ! 00200000 
   F_qi = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( moist_index_table( PARAM_qg , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qg = moist_num_table(idomain)
     moist_index_table( PARAM_qg , idomain ) = P_qg
   ELSE
     P_qg = moist_index_table( PARAM_qg , idomain )
   END IF
   moist_boundary_table( idomain, P_qg ) = .TRUE.
   moist_dname_table( idomain, P_qg ) = 'QGRAUP'
   moist_desc_table( idomain, P_qg ) = 'Graupel mixing ratio'
   moist_units_table( idomain, P_qg ) = 'kg kg-1'
   moist_streams_table( idomain, P_qg )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qg )%stream(2) = 2097152 ! 00200000 
   F_qg = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qi , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qi = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qi , idomain ) = P_g_qi
   ELSE
     P_g_qi = g_moist_index_table( PARAM_g_qi , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qi ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qi ) = 'G_QICE'
   g_moist_desc_table( idomain, P_g_qi ) = 'Ice mixing ratio'
   g_moist_units_table( idomain, P_g_qi ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qi )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qi )%stream(2) = 2097152 ! 00200000 
   F_g_qi = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qg , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qg = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qg , idomain ) = P_g_qg
   ELSE
     P_g_qg = g_moist_index_table( PARAM_g_qg , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qg ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qg ) = 'G_QGRAUP'
   g_moist_desc_table( idomain, P_g_qg ) = 'Graupel mixing ratio'
   g_moist_units_table( idomain, P_g_qg ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qg )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qg )%stream(2) = 2097152 ! 00200000 
   F_g_qg = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qi , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qi = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qi , idomain ) = P_a_qi
   ELSE
     P_a_qi = a_moist_index_table( PARAM_a_qi , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qi ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qi ) = 'A_QICE'
   a_moist_desc_table( idomain, P_a_qi ) = 'Ice mixing ratio'
   a_moist_units_table( idomain, P_a_qi ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qi )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qi )%stream(2) = 2097152 ! 00200000 
   F_a_qi = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qg , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qg = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qg , idomain ) = P_a_qg
   ELSE
     P_a_qg = a_moist_index_table( PARAM_a_qg , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qg ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qg ) = 'A_QGRAUP'
   a_moist_desc_table( idomain, P_a_qg ) = 'Graupel mixing ratio'
   a_moist_units_table( idomain, P_a_qg ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qg )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qg )%stream(2) = 2097152 ! 00200000 
   F_a_qg = .TRUE.
   IF ( scalar_index_table( PARAM_qni , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qni = scalar_num_table(idomain)
     scalar_index_table( PARAM_qni , idomain ) = P_qni
   ELSE
     P_qni = scalar_index_table( PARAM_qni , idomain )
   END IF
   scalar_boundary_table( idomain, P_qni ) = .TRUE.
   scalar_dname_table( idomain, P_qni ) = 'QNICE'
   scalar_desc_table( idomain, P_qni ) = 'Ice Number concentration'
   scalar_units_table( idomain, P_qni ) = '  kg-1'
   scalar_streams_table( idomain, P_qni )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qni )%stream(2) = 2097152 ! 00200000 
   F_qni = .TRUE.
   IF ( scalar_index_table( PARAM_qnr , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnr = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnr , idomain ) = P_qnr
   ELSE
     P_qnr = scalar_index_table( PARAM_qnr , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnr ) = .TRUE.
   scalar_dname_table( idomain, P_qnr ) = 'QNRAIN'
   scalar_desc_table( idomain, P_qnr ) = 'Rain Number concentration'
   scalar_units_table( idomain, P_qnr ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnr )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnr )%stream(2) = 2097152 ! 00200000 
   F_qnr = .TRUE.
   IF ( scalar_index_table( PARAM_qnc , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qnc = scalar_num_table(idomain)
     scalar_index_table( PARAM_qnc , idomain ) = P_qnc
   ELSE
     P_qnc = scalar_index_table( PARAM_qnc , idomain )
   END IF
   scalar_boundary_table( idomain, P_qnc ) = .TRUE.
   scalar_dname_table( idomain, P_qnc ) = 'QNCLOUD'
   scalar_desc_table( idomain, P_qnc ) = 'cloud water Number concentration'
   scalar_units_table( idomain, P_qnc ) = '  kg(-1)'
   scalar_streams_table( idomain, P_qnc )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qnc )%stream(2) = 2097152 ! 00200000 
   F_qnc = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==95)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( moist_index_table( PARAM_qs , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qs = moist_num_table(idomain)
     moist_index_table( PARAM_qs , idomain ) = P_qs
   ELSE
     P_qs = moist_index_table( PARAM_qs , idomain )
   END IF
   moist_boundary_table( idomain, P_qs ) = .TRUE.
   moist_dname_table( idomain, P_qs ) = 'QSNOW'
   moist_desc_table( idomain, P_qs ) = 'Snow mixing ratio'
   moist_units_table( idomain, P_qs ) = 'kg kg-1'
   moist_streams_table( idomain, P_qs )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qs )%stream(2) = 2097152 ! 00200000 
   F_qs = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qs , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qs = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qs , idomain ) = P_g_qs
   ELSE
     P_g_qs = g_moist_index_table( PARAM_g_qs , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qs ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qs ) = 'G_QSNOW'
   g_moist_desc_table( idomain, P_g_qs ) = 'Snow mixing ratio'
   g_moist_units_table( idomain, P_g_qs ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qs )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qs )%stream(2) = 2097152 ! 00200000 
   F_g_qs = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qs , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qs = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qs , idomain ) = P_a_qs
   ELSE
     P_a_qs = a_moist_index_table( PARAM_a_qs , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qs ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qs ) = 'A_QSNOW'
   a_moist_desc_table( idomain, P_a_qs ) = 'Snow mixing ratio'
   a_moist_units_table( idomain, P_a_qs ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qs )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qs )%stream(2) = 2097152 ! 00200000 
   F_a_qs = .TRUE.
   IF ( scalar_index_table( PARAM_qt , idomain ) .lt. 1 ) THEN
     scalar_num_table(idomain) = scalar_num_table(idomain) + 1
     P_qt = scalar_num_table(idomain)
     scalar_index_table( PARAM_qt , idomain ) = P_qt
   ELSE
     P_qt = scalar_index_table( PARAM_qt , idomain )
   END IF
   scalar_boundary_table( idomain, P_qt ) = .TRUE.
   scalar_dname_table( idomain, P_qt ) = 'CWM'
   scalar_desc_table( idomain, P_qt ) = 'Total condensate mixing ratio'
   scalar_units_table( idomain, P_qt ) = 'kg kg-1'
   scalar_streams_table( idomain, P_qt )%stream(1) = 33554433 ! 02000001 
   scalar_streams_table( idomain, P_qt )%stream(2) = 2097152 ! 00200000 
   F_qt = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==98)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
  END IF
  IF (model_config_rec%mp_physics(idomain)==99)THEN
   IF ( moist_index_table( PARAM_qv , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qv = moist_num_table(idomain)
     moist_index_table( PARAM_qv , idomain ) = P_qv
   ELSE
     P_qv = moist_index_table( PARAM_qv , idomain )
   END IF
   moist_boundary_table( idomain, P_qv ) = .TRUE.
   moist_dname_table( idomain, P_qv ) = 'QVAPOR'
   moist_desc_table( idomain, P_qv ) = 'Water vapor mixing ratio'
   moist_units_table( idomain, P_qv ) = 'kg kg-1'
   moist_streams_table( idomain, P_qv )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qv )%stream(2) = 2097152 ! 00200000 
   F_qv = .TRUE.
   IF ( moist_index_table( PARAM_qc , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qc = moist_num_table(idomain)
     moist_index_table( PARAM_qc , idomain ) = P_qc
   ELSE
     P_qc = moist_index_table( PARAM_qc , idomain )
   END IF
   moist_boundary_table( idomain, P_qc ) = .TRUE.
   moist_dname_table( idomain, P_qc ) = 'QCLOUD'
   moist_desc_table( idomain, P_qc ) = 'Cloud water mixing ratio'
   moist_units_table( idomain, P_qc ) = 'kg kg-1'
   moist_streams_table( idomain, P_qc )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qc )%stream(2) = 2097152 ! 00200000 
   F_qc = .TRUE.
   IF ( moist_index_table( PARAM_qr , idomain ) .lt. 1 ) THEN
     moist_num_table(idomain) = moist_num_table(idomain) + 1
     P_qr = moist_num_table(idomain)
     moist_index_table( PARAM_qr , idomain ) = P_qr
   ELSE
     P_qr = moist_index_table( PARAM_qr , idomain )
   END IF
   moist_boundary_table( idomain, P_qr ) = .TRUE.
   moist_dname_table( idomain, P_qr ) = 'QRAIN'
   moist_desc_table( idomain, P_qr ) = 'Rain water mixing ratio'
   moist_units_table( idomain, P_qr ) = 'kg kg-1'
   moist_streams_table( idomain, P_qr )%stream(1) = 33554433 ! 02000001 
   moist_streams_table( idomain, P_qr )%stream(2) = 2097152 ! 00200000 
   F_qr = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qv , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qv = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qv , idomain ) = P_g_qv
   ELSE
     P_g_qv = g_moist_index_table( PARAM_g_qv , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qv ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qv ) = 'G_QVAPOR'
   g_moist_desc_table( idomain, P_g_qv ) = 'Water vapor mixing ratio'
   g_moist_units_table( idomain, P_g_qv ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qv )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qv )%stream(2) = 2097152 ! 00200000 
   F_g_qv = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qc , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qc = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qc , idomain ) = P_g_qc
   ELSE
     P_g_qc = g_moist_index_table( PARAM_g_qc , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qc ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qc ) = 'G_QCLOUD'
   g_moist_desc_table( idomain, P_g_qc ) = 'Cloud water mixing ratio'
   g_moist_units_table( idomain, P_g_qc ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qc )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qc )%stream(2) = 2097152 ! 00200000 
   F_g_qc = .TRUE.
   IF ( g_moist_index_table( PARAM_g_qr , idomain ) .lt. 1 ) THEN
     g_moist_num_table(idomain) = g_moist_num_table(idomain) + 1
     P_g_qr = g_moist_num_table(idomain)
     g_moist_index_table( PARAM_g_qr , idomain ) = P_g_qr
   ELSE
     P_g_qr = g_moist_index_table( PARAM_g_qr , idomain )
   END IF
   g_moist_boundary_table( idomain, P_g_qr ) = .FALSE.
   g_moist_dname_table( idomain, P_g_qr ) = 'G_QRAIN'
   g_moist_desc_table( idomain, P_g_qr ) = 'Rain water mixing ratio'
   g_moist_units_table( idomain, P_g_qr ) = 'kg kg-1'
   g_moist_streams_table( idomain, P_g_qr )%stream(1) = 1 ! 00000001 
   g_moist_streams_table( idomain, P_g_qr )%stream(2) = 2097152 ! 00200000 
   F_g_qr = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qv , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qv = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qv , idomain ) = P_a_qv
   ELSE
     P_a_qv = a_moist_index_table( PARAM_a_qv , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qv ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qv ) = 'A_QVAPOR'
   a_moist_desc_table( idomain, P_a_qv ) = 'Water vapor mixing ratio'
   a_moist_units_table( idomain, P_a_qv ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qv )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qv )%stream(2) = 2097152 ! 00200000 
   F_a_qv = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qc , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qc = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qc , idomain ) = P_a_qc
   ELSE
     P_a_qc = a_moist_index_table( PARAM_a_qc , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qc ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qc ) = 'A_QCLOUD'
   a_moist_desc_table( idomain, P_a_qc ) = 'Cloud water mixing ratio'
   a_moist_units_table( idomain, P_a_qc ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qc )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qc )%stream(2) = 2097152 ! 00200000 
   F_a_qc = .TRUE.
   IF ( a_moist_index_table( PARAM_a_qr , idomain ) .lt. 1 ) THEN
     a_moist_num_table(idomain) = a_moist_num_table(idomain) + 1
     P_a_qr = a_moist_num_table(idomain)
     a_moist_index_table( PARAM_a_qr , idomain ) = P_a_qr
   ELSE
     P_a_qr = a_moist_index_table( PARAM_a_qr , idomain )
   END IF
   a_moist_boundary_table( idomain, P_a_qr ) = .FALSE.
   a_moist_dname_table( idomain, P_a_qr ) = 'A_QRAIN'
   a_moist_desc_table( idomain, P_a_qr ) = 'Rain water mixing ratio'
   a_moist_units_table( idomain, P_a_qr ) = 'kg kg-1'
   a_moist_streams_table( idomain, P_a_qr )%stream(1) = 1 ! 00000001 
   a_moist_streams_table( idomain, P_a_qr )%stream(2) = 2097152 ! 00200000 
   F_a_qr = .TRUE.
  END IF
  IF (model_config_rec%bl_pbl_physics(idomain)==98)THEN
  END IF
  IF (model_config_rec%cu_physics(idomain)==98)THEN
  END IF
  IF (model_config_rec%tracer_opt(idomain)==2)THEN
   IF ( tracer_index_table( PARAM_tr17_1 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_1 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_1 , idomain ) = P_tr17_1
   ELSE
     P_tr17_1 = tracer_index_table( PARAM_tr17_1 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_1 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_1 ) = 'tr17_1'
   tracer_desc_table( idomain, P_tr17_1 ) = 'tr17_1'
   tracer_units_table( idomain, P_tr17_1 ) = '-'
   tracer_streams_table( idomain, P_tr17_1 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_1 )%stream(2) = 2097152 ! 00200000 
   F_tr17_1 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_2 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_2 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_2 , idomain ) = P_tr17_2
   ELSE
     P_tr17_2 = tracer_index_table( PARAM_tr17_2 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_2 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_2 ) = 'tr17_2'
   tracer_desc_table( idomain, P_tr17_2 ) = 'tr17_2'
   tracer_units_table( idomain, P_tr17_2 ) = '-'
   tracer_streams_table( idomain, P_tr17_2 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_2 )%stream(2) = 2097152 ! 00200000 
   F_tr17_2 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_3 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_3 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_3 , idomain ) = P_tr17_3
   ELSE
     P_tr17_3 = tracer_index_table( PARAM_tr17_3 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_3 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_3 ) = 'tr17_3'
   tracer_desc_table( idomain, P_tr17_3 ) = 'tr17_3'
   tracer_units_table( idomain, P_tr17_3 ) = '-'
   tracer_streams_table( idomain, P_tr17_3 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_3 )%stream(2) = 2097152 ! 00200000 
   F_tr17_3 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_4 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_4 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_4 , idomain ) = P_tr17_4
   ELSE
     P_tr17_4 = tracer_index_table( PARAM_tr17_4 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_4 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_4 ) = 'tr17_4'
   tracer_desc_table( idomain, P_tr17_4 ) = 'tr17_4'
   tracer_units_table( idomain, P_tr17_4 ) = '-'
   tracer_streams_table( idomain, P_tr17_4 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_4 )%stream(2) = 2097152 ! 00200000 
   F_tr17_4 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_5 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_5 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_5 , idomain ) = P_tr17_5
   ELSE
     P_tr17_5 = tracer_index_table( PARAM_tr17_5 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_5 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_5 ) = 'tr17_5'
   tracer_desc_table( idomain, P_tr17_5 ) = 'tr17_5'
   tracer_units_table( idomain, P_tr17_5 ) = '-'
   tracer_streams_table( idomain, P_tr17_5 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_5 )%stream(2) = 2097152 ! 00200000 
   F_tr17_5 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_6 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_6 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_6 , idomain ) = P_tr17_6
   ELSE
     P_tr17_6 = tracer_index_table( PARAM_tr17_6 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_6 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_6 ) = 'tr17_6'
   tracer_desc_table( idomain, P_tr17_6 ) = 'tr17_6'
   tracer_units_table( idomain, P_tr17_6 ) = '-'
   tracer_streams_table( idomain, P_tr17_6 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_6 )%stream(2) = 2097152 ! 00200000 
   F_tr17_6 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_7 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_7 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_7 , idomain ) = P_tr17_7
   ELSE
     P_tr17_7 = tracer_index_table( PARAM_tr17_7 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_7 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_7 ) = 'tr17_7'
   tracer_desc_table( idomain, P_tr17_7 ) = 'tr17_7'
   tracer_units_table( idomain, P_tr17_7 ) = '-'
   tracer_streams_table( idomain, P_tr17_7 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_7 )%stream(2) = 2097152 ! 00200000 
   F_tr17_7 = .TRUE.
   IF ( tracer_index_table( PARAM_tr17_8 , idomain ) .lt. 1 ) THEN
     tracer_num_table(idomain) = tracer_num_table(idomain) + 1
     P_tr17_8 = tracer_num_table(idomain)
     tracer_index_table( PARAM_tr17_8 , idomain ) = P_tr17_8
   ELSE
     P_tr17_8 = tracer_index_table( PARAM_tr17_8 , idomain )
   END IF
   tracer_boundary_table( idomain, P_tr17_8 ) = .TRUE.
   tracer_dname_table( idomain, P_tr17_8 ) = 'tr17_8'
   tracer_desc_table( idomain, P_tr17_8 ) = 'tr17_8'
   tracer_units_table( idomain, P_tr17_8 ) = '-'
   tracer_streams_table( idomain, P_tr17_8 )%stream(1) = 33554433 ! 02000001 
   tracer_streams_table( idomain, P_tr17_8 )%stream(2) = 2097152 ! 00200000 
   F_tr17_8 = .TRUE.
  END IF
!ENDOFREGISTRYGENERATEDINCLUDE
