































MODULE module_state_description
  
  INTEGER, PARAMETER :: slabscheme = 1
  INTEGER, PARAMETER :: lsmscheme = 2
  INTEGER, PARAMETER :: ruclsmscheme = 3
  INTEGER, PARAMETER :: noahmpscheme = 4
  INTEGER, PARAMETER :: dfi_setup = 0
  INTEGER, PARAMETER :: dfi_bck = 1
  INTEGER, PARAMETER :: dfi_fwd = 2
  INTEGER, PARAMETER :: dfi_fst = 3
  INTEGER, PARAMETER :: dfi_startfwd = 4
  INTEGER, PARAMETER :: dfi_startbck = 5
  INTEGER, PARAMETER :: dfi_nodfi = 0
  INTEGER, PARAMETER :: dfi_dfl = 1
  INTEGER, PARAMETER :: dfi_ddfi = 2
  INTEGER, PARAMETER :: dfi_tdfi = 3
  INTEGER, PARAMETER :: io_intio = 1
  INTEGER, PARAMETER :: io_netcdf = 2
  INTEGER, PARAMETER :: io_hdf = 3
  INTEGER, PARAMETER :: io_phdf5 = 4
  INTEGER, PARAMETER :: io_grib1 = 5
  INTEGER, PARAMETER :: io_mcel = 6
  INTEGER, PARAMETER :: io_esmf = 7
  INTEGER, PARAMETER :: io_yyy = 8
  INTEGER, PARAMETER :: io_zzz = 9
  INTEGER, PARAMETER :: io_grib2 = 10
  INTEGER, PARAMETER :: io_pnetcdf = 11
  INTEGER, PARAMETER :: io_pio = 12
  INTEGER, PARAMETER :: no_wrfhydro = 0
  INTEGER, PARAMETER :: wrfhydro = 1
  INTEGER, PARAMETER :: dyn_nodyn = 0
  INTEGER, PARAMETER :: dyn_em = 2
  INTEGER, PARAMETER :: dyn_em_sn = 102
  INTEGER, PARAMETER :: dyn_em_tl = 202
  INTEGER, PARAMETER :: dyn_em_ad = 302
  INTEGER, PARAMETER :: dyn_em_tst = 402
  INTEGER, PARAMETER :: dyn_em_var = 502
  INTEGER, PARAMETER :: passiveqv = 0
  INTEGER, PARAMETER :: kesslerscheme = 1
  INTEGER, PARAMETER :: linscheme = 2
  INTEGER, PARAMETER :: wsm3scheme = 3
  INTEGER, PARAMETER :: wsm5scheme = 4
  INTEGER, PARAMETER :: fer_mp_hires = 5
  INTEGER, PARAMETER :: fer_mp_hires_advect = 15
  INTEGER, PARAMETER :: wsm6scheme = 6
  INTEGER, PARAMETER :: gsfcgcescheme = 7
  INTEGER, PARAMETER :: thompson = 8
  INTEGER, PARAMETER :: milbrandt2mom = 9
  INTEGER, PARAMETER :: morr_two_moment = 10
  INTEGER, PARAMETER :: cammgmpscheme = 11
  INTEGER, PARAMETER :: sbu_ylinscheme = 13
  INTEGER, PARAMETER :: wdm5scheme = 14
  INTEGER, PARAMETER :: wdm6scheme = 16
  INTEGER, PARAMETER :: nssl_2mom = 17
  INTEGER, PARAMETER :: nssl_2momccn = 18
  INTEGER, PARAMETER :: nssl_1mom = 19
  INTEGER, PARAMETER :: nssl_1momlfo = 21
  INTEGER, PARAMETER :: nssl_2momg = 22
  INTEGER, PARAMETER :: thompsonaero = 28
  INTEGER, PARAMETER :: etampnew = 95
  INTEGER, PARAMETER :: lscondscheme = 98
  INTEGER, PARAMETER :: mkesslerscheme = 99
  INTEGER, PARAMETER :: surfdragscheme = 98
  INTEGER, PARAMETER :: ducuscheme = 98
  INTEGER, PARAMETER :: tracer_test1 = 2
  
  INTEGER, PARAMETER :: PARAM_qv = 1
  INTEGER            ::     P_qv = 1
  LOGICAL            ::     F_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_qc = 2
  INTEGER            ::     P_qc = 1
  LOGICAL            ::     F_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_qr = 3
  INTEGER            ::     P_qr = 1
  LOGICAL            ::     F_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_qi = 4
  INTEGER            ::     P_qi = 1
  LOGICAL            ::     F_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_qs = 5
  INTEGER            ::     P_qs = 1
  LOGICAL            ::     F_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_qg = 6
  INTEGER            ::     P_qg = 1
  LOGICAL            ::     F_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_qh = 7
  INTEGER            ::     P_qh = 1
  LOGICAL            ::     F_qh = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_moist = 8
  INTEGER            ::       NUM_moist = 1
  INTEGER, PARAMETER :: PARAM_qndrop = 1
  INTEGER            ::     P_qndrop = 1
  LOGICAL            ::     F_qndrop = .FALSE.
  INTEGER, PARAMETER :: PARAM_qni = 2
  INTEGER            ::     P_qni = 1
  LOGICAL            ::     F_qni = .FALSE.
  INTEGER, PARAMETER :: PARAM_qt = 3
  INTEGER            ::     P_qt = 1
  LOGICAL            ::     F_qt = .FALSE.
  INTEGER, PARAMETER :: PARAM_qns = 4
  INTEGER            ::     P_qns = 1
  LOGICAL            ::     F_qns = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnr = 5
  INTEGER            ::     P_qnr = 1
  LOGICAL            ::     F_qnr = .FALSE.
  INTEGER, PARAMETER :: PARAM_qng = 6
  INTEGER            ::     P_qng = 1
  LOGICAL            ::     F_qng = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnh = 7
  INTEGER            ::     P_qnh = 1
  LOGICAL            ::     F_qnh = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnn = 8
  INTEGER            ::     P_qnn = 1
  LOGICAL            ::     F_qnn = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnc = 9
  INTEGER            ::     P_qnc = 1
  LOGICAL            ::     F_qnc = .FALSE.
  INTEGER, PARAMETER :: PARAM_qvolg = 10
  INTEGER            ::     P_qvolg = 1
  LOGICAL            ::     F_qvolg = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_scalar = 11
  INTEGER            ::       NUM_scalar = 1
  INTEGER, PARAMETER :: PARAM_a_qv = 1
  INTEGER            ::     P_a_qv = 1
  LOGICAL            ::     F_a_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_a_qc = 2
  INTEGER            ::     P_a_qc = 1
  LOGICAL            ::     F_a_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_a_qr = 3
  INTEGER            ::     P_a_qr = 1
  LOGICAL            ::     F_a_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_a_qi = 4
  INTEGER            ::     P_a_qi = 1
  LOGICAL            ::     F_a_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_a_qs = 5
  INTEGER            ::     P_a_qs = 1
  LOGICAL            ::     F_a_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_a_qg = 6
  INTEGER            ::     P_a_qg = 1
  LOGICAL            ::     F_a_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_a_qh = 7
  INTEGER            ::     P_a_qh = 1
  LOGICAL            ::     F_a_qh = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_a_moist = 8
  INTEGER            ::       NUM_a_moist = 1
  INTEGER, PARAMETER :: PARAM_g_qv = 1
  INTEGER            ::     P_g_qv = 1
  LOGICAL            ::     F_g_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_g_qc = 2
  INTEGER            ::     P_g_qc = 1
  LOGICAL            ::     F_g_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_g_qr = 3
  INTEGER            ::     P_g_qr = 1
  LOGICAL            ::     F_g_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_g_qi = 4
  INTEGER            ::     P_g_qi = 1
  LOGICAL            ::     F_g_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_g_qs = 5
  INTEGER            ::     P_g_qs = 1
  LOGICAL            ::     F_g_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_g_qg = 6
  INTEGER            ::     P_g_qg = 1
  LOGICAL            ::     F_g_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_g_qh = 7
  INTEGER            ::     P_g_qh = 1
  LOGICAL            ::     F_g_qh = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_g_moist = 8
  INTEGER            ::       NUM_g_moist = 1
  INTEGER, PARAMETER :: PARAM_NUM_a_scalar = 1
  INTEGER            ::       NUM_a_scalar = 1
  INTEGER, PARAMETER :: PARAM_NUM_g_scalar = 1
  INTEGER            ::       NUM_g_scalar = 1
  INTEGER, PARAMETER :: PARAM_NUM_chem = 1
  INTEGER            ::       NUM_chem = 1
  INTEGER, PARAMETER :: PARAM_tr17_1 = 1
  INTEGER            ::     P_tr17_1 = 1
  LOGICAL            ::     F_tr17_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_2 = 2
  INTEGER            ::     P_tr17_2 = 1
  LOGICAL            ::     F_tr17_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_3 = 3
  INTEGER            ::     P_tr17_3 = 1
  LOGICAL            ::     F_tr17_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_4 = 4
  INTEGER            ::     P_tr17_4 = 1
  LOGICAL            ::     F_tr17_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_5 = 5
  INTEGER            ::     P_tr17_5 = 1
  LOGICAL            ::     F_tr17_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_6 = 6
  INTEGER            ::     P_tr17_6 = 1
  LOGICAL            ::     F_tr17_6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_7 = 7
  INTEGER            ::     P_tr17_7 = 1
  LOGICAL            ::     F_tr17_7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_8 = 8
  INTEGER            ::     P_tr17_8 = 1
  LOGICAL            ::     F_tr17_8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_tracer = 9
  INTEGER            ::       NUM_tracer = 1
  INTEGER, PARAMETER :: P_XSB                          = 1
  INTEGER, PARAMETER :: P_XEB                          = 2
  INTEGER, PARAMETER :: P_YSB                          = 3
  INTEGER, PARAMETER :: P_YEB                          = 4
  INTEGER, PARAMETER :: NUM_TIME_LEVELS = 2
  INTEGER , PARAMETER :: PARAM_FIRST_SCALAR = 2
CONTAINS
SUBROUTINE init_module_state_description
END SUBROUTINE init_module_state_description
END MODULE module_state_description

