!STARTOFREGISTRYGENERATEDINCLUDE 'inc/dummy_decl.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
real                                     :: cfn
real                                     :: cfn1
integer                                  :: step_number
logical                                  :: this_is_an_ideal_run
logical                                  :: stepping_to_time
integer                                  :: last_step_updated
real                                     :: rdx
real                                     :: rdy
real                                     :: dts
real                                     :: dtseps
real                                     :: resm
real                                     :: zetatop
real                                     :: cf1
real                                     :: cf2
real                                     :: cf3
integer                                  :: itimestep
real                                     :: xtime
real                                     :: julian
integer                                  :: lbc_fid
logical                                  :: tiled
logical                                  :: patched
logical                                  :: press_adj
real                                     :: dtbc
real                                     :: u_frame
real                                     :: v_frame
real                                     :: p_top
real                                     :: t00
real                                     :: p00
real                                     :: tlp
real                                     :: tiso
real                                     :: tlp_strat
real                                     :: p_strat
real                                     :: max_msftx
real                                     :: max_msfty
integer                                  :: stepra
logical                                  :: warm_rain
logical                                  :: adv_moist_cond
integer                                  :: save_topo_from_real
integer                                  :: auxhist1_oid
integer                                  :: auxhist2_oid
integer                                  :: auxhist3_oid
integer                                  :: auxhist4_oid
integer                                  :: auxhist5_oid
integer                                  :: auxhist6_oid
integer                                  :: auxhist7_oid
integer                                  :: auxhist8_oid
integer                                  :: auxhist9_oid
integer                                  :: auxhist10_oid
integer                                  :: auxhist11_oid
integer                                  :: auxhist12_oid
integer                                  :: auxhist13_oid
integer                                  :: auxhist14_oid
integer                                  :: auxhist15_oid
integer                                  :: auxhist16_oid
integer                                  :: auxhist17_oid
integer                                  :: auxhist18_oid
integer                                  :: auxhist19_oid
integer                                  :: auxhist20_oid
integer                                  :: auxhist21_oid
integer                                  :: auxhist22_oid
integer                                  :: auxhist23_oid
integer                                  :: auxhist24_oid
integer                                  :: auxinput1_oid
integer                                  :: auxinput2_oid
integer                                  :: auxinput3_oid
integer                                  :: auxinput4_oid
integer                                  :: auxinput5_oid
integer                                  :: auxinput6_oid
integer                                  :: auxinput7_oid
integer                                  :: auxinput8_oid
integer                                  :: auxinput9_oid
integer                                  :: auxinput10_oid
integer                                  :: auxinput11_oid
integer                                  :: auxinput12_oid
integer                                  :: auxinput13_oid
integer                                  :: auxinput14_oid
integer                                  :: auxinput15_oid
integer                                  :: auxinput16_oid
integer                                  :: auxinput17_oid
integer                                  :: auxinput18_oid
integer                                  :: auxinput19_oid
integer                                  :: auxinput20_oid
integer                                  :: auxinput21_oid
integer                                  :: auxinput22_oid
integer                                  :: auxinput23_oid
integer                                  :: auxinput24_oid
integer                                  :: oid
TYPE(vp_type)                               :: vv
TYPE(vp_type)                               :: vv6
TYPE(vp_type)                               :: vp
TYPE(vp_type)                               :: vp6
TYPE(ep_type)                               :: ep
TYPE(xb_type)                               :: xb
TYPE(x_type)                               :: xa
TYPE(x_type)                               :: x6a
TYPE(xpose_type)                               :: xp
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: lu_index
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: lu_mask
real      ,DIMENSION(grid%sm33:grid%em33)           :: znu
real      ,DIMENSION(grid%sm33:grid%em33)           :: znw
real      ,DIMENSION(1:grid%num_soil_layers)           :: zs
real      ,DIMENSION(1:grid%num_soil_layers)           :: dzs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: u_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: u_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ru
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ru_m
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ru_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: v_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: v_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rv_m
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rv_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: phb
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: php
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_init
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mu_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mu_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mub
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mu0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mudf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: muu
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: muv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mut
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: muts
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: al
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: alt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: alb
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: zx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: zy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rdz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rdzw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pb
real      ,DIMENSION(grid%sm33:grid%em33)           :: fnm
real      ,DIMENSION(grid%sm33:grid%em33)           :: fnp
real      ,DIMENSION(grid%sm33:grid%em33)           :: rdnw
real      ,DIMENSION(grid%sm33:grid%em33)           :: rdn
real      ,DIMENSION(grid%sm33:grid%em33)           :: dnw
real      ,DIMENSION(grid%sm33:grid%em33)           :: dn
real      ,DIMENSION(grid%sm33:grid%em33)           :: t_base
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: z_at_w
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p_hyd
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p_hyd_w
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: q2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: t2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: th2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: psfc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: qsfc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: u10
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: v10
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: imask_nostag
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: imask_xstag
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: imask_ystag
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: imask_xystag
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(1:grid%spec_bdy_width)           :: fcx
real      ,DIMENSION(1:grid%spec_bdy_width)           :: gcx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,1:grid%num_soil_layers)           :: smois
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,1:grid%num_soil_layers)           :: sh2o
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,1:grid%num_soil_layers)           :: smcrel
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xice
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: ivgtyp
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: isltyp
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: vegfra
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: snow
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: snowh
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: htop
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: hbot
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,1:grid%num_soil_layers)           :: tslb
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msft
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfu
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msftx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfty
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfux
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfuy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfvx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfvx_inv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: msfvy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: f
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: e
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: sina
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: cosa
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: ht
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: ht_fine
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: ht_int
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: ht_input
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: tsk
real      ,DIMENSION(grid%sm33:grid%em33)           :: u_base
real      ,DIMENSION(grid%sm33:grid%em33)           :: v_base
real      ,DIMENSION(grid%sm33:grid%em33)           :: qv_base
real      ,DIMENSION(grid%sm33:grid%em33)           :: z_base
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: rainc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: rainnc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: raincv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: rainncv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xlat
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xlong
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xlat_u
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xlong_u
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xlat_v
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xlong_v
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: albedo
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: clat
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: albbck
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: embck
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: emiss
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: snotime
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: noahres
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: cldefi
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: tmn
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: xland
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: snowc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: znt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: ust
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: mol
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: pblh
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: hfx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: qfx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: regime
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_u_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_u_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_u_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_u_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_v_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_v_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_v_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_v_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_w_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_w_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_w_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_w_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_ph_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_ph_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_ph_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_ph_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_t_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_t_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_t_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_t_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: a_mu_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: a_mu_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: g_mu_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: g_mu_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_p
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_p
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_z
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_z
real      ,DIMENSION(grid%sm33:grid%em33)           :: kma_a
real      ,DIMENSION(grid%sm33:grid%em33)           :: kma_b
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_a_moist)           :: a_moist
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_g_moist)           :: g_moist
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: g_h_diabatic
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_h_diabatic
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: g_rainc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: a_rainc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: g_rainnc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: a_rainnc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: g_raincv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: a_raincv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: g_rainncv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: a_rainncv
real      ,DIMENSION(grid%sm31:grid%em31)           :: dummy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_a_scalar)           :: a_scalar
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_g_scalar)           :: g_scalar
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: landmask
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32)           :: sst
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%sm32:grid%em32,num_tracer)           :: tracer
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_bye
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btxs
real      ,DIMENSION(grid%sm32:grid%em32,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%spec_bdy_width,num_tracer)           :: tracer_btye
#undef COPY_IN
#undef COPY_OUT
!ENDOFREGISTRYGENERATEDINCLUDE
