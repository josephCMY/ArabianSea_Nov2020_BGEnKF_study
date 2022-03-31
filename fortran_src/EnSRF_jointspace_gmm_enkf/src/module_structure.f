module constants

!-----------------------------------------------------------------------------
!  PURPOSE: Common reference point for constants.
!
!  METHOD:  Straightforward definitions.
!------------------------------------------------------------------------------

   IMPLICIT NONE

!------------------------------------------------------------------------------
!  [1.0] Physical parameter constants (all NIST standard values):
!------------------------------------------------------------------------------
   real, parameter :: g=9.81,                       &! gravitational constant
                      to=300.,                      &! base pot. temp.
                      cp=1004.5,                    &! specific heat
                      rd=287.,                      &! dry air gas constant
                      rv=461.6,                     &! moist air gas constant
                      lhv=2.50*10**6,               &! latent heat of vap.
                      R_earth=6374.,                &! radius of the earth
                      Pr=100000.,                   &! reference pressure
                      Po=101325.,                   &! altimeter reference P
                      Talt=288.15,                  &! attimeter reference T
                      pbrh=611.,                    &! vapor pressure at 0 C (Pa)
                      T_frez=273.15,                &! water freezing point
                      RvRd=Rv/Rd,                   &! Rd/Rv
                      kappa=2./7.,                  &! kappa for pot. temp
                      pid=3.1415/180.,              &! radians to degrees
                      om_earth=(2.*3.1415)/(24.*3600.),&! earth rot
                      lap_std_atmos=0.0065,         &! standard atmosphere
                      th_denom=0.03727594,          &! Pr**(-kappa)
                      missing=-888888.,             &! missing value
                      es_alpha = 611.2,             &! saturation vapor pressure
                      es_beta = 17.67,              &!
                      es_gamma = 243.5
  real, parameter ::  pi = 3.1415927,               &!
                      deg_per_rad = 180./pi,        &!
                      rad_per_deg = pi / 180.,      &!
                      earth_radius_m = 6371200.
  integer,parameter :: PROJ_LATLON = 0,             &!
                       PROJ_MERC   = 1,             &!
                       PROJ_LC     = 3,             &!
                       PROJ_PS     = 5


end module constants

!==============================================================================

module namelist_define

   implicit none

!------------------------------------------------------------------------------
!  [1.0] Namelist parameters:
!------------------------------------------------------------------------------
!-- enkf_parameter
   integer       :: numbers_en          ! ensemble size
   character(len=10) :: expername       ! experiment name, such as hurricane, clear-air
   character(len=10), dimension(30) :: enkfvar   ! include the variables which will 
                                                 ! be updated and be used to 
                                                 ! calculate the XB
   character(len=10), dimension(20) :: updatevar ! updated variables
   integer       :: update_is, update_ie, update_js, update_je, update_ks, update_ke  ! domain to be updated
   real          :: inflate             !
   real          :: mixing              !
   logical       :: use_vortex_RTPP     !.true.: use different mixing coeff for 
                                        !        inner-core and env
   real          :: mixing_core
   integer       :: print_detail        ! used to debug the code
   logical       :: use_sfc_theta, use_sfc_td ! Using theta and Td instead of T and 
                                              ! Qv for surface obs when available
   logical       :: use_nonlinear_enkf  ! .true. : Assimilate IR obs, with CRTM 
                                        !          executed in a batchwise fashion.
   logical       :: batchwise_order     ! .true.: Arrange obs order to be
                                        !         consistent with nonlinear enkf
   logical       :: use_gmm_enkf        ! .true.: Assimilate IR obs using GMM EnKF
   logical       :: use_jointspace_gmm_enkf  ! .true.: Assimialte IR obs using joint-space GMM-EnKF
   integer       :: max_kernel_num      ! Number of Gaussian kernels for GMM EnKF
   integer       :: min_expanding_kernel_prior_size ! Minimum number of members within expanding cluster.
                                                    ! Otherwise, default to single-kernel EnKF.

!-- parallel
   logical       :: manual_parallel, random_order
   integer       :: nicpu, njcpu, nmcpu




!-- use_osse
   logical       :: use_ideal_obs       ! .true. : use ideal observation, no real obs. position. Just for use_sounding or use_groundbase_radar
   integer       :: gridobs_is, gridobs_ie, gridobs_js, gridobs_je, gridobs_ks, gridobs_ke  ! domain to pick up data as gridobservation, for use_ideal_obs
   integer       :: gridobs_int_x, gridobs_int_k     ! horizontal and vertical interval grid for picking up gridobs, for use_ideal_obs
   logical       :: use_simulated       ! .true. : use simulated radar data

!-- use_hurricane_position_intensity
   logical       :: use_hurricane_PI    ! .true. : assimilated hurricane position and intensity
   integer       :: hroi_hurricane_PI   ! horizontal radius of influence for hurricane PI
   integer       :: vroi_hurricane_PI   ! vertical radius of influence for hurricane PI

!-- use_surface_obs
   logical       :: use_surface         ! .true. : assimilated SURFACE LAND (SYNOPTIC, METAR) REPORTS
   integer       :: datathin_surface    ! 0=all data, 2=1/2 data, 10=1/10 data
                                        ! 2: get the 1st, 3rd, 5th ... data
                                        !-2: get the 2nd, 4th, 6th ... data
   integer       :: hroi_surface        ! horizontal radius of influence for surface
   integer       :: vroi_surface        ! vertical radius of influence for surface

!-- use_sounding_obs
   logical       :: use_sounding        ! .true. : assimilated UPPER-AIR (RAOB, PIBAL, RECCO, DROPS) REPORTS
   integer       :: datathin_sounding   ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_sounding       ! horizontal radius of influence for sounding
   integer       :: vroi_sounding       ! vertical radius of influence for sounding

!-- use_profiler_obs
   logical       :: use_profiler        ! .true. : assimilated WIND PROFILER REPORTS
   integer       :: datathin_profiler   ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_profiler       ! horizontal radius of influence for profiler
   integer       :: vroi_profiler       ! vertical radius of influence for profiler

!-- use_aircft_obs
   logical       :: use_aircft          ! .true. : assimilated AIREP/PIREP, AMDAR (ASDAR/ACARS), E-ADAS (AMDAR BUFR) AIRCRAFT REPORTS
   integer       :: datathin_aircft     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_aircft         ! horizontal radius of influence for aircft
   integer       :: vroi_aircft         ! vertical radius of influence for aircft

!-- use_metar_obs
   logical       :: use_metar           ! .true. : assimilated auto-meteo-station report
   integer       :: datathin_metar      ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_metar          ! horizontal radius of influence for sfcshp
   integer       :: vroi_metar          ! vertical radius of influence for sfcshp

!-- use_sfcshp_obs
   logical       :: use_sfcshp          ! .true. : assimilated SURFACE MARINE (SHIP, BUOY, C-MAN PLATFORM) REPORTS
   integer       :: datathin_sfcshp     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_sfcshp         ! horizontal radius of influence for sfcshp
   integer       :: vroi_sfcshp         ! vertical radius of influence for sfcshp

!-- use_spssmi_obs
   logical       :: use_spssmi          ! .true. : assimilated DMSP SSM/I RETRIEVAL PRODUCTS (REPROCESSED WIND SPEED, TPW)
   integer       :: datathin_spssmi     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_spssmi         ! horizontal radius of influence for spssmi
   integer       :: vroi_spssmi         ! vertical radius of influence for spssmi

!-- use_atovs_obs
   logical       :: use_atovs           ! .true. : assimilated WIND PROFILER REPORTS
   integer       :: datathin_atovs      ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_atovs          ! horizontal radius of influence for atovs
   integer       :: vroi_atovs          ! vertical radius of influence for atovs

!-- use_satwnd_obs
   logical       :: use_satwnd          ! .true. : assimilated SATELLITE-DERIVED WIND REPORTS
   integer       :: datathin_satwnd     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_satwnd         ! horizontal radius of influence for satwnd
   integer       :: vroi_satwnd         ! vertical radius of influence for satwnd

!-- use_gpspw_obs
   logical       :: use_gpspw          ! .true. : assimilated SATELLITE-DERIVED WIND REPORTS
   integer       :: datathin_gpspw     ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_gpspw         ! horizontal radius of influence for satwnd
   integer       :: vroi_gpspw         ! vertical radius of influence for satwnd

!-- use_groundbase_radar_obs
   integer       :: radar_number        ! radar number
   logical       :: use_radar_rf        ! .true. : assimilated radar_rf
   logical       :: use_radar_rv        ! .true. : assimilated radar_rv
   integer       :: datathin_radar      ! 0=all data, 2=1/2 data, 10=1/10 data
   integer       :: hroi_radar          ! horizontal radius of influence for radar
   integer       :: vroi_radar          ! vertical radius of influence for radar

!-- use_airborne_radar
   logical       :: use_airborne_rf     ! .true. : assimilated airborne_rf
   logical       :: use_airborne_rv     ! .true. : assimilated airborne_rv
   integer       :: datathin_airborne   ! 0=all data, 2=1/2 data, 10=1/10 airborne
   integer       :: hroi_airborne       ! horizontal radius of influence for airborne
   integer       :: vroi_airborne       ! vertical radius of influence for airborne

!-- use_radiance
   logical       :: use_radiance        ! .true. : assimilated radiance
   integer       :: datathin_radiance   ! 0=all data, 2=1/2 data, 10=1/10 radiance
   integer       :: hroi_radiance       ! horizontal radius of influence for radiance
   integer       :: vroi_radiance       ! vertical radius of influence for radiance if use_vroi_radiance_halfsfc = .false.
   logical       :: use_vroi_radiance_halfsfc ! .true. : VROI is set to (obs-level)*2, .false.:vroi_radiance is used
   logical       :: use_aoei            ! .true. : apply Adaptive Observation Error Inflation
   logical       :: use_elf             ! .true. : apply Empirical Localization Function
   logical       :: use_abei          ! .true. : apply Empirical BT-based inflation

!-- use_microwave
   logical       :: use_microwave        ! .true. : assimilated microwave
   integer       :: datathin_microwave   ! 0=all data, 2=1/2 data, 10=1/10 microwave
   integer       :: hroi_microwave       ! horizontal radius of influence for microwave
   integer       :: vroi_microwave       ! vertical radius of influence for microwave
   logical       :: aoei_microwave       ! .true. : apply Adaptive Observation Error Inflation
   logical       :: use_slant_path       ! .true. : use slant path for CRTM calculation

!-- Namelist contents :

   namelist /enkf_parameter / numbers_en, expername, enkfvar, updatevar,&
                                 update_is, update_ie, update_js, update_je, update_ks, update_ke,   &
                                 inflate, mixing, use_vortex_RTPP, mixing_core,  &
                                 random_order, print_detail, use_sfc_theta, use_sfc_td, &
                                 batchwise_order, use_nonlinear_enkf, &
                                 use_gmm_enkf, use_jointspace_gmm_enkf, max_kernel_num, &
                                 min_expanding_kernel_prior_size 
   namelist /parallel       / manual_parallel, nicpu, njcpu, nmcpu
   namelist /osse           / use_ideal_obs, gridobs_is, gridobs_ie, gridobs_js, gridobs_je,      &
                                 gridobs_ks, gridobs_ke, gridobs_int_x, gridobs_int_k, use_simulated
   namelist /hurricane_PI   / use_hurricane_PI, hroi_hurricane_PI, vroi_hurricane_PI
   namelist /surface_obs    / use_surface, datathin_surface, hroi_surface, vroi_surface
   namelist /sounding_obs   / use_sounding, datathin_sounding, hroi_sounding, vroi_sounding
   namelist /profiler_obs   / use_profiler, datathin_profiler, hroi_profiler, vroi_profiler
   namelist /aircft_obs     / use_aircft, datathin_aircft, hroi_aircft, vroi_aircft
   namelist /metar_obs      / use_metar , datathin_metar , hroi_metar , vroi_metar
   namelist /sfcshp_obs     / use_sfcshp, datathin_sfcshp, hroi_sfcshp, vroi_sfcshp
   namelist /spssmi_obs     / use_spssmi, datathin_spssmi, hroi_spssmi, vroi_spssmi
   namelist /atovs_obs      / use_atovs, datathin_atovs, hroi_atovs, vroi_atovs
   namelist /satwnd_obs     / use_satwnd, datathin_satwnd, hroi_satwnd, vroi_satwnd
   namelist /gpspw_obs      / use_gpspw, datathin_gpspw, hroi_gpspw, vroi_gpspw
   namelist /radar_obs      / radar_number, use_radar_rf, use_radar_rv, datathin_radar, hroi_radar, vroi_radar
   namelist /airborne_radar / use_airborne_rf, use_airborne_rv, datathin_airborne, hroi_airborne, vroi_airborne
   namelist /radiance/use_radiance,datathin_radiance,hroi_radiance,vroi_radiance,use_vroi_radiance_halfsfc,use_aoei,use_elf,use_abei
   namelist /microwave      / use_microwave, datathin_microwave, hroi_microwave,vroi_microwave, aoei_microwave, use_slant_path

end module namelist_define

!==============================================================================

module mapinfo_define

 IMPLICIT NONE

 ! Define data structures to define various projections

      TYPE proj_info

      INTEGER    :: code     ! Integer code for projection type
      REAL       :: lat1    ! SW latitude (1,1) in degrees (-90->90N)
      REAL       :: lon1    ! SW longitude (1,1) in degrees (-180->180E)
      REAL       :: dx       ! Grid spacing in meters at truelats, used
                             ! only for ps, lc, and merc projections
      REAL       :: dlat     ! Lat increment for lat/lon grids
      REAL       :: dlon     ! Lon increment for lat/lon grids
      REAL       :: stdlon   ! Longitude parallel to y-axis (-180->180E)
      REAL       :: truelat1 ! First true latitude (all projections)
      REAL       :: truelat2 ! Second true lat (LC only)
      REAL       :: hemi     ! 1 for NH, -1 for SH
      REAL       :: cone     ! Cone factor for LC projections
      REAL       :: polei    ! Computed i-location of pole point
      REAL       :: polej    ! Computed j-location of pole point
      REAL       :: rsw      ! Computed radius to SW corner
      REAL       :: rebydx   ! Earth radius divided by dx
      LOGICAL    :: init     ! Flag to indicate if this struct is
                                 ! ready for use
      INTEGER    :: nx
      INTEGER    :: ny
      END TYPE proj_info
end module mapinfo_define




!=============================================================================
module wrf_field_indices_define

!-----------------------------------------------------------------------------
!  BACKGROUND:
!  Variable x in enkf.mpi has dimensions (ni, nj, nk, nv, nm).
!  The nv dimension corresponds to the wrf fields listed in enkfvar (defined in
!  namelist.enkf).
!  We will allocate the indices corresponding to each field here.
!------------------------------------------------------------------------------

   IMPLICIT NONE

   ! Hydrometeor and qvapor fields (assuming 6-class)
   integer :: ind_qvapor, ind_qcloud, ind_qrain
   integer :: ind_qice, ind_qsnow, ind_qgraup

   ! Thermo fields
   integer :: ind_p, ind_pb, ind_tsk, ind_t
   integer :: ind_t2,ind_th2, ind_q2, ind_psfc

   ! Dynamic fields
   integer :: ind_ph, ind_phb, ind_hgt
   integer :: ind_u, ind_v, ind_w
   integer :: ind_u10, ind_v10, ind_mu, ind_mub


end module wrf_field_indices_define


!=============================================================================
module common_variables

!-----------------------------------------------------------------------------
!  BACKGROUND:
!  There's a boatload of commonly used variables in each subroutine.
!  To reduce the number of subroutine arguments and variable definitions, I am
!  putting a bunch of commonly used variables in this mod.
!  That way, there's no need to use arguments to pass info around.
!  Might also save a bunch of memory.
!------------------------------------------------------------------------------

   IMPLICIT NONE

   ! Constants
   integer :: ni, nj, nk, nm, nv, ix, jx, kx, kxs
   character(len=80) :: times
   real :: p_top

   ! Constant arrays (lat, lon, landuse index, landsea mask)
   real, allocatable, dimension(:,:  ):: xlat, xlong, lu_index, xland
   real, allocatable, dimension(:    ):: znw, znu

   ! Variable arrays
   real,    allocatable, dimension(:,:,:,:,:) :: x
   real,    allocatable, dimension(:,:,:,:  ) :: xm
   real,    allocatable, dimension(:,:      ) :: ya, yf, yasend
   real,    allocatable, dimension(:        ) :: yfm_radiance, yam_radiance
   integer, allocatable, dimension(:        ) :: ind

   ! Size of variables
   integer, dimension( 50, 3 ) :: var_dims

end module common_variables




! ===========================================================================
module boundary_variables

   ! Variables used to extend the boundary of the slabs.
   ! Very important for compute_yf


   implicit none
   ! slab boundary extending arrays
   real, allocatable, dimension(:,:,:,:,:)::  x_we_bdy,  x_sn_bdy,  x_crnrs
   real, allocatable, dimension(:,:,:,:  ):: xm_we_bdy, xm_sn_bdy, xm_crnrs
end module boundary_variables






!==============================================================================
module obs_define

!                    _
!                    | radar_stn   % numObs, levels, lat, lon, elv, i_radar, j_radar
! obs % radar(200) % | radar_point % ii, jj, hgt
!                    | rf
!                    | rv
!                    -

   implicit none

   logical                     :: flag_level_assignment

   ! Variables relating to determining unique obs names
   character( len=10 ), dimension( 200) :: uniq_obs_types
   integer :: n_uniq_obs_types


   TYPE Radar_stn_type
        INTEGER                :: idn           ! station id number
        CHARACTER (LEN = 4)    :: idc           ! station id character
        CHARACTER (LEN = 50)   :: name          ! Station name
        CHARACTER (LEN = 19)   :: date_char     ! CCYY-MM-DD_HH:MM:SS date
        INTEGER                :: numObs        ! number of Obs
        INTEGER                :: levels        ! number of levels
        REAL                   :: lat           ! Latitude in degree
        REAL                   :: lon           ! Longitude in degree
        REAL                   :: elv           ! Elevation in m
        real, dimension(20)    :: elv_ang       ! elevation angles
        real                   :: i_radar       ! radar base position according to wrf model domain
        real                   :: j_radar       ! i_radar: west-east; j_radar: south-north
        real                   :: min_dis       ! distance between 1st data and radar
        real                   :: max_dis       ! valid data distance
        real                   :: d_dis         ! resolution in radial
        real                   :: d_azim        ! resolution in azimuth
        real                   :: err_rf        ! reflective observation error
        real                   :: err_rv        ! radial velocity observation error
   END TYPE Radar_stn_type

   type radar_data_point_type
        real                    :: ii          ! each valid radar data point according to wrf domain
        real                    :: jj          ! ii: west-east; jj: south-north
        real                    :: kk          ! kk: half eta level
        real                    :: hgt         ! each valid radar data height above sea level
        real                    :: rdis        ! radial distance   ! just output for plot
        real                    :: azim        ! azimuth           ! just output for plot
        real                    :: elev        ! elevation angles  ! just output for plot
   end type radar_data_point_type

   type radar_data_type
        type ( Radar_stn_type )        :: radar_stn
        type ( radar_data_point_type ), allocatable,dimension(:) :: radar_point
        real, allocatable,dimension(:)             :: rf    ! radar reflective
        real, allocatable,dimension(:)             :: rv    ! radar radial velocity
   end type radar_data_type

   type airborne_data_type
        integer                                    :: num
        real, allocatable,dimension(:)             :: lat, lon, height, radar_ii, radar_jj   !for radar location
        real, allocatable,dimension(:)             :: azim, elev, range, rf, rv, ii, jj, hh
   end type airborne_data_type

   type gts_data_type
        integer                                    :: num
        character(len=16),allocatable,dimension(:) :: platform
        character(len=19),allocatable,dimension(:) :: date
        real, allocatable,dimension(:)             :: latitude, longitude, elevation
        real, allocatable,dimension(:,:)           :: slp                              !slp(:,3), 3: data, qc, error
        real, allocatable,dimension(:,:)           :: pw                               !pw(:,3), 3: data, qc, error
        integer, allocatable,dimension(:)          :: levels
        real, allocatable,dimension(:,:,:)         :: pres, spd, wd, height, t, td, rh !(:,:,3), 3: data, qc, error
   end type gts_data_type

   type radiance_data_type
        integer                                    :: num
        character(len=16),allocatable,dimension(:) :: platform
        real, allocatable,dimension(:)             :: lat, lon,ii, jj, kk, tb, err
        integer, allocatable,dimension(:)          :: ch, hroi, hroi_d
   end type radiance_data_type

   type microwave_data_type
        integer                                    :: num
        character(len=16),allocatable,dimension(:) :: platform
        real, allocatable,dimension(:)             :: lat, lon, ii, jj, tb, err
        integer, allocatable,dimension(:)          :: ch, hroi, hroi_d
        real, allocatable,dimension(:)             :: efov_aScan, efov_cScan, scan_angle, zenith_angle, azimuth_angle
        real, allocatable,dimension(:)             :: sat_lat, sat_lon, sat_alt
   end type microwave_data_type

   type raw_type
        integer                                  :: radar_stn_num
        type ( Radar_data_type ), allocatable,dimension( : )  :: radar
        type ( airborne_data_type  )             :: airborne
        type ( gts_data_type      )              :: gts
        type ( radiance_data_type      )         :: radiance
        type ( microwave_data_type     )         :: microwave
   end type raw_type

   type obs_type
       integer                                   :: num          !! observation number
       real, allocatable, dimension(:)           :: dat              !! observation data
       character(len=10), allocatable, dimension(:) :: type  !! observation type
       real, allocatable, dimension(:)           :: err         !! observation error
       real, allocatable, dimension(:,:)         :: position    !! (ob_num,4)
                                                                    !! observation data position in wrf domain (ii,jj,kk,hh)
                                                                    !! for radar, it is km; for sounding, it's mb;
       real, allocatable, dimension(:,:)         :: sta         !! (ob_num,4)
                                                                !! station attribution. For radar, it includes
                                                                !! radar station's (ii,jj,kk,hh) in wrf domain;
                                                                !! for surface obs, it include (elevation, station pressure, t, q)
       integer, allocatable, dimension(:,:)      :: roi         !! (ob_num,3) : 1=horizontal, 2=vertical, 3=h.(non-Q in radiance)
       character(len=16), allocatable, dimension(:) :: sat      !! Name of the Satellite
       integer, allocatable, dimension(:)        :: ch          !! channel of the satellite
       real, allocatable, dimension(:)           :: efov_aScan    !! effective field of view in the along-scan direction
       real, allocatable, dimension(:)           :: efov_cScan    !! effective field of view in the cross-scan direction
       real, allocatable, dimension(:)           :: scan_angle    !! scan angle
       real, allocatable, dimension(:)           :: zenith_angle  !! zenith (incidence) angle
       real, allocatable, dimension(:)           :: azimuth_angle !! azimuth angle (degrees clockwise from north)
       integer, allocatable, dimension(:)        :: batch_id      !! Assimilation batch id
       integer  :: num_batches                                    !! Number of batches
       integer, allocatable, dimension(:) :: kick_flag

   end type obs_type

   type ( raw_type )   :: raw
   type ( obs_type )   :: obs, master_obs

!---------------------------------------------------------------------------------------------------------------------


  CONTAINS


  ! ==============================================================================
  ! Subroutine to allocate an observation-associated structured type
  ! Man-Yau Chan
  subroutine allocate_obs_type( obs_struct, num_obs )

    implicit none
    type( obs_type ), intent(inout) :: obs_struct
    integer,          intent(in) :: num_obs

    allocate( obs_struct%dat           ( num_obs    ) )
    allocate( obs_struct%type          ( num_obs    ) )
    allocate( obs_struct%err           ( num_obs    ) )
    allocate( obs_struct%position      ( num_obs, 4 ) )
    allocate( obs_struct%sta           ( num_obs, 4 ) )
    allocate( obs_struct%roi           ( num_obs, 3 ) )
    allocate( obs_struct%sat           ( num_obs    ) )
    allocate( obs_struct%ch            ( num_obs    ) )
    allocate( obs_struct%batch_id      ( num_obs    ) )
    allocate( obs_struct%kick_flag     ( num_obs    ) )
    allocate( obs_struct%efov_aScan    ( num_obs    ) )
    allocate( obs_struct%efov_cScan    ( num_obs    ) )
    allocate( obs_struct%scan_angle    ( num_obs    ) )
    allocate( obs_struct%zenith_angle  ( num_obs    ) )
    allocate( obs_struct%azimuth_angle ( num_obs    ) )


    obs_struct%dat           = -888888.
    obs_struct%type          = '          '
    obs_struct%err           = -888888.
    obs_struct%position      = -888888.
    obs_struct%sta           = -888888.
    obs_struct%roi           = -888888
    obs_struct%sat           = '                '
    obs_struct%ch            = -888888
    obs_struct%batch_id      = -888888
    obs_struct%efov_aScan    = -888888
    obs_struct%efov_cScan    = -888888
    obs_struct%scan_angle    = -888888
    obs_struct%zenith_angle  = -888888
    obs_struct%azimuth_angle = -888888
    obs_struct%kick_flag = 0


  end subroutine allocate_obs_type




  ! ==============================================================================
  ! Subroutine to deallocate an observation-associated structured type variable
  ! Man-Yau Chan
  subroutine deallocate_obs_type( obs_struct )


    implicit none
    type( obs_type ), intent( inout ) :: obs_struct

    deallocate( obs_struct%dat           )
    deallocate( obs_struct%type          )
    deallocate( obs_struct%err           )
    deallocate( obs_struct%position      )
    deallocate( obs_struct%sta           )
    deallocate( obs_struct%roi           )
    deallocate( obs_struct%sat           )
    deallocate( obs_struct%ch            )
    deallocate( obs_struct%batch_id      )
    deallocate( obs_struct%kick_flag     )
    deallocate( obs_struct%efov_aScan    )
    deallocate( obs_struct%efov_cScan    )
    deallocate( obs_struct%scan_angle    )
    deallocate( obs_struct%zenith_angle  )
    deallocate( obs_struct%azimuth_angle )

  end subroutine deallocate_obs_type



  ! =========================================================================
  ! Subroutine to copy obs_type variable (copy obsA to obsB)
  ! Man-Yau Chan
  subroutine copy_obs_type_variable( obsA, obsB )

    ! Variable definitions (if any)
    implicit none
    type( obs_type ), intent(inout) :: obsA, obsB
    logical :: flag_filled


    ! 1) Prepare obsB to receive info from obsA.
    ! ------------------------------------------
    ! If obsB stuff has already been allocated, deallocate obsB
    if ( allocated( obsB%dat ) ) call deallocate_obs_type( obsB )

    ! Allocate obsB based on obsA number
    call allocate_obs_type( obsB, obsA%num )


    ! 2) Copy info from obsA to obsB
    ! ------------------------------
    obsB%dat           ( :    ) = obsA%dat           ( :    )
    obsB%type          ( :    ) = obsA%type          ( :    )
    obsB%err           ( :    ) = obsA%err           ( :    )
    obsB%position      ( :, : ) = obsA%position      ( :, : )
    obsB%sta           ( :, : ) = obsA%sta           ( :, : )
    obsB%kick_flag     ( :    ) = obsA%kick_flag     ( :    )
    obsB%roi           ( :, : ) = obsA%roi           ( :, : )
    obsB%sat           ( :    ) = obsA%sat           ( :    )
    obsB%ch            ( :    ) = obsA%ch            ( :    )
    obsB%batch_id      ( :    ) = obsA%batch_id      ( :    )
    obsB%efov_aScan    ( :    ) = obsA%efov_aScan    ( :    )
    obsB%efov_cScan    ( :    ) = obsA%efov_cScan    ( :    )
    obsB%scan_angle    ( :    ) = obsA%scan_angle    ( :    )
    obsB%zenith_angle  ( :    ) = obsA%zenith_angle  ( :    )
    obsB%azimuth_angle ( :    ) = obsA%azimuth_angle ( :    )
    obsB%num                    = obsA%num + 0
    obsB%num_batches            = obsA%num_batches + 0


  end subroutine copy_obs_type_variable



  ! ===============================================================
  ! Subroutine to copy indA-th obs from obsA to indB-th obs in obsB
  ! Man-Yau Chan
  subroutine copy_obs_type_entry( obsA, indA, obsB, indB )


    ! Variable definitions (if any)
    implicit none
    type( obs_type ), intent(inout) :: obsA, obsB
    integer, intent(in) :: indA, indB
    logical :: flag_filled

    ! Copy item indA from obsA to item indB in obsB
    obsB%dat           ( indB    ) = obsA%dat           ( indA   )
    obsB%type          ( indB    ) = obsA%type          ( indA    )
    obsB%err           ( indB    ) = obsA%err           ( indA    )
    obsB%position      ( indB, : ) = obsA%position      ( indA, : )
    obsB%sta           ( indB, : ) = obsA%sta           ( indA, : )
    obsB%roi           ( indB, : ) = obsA%roi           ( indA, : )
    obsB%sat           ( indB    ) = obsA%sat           ( indA    )
    obsB%ch            ( indB    ) = obsA%ch            ( indA    )
    obsB%kick_flag     ( indB    ) = obsA%kick_flag     ( indA    )
    obsB%batch_id      ( indB    ) = obsA%batch_id      ( indA    )
    obsB%efov_aScan    ( indB    ) = obsA%efov_aScan    ( indA    )
    obsB%efov_cScan    ( indB    ) = obsA%efov_cScan    ( indA    )
    obsB%scan_angle    ( indB    ) = obsA%scan_angle    ( indA    )
    obsB%zenith_angle  ( indB    ) = obsA%zenith_angle  ( indA    )
    obsB%azimuth_angle ( indB    ) = obsA%azimuth_angle ( indA    )

  end subroutine copy_obs_type_entry



  
  ! =================================================================
  ! Subroutine to discern all unique observation types
  subroutine seek_unique_obs_types( obsA )

    ! Variable definitions
    implicit none 
    type( obs_type ), intent(inout) :: obsA
    integer :: iob, i_uniq
    logical :: flag_found

    ! Initialize variables
    uniq_obs_types = '          '
    n_uniq_obs_types=0


    ! Iterate thru each observation in list
    read_all_obs_types: do iob = 1, obsA%num


      ! Check uniq_types. 
      ! ----------------
      ! If type has been recorded, indicate with flag_found
      flag_found = .false.
      check_recorded_types: do i_uniq = 1, n_uniq_obs_types + 1

        ! Comparing recorded types
        if ( trim(uniq_obs_types(i_uniq)) == trim(obs%type(iob)) ) then
          flag_found = .true.
          exit check_recorded_types
        endif
     
      enddo check_recorded_types

      
      ! If obs%type is unrecorded, record as a unique value
      ! ---------------------------------------------------
      record_obs_type: if (flag_found == .false.) then
        uniq_obs_types( n_uniq_obs_types +1) = obs%type(iob)
        n_uniq_obs_types = n_uniq_obs_types + 1
      endif record_obs_type


    enddo read_all_obs_types

  end subroutine seek_unique_obs_types





  ! ==================================================================
  ! Subroutine to split obs data into non-overlapping batches
  ! Currently, all non-IR obs are lumped together in a single batch.
  subroutine assign_obs_batch_id( obsA )


    use common_variables

    ! Variable defintions
    implicit none
    type( obs_type ), intent(inout) :: obsA
    integer       :: imax, imin, jmax, jmin, batch_num, seed, zone_rad
    logical       :: flag_sort
    integer       :: iob, iob2, itype

    ! Array to figure out obs overlaps
    integer, dimension(ix+1, jx+1) :: obs_zone_map


    ! First figure out all unique obs types
    call seek_unique_obs_types( obsA )


    ! Allocate obs batch id
    if ( allocated( obsA%batch_id ) ) deallocate( obsA%batch_id )
    allocate( obsA%batch_id ( obsA%num    ) )
    obsA%batch_id = -888888


    ! Preparing to split obs into batches
    flag_sort = .true.
    batch_num = 1

    ! Iterate thru all unique obs types
    loop_sort_obs_types: do itype = 1, n_uniq_obs_types

      flag_sort = .true.

      ! Iterate thru all obs of indicated type
      sort_batches: do while ( flag_sort )

        ! 1) Seeking an unsorted obs and initialize exclusion zone
        ! --------------------------------------------------------
        seek_unsorted: do iob = 1, obsA%num

        ! Exclude non-selected obs type
          if (trim(obsA%type(iob)) .ne. trim(uniq_obs_types(itype)) ) &
            cycle seek_unsorted

          ! Stop when we reach an unsorted obs
          if (obsA%batch_id(iob) < 0) exit seek_unsorted

        end do seek_unsorted

        ! Initialize the exclusion zone
        obs_zone_map = 0

        ! Input observation into the zone
        zone_rad = max( obsA%roi(iob,1), obsA%roi(iob,3) ) + 3
        imax = min( ix+1, ceiling( obsA%position( iob, 1 ) + zone_rad ) )
        imin = max(    1,   floor( obsA%position( iob, 1 ) - zone_rad ) )
        jmax = min( jx+1, ceiling( obsA%position( iob, 2 ) + zone_rad ) )
        jmin = max(    1,   floor( obsA%position( iob, 2 ) - zone_rad ) )
        obs_zone_map( imin:imax, jmin:jmax ) = batch_num
        !write(*,*) 'cleared first zone map'

        ! Filling in the value of the batch for the first obs
        obsA%batch_id( iob ) = batch_num



        ! 2) Loop over all unsorted observations
        ! --------------------------------------
        seek_batch_obs: do iob2 = 1, obsA%num

          ! Exclude non-selected obs type
          if (trim(obsA%type(iob2)) .ne. trim(uniq_obs_types(itype)) ) &
            cycle seek_batch_obs

          ! Exclude sorted obs
          if ( obsA%batch_id(iob2) > 0 ) cycle seek_batch_obs

          ! Exclude obs that falls within the ROI of other obs
          imin = int(obsA%position(iob2,1)); imax = imin+1
          jmin = int(obsA%position(iob2,2)); jmax = jmin+1
          if ( sum(obs_zone_map(imin:imax, jmin:jmax)) > 0 ) &
            cycle seek_batch_obs

          ! If obs survives till here, then the obs should be filed into
          ! the batch
          obsA%batch_id( iob2 ) = batch_num

          ! Expanding the exclusion zone to account for the new obs
          zone_rad = max( obsA%roi(iob2,1), obsA%roi(iob2,3) ) + 3
          imax = min( ceiling( obsA%position(iob2,1) + zone_rad ), ix+1 )
          imin = max(   floor( obsA%position(iob2,1) - zone_rad ),    1 )
          jmax = min( ceiling( obsA%position(iob2,2) + zone_rad ), jx+1 )
          jmin = max(   floor( obsA%position(iob2,2) - zone_rad ),    1 )
          obs_zone_map( imin:imax, jmin:jmax ) = batch_num
          !write(*,*) 'cleared ',iob2,' zone map'

        end do seek_batch_obs



        ! 3) Check for unsorted obs
        ! -------------------------
        flag_sort = .false.
        batch_num = batch_num + 1
        loop_check_unsorted: do iob2 = 1, obsA%num

          ! Skip over non-selected obs type
          if (trim(obsA%type(iob2)) .ne. trim(uniq_obs_types(itype)) ) &
            cycle loop_check_unsorted

          ! If unsorted obs detected, will do another round of batch
          ! sorting
          if ( obsA%batch_id(iob2) < 0 ) then
            flag_sort = .true.
            exit loop_check_unsorted
          endif

        end do loop_check_unsorted

      end do sort_batches
    enddo loop_sort_obs_types

    ! Record the maximum batch num
    obsA%num_batches = batch_num-1

  end subroutine assign_obs_batch_id
  ! ==================================================================






  ! ====================================================================
  ! Special subroutine to rearrange obs in order of ascending batch ids
  subroutine sort_obs_ascending_batch_id( obsA )


    ! Variable definitions
    implicit none
    type( obs_type ), intent(inout) :: obsA
    type( obs_type ) :: holder
    integer :: b_id0, b_id, iob1, iob2


    ! Make copy of the input_obs
    call copy_obs_type_variable( obsA, holder )


    ! Reinitialize obs variable
    call deallocate_obs_type( obsA )
    call allocate_obs_type( obsA, holder%num )
    obsA%num = holder%num


    ! Initialize counting variables
    iob1 = 1; iob2 = 1


    ! Iteratively sort through batch ids
    batch_order_loop: do b_id0 = 0, holder%num_batches

      ! Special handling for obs not handled by the batch assignment
      if (b_id0 == 0) then
        b_id = -888888
      else
        b_id = b_id0
      endif

      ! Iterate thru the observations to seek obs with the matching
      ! batch_id
      inner_rearrange_loop: do iob1 = 1, holder%num

        ! If batch id matches, copy in
        if ( holder%batch_id(iob1) == b_id ) then

          ! Copy matching batch item
          call copy_obs_type_entry( holder, iob1, obsA, iob2 )

          ! Increment counter
          iob2 = iob2 +1

        endif

      enddo inner_rearrange_loop

    enddo batch_order_loop


    ! Sanity check: batch_ids should be in ascending order now
    !write(*,*) obsA%batch_id
    loop_check_b_id_order: do iob1 = 2, holder%num

      ! Failure condition: batch_id decreases
      if ( obsA%batch_id(iob1-1)  .gt. obsA%batch_id(iob1) ) then
        write(*,'(a,i8,a,i8,a,i8,a,i8)') &
           'batch_id ',iob1-1, ' = ', &
           obsA%batch_id(iob1-1), ', but batch_id ',iob1, ' = ', &
           obsA%batch_id(iob1)
        write(*,*) ' batch ids must be ascending! Exiting.'
        call exit
      endif

    enddo loop_check_b_id_order


  end subroutine sort_obs_ascending_batch_id
  ! ===========================================================================





end module obs_define
!==============================================================================
