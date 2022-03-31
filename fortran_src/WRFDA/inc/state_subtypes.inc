!STARTOFREGISTRYGENERATEDINCLUDE 'inc/state_subtypes.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
TYPE vp_type
real      ,DIMENSION(:,:,:)   ,POINTER   :: v1
real      ,DIMENSION(:,:,:)   ,POINTER   :: v2
real      ,DIMENSION(:,:,:)   ,POINTER   :: v3
real      ,DIMENSION(:,:,:)   ,POINTER   :: v4
real      ,DIMENSION(:,:,:)   ,POINTER   :: v5
real      ,DIMENSION(:,:,:)   ,POINTER   :: v6
real      ,DIMENSION(:,:,:)   ,POINTER   :: v7
real      ,DIMENSION(:,:,:)   ,POINTER   :: v8
real      ,DIMENSION(:,:,:)   ,POINTER   :: v9
real      ,DIMENSION(:,:,:)   ,POINTER   :: v10
real      ,DIMENSION(:,:,:)   ,POINTER   :: v11
real      ,DIMENSION(:,:,:,:) ,POINTER   :: alpha
END TYPE vp_type
TYPE ep_type
integer                                  :: ne
real      ,DIMENSION(:,:,:,:) ,POINTER   :: v1
real      ,DIMENSION(:,:,:,:) ,POINTER   :: v2
real      ,DIMENSION(:,:,:,:) ,POINTER   :: v3
real      ,DIMENSION(:,:,:,:) ,POINTER   :: v4
real      ,DIMENSION(:,:,:,:) ,POINTER   :: v5
END TYPE ep_type
TYPE xb_type
integer                                  :: map
real                                     :: ptop
real                                     :: ztop
real                                     :: ds
real                                     :: sigmaf_0
integer                                  :: mix
integer                                  :: mjy
integer                                  :: mkz
integer                                  :: num_of_var
integer                                  :: ids
integer                                  :: ide
integer                                  :: jds
integer                                  :: jde
integer                                  :: kds
integer                                  :: kde
integer                                  :: ims
integer                                  :: ime
integer                                  :: jms
integer                                  :: jme
integer                                  :: kms
integer                                  :: kme
integer                                  :: its
integer                                  :: ite
integer                                  :: jts
integer                                  :: jte
integer                                  :: kts
integer                                  :: kte
real      ,DIMENSION(:,:)     ,POINTER   :: grid_box_area
real      ,DIMENSION(:,:)     ,POINTER   :: map_factor
real      ,DIMENSION(:,:,:)   ,POINTER   :: u
real      ,DIMENSION(:,:,:)   ,POINTER   :: v
real      ,DIMENSION(:,:,:)   ,POINTER   :: w
real      ,DIMENSION(:,:,:)   ,POINTER   :: t
real      ,DIMENSION(:,:,:)   ,POINTER   :: p
real      ,DIMENSION(:,:,:)   ,POINTER   :: q
real      ,DIMENSION(:,:,:)   ,POINTER   :: qrn
real      ,DIMENSION(:,:,:)   ,POINTER   :: qcw
real      ,DIMENSION(:,:,:)   ,POINTER   :: qci
real      ,DIMENSION(:,:,:)   ,POINTER   :: qsn
real      ,DIMENSION(:,:,:)   ,POINTER   :: qgr
real      ,DIMENSION(:,:,:)   ,POINTER   :: qt
real      ,DIMENSION(:,:,:)   ,POINTER   :: delt
real      ,DIMENSION(:,:,:)   ,POINTER   :: xb_p_x
real      ,DIMENSION(:,:,:)   ,POINTER   :: xb_p_y
real      ,DIMENSION(:,:,:)   ,POINTER   :: h
real      ,DIMENSION(:,:,:)   ,POINTER   :: hf
real      ,DIMENSION(:,:,:)   ,POINTER   :: wh
real      ,DIMENSION(:,:,:)   ,POINTER   :: rho
real      ,DIMENSION(:,:,:)   ,POINTER   :: rh
real      ,DIMENSION(:,:,:)   ,POINTER   :: es
real      ,DIMENSION(:,:,:)   ,POINTER   :: qs
real      ,DIMENSION(:,:,:)   ,POINTER   :: td
real      ,DIMENSION(:,:,:)   ,POINTER   :: vertical_inner_product
real      ,DIMENSION(:,:,:)   ,POINTER   :: ref
real      ,DIMENSION(:,:,:)   ,POINTER   :: reflog
real      ,DIMENSION(:,:)     ,POINTER   :: rough
real      ,DIMENSION(:,:)     ,POINTER   :: ztd
real      ,DIMENSION(:,:)     ,POINTER   :: tpw
real      ,DIMENSION(:,:)     ,POINTER   :: slp
real      ,DIMENSION(:,:)     ,POINTER   :: speed
real      ,DIMENSION(:,:)     ,POINTER   :: tb19v
real      ,DIMENSION(:,:)     ,POINTER   :: tb19h
real      ,DIMENSION(:,:)     ,POINTER   :: tb22v
real      ,DIMENSION(:,:)     ,POINTER   :: tb37v
real      ,DIMENSION(:,:)     ,POINTER   :: tb37h
real      ,DIMENSION(:,:)     ,POINTER   :: tb85v
real      ,DIMENSION(:,:)     ,POINTER   :: tb85h
real      ,DIMENSION(:,:)     ,POINTER   :: cori
real      ,DIMENSION(:,:)     ,POINTER   :: tgrn
real      ,DIMENSION(:,:)     ,POINTER   :: lat
real      ,DIMENSION(:,:)     ,POINTER   :: lon
real      ,DIMENSION(:,:)     ,POINTER   :: terr
real      ,DIMENSION(:,:)     ,POINTER   :: psac
real      ,DIMENSION(:,:)     ,POINTER   :: lanu
real      ,DIMENSION(:,:)     ,POINTER   :: snow
real      ,DIMENSION(:,:)     ,POINTER   :: xland
real      ,DIMENSION(:,:)     ,POINTER   :: landmask
real      ,DIMENSION(:,:)     ,POINTER   :: u10
real      ,DIMENSION(:,:)     ,POINTER   :: v10
real      ,DIMENSION(:,:)     ,POINTER   :: t2
real      ,DIMENSION(:,:)     ,POINTER   :: q2
real      ,DIMENSION(:,:)     ,POINTER   :: psfc
real      ,DIMENSION(:,:)     ,POINTER   :: tsk
real      ,DIMENSION(:,:)     ,POINTER   :: smois
real      ,DIMENSION(:,:)     ,POINTER   :: tslb
real      ,DIMENSION(:,:)     ,POINTER   :: xice
integer   ,DIMENSION(:,:)     ,POINTER   :: ivgtyp
integer   ,DIMENSION(:,:)     ,POINTER   :: isltyp
real      ,DIMENSION(:,:)     ,POINTER   :: vegfra
real      ,DIMENSION(:,:)     ,POINTER   :: snowh
real      ,DIMENSION(:,:)     ,POINTER   :: regime
real      ,DIMENSION(:)       ,POINTER   :: znu
real      ,DIMENSION(:)       ,POINTER   :: znw
real      ,DIMENSION(:)       ,POINTER   :: dn
real      ,DIMENSION(:)       ,POINTER   :: dnw
real      ,DIMENSION(:)       ,POINTER   :: sigmah
real      ,DIMENSION(:)       ,POINTER   :: sigmaf
real      ,DIMENSION(:,:)     ,POINTER   :: coefx
real      ,DIMENSION(:,:)     ,POINTER   :: coefy
real      ,DIMENSION(:,:)     ,POINTER   :: coefz
real      ,DIMENSION(:)       ,POINTER   :: kma_a
real      ,DIMENSION(:)       ,POINTER   :: kma_b
END TYPE xb_type
TYPE x_type
real      ,DIMENSION(:,:,:)   ,POINTER   :: u
real      ,DIMENSION(:,:,:)   ,POINTER   :: v
real      ,DIMENSION(:,:,:)   ,POINTER   :: w
real      ,DIMENSION(:,:,:)   ,POINTER   :: t
real      ,DIMENSION(:,:,:)   ,POINTER   :: q
real      ,DIMENSION(:,:,:)   ,POINTER   :: p
real      ,DIMENSION(:,:,:)   ,POINTER   :: geoh
real      ,DIMENSION(:,:,:)   ,POINTER   :: rh
real      ,DIMENSION(:,:,:)   ,POINTER   :: rho
real      ,DIMENSION(:,:,:)   ,POINTER   :: wh
real      ,DIMENSION(:,:,:)   ,POINTER   :: grad_p_x
real      ,DIMENSION(:,:,:)   ,POINTER   :: grad_p_y
real      ,DIMENSION(:,:,:)   ,POINTER   :: qrn
real      ,DIMENSION(:,:,:)   ,POINTER   :: qcw
real      ,DIMENSION(:,:,:)   ,POINTER   :: qci
real      ,DIMENSION(:,:,:)   ,POINTER   :: qsn
real      ,DIMENSION(:,:,:)   ,POINTER   :: qgr
real      ,DIMENSION(:,:,:)   ,POINTER   :: qt
real      ,DIMENSION(:,:)     ,POINTER   :: tgrn
real      ,DIMENSION(:,:)     ,POINTER   :: psfc
real      ,DIMENSION(:,:)     ,POINTER   :: mu
real      ,DIMENSION(:,:)     ,POINTER   :: u10
real      ,DIMENSION(:,:)     ,POINTER   :: v10
real      ,DIMENSION(:,:)     ,POINTER   :: t2
real      ,DIMENSION(:,:)     ,POINTER   :: q2
real      ,DIMENSION(:,:,:)   ,POINTER   :: ref
real      ,DIMENSION(:,:)     ,POINTER   :: ztd
real      ,DIMENSION(:,:)     ,POINTER   :: tpw
real      ,DIMENSION(:,:)     ,POINTER   :: speed
real      ,DIMENSION(:,:)     ,POINTER   :: tb19v
real      ,DIMENSION(:,:)     ,POINTER   :: tb19h
real      ,DIMENSION(:,:)     ,POINTER   :: tb22v
real      ,DIMENSION(:,:)     ,POINTER   :: tb37v
real      ,DIMENSION(:,:)     ,POINTER   :: tb37h
real      ,DIMENSION(:,:)     ,POINTER   :: tb85v
real      ,DIMENSION(:,:)     ,POINTER   :: tb85h
END TYPE x_type
TYPE xpose_type
integer                                  :: domdesc
integer                                  :: comms(max_comms)
integer                                  :: ids
integer                                  :: ide
integer                                  :: jds
integer                                  :: jde
integer                                  :: kds
integer                                  :: kde
integer                                  :: ims
integer                                  :: ime
integer                                  :: jms
integer                                  :: jme
integer                                  :: kms
integer                                  :: kme
integer                                  :: its
integer                                  :: ite
integer                                  :: jts
integer                                  :: jte
integer                                  :: kts
integer                                  :: kte
integer                                  :: ips
integer                                  :: ipe
integer                                  :: jps
integer                                  :: jpe
integer                                  :: kps
integer                                  :: kpe
integer                                  :: idsx
integer                                  :: idex
integer                                  :: jdsx
integer                                  :: jdex
integer                                  :: kdsx
integer                                  :: kdex
integer                                  :: imsx
integer                                  :: imex
integer                                  :: jmsx
integer                                  :: jmex
integer                                  :: kmsx
integer                                  :: kmex
integer                                  :: itsx
integer                                  :: itex
integer                                  :: jtsx
integer                                  :: jtex
integer                                  :: ktsx
integer                                  :: ktex
integer                                  :: ipsx
integer                                  :: ipex
integer                                  :: jpsx
integer                                  :: jpex
integer                                  :: kpsx
integer                                  :: kpex
integer                                  :: idsy
integer                                  :: idey
integer                                  :: jdsy
integer                                  :: jdey
integer                                  :: kdsy
integer                                  :: kdey
integer                                  :: imsy
integer                                  :: imey
integer                                  :: jmsy
integer                                  :: jmey
integer                                  :: kmsy
integer                                  :: kmey
integer                                  :: itsy
integer                                  :: itey
integer                                  :: jtsy
integer                                  :: jtey
integer                                  :: ktsy
integer                                  :: ktey
integer                                  :: ipsy
integer                                  :: ipey
integer                                  :: jpsy
integer                                  :: jpey
integer                                  :: kpsy
integer                                  :: kpey
real      ,DIMENSION(:,:)     ,POINTER   :: vxy
real      ,DIMENSION(:,:,:)   ,POINTER   :: v1z
real      ,DIMENSION(:,:,:)   ,POINTER   :: v1x
real      ,DIMENSION(:,:,:)   ,POINTER   :: v1y
real      ,DIMENSION(:,:,:)   ,POINTER   :: v2z
real      ,DIMENSION(:,:,:)   ,POINTER   :: v2x
real      ,DIMENSION(:,:,:)   ,POINTER   :: v2y
END TYPE xpose_type
!ENDOFREGISTRYGENERATEDINCLUDE
